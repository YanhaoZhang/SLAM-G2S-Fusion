classdef NLSFormulation < PoseGraph

    %
    
    properties
        
        % variable sequence given by node IDs
        NodesIDSeq = []

        % measurement sequence given by edge IDs
        EdgesIDSeq = []
        
        % maximum iterations
        MaxIters = 10; 
        
        % consumed iteration by Gauss-Newton
        UsedIters = uint32(0)
        
        % number of loop-closures
        NumLpCls = uint32(0) 
        
        % flag for debug mode
        debug = true
        
        % Fixed nodes array
        FixNodes = [ 1 ]
        
        % condition number of information matrix
        InfoCond = 0
        
        % use isotropic noise
        IsotropicNoise = 0 % 0 - disabled,  otherwise, enabled % set coefficient 1 - 1000
        
       
    end
    
    properties(Access = public)
        
        % A brief cache of edge information
        % edge and variable indices and measurement data
        EdgeInfoCache = struct( ...
            'mvIndex', [], ...
            'mdata', [] )
        
        % structure / Laplacian matrix of the graph
        P = struct( ...
            'i', [], ...
            'j', [], ...
            'v', [] )        
        Q = struct( ...
            'i', [], ...
            'j', [], ...
            'v', [] )    
        
        
        % matrices for linearization
        A = struct( ...
            'i', [], ...
            'j', [], ...
            'v', [] )
       
         
        
        b = []
        
        % inverse-covariance matrix
        CovMatrInv = struct( ...
            'i', [], ...
            'j', [], ...
            'v', [] )
        
        % preserve loop-closure structure to speed up search
        LpCls = struct( ...
            'LpClEdgeID', [], ...
            'LpClMetric', [], ...
            'Enabled', [] )
        
        % mapping LpCls ID to its storage idx
        LpClsID2Idx = []
        
        % mapping node ID to its variable index
        NodesID2VarIdx = []
        
        % mapping edge ID to its measurement index 
        EdgesID2MstIdx = []
        
    end
    
    
    properties(Access = public)
        
        %
        p_b = []  
        p_A= []
        p_x = []
        
    end
    
    
    methods (Access = public)
        
        % constructor
        function this = NLSFormulation ()
            % call superclass constructor
            this@PoseGraph();
            % add field 'Error' for each edge
            % Error = Log(inv(measure) * inv(pose1) * pose2)
            [ this.Edges.Error ] = deal({});
            if (nargin == 0)
                % initialization
            end
        end
        
        % overload addNode
        function addNode(this, id, pose)
            % call superclass addNode
            this.addNode@PoseGraph(id, pose);
            % Initial Node ID sequence is given by the node adding order
            this.NodesIDSeq = [ this.NodesIDSeq, id ];
        end
        
        % overload addEdge
        function addEdge(this, id, sid, tid, rlvpose, covinv) 
            % call superclass addEdge
            this.addEdge@PoseGraph(id, sid, tid, rlvpose, covinv);
            % Initial Edge ID sequence is given by the edge adding order
            this.EdgesIDSeq = [ this.EdgesIDSeq, id ];
            % add edge ID to LpClEdgeID if this edge is not odometry edge
            if (tid - sid > 1 || sid - tid > 1) 
                this.NumLpCls = this.NumLpCls + 1;                
                % put edge id into LpClEdgeID
                this.LpCls.LpClEdgeID = [ this.LpCls.LpClEdgeID, id ];
                this.LpCls.LpClMetric = [ this.LpCls.LpClMetric, inf ];
                this.LpCls.Enabled = [ this.LpCls.Enabled, 1 ];
                % mapping LpCls ID to its storage index (idx)
                this.LpClsID2Idx(id) = this.NumLpCls;
            end
        end
        
        % check if an added edge is enabled(1) or not(0)
        function flag = isEnabled (this, id)
            idx = this.EdgeID2Idx(id);
            sid = this.Edges(idx).snID;
            tid = this.Edges(idx).tnID;
            if( tid - sid == 1 || sid - tid == 1 ) 
                flag = 1; % always enable odometry edge
            else
                flag = this.LpCls.Enabled(this.LpClsID2Idx(id));
            end
        end
        
        % set loop-closure edge falg = enabled(1)/disabled(0)
        function setEnabled (this, id, flag)
            idx = this.EdgeID2Idx(id);
            sid = this.Edges(idx).snID;
            tid = this.Edges(idx).tnID;
            if( tid - sid == 1 || sid - tid == 1 ) 
                error('do not enabled/disable odometry edges!');
            end
            this.LpCls.Enabled(this.LpClsID2Idx(id)) = flag;
        end
        
        %-+-%
        function objFunc = computeObjFunc (this)
            objFunc = 0;
            for idx = 1 : this.NumEdges
                id = this.Edges(idx).EdgeID;
                if(this.EdgesID2MstIdx(id)>0) % enabled
                    vecError = this.Edges(idx).Error;
                    if(this.IsotropicNoise > 0)
                        objFunc = objFunc + ...
                            vecError' * ( this.IsotropicNoise * eye(6) ) * vecError;
                    else
                        objFunc = objFunc + ...
                            vecError' * this.Edges(idx).CovInv * vecError;
                    end
                end
            end
        end
       
        %-+-%
        function RlvErrorVec = obtainRelativeError (this)
            RlvErrorVec = [];
            for idx = 1 : this.NumEdges
                id = this.Edges(idx).EdgeID;
                if(this.EdgesID2MstIdx(id)>0) % enabled
                    RlvErrorVec = [ RlvErrorVec; this.Edges(idx).Error;];
                end
            end
        end
        
        
        
        %
        function solveByGaussNewton (this)
            this.initializeProblem();
            this.obtainCovMatrInv();
            %
            % this.createStructureLaplacianMatrix ();
            % Pmatr = sparse(this.P.i, this.P.j, this.P.v);
            % Qmatr = sparse(this.Q.i, this.Q.j, this.Q.v);
            % 
            covinv = sparse(this.CovMatrInv.i, ...
                                            this.CovMatrInv.j, ...
                                            this.CovMatrInv.v);
            %
            this.computeEdgeError();
            % 

            %
            iter = 1;
            
            %---- condition number of information matrix
%             this.InfoCond = cond(covinv);
%             disp (['Condition Number of Information Matrix: ', num2str(this.InfoCond)]);

            while( iter < this.MaxIters )
                
                
                % ++++++++ Stage 1: Linearization Process ++++++++ %
                %-- used to obtain new linearized matrix [ A   b ]
                this.linearizeProblem();
                Amatr = sparse(this.A.i, this.A.j, this.A.v);
                bvec = this.b;
                
                
                % ++++++++ Stage 2: Linear Solve Process ++++++++ %
                %-- solve normal equation to get x_star
                x_star_vec = this.solveNormalEquation (Amatr, bvec, covinv);
     
                
                % ++++++++ Stage 3: Update Process ++++++++ %
                %-- update pose estimate
                this.updatePoseEstimate (x_star_vec);
                %-- update edge error, namely linearization point
                this.computeEdgeError();
                
                
                %## Optimized linearization point ##%
                LinearizationPoint_Optimized = this.obtainRelativeError();
                LinearizationPoint_LinearUpdate = bvec - Amatr * x_star_vec;
                
                
                %## Record Optimized Linear and Nonlinear Objective Function ##% 
                objFuncLinear = (bvec - Amatr * x_star_vec)' * covinv * (bvec - Amatr * x_star_vec);
                
                objFunc_Optimized = LinearizationPoint_Optimized' * covinv * LinearizationPoint_Optimized;
                
                objFunc_LinearUpdate = LinearizationPoint_LinearUpdate' * covinv * LinearizationPoint_LinearUpdate;

                Error_LinearUpdate = LinearizationPoint_Optimized - LinearizationPoint_LinearUpdate;
               
                
                fprintf(1, [ '\n#  Iter ', num2str(iter), ' Optimized  #', '\n\n']);
                
                
                 fprintf(2, [ '     objFuncs:     nonlinear_Optimized = ', num2str(objFunc_Optimized),...
                    '     nonlinear_LinearUpdated = ', num2str(objFunc_LinearUpdate), ...
                    '     linear = ', num2str(objFuncLinear), ...
                    '\n\n']);               
                
                fprintf(1, ['      x_star_norm:     ', num2str( sqrt(x_star_vec' * x_star_vec) ), ...
                               '      x_star_avg:     ', num2str( sqrt(x_star_vec' * x_star_vec)/this.NumNodes ), ...
                               '      x_star_max:     ', num2str( max(abs(x_star_vec)) ), ...
                               '\n\n' ]);
                
                fprintf(1, ['      b_norm:     ', num2str( sqrt(bvec' * bvec) ), ...
                               '      b_avg:     ', num2str( sqrt(bvec' * bvec)/this.NumEdges ), ...
                               '      b_max:     ', num2str( max(abs(bvec)) ), ...
                               '\n\n' ]);
                           
                 fprintf(1, ['      Error_LinearUpdate_norm:     ', num2str( sqrt(Error_LinearUpdate' * Error_LinearUpdate) ), ...
                                '      Error_LinearUpdate_avg:     ', num2str( sqrt(Error_LinearUpdate' * Error_LinearUpdate)/this.NumEdges ), ...
                                '      Error_LinearUpdate_max:     ', num2str( max(abs(Error_LinearUpdate)) ), ...
                               '\n\n' ]);                          
                           



                % save the linear system for the next iteration
                this.p_A = Amatr;
                this.p_b = bvec;
                this.p_x = x_star_vec;
                %---------------------------------------%
                if(norm(x_star_vec) < 1e-6)
                    this.UsedIters = iter;
                    break
                end
                iter = iter + 1;
            end
        end
        
        %
        function changeVariableOrdering (this)
            if (nargin > 1)
                % for example, sort in ascending order
                this.NodesIDSeq = sort(this.NodesIDSeq, 'ascend');
            end
        end
        
        %
        function changeMeasurementOrdering (this)
            if (nargin > 1)
                % for example, sort in ascending order
                this.EdgesIDSeq = sort(this.EdgesIDSeq, 'ascend');
            end
        end
        
        %
        function initializeProblem (this)
            this.changeVariableOrdering();
            this.changeMeasurementOrdering();
            this.initializeVariableIndex();
            this.initializeMeasurementIndex();
            this.initializeVectorMemory();
        end
        
        % add g2o data into the graph
        add_data_g2o (this, FilePathName);
        % plot trajectory of odometry
        [ h1, h2 ] = plotTrajectory(this, fid);
        % dump 3d data
        DumpData3D (this, filename);
    end
    

    % private methods
    methods (Access = protected)

        %
        function initializeVariableIndex (this)
            vidx = 1;
            for idx = 1 : this.NumNodes
                id = this.NodesIDSeq(idx);
                if(~ismember(id, this.FixNodes))
                    this.NodesID2VarIdx(id) = vidx;
                    vidx = vidx + 1;
                else
                    % default value for fixNodes
                    this.NodesID2VarIdx(id) = 0;
                end                
            end
            this.NodesID2VarIdx = 6 * this.NodesID2VarIdx - 5;
        end
        
        %
        function initializeMeasurementIndex (this)
            midx = 1;
            for idx = 1 : this.NumEdges
                id = this.EdgesIDSeq(idx);
                if(this.isEnabled(id) > 0) % enabled
                    this.EdgesID2MstIdx(id) = midx;
                    midx = midx + 1;
                else
                    % default value for disabled edges
                    this.EdgesID2MstIdx(id) = 0;
                end
            end
            this.EdgesID2MstIdx = 6 * this.EdgesID2MstIdx - 5;  
        end
        
        %
        function initializeVectorMemory (this)
            n_edges = 0;
            for idx = 1 : this.NumEdges
                id = this.Edges(idx).EdgeID;
                if(this.EdgesID2MstIdx(id)>0) % enabled
                    n_edges = n_edges + 1;
                end
            end
            % i,j,v vectors for inverse-covariance matrix;
            this.CovMatrInv.i = ones(1, n_edges*36);
            this.CovMatrInv.j = ones(1, n_edges*36);
            this.CovMatrInv.v = zeros(1, n_edges*36);
            % i,j,v vector for objective function
            this.A.i = ones(1, n_edges*36*2);
            this.A.j = ones(1, n_edges*36*2);
            this.A.v = zeros(1, n_edges*36*2);
            this.b = zeros(n_edges*6, 1);
        end
        
                
    end
    
    
    % private methods
    methods(Access = protected)
        
        %
        function obtainCovMatrInv( this )
            ptcov = 1;
            for idx = 1 : this.NumEdges
                id = this.Edges(idx).EdgeID;
                mstIdx = this.EdgesID2MstIdx(id);
                if(mstIdx > 0) % enabled
                    if(this.IsotropicNoise > 0)
                        this.CovMatrInv.v(ptcov:ptcov+35) = ...
                            reshape( this.IsotropicNoise * eye(6), 1, 36);
                    else
                        this.CovMatrInv.v(ptcov:ptcov+35) = ...
                            reshape(this.Edges(idx).CovInv, 1, 36);
                    end
                    this.CovMatrInv.i(ptcov:ptcov+35) = ...
                        kron(ones(1,6), mstIdx:mstIdx+5);
                    this.CovMatrInv.j(ptcov:ptcov+35) = ...
                        kron(mstIdx:mstIdx+5, ones(1,6));              
                    ptcov = ptcov + 36;
                end
            end  
        end
        
        %
        function linearizeProblem (this)
           pta = 1;
           for idx = 1 : this.NumEdges
               % extract edge information
               id = this.Edges(idx).EdgeID;
               mstIdx = this.EdgesID2MstIdx(id);
               if(mstIdx > 0) % enabled
                   rlvpose = this.Edges(idx).RlvPose;
                   error = this.Edges(idx).Error;
                   % compute linearization Jacobian
                   J_sn = SE3.InvJl(error) * SE3.Adj(SE3.Inv(rlvpose));
                   J_tn = -SE3.InvJr(error);
                   % variable indices
                   svarIdx = this.NodesID2VarIdx(this.Edges(idx).snID);
                   tvarIdx = this.NodesID2VarIdx(this.Edges(idx).tnID);
                   % put Jacobian into A matrix
                   if(svarIdx>0) % if not a fixed node
                       this.A.v(pta:pta+35) = reshape(J_sn, 1, 36);
                       this.A.i(pta:pta+35) = kron(ones(1,6), mstIdx:mstIdx+5);
                       this.A.j(pta:pta+35) = kron(svarIdx:svarIdx+5, ones(1,6)); 
                       pta = pta + 36;
                   end
                   if(tvarIdx>0) % if not a fixed node
                       this.A.v(pta:pta+35) = reshape(J_tn, 1, 36);
                       this.A.i(pta:pta+35) = kron(ones(1,6), mstIdx:mstIdx+5);
                       this.A.j(pta:pta+35) = kron(tvarIdx:tvarIdx+5, ones(1,6)); 
                       pta = pta + 36;
                   end
                   this.b(mstIdx:mstIdx+5, 1) = error;
               end
           end
        end
        
        
        % obtain linearization point in terms of error
        %
        function eta = obtainLinearizationPoint (this)
           for idx = 1 : this.NumEdges
               % extract edge information
               id = this.Edges(idx).EdgeID;
               mstIdx = this.EdgesID2MstIdx(id);
               if(mstIdx > 0) % enabled
                   error = this.Edges(idx).Error;
                   eta(mstIdx:mstIdx+5, 1) = error;
               end
           end
        end
        
        
        % create structure / Laplacian matrix of the graph
        function createStructureLaplacianMatrix (this)
           pta = 1;
           for idx = 1 : this.NumEdges
               % extract edge information
               id = this.Edges(idx).EdgeID;
               mstIdx = this.EdgesID2MstIdx(id);
               if(mstIdx > 0) % enabled
                   rlvpose = this.Edges(idx).RlvPose;
                   % compute linearization Jacobian
                   J_sn =  SE3.Adj(SE3.Inv(rlvpose));
                   J_tn =  eye(6);
                   % variable indices
                   svarIdx = this.NodesID2VarIdx(this.Edges(idx).snID);
                   tvarIdx = this.NodesID2VarIdx(this.Edges(idx).tnID);
                   % A brief cache of the edge information
                   this.EdgeInfoCache(idx).mvIndex = [svarIdx, tvarIdx];
                   this.EdgeInfoCache(idx).mdata = rlvpose;
                   % construct constant structure matrix
                   %--- P matrix  [.... Ad(T^{-1}) .... eye(6) ....]
                   %--- Q matrix  [.... Ad(T^{-1}) .... - eye(6) ....]
                   if(svarIdx>0) % if not a fixed node
                       this.P.v(pta:pta+35) = reshape(J_sn, 1, 36);
                       this.P.i(pta:pta+35) = kron(ones(1,6), mstIdx:mstIdx+5);
                       this.P.j(pta:pta+35) = kron(svarIdx:svarIdx+5, ones(1,6));
                       this.Q.v(pta:pta+35) = reshape(J_sn, 1, 36);
                       this.Q.i(pta:pta+35) = kron(ones(1,6), mstIdx:mstIdx+5);
                       this.Q.j(pta:pta+35) = kron(svarIdx:svarIdx+5, ones(1,6)); 
                       pta = pta + 36;
                   end
                   if(tvarIdx>0) % if not a fixed node
                       this.P.v(pta:pta+35) = reshape(J_tn, 1, 36);
                       this.P.i(pta:pta+35) = kron(ones(1,6), mstIdx:mstIdx+5);
                       this.P.j(pta:pta+35) = kron(tvarIdx:tvarIdx+5, ones(1,6)); 
                       this.Q.v(pta:pta+35) = reshape( - J_tn, 1, 36);
                       this.Q.i(pta:pta+35) = kron(ones(1,6), mstIdx:mstIdx+5);
                       this.Q.j(pta:pta+35) = kron(tvarIdx:tvarIdx+5, ones(1,6)); 
                       pta = pta + 36;
                   end
               end
           end
        end
        
        
        %
        function updatePoseEstimate (this, updatevec)
            for idx = 1 : this.NumNodes
                id = this.Nodes(idx).NodeID;
                varIdx = this.NodesID2VarIdx(id);
                if(varIdx > 0) % if not a fixed node
                    vec = updatevec(varIdx:varIdx+5);
                    this.Nodes(idx).Pose = ...
                        this.Nodes(idx).Pose * SE3.Exp(vec);
                end
            end
        end
        
 
        %
        function computeEdgeError (this)
            for idx = 1 : this.NumEdges
                % uncoment if compute Error for enabled edges only
                % id = obj.Edges(idx).EdgeID;
                % if(obj.EdgesID2MstIdx(id) > 0) % enabled
                sid = this.Edges(idx).snID;
                tid = this.Edges(idx).tnID;
                rlvpose = this.Edges(idx).RlvPose;
                this.Edges(idx).Error = SE3.Log( SE3.Inv(rlvpose) * ...
                    SE3.Inv(this.Node(sid)) * this.Node(tid) );
                % end
            end
        end
       
        
        % ----- used to obtain poses for debug
        function Matr = obtainPoseMatrix (this)
            for idx = 1 : this.NumNodes
                id = this.Nodes(idx).NodeID;
                varIdx = this.NodesID2VarIdx(id);
                if(varIdx > 0) % if not a fixed node
                    Matr{varIdx} = this.Nodes(idx).Pose;
                end
            end
        end      
        
       
    end
    
    % private and static methods
    methods(Access = private, Static = true)    
        
        %
        function X = solveNormalEquation (cA, cb, cPinv)
            X = (cA'*cPinv*cA)\(cA'*cPinv*cb);
        end
        
        %
        function matr = adHat (vec)
            matr = [ SO3.Hat(vec (4:6)),  SO3.Hat(vec (1:3));
                          zeros(3,3),               SO3.Hat(vec (4:6)) ];
        end

        %
        function matr = adtHat (vec)
            matr = [ zeros(3,3),                SO3.Hat(vec (1:3));
                          SO3.Hat(vec (1:3)),   SO3.Hat(vec (4:6)) ];
        end            
        
        
        
        %
        function matr = blkdiag_ad_2 (px, qx, xi)
            dim = size(px, 1);
            pt = 1;
            for i = 1 : 6 : dim
                % term1
                term1 = NLSFormulation.adHat (px(i : i+5, 1));
                qx_k = qx(i : i+5, 1);
                ad_xi_k = NLSFormulation.adHat( xi(i : i+5, 1) );
                % term2
                term2 = - NLSFormulation.adHat ( ad_xi_k  * qx_k ) / 6;
                % term3
                term3 = - ad_xi_k * NLSFormulation.adHat(qx_k) / 6;
                % summation
                tmp = term1 + term2 + term3;
                vv(pt:pt+35) = reshape(tmp, 1, 36);
                vi(pt:pt+35) = kron(ones(1,6), i:i+5);
                vj(pt:pt+35) = kron(i:i+5, ones(1,6));
                pt = pt + 36; 
            end
            matr = sparse(vi, vj, vv);   
        end
        
        
        %
        function matr = blkdiag_ad (vec)
            dim = size(vec, 1);
            pt = 1;
            for i = 1 : 6 : dim
                tmp = NLSFormulation.adHat (vec(i : i+5, 1));
                vv(pt:pt+35) = reshape(tmp, 1, 36);
                vi(pt:pt+35) = kron(ones(1,6), i:i+5);
                vj(pt:pt+35) = kron(i:i+5, ones(1,6));
                pt = pt + 36; 
            end
            matr = sparse(vi, vj, vv);
        end
        
        
        %
        function matr = blkdiag_adt (vec)
            dim = size(vec, 1);
            pt = 1;
            for i = 1 : 6 : dim
                tmp = NLSFormulation.adtHat (vec(i : i+5, 1));
                vv(pt:pt+35) = reshape(tmp, 1, 36);
                vi(pt:pt+35) = kron(ones(1,6), i:i+5);
                vj(pt:pt+35) = kron(i:i+5, ones(1,6));
                pt = pt + 36; 
            end
            matr = sparse(vi, vj, vv);
        end
        
        
        % vec2 - vec1
        function diff = computeSE3Diff (vec1, vec2)
            for i =  1 : 6 : size(vec1, 1)
                diff(i:i+5,1) = SE3.Log(SE3.Inv( SE3.Exp(vec1(i : i+5, 1)) ) * SE3.Exp(vec2(i : i+5, 1)) );
            end
        end
        
        
        % row wise check
        function diff = FrobeniusNormDiffRowBasedCheck(Matr1, Matr2, stepSize)
            if( size(Matr1,1) ~= size(Matr2,1) )
                error('Matrix row dimensional inconsistent:  FrobeniusNormDiffRowBasedCheck@NLSFormulation.');
            end
            diff = 0;
            for idx = 1 : stepSize : size(Matr1, 1)
                n1 = norm(Matr1(idx : idx+stepSize-1, :), 'fro');
                n2 = norm(Matr2(idx : idx+stepSize-1, :), 'fro');
                diff = diff + (n1 - n2) * (n1 - n2);
            end
            diff = sqrt(diff);
        end
        
        
    end
    
    
end


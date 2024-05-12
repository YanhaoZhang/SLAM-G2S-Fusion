classdef ConstrainedFormulation < PoseGraph

    %
    
    properties (Access = public)

        % Edges Orders
        EdgesIDSeq = [] 
        
        % maximum iterations
        MaxIters = 20;
        
        % consumed iteration by SQP
        UsedIters = uint32(0)
        
        % number of loop-closures
        NumLpCls = uint32(0)
        
        % options for choosing linear solver
        % 1 - direct methods
        % 2 - pcg
        LinearSolverOpts = 1
        
        % flag for display debug info
        debug = true
        
    end
    
    properties (Access = protected)
        
        % measurement covariance matrix
        CovMatrix = struct( ...
            'i', [], ...
            'j', [], ...
            'v', [] )
        
        % matrices for linearization
        A = struct( ...
            'i', [], ...
            'j', [], ...
            'v', [] )
        
        C = struct( ...
            'i', [], ...
            'j', [], ...
            'v', [] )
                    
        b = []
        
        d = []
        
        % preserve loop-closure structure to speed up search
        LpCls = struct( ...
            'LpClEdgeID', [], ...
            'LpClMetric', [], ...
            'Enabled', [])
        
        % preserve odometry
        Odoms = struct( ...
            'OdomEdgeID', [])
        
        % mapping LpCls ID to its storage idx
        LpClsID2Idx = []
        
        % Mapping edge ID to its variable index
        EdgesID2VarIdx = []

        % compute preconditioner in first_iteration
        first_iter = true
        
    end
    
    properties (Access = protected)
        
        % System Matrix
        SysMatr
        
        % System Vector
        SysVec
        
        % preconditioner for PCG : InvM = inv(SysMatr)
        InvM
        
    end
    
    
    
    % public methods
    methods (Access = public)
        
        % constructor
        function this = ConstrainedFormulation ()
            % call superclass constructor
            this@PoseGraph();
            % add field 'estRlvPose' 
            % 'estRlvPose' is the variable of constraine formulation
            [ this.Edges.estRlvPose ] = deal({});
            if (nargin == 0)
                % intialization
            end
        end
        
        % overload addEdge to add values to field 'estRlvPose'
        function addEdge(this, id, sid, tid, rlvpose, covinv) 
            % call superclass addEdge
            this.addEdge@PoseGraph(id, sid, tid, rlvpose, covinv);
            % extend the method by adding value to 'Edges.estRlvPose'
            this.Edges(this.NumEdges).estRlvPose = rlvpose;
            % add edge ID to LpClEdgeID if this edge is not odometry edge
            if (tid - sid > 1 || sid - tid > 1) 
                this.NumLpCls = this.NumLpCls + 1;                
                % put edge id into LpClEdgeID
                this.LpCls.LpClEdgeID = [ this.LpCls.LpClEdgeID, id ];
                this.LpCls.LpClMetric = [ this.LpCls.LpClMetric, inf ];
                this.LpCls.Enabled = [ this.LpCls.Enabled, 1];
                % mapping LpCls ID to its storage index (idx)
                this.LpClsID2Idx(id) = this.NumLpCls;
            else
                this.Odoms.OdomEdgeID = [ this.Odoms.OdomEdgeID, id ];
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
        
        % get estEdge(id)
        function [sid, tid, estrlvpose] = estEdge(this, id)
            idx = this.EdgeID2Idx(id);
            sid = this.Edges(idx).snID;
            tid = this.Edges(idx).tnID;
            estrlvpose = this.Edges(idx).estRlvPose;
        end

        % get estEdgeST(sid, tid)
        function [id, estrlvpose] = estEdgeST(this, sid, tid)
            if( tid - sid == 1 || sid - tid == 1 )
                % use OdomEdgeIdx to speed up search
                idx = this.OdomEdgeIdx(sid);
            else
                for ii = this.EdgeID2Idx(this.AdjacentEdges(sid).Source)
                    if(this.Edges(ii).tnID == tid)
                        idx = ii;
                        break
                    end                    
                end
            end
            id = this.Edges(idx).EdgeID;
            estrlvpose = this.Edges(idx).estRlvPose;
        end        
        
        %
        function error = computeObjFunc (this)
            error = 0;
            for i = 1 : this.NumEdges
                vecError = SE3.Log( SE3.Inv(this.Edges(i).RlvPose) * ...
                    this.Edges(i).estRlvPose );
                error = error + ...
                    vecError' * this.Edges(i).CovInv * vecError;
            end
        end
        
        %
        function computePoseEstimate(this)
            for idx = 1 : this.NumNodes
                nodeid = this.Nodes(idx).NodeID;
                if(nodeid > 1)
                    [~, estrlvpose] = this.estEdgeST(nodeid-1, nodeid);
                    this.Nodes(idx).Pose = this.Node(nodeid - 1) * estrlvpose;
                end
            end
        end
        
        %
        function solveBySQP (this)
            % initialize variable ordering and index
            % initialize vector memeory and covariance
            this.initializeProblem();
            % body of SQP solver
            iter = 1;
            while (  iter < this.MaxIters )
                % linearize objective function
                this.linearizeObjFunc ();
                % linearize all enabled constraints
                this.linearizeConsEnabled ();
                % compute Euclidean increment
                vec_inc = this.solveEquilibriumEquation ();
                % update Euclidean increment to SE3 relative pose
                this.updateRlvPoseEstimate (vec_inc);
                % not necessary but cheap
                % this.computePoseEstimate(); 
                norm_inc = norm(vec_inc);
                norm_cons = this.computeConstraintError();
                % * debug *%
                if(this.debug == true)
                    objFunc = this.computeObjFunc ();
                    disp( [ 'norm_increments = ', num2str(norm_inc), ...
                        '      max_norm_constraints = ', num2str(norm_cons), ...
                        '      objFunc = ', num2str(objFunc) ] );
                end
                % *
                %obj.plotTrajectory( iter );  
                % stop criteria
                if(norm_inc < 1e-3 && norm_cons < 1e-3 )
                    break;
                end
                iter = iter + 1;
            end
            % used iterations by SQP
            this.UsedIters = iter;
            % do this each iteration or overall once
            this.computePoseEstimate();      
        end
      
        %
        function initializeVariableOrdering (this)
            % Initialize Edge ID sequence
            this.EdgesIDSeq = [ this.Odoms.OdomEdgeID, ...
                this.LpCls.LpClEdgeID ];
            if (nargin > 1)
                % for example, sort in ascending order
                this.EdgesIDSeq = sort(this.EdgesIDSeq, 'ascend');
            end
        end      
        
        %
        function initializeProblem (this)
            this.initializeVariableOrdering();
            this.initializeVariableIndex();
            this.initializeVectorMemory();
            % obtain relative pose measurement covariance matrix
            this.obtainMeasureCov ();
        end
        
        % add g2o data into the graph
        add_data_g2o (this, FilePathName);
        % plot trajectory of odometry
        [ h1, h2 ] = plotTrajectory(this, fid);
        % dum data
        DumpData3D(this, filename);
    end
    
    
    % private methods
    methods(Access = protected)
        
        %
        function obtainMeasureCov (this)
            ptcov = 1;
            for idx = 1 : this.NumEdges
                % obtain edge information
                id = this.Edges(idx).EdgeID;
                % inverse to obtain covariance
                blk = inv(this.Edges(idx).CovInv);
                % obtain variable index
                varIdx = this.EdgesID2VarIdx(id);
                % put blk into CovMatrix
                this.CovMatrix.v(ptcov:ptcov+35) = reshape(blk, 1, 36);
                this.CovMatrix.i(ptcov:ptcov+35) = ...
                    kron(ones(1,6), varIdx:varIdx+5);
                this.CovMatrix.j(ptcov:ptcov+35) = ...
                    kron(varIdx:varIdx+5, ones(1,6));
                ptcov = ptcov + 36;            
            end  
        end
        
        %
        function linearizeObjFunc (this)
            pta = 1;
            for idx = 1 : this.NumEdges
                % obtain edge
                id = this.Edges(idx).EdgeID;
                rlvpose = this.Edges(idx).RlvPose;
                estrlvpose = this.Edges(idx).estRlvPose;
                % linearization results
                eta = SE3.Log( SE3.Inv(rlvpose) * estrlvpose );
                %blk = SE3.InvJr( eta );
                blk = SE3.Jr( eta );
                % put result into linearization matrix
                % obtain variable index with reagrd to edge ID
                varIdx = this.EdgesID2VarIdx(id);
                % set A and b matrix
                this.A.v(pta:pta+35) = reshape(blk, 1, 36);
                this.A.i(pta:pta+35) = kron(ones(1,6), varIdx:varIdx+5);
                this.A.j(pta:pta+35) = kron(varIdx:varIdx+5, ones(1,6));
                pta = pta + 36;
                this.b(varIdx:varIdx+5, 1) = -eta;
            end
        end
        
        %
        function linearizeConsEnabled (this)
            ptc = 1;
            ptd = 1;
            num_lps = 0;
            for ii = 1 : this.NumLpCls
                if(this.LpCls.Enabled(ii) > 0)
                    num_lps = num_lps + 1;
                    % linearize one constraint with respect to edge id
                    id = this.LpCls.LpClEdgeID(ii);
                    [ Ctmp, dtmp, n_edges] = this.linearizeCons (id);
                    % stack Ctmp to C with respect to rows
                    Ctmp.i = Ctmp.i + 6*(num_lps - 1);
                    this.C.i(ptc:ptc+36*n_edges-1) = Ctmp.i;
                    this.C.j(ptc:ptc+36*n_edges-1) = Ctmp.j;
                    this.C.v(ptc:ptc+36*n_edges-1) = Ctmp.v;
                    ptc = ptc + 36*n_edges;
                    this.d(ptd:ptd+5) = dtmp;
                    ptd = ptd + 6;
                end
            end
        end
        
        %
        function [ Ctmp, dtmp, n_edges] = linearizeCons (this, id)
            % loop-closure edge
            [sid, tid, estlp] = this.estEdge(id);
            n_edges = abs(tid-sid)+1;
            % preallocate space for Ctmp
            Ctmp = struct( ...
                'i', ones(1, 36*n_edges), ...
                'j', ones(1, 36*n_edges), ...
                'v', zeros(1, 36*n_edges));  
            % pointer to current Ctmp index
            ptc = 1;
            % linearize loop-closure edge
            blk = -SE3.Adj(estlp);   
            varIdx = this.EdgesID2VarIdx(id);
            %-% put blk into Ctmp
            Ctmp.v(ptc:ptc+35) = reshape(blk, 1, 36);
            Ctmp.i(ptc:ptc+35) = kron(ones(1,6), 1:6);
            Ctmp.j(ptc:ptc+35) = kron(varIdx:varIdx+5, ones(1,6));
            ptc = ptc + 36;
            %-% 
            % linearize odometry edges
            estodom = eye(4);                
            for idd = sid : tid-1
                [edgeid, estrlvpose] = this.estEdgeST(idd, idd+1);    
                estodom = estodom * estrlvpose;
                blk = SE3.Adj(estodom);
                varIdx = this.EdgesID2VarIdx(edgeid);
                %-% put blk into Ctmp
                Ctmp.v(ptc:ptc+35) = reshape(blk, 1, 36);
                Ctmp.i(ptc:ptc+35) = kron(ones(1,6), 1:6);
                Ctmp.j(ptc:ptc+35) = kron(varIdx:varIdx+5, ones(1,6));
                ptc = ptc + 36;
                %-%  
            end
            eta = SE3.Log( estodom * SE3.Inv(estlp) );
            dtmp = -SE3.Jl( eta ) * eta;
        end
        

        %
        function updateRlvPoseEstimate (this, updatevec)
            for idx = 1 : this.NumEdges
                id = this.Edges(idx).EdgeID;
                varIdx = this.EdgesID2VarIdx(id);
                vec = updatevec(varIdx : varIdx+5);   
                this.Edges(idx).estRlvPose = ...
                    this.Edges(idx).estRlvPose * SE3.Exp(vec);   
            end
        end
        
        %
        function error = computeConstraintError( this )
            error = 0;
            for ii = 1 : this.NumLpCls
                if(this.LpCls.Enabled(ii) > 0)
                    id = this.LpCls.LpClEdgeID(ii);
                    [sid, tid, estlp] = this.estEdge(id);
                    estodom = eye(4);
                    for idd = sid : tid-1
                        [~, estrlvpose] = this.estEdgeST(idd, idd+1);                
                        estodom = estodom * estrlvpose;
                    end
                    error = max(error, norm(SE3.Log(SE3.Inv(estlp) * estodom)));
                end
            end
        end
        
        %
        function initializeVariableIndex( this )
            % intialize edge id order to variable idx
            for idx = 1 : this.NumEdges
                id = this.EdgesIDSeq (idx); 
                this.EdgesID2VarIdx(id) = idx;
            end
            this.EdgesID2VarIdx = 6 * this.EdgesID2VarIdx - 5;
        end
        
        %
        function initializeVectorMemory( this )
            % i,j,v vectors for covariance matrix;
            this.CovMatrix.i = ones(1, this.NumEdges*36);
            this.CovMatrix.j = ones(1, this.NumEdges*36);
            this.CovMatrix.v = zeros(1, this.NumEdges*36);
            % i,j,v vector for objective function
            this.A.i = ones(1, this.NumEdges*36);
            this.A.j = ones(1, this.NumEdges*36);
            this.A.v = zeros(1, this.NumEdges*36);
            this.b = zeros(this.NumEdges*6, 1);
            % i,j,v vector for constraints
            n_edges = 0;
            lp_num = 0;
            for ii = 1 : this.NumLpCls
                if(this.LpCls.Enabled(ii)>0)
                    id = this.LpCls.LpClEdgeID(ii);
                    idx = this.EdgeID2Idx(id);
                    sid = this.Edges(idx).snID;
                    tid = this.Edges(idx).tnID;
                    lp_num = lp_num + 1;
                    n_edges = n_edges +  abs(tid - sid) + 1;
                end
            end            
            this.C.i = ones(1, n_edges*36);
            this.C.j = ones(1, n_edges*36);
            this.C.v = zeros(1, n_edges*36);
            this.d = zeros(lp_num*6, 1);
        end 
        
    end
    
    % private methods
    methods(Access = private)

        %
        function X = solveEquilibriumEquation (this)
            % construct sparse matrices
            CovMatr = sparse(this.CovMatrix.i, ...
                this.CovMatrix.j, ...
                this.CovMatrix.v);                  
            Amatr = sparse(this.A.i, this.A.j, this.A.v);
            Cmatr = sparse(this.C.i, this.C.j, this.C.v, ...
                size(this.d, 1), double(this.NumEdges*6));
            % convert to a minimum norm problem regarding Y
            Q = Amatr * CovMatr * Amatr';
            this.SysMatr = Cmatr * Q * Cmatr';
            ubx = Amatr * this.b;
            this.SysVec = Cmatr * ubx - this.d;
            % Solver for Lambda
            switch this.LinearSolverOpts
                case 1
                    Lambda = this.LinearSolverDirect ();
                case 2
                    if (this.first_iter)
                        this.computePreconditioner ();
                        this.first_iter = false;
                    end
                    Lambda = this.LinearSolverPCG ();
                otherwise
                    error('undefined linear solver type');
            end
            % recover X from Lambda
            X = ubx - Q * Cmatr' * Lambda;
        end
        
        %
        function computePreconditioner (this)
            tinv = tic();
            this.InvM = inv(this.SysMatr);
            TimeSysInversion = toc(tinv)
        end
            
            
    end
    
    methods (Access = private)  
        % solvers for AX = b
        
        %
        function X = LinearSolverDirect (this)
            X = this.SysMatr \  this.SysVec;
        end 
        
        %
        function X = LinearSolverPCG (this)
            [X, iters] = ConstrainedFormulation.pcg_pr( ...
                this.SysMatr, ...
                this.SysVec, ...
                zeros(size(this.SysVec,1),1), ...
                this.InvM, ...
                1000, ...
                1e-8);
            disp(['pcg iters = ', num2str(iters)]);          
        end
                
    end
    
    methods(Static = true)
        
        %
        function [x, i] = pcg_pr(A, b, x, InvM, mit, tol)
            i = 1;
            rp = b - A * x;
            d = InvM * rp;
            rsold = rp' * d;
            rs0 = rsold;
            disp(['tol*rso = ', num2str(tol*rs0)])
            while (i < mit)
                Ad = A * d;
                alpha = rsold / (d' * Ad);
                x = x + alpha * d;
                if( mod(i+1,50) == 0 )
                    rn = b - A * x;
                else
                    rn = rp - alpha * Ad;
                end
                s = InvM * rn;
                rsnew = rn' * s;
                rsmid = rp' * s;
                if sqrt(rsnew) < tol*rs0
                %if sqrt(rsnew) < tol
                    break;
                end                
                beta = (rsnew - rsmid)/rsold;
                d = s + beta * d;
                rsold = rsnew;
                rp = rn;
                i = i + 1;                
            end
        end
        
        
    end
    
    
end


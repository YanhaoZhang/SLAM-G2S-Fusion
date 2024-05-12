classdef IncrementalNLSFramework < NLSFormulation

    % 
    
    properties
        
        ConsMetric = []
        
        ObjGrowth = []
        
        NoiseSchme = 1
        
    end
    
    properties (Access = public)
        
        PoseCov
        
    end
   
    
    methods (Access = public)
        
        % constructor
        function this = IncrementalNLSFramework ()
            % call superclass constructor
            this@NLSFormulation();
            if (nargin == 0)
                % initialization
            end
        end
        
        % addNode
        function addNode(this, id, pose)
            % call superclass addNode
            this.addNode@NLSFormulation(id, pose);
        end
        
        % overload addEdge to initialize LpCls.Enabled = false
        function addEdge(this, id, sid, tid, rlvpose, covinv)
            % call superclass addEdge
            this.addEdge@NLSFormulation(id, sid, tid, rlvpose, covinv);
            if (tid - sid > 1 || sid - tid > 1)
                this.LpCls.Enabled(this.NumLpCls) = 0; % false
            end
        end
        
    end
    
    
    methods (Access = public)
        
        % batch incremental GaussNewton
        function solveByBatchIncrementalGaussNewton (this)
            % compute the problem with odometry only
            this.solveByGaussNewton();
            this.computeMarginal ();
            fprintf(1, '\nFinished initialization, the original poses and its marginal computed!\n\n' );
            while(true)
                % unoptimized loop-closures
                idxs = find(this.LpCls.Enabled == 0);
                if(size(idxs,2)==0)
                    disp('all loop-closures have been optimized');
                    break;
                end
                % compute constraint metrics for un-optimized loop-closures
                for idx = idxs
                    id = this.LpCls.LpClEdgeID(idx);
                    this.LpCls.LpClMetric(idx) = ...
                        this.computeLpClMetric (id);
                end
                 % Enable loop-closure with lowest metric value
                [ v, i ] = min(this.LpCls.LpClMetric(idxs));
                idx = idxs(i);
                this.LpCls.Enabled(idx) = 1;
                % solve by GN
                objFunc1 = this.computeObjFunc();
                this.solveByGaussNewton();
                this.computeMarginal();
                objFunc2 = this.computeObjFunc();
                fprintf(1, ['Summary >>   ConstraintMetric = ', num2str(v), ...
                    '     ObjectiveGrowth = ', num2str(objFunc2 - objFunc1), '\n\n']);
                %
                this.ConsMetric = [ this.ConsMetric, v ];
                this.ObjGrowth = [ this.ObjGrowth, objFunc2-objFunc1 ];
                %
               % this.plotTrajectory();
            end
        end
        
        % online incremental GaussNewton
        function solveByOnlineIncrementalGaussNewton (this)
            % compute the problem with odometry only
            this.solveByGaussNewton();
            this.computeMarginal ();
            fprintf(1, '\nFinished initialization, the original poses and its marginal computed!\n\n' );
            while(true)
                % unoptimized loop-closures
                idx = find(this.LpCls.Enabled == 0, 1, 'first');
                if(size(idx,2)==0)
                    disp('all loop-closures have been optimized');
                    break;
                end
                % compute constraint metrics for un-optimized loop-closures
                id = this.LpCls.LpClEdgeID(idx);
                v = this.computeLpClMetric (id);
                if(v < 12)
                    this.LpCls.LpClMetric(idx) = v;
                    this.LpCls.Enabled(idx) = 1;
                else
                    % do something else
                end
                % solve by GN
                objFunc1 = this.computeObjFunc();
                this.solveByGaussNewton();
                this.computeMarginal();
                objFunc2 = this.computeObjFunc();
                fprintf(1, ['Summary >>   ConstraintMetric = ', num2str(v), ...
                    '     ObjectiveGrowth = ', num2str(objFunc2 - objFunc1), '\n\n']);
                %
                this.ConsMetric = [ this.ConsMetric, v ];
                this.ObjGrowth = [ this.ObjGrowth, objFunc2-objFunc1 ];
                %
               % this.plotTrajectory();
            end   
        end
        
        %-%
        function metr = computeLpClMetric (this, id)
            idx = this.EdgeID2Idx(id);
            % edge error vector by logarithm mapping
            edgeError = this.Edges(idx).Error;
            % covariane of the edge error vector
            cov = this.computeEdgeErrorCov (id);
            % metric
            metr = edgeError' *(cov\edgeError);
        end
        
    end
    
    methods (Access = public)

        % 
        function cov = computeEdgeErrorCov (this, id)
            idx = this.EdgeID2Idx(id);
            sid = this.Edges(idx).snID;
            tid = this.Edges(idx).tnID;
            rlvpose = this.Edges(idx).RlvPose;
            covinv = this.Edges(idx).CovInv;
            if(this.NoiseSchme == 1)
                pose1 = this.Node(sid);
                pose2 = this.Node(tid);
                estrlvpose = SE3.Inv(pose1) * pose2;
                J1 = -SE3.InvAdj( SE3.Inv(rlvpose)*estrlvpose ); % mst
                J2 = -SE3.InvAdj( estrlvpose ); % snode
                J3 = eye(6); % tnode
            end
            % obtain edge error covariance
            J = [J1, J2, J3];
            stjcov = this.obtainJointMarginal (sid, tid);
            cov = J * blkdiag(inv(covinv), stjcov) * J';  
        end
        
        %
        function cov = obtainJointMarginal (this, id1, id2)
            varIdx1 = this.NodesID2VarIdx(id1);
            varIdx2 = this.NodesID2VarIdx(id2);
            if(varIdx1>0 && varIdx2>0) % no fixed nodes
                cov = this.PoseCov( ...
                    [varIdx1:varIdx1+5, varIdx2:varIdx2+5], ...
                    [varIdx1:varIdx1+5, varIdx2:varIdx2+5] );
            elseif(varIdx1<0 && varIdx2>0) % fix node 1
                cov = blkdiag( zeros(6,6), ...
                    this.PoseCov(varIdx2:varIdx2+5, varIdx2:varIdx2+5) );
            elseif(varIdx1>0 && varIdx2<0) % fix node 2
                cov = blkdiag( ...
                    this.PoseCov(varIdx1:varIdx1+5, varIdx1:varIdx1+5), ...
                    zeros(6,6) );
            elseif(varIdx1<0 && varIdx2<0) % fix node1 && node2
                cov = zeros(12,12);
            end
        end
        
        %
        function computeMarginal (this)
            cPinv = sparse(this.CovMatrInv.i, ...
                                   this.CovMatrInv.j, ...
                                   this.CovMatrInv.v);
            cA = sparse(this.A.i, this.A.j, this.A.v);
            posecovinv = cA' * cPinv * cA;
            % compute pose covariance directly
            this.PoseCov = inv(posecovinv); % 
        end
        
    end
    
end


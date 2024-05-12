classdef IncrementalConstrainedFramework < ConstrainedFormulation

    % 
    
    properties
        
        % estimated edge measurement covariance matrix
        % corrected by objective and constraint coefficients
        estCovMatrix
        
        ConsMetric = []
        
        ObjGrowth = []
        
        test = true
        
    end
    
    methods (Access = public)
        
        % constructor
        function this = IncrementalConstrainedFramework ()
            % call superclass constructor
            this@ConstrainedFormulation();
            if (nargin == 0)
                % intialization
            end
        end
        
        % overload addEdge to initialize LpCls.Enabled = false
        function addEdge(this, id, sid, tid, rlvpose, covinv) 
            % call superclass addEdge
            this.addEdge@ConstrainedFormulation( ...
                id, sid, tid, rlvpose, covinv);
            if (tid - sid > 1 || sid - tid > 1) 
                this.LpCls.Enabled(this.NumLpCls) = 0; % false
            end
        end
        
        % update edge covariance
        function updateEdgeCovariance (this)
            % construct sparse matrices
            SqrtCovMatr = sparse(this.SqrtCovMatrix.i, ...
                this.SqrtCovMatrix.j, ...
                this.SqrtCovMatrix.v);                  
            Amatr = sparse(this.A.i, this.A.j, this.A.v);
            Cmatr = sparse(this.C.i, this.C.j, this.C.v, ...
                size(this.d, 1), double(this.NumEdges*6));
            % convert to a minimum norm problem regarding Y
            H = Cmatr * Amatr;
            h = this.d - H * this.b;
            H = H * SqrtCovMatr;
            % covariance correction by loop-closure coefficient
            cov_rh = - H' * inv(H*H') * H;
            J = Amatr * SqrtCovMatr;
            this.estCovMatrix = J * J' + J * cov_rh * J';
        end
        
        % compute loop-closure metric with id
        function vmetr = computeConstraintMetric (this, id)
            [ A, ~, ~] = this.linearizeCons (id);
            Amatr = sparse(A.i, A.j, A.v, 6, double(this.NumEdges*6));
            LpClCovMatr = Amatr * this.estCovMatrix * Amatr';
            % estRlvPose is intialized as RlvPose for new data @addEdge
            [sid, tid, estlp] = this.estEdge(id);
            estodom = eye(4);
            for idd = sid : tid-1
                [~, estrlvpose] = this.estEdgeST(idd, idd+1);
                estodom = estodom * estrlvpose;
            end
            vecError = SE3.Log(SE3.Inv(estodom) * estlp);
            vmetr = vecError' * (LpClCovMatr\vecError);
            %vmetr = vecError' * vecError;
        end
        
        % solve a dataset by batch incremental SQP
        function solveByBatchIncrementalSQP (this)
            this.solveBySQP();
            this.updateEdgeCovariance();
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
                        this.computeConstraintMetric(id);
                end
                % Enable loop-closure with lowest metric value
                [ v, i ] = min(this.LpCls.LpClMetric(idxs));            
                idx = idxs(i);
                this.LpCls.Enabled(idx) = 1;
                % solve by SQP
                objFunc1 = this.computeObjFunc();
                this.solveBySQP();
                objFunc2 = this.computeObjFunc();
                disp(['ConstraintMetric = ', num2str(v), ...
                    '     ObjectiveGrowth = ', num2str(objFunc2 - objFunc1)]);
                %
                this.ConsMetric = [ this.ConsMetric, v ];
                this.ObjGrowth = [ this.ObjGrowth, objFunc2-objFunc1 ];
                % compute Edge covariance : expensive
                this.updateEdgeCovariance();
                this.plotTrajectory();
            end 
        end
        
    end
    
end


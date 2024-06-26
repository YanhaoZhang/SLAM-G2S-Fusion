classdef PoseGraph < handle
    
    % graph structure for a specific type of optimization problem
    % --------------------- methods ---------------------
    % function [pose, estpose] = Node(obj, id)
    % function [rlvpose, covinv] = Edge(obj, id)
    % function [rlvpose, covinv] = EdgeST(obj, sid, tid)
    % function [ id ] = locateEdge(obj, sid, tid)
    % function [ flag ] = existEdge(obj, sid, tid)
    % function [ nids ] = predecessors(obj, nid)
    % function [ nids ] = successors(obj, nid)
    % function obj = addNode(obj, id, pose)
    % function obj = addEdge(obj, id, sid, tid, rlvpose, covinv)
    
    properties (Access = public)
        
        % define Nodes and Edges with structure     
        Nodes = struct( ...
            'NodeID', {}, ...
            'Pose', {} )
        
        Edges = struct( ...
            'EdgeID', {}, ...
            'snID', {}, ...
            'tnID', {}, ...
            'RlvPose', {}, ...
            'CovInv', {} )
        
        % Num of Edges
        NumEdges = 0

        % Num of Nodes
        NumNodes = 0
                      
    end
    
    properties (Access = protected)
        % nodes and edges connecting to a specific node
        
        % nodes ID in adjacent list
        AdjacentNodes = struct( 'Source', {}, 'Target', {} )
        
        % edges ID in adjacent list
        AdjacentEdges = struct( 'Source', {}, 'Target', {} )            
        
        % several idx mapping vectors, reserve 10000
        
        % store odometry edge idx
        OdomEdgeIdx %= zeros(1, 10000, 'uint32')

        % map node ID to its storage index (idx)
        NodeID2Idx %= zeros(1, 10000, 'uint32')
        
        % map edge ID to its storage index (idx)
        EdgeID2Idx %= zeros(1, 10000, 'uint32')

    end
    
    
    
    methods
        
        % constructor
        function this = PoseGraph(nodes, edges)
            if(nargin > 0)
                if(isfield(node, {'NodeID', 'Pose'}))
                    for ii = 1 : size(nodes, 2)
                        this.addNode( ...
                            nodes(ii).NodeID, ...
                            nodes(ii).Pose );
                    end
                end
                if(isfield(edge, {'EdgeID', 'snID', 'tnID', 'RlvPose', 'CovInv'}))
                    for ii = 1 : size(edges, 2)
                        this.addEdge( ...
                            edges(ii).EdgeID, ...
                            edges(ii).snID, ...
                            edges(ii).tnID, ...
                            edges(ii).RlvPose, ...
                            edges(ii).CovInv );
                    end
                end
            end
        end
        
        % get node(id)
        function [pose] = Node(this, id)
            idx = this.NodeID2Idx(id);
            pose = this.Nodes(idx).Pose;
        end
        
        % get edge(id)
        function [sid, tid, rlvpose, covinv] = Edge(this, id)
            idx = this.EdgeID2Idx(id);
            sid = this.Edges(idx).snID;
            tid = this.Edges(idx).tnID;
            rlvpose = this.Edges(idx).RlvPose;
            covinv = this.Edges(idx).CovInv;
        end
        
        % get edge(sid, tid)
        function [id, rlvpose, covinv] = EdgeST(this, sid, tid)
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
            rlvpose = this.Edges(idx).RlvPose;
            covinv = this.Edges(idx).CovInv;
        end
        
        % locate edge id for edge(sid, tid)
        function [ id ] = locateEdge(this, sid, tid)
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
        end
        
        % check if edge(sid, tid) exist in the graph
        function [ flag ] = existEdge(this, sid, tid)
            flag = ismember(sid, ...
                this.AdjacentNodes(tid).Source);
        end
        
        % predecessors of a node in the graph
        function [ nids ] = predecessors(this, nid)
            nids = this.AdjacentNodes(nid).Source;
        end
        
        % successors of a node in the graph
        function [ nids ] = successors(this, nid)
            nids = this.AdjacentNodes(nid).Target;
        end
        
        % add node
        function addNode(this, id, pose)          
                count = this.NumNodes + 1;
                this.NumNodes = count; 
                % preserve node
                this.Nodes(count).NodeID = id;
                this.Nodes(count).Pose = pose;
                % mapping node id to its storage index (idx)
                this.NodeID2Idx(id) = count;
                % initialize AdjacentNodes
                this.AdjacentNodes(id)= ...
                    struct( 'Source', [], 'Target', [] );
                % initialize AdjacentEdges
                this.AdjacentEdges(id) = ...
                    struct( 'Source', [], 'Target', [] );               
        end
        
        % add edge
        function addEdge(this, id, sid, tid, rlvpose, covinv)         
                count = this.NumEdges + 1;
                this.NumEdges = count;
                % preserve edge
                this.Edges(count).EdgeID = id;
                this.Edges(count).snID = sid;
                this.Edges(count).tnID = tid;
                this.Edges(count).RlvPose = rlvpose;
                this.Edges(count).CovInv = covinv;
                % mapping edge id to its storage index (idx)
                this.EdgeID2Idx(id) = count;
                % add to AdjacentNodes
                % # sources nodes of a given node
                this.AdjacentNodes(tid).Source = ...
                    [ this.AdjacentNodes(tid).Source, sid ];
                % # target nodes of a given node
                this.AdjacentNodes(sid).Target = ...
                    [ this.AdjacentNodes(sid).Target, tid ];
                % add to AdjacentEdges
                % # edges using a given node as source
                this.AdjacentEdges(sid).Source = ...
                    [ this.AdjacentEdges(sid).Source, id ];
                % # edges using a given node as target
                this.AdjacentEdges(tid).Target = ...
                    [ this.AdjacentEdges(tid).Target, id ];
                % storage idx for odometry edges
                if( tid - sid == 1 || sid - tid == 1 )
                    this.OdomEdgeIdx(sid) = count;
                end                
        end
        
    end

    
end


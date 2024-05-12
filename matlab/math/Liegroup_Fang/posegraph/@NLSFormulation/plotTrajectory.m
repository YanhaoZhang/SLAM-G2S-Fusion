function [ h1, h2 ] = plotTrajectory(this, fid)

    if(nargin > 1)
        figure(fid)
    else
        figure(100)
        clf
    end
    
    Matr = [];
    % ASSUME : NodeID = NodeIdx
    for idx = 1 : this.NumNodes
        pose = this.Nodes(idx).Pose;
        t = pose((1:3), 4);
        Matr = [ Matr,  t];
    end
    
    if(~isreal(Matr))
        fprintf(2, 'Complex matrix: Matr @plotTrajectory \n');
        Matr = real(Matr);
    end
    if(Matr(3,:) == 0)
        PlanarTrajectory = true;    
    else
        PlanarTrajectory = false;
    end
    
    if(PlanarTrajectory)
        h1 = plot(Matr(1,:), Matr(2,:), '-r.');
    else
        h1 = plot3(Matr(1,:), Matr(2,:), Matr(3,:), '-r.');
    end
    axis equal
    hold on 
    % plot loop-closure edges
    for ii = 1 : this.NumLpCls
        if ( this.LpCls.Enabled(ii) )
            id = this.LpCls.LpClEdgeID(ii);
            idx = this.EdgeID2Idx(id);
            sid = this.Edges(idx).snID;
            tid = this.Edges(idx).tnID;
            edge = Matr(:, [sid, tid]);
            if(PlanarTrajectory)
                h2 = plot(edge(1,:), edge(2,:), 'b.--');
            else
                h2 = plot3(edge(1,:), edge(2,:), edge(3,:), 'b.--');
            end
        end
    end

    hold off
    pause(0.01);
end


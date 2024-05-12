function [ h ] = plotTrajectory(this, fid)

    if(nargin > 1)
        figure(fid)
    else
        figure(100)
        clf
    end
    
    hold on
    
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
        h = plot(Matr(1,:), Matr(2,:), 'ro');
    else
        h = plot3(Matr(1,:), Matr(2,:), Matr(3,:), 'ro');
    end
    
    
%     for idx = 1 : this.NumEdges
%         if ( this.isEnabled(this.Edges(idx).EdgeID) )
%             sid = this.Edges(idx).snID;
%             tid = this.Edges(idx).tnID;
%             edge = Matr(:, [sid, tid]);
%             if(tid - sid == 1)
%                 if(PlanarTrajectory)
%                     h = plot(edge(1,:), edge(2,:), 'g.--');
%                 else
%                     h = plot3(edge(1,:), edge(2,:), edge(3,:), 'g.--');
%                 end
%             else
%                 if(PlanarTrajectory)
%                     h = plot(edge(1,:), edge(2,:), 'b.--');           
%                 else
%                     h = plot3(edge(1,:), edge(2,:), edge(3,:), 'b.--');
%                 end   
%             end
%         end
%     end    
    
    hold off
    
    pause(0.01);
    
end


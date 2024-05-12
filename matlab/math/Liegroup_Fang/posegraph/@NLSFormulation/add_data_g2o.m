% add g2o data into graph
function add_data_g2o (this, FilePathName)

EdgeID = 1;
fid = fopen(FilePathName,'r');
g2oLine = fgetl(fid);

while g2oLine(1) ~=-1
    
    TypeName = strtok(g2oLine);
    
    %-% 2D node
    if strcmp(TypeName, 'VERTEX_SE2')
        
        lineInfo = textscan(g2oLine, '%*s %d %f %f %f');
        
        % compute pose
        PoseIndex = lineInfo{1};
        angle = lineInfo{4};
        axis = [ 0; 0; 1; ];        
        R = SO3.Exp(angle*axis);        
        t = [ lineInfo{2};
                lineInfo{3};
                0; ];            
        pose = [ R,            t;
                     zeros(1,3),  1]; 
        % add node into graph
        NodeID = PoseIndex + 1;           
        this.addNode(NodeID, pose);      
        
    %-% 2D edge    
    elseif strcmp(TypeName, 'EDGE_SE2')
        
        lineInfo = textscan(g2oLine, '%*s %d %d %f %f %f %f %f %f %f %f %f');

        % extract R, t
        StartingPose = lineInfo{1};
        EndingPose = lineInfo{2};  
        angle = lineInfo{5};
        axis = [ 0; 0; 1; ];
        R = SO3.Exp(angle*axis);
        t = [ lineInfo{3};
                lineInfo{4};
                0; ];
        % construct relative-pose
        % transformation to ensure sid < tid
        if(EndingPose>StartingPose)
            sid = StartingPose+1;
            tid = EndingPose+1;
            rlvpose = [ R,    t;
                            zeros(1,3), 1 ];
        elseif(EndingPose<StartingPose)
            sid = EndingPose+1;
            tid = StartingPose+1;
            rlvpose = [ R',   -R' * t;
                            zeros(1,3), 1 ];            
        else
            error('self-cycle in the graph: startnode = targetnode!');
        end
        
        % initialize pose estimates
        if(this.NumNodes == 0)
            this.addNode(1, eye(4));
        end
        if(this.NumNodes < sid)
            this.addNode(sid, this.Nodes(this.NumNodes).Pose * rlvpose);
        end
        if(this.NumNodes < tid)
            this.addNode(tid, this.Nodes(this.NumNodes).Pose * rlvpose);
        end
        
        % extract invere-covariance
        covinv2d = ... 
            [lineInfo{6}, lineInfo{7}, lineInfo{8};
             lineInfo{7}, lineInfo{9}, lineInfo{10};
             lineInfo{8}, lineInfo{10}, lineInfo{11}; ];
        % construct 3D inverse-covariance by setting redundant variables' inverse-covariance to identity
        covinv = eye(6);           
        covinv([1,2,6], [1,2,6]) = covinv2d;
        % add edge into graph
        this.addEdge(EdgeID, sid, tid, rlvpose, covinv);
        EdgeID = EdgeID + 1;   
        
    %-% 3D node
    elseif strcmp(TypeName, 'VERTEX_SE3:QUAT')
        
        lineInfo = textscan(g2oLine, '%*s %d %f %f %f %f %f %f %f');
        
        % extract pose
        PoseIndex = lineInfo{1};
        q_x = lineInfo{5};
        q_y = lineInfo{6};
        q_z = lineInfo{7};
        q_w = lineInfo{8}; 
        R = R_q([q_w; q_x; q_y; q_z]);
        t = [ lineInfo{2};
                lineInfo{3};
                lineInfo{4}; ];
        pose = [ R,            t;
                     zeros(1,3),  1]; 
        % add node into graph
        NodeID = PoseIndex + 1;
        this.addNode(NodeID, pose);       
        
    %-% 3D edge
    elseif strcmp(TypeName, 'EDGE_SE3:QUAT')
        
        lineInfo = textscan(g2oLine, '%*s %d %d %f %f %f %f %f %f %f    %f %f %f %f %f %f    %f %f %f %f %f   %f %f %f %f   %f %f %f    %f %f    %f'); 

        % extract R, t
        StartingPose = lineInfo{1};
        EndingPose = lineInfo{2};
        q_x = lineInfo{6};
        q_y = lineInfo{7};
        q_z = lineInfo{8};
        q_w = lineInfo{9};
        R = R_q([q_w; q_x; q_y; q_z]);
        t = [ lineInfo{3};
                lineInfo{4};
                lineInfo{5}; ];
        % construct relative-pose
        % transformation to ensure sid < tid
        if(EndingPose>StartingPose) 
            sid = StartingPose+1;
            tid = EndingPose+1;
            rlvpose = [ R,    t;
                            zeros(1,3), 1 ];
        elseif(EndingPose<StartingPose)
            sid = EndingPose+1;
            tid = StartingPose+1;
            rlvpose = [ R',   -R' * t;
                            zeros(1,3), 1 ];
        else
            error('self-cycle in the graph: startnode = targetnode!');
        end
        
        % initialize pose estimates
        if(this.NumNodes == 0)
            this.addNode(1, eye(4));
        end
        if(this.NumNodes < sid)
            this.addNode(sid, this.Nodes(this.NumNodes).Pose * rlvpose);
        end
        if(this.NumNodes < tid)
            this.addNode(tid, this.Nodes(this.NumNodes).Pose * rlvpose);
        end
        
        % extract invere-covariance
        covinv = ...
            [lineInfo{10}, lineInfo{11}, lineInfo{12}, lineInfo{13}, lineInfo{14}, lineInfo{15};
            lineInfo{11}, lineInfo{16}, lineInfo{17}, lineInfo{18}, lineInfo{19}, lineInfo{20};
            lineInfo{12}, lineInfo{17}, lineInfo{21}, lineInfo{22}, lineInfo{23}, lineInfo{24};
            lineInfo{13}, lineInfo{18}, lineInfo{22}, lineInfo{25}, lineInfo{26}, lineInfo{27};
            lineInfo{14}, lineInfo{19}, lineInfo{23}, lineInfo{26}, lineInfo{28}, lineInfo{29};
            lineInfo{15}, lineInfo{20}, lineInfo{24}, lineInfo{27}, lineInfo{29}, lineInfo{30};];
        % compute axis-angle of rotation
        angle = 2*acos(q_w);
        if(angle<1e-10)
            axis = [0;0;1];
        else
            axis = [q_x; q_y; q_z]/sin(angle/2);
        end
        if ( angle > pi)
            angle = 2*pi - angle;
            axis = -axis;
        end
        AngleAxis = angle*axis;
        % compute Jacobian: quaternion with respect to axis-angle
        J = J_dq_daa( AngleAxis );
        % inverse-covariance of axis-angle rotation
        covinv((4:6), (4:6))= J((2:4),:)' * covinv((4:6), (4:6)) * J((2:4),:);
        % add edge into graph
        this.addEdge(EdgeID, sid, tid, rlvpose, covinv);
        EdgeID = EdgeID + 1;
        
    end
    
    g2oLine = fgetl(fid);
    
end

fclose(fid);

fprintf(1, '\nSucceed in adding data: "%s"\n\n', FilePathName);

end

% transform a quaternion into a rotation matrix
function [ RotationMatrixQuaternion ] = R_q( q )
q_w = q(1);
q_x = q(2);
q_y = q(3);
q_z = q(4);

RotationMatrixQuaternion = [q_w^2+q_x^2-q_y^2-q_z^2,  2*(q_x*q_y-q_w*q_z),    2*(q_x*q_z+q_w*q_y);
                            2*(q_x*q_y+q_w*q_z),    q_w^2-q_x^2+q_y^2-q_z^2,  2*(q_y*q_z-q_w*q_x);
                            2*(q_x*q_z-q_w*q_y),      2*(q_y*q_z+q_w*q_x),  q_w^2-q_x^2-q_y^2+q_z^2;];

end

% Jacobian \frac{dq}{daa}
function [ J ] = J_dq_daa( AxisAngle )

t = norm(AxisAngle);
c = cos(t/2);
s = sin(t/2);

if(t<1e-10)
    J = [ 0 0 0;
          1/2 0 0;
          0 1/2 0;
          0 0 1/2;];
else
    J = [ -s/(2*t)*AxisAngle(1),  -s/(2*t)*AxisAngle(2),   -s/(2*t)*AxisAngle(3);
          c/(2*t^2)*AxisAngle(1)^2+s/(t^3)*(AxisAngle(2)^2+AxisAngle(3)^2),  (c/(2*t^2)-s/(t^3))*AxisAngle(1)*AxisAngle(2),  (c/(2*t^2)-s/(t^3))*AxisAngle(1)*AxisAngle(3);
          (c/(2*t^2)-s/(t^3))*AxisAngle(1)*AxisAngle(2),  c/(2*t^2)*AxisAngle(2)^2+s/(t^3)*(AxisAngle(1)^2+AxisAngle(3)^2),  (c/(2*t^2)-s/(t^3))*AxisAngle(2)*AxisAngle(3);
          (c/(2*t^2)-s/(t^3))*AxisAngle(1)*AxisAngle(3),  (c/(2*t^2)-s/(t^3))*AxisAngle(2)*AxisAngle(3),   c/(2*t^2)*AxisAngle(3)^2+s/(t^3)*(AxisAngle(1)^2+AxisAngle(2)^2);];
end

end



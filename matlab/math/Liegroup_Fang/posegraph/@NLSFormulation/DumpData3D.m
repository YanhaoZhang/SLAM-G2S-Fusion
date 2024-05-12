function DumpData3D (this, filename)


edge_str3d = 'Edge3 %d %d';
pose_str3d =  '%f %f %f %f %f %f';
cov_str3d = ['%f %f %f %f %f %f ', ...
                        '%f %f %f %f %f ', ...
                              '%f %f %f %f ', ...
                                   '%f %f %f ', ...
                                        '%f %f ', ...
                                             '%f'];
   
format_str3d = [edge_str3d, ' ', pose_str3d, ' ', cov_str3d, '\n'];
    
% format_str2d = 'Edge2 %d %d %f %f %f %f %f %f %f %f %f';

fid = fopen(filename, 'w');

for idx = 1 : this.NumEdges
  
    id1 = this.Edges(idx).snID-1;  %convert to zero based id
    
    id2 = this.Edges(idx).tnID-1;   %convert to zero based id
    
    rlvpose = this.Edges(idx).RlvPose;
    
    covinv = this.Edges(idx).CovInv;
    
   
    tvec = rlvpose([1:3], 4);
    rvec = SO3.Log( rlvpose([1:3], [1:3]) );
    
    
    
    fprintf(fid, format_str3d, id1, id2, tvec(1), tvec(2), tvec(3), rvec(1), rvec(2), rvec(3), ...
        covinv(1,1), covinv(1,2), covinv(1,3), covinv(1,4), covinv(1,5), covinv(1,6), ...
                           covinv(2,2), covinv(2,3), covinv(2,4), covinv(2,5), covinv(2,6), ...
                                              covinv(3,3), covinv(3,4), covinv(3,5), covinv(3,6), ...
                                                                 covinv(4,4), covinv(4,5), covinv(4,6), ...
                                                                                    covinv(5,5), covinv(5,6), ...
                                                                                                       covinv(6,6) );
    
end

end
function [interId_] = intersection(A, B)

% find the common elements of A and B. Without remove the redundant
% elements. B must be a subset of A

numB = size(B,1);

interId = zeros(numB*2,1);   % initialize using numB*2 elements just in case there are duplicate elements in A

idx=0;
for i=1:numB

id_tmp = find(A == B(i));


for l=1:length(id_tmp)
    idx = idx+1;
    interId(idx) = id_tmp(l);
end
end


% interId_ = find(interId);

% remove the rest zero elements
interId_ = interId(1:idx);


end
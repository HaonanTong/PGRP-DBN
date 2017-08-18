function [ C ] = f_alignment( A, B, cntr )
% A, B two structure contains matrix, at
% A.matrix = [1:9;11:19;21:29];
% A.at = 4;
% B.matrix = [31:39;41:49];
% B.at = 4;
% C = f_alignment(A,B,0);

T1 = A.at;
T2 = B.at;
lA = size(A.matrix,2);
lB = size(B.matrix,2);
l2 = lA-T1;
l4 = lB-T2;
if T1 == T2 && lA == lB
    if ~mod(cntr,2)
       C.matrix = [A.matrix;B.matrix];
       C.at = A.at;
    else
       C.matrix = [B.matrix;A.matrix];
       C.at = A.at; 
    end
    
   return
end

if T1 > T2
    C = f_alignment(B,A,cntr+1);
    return 
end

if T1 < T2
    B.matrix = B.matrix(:,(T2-T1+1):end);
    B.at = A.at;
    C = f_alignment(A,B,cntr);
    return
end

% l1 = l2
if l2 < l4
    C = f_alignment(B,A,cntr+1);
    return
end

if l2 > l4
    A.matrix = A.matrix(:,1:lB);
    C = f_alignment(A,B,cntr);
    return
end

end


function mat = cl2mat(cl)
% mat = cl2mat(cl)
% Convert cluster matrix to matrix
% Input: cl: n-vector for class labels
% Output: mat: n*k binary matrix

% Author: Bowei Yan
% Last modified: June 13, 2017

n = length(cl);
mat = zeros(n);
for i = 1:n,
    for j = 1:n,
        mat(i,j) = 1*(cl(i)==cl(j));
    end
end

function mat = cl2mat(cl)
% mat = cl2mat(cl)
% Convert cluster matrix to matrix
% Input: cl: n-vector for class labels
% Output: mat: n*k binary matrix

% Author: Bowei Yan
% Last modified: June 14, 2017

n = length(cl);
k = length(unique(cl));
mat = zeros(n,k);

for j = 1:k,
    mat(:,j) = 1*(cl==j);
end


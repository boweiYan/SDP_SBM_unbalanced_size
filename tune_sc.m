function [l,C_best] = tune_sc(A,Y,k)
% [l,C_best] = tune_sc(A,Y,k)
% ACASC with tuning; find the tuning parameter by k-means loss (as suggested in the paper)
% Ref:
% N. Binkiewicz, J. T. Vogelstein, K. Rohe; Covariate-assisted spectral clustering. 
% Biometrika 2017; 104 (2): 361-377. doi: 10.1093/biomet/asx008

% Input: A: network
%        Y: feature matrix
%        k: number of clusters
% Output: l: chosen lambda
%         C_best: final clustering

% Author: Bowei Yan
% Last modified: June 13, 2017

[n,~] = size(A);
tau = mean(sum(A,1));
% regularize on D not A
D = sqrt(sum(A,1)+tau);
L = inv(diag(D))*A*inv(diag(D));
eig_L = eigs(L*L,k+1);
eig_X = eigs(Y*Y',k+1);
lower = (eig_L(k)-eig_L(k+1))/eig_X(1);
upper = eig_L(1)/(eig_X(1));
step = 10;
stepsize = (upper-lower)/step;

lambda = lower:stepsize:upper;
O = zeros(length(lambda),1);
for i=1:length(lambda),
    L_tau = inv(diag(D))*A*inv(diag(D))+lambda(i)*Y*Y';
    cl_acasc_reg(i,:) = rsc(L_tau,k,'adj');
    [u,~] = eigs(L_tau,k);
    for class = 1:k,
        idx = find(cl_acasc_reg(i,:)==class);
        center = mean(u(idx,:),1);
        for j = 1:length(idx),
            O(i) = O(i)+norm(u(idx(j),:)-center)^2;
        end
    end
end
[~,min_i] = min(O);
l = lambda(min_i);
C_best = cl_acasc_reg(min_i,:);
end
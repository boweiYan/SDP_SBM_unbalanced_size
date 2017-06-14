function [class,U1] = rsc(A, K, method, prior)
% Regularized Spectral Clustering
% [class,U1] = rsc(A, K, method, prior)
% Input: A: Adjacency matrix
%        K: number of clusters
%        method: 'pos' - Regularized Spectral Clustering
%                'lap' - Spectral Clustering (use normalized graph laplacian)
%                'adj' - Use adjacency matrix
% Output: class labels and eigenvectors
% Author: Bowei Yan
% Last Updated: Mar 14, 2017

    if nargin==3,
        prior = 1/K*ones(K,1);
    end
    [nv,~] = size(A);
    tau = mean(A(:));
    if strcmp(method,'pos')
        A_tau = A + tau*ones(nv,nv);
        L_tau = normalizeSym(A_tau);
    elseif strcmp(method,'lap')
        L_tau = normalizeSym(A);
    elseif strcmp(method,'adj')      
        L_tau = (A+A')/2;
    end
    [U1,~] = eigs(L_tau,K);
    % Row renormalization for spectral clustering
    %if strcmp(method,'noreg'),
        U1 = normr(U1);
    %end
    % K-means with restarts
    maxsum = inf;
    nrestart = 100;
    for i=1:nrestart,
        [class0,~,sumD] = kmeans(U1(:,1:K),K);
	if maxsum > sum(sumD),
		maxsum = sum(sumD);
		class=class0;
	end
    end
end

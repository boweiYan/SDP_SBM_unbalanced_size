% Example Community Detection for network

% Author: Bowei Yan
% Last Updated: June 14, 2017

%% Generate data
n = 300; k = 3;
clsize = ones(k,1);
clsize = clsize/sum(clsize);
p = 0.3; q = 0.05;
prob = q*ones(k)+(p-q)*eye(k);
rho = 1;
[A,Z,~,~] = create_block_model(n,rho,prob,clsize,1);

alpha = k/n;
opts = struct('rho',1,'T',10000,'tol',0.1,'report_interval',20,'quiet',0);

Xhat = admm_imb(A,k,alpha,opts);

figure,subplot(121),imagesc(A),title('adjacency')
subplot(122),imagesc(Xhat);title('Solution of SDP');

%% To get labels and NMI
cl = rsc(Xhat, k, 'adj');
nmi(cl,Z)
function [X,itr,rs,ier]= admm_clustering1(W,k_cluster,opts1)
%
% ADMM for solving the following SDP
% Used for clustering Stochastic Block Model with unknown cluster number
% max  trace(WX)
% s.t. X\succeq 0, X1=1, 0\le X\le 1
% Input: 
%     W: A+\lambda I, where A is adjacency matrix, and lambda is a tuning
%     parameter
%     opts1: options
%         max.ite:   max iteration        
%         tol: tolerance for stopping criterion
%         report_interval: frequency to print intermediate result
% 
% Output:
%     X: optmal solution
%     delta: relative error when converge
%     T_term: number of iteration taken to converge
% Author: Xiuyuan Cheng
% Modified by: Bowei Yan
% Last Updated: Jul 18, 2017
%

n=size(W,1);
assert( norm(W-W','fro') < eps*10); 

%% extended problem

%
C=-W;
%nrmC = sum(abs(C(:)))+1;
nrmC = norm(C,1)+1;

%
m=n+1;

funcAX1=@(X) sum(X,2);
b1=ones(m-1,1);

funcAX2=@(X) trace(X);
b2=k_cluster;

funcAX=@(X) [funcAX1(X);funcAX2(X)];
b=[b1;b2];


% As=zeros(n,n,m);%[n,n,m]
% for i=1:m-1
%     Ai=zeros(n,n);
%     Ai(:,i)=1;
%     Ai=(Ai+Ai')/2;%Ai needs to be symmetric
%     As(:,:,i)=Ai;
% end
% As(:,:,m)=eye(n);
% vecA =reshape(As,n^2,m); %[n^2,m]
% 
% funcATy=@(y) reshape( sum(bsxfun(@times,vecA,y'),2), ...
%                      n,n);
                 
% A*AT
% AAT=vecA'*vecA; %eye(m)*m + ??
% invAAT=inv(AAT);
AAT = zeros(n+1);
AAT(1:n,1:n) = (n*eye(n)+ones(n))/2;
AAT(n+1,n+1) = n;
AAT(1:n,n+1) = ones(n,1);
AAT(n+1,1:n) = ones(1,n);
invAAT = inv(AAT);


%% parameter of AltSDP
tol=opts1.tol;

opts.mxitr = opts1.max_ite;
opts.record = 0;

opts.mu = 5; %initial mu
opts.rmu = 0.5;
opts.mu_min = 1e-4;
opts.mu_max = 1e4;

opts.min_mu_itr= 5;

opts.max_mu_itr= 100; %h4

opts.max_itr_stag_lev1 = 20; %h1
opts.max_itr_stag_lev2 = 50; %h2
opts.max_itr_stag_lev3 = 150; %h3

opts.rho = 1.6; %rho to update X

opts.delta_mu_l = 1;
opts.delta_mu_u = 1;

opts.ftol = tol;

%% initialize

% record of residual
rs=zeros(opts.mxitr ,4);

%
X = eye(n); %

y=zeros(m,1);
S=zeros(n,n);
Z=zeros(n,n);

mu=opts.mu;

%%

%gap=|b.y-<C,X>|/(1+|b.y|+<C,X>)
ppobj=trace(C'*X);
pdobj = 0;
gap = abs(ppobj-pdobj)/(1.+abs(ppobj)+abs(pdobj)); 

% pinf = |A(X)-b|_2/(1+|b|)
resi=funcAX(X)-b;
pinf = norm(resi,2)/(1+norm(b,2));

% dinf=|A^*(y)+S+Z-C|/(1+|C|)
dinf = 1;

itmu_pinf = 0;
itmu_dinf = 0;
itr_stag = 0;

ref_inf = max([pinf,dinf,gap]);

%%
for itr=1:opts.mxitr
    
%     %%
%     pobj_old = ppobj;  
%     dobj_old = pdobj;
        
    %% step 1: compute y
    
    % y=-(AA*)^(-1)(mu*(A(X)-y)+A(S+Z-C))
    y=-invAAT*(mu*resi+funcAX(S+Z-C));
    
    %% step2: compute Z
    %W =  C - ATy -S- mu*X;  WPos = (W>0);       
    %Z = zeros(n); Z(WPos) = W(WPos);
    ATy=funcATy(y);
    
    W=C-ATy-S-mu*X;
    Z=max(W,0);
    
    %% step 3: compute S by projection
    % V =  C - ATy-Z - mu*X
    V=C-ATy-Z-mu*X;
   
    % 
    [uV,dV]=eig(V);
    [evalsV,tmp]=sort(diag(dV),'descend');
    evecsV=uV(:,tmp);
    
    % S=V pos part
    idpos=find(evalsV>0);
    S=evecsV(:,idpos)*diag(evalsV(idpos))*evecsV(:,idpos)';
    
    S=(S+S')/2;
    
    %% step 3: compute X
    %X=(1/mu)*(S-V);
    Xnew=(1/mu)*(S-V);
    
    X=(1-opts.rho)*X+opts.rho*Xnew;
    
        
    %% check optimality
    
    %gap=|b.y-<C,X>|/(1+|b.y|+<C,X>)
    ppobj= trace(C'*X);
    pdobj = b'*y;
    gap = abs(ppobj-pdobj)/(1.+abs(ppobj)+abs(pdobj));
    
    
    % pinf = |A(X)-b|_2/(1+|b|)
    resi=funcAX(X)-b;
    pinf = norm(resi,2)/(1+norm(b,2));
    
    
    % dinf=|ATy+S+Z-C|/(1+|C|)
    dinf = norm(funcATy(y)+S+Z-C,'fro')/nrmC;
    
    %
    dtmp = max([pinf,dinf,gap]);
    
    % optimal solution
    if( dtmp <= opts.ftol )
        ier=1;
        break;
    end
    
    %
    if mod(itr,opts1.report_interval)==1
        fprintf('%d | %4.2e\t%4.2e\t gap=%4.2e\t %4.2e\n',itr,pinf,dinf,gap,mu);
    end
    
    rs(itr,:)=[pinf,dinf,gap,mu];
    
    %% detection of stagnation
    if dtmp < ref_inf
        ref_inf = dtmp;
        itr_stag = 0;
    else
        itr_stag = itr_stag + 1;
    end
    
    if (  (itr_stag > opts.max_itr_stag_lev1 &&  dtmp<= opts.ftol*10)...
            || ( itr_stag > opts.max_itr_stag_lev2 &&  dtmp<= opts.ftol*100)...
            || ( itr_stag > opts.max_itr_stag_lev3 && dtmp <= opts.ftol*1000) )
        ier=2;
        break;
    end
    
    %% update mu
    
    % line 417 of mexAltSDP_ThetaPlus
    if pinf/dinf <= opts.delta_mu_l
        itmu_pinf = itmu_pinf + 1;
        itmu_dinf = 0;
        if itmu_pinf > opts.max_mu_itr
            mu=max(mu*opts.rmu, opts.mu_min);
            itmu_pinf=0;
        end
    elseif pinf/dinf >opts.delta_mu_u
        itmu_dinf=itmu_dinf+1;
        itmu_pinf=0;
        if itmu_dinf > opts.max_mu_itr
            mu = min(mu/opts.rmu, opts.mu_max);
            itmu_dinf=0;
        end
    end
    
    
end

ier=3;
end

function aty = funcATy(y)
n=length(y)-1;
aty = zeros(n);
idmat = eye(n);
for i=1:n,
    aty = aty+y(i)/2*(ones(n,1)*idmat(i,:)+idmat(:,i)*ones(1,n));
end
aty = aty+y(end)*eye(n);
end
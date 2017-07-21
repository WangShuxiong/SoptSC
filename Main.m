function [P,nC,W,idx,eigenvalues] = Main(nC,data)
alpha = 0.2;

realdata = data;
realdata = realdata-min(realdata(:));
realdata = realdata./max(realdata(:));


[~,n] = size(realdata);
for i = 1:n
    realdata(:,i) = realdata(:,i)/norm(realdata(:,i));
end

lambda = 0.5;
K = ceil(alpha*n);
if K<=10
    K = 10;
elseif K >= 20
    K = 20;
end

[W,P] = SimilarityM(realdata,lambda,K);
W(W<=eps) = 0;

%% Determinning the number of clusters
eigenvalues = [];
if isempty(nC)
    [eigenvalues,No_cluster] = Num_cluster(W);
    nC = No_cluster;
end

% [decomW,~] = eigs(W,nC);
% params.Hinit = 2*full(sqrt(mean(mean(W))/nC))*decomW;
flag = 1;
[chuzhiA,~] = nndsvd(W,nC,flag);
params.tol = 10^(-6);
params.Hinit = chuzhiA;
[HH,~,~] = symnmf_newton(W, nC, params);

[~,idx] = max(HH,[],2);
end
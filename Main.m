function [P,No_cluster,W,idx,eigenvalues,H] = Main(nC,data)
% This is the main function to run SoptSC.
%
% Input:
%          nC: Number of cluster.
%        data: A m*n matrix with m rows (genes) and n columns (cells).
%
% Output:
%           W: Cell-to-cell similarity matrix
%           P: Transition matrix
%  No_cluster: Number of cluster computed by SoptSC if nC = [];
%              otherwise, No_cluster = nC
%         idx: Cluster labels
% eigenvalues: Eigenvalues of graph Laplacian of the consensus matrix
%
%


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

% Determine the number of clusters
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
H = HH;
[~,idx] = max(HH,[],2);
No_cluster  = nC;
end
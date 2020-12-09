function obj = run_SoptSC(obj,opt)
% G_filter_idx
% SOptSC identifies clusters from single cell data
%
% Input
%   -data:      An m*n single cell data matrix with m rows(genes) and n columns(cells)
%   -NC:        Number of Cluster specified by user, if NC = [] (default), then the
%               algorithm will compute the number.
%
%   -No_exc_cell: Gene selection parameter range from [0,No_cells] (No_cells represents
%                 the number of cells), we remove genes that are expressed less than
%                 No_exc_cell cells and more than (No_cells - No_exc_cell)
%                 cells (Default value: No_exc_cell = 6)
%   -No_features: Gene selection parameter, which represent the number of highly expressed
%                 genes selected for clustering (Default: No_features = 2000)
%
%
%
% Output
%   --  W: Cell-to-cell similarity matrix.
%   --  No_cluster: Number of clusters
%   --  cluster_label: cluster labels for all cells.
%   --  latent: An nx2 matrix representing low dimensional latent space for all cells
%   --  H: Non-negative matrix such that W = H*H^T
%   --  eigenvalues: eigenvalues of the graph Laplacian of the truncated
%       consensus matrix which is used to estimate the number of clusters

data = full(obj.data);

NC = opt.No_cluster;
No_exc_cell = opt.No_exc_cell;
No_features = opt.No_features;



if nargin==1
    NC = [];
    No_exc_cell = 6;
    No_features = 2000;
elseif nargin == 2
    No_exc_cell = 6;
    No_features = 2000;
elseif nargin == 3
    No_features = 2000;
end
[No_genes,No_cells] = size(data);

% Data preprocess: Gene filtering
if No_exc_cell > 0
    gene_nnz = zeros(No_genes,1);
    alpha_filter = No_exc_cell./No_cells;
    for i = 1:No_genes
        gene_nnz(i) = nnz(data(i,:))./No_cells;
    end
    G_filter_idx1 = union(find(gene_nnz<=alpha_filter),find(gene_nnz>=1-alpha_filter));
    G_filter_idx = setdiff(1:No_genes,G_filter_idx1);
else
    G_filter_idx = 1:size(data,1);
end

data1 = data(G_filter_idx,:);

[coeff, score, pca_eigvalue] = pca(data1');
[~,No_Comps] = max(abs(pca_eigvalue(2:end-1) - pca_eigvalue(3:end)));

aa = max(coeff(:,1:No_Comps+1)');
bb = sort(aa,'descend');

if size(data1,1) <=1000
    No_sel_genes = size(data1,1);
else
    No_sel_genes = min([No_features round(size(data1,1))]);
end

gene_selection = find(aa>=bb(No_sel_genes));
input_data = data1(gene_selection,:);


eigenvalues = [];
nC = NC;
[No_cluster,W,idx,eigenvalues,H] = Main(nC,input_data,opt);

cluster_label = idx;

cluster_name = cell(No_cluster,1);
for i = 1:No_cluster
    cluster_name{i} = ['C' num2str(i)];
end

obj.W = sparse(W);
obj.No_cluster = No_cluster;
obj.cluster_label = cluster_label;
obj.H = H;
obj.eigenvalues = eigenvalues;
obj.cluster_name = cluster_name;

end



%%%%%%%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%
function [No_cluster,W,idx,eigenvalues,H] = Main(nC,data,opt)
% This function perform the SoptSC algrithm.
%
% Input:
%          nC: Number of cluster.
%        data: A m*n matrix with m rows(genes) and n columns(cells).
%
% Output:
%           W: Cell-to-cell similarity matrix.
%  No_cluster: Number of cluster computed by SoptSC if nC = [];
%             otherwise, No_cluster = nC
%         idx: Cluster label
% eigenvalues: Eigenvalues of graph Laplacian of the consensus matrix
%           H: Non-negative matrix such that W = H*H^T

realdata = data;
realdata = realdata-min(realdata(:));
realdata = realdata./max(realdata(:));

[~,n] = size(realdata);
for i = 1:n
    realdata(:,i) = realdata(:,i)/norm(realdata(:,i));
end

lambda = 0.5;

W = SimilarityM(realdata,lambda,data,opt);

WB = W;
n = size(W,1);
D = diag(WB*ones(n,1));
Prw = eye(size(W)) - D^(-1/2)*WB*D^(-1/2);
if n>=1000
    No_eigs = 100;
    all_eigs = real(eigs(Prw,No_eigs,'sm'));
else
    all_eigs = real(eig(Prw));
end

ZZ = sort(abs(real(all_eigs)));
No_cluster1 = length(find(ZZ<=0.01));

%% Determinning the number of clusters
eigenvalues = [];
if isempty(nC)
    [eigenvalues,No_cluster] = Num_cluster(W,No_cluster1);
    nC = No_cluster;
end


flag = 1;
[chuzhiA,~] = nndsvd(W,nC,flag);
params.tol = 10^(-6);
params.Hinit = chuzhiA;
[H,~,~] = symnmf_newton(W, nC, params);
[~,idx] = max(H,[],2);
No_cluster  = nC;
end


%%%%%%%%%%%%%%%%%%%%% Similarity matrix Function %%%%%%%%%%%%%%%%%%%%%%%
function W = SimilarityM(X,lambda,data,opt)
%
% Computing cell-to-cell similarity matrix by solving the following
% optimization problem via ADMM
%
%    min_{Z,E}  ||Z||_* + lambda ||E||_{2,1}
%    s.t.       X = XZ + E;
%               Z'1 = 1;
%               Z_{i,j} = 0 for (i,j)\in Omega
%
% Input
%   1) X: Single cell data, a m*n matrix, with m rows(genes) and n columns(cells);
%   2) lambda: the default value is 0.5;
%
% Output
%   1) W: Cell-to-cell similarity matrix


[m,n] = size(X);
% KNN Search for finding nearest neighbors
if m>=60
    [coeff1,X1,pca_eigvalue1] = pca(X','NumComponents',60);
else
    [coeff1,X1,pca_eigvalue1] = pca(X','NumComponents',m);
end


[~,No_Comps1] = max(abs(pca_eigvalue1(2:end-1) - pca_eigvalue1(3:end)));

cc = cumsum(pca_eigvalue1(2:end));
dd = cc(2:end)./sum(pca_eigvalue1(2:end));

K1 = length(find(dd<=0.3));

if n > 5000
    K = 100;
else
    if K1 <= 10
        K = 10;
    elseif K1 >=30
        K = 30;
    else
        K = K1+1;
    end
end

if isfield(opt,'K')
    K = opt.K;
end

display(K);
dim_init = 3;

% InitY = pca(data,'NumComponents',dim_init);
% X2 = tsne(X','Standardize',true,'Perplexity',20,'NumDimensions',dim_init,'InitialY',InitY);

ppx = 30;
% dim_init = 2;
[X2, ~, ~]=run_umap(X','n_neighbors',ppx,'n_components',dim_init,...
    'verbose','text','method','Java','metric','correlation'); % row -- data point; column -- feature; correlation


D = ones(n,n);
if No_Comps1>=1
    No_Comps1 = 1;
end

[IDX,~] = knnsearch(X2(:,1:No_Comps1+2),X2(:,1:No_Comps1+2),'k',K);

for jj = 1:n
    D(jj,IDX(jj,:)) = 0;
end

Z = computM(D,X,lambda);
Z(Z<=eps) = 0;
W = 0.5.*(abs(Z)+abs(Z'));



end


function Z = computM(D,X,lambda)
%% ADMM iteration
maxiter = 100;
Err = zeros(maxiter,2);
rho = 5;        % 5
mu = 10^(-6);   % 10^(-6)
mumax = 10^(6);
epsilon = 10^(-5);
[m,n] = size(X);
Z = zeros(n); E = zeros(m,n);
Y1 = zeros(m,n); Y2 = zeros(1,n); Y3 = zeros(n,n);
iter = 0;

% display('Iter  Err');
while 1
    iter = iter + 1;
    if iter >= maxiter
        break;
    end
    
    % step 1: Update J
    mu = min(rho*mu,mumax);
    [U,S,V] = svd(Z-Y3/mu);
    
    R = length(diag(S));
    Dmu = zeros(R,R);
    MM = max(diag(S)-(1/mu)*ones(R,1),zeros(R,1));
    %     for i = 1:R
    %         Dmu(i,i) = MM(i);
    %     end
    Dmu(1:R+1:end) = MM;
    
    J = U*Dmu*V';
    
    if iter >= 3
        if Err(iter-1,1) >= Err(iter-2,1) || max(max(Err(iter,:)),norm(X-X*Z)) <= epsilon
            break;
        end
    end
    
    
    % Update E
    Q = X-X*Z+Y1/mu;
    %     for j = 1:n
    %         if norm(Q(:,j)) > lambda/mu
    %             E(:,j) = ((norm(Q(:,j))-lambda/mu)/norm(Q(:,j)))*Q(:,j);
    %         else
    %             E(:,j) = zeros(length(Q(:,j)),1);
    %         end
    %     end
    
    aa = max(sqrt(sum(Q.^2)) - lambda/mu,0)./sqrt(sum(Q.^2));
    E = Q.*aa;
    
%     eta = norm(X)^2 + norm(ones(n,1))^2+1;
    eta = sum( X(:).^2)  + n +1;
    H = -X'*(X - X*Z - E + (1/mu)*Y1) - ones(n,1)*(ones(n,1)'- ones(n,1)'*Z + (1/mu)*Y2) + ...
        (Z - J + (1/mu)*Y3);
    Z = Z - (1/eta)*H;
    Z(D>0) = 0;
    
    % Update Dual variable
    Y1 = Y1 + mu*(X-X*Z-E);
    Y2 = Y2 + mu*(ones(1,n) - ones(1,n)*Z);
    Y3 = Y3 + mu*(Z - J);
    
%     fprintf('%d, %8.6f\n',iter,norm(X-X*Z-E));
%     Err(iter,:) = [norm(X-X*Z-E) norm(Z-J)];
    Err(iter,:) = [sqrt( sum(sum((X-X*Z-E).^2)) ) sqrt(sum(sum((Z-J).^2))) ];
%     if max(Err(iter,:)) <= epsilon
%         break;
%     end
end
end

function [eigenvalues,No_cluster] = Num_cluster(W,NC1)
% Computing number of cluster
nno_cluster = 2:20;  % 2:20

if NC1 <= 5
    tau = 0.3; % 0.4
elseif NC1 <= 10
    tau = 0.4;
else
    tau = 0.5;
end

tol = 0.01;
[all_eigs,M_all] = consen(W,nno_cluster,tau);

zz = sort(abs(real(all_eigs)));

if length(zz)>=2
    gap = zz(2:end) - zz(1:end-1);
    [~,No_cluster1] = max(gap);
end

No_cluster2 = length(find(zz<=tol));
No_cluster = No_cluster1;
display('Number of cluster based on zero eigenvalues & Largest gap ');
display([No_cluster2 No_cluster]);

eigenvalues = zz;
end


function [all_eigs,M_all] = consen(S,idx,tau)

flag = 1;
n = size(S,1);

M_all = zeros(size(S));
[LS,~] = pca(S,'NumComponents',3);
for i = 1:length(idx)
%     [chuzhiA,~] = nndsvd(S,i,flag);
%     params.tol = 10^(-6);
%     params.Hinit = chuzhiA;
%     
%     [H_all,~,~] = symnmf_newton(S,i,params);
%     [~,indx] = max(H_all, [], 2);
    indx = kmeans(LS,i,'Start',zeros(i,3));
    M = consen(indx);
    M_all = M_all + M;
end

M_all(M_all <= tau*length(idx)) = 0;
M_all = (1/2)*(M_all + M_all');

D = diag(M_all*ones(n,1));
Prw = eye(n) - D^(-1/2)*M_all*D^(-1/2);
% Prw = eye(n) - D^(-1)*M_all;
% P = (eye(n)/D)*M_all;
% all_eigs = eig(P);

if n>=1000
    No_eigs = 100;
    all_eigs = real(eigs(Prw,No_eigs,'smallestreal'));%'smallestreal''smallestabs'
else
    all_eigs = real(eig(Prw));
end

% all_eigs = real(eig(Prw));

    function M = consen(indx)
        
        % [m,~] = size(S);
        m = length(indx);
        
        M = zeros(m);
        
        for ii = 1:m
            for j = 1:m
                if indx(ii)==indx(j)
                    M(ii,j) = 1;
                end
            end
        end
    end
end

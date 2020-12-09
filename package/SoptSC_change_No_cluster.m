function obj = SoptSC_change_No_cluster(obj,No_cluster)

W = obj.W;
nC = No_cluster;

flag = 1;
[chuzhiA,~] = nndsvd(W,nC,flag);
params.tol = 10^(-6);
params.Hinit = chuzhiA;
[H,~,~] = symnmf_newton(W, nC, params);
[~,idx] = max(H,[],2);

cluster_name = cell(No_cluster,1);
for i = 1:No_cluster
    cluster_name{i} = ['C' num2str(i)];
end

obj.H = H;
obj.cluster_label = idx;
obj.No_cluster = No_cluster;

obj.cluster_name = cluster_name;











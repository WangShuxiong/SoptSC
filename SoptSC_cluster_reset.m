function [No_cluster,cluster_label,H] = SoptSC_cluster_reset(W,No_cluster_by_usr)
% Reclustering cells based on given No_of_clusters by user and cell-cell
% similarity matrix W

No_cluster = No_cluster_by_usr;

flag = 1;
[chuzhiA,~] = nndsvd(W,No_cluster,flag);
params.tol = 10^(-6);
params.Hinit = chuzhiA;
[H,~,~] = symnmf_newton(W, No_cluster, params);
[~,cluster_label] = max(H,[],2);


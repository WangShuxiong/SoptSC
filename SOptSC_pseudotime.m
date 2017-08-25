function cell_order = SOptSC_pseudotime(init_cluster,P,cluster_label,latent,No_cluster)
%SoptSC infers pseudotime from single cell data
%   This function perform pseudotime estimation by SoptSC.
%
%   Input:
%   -- init_cluster: starting cluster specified by user
%   -- P: transition matrix
%   -- cluster_label: cluster labels for all cells.
%   -- No_cluster: Number of clusters
%   -- latent: low dimension space induced from cell-cell transition matrix
%   
%   Output:
%   --  cell_order: cell order inferred by SOptSC
%
%% Reconstruct pseudotime
nC = No_cluster;
%[dvis1,~] = eigs((P+P')/2,5);
% [dvis1,~] = eigs((P+P')/2);
dvis1 = pca(P,6);
[~,~,W1,~] = Main(nC,dvis1');

%% Rank-1 NMF
flag = 1;
[chuzhiA,~] = nndsvd(W1,1,flag);
params.tol = 10^(-6);
params.Hinit = chuzhiA;
[H_all,~,~] = symnmf_newton(W1, 1, params);

ZZ = find(cluster_label==init_cluster); 
[~,init_idx] = max(H_all(ZZ));

init_point = ZZ(init_idx);
H1_all = abs(H_all - H_all(init_point));
[~,cell_order] = sort(H1_all);

% Pseudotime visualization
c = linspace(0,1,size(W1,1));
colormap parula;
cmap = colormap;
mymap = cmap(1:58,:);
colormap(mymap);
figure(1);
scatter(latent(cell_order,2),latent(cell_order,3),40,c,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);

box on;
cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;
set(gca,'LineWidth',1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
print(1,'-dtiff','Results\pseudotime.tiff');

% Cell lineage inference
root_cell = init_point;
[Tree,pred,cluster_center] = lineage(cluster_label,P,root_cell);
end


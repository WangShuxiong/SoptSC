function [W,P,No_cluster,cluster_label,cell_order] = SOptSC(data,init_point,NC)
% SOptSC: Similarity-based Optimization for Single Cell data
%
% Input
%   -- data:
%       a m*n single cell data matrix with m rows (genes) and n columns (cells)
%   -- init_point:
%       index of the initial point for pseudotemporal ordering and cell lineage inference
%   -- NC:
%       Number of Cluster specified by user, if NC = [] (default),
%       then SOptSC will predict the number of clusters.
%   -- plot: if plot = 1;  plots the result (default).
%            if plot = 0;  plots off.
%
% Output
%   --  W: Cell-to-cell similarity matrix.
%   --  P: Transition matrix
%   --  No_cluster: Number of clusters, if NC is specified by user,
%       No_cluster = NC; otherwise, SOptSC will compute NC.
%   --  cluster_label: cluster labels for all cells.
%   --  cell_order: cell order inferred by SOptSC
%

if nargin==2
    NC = [];
    %     plot = 1;
end

if nargin==3
    %     plot = 1;
end


nC = NC;
[P,No_cluster,W,idx,eigenvalues] = Main(nC,data);
[dvis,~] = eigs(P);

%% Reconstruct pseudotime
nC = No_cluster;
% [dvis1,~] = eigs((P+P')/2);
% [dvis1,~] = eigs((P+P')/2);
dvis1 = pca(P,6);
[~,~,W1,~] = Main(nC,dvis1');

%% Rank-1 NMF
flag = 1;
[chuzhiA,~] = nndsvd(W1,1,flag);
params.tol = 10^(-6);
params.Hinit = chuzhiA;
[H_all,~,~] = symnmf_newton(W1, 1, params);
H1_all = abs(H_all - H_all(init_point));
[~,cell_order] = sort(H1_all);

% Pseudotime visualization
c = linspace(0,1,size(W,1));
colormap parula;
cmap = colormap;
mymap = cmap(1:58,:);
colormap(mymap);
figure(1);
%scatter(dvis(cell_order,2),dvis(cell_order,3),40,c,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
scatter(dvis(cell_order,2),dvis(cell_order,3),40,c,'filled');

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

%% Visualization of subpopulations
figure(2);
for ik = 1:No_cluster
    %scatter(dvis(find(idx==ik),2),dvis(find(idx==ik),3),40,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    scatter(dvis(find(idx==ik),2),dvis(find(idx==ik),3),40,'filled');
    hold on;
end
box on;
set(gca,'LineWidth',1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ytick',[]);

lgd = cell(1,No_cluster);
for i = 1:No_cluster
    if i<10
        vv = 'ClusterC';
        vv(8:8) = num2str(i);
        lgd{i} = vv;
    else
        vv = 'ClusterCC';
        vv(8:9) = num2str(i);
        lgd{i} = vv;
    end
end
legend(lgd,'FontSize',10,'Location','best');%,'Orientation','horizontal');
set(gca,'FontName','Arial');
set(gca,'FontSize',12);

print(2,'-dtiff','Results\subpopulation.tiff');

% Cell lineage inference
root_cell = init_point;
[Tree,pred,cluster_center] = lineage(idx,P,root_cell);

% Display eigengap of the graph Laplacian
if isempty(NC)
    figure(3);
    scatter(1:min([30 size(W,1)]),eigenvalues(1:min([30 size(W,1)])),20,'filled');
    box on;
    set(gca,'LineWidth',1.5);
    xlabel('i');
    ylabel('Eigenvalues of graph Laplacian \lambda_i');
    set(gca,'FontName','Arial');
    set(gca,'FontSize',12);
    print(3,'-dtiff','Results\EigenGap.tiff');
end
cluster_label = idx;
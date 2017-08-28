function [W,P,No_cluster,cluster_label,latent] = SOptSC_cluster(data,NC)
% SOptSC identifies clusters from single cell data
%
% Input
%   -- data:
%       a m*n single cell data matrix with m rows(genes) and n columns(cells)
%   -- NC:
%       Number of Cluster specified by user, if NC = [] (default), then the algorithm will compute the number.
%   -- plot: if plot = 1;  plot the result (default).
%            if plot = 0;  plot of  f.
%
% Output
%   --  W: Cell-to-cell similarity matrix.
%   --  P: Transition matrix
%   --  No_cluster: Number of clusters, if NC is specified by user,
%       No_cluster = NC; otherwise, the method will compute one.
%   --  cluster_label: cluster labels for all cells.
%   --  latent: low dimension space induced from cell-cell transition
%   matrix

if nargin==1
    NC = [];
end


nC = NC;
[P,No_cluster,W,idx,eigenvalues] = Main(nC,data);
[dvis,~] = eigs(P);


%% Subpopulations visualization
figure(2);
for ik = 1:No_cluster
    scatter(dvis(find(idx==ik),2),dvis(find(idx==ik),3),40,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
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



% display eigen-gap of graph Laplacian
if isempty(NC)
    figure(3);
    scatter(1:min([30 size(W,1)]),eigenvalues(1:min([30 size(W,1)])),20,'filled');
    box on;
    set(gca,'LineWidth',1.5);
    xlabel('i');
    ylabel('Eigenvalue of graph Laplacian \lambda_i');
    set(gca,'FontName','Arial');
    set(gca,'FontSize',12);
    print(3,'-dtiff','Results\EigenGap.tiff');
end
cluster_label = idx;
latent = dvis;
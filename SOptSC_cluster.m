function [W,No_cluster,cluster_label,latent,H,Gene_sel_idx] = SOptSC_cluster(data,NC,No_exc_cell,No_features)
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
%   -W:             Cell-to-cell similarity matrix.
%   -No_cluster:    Number of clusters inferred by SoptSC
%   -cluster_label: Cluster labels identified by SoptSC.
%   -latent:        Low dimension projection of W by SoptSC
%   -H:             Non-negative low-rank matrix such that W = H*H^T
%   -Gene_sel_idx:  Gene indices indentified by SoptSC 


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

% Data preprocess
% Gene filtering
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

[coeff, ~, pca_eigvalue] = pca(data1');
[~,No_Comps] = max(abs(pca_eigvalue(2:end-1) - pca_eigvalue(3:end)));
display(No_Comps);

aa = max(coeff(:,1:No_Comps+1)');
bb = sort(aa,'descend');

if size(data1,1) <=1000
    No_sel_genes = size(data1,1);
else
    No_sel_genes = min([No_features round(size(data1,1))]);
end
% No_sel_genes = round(1*size(data1,1));

% display(No_sel_genes);
gene_selection = find(aa>=bb(No_sel_genes));
input_data = data1(gene_selection,:);


Gene_sel_idx = G_filter_idx(gene_selection);
display(length(Gene_sel_idx));

nC = NC;
[No_cluster,W,idx,eigenvalues,H] = Main(nC,input_data);

W1 = W./(ones(1,size(W,1))*W*ones(size(W,1),1));
dvis = pca1(W1,2);

%% Subpopulations visualization
switch1 = 1;
if switch1==1
    figure(1);
    for ik = 1:No_cluster
        scatter(dvis(find(idx==ik),1),dvis(find(idx==ik),2),40,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
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
    
    print(1,'-dtiff', strcat(ResFolder,'/Subpopulation.tiff'));
     
    % display eigen-gap of graph Laplacian
    if isempty(NC)
        figure(2);
        scatter(1:min([30 size(W,1)]),eigenvalues(1:min([30 size(W,1)])),20,'filled');
        box on;
        set(gca,'LineWidth',1.5);
        xlabel('i');
        ylabel('Eigenvalue of graph Laplacian \lambda_i');
        set(gca,'FontName','Arial');
        set(gca,'FontSize',12);
        print(2,'-dtiff', strcat(ResFolder,'/EigenGap.tiff'));  
    end
end
cluster_label = idx;
latent = dvis;
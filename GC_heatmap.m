function Gene_labels = GC_heatmap(data,cluster_label,H,No_exc_cell,No_select_genes)
% This function assign genes to each cluster by SoptSC
%
%   Input:
%       data: full gene-cell data matrix
%       cluster_label: cluster labels for all cells
%       H: non-negative matrix such that W = H*H^T
%       No_exc_cell: Gene selection parameter range from [0,No_cells] (No_cells represents
%                   the number of cells), we remove genes that are expressed less than 
%                   No_exc_cell cells and more than (No_cells - No_exc_cell)
%                   cells (Default value: No_exc_cell = 6)
%       No_select_genes: maximal number of genes to be ploted
%
%   Output:
%       Gene_labels: gene label information for each gene associated with a specific
%       cluster.
%       -- 1st columns of Gene_labels represents gene indices;
%       -- 2nd column of Gene_labels represents the cluster index that the gene belongs to; 
%       -- 3rd column represent gene score associated with corresponding cluster.

[~,gene_idx] = Data_selection(data,No_exc_cell,No_select_genes);

data = data(gene_idx,:);

NC = max(cluster_label);
m = size(data,1);
Gene_labels = zeros(m,3);

Gene_labels(:,1) = gene_idx;

%% data normalization
for i = 1:size(data,2)
    data(:,i) = data(:,i)./norm(data(:,i),2);
end

% for i = 1:m
%     data(i,:) = data(i,:)./norm(data(i,:),2);
% end


G_latent = data*H;
[Gene_value,Gene_label] = max(G_latent,[],2);

Gene_labels(:,2) = Gene_label;
Gene_labels(:,3) = Gene_value;


OGI = [];
OGV = [];

CGI = [];
for i = 1:NC
    % order genes within each cluster
    Z = find(Gene_label==i);
    Z1 = Gene_value(Z);
    [Z1V,I] = sort(Z1,'descend');
    
    Z2 = Z(I);
    OGI = [OGI; Z2];
    OGV = [OGV; Z1V];
    
    %order cells
    Y = find(cluster_label==i);
    CGI = [CGI; Y];
end

%% allgenes string to cell
% allgs = cell(size(allgenes));
% for i = 1:length(allgs)
%     allgs{i} = allgenes{i};
% end
% %% plot gene-cell heatmap
% RowLabelsValue = allgs(gene_idx);
% % HeatMap(data(OGI,CGI),'RowLabels', ColumnLabelsValue)
% HMdata = data(OGI,CGI);
% 
% HeatMap(HMdata,'RowLabels', RowLabelsValue,'Standardize',2,'DisplayRange',2,'Colormap',redbluecmap)


%% data normalization and zscroe
% kk = 1 row; kk = 2 column
figure;
idata = data(OGI,CGI);
kk = 2;
center = mean(idata,kk);
scale = std(idata, 0,kk);

tscale = scale;
%=Check for zeros and set them to 1 so not to scale them.
scale(tscale == 0) = 1;
%== Center and scale the data
idata = bsxfun(@minus, idata, center);
sdata = bsxfun(@rdivide, idata, scale);

colormap redbluecmap;
clims = [-3 3];
imagesc(sdata,clims);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

% yticks(1:size(data,1));
% yticklabels(allgenes(OGI));

cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;

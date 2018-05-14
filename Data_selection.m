function [input_data,Gene_sel_idx] = Data_selection(data,No_exc_cell,No_features)
% G_filter_idx
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
%   --  H: Non-negative matrix such that W = H*H^T

% percent_features = 0.5;

if nargin==1
    NC = [];
end
[No_genes,No_cells] = size(data);


% Data preprocess
% 1st filtering genes
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
display(No_Comps);

aa = max(coeff(:,1:No_Comps+1)');
bb = sort(aa,'descend');

if size(data1,1) <=1000
    No_sel_genes = size(data1,1);
else
    No_sel_genes = min([No_features round(size(data1,1))]);
end
% No_sel_genes = round(1*size(data1,1));

display(No_sel_genes);

gene_selection = find(aa>=bb(No_sel_genes));

input_data = data1(gene_selection,:);


Gene_sel_idx = G_filter_idx(gene_selection);
display(length(Gene_sel_idx));
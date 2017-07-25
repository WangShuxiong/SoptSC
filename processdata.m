function newdata = processdata(data,alpha,beta)
% This function preprocesses single cell data
% by selecting a subset of genes based on
%   1) the variance of each gene, tolerance is alpha 
%   2) the number of cells that a particular gene is expressed
%      beta is the ratio where 0<= beta <=1
%
%

[m,n] = size(data);
V_gene1 = var(data');
gene_sel_idx1 = find(V_gene1>alpha);

Z = zeros(1,m);
for i = 1:m
    Z(i) = length(find(data(i,:)>1));
end
gene_sel_idx2 = find(Z>n*beta);
gene_sel_idx = intersect(gene_sel_idx1,gene_sel_idx2);
newdata = data(gene_sel_idx,:);

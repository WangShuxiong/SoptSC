%% Setup paths and load data
clear;
clc;
echo on;
addpath('Data');
addpath('NNDSVD');
addpath('symnmf2');

% Select and load data matrix in format: Genes x Cells 
load HEE_matrix.mat;
data = HEE_matrix';

% Load all gene names as string
load allgenes.mat;
% allgenes = allgenes;
% load Test_3_Pollen.mat;
% data = in_X';
%% Optional step: preprocess data by selecting a subset of genes
alpha = 0.5; % Variance in gene expression (threshold)
beta = 0.5;  % Number of cells in which a gene is expressed
             % (threshold ratio)

[select_genes,data_processed] = processdata(data,alpha,beta);


%% Step 1: Run SoptSC to identify clusters and subpopulation composition
NC = [];    % NC is the number of clusters: can be specified by user, or
            % if not given (NC = []), it will be inferred

% [W,P,No_cluster,cluster_label,latent] = SOptSC_cluster(input',NC);

[W,P,No_cluster,cluster_label,latent,H] = SOptSC_cluster(data_processed,NC);
% Output
%   --  W: Cell-to-cell similarity matrix.
%   --  P: Transition matrix.
%   --  No_cluster: Number of clusters.
%   --  cluster_label: cluster labels for all cells.
%   --  H: decomposed matrix for W.
%   --  latent: low dimensional space (first three eigenvectors) of
%   transition matrix P.

%% Step 2: Plot gene expression to compare with clusters
gene_set = {'SALL2','FGFR4'};   % Marker genes selected by user
plot_genes(gene_set,allgenes,data,latent)

%% Step 3: Plot gene-cell heatmap

topn = 20;  % Number of top genes selected for each cluster to plot.
Gene_labels = GC_heatmapTopn(data_processed,cluster_label,H,allgenes,select_genes,topn);
%% Step 4: Run SoptSC to infer pseudotime and lineage tree
init_cluster = 2;               % Starting cluster specified by user
                                % based on the analysis from step 2
cell_order = SOptSC_pseudotime(init_cluster,P,cluster_label,latent,No_cluster);


%% Setup paths and load data
clear;
clc;
echo on;
addpath('Data');
addpath('NNDSVD');
addpath('symnmf2');

% Load data
load Guo2010.mat;

%% Step 1: Run SoptSC to identify clusters and subpopulation composition
NC = [];    
No_exc_cell = 0;
No_features = 2000;
[W,No_cluster,cluster_label,latent,H,Gene_sel_idx] = SOptSC_cluster(data,NC,No_exc_cell,No_features);

%% Pseudotime and lineage inference
root_cluster = 0;
root_cell = 0;
reverse = 0;
[Lineage, Ptime] = Lineage_Ptime(W,No_cluster,cluster_label,root_cluster,root_cell,latent,reverse);



%% Plot gene expression to compare with clusters
gene_set = {'Nanog','Cdx2'};   % Marker genes selected by user
plot_genes(gene_set,allgenes,data,latent)
%% Setup paths and load data
clear;
clc;
echo on;

addpath('Data');
addpath('NNDSVD');
addpath('symnmf2');

% REQUIRED INPUTS: 'data' and 'allgenes (gene names)'
% Load data
load Guo2010.mat;

% To save results 
ResFolder = strcat('Results_',string(datetime('now','Format','yyyyMMdd_HHmmss')));
mkdir(ResFolder);

%% Step 1: Run SoptSC to identify clusters and subpopulation composition
NC = [];    
No_exc_cell = 0;
No_features = 2000;
[W,No_cluster,cluster_label,latent,H,Gene_sel_idx] = SOptSC_cluster(data,NC,No_exc_cell,No_features,ResFolder);

%% Pseudotime and lineage inference
root_cluster = 0;
root_cell = 0;
reverse = 0;
[Lineage, Ptime] = Lineage_Ptime(W,No_cluster,cluster_label,root_cluster,root_cell,latent,reverse,ResFolder);


%% Plot gene expression to compare with clusters
gene_set = {'Nanog','Cdx2'};   % Marker genes selected by user
plot_genes(gene_set,allgenes,data,latent)
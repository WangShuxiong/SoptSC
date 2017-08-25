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

%% Optional step: preprocess data by selecting a subset of genes
alpha = 0.5; % Variance in gene expression (threshold)
beta = 0.5;  % Number of cells in which a gene is expressed
             % (threshold ratio)
data = processdata(data,alpha,beta);

%% Run SOptSC
NC = [];
init_point = 1;
[W,P,No_cluster,cluster_label,cell_order] = SOptSC(data,init_point,NC);

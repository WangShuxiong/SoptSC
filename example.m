clear;
clc;
addpath('Data');
addpath('NNDSVD');
addpath('symnmf2');
% data = importdata('data.txt');
load HEE_matrix.mat;
data = HEE_matrix';

%% Optional step: preprocess data by selecting a subset of genes
alpha = 0.5;
beta = 0.5;
data = processdata(data,alpha,beta);

%% Run SOptSC
NC = 3;
init_point = 1;
[W,P,No_cluster,cluster_label,cell_order] = SOptSC(data,init_point,NC);

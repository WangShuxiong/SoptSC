% add path
clc;
addpath('Data');
addpath('package/symnmf');
addpath('package');
%% Read 10X data
data_dir = 'Data/GRCh38'; % path where 10x data located
[count, gene_annotation, gene_ids, cell_annotation] = load_10x(data_dir,'sparse',1);

%% load the processed data
load Data/SoptSC_object.mat;
%%
opt.alpha = 10;     % remove cells that express less than opt.alpha genes (default: 10)
opt.beta = 0.1;     % remove cells that express more than opt.beta*100% genes (default: 0.1)
SoptSC_object = creat_soptsc_object(count,gene_annotation,gene_ids,cell_annotation,opt);

%% save all the figures in to folder
folder = 'Results';

%% Step 1: Run SoptSC to identify clusters and subpopulation composition
opt.No_cluster = [];    
opt.No_exc_cell = 10;
opt.No_features = 2000;
% opt.K = 30;
tic;
SoptSC_object = run_SoptSC(SoptSC_object,opt); % run_SoptSC
time = toc
%% Change the number of cluster to be 10
No_cluster = 7;
SoptSC_object = SoptSC_change_No_cluster(SoptSC_object,No_cluster);

%% Run dimension reduction
opt.ppx = 50;   % ppx: peplexicty for umap
opt.topn = 300;
opt.input = 0; % 0: used topn gene expression as input; 1: use similarity as input
SoptSC_object = run_dimensionreduction(SoptSC_object,opt);

%%
plot_eigengap(SoptSC_object,folder);

%% plot heatmap of topn DEGs
opt1.figsize = [0 0 6 10];
opt1.figname = 'GC_heatmap';
topn = 10;
plot_DEGs_heatmap(SoptSC_object,topn,folder,opt1);

%% user-defined cluster colors: No_cluster x 3 matrix (each row represent rgb)
% zzz = hsv;
% cluster_color = zzz(1:round(256./No_cluster):end,:);
cluster_color = get(gca,'colororder');
%% plot cluster
opt.figname = 'Cluster';
opt.cluster_color = cluster_color;
opt.figsize = [0 0 5 3]; % figure size
plot_cluster(SoptSC_object,folder,opt);

%%
opt.figname = 'Cluster_Cell_Percantage';
opt.figsize = [0 0 4 1.5]; % figure size
plot_cell_percent(SoptSC_object,folder,opt);

%%
Marker = {'KRT14','KRT10','PTTG1','RRM2','GJB2','MLANA','ASS1','KRT6A'};

%% plot mean expression for selected genes
opt.threshold = 0.1;
opt.figname = 'dot_gene_map';
opt.figsize = [0 0 4 2];
plot_dot_mean_genes(SoptSC_object,Marker,folder,opt);

%% Violin
opt.figsize = [0 0 4 1.5];
plot_violin(SoptSC_object,Marker,folder,opt);

%% feature plot
opt.figsize = [0 0 4 3];
plot_marker(SoptSC_object,Marker,folder,opt);

%% user defined cluster name
SoptSC_object.cluster_name = {'BAS-II','BAS-IV','SPN-II','BAS-III','GRN','SPN-I','BAS-I'};

%% plot cluster with new cluster names
opt.figname = 'Cluster_named';
opt.cluster_color = cluster_color;
opt.figsize = [0 0 5 3]; % figure size
plot_cluster(SoptSC_object,folder,opt);

%% run lineage and pseudotime inference
opt.root_cluster = 7; 
opt.root_cell = 0;
opt.reverse = 0;
SoptSC_object = run_lineage_and_pseudotime(SoptSC_object,opt);

%% plot pseudotime for cells in UMAP space
opt.figname = 'Pseudotime';
opt.figsize = [0 0 4 3];
plot_pseudotime(SoptSC_object,folder,opt);

%% plot lineage tree
opt.figname = 'Lineage';
opt.figsize = [0 0 3 5];
plot_lineage(SoptSC_object,folder,opt);

%% plot marker on lineage treee
opt.figsize = [0 0 3 5];
plot_lineage_marker(SoptSC_object,Marker,folder,opt);


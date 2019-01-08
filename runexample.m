% This is an example for applying SoptSC to scRNA-seq data (Joost et al. 2016) to perform
%   - Clustering
%   - Marker genes identified for each cluster
%   - Cell lineage
%   - Pseudotiime
%   - Cell-to-cell signaling network for given ligand-receptor pairs and
%     downstream target genes
%   - Cluster-to-cluster signaling network for given ligand-receptor pairs
%     and downstream target genes
%
%   All details of the function can be found in the corresponding function
%   files. Plsease refer to each M file of functions for the discription of
%   input and output.
%
%   Results are save in the folder: Results
%
%   Contact: Shuxiong Wang (Email: shuxionw@uci.edu) if you have any
%   question.

clear;
clc;
echo on;
addpath('Data');
addpath('NNDSVD');
addpath('symnmf2');
addpath('Signaling');
addpath('vinlinplot');
addpath('Results')

Data_all = importdata('JoostData.txt');
data_matrix = Data_all.data;
allgenes = Data_all.textdata(2:end,1);


%% Log normalization
data = log10(data_matrix +1);

%% Step 1: Run SoptSC to identify clusters and subpopulation composition
resfolder = 'Results';     % all results are saved here
NC = [];    % NC is the number of clusters: can be specified by user, or
            % if not given (NC = []), it will be inferred
No_cells = size(data,2);
No_exc_cell = 0.03*No_cells;
No_features = 3000;

[W,No_cluster,cluster_label,H,eigenvalues] = SoptSC_cluster(data,NC,No_exc_cell,No_features,resfolder);

%% Plot Eigen-gap of truncated graph Laplacian of the consensus matrix (if NC = [])
plot_eigengap(eigenvalues,resfolder);


%% Plot cluster on 2-dimensional space
method = 'tsne';        % set method as 'pca' or 'tsne'
latent = plot_cluster(W,cluster_label,No_cluster,method,resfolder);

%% Gene-cell heatmap for all genes
No_exc_cell1 =6;
No_select_genes = 20000;

Gene_labels_all = GC_heatmap(data,cluster_label,H,No_exc_cell1,No_select_genes);

%% Gene-cell heatmap for top n markers w.r.t each cluster
topn = 10;
Gene_labels_topn = GC_heatmapTopn(data,cluster_label,H,allgenes,No_exc_cell1,No_select_genes,topn,resfolder);



%% Plot gene expression on the low-dimensional projection of cells
Marker = {'Krt14','Mt2','Krt10','Ptgs1','Lor','Flg2'};
plot_marker(data,Marker,allgenes,latent,resfolder) 


%% Violin plot of marker genes along clusters
plot_marker_violin(data,allgenes,Marker,cluster_label,No_cluster,resfolder)


%% Bar plot of marker genes along clusters
boxplot_marker(data,allgenes,Marker,cluster_label,No_cluster,resfolder);


%% Pseudotime and lineage inference

root_cluster = 0;
root_cell = 0;
reverse = 1;
[Lineage, Ptime,Cell_dist] = Lineage_Ptime(W,No_cluster,cluster_label,root_cluster,root_cell,latent,reverse);


%% Plot cluster color and pseudotime color on the lineage tree
plot_lineage(Lineage,No_cluster,cluster_label,Cell_dist,resfolder);


%% Plot pseudotime on latent space
plot_pseudotime(latent,Ptime,resfolder);

%% Plot marker gene on the lineage tree

plot_lineage_marker(data,Lineage,allgenes,No_cluster,cluster_label,Marker,resfolder)




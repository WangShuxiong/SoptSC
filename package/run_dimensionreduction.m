function obj = run_dimensionreduction(obj,opt)

W = obj.W;
% data_degs = data1(DE_idx,:);
% ppx = 30;
% dim_init = 2;
% [umap, ~, cluster_lab_umap]=run_umap(data_degs','n_neighbors',ppx,'n_components',dim_init,...
%     'verbose','text','method','Java','metric','correlation'); % row -- data point; column -- feature; correlation

% UMAP based on W
ppx = 30;
dim_init = 2;
if isfield(opt,'ppx')
    ppx = opt.ppx;
end
if isfield(opt,'dim_init')
    dim_init = opt.dim_init;
end

[umap,~,~]=run_umap(W,'n_neighbors',ppx,'n_components',dim_init,...
    'verbose','text','method','Java','metric','correlation'); % row -- data point; column -- feature

obj.umap = umap;


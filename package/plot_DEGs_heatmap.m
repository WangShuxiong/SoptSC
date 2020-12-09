function gene_idxv0 = plot_DEGs_heatmap(obj,topn,folder,folder1,opt)
data = obj.data;
allgenes = obj.gene_annotation;
cluster_labs = obj.cluster_label;
H = obj.H;
lgd = obj.cluster_name;

H = H./max(H(:));

zz = mean(data,2);
idx0 = [1:length(allgenes)]';
idx = find(zz==0);
data(idx,:) = [];
allgenes(idx) = [];
idx0(idx) = [];
No_clusterr = length(unique(cluster_labs));

[No_gene,~] = size(data);
gene_mean = zeros(No_gene,No_clusterr);
gene_DE_score = zeros(No_gene,No_clusterr);


[gene_mean,gene_DE_score] = corr(data',H);

[~,gene_value_idx] = max(gene_mean,[],2);

% topn markers for each cluster based on DE score
gclusters = [];
gscore = [];

gene_idxv = [];
cluster_order = [];

for i = 1:No_clusterr
    zz_idx = find(gene_value_idx == i);
    zz_DEscore = gene_DE_score(zz_idx,i);
    [zzvalue,zz1] = sort(zz_DEscore);
    topn1 = min([topn length(zz_idx)]);
    gene_idxv = [gene_idxv; zz_idx(zz1(1:topn1))];
    gclusters = [gclusters;i.*ones(topn1,1)];
    gscore = [gscore;zzvalue(1:topn1)];
    cluster_order = [cluster_order;find(cluster_labs==i)];
end



GL500 = [idx0(gene_idxv) gclusters gscore];
T = array2table(GL500,'RowNames',allgenes(gene_idxv),'VariableNames',{'Gene_indx','Cluster','P_val'});
writetable(T,[folder1 '/' opt.figname '-' num2str(topn) '.csv'],'WriteRowNames',true,'WriteVariableNames',true,'Delimiter',',');  

gene_idxv0 = idx0(gene_idxv);

datav = data(gene_idxv,cluster_order);
idata = datav;
kk = 2;
center = mean(idata,kk);
scale = std(idata, 0,kk);
tscale = scale;
%=Check for zeros and set them to 1 so not to scale them.
scale(tscale == 0) = 1;
%== Center and scale the data
idata = bsxfun(@minus, idata, center);
sdata = bsxfun(@rdivide, idata, scale);
thresh = 3;

colormap redbluecmap;
clims = [-thresh thresh];
imagesc(sdata,clims);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

if size(datav,1) < 200
yticks(1:size(datav,1));
yticklabels(allgenes(gene_idxv));
end


No_cells_inC = [];
for i = 1:No_clusterr
    No_cells_inC = [No_cells_inC; length(find(cluster_labs==i))];
end
xtkval = cumsum(No_cells_inC);
xtkval1 = zeros(size(xtkval));
for i = 1:No_clusterr
    if i==1
        xtkval1(i) = 0.5.*No_cells_inC(i);
    else
        xtkval1(i) = xtkval(i-1) + 0.5.*No_cells_inC(i);
    end
end

hold on;
% draw lines between clusters
yy = length(gene_idxv);
for i = 1:No_clusterr-1
     line([xtkval(i) xtkval(i)],[0 yy+0.5],'color','w',...
        'LineStyle','-','LineWidth',2);
    hold on;
end

xticks(xtkval1);
xticklabels(lgd);
xtickangle(45);
cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = opt.figsize;
box off;

print([folder '\' opt.figname '_' num2str(topn)],'-depsc','-r300'); %'-dpdf', ,'-fillpage'


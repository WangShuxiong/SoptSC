function plot_dot_mean_genes(obj,mk,folder,opt)
data = obj.data;
allgenes = obj.gene_annotation;
lgd = obj.cluster_name;
cluster_label = obj.cluster_label;
No_cluster = length(unique(cluster_label));
reorder = 1:No_cluster;
figname = opt.figname;

if isfield(opt,'reorder')
    reorder = opt.reorder;
end


figsize =  [0 0 21 8];
if isfield(opt,'figsize')
    figsize = opt.figsize;
end

threshold = 0.1;
if isfield(opt,'figsize')
    threshold = opt.threshold;
end


No_LR = length(mk);
No_cluster = length(lgd);

aa = [1:No_LR]';
bb = repmat(aa,No_cluster,1);
x = bb(:);

y = zeros(size(x));
for i = 1:No_cluster
    y((i-1)*No_LR + 1:i*No_LR) = No_cluster + 1 - i;
end

[~,~,gene_idx] = intersect(mk,allgenes,'stable');
mk = allgenes(gene_idx);

CC_single_M = zeros(length(mk),No_cluster);
CC_single_M1 = zeros(size(CC_single_M));   % dot size

data1 = data(gene_idx,:);

for i = 1:No_cluster
    CC_single_M(:,i) = mean(data1(:,cluster_label==i),2);
    for j = 1:No_LR
        CC_single_M1(j,i) = nnz(data1(j,cluster_label==i))./nnz(data1(j,:));
    end
end
CC_single_M = CC_single_M./max(CC_single_M,[],2);


% dot size
zz = CC_single_M1;
dot_size = zz(:);
dot_size(dot_size <= threshold*max(zz(:))) = 0;
dot_size = 80.*dot_size./max(zz(:)) +1;


CC_single_M = CC_single_M./max(CC_single_M(:));
% dot color scale
zz = CC_single_M;
dot_color = zz(:);

% colormap(redbluecmap);
% colormap(redbluemp_color(30));

a = redbluecmap;
a(6,:) = [];
colormap(a(2:9,:));

% colormap(cool);
scatter(x,y,dot_size,dot_color,'filled');
ylim([0 No_cluster+1]);
xlim([0 No_LR+1]);


% ynames = cell(No_cluster,1);
% for i = 1:No_cluster
%     ynames{i} = [mk{No_cluster-i+1}];
% end
ynames = flip(lgd(reorder));


yticks(1:length(ynames));
yticklabels(ynames);

xticks(1:length(mk));
xticklabels(mk);
xtickangle(90);

% a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',10);
% 
% a = get(gca,'YTickLabel');  
% set(gca,'YTickLabel',a,'fontsize',10);
set(gca,'FontName','Arial');
set(gca,'FontSize',10);
 

ax = gca;
ax.TickDir = 'out';
% ax.LineWidth = 1;
box off;
% box on;
grid on;

% cb = colorbar('Ticks',tickval,'TickLabels',lgd,'FontSize',12);
cb = colorbar;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.3*cpos(3);
cb.Position = cpos;
ax.Position = axpos;

% Set the size of output fig
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 4 3];

fig.Units = 'Inches';
fig.Position =figsize; % 8 6 for mock edge;
                          % 6 6 dnlef
print([folder '\' figname '_' num2str(threshold*10)],'-depsc','-r300');



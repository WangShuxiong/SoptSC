function plot_violin(obj,Marker,folder,opt)
data = full(obj.data);
allgenes = obj.gene_annotation;
cluster_label = obj.cluster_label;
No_cluster = length(unique(cluster_label));
lgd = obj.cluster_name;

zzz = hsv;
cluster_color = zzz(1:round(255./No_cluster):end,:);
if isfield(opt,'cluster_color')
    cluster_color = opt.cluster_color;
end
mycolor = cluster_color;

% reorder = 1:No_cluster;
% if isfield(opt1,'reorder')
%     reorder = opt1.reorder;
% end

figuresize = [0 0 4 2];
if isfield(opt,'figsize')
    figuresize = opt.figsize;
end


for i = 1:length(Marker)
    [~,ia,~] = intersect(allgenes,Marker{i},'stable');
    figure(i)
    h = violinplot(data(ia,:),cluster_label',1:No_cluster);
    for j = 1:length(h)
        h(j).ViolinColor = mycolor(j,:);
    end
    set(gca,'Xtick',1:No_cluster)
    set(gca,'Xticklabel',lgd,'FontSize',12);
    title(Marker{i});
    
    set(gca,'FontName','Arial');
    set(gca,'FontSize',12);
    ax = gca;
    ax.TickDir = 'out';
    ax.LineWidth = 1.5;
    
    % Set the size of output fig
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    fig.Units = 'Inches';
    fig.Position = figuresize;
    box off;
    
    print([folder '\violin_' Marker{i}],'-depsc','-r300'); %'-dpdf',
end

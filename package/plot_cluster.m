function plot_cluster(obj,folder,opt1)
latent = obj.umap;
cluster_label = obj.cluster_label;
No_cluster = length(unique(cluster_label));
lgd = obj.cluster_name;
reorder = 1:No_cluster;
figname = opt1.figname;

zzz = hsv;
cluster_color = zzz(1:round(256./No_cluster):end,:);
if isfield(opt1,'cluster_color')
    cluster_color = opt1.cluster_color;
end

if isfield(opt1,'reorder')
    reorder = opt1.reorder;
end


for ik = 1:No_cluster
    scatter(latent(cluster_label==reorder(ik),1),latent(cluster_label==reorder(ik),2),10,cluster_color(ik,:),'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    hold on;
end
box off;
xlabel('UMAP1');
ylabel('UMAP2');
legend(lgd(reorder),'FontSize',12,'Location','eastoutside');%,'Orientation','horizontal');
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = opt1.figsize;
box off;

print([folder '\' figname],'-depsc','-r300'); %'-dpdf','-fillpage'

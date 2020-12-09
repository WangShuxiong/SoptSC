function plot_lineage(obj,folder,opt)
Lineage = obj.lineage;
cluster_label = obj.cluster_label;
Cell_dist = obj.pseudotime;
CC_adjacent = obj.CC_adj;
lgd = obj.cluster_name;
No_cluster = length(unique(cluster_label));


figname = 'Lineage';
if isfield(opt,'figname')
    figname = opt.figname;
end

figuresize = [0 0 4 3];
if isfield(opt,'figsize')
    figuresize = opt.figsize;
end

zzz = hsv;
cluster_color = zzz(1:round(255./No_cluster):end,:);
if isfield(opt,'cluster_color')
    cluster_color = opt.cluster_color;
end
mycolor = cluster_color;

Nodesize = zeros(No_cluster,1);
for i = 1:No_cluster
    Nodesize(i) = length(find(cluster_label==i));
end
Nodesize = 30*Nodesize./max(Nodesize) + 3;


WM = CC_adjacent + CC_adjacent';
pred = Lineage;
aa = pred(pred~=0);
bb = find(pred~=0);
Tw = zeros(length(aa),1);
for i = 1:length(aa)
    Tw(i) = WM(aa(i),bb(i));
end
rootedTree = digraph(pred(pred~=0),find(pred~=0),Tw);
Lwidth = 5*rootedTree.Edges.Weight./max(rootedTree.Edges.Weight);



pred = Lineage;
rootedTree = digraph(pred(pred~=0),find(pred~=0));

plot(rootedTree,'Marker','o','MarkerSize',Nodesize,'NodeColor',mycolor(1:No_cluster,:),...
    'LineWidth',Lwidth,'NodeLabel',lgd,'EdgeColor',[0.69 0.77 0.87]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
% title('Cluster on lineage')
box off;
axis off;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = figuresize; % [0 0 2 4];

print([folder '\' figname '-cluster'],'-depsc','-r300'); 
% print([folder '\Lineage_Cluster_Color_' figname],'-dpdf','-r300'); 



cmap = parula;
mymap = cmap(1:round(256*58./64),:);

ptimecolor = zeros(No_cluster,1);

for i = 1:No_cluster
     ptimecolor(i) = mean(Cell_dist(cluster_label==i));
%     ptimecolor(i) = mean(Ptime(find(cluster_label==i)));
end

pred = Lineage;
rootedTree = digraph(pred(pred~=0),find(pred~=0));
figure;
colormap(mymap);
plot(rootedTree,'Marker','o','MarkerSize',Nodesize,'NodeCData',ptimecolor, 'NodeColor','flat','NodeLabel',lgd,...
    'LineWidth',Lwidth,'EdgeColor',[0.69 0.77 0.87]);

% cb = colorbar;
% ax = gca;
% axpos = ax.Position;
% cpos = cb.Position;
% cpos(3) = 0.5*cpos(3);
% cb.Position = cpos;
% ax.Position = axpos;
% % lim = caxis
% % cb.Limits = lim;
% aa = cell(1,2);
% aa{1} = 'low';
% aa{2} = 'high';
% cb.TickLabels{1} = aa{1};
% cb.TickLabels{end} = aa{2};
% 
% for ii = 2:length(cb.TickLabels)-1
%     cb.TickLabels{ii} = [];
% end
% box off;

axis off;
% set(gca,'LineWidth',1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'FontName','Arial');
set(gca,'FontSize',12);


fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

fig.Units = 'Inches';
fig.Position = figuresize; %[0 0 2 4];
box off;

print([folder '\' figname '-ptime'],'-depsc','-r300'); 

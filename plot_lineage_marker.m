function plot_lineage_marker(data,Lineage,allgenes,No_cluster,cluster_label,marker)
% Plot marker on lineage

mk = marker;
[~,ia] = intersect(allgenes,mk);
mk_val = data(ia,:);

mkcolor = zeros(No_cluster,1);

for i = 1:No_cluster
     mkcolor(i) = mean(mk_val(cluster_label==i));
%     ptimecolor(i) = mean(Ptime(find(cluster_label==i)));
end

pred = Lineage;
rootedTree = digraph(pred(pred~=0),find(pred~=0));

load newmap.mat;
figure;
mymap = newmap;
colormap(mymap);
% plot(rootedTree,'Marker','o','MarkerSize',20,'NodeCData',ptimecolor, 'NodeColor','flat','NodeLabel',[]);
plot(rootedTree,'Marker','o','MarkerSize',20,'NodeCData',mkcolor, 'NodeColor','flat','NodeLabel',[]);

set(gca,'xtick',[]);
set(gca,'ytick',[]);

cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;
% lim = caxis
% cb.Limits = lim;
aa = cell(1,2);
aa{1} = 'low';
aa{2} = 'high';
cb.TickLabels{1} = aa{1};
cb.TickLabels{end} = aa{2};

for ii = 2:length(cb.TickLabels)-1
    cb.TickLabels{ii} = [];
end
box on;

% set(gca,'LineWidth',1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'FontName','Arial');
set(gca,'FontSize',12);


title(mk)
print(['Results\Lineage_Marker_' marker],'-dpdf','-r300'); 

function plot_cell_percent(obj,folder,opt)
cluster_labs = obj.cluster_label;
lgd = obj.cluster_name;
reorder = 1:length(lgd);
figname = opt.figname;

No_cluster = length(unique(cluster_labs));

zzz = hsv;
cluster_color = zzz(1:round(255./No_cluster):end,:);
if isfield(opt,'cluster_color')
    cluster_color = opt.cluster_color;
end
mycolor = cluster_color;

if isfield(opt,'reorder')
    reorder = opt.reorder;
end


No_cells = [];
for i = 1:No_cluster
    No_cells = [No_cells; length(find(cluster_labs==i))];
end



figuresize = [0 0 2.5 1.2];
if isfield(opt,'figsize')
    figuresize = opt.figsize;
end

% 
% b = bar(No_cells(reorder));
% b.FaceColor = 'flat';
% 
% for i = 1:No_cluster
%     b.CData(i,:) = mycolor(i,:);
% end
% % box off;
% grid on;
% box off;
% 
% xtips1 = b(1).XEndPoints;
% ytips1 = b(1).YEndPoints;
% labels2 = string(lgd(1:No_cluster));
% 
% for i = 1:No_cluster
%     labels2(i) = num2str(No_cells(i));
% end
% 
% text(xtips1,ytips1,labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% 
% 
% 
% % title('Number of Cells');
% xticks(1:No_cluster);
% xticklabels(lgd(reorder));
% xtickangle(45);
% ax = gca;
% ax.TickDir = 'out';
% ax.LineWidth = 1.5;
% 
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% 
% % fig.PaperUnits = 'inches';
% % fig.PaperPosition = [0 0 4 3];
% 
% fig.Units = 'Inches';
% fig.Position = figuresize;
% print([folder '\No_Cells_' figname],'-depsc','-r300'); %'-dpdf',,'-fillpage'

figure;
b = bar(No_cells(reorder)./sum(No_cells));
b.FaceColor = 'flat';

% labels1 = string(b(1).YData);

zz = No_cells(reorder)./sum(No_cells);

xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels2 = string(lgd);
for i = 1:length(lgd)
    labels2(i) = [num2str(round(zz(i)*100)) '%'];
end

text(xtips1,ytips1,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')


for i = 1:No_cluster
    b.CData(i,:) = mycolor(i,:);
end
grid on;
box off;
% title('Percentage of Cells');
xticks(1:No_cluster);
xticklabels(lgd(reorder));
xtickangle(45);

ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 1.5;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
fig.Units = 'Inches';
fig.Position = figuresize;


% print([folder '\No_Cells_Percent_' figname],'-dpdf','-r300'); %'-dpdf',,'-fillpage'
print([folder '\' figname],'-depsc','-r300'); %'-dpdf',,'-fillpage'


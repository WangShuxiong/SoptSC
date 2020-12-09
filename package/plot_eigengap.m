function plot_eigengap(obj,folder)

eigenvalues = obj.eigenvalues;

scatter(1:min([30 size(eigenvalues,1)]),eigenvalues(1:min([30 size(eigenvalues,1)])),20,'filled');
grid on;
box off;
yticks(0:0.2:1);
xticks(1:3:30)
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
fig.Position = [0 0 4 2];

print([folder '\EigenGap'],'-depsc','-r300');


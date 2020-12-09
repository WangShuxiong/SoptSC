function plot_pseudotime(obj,folder,opt)
latent = obj.umap;
Cell_dist = obj.pseudotime;

figname = 'Pseudotime';
if isfield(opt,'figname')
    figname = opt.figname;
end

scatter(latent(:,1),latent(:,2),15,Cell_dist,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
figuresize = [0 0 5.5 3];
if isfield(opt,'figsize')
    figuresize = opt.figsize;
end


cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;
% box on;
aa = cell(1,2);
aa{1} = 'low';
aa{2} = 'high';
cb.TickLabels{1} = aa{1};
cb.TickLabels{end} = aa{2};

for ii = 2:length(cb.TickLabels)-1
    cb.TickLabels{ii} = [];
end

xlabel('UMAP1');
ylabel('UMAP2');
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
box off;
% axis off;

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

print([folder '\' figname],'-depsc','-r300'); 

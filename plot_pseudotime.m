function plot_pseudotime(latent,Ptime,resfolder)
%% plot pseudotime on the latent space
No_cell = size(latent,1);
c = linspace(0,1,No_cell);
figure;
colormap parula;
cmap = colormap;
mymap = cmap(1:58,:);
colormap(mymap);
scatter(latent(:,1),latent(:,2),30,Ptime./max(Ptime(:)),'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);


cb = colorbar;
ax = gca;
axpos = ax.Position;
cpos = cb.Position;
cpos(3) = 0.5*cpos(3);
cb.Position = cpos;
ax.Position = axpos;
box on;

set(gca,'LineWidth',1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'FontName','Arial');
set(gca,'FontSize',12);

print([resfolder '\Pseudotime'],'-dpdf','-r300'); 

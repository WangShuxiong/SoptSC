function plot_marker(obj,gene_set,folder,opt)
data = obj.data;
allgenes = obj.gene_annotation;
latent = obj.umap;

[~,~,gene_set_idx] = intersect(gene_set,allgenes,'stable');
gene_set_no = length(gene_set_idx);

figuresize = [0 0 4 3];
if isfield(opt,'figsize')
    figuresize = opt.figsize;
end


MM = data(gene_set_idx,:);

%% Marker genes expression on each subpopulation

% aa1 = hot;
% bb = bone;
% zz = bb(176:end,:);
% mymap = [zz;flip(aa1(96:2:end,:))];
% 

for ik = 1:gene_set_no
    figure(ik);
    colormap redbluecmap; % jet,
    scatter(latent(:,1),latent(:,2),5,MM(ik,:),'filled','MarkerEdgeAlpha',0.8,'MarkerFaceAlpha',0.8);
    
    %    box on;
    box off;
    set(gca,'FontName','Arial');
    set(gca,'FontSize',12);
    xlabel('UMAP1');
    ylabel('UMAP2');
    ax = gca;
    ax.TickDir = 'out';
    ax.LineWidth = 1.5;
    title(gene_set{ik});
    
    ax = gca;
    cb = colorbar;
    axpos = ax.Position;
    cpos = cb.Position;
    cpos(3) = 0.5*cpos(3);
    cb.Position = cpos;
    ax.Position = axpos;
    %     cb.TickLabels{1} = 0;
    %     cb.TickLabels{end} = 1;
    %
    %     for ii = 2:length(cb.TickLabels)-1
    %         cb.TickLabels{ii} = [];
    %     end
    %      lim = caxis
    %     cb.Limits = [0 1];
    %     aa = cell(1,2);
    %     aa{1} = 'low';
    %     aa{2} = 'high';
    % Set the size of output fig
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
        
    fig.Units = 'Inches';
    fig.Position = figuresize;
    box off;
    print([folder '\FeaturePlot-' gene_set{ik}],'-depsc','-r300'); %'-dpdf',
end



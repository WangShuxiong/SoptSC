function plot_lineage_marker(obj,marker,folder,opt)
Lineage = obj.lineage;
cluster_label = obj.cluster_label;
CC_adjacent = obj.CC_adj;
lgd = obj.cluster_name;
No_cluster = length(unique(cluster_label));

data = obj.data;
allgenes = obj.gene_annotation;

figuresize = [0 0 4 3];
if isfield(opt,'figsize')
    figuresize = opt.figsize;
end


Nodesize = zeros(No_cluster,1);
for i = 1:No_cluster
    Nodesize(i) = length(find(cluster_label==i));
end
Nodesize = 20*Nodesize./max(Nodesize);

marker = intersect(allgenes,marker);
No_genes = length(marker);

% plot cluster color on lineage tree
WM = CC_adjacent + CC_adjacent';
pred = Lineage;
aa = pred(pred~=0);
bb = find(pred~=0);
Tw = zeros(length(aa),1);
for i = 1:length(aa)
    Tw(i) = WM(aa(i),bb(i));
end
rootedTree = digraph(pred(pred~=0),find(pred~=0),Tw);
Lwidth = 4*rootedTree.Edges.Weight./max(rootedTree.Edges.Weight);

for i = 1:No_genes
    [~,~,gene_idx] = intersect(marker{i},allgenes,'stable');
    mycolor = zeros(No_cluster,1);
    for j = 1:No_cluster
        mycolor(j) = mean(data(gene_idx,cluster_label==j));
    end
    figure;
    colormap redbluecmap;
    % figure;
    % plot(rootedTree,'Marker','o','MarkerSize',20,'NodeColor',mycolor(1:No_cluster,:),'NodeLabel',[]);
    plot(rootedTree,'Marker','o','MarkerSize',Nodesize,'NodeCData',mycolor,'NodeColor','flat',...
        'LineWidth',Lwidth,'ArrowSize',12,'EdgeColor',[0.69 0.77 0.87],...
        'NodeLabel',lgd); %,'WeightEffect','direct'? ,'layout','force'
    
    cb = colorbar;
    ax = gca;
    axpos = ax.Position;
    cpos = cb.Position;
    cpos(3) = 0.5*cpos(3);
    cb.Position = cpos;
    ax.Position = axpos;
    
    box off;
    axis off;
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    title(marker{i})
    set(gca,'FontName','Arial');
    set(gca,'FontSize',12);
    
    % Set the size of output fig
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    
    fig.Units = 'Inches';
    fig.Position = figuresize;
    
    print([folder '\L' '_' marker{i}],'-depsc','-r300');
end
function plot_sig_network(Pidv,Pall,cluster_label,Lig,Rec,threshold,folder)
% plot cell-cell and cluster-cluster signaling network based on the probability matrix
% Input:
%   -- Pidv: Cell-to-cell interaction probability of individual
%      ligand-receptor pair where Pidv{i} is the probability corresponds to
%      the ith ligand-receptor pair.
%   -- Pall: cell-to-cell interaction probability based on all
%      ligand-receptor pairs and their target genes.
%   -- threshold: restricting probabiltiy between cells less than threshold
%      to be zero.
%   -- folder: folder name where the results will be saved to.
%
%   Output:
%   -- cell-cell signaling network for individual ligand-receptor pair and
%   all ligand-receptor pairs.
No_LR = length(Lig);


max_prob_idv = 0;
for ii = 1:No_LR
    P1 = Pidv{ii};
    if max_prob_idv <= max(P1(:))
        max_prob_idv = max(P1(:));
    end
end

if max_prob_idv <= max(Pall(:))
    max_prob_idv = max(Pall(:));
end

if threshold >= max_prob_idv
    display('Error: Threshold Is Too Large !!!');
    return;
end
    


No_cluster = length(unique(cluster_label));
No_cells = length(cluster_label);

cmap1 = jet;
mymap1 = cmap1(1:end,:);
ncolor = size(mymap1,1);
mycolor = mymap1(1:round(ncolor./No_cluster):ncolor,:);

mycolor_cells = zeros(No_cells,3);



%% Order cells based on cluster labels
true_labs = [];
new_order = [];
No_cells_cluster = [];

for i = 1:No_cluster
    true_labs = [true_labs; cluster_label(cluster_label ==i)];
    new_order = [new_order; find(cluster_label ==i)];
    No_cells_cluster = [No_cells_cluster; length(find(cluster_label ==i))];
end

for ii = 1:No_cells
    mycolor_cells(ii,:) = mycolor(true_labs(ii),:);
end

%% plot cell-cell signaling network based on cluster labels

% color5 = zeros(size(W,1),3);
% for i = 1:5
%     color5(find(true_labs1==i),:) = ones(length(find(true_labs1==i)),1).*zzz(i,:);
% end
lgd = cell(1,No_cluster);
for i = 1:No_cluster
    if i<10
        vv = 'CC';
        vv(2:2) = num2str(i);
        lgd{i} = vv;
    else
        vv = 'CCC';
        vv(2:3) = num2str(i);
        lgd{i} = vv;
    end
end

zz = (0:No_cluster)./(No_cluster);
tickval = 0.5.*(zz(2:end) + zz(1:end-1));

% plot cell-cell signlaing network based on each individual ligand-receptor pair
for j = 1:No_LR
    P = Pidv{j};
    P(P<=threshold) = 0;
    adjacentM =P(new_order,new_order);
    adjacentM(1:No_cells+1:end) = 0;
    if max(adjacentM(:)) > 0
        adjacentM = adjacentM./max(adjacentM(:));
    end
    bg = digraph(adjacentM);
    bg.Edges.LWidths = 3*bg.Edges.Weight/max(bg.Edges.Weight);
    figure;
    
    Gh = plot(bg,'Marker','o','MarkerSize',5,'Layout','circle','NodeLabel',[],...
        'NodeColor',mycolor_cells,'EdgeColor',[0.690196 0.768627 0.870588],'LineWidth',...
        bg.Edges.LWidths,'ArrowSize',8);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    a = Lig{j};
    b= Rec{j};
    title([a{1} '\_' b{1}],'fontsize',12);
    colormap(mycolor);
    colorbar('Ticks',tickval,'TickLabels',lgd,'FontSize',12);
    
    print([folder '\Cell_LR_' a{1} '_' b{1}],'-dpdf','-r300');
end


%% plot cell-cell signlaing network based on all ligand-receptor pairs
P = Pall;
P(P<=threshold) = 0;
adjacentM =P(new_order,new_order);
adjacentM(1:No_cells+1:end) = 0;
if max(adjacentM(:)) > 0
    adjacentM = adjacentM./max(adjacentM(:));
end
bg = digraph(adjacentM);
bg.Edges.LWidths = 3*bg.Edges.Weight/max(bg.Edges.Weight);
figure;

Gh = plot(bg,'Marker','o','MarkerSize',5,'Layout','circle','NodeLabel',[],...
    'NodeColor',mycolor_cells,'EdgeColor',[0.690196 0.768627 0.870588],'LineWidth',...
    bg.Edges.LWidths,'ArrowSize',8);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('cell\_LR\_all\_pairs','fontsize',12);
colormap(mycolor);
colorbar('Ticks',tickval,'TickLabels',lgd,'FontSize',12);

print([folder '\Cell_LR_all_pairs'],'-dpdf','-r300');



%% plot cluster-cluster signaling network based on each individual ligand-receptor pair

for j = 1:No_LR
    P = Pidv{j};
    P(P<=threshold) = 0;
    P_cluster = zeros(No_cluster);
    for i1 = 1:No_cluster
        for j1 = 1:No_cluster
            P_cluster(i1,j1) = sum(sum(P(cluster_label==i1,cluster_label ==j1)));
        end
    end
    
    adjacentM =P_cluster;
    adjacentM(1:No_cluster+1:end) = 0;
    if max(adjacentM(:)) > 0
        adjacentM = adjacentM./max(adjacentM(:));
    end
    bg = digraph(adjacentM);
    bg.Edges.LWidths = 10*bg.Edges.Weight/max(bg.Edges.Weight);
    figure;
    Gh = plot(bg,'Marker','o','MarkerSize',25,'Layout','circle','NodeLabel',[],...
        'NodeColor',mycolor,'EdgeColor',[0.690196 0.768627 0.870588],'LineWidth',...
        bg.Edges.LWidths,'ArrowSize',10);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    a = Lig{j};
    b= Rec{j};
    title([a{1} '\_' b{1}],'fontsize',12);
    colormap(mycolor);
    colorbar('Ticks',tickval,'TickLabels',lgd,'FontSize',12);
    
    print([folder '\Cluster_LR_' a{1} '_' b{1}],'-dpdf','-r300');
end


%% plot cluster-cluster signlaing network based on all ligand-receptor pairs


    P = Pall;
    P(P<=threshold) = 0;
    P_cluster = zeros(No_cluster);
    for i1 = 1:No_cluster
        for j1 = 1:No_cluster
            P_cluster(i1,j1) = sum(sum(P(cluster_label==i1,cluster_label ==j1)));
        end
    end
    
    adjacentM =P_cluster;
    adjacentM(1:No_cluster+1:end) = 0;
    if max(adjacentM(:)) > 0
        adjacentM = adjacentM./max(adjacentM(:));
    end
    bg = digraph(adjacentM);
    bg.Edges.LWidths = 10*bg.Edges.Weight/max(bg.Edges.Weight);
    figure;
    Gh = plot(bg,'Marker','o','MarkerSize',25,'Layout','circle','NodeLabel',[],...
        'NodeColor',mycolor,'EdgeColor',[0.690196 0.768627 0.870588],'LineWidth',...
        bg.Edges.LWidths,'ArrowSize',10);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    title('Cluster\_LR\_all\_pairs','fontsize',12);
    colormap(mycolor);
    colorbar('Ticks',tickval,'TickLabels',lgd,'FontSize',12);
    
    print([folder '\Cluster_LR_all_pairs'],'-dpdf','-r300');







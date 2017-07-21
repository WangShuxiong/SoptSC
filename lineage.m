function [Tree,pred,cluster_center] = lineage(idx,W,root_cell)
% [dvis,~] = eigs(W,3);
dvis = pca(W,3);
% H = H_all;
cluster_idx = unique(idx,'stable');
NN = max(idx);
DisM = zeros(NN);

cluster_center = zeros(3,NN);

for i = 1:NN
    zz = find(idx==cluster_idx(i));
    cluster_center(:,i) = mean(dvis(zz,:));
    
    if ismember(root_cell,zz)
        root_cluster = cluster_idx(i);
    end
end

% computing cluster-to-cluster distance
for i = 1:NN
    for j = 1:NN
        DisM(i,j) = norm(cluster_center(:,i)-cluster_center(:,j),1);
    end
end

bg = graph(DisM);
% Gh = plot(bg);
% Gh.LineWidth = bg.Edges.Weight;
% d = distances(bg);
% [~,rootvalue] = max(sum(d));
rootvalue = find(cluster_idx==root_cluster);
% display(rootvalue);

[Tree,pred] = minspantree(bg,'Root',rootvalue); % ,'Method','Kruskal'

rootedTree = digraph(pred(pred~=0),find(pred~=0));
% eLabels = Tree.Edges.Weight;
% nLabels = ;
MS = 10*ones(NN,1);
figure(4);
% Gtree = plot(rootedTree,'EdgeLabel',eLabels,'Markersize',MS); % 'Layout','force','EdgeLabel',eLabels,'NodeLabel',nLabels
plot(rootedTree,'Markersize',MS,'NodeLabel',cluster_idx); % 'Layout','force','EdgeLabel',eLabels,'NodeLabel',nLabels

set(gca,'xtick',[]); 
set(gca,'ytick',[]);
print(4,'-dtiff','Results\Lineage.tiff');
%Gtree.LineWidth = Tree.Edges.Weight;
% display(DisM);
% 
% disp('============================');
% disp('Tree predessensor:');
% fprintf('%d\n',pred);
% disp('============================');
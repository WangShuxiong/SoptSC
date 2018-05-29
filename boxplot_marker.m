function boxplot_marker(data,allgenes,marker,cluster_labs,No_cluster,folder)
% Box plot for each gene along all clusters

% colormap jet;
cmap1 = jet;
mymap1 = cmap1(1:end,:);
ncolor = size(mymap1,1);
mycolor = mymap1(1:round(ncolor./No_cluster):1+ncolor,:);

group = cell(size(cluster_labs));
for i = 1:length(cluster_labs)
    group{i} = ['C' num2str(cluster_labs(i))];
end


ia = zeros(1,length(marker));
for i = 1:length(marker)
    for j = 1:length(allgenes)
        if strcmp(upper(marker{i}),upper(allgenes{j}))
            ia(i) = j;
        end
    end
end


gname = allgenes(ia);

display(allgenes(ia));

MM = data(ia,:);

for i = 1:No_cluster
    cluster_notation{i} = ['C' num2str(i)];
end

n = length(ia);
for i = 1:n
    figure;
    boxplot(MM(i,:),group,'GroupOrder',cluster_notation,'Notch','on','PlotStyle','compact','Widths',0.9,'Colors',mycolor,...
        'LabelOrientation','horizontal'); 
    title(gname{i});
    print([folder '\boxplot_mk_' gname{i}],'-dpdf','-r300');
end
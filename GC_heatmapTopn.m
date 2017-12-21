function Gene_labels = GC_heatmapTopn(data,cluster_label,H,allgenes,gene_idx,topn)
% This function assign genes to each cluster by SoptSC
%
%   Input:
%       data: gene-cell matrix
%       cluster_label: cluster labels for all cells
%       H: non-negative matrix such that W = H*H^T
%
%   Output:
%       Gene_labels: gene label information for each gene associated with a specific
%       cluster, the first columns of Gene_labels represents gene indices;
%       the second column of Gene_labels represents the cluster index that 
%       the gene belongs to; the third column represent gene score associated
%       with corresponding cluster.
Dataplot = data;

NC = max(cluster_label);
m = size(data,1);
Gene_labels = zeros(m,3);

Gene_labels(:,1) = gene_idx;


%% data normalization
for i = 1:size(data,2)
    data(:,i) = data(:,i)./norm(data(:,i),2);
end



G_latent = data*H;
[Gene_value,Gene_label] = max(G_latent,[],2);

Gene_labels(:,2) = Gene_label;
Gene_labels(:,3) = Gene_value;


OGI = [];
OGV = [];

CGI = [];
topn1 = topn;
for i = 1:NC
    % order genes within each cluster
    Z = find(Gene_label==i);
    Z1 = Gene_value(Z);
    [Z1V,I] = sort(Z1,'descend');
    
    if topn1 > length(Z)
        topn1 = length(Z);
    end
    
    display(topn1);
    Z2 = Z(I(1:topn1));
    OGI = [OGI; Z2];
    OGV = [OGV; Z1V(1:topn1)];
    
    %order cells
    Y = find(cluster_label==i);
    CGI = [CGI; Y];
    topn1 = topn;
end

%% allgenes string to cell
% allgs = cell(size(allgenes));
% for i = 1:length(allgs)
%     allgs{i} = allgenes{i};
% end
% %% plot gene-cell heatmap
% RowLabelsValue = allgs(gene_idx);
% % HeatMap(data(OGI,CGI),'RowLabels', ColumnLabelsValue)
% HMdata = data(OGI,CGI);
% 
% HeatMap(HMdata,'RowLabels', RowLabelsValue,'Standardize',2,'DisplayRange',2,'Colormap',redbluecmap)


%% data normalization and zscroe
% kk = 1 row; kk = 2 column
figure;
idata = Dataplot(OGI,CGI);
kk = 2;
center = mean(idata,kk);
scale = std(idata, 0,kk);

tscale = scale;
%=Check for zeros and set them to 1 so not to scale them.
scale(tscale == 0) = 1;
%== Center and scale the data
idata = bsxfun(@minus, idata, center);
sdata = bsxfun(@rdivide, idata, scale);

colormap redbluecmap;
clims = [-3 3];
imagesc(sdata,clims);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

yticks(1:length(OGI));
yticklabels(allgenes(OGI));
% xticks(1:length(CGI));
% xticklabels(CGI);

% cb = colorbar;
% ax = gca;
% axpos = ax.Position;
% cpos = cb.Position;
% cpos(3) = 0.5*cpos(3);
% cb.Position = cpos;
% ax.Position = axpos;


% % 
% % hold on;
% % a = get(gca,'Xlim');
% % b = get(gca,'Ylim');
% % 
% % axis([a(1) a(2) b(1) b(2)+1]);
% % 
% % 
% % for i = 1:NC
% %     bili(i) = length(find(cluster_label ==i));
% % end
% % bili = bili./sum(bili);
% % bili = cumsum(bili);
% % 
% % xpoint = a(1) + (a(2)-a(1)).*bili;
% % ypoint = b(2)+1;
% % 
% % point = [xpoint; ypoint.*ones(size(xpoint))];
% % 
% % for i = 1:NC-1
% %     line(point(:,i),point(:,i+1),'Color','r');
% %     hold on;
% % end

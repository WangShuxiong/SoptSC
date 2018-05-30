function P = LRact_Prob(data,allgenes,L,R,TG)

[~,n] = size(data);
No_TGs = length(TG);
TG_data = zeros(No_TGs,n);

[~,L_idx,~] = intersect(allgenes,L,'stable');
[~,R_idx,~] = intersect(allgenes,R,'stable');

% Normalize L_data and R_data along gene
% Normalize L_data and R_data along cell
if max(data(L_idx,:)) <=0 
    L_data = data(L_idx,:);
    display(['Error: the selected ligand is not expressed: ' allgenes{L_idx}]);
else
    L_data = data(L_idx,:)./max(data(L_idx,:));
end

if max(data(R_idx,:)) <=0 
    R_data = data(R_idx,:);
    display(['Error: the selected receptor is not expressed: ' allgenes{R_idx}]);
else
    R_data = data(R_idx,:)./max(data(R_idx,:));
end
% L_data = data(L_idx,:);
% R_data = data(R_idx,:);

display(allgenes{L_idx});
display(allgenes{R_idx});

for i = 1:No_TGs
    aa = TG{i};
    [~,TG_idx,~] = intersect(allgenes,aa,'stable');
    if max(data(TG_idx,:)) > 0
        TG_data(i,:) = data(TG_idx,:)./max(data(TG_idx,:));
        display(allgenes{TG_idx});
    else
         display(['Note: the unexpressed target gene is excluded: ' allgenes{TG_idx}]);
        TG_data(i,:) = data(TG_idx,:);
    end
end

% normal_c = max([max(L_data) max(R_data) max(TG_data(:))]);
% L_data = data(L_idx,:)./normal_c;
% R_data = data(R_idx,:)./normal_c;
% TG_data = TG_data./normal_c;

if No_TGs > 1
    TG_all = mean(TG_data);
else
    TG_all = TG_data;
end

Coef =zeros(n);
for i = 1:n
    for j = 1:n
        alpha = (exp(-1./(L_data(i).*R_data(j))));
        beta = exp(-1./TG_all(j));
        if beta == 0
            Coef(i,j) = 0;
        else
            Coef(i,j) = alpha./(alpha + beta);
%         Coef(i,j) = (exp(-1./(L_data(i)+R_data(j))))./(exp(-1./(L_data(i)+R_data(j))) + ...
%             exp(-1./TG_all(j)));
        end
    end
end

P = zeros(n);
for i = 1:n
    b = 0;
    for kk = 1:n
        b = b + (exp(-1./(L_data(i).*R_data(kk)))).*( exp(-1./TG_all(kk)) ).*Coef(i,kk);
    end
    
    for j = 1:n
        a = (exp(-1./(L_data(i).*R_data(j)))).*( exp(-1./(TG_all(j)) ) ).*(Coef(i,j));
        if b == 0
            P(i,j) = 0;
        else
            P(i,j) = a./b;
        end
    end
end



% network visulization
% %% select top edges
% topn = 0.01;
% ntop = round(topn*n*(n-1)*0.5);
% Wtop = P;
% ZZ = abs(Wtop(:));
% B = sort(ZZ,'descend');
%
% Wtop(abs(Wtop)<=B(ntop)) = 0;
% imagesc(Wtop);
%
% %% plot graph
% % nLabels = true_time;
% nLabels = cluster_label;
% adjacentM = max(abs(Wtop),abs(Wtop'));
% adjacentM(1:n+1:end) = 0;
% bg = graph(adjacentM);
% bg.Edges.LWidths = 5*bg.Edges.Weight/max(bg.Edges.Weight);
% % Gh.NodeCData = true_time;
% Gh = plot(bg,'Marker','s','MarkerSize',5,'Layout','force','NodeLabel',[],...
%     'NodeCData',nLabels,'NodeColor','flat','LineWidth',bg.Edges.LWidths,'MarkerSize',5);
%
% for i=1:length(nLabels)
%    text(Gh.XData(i)+0.01,Gh.YData(i),num2str(nLabels(i)),'fontsize',12);
% end
%
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);


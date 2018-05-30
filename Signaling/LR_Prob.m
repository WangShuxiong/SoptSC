function P = LR_Prob(data,allgenes,L,R)

[~,n] = size(data);

[~,L_idx,~] = intersect(allgenes,L,'stable');
[~,R_idx,~] = intersect(allgenes,R,'stable');

% Normalize L_data and R_data along gene
% Normalize L_data and R_data along cell
if max(data(L_idx,:)) <=0 
   L_data = data(L_idx,:);
else
    L_data = data(L_idx,:)./max(data(L_idx,:));
end

if max(data(R_idx,:)) <=0 
    R_data = data(R_idx,:);
else
    R_data = data(R_idx,:)./max(data(R_idx,:));
end


display([nnz(L_data) nnz(R_data)]);

display(allgenes{L_idx});
display(allgenes{R_idx});




P = zeros(n);
for i = 1:n
    b = 0;
    for kk = 1:n
        b = b + (exp(-1./(L_data(i).*R_data(kk))));
    end
    
    for j = 1:n
        a = (exp(-1./(L_data(i).*R_data(j))));
        if b == 0 || a == 0
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


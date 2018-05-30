function P = LRboth_Prob(data,allgenes,L,R,TG_act,TG_inh)


[~,n] = size(data);

TG = TG_act;
TG1 = TG_inh;
No_TGs = length(TG);
No_TGs1 = length(TG1);

TG_data = zeros(No_TGs,n);
TG_data1 = zeros(No_TGs1,n);

[~,L_idx,~] = intersect(allgenes,L,'stable');
[~,R_idx,~] = intersect(allgenes,R,'stable');

% Normalize L_data and R_data along gene
% Normalize L_data and R_data along cell
if max(data(L_idx,:)) <=0 
    L_data = data(L_idx,:);
%    display(['Error: the selected ligand is not expressed: ' allgenes{L_idx}]);
else
    L_data = data(L_idx,:)./max(data(L_idx,:));
end

if max(data(R_idx,:)) <=0 
    R_data = data(R_idx,:);
%    display(['Error: the selected receptor is not expressed: ' allgenes{R_idx}]);
else
    R_data = data(R_idx,:)./max(data(R_idx,:));
end

display(allgenes{L_idx});
display(allgenes{R_idx});

for i = 1:No_TGs
    aa = TG{i};
    [~,TG_idx,~] = intersect(allgenes,aa,'stable');
    if max(data(TG_idx,:)) > 0
        TG_data(i,:) = data(TG_idx,:)./max(data(TG_idx,:));
        display(allgenes{TG_idx});
    else
%         display(['Note: the unexpressed target gene is excluded: ' allgenes{TG_idx}]);
        TG_data(i,:) = data(TG_idx,:);
    end
end


for i = 1:No_TGs1
    aa = TG1{i};
    [~,TG_idx1,~] = intersect(allgenes,aa,'stable');
    if max(data(TG_idx1,:)) > 0
        TG_data1(i,:) = data(TG_idx1,:)./max(data(TG_idx1,:));
        display(allgenes{TG_idx1});
    else
%         display(['Note: the unexpressed target gene is excluded: ' allgenes{TG_idx1}]);
        TG_data1(i,:) = data(TG_idx1,:);
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


if No_TGs1 > 1
    TG_all1 = mean(TG_data1);
else
    TG_all1 = TG_data1;
end

% Comput C_i,j and E_i,j
Eoef = zeros(n);
Coef =zeros(n);
for i = 1:n
    for j = 1:n
        alpha = (exp(-1./(L_data(i).*R_data(j))));
        beta = exp(-1./TG_all(j));
        gamma = exp(-TG_all1(j));

        if beta == 0
            Coef(i,j) = 0;
        else
            Coef(i,j) = alpha./(alpha + beta);
        end
        
        if gamma == 0
            Eoef(i,j) = 0;
        else
            Eoef(i,j) = alpha./(alpha + gamma);
        end
        
    end
end

P = zeros(n);
for i = 1:n
    b = 0;
    for kk = 1:n
        b = b + (exp(-1./(L_data(i).*R_data(kk)))).*( exp(-1./TG_all(kk)) ).*(Coef(i,kk)).*(exp(-TG_all1(kk))).*(Eoef(i,kk));
    end
    
    for j = 1:n
        a = (exp(-1./(L_data(i).*R_data(j)))).*( exp(-1./(TG_all(j)) ) ).*(Coef(i,j)).*(exp(-TG_all1(j))).*(Eoef(i,j));
        if b == 0
            P(i,j) = 0;
        else
            P(i,j) = a./b;
        end
    end
end


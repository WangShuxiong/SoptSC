function [eigenvalues,No_cluster] = Num_cluster(W)
% Computing number of cluster
nno_cluster = 1:25; tau = 0.3; tol = 0.01;
[all_eigs,M_all] = consen(W,nno_cluster,tau);

zz = sort(abs(real(all_eigs)));

% zz_idx = find(zz<=1);
% ZZ = zz(zz_idx);
ZZ = zz;

if length(ZZ)>=2
    gap = ZZ(2:end) - ZZ(1:end-1);
    [~,No_cluster1] = max(gap);
%     display(gapval);
%     display(No_cluster1);
end

No_cluster2 = length(find(zz<=tol));
No_cluster = No_cluster1;
display('Number of cluster based on zero eigenvalues & Largest gap ');
display([No_cluster2 No_cluster]);

% display('First 50 eigenvalues of graph Laplacian:');
% display(zz(1:50));
% 
% figure;
% scatter(1:30,zz(1:30),20,'filled');
% box on;
% set(gca,'LineWidth',1.5);
% xlabel('i');
% ylabel('Eigenvalue of graph Laplacian \lambda_i');
% set(gca,'FontName','Arial');
% set(gca,'FontSize',12);
eigenvalues = zz;
end
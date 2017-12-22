function [eigenvalues,No_cluster] = Num_cluster(W)
% Computing the number of clusters

nno_cluster = 1:25; tau = 0.3; tol = 0.01;
[all_eigs,M_all] = consen(W,nno_cluster,tau);

zz = sort(abs(real(all_eigs)));

ZZ = zz;

if length(ZZ)>=2
    gap = ZZ(2:end) - ZZ(1:end-1);
    [~,No_cluster1] = max(gap);
end

% No_cluster2 = length(find(zz<=tol));
No_cluster = No_cluster1;
display('Number of clusters inferred by SoptSC is:');
display(No_cluster);

eigenvalues = zz;
end
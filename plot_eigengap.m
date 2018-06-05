function plot_eigengap(eigenvalues,folder)
figure;
scatter(1:min([30 size(eigenvalues,1)]),eigenvalues(1:min([30 size(eigenvalues,1)])),20,'filled');
box on;
set(gca,'LineWidth',1.5);
xlabel('i');
ylabel('Eigenvalue of graph Laplacian \lambda_i');
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
print([folder '\EigenGap'],'-dpdf','-r300'); 
end

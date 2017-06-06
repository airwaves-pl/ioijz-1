
x = [2.5 0.5 2.2 1.9 3.1 2.3 2 1 1.5 1.1; 2.4 0.7 2.9 2.2 3 2.7 1.6 1.1 1.6 0.9];
s = cov(x');
[eigenvectors, eigenvalues] = eigs(s, 2);
pcs = (x') * eigenvectors;
figure;
hold on
plot(x(1:1,:), x(2:2,:), 'or', pcs(:,1:1), pcs(:,2:2), 'ob')
plotv(eigenvectors(:,1), '-r');
plotv(eigenvectors(:,2), '-b');
legend('dane oryginalne','sk³adowe g³ówne','wektor cechy 1','wektor cechy 2','Location','northwest','Orientation','vertical')
hold off
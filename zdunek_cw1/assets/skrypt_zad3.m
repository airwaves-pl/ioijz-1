
x = [2.5 0.5 2.2 1.9 3.1 2.3 2 1 1.5 1.1; 2.4 0.7 2.9 2.2 3 2.7 1.6 1.1 1.6 0.9];
s = cov(x');

[eigenvectors, eigenvalues] = eig(s);
eigenvector_power = Power_f(s,1e-10);
T=lanczos(s, [1; 1], 200);

figure;
hold on
plot(x(1:1,:), x(2:2,:), 'ok')
plotv(eigenvector_power(:,1), '-r');
plotv(eigenvectors(:,2), '-g');

legend('dane oryginalne','wektor w³asny (power)','wektor w³asny (eig)','Location','northwest','Orientation','vertical')
hold off


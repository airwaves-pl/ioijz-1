% dane wejsciowe
x = [2.5 0.5 2.2 1.9 3.1 2.3 2 1 1.5 1.1;
    2.4 0.7 2.9 2.2 3 2.7 1.6 1.1 1.6 0.9];
% obliczenie wartosci srednich
% odjecie ich od macierzy
x_norm = bsxfun(@minus, x, mean(x, 2));
% wyznaczenie macierzy kowariancji
s = cov(x_norm');
% wyznaczenie wartosci i wektorow wlasnych
[eigenvectors, eigenvalues] = eigs(s, 2);
% rzutowanie
pcs = eigenvectors' * x;

figure;
hold on
plot(x(1,:), x(2,:), 'or', pcs(1,:), pcs(2,:), 'ob')
plotv(eigenvectors(:,1), '-');
plotv(eigenvectors(:,2), '-');
legend('dane oryginalne', 'skladowe glowne', 'wektor cechy 1','wektor cechy 2','Location','northeast','Orientation','vertical')
hold off

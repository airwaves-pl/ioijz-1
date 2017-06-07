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
% porzadkujemy wartosci wlasne w kolejnosci malejacej
[eigenvalues order] = sort(diag(eigenvalues), 'descend');
eigenvectors = eigenvectors(:,order);
% rzutowanie
pcs = eigenvectors' * x;

figure;
hold on
plot(x(2,:), x(1,:), 'or', pcs(2,:), pcs(1,:), 'ob')
plotv(eigenvectors(:,1), '-r');
plotv(eigenvectors(:,2), '-b');
legend('dane oryginalne','skladowe glowne','wektor cechy 1','wektor cechy 2','Location','northeast','Orientation','vertical')
hold off

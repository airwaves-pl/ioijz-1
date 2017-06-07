clear all;
clc;

%macierz pocz¹tkowa (wiersze - pomiary, kolumny - rodzaje badañ(wymiary))
X = [2.5 0.5 2.2 1.9 3.1 2.3 2 1 1.5 1.1; 
    2.4 0.7 2.9 2.2 3 2.7 1.6 1.1 1.6 0.9];

%okreœlenie iloœci wierszy
n = size(X);
n = n(1);

%obliczenie macierzy kowariancji
R = cov(X');

%analiza wektorów w³asnych
[eigenvectors, eigenvalues] = eigs(R, n);

%sortowanie wektorów
eigenvectors = sort(eigenvectors, 'ascend')

%utworzenie wykresu
hold on;
plot(X(1:1,:), X(2:2,:), 'ob');

plotv(eigenvectors(:,1), '-r')
plotv(eigenvectors(:,2), '-b')
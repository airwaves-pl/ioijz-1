clear all;
clc;

%macierz pocz�tkowa (wiersze - pomiary, kolumny - rodzaje bada�(wymiary))
X = [2.5 0.5 2.2 1.9 3.1 2.3 2 1 1.5 1.1; 
    2.4 0.7 2.9 2.2 3 2.7 1.6 1.1 1.6 0.9];

%okre�lenie ilo�ci wierszy
n = size(X);
n = n(1);

%obliczenie macierzy kowariancji
R = cov(X');

%analiza wektor�w w�asnych
[eigenvectors, eigenvalues] = eigs(R, n);

%sortowanie wektor�w
eigenvectors = sort(eigenvectors, 'ascend')

%utworzenie wykresu
hold on;
plot(X(1:1,:), X(2:2,:), 'ob');

plotv(eigenvectors(:,1), '-r')
plotv(eigenvectors(:,2), '-b')
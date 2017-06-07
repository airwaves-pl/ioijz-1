function [eigenvectors] = PCA_f(X)
    %okreœlenie iloœci wierszy
    n = size(X);
    n = n(1);

    %obliczenie macierzy kowariancji
    R = cov(X');

    %analiza wektorów w³asnych
    [eigenvectors, eigenvalues] = eigs(R, n);
    
    %sortowanie wektorów
    eigenvectors = sort(eigenvectors, 'ascend')

end
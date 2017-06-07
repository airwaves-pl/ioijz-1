function [eigenvectors] = PCA_f(X)
    %okre�lenie ilo�ci wierszy
    n = size(X);
    n = n(1);

    %obliczenie macierzy kowariancji
    R = cov(X');

    %analiza wektor�w w�asnych
    [eigenvectors, eigenvalues] = eigs(R, n);
    
    %sortowanie wektor�w
    eigenvectors = sort(eigenvectors, 'ascend')

end
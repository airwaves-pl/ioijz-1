function [A,X,res,MSE] = skrypt_zad1_nmf_mue(A,X,Y,N)
    MSE = [];
    for k = 1:N
        % obliczanie A
        s1 = A.*(Y*X');
        s2 = (A*X*X');
        A = max(0, s1./s2);  
        A = A*diag(1./sum(A,1));

        % obliczanie X
        s1 = X.*(A'*Y);
        s2 =A'*A*X;
        X = max(0, s1./s2);

        % blad residualny
        res(k) = norm(Y - A*X, 'fro')/norm(Y,'fro');
        
        % blad srednio-kwadratowy
        MSE = [MSE immse(Y,A*X)];
    end
end
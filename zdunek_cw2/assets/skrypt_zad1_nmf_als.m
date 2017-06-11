function[A,X,res] = skrypt_zad1_nmf_als(A,X,Y,N)
    for k = 1:N
        % obliczanie A
        A = max(0,Y*X'*inv(X*X'));
        A = A*diag(1./sum(A,1));

        % obliczanie X
        X = max(0,inv(A'*A)*A'*Y);
        
        % blad residualny
        res(k) = norm(Y - A*X,'fro')/norm(Y,'fro');
    end
end
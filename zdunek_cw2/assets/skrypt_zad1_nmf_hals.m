function [A,X,res] = skrypt_zad1_nmf_hals(A,X,Y,J,MaxIter)
    for k = 1:MaxIter
        for j = 1:J % obliczanie A - po kolumnach
            YXp = Y*X';
            XXp = X*X';
            A(:,j) = max(0, A(:,j)+(YXp(:,j)-A*XXp(:,j))/(XXp(j,j)+eps));  
            A = A*diag(1./sum(A,1));
        end
        for j = 1:J % obliczanie X - po wierszach
            ApY = A'*Y;
            ApA = A'*A;
            X(j,:) = max(0, X(j,:)+(ApY(j,:)-ApA(j,:)*X)/(ApA(j,j)+eps));
        end
        res(k) = norm(Y - A*X, 'fro')/norm(Y,'fro'); % blad residualny
    end
end
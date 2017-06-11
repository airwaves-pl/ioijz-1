function [A,X,res] = NMF_HALS(A,X,Y,J,MaxIter)
  for k = 1:MaxIter
    for j = 1:J
      % obliczanie A - po kolumnach
      aj = A(:,j);
      
      s1 = (Y*X');
      s1 = s1(:,j);
      
      s2 = A*(X*X');
      s2 = s2(:,j);
      
      s3 = (X*X');
      s3 = s3(j:j);
      
      combined = (aj + ((s1 - s2) / s3));
      
      A(:,j) = max(0, combined);  
      A = A*diag(1./sum(A,1));
    end
    for j = 1:J
        % obliczanie X - po wierszach
        xj = X(j,:);
        
        s1 = (A'*Y);
        s1 = s1(j,:);
        
        s2 = (A'*A);
        s2 = s2(j,:);
        s2 = s2*X;
        
        s3 = (A'*A);
        s3 = s3(j:j);
        
        combined(xj + ((s1-s2)/s3));
        
        X(j,:) = max(0, combined);
    end
    % blad residualny
    res(k) = norm(Y - A*X, 'fro')/norm(Y,'fro');
  end
end
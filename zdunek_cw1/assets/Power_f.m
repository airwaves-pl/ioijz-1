function [new_result] = Power_f(X, tolerance)

    n = size(X);
    n = n(1);
    
    new_result = ones(n,1);
    m = 0;
    
    while(1)
        m_old = m;
        old_result = new_result;    
        new_result = X * new_result;

        m = max(new_result);
        new_result = new_result/m;
        
        if abs(m-m_old) < tolerance && norm(new_result-old_result,2) < tolerance
            break;
        end
    end
end
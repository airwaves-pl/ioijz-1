clear

% Dane oryginalne
I = 100; T = 1000; J = 10;
Aw = max(0,randn(I,J));
Xw = max(0,randn(J,T));

Y = Aw*Xw;

% Inicjalizacja
A = rand(size(Y,1),J);
X = rand(J,size(Y,2));

MaxIter = 100;

for k = 1:MaxIter
    
    A = max(0,Y*X'*inv(X*X'));
    A = A*diag(1./sum(A,1));
    
    X = max(0,inv(A'*A)*A'*Y);
    res(k) = norm(Y - A*X,'fro')/norm(Y,'fro');
    
end

% to wy¿ej to algorytm optymalizacji naprzemiennej

figure
semilogy(res)

% to samo ale z wykorzystaniem wbudowanej funkcji

rng(1) % For reproducibility
[W,H] = nnmf(Y,10);
D = norm(Y-W*H,'fro')/sqrt(I*T); % b³¹d residualny


SIR = CalcSIR(Aw,A);

% ---- to co napisal na konsoli ----

Aws = Aw*diag(1./sum(Aw,1));
grid on
eps
Aws(1:5,:)


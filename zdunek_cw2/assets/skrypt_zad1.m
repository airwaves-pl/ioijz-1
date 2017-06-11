clc;
clear all;
close all;

% dane oryginalne
I = 100; T = 1000; J = 10;
Aw = max(0,randn(I,J));
Xw = max(0,randn(J,T));

Y = Aw*Xw;

% inicjalizacja
A = rand(size(Y,1),J);
X = rand(J,size(Y,2));

A1 = A;
A2 = A;
A3 = A;

X1 = X;
X2 = X;
X3 = X;

MaxIter = 100;

[A1,X1,res1] = NMF_ALS(A1,X1,Y,MaxIter);
[A2,X2,res2] = NMF_MUE(A2,X2,Y,MaxIter);
%[A3,X3,res3] = NMF_HALS(A3,X3,Y,J,MaxIter);

figure
hold on;
semilogy(res1)
semilogy(res2)
%semilogy(res3)
legend('ALS', 'MUE')

% to samo ale z wykorzystaniem wbudowanej funkcji

rng(1) % for reproducibility
[W,H] = nnmf(Y,10);
D = norm(Y-W*H,'fro')/sqrt(I*T); % blad residualny (resztowy)


SIR = CalcSIR(Aw,A);

% ---- to co pisal na konsoli ----

Aws = Aw*diag(1./sum(Aw,1));
grid on
eps
Aws(1:5,:)


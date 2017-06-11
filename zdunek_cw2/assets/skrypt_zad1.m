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

A1 = A; A2 = A; A3 = A;
X1 = X; X2 = X; X3 = X;

MaxIter = 300;

[A1,X1,res1] = skrypt_zad1_nmf_als(A1,X1,Y,MaxIter);
[A2,X2,res2] = skrypt_zad1_nmf_mue(A2,X2,Y,MaxIter);
[A3,X3,res3] = skrypt_zad1_nmf_hals(A3,X3,Y,J,MaxIter);

figure
semilogy(res1);
hold on;
semilogy(res2);
semilogy(res3);
legend('ALS', 'MUE', 'HALS');
hold off;
grid on;

SIR_A1 = CalcSIR(A',A1');
SIR_A2 = CalcSIR(A',A2');
SIR_A3 = CalcSIR(A',A3');
%%
figure;
hold on
title('Wartosci SIR');
plot(1:100, SIR_A1, 1:100, SIR_A2, 1:100, SIR_A3);
xlim([1 100]);
legend('ALS','MUE','HALS', 'Location','southeast','Orientation','vertical')
ylabel('SIR')
hold off

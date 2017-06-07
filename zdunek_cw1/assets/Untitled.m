
testb = [1 2 3 4 5];
testa = [1 2 2 3 4];

[C,order] = confusionmat(testb, testa);
x = diag(C);
y = [1 0 1 1 0];

plotconfusion(x', y);

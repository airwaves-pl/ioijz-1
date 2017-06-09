function [ delta_hosvd,elapsed_time_hosvd ] = Tucker_orth_3D( Y_train,Y_test,Class_train_idx,Class_test_idx,J )

% TRENOWANIE

%% HOSVD

tic

DimY = size(Y_train);

% unfolding
Y1 = reshape(Y_train,DimY(1),DimY(2)*DimY(3));
Y2 = reshape(permute(Y_train,[2 1 3]),DimY(2),DimY(1)*DimY(3));
Y3 = reshape(permute(Y_train,[3 1 2]),DimY(3),DimY(1)*DimY(2));

% dekompozycja tensorów
[E1,D1] = eig(Y1*Y1');
Ar = fliplr(E1(:,DimY(1)-J(1)+1:DimY(1)));

[E2,D2] = eig(Y2*Y2');
Br = fliplr(E2(:,DimY(2)-J(2)+1:DimY(2)));

[E3,D3] = eig(Y3*Y3');
Cr = fliplr(E3(:,DimY(3)-J(3)+1:DimY(3)));

Gr = ntimes(ntimes(ntimes(Y_train,Ar',1,2),Br',1,2),Cr',1,2); % core tensor
Y_hat = ntimes(ntimes(ntimes(Gr,Ar,1,2),Br,1,2),Cr,1,2); % tensor 3-way

Cr = Cr.*repmat(1./sqrt(sum(Cr.^2,2)+eps),1,size(Cr,2));



elapsed_time_hosvd = toc;
%%

% TESTOWANIE

Y3 = reshape(permute(Y_test,[3 1 2]),size(Y_test,3),size(Y_test,1)*size(Y_test,2));
G3 = reshape(permute(Gr,[3 1 2]),[J(3),J(1)*J(2)]);

% matrycyzacja
Ct = Y3*pinv(double(G3)*(kron(Br,Ar))');
Ct = Ct.*repmat(1./sqrt(sum(Ct.^2,2)+eps),1,size(Ct,2));

% klasyfikacja 1-NN
mdl_class = fitcknn(Cr,Class_train_idx,'NumNeighbors',1);
prediction = predict(mdl_class, Ct); 

% dok³adnoœæ klasyfikacji
delta_hosvd = 100*(length(find((prediction - Class_test_idx)==0))/length(Class_test_idx));

end

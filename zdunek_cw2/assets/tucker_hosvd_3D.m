function [ A,B,C,G,Y_hat ] = tucker_hosvd( Y,J )

% Size
DimY = size(Y);

% Unfolding
Y1 = reshape(Y,DimY(1),DimY(2)*DimY(3));
Y2 = reshape(permute(Y,[2 1 3]),DimY(2),DimY(1)*DimY(3));
Y3 = reshape(permute(Y,[3 1 2]),DimY(3),DimY(1)*DimY(2));

% Tensor decomposition
[E1,D1] = eig(Y1*Y1');
A = fliplr(E1(:,DimY(1)-J(1)+1:DimY(1)));

[E2,D2] = eig(Y2*Y2');
B = fliplr(E2(:,DimY(2)-J(2)+1:DimY(2)));

[E3,D3] = eig(Y3*Y3');
C = fliplr(E3(:,DimY(3)-J(3)+1:DimY(3)));

G = ntimes(ntimes(ntimes(Y,A',1,2),B',1,2),C',1,2); % core tensor
Y_hat = ntimes(ntimes(ntimes(G,A,1,2),B,1,2),C,1,2); % tensor 3-way

C = C.*repmat(1./sqrt(sum(C.^2,2)+eps),1,size(C,2));

end

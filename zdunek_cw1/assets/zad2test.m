clc;
Persons = 3;
ImagesPerPerson = 10;
J = 10;
M = [];
for p=(1:Persons)
    for i=(1:ImagesPerPerson)
        x = FaceData(p, i).Image;
        [irow, icol] = size(x);
        x = double(x);
        temp = reshape(x', irow * icol, 1);
        M = [M temp];
    end
end

%subtract mean
x=bsxfun(@minus, M, mean(M,2));

% calculate covariance
s = cov(x');

% obtain eigenvalue & eigenvector
[V,D] = eigs(s,J);

Z = (x') * V;

figure;
for i=(1:J)
    C = V(:,i);
    CC = reshape(C, [46, 56]);
    CC = CC';
    subplot(round(sqrt(J)), round(sqrt(J)) + 1, i);
    imagesc(CC);
    title(i)
    colormap gray;
end

id = (1:J);
result = kmeans(Z', Persons);
groups = sortrows([id', result], 2);

T=[1 1 1 1; 2 2 2 2; 3 3 3 3];
[Acc,rand_index,match]=AccMeasure(T,result);


%Class = knnclassify(V, M, Group);



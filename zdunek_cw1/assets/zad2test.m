clc;
Persons = 10
ImagesPerPerson = 10
J = 20
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

[eigenvectors, eigenvalues] = eigs(M*(M'), J);
V = eigenvectors;
figure;
for i=(1:J)
    C = V(:,i);
    CC = reshape(C, [46, 56]);
    CC = CC';
    subplot(round(sqrt(J)) + 1, round(sqrt(J)) + 1, i);
    imagesc(CC);
    colormap gray;
end

close all;
clc;
clear;
load('FaceData_56_46.mat');

Persons = 3
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
    subplot(round(sqrt(J)), round(sqrt(J)) + 1, i);
    imagesc(CC);
    title(i)
    colormap gray;
end

id = (1:J);
result = kmeans(V', Persons);
result = [id', result];
%sortrows(result, 2);

%-----------------------
% ORYGINA�Y
%-----------------------

n = Persons*ImagesPerPerson;

id2 = (1:n);
result2 = kmeans(M', Persons);
result2 = [id2', result2];
%sortrows(result2, 2);

figure;

for i=(1:n)
    %C = V(:,i);
    %CC = reshape(C, [46, 56]);
    %CC = CC';
    subplot(round(sqrt(n)), round(sqrt(n)) + 1, i);
    imagesc(reshape(M(:,i), [46, 56])');
    title(i)
    colormap gray;
end
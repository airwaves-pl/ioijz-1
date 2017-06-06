clc;
clear;
close all;

load('yalefaces.mat'); 
load('FaceData_56_46.mat');

%imagesc(yalefaces(:,:,1))

% -- FRAGMENT Z WIKIPEDII --

% J = 30;
% [h,w,n] = size(yalefaces);
% d = h*w;
% % vectorize images
% x = reshape(yalefaces,[d n]);
% x = double(x);
% %subtract mean
% x=bsxfun(@minus, x, mean(x,2));
% % calculate covariance
% s = cov(x');
% % obtain eigenvalue & eigenvector
% [V,D] = eigs(s,J);
% eigval = diag(D);
% % sort eigenvalues in descending order
% eigval = eigval(end:-1:1);
% V = fliplr(V);
% % show 0th through 15th principal eigenvectors
% eig0 = reshape(mean(x,2), [h,w]);
% figure,subplot(6,6,1)
% imagesc(eig0)
% colormap gray
% for i = 1:J
% 	subplot(6,6,i+1)
% 	imagesc(reshape(V(:,i),h,w))
% end


% -- NASZ KOD --

Persons = 3;
ImagesPerPerson = 10;
J = 10;

Group = [];
M = [];
for p=(1:Persons)
    for i=(1:ImagesPerPerson)
        x = FaceData(p, i).Image;
        [irow, icol] = size(x);
        x = double(x);
        temp = reshape(x', irow * icol, 1);
        M = [M temp];
        Group = [Group p];
    end
end


%[eigenvectors, eigenvalues] = eigs(M*(M'), J);
%V = eigenvectors;


%subtract mean
x=bsxfun(@minus, M, mean(M,2));
% calculate covariance
s = cov(x');
% obtain eigenvalue & eigenvector
[V,D] = eigs(s,J);
Z = (x') * V;


figure; %oryginalne
nOfImages = Persons*ImagesPerPerson;
for i=(1:nOfImages)
    C = M(:,i);
    CC = reshape(C, [46, 56]);
    subplot(round(sqrt(nOfImages)), round(sqrt(nOfImages)) + 1, i);
    imagesc(CC');
    title(i)
    colormap gray;
end

figure; %eigenfaces
for i=(1:J)
    C = V(:,i);
    CC = reshape(C, [46, 56]);
    subplot(round(sqrt(J)), round(sqrt(J)) + 1, i);
    imagesc(CC');
    title(i)
    colormap gray;
end

id = (1:J);

kmeans_result = kmeans(Z', Persons);
groups = sortrows([id', kmeans_result], 2);

known_groups = [1 2 3 2 1 2 3 2 3 1];

confusionmat_result = confusionmat(known_groups, kmeans_result);

figure;
plotconfusion(known_groups, kmeans_result')

[acc_eigen, rand_index_eigen, match_eigen] = AccMeasure(known_groups, kmeans_result);

kmeans_result_original = kmeans(M', Persons);
[acc_orig, rand_index_orig, match_orig] = AccMeasure(Group, kmeans_result_original);

Class = knnclassify(V', M', Group);

figure;
plotconfusion(known_groups, Class')


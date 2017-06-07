clc;
clear;
close all;

load('FaceData_56_46.mat');

Persons = 3;
ImagesPerPerson = 10;

%wczytywanie danych
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

figure;
suptitle('Twarze oryginalne');
nOfImages = Persons*ImagesPerPerson;
for i=(1:nOfImages)
    C = M(:,i);
    CC = reshape(C, [46, 56]);
    subplot(round(sqrt(nOfImages)), round(sqrt(nOfImages)) + 1, i);
    imagesc(CC');
    title(i)
    colormap gray;
end

for group_size=(1:3)
    tic;
    
    % ksrednich dla obrazow oryginalnych
    kmeans_result_orig = kmeans(M', group_size);
    [acc_orig, rand_index_orig, match_orig] = AccMeasure(Group, kmeans_result_orig');

    disp(sprintf('Czas grupowania dla obrazow oryginalnych n=%d: %2.3fs',group_size, toc))
    disp(sprintf('Dokladnosc grupowania dla obrazow oryginalnych n=%d: %2.2f',group_size, acc_orig))
end

% calculate covariance
s = cov(M');

J_serie = [4 10 20 30];
for J_current=(1:length(J_serie))
    J = J_serie(J_current);

    % obtain eigenvalue & eigenvector
    [V,D] = eigs(s,J);

    Z = (M') * V;

    figure;
    suptitle('Twarze w³asne');
    for i=(1:J)
        C = V(:,i);
        CC = reshape(C, [46, 56]);
        subplot(round(sqrt(J)), round(sqrt(J)) + 1, i);
        imagesc(CC');
        title(i)
        colormap gray;
    end

    known_groups = ones(1,J);

    tic;
    
    % ksrednich dla obrazow redukowanych
    kmeans_result_eigen = kmeans(Z', Persons);
    [acc_eigen, rand_index_eigen, match_eigen] = AccMeasure(known_groups, kmeans_result_eigen');
    
    disp(sprintf('Czas grupowania dla obrazow redukowanych J=%d: %2.3fs',J, toc))
    disp(sprintf('Dokladnosc grupowania dla obrazow redukowanych J=%d: %2.2f',J, acc_eigen))


    Class = knnclassify(V', M', Group);
    [acc_classify, rand_index_classify, match_classify] = AccMeasure(known_groups, Class');

    figure;
    plotconfusion(known_groups, Class')

end

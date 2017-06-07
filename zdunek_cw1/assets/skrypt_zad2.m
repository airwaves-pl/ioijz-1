clc;
clear;
close all;

load('FaceData_56_46.mat');

Persons = 3;
ImagesPerPerson = 10;
J_serie = [4 10 20 30];

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

for J_current=(1:length(J_serie))
    J = J_serie(J_current);

    %subtract mean
    x=bsxfun(@minus, M, mean(M,2));
    % calculate covariance
    s = cov(x');
    % obtain eigenvalue & eigenvector
    [V,D] = eigs(s,J);
    Z = (x') * V;

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
    
    % ksrednich dla obrazow oryginalnych
    kmeans_result_orig = kmeans(M', Persons);
    [acc_orig, rand_index_orig, match_orig] = AccMeasure(Group, kmeans_result_orig');
    
    disp(sprintf('Czas grupowania dla obrazow oryginalnych J=%d: %2.3fs',J, toc))
    disp(sprintf('Dokladnosc grupowania dla obrazow oryginalnych J=%d: %2.2f',J, acc_orig))
    
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

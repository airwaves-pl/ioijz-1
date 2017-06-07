clc;
clear;
close all;

% parametry
Persons = 3;
ImagesPerPerson = 10;
NumberOfGroups = 3;
J_serie = [4 10 20 30];

% wczytywanie danych do macierzy wejsciowej
load('FaceData_56_46.mat');
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

M2 = [];
for p=(1:Persons)
    for i=(1:ImagesPerPerson)
        x = FaceData(p, i).Image;
        [irow, icol] = size(x);
        x = double(x);
        temp = reshape(x', irow * icol, 1);
        % obliczenie wartosci srednich
        % odjecie ich od macierzy
        temp = bsxfun(@minus, temp, mean(temp, 1));
        M2 = [M2 temp];
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

% M2 - zawiera usrednione obrazy w formie macierzy
% wyznaczenie macierzy kowariancji
s = cov(M2');

for g=(3:NumberOfGroups)
    tic;
    kmeans_original = kmeans(M', NumberOfGroups);
    disp(sprintf('Czas grupowania dla obrazow oryginalnych G=%d: %2.3fs', g, toc))
    [acc_orig, rand_index_orig, match_orig] = AccMeasure(Group, kmeans_original');
    disp(sprintf('Dokladnosc grupowania do obrazow oryginalnych G=%d: %2.2f', g, acc_orig))
end


knnclassify_org = knnclassify(M', M', Group);
[acc_redknno, rand_index_redknno, match_redknno] = AccMeasure(Group, knnclassify_org');
disp(sprintf('Dokladnosc klasyfikacji dla obrazow oryginalnych: %2.2f', acc_redknno))

for J_current=(1:length(J_serie))
    J = J_serie(J_current);
    % wyznaczenie wartosci i wektorow wlasnych
    [V, D] = eigs(s, J);
    % porzadkujemy wartosci wlasne w kolejnosci malejacej
    [D order] = sort(diag(D), 'descend');
    V = V(:,order);
    % rzutowanie
    Z = V' * M;

    %figure;
    %suptitle('Twarze w³asne');
    %for i=(1:J)
    %    C = V(:,i);
    %    CC = reshape(C, [46, 56]);
    %    subplot(round(sqrt(J)), round(sqrt(J)) + 1, i);
    %    imagesc(CC');
    %    title(i)
    %    colormap gray;
    %end

    figure;
    suptitle(sprintf('Twarze zredukowane J=%d', J));
    nOfImages = Persons*ImagesPerPerson;
    for i=(1:nOfImages)
        C = Z(:,i)' * V(:,1:J)';
        CC = reshape(C, [46, 56]);
        subplot(round(sqrt(nOfImages)), round(sqrt(nOfImages)) + 1, i);
        imagesc(CC');
        title(i)
        colormap gray;
    end
    % Porownanie wynikow zredukowanych
    for g=(3:NumberOfGroups)
        tic;
        kmeans_reduced = kmeans(Z', NumberOfGroups);
        disp(sprintf('Czas grupowania dla obrazow zredukowanych G=%d, J=%d: %2.3fs', g, J, toc))
        [acc_red, rand_index_red, match_red] = AccMeasure(Group, kmeans_reduced');
        disp(sprintf('Dokladnosc grupowania do obrazow zredukowanych G=%d, J=%d: %2.2f', g, J, acc_red))
    end
    
    knnclassify_red = knnclassify(Z' * V', M', Group);
    [acc_redknn, rand_index_redknn, match_redknn] = AccMeasure(Group, knnclassify_red');
    disp(sprintf('Dokladnosc klasyfikacji dla obrazow zredukowanych J=%d: %2.2f', J, acc_redknn))
    
end

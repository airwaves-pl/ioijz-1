clc;
clear;
close all;
load('FaceData_56_46.mat');

Persons = 8;
ImagesPerPerson = 10;
nOfImages = Persons*ImagesPerPerson;

% wczytanie danych do tensora
P = zeros(1,nOfImages);
Y=zeros(56,46,nOfImages);
img_index = 1;
for p=(1:Persons)
    for i=(1:ImagesPerPerson)
        P(img_index) = p;
        Y(:,:,img_index) = FaceData(p, i).Image;
        img_index = img_index + 1;
    end
end
P = P';

figure;
suptitle('Twarze oryginalne');
for i=(1:nOfImages)
    subplot(Persons, ImagesPerPerson, i);
    imagesc(Y(:,:,i));
    title(i)
    colormap gray;
    set(gca,'XtickLabel',[],'YtickLabel',[]);
end

% rozdzielenie na dwa zbiory (5-folds CV)
CV = cvpartition(P,'kfold',5);
train_idx = CV.training(1);
test_idx = CV.test(1);

% utworzenie tensorow trenujacego i testowego
Y_train = Y(:,:,train_idx);
Y_test = Y(:,:,test_idx);
Class_train_idx = P(train_idx);
Class_test_idx = P(test_idx);

J_serie = [4 10 20 30];

res_rands = zeros(1,length(J_serie));
res_time_all = zeros(1,length(J_serie));
res_time_train = zeros(1,length(J_serie));
res_acc = zeros(1,length(J_serie));
res_delta = zeros(1,length(J_serie));
res_groups_kmeans = [];

for J_current=(1:length(J_serie))
    J(1:3) = J_serie(J_current);

    % dekompozycja hosvd (pod kmeansa)
    tic
    [A, B, C, G, Y_hat] = skrypt_zad3_hosvd(Y, J);
    res_time_all(J_current) = toc;

    figure;
    suptitle(sprintf('Twarze zredukowane J=%d (HOSVD)', J_serie(J_current)));
    for i=(1:nOfImages)
        subplot(Persons, ImagesPerPerson, i);
        imagesc(Y_hat(:,:,i));
        title(i)
        colormap gray;
        set(gca,'XtickLabel',[],'YtickLabel',[]);
    end

    % grupowanie metoda ksrednich dla faktora U^(3) - stala liczba grup (ilosc osob)
    kmeans_result = kmeans(C, Persons);
    res_groups_kmeans = [res_groups_kmeans kmeans_result];
    [res_acc(J_current), res_rands(J_current), ~] = AccMeasure(P, kmeans_result');
    
    % dekompozycja hosvd
    tic
    [Ar, Br, Cr, Gr, Yr_hat] = skrypt_zad3_hosvd(Y_train, J);
    res_time_train(J_current) = toc;

    % projekcja
    Y3 = reshape(permute(Y_test,[3 1 2]),size(Y_test,3),size(Y_test,1)*size(Y_test,2));
    G3 = reshape(permute(Gr,[3 1 2]),[J(3),J(1)*J(2)]);
    Ct = Y3*pinv(double(G3)*(kron(Br,Ar))');
    Ct = Ct.*repmat(1./sqrt(sum(Ct.^2,2)+eps),1,size(Ct,2));

    % klasyfikacja w przestrzeni cech U^(3)
    mdl_class = fitcknn(Cr,Class_train_idx,'NumNeighbors',1);
    prediction = predict(mdl_class, Ct); 
    
    % dokladnosc klasyfikacji
    res_delta(J_current) = 100*(length(find((prediction - Class_test_idx)==0))/length(Class_test_idx));

end


load('FaceData_56_46.mat');

Persons = 10;
ImagesPerPerson = 10;

%wczytywanie danych
Group = [];
Tensor=zeros(56,46,Persons*ImagesPerPerson);
img_index = 1;
for p=(1:Persons)
    for i=(1:ImagesPerPerson)
        x = FaceData(p, i).Image;
        [irow, icol] = size(x);
        x = double(x);
        temp = reshape(x', irow * icol, 1);
        Group = [Group p];
        Tensor(:,:,img_index)=x;
        img_index = img_index + 1;
    end
end


% 5-folds CV

Y = Group';

CV = cvpartition(Y,'kfold',5);
train_idx = CV.training(1);
test_idx = CV.test(1);

Y_train = Tensor(:,:,train_idx);
Y_test = Tensor(:,:,test_idx);

Class_train_inx = Y(train_idx);
Class_test_inx = Y(test_idx);

Jt = [4 4 4];

% Orth-Tucker
[delta_hosvd,elapsed_time_hosvd] = Tucker_orth_3D(Y_train,Y_test,Class_train_inx,Class_test_inx,Jt);


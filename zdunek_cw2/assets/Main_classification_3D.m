
load('FaceData_56_46.mat');

Persons = 6;
ImagesPerPerson = 10;

%wczytywanie danych
Group = [];
Tensor=zeros(56,46,Persons*ImagesPerPerson);
for p=(1:Persons)
    for i=(1:ImagesPerPerson)
        x = FaceData(p, i).Image;
        [irow, icol] = size(x);
        x = double(x);
        temp = reshape(x', irow * icol, 1);
        Group = [Group p];
        Tensor(:,:,p)=x;
    end
end


% tutaj trzeba foldowac Tensor na 2 zbiory

% tymczasowo
Y_train = Tensor;
Y_test = Tensor;

Class_train_inx = Group';
Class_test_inx = Group';

Jt = [4 4 4];

% Orth-Tucker
[C_hosvd,delta_hosvd,elapsed_time_hosvd] = Tucker_orth_3D(Y_train,Y_test,Class_train_inx,Class_test_inx,Jt);

% Visualization
figure
hintonw((squeeze(mean(C_hosvd,3)))')
title(['OrthTucker: P = ',num2str(mean(delta_hosvd) ), ' %'])
ylabel('Output')


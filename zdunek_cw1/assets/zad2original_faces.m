clear;
clc;
load FaceData_56_46;
Persons = 3
ImagesPerPerson = 10
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
figure;
for i=(1:Persons*ImagesPerPerson)
    C = M(:,i);
    CC = reshape(C, [46, 56]);
    subplot(round(sqrt(Persons*ImagesPerPerson)), round(sqrt(Persons*ImagesPerPerson)) + 1, i);
    imagesc(CC');
    title(i)
    colormap gray;
end

id = (1:Persons*ImagesPerPerson);
result = kmeans(M', Persons);
result = [id', result]
%sortrows(result, 2)
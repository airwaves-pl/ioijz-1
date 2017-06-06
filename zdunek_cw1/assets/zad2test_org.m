
clear;
clc;
load FaceData_56_46
Persons = 3;
ImagesPerPerson = 10;
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
n=Persons*ImagesPerPerson;
for i=(1:n)
    subplot(round(sqrt(n)), round(sqrt(n)) + 1, i);
    imagesc(reshape(M(:,i), [46, 56])');
    title(i)
    colormap gray;
end

id = (1:n);
result = kmeans(M', Persons);
result = [id', result];
sortrows(result, 2)


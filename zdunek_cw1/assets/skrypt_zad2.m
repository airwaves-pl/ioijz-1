
image(FaceData(2,2).Image)

%C = confusionmat(group,grouphat);

X=[randn(200,2);randn(200,2)+6,;[randn(200,1)+12,randn(200,1)]]; T=[ones(200,1);ones(200,1).*2;ones(200,1).*3]; 
 idx=kmeans(X,3,'emptyaction','singleton','Replicates',5); 
  [Acc,rand_index,match]=AccMeasure(T,idx)
  
close all;
load yalefaces
[h,w,n] = size(yalefaces);
d = h*w;
% vectorize images
x = reshape(yalefaces,[d n]);
x = double(x);
%subtract mean
x=bsxfun(@minus, x, mean(x,2));
% calculate covariance
s = cov(x');
% obtain eigenvalue & eigenvector
[V,D] = eig(s);
eigval = diag(D);
% sort eigenvalues in descending order
eigval = eigval(end:-1:1);
V = fliplr(V);
% show 0th through 15th principal eigenvectors
eig0 = reshape(mean(x,2), [h,w]);
figure,subplot(4,4,1)
imagesc(eig0)
colormap gray
for i = 1:15
subplot(4,4,i+1)
imagesc(reshape(V(:,i),h,w))
end
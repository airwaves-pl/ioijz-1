
image(FaceData(2,2).Image)

%C = confusionmat(group,grouphat);

X=[randn(200,2);randn(200,2)+6,;[randn(200,1)+12,randn(200,1)]]; T=[ones(200,1);ones(200,1).*2;ones(200,1).*3]; 
 idx=kmeans(X,3,'emptyaction','singleton','Replicates',5); 
  [Acc,rand_index,match]=AccMeasure(T,idx)
  
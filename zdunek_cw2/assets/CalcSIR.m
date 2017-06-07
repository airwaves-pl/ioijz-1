function SIR = CalcSIR(A,Aest)
% Sergio Cruces & Andrzej Cichocki 
A=A*diag(1./(sqrt(sum(A.^2))+eps));
Aest=Aest*diag(1./(sqrt(sum(Aest.^2))+eps));
for i=1:size(Aest,2)
  [MSE(i),ind]=min(sum(bsxfun(@minus,Aest(:,i),A).^2,1));
%  A(:,ind)=[];
end
SIR=-10*log10(MSE); 

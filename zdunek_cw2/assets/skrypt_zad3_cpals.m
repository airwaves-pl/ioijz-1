function [accuracy, time, C] = skrypt_zad3_cpals(Ytrain, Ytest, TrainIndex, TestIndex, J)
% Otrzymujemy faktory, u nas beda tylko A, B, C, a G to tensor rdzeniowy
tic
%% 
%Na podstawie tensor_toolbox
%norma tensora
X = tensor(Ytrain);
tensorNorm = norm(X);
%maksymalna liczba iteracji
maxiters = 50;

%tworzymy poczatkowe losowe wartosci dla macierzy
U = cell(3,1);
U{1} = rand(size(X, 1), J);
U{2} = rand(size(X, 2), J);
U{3} = rand(size(X, 3), J);

fit = 0;

% Iteracje do momentu stagnacji
for k = 1:maxiters

    oldFit = fit;

    % iteruje wszystkie 3 mody tensora
    for n = 1:3

        % Calculate newU = X_(n) * khatrirao(all U except n, 'r').
        newU = mttkrp(X, U, n);

        % macierz wspolczynnikow
        Y = ones(J,J);
        for i = [1:n-1, n+1:3]
            Y = Y .* (U{i}'*U{i});
        end

        newU = (Y \ newU')'; 
        lambda = max( max(newU, [], 1), 1 )'; %norma
        U{n} = newU * spdiags(1./lambda,0 ,J ,J);

    end

    Tensor = ktensor(lambda, U); %tensor kruskala (na podstawie tensor toolbox)
    residual = sqrt( tensorNorm^2 + norm(Tensor)^2 - 2 * innerprod(X,Tensor) );
    fit = 1 - (residual / tensorNorm); %ze wzoru

    %sprawdzam czy nastapila juz zbieznosc (mala zmiana miedzy iteracjami)
    if (abs(oldFit - fit) < eps)
        break;
    end        
end 
   
   A = Tensor.U{1};
   B = Tensor.U{2};
   C = Tensor.U{3};
   
   C = C.*repmat(1./sqrt(sum(C.^2,2)+eps), 1, size(C,2));
%%
time = toc;

Y3 = reshape(permute(Ytest, [3 1 2]), size(Ytest,3), size(Ytest,1) * size(Ytest,2));
Ct = Y3 * pinv((khatrirao(B,A))');

Ct = Ct.*repmat(1./sqrt(sum(Ct.^2,2) + eps), 1, size(Ct, 2));

% Ytest = Ytrain;
% TestIndex = TrainIndex;

mdl = fitcknn(C, TrainIndex, 'NumNeighbors', 1);

%predykcja
prediction = predict(mdl, Ct);

predictions = [];
predictions = [predictions; [prediction, TestIndex]];
count = squeeze(sum(predictions(:,1) == predictions(:,2))).';
accuracy = count * 100 / numel(predictions(:,1));

end
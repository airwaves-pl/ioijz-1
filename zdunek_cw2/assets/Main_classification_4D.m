
% tutaj trzeba wczytac obrazy do tensorow i je foldowac na 2 zbiory

Jt = 4;



% NMF
% 
% C(:,:,n,2) = C_cp; delta(n,2) = delta_cp; ET(n,2) = elapsed_time_cp; res.res_nmf = res_nmf;

%  % Nonnegative ALS-CP
%   [C_cp,delta_cp,elapsed_time_cp,res_cp_als] = CP_4D(Y_train,Y_test,Class_train_inx,Class_test_inx,J,MCRuns,Toi,MaxIter);
%   C(:,:,n,3) = C_cp; delta(n,3) = delta_cp; ET(n,3) = elapsed_time_cp; res.res_cp_als = res_cp_hals;

%  % Nonnegative HALS-CP
%   [C_cp,delta_cp,elapsed_time_cp,res_cp_hals] = HALS_CP_4D(Y_train,Y_test,Class_train_inx,Class_test_inx,J,MCRuns,Toi,MaxIter);
%   C(:,:,n,4) = C_cp; delta(n,4) = delta_cp; ET(n,4) = elapsed_time_cp; res.res_cp_hals = res_cp_hals;

  % Orth-Tucker
    [C_hosvd,delta_hosvd,elapsed_time_hosvd] = Tucker_orth_4D(Y_train,Y_test,Class_train_inx,Class_test_inx,Jt);
    C(:,:,n,5) = C_hosvd; delta(n,5) = delta_hosvd; ET(n,5) = elapsed_time_hosvd; 

%  %  Orth-Tucker(core)
%     [C_hosvd,delta_hosvd,elapsed_time_hosvd] = Tucker_orth_4D_core(Y_train,Y_test,Class_train_inx,Class_test_inx,Jtc);
%     C(:,:,n,6) = C_hosvd; delta(n,6) = delta_hosvd; ET(n,6) = elapsed_time_hosvd;

% % Multiclass-NMP
%     [C_cp,delta_cp, elapsed_time_cp ] = Multiclass_NMP_4D(Y_train,Y_test,Class_train_inx,Class_train_inx,J,MCRuns,MaxIter);
%     C(:,:,n,6) = C_cp; delta(n,6) = delta_cp; ET(n,6) = elapsed_time_cp;

%   end

     % Visualization
    if visualization==1
        
        figure
        subplot(2,3,1)
        hintonw((squeeze(mean(C(:,:,:,1),3)))')
        title(['PCA: P = ',num2str(mean(delta(:,1)) ), ' %'])
        ylabel('Output')
        subplot(2,3,2)
        hintonw((squeeze(mean(C(:,:,:,2),3)))')
        title(['NMP: P = ',num2str(mean(delta(:,2)) ), ' %'])
        ylabel('Output')
        subplot(2,3,3)
        hintonw((squeeze(mean(C(:,:,:,3),3)))')
        title(['CP(ADS): P = ',num2str(mean(delta(:,3)) ), ' %'])
        ylabel('Output')
        subplot(2,3,4)
        hintonw((squeeze(mean(C(:,:,:,4),3)))')
        title(['CP(HALS): P = ',num2str(mean(delta(:,4)) ), ' %'])
        ylabel('Output')
        
        hintonw((squeeze(mean(C(:,:,:,5),3)))')
        title(['OrthTucker: P = ',num2str(mean(delta(:,5)) ), ' %'])
        ylabel('Output')
        
    end

% pokaza³ diagramy hintona

     
     
% ---- to co napisal na konsoli ----

size(Y_train)
length(Class_train_inx) %ans=1760


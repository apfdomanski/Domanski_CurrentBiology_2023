function FracCorr=DecodePredErrCrossTemporal_CV_GPU_old(FR,evt0,reg)
%
% (c) 2020 Aleksander PF Domanski PhD UoB
% aleks.domanski@bristol.ac.uk
% Input Variables:
% FR:   Firing rate matrix [units,time x no. trials]
% evt0: Event category labels (either 1 or 2) [no. trials]
% reg:  Regularization parameter
%
% Output Variables:
% PE:  Decoding result
%
if nargin<3, reg=0; end;
FoldCVE = 2;
nDraws = 10;
Ltr=round(size(FR,1)/length(evt0));
PE = nan(Ltr,Ltr,nDraws);
n_(1) = sum(evt0==1);
n_(2) = sum(evt0==2);

DrawPool = min(ceil(n_/FoldCVE)); % how many trials to draw to train on?
for s =1:2
    for iDraw = 1:nDraws
        idTrain{s}{iDraw} = randsample(n_(s),DrawPool);
        idTest{s}{iDraw}  = setdiff(1:n_(s), idTrain{s}{iDraw});
    end
end



for iDraw =1:nDraws
    Creg_ = cell(1,Ltr);
    m1_ = Creg_;
    m2_ = Creg_;

    fprintf('Cross-temporal decoding (draw %d/%d)...\n',iDraw,nDraws)
    PE_ = zeros(Ltr,Ltr); 
    N_=0;
    
    
    for t_Train=1:Ltr % Loop across each training time bin
        X=FR(t_Train:Ltr:end,:);
        Xtrain1=X(evt0==1,:);
        Xtrain2=X(evt0==2,:);

        Xtrain1=Xtrain1(idTrain{1}{iDraw},:);
        Xtrain2=Xtrain2(idTrain{2}{iDraw},:);
        
        Xtrain1 = Xtrain1(sum(isnan(Xtrain1),2)==0,:);
        Xtrain2 = Xtrain2(sum(isnan(Xtrain2),2)==0,:);
        
        Cpool=((DrawPool-1)*cov(Xtrain1)+(DrawPool-1)*cov(Xtrain2))./(2*DrawPool-2);            % Pooled covariance matrix of 2 response sets
        Creg=Cpool+reg*eye(size(Cpool));                                                        % Regularized (pooled) covariance matrix
        m1=mean(Xtrain1); m2=mean(Xtrain2);                                                     % Mean at this time point
        
        Creg_{t_Train} = Creg^-1;
        m1_{t_Train}=m1;
        m2_{t_Train}=m2;
    end
    Creg_ = repmat(Creg_,Ltr,1);
    m1_   = repmat(m1_,Ltr,1);
    m2_   = repmat(m2_,Ltr,1);
    
    
  
       
        
    
    for s=1:2 % Loop across each outcome type
        Xtest_ = cell(Ltr,1);
        n_     = cell(Ltr,1);
        for t_Test=1:Ltr % Loop across each time bin
            
            X_ = FR(t_Test:Ltr:end,:);
            Xtest=X_(evt0==s,:);
            Xtest=Xtest(idTest{s}{iDraw},:);
            Xtest = Xtest(sum(isnan(Xtest),2)==0,:);
            n=size(Xtest,1);
            
            Xtest_{t_Test} = Xtest;
            n_{t_Test} = n;
        end
        Xtest_ = repmat(Xtest_,1,Ltr);
        n_ = repmat(n_,1,Ltr);
        
        % Mah. distance from response 1 mean at this time point
        A  = cellfun(@(Xtest,n,m1) (Xtest-ones(n,1)*m1), Xtest_,n_, m1_,'UniformOutput',false);
        B  = cellfun(@(Xtest,n,m1) (Xtest-ones(n,1)*m1)', Xtest_,n_, m1_,'UniformOutput',false);
        D1 = cellfun(@(A,Creg,B) diag(A*Creg*B), A,Creg_,B,'UniformOutput',false);
        
        % Mah. distance from response 2 mean at this time point
        A  = cellfun(@(Xtest,n,m1) (Xtest-ones(n,1)*m1), Xtest_,n_, m2_,'UniformOutput',false);
        B  = cellfun(@(Xtest,n,m1) (Xtest-ones(n,1)*m1)', Xtest_,n_, m2_,'UniformOutput',false);
        D2 = cellfun(@(A,Creg,B) diag(A*Creg*B), A,Creg_,B,'UniformOutput',false);
        
        % Predicted class at this time point
        cpred = cellfun(@(D1,D2)(D1>D2)+1,D1,D2,'UniformOutput',false);
        % Classifcation error at this time point
        PE_ = PE_ + cell2mat(cellfun(@(cpred) sum(cpred~=s),cpred,'UniformOutput',false));
        N_ = N_+n;
   
    end
    
    PE(:,:,iDraw)=PE_./N_;                                                         % Normalise by length
end
FracCorr=1-nanmean(PE,3);

















































function [FracCorr FracCorrShuffCIL FracCorrShuffCIH]=DecodePredErrCrossTemporal_CV(FR,evt0,reg)
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
if nargout>1
    ShuffleYN = true;
else
    ShuffleYN = false;
end
   
if nargin<3, reg=0; end
FoldCVE = 2;
nDraws = 10;
Ltr=round(size(FR,1)/length(evt0));
PE = nan(Ltr,Ltr,nDraws);
n_(1) = sum(evt0==1);
n_(2) = sum(evt0==2);

DrawPool = min(ceil(n_/FoldCVE)); % how many trials to draw to train on?

%Pre-draw Training and Test sets of trials
for s =1:2
    for iDraw = 1:nDraws
        idTrain{s}{iDraw} = randsample(n_(s),DrawPool);
        idTest{s}{iDraw}  = setdiff(1:n_(s), idTrain{s}{iDraw});
    end
end




for iDraw =1:nDraws
    fprintf('Cross-temproal decoding (draw %d/%d)...\n',iDraw,nDraws)
    PE_=zeros(Ltr,Ltr); N_=PE_;
    
    
    parfor tTrain=1:Ltr % Loop across each time bin

        X=FR(tTrain:Ltr:end,:);
        Xtrain1=X(evt0==1,:);
        Xtrain2=X(evt0==2,:);

        Xtrain1=Xtrain1(idTrain{1}{iDraw},:);
        Xtrain2=Xtrain2(idTrain{2}{iDraw},:);
        
        Xtrain1 = Xtrain1(sum(isnan(Xtrain1),2)==0,:);
        Xtrain2 = Xtrain2(sum(isnan(Xtrain2),2)==0,:);
        
        Cpool=((DrawPool-1)*cov(Xtrain1)+(DrawPool-1)*cov(Xtrain2))./(2*DrawPool-2);% Pooled covariance matrix of 2 response sets
        Creg=Cpool+reg*eye(size(Cpool));                                       % Regularized (pooled) covariance matrix
        m1=mean(Xtrain1); m2=mean(Xtrain2);                                    % Mean at this time point
        
        
        for tTest=1:Ltr % Loop across each time bin
            X_ = FR(tTest:Ltr:end,:);
            for s=1:2 % Loop across each outcome type
                
                Xtest=X_(evt0==s,:);
                Xtest=Xtest(idTest{s}{iDraw},:);
                Xtest = Xtest(sum(isnan(Xtest),2)==0,:);
                
                n=size(Xtest,1);
                
                D1=diag((Xtest-ones(n,1)*m1)*Creg^-1*(Xtest-ones(n,1)*m1)');       % Mah. distance from response 1 mean at this time point
                D2=diag((Xtest-ones(n,1)*m2)*Creg^-1*(Xtest-ones(n,1)*m2)');       % Mah. distance from response 2 mean at this time point
                cpred=(D1>D2)+1;                                                   % Predicted class at this time point
                PE_(tTrain,tTest)=PE_(tTrain,tTest)+sum(cpred~=s);                 % Classifcation error at this time point
                N_(tTrain,tTest)=N_(tTrain,tTest)+n;
                
                
            end
            
        end
    end
    
    PE(:,:,iDraw)=PE_./N_;                                                         % Normalise by length
end
FracCorr=1-nanmean(PE,3);

if ShuffleYN
    fprintf('Cross-temproal decoding: TRial shuffled data (draw %d/%d)...\n',iDraw,nDraws)
    PE_=zeros(Ltr,Ltr); N_=PE_;
    
    
    
end





































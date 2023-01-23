function [FracCorr,FracCorr_CIshuffLow,FracCorr_CIshuffHigh]=DecodePredErrCrossTemporal_CV_GPU(FR,evt0,reg,runShuffle)
%
% (c) 2020 Aleksander PF Domanski PhD UoB
% aleks.domanski@bristol.ac.uk
%
% Cross-temporal regularised leave-one-out cross-validated decoder
% - Trains and tests at each time-point using separate training  and test sets
% - Uses count-matched subsets of trials to normalise decoding performance
% - Optionally calculates bootstrap permutated confidence bounds on decoding performance
%
% Input Variables:
% FR:   Firing rate matrix [units,time x no. trials]
% evt0: Event category labels (either 1 or 2) [no. trials]
% reg:  Regularization parameter
% runShuffle: Trial permutation for stats calculation (true or false)
%
% Output Variables:
% FracCorr:  Decoding result
% FracCorr_CIshuffLow:  Trial-label shuffled decoding result (lower bound)
% FracCorr_CIshuffHigh:  Trial-label shuffled decoding result (upper bound)
%
if nargin<3 || isempty(reg), reg=0.05; end;
if nargin<4 || isempty(runShuffle), runShuffle = false; end
FoldCVE = 2;                            % witheld fraction of trials for cross-validation testing
nDraws = 10;                            % how many trial draws
Ltr=round(size(FR,1)/length(evt0));     % Trial length (samples)
n_(1) = sum(evt0==1);                   % No. trials (outcome 1)
n_(2) = sum(evt0==2);                   % No. trials (outcome 2)

PE = nan(Ltr,Ltr,nDraws);               % Prediction error result at each train/test timestep pair
if runShuffle
    PEshufHigh = nan(Ltr,Ltr,nDraws);   % Prediction error result: 95% decoding performance
    PEshufLow  = nan(Ltr,Ltr,nDraws);   % Prediction error result: 95% decoding performance
end

DrawPool = min(ceil(n_/FoldCVE)); % how many trials to draw to train on
nBS = round(100/nDraws); % 11/21: was 1000
for s =1:2
    for iDraw = 1:nDraws
        idTrain{s}{iDraw} = randsample(n_(s),DrawPool);
        idTest{s}{iDraw}  = setdiff(1:n_(s), idTrain{s}{iDraw});
    end
end

for iDraw =1:nDraws
    fprintf('Cross-temporal decoding (trial draw %d/%d)...\n',iDraw,nDraws)
    % Train the decoder
    % m1_, m2_ are mean firing rates of each unit for conditions 1 and 2 at each training time-step.
    % Creg_ is the regularised pooled covariance of unit firing rates of units at each training time-step.
    [m1_,m2_,Creg_] = TrainFun(FR,reg,Ltr,evt0,DrawPool,false,idTrain,iDraw);
    PE(:,:,iDraw)   = TestFun(FR,Ltr,evt0,false,m1_,m2_,Creg_,idTest,iDraw);
    
    
    if runShuffle
        PEshuf = nan(Ltr,Ltr,nBS);               % Prediction error result: 
        parfor bs =1:nBS
            fprintf('Cross-temporal decoding on shuffled labels (trial draw %d/%d, BS no. %d/%d)...\n',iDraw,nDraws,bs,nBS)
            
            [m1_,m2_,Creg_] = TrainFun(FR,reg,Ltr,evt0,DrawPool,true,idTrain,iDraw);
            PEshuf(:,:,bs)   = TestFun(FR,Ltr,evt0,true,m1_,m2_,Creg_,idTest,iDraw);
        end
            PEshuf = sort(PEshuf,3,'ascend');
            PEshufHigh(:,:,iDraw) = PEshuf(:,:,round(0.95*nBS));
            PEshufLow(:,:,iDraw)  = PEshuf(:,:,round(0.05*nBS));
    end
end
FracCorr             = 1-nanmean(PE,3);
if runShuffle
    FracCorr_CIshuffHigh = 1-nanmean(PEshufHigh,3);
    FracCorr_CIshuffLow  = 1-nanmean(PEshufLow,3);
else 
    FracCorr_CIshuffHigh = nan(size(FracCorr));
    FracCorr_CIshuffLow  = nan(size(FracCorr));
end
end


function [m1_,m2_,Creg_] = TrainFun(FR,reg,Ltr,evt0,DrawPool,shuffleYN,idTrain,iDraw)
if shuffleYN, evt0 = evt0(randperm(length(evt0))); end
Creg_ = cell(1,Ltr);
m1_ = Creg_;
m2_ = Creg_;

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
end

function PE_ = TestFun(FR,Ltr,evt0,shuffleYN,m1_,m2_,Creg_,idTest,iDraw)
if shuffleYN, evt0 = evt0(randperm(length(evt0))); end
PE_ = zeros(Ltr,Ltr);
N_= 0;
for s=1:2 % Loop across each outcome type
    Xtest_ = cell(Ltr,1);
    n_     = cell(Ltr,1);
    for t_Test=1:Ltr % Loop across each time bin
        
        X_    = FR(t_Test:Ltr:end,:);
        Xtest = X_(evt0==s,:);
        Xtest = Xtest(idTest{s}{iDraw},:);
        Xtest = Xtest(sum(isnan(Xtest),2)==0,:);
        n     = size(Xtest,1);
        
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
PE_ = PE_./N_; % Nomalise by length;
end


function [CVE,...                time-dependent mean decoding score [No. Timepoints x 1]
    Score,...              Time-averaged/Peak decoding score, depending on InfoCriteria [1 x 1]
    Score_Shuffled,...     Median Time-averaged decoding score from shuffled trials, [1 x 1]
    Score_ShuffledCIL,...  5 percentile Time-averaged decoding score from shuffled trials, [1 x 1]
    Score_ShuffledCIH,...  95 percentile Time-averaged decoding score from shuffled trials, [1 x 1]
    CVE_shuffledHL] ...    mean time-dependent decoding score [No. Timepoints x 2]
        = RunDecodeGroupsDestroyNoiseCorrs( unitIds , nReps , Ltrials_ , Rtrials_ , TrialsSubsample , InfoCriteria , timeRange)
%%
% if isinf(TrialsSubsample)
%     nReps=1;
% end
    
plotOnline = false;
LTR              = size(Ltrials_{1},1);
CVE              = nan(LTR,nReps);
Score            = nan(nReps,1);
if nargin<7
    timeRange=1:LTR;
end
if nargout>2
    PermuteTrialLabels=true;
    Score_Shuffled    = nan(nReps,1);
    Score_ShuffledCIL = nan(nReps,1);
    Score_ShuffledCIH = nan(nReps,1);
    CVE_shuffledHL=[];
    CVE_shuffledH = nan(LTR,nReps);
    CVE_shuffledL = nan(LTR,nReps);
else
    PermuteTrialLabels=false;
    Score_Shuffled    = NaN;
    Score_ShuffledCIL = NaN;
    Score_ShuffledCIH = NaN;
    CVE_shuffledHL=[];
    CVE_shuffledH = NaN;
    CVE_shuffledL = NaN;
end
if nReps > 1
    nBS = 5;
else
    nBS = 50;
end

% Draw sub-groups of trials to maintain DoF between experiments

parfor nrep = 1:nReps
    CVE_=[]; CVEshuf_=[];
    
    if isinf(TrialsSubsample)
        rsL = 1:length(Ltrials_);
        rsR = 1:length(Rtrials_);
        evt0  = [ones(length(Ltrials_),1);2*ones(length(Rtrials_),1)];
    else
        rsL = randsample(length(Ltrials_),TrialsSubsample);
        rsR = randsample(length(Rtrials_),TrialsSubsample);
        evt0  = [ones(TrialsSubsample,1);2*ones(TrialsSubsample,1)];
    end
    
    % Trial Label permutation step
    FR_L = [];  FR_R = [];
    for iUnit =1:length(unitIds)
        unitID = unitIds(iUnit);
        
            rsL_ = rsL(randperm(length(rsL)));
            rsR_ = rsR(randperm(length(rsR)));
       
        for iTrial = 1:length(rsL)
            FR_L{iTrial}(:,iUnit) = Ltrials_{rsL_(iTrial)}(:,unitID);
        end
        for iTrial = 1:length(rsR)
            FR_R{iTrial}(:,iUnit) = Rtrials_{rsR_(iTrial)}(:,unitID);
        end
        
    end
    FR_   = [cell2mat(FR_L');cell2mat(FR_R')];
    
    
    % (1) Decode real
    if exist('DecodeCVE_mex','file')
        CVE_ = 1-DecodeCVE_mex(FR_,evt0,0.05);
    else
        CVE_ = 1-DecodeCVE(FR_,evt0,0.05);
    end
    CVE(:,nrep) = CVE_;
    switch InfoCriteria
        case 'mean'
            Score(nrep) = mean(CVE_(timeRange),2);
        case 'max'
            Score(nrep) = max(CVE_(timeRange),[],2);
    end
    % (2) Decode shuffled versions - aim for 50 total random draws
    if PermuteTrialLabels
        for iDraw = 1:nBS
            if exist('DecodeCVE_mex','file')
                CVEshuf_(iDraw,:) = 1-DecodeCVE_mex(FR_,evt0(randperm(length(evt0))),0.05);
            else
                CVEshuf_(iDraw,:) = 1-DecodeCVE(FR_,evt0(randperm(length(evt0))),0.05);
            end
        end
        Score_Shuffled(nrep)      =  nanmean(median(CVEshuf_(:,timeRange),1));
        Score_ShuffledCIL(nrep)   =  prctile(median(CVEshuf_(:,timeRange),2),5);
        Score_ShuffledCIH(nrep)   =  prctile(median(CVEshuf_(:,timeRange),2),95);
        CVE_shuffledL(:,nrep) =  prctile(CVEshuf_,5,1)';
        CVE_shuffledH(:,nrep) =  prctile(CVEshuf_,95,1)';
    end
end

Score    = nanmean(Score);
CVE      = nanmean(CVE,2);
if PermuteTrialLabels
    Score_Shuffled    =  nanmean(Score_Shuffled);
    Score_ShuffledCIL =  nanmean(Score_ShuffledCIL);
    Score_ShuffledCIH =  nanmean(Score_ShuffledCIH);
    CVE_shuffledHL    =  [nanmean(CVE_shuffledL,2), nanmean(CVE_shuffledH,2)];
else
    
end

if plotOnline
    figure; hold on
    plot(CVE,'k')
    plot(CVE_shuffledHL,'r')
    plot([1 LTR],[Score Score],'k')
    plot([1 LTR],[Score_Shuffled Score_Shuffled],'r')
    plot([1 LTR],[Score_ShuffledCIH Score_ShuffledCIH],':r')
    plot([1 LTR],[Score_ShuffledCIL Score_ShuffledCIL],':r')
end

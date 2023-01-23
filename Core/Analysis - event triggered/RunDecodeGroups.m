function [CVE,...                
    Score,...              
    Score_Shuffled,...     
    Score_ShuffledCIL,...  
    Score_ShuffledCIH,...  
    CVE_shuffledHL] ...    
        = RunDecodeGroups( unitIds , nReps , Ltrials_ , Rtrials_ , TrialsSubsample , InfoCriteria , timeRange)
% Wrapper function for DecodeCVE.m: Handles trial subsampling and summary statistics
%
%  ***Inputs:
%  unitIds              Indices of each to include
%  nReps                How many jackknife samples to run for shuffled outcome distributions
%  Ltrials_             Input Data for Left trials (i.e. outcome one)  {nTrials}[TrialLength , no Units]
%  Rtrials_             Input Data for Right trials (i.e. outcome two) {nTrials}[TrialLength , no Units]
%  TrialsSubsample      How many trials to draw per subsample (Or Inf to use all trials)
%  InfoCriteria         Score Summary: 'max' or 'mean'
%  TimeRange            Optional: Either [Min Sample Max Sample] or 'Inf' to run on whole trial duration 
%
%  ***Outputs:
% CVE                   time-dependent mean decoding score [No. Timepoints x 1]
% Score                 time-dependent mean decoding score [No. Timepoints x 1]
% Score_Shuffled        Median Time-averaged decoding score from shuffled trials, [1 x 1]
% Score_ShuffledCIL     5 percentile Time-averaged decoding score from shuffled trials, [1 x 1]
% Score_ShuffledCIH     95 percentile Time-averaged decoding score from shuffled trials, [1 x 1]
% CVE_shuffledHL        mean time-dependent decoding score [No. Timepoints x 2]

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
    shuffleYN=true;
    Score_Shuffled    = nan(nReps,1);
    Score_ShuffledCIL = nan(nReps,1);
    Score_ShuffledCIH = nan(nReps,1);
    CVE_shuffledHL=[];
    CVE_shuffledH = nan(LTR,nReps);
    CVE_shuffledL = nan(LTR,nReps);
else
    shuffleYN=false;
end

% Draw sub-groups of trials to maintain DoF between experiments

%  CVE_shuffledL_ = nan(LTR,nReps);
%     CVE_shuffledH  = nan(LTR,nReps);

nBS = 1000;

    
for nrep = 1:nReps % parfor
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
    
    FR_   = [cell2mat(Ltrials_(rsL)');cell2mat(Rtrials_(rsR)')];
    
    % (1) Decode real
    if exist('DecodeCVE_mex','file')
        CVE_ = 1-DecodeCVE_mex(FR_(:,unitIds),evt0,0.05);
    else
        CVE_ = 1-DecodeCVE(FR_(:,unitIds),evt0,0.05);        
    end
    CVE(:,nrep) = CVE_;
    switch InfoCriteria
        case 'mean'
            Score(nrep) = mean(CVE_(timeRange),2);
        case 'max'
            Score(nrep) = max(CVE_(timeRange),[],2);
    end
    % (2) Decode shuffled versions - aim for 50 total random draws
    if shuffleYN
        parfor iDraw = 1:nBS
%             [nrep, iDraw]
            if exist('DecodeCVE_mex','file')
                CVEshuf_(iDraw,:) = 1-DecodeCVE_mex(FR_(:,unitIds),evt0(randperm(length(evt0))),0.05);
            else
                CVEshuf_(iDraw,:) = 1-DecodeCVE(FR_(:,unitIds),evt0(randperm(length(evt0))),0.05);
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
if shuffleYN
    Score_Shuffled    =  mean(Score_Shuffled);
    Score_ShuffledCIL =  mean(Score_ShuffledCIL);
    Score_ShuffledCIH =  mean(Score_ShuffledCIH);
    CVE_shuffledHL    =  [nanmean(CVE_shuffledL,2), nanmean(CVE_shuffledH,2)];
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

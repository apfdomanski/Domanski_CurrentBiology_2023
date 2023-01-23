function trialCounts = parseNLX_TTLsForRatios(inputData)
times = inputData (:,1);
TTLs =  inputData (:,2);

hasDelay= ismember(4128,TTLs);

noTrials  = numel(find(ismember(TTLs,[9,18])));
noLTrials = numel(find(TTLs==9));
noRTrials = numel(find(TTLs==18));

TrialStartIdx = find(ismember(TTLs,[9,18]));
TrialArray = cell(1,length(TrialStartIdx)); % [Timestamps,TTLs] associated with each trial
TrialTypeArray = zeros(noTrials,1);         % Left Sample = 1, Right Sample = 2
TrialDelayArray = zeros(noTrials,1);        % Delay for each trial
TrialOutcomeArray = zeros(noTrials,1);      % Outcome: 0 = Error, 1 = Correct
AbortedTrialsIdx =[];                       % Indices of trials lacking aborted before sample press
%%
for iTrial =1:length(TrialStartIdx)
    try
        % find events belonging to each trial epoch
        if iTrial<length(TrialStartIdx)
            idx = TrialStartIdx(iTrial):TrialStartIdx(iTrial+1)-1;        
            TrialArray{iTrial} = [times(idx),TTLs(idx)];
        else
            idx = TrialStartIdx(iTrial):length(TTLs);
            TrialArray{iTrial} = [times(idx),TTLs(idx)];
            i = find(TrialArray{iTrial}(:,2)==0,1,'first');
            TrialArray{iTrial} = TrialArray{iTrial}(1:i,:);
        end
        TTLs_ = TTLs(idx);
        times_= times(idx);
        % Left or Right?
        if sum(ismember(TTLs_,9))
            TrialTypeArray(iTrial) = 1;
        elseif sum(ismember(TTLs_,18))
            TrialTypeArray(iTrial) = 2;
        end
        
        if hasDelay
            %Delay length?
            tStart = times_(TTLs_==2080); % press sample
            tStop  = times_(TTLs_==4128); % end of delay
            TrialDelayArray(iTrial) = round((tStop  - tStart) /1e6);
        else
            tStart = times_(ismember(TTLs_,[6144,6176]));% press sample
            %tStart = times_(TTLs_==6176); % press sample
            tStop  = times_(TTLs_==8219); % end of delay
            TrialDelayArray(iTrial) = 0;
        end
        % Correct or Error?
        if sum(ismember(TTLs_,-32768))
            TrialOutcomeArray(iTrial) = 0;
        elseif sum(ismember(TTLs_,16420))
            TrialOutcomeArray(iTrial) = 1;
        elseif ~sum(ismember(TTLs_,16420)) | ~sum(ismember(TTLs_,-32768))
            AbortedTrialsIdx=[AbortedTrialsIdx;iTrial];
        end
        
    catch
        AbortedTrialsIdx=[AbortedTrialsIdx;iTrial];
    end
end
TrialDelayArray(AbortedTrialsIdx)=[];
TrialArray(AbortedTrialsIdx)=[];
TrialTypeArray(AbortedTrialsIdx)=[];
TrialOutcomeArray(AbortedTrialsIdx)=[];
RewardConsumed=[];
trialCounts = [sum(TrialOutcomeArray), sum(~TrialOutcomeArray),length(AbortedTrialsIdx)];
trialCounts = trialCounts./sum(trialCounts);
% trialCounts = trialCounts./sum(trialCounts(1:2));

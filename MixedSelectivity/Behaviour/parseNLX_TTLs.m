% no delay
%
% TTL       code
% 9         sample left
% 18        sample right
% 6144      press sample	or 6176?
% 32        maglight
% 8219      initiate poke
% 16420     press correct
% -32768	press incorrect
% 32 =       reward delivery signal 
% 0         reward consumption (Beambreak, triggers lights off)

% with delay
%
% 9         sample left
% 18        sample right
% 2079      press sample % Actually, I think this should be 2080
% 4128      end of delay
% 32        maglight
% 8219      initiate poke
% 16420     press correct
% -32768	press incorrect
% 32        reward delivery signal 
% 0         reward consumption (Beambreak, triggers lights off)

% inputData = [];
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
for iTrial =1:length(TrialArray)
    RewardConsumed(iTrial) = TrialArray{iTrial}(end,2)==0;
    if ~RewardConsumed(iTrial)
        TrialArray{iTrial} = [TrialArray{iTrial}; [NaN,0]];
    end
end
%%
delays = unique(TrialDelayArray);   
if numel(delays)==3
    delays_ = {'Short','Medium','Long'};
else
    delays_ = cellstr(num2str(delays));
    delays_ = cellfun(@(c)['Delay_' c],delays_,'uni',false);
end
TrialTypes = [9,18];                
TrialTypes_={'Left','Right'};
% Code_ = [6176,2080]; %sprintf will select the correct code based on hasDelay
Code_ = [6144,2080]; %sprintf will select the correct code based on hasDelay
for iDelay = 1:length(delays)
    for iTrialType=1:2
        
        % Correct Trials
        Trials_ = cell2mat({TrialArray{TrialDelayArray==delays(iDelay) & TrialOutcomeArray & TrialTypeArray==iTrialType}}');
        if ~isempty(Trials_)
            eval(sprintf('t.%s.CueLight_%sCorrect = Trials_(Trials_(:,2)==%d,1);',      delays_{iDelay},TrialTypes_{iTrialType},TrialTypes(iTrialType)));
            eval(sprintf('t.%s.SamplePress_%sCorrect = Trials_(Trials_(:,2)==%d,1);',   delays_{iDelay},TrialTypes_{iTrialType},Code_(hasDelay+1)));
            eval(sprintf('t.%s.DelayEnd_%sCorrect = Trials_(Trials_(:,2)==4128,1);',    delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.NosePoke_%sCorrect = Trials_(Trials_(:,2)==8219,1);',    delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.ChoicePress_%sCorrect = Trials_(Trials_(:,2)==16420,1);',delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.RewardConsume_%sCorrect = Trials_(Trials_(:,2)==0,1);',  delays_{iDelay},TrialTypes_{iTrialType}));
        else
            
            eval(sprintf('t.%s.CueLight_%sCorrect = [];',   delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.SamplePress_%sCorrect = [];',delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.DelayEnd_%sCorrect = [];',   delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.NosePoke_%sCorrect = [];',   delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.ChoicePress_%sCorrect = [];',delays_{iDelay},TrialTypes_{iTrialType}));
            
        end
        
        % This implentation leads to misalinged timestamps within trials
        
%         % Error Trials
%         Trials_ = cell2mat({TrialArray{TrialDelayArray == delays(iDelay) & ~TrialOutcomeArray & TrialTypeArray==iTrialType}}');
%         if ~isempty(Trials_)
%             eval(sprintf('t.%s.CueLight_%sError = Trials_(Trials_(:,2)==%d,1);',  delays_{iDelay},TrialTypes_{iTrialType},TrialTypes(iTrialType)));
%             eval(sprintf('t.%s.SamplePress_%sError = Trials_(Trials_(:,2)==%d,1);',  delays_{iDelay},TrialTypes_{iTrialType},Code_(hasDelay+1)));
%             eval(sprintf('t.%s.DelayEnd_%sError = Trials_(Trials_(:,2)==4128,1);',  delays_{iDelay},TrialTypes_{iTrialType}));
%             eval(sprintf('t.%s.NosePoke_%sError = Trials_(Trials_(:,2)==8219,1);',  delays_{iDelay},TrialTypes_{iTrialType}));
%             eval(sprintf('t.%s.ChoicePress_%sError = Trials_(Trials_(:,2)==-32768,1);',  delays_{iDelay},TrialTypes_{iTrialType}));
%         else
%             eval(sprintf('t.%s.CueLight_%sError = [];',   delays_{iDelay},TrialTypes_{iTrialType}));
%             eval(sprintf('t.%s.SamplePress_%sError = [];',delays_{iDelay},TrialTypes_{iTrialType}));
%             eval(sprintf('t.%s.DelayEnd_%sError = [];',   delays_{iDelay},TrialTypes_{iTrialType}));
%             eval(sprintf('t.%s.NosePoke_%sError = [];',   delays_{iDelay},TrialTypes_{iTrialType}));
%             eval(sprintf('t.%s.ChoicePress_%sError = [];',delays_{iDelay},TrialTypes_{iTrialType}));    
%         end
       
    end    
end
%%
% Error Trials

for iDelay = 1:length(delays)
    for iTrialType=1:2
        trialIdx = find(TrialDelayArray == delays(iDelay) & ~TrialOutcomeArray & TrialTypeArray==iTrialType);
        if ~isempty(trialIdx)
            
            for iTrial = 1:length(trialIdx)
                Trials_ = TrialArray{trialIdx(iTrial)};
                try
                    eval(sprintf('t.%s.CueLight_%sError(iTrial) = Trials_(Trials_(:,2)==%d,1);',     delays_{iDelay},TrialTypes_{iTrialType},TrialTypes(iTrialType)));
                catch
                    eval(sprintf('t.%s.CueLight_%sError(iTrial) = NaN;',     delays_{iDelay},TrialTypes_{iTrialType},TrialTypes(iTrialType)));
                end
                try
                    eval(sprintf('t.%s.SamplePress_%sError(iTrial) = Trials_(Trials_(:,2)==%d,1);',  delays_{iDelay},TrialTypes_{iTrialType},Code_(hasDelay+1)));
                catch
                    eval(sprintf('t.%s.SamplePress_%sError(iTrial) = NaN;',  delays_{iDelay},TrialTypes_{iTrialType}));
                end
                try
                    eval(sprintf('t.%s.DelayEnd_%sError(iTrial) = Trials_(Trials_(:,2)==4128,1);',   delays_{iDelay},TrialTypes_{iTrialType}));
                catch
                    eval(sprintf('t.%s.DelayEnd_%sError(iTrial) = NaN;',   delays_{iDelay},TrialTypes_{iTrialType}));
                end
                try
                    eval(sprintf('t.%s.NosePoke_%sError(iTrial) = Trials_(Trials_(:,2)==8219,1);',   delays_{iDelay},TrialTypes_{iTrialType}));
                catch
                    eval(sprintf('t.%s.NosePoke_%sError(iTrial) = NaN;',   delays_{iDelay},TrialTypes_{iTrialType}));
                end
                try
                    eval(sprintf('t.%s.ChoicePress_%sError(iTrial) = Trials_(Trials_(:,2)==-32768,1);',  delays_{iDelay},TrialTypes_{iTrialType}));
                catch
                    eval(sprintf('t.%s.ChoicePress_%sError(iTrial) = NaN;',  delays_{iDelay},TrialTypes_{iTrialType}));
                end
            end
        else
            eval(sprintf('t.%s.CueLight_%sError = [];',   delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.SamplePress_%sError = [];',delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.DelayEnd_%sError = [];',   delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.NosePoke_%sError = [];',   delays_{iDelay},TrialTypes_{iTrialType}));
            eval(sprintf('t.%s.ChoicePress_%sError = [];',delays_{iDelay},TrialTypes_{iTrialType}));
        end
    end
end
%%
clearvars -except t inputData
%%
% clear A
% B = t.Short.ChoicePress_LeftCorrect
% for iTrial=1:length(t.Short.RewardConsume_LeftCorrect)
%     A = t.Short.RewardConsume_LeftCorrect; A(iTrial)=[];
%     
%     A_(iTrial) = var(B-A);
% end
    
    
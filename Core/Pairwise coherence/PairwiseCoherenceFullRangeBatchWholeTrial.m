%%
% Compares continuous L/R discrimination during the delay period,
% calcluated on short/medium/long delay trials independently
%% %%%%%% PREAMBLE %%%%%%

clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
DelaysDur = [4,8,16];

Target = 'LONG';

tlimsAll = [-30 30];
tlimsTrials = [-5 5];
tlimsShort=[-5 6];
tlimsMedium=[-5 10];
tlimsLong=[-5 20];

plotOnline = false;
bw=0.05;
Nbs = 500;

minFR = 0.1;            % Minimal lifetime mean spike rate between pairs of cells
LimitBySpikes = true;   % Discard trials with low spike counts
minSpikesPerTrial = 30; % Assuming 10s trial duration Run3:10, Run5:30

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
    home = getenv('HOME');
else
    
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    addpath(genpath('/panfs/panasas01/phph/ad15419/MATLAB'))
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);

fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileListAss(find(name_flag))=[];% fnames_(name_flag)=[];
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

UnitSelection = 'pairs'; %{'pairs','groups'}
groupSize = 8;
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};

AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time
switch AssemblyChoice
    case 1
        pat2 = [pat 'KDE_bins' filesep Target filesep];
    case 2
        pat2 = [pat 'KDE_binsTaskonly' filesep 'LONGTaskonly' filesep];
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end
mkdir([pat 'MixedSelectivity'])

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
EventList = {'CueLight','SamplePress','DelayEnd','NosePoke','ChoicePress','RewardConsume'};
Direction = {'Left','Right'};
Outcome   = {'Correct','Error'};
%% Batch process units
for iFile =1:length(fileList)
    tic
    %% Get the files
    fname=strtok(fileList(iFile).name,'_');
    fprintf('Analysing run %d/%d %s...\n',iFile,length(fileList),fname)

    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname),'PFCcells','HPcells');
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
    spikes = [PFCcells;HPcells];
    
    evtTimes = ...
        [cell2mat(struct2cell(structfun(@(x) x(:),t.Short,'UniformOutput',false)));...
        cell2mat(struct2cell(structfun(@(x) x(:),t.Medium,'UniformOutput',false)));...
        cell2mat(struct2cell(structfun(@(x) x(:),t.Long,'UniformOutput',false)))];
    
    tlims_  = [min(evtTimes) max(evtTimes)] /1e6+tlimsAll;   
    %% Batch calculate coherence
    defaultProfile = parallel.defaultClusterProfile;
    myCluster = parcluster(defaultProfile);
    try
        parpool(myCluster);
    end
    
    % First run
%     params.pad       = -1;%0
%     params.Fs        = 1000;%100
%     params.tapers    = [4 9];% 4 9 % 10 19
%     params.movingwin = [5 1];% 10, 8
%     params.fpass     = [0.5 100];
%     params.err       = 0;
%     params.trialave  = 0;
%     params.fscorr    = 0;
%     % 2nd run (and 3rd,5th: 3rd run has min trial count for spikes =10,
%     5th=30
    params.pad       = -1;%0
    params.Fs        = 100;%100
    params.tapers    = [3 5];% 4 9 % 10 19
    params.movingwin = [5 1];% 10, 8
    params.fpass     = [0.5 100];
    params.err       = 0;
    params.trialave  = 0;
    params.fscorr    = 0;
    
 % 4th run: has min trial count for averaging
%     params.Fs        = 100;%100
%     params.tapers    = [3 5];% 4 9 % 10 19
%     params.movingwin = [0.1 0.05];% 10, 8
%     params.fpass     = [0.5 100];
%     params.err       = 0;
%     params.trialave  = 0;
%     params.fscorr    = 0;
    
    coherence_ = cell(length(spikes));
    coherence_mean = cell(length(spikes));
    time_      = cell(length(spikes));
    freq_      = cell(length(spikes));
    
    %     cache minimum firing rate of pair
    MeanRate = cell(length(spikes));
    MinRate=zeros(length(spikes));
    for i = 1:length(spikes)
        for j = i+1:length(spikes)
                data1 = extractdatapt(spikes{i}.t*1e-4, tlims_,0);
                data2 = extractdatapt(spikes{j}.t*1e-4, tlims_,0);
                MeanRate{i,j} = [length(data1.times)./range(data1.times);length(data2.times)./range(data2.times)];
                MinRate(i,j) = min(MeanRate{i,j});
        end
    end
            
    for i = 1:length(spikes)-1
        fprintf('Evaulating spike train coherence from unit seed [%d/%d])\n',i,length(spikes))
        data1 = extractdatapt(spikes{i}.t*1e-4, tlims_,0);
        pairsToRun = find(MinRate(i,:)>minFR);
        % cache spectrograms
        if ~isempty(pairsToRun)
        C_ = cell(length(pairsToRun),1);
        t_ = cell(length(pairsToRun),1);
        f_ = cell(length(pairsToRun),1);
        parfor j = 1:length(pairsToRun)
            fprintf('Evaulating spike train coherence between unit %d and %d (of %d units)\n', i,j,length(pairsToRun))
            data2 = extractdatapt(spikes{pairsToRun(j)}.t*1e-4, tlims_,0);
            [C_{j},~,~,~,~,t_{j},f_{j}]=cohgrampt(data1,data2,params.movingwin,params,params.fscorr);
        end
        t_ = t_{1};
        f_ = f_{1};
        % cutout times of interest
        dt = mean(diff(t_));
        nPoints = sum(abs(tlimsTrials))/dt;
        for j = 1:length(pairsToRun)
            data2 = extractdatapt(spikes{pairsToRun(j)}.t*1e-4, tlims_,0);
            for iDelay = 1:length(Delays_)
                for iOutcome = 1:2
                    if iOutcome==1
                        nEvents = length(EventList);
                    else
                        nEvents = length(EventList)-1;
                    end
                    
                    % Loop across events
                    for iEvent = 1:nEvents
                        Ccutout = cell(2,1);
                        CcutoutMean = cell(2,1);
                        for iDirection =1:2
                            evt_ = sprintf('t.%s.%s_%s%s',Delays_{iDelay},EventList{iEvent},Direction{iDirection},Outcome{iOutcome});
                            Times_ = eval([evt_ '(:);']);Times_(isnan(Times_))=[];
                            
                            Ccutout{iDirection,1} = [];
                            CcutoutMean{iDirection,1} = [];
                            for iTrial =1:size(Times_,1)
                                try
                                    tlim_  = Times_(iTrial,1)/1e6+tlimsTrials(1);
                                    
                                    tlim_  = closest(t_,tlim_);
                                    tlim_  = tlim_(1) : (tlim_(1)+nPoints-1);
                                    
                                    retainTrial= true;
                                    if LimitBySpikes 
                                        data1_ = extractdatapt(data1, t_([tlim_(1),tlim_(end)]),0);
                                        data2_ = extractdatapt(data2, t_([tlim_(1),tlim_(end)]),0);

                                        if min([length(data1_.times),length(data2_.times)])<=minSpikesPerTrial
                                           retainTrial = false;
                                        end
                                    end
                                    if retainTrial
                                        Ccutout{iDirection,1} = cat(3,Ccutout{iDirection,1},C_{j}(tlim_,:));  
                                    end
                                end
                            end
                            CcutoutMean{iDirection,1} = nanmean(Ccutout{iDirection,1},3); 
                        end
                        eval(sprintf('coherence_{i,pairsToRun(j)}.%s.%s.%s = Ccutout;',Delays_{iDelay},EventList{iEvent},Outcome{iOutcome}))
                        eval(sprintf('coherence_mean{i,pairsToRun(j)}.%s.%s.%s = CcutoutMean;',Delays_{iDelay},EventList{iEvent},Outcome{iOutcome}))
                        
                        
                    end
                    
                    % Special case to calculate delay period
                    Ccutout = cell(2,1);
                    CcutoutMean{iDirection,1} = [];
                    for iDirection = 1:2
                        evt_ = sprintf('t.%s.%s_%s%s',Delays_{iDelay},'SamplePress',Direction{iDirection},Outcome{iOutcome});
                        Times_ = eval([evt_ '(:);']);Times_(isnan(Times_))=[];
                        Ccutout{iDirection,1} = [];
                        CcutoutMean{iDirection,1} = [];
                        for iTrial =1:size(Times_,1)
                            try
                                tlim_  = Times_(iTrial,1)/1e6+tlimsTrials(1);
                                tlim_  = closest(t_,tlim_);
                                tlim_  = tlim_(1) : (tlim_(1) + DelaysDur(iDelay)/dt + nPoints-1);
                                
                                retainTrial= true;
                                if LimitBySpikes
                                    data1_ = extractdatapt(data1, t_([tlim_(1),tlim_(end)]),0);
                                    data2_ = extractdatapt(data2, t_([tlim_(1),tlim_(end)]),0);
                                    
                                    if min([length(data1_.times),length(data2_.times)])<=minSpikesPerTrial
                                        retainTrial = false;
                                    end
                                end
                                if retainTrial
                                    Ccutout{iDirection,1} = cat(3,Ccutout{iDirection,1},C_{j}(tlim_,:));
                                end
                                
                                
                            end
                        end
                        CcutoutMean{iDirection,1} = nanmean(Ccutout{iDirection,1},3); 
                    end
                    eval(sprintf('coherence_{i,pairsToRun(j)}.%s.%s.%s = Ccutout;',Delays_{iDelay},'DelayPeriod',Outcome{iOutcome}))
                    eval(sprintf('coherence_mean{i,pairsToRun(j)}.%s.%s.%s = CcutoutMean;',Delays_{iDelay},'DelayPeriod',Outcome{iOutcome}))
                end
                
            end
            coherence_{i,pairsToRun(j)}.f = f_;
            coherence_{i,pairsToRun(j)}.MeanRate = MeanRate{i,pairsToRun(j)};
        end
        end
    end
    %% Save
%     save(fullfile(pat,'SpikeTrainCoherence', [fname,'_coherence.mat']),'coherence_','params','t_','f_','MeanRate','-v7.3')
    save(fullfile(pat,'SpikeTrainCoherence', [fname,'_coherenceMean.mat']),'coherence_mean','params','t_','f_','MeanRate','-v7.3')
    disp([sprintf('Finished  file %d of %d in %.2g hours',iFile, length(fileList), toc/3600), ' (' fname ')'])
end

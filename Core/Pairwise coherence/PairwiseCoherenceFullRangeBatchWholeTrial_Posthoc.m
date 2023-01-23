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

minFR = 0.5;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

clear Av

warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
%     pat = '/Volumes/HDD2/DNMTP/raw/';
    pat = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw';
    home = getenv('HOME');
else
    
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    addpath(genpath('/panfs/panasas01/phph/ad15419/MATLAB'))
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'SpikeTrainCoherence' filesep '*' Target '*coherence.mat']);

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
for iFile = 1:length(fileList)
    tic
    %% Get the files
    fname=strtok(fileList(iFile).name,'_');
    fprintf('Analysing run %d/%d %s...\n',iFile,length(fileList),fname)

    load(fullfile(pat,'SpikeTrainCoherence',sprintf('%s_coherence.mat',fname)));
   
    
    %% Batch calculate coherence
   
    coherence_mean = cell(size(coherence_));
    
    for i = 1:size(coherence_,1)
        for j = 1:size(coherence_,1)
            
            if ~isempty(coherence_{i,j})
                for iDelay = 1:length(Delays_)
                    for iOutcome = 1:2
                        if iOutcome==1
                            nEvents = length(EventList);
                        else
                            nEvents = length(EventList)-1;
                        end
                        
                        % Loop across events
                        for iEvent = 1:nEvents
                            CcutoutMean = cell(2,1);
                            evt_ = sprintf('%s.%s.%s',Delays_{iDelay},EventList{iEvent},Outcome{iOutcome});
                            eval(sprintf('Ccutout = coherence_{i,j}.%s;',evt_));
                            CcutoutMean = cellfun(@(x) nanmean(x,3),Ccutout,'UniformOutput',false);
                            eval(sprintf('coherence_mean{i,j}.%s = CcutoutMean;',evt_))
                        end
                        
                        % Special case to calculate delay period
                        evt_ = sprintf('%s.%s.%s',Delays_{iDelay},'DelayPeriod',Outcome{iOutcome});
                        eval(sprintf('Ccutout = coherence_{i,j}.%s;',evt_));
                        CcutoutMean = cellfun(@(x) nanmean(x,3),Ccutout,'UniformOutput',false);
                        eval(sprintf('coherence_mean{i,j}.%s = CcutoutMean;',evt_))
                    end
                end
            end
        end
    end
    %% Save
    save(fullfile(pat,'SpikeTrainCoherence', [fname,'_coherenceMean.mat']),'coherence_mean','params','t_','f_','MeanRate','-v7.3')
    disp([sprintf('Finished  file %d of %d in %.2g hours',iFile, length(fileList), toc/3600), ' (' fname ')'])
end



















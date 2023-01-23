%% Preamble, get file list
clear all; close all
set(0,'DefaultFigureWindowStyle','docked')
flags.plotIndividualChanges = false;
% pat{1} -  % location of the processed assembly memberships
% pat{2} -  % location of raw spike times
if ispc
    %     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - epoch specific FRs\';
    pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - use task mean FRs\';
    pat{2} = 'C:\Analysis\AssemblyAnalysis\raw';          
elseif ismac
    pat{1} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Decoding vs sleep/';
    pat{2} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Task';
else
%     path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    addpath(genpath('/panfs/panasas01/phph/ad15419/MATLAB/SpikeTrainAna'))
    p.home = getenv('HOME');
    cd([p.home '/MATLAB/SpikeTrainAna/'])
    pat{1} = [p.home '/Sleep/DecodingVsSleep/'];
    pat{2} = [p.home '/raw/'];
end
cd(pat{1})
RatList = dir([pat{1} filesep '*_DecodingVSsleep.mat']);

% Uncomment to ignore some specific files...

% Ignore= {'JaroslawLONG1'; ...
%          'MiroslawLONG2'; ...
%          'NorbertLONG2' ; ...
%          'OnufryLONG2'};
% for iList = 1:length(Ignore)
%     RatList = RatList(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
%            {RatList.name},repmat({Ignore(iList)},1,length(RatList)),'UniformOutput',false)));
% end

%% Runtime
defaultProfile = parallel.defaultClusterProfile;
myCluster = parcluster(defaultProfile);
parpool(myCluster);

params.pad       = -1;%0
params.Fs        = 1000;%100
params.tapers    = [4 9];% 10 19
params.movingwin = [5 1];% 10, 8
params.fpass     = [0.5 200];
params.err       = 0;
params.trialave  = 0;
params.fscorr    = 0;


for iRat =1:length(RatList)
    tic
    disp([sprintf('Working on rat %d of %d',iRat, length(RatList)), ' (' RatList(iRat).name ')'])
    load([pat{1} RatList(iRat).name],'P','D','FRtrials','FAtrials');
    load([pat{2} filesep strtok(RatList(iRat).name,'_')],'HPcells','PFCcells')
    
    spikes = [PFCcells(FAtrials.unitIDs{1});...
              HPcells(FAtrials.unitIDs{2})];
    
	coherence = cell(length(spikes),...
                     length(spikes));
    t     = cell(length(spikes),...
                     length(spikes));

	%%% Prepare a list of cell pairs of interest of which to restrict
	%%% analysis to
    index = cell(length(spikes),...
                length(spikes));
            
    % exclusively analyse HP-PFC pairs
    for i = 1:length(HPcells) %HP
        for j = length(PFCcells)+1:length(spikes) %PFC
            index{i,j} = [i,j];
        end
    end
    runList = find(~cellfun(@isempty,index));
    outputList = cell(length(runList),1);   
    outputList_t = cell(length(runList),1);   
    parfor i = 1:length(runList)
            C_ = [];
            fprintf('Evaulating spike train coherence between unit %d and %d (of %d unit pairs)\n', index{runList(i)},length(runList))
            
            try
                
                data1 = extractdatapt(spikes{index{runList(i)}(1)}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
                data2 = extractdatapt(spikes{index{runList(i)}(2)}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
                
                [C,~,~,~,~,t_,f_]=cohgrampt(data1,data2,params.movingwin,params,params.fscorr);
                % fprintf('Time resolution = %.1d, frequency resolution = %.1d',min(diff(t)), unique(diff(f)))
                
                % Fudge the frequency resolution to 1Hz
                for band =1:params.fpass(2)
                    C_(:,band)= nanmean(C(:,find(f_>(band-1) & f_<=band)),2);
                end
                %coherence{i,j} = nanmean(C(:,find(f>=3.5 & f<=5.5)),2);
%                 coherence{index{runList(i)}(1),index{runList(i)}(2)} = C_;
                outputList{i} = C_
                outputList_t{i} = t_;
            catch
                fprintf('Error on coherence calculation between unit %d and %d\n', index{runList(i)})
            end
    end
    
    % reshape outputs
    for i = 1:length(runList)
        coherence{runList(i)} = outputList{i};
        t    {runList(i)} = outputList_t{i};
        outputList{i} = [];
        outputList_t{i} = [];
    end
    f = 1:params.fpass(2);
    save([pat{1} strtok(RatList(iRat).name,'_'),'_coherenceFull.mat'],'coherence','params','t','f','-v7.3')
    disp([sprintf('Finished  rat %d of %d in %.2g hours',iRat, length(RatList), toc/3600), ' (' RatList(iRat).name ')'])
end


% delete(myCluster);
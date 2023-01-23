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
    %    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
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

params.band = [0 1.5;...    % alpha
               3   5;...    % 4Hz
               7   12;...   % theta
               15  25;...   % beta
               30  50;...   % low gamma
               51  80;...   % high gamma
               81  100;...  % very high gamma
               120 180];    % spindle band    
for iRat =1:length(RatList)
    tic

    disp([sprintf('Working on rat %d of %d',iRat, length(RatList)), ' (' RatList(iRat).name ')'])
    load([pat{1} RatList(iRat).name],'P','D','FRtrials','FAtrials');
    load([pat{2} filesep strtok(RatList(iRat).name,'_')],'HPcells','PFCcells')
    
    spikes = [PFCcells(FAtrials.unitIDs{1});...
              HPcells(FAtrials.unitIDs{2})];
    
	coherence = cell(length(spikes),length(spikes));
    t_ = cell(length(spikes),length(spikes));
    %f_  = cell(length(spikes),length(spikes),params.fpass(2),1);
    for i = 1:length(spikes)
        parfor j = 1:length(spikes)
            
            fprintf('Evaulating spike train coherence between unit %d and %d (of %d units)\n', i, j,length(spikes))
            
            try
                
                data1 = extractdatapt(spikes{i}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
                data2 = extractdatapt(spikes{j}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
                
                [C,~,~,~,~,t,f] = cohgrampt(data1,data2,params.movingwin,params,params.fscorr);
                % fprintf('Time resolution = %.1d, frequency resolution = %.1d',min(diff(t)), unique(diff(f)))

                for band = 1:size(params.band,1)
                    coherence{i,j}(:,band) = nanmean(C(:, f > params.band(band,1) & f <= params.band(band,2)),2);
                end
                %coherence{i,j} = nanmean(C(:,find(f>=3.5 & f<=5.5)),2); % 4Hz only
                t_{i,j} = t;
                %f_{i,j} = f;
                
            catch
                fprintf('Error on coherence calculation between unit %d and %d\n', i, j)
            end
            
        end
    end
    save([pat{1} strtok(RatList(iRat).name,'_'),'_coherence.mat'],'coherence','params','t_','-v7.3')%,'f_'
    disp([sprintf('Finished  rat %d of %d in %.2g hours',iRat, length(RatList), toc/3600), ' (' RatList(iRat).name ')'])
end


%delete(myCluster);
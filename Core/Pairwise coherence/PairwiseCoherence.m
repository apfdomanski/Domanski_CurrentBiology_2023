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
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    pat{1} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Decoding vs sleep/';
    pat{2} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Task';
end
pat{3} = [pat{2} filesep 'KDE_bins'];                  % location of calculated firing rates for Task period
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
%% Demo version

% Group = struct;
% params.pad       = 0;
% params.Fs        = 1;
% params.tapers    = [3 5];
% params.movingwin = [5 1];
% params.fpass     = [0.5 200];
% params.err       = 0;
% params.trialave  = 0;
% params.fscorr    = 0;
% 
% for iRat =1:length(RatList)
%     disp([sprintf('Working on rat %d of %d',iRat, length(RatList)), ' (' RatList(iRat).name ')'])
%     load([pat{1} RatList(iRat).name],'P','D','FRtrials','FAtrials');
%     load([pat{2} filesep strtok(RatList(iRat).name,'_')],'HPcells','PFCcells')
%     
%     spikes{1} = PFCcells(FAtrials.unitIDs{1});
%     spikes{2} = HPcells(FAtrials.unitIDs{2});
%     spikes{3} = [spikes{1};spikes{2}];
%     
%     data1 = extractdatapt(spikes{1}{1}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
%     data2 = extractdatapt(spikes{1}{2}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
%     [C,phi,S12,S1,S2,t,f]=cohgrampt(data1,data2,params.movingwin,params,params.fscorr);
%     
% end



params.pad       = -1;
params.Fs        = 1000;
params.tapers    = [4 9];% 10 19
params.movingwin = [5 1];% 10, 8
params.fpass     = [0.5 100];
params.err       = 0;
params.trialave  = 0;
params.fscorr    = 0;
tic

% for i=1:length(spikes)
%     a(i)=length(spikes{i}.t)
% end
% bins  = 1:500:15000
% figure, plot(bins,histc(a,bins))

i = 39, j = 44
data1 = extractdatapt(spikes{1}{i}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
data2 = extractdatapt(spikes{1}{j}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
fprintf('%d and %d spikes\n',length(data1.times),length(data2.times))
[C,~,~,~,~,t,f]=cohgrampt(data1,data2,params.movingwin,params,params.fscorr);
fprintf('Took %.0d seconds. Time resolution = %.1d, frequency resolution = %.1d\n',...
        round(toc),min(diff(t)), min(unique(diff(f))))
clf, imagesc(t-min(t),f,smooth2a(C,50,0)); set(gca,'YDir','normal','YScale','linear'),hold on
% plot(t-min(t),40+100*nanmean(C(:,find(f>=50 & f<=60)),2),'r')
plot(t-min(t),100*nanmean(C(:,find(f>=3.5 & f<=5.5)),2),'k')

xlim([0 600])
caxis([0 1])
%% Runtime
defaultProfile = parallel.defaultClusterProfile;
myCluster = parcluster(defaultProfile);
parpool(myCluster);



params.pad       = -1;%0
params.Fs        = 1000;%100
params.tapers    = [4 9];% 10 19
params.movingwin = [5 1];% 10, 8
params.fpass     = [0.5 100];
params.err       = 0;
params.trialave  = 0;
params.fscorr    = 0;
tic

[C,~,~,~,~,t,f]=cohgrampt(data1,data2,params.movingwin,params,params.fscorr);
fprintf('Time resolution = %.1d, frequency resolution = %.1d\n',min(diff(t)), unique(diff(f)))
clf, imagesc(t-min(t),f,C); set(gca,'YDir','normal','YScale','linear')    
hold on, plot(t,100*nanmean(C(:,find(f>=3.5 & f<=5.5)),2),'k')
xlim([min(t) min(t)+600])
caxis([0 1])

for iRat =1%:length(RatList)
    disp([sprintf('Working on rat %d of %d',iRat, length(RatList)), ' (' RatList(iRat).name ')'])
    load([pat{1} RatList(iRat).name],'P','D','FRtrials','FAtrials');
    load([pat{2} filesep strtok(RatList(iRat).name,'_')],'HPcells','PFCcells')
    
    spikes = [PFCcells(FAtrials.unitIDs{1});...
              HPcells(FAtrials.unitIDs{2})];

    coherence =cell(length(spikes),length(spikes));
    for i = 1:length(spikes)
        parfor j = 1:length(spikes)
            fprintf('Evaulating spike train coherence between unit %d and %d (of %d units)\n', i, j,length(spikes))
            try
                data1 = extractdatapt(spikes{i}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
                data2 = extractdatapt(spikes{j}.t*1e-4, [min(FRtrials.Tmtx), max(FRtrials.Tmtx)],0);
                [C,~,~,~,~,t,f]=cohgrampt(data1,data2,params.movingwin,params,params.fscorr);
                % fprintf('Time resolution = %.1d, frequency resolution = %.1d',min(diff(t)), unique(diff(f)))
                coherence{i,j} = nanmean(C(:,find(f>=3.5 & f<=5.5)),2);
                %imagesc(t/3600,f,C); set(gca,'YDir','normal')
            catch
                fprintf('Error on coherence calculation between unit %d and %d\n', i, j)
            end
            
            
        end
    end
end
toc
coher.pairwiseUnits{iRat}= coherence;
save('JarL1Coherence.mat','coherence')
%% Post process
% len_= length(t);    
len_= max(unique(cellfun(@length,coherence(1:end))));
temp = cell(size(coherence));
SameAssCoherence = cell(length(FAtrials.units{3}),1);     % Same assembly
DiffAssCoherence = cell(length(FAtrials.units{3}),1);     % Different assembly
AllAssMembers    = unique(cell2mat(FAtrials.units{3}));   % All member assembly units
noUnits          = cellfun(@length,FAtrials.unitIDs);
for Ass_idx = 1:length(FAtrials.units{3})                 % loop across assemblies
    
    SameAssCoherence{Ass_idx}=[];
    AssMembers      = FAtrials.units{3}{Ass_idx};         % Member units for this assembly
    OtherAssMembers = setdiff(AllAssMembers,AssMembers);
   
    for i = 1:noUnits(3)
        for j = 1:noUnits(3)
            
            % [i j]
            
            %%%% Conditional sorting of unit pairs  %%%%
            
            % (1) Members of same assembly
            if      i~=j && ...                                                    % ...Not the same unit
                    ismember(i,AssMembers) && ismember(j,AssMembers) && ...        % ...Both members of the same assembly
                    i<=noUnits(1) && j > noUnits(1)                                % ...Live in different areas
           
                    temp{i,j} = [i,j];
                    SameAssCoherence{Ass_idx} = [SameAssCoherence{Ass_idx}, [nan(len_-length(coherence{i,j}),1);coherence{i,j}]];
            
            % (2) Members of different assembly
            elseif  i~=j && ...                                                    % ...Not the same unit
                    ismember(i,AssMembers) && ismember(j,OtherAssMembers) && ...   % ...Members of different assemblies
                    i<=noUnits(1) && j > noUnits(1)                                % ...Live in different areas
           
                    DiffAssCoherence{Ass_idx} = [DiffAssCoherence{Ass_idx}, [nan(len_-length(coherence{i,j}),1);coherence{i,j}]];
                
            else
            end
        end
    end
end
idx = cellfun(@isempty,temp);
idx(:,1:noUnits(1)) = 0; 
idx(noUnits(2):noUnits(3),:)=0;
temp = coherence(idx);    
NonAssCoherence = [];
for i=1:length(temp)
    NonAssCoherence= [NonAssCoherence, [nan(len_-length(temp{i}),1);temp{i}]];
end
    NonAssCoherence(:,nanvar(NonAssCoherence)<0.02)=[];

%%
bins =0:0.01:1;
AssCoherenceBinned=[];
Ass2CoherenceBinned=[];
NonAssCoherenceBinned=[];
temp = cell2mat(SameAssCoherence');
temp2 = cell2mat(DiffAssCoherence');

for iAss = 1:size(temp,2)
     AssCoherenceBinned(:,iAss) = histc(temp(:,iAss),bins);
end
for iAss = 1:size(temp2,2)
     Ass2CoherenceBinned(:,iAss) = histc(temp2(:,iAss),bins);
end

for iAss = 1:size(NonAssCoherence,2)    
    NonAssCoherenceBinned(:,iAss) = histc(NonAssCoherence(:,iAss),bins);
end

length(bins)
nBS = 1000;
NonAssCoherenceBinnedBS=[];
AssCoherenceBinnedBS=[];
OtherAssCoherenceBinnedBS=[];

k = round(0.5*(min([size(NonAssCoherenceBinned,2),size(Ass2CoherenceBinned,2),size(AssCoherenceBinned,2)])));
for iBS = 1:nBS
   fprintf('Bootsrap draw no. %d of %d\n',iBS,nBS)
   NonAssCoherenceBinnedBS  (:,iBS) = nanmean(NonAssCoherenceBinned(:,randperm(size(NonAssCoherenceBinned,2),k)),2);
   AssCoherenceBinnedBS     (:,iBS) = nanmean(AssCoherenceBinned(:,randperm(size(AssCoherenceBinned,2),k)),2);
   OtherAssCoherenceBinnedBS(:,iBS) = nanmean(Ass2CoherenceBinned(:,randperm(size(Ass2CoherenceBinned,2),k)),2);
end
AssCoherenceBinnedBS_mean = nanmean(AssCoherenceBinnedBS');
NonAssCoherenceBinnedBS_mean= nanmean(NonAssCoherenceBinnedBS');
OtherAssCoherenceBinnedBS_mean= nanmean(OtherAssCoherenceBinnedBS');

AssCoherenceBinnedBS_CI = nanstd(AssCoherenceBinnedBS');
NonAssCoherenceBinnedBS_CI = nanstd(NonAssCoherenceBinnedBS');
OtherAssCoherenceBinnedBS_CI = nanstd(OtherAssCoherenceBinnedBS');

AssCoherenceBinnedBS_mean = AssCoherenceBinnedBS_mean./length(bins);
OtherAssCoherenceBinnedBS_mean = OtherAssCoherenceBinnedBS_mean./length(bins);
NonAssCoherenceBinnedBS_mean = NonAssCoherenceBinnedBS_mean./length(bins);

AssCoherenceBinnedBS_CI    = AssCoherenceBinnedBS_CI./length(bins);
NonAssCoherenceBinnedBS_CI = NonAssCoherenceBinnedBS_CI./length(bins);
OtherAssCoherenceBinnedBS_CI = OtherAssCoherenceBinnedBS_CI./length(bins);

figure; hold on
plot(bins,AssCoherenceBinnedBS_mean,'b','LineWidth',2)
plot(bins,OtherAssCoherenceBinnedBS_mean,'g','LineWidth',2)
plot(bins,NonAssCoherenceBinnedBS_mean,'k','LineWidth',2)
legend('Same CA1-mPFC assembly member pairs ',...
       'Different CA1-mPFC assembly member pairs ',...
       'Non-member CA1-mPFC pairs'); legend boxoff
ciplot(AssCoherenceBinnedBS_mean+AssCoherenceBinnedBS_CI,...
       AssCoherenceBinnedBS_mean-AssCoherenceBinnedBS_CI,...
       bins,'b')
ciplot(OtherAssCoherenceBinnedBS_mean+OtherAssCoherenceBinnedBS_CI,...
       OtherAssCoherenceBinnedBS_mean-OtherAssCoherenceBinnedBS_CI,...
       bins,'g')   
ciplot(NonAssCoherenceBinnedBS_mean+NonAssCoherenceBinnedBS_CI,...
       NonAssCoherenceBinnedBS_mean-NonAssCoherenceBinnedBS_CI,...
       bins,'k')
axis([0 1 0 Inf])
view(90, -90)
xlabel('4-5Hz Spike train Coherence')
ylabel('Distribution of experiment time')

%% Preamble, get file list
clear all; close all
set(0,'DefaultFigureWindowStyle','docked')

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
RatList = dir([pat{1} filesep '*_Decoding*.mat']);

% Uncomment to ignore some specific files...

% Ignore= {'JaroslawLONG1'; ...
%          'MiroslawLONG2'; ...
%          'NorbertLONG2' ; ...
%          'OnufryLONG2'};

Ignore ={};
% Ignore ={'MiroslawLONG0';...
%          'KrzysztofLONG2';...
%          'KrzysztofLONG1';...
%          'KrzesimirLONG2';...
%          'KrzesimirLONG1';...
%          'JaroslawLONG2';...
%          'JaroslawLONG1'};

for iList = 1:length(Ignore)
    RatList = RatList(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
        {RatList.name},repmat({Ignore(iList)},1,length(RatList)),'UniformOutput',false)));
end
clear iList
%%
iRat = 2


disp([sprintf('Working on rat %d of %d',iRat, length(RatList)), ' (' RatList(iRat).name ')'])

load([pat{1} strtok(RatList(iRat).name,'_') '_Task_DecodingVSsleep.mat'],'P','D','FRtrials','FAtrials','FAcont','FRcont');
load([pat{2} filesep strtok(RatList(iRat).name,'_')],'HPcells','PFCcells')

%% Plot


 spikes = {PFCcells(FAtrials.unitIDs{1}),HPcells(FAtrials.unitIDs{2}),...
          [PFCcells(FAtrials.unitIDs{1});HPcells(FAtrials.unitIDs{2})]};
          
          
% AssemblyRasterPlot(FAcont.Task,FRcont.Task,spikes,FRtrials,[4760 5000],3)
% rat 4



AssemblyRasterPlot(FAcont.Task,FRcont.Task,spikes,FRtrials,[5600 5820],3) 
% AssemblyRasterPlotAll(FAcont.Task,FRcont.Task,spikes,FRtrials,[5600 5820],3) 

%%
figure ; hold on 
col_ = [0.635,0.078,0.184];% %[0,0.447,0.741]
temp = [FAcont.Pre.FSC{3}(:,1)-mean(FAcont.Pre.FSC{3}(:,1));...
        FAcont.Task.FSC{3}(:,1)-mean(FAcont.Task.FSC{3}(:,1));...
        FAcont.Post.FSC{3}(:,1)-mean(FAcont.Post.FSC{3}(:,1))];

tb = [FRcont.Pre{1, 1}.Tmtx(1:end-1)';FAcont.Task.TmtxS';FRcont.Post{1, 1}.Tmtx(1:end-1)'];
% plot(tb,temp,'Color',col_)

plot(FRcont.Pre{1, 1}.Tmtx(1:end-1)/3600,FAcont.Pre.FSC{3}(:,1)-mean(FAcont.Pre.FSC{3}(:,1)),'Color',col_)
plot(FAcont.Task.TmtxS/3600, FAcont.Task.FSC{3}(:,1)-mean(FAcont.Task.FSC{3}(:,1)),'Color',col_)
plot(FRcont.Post{1, 1}.Tmtx(1:end-1)/3600,FAcont.Post.FSC{3}(:,1)-mean(FAcont.Post.FSC{3}(:,1)),'Color',col_)

plot([3 3.25],[6 6],'LineWidth',1.5,'Color','k')
plot([3 3],[6 10],'LineWidth',1.5,'Color','k')
text(3.125,7,'15mins','HorizontalAlignment','center')
text(2.95,8,{'Factor';'Score';'(A.U.)'},'HorizontalAlignment','right')
axis off

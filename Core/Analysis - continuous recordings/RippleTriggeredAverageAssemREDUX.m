% Group analysis for single unit / assembly decoding capability, activation
% patterns for assemblies
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
%% Preamble, get file list
clear all; close all
set(0,'DefaultFigureWindowStyle','docked')
flags.plotIndividualChanges = false;
if ispc
%     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\';       % location of the processed sleep assembly activations
%     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - epoch specific FRs\';
    pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - use task mean FRs\';
    pat{2} = 'C:\Analysis\AssemblyAnalysis\raw';           % location of raw spike times
else
    pat{1} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Decoding vs sleep/';
    pat{2} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Task';
end
pat{3} = [pat{2} filesep 'KDE_bins'];                  % location of calculated firing P.rates for Task period
cd(pat{1})
Rats = dir([pat{1} filesep '*_DecodingVSsleep.mat']);

% Ignore= {'JaroslawLONG1'; ...
%          'MiroslawLONG2'; ...
%          'NorbertLONG2' ; ...
% %          'OnufryLONG2'};
% for iList = 1:length(Ignore)
%     Rats = Rats(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
%            {Rats.name},repmat({Ignore(iList)},1,length(Rats)),'UniformOutput',false)));
% end
iRat = 2

%% Analyse interplay between SWRs and assembly activation times
    load([pat{1} Rats(iRat).name]);

% Run Project_Sleep_into_Task first (... including DecodeAssemFromSleepProjections)
temp = load('C:\Analysis\AssemblyAnalysis\Sleep\Ripples\JaroslawLONG2Ripp.mat');
temp_ = reshape([temp.Jlong2rip5sd_R1.tsabsmax],size(temp.Jlong2rip5sd_R1)); 
rip.PreTimes  = reshape([temp.Jlong2rip5sd_R1.tsabsmax],size(temp.Jlong2rip5sd_R1)); 
rip.TaskTimes =  reshape([temp.Jlong2rip5sd_R1.tsabsmax],size(temp.Jlong2rip5sd_R1)); 
rip.PostTimes = reshape([temp.Jlong2rip5sd_R1.tsabsmax],size(temp.Jlong2rip5sd_R1)); 

rip.PreTimes(rip.PreTimes>min(FAcont.Task.TmtxS)) = [];
rip.TaskTimes (rip.TaskTimes<min(FAcont.Task.TmtxS) | rip.TaskTimes>max(FAcont.Task.TmtxS))=[];
rip.PostTimes (rip.PostTimes<max(FAcont.Task.TmtxS)) = [];

%% Plot Task
tb =0:0.05:20;
for s= 1:3
    for iAss = 1:size(FAcont.Task.FSC{s},2)
        for iRip = 1:length(rip.TaskTimes)
            time_ = find(abs(FAcont.Task.TmtxS-rip.TaskTimes(iRip))==min(abs(FAcont.Task.TmtxS-rip.TaskTimes(iRip))));
            try 
            
                rip_FSC{s}{iAss}(:,iRip) = FAcont.Task.FSC{s}([time_-200:time_+200],iAss);
                
            end
        end
        %rip_FSC_mean{s}(:,iAss)= nanmean(rip_FSC{s}{iAss}-nanmean(rip_FSC{s}{iAss}(1:20,:)),2);
        %rip_FSC_SEM{s}(:,iAss)= nansem(rip_FSC{s}{iAss}-nanmean(rip_FSC{s}{iAss}(1:20,:)),2);
        rip_FSC_mean{s}(:,iAss)= nanmean(rip_FSC{s}{iAss},2);
        rip_FSC_SEM{s}(:,iAss)= nansem(rip_FSC{s}{iAss},2);
    end
end
                
figure;
for s= 1:3
                    subplot(1,3,s); hold on

    for iAss = 1:size(FAcont.Task.FSC{s},2)
                
                plot(tb-10, rip_FSC_mean{s}(:,iAss),'k')
                ciplot(rip_FSC_mean{s}(:,iAss)+rip_FSC_SEM{s}(:,iAss),...
                       rip_FSC_mean{s}(:,iAss)-rip_FSC_SEM{s}(:,iAss),...
                       tb-10,'k')
    end
end        
%% plot Pre
t_temp = FRcont.Pre{1, 1}.Tmtx(1:end-1);

for s= 1:3
    for iAss = 1:size(FAcont.Pre.FSC{s},2)
        for iRip = 1:length(rip.PreTimes)
            time_ = find(abs(t_temp-rip.PreTimes(iRip))==min(abs(t_temp-rip.PreTimes(iRip))));
            try 
            
                rip_FSC{s}{iAss}(:,iRip) = FAcont.Pre.FSC{s}([time_-200:time_+200],iAss);
                
            end
        end
        %rip_FSC_mean{s}(:,iAss)= nanmean(rip_FSC{s}{iAss}-nanmean(rip_FSC{s}{iAss}(1:20,:)),2);
        %rip_FSC_SEM{s}(:,iAss)= nansem(rip_FSC{s}{iAss}-nanmean(rip_FSC{s}{iAss}(1:20,:)),2);
        rip_FSC_mean{s}(:,iAss)= nanmean(rip_FSC{s}{iAss},2);
        rip_FSC_SEM{s}(:,iAss)= nansem(rip_FSC{s}{iAss},2);
    end
end
                
figure;
for s= 1:3
                    subplot(1,3,s); hold on

    for iAss = 1:size(FAcont.Pre.FSC{s},2)
                
                plot(tb-10, rip_FSC_mean{s}(:,iAss),'k')
                ciplot(rip_FSC_mean{s}(:,iAss)+rip_FSC_SEM{s}(:,iAss),...
                       rip_FSC_mean{s}(:,iAss)-rip_FSC_SEM{s}(:,iAss),...
                       tb-10,'k')
    end
end
%% Plot Post
t_temp = FRcont.Post{1, 1}.Tmtx(1:end-1);

for s= 1:3
    for iAss = 1:size(FAcont.Post.FSC{s},2)
        for iRip = 1:length(rip.PostTimes)
            time_ = find(abs(t_temp-rip.PostTimes(iRip))==min(abs(t_temp-rip.PostTimes(iRip))));
            try 
            
                rip_FSC{s}{iAss}(:,iRip) = FAcont.Post.FSC{s}([time_-200:time_+200],iAss);
                
            end
        end
        %rip_FSC_mean{s}(:,iAss)= nanmean(rip_FSC{s}{iAss}-nanmean(rip_FSC{s}{iAss}(1:20,:)),2);
        %rip_FSC_SEM{s}(:,iAss)= nansem(rip_FSC{s}{iAss}-nanmean(rip_FSC{s}{iAss}(1:20,:)),2);
        rip_FSC_mean{s}(:,iAss)= nanmean(rip_FSC{s}{iAss},2);
        rip_FSC_SEM{s}(:,iAss)= nansem(rip_FSC{s}{iAss},2);
    end
end
                
figure;
for s= 1:3
                    subplot(1,3,s); hold on

    for iAss = 1:size(FAcont.Pre.FSC{s},2)
                
                plot(tb-10, rip_FSC_mean{s}(:,iAss),'k')
                ciplot(rip_FSC_mean{s}(:,iAss)+rip_FSC_SEM{s}(:,iAss),...
                       rip_FSC_mean{s}(:,iAss)-rip_FSC_SEM{s}(:,iAss),...
                       tb-10,'k')
    end
end
            
    
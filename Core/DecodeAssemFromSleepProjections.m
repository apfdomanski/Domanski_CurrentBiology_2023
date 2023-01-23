% Takes factor model computed for continuous data for the task epoch of experiment, then
% cast model onto cutout segments for decoding analysis
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
%% Preample, load 
tic
clear all 
close all
P.target = 'Task';
% P.rat    = 'JaroslawLONG1';  %*   Post sleep import error
% P.rat    = 'JaroslawLONG2';  %*
% P.rat    = 'KrzesimirLONG1'; %*
% P.rat    = 'KrzesimirLONG2'; %*   No factors!
% P.rat    = 'KrzysztofLONG1'; %*
% P.rat    = 'KrzysztofLONG2'; %*
% P.rat    = 'MiroslawLONG0';  %*
% P.rat    = 'MiroslawLONG2';  %*   No Post sleep recorded, No factors!
% P.rat    = 'NorbertLONG1';   %*
% P.rat    = 'NorbertLONG2';   %*   Post-sleep error?
% P.rat    = 'OnufryLONG1';    %*
P.rat    = 'OnufryLONG2';    %*   Post sleep error

P.flag.noPostSleep          = true;

P.flag.CutoutLeverPresses   = true;
P.flag.PlotVerbose          = false;
P.flag.useTaskMeanFR        = true;
P.flag.sideSpecificSCdecode = false;
P.flag.runSingleTrialDecode = false;

P.expt   = 'LONG';
P.sep='\';

pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\';        % location of the processed sleep assembly activations
pat{2}  = 'C:\Analysis\AssemblyAnalysis\raw';           % location of raw spike times
pat{3} = [pat{2} P.sep 'KDE_bins'];                     % location of calculated firing P.rates for Task period
cd(pat{1})

FileList.A(1) = dir([pat{1} '*' P.rat '*PFC*' 'PreSleep.mat']);
FileList.A(2) = dir([pat{1} '*' P.rat '*HP*'  'PreSleep.mat']);
FileList.B(1) = dir([pat{1} '*' P.rat '*PFC*'  'Task.mat']);
FileList.B(2) = dir([pat{1} '*' P.rat '*HP*'   'Task.mat']);
try 
    FileList.C(1) = dir([pat{1} '*' P.rat '*PFC*'  'PostSleep.mat']);
    FileList.C(2) = dir([pat{1} '*' P.rat '*HP*'   'PostSleep.mat']);
catch 
    FileList.C(1) = cell(1,1);
    FileList.C(2) = cell(1,1);
    P.flag.noPostSleep = true;
end

P.names={'PFC','HP','Joint HP~PFC'};
set(0,'DefaultFigureWindowStyle','docked')
%% Load FA results from the FA continuous task dataset
fnOut  = dir ([pat{1} P.target P.sep '*' P.rat '*FSC.mat*']); % Factor analysis results
load([pat{1} P.target P.sep fnOut(1).name],'FSCsel','ciHld','ciHsc','units','unit_IDs','keepers');
% load([pat{1} P.target P.sep fnOut(3).name],'FSC','ciHld','ciHsc');
fnOut  = dir ([pat{1} P.target P.sep '*' P.rat '*AssemRes2.mat*']); % Factor analysis results
FAcont.Task = load([pat{1} P.target P.sep fnOut(1).name],'FL','psix','nassem','TmtxS');
FAcont.Task.FSC = FSCsel; 
FAcont.Task.ciHld = ciHld; 
FAcont.Task.ciHsc = ciHsc; 
FAcont.Task.units = units; 
FAcont.Task.unitIDs = unit_IDs; 
FAcont.Task.keepers = keepers;
FAcont.Task.TmtxS = cell2mat(FAcont.Task.TmtxS{1});

clear FSCsel ciHld ciHsc fnOut units unit_IDs keepers

% check joint assemblies for validity
reject_ = [];
FAcont.Task.unitIDs{3} = [FAcont.Task.unitIDs{1},FAcont.Task.unitIDs{2}+max(FAcont.Task.unitIDs{1})];
for iAss = 1:length(FAcont.Task.units{3})
    units_ = FAcont.Task.units{3}{iAss};
    %[sum(units_<=length(FAcont.Task.unitIDs{1})), sum(units_>length(FAcont.Task.unitIDs{1}))]
    if sum(units_<=length(FAcont.Task.unitIDs{1}))<1 || sum(units_>length(FAcont.Task.unitIDs{1}))<1
     reject_ = [reject_ , iAss];
     disp(sprintf('Removing factor %d\n',iAss))
    end
end
FAcont.Task.units{3}(reject_) = [];
FAcont.Task.FSC{3}(:,reject_)=[];
FAcont.Task.keepers{3}(reject_) = [];
% FAcont.Task.FL{3}{FAcont.Task.nassem{3}(3)}(:,reject_) = [];
clear iAss reject_ units_
%% Load firing rate matrices from cutout task dataset 
P.twin     = 10;   % Time window
P.minFR    = 0.1;  % minimal acceptable firing P.rate
P.critCvBW = 1e6;  % critical max bound on kernel width drift over time (variance)
P.bw       = 0.05;

% Enforce that Unit IDs are consistent between two analyses
% (1) get unit IDs from continuous dataset
% fn = [pat{1} P.sep P.rat '_PFC_iFR' num2str(P.bw*1e3) '_' P.target '.mat'];
% [~,~,FAcont.Task.unitIDs] = SelCellsSleep(fn,P.twin,P.minFR,P.critCvBW,'iFR');
% 18/10/16 - FAcont.Task.unitIDs now comes through with ../Task/[Name]_Task_FSC.mat

% (2) load pre-specified units from cutout dataset
fn = [pat{3} P.sep P.rat '_PFC_iFR' num2str(P.bw*1e3) '.mat'];
[FRtrials.TmtxS,...
 FRtrials.iFR0, ...
 FRtrials.EvtLs, ...
 FRtrials.EvtTs, ...
 FRtrials.unitIDs] = SelTrials(fn,P.twin,FAcont.Task.unitIDs,'iFR','cont');
% plot(FRtrials.TmtxS{1}{1}-FRtrials.TmtxS{1}{1}(1),FRtrials.iFR0{1}{1})

% (3) Go again on error trials, keeping them separate JIC
[E.TmtxS,...
 E.iFR0, ...
 E.EvtLs, ...
 E.EvtTs, ...
 E.unitIDs] = SelErrorTrials(fn,P.twin,FAcont.Task.unitIDs,'iFR','cont');

FRtrials.Err = E;
clear fn E

% reconstruct parameters for joint assems
FAcont.Task.unitIDs{3}    = [FAcont.Task.unitIDs{1},FAcont.Task.unitIDs{2}+max(FAcont.Task.unitIDs{1})];
FRtrials.unitIDs{3}       = FAcont.Task.unitIDs{3};
FRtrials.iFR0{3}          = cellfun(@horzcat,FRtrials.iFR0{1},FRtrials.iFR0{2},'UniformOutput',false);
FRtrials.Err.iFR0{3}      = cellfun(@horzcat,FRtrials.Err.iFR0{1},FRtrials.Err.iFR0{2},'UniformOutput',false);

P.nUnits   = cellfun(@length,FAcont.Task.unitIDs);
% P.nFactors =  [FAcont.Task.nassem{1}(3), FAcont.Task.nassem{2}(3), FAcont.Task.nassem{3}(3)];
% 18/10/16 - Corrected to exclude single-unit factors
% 25/1/17  - Corrected to exclude incorrect joint-assemblies
P.nFactors = cell2mat(cellfun(@numel,FAcont.Task.units,'UniformOutput',false));
disp(['No. Units = ' num2str(P.nUnits)])
disp(['No. Factors = ' num2str(P.nFactors)])
%% Load {Pre/Task/Post] epoch Firing Rates
disp(['*** loading Continous FR dataset...'])
for s=1:2
    FRcont.Pre{s}   = load([pat{1} P.sep FileList.A(s).name]);
    FRcont.Task{s}  = load([pat{1} P.sep FileList.B(s).name]);
    if P.flag.noPostSleep
        FRcont.Post{s}  =[];
    else
        FRcont.Post{s}  = load([pat{1} P.sep FileList.C(s).name]);
    end
end
FRcont.Pre{3}.iFR  = [FRcont.Pre{1}.iFR,FRcont.Pre{2}.iFR];
FRcont.Task{3}.iFR = [FRcont.Task{1}.iFR,FRcont.Task{2}.iFR];
if P.flag.noPostSleep
    FRcont.Post{3}.iFR =[];
else
    FRcont.Post{3}.iFR = [FRcont.Post{1}.iFR,FRcont.Post{2}.iFR];
end
for s = 1:3
    FRcont.Pre{s}.iFR_   = FRcont.Pre{s}.iFR(:,FAcont.Task.unitIDs{s});
    FRcont.Task{s}.iFR_  = FRcont.Task{s}.iFR(:,FAcont.Task.unitIDs{s});
    if P.flag.noPostSleep
        FRcont.Post{s}.iFR_ =[];
    else
        FRcont.Post{s}.iFR_   = FRcont.Post{s}.iFR(:,FAcont.Task.unitIDs{s});
    end
end
%% Reconstruct continuous FA model (FAcont) for [Pre/Post] from Task epoch
% i.e. take model computed for continous Task epoch, force onto unit firing
% rate cutouts in the Task dataset

% Factor model:  FR = meanFR + FL * FSC + Psi     ...Or: FSC = (FR-meanFR-Psi)/FL

% P.flag.useTaskMeanFR : Flag to 

if ~P.flag.useTaskMeanFR % Use mean FR in each period independently
    for s = 1:3
        FAcont.Pre.FSC{s}  =  (FRcont.Pre{s}.iFR_ - ...                                                % FR
            repmat((mean(FRcont.Pre{s}.iFR_) - ...                                                     % meanFR
            FAcont.Task.psix{s}{FAcont.Task.nassem{s}(3)}'),size(FRcont.Pre{s}.iFR_,1),1))...          % Psi
            / FAcont.Task.FL{s}{FAcont.Task.nassem{s}(3)}';                                            % FL
        if P.flag.noPostSleep
            FAcont.Post.FSC{s} =[];
        else
            FAcont.Post.FSC{s} =  (FRcont.Post{s}.iFR_ - ...                                           % FR
                repmat((mean(FRcont.Post{s}.iFR_) - ...                                                % meanFR
                FAcont.Task.psix{s}{FAcont.Task.nassem{s}(3)}'),size(FRcont.Post{s}.iFR_,1),1))...     % Psi
                / FAcont.Task.FL{s}{FAcont.Task.nassem{s}(3)}';
        end
        FAcont.Pre.FSC{s}  = FAcont.Pre.FSC{s}(:,FAcont.Task.keepers{s});
        if P.flag.noPostSleep
            FAcont.Post.FSC{s} = [];
        else
            FAcont.Post.FSC{s} = FAcont.Post.FSC{s}(:,FAcont.Task.keepers{s});
        end
    end
else %use mean FRs from Task epoch for each condition
    for s = 1:3
        FAcont.Pre.FSC{s}  =  (FRcont.Pre{s}.iFR_ - ...                                                % FR
            repmat((mean(FRcont.Task{s}.iFR_) - ...                               % meanFR
            FAcont.Task.psix{s}{FAcont.Task.nassem{s}(3)}'),size(FRcont.Pre{s}.iFR_,1),1))... % Psi
            / FAcont.Task.FL{s}{FAcont.Task.nassem{s}(3)}';                                        % FL
        if P.flag.noPostSleep
            FAcont.Post.FSC{s} =[];
        else
            FAcont.Post.FSC{s} =  (FRcont.Post{s}.iFR_ - ...                                                % FR
                repmat((mean(FRcont.Task{s}.iFR_) - ...                               % meanFR
                FAcont.Task.psix{s}{FAcont.Task.nassem{s}(3)}'),size(FRcont.Post{s}.iFR_,1),1))... % Psi
                / FAcont.Task.FL{s}{FAcont.Task.nassem{s}(3)}';
        end
        FAcont.Pre.FSC{s}  = FAcont.Pre.FSC{s}(:,FAcont.Task.keepers{s});
        if P.flag.noPostSleep
            FAcont.Post.FSC{s} = [];
        else
            FAcont.Post.FSC{s} = FAcont.Post.FSC{s}(:,FAcont.Task.keepers{s});
        end
    end
end
%% Reconstruct cutout FA model (FAtrials) on based on continuous data
% i.e. take model computed for continous Task epoch, force onto unit firing
% rate cutouts in the Task dataset

% Factor model:  FR = meanFR + FL * FSC + Psi     ...Or: FSC = (FR-meanFR-Psi)/FL

for s = 1:3
        for iTrial=1:length(FRtrials.iFR0{s})
            % FSC = (FR-meanFR-Psi)/FL
            FAtrials.FSC{s}{iTrial} =    (FRtrials.iFR0{s}{iTrial} - ...                                                % FR
                                          repmat((mean(cell2mat(FRtrials.iFR0{s}')) - ...                               % meanFR
                                          FAcont.Task.psix{s}{FAcont.Task.nassem{s}(3)}'),size(FRtrials.iFR0{s}{iTrial},1),1))... % Psi
                                           / FAcont.Task.FL{s}{FAcont.Task.nassem{s}(3)}';                                        % FL    
        end
end
clear s iTrial

% Excise an equal duration section for each trial, keep only factors that
% meet assembly criteria (must have >1 unit)
FRtrials.TmtxS_ = cell(1,2); 
FRtrials.iFR0_  = cell(1,2); 
FAtrials.FSC_   = cell(1,2); 
cutout=round((P.twin+5)/P.bw); 
P.Ltr=cutout;
for s=1:3
    for iTrial=1:length(FRtrials.TmtxS{1})
            %disp('...Excising delay period')
            FRtrials.TmtxS_{s}{iTrial} = FRtrials.TmtxS{1}{iTrial}([1:cutout end-cutout+1:end]);   
            FRtrials.iFR0_{s}{iTrial}  = FRtrials.iFR0{s}{iTrial}([1:cutout end-cutout+1:end],:);
            FAtrials.FSC_{s}{iTrial}   = FAtrials.FSC{s}{iTrial}([1:cutout end-cutout+1:end],FAcont.Task.keepers{s});
            
            %disp('...Keeping delay period')
            FRtrials.DelayInc.TmtxS_{s}{iTrial}  = FRtrials.TmtxS{1}{iTrial};
            FRtrials.DelayInc.iFR0_{s}{iTrial}   = FRtrials.iFR0{s}{iTrial}; 
            FAtrials.DelayInc.FSC_{s}{iTrial}    = FAtrials.FSC{s}{iTrial}(:,FAcont.Task.keepers{s});  
            
            FRtrials.DelayOnly.TmtxS_{s}{iTrial} = FRtrials.TmtxS{1}{iTrial}([(round(cutout/2)+1):(end-round(cutout/2))]);   
            FRtrials.DelayOnly.iFR0_{s}{iTrial}  = FRtrials.iFR0{s}{iTrial}([(round(cutout/2)+1):(end-round(cutout/2))],:);
            FAtrials.DelayOnly.FSC_{s}{iTrial}   = FAtrials.FSC{s}{iTrial}([(round(cutout/2)+1):(end-round(cutout/2))],FAcont.Task.keepers{s});
    end;    
end;
FRtrials.Tmtx=cell2mat(FRtrials.TmtxS_{1})';
FAtrials.units = FAcont.Task.units;
FAtrials.unitIDs = FAcont.Task.unitIDs;
FAtrials.keepers = FAcont.Task.keepers;

% iUnit = 8
% figure
% for s=1:2
%     for iTrial=1:length(FRtrials.TmtxS{1})
%             subplot(1,2,s); hold on
%             tempx=FRtrials.DelayInc.TmtxS_{s}{iTrial};
%             tempx=tempx-tempx(cutout/2);
%             plot(tempx,iTrial*40+FRtrials.DelayInc.iFR0_{s}{iTrial}(:,iUnit))
%             scatter([0, tempx(end-cutout/2)],[iTrial*40 iTrial*40], '.k')
%     end
% end
% clear s iTrial tempx
%% Reconstruct cutout FA model (FAtrials) on based on continuous data - Error trials

for s = 1:3
        for iTrial=1:length(FRtrials.Err.iFR0{s})
            % FSC = (FR-meanFR-Psi)/FL
            FAtrials.Err.FSC{s}{iTrial} =    (FRtrials.Err.iFR0{s}{iTrial} - ...                                                % FR
                                          repmat((mean(cell2mat(FRtrials.Err.iFR0{s}')) - ...                               % meanFR
                                          FAcont.Task.psix{s}{FAcont.Task.nassem{s}(3)}'),size(FRtrials.Err.iFR0{s}{iTrial},1),1))... % Psi
                                           / FAcont.Task.FL{s}{FAcont.Task.nassem{s}(3)}';                                        % FL    
        end
end
clear s iTrial

% Excise an equal duration section for each trial, keep only factors that
% meet assembly criteria (must have >1 unit)
FRtrials.Err.TmtxS_ = cell(1,2); 
FRtrials.Err.iFR0_  = cell(1,2); 
FAtrials.Err.FSC_   = cell(1,2); 
cutout=round((P.twin+5)/P.bw); 
P.Ltr=cutout;
for s=1:3
    for iTrial=1:length(FRtrials.Err.TmtxS{1})
            %disp('...Excising delay period')
            FRtrials.Err.TmtxS_{s}{iTrial} = FRtrials.Err.TmtxS{1}{iTrial}([1:cutout end-cutout+1:end]);   
            FRtrials.Err.iFR0_{s}{iTrial}  = FRtrials.Err.iFR0{s}{iTrial}([1:cutout end-cutout+1:end],:);
            FAtrials.Err.FSC_{s}{iTrial}   = FAtrials.Err.FSC{s}{iTrial}([1:cutout end-cutout+1:end],FAcont.Task.keepers{s});
            
            %disp('...Keeping delay period')
            FRtrials.Err.DelayInc.TmtxS_{s}{iTrial}  = FRtrials.Err.TmtxS{1}{iTrial};
            FRtrials.Err.DelayInc.iFR0_{s}{iTrial}   = FRtrials.Err.iFR0{s}{iTrial}; 
            FAtrials.Err.DelayInc.FSC_{s}{iTrial}    = FAtrials.Err.FSC{s}{iTrial}(:,FAcont.Task.keepers{s});  
            
            FRtrials.Err.DelayOnly.TmtxS_{s}{iTrial} = FRtrials.Err.TmtxS{1}{iTrial}([(round(cutout/2)+1):(end-round(cutout/2))]);   
            FRtrials.Err.DelayOnly.iFR0_{s}{iTrial}  = FRtrials.Err.iFR0{s}{iTrial}([(round(cutout/2)+1):(end-round(cutout/2))],:);
            FAtrials.Err.DelayOnly.FSC_{s}{iTrial}   = FAtrials.Err.FSC{s}{iTrial}([(round(cutout/2)+1):(end-round(cutout/2))],FAcont.Task.keepers{s});
    end;    
end;
FRtrials.Err.Tmtx=cell2mat(FRtrials.Err.TmtxS_{1})';
FAtrials.Err.units = FAcont.Task.units;
FAtrials.Err.unitIDs = FAcont.Task.unitIDs;
FAtrials.Err.keepers = FAcont.Task.keepers;
%% Trial-outcome sorted activations - correct

% remap trial outcome labels
P.ntr = length(FRtrials.EvtTs)/2;
FRtrials.evt0 = FRtrials.EvtLs(((1:P.ntr)-1)*2+1)-2;

FRtrials.iFR0_Sorted = cell(3,1);
FRtrials.iFR0_Mean   = cell(3,1);
FRtrials.iFR0_SEM    = cell(3,1);
FAtrials.FSC_Sorted  = cell(3,1);
FAtrials.FSC_Mean    = cell(3,1);
FAtrials.FSC_SEM     = cell(3,1);

% Lever-press cutouts
tempSize = cell2mat(cellfun(@size,FRtrials.iFR0_{1},'UniformOutput',false)');
tempSize_= min(tempSize(:,1));
P.Ltr_cont = tempSize;
P.Ltr_shortest = tempSize_;

% for trials with delay period still in, align the timebase to that of the shortest trial length of shortest trial
tempSize = cell2mat(cellfun(@size,FRtrials.DelayInc.iFR0_{1},'UniformOutput',false)');
tempSize_= max(tempSize(:,1));
FRtrials.DelayInc.Ltr_cont = tempSize;
FRtrials.DelayInc.Ltr_longest = tempSize_;

tempSize = cell2mat(cellfun(@size,FRtrials.DelayOnly.iFR0_{1},'UniformOutput',false)');
tempSize_= max(tempSize(:,1));
FRtrials.DelayOnly.Ltr_cont = tempSize;
FRtrials.DelayOnly.Ltr_longest = tempSize_;

for iOutcome = 1:2
    trials = find(FRtrials.evt0 == iOutcome);
    for s=1:3
        
        % Mean unit firing rates
        tempFR   = FRtrials.iFR0_{s}(trials);
        tempFR_  = FRtrials.DelayInc.iFR0_{s}(trials);
        tempFR__ = FRtrials.DelayOnly.iFR0_{s}(trials);
        for iUnit=1:size(tempFR{1},2)
            for iTrial=1:length(trials)
                     FRtrials.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial)=tempFR{iTrial}(:,iUnit);
                     
                     % Right-pad with NaNs for any dataset with varialble trial lengths
                     tempCutout = tempFR_{iTrial}(:,iUnit);
                     FRtrials.DelayInc.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial)= ...
                        [tempCutout;nan(FRtrials.DelayInc.Ltr_longest-length(tempCutout),1)]; 
                    
                     tempCutout = tempFR__{iTrial}(:,iUnit);
                     FRtrials.DelayOnly.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial)= ...
                        [tempCutout;nan(FRtrials.DelayOnly.Ltr_longest-length(tempCutout),1)];
                    
                    
            end
            FRtrials.iFR0_Mean{s}{iOutcome}(:,iUnit) = mean (FRtrials.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            FRtrials.iFR0_SEM{s}{iOutcome} (:,iUnit) = nansem(FRtrials.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            
            FRtrials.DelayInc.iFR0_Mean{s}{iOutcome}(:,iUnit) = nanmean (FRtrials.DelayInc.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            FRtrials.DelayInc.iFR0_SEM{s}{iOutcome} (:,iUnit) = nansem(FRtrials.DelayInc.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            
            FRtrials.DelayOnly.iFR0_Mean{s}{iOutcome}(:,iUnit) = nanmean (FRtrials.DelayOnly.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            FRtrials.DelayOnly.iFR0_SEM{s}{iOutcome} (:,iUnit) = nansem(FRtrials.DelayOnly.iFR0_Sorted{s}{iOutcome}{iUnit},2);
        end
        
        % Mean factor scores
        tempFA   = FAtrials.FSC_{s}(trials);
        tempFA_  = FAtrials.DelayInc.FSC_{s}(trials);
        tempFA__ = FAtrials.DelayOnly.FSC_{s}(trials);
        for iFactor=1:size(tempFA{1},2)
            for iTrial=1:length(trials)
                
                FAtrials.FSC_Sorted{s}{iOutcome}{iFactor}(:,iTrial)=tempFA{iTrial}(:,iFactor);
                
                % Right-pad with NaNs for any dataset with varialble trial lengths
                tempCutout = tempFA_{iTrial}(:,iFactor);
                FAtrials.DelayInc.FSC_Sorted{s}{iOutcome}{iFactor}(:,iTrial)  =   [tempCutout;nan(FRtrials.DelayInc.Ltr_longest-length(tempCutout),1)];
                
                tempCutout = tempFA__{iTrial}(:,iFactor);
                FAtrials.DelayOnly.FSC_Sorted{s}{iOutcome}{iFactor}(:,iTrial) =   [tempCutout;nan(FRtrials.DelayOnly.Ltr_longest-length(tempCutout),1)];
                    
            end
            FAtrials.FSC_Mean{s}{iOutcome}(:,iFactor) = mean  (FAtrials.FSC_Sorted{s}{iOutcome}{iFactor},2);
            FAtrials.FSC_SEM{s}{iOutcome} (:,iFactor) = nansem(FAtrials.FSC_Sorted{s}{iOutcome}{iFactor},2);
            
            FAtrials.DelayInc.FSC_Mean{s}{iOutcome}(:,iFactor) = nanmean(FAtrials.DelayInc.FSC_Sorted{s}{iOutcome}{iFactor},2);
            FAtrials.DelayInc.FSC_SEM{s}{iOutcome} (:,iFactor) = nansem (FAtrials.DelayInc.FSC_Sorted{s}{iOutcome}{iFactor},2);
            
            FAtrials.DelayOnly.FSC_Mean{s}{iOutcome}(:,iFactor) = nanmean(FAtrials.DelayOnly.FSC_Sorted{s}{iOutcome}{iFactor},2);
            FAtrials.DelayOnly.FSC_SEM{s}{iOutcome} (:,iFactor) = nansem (FAtrials.DelayOnly.FSC_Sorted{s}{iOutcome}{iFactor},2);

        end
        
    end
end

% figure; hold on
% plot(staggerplot(FAtrials.DelayInc.FSC_Mean{1}{1},0,10),'b')
% plot(staggerplot(FAtrials.DelayInc.FSC_Mean{1}{2},0,10),'r')
%% Trial-outcome sorted activations - errors

% remap trial outcome labels
P.ntrErr = length(FRtrials.Err.EvtTs)/2;
FRtrials.Err.evt0 = FRtrials.Err.EvtLs(((1:P.ntrErr)-1)*2+1)/10-2;

FRtrials.Err.iFR0_Sorted = cell(3,1);
FRtrials.Err.iFR0_Mean   = cell(3,1);
FRtrials.Err.iFR0_SEM    = cell(3,1);
FAtrials.Err.FSC_Sorted  = cell(3,1);
FAtrials.Err.FSC_Mean    = cell(3,1);
FAtrials.Err.FSC_SEM     = cell(3,1);

% Lever-press cutouts
tempSize = cell2mat(cellfun(@size,FRtrials.Err.iFR0_{1},'UniformOutput',false)');
tempSize_= min(tempSize(:,1));
P.Ltr_cont = tempSize;
P.Ltr_shortest = tempSize_;

% for trials with delay period still in, align the timebase to that of the shortest trial length of shortest trial
tempSize = cell2mat(cellfun(@size,FRtrials.Err.DelayInc.iFR0_{1},'UniformOutput',false)');
tempSize_= max(tempSize(:,1));
FRtrials.Err.DelayInc.Ltr_cont = tempSize;
FRtrials.Err.DelayInc.Ltr_longest = tempSize_;

tempSize = cell2mat(cellfun(@size,FRtrials.Err.DelayOnly.iFR0_{1},'UniformOutput',false)');
tempSize_= max(tempSize(:,1));
FRtrials.Err.DelayOnly.Ltr_cont = tempSize;
FRtrials.Err.DelayOnly.Ltr_longest = tempSize_;

for iOutcome = 1:2
    trials = find(FRtrials.Err.evt0 == iOutcome);
    for s=1:3
        
        % Mean unit firing rates
        tempFR   = FRtrials.Err.iFR0_{s}(trials);
        tempFR_  = FRtrials.Err.DelayInc.iFR0_{s}(trials);
        tempFR__ = FRtrials.Err.DelayOnly.iFR0_{s}(trials);
        for iUnit=1:size(tempFR{1},2)
            for iTrial=1:length(trials)
                     FRtrials.Err.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial)=tempFR{iTrial}(:,iUnit);
                     
                     % Right-pad with NaNs for any dataset with varialble trial lengths
                     tempCutout = tempFR_{iTrial}(:,iUnit);
                     FRtrials.Err.DelayInc.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial)= ...
                        [tempCutout;nan(FRtrials.Err.DelayInc.Ltr_longest-length(tempCutout),1)]; 
                    
                     tempCutout = tempFR__{iTrial}(:,iUnit);
                     FRtrials.Err.DelayOnly.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial)= ...
                        [tempCutout;nan(FRtrials.Err.DelayOnly.Ltr_longest-length(tempCutout),1)];
                    
                    
            end
            FRtrials.Err.iFR0_Mean{s}{iOutcome}(:,iUnit) = mean (FRtrials.Err.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            FRtrials.Err.iFR0_SEM{s}{iOutcome} (:,iUnit) = nansem(FRtrials.Err.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            
            FRtrials.Err.DelayInc.iFR0_Mean{s}{iOutcome}(:,iUnit) = nanmean (FRtrials.Err.DelayInc.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            FRtrials.Err.DelayInc.iFR0_SEM{s}{iOutcome} (:,iUnit) = nansem(FRtrials.Err.DelayInc.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            
            FRtrials.Err.DelayOnly.iFR0_Mean{s}{iOutcome}(:,iUnit) = nanmean (FRtrials.Err.DelayOnly.iFR0_Sorted{s}{iOutcome}{iUnit},2);
            FRtrials.Err.DelayOnly.iFR0_SEM{s}{iOutcome} (:,iUnit) = nansem(FRtrials.Err.DelayOnly.iFR0_Sorted{s}{iOutcome}{iUnit},2);
        end
        
        % Mean factor scores
        tempFA   = FAtrials.Err.FSC_{s}(trials);
        tempFA_  = FAtrials.Err.DelayInc.FSC_{s}(trials);
        tempFA__ = FAtrials.Err.DelayOnly.FSC_{s}(trials);
        for iFactor=1:size(tempFA{1},2)
            for iTrial=1:length(trials)
                
                FAtrials.Err.FSC_Sorted{s}{iOutcome}{iFactor}(:,iTrial)=tempFA{iTrial}(:,iFactor);
                
                % Right-pad with NaNs for any dataset with varialble trial lengths
                tempCutout = tempFA_{iTrial}(:,iFactor);
                FAtrials.Err.DelayInc.FSC_Sorted{s}{iOutcome}{iFactor}(:,iTrial)  =   [tempCutout;nan(FRtrials.Err.DelayInc.Ltr_longest-length(tempCutout),1)];
                
                tempCutout = tempFA__{iTrial}(:,iFactor);
                FAtrials.Err.DelayOnly.FSC_Sorted{s}{iOutcome}{iFactor}(:,iTrial) =   [tempCutout;nan(FRtrials.Err.DelayOnly.Ltr_longest-length(tempCutout),1)];
                    
            end
            FAtrials.Err.FSC_Mean{s}{iOutcome}(:,iFactor) = mean  (FAtrials.Err.FSC_Sorted{s}{iOutcome}{iFactor},2);
            FAtrials.Err.FSC_SEM{s}{iOutcome} (:,iFactor) = nansem(FAtrials.Err.FSC_Sorted{s}{iOutcome}{iFactor},2);
            
            FAtrials.Err.DelayInc.FSC_Mean{s}{iOutcome}(:,iFactor) = nanmean(FAtrials.Err.DelayInc.FSC_Sorted{s}{iOutcome}{iFactor},2);
            FAtrials.Err.DelayInc.FSC_SEM{s}{iOutcome} (:,iFactor) = nansem (FAtrials.Err.DelayInc.FSC_Sorted{s}{iOutcome}{iFactor},2);
            
            FAtrials.Err.DelayOnly.FSC_Mean{s}{iOutcome}(:,iFactor) = nanmean(FAtrials.Err.DelayOnly.FSC_Sorted{s}{iOutcome}{iFactor},2);
            FAtrials.Err.DelayOnly.FSC_SEM{s}{iOutcome} (:,iFactor) = nansem (FAtrials.Err.DelayOnly.FSC_Sorted{s}{iOutcome}{iFactor},2);

        end
        
    end
end
% 
% figure; hold on
% plot(staggerplot(FAtrials.Err.DelayInc.FSC_Mean{1}{1},0,10),'b')
% plot(staggerplot(FAtrials.Err.DelayInc.FSC_Mean{1}{2},0,10),'r')
%% Decoding on delay period (experimental...)
% sum and mean FR in delay period
for iOutcome = 1:2
    trials = find(FRtrials.evt0 == iOutcome);
    for s=1:3
        for iTrial=1:length(trials)
            
            for iUnit =1:numel(FRtrials.DelayOnly.iFR0_Sorted{s}{iOutcome})
                FRtrials.DelayOnly.iFR0_meanTrials{s}{iOutcome}(iUnit,iTrial)  = nanmean(FRtrials.DelayOnly.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial));
                FRtrials.DelayOnly.iFR0_sumTrials{s} {iOutcome}(iUnit,iTrial)  = nansum(FRtrials.DelayOnly.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial))./ ...
                                                                                   sum(~isnan(FRtrials.DelayOnly.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial)));
            end
        end
    end
end
rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
for s=1:3
	for iUnit =1:numel(FRtrials.DelayOnly.iFR0_Sorted{s}{1})

       [~,FRtrials.DelayOnly.iFR0_meanTrials_ts{s}(iUnit)] = ttest2( FRtrials.DelayOnly.iFR0_meanTrials{s}{1}(iUnit,:),FRtrials.DelayOnly.iFR0_meanTrials{s}{2}(iUnit,:));
       [~,FRtrials.DelayOnly.iFR0_sumTrials_ts{s}(iUnit)]  = ttest2( FRtrials.DelayOnly.iFR0_sumTrials{s}{1}(iUnit,:),FRtrials.DelayOnly.iFR0_meanTrials{s}{2}(iUnit,:));
    end
end
addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
% clear iOutcome trials s iUnit iFactor iTrial temp tempSize tempSize_ cutout tempCutout tempFA tempFA_ tempFA__ tempFR tempFR_ tempFR__
%  figure;
%  scatter (FRtrials.DelayOnly.iFR0_meanTrials_ts{1},...
%           FRtrials.DelayOnly.iFR0_sumTrials_ts{1})
%% plotting optional
if P.flag.PlotVerbose
%% plot unit L/R profiles - subplots
if P.flag.CutoutLeverPresses
    for s=1:2
        figure('Name',P.names{s});
        for iUnit=1:P.nUnits(s)
            subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
            axis off; hold on
            %{area}{outcome}[time,Unit]
            ciplot(FRtrials.iFR0_Mean{s}{1}(:,iUnit)+FRtrials.iFR0_SEM{s}{1}(:,iUnit),...
                   FRtrials.iFR0_Mean{s}{1}(:,iUnit)-FRtrials.iFR0_SEM{s}{1}(:,iUnit),...
                   (1:P.Ltr*2)*P.bw,'r')

            ciplot(FRtrials.iFR0_Mean{s}{2}(:,iUnit)+FRtrials.iFR0_SEM{s}{2}(:,iUnit),...
                   FRtrials.iFR0_Mean{s}{2}(:,iUnit)-FRtrials.iFR0_SEM{s}{2}(:,iUnit),...
                   (1:P.Ltr*2)*P.bw,'g')

            plot((1:P.Ltr*2)*P.bw,FRtrials.iFR0_Mean{s}{1}(:,iUnit),'r','LineWidth',1.5)
            plot((1:P.Ltr*2)*P.bw,FRtrials.iFR0_Mean{s}{2}(:,iUnit),'g','LineWidth',1.5)

        end
    end
else
    for s=1:2
        figure('Name',P.names{s});
        for iUnit=1:P.nUnits(s)
            subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
%             axis off; 
            hold on
            %{area}{outcome}[time,Unit]
            ciplot(FRtrials.iFR0_Mean{s}{1}(:,iUnit)+FRtrials.iFR0_SEM{s}{1}(:,iUnit),...
                   FRtrials.iFR0_Mean{s}{1}(:,iUnit)-FRtrials.iFR0_SEM{s}{1}(:,iUnit),...
                   (1:P.Ltr_shortest)*P.bw,'r')

            ciplot(FRtrials.iFR0_Mean{s}{2}(:,iUnit)+FRtrials.iFR0_SEM{s}{2}(:,iUnit),...
                   FRtrials.iFR0_Mean{s}{2}(:,iUnit)-FRtrials.iFR0_SEM{s}{2}(:,iUnit),...
                   (1:P.Ltr_shortest)*P.bw,'g')

            plot((1:P.Ltr_shortest)*P.bw,FRtrials.iFR0_Mean{s}{1}(:,iUnit),'r','LineWidth',1.5)
            plot((1:P.Ltr_shortest)*P.bw,FRtrials.iFR0_Mean{s}{2}(:,iUnit),'g','LineWidth',1.5)

        end
    end    
end
%% plot unit L/R profiles - overlaid
stepSize=5;
if P.flag.CutoutLeverPresses
    for s=1:2
        figure('Name',P.names{s});
        for iUnit=1:P.nUnits(s)
            subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
            axis off; hold on
            %{area}{outcome}[time,Unit]
            ciplot(FRtrials.iFR0_Mean{s}{1}(:,iUnit)+FRtrials.iFR0_SEM{s}{1}(:,iUnit),...
                   FRtrials.iFR0_Mean{s}{1}(:,iUnit)-FRtrials.iFR0_SEM{s}{1}(:,iUnit),...
                   (1:P.Ltr*2)*P.bw,'r')

            ciplot(FRtrials.iFR0_Mean{s}{2}(:,iUnit)+FRtrials.iFR0_SEM{s}{2}(:,iUnit),...
                   FRtrials.iFR0_Mean{s}{2}(:,iUnit)-FRtrials.iFR0_SEM{s}{2}(:,iUnit),...
                   (1:P.Ltr*2)*P.bw,'g')

            plot((1:P.Ltr*2)*P.bw,FRtrials.iFR0_Mean{s}{1}(:,iUnit),'r','LineWidth',1.5)
            plot((1:P.Ltr*2)*P.bw,FRtrials.iFR0_Mean{s}{2}(:,iUnit),'g','LineWidth',1.5)

        end
    end
else
    for s=1:2
        cmap = jet(P.nUnits(s));
        figure('Name',P.names{s});
        for iUnit=1:P.nUnits(s)
%             subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
%             axis off; 
            hold on
            %{area}{outcome}[time,Unit]
            ciplot(stepSize*iUnit+zscore(FRtrials.iFR0_Mean{s}{1}(:,iUnit)+FRtrials.iFR0_SEM{s}{1}(:,iUnit)),...
                   stepSize*iUnit+zscore(FRtrials.iFR0_Mean{s}{1}(:,iUnit)-FRtrials.iFR0_SEM{s}{1}(:,iUnit)),...
                   (1:P.Ltr_shortest)*P.bw,'r')

            ciplot(stepSize*iUnit+zscore(FRtrials.iFR0_Mean{s}{2}(:,iUnit)+FRtrials.iFR0_SEM{s}{2}(:,iUnit)),...
                   stepSize*iUnit+zscore(FRtrials.iFR0_Mean{s}{2}(:,iUnit)-FRtrials.iFR0_SEM{s}{2}(:,iUnit)),...
                   (1:P.Ltr_shortest)*P.bw,'g')

            plot((1:P.Ltr_shortest)*P.bw,stepSize*iUnit+zscore(FRtrials.iFR0_Mean{s}{1}(:,iUnit)),'r','LineWidth',1.5)
            plot((1:P.Ltr_shortest)*P.bw,stepSize*iUnit+zscore(FRtrials.iFR0_Mean{s}{2}(:,iUnit)),'g','LineWidth',1.5)

        end
    end    
end
%% plot factor L/R profiles - subplots
if P.flag.CutoutLeverPresses

    for s=1:3
        figure('Name',P.names{s});
        for iFactor=1:P.nFactors(s)
            subplot(ceil(P.nFactors(s)^0.5),ceil(P.nFactors(s)^0.5),iFactor)
            axis off; hold on
            %{area}{outcome}[time,Unit]
            ciplot(FAtrials.FSC_Mean{s}{1}(:,iFactor)+FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   FAtrials.FSC_Mean{s}{1}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   (1:P.Ltr*2)*P.bw,'r')

            ciplot(FAtrials.FSC_Mean{s}{2}(:,iFactor)+FAtrials.FSC_SEM{s}{2}(:,iFactor),...
                   FAtrials.FSC_Mean{s}{2}(:,iFactor)-FAtrials.FSC_SEM{s}{2}(:,iFactor),...
                   (1:P.Ltr*2)*P.bw,'g')

            plot((1:P.Ltr*2)*P.bw,FAtrials.FSC_Mean{s}{1}(:,iFactor),'r','LineWidth',1.5)
            plot((1:P.Ltr*2)*P.bw,FAtrials.FSC_Mean{s}{2}(:,iFactor),'g','LineWidth',1.5)

        end
    end
else
    for s=1:3
        figure('Name',P.names{s});
        for iFactor=1:P.nFactors(s)
            subplot(ceil(P.nFactors(s)^0.5),ceil(P.nFactors(s)^0.5),iFactor)
            axis off; hold on
            %{area}{outcome}[time,Unit]
            ciplot(FAtrials.FSC_Mean{s}{1}(:,iFactor)+FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   FAtrials.FSC_Mean{s}{1}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   (1:P.Ltr_shortest)*P.bw,'r')

            ciplot(FAtrials.FSC_Mean{s}{2}(:,iFactor)+FAtrials.FSC_SEM{s}{2}(:,iFactor),...
                   FAtrials.FSC_Mean{s}{2}(:,iFactor)-FAtrials.FSC_SEM{s}{2}(:,iFactor),...
                   (1:P.Ltr_shortest)*P.bw,'g')

            plot((1:P.Ltr_shortest)*P.bw,FAtrials.FSC_Mean{s}{1}(:,iFactor),'r','LineWidth',1.5)
            plot((1:P.Ltr_shortest)*P.bw,FAtrials.FSC_Mean{s}{2}(:,iFactor),'g','LineWidth',1.5)

        end
    end
    
end

clear iFactor iUnit
%% plot factor L/R profiles - overlaid
stepSize=5;
if P.flag.CutoutLeverPresses
    for s=1:3
        figure('Name',P.names{s}); hold on
        for iFactor=1:P.nFactors(s)
            % subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
            % axis off
            % {area}{outcome}[time,Unit]
            ciplot(stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor)+FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   (1:P.Ltr*2)*P.bw,'r')

           ciplot(stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor)+FAtrials.FSC_SEM{s}{2}(:,iFactor),...
                   stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   (1:P.Ltr*2)*P.bw,'g')

      
             plot((1:P.Ltr*2)*P.bw,stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor),'r','LineWidth',1.5)
            plot((1:P.Ltr*2)*P.bw,stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor),'g','LineWidth',1.5)

        end
    end
else
    for s=1:3
        cmap = jet(P.nUnits(s));
        figure('Name',P.names{s});hold on
        for iFactor=1:P.nFactors(s)
            % subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
            % axis off
            % {area}{outcome}[time,Unit]
            ciplot(stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor)+FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   (1:P.Ltr_shortest)*P.bw,'r')

            ciplot(stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor)+FAtrials.FSC_SEM{s}{2}(:,iFactor),...
                   stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   (1:P.Ltr_shortest)*P.bw,'g')

      
            plot((1:P.Ltr_shortest)*P.bw,stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor),'r','LineWidth',1.5)
            plot((1:P.Ltr_shortest)*P.bw,stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor),'g','LineWidth',1.5)
        end
    end    
end
end
%% Left/Right decoding on FRs and Assemblies - correct trials
P.reg = 0.05;
for s=1:3
    
    %   Firing rate decoder - F test
    
    if s ==3 
        input = cell2mat([FRtrials.iFR0_{1}' , FRtrials.iFR0_{2}']);
    else
        input = cell2mat(FRtrials.iFR0_{s}');
    end
    [~,~,...
     D.units.Ft2x{s},...
     D.units.Rt2x{s},...
     D.units.Ft2ciL0{s},...
     D.units.Ft2ciH0{s},...
     D.units.TS{s},...
     D.units.dfnum{s},...
     D.units.dfd{s}] = DecodeStats(input,FRtrials.evt0,P.reg);
     % NB, for  Significant decoders: (tpdf(D.units.TS{s},numel(FRtrials.evt0)-2))<0.05
     
     %   Assembly decoder - F test
    if ~isempty(FAtrials.keepers{s}) %... if no assemblies in this area
    [~,~,...
     D.Assem.Ft2x{s},...
     D.Assem.Rt2x{s},...
     D.Assem.Ft2ciL0{s},...
     D.Assem.Ft2ciH0{s},...
     D.Assem.TS{s},...
     D.Assem.dfnum{s},...
     D.Assem.dfd{s}] = DecodeStats(cell2mat(FAtrials.FSC_{s}'),FRtrials.evt0,P.reg);
    else
     D.Assem.Ft2x{s}    = nan(1,P.Ltr*2);   
     D.Assem.Rt2x{s}    = nan(1,P.Ltr*2);  
     D.Assem.Ft2ciL0{s} = nan(1,P.Ltr*2);  
     D.Assem.Ft2ciH0{s} = nan(1,P.Ltr*2);  
     D.Assem.Rt2x{s}    = nan(1,P.Ltr*2);  
     D.Assem.TS{s}      = nan(P.Ltr*2,1);  
     D.Assem.dfnum{s}   = NaN;
     D.Assem.dfd{s}     = NaN;
    end
     % Peak decoding
     D.Assem.TSPeak{s} = max(D.Assem.TS{s},[],1);
     D.Assem.TSIntegral{s} = sum(D.Assem.TS{s},1);

end
%% Left/Right decoding on FRs and Assemblies - error trials
P.reg = 0.05;
for s=1:3
    
    %   Firing rate decoder - F test
    
    if s ==3 
        input = cell2mat([FRtrials.Err.iFR0_{1}' , FRtrials.Err.iFR0_{2}']);
    else
        input = cell2mat(FRtrials.Err.iFR0_{s}');
    end
    [~,~,...
     D.units.Err.Ft2x{s},...
     D.units.Err.Rt2x{s},...
     D.units.Err.Ft2ciL0{s},...
     D.units.Err.Ft2ciH0{s},...
     D.units.Err.TS{s},...
     D.units.Err.dfnum{s},...
     D.units.Err.dfd{s}] = DecodeStats(input,FRtrials.Err.evt0,P.reg);
     % NB, for  Significant decoders: (tpdf(D.units.TS{s},numel(FRtrials.evt0)-2))<0.05
     
     %   Assembly decoder - F test
    if ~isempty(FAtrials.Err.keepers{s}) %... if no assemblies in this area
    [~,~,...
     D.Assem.Err.Ft2x{s},...
     D.Assem.Err.Rt2x{s},...
     D.Assem.Err.Ft2ciL0{s},...
     D.Assem.Err.Ft2ciH0{s},...
     D.Assem.Err.TS{s},...
     D.Assem.Err.dfnum{s},...
     D.Assem.Err.dfd{s}] = DecodeStats(cell2mat(FAtrials.Err.FSC_{s}'),FRtrials.Err.evt0,P.reg);
    else
     D.Assem.Err.Ft2x{s}    = nan(1,P.Ltr*2);   
     D.Assem.Err.Rt2x{s}    = nan(1,P.Ltr*2);  
     D.Assem.Err.Ft2ciL0{s} = nan(1,P.Ltr*2);  
     D.Assem.Err.Ft2ciH0{s} = nan(1,P.Ltr*2);  
     D.Assem.Err.Rt2x{s}    = nan(1,P.Ltr*2);  
     D.Assem.Err.TS{s}      = nan(P.Ltr*2,1);  
     D.Assem.Err.dfnum{s}   = NaN;
     D.Assem.Err.dfd{s}     = NaN;
    end
     % Peak decoding
     D.Assem.Err.TSPeak{s} = max(D.Assem.Err.TS{s},[],1);
     D.Assem.Err.TSIntegral{s} = sum(D.Assem.Err.TS{s},1);

end
%% Sample/Choice decoding on FRs and Assemblies - correct trials
side = 2;
Dsc = struct;


% NB need to align lever presses first....
offset = 5/P.bw;
cutout = (P.Ltr - offset)*2;
P.Ltr_SC = cutout/2; % Length of time cutout around each sample or choice event
cutout = (1:cutout) + offset;

tempFR  = FRtrials.iFR0_;
for s = 1:3
    for iTrial = 1:length(tempFR{s})
        tempFR{s}{iTrial} = tempFR{s}{iTrial}(cutout,:);
    end
    if P.flag.sideSpecificSCdecode
        tempFR{s}(FRtrials.evt0~=side)=[];
        outcome = repmat([1 2],1,length(tempFR{s}));
    else
        outcome = repmat([1 2],1,length(FRtrials.evt0));
    end
end


for s=1:3
    
    %   Firing rate decoder 
    
    if s ==3 
        input = cell2mat([tempFR{1}' , tempFR{2}']);
    else
        input = cell2mat(tempFR{s}');
    end
    [~,~,...
     Dsc.units.Ft2x{s},...
     Dsc.units.Rt2x{s},...
     Dsc.units.Ft2ciL0{s},...
     Dsc.units.Ft2ciH0{s},...
     Dsc.units.TS{s},...
     Dsc.units.dfnum{s},...
     Dsc.units.dfd{s}] = DecodeStats(input,outcome,P.reg);
     
     %   Assembly decoder 
    if ~isempty(FAtrials.keepers{s}) %... if no assemblies in this area
    tempFSC  = FAtrials.FSC_{s};
    for iTrial = 1:length(tempFSC)
        tempFSC{iTrial} = tempFSC{iTrial}(cutout,:);
    end
    if P.flag.sideSpecificSCdecode
        tempFSC(FRtrials.evt0~=side)=[];
    end
    
    [~,~,...
     Dsc.Assem.Ft2x{s},...
     Dsc.Assem.Rt2x{s},...
     Dsc.Assem.Ft2ciL0{s},...
     Dsc.Assem.Ft2ciH0{s},...
     Dsc.Assem.TS{s},...
     Dsc.Assem.dfnum{s},...
     Dsc.Assem.dfd{s}] = DecodeStats(cell2mat(tempFSC'),outcome,P.reg);
    else
     Dsc.Assem.Ft2x{s}    = nan(1,P.Ltr_SC);   
     Dsc.Assem.Rt2x{s}    = nan(1,P.Ltr_SC);  
     Dsc.Assem.Ft2ciL0{s} = nan(1,P.Ltr_SC);  
     Dsc.Assem.Ft2ciH0{s} = nan(1,P.Ltr_SC);  
     Dsc.Assem.Rt2x{s}    = nan(1,P.Ltr_SC);  
     Dsc.Assem.TS{s}      = nan(P.Ltr_SC,1);  
     Dsc.Assem.dfnum{s}   = NaN;
     Dsc.Assem.dfd{s}     = NaN;
    end
     % Peak decoding
     Dsc.Assem.TSPeak{s}     = max(Dsc.Assem.TS{s},[],1);
     Dsc.Assem.TSIntegral{s} = sum(Dsc.Assem.TS{s},1);
end


x_temp = (1:P.Ltr_SC)*P.bw-P.Ltr_SC/2*P.bw;

if P.flag.PlotVerbose
    figure('color','w','name','Sample/Choice decoding');
    subplot(1,2,1); hold on
    ciplot(Dsc.units.Ft2ciH0{1},Dsc.units.Ft2ciL0{1},x_temp,'r')
    ciplot(Dsc.units.Ft2ciH0{2},Dsc.units.Ft2ciL0{2},x_temp,'g')
    plot(x_temp, Dsc.units.Ft2x{1},'r')
    plot(x_temp, Dsc.units.Ft2x{2},'g')
    title('Unit Sample/Choice decoding')
    axis([-Inf Inf 0 30])
    
    subplot(1,2,2); hold on
    ciplot(Dsc.Assem.Ft2ciH0{1},Dsc.Assem.Ft2ciL0{1},x_temp,'r')
    ciplot(Dsc.Assem.Ft2ciH0{2},Dsc.Assem.Ft2ciL0{2},x_temp,'g')
    ciplot(Dsc.Assem.Ft2ciH0{3},Dsc.Assem.Ft2ciL0{3},x_temp,'b')
    plot(x_temp, Dsc.Assem.Ft2x{1},'r')
    plot(x_temp, Dsc.Assem.Ft2x{2},'g')
    plot(x_temp, Dsc.Assem.Ft2x{3},'b')
    title('Assembly Sample/Choice decoding')
    axis([-Inf Inf 0 30])
end
clear s input iTrial cutout tempFSC  x_temp side
%% Sample/Choice decoding on FRs and Assemblies - correct trials
side = 2;

% NB need to align lever presses first....
offset = 5/P.bw;
cutout = (P.Ltr - offset)*2;
P.Ltr_SC = cutout/2; % Length of time cutout around each sample or choice event
cutout = (1:cutout) + offset;

tempFR  = FRtrials.Err.iFR0_;
for s = 1:3
    for iTrial = 1:length(tempFR{s})
        tempFR{s}{iTrial} = tempFR{s}{iTrial}(cutout,:);
    end
    if P.flag.sideSpecificSCdecode
        tempFR{s}(FRtrials.Err.evt0~=side)=[];
        outcome = repmat([1 2],1,length(tempFR{s}));
    else
        outcome = repmat([1 2],1,length(FRtrials.Err.evt0));
    end
end


for s=1:3
    
    %   Firing rate decoder 
    
    if s ==3 
        input = cell2mat([tempFR{1}' , tempFR{2}']);
    else
        input = cell2mat(tempFR{s}');
    end
    [~,~,...
     Dsc.units.Err.Ft2x{s},...
     Dsc.units.Err.Rt2x{s},...
     Dsc.units.Err.Ft2ciL0{s},...
     Dsc.units.Err.Ft2ciH0{s},...
     Dsc.units.Err.TS{s},...
     Dsc.units.Err.dfnum{s},...
     Dsc.units.Err.dfd{s}] = DecodeStats(input,outcome,P.reg);
     
     %   Assembly decoder 
    if ~isempty(FAtrials.Err.keepers{s}) %... if no assemblies in this area
    tempFSC  = FAtrials.Err.FSC_{s};
    for iTrial = 1:length(tempFSC)
        tempFSC{iTrial} = tempFSC{iTrial}(cutout,:);
    end
    if P.flag.sideSpecificSCdecode
        tempFSC(FRtrials.Err.evt0~=side)=[];
    end
    
    [~,~,...
     Dsc.Assem.Err.Ft2x{s},...
     Dsc.Assem.Err.Rt2x{s},...
     Dsc.Assem.Err.Ft2ciL0{s},...
     Dsc.Assem.Err.Ft2ciH0{s},...
     Dsc.Assem.Err.TS{s},...
     Dsc.Assem.Err.dfnum{s},...
     Dsc.Assem.Err.dfd{s}] = DecodeStats(cell2mat(tempFSC'),outcome,P.reg);
    else
     Dsc.Assem.Err.Ft2x{s}    = nan(1,P.Ltr_SC);   
     Dsc.Assem.Err.Rt2x{s}    = nan(1,P.Ltr_SC);  
     Dsc.Assem.Err.Ft2ciL0{s} = nan(1,P.Ltr_SC);  
     Dsc.Assem.Err.Ft2ciH0{s} = nan(1,P.Ltr_SC);  
     Dsc.Assem.Err.Rt2x{s}    = nan(1,P.Ltr_SC);  
     Dsc.Assem.Err.TS{s}      = nan(P.Ltr_SC,1);  
     Dsc.Assem.Err.dfnum{s}   = NaN;
     Dsc.Assem.Err.dfd{s}     = NaN;
    end
     % Peak decoding
     Dsc.Assem.Err.TSPeak{s}     = max(Dsc.Assem.Err.TS{s},[],1);
     Dsc.Assem.Err.TSIntegral{s} = sum(Dsc.Assem.Err.TS{s},1);
end


x_temp = (1:P.Ltr_SC)*P.bw-P.Ltr_SC/2*P.bw;

if P.flag.PlotVerbose
    figure('color','w','name','Sample/Choice decoding');
    subplot(1,2,1); hold on
    ciplot(Dsc.units.Err.Ft2ciH0{1},Dsc.units.Err.Ft2ciL0{1},x_temp,'r')
    ciplot(Dsc.units.Err.Ft2ciH0{2},Dsc.units.Err.Ft2ciL0{2},x_temp,'g')
    plot(x_temp, Dsc.units.Err.Ft2x{1},'r')
    plot(x_temp, Dsc.units.Err.Ft2x{2},'g')
    title('Unit Sample/Choice decoding')
    axis([-Inf Inf 0 30])
    
    subplot(1,2,2); hold on
    ciplot(Dsc.Assem.Err.Ft2ciH0{1},Dsc.Assem.Err.Ft2ciL0{1},x_temp,'r')
    ciplot(Dsc.Assem.Err.Ft2ciH0{2},Dsc.Assem.Err.Ft2ciL0{2},x_temp,'g')
    ciplot(Dsc.Assem.Err.Ft2ciH0{3},Dsc.Assem.Err.Ft2ciL0{3},x_temp,'b')
    plot(x_temp, Dsc.Assem.Err.Ft2x{1},'r')
    plot(x_temp, Dsc.Assem.Err.Ft2x{2},'g')
    plot(x_temp, Dsc.Assem.Err.Ft2x{3},'b')
    title('Assembly Sample/Choice decoding')
    axis([-Inf Inf 0 30])
end
clear s input iTrial cutout tempFSC  x_temp side
%% Firing rate decoder - Cross-validation
if P.flag.runSingleTrialDecode
Nbs = 10;
x_temp = (1:P.Ltr*2)*P.bw;

for s=1:3
    % 1 Make input
    if s ==3
        input = cell2mat([FRtrials.iFR0_{1}' , FRtrials.iFR0_{2}']);
    else
        input = cell2mat(FRtrials.iFR0_{s}');
    end
    
    if  P.nUnits(s)==min(P.nUnits)
        nrep=1;
    else
        nrep=10;
    end;
    rs=cell(1,nrep);
    pe0=zeros(nrep,P.Ltr*2);
    % 2 Calcuate decoder accuracy using CVE
    clear temp pClass1 pe0 pe0units pClass1units

    for i=1:nrep
        
        [pe0(i,:)]=DecodeCVE(input(:,rs{i}),FRtrials.evt0,P.reg); % run leave-one out CVE on this draw
        
        % pe0             Cross-validation error (1-accuracy) predicted from the population
        % pClass1         Response at each timestep predicted from the population
        % pe0units        Cross-validation error (1-accuracy) predicted from each unit
        % pClass1units    Response at each timestep predicted from each unit
        
        rs{i}=randsample(P.nUnits(s),min(P.nUnits)); % random draw of units
        % rs{i}=1:P.nUnits(s);
        [pe0(i,:),pClass1{i},pe0units{i},pClass1units{i}] = DecodeCVEdistance(input(:,rs{i}),FRtrials.evt0,P.reg); 
        
        for iUnit = 1:P.nUnits(s)
            figure; hold on
            y_temp1 = pClass1units{i}{iUnit}(:,find(FRtrials.evt0==1))
            plot(x_temp,y_temp1,'b')
            y_temp2 = pClass1units{i}{iUnit}(:,find(FRtrials.evt0==2));
            plot(x_temp,1-y_temp2,'r')
            for iTrial = 1:length(FRtrials.evt0)
                 temp(:,iTrial) = xcorr(pClass1units{i}{iUnit}(:,iTrial));
            end
                    figure; hold on
                    plot(mean(temp(:,FRtrials.evt0==1),2))
                    plot(mean(temp(:,FRtrials.evt0==2),2))
        end


        
         
    plot(x_temp,pClass1(:,3))

    CVE{s}(1,:)=mean(pe0,1); % Collapse results across all draws
    
    cveBSciH{s}=zeros(1,2*P.Ltr); cveBSciL{s}=zeros(1,2*P.Ltr);
    CVEbs{s}=[];
    for b=1:Nbs
        
        k=randperm(length(FRtrials.evt0));
        evt1=FRtrials.evt0(k);
        cve0=zeros(nrep,2*P.Ltr);
        for i=1:nrep
            fprintf('Bootstrap %d of %d, permutation %d of %d...\n',b,Nbs,i,nrep)
            cve0(i,:)=DecodeCVE(input(:,rs{i}),evt1,P.reg);
        end;
        CVEbs{s}(b,:)=mean(cve0,1);
    end;
    for t=1:2*P.Ltr
        cves=sort(CVEbs{s}(:,t),'ascend');
        cveBSciH{s}(t)=cves(round(0.95*Nbs));
        cveBSciL{s}(t)=cves(round(0.05*Nbs));
    end;
    
%     x_temp = (1:P.Ltr*2)*P.bw;
%     figure('name',P.names{s}); hold on
%     plot([10 10],[0 1],'color',[0.89 0.95 0.89],'LineWidth',20)
%     plot([20 20],[0 1],'color',[0.95 0.89 0.89],'LineWidth',20)
%     plot([0 P.Ltr*2*P.bw],[0.5 0.5],':b');
%     plot([15 15],[0 1],':k');
%   	plot(x_temp,1-CVE{s},'color','b','LineWidth',1.5);
%     h=area(x_temp',1-cveBSciL{s}');
%     set(h(1),'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6],'FaceAlpha',1);
% %     plot(x_temp,1-cveBSciL{s},'color','b','LineWidth',1);
% axis([0 Inf 0 1]); xlabel('Time (s)'); ylabel('Decoder Accuracy')
    end
end
if P.flag.PlotVerbose
figure;
for s=1:3
subplot(3,1,s); hold on
    title(P.names{s});
    plot([10 10],[0 1],'color',[0.89 0.95 0.89],'LineWidth',20)
    plot([20 20],[0 1],'color',[0.95 0.89 0.89],'LineWidth',20)
    plot([0 P.Ltr*2*P.bw],[0.5 0.5],':b');
    plot([15 15],[0 1],':k');
    plot(x_temp,1-CVE{s},'color','b','LineWidth',1.5);
    h=area(x_temp',1-cveBSciL{s}');
    set(h(1),'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6],'FaceAlpha',1);
    %   plot(x_temp,1-cveBSciL{s},'color','b','LineWidth',1);
    axis([0 Inf 0 1]); xlabel('Time (s)'); ylabel('Decoder Accuracy')
end
figure; hold on
    plot(x_temp,1-CVE{1},'color','r','LineWidth',1.5);
    plot(x_temp,1-CVE{2},'color','b','LineWidth',1.5);
    plot(x_temp,1-CVE{3},'color','g','LineWidth',1.5);
    plot([10 10],[0 1],'color',[0.89 0.95 0.89],'LineWidth',2)
    plot([20 20],[0 1],'color',[0.95 0.89 0.89],'LineWidth',2)
    plot([0 P.Ltr*2*P.bw],[0.5 0.5],':b');
    plot([15 15],[0 1],':k');
    legend(P.names)
end
clear s input
end
%% Firing rate decoder - Cross-validation - single trial posterior probabilities
if P.flag.runSingleTrialDecode
Nbs = 10;
x_temp = (1:P.Ltr*2)*P.bw;

% Single trial Unit decoding
for s=1:2
    % 1 Make input
    if s ==3
        input = cell2mat([FRtrials.iFR0_{1}' , FRtrials.iFR0_{2}']);
    else
        input = cell2mat(FRtrials.iFR0_{s}');
    end
    
    if  P.nUnits(s)==min(P.nUnits)
        nrep=1;
    else
        nrep=10;
    end;
    rs=cell(1,nrep);
    % 2 Calcuate decoder accuracy using CVE
    clear temp pClass1 pe0 pe0units pClass1units

    for i=1:nrep
        fprintf('...Draw %d of %d\n',i,nrep)
        rs{i}=randsample(P.nUnits(s),min(P.nUnits)); % random draw of units
        [Accuracy{s}(i,:),pCorrect{s}{i}] = GetBinaryPosteriorProbs(input(:,rs{i}),FRtrials.evt0); 
    end
    if  P.nUnits(s)~=min(P.nUnits)
        Accuracy{s} =mean(Accuracy{s});
        pCorrect{s}{1} =mean(cell2mat(cellfun(@nanmean,cellfun(@transpose,pCorrect{s},'UniformOutput',false),'UniformOutput',false)'));
    end
    [~,~,Accuracyunits{s},pCorrectunits{s}] = GetBinaryPosteriorProbs(input,FRtrials.evt0); 

         % Plot each single-trial decoder for each unit
         figure;
         for iUnit = 1:P.nUnits(s)
             subplot(ceil(P.nUnits(s)^0.5),ceil(P.nUnits(s)^0.5),iUnit);hold on
             y_temp1 = pCorrectunits{s}{iUnit}(:,:);    
             plot(x_temp,0.5*ones(length(x_temp),1),':k')
             plot(x_temp,y_temp1(:,1),'b');
             axis([0,2*P.Ltr*P.bw,0,1])
             axis off
         end
end

figure; hold on
    col_={'r','b'};
    for s=1:2
        area(x_temp,D.SingleTrial.Units.Accuracy{s},'edgecolor',col_{s},'facecolor',col_{s},'facealpha',0.8)
    end
    plot(x_temp,0.5*ones(length(x_temp),1),':k')

if P.flag.PlotVerbose
         % Plot each trial for each unit - Autocorrelation
         figure;
         for iUnit = 1:P.nUnits(s)
             subplot(ceil(P.nUnits(s)^0.5),ceil(P.nUnits(s)^0.5),iUnit);hold on
             for iTrial = 1:length(FRtrials.evt0)
                  [temp(:,iTrial),lags] = xcorr(pCorrectunits{i}{iUnit}(:,iTrial),'coeff');
             end
             plot(lags*P.bw,temp)
             % plot(mean(temp,2))
         end
         
         % Grid of Unit-Unit x-correlation
         temp=[];
         counter=1;
         for j = 1:P.nUnits(s)
             for k = 1:P.nUnits(s)
                 Xj = pCorrectunits{i}{j}(:,:);
                 Xk = pCorrectunits{i}{k}(:,:);
                 for iTrial=1:length(FRtrials.evt0)
                     [temp{j,k}(iTrial,:),lags] = xcorr(Xj(:,iTrial),Xk(:,iTrial));
                 end
                 lookup_(j,k)=counter;
                 counter=counter+1;
             end
         end
         
         figure('color','w'); hold on
         for l = 1:P.nUnits(s)*2
             if ismember(l,find(tril(lookup_)))
                 subplot(P.nUnits(s),P.nUnits(s),l);hold on
                 y_temp1 =temp{l}(:,:);
                 plot(lags*P.bw,mean(y_temp1),'b');
                 axis off
             end
         end
         
         
         
         % Plot quick and dirty FFT - each unit
         NFFT=1028;	 	 
         L =P.Ltr*2;
         fVals=(1/P.bw)*(0:NFFT/2-1)/NFFT;
         Px=[];
         figure;
         for iUnit = 1:P.nUnits(s)
             subplot(ceil(P.nUnits(s)^0.5),ceil(P.nUnits(s)^0.5),iUnit);hold on
             for iTrial = 1:length(FRtrials.evt0)
                x= pCorrectunits{i}{iUnit}(:,iTrial); 
                X=fft(x,NFFT);	 	 
                %Px(iTrial,:)=X.*conj(X)/(NFFT*L); 
                Px(iTrial,:)=10*log10(abs(X));
             end
                plot(fVals,mean(Px(:,1:NFFT/2)),'b');	 	 
%                 xlim([0 10])
                set(gca,'XScale','log')
                % axis off
         end
end  % Some other stuff
D.SingleTrial.Units.Accuracy = Accuracy;
D.SingleTrial.Units.pCorrect = pCorrect;
D.SingleTrial.Units.Accuracyunits = Accuracyunits;
D.SingleTrial.Units.pCorrectunits = pCorrectunits;

% Single trial Assembly decoding
nAssems = cellfun(@length,FAcont.Task.keepers);
for s=1:2
    % 1 Make input
   
        input = cell2mat(FAtrials.FSC_{s}');
    
    if  nAssems(s)==min(nAssems)
        nrep=1;
    else
        nrep=10;
    end;
   
    % 2 Calcuate decoder accuracy using CVE
    clear temp pClass1 pe0 pe0units pClass1units

    for i=1:nrep
        fprintf('...Draw %d of %d\n',i,nrep)
        rs{i}=randsample(nAssems(s),min(nAssems)); % random draw of units
        [Accuracy{s}(i,:),pCorrect{s}{i}] = GetBinaryPosteriorProbs(input(:,rs{i}),FRtrials.evt0); 
    end
    if  nAssems(s)~=min(nAssems)
        Accuracy{s} =mean(Accuracy{s});
        pCorrect{s}{1} =mean(cell2mat(cellfun(@nanmean,cellfun(@transpose,pCorrect{s},'UniformOutput',false),'UniformOutput',false)'));
    end
    [~,~,AccuracyAssems{s},pCorrectAssems{s}] = GetBinaryPosteriorProbs(input,FRtrials.evt0); 

         % Plot each single-trial decoder for each unit
         figure;
         for iAssem = 1:nAssems(s)
             subplot(ceil(nAssems(s)^0.5),ceil(nAssems(s)^0.5),iAssem);hold on
             y_temp1 = pCorrectAssems{s}{iAssem}(:,:);    
             plot(x_temp,0.5*ones(length(x_temp),1),':k')
             plot(x_temp,y_temp1(:,1),'b');
             axis([0,2*P.Ltr*P.bw,0,1])
             axis off
         end
end
     
figure; hold on
    col_={'r','b'};
    for s=1:2
        area(x_temp,D.SingleTrial.Assems.Accuracy{s},'edgecolor',col_{s},'facecolor',col_{s},'facealpha',0.8)
    end
    plot(x_temp,0.5*ones(length(x_temp),1),':k')

D.SingleTrial.Assems.Accuracy = Accuracy;
D.SingleTrial.Assems.pCorrect = pCorrect;
D.SingleTrial.Assems.AccuracyAssems = AccuracyAssems;
D.SingleTrial.Assems.pCorrectAssems = pCorrectAssems;



for s =1
    for iAss = 1:nAssems(s)
        figure
        for iTrial=1:P.ntr
        temp = [];
        subplot(ceil(P.ntr^0.5),ceil(P.ntr^0.5),iTrial);hold on
        memberUnits_ = FAcont.Task.units{s}{iAss};
        for iUnit=1:length(memberUnits_)
            temp(:,iUnit)= D.SingleTrial.Units.pCorrectunits{s}{memberUnits_(iUnit)}(:,iTrial);
        end
        plot(x_temp,temp,'r','LineWidth',1)
        plot(x_temp,D.SingleTrial.Assems.pCorrectAssems{s}{iAss}(:,iTrial),'k','LineWidth',2)
        axis([0,2*P.Ltr*P.bw,0,1])
        plot(x_temp,0.5*ones(length(x_temp),1),':k')
        axis off
        end
    end
end




for s =1        
    figure
    for iAss = 1:nAssems(s)
        temp = [];
        subplot(ceil(nAssems(s)^0.5),ceil(nAssems(s)^0.5),iAss);hold on
        memberUnits_ = FAcont.Task.units{s}{iAss};
        for iUnit=1:length(memberUnits_)
            temp(:,iUnit)= mean(D.SingleTrial.Units.pCorrectunits{s}{memberUnits_(iUnit)},2);
        end
        plot(x_temp,temp,'r','LineWidth',1)
        plot(x_temp,mean(D.SingleTrial.Assems.pCorrectAssems{s}{iAss},2),'k','LineWidth',2)
        axis([0,2*P.Ltr*P.bw,0,1])
        plot(x_temp,0.5*ones(length(x_temp),1),':k')
%         axis off
    end
end

clear s input
end
%% LR Decoding on Assembly members vs non-members
P.reg = 0.05;
%   Firing rate decoder 
for s=1:3
    if s==3
        temp = cell2mat([FRtrials.iFR0_{1}' , FRtrials.iFR0_{2}']);
    else
        temp = cell2mat(FRtrials.iFR0_{s}');
    end
    memberUnits_=[];
    if numel(FAtrials.units{s}) ~= 0
        % Assembly member units (separated by assembly)
        for iAss = 1:numel(FAtrials.units{s}) 

            memberUnits_ = FAtrials.units{s}{iAss};
            [~,~,...
             D.MemberUnits.Ft2x{s}{iAss},...
             D.MemberUnits.Rt2x{s}{iAss},...
             D.MemberUnits.Ft2ciL0{s}{iAss},...
             D.MemberUnits.Ft2ciH0{s}{iAss},...
             D.MemberUnits.TS{s}{iAss},...
             D.MemberUnits.dfnum{s}{iAss},...
             D.MemberUnits.dfd{s}{iAss}] = DecodeStats(temp(:,memberUnits_),FRtrials.evt0,P.reg);
        end

        % Assembly member units (all assemblies collapsed)
        memberUnits_ = unique(cell2mat(FAtrials.units{s}));

        [~,~,...
         D.MemberUnitsCollapsed.Ft2x{s},...
         D.MemberUnitsCollapsed.Rt2x{s},...
         D.MemberUnitsCollapsed.Ft2ciL0{s},...
         D.MemberUnitsCollapsed.Ft2ciH0{s},...
         D.MemberUnitsCollapsed.TS{s},...
         D.MemberUnitsCollapsed.dfnum{s},...
         D.MemberUnitsCollapsed.dfd{s}] = DecodeStats(temp(:,memberUnits_),FRtrials.evt0,P.reg);
    else
    end
    % Assembly non-member units
    temp(:,memberUnits_)=[];
    [~,~,...
     D.NonMemberUnits.Ft2x{s},...
     D.NonMemberUnits.Rt2x{s},...
     D.NonMemberUnits.Ft2ciL0{s},...
     D.NonMemberUnits.Ft2ciH0{s},...
     D.NonMemberUnits.TS{s},...
     D.NonMemberUnits.dfnum{s},...
     D.NonMemberUnits.dfd{s}] = DecodeStats(temp,FRtrials.evt0,P.reg);
end
%% SC Decoding on Assembly members vs non-members
% tempFR rolls over from previous block, just need to prune units...
%   Firing rate decoder 
for s=1:3
    temp = cell2mat(tempFR{s}');
    memberUnits_=[];
    if numel(FAtrials.units{s}) ~= 0
        % Assembly member units (separated by assembly)
        for iAss = 1:numel(FAtrials.units{s}) 
            
            
            memberUnits_    = FAtrials.units{s}{iAss};
            [~,~,...
             Dsc.MemberUnits.Ft2x{s}{iAss},...
             Dsc.MemberUnits.Rt2x{s}{iAss},...
             Dsc.MemberUnits.Ft2ciL0{s}{iAss},...
             Dsc.MemberUnits.Ft2ciH0{s}{iAss},...
             Dsc.MemberUnits.TS{s}{iAss},...
             Dsc.MemberUnits.dfnum{s}{iAss},...
             Dsc.MemberUnits.dfd{s}{iAss}] = DecodeStats(temp(:,memberUnits_),outcome,P.reg);
        end

        % Assembly member units (all assemblies collapsed)
        memberUnits_ = unique(cell2mat(FAtrials.units{s}));
        

        [~,~,...
         Dsc.MemberUnitsCollapsed.Ft2x{s},...
         Dsc.MemberUnitsCollapsed.Rt2x{s},...
         Dsc.MemberUnitsCollapsed.Ft2ciL0{s},...
         Dsc.MemberUnitsCollapsed.Ft2ciH0{s},...
         Dsc.MemberUnitsCollapsed.TS{s},...
         Dsc.MemberUnitsCollapsed.dfnum{s},...
         Dsc.MemberUnitsCollapsed.dfd{s}] = DecodeStats(temp(:,memberUnits_),outcome,P.reg);
    else
    end
    % Assembly non-member units
    %nonmemberUnits_ = setdiff(1:size(FRtrials.iFR0_{s}{1},2), memberUnits_);
    temp(:,memberUnits_)=[];
    [~,~,...
     Dsc.NonMemberUnits.Ft2x{s},...
     Dsc.NonMemberUnits.Rt2x{s},...
     Dsc.NonMemberUnits.Ft2ciL0{s},...
     Dsc.NonMemberUnits.Ft2ciH0{s},...
     Dsc.NonMemberUnits.TS{s},...
     Dsc.NonMemberUnits.dfnum{s},...
     Dsc.NonMemberUnits.dfd{s}] = DecodeStats(temp,outcome,P.reg);
end
%% Plotting optional
if P.flag.PlotVerbose
%% Plot LR decoding - units 
for s=1:2
%     stepSize=0;
%     figure('Name',[P.names{s},' Units (individual)'],'color','w'); hold on
%     for iUnit=1:P.nUnits(s)
%         subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
%         axis off; hold on       
%         %         ciplot(stepSize*iUnit+(FRtrials.iFR0_Mean{s}{1}(:,iUnit)+FRtrials.iFR0_SEM{s}{1}(:,iUnit)),...
%         %               stepSize*iUnit+(FRtrials.iFR0_Mean{s}{1}(:,iUnit)-FRtrials.iFR0_SEM{s}{1}(:,iUnit)),...
%         %               (1:P.Ltr*2)*P.bw,'r')
%         % % 
%         %        ciplot(stepSize*iUnit+(FRtrials.iFR0_Mean{s}{2}(:,iUnit)+FRtrials.iFR0_SEM{s}{2}(:,iUnit)),...
%         %               stepSize*iUnit+(FRtrials.iFR0_Mean{s}{2}(:,iUnit)-FRtrials.iFR0_SEM{s}{2}(:,iUnit)),...
%         %               (1:P.Ltr*2)*P.bw,'g')
%         %        plot((1:P.Ltr*2)*P.bw,FRtrials.iFR0_Mean{s}{1}(:,iUnit),'r','LineWidth',1.5)
%         %        plot((1:P.Ltr*2)*P.bw,FRtrials.iFR0_Mean{s}{2}(:,iUnit),'g','LineWidth',1.5)
%        plot((1:P.Ltr*2)*P.bw,stepSize*iUnit+(D.units.TS{s}(:,iUnit)),'color',[0.6 0.6 0.6],'LineWidth',1.5)
%    end
    
    figure('Name',[P.names{s}, ' Units (group)'],'color','w'); hold on
    subplot(2,1,1)
        Y1=[D.units.Ft2ciL0{s};D.units.Ft2ciH0{s}-D.units.Ft2ciL0{s}];
        h=area((1:P.Ltr*2)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','k','FaceAlpha',0.6); hold on
        plot((1:P.Ltr*2)*P.bw,D.units.Ft2x{s},'k','LineWidth',1.5)
	    axis([0 P.Ltr*P.bw*2 0 10])
        title('Unit group decoding: Correct trials')
    subplot(2,1,2)
        Y1=[D.units.Err.Ft2ciL0{s};D.units.Err.Ft2ciH0{s}-D.units.Err.Ft2ciL0{s}];
        h=area((1:P.Ltr*2)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','r','FaceAlpha',0.6); hold on
        plot((1:P.Ltr*2)*P.bw,D.units.Err.Ft2x{s},'r','LineWidth',1.5)
        axis([0 P.Ltr*P.bw*2 0 10])
        title('Unit group decoding: Error trials')
end
%% Plot LR decoding - factors
for s=1:3
    figure('Name',[P.names{s},' Assemblies (group)'],'color','w'); hold on
    for iFactor=1:P.nFactors(s)
            % subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
            axis off
            % {area}{outcome}[time,Unit]
%             ciplot(stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor)+FAtrials.FSC_SEM{s}{1}(:,iFactor),...
%                    stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
%                    (1:P.Ltr*2)*P.bw,'r')
% 
%            ciplot(stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor)+FAtrials.FSC_SEM{s}{2}(:,iFactor),...
%                    stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
%                    (1:P.Ltr*2)*P.bw,'g')

%              plot((1:P.Ltr*2)*P.bw,stepSize*iUnit+(D.Assem.TS{s}(:,iFactor)),'k','LineWidth',1.5)


    end
    subplot(2,1,1)
        Y1=[D.Assem.Ft2ciL0{s};D.Assem.Ft2ciH0{s}-D.Assem.Ft2ciL0{s}];
        h=area((1:P.Ltr*2)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','k','FaceAlpha',0.6); hold on
        plot((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s},'k','LineWidth',1.5)
	    axis([0 P.Ltr*P.bw*2 0 10])
        title('Assembly group decoding: Correct trials')
    subplot(2,1,2)
        Y1=[D.Assem.Err.Ft2ciL0{s};D.Assem.Err.Ft2ciH0{s}-D.Assem.Err.Ft2ciL0{s}];
        h=area((1:P.Ltr*2)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','r','FaceAlpha',0.6); hold on
        plot((1:P.Ltr*2)*P.bw,D.Assem.Err.Ft2x{s},'r','LineWidth',1.5)
        axis([0 P.Ltr*P.bw*2 0 10])
        title('Assembly group decoding: Error trials')
end
%% Plot SC decoding - units 
for s=1:2
    figure('Name',[P.names{s}, ' Unit (group)'],'color','w'); hold on
    subplot(2,1,1)
        Y1=[Dsc.units.Ft2ciL0{s};Dsc.units.Ft2ciH0{s}-Dsc.units.Ft2ciL0{s}];
        h=area((1:P.Ltr_SC)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','k','FaceAlpha',0.6); hold on
        plot((1:P.Ltr_SC)*P.bw,Dsc.units.Ft2x{s},'k','LineWidth',1.5)
	    axis([0 P.Ltr_SC*P.bw 0 10])
        title('Unit group decoding: Correct trials')
    subplot(2,1,2)
        Y1=[Dsc.units.Err.Ft2ciL0{s};Dsc.units.Err.Ft2ciH0{s}-Dsc.units.Err.Ft2ciL0{s}];
        h=area((1:P.Ltr_SC)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','r','FaceAlpha',0.6); hold on
        plot((1:P.Ltr_SC)*P.bw,Dsc.units.Err.Ft2x{s},'r','LineWidth',1.5)
	    axis([0 P.Ltr_SC*P.bw 0 10])
        title('Unit group decoding: Error trials')
end
%% Plot SC decoding - Assemblies
for s=1:3
    figure('Name',[P.names{s}, ' Assembly (group)'],'color','w'); hold on
    subplot(2,1,1)
        Y1=[Dsc.Assem.Ft2ciL0{s};Dsc.Assem.Ft2ciH0{s}-Dsc.Assem.Ft2ciL0{s}];
        h=area((1:P.Ltr_SC)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','k','FaceAlpha',0.6); hold on
        plot((1:P.Ltr_SC)*P.bw,Dsc.Assem.Ft2x{s},'k','LineWidth',1.5)
	    axis([0 P.Ltr_SC*P.bw 0 10])
        title('Assembly group decoding: Correct trials')
    subplot(2,1,2)
        Y1=[Dsc.Assem.Err.Ft2ciL0{s};Dsc.Assem.Err.Ft2ciH0{s}-Dsc.Assem.Err.Ft2ciL0{s}];
        h=area((1:P.Ltr_SC)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','r','FaceAlpha',0.6); hold on
        plot((1:P.Ltr_SC)*P.bw,Dsc.Assem.Err.Ft2x{s},'r','LineWidth',1.5)
	    axis([0 P.Ltr_SC*P.bw 0 10])
        title('Assembly group decoding: Error trials')
end
end
%% Plot decoding Assems vs single units
if P.flag.PlotVerbose
figure('name',P.rat,'color','w')
% Plot decoding
maxH = 10;

for s=1:3
    subplot(2,3,s); hold on
    title(P.names{s})

%     plot([5 5],[0 10],'color',[0.89 0.95 0.89],'LineWidth',20)
%     plot([25 25],[0 10],'color',[0.95 0.89 0.89],'LineWidth',20) 
    plot([10 10],[0 10*maxH],'color',[0.89 0.95 0.89],'LineWidth',20)
    plot([20 20],[0 10*maxH],'color',[0.95 0.89 0.89],'LineWidth',20)
    
%     Y1=[D.Assem.Ft2ciL0{s};D.Assem.Ft2ciH0{s}-D.Assem.Ft2ciL0{s}];
%     h=area((1:P.Ltr*2)*P.bw,Y1');
%     set(h(1),'FaceColor','w'); set(h(2),'FaceColor','w'); hold on
    plot((1:P.Ltr*2)*P.bw,D.Assem.Ft2ciH0{s},'k')
    plot((1:P.Ltr*2)*P.bw,D.units.Ft2ciH0{s},'b')

	plot((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s},'k','LineWidth',1.5)
    plot((1:P.Ltr*2)*P.bw,D.units.Ft2x{s},'b','LineWidth',1.5)
%   area((1:P.Ltr*2)*P.bw,D.units.Ft2x{s},'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
% 	area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s},'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)

% area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s}-D.units.Ft2x{s},'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    plot([15 15],[0 10*maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 10*maxH],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    set(gca,'XTick',[])
    
    if s==1
        ylabel({'Decoding score';'(F-value)'})
        axis([0 Inf 0 maxH])

    elseif s==2
%         xlabel('Time (s)')
        text(15,-1,'Single unit decoding','Color',[0 0 1],'HorizontalAlignment','center' )
        text(15,-3,'Assembly decoding','Color',[0 0 0],'HorizontalAlignment','center' )
        axis([0 Inf 0 2*maxH])
    elseif s==3      
        axis([0 Inf 0 2*maxH])
    end
end
% Plot assembly advantage
for s=1:3
    subplot(2,3,s+3); hold on
%     title(P.names{s})

    plot([9 9],[-10 10],'color',[0.89 0.95 0.89],'LineWidth',20)
    plot([20 20],[-10 10],'color',[0.95 0.89 0.89],'LineWidth',20)
    
    temp = D.Assem.Ft2x{s}-D.units.Ft2x{s}; temp(temp>0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    temp = D.Assem.Ft2x{s}-D.units.Ft2x{s}; temp(temp<0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'FaceColor',[0.6 0.6 0.6],'EdgeColor','k','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    % area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s}-D.units.Ft2x{s},'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
%     plot([0 P.Ltr*2*P.bw],[0 0],'color','k','LineWidth',2)
    plot([15 15],[-10 10],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[-10 10],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    axis([0 Inf -4 4])
    if s==1
%         ylabel({'Assembly decodimg advantage';'Assembly - Unit'})
        ylabel({'Assembly decodimg advantage'})
    elseif s==2
        xlabel('Time (s)')
    end
end
 fnam=[pat{1} P.rat '_' P.target 'AssVsUnits']; % your file name
     snam='8x5';
     s=hgexport('readstyle',snam);
	 s.Format = 'jpeg'; 
     hgexport(gcf,[fnam, '.jpeg'],s);
     sdf('8x5')
     savefig(gcf,[fnam, '.fig'])
     
end
%% Process activation patterns in reconstructed continuous assemblies

% corrected factor model loadings (rejecting assemblies with only one unit)
for s= 1:3
    FAcont.Task.FL_{s} = FAcont.Task.FL{s}{FAcont.Task.nassem{s}(3)};
    FAcont.Task.FL_{s} = FAcont.Task.FL_{s}(:,FAcont.Task.keepers{s});
end

FAcont.Pre.patterns   = getAssemPatterns('Pre',P.bw,FAcont.Pre.FSC,...
                                         FAcont.Task.ciHld,FAcont.Task.ciHsc,FAcont.Task.FL_);
                                     
FAcont.Task.patterns  = getAssemPatterns('Task',P.bw,FAcont.Task.FSC,...
                                         FAcont.Task.ciHld,FAcont.Task.ciHsc,FAcont.Task.FL_);
     if P.flag.noPostSleep
         FAcont.Post.patterns   = struct;
     else
        FAcont.Post.patterns   = getAssemPatterns('Post',P.bw,FAcont.Post.FSC,...
                                         FAcont.Task.ciHld,FAcont.Task.ciHsc,FAcont.Task.FL_);                               
     end
%% Process activation patterns in delay period only
FAtrials.DelayOnly.patterns  = getCutoutAssemPatterns('Delay Period',P.bw,FAtrials.DelayOnly.FSC_,...
                                                FAcont.Task.ciHld,FAcont.Task.ciHsc,FAcont.Task.FL_,FRtrials.evt0,P.flag.PlotVerbose); %
%% Plot patterns before and after task, and during sleep
if P.flag.PlotVerbose
    plotAssemPatterns(FAcont.Pre.patterns)
    if ~P.flag.noPostSleep
        plotAssemPatterns(FAcont.Post.patterns)
    end
end
%% Example plot...

% figure; hold on
% s= 3
% iAss = 1
% norm_ = numel(FAcont.Task.patterns.pks.repHistBins);
% plot(1./FAcont.Task.patterns.pks.repHistBins,...
%      smooth_hist(FAcont.Task.patterns.pks.repHist{s}(iAss,:)),'LineWidth',1.5,'Color','g')
% plot(1./FAcont.Pre.patterns.pks.repHistBins,...
%      smooth_hist(FAcont.Pre.patterns.pks.repHist{s}(iAss,:)),'LineWidth',1.5,'Color','r')
% plot(1./FAcont.Post.patterns.pks.repHistBins,...
%      smooth_hist(FAcont.Post.patterns.pks.repHist{s}(iAss,:)),'LineWidth',1.5,'Color','b')
%  set(gca,'XScale','log');
% ylabel('Fraction of repeats')
% xlabel('Assembly Activation rate (Hz)')
%% Compare assembly activity between epochs
% P1vsP2 = compareAssemPatterns(P1,P2,plotYN,PermittedTransitions,RankOrder)
FAcont.PreVsTask  = compareAssemPatterns(FAcont.Pre.patterns, FAcont.Task.patterns,true,3,true);
if ~P.flag.noPostSleep
    FAcont.PreVsPost  = compareAssemPatterns(FAcont.Pre.patterns, FAcont.Post.patterns,true,3,true);
    FAcont.TaskVsPost = compareAssemPatterns(FAcont.Task.patterns,FAcont.Post.patterns,true,3,true);
else
    FAcont.PreVsPost  = struct;
    FAcont.TaskVsPost = struct;
end
FAcont.ComparisonNames={'Pre vs. Post','Pre vs. Task','Task vs. Post'};
if P.flag.PlotVerbose
    colors = [0.2 0.5 0.9 ; ...
          0.9 0.2 0.5 ; ...
          0.9 0.9 0.5] ; ...
    minmax = 0.2;
    bins = -minmax:0.01:minmax; 
    figure('color', 'w'); 
        subplot('Position',[0.1 0.15 0.55 0.8]);hold on
        plot(sort(FAcont.PreVsPost.ptn.delta_transMatrix(1:end)),'LineWidth',1.5,'color',colors(1,:))
        plot(sort(FAcont.PreVsTask.ptn.delta_transMatrix(1:end)),'LineWidth',1.5,'color',colors(2,:))
        plot(sort(FAcont.TaskVsPost.ptn.delta_transMatrix(1:end)),'LineWidth',1.5,'color',colors(3,:))
        plot([1 sum(~isnan(FAcont.PreVsTask.ptn.delta_transMatrix(1:end)))],[0 0],':k')

        axis([1 Inf -minmax minmax])
        ylabel('Change in p(Pattern)')
        xlabel('Rank ordered assembly pattern no.')
        subplot('Position',[0.7 0.15 0.25 0.8]);hold on
        x{1} = histc(sort(FAcont.PreVsPost.ptn.delta_transMatrix(1:end)),bins);
        x{2} = histc(sort(FAcont.PreVsTask.ptn.delta_transMatrix(1:end)),bins);
        x{3} = histc(sort(FAcont.TaskVsPost.ptn.delta_transMatrix(1:end)),bins);
        plot(x{1}./max(x{1}),bins,'LineWidth',1.5,'color',colors(1,:))
        plot(x{2}./max(x{2}),bins,'LineWidth',1.5,'color',colors(2,:))
        plot(x{3}./max(x{3}),bins,'LineWidth',1.5,'color',colors(3,:))
        legend(FAcont.ComparisonNames)
        set(gca,'YTick',[])
        % set(gca,'XScale','log')
        xlabel('Normalized distibution')
        plot([0 1],[0 0],':k')
        axis([0 1 -.4 .4])
end
clear x bins
if ~P.flag.noPostSleep && P.flag.PlotVerbose
    figure('color','w'); hold on
    plot(1./FAcont.Task.patterns.pks.repHistBins,(smooth_hist(FAcont.Pre.patterns.pks.repHist{1}(1,:))),'color','r','LineWidth',1.5)
    plot(1./FAcont.Task.patterns.pks.repHistBins,(smooth_hist(FAcont.Task.patterns.pks.repHist{1}(1,:))),'color','g','LineWidth',1.5)
    plot(1./FAcont.Task.patterns.pks.repHistBins,(smooth_hist(FAcont.Post.patterns.pks.repHist{1}(1,:))),'color','b','LineWidth',1.5)
    set(gca,'XScale','log')
    legend(P.names); legend('boxoff')
    xlabel('Activation rate (Hz)')
    ylabel('Distribution')
end
%% Delay sequence probability vs sum decoding power
% Average p(sequence) over both left and right outcomes
Seq = (FAtrials.DelayOnly.patterns.ptn.transMatrix_sortedMean{4}{1} + FAtrials.DelayOnly.patterns.ptn.transMatrix_sortedMean{4}{2})/2;
SumDecode     = [];
SumDecodeSort = [];
SeqSort       = [];
SeqSortidx    = [];

% Decode = cell2mat(D.Assem.TSPeak);
Decode = cell2mat(D.Assem.TSIntegral);

Decode(isnan(Decode))=[]; % NB no assembly case has inserted NaN;
for iAss = 1:length(Decode)
    for iAss_ = 1:length(Decode)
        SumDecode(iAss,iAss_) = Decode(iAss)+Decode(iAss_);
    end
end
SumDecode(find(eye(length(Decode)))) = Decode;  % Prevent doubling of decoding capacity along identity diagonal
% Optional: Prune zero probabiity patterns
SumDecode(Seq == 0)= NaN;
Seq(Seq==0)=NaN;

% plot sorted by (change in) sequence probability
[SeqSort,SeqSortidx] = sort(Seq(1:end));
SumDecodeSort = SumDecode(SeqSortidx);

if P.flag.PlotVerbose
figure; hold on
    x= SeqSort;
    y = SumDecodeSort;
    scatter(x,y)
    f = ezfit(x,y,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.85*max(x), 0.8*max(y),strcat('R^2=', sprintf('%3.1g',f.r^2)),'color','r')
    xlabel('Sequence probability in delay period')
    ylabel('Sum Decoding capacity')
end
%% Pattern Stability vs decoding power
if P.flag.PlotVerbose
Seq = FAcont.PreVsPost.ptn.delta_transMatrix; % Change in sequence probability
% Seq = FAcont.Post.patterns.ptn.transMatrix{4}; % Pooled assemly sequence
SumDecode     = [];
SumDecodeSort = [];
SeqSort       = [];
SeqSortidx    = [];

Decode = cell2mat(D.Assem.TSPeak);
% Decode = cell2mat(D.Assem.TSIntegral);

Decode(isnan(Decode))=[]; % NB no assembly case has inserted NaN;
for iAss = 1:length(Decode)
    for iAss_ = 1:length(Decode)
        SumDecode(iAss,iAss_) = Decode(iAss)+Decode(iAss_);
    end
end
SumDecode(find(eye(length(Decode)))) = Decode;  % Prevent doubling of decoding capcity along identity diagonal
% Optional: Prune zero probabiity patterns

SumDecode(Seq == 0)= NaN;
Seq(Seq==0)=NaN;


% plot sorted by (change in) sequence probability
[SeqSort,SeqSortidx] = sort(Seq(1:end));
SumDecodeSort = SumDecode(SeqSortidx);

if P.flag.PlotVerbose
    figure; hold on
    scatter(1:length(Decode)^2,SeqSort,50,SumDecodeSort,'filled','LineWidth',1)
    plot([1 sum(~isnan(SeqSort))],[0 0],':k')
    xlabel('Pattern rank-ordered by probability')
    ylabel('\Delta p(Pattern)')
    colormap parula
    c = colorbar('Location','manual','Position',[0.75 0.2 0.05 0.2]);
    c.Label.String = {'Sum decoding capacity';'of assembly sequence'};

    axis tight
end
% plot sorted by total decoding power
% [DecodeSort,DecodeSortidx] = sort(SumDecode(1:end));
% SeqSort = Seq(DecodeSortidx);
% % figure; scatter(1:length(Decode)^2,DecodeSort,10+10*SeqSort,SeqSort,'LineWidth',1)
% figure; scatter(1:length(Decode)^2,DecodeSort,50,SeqSort,'filled','LineWidth',1)
% xlabel('Pattern rank-ordered by sum decoding power')
% ylabel('Sum decoding power')
% title('Colour: p(Sequence)')
% title('Colour: Sum decoding power of assemblies in sequence')
% colormap gray
% axis tight
end
%% Plot strength histograms for each epoch/assembly
if P.flag.PlotVerbose
    for s = 1:3
        figure
        nAss = size(FAcont.Pre.patterns.pks.repHist{s},1);
        for iAss = 1:nAss
                subplot(1,nAss,iAss);
                hold on
                plot(FAcont.Pre.patterns.pks.LogStrengthHist{s}(:,iAss),FAcont.Pre.patterns.pks.LogBins,'r')
                plot(FAcont.Task.patterns.pks.LogStrengthHist{s}(:,iAss),FAcont.Task.patterns.pks.LogBins,'g')
                plot(FAcont.Post.patterns.pks.LogStrengthHist{s}(:,iAss),FAcont.Post.patterns.pks.LogBins,'b')
    %             legend('Pre','Task','Post')
                scatter(FAcont.Pre.patterns.pks.LogStrengthQ3binValue{s}(iAss),FAcont.Pre.patterns.pks.LogStrengthQ3bin{s}(iAss),'r')
                scatter(FAcont.Task.patterns.pks.LogStrengthQ3binValue{s}(iAss),FAcont.Task.patterns.pks.LogStrengthQ3bin{s}(iAss),'g')
                scatter(FAcont.Post.patterns.pks.LogStrengthQ3binValue{s}(iAss),FAcont.Post.patterns.pks.LogStrengthQ3bin{s}(iAss),'b')

                scatter(FAcont.Pre.patterns.pks.ModeLogStrengthPos{s}(iAss),FAcont.Pre.patterns.pks.ModeLogStrength{s}(iAss),'r')
                scatter(FAcont.Task.patterns.pks.ModeLogStrengthPos{s}(iAss),FAcont.Task.patterns.pks.ModeLogStrength{s}(iAss),'g')
                scatter(FAcont.Post.patterns.pks.ModeLogStrengthPos{s}(iAss),FAcont.Post.patterns.pks.ModeLogStrength{s}(iAss),'b')
    %                 title([FAcont.Task.patterns.titles{s} ' assembly no. ' num2str(iAss)])
    %                 ylabel('Activation strength')
    %                 xlabel('Distribution')
        end
    end
end
%% Plot assembly activation in each epoch
if P.flag.PlotVerbose
    names_= {'Pre','Task','Post'};
    % cmap ={'r','g','b'}
    xmax=20;
    figure('name','Assembly activation strength')
    for s = 1:3
        subplot(3,1,s);hold on
        nAss = size(FAcont.Pre.patterns.pks.repHist{s},1);
        for iAss = 1:nAss
            plot([1 2 3],...
                [FAcont.Pre.patterns.pks.LogStrengthQ3bin{s}(iAss),...
                 FAcont.Task.patterns.pks.LogStrengthQ3bin{s}(iAss),...
                 FAcont.Post.patterns.pks.LogStrengthQ3bin{s}(iAss)],...
                '-ok')

    %         plot([1 2 3],...
    %             [FAcont.Pre.patterns.pks.ModeLogStrength{s}(iAss),...
    %              FAcont.Task.patterns.pks.ModeLogStrength{s}(iAss),...
    %              FAcont.Post.patterns.pks.ModeLogStrength{s}(iAss)],...
    %             '-ok')     
            if s ==3
                set(gca,'XTick',[1 2 3],'XTickLabels',names_)
            else
                set(gca,'XTick',[])
            end 
        end
        axis([0.8 3.2 0 xmax])
        title([P.names{s},' Assemblies'])
        ylabel('Activation strength')
    end
end
%% Plot change in assembly activation [rate,strength] vs peak decoding score in task
if P.flag.PlotVerbose
    xmax=20;
    figure('name','Change in assembly activation rate vs. decoding score');
    for s = 1:3
        subplot(1,3,s)
        if ~isempty(FAtrials.keepers{s})
            title(FAcont.Pre.patterns.titles{s})
            scatter(D.Assem.TSPeak{s}, FAcont.PreVsPost.pks.delta_repMode{s})
            %    axis([-5 5 0 1])
            title(P.names{s})
        else
            title(['No ',FAcont.Pre.patterns.titles{s},' assemblies '])
        end
        if s == 1
            ylabel({'\Delta Assembly activation rate';FAcont.ComparisonNames{1}})
        elseif s ==2
            xlabel('Peak decoding score during task')
        end
    end

    figure('name','Change in assembly activation strength vs. decoding score');
    for s = 1:3
        subplot(1,3,s)
        if ~isempty(FAtrials.keepers{s})
            title(FAcont.Pre.patterns.titles{s})
            scatter(D.Assem.TSPeak{s}, FAcont.PreVsPost.pks.delta_LogStrengthQ3bin{s})
            %    axis([-5 5 0 1])
            title(P.names{s})
        else
            title(['No ',FAcont.Pre.patterns.titles{s},' assemblies '])

        end

        if s == 1
            ylabel({'\Delta Assembly activation strength';FAcont.ComparisonNames{1}})
        elseif s ==2
            xlabel('Peak decoding score during task')
        end
    end
end
%% clean up and save
clear cutout Decode f f_ iAss_ iAss iFactor iOutcome iTrial iUnit lastfit maxH memberUnits_ s Seq SeqSort SeqSortidx SumDecode SumDecodeSort
clear temp tempCutout tempFA tempFA_ tempFA__ tempFR tempFR_ tempFR__ tempSize tempSize_ trials x y 
save([pat{1} P.rat '_' P.target '_DecodingVSsleep.mat'])
set(0,'DefaultFigureWindowStyle','normal')
disp(sprintf('>>Done in %4.1f seconds.', toc))
% Takes factor model computed for continuous data for the task epoch of experiment, then
% cast model onto cutout segments for decoding analysis
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
%% Preample, load 
clear all 
% close all
P.target = 'Task';
P.rat    = 'JaroslawLONG1';  %    Post sleep import error
% P.rat    = 'JaroslawLONG2';  %*
% P.rat    = 'KrzesimirLONG1'; %*
% P.rat    = 'KrzesimirLONG2'; %*
% P.rat    = 'KrzysztofLONG1'; %*
% P.rat    = 'KrzysztofLONG2'; %*
% P.rat    = 'MiroslawLONG0';  %*
% P.rat    = 'MiroslawLONG2';  %*   No Post sleep recorded
% P.rat    = 'NorbertLONG1';   %*
% P.rat    = 'NorbertLONG2';   %
% P.rat    = 'OnufryLONG1';    %*
% P.rat    = 'OnufryLONG2';    %      Post sleep error
P.expt    = 'LONG';
if ispc
   pat{1} = 'C:\Analysis\AssemblyAnalysis\raw';           % location of raw spike times
else
   pat{1} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Task';
end
pat{2} = [pat{1} filesep 'KDE_binsTaskonly'];                   % location of calculated firing P.rates for Task period
pat{3} = [pat{2} filesep P.expt 'Taskonly'];                       % location of the processed sleep assembly activations
cd(pat{1})

P.names={'PFC','HP','Joint HP~PFC'};
P.twin     = 10;   % Time window
P.minFR    = 0.1;  % minimal acceptable firing P.rate
P.critCvBW = 1e6;  % critical max bound on kernel width drift over time (variance)
P.bw       = 0.05; % Binwidth to operate on (time resolution)
P.flag.CutoutLeverPresses = true;
P.flag.PlotVerbose = true;
%% Specify time ranges 

P.tRangesNames= {'Pre-sample'; 'Sample Press'; 'Delay'; 'Choice Press'; 'Post-sample'};
   

% % P.tRanges = [2 8;...
% %              9 11;...
% %              13 17;...
% %              19 21;...
% %              24 30];
            
% P.tRanges = [0 8.99;...
%              9 11;...
%              11.01 18.99;...
%              19 21;...
%              21.01 30];
%             
% P.tRangesColors= [0.9 0.9 0.9;...
%                   0.7 0.9 0.7;...
%                 0.9 0.9 0.9;...
%                 0.9 0.7 0.7;...
%                 0.9 0.9 0.9];            
% %             
            
P.tRanges = [9 11;...
             19 21];
P.tRangesColors= [0.7 0.9 0.7;...
                  0.9 0.7 0.7];        

P.tRangesAdjusted= P.tRanges/P.bw;
%% Load data from cutout task dataset 
fn = {[pat{2} filesep P.rat '_HP_iFR50.mat'] ; ...
      [pat{2} filesep P.rat '_PFC_iFR50.mat']; ...
      [pat{3} filesep P.rat '_iFR50_FSC.mat']};

FRtrials.HP  = load(fn{1}); % spike rate matrices
FRtrials.PFC = load(fn{2});
FAtrials = load(fn{3});     % 2- Factor analysis results
    
% Import trial cutouts

% Enforce that input unit IDs are consistent between two analyses
SelectedUnits = {find(FRtrials.PFC.CvBW<=P.critCvBW); find(FRtrials.HP.CvBW<=P.critCvBW)};

%  load FR cutouts (lever press +/- P.twin)
[FRtrials.TmtxS,...
 FRtrials.iFR0, ...
 FRtrials.EvtLs,...
 FRtrials.EvtTs,...
 FRtrials.unitIDs] = SelTrials(fn{2},P.twin,SelectedUnits,'iFR','corr');

%  load Assem cutouts (lever press +/- P.twin)
[FAtrials.TmtxS,...
 FAtrials.FSCs, ...
 FAtrials.EvtLs,...
 FAtrials.EvtTs]   = SelTrialsFSCs(fn{3},P.twin,'corr');

clear fn

if P.flag.PlotVerbose
    
    figure('name','Sanity: Assem cutouts vs. FR continuous'); hold on
    plot(FRtrials.HP.Tmtx(1:end-1),FRtrials.HP.iFR(:,1))
    plot(FAtrials.Tmtx,FAtrials.FSCsel{1})
    legend('FRTrials.iFR','FAtrials.FSCsel')
    
    figure('name','Sanity: Assem cutouts vs. FR cutout'); hold on
    plot(FRtrials.TmtxS{1}{2},FRtrials.iFR0{1}{2}(:,1))
    plot(FAtrials.TmtxS{1}{2},FAtrials.FSCs{1}{2}(:,1))
    legend('FR','Assem')
    
    figure('name','Sanity: All assem cutouts vs. FR cutout'); hold on
    temp = cell2mat(FRtrials.iFR0{1}');
    plot(cell2mat(FRtrials.TmtxS{1}),temp(:,1)); clear temp
    plot(cell2mat(FAtrials.TmtxS{1}'),cell2mat(FAtrials.FSCs{1}'))
    legend('FR','Assem')
    
    figure;
    plot(cellfun(@length,FRtrials.iFR0{1})-cellfun(@length,FAtrials.TmtxS{1}))
    xlabel('Trial no.')
    ylabel('FR Trial Length - Assem Trial Length')
    
    figure; hold on
    scatter(FRtrials.EvtTs,FRtrials.EvtLs)
    scatter(FAtrials.EvtTs,FAtrials.EvtLs+0.01)
    
end

% reconstruct parameters for joint assems
FAtrials.unitIDs{1} = SelectedUnits{1};
FAtrials.unitIDs{2} = SelectedUnits{2};
FAtrials.unitIDs{3} = [SelectedUnits{1},SelectedUnits{2}+max(SelectedUnits{1})];
FRtrials.iFR0{3}    = cellfun(@horzcat,FRtrials.iFR0{1},FRtrials.iFR0{2},'UniformOutput',false);

P.nUnits   = cellfun(@length,FAtrials.unitIDs);
P.nFactors = cellfun(@length,FAtrials.units);
% P.nFactors =  [FAtrials.nassem{1}(3), FAtrials.nassem{2}(3), FAtrials.nassem{3}(3)];
disp(['No. Units = ' num2str(P.nUnits)])
disp(['No. Factors = ' num2str(P.nFactors)])
% Excise an equal duration section for each trial
FRtrials.TmtxS_= cell(1,2); FRtrials.iFR0_= FRtrials.TmtxS_; FAtrials.FSC_ = FRtrials.TmtxS_;
cutout=round((P.twin+5)/P.bw); 
P.Ltr=cutout;
for s=1:3
    for iTrial=1:length(FRtrials.TmtxS{1})
           % disp('Excising delay period')
            FRtrials.TmtxS_{s}{iTrial} = FRtrials.TmtxS{1}{iTrial}([1:cutout end-cutout+1:end]);   %commented out to include delay period activity
            FRtrials.iFR0_{s}{iTrial}  = FRtrials.iFR0{s}{iTrial}([1:cutout end-cutout+1:end],:);
            FAtrials.TmtxS_{s}{iTrial} = FAtrials.TmtxS{1}{iTrial}([1:cutout end-cutout+1:end]);   %commented out to include delay period activity
            FAtrials.FSC_{s}{iTrial}   = FAtrials.FSCs{s}{iTrial}([1:cutout end-cutout+1:end],:);
     end;    
end;
% FRtrials.Tmtx=cell2mat(FRtrials.TmtxS_{1})';


if P.flag.PlotVerbose
    vstr= {'L Choice','R Choice','L Sample','R Sample'};
    figure; hold on
    for i=1:length(FRtrials.TmtxS{1})
        plot(FRtrials.TmtxS{1}{i},i+0*FRtrials.TmtxS{1}{i},'r')
        plot(FRtrials.TmtxS_{1}{i},0.02+i+0*FRtrials.TmtxS_{1}{i},'r','LineWidth',1.5)
    end
%     for i=1:length(FAtrials.TmtxS{1})
%         plot(FAtrials.TmtxS{1}{i},0.05+i+0*FAtrials.TmtxS{1}{i},'k')
%         plot(FAtrials.TmtxS_{1}{i},0.07+i+0*FAtrials.TmtxS_{1}{i},'k','LineWidth',1.5)
%     end
    for i=1:length(FRtrials.EvtTs)
        text(FRtrials.EvtTs(i),round(i/2),vstr{FRtrials.EvtLs(i)},'Rotation',90)
    end
    xlabel('Time (s)')
    ylabel('Trial')
    
    
    
    figure('name','Assembly vs. Unit timebase alignment sanity check...'); hold on
    plot(FRtrials.TmtxS_{1}{1},FRtrials.iFR0_{1}{1}(:,3))
    
    temp = cell2mat(FRtrials.iFR0_{1}');
    plot(cell2mat(FRtrials.TmtxS_{1}),temp(:,1))
    plot(cell2mat(FAtrials.TmtxS_{1}'),cell2mat(FAtrials.FSC_{1}'))
    
    
    for iTrial = 1:10%length(FAtrials.TmtxS{1})
        subplot(5,5,iTrial); hold on
        plot(FAtrials.TmtxS_{1}{iTrial},FAtrials.FSC_{1}{iTrial})
        plot(FAtrials.TmtxS{1}{iTrial},FAtrials.FSCs{1}{iTrial}+0.1)
    end    
    
    
    
end

clear s iTrial 
%% Trial-outcome sorted activations

% remap trial outcome labels
P.ntr = length(FRtrials.EvtTs)/2;
FRtrials.evt0 = FRtrials.EvtLs(((1:P.ntr)-1)*2+1)-2;

FRtrials.iFR0_Sorted   = cell(3,1);    
FAtrials.FSC_Sorted = cell(3,1);  
FRtrials.iFR0_Mean = cell(3,1);  
FRtrials.iFR0_SEM = cell(3,1);  
FRtrials.FSC_Mean = cell(3,1);  
FRtrials.FSC_SEM = cell(3,1);  

% length of shortest trial 
% tempSize = cell2mat(cellfun(@size,FRtrials.iFR0_{1},'UniformOutput',false)'); 
% tempSize_= min(tempSize(:,1));    
% P.Ltr_cont = tempSize;
% P.Ltr_shortest = tempSize_;
  

% Sort by L/R trial outcome, take mean and SEM
for iOutcome = 1:2
    trials = find(FRtrials.evt0 == iOutcome);
        for s=1:3

            % Unit firing rates
            temp = FRtrials.iFR0_{s}(trials);
            for iUnit=1:size(temp{1},2)

                for iTrial=1:length(trials)
                    FRtrials.iFR0_Sorted{s}{iOutcome}{iUnit}(:,iTrial)=temp{iTrial}(:,iUnit);
                end

                FRtrials.iFR0_Mean{s}{iOutcome}(:,iUnit) = mean (FRtrials.iFR0_Sorted{s}{iOutcome}{iUnit},2);
                FRtrials.iFR0_SEM{s}{iOutcome} (:,iUnit) = nansem(FRtrials.iFR0_Sorted{s}{iOutcome}{iUnit},2);

            end

            % Factor scores
            temp = FAtrials.FSC_{s}(trials);
            for iFactor=1:size(temp{1},2)
                for iTrial=1:length(trials)
                    FAtrials.FSC_Sorted{s}{iOutcome}{iFactor}(:,iTrial)=temp{iTrial}(:,iFactor);
                end
                    FAtrials.FSC_Mean{s}{iOutcome}(:,iFactor) = mean (FAtrials.FSC_Sorted{s}{iOutcome}{iFactor},2);
                    FAtrials.FSC_SEM{s}{iOutcome} (:,iFactor) = nansem(FAtrials.FSC_Sorted{s}{iOutcome}{iFactor},2);
            end
        end
    
        
end
clear iOutcome trials s iUnit iFactor iTrial temp tempSize tempSize_
%% Run decoding

P.reg = 0.05;
%   Firing rate decoder 
for s=1:2
        [~,~,...
         D.units.Ft2x{s},...
         D.units.Rt2x{s},...
         D.units.Ft2ciL0{s},...
         D.units.Ft2ciH0{s},...
         D.units.TS{s},...
         D.units.dfnum{s},...
         D.units.dfd{s}] = DecodeStats(cell2mat(FRtrials.iFR0_{s}'),FRtrials.evt0,P.reg);
end
s=3;
[~,~,...
 D.units.Ft2x{s},...
 D.units.Rt2x{s},...
 D.units.Ft2ciL0{s},...
 D.units.Ft2ciH0{s},...
 D.units.TS{s},...
 D.units.dfnum{s},...
 D.units.dfd{s}] = DecodeStats(cell2mat([FRtrials.iFR0_{1}' , FRtrials.iFR0_{2}']),FRtrials.evt0,P.reg);
%   Assembly decoder 
for s=1:3
        try 
            [~,~,...
         D.Assem.Ft2x{s},...
         D.Assem.Rt2x{s},...
         D.Assem.Ft2ciL0{s},...
         D.Assem.Ft2ciH0{s},...
         D.Assem.TS{s},...
         D.Assem.dfnum{s},...
         D.Assem.dfd{s}] = DecodeStats(cell2mat(FAtrials.FSC_{s}'),FRtrials.evt0,P.reg);
        catch
        end

end
clear s
%% Decoding on Assembly members vs non-members


P.reg = 0.05;
%   Firing rate decoder 
for s=1:3
    if s==3
        temp = cell2mat([FRtrials.iFR0_{1}' , FRtrials.iFR0_{2}']);
    else
        temp = cell2mat(FRtrials.iFR0_{s}');
    end
    
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
 
    % Assembly non-	 units
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
%% plot Units and factors' LR preferences
% plot unit L/R profiles - subplots
if P.flag.PlotVerbose
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

%% plot unit L/R profiles - overlaid
stepSize=2;
    for s=1:2
        figure('Name',P.names{s});
        for iUnit=1:P.nUnits(s)
%             subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
            axis off; hold on
            %{area}{outcome}[time,Unit]
            ciplot(stepSize*iUnit+FRtrials.iFR0_Mean{s}{1}(:,iUnit)+FRtrials.iFR0_SEM{s}{1}(:,iUnit),...
                   stepSize*iUnit+FRtrials.iFR0_Mean{s}{1}(:,iUnit)-FRtrials.iFR0_SEM{s}{1}(:,iUnit),...
                   (1:P.Ltr*2)*P.bw,'r')

            ciplot(stepSize*iUnit+FRtrials.iFR0_Mean{s}{2}(:,iUnit)+FRtrials.iFR0_SEM{s}{2}(:,iUnit),...
                   stepSize*iUnit+FRtrials.iFR0_Mean{s}{2}(:,iUnit)-FRtrials.iFR0_SEM{s}{2}(:,iUnit),...
                   (1:P.Ltr*2)*P.bw,'g')

            plot((1:P.Ltr*2)*P.bw,stepSize*iUnit+FRtrials.iFR0_Mean{s}{1}(:,iUnit),'r','LineWidth',1.5)
            plot((1:P.Ltr*2)*P.bw,stepSize*iUnit+FRtrials.iFR0_Mean{s}{2}(:,iUnit),'g','LineWidth',1.5)

        end
    end  
%% plot factor L/R profiles - subplots

    for s=1:3
        figure('Name',P.names{s});
        for iFactor=1:size(FAtrials.FSC_Mean{s}{1},2)% P.nFactors(s)
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

clear iFactor iUnit
%% plot factor L/R profiles - overlaid
stepSize=5;
    for s=1:3
        figure('Name',P.names{s}); hold on
        for iFactor=1:size(FAtrials.FSC_Mean{s}{1},2)% P.nFactors(s)
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
end
%% Plot decoding 
if P.flag.PlotVerbose
% plot units
stepSize=0;
for s=1:2   
    figure('Name',P.names{s}); hold on
    for iUnit=1:P.nUnits(s)

       axis off
       ciplot(stepSize*iUnit+(FRtrials.iFR0_Mean{s}{1}(:,iUnit)+FRtrials.iFR0_SEM{s}{1}(:,iUnit)),...
              stepSize*iUnit+(FRtrials.iFR0_Mean{s}{1}(:,iUnit)-FRtrials.iFR0_SEM{s}{1}(:,iUnit)),...
              (1:P.Ltr*2)*P.bw,'r')
% 
       ciplot(stepSize*iUnit+(FRtrials.iFR0_Mean{s}{2}(:,iUnit)+FRtrials.iFR0_SEM{s}{2}(:,iUnit)),...
              stepSize*iUnit+(FRtrials.iFR0_Mean{s}{2}(:,iUnit)-FRtrials.iFR0_SEM{s}{2}(:,iUnit)),...
              (1:P.Ltr*2)*P.bw,'g')

       plot((1:P.Ltr*2)*P.bw,stepSize*iUnit+(D.units.TS{s}(:,iUnit)),'k','LineWidth',1.5)

       
    end

    Y1=[D.units.Ft2ciL0{s};D.units.Ft2ciH0{s}-D.units.Ft2ciL0{s}];
    h=area((1:P.Ltr*2)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','y'); hold on
	plot((1:P.Ltr*2)*P.bw,D.units.Ft2x{s},'k','LineWidth',1.5)
end
% Plot  factors
stepSize=0;
for s=1:3
    figure('Name',P.names{s})
    for iFactor=1:size(FAtrials.FSC_Mean{s}{1},2)% P.nFactors(s)
%             subaxis(ceil(P.nUnits(s)^0.5),floor(P.nUnits(s)^0.5),iUnit)
            axis off; hold on
            % {area}{outcome}[time,Unit]
            ciplot(stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor)+FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   stepSize*iFactor+FAtrials.FSC_Mean{s}{1}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   (1:P.Ltr*2)*P.bw,'r')
% 
           ciplot(stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor)+FAtrials.FSC_SEM{s}{2}(:,iFactor),...
                   stepSize*iFactor+FAtrials.FSC_Mean{s}{2}(:,iFactor)-FAtrials.FSC_SEM{s}{1}(:,iFactor),...
                   (1:P.Ltr*2)*P.bw,'g')

             plot((1:P.Ltr*2)*P.bw,stepSize*iFactor+(D.Assem.TS{s}(:,iFactor)),'k','LineWidth',1.5)


	end

    Y1=[D.Assem.Ft2ciL0{s};D.Assem.Ft2ciH0{s}-D.Assem.Ft2ciL0{s}];
    h=area((1:P.Ltr*2)*P.bw,Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','y'); hold on
	plot((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s},'k','LineWidth',1.5)
end
end
clear Y1
%% Plot relative strength of decoding
if P.flag.PlotVerbose

    Assem_Color = [1 0 0];
Unit_Color  = [0 0 0];
boxHeight = 32;
figure('color','w','name',['Decoding: ' P.rat])
% Plot decoding
for s=1:3
    subplot(2,3,s); hold on
    title(P.names{s})
    for iTrange = 1:size(P.tRanges,1)
          area(P.tRanges(iTrange,:),[boxHeight boxHeight],...
          'FaceColor',P.tRangesColors(iTrange,:),...
          'EdgeColor',P.tRangesColors(iTrange,:),...
          'LineWidth',1.5,...
          'FaceAlpha',0.5,...
          'EdgeAlpha',0.)
    end
%     plot([5 5],[0 10],'color',[0.89 0.95 0.89],'LineWidth',20)
%     plot([25 25],[0 10],'color',[0.95 0.89 0.89],'LineWidth',20)
    
%     Y1=[D.Assem.Ft2ciL0{s};D.Assem.Ft2ciH0{s}-D.Assem.Ft2ciL0{s}];
%     h=area((1:P.Ltr*2)*P.bw,Y1');
%     set(h(1),'FaceColor','w'); set(h(2),'FaceColor','w'); hold on
    plot((1:P.Ltr*2)*P.bw,D.Assem.Ft2ciH0{s},'color',Assem_Color,'LineWidth',1.2,'LineStyle','-')
    plot((1:P.Ltr*2)*P.bw,D.units.Ft2ciH0{s},'color',Unit_Color,'LineWidth',1.2,'LineStyle','-')

	plot((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s},'color',Assem_Color,'LineWidth',1.5)
    plot((1:P.Ltr*2)*P.bw,D.units.Ft2x{s},'color',Unit_Color,'LineWidth',1.5)
% 	area((1:P.Ltr*2)*P.bw,D.units.Ft2x{s},'FaceColor',Unit_C,'EdgeColor',Unit_C,'LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
% 	area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s},'FaceColor',Assem_C,'EdgeColor',Assem_C,'LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)

% area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s}-D.units.Ft2x{s},'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)

    plot([15 15],[0 10],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    set(gca,'XTick',[])
    axis([0 Inf 0 37])
    if s==1
        ylabel({'Decoding score';'(F-value)'})
%     elseif s==2
%         xlabel('Time (s)')
        text(1,boxHeight+2,'Single unit decoding','Color',Unit_Color,'HorizontalAlignment','left' )
        text(1,boxHeight+4,'Assembly decoding','Color',Assem_Color,'HorizontalAlignment','left' )
    end
end
% Plot assembly advantage
for s=1:3
    subplot(2,3,s+3); hold on
    title(P.names{s})
    for iTrange = 1:size(P.tRanges,1)
          area(P.tRanges(iTrange,:),[boxHeight boxHeight],...
          'FaceColor',P.tRangesColors(iTrange,:),...
          'EdgeColor',P.tRangesColors(iTrange,:),...
          'LineWidth',1.5,...
          'FaceAlpha',0.5,...
          'EdgeAlpha',0.)
      
          area(P.tRanges(iTrange,:),[-boxHeight -boxHeight],...
          'FaceColor',P.tRangesColors(iTrange,:),...
          'EdgeColor',P.tRangesColors(iTrange,:),...
          'LineWidth',1.5,...
          'FaceAlpha',0.5,...
          'EdgeAlpha',0.)
    end
    
    temp = (D.Assem.Ft2x{s}-D.units.Ft2x{s}); temp(temp<0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'BaseValue',0,'FaceColor',Assem_Color,'EdgeColor',Assem_Color,'LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    temp = (D.Assem.Ft2x{s}-D.units.Ft2x{s}); temp(temp>0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'BaseValue',0,'FaceColor',Unit_Color,'EdgeColor',Unit_Color,'LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    % area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s}-D.units.Ft2x{s},'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    
%     plot([0 P.Ltr*2*P.bw],[0 0],'color','k','LineWidth',2)

    plot([15 15],[-10 10],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    axis([0 Inf -25 25])
    if s==1
%         ylabel({'Assembly decodimg advantage';'Assembly - Unit'})
        ylabel({'Relative Advantage'})
        text(0.5,2,'Assembly advangage','Color',Assem_Color,'HorizontalAlignment','left','Rotation',90)
        text(0.5,-23,'Single unit advantage','Color',Unit_Color,'HorizontalAlignment','left','Rotation',90 )
 
    elseif s==2
        xlabel('Time (s)')
    end
end
end
%% Sort units/assems profiles by peak time / amp

% unit firing rates
for s=1:3
        
    [FRtrials.iFRpeak{s}(1,:), FRtrials.iFRpeakLoc{s}(1,:)] = max(FRtrials.iFR0_Mean{s}{1});
    [FRtrials.iFRpeak{s}(2,:), FRtrials.iFRpeakLoc{s}(2,:)] = max(FRtrials.iFR0_Mean{s}{2});
    
     FRtrials.iFRpeakGlobal{s} = max(FRtrials.iFRpeak{s});
     FRtrials.iFRpref{s}       = (FRtrials.iFRpeak{s}(2,:)>=FRtrials.iFRpeak{s}(1,:))+1;
     
     for uID = 1:length(FRtrials.iFRpref{s})
         FRtrials.iFRlocGlobal{s}(uID)    =  FRtrials.iFRpeakLoc{s}(FRtrials.iFRpref{s}(uID),uID) ;
     end
     
    [~,FRtrials.iFRPeaktimeSorted{s}] = sort(FRtrials.iFRlocGlobal{s});
    [~,FRtrials.iFRPeakRateSorted{s}] = sort(FRtrials.iFRpeakGlobal{s});
    
    % Sort by decoding power 
    [FRtrials.DPeak{s}(1,:), FRtrials.DPeakLoc{s}(1,:)] = max(D.units.TS{s});
    
    [~,FRtrials.DPeakTimeSorted{s}] = sort(FRtrials.DPeakLoc{s});
    [~,FRtrials.DPeakRateSorted{s}] = sort(FRtrials.DPeak{s});

end

% Sort assembly activation
for s=1:3
    
    [FAtrials.AssemPeak{s}(1,:), FAtrials.AssemPeakLoc{s}(1,:)] = max(FAtrials.FSC_Mean{s}{1});
    [FAtrials.AssemPeak{s}(2,:), FAtrials.AssemPeakLoc{s}(2,:)] = max(FAtrials.FSC_Mean{s}{2});
    
     FAtrials.AssemPeakGlobal{s} = max(FAtrials.AssemPeak{s});
     FAtrials.Assempref{s}       = (FAtrials.AssemPeak{s}(2,:)>=FAtrials.AssemPeak{s}(1,:))+1;
     
     for AssemID = 1:length(FAtrials.Assempref{s})
         FAtrials.AssemlocGlobal{s}(AssemID)    =  FAtrials.AssemPeakLoc{s}(FAtrials.Assempref{s}(AssemID),AssemID) ;
     end
     
    [~,FAtrials.AssemPeakTimeSorted{s}] = sort(FAtrials.AssemlocGlobal{s});
    [~,FAtrials.AssemPeakRateSorted{s}] = sort(FAtrials.AssemPeakGlobal{s});

    % Sort by decoding power 
    [FAtrials.DPeak{s}(1,:), FAtrials.DPeakLoc{s}(1,:)] = max(D.Assem.TS{s});
    
    [~,FAtrials.DPeakTimeSorted{s}] = sort(FAtrials.DPeakLoc{s});
    [~,FAtrials.DPeakRateSorted{s}] = sort(FAtrials.DPeak{s});
end

clear uID AssemID
%% classify units/assems by peak time of activation 
ClassTemp= 1:size(P.tRanges,1);
for s=1:3
    FRtrials.unitClass{s} =zeros(size(FRtrials.DPeakLoc{s}));
    FAtrials.unitClass{s} =zeros(size(FAtrials.DPeakLoc{s}));
    for classID = 1:size(P.tRanges,1)
        % Maximum Likelihood classification of units by peak decoding time
        FRtrials.unitClass{s}(FRtrials.DPeakLoc{s}>P.tRangesAdjusted(classID,1) & FRtrials.DPeakLoc{s}<=P.tRangesAdjusted(classID,2)) ...
            = classID;
        % Maximum Likelihood classification of factors by peak decoding time
        FAtrials.unitClass{s}(FAtrials.DPeakLoc{s}>P.tRangesAdjusted(classID,1) & FAtrials.DPeakLoc{s}<=P.tRangesAdjusted(classID,2)) ...
            = classID;
    end
end    
%% Plot sorted by decoding / FR
if P.flag.PlotVerbose
    figure; hold on
    for s=1:3
        subplot(1,3,s); hold on
        plot(staggerplot((FRtrials.iFR0_Mean{s}{1}(:,FRtrials.iFRPeakRateSorted{s})),0,10),'g');
        plot(staggerplot((FRtrials.iFR0_Mean{s}{2}(:,FRtrials.iFRPeakRateSorted{s})),0,10),'r');
    end

    figure; hold on
    for s=1:3
        subplot(1,3,s); hold on
        plot(staggerplot((FAtrials.FSC_Mean{s}{1}(:,FAtrials.AssemPeakRateSorted{s})),0,5),'g');
        plot(staggerplot((FAtrials.FSC_Mean{s}{2}(:,FAtrials.AssemPeakRateSorted{s})),0,5),'r');
    end
end
%% plot activation scores sorted by peak time

figure('color','w','name',P.rat); hold on
for s=1:3
    subplot(2,3,s); 
    hold on
%     plot(staggerplot((D.units.TS{s}(:,FRtrials.DPeakTimeSorted{s})),0,3),'k');
    x_= repmat((1:P.Ltr*2)*P.bw,size(D.units.TS{s},2),1);
    y_ = repmat(1:size(D.units.TS{s},2),P.Ltr*2,1)';
    z_= zscore(D.units.TS{s}(:,FRtrials.DPeakTimeSorted{s}))';
    imagesc(x_(1:end),y_(1:end),z_)
    colormap(flipud(gray))
    caxis([0 3])
    boxHeight =size(D.units.TS{s},2);
    for iTrange = 1:size(P.tRanges,1)
          area(P.tRanges(iTrange,:),[boxHeight boxHeight],...
                  'FaceColor',P.tRangesColors(iTrange,:),...
                  'EdgeColor',P.tRangesColors(iTrange,:),...
                  'LineWidth',1.5,...
                  'FaceAlpha',0.5,...
                  'EdgeAlpha',0)
    end
    plot([15 15],[1 boxHeight],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    axis([0 30 1 boxHeight]) 
    if s==1
        ylabel('Unit no.')
    elseif s==2
        
        title('Decoding Strength (z-scored)')
    end
end


% figure; hold on
for s=1:3
    try
    subplot(2,3,s+3); 
    hold on
%     plot(staggerplot((D.Assem.TS{s}(:,FAtrials.DPeakTimeSorted{s})),0,3),'k');
    x_= repmat((1:P.Ltr*2)*P.bw,size(D.Assem.TS{s},2),1);
    y_ = repmat(1:size(D.Assem.TS{s},2),P.Ltr*2,1)';
        
    z_= zscore(D.Assem.TS{s}(:,FAtrials.DPeakTimeSorted{s}))';
    imagesc(x_(1:end),y_(1:end),z_)
    % colormap((gray))
    caxis([0 3])
    boxHeight =size(D.Assem.TS{s},2);

    for iTrange = 1:size(P.tRanges,1)
          area(P.tRanges(iTrange,:),[2*boxHeight 2*boxHeight],...
                  'FaceColor',P.tRangesColors(iTrange,:),...
                  'EdgeColor',P.tRangesColors(iTrange,:),...
                  'LineWidth',1.5,...
                  'FaceAlpha',0.5,...
                  'EdgeAlpha',0.)
    end
    plot([15 15],[1 boxHeight],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    axis([0 30 1 boxHeight])
    end
    title(P.names{s})
    if s==1
        ylabel('Assem. no.')
    elseif s==2
        xlabel('Time (s)')
    end
end
clear x_ y_ z_

 fnam=[pat{3}  filesep  'figures' filesep P.rat '_' P.target '_DecodingTime']; % your file name
     snam='8x5';
     s=hgexport('readstyle',snam);
	 s.Format = 'jpeg'; 
     hgexport(gcf,[fnam, '.jpeg'],s);
     sdf('8x5')
     savefig(gcf,[fnam, '.fig'])
%% Plot decoding - both
figure('name',P.rat,'color','w')
maxH = 10;
% Plot assembly and unit decoding
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
 fnam=[pat{3}  filesep  'figures' filesep P.rat '_' P.target '_AssVsUnits']; % your file name
     snam='8x5';
     s=hgexport('readstyle',snam);
	 s.Format = 'jpeg'; 
     hgexport(gcf,[fnam, '.jpeg'],s);
     sdf('8x5')
     savefig(gcf,[fnam, '.fig'])
%% Plot decoding - members vs. non-members
figure('name',P.rat,'color','w')
% Plot decoding
maxH = 10;
LeverLineWidth = 5; 
for s=1:3
    subplot(2,3,s); hold on
    title(P.names{s})

%     plot([5 5],[0 10],'color',[0.89 0.95 0.89],'LineWidth',20)
%     plot([25 25],[0 10],'color',[0.95 0.89 0.89],'LineWidth',20) 
    plot([10 10],[0 10*maxH],'color',[0.89 0.95 0.89],'LineWidth',LeverLineWidth)
    plot([20 20],[0 10*maxH],'color',[0.95 0.89 0.89],'LineWidth',LeverLineWidth)
    
%     Y1=[D.Assem.Ft2ciL0{s};D.Assem.Ft2ciH0{s}-D.Assem.Ft2ciL0{s}];
%     h=area((1:P.Ltr*2)*P.bw,Y1');
%     set(h(1),'FaceColor','w'); set(h(2),'FaceColor','w'); hold on
    plot((1:P.Ltr*2)*P.bw,D.MemberUnitsCollapsed.Ft2ciH0{s},'k')
    plot((1:P.Ltr*2)*P.bw,D.NonMemberUnits.Ft2ciH0{s},'r')

	plot((1:P.Ltr*2)*P.bw,D.MemberUnitsCollapsed.Ft2x{s},'k','LineWidth',1.5)
    plot((1:P.Ltr*2)*P.bw,D.NonMemberUnits.Ft2x{s},'r','LineWidth',1.5)
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
        text(15,-1,'Assembly members','Color',[0 0 0],'HorizontalAlignment','center' )
        text(15,-3,'Assembly non-members','Color',[1 0 0],'HorizontalAlignment','center' )
        axis([0 Inf 0 2*maxH])
    elseif s==3      
        axis([0 Inf 0 2*maxH])
    end
end
% Plot assembly advantage
for s=1:3
    subplot(2,3,s+3); hold on
%     title(P.names{s})

    plot([9 9],[-10 10],'color',[0.89 0.95 0.89],'LineWidth',LeverLineWidth)
    plot([20 20],[-10 10],'color',[0.95 0.89 0.89],'LineWidth',LeverLineWidth)
    
    temp = D.MemberUnitsCollapsed.Ft2x{s}-D.NonMemberUnits.Ft2x{s}; temp(temp>0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'FaceColor',[0.9 0.6 0.6],'EdgeColor','r','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    temp = D.MemberUnitsCollapsed.Ft2x{s}-D.NonMemberUnits.Ft2x{s}; temp(temp<0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'FaceColor',[0.6 0.6 0.6],'EdgeColor','k','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    % area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s}-D.units.Ft2x{s},'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    
%     plot([0 P.Ltr*2*P.bw],[0 0],'color','k','LineWidth',2)
    plot([15 15],[-10 10],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[-10 10],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    axis([0 Inf -10 10])
    if s==1
%         ylabel({'Assembly decodimg advantage';'Assembly - Unit'})
        ylabel({'Assembly decodimg advantage'})
    elseif s==2
        xlabel('Time (s)')
    end
end
 fnam=[pat{3}  filesep  'figures' filesep P.rat '_' P.target '_membersVsnonmembers']; % your file name
     snam='8x5';
     s=hgexport('readstyle',snam);
	 s.Format = 'jpeg'; 
     hgexport(gcf,[fnam, '.jpeg'],s);
     sdf('8x5')
     savefig(gcf,[fnam, '.fig'])
%% Save 
save([pat{3} filesep P.rat '_' P.target '_decoding.mat'])
% Group analysis for single unit / assembly decoding capability, activation
% patterns for assemblies
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
%% Preamble, get file list
clear all; close all
set(0,'DefaultFigureWindowStyle','docked')
flags.plotIndividualChanges = true;
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

Ignore= {'JaroslawLONG1'; ...
         'MiroslawLONG2'; ...
         'NorbertLONG2' ; ...
         'OnufryLONG2'};
     
% Ignore= {};

for iList = 1:length(Ignore)
    Rats = Rats(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
           {Rats.name},repmat({Ignore(iList)},1,length(Rats)),'UniformOutput',false)));
end
clear Ignore iList 
%% Process assembly activation averages
Group = struct;
Group.Assem.FSC = {};
Group.Assem.FSC_Mean = cell(3,1);
Group.Assem.FSC_SEM = cell(3,1);
for iRat= 1:length(Rats)
    load([pat{1} Rats(iRat).name],'P','D','FAtrials');
    disp(['Working on ' Rats(iRat).name '...'])
    for s=1:length(P.names) 
        if ~isempty(FAtrials.FSC_{s}(1))
    temp = cell2mat(FAtrials.FSC_{s}'); temp_m=[]; temp_s=[];
    for i = 1:size(temp,2)
        temp_m (:,i) = nanmean(reshape(temp(:,i),[600,length(temp(:,i))/600])');
        temp_s (:,i) = nansem(reshape(temp(:,i),[600,length(temp(:,i))/600])');
    end
    Group.Assem.FSC_Mean{s} = [Group.Assem.FSC_Mean{s}; temp_m'];
    Group.Assem.FSC_SEM{s}  = [Group.Assem.FSC_SEM{s}; temp_m'];
    
        end
    end
end
clear iRat s i temp_m temp_s temp
%%

figure('color','w'); hold on
tb = (1:P.Ltr*2)*P.bw;
clr2={'b','g','r'};

for s= 1:3
    temp_m  = Group.Assem.FSC_Mean{s};
    temp_s  = Group.Assem.FSC_SEM{s};
    clr_    = rand(size(temp_m,1),3);
    for i=1:size(temp_m,1)
        subplot(1,3,s); hold on
%         ciplot(5*i+temp_m(i,:)-temp_s(i,:),...
%             5*i+temp_m(i,:)+temp_s(i,:),...
%             tb,clr_(i,:));
        plot(tb,5*i+temp_m(i,:),'color',clr_(i,:),'LineWidth',1.5)
        %     plot(tb,staggerplot(temp_m',0, 5))
        %     axis([min(tb) max(tb), -10 10])
    end
end
clear clr2 clr_ temp_m temp_s tb i s
%%

figure('color','w'); hold on
tb = (1:P.Ltr*2)*P.bw;
clr2={'b','g','r'};

for s= 1:3
   temp_m  = Group.Assem.FSC_Mean{s};
   temp_s  = Group.Assem.FSC_SEM{s};


   ciplot(mean(temp_m)+nansem(temp_m),...
          mean(temp_m)-nansem(temp_m),...
          tb,clr2{s},0.6); 
%    plot(tb,mean(temp_m),'color',clr2{s},'LineWidth',1.5)
    axis([min(tb) max(tb), -2 5])
   
end
legend(P.names)
xlabel('Time (s)')
ylabel('Mean Assembly activation strength')
clear clr2 clr_ temp_m temp_s tb i s
%%
   ylims1 = get(gca,'Ylim');
   subplot(1,2,2); hold on
   ciplot(v.spikes.IncorrectResponseiFRmean(iUnit,:)+v.spikes.IncorrectResponseiFRSEM(iUnit,:),...
          v.spikes.IncorrectResponseiFRmean(iUnit,:)-v.spikes.IncorrectResponseiFRSEM(iUnit,:),...
          v.spikes.ResponsePSTH_tb',clr_(iUnit,:));
   ylims2 = get(gca,'Ylim');

end
subplot(1,2,1);
plot([0 0],max([ylims1;ylims2]),':k')
axis([min(v.spikes.ResponsePSTH_tb) max(v.spikes.ResponsePSTH_tb) max([ylims1;ylims2])])
plot([1 1],0.9*max([ylims1;ylims2]),'b')
text(1.2,max(0.9*max([ylims1;ylims2])),{'Reward';'Delivery'},'Color','b')

title('Correct arm entry')
ylabel('z-scored Firing Rate (S.D.)')
xlabel('Peri-arm entry time (s)')
subplot(1,2,2);
plot([0 0],max([ylims1;ylims2]),':k')
axis([min(v.spikes.ResponsePSTH_tb) max(v.spikes.ResponsePSTH_tb) max([ylims1;ylims2])])
title('Incorrect arm entry')
xlabel('Peri-arm entry time (s)')



%% Process decoding results
for iRat= 1:length(Rats)
    load([pat{1} Rats(iRat).name],'P','D','FAtrials');
    disp(['Working on ' Rats(iRat).name '...'])
    for s=1:length(P.names) 
        % Collate decoding results
        Group.Units.Ft2ciH0{s} (iRat,:)      = D.units.Ft2ciH0{s};
        Group.Units.Ft2x{s}    (iRat,:)      = D.units.Ft2x{s};
        Group.Assem.Ft2ciH0{s} (iRat,:)      = D.Assem.Ft2ciH0{s};
        Group.Assem.Ft2x{s}    (iRat,:)      = D.Assem.Ft2x{s};
        Group.Assem.TSPeak{s}  {iRat,:}      = D.Assem.TSPeak{s};
        Group.Assem.TSIntegral{s}{iRat,:}    = D.Assem.TSIntegral{s};
%         Group.AssemAdvantage{s}(iRat,:)      = (D.Assem.Ft2x{s}-D.units.Ft2x{s})./(D.Assem.Ft2x{s}+D.units.Ft2x{s});
        Group.AssemAdvantage{s}(iRat,:)      = D.Assem.Ft2x{s}-D.units.Ft2x{s};
    end
    % Get transition matrix for delay period
    Group.Assem.Task.DelayTransMatrix{iRat} = (FAtrials.DelayOnly.patterns.ptn.transMatrix_sortedMean{4}{1} + FAtrials.DelayOnly.patterns.ptn.transMatrix_sortedMean{4}{2})/2;
end
% Collapse decoding
Group.Units.Ft2ciH0_Mean  = cellfun(@nanmean,Group.Units.Ft2ciH0,'UniformOutput',false);
Group.Units.Ft2x_Mean     = cellfun(@nanmean,Group.Units.Ft2x,'UniformOutput',false);
Group.Units.Ft2x_SEM      = cellfun(@nansem,Group.Units.Ft2x,'UniformOutput',false);

Group.Assem.Ft2ciH0_Mean  = cellfun(@nanmean,Group.Assem.Ft2ciH0,'UniformOutput',false);
Group.Assem.Ft2x_Mean     = cellfun(@nanmean,Group.Assem.Ft2x,'UniformOutput',false);
Group.Assem.Ft2x_SEM      = cellfun(@nansem,Group.Assem.Ft2x,'UniformOutput',false);

Group.AssemAdvantage_Mean = cellfun(@nanmean,Group.AssemAdvantage,'UniformOutput',false);
Group.AssemAdvantage_SEM  = cellfun(@nansem,Group.AssemAdvantage,'UniformOutput',false);

clear iRat s D FAtrials
%% Plot group decoding (units/assemblies')
figure('name','Group decoding','color','w')
% Plot decoding
maxH = 20;              % yMax of decoding plots
LeverPressWidth = 5;    % Width of lever press bar (s)
MinMaxAdvantage = 10;   % yMin and yMax of assembly advantage plot 
SAMPLEln  = [0.89 0.95 0.89];
CHOICEln  = [0.95 0.89 0.89];
alpha = 0.1;
for s=1:3
    subplot(2,3,s); hold on
    title(P.names{s})

    plot([9 9],[0 10*maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 10*maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    
    plot((1:P.Ltr*2)*P.bw,Group.Assem.Ft2ciH0_Mean{s},'k')
    plot((1:P.Ltr*2)*P.bw,Group.Units.Ft2ciH0_Mean{s},'b')

    ciplot(Group.Assem.Ft2x_Mean{s}+Group.Assem.Ft2x_SEM{s},...
           Group.Assem.Ft2x_Mean{s}-Group.Assem.Ft2x_SEM{s},...
           (1:P.Ltr*2)*P.bw,[0.6 0.6 0.6],alpha)  
    ciplot(Group.Units.Ft2x_Mean{s}+Group.Units.Ft2x_SEM{s},...
           Group.Units.Ft2x_Mean{s}-Group.Units.Ft2x_SEM{s},...
           (1:P.Ltr*2)*P.bw,[0.6 0.6 0.9],alpha)
	plot((1:P.Ltr*2)*P.bw,Group.Assem.Ft2x_Mean{s},'k','LineWidth',1.5)
    plot((1:P.Ltr*2)*P.bw,Group.Units.Ft2x_Mean{s},'b','LineWidth',1.5)
    plot([15 15],[0 10*maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 10*maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    set(gca,'XTick',[])
    
    axis([0 Inf 0 maxH ])
    if s==1
        ylabel({'Decoding score';'(F-value)'})
        axis([0 Inf 0 maxH/2 ])

    elseif s==2
%         xlabel('Time (s)')
        text(15,-1,'Single unit decoding','Color',[0 0 1],'HorizontalAlignment','center' )
        text(15,-3,'Assembly decoding','Color',[0 0 0],'HorizontalAlignment','center' )
        axis([0 Inf 0 maxH ])
    elseif s==3
        axis([0 Inf 0 maxH ])

    end
end
% Plot assembly advantage
for s=1:3
    subplot(2,3,s+3); hold on

    plot([9 9],[-10 10],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[-10 10],'color',CHOICEln,'LineWidth',LeverPressWidth)
    
    temp = Group.AssemAdvantage_Mean{s}; temp(temp>0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',alpha,'EdgeAlpha',1)
    
    temp = Group.AssemAdvantage_Mean{s}; temp(temp<0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'FaceColor',[0.6 0.6 0.6],'EdgeColor','k','LineWidth',1.5,'FaceAlpha',alpha,'EdgeAlpha',1)
        plot([15 15],[-10 10],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[-10 10],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    axis([0 Inf -MinMaxAdvantage MinMaxAdvantage])
    if s==1
        ylabel({'Assembly decoding advantage'})
    elseif s==2
        xlabel('Time (s)')
    end
end

%  fnam=[pat{1}  'Group_AssVsUnits']; % your file name
%      snam='8x5';
%      s=hgexport('readstyle',snam);
% 	 s.Format = 'jpeg'; 
%      hgexport(gcf,[fnam, '.jpeg'],s);
%      sdf('8x5')
%      savefig(gcf,[fnam, '.fig'])
%% Process epoch change variables for units and factors
hasNoPostSleep = false;
epochs  = {'Pre','Task','Post'}; 
epochsD = {'PreVsTask','TaskVsPost','PreVsPost'};
for iRat= 1:length(Rats)
    load([pat{1} Rats(iRat).name],'P','D','FAcont','FRcont','FAtrials');
   
    hasNoPostSleep =  ~isfield(FAcont.Post.patterns,'filename');
    if hasNoPostSleep 
        disp(['Working on ' Rats(iRat).name '... no Post sleep'])
    else
        disp(['Working on ' Rats(iRat).name '... '])
    end
    
    
    for s=1:length(P.names)
       
        
        % Change between epochs
        for iEpoch = 1:length(epochsD)
            % Collate epoch to epoch changes in assembly strength and repeat rate
            if  hasNoPostSleep && ~isempty(strfind(epochsD{iEpoch},'Post'))
%                 has Post sleep and comparison doesn't include Post time
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_repMode{s}{iRat,:}          = [];'])
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_LogStrengthQ3bin{s}{iRat,:} = [];'])
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_ModeLogStrength{s}{iRat,:}  = [];'])
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_Eccentric{s}{iRat,:}        = [];'])
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_transMatrix{iRat,:}         = [];'])
            else
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_repMode{s}{iRat,:}          = FAcont.' epochsD{iEpoch} '.pks.delta_repMode{s};'])
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_LogStrengthQ3bin{s}{iRat,:} = FAcont.' epochsD{iEpoch} '.pks.delta_LogStrengthQ3bin{s};'])
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_ModeLogStrength{s}{iRat,:}  = FAcont.' epochsD{iEpoch} '.pks.delta_ModeLogStrength{s};'])
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_Eccentric{s}{iRat,:}        = FAcont.' epochsD{iEpoch} '.pks.delta_Eccentric{s};'])
                eval(['Group.Assem.' epochsD{iEpoch} '.delta_transMatrix{iRat,:}         = FAcont.' epochsD{iEpoch} '.ptn.delta_transMatrix(1:end);'])
            end
        end
        

        
        
        % Each epoch
        for iEpoch = 1:length(epochs)
            if  hasNoPostSleep && ~isempty(strfind(epochs{iEpoch},'Post'))
                 % Collate assembly activation rates
                eval(['Group.Assem.'  epochs{iEpoch} '.mean_repMedian{s}{iRat,:}         = [];'])
                eval(['Group.Assem.'  epochs{iEpoch} '.mean_repMode{s}{iRat,:}           = [];'])
                % Collate assembly activation rhytmicity
                eval(['Group.Assem.'  epochs{iEpoch} '.mean_repSEM{s}{iRat,:}            = [];'])
                eval(['Group.Assem.'  epochs{iEpoch} '.mean_repSTD{s}{iRat,:}            = [];'])
                % Collate assembly activation strength
                eval(['Group.Assem.'  epochs{iEpoch} '.LogStrengthQ3bin{s}{iRat,:}       = [];'])
                eval(['Group.Assem.'  epochs{iEpoch} '.ModeLogStrength{s}{iRat,:}        = [];'])

                % Collate average unit firing rate
                unitIDs = FAcont.Task.unitIDs{s};
                eval(['Group.Units.' epochs{iEpoch} '.meanFR{s}{iRat,:}  = [];'])

                % Collate assembly membership information (no. Assems each unit is member of)
                eval(['Group.Units.' epochs{iEpoch} '.AssOverlap{s}{iRat,:} = [];'])
                
            else
                % Collate assembly activation rates
                eval(['Group.Assem.'  epochs{iEpoch} '.mean_repMedian{s}{iRat,:}         = FAcont.'  epochs{iEpoch} '.patterns.pks.repMedian{s};'])
                eval(['Group.Assem.'  epochs{iEpoch} '.mean_repMode{s}{iRat,:}           = FAcont.'  epochs{iEpoch} '.patterns.pks.repMode{s};'])
                % Collate assembly activation rhytmicity
                eval(['Group.Assem.'  epochs{iEpoch} '.mean_repSEM{s}{iRat,:}            = FAcont.'  epochs{iEpoch} '.patterns.pks.repSEM{s};'])
                eval(['Group.Assem.'  epochs{iEpoch} '.mean_repSTD{s}{iRat,:}            = FAcont.'  epochs{iEpoch} '.patterns.pks.repSTD{s};'])
                % Collate assembly activation strength
                eval(['Group.Assem.'  epochs{iEpoch} '.LogStrengthQ3bin{s}{iRat,:}       = FAcont.'  epochs{iEpoch} '.patterns.pks.LogStrengthQ3bin{s};'])
                eval(['Group.Assem.'  epochs{iEpoch} '.ModeLogStrength{s}{iRat,:}        = FAcont.'  epochs{iEpoch} '.patterns.pks.ModeLogStrength{s};'])

                % Collate average unit firing rate
                unitIDs = FAcont.Task.unitIDs{s};
                if    s<length(P.names)
                    eval(['Group.Units.' epochs{iEpoch} '.meanFR{s}{iRat,:}  = FRcont.' epochs{iEpoch} '{s}.avgFR(unitIDs);'])
                elseif s==length(P.names)
                    eval(['Group.Units.' epochs{iEpoch} '.meanFR{3}{iRat,:} = [Group.Units.' epochs{iEpoch} '.meanFR{1}{iRat},Group.Units.' epochs{iEpoch} '.meanFR{2}{iRat}];'])
                end

                % Collate assembly membership information (no. Assems each unit is member of)
                eval(['Group.Units.' epochs{iEpoch} '.AssOverlap{s}{iRat,:} =  FAcont.' epochs{iEpoch} '.patterns.units.AssOverlap{s};'])
            end
        end
       
    end
    
       
end


% Collapse assembly membership vs excitability
for s=1:length(P.names)
    Rats = find(~cellfun(@isempty,Group.Units.Pre.AssOverlap{s}));
    Group.Units.Pre.meanFR_{s}     = cell2mat(Group.Units.Pre.meanFR{s}(Rats)');
    Group.Units.Pre.AssOverlap_{s} = cell2mat(Group.Units.Pre.AssOverlap{s}(Rats));
    
    Rats = find(~cellfun(@isempty,Group.Units.Task.AssOverlap{s}));
    Group.Units.Task.meanFR_{s}     = cell2mat(Group.Units.Task.meanFR{s}(Rats)');
    Group.Units.Task.AssOverlap_{s} = cell2mat(Group.Units.Task.AssOverlap{s}(Rats));
    
    Rats = find(~cellfun(@isempty,Group.Units.Post.AssOverlap{s}));
    Group.Units.Post.meanFR_{s}     = cell2mat(Group.Units.Post.meanFR{s}(Rats)');
    Group.Units.Post.AssOverlap_{s} = cell2mat(Group.Units.Post.AssOverlap{s}(Rats));
end
% Averages
for iEpoch = 1:length(epochs)
    for s=1:length(P.names)
        eval(['U_ = unique(Group.Units.' epochs{iEpoch} '.AssOverlap_{s});'])
        for Uid = 1:length(U_)
            eval(['units = find(Group.Units.' epochs{iEpoch} '.AssOverlap_{s} == U_(Uid));'])
            % [No assems , mean FR, SEM FR]
             eval(['Group.Units.' epochs{iEpoch} '.nAssVsFR{s}(1,Uid) = U_(Uid);'])
             eval(['Group.Units.' epochs{iEpoch} '.nAssVsFR{s}(2,Uid) = nanmean(Group.Units.' epochs{iEpoch} '.meanFR_{s}(units));'])
             eval(['Group.Units.' epochs{iEpoch} '.nAssVsFR{s}(3,Uid) = nansem(Group.Units.' epochs{iEpoch} '.meanFR_{s}(units));'])
            
        end
    end
end



clear FAcont FAtrials FRcont D iEpoch iRat s unitIDs U_ Uid units
%% optional plotting 
if flags.plotIndividualChanges
%% Plot change in assembly activation RATE vs peak decoding score in task (Pre-Post)
figure('name','Change in assembly activation rate vs. decoding score');
ymax = 20;
xminmax = 1.5;
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 0],[0 ymax],':k')
    tempx = 1./cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_repMode{s}','UniformOutput',false)');
    tempy = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
%     try
        scatter(tempx,tempy,'LineWidth',1,'MarkerEdgeColor','k')
%     catch
%     end
    title(P.names{s})
    if s==1
        ylabel('Assembly peak decoding score')
    elseif s==2
        xlabel('Pre-task to Post-task change in activation rate (Hz)')
    end
    
    axis([-xminmax  xminmax  0 ymax])
    axis square
end
%% Plot change in assembly activation RATE vs peak decoding score in task (Task-Post)
figure('name','Change in assembly activation rate vs. decoding score');
ymax = 20;
xminmax = 1.5;
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 0],[0 ymax],':k')
    tempx = 1./cell2mat(cellfun(@transpose,Group.Assem.TaskVsPost.delta_repMode{s}','UniformOutput',false)');
    tempy = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
%     try
        scatter(tempx,tempy,'LineWidth',1,'MarkerEdgeColor','k')
%     catch
%     end
    title(P.names{s})
    if s==1
        ylabel('Assembly peak decoding score')
    elseif s==2
        xlabel('Task to Post-task change in activation rate (Hz)')
    end
    
    axis([-xminmax  xminmax  0 ymax])
        axis square

end
%% Plot change in assembly activation STRENGTH vs peak decoding score in task (Pre-Post)
ymax = 20;
xminmax = 1.5;
figure('name','Change in assembly activation strength vs. decoding score');
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 0],[0 ymax],':k')
    tempx = cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_LogStrengthQ3bin{s}','UniformOutput',false)');
    tempy = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    
    scatter (tempx,tempy,'LineWidth',1,'MarkerEdgeColor','k')

    title(P.names{s})
    if s==1
        ylabel('Assembly peak decoding score')
    elseif s==2
        xlabel('Pre-Task to Post-task change in activation strength')
    end   
    axis([-xminmax  xminmax  0 ymax])
    axis square

end
%% Plot change in assembly activation STRENGTH vs peak decoding score in task (Task-Post)
ymax = 20;
xminmax = 1.5;
figure('name','Change in assembly activation strength vs. decoding score');
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 0],[0 ymax],':k')
    tempx = cell2mat(cellfun(@transpose,Group.Assem.TaskVsPost.delta_LogStrengthQ3bin{s}','UniformOutput',false)');
    tempy = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    
    scatter (tempx,tempy,'LineWidth',1,'MarkerEdgeColor','k')

    title(P.names{s})
    if s==1
        ylabel('Assembly peak decoding score')
    elseif s==2
        xlabel('Task to Post-task change in activation strength')
    end   
    axis([-xminmax  xminmax  0 ymax])
    axis square
end
end
%% Plot repeat rate stability vs decoding score
ymax = 20;
xminmax = 1.5;
figure('name','Task assembly activation rhytmicity vs. decoding score');
for s = 1:3
    subplot(1,3,s); hold on
    tempx = 1./cell2mat(cellfun(@transpose, Group.Assem.Post.mean_repSTD {s}','UniformOutput',false)');
	tempy = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    scatter (tempx,tempy,'LineWidth',1,'MarkerEdgeColor','k')
    

    x = [min(tempx) max(tempx)];
    f = ezfit(tempx,tempy,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.6*max(x), 0.8*max(ymax),strcat('R^2=', sprintf('%3.1g',f.r^2)),'color','r')
    
%     set(gca,'XScale','log')
    
    title(P.names{s})
    if s==1
        ylabel('Assembly peak decoding score')
    elseif s==2
        xlabel('Post-task assembly activation rate stability')
    end   
    axis([0 0.1 0 ymax])
	axis square
end
%% Plot repeat rate and stability as function of decoding score
ymax = 0.15;

figure('name','Change in assembly activation rate vs. decoding score');
for s = 1:3
    subplot(1,3,s); hold on
    tempx = 1./cell2mat(cellfun(@transpose, Group.Assem.Post.mean_repMode{s}','UniformOutput',false)');
    tempy = 1./cell2mat(cellfun(@transpose, Group.Assem.Post.mean_repSTD{s}','UniformOutput',false)');
%     tempy = cell2mat(cellfun(@transpose, Group.Assem.Post.LogStrengthQ3bin{s}','UniformOutput',false)');
    tempz = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');

    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];   
    tempz(isnan(tempz)) = [];
    
    scatter(tempx,tempy,30,tempz,'filled','LineWidth',0.8,'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.9)
   
    x = [min(tempx) max(tempx)];
    f = ezfit(tempx,tempy,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.6*max(x), 0.8*max(ymax),strcat('R^2=', sprintf('%3.1g',f.r^2)),'color','r')
    
    caxis([0 10])
    axis([0 0.15 0 ymax])
    title([P.names{s} ' Assemblies'])
    if s==1
        ylabel('Repeat rate stability')
%         ylabel('Assembly peak decoding score')
    elseif s==2
        xlabel('Average repeat rate (Hz) ')
    end   
    axis square
%    colormap gray
%     set(gca,'XScale','log')
end
colormap parula
c = colorbar('Location','manual','Position',[0.92 0.3 0.02 0.2]);
c.Label.String = {'Decoding power'};

%% Plot repeat Pre  vs. Post strength vs. decoding score
xymax = 15;
P.names = {'mPFC','dCA1','Inter-Regional'};
figure('name','Change in assembly activation strength vs. decoding score');
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 xymax], [0 xymax],':k')
    tempx = cell2mat(cellfun(@transpose,  Group.Assem.Pre.LogStrengthQ3bin{s}','UniformOutput',false)');
%     tempy = 1./cell2mat(cellfun(@transpose,  Group.Assem.Post.mean_repMode{s}','UniformOutput',false)');
        tempy = cell2mat(cellfun(@transpose, Group.Assem.Post.LogStrengthQ3bin{s}','UniformOutput',false)');
    tempz = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
    scatter (tempx,tempy,20+5*tempz,tempz,'filled','LineWidth',1,'MarkerEdgeColor','k')
    [R_,P_] = corrcoef(tempx,tempy);

    x = [min(tempx) max(tempx)];
    f = ezfit(tempx,tempy,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.1*max(x), 0.8*max(xymax),{strcat('R^2=', sprintf('%3.1g',R_(1,2)^2));strcat('p=', sprintf('%3.1g',P_(1,2)))},'color','r')%f.r^2)
  
    
    
    caxis([0 10])
    axis([0 xymax 0 xymax])
%     set(gca,'XScale','log','YScale','log')
     title(P.names{s})
    if s==1
        ylabel('Post-Task rate strength')
    elseif s==2
        xlabel('Pre-Task assembly activation strength')    
    end   
    axis square
end
colormap parula
c = colorbar('Location','East','Position',[0.75 0.28 0.03 0.15]);
c.Label.String = {'Decoding power'};
%% Plot repeat Pre  vs. Post rate vs. decoding score
xymax = 0.55;
P.names = {'mPFC','dCA1','Inter-Regional'};
figure('name','Change in assembly activation rate vs. decoding score');
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 xymax], [0 xymax],':k')
    tempx = 1./cell2mat(cellfun(@transpose,  Group.Assem.Pre.mean_repMode{s}','UniformOutput',false)');
    tempy = 1./cell2mat(cellfun(@transpose,  Group.Assem.Post.mean_repMode{s}','UniformOutput',false)');
    %     tempy = cell2mat(cellfun(@transpose, Group.Assem.Post.LogStrengthQ3bin{s}','UniformOutput',false)');
    tempz = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
    scatter (tempx,tempy,20+5*tempz,tempz,'filled','LineWidth',1,'MarkerEdgeColor','k')
    [R_,P_] = corrcoef(tempx,tempy);

    x = [min(tempx) max(tempx)];
    f = ezfit(tempx,tempy,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.1*max(x), 0.8*max(xymax),{strcat('R^2=', sprintf('%3.1g',R_(1,2)^2));strcat('p=', sprintf('%3.1g',P_(1,2)))},'color','r')%f.r^2)
  
    
    
    caxis([0 10])
%     axis([0 xymax 0 xymax])
%     set(gca,'XScale','log','YScale','log')
     title(P.names{s})
    if s==1
        ylabel('Post-Task rate (Hz)')
    elseif s==2
        xlabel('Pre-Task assembly activation rate (Hz)')    
    end   
    axis square
end
colormap parula
c = colorbar('Location','East','Position',[0.75 0.28 0.03 0.15]);
c.Label.String = {'Decoding power'};
%% Plot change in assembly activation STRENGTH and rate (Pre-Post)
figure('name','PreVsPost rate/strength changes','NumberTitle','off');
xminmax = 2;
yminmax = 10;
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 0],[-yminmax yminmax],':k')
    plot([-xminmax xminmax],[0 0],':k')
    tempx = 1./cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_repMode{s}','UniformOutput',false)');
    tempy =    cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_LogStrengthQ3bin{s}','UniformOutput',false)');
    tempz =    cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
    
    scatter (tempx,tempy,40,tempz,'filled','LineWidth',1,'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.9)
    title(P.names{s})
    if s==1
        ylabel('\Delta strength (A.U.)')
    elseif s==2
        xlabel('\Delta rate (Hz)')
    end   
    axis square
    caxis([0 10])
end

    colormap parula
    c = colorbar('Location','manual','Position',[0.92 0.3 0.02 0.2]);
    c.Label.String = {'Decoding capacity (Peak t-score)'};

%% Pre activation rate vs. post-Task activation rate
figure('name','Task activation rate vs. post-Task activation rate'); hold on
for s=1:3
    subplot(1,3,s)
 hold on
    plot([0 xymax], [0 xymax],':k')
    tempx = 1./cell2mat(cellfun(@transpose,  Group.Assem.Pre.mean_repMode{s}','UniformOutput',false)');
    tempy = 1./cell2mat(cellfun(@transpose,  Group.Assem.Post.mean_repMode{s}','UniformOutput',false)');
    %     tempy = cell2mat(cellfun(@transpose, Group.Assem.Post.LogStrengthQ3bin{s}','UniformOutput',false)');
    tempz = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
    scatter (tempx,tempy,20+5*tempz,tempz,'filled','LineWidth',1,'MarkerEdgeColor','k')
    
    caxis([0 10])
    axis([0 xymax 0 xymax])
%     set(gca,'XScale','log','YScale','log')
     title(P.names{s})
     if s==1
        ylabel('Assembly activation rate: Post-Task (Hz)')
     elseif s==2
        xlabel('Assembly activation rate: Task (Hz)')    
     end
    axis square    
    colormap parula
end
c = colorbar('Location','East','Position',[0.75 0.28 0.03 0.15]);
c.Label.String = {'Decoding power'};
%% Plot repeat Task vs. Post strength vs. decoding score
xymax = 3;
figure('name','Change in assembly activation strength vs. decoding score');
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 xymax], [0 xymax],':k')
    tempx = 1./cell2mat(cellfun(@transpose, Group.Assem.Pre.ModeLogStrength{s}','UniformOutput',false)');
    tempy = 1./cell2mat(cellfun(@transpose, Group.Assem.Post.ModeLogStrength{s}','UniformOutput',false)');
    %     tempy = cell2mat(cellfun(@transpose, Group.Assem.Post.LogStrengthQ3bin{s}','UniformOutput',false)');
    tempz = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
        
    tempy(tempx>8) = [];
    tempz(tempx>8) = [];
    tempx(tempx>8) = [];
    scatter (tempx,tempy,20+5*tempz,tempz,'filled','LineWidth',1,'MarkerEdgeColor','k')
    [R_,P_] = corrcoef(tempx,tempy);

    x = [min(tempx) max(tempx)];
    f = ezfit(tempx,tempy,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.1*max(x), 0.8*max(xymax),{strcat('R^2=', sprintf('%3.1g',R_(1,2)^2));strcat('p=', sprintf('%3.1g',P_(1,2)))},'color','r')%f.r^2)
  
    caxis([0 10])
    axis([0 xymax 0 xymax])
%     set(gca,'XScale','log','YScale','log')
     title(P.names{s})
    if s==1
        ylabel('Post-Task strength')    
    elseif s==2
        xlabel('Pre-Task assembly activation strength')    
    end   
    axis square
    
end
c = colorbar('Location','East','Position',[0.75 0.28 0.03 0.15]);
c.Label.String = {'Decoding power'};
%% Plot change in assembly activation STRENGTH and rate (Task-Post)
figure('name','TaskVsPost rate/strength changes','NumberTitle','off');
xminmax = 1.5;
yminmax = 10;
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 0],[-yminmax yminmax],':k')
    plot([-xminmax xminmax],[0 0],':k')
    tempx = 1./cell2mat(cellfun(@transpose,Group.Assem.TaskVsPost.delta_repMode{s}','UniformOutput',false)');
    tempy = cell2mat(cellfun(@transpose,Group.Assem.TaskVsPost.delta_LogStrengthQ3bin{s}','UniformOutput',false)');
    tempz = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
    
    scatter (tempx,tempy,50,tempz,'filled','LineWidth',1,'MarkerEdgeColor','k')
    title(P.names{s})
    if s==1
        ylabel('\Delta strength')
    elseif s==2
        xlabel('\Delta rate')
    end   
    axis([-xminmax xminmax -yminmax yminmax])
    axis square
    caxis([0 10])
end
 colormap parula
    c = colorbar('Location','manual','Position',[0.92 0.3 0.02 0.2]);
    c.Label.String = {'Decoding power'};
%% Plot repeat Task vs. Post rate vs. decoding score
xymax = 0.2;
P.names = {'mPFC','dCA1','Inter-Regional'};
figure('name','Change in assembly activation rate vs. decoding score');
for s = 1:3
    subplot(1,3,s); hold on
    plot([0 xymax], [0 xymax],':k')
    tempx = 1./cell2mat(cellfun(@transpose,  Group.Assem.Task.mean_repMode{s}','UniformOutput',false)');
    tempy = 1./cell2mat(cellfun(@transpose,  Group.Assem.Post.mean_repMode{s}','UniformOutput',false)');
    %     tempy = cell2mat(cellfun(@transpose, Group.Assem.Post.LogStrengthQ3bin{s}','UniformOutput',false)');
    tempz = cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
    scatter (tempx,tempy,20+5*tempz,tempz,'filled','LineWidth',1,'MarkerEdgeColor','k')
    [R_,P_] = corrcoef(tempx,tempy);

    x = [min(tempx) max(tempx)];
    f = ezfit(tempx,tempy,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.1*max(x), 0.8*max(xymax),{strcat('R^2=', sprintf('%3.1g',R_(1,2)^2));strcat('p=', sprintf('%3.1g',P_(1,2)))},'color','r')%f.r^2)
  
    
    
    caxis([0 10])
    axis([0 xymax 0 xymax])
%     set(gca,'XScale','log','YScale','log')
     title(P.names{s})
    if s==1
        ylabel('Post-Task rate (Hz)')
    elseif s==2
        xlabel('Task assembly activation rate (Hz)')    
    end   
    axis square
end
colormap parula
c = colorbar('Location','East','Position',[0.75 0.28 0.03 0.15]);
c.Label.String = {'Decoding power'};

%% Plot excitability vs no factors each unit assigned to
figure('name','Pre: Excitabiity vs assembly membership','color','w');
for s = 1:3
    subplot(1,3,s); hold on
    bar(Group.Units.Task.nAssVsFR{s}(1,:),Group.Units.Task.nAssVsFR{s}(2,:),...
        'FaceColor','w','EdgeColor','b','LineWidth',1.5)

    
    x=Group.Units.Task.AssOverlap_{s};
    y=Group.Units.Task.meanFR_{s};
    scatter(x+0.05*randn(size(x)),y,10,[0.8 0.8 0.99])
    errorbar(Group.Units.Task.nAssVsFR{s}(1,:),...
             Group.Units.Task.nAssVsFR{s}(2,:),...
             Group.Units.Task.nAssVsFR{s}(3,:),'LineWidth',2)
%     axis([-1 3 0 10])
    title(P.names{s})
%     set(gca,'XTick',[0 1 2])
    if s==1
        ylabel('Mean firing rate (Spikes/s)')
        set(gca,'XTick',[0 1 2])
    elseif s==2
        xlabel('No. assemblies each unit assigned to')
        set(gca,'XTick',[0 1 2])
    elseif s==3
        set(gca,'XTick',0:4)
%         axis([-1 5 0 10])

    end
end

names = {'0','1','2','3','4'};
Colors = jet(length(names));
figure('name','Pre: assembly overlap','color','w');
for s = 1:3
    subplot(3,1,s)
        temp = Group.Units.Task.AssOverlap_{s};
        temp = histc(temp(temp>0),0:length(names)-1);
        tempLabels = names(temp~=0);
        tempColors = Colors(temp~=0,:);
        Hpie= pie(temp,names);
        hp = findobj(Hpie,'Type','patch');
        for iClass= 1:numel(hp)
            set(hp(iClass), 'FaceColor', Colors(iClass,:),...
                            'EdgeColor', [0 0 0],...
                            'FaceAlpha', 1,...
                            'LineWidth', 1);
        end
        title([P.names{s},' units'])
        
end
%% Plot change in sequence probability between epochs
x = cell(3,1);
bins = -0.2:0.01:0.2;
bins_D = 0.01;
colors = [0.2 0.5 0.9 ; ...
          0.9 0.2 0.5 ; ...
          0.9 0.9 0.5] ; ...
alpha = 0.4;
          
for iEpoch = 1:length(epochsD)
    eval(['tempdata = Group.Assem.' epochsD{iEpoch} '.delta_transMatrix;'])
    eval(['Group.Assem.' epochsD{iEpoch} '.delta_transMatrixCumHist = [];'])

    for iRat = 1:length(tempdata)
        temp = tempdata{iRat};
        temp = cumsum(histc(temp,bins)); 
        temp = temp./max(temp);
        if length(tempdata{iRat})> 10
            eval(['Group.Assem.' epochsD{iEpoch} '.delta_transMatrixCumHist(iRat,1:length(bins)) = temp;'])
        else
            eval(['Group.Assem.' epochsD{iEpoch} '.delta_transMatrixCumHist(iRat,1:length(bins)) = NaN;'])
        end
        eval(['Group.Assem.' epochsD{iEpoch} '.delta_transMatrixCumHistBins     = bins;'])
    end
end    

figure('color', 'w'); 
subplot('Position',[0.1 0.15 0.55 0.8]);
    hold on
    plot([0 0],[0 1],':k')
    plot([min(bins) max(bins)],[0.5 0.5],':k')
    ciplot(nanmean(Group.Assem.PreVsTask.delta_transMatrixCumHist)+nansem(Group.Assem.PreVsTask.delta_transMatrixCumHist),...
           nanmean(Group.Assem.PreVsTask.delta_transMatrixCumHist)-nansem(Group.Assem.PreVsTask.delta_transMatrixCumHist),...
           bins,colors(1,:),alpha)
       
    ciplot(nanmean(Group.Assem.TaskVsPost.delta_transMatrixCumHist)+nansem(Group.Assem.TaskVsPost.delta_transMatrixCumHist),...
           nanmean(Group.Assem.TaskVsPost.delta_transMatrixCumHist)-nansem(Group.Assem.TaskVsPost.delta_transMatrixCumHist),...
           bins,colors(2,:),alpha)
       
    ciplot(nanmean(Group.Assem.PreVsPost.delta_transMatrixCumHist)+nansem(Group.Assem.PreVsPost.delta_transMatrixCumHist),...
           nanmean(Group.Assem.PreVsPost.delta_transMatrixCumHist)-nansem(Group.Assem.PreVsPost.delta_transMatrixCumHist),...
           bins,colors(3,:),alpha)    
    plot(bins,nanmean(Group.Assem.PreVsTask.delta_transMatrixCumHist), 'LineWidth',1.5,'color',colors(1,:))
    plot(bins,nanmean(Group.Assem.TaskVsPost.delta_transMatrixCumHist),'LineWidth',1.5,'color',colors(2,:))
    plot(bins,nanmean(Group.Assem.PreVsPost.delta_transMatrixCumHist), 'LineWidth',1.5,'color',colors(3,:))

    axis([min(bins) max(bins) 0 1])
    view(90,-90)
    xlabel('Change in p(Pattern)')
    ylabel('Fraction of rank ordered assembly patterns.')
    
subplot('Position',[0.7 0.15 0.25 0.8]);hold on
    x{1} =nanmean(Group.Assem.PreVsTask.delta_transMatrixCumHist);
    x{2} =nanmean(Group.Assem.TaskVsPost.delta_transMatrixCumHist);
    x{3} =nanmean(Group.Assem.PreVsPost.delta_transMatrixCumHist);
    for s= 1:3
        x{s} = diff(x{s});
    end
    plot(x{1}./max(x{1}),bins(1:end-1)+bins_D,'LineWidth',1.5,'color',colors(1,:))
    plot(x{2}./max(x{2}),bins(1:end-1)+bins_D,'LineWidth',1.5,'color',colors(2,:))
    plot(x{3}./max(x{3}),bins(1:end-1)+bins_D,'LineWidth',1.5,'color',colors(3,:))
    
    legend(epochsD)
    set(gca,'YTick',[])
    % set(gca,'XScale','log')
    xlabel('Normalized distibution')
    plot([0 1],[0 0],':k')
    axis([0 1 min(bins) max(bins)])
    
    
clear iEpoch tempdata temp iRat bins s x colors ans
%% Delay period repeat rate vs sum decoding score
% Average p(sequence) over both left and right outcomes
noRats = length(Group.Assem.Task.DelayTransMatrix);
x= []; y=[];

for iRat = 1:noRats
    SumDecode     = [];
    SumDecodeSort = [];
    SeqSort       = [];
    SeqSortidx    = [];
    Seq = Group.Assem.Task.DelayTransMatrix{iRat};

    Decode = [Group.Assem.TSPeak{1}{iRat} Group.Assem.TSPeak{2}{iRat} Group.Assem.TSPeak{3}{iRat}];
%     Decode = [Group.Assem.TSIntegral{1}{iRat} Group.Assem.TSIntegral{2}{iRat} Group.Assem.TSIntegral{3}{iRat}];

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

x= [x,SeqSort];
y= [y,SumDecodeSort];
y(x>0.3)=[];
x(x>0.3)=[];
end
figure; hold on
    scatter(x,y)
    f = ezfit(x,y,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.85*max(x), 0.8*max(y),strcat('R^2=', sprintf('%3.1g',f.r^2)),'color','r')
    xlabel('Ass(n) > Ass(n+1) sequence probability in delay period')
    ylabel('Sum decoding capacity of Ass(n) + Ass(n+1)')

%% Task vs. Pre-Post change
% (1) Rate
figure('name',' Task rate vs PreVsPost rate changes','NumberTitle','off');
xmax = 0.15;
yminmax = 2;
for s = 1:3
    subplot(1,3,s); hold on
%     plot([0 0],[0 yminmax],':k')
    plot([0 xmax],[0 0],':k')
    tempx = 1./cell2mat(cellfun(@transpose,Group.Assem.Task.mean_repMode{s}','UniformOutput',false)');
    tempy = 1./cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_repMode{s}','UniformOutput',false)');

%     tempx = cell2mat(cellfun(@transpose,Group.Assem.Task.delta_LogStrengthQ3bin{s}','UniformOutput',false)');
%     tempy = cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_LogStrengthQ3bin{s}','UniformOutput',false)');
    
    tempz =  cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
    tempx(tempy==Inf) = [];
    tempz(tempy==Inf) = [];
    tempy(tempy==Inf) = [];
    
    scatter (tempx,tempy,40,tempz,'filled','LineWidth',1,'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.9)
    
    [R_,P_] = corrcoef(tempx,tempy);
    x = [min(tempx) max(tempx)];
    f = ezfit(tempx,tempy,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.1*xmax, 0.8*yminmax,{strcat('R^2=', sprintf('%3.1g',R_(1,2)^2));strcat('p=', sprintf('%3.1g',P_(1,2)))},'color','r')%f.r^2)
  
    
    
    title(P.names{s})
    if s==1
        ylabel('Pre-Post \Delta rate')
    elseif s==2
        xlabel('Task activation rate (Hz)')
    end   
    axis square ; axis([0 xmax -yminmax yminmax])
    caxis([0 10])
end

colormap parula
c = colorbar('Location','manual','Position',[0.92 0.3 0.02 0.2]);
c.Label.String = {'Assembly decoding capacity (Peak t-score)'};

%% (2) Strength
figure('name',' Task rate vs PreVsPost Strength changes','NumberTitle','off');
xmax = 15;
yminmax = 10;
for s = 1:3
    subplot(1,3,s); hold on
%     plot([0 0],[0 yminmax],':k')
    plot([0 xmax],[0 0],':k')
%     tempx = 1./cell2mat(cellfun(@transpose,Group.Assem.Task.mean_repMode{s}','UniformOutput',false)');
%     tempy = 1./cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_repMode{s}','UniformOutput',false)');

    tempx = cell2mat(cellfun(@transpose,Group.Assem.Task.LogStrengthQ3bin{s}','UniformOutput',false)');
    tempy = cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_LogStrengthQ3bin{s}','UniformOutput',false)');
    
    tempz =  cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
    
    scatter (tempx,tempy,40,tempz,'filled','LineWidth',1,'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.9)
    [R_,P_] = corrcoef(tempx,tempy);
    x = [min(tempx) max(tempx)];
    f = ezfit(tempx,tempy,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.1*xmax, 0.8*yminmax,{strcat('R^2=', sprintf('%3.1g',R_(1,2)^2));strcat('p=', sprintf('%3.1g',P_(1,2)))},'color','r')%f.r^2)
  
    
    
    
    title(P.names{s})
    if s==1
        ylabel('Pre to Post Task \Delta strength')
    elseif s==2
        xlabel(' Task activation strength (A.U.)')
    end   
    axis square ; axis([0 xmax -yminmax yminmax])
    caxis([0 10])
end

colormap parula
c = colorbar('Location','manual','Position',[0.92 0.3 0.02 0.2]);
c.Label.String = {'Assembly decoding capacity (Peak t-score)'};

%% (3) Task rate vs. pre-to-post strength change

figure('name',' Task rate vs PreVsPost Strength changes','NumberTitle','off');
xmax = 0.15;
yminmax = 10;
for s = 1:3
    subplot(1,3,s); hold on
%     plot([0 0],[0 yminmax],':k')
    plot([0 xmax],[0 0],':k')
    tempx = 1./cell2mat(cellfun(@transpose,Group.Assem.Task.mean_repMode{s}','UniformOutput',false)');
%     tempy = 1./cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_repMode{s}','UniformOutput',false)');

%     tempx = cell2mat(cellfun(@transpose,Group.Assem.Task.LogStrengthQ3bin{s}','UniformOutput',false)');
    tempy = cell2mat(cellfun(@transpose,Group.Assem.PreVsPost.delta_LogStrengthQ3bin{s}','UniformOutput',false)');
    
    tempz =  cell2mat(cellfun(@transpose,Group.Assem.TSPeak{s}','UniformOutput',false)');
    tempx(isnan(tempx)) = [];
    tempy(isnan(tempy)) = [];
    tempz(isnan(tempz)) = [];
    
    scatter (tempx,tempy,40,tempz,'filled','LineWidth',1,'MarkerEdgeColor',[0 0 0],'MarkerFaceAlpha',0.9)
    [R_,P_] = corrcoef(tempx,tempy);
    x = [min(tempx) max(tempx)];
    f = ezfit(tempx,tempy,'a*x+b');
    plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color','r')
    text(0.1*xmax, 0.8*yminmax,{strcat('R^2=', sprintf('%3.1g',R_(1,2)^2));strcat('p=', sprintf('%3.1g',P_(1,2)))},'color','r')%f.r^2)
  
    
    title(P.names{s})
    if s==1
        ylabel('Pre to Post Task \Delta strength')
    elseif s==2
        xlabel(' Task activation rate (Hz)')
    end   
    axis square ; axis([0 xmax -yminmax yminmax])
    caxis([0 10])
end

colormap parula
c = colorbar('Location','manual','Position',[0.92 0.3 0.02 0.2]);
c.Label.String = {'Assembly decoding capacity (Peak t-score)'};
    
%% Pre-vs-Post change vs. encoding type
% Load class assignments
load([  pat{1} 'UnitAssemClassified.mat'])
for iList = 1:length(Ignore)
    listout= listout(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
           {listout.name_},repmat({Ignore(iList)},1,length(listout)),'UniformOutput',false)));
end

% Pre-vs-Post change vs. encoding type
for s = 1:3
    Group.Assem.Class{s}=[]; 
    for iRat = 1:noRats
        Group.Assem.Class{s} = [Group.Assem.Class{s};listout(iRat).AssemClass{s}];
    end
end

for s = 1:3
    temp1 = cell2mat(Group.Assem.PreVsPost.delta_repMode{s}')';
    temp2 = cell2mat(Group.Assem.PreVsPost.delta_LogStrengthQ3bin{s}')';
    temp3 = cell2mat(Group.Assem.Task.mean_repMode{s}')';
    temp4 = cell2mat(Group.Assem.Task.LogStrengthQ3bin{s}')';
    for iClass = 1:size(kmeans_.Centroids,1)
        Group.Assem.PreVsPost.delta_repModeByClass{s}{iClass}          =  temp1(find(Group.Assem.Class{s}==iClass));
        Group.Assem.PreVsPost.delta_LogStrengthQ3binByClass{s}{iClass} =  temp2(find(Group.Assem.Class{s}==iClass));
        Group.Assem.Task.mean_repModeByClass{s}{iClass}                =  temp3(find(Group.Assem.Class{s}==iClass));
        Group.Assem.Task.LogStrengthQ3binByClass{s}{iClass}            =  temp4(find(Group.Assem.Class{s}==iClass));
        
    end
    
    Group.Assem.PreVsPost.delta_repModeByClassMean{s} = cellfun(@nanmean,Group.Assem.PreVsPost.delta_repModeByClass{s});
    Group.Assem.PreVsPost.delta_repModeByClassSEM{s}  = cellfun(@nansem,Group.Assem.PreVsPost.delta_repModeByClass{s});

    Group.Assem.PreVsPost.delta_LogStrengthQ3binByClassMean{s} = cellfun(@nanmean,Group.Assem.PreVsPost.delta_LogStrengthQ3binByClass{s});
    Group.Assem.PreVsPost.delta_LogStrengthQ3binByClassSEM{s}  = cellfun(@nansem,Group.Assem.PreVsPost.delta_LogStrengthQ3binByClass{s});
    
    Group.Assem.Task.mean_repModeByClassMean{s} = cellfun(@nanmean,Group.Assem.Task.mean_repModeByClass{s});
    Group.Assem.Task.mean_repModeByClassSEM{s}  = cellfun(@nansem,Group.Assem.Task.mean_repModeByClass{s});
    
    Group.Assem.Task.LogStrengthQ3binByClassMean{s} = cellfun(@nanmean,Group.Assem.Task.LogStrengthQ3binByClass{s});
    Group.Assem.Task.LogStrengthQ3binByClassSEM{s}  = cellfun(@nansem,Group.Assem.Task.LogStrengthQ3binByClass{s});
end

ClassNames = {'Pre-Sample',...
              'Sample',...
              'Early Delay',...
              'Late Delay',...
              'Sample/Choice',...
              'Post-Choice'};
ClassColors = {[0.7 0.9 0.9],...
               [0.1 0.9 0.1],...
               [0.9 0.7 0.9],...
               [0.9 0.7 0.9],...
               [0.9 0.6 0.1],...
               [0.9 0.1 0.1]};
           
classes = 1:size(kmeans_.Centroids,1);

clear temp1 temp2 temp3 temp4
%% Plot change in activation rate
yminmax = 35;
figure('color','w')
for s=1:3
   subplot(3,1,s); hold on 
   for iClass = 1:size(kmeans_.Centroids,1)
        B = bar(iClass,Group.Assem.PreVsPost.delta_repModeByClassMean{s}(iClass));
        set(B,'FaceColor','w',...
              'EdgeColor',ClassColors{iClass},...
              'LineWidth',2);     
          
        errorbar(iClass,Group.Assem.PreVsPost.delta_repModeByClassMean{s}(iClass),...
                         Group.Assem.PreVsPost.delta_repModeByClassSEM{s}(iClass),...
                         'Color',ClassColors{iClass},...
                         'LineWidth',2)
   end
   axis tight
   ylim_ = max(abs(get(gca,'ylim')));
   ylim([-ylim_ ylim_])
%    axis([0 size(kmeans_.Centroids,1)+1 -yminmax yminmax])
   set(gca,'xtick',classes,'xTicklabel',ClassNames)
   title(P.names{s})
   set(gca,'xtick',[])
    if s==2
        ylabel('Pre to Post Task \Delta activation rate')
    elseif s==3
        xlabel('Assembly category')
        set(gca,'xtick',classes,'xTicklabel',ClassNames)
    end   
end
%% Plot change in activation strength
yminmax = 2;
figure('color','w')
for s=1:3
   subplot(3,1,s); hold on 
   for iClass = 1:size(kmeans_.Centroids,1)
        B = bar(iClass,Group.Assem.PreVsPost.delta_LogStrengthQ3binByClassMean{s}(iClass));
        set(B,'FaceColor','w',...
              'EdgeColor',ClassColors{iClass},...
              'LineWidth',2);     
        errorbar(iClass,Group.Assem.PreVsPost.delta_LogStrengthQ3binByClassMean{s}(iClass),...
                         Group.Assem.PreVsPost.delta_LogStrengthQ3binByClassSEM{s}(iClass),...
                         'Color',ClassColors{iClass},...
                         'LineWidth',2)
   end
      axis tight
      ylim_ = max(abs(get(gca,'ylim')));
   ylim([-ylim_ ylim_])
% %    axis([0 size(kmeans_.Centroids,1)+1 -yminmax yminmax])
   set(gca,'xtick',classes,'xTicklabel',ClassNames)
   title(P.names{s})
   set(gca,'xtick',[])
    if s==2
        ylabel('Pre to Post Task \Delta activation strength')
    elseif s==3
        xlabel('Assembly category')
                set(gca,'xtick',classes,'xTicklabel',ClassNames)

    end   
    
end
%% Task activation strength
yminmax = 2;
figure('color','w')
for s=1:3
   subplot(3,1,s); hold on 
   for iClass = 1:size(kmeans_.Centroids,1)
        B = bar(iClass,Group.Assem.Task.LogStrengthQ3binByClassMean{s}(iClass));
        set(B,'FaceColor','w',...
              'EdgeColor',ClassColors{iClass},...
              'LineWidth',2);     
        errorbar(iClass,Group.Assem.Task.LogStrengthQ3binByClassMean{s}(iClass),...
                         Group.Assem.Task.LogStrengthQ3binByClassSEM{s}(iClass),...
                         'Color',ClassColors{iClass},...
                         'LineWidth',2)
   end
      axis tight
      ylim_ = max(abs(get(gca,'ylim')));
   ylim([-ylim_ ylim_])
% %    axis([0 size(kmeans_.Centroids,1)+1 -yminmax yminmax])
   set(gca,'xtick',classes,'xTicklabel',ClassNames)
   title(P.names{s})
   set(gca,'xtick',[])
    if s==2
        ylabel('Task activation strength')
    elseif s==3
        xlabel('Assembly category')
                set(gca,'xtick',classes,'xTicklabel',ClassNames)

    end   
    
end

%% Task activation rate
yminmax = 35;
figure('color','w')
for s=1:3
   subplot(3,1,s); hold on 
   for iClass = 1:size(kmeans_.Centroids,1)
        B = bar(iClass,Group.Assem.Task.mean_repModeByClassMean{s}(iClass));
        set(B,'FaceColor','w',...
              'EdgeColor',ClassColors{iClass},...
              'LineWidth',2);     
          
        errorbar(iClass,Group.Assem.Task.mean_repModeByClassMean{s}(iClass),...
                         Group.Assem.Task.mean_repModeByClassSEM{s}(iClass),...
                         'Color',ClassColors{iClass},...
                         'LineWidth',2)
   end
   axis tight
   ylim_ = max(abs(get(gca,'ylim')));
   ylim([-ylim_ ylim_])
%    axis([0 size(kmeans_.Centroids,1)+1 -yminmax yminmax])
   set(gca,'xtick',classes,'xTicklabel',ClassNames)
   title(P.names{s})
   set(gca,'xtick',[])
    if s==2
        ylabel('Mean Task Activation rate')
    elseif s==3
        xlabel('Assembly category')
        set(gca,'xtick',classes,'xTicklabel',ClassNames)
    end   
end

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
%          'OnufryLONG2'};

Ignore = {};
for iList = 1:length(Ignore)
    Rats = Rats(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
           {Rats.name},repmat({Ignore(iList)},1,length(Rats)),'UniformOutput',false)));
end
clear iList
%% Process decoding results
Group =struct;
for iRat= 1:length(Rats)
    load([pat{1} Rats(iRat).name],'P','D','Dsc','FAtrials','FRtrials');
    disp(['Working on ' Rats(iRat).name '...'])
    for s=1:length(P.names)       
        %% Collate Assembly activation scores
        % Correct trials %%%%%%%%%%%%%%%%%%
        Group.Assem.FSC_{s}{iRat,:}   = FAtrials.FSC_{s};
        Group.Assem.FSCmeanL{s}{iRat} = zeros(size(FAtrials.FSC_{s}{1}));
        Group.Assem.FSCmeanR{s}{iRat} = zeros(size(FAtrials.FSC_{s}{1}));
        Group.Assem.FSCmeanBoth{s}{iRat} = zeros(size(FAtrials.FSC_{s}{1}));
        Llist = find(FRtrials.evt0==1); Rlist = find(FRtrials.evt0==2);
        for i = 1:length(Llist)
             Group.Assem.FSCmeanL{s}{iRat} =  Group.Assem.FSCmeanL{s}{iRat} + FAtrials.FSC_{s}{Llist(i)};
        end
        for i = 1:length(Rlist)
             Group.Assem.FSCmeanR{s}{iRat} =  Group.Assem.FSCmeanR{s}{iRat} + FAtrials.FSC_{s}{Rlist(i)};
        end
        Group.Assem.FSCmeanL{s}{iRat} = Group.Assem.FSCmeanL{s}{iRat}./length(Llist);    
        Group.Assem.FSCmeanR{s}{iRat} = Group.Assem.FSCmeanR{s}{iRat}./length(Rlist);  
        Group.Assem.FSCmeanBoth{s}{iRat} = (Group.Assem.FSCmeanL{s}{iRat}+Group.Assem.FSCmeanR{s}{iRat})./2;
        
        % Error trials %%%%%%%%%%%%%%%%%%
        Group.Assem.Err.FSC_{s}{iRat,:}      = FAtrials.Err.FSC_{s};
        Group.Assem.Err.FSCmeanL{s}{iRat}    = zeros(size(FAtrials.Err.FSC_{s}{1}));
        Group.Assem.Err.FSCmeanR{s}{iRat}    = zeros(size(FAtrials.Err.FSC_{s}{1}));
        Group.Assem.Err.FSCmeanBoth{s}{iRat} = zeros(size(FAtrials.Err.FSC_{s}{1}));
        Llist = find(FRtrials.Err.evt0==1); Rlist = find(FRtrials.Err.evt0==2);
        for i = 1:length(Llist)
             Group.Assem.Err.FSCmeanL{s}{iRat} =  Group.Assem.Err.FSCmeanL{s}{iRat} + FAtrials.Err.FSC_{s}{Llist(i)};
        end
        for i = 1:length(Rlist)
             Group.Assem.Err.FSCmeanR{s}{iRat} =  Group.Assem.Err.FSCmeanR{s}{iRat} + FAtrials.Err.FSC_{s}{Rlist(i)};
        end
        Group.Assem.Err.FSCmeanL{s}{iRat} = Group.Assem.Err.FSCmeanL{s}{iRat}./length(Llist);    
        Group.Assem.Err.FSCmeanR{s}{iRat} = Group.Assem.Err.FSCmeanR{s}{iRat}./length(Rlist);  
        Group.Assem.Err.FSCmeanBoth{s}{iRat} = (Group.Assem.Err.FSCmeanL{s}{iRat}+Group.Assem.Err.FSCmeanR{s}{iRat})./2;
        %% Collate LR decoding results
        % Correct trials %%%%%%%%%%%%%%%%%%
        Group.Units.Ft2ciH0{s} (iRat,:)      = D.units.Ft2ciH0{s};
        Group.Units.Ft2x{s}    (iRat,:)      = D.units.Ft2x{s};
        Group.Units.TSPeak{s}  {iRat,:}      = nanmax(D.units.TS{s},1);
        Group.Units.TSIntegral{s}{iRat,:}    = nansum(D.units.TS{s},1);
        Group.Units.LR_TS{s}{iRat,:}         = D.units.TS{s};
        
        Group.MemberUnits.LR_TS{s}{iRat,:}         = D.MemberUnitsCollapsed.TS{s};
        if ~isempty(D.MemberUnitsCollapsed.Ft2x{s})
            Group.MemberUnits.LR_Ft2x{s}(iRat,:)        = D.MemberUnitsCollapsed.Ft2x{s};
            Group.MemberUnits.LR_Ft2ciH0{s}(iRat,:)     = D.MemberUnitsCollapsed.Ft2ciH0{s};
        else
            Group.MemberUnits.LR_Ft2x{s}(iRat,:)        = nan(P.Ltr*2,1);
            Group.MemberUnits.LR_Ft2ciH0{s}(iRat,:)     = nan(P.Ltr*2,1);
        end
        %Group.NonMemberUnits.LR_TS{s}{iRat,:}      = D.NonMemberUnits.TS{s};
        Group.NonMemberUnits.LR_Ft2x{s}(iRat,:)     = D.NonMemberUnits.Ft2x{s};
        Group.NonMemberUnits.LR_Ft2ciH0{s}(iRat,:)  = D.NonMemberUnits.Ft2ciH0{s};
              
        Group.Assem.Ft2ciH0{s} (iRat,:)      = D.Assem.Ft2ciH0{s};
        Group.Assem.Ft2x{s}    (iRat,:)      = D.Assem.Ft2x{s};
        Group.Assem.TSPeak{s}  {iRat,:}      = D.Assem.TSPeak{s};
        Group.Assem.TSIntegral{s}{iRat,:}    = D.Assem.TSIntegral{s};
        Group.Assem.LR_TS{s}{iRat,:}           = D.Assem.TS{s};

        % Group.AssemAdvantage{s}(iRat,:)      = (D.Assem.Ft2x{s}-D.units.Ft2x{s})./(D.Assem.Ft2x{s}+D.units.Ft2x{s});
        Group.AssemAdvantage{s}(iRat,:)      = D.Assem.Ft2x{s}-D.units.Ft2x{s};
        
        
        % Error trials %%%%%%%%%%%%%%%%%%
        Group.Units.Err.Ft2ciH0{s} (iRat,:)      = D.units.Err.Ft2ciH0{s};
        Group.Units.Err.Ft2x{s}    (iRat,:)      = D.units.Err.Ft2x{s};
        Group.Units.Err.TSPeak{s}  {iRat,:}      = nanmax(D.units.Err.TS{s},1);
        Group.Units.Err.TSIntegral{s}{iRat,:}    = nansum(D.units.Err.TS{s},1);
        Group.Units.Err.LR_TS{s}{iRat,:}         = D.units.Err.TS{s};
        
        
              
        Group.Assem.Err.Ft2ciH0{s} (iRat,:)      = D.Assem.Err.Ft2ciH0{s};
        Group.Assem.Err.Ft2x{s}    (iRat,:)      = D.Assem.Err.Ft2x{s};
        Group.Assem.Err.TSPeak{s}  {iRat,:}      = D.Assem.Err.TSPeak{s};
        Group.Assem.Err.TSIntegral{s}{iRat,:}    = D.Assem.Err.TSIntegral{s};
        Group.Assem.Err.LR_TS{s}{iRat,:}         = D.Assem.Err.TS{s};

%         % Group.AssemAdvantage{s}(iRat,:)      = (D.Assem.Ft2x{s}-D.units.Ft2x{s})./(D.Assem.Ft2x{s}+D.units.Ft2x{s});
%         Group.AssemAdvantage{s}(iRat,:)      = D.Assem.Ft2x{s}-D.units.Ft2x{s};
        %% Collate SC decoding results
        % Correct trials %%%%%%%%%%%%%%%%%%
        Group.Units.SC_Ft2ciH0{s} (iRat,:)      = Dsc.units.Ft2ciH0{s};
        Group.Units.SC_Ft2x{s}    (iRat,:)      = Dsc.units.Ft2x{s};
        Group.Units.SC_TSPeak{s}  {iRat,:}      = nanmax(Dsc.units.TS{s},1);
        Group.Units.SC_TSIntegral{s}{iRat,:}    = nansum(Dsc.units.TS{s},1);
        
        Group.MemberUnits.SC_TS{s}{iRat,:}        = Dsc.MemberUnitsCollapsed.TS{s};
        if ~isempty(Dsc.MemberUnitsCollapsed.Ft2x{s})
            Group.MemberUnits.SC_Ft2x{s}(iRat,:)        = Dsc.MemberUnitsCollapsed.Ft2x{s};
            Group.MemberUnits.SC_Ft2ciH0{s}(iRat,:)     = Dsc.MemberUnitsCollapsed.Ft2ciH0{s};
        else
            Group.MemberUnits.SC_Ft2x{s}(iRat,:)        = nan(P.Ltr_SC,1);
            Group.MemberUnits.SC_Ft2ciH0{s}(iRat,:)     = nan(P.Ltr_SC,1);
        end
        
        % Group.NonMemberUnits.SC_TS{s}{iRat,:}     = Dsc.NonMemberUnits.TS{s};
        Group.NonMemberUnits.SC_Ft2x{s}(iRat,:)     = Dsc.NonMemberUnits.Ft2x{s};
        Group.NonMemberUnits.SC_Ft2ciH0{s}(iRat,:)  = Dsc.NonMemberUnits.Ft2ciH0{s};
        
        Group.Assem.SC_Ft2ciH0{s} (iRat,:)      = Dsc.Assem.Ft2ciH0{s};
        Group.Assem.SC_Ft2x{s}    (iRat,:)      = Dsc.Assem.Ft2x{s};
        Group.Assem.SC_TSPeak{s}  {iRat,:}      = Dsc.Assem.TSPeak{s};
        Group.Assem.SC_TSIntegral{s}{iRat,:}    = Dsc.Assem.TSIntegral{s};
        Group.Assem.SC_TS{s}{iRat,:}            = Dsc.Assem.TS{s};
        
        % Group.AssemAdvantageSC{s}(iRat,:)      = (D.Assem.Ft2x{s}-D.units.Ft2x{s})./(D.Assem.Ft2x{s}+D.units.Ft2x{s});
        Group.AssemAdvantageSC{s}(iRat,:)      = Dsc.Assem.Ft2x{s}-Dsc.units.Ft2x{s};
        
        % Error trials %%%%%%%%%%%%%%%%%%
        Group.Units.Err.SC_Ft2ciH0{s} (iRat,:)      = Dsc.units.Err.Ft2ciH0{s};
        Group.Units.Err.SC_Ft2x{s}    (iRat,:)      = Dsc.units.Err.Ft2x{s};
        Group.Units.Err.SC_TSPeak{s}  {iRat,:}      = nanmax(Dsc.units.Err.TS{s},1);
        Group.Units.Err.SC_TSIntegral{s}{iRat,:}    = nansum(Dsc.units.Err.TS{s},1);
        
        Group.Assem.Err.SC_Ft2ciH0{s} (iRat,:)      = Dsc.Assem.Err.Ft2ciH0{s};
        Group.Assem.Err.SC_Ft2x{s}    (iRat,:)      = Dsc.Assem.Err.Ft2x{s};
        Group.Assem.Err.SC_TSPeak{s}  {iRat,:}      = Dsc.Assem.Err.TSPeak{s};
        Group.Assem.Err.SC_TSIntegral{s}{iRat,:}    = Dsc.Assem.Err.TSIntegral{s};
        Group.Assem.Err.SC_TS{s}{iRat,:}            = Dsc.Assem.Err.TS{s};
        %% Find significant decoders 
        % Correct trials %%%%%%%%%%%%%%%%%%
        Group.Units.TSsigLR{s}{iRat,:} = (tpdf(D.units.TS{s},numel(FRtrials.evt0)-2))<0.05;
        Group.Units.TSsigSC{s}{iRat,:} = (tpdf(Dsc.units.TS{s},2*numel(FRtrials.evt0)-2))<0.05; 
        
        Group.MemberUnits.TSsigLR{s}{iRat,:} = (tpdf(D.MemberUnitsCollapsed.TS{s},numel(FRtrials.evt0)-2))<0.05;
        Group.MemberUnits.TSsigSC{s}{iRat,:} = (tpdf(Dsc.MemberUnitsCollapsed.TS{s},2*numel(FRtrials.evt0)-2))<0.05; 
        
        if ~isempty(FAtrials.units{s})
            Group.Assem.TSsigLR{s}{iRat,:} = (tpdf(D.Assem.TS{s},numel(FRtrials.evt0)-2))<0.05;
            Group.Assem.TSsigSC{s}{iRat,:} = (tpdf(Dsc.Assem.TS{s},2*numel(FRtrials.evt0)-2))<0.05;
        else
            Group.Assem.TSsigLR{s}{iRat,:} = [];
            Group.Assem.TSsigSC{s}{iRat,:} = [];
        end
        % Error trials %%%%%%%%%%%%%%%%%%
        Group.Units.Err.TSsigLR{s}{iRat,:} = (tpdf(D.units.Err.TS{s},numel(FRtrials.Err.evt0)-2))<0.05;
        Group.Units.Err.TSsigSC{s}{iRat,:} = (tpdf(Dsc.units.Err.TS{s},2*numel(FRtrials.Err.evt0)-2))<0.05; 
                
        if ~isempty(FAtrials.Err.units{s})
            Group.Assem.Err.TSsigLR{s}{iRat,:} = (tpdf(D.Assem.Err.TS{s},numel(FRtrials.Err.evt0)-2))<0.05;
            Group.Assem.Err.TSsigSC{s}{iRat,:} = (tpdf(Dsc.Assem.Err.TS{s},2*numel(FRtrials.Err.evt0)-2))<0.05;
        else
            Group.Assem.Err.TSsigLR{s}{iRat,:} = [];
            Group.Assem.Err.TSsigSC{s}{iRat,:} = [];
        end
    end 
end
clear Llist Rlist D Dsc FAtrials FRtrials iRat s i
%% Collapse imported data - Correct trials

% Collapse all factor scores
 for s=1:length(P.names)
      Group.Assem.FSCmeanBothCollapse{s}      = cell2mat(Group.Assem.FSCmeanBoth{s});
      Group.Assem.FSCmeanBothCollapse{s}      = Group.Assem.FSCmeanBothCollapse{s}-mean(Group.Assem.FSCmeanBothCollapse{s}(1:20,:));
      Group.Assem.FSCmeanBothCollapseMean{s} = nanmean(Group.Assem.FSCmeanBothCollapse{s},2);
      Group.Assem.FSCmeanBothCollapseSEM{s}  = nansem(Group.Assem.FSCmeanBothCollapse{s},2);
 end
% Collapse LR decoding
Group.Units.Ft2ciH0_Mean  = cellfun(@nanmean,Group.Units.Ft2ciH0,'UniformOutput',false);
Group.Units.Ft2x_Mean     = cellfun(@nanmean,Group.Units.Ft2x,'UniformOutput',false);
Group.Units.Ft2x_SEM      = cellfun(@nansem,Group.Units.Ft2x,'UniformOutput',false);

Group.MemberUnits.LR_Ft2ciH0_mean = cellfun(@nanmean,Group.MemberUnits.LR_Ft2ciH0,'UniformOutput',false);
Group.MemberUnits.LR_Ft2x_mean = cellfun(@nanmean,Group.MemberUnits.LR_Ft2x,'UniformOutput',false);
Group.MemberUnits.LR_Ft2x_SEM  = cellfun(@nansem,Group.MemberUnits.LR_Ft2x,'UniformOutput',false);

Group.NonMemberUnits.LR_Ft2ciH0_mean = cellfun(@nanmean,Group.NonMemberUnits.LR_Ft2ciH0,'UniformOutput',false);
Group.NonMemberUnits.LR_Ft2x_mean = cellfun(@nanmean,Group.NonMemberUnits.LR_Ft2x,'UniformOutput',false);
Group.NonMemberUnits.LR_Ft2x_SEM  = cellfun(@nansem,Group.NonMemberUnits.LR_Ft2x,'UniformOutput',false);

Group.Assem.Ft2ciH0_Mean  = cellfun(@nanmean,Group.Assem.Ft2ciH0,'UniformOutput',false);
Group.Assem.Ft2x_Mean     = cellfun(@nanmean,Group.Assem.Ft2x,'UniformOutput',false);
Group.Assem.Ft2x_SEM      = cellfun(@nansem,Group.Assem.Ft2x,'UniformOutput',false);

Group.AssemAdvantage_Mean = cellfun(@nanmean,Group.AssemAdvantage,'UniformOutput',false);
Group.AssemAdvantage_SEM  = cellfun(@nansem,Group.AssemAdvantage,'UniformOutput',false);

% for s= 1:3
%     Group.MemberUnits.LR_TS_mean{s} = nanmean(cell2mat(Group.MemberUnits.LR_TS{s}'));
%     Group.MemberUnits.LR_TS_SEM{s}  = nansem(cell2mat(Group.MemberUnits.LR_TS{s}'));
%     Group.NonMemberUnits.LR_TS_mean{s} = nanmean(cell2mat(Group.NonMemberUnits.LR_TS{s}'));
%     Group.NonMemberUnits.LR_TS_SEM{s}  = nansem(cell2mat(Group.NonMemberUnits.LR_TS{s}'));
% end

% Collapse SC decoding
Group.Units.SC_Ft2ciH0_Mean  = cellfun(@nanmean,Group.Units.SC_Ft2ciH0,'UniformOutput',false);
Group.Units.SC_Ft2x_Mean     = cellfun(@nanmean,Group.Units.SC_Ft2x,'UniformOutput',false);
Group.Units.SC_Ft2x_SEM      = cellfun(@nansem,Group.Units.SC_Ft2x,'UniformOutput',false);

Group.MemberUnits.SC_Ft2ciH0_mean = cellfun(@nanmean,Group.MemberUnits.SC_Ft2ciH0,'UniformOutput',false);
Group.MemberUnits.SC_Ft2x_mean    = cellfun(@nanmean,Group.MemberUnits.SC_Ft2x,'UniformOutput',false);
Group.MemberUnits.SC_Ft2x_SEM     = cellfun(@nansem,Group.MemberUnits.SC_Ft2x,'UniformOutput',false);

Group.NonMemberUnits.SC_Ft2ciH0_mean = cellfun(@nanmean,Group.NonMemberUnits.SC_Ft2ciH0,'UniformOutput',false);
Group.NonMemberUnits.SC_Ft2x_mean    = cellfun(@nanmean,Group.NonMemberUnits.SC_Ft2x,'UniformOutput',false);
Group.NonMemberUnits.SC_Ft2x_SEM     = cellfun(@nansem,Group.NonMemberUnits.SC_Ft2x,'UniformOutput',false);

Group.Assem.SC_Ft2ciH0_Mean  = cellfun(@nanmean,Group.Assem.SC_Ft2ciH0,'UniformOutput',false);
Group.Assem.SC_Ft2x_Mean     = cellfun(@nanmean,Group.Assem.SC_Ft2x,'UniformOutput',false);
Group.Assem.SC_Ft2x_SEM      = cellfun(@nansem,Group.Assem.SC_Ft2x,'UniformOutput',false);
Group.AssemAdvantageSC_Mean  = cellfun(@nanmean,Group.AssemAdvantageSC,'UniformOutput',false);
Group.AssemAdvantageSC_SEM   = cellfun(@nansem,Group.AssemAdvantageSC,'UniformOutput',false);
% for s= 1:3
%     Group.MemberUnits.SC_TS_mean{s} = nanmean(cell2mat(Group.MemberUnits.SC_TS{s}'));
%     Group.MemberUnits.SC_TS_SEM{s}  = nansem(cell2mat(Group.MemberUnits.SC_TS{s}'));
%     Group.NonMemberUnits.SC_TS_mean{s} = nanmean(cell2mat(Group.NonMemberUnits.SC_TS{s}'));
%     Group.NonMemberUnits.SC_TS_SEM{s}  = nansem(cell2mat(Group.NonMemberUnits.SC_TS{s}'));
% end
clear iRat s D FAtrials FRtrials i 
%% Collapse imported data - Error trials

% Collapse all factor scores
 for s=1:length(P.names)
      Group.Assem.Err.FSCmeanBothCollapse{s}      = cell2mat(Group.Assem.Err.FSCmeanBoth{s});
      Group.Assem.Err.FSCmeanBothCollapse{s}      = Group.Assem.Err.FSCmeanBothCollapse{s}-mean(Group.Assem.Err.FSCmeanBothCollapse{s}(1:20,:));
      Group.Assem.Err.FSCmeanBothCollapseMean{s} = nanmean(Group.Assem.Err.FSCmeanBothCollapse{s},2);
      Group.Assem.Err.FSCmeanBothCollapseSEM{s}  = nansem(Group.Assem.Err.FSCmeanBothCollapse{s},2);
 end
% Collapse LR decoding
Group.Units.Err.Ft2ciH0_Mean  = cellfun(@nanmean,Group.Units.Err.Ft2ciH0,'UniformOutput',false);
Group.Units.Err.Ft2x_Mean     = cellfun(@nanmean,Group.Units.Err.Ft2x,'UniformOutput',false);
Group.Units.Err.Ft2x_SEM      = cellfun(@nansem,Group.Units.Err.Ft2x,'UniformOutput',false);

Group.Assem.Err.Ft2ciH0_Mean  = cellfun(@nanmean,Group.Assem.Err.Ft2ciH0,'UniformOutput',false);
Group.Assem.Err.Ft2x_Mean     = cellfun(@nanmean,Group.Assem.Err.Ft2x,'UniformOutput',false);
Group.Assem.Err.Ft2x_SEM      = cellfun(@nansem,Group.Assem.Err.Ft2x,'UniformOutput',false);

% Group.AssemAdvantage_Mean = cellfun(@nanmean,Group.Err.AssemAdvantage,'UniformOutput',false);
% Group.AssemAdvantage_SEM  = cellfun(@nansem,Group.Err.AssemAdvantage,'UniformOutput',false);

% Collapse SC decoding
Group.Units.Err.SC_Ft2ciH0_Mean  = cellfun(@nanmean,Group.Units.Err.SC_Ft2ciH0,'UniformOutput',false);
Group.Units.Err.SC_Ft2x_Mean     = cellfun(@nanmean,Group.Units.Err.SC_Ft2x,'UniformOutput',false);
Group.Units.Err.SC_Ft2x_SEM      = cellfun(@nansem,Group.Units.Err.SC_Ft2x,'UniformOutput',false);

Group.Assem.Err.SC_Ft2ciH0_Mean  = cellfun(@nanmean,Group.Assem.Err.SC_Ft2ciH0,'UniformOutput',false);
Group.Assem.Err.SC_Ft2x_Mean     = cellfun(@nanmean,Group.Assem.Err.SC_Ft2x,'UniformOutput',false);
Group.Assem.Err.SC_Ft2x_SEM      = cellfun(@nansem,Group.Assem.Err.SC_Ft2x,'UniformOutput',false);

clear iRat s D FAtrials FRtrials i 
%% Plot all factor scores
figure('name','Grand average factor scores','color','w');hold on

% Plot decoding
maxH = 10;               % yMax of decoding plots
minH = -5;
LeverPressWidth = 5;    % Width of lever press bar (s)
MinMaxAdvantage = 10;   % yMin and yMax of assembly advantage plot 
SAMPLEln  = [0.89 0.95 0.89];
CHOICEln  = [0.95 0.89 0.89];
col_ = {'r','g','b'};
alpha = 0.1;


    
for s=1:3
    subplot(1,3,s);hold on
    title(P.names{s})
    plot([10 10],[minH maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[minH maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
%     ciplot(Group.Assem.FSCmeanBothCollapseMean{s}+Group.Assem.FSCmeanBothCollapseSEM{s},...
%            Group.Assem.FSCmeanBothCollapseMean{s}-Group.Assem.FSCmeanBothCollapseSEM{s},...
%            (1:P.Ltr*2)*P.bw,col_{s},alpha) 
       
    plot((1:P.Ltr*2)*P.bw,Group.Assem.FSCmeanBothCollapseMean{s},col_{s},'LineWidth',1.5)
    plot((1:P.Ltr*2)*P.bw,Group.Assem.FSCmeanBothCollapse{s},col_{s})



    plot([15 15],[minH maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[minH maxH],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')

%     set(gca,'XTick',[])
    
    axis([0 Inf minH Inf ])
    if s==1
        ylabel({'Assembly activation score'})

    elseif s==2
        xlabel('Time (s)')
    end
end
%% Plot all factor scores mean
figure('name','Grand average factor scores','color','w');hold on

maxH = 2;               % yMax of decoding plots
minH = -1;
LeverPressWidth = 5;    % Width of lever press bar (s)
MinMaxAdvantage = 10;   % yMin and yMax of assembly advantage plot 
SAMPLEln  = [0.89 0.95 0.89];
CHOICEln  = [0.95 0.89 0.89];
col_ = {'r','g','b'};
alpha = 0.1;
title('Grand mean Assembly activation scores')
plot([10 10],[minH maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
plot([20 20],[minH maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    
for s=1:3
    ciplot(Group.Assem.FSCmeanBothCollapseMean{s}+Group.Assem.FSCmeanBothCollapseSEM{s},...
           Group.Assem.FSCmeanBothCollapseMean{s}-Group.Assem.FSCmeanBothCollapseSEM{s},...
           (1:P.Ltr*2)*P.bw,col_{s},alpha) 
       
    plot((1:P.Ltr*2)*P.bw,Group.Assem.FSCmeanBothCollapseMean{s},col_{s})
%     plot((1:P.Ltr*2)*P.bw,Group.Assem.FSCmeanBothCollapse{s},col_{s})

end
plot([15 15],[minH maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
plot([15 15],[minH maxH],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
axis([0 Inf minH maxH ])
ylabel({'Assembly activation score'})
xlabel('Time (s)')
text(25,-0.2,'PFC','color','r')
text(25,-0.3,'HP','color','g')
text(25,-0.4,'Joint','color','b')
%% Plot decoding of PFC/HP units
clear data_norm data_mean data_SEM

figure('name','Group LR decoding','color','w')
textAlign_ = 'middle';
normFlag = true;
normWin = [0 7];
LeverPressWidth = 2; % Width of lever press bar (s)
NBS = 1000;           % No draws for permutation tests
pValue = 0.05;       % Bonferroni corrected inside the function


if normFlag
    maxH = 10;              % yMax of decoding plots
else
    maxH = 10;              % yMax of decoding plots
end


SAMPLEln  = [0.6 0.95 0.6];
CHOICEln  = [0.95 0.6 0.6];
alpha = 0.8;
x_temp = (1:P.Ltr*2)*P.bw;
bp=find(x_temp >=normWin(1) & x_temp <=normWin(2));
col_ = {'k','b','r'}; %[0.9 0.6 0.3]};
label_str_ = {'',... % All single units
              'Error trials',...
              'Correct trials'};
label_pos_ = [0,9,10];
figure;hold on
title('Units decoding compared')


    plot([10 10],[0 maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    


    data_={nan*Group.Units.Ft2x{1},...
               Group.Units.Ft2x{1},...  
               Group.Units.Ft2x{2}};
  
    for a =1:length(data_)
       if normFlag
           B = nanmean(data_{a}(:,bp),2);
           data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
           data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
           data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));
       else
           data_norm{a} = data_{a};
           data_mean{a} =nanmean(data_{a});
           data_SEM{a}  =nansem(data_{a});
       end
       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end
    
    plot([15 15],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(9,maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
    text(19,maxH,'Choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

%     rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%     for t= 1:length(x_temp)
%         p(t) = ttest2(data_norm{2}(:,t),data_norm{3}(:,t));
%     end
%     addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
    
    a=smooth_hist(data_norm{2}'); b=smooth_hist(data_norm{3}');
    a(:,isnan(a(1,:)))=[];     b(:,isnan(b(1,:)))=[];

    [p,~] = permtest2vec(a,b,NBS,pValue);
%     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
%     a = nan(size(p));a(p)=1;
%     plot(x_temp,maxH*a+1,'-k','LineWidth',2.5)
%     set(gca,'XTick',[])
        
        sig1 =  (data_mean{3}>data_mean{2} & p'); 
        a = nan(size(sig1));a(sig1)=1;
    plot(x_temp,maxH*a+1,'-r','LineWidth',2.5)
    
        sig1 =  (data_mean{3}<data_mean{2} & p'); 
        a = nan(size(sig1));a(sig1)=1;
    plot(x_temp,maxH*a+1,'-b','LineWidth',2.5)
    axis([0 Inf 0 maxH+2 ])
    if s==1
        ylabel({'Normalized decoding score'})
    for a=1:length(label_str_)
        text(30,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','right' )
    end
    elseif s==2
        xlabel('Time (s)')
        text(30,maxH+1,'p<0.05')
    end
%% Plot group LR decoding (units/assemblies)
figure('name','Group LR decoding','color','w')
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
%% Plot group LR decoding on Correct vs. error trials...normalised F-scores
clear data_norm data_mean data_SEM

figure('name','Group LR decoding','color','w')
textAlign_ = 'middle';
normFlag = true;
normWin = [0 7];
LeverPressWidth = 2; % Width of lever press bar (s)
NBS = 100;           % No draws for permutation tests
pValue = 0.05;       % Bonferroni corrected inside the function


if normFlag
    maxH = 20;              % yMax of decoding plots
else
    maxH = 10;              % yMax of decoding plots
end


SAMPLEln  = [0.6 0.95 0.6];
CHOICEln  = [0.95 0.6 0.6];
alpha = 0.8;
x_temp = (1:P.Ltr*2)*P.bw;
bp=find(x_temp >=normWin(1) & x_temp <=normWin(2));
col_ = {'k',[0.9 0.2 0.2],'k'}; %[0.9 0.6 0.3]};
label_str_ = {'',... % All single units
              'Error trials',...
              'Correct trials'};
label_pos_ = [0,9,10];

for s=1:2
    subplot(2,3,s); hold on
    title([P.names{s} , ' Units'])

    plot([10 10],[0 maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    


    data_={nan*Group.Units.Ft2x{s},...
               Group.Units.Err.Ft2x{s},...  
               Group.Units.Ft2x{s}};
  
    for a =1:length(data_)
       if normFlag
           B = nanmean(data_{a}(:,bp),2);
           data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
           data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
           data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));
       else
           data_norm{a} = data_{a};
           data_mean{a} =nanmean(data_{a});
           data_SEM{a}  =nansem(data_{a});
       end
       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end
    
    plot([15 15],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(9,maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
    text(19,maxH,'Choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

%     rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%     for t= 1:length(x_temp)
%         p(t) = ttest2(data_norm{2}(:,t),data_norm{3}(:,t));
%     end
%     addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
    a=smooth_hist(data_norm{2}'); b=smooth_hist(data_norm{3}');
    a(:,isnan(a(1,:)))=[];     b(:,isnan(b(1,:)))=[];

    [p,~] = permtest2vec(a,b,NBS,pValue);
%     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
    a = nan(size(p));a(p)=1;
    plot(x_temp,maxH*a+1,'-k','LineWidth',2.5)
%     set(gca,'XTick',[])
    
    axis([0 Inf 0 maxH+2 ])
    if s==1
        ylabel({'Normalized decoding score'})
    for a=1:length(label_str_)
        text(30,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','right' )
    end
    elseif s==2
        xlabel('Time (s)')
        text(30,maxH+1,'p<0.05')
    end
end


if normFlag
    maxH = 20;              % yMax of decoding plots
else
    maxH = 10;              % yMax of decoding plots
end


for s=1:3
    subplot(2,3,s+3); hold on
    title([P.names{s} , ' Assemblies'])

    plot([10 10],[0 maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    


    data_={nan*Group.Units.Ft2x{s},...
       Group.Assem.Err.Ft2x{s},...  
       Group.Assem.Ft2x{s}};
  
    for a =1:length(data_)
       if normFlag
           B = nanmean(data_{a}(:,bp),2);
           data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
           data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
           data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));
       else
           data_norm{a} =data_{a};
           data_mean{a} =nanmean(data_{a});
           data_SEM{a}  =nansem(data_{a});
       end
       
       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end
    
    plot([15 15],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(9,maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
    text(19,maxH,'Choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

%     rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%     for t= 1:length(x_temp)
%         p(t) = ttest2(data_norm{2}(:,t),data_norm{3}(:,t));
%     end
%     addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
    a=smooth_hist(data_norm{2}'); b=smooth_hist(data_norm{3}');
    a(:,isnan(a(1,:)))=[];     b(:,isnan(b(1,:)))=[];

    [p,~] = permtest2vec(a,b,NBS,pValue);
%     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
    a = nan(size(p));a(p)=1;
    plot(x_temp,maxH*a+1,'-k','LineWidth',2.5)
%     set(gca,'XTick',[])
    
    axis([0 Inf 0 maxH+2 ])
    if s==1
        ylabel({'Normalized decoding score'})
    elseif s==2
        xlabel('Time (s)')
    elseif s==3
    text(30,maxH+1,'p<0.05')
    end
end
%% Plot group SC decoding on Correct vs. error trials...normalised F-scores
clear data_norm data_mean data_SEM
normFlag = true;
figure('name','Group SC decoding ','color','w')
% Plot decoding
textAlign_ = 'middle';
LeverPressWidth = 2;    % Width of lever press bar (s)
% SAMPLEln  = [0.89 0.95 0.89];
% CHOICEln  = [0.95 0.89 0.89];

NBS = 5000; % No draws for permutation tests

PRESSln  = [0.6 0.6 0.6];

alpha = 0.8;
x_temp = (1:P.Ltr_SC)*P.bw; x_temp = x_temp - 5;
bp=find(x_temp >=-5 & x_temp <=-4);
col_ = {'k',[0.9 0.2 0.2],'k'}; %[0.9 0.6 0.3]};
label_str_ = {'',... % All single units
              'Error trials',...
              'Correct trials'};
label_pos_ = [0,0.9,1]*maxH;

if normFlag
    maxH = 5;              % yMax of decoding plots
else
    maxH = 10;              % yMax of decoding plots
end

for s=1:2
    subplot(2,3,s); hold on
    title([P.names{s} , ' Units'])
    plot([0 0],[0 maxH],'color',PRESSln,'LineWidth',LeverPressWidth)
    data_={nan*Group.Units.SC_Ft2x{s},...
               Group.Units.Err.SC_Ft2x{s},...  
               Group.Units.SC_Ft2x{s}};
    for a =1:length(data_)
        if normFlag
            B = nanmean(data_{a}(:,bp),2);
            data_norm{a} = data_{a}./(B*ones(1,length(x_temp)));
            data_mean{a} = nanmean(data_{a}./(B*ones(1,length(x_temp))));
            data_SEM{a}  = nansem(data_{a}./(B*ones(1,length(x_temp))));
        else 
            B = nanmean(data_{a}(:,bp),2);
            data_norm{a} = data_{a};
            data_mean{a} = nanmean(data_{a});
            data_SEM{a}  = nansem(data_{a});
        end
       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end
    
    plot([0 0],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([0 0],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(-0.5,maxH,'Sample/Choice','Rotation',90,'Color',PRESSln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

%     rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%     for t= 1:length(x_temp)
%         p(t) = ttest2(data_norm{2}(:,t),data_norm{3}(:,t));
%     end
%     addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
    a=smooth_hist(data_norm{2}'); b=smooth_hist(data_norm{3}');
    a(:,isnan(a(1,:)))=[];     b(:,isnan(b(1,:)))=[];

    [p,~] = permtest2vec(a,b,NBS);
%     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
    a = nan(size(p));a(p)=1;
    plot(x_temp,maxH*a+1,'-k','LineWidth',2.5)
    
%     set(gca,'XTick',[])
    
    axis([min(x_temp) max(x_temp) 0 maxH+1 ])
    if s==1
        ylabel({'Normalized decoding score'})
    for a=1:length(label_str_)
        text(5,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','right' )
    end
    elseif s==2
        xlabel('Time (s)')
        text(max(x_temp),maxH+1,'p<0.05')
    end
end

if normFlag
    maxH = 25;              % yMax of decoding plots
else
    maxH = 35;              % yMax of decoding plots
end
for s=1:3
    subplot(2,3,s+3); hold on
    title([P.names{s} , ' Assemblies'])
    plot([0 0],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([0 0],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(-0.5,maxH,'Sample/Choice','Rotation',90,'Color',PRESSln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
   


    data_={nan*Group.Units.SC_Ft2x{s},...
       Group.Assem.Err.SC_Ft2x{s},...  
       Group.Assem.SC_Ft2x{s}};
  
    for a =1:length(data_)
        if normFlag
            B = nanmean(data_{a}(:,bp),2);
            data_norm{a} = data_{a}./(B*ones(1,length(x_temp)));
            data_mean{a} = nanmean(data_{a}./(B*ones(1,length(x_temp))));
            data_SEM{a}  = nansem(data_{a}./(B*ones(1,length(x_temp))));
        else 
            B = nanmean(data_{a}(:,bp),2);
            data_norm{a} = data_{a};
            data_mean{a} = nanmean(data_{a});
            data_SEM{a}  = nansem(data_{a});
        end

       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end
    
    plot([0 0],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([0 0],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(-0.5,maxH,'Sample/Choice','Rotation',90,'Color',PRESSln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

%     rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%     for t= 1:length(x_temp)
%         p(t) = ttest2(data_norm{2}(:,t),data_norm{3}(:,t));
%     end
%     addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
    a=smooth_hist(data_norm{2}'); b=smooth_hist(data_norm{3}');
    a(:,isnan(a(1,:)))=[];     b(:,isnan(b(1,:)))=[];

    [p,~] = permtest2vec(a,b,NBS);
%     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
    a = nan(size(p));a(p)=1;
    plot(x_temp,maxH*a+1,'-k','LineWidth',2.5)
%     set(gca,'XTick',[])
    
    axis([min(x_temp) max(x_temp) 0 maxH+1 ])
    if s==1
        ylabel({'Normalized decoding score'})
    elseif s==2
        xlabel('Time (s)')
    elseif s==3
    text(30,maxH+1,'p<0.05')
    end
end
%% Plot group LR decoding (members/assemblies) ...normalised F-scores
clear data_norm data_mean data_SEM

figure('name','Group LR decoding','color','w')
% Plot decoding
textAlign_ = 'middle';
maxH = 25;              % yMax of decoding plots
LeverPressWidth = 2;    % Width of lever press bar (s)
MinMaxAdvantage = 10;   % yMin and yMax of assembly advantage plot 
% SAMPLEln  = [0.89 0.95 0.89];
% CHOICEln  = [0.95 0.89 0.89];

SAMPLEln  = [0.6 0.95 0.6];
CHOICEln  = [0.95 0.6 0.6];
alpha = 0.8;
x_temp = (1:P.Ltr*2)*P.bw;
bp=find(x_temp >=0 & x_temp <=7);
col_ = {'k','k',[0.9 0.6 0.3]};
label_str_ = {'',... % All single units
              'Assemblies',...
              {'Assembly';'member';'units'}};
label_pos_ = [0,22,18];

for s=1:3
    subplot(1,3,s); hold on
    title(P.names{s})

    plot([10 10],[0 maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    


    data_={nan*Group.Units.Ft2x{s},...
       Group.Assem.Ft2x{s},...  
       Group.MemberUnits.LR_Ft2x{s}};
  
    for a =1:length(data_)
       
       B = nanmean(data_{a}(:,bp),2);
       data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
       data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
       data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));

       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end
    
    plot([15 15],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(9,maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
    text(19,maxH,'Choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

%     rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%     for t= 1:length(x_temp)
%         p(t) = ttest2(data_norm{2}(:,t),data_norm{3}(:,t));
%     end
%     addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
    a=smooth_hist(data_norm{2}'); b=smooth_hist(data_norm{3}');
    a(:,isnan(a(1,:)))=[];     b(:,isnan(b(1,:)))=[];

    [p,~] = permtest2vec(a,b,500);
%     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
    a = nan(size(p));a(p)=1;
    plot(x_temp,maxH*a+1,'-k','LineWidth',2.5)
%     set(gca,'XTick',[])
    
    axis([0 Inf 0 maxH+2 ])
    if s==1
        ylabel({'Normalized decoding score'})
%         axis([0 Inf 0 maxH ])
    for a=1:length(label_str_)
        text(30,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','right' )
    end
    elseif s==2
        xlabel('Time (s)')

%         axis([0 Inf 0 maxH ])
    elseif s==3
%         axis([0 Inf 0 maxH ])
    text(30,maxH+1,'p<0.05')
    end
end
%% Plot group LR decoding (members/assemblies) ...normalised F-scores vertical layout
clear data_norm data_mean data_SEM

figure('name','Group LR decoding','color','w')
% Plot decoding
textAlign_ = 'middle';
maxH = 25;              % yMax of decoding plots
LeverPressWidth = 2;    % Width of lever press bar (s)
MinMaxAdvantage = 10;   % yMin and yMax of assembly advantage plot 
% SAMPLEln  = [0.89 0.95 0.89];
% CHOICEln  = [0.95 0.89 0.89];

SAMPLEln  = [0.6 0.95 0.6];
CHOICEln  = [0.95 0.6 0.6];
alpha = 0.8;
x_temp = (1:P.Ltr*2)*P.bw;
bp=find(x_temp >=0 & x_temp <=7);
col_ = {'k','k',[0.9 0.6 0.3]};
label_str_ = {'',... % All single units
              'Assemblies',...
              {'Assembly';'member';'units'}};
label_pos_ = [0,22,18];

for s=1:3
    subplot(3,1,s); hold on
    title(P.names{s})

    plot([10 10],[0 maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    


    data_={nan*Group.Units.Ft2x{s},...
       Group.Assem.Ft2x{s},...  
       Group.MemberUnits.LR_Ft2x{s}};
  
    for a =1:length(data_)
       
       B = nanmean(data_{a}(:,bp),2);
       data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
       data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
       data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));

       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end
    
    plot([15 15],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(9,maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
    text(19,maxH,'Choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

%     rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%     for t= 1:length(x_temp)
%         p(t) = ttest2(data_norm{2}(:,t),data_norm{3}(:,t));
%     end
%     addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
    a=smooth_hist(data_norm{2}'); b=smooth_hist(data_norm{3}');
    a(:,isnan(a(1,:)))=[];     b(:,isnan(b(1,:)))=[];

    [p,~] = permtest2vec(a,b,500);
%     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
    a = nan(size(p));a(p)=1;
    plot(x_temp,maxH*a+1,'-k','LineWidth',2.5)
%     set(gca,'XTick',[])
    
    axis([0 Inf 0 maxH+2 ])
    if s==1
%         ylabel({'Normalized decoding score'})
%         axis([0 Inf 0 maxH ])
%     for a=1:length(label_str_)
%         text(30,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','right' )
%     end
%     elseif s==2
%         xlabel('Time (s)')

%         axis([0 Inf 0 maxH ])
    elseif s==3
%         axis([0 Inf 0 maxH ])
    text(30,maxH+1,'p<0.05')
    end
end
%% Plot group LR decoding on Correct vs. error trials...normalised F-scores
clear data_norm data_mean data_SEM

figure('name','Group LR decoding','color','w')
textAlign_ = 'middle';
normFlag = true;
normWin = [0 7];
LeverPressWidth = 2; % Width of lever press bar (s)
NBS = 100;           % No draws for permutation tests
pValue = 0.05;       % Bonferroni corrected inside the function


if normFlag
    maxH = 20;              % yMax of decoding plots
else
    maxH = 10;              % yMax of decoding plots
end


SAMPLEln  = [0.6 0.95 0.6];
CHOICEln  = [0.95 0.6 0.6];
alpha = 0.8;
x_temp = (1:P.Ltr*2)*P.bw;
bp=find(x_temp >=normWin(1) & x_temp <=normWin(2));
% col_ = {[0.9 0.6 0.3],[0.9 0.2 0.2],[0.9 0.6 0.3]}; %[0.9 0.60.3]};
col_ = {'k',[0.9 0.2 0.2],'k'}; %[0.9 0.6 0.3]};

label_str_ = {'',... % All single units
              'Error trials',...
              'Correct trials'};
label_pos_ = [0,9,10];

for s=1:2
    subplot(1,3,s); hold on
    title([P.names{s} , ' Units'])

    plot([10 10],[0 maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    


    data_={nan*Group.Units.Ft2x{s},...
               Group.Units.Err.Ft2x{s},...  
               Group.Units.Ft2x{s}};
  
    for a =1:length(data_)
       if normFlag
           B = nanmean(data_{a}(:,bp),2);
           data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
           data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
           data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));
       else
           data_norm{a} = data_{a};
           data_mean{a} =nanmean(data_{a});
           data_SEM{a}  =nansem(data_{a});
       end
       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end
    
    plot([15 15],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(9,maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
    text(19,maxH,'Choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

%     rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%     for t= 1:length(x_temp)
%         p(t) = ttest2(data_norm{2}(:,t),data_norm{3}(:,t));
%     end
%     addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
    a=smooth_hist(data_norm{2}'); b=smooth_hist(data_norm{3}');
    a(:,isnan(a(1,:)))=[];     b(:,isnan(b(1,:)))=[];

    [p,~] = permtest2vec(a,b,NBS,pValue);
%     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
    a = nan(size(p));a(p)=1;
    plot(x_temp,maxH*a+1,'-k','LineWidth',2.5)
%     set(gca,'XTick',[])
    
    axis([0 Inf 0 maxH+2 ])
    if s==1
        ylabel({'Normalized decoding score'})
    for a=1:length(label_str_)
        text(30,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','right' )
    end
    elseif s==2
        xlabel('Time (s)')
        text(30,maxH+1,'p<0.05')
    end
end


col_ = {'k',[0.9 0.2 0.2],'k'}; %[0.9 0.6 0.3]};

if normFlag
    maxH = 20;              % yMax of decoding plots
else
    maxH = 10;              % yMax of decoding plots
end


for s=3
    subplot(1,3,s); hold on
    title([P.names{s} , ' Assemblies'])

    plot([10 10],[0 maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    


    data_={nan*Group.Units.Ft2x{s},...
       Group.Assem.Err.Ft2x{s},...  
       Group.Assem.Ft2x{s}};
  
    for a =1:length(data_)
       if normFlag
           B = nanmean(data_{a}(:,bp),2);
           data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
           data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
           data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));
       else
           data_norm{a} =data_{a};
           data_mean{a} =nanmean(data_{a});
           data_SEM{a}  =nansem(data_{a});
       end
       
       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end
    
    plot([15 15],[0 maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    text(9,maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
    text(19,maxH,'Choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

%     rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%     for t= 1:length(x_temp)
%         p(t) = ttest2(data_norm{2}(:,t),data_norm{3}(:,t));
%     end
%     addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
    a=smooth_hist(data_norm{2}'); b=smooth_hist(data_norm{3}');
    a(:,isnan(a(1,:)))=[];     b(:,isnan(b(1,:)))=[];

    [p,~] = permtest2vec(a,b,NBS,pValue);
%     plot(x_temp(find(p)),0.9*maxH*p(find(p)),'*k')
    a = nan(size(p));a(p)=1;
    plot(x_temp,maxH*a+1,'-k','LineWidth',2.5)
%     set(gca,'XTick',[])
    
    axis([0 Inf 0 maxH+2 ])
    if s==1
        ylabel({'Normalized decoding score'})
    elseif s==2
        xlabel('Time (s)')
    elseif s==3
    text(30,maxH+1,'p<0.05')
    end
end
%% Plot group LR decoding (members/non-members/assemblies) ...normalised F-scores
figure('name','Group LR decoding','color','w')
% Plot decoding
maxH = 20;              % yMax of decoding plots
LeverPressWidth = 5;    % Width of lever press bar (s)
MinMaxAdvantage = 10;   % yMin and yMax of assembly advantage plot 
SAMPLEln  = [0.89 0.95 0.89];
CHOICEln  = [0.95 0.89 0.89];
alpha = 0.5;
x_temp = (1:P.Ltr*2)*P.bw;
bp=find(x_temp >=0 & x_temp <=7);
col_ = {'k','b','g','r'};
label_str_ = {'',... %All single units
          'Assemblies',...
          'Assembly member units',...
          'Assembly non-member units'};
label_pos_ = [-1,-3,-5,-7];

for s=1:3
    subplot(2,3,s); hold on
    title(P.names{s})

    plot([9 9],[0 10*maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 10*maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    


    data_={nan*Group.Units.Ft2x{s},...
       Group.Assem.Ft2x{s},...  
       Group.MemberUnits.LR_Ft2x{s},...
       Group.NonMemberUnits.LR_Ft2x{s}};
   
  
    for a =1:length(data_)
       B = nanmean(data_{a}(:,bp),2);
       data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
       data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
       data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));

       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end

    plot([15 15],[0 10*maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 10*maxH ],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    set(gca,'XTick',[])
    
    axis([0 Inf 0 maxH ])
    if s==1
        ylabel({'Decoding (norm.)'})
        axis([0 Inf 0 maxH ])

    elseif s==2
%         xlabel('Time (s)')
    for a=1:length(label_str_)
        text(15,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','center' )
    end
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
    area(x_temp,temp,'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',alpha,'EdgeAlpha',1)
    
    temp = Group.AssemAdvantage_Mean{s}; temp(temp<0)= 0;
    area(x_temp,temp,'FaceColor',[0.6 0.6 0.6],'EdgeColor','k','LineWidth',1.5,'FaceAlpha',alpha,'EdgeAlpha',1)
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
%% Plot group SC decoding (units/assemblies)
figure('name','Group SC decoding','color','w')
% Plot decoding
maxH = 40;              % yMax of decoding plots
LeverPressWidth = 2;    % Width of lever press bar (s)
MinMaxAdvantage = 20;   % yMin and yMax of assembly advantage plot 
alpha = 0.5;
x_temp = (1:P.Ltr_SC)*P.bw; x_temp = x_temp - 5;

for s=1:3
    subplot(2,3,s); hold on
    title(P.names{s})

    plot([0 0],[0 1*maxH],'color',[0.6 0.6 0.6],'LineWidth',LeverPressWidth)
    
    plot(x_temp,Group.Assem.SC_Ft2ciH0_Mean{s},'k')
    plot(x_temp,Group.Units.SC_Ft2ciH0_Mean{s},'b')

    ciplot(Group.Assem.SC_Ft2x_Mean{s}+Group.Assem.SC_Ft2x_SEM{s},...
           Group.Assem.SC_Ft2x_Mean{s}-Group.Assem.SC_Ft2x_SEM{s},...
           x_temp,[0.6 0.6 0.6],alpha)  
    ciplot(Group.Units.SC_Ft2x_Mean{s}+Group.Units.SC_Ft2x_SEM{s},...
           Group.Units.SC_Ft2x_Mean{s}-Group.Units.SC_Ft2x_SEM{s},...
           x_temp,[0.6 0.6 0.9],alpha)
	plot(x_temp,Group.Assem.SC_Ft2x_Mean{s},'k','LineWidth',1.5)
    plot(x_temp,Group.Units.SC_Ft2x_Mean{s},'b','LineWidth',1.5)
    set(gca,'XTick',[])
    
    if s==1
        ylabel({'Decoding score';'(F-value)'})

    elseif s==2
%         xlabel('Time (s)')
        text(0,-1,'Single unit decoding','Color',[0 0 1],'HorizontalAlignment','center' )
        text(0,-3,'Assembly decoding','Color',[0 0 0],'HorizontalAlignment','center' )
    axis([min(x_temp) max(x_temp) 0 maxH ])
    elseif s==3
    axis([min(x_temp) max(x_temp) 0 maxH ])

    end
end

% Plot assembly advantage
for s=1:3
    subplot(2,3,s+3); hold on

    plot([0 0],[0 1*maxH],'color',[0.6 0.6 0.6],'LineWidth',LeverPressWidth)
    
    temp = Group.AssemAdvantageSC_Mean{s}; temp(temp>0)= 0;
    area(x_temp,temp,'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',alpha,'EdgeAlpha',1)
    
    temp = Group.AssemAdvantageSC_Mean{s}; temp(temp<0)= 0;
    area(x_temp,temp,'FaceColor',[0.6 0.6 0.6],'EdgeColor','k','LineWidth',1.5,'FaceAlpha',alpha,'EdgeAlpha',1)
    axis([min(x_temp) max(x_temp) -5 MinMaxAdvantage])
    if s==1
        ylabel({'Assembly decoding advantage'})
    elseif s==2
        xlabel('Time (s)')
    end
end
%% Plot group SC decoding (members/non-members/assemblies) ...normalised F-scores
figure('name','Group SC decoding','color','w')
% Plot decoding
maxH = 20;              % yMax of decoding plots
LeverPressWidth = 2;    % Width of lever press bar (s)
MinMaxAdvantage = 20;   % yMin and yMax of assembly advantage plot 
alpha = 0.5;

col_ = {'k','b','g','r'};
label_str_ = {'',... %All single units
          'Assemblies',...
          'Assembly member units',...
          'Assembly non-member units'};
label_pos_ = [-1,-3,-5,-7];

x_temp = (1:P.Ltr_SC)*P.bw; x_temp = x_temp - 5;
bp=find(x_temp >=-5 & x_temp <=-4);

for s=1:3
    subplot(2,3,s); hold on
    title(P.names{s})

    plot([0 0],[0 1*maxH],'color',[0.6 0.6 0.6],'LineWidth',LeverPressWidth)
    


    
     data_={nan*Group.Units.SC_Ft2x{s},...
            Group.Assem.SC_Ft2x{s},...  
            Group.MemberUnits.SC_Ft2x{s},...
            Group.NonMemberUnits.SC_Ft2x{s}};
    for a =1:length(data_)
       B = nanmean(data_{a}(:,bp),2);
       data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
       data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
       data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));

       plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
       ciplot(data_mean{a}-data_SEM{a},...
              data_mean{a}+data_SEM{a},...
              x_temp,col_{a},alpha)  
    end

    set(gca,'XTick',[])
    
    if s==1
        ylabel({'Decoding (norm.)'})
    elseif s==2
        % xlabel('Time (s)')
        for a=1:length(label_str_)
            text(0,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','center' )
        end
        axis([min(x_temp) max(x_temp) 0 maxH ])
    elseif s==3
        axis([min(x_temp) max(x_temp) 0 maxH ])
    end
end

% Plot assembly advantage
for s=1:3
    subplot(2,3,s+3); hold on

    plot([0 0],[0 1*maxH],'color',[0.6 0.6 0.6],'LineWidth',LeverPressWidth)
    
    temp = Group.AssemAdvantageSC_Mean{s}; temp(temp>0)= 0;
    area(x_temp,temp,'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',alpha,'EdgeAlpha',1)
    
    temp = Group.AssemAdvantageSC_Mean{s}; temp(temp<0)= 0;
    area(x_temp,temp,'FaceColor',[0.6 0.6 0.6],'EdgeColor','k','LineWidth',1.5,'FaceAlpha',alpha,'EdgeAlpha',1)
    axis([min(x_temp) max(x_temp) -5 MinMaxAdvantage])
    if s==1
        ylabel({'Assembly decoding advantage'})
    elseif s==2
        xlabel('Time (s)')
    end
end
%% LR decoding SC decoding
figure('name','SC vs LR decoding','color','w')
for s=1:2
    subplot(2,3,s); hold on
    title([P.names{s} ' Units'])
%     tempLR = cell2mat(Group.Units.TSPeak{s}');
%     tempSC = cell2mat(Group.Units.SC_TSPeak{s}');
%     plot([0 30],[0 30],':k')
%     scatter(tempLR,tempSC)
%     xlabel('Peak spatial decoding')
%     ylabel('Peak contextual decoding')
%     axis square
%     axis([0 30 0 30])
    
    tempLR = cell2mat(Group.Units.TSIntegral{s}')./P.Ltr*2;
    tempSC = cell2mat(Group.Units.SC_TSIntegral{s}')./P.Ltr_SC;
    plot([0 25],[0 25],':k')
    scatter(tempLR,tempSC)
    xlabel('Spatial decoding integral')
    axis square
    axis([0 25 0 25])
    if s==1
        ylabel('Contextual decoding integral ')
    end
end

for s=1:3
    subplot(2,3,s+3); hold on
    title([P.names{s} ' Assemblies'])
%     tempLR = cell2mat(Group.Assem.TSPeak{s}');
%     tempSC = cell2mat(Group.Assem.SC_TSPeak{s}');
%     plot([0 30],[0 30],':k')
%     scatter(tempLR,tempSC)
%     xlabel('Peak spatial decoding')
%     ylabel('Peak contextual decoding')
%     axis square
%     axis([0 30 0 30])
    
    tempLR = cell2mat(Group.Assem.TSIntegral{s}')./P.Ltr*2;
    tempSC = cell2mat(Group.Assem.SC_TSIntegral{s}')./P.Ltr_SC;
    plot([0 25],[0 25],':k')
    scatter(tempLR,tempSC)
    xlabel('Spatial decoding integral')
    if s==1
        ylabel('Contextual decoding integral')
    end
    axis square
    axis([0 25 0 25])
    
end
%% Pies - Significant decoders - none / spatial only / contextual only / spatial and contextual
Group.Units.TSsigLRsum = [];
Group.Units.TSsigSCsum = [];
% Fraction of window to that unit/assembly must be a significant decoder
% for in order to be considered informative (0.05 corresponds to 1s) 
TimeCutoff = 0.05; 
ClassNames = {'Uninformative', 'Spatial only', 'Contextual only ', 'Spatial and Contextual'};
ClassColors = {[0.9 0.9 0.9],[1 0 0],[0 1 0],[0.3 0.3 0.9]};
for s= 1:3
    for iRat = 1:length(Rats)
        
        Group.Units.TSsigLRsum{s}{iRat}      =  sum(Group.Units.TSsigLR{s}{iRat})./size(Group.Units.TSsigLR{s}{iRat},1);
        Group.Units.TSsigSCsum{s}{iRat}      =  sum(Group.Units.TSsigSC{s}{iRat})./size(Group.Units.TSsigSC{s}{iRat},1);
%         LR = cellfun(@isempty,Group.MemberUnits.TSsigLR{s})
        Group.MemberUnits.TSsigLRsum{s}{iRat}      =  sum(Group.MemberUnits.TSsigLR{s}{iRat})./size(Group.MemberUnits.TSsigLR{s}{iRat},1);
        Group.MemberUnits.TSsigSCsum{s}{iRat}      =  sum(Group.MemberUnits.TSsigSC{s}{iRat})./size(Group.MemberUnits.TSsigSC{s}{iRat},1);
        
        if ~isempty(Group.Assem.TSsigLR{s}{iRat})
            Group.Assem.TSsigLRsum{s}{iRat}=sum(Group.Assem.TSsigLR{s}{iRat})./size(Group.Assem.TSsigLR{s}{iRat},1);
            Group.Assem.TSsigSCsum{s}{iRat}=sum(Group.Assem.TSsigSC{s}{iRat})./size(Group.Assem.TSsigSC{s}{iRat},1);
        else
            Group.Assem.TSsigLRsum{s}{iRat}=[];
            Group.Assem.TSsigSCsum{s}{iRat}=[];
        end
    end

    
    Group.Units.TSsigLRSC{s} = [cell2mat(Group.Units.TSsigLRsum{s})', cell2mat(Group.Units.TSsigSCsum{s})'];
    Group.MemberUnits.TSsigLRSC{s} = [cell2mat(Group.MemberUnits.TSsigLRsum{s})', cell2mat(Group.MemberUnits.TSsigSCsum{s})'];
    
    Group.Units.TSsigLRSCfrac{s}(1,1) =  sum(Group.Units.TSsigLRSC{s}(:,1) <= TimeCutoff & Group.Units.TSsigLRSC{s}(:,2) <= TimeCutoff);
    Group.Units.TSsigLRSCfrac{s}(1,2) =  sum(Group.Units.TSsigLRSC{s}(:,1)  > TimeCutoff & Group.Units.TSsigLRSC{s}(:,2) <= TimeCutoff);
    Group.Units.TSsigLRSCfrac{s}(1,3) =  sum(Group.Units.TSsigLRSC{s}(:,1) <= TimeCutoff & Group.Units.TSsigLRSC{s}(:,2)  > TimeCutoff);
    Group.Units.TSsigLRSCfrac{s}(1,4) =  sum(Group.Units.TSsigLRSC{s}(:,1)  > TimeCutoff & Group.Units.TSsigLRSC{s}(:,2)  > TimeCutoff);
    
    Group.Units.TSsigLRSCfrac{s}      = Group.Units.TSsigLRSCfrac{s}./sum(Group.Units.TSsigLRSCfrac{s});
    
    
    Group.Assem.TSsigLRSC{s} = [cell2mat(Group.Assem.TSsigLRsum{s})', cell2mat(Group.Assem.TSsigSCsum{s})'];
    
    Group.Assem.TSsigLRSCfrac{s}(1,1) =  sum(Group.Assem.TSsigLRSC{s}(:,1) <= TimeCutoff & Group.Assem.TSsigLRSC{s}(:,2) <= TimeCutoff);
    Group.Assem.TSsigLRSCfrac{s}(1,2) =  sum(Group.Assem.TSsigLRSC{s}(:,1)  > TimeCutoff & Group.Assem.TSsigLRSC{s}(:,2) <= TimeCutoff);
    Group.Assem.TSsigLRSCfrac{s}(1,3) =  sum(Group.Assem.TSsigLRSC{s}(:,1) <= TimeCutoff & Group.Assem.TSsigLRSC{s}(:,2)  > TimeCutoff);
    Group.Assem.TSsigLRSCfrac{s}(1,4) =  sum(Group.Assem.TSsigLRSC{s}(:,1)  > TimeCutoff & Group.Assem.TSsigLRSC{s}(:,2)  > TimeCutoff);
    
    Group.Assem.TSsigLRSCfrac{s}      = Group.Assem.TSsigLRSCfrac{s}./sum(Group.Assem.TSsigLRSCfrac{s});
 
end


figure('color','w');
hold on
x = 100*[cell2mat(Group.Units.TSsigLRsum{1}),cell2mat(Group.Units.TSsigLRsum{2})];
y = 100*[cell2mat(Group.Units.TSsigSCsum{1}),cell2mat(Group.Units.TSsigSCsum{2})];
group = [ones(length(cell2mat(Group.Units.TSsigLRsum{1})),1);2*ones(length(cell2mat(Group.Units.TSsigLRsum{2})),1)];
scatterhist(x,y,'Group',group,'NBins',50,...
    'Location','SouthWest',...
    'Direction','out',...
    'Kernel','on',...
    'Bandwidth',2*[1,1; 1,1],...
    'Color','rg','LineStyle',{'-','-'},...
    'LineWidth',[2,2,2],'Marker','.','MarkerSize',[5 5]);
    axis(100*[0 1 0 1]);  
    legend off
    axis square
    xlabel('Significant spatial decoding (% time)')
    ylabel('Significant contextual decoding (% time)')   
    %title('Units')

figure('color','w');
hold on
x = 100*[cell2mat(Group.Assem.TSsigLRsum{1}),cell2mat(Group.Assem.TSsigLRsum{2}),cell2mat(Group.Assem.TSsigLRsum{3})];
y = 100*[cell2mat(Group.Assem.TSsigSCsum{1}),cell2mat(Group.Assem.TSsigSCsum{2}),cell2mat(Group.Assem.TSsigSCsum{3})];
group = [ones(length(cell2mat(Group.Assem.TSsigLRsum{1})),1);...
         2*ones(length(cell2mat(Group.Assem.TSsigLRsum{2})),1);...
         3*ones(length(cell2mat(Group.Assem.TSsigLRsum{3})),1)];
scatterhist(x,y,'Group',group,'NBins',50,...
    'Location','SouthWest',...
    'Direction','out',...
    'Kernel','on',...
    'Bandwidth',2*[1,1,1; 1,1,1],...
    'Color','rgb','LineStyle',{'-','-','-'},...
    'LineWidth',[2,2,2],'Marker','.','MarkerSize',[5 5 5]);
    axis(100*[0 1 0 1]);  
    legend off
    axis square
    xlabel('Significant spatial decoding (% time)')
    ylabel('Significant contextual decoding (% time)')    
%     title('Assemblies')



figure('color','w');
for s= 1:2
    subplot(2,3,s); hold on
    scatter(cell2mat(Group.Units.TSsigLRsum{s}),...
            cell2mat(Group.Units.TSsigSCsum{s}))
%     scatterhist(cell2mat(Group.Units.TSsigLRsum{s}),cell2mat(Group.Units.TSsigSCsum{s}),'NBins',50,...
%         'Location','SouthWest',...
%         'Direction','out',...
%         'Kernel','on',...
%         'Bandwidth',[0.02; 0.02])
    axis([0 1 0 1]);  
    axis square
    xlabel('Fraction significant spatial decoders')
    ylabel('Fraction significant contextual decoders')
end
for s= 1:3
    subplot(2,3,s+3)
    scatter(cell2mat(Group.Assem.TSsigLRsum{s}),...
            cell2mat(Group.Assem.TSsigSCsum{s}))
    axis([0 1 0 1]);  
    axis square
end

figure('color','w');
for s = 1:2
    subplot(2,3,s); hold on
    temp = Group.Units.TSsigLRSCfrac{s};
    tempLabels = ClassNames(temp~=0);
    tempColors = cell2mat(ClassColors(temp~=0)');
    Hpie = pie(temp,ClassNames);
    scatter(0,0,10000,'w','filled')
    hp = findobj(Hpie,'Type','patch');
    for iClass= 1:numel(hp)
        set(hp(iClass), 'FaceColor', tempColors(iClass,:),...
            'EdgeColor', [0.6 0.6 0.6],...
            'FaceAlpha', 1,...
            'LineWidth', 1,...
            'LineStyle','none');
    end
    title({[P.names{s},' Units'];''})
    text(0,0,{sprintf('n = %d',size(Group.Units.TSsigLRSC{s},1));'Units'},'HorizontalAlignment','center')
    axis square
    axis off
end
for s = 1:3
    subplot(2,3,s+3); hold on
    
    temp = Group.Assem.TSsigLRSCfrac{s};
    tempLabels = ClassNames(temp~=0);
    tempColors = cell2mat(ClassColors(temp~=0)');
    Hpie = pie(temp,ClassNames);
    scatter(0,0,10000,'w','filled')
    hp = findobj(Hpie,'Type','patch');
    for iClass= 1:numel(hp)
        set(hp(iClass), 'FaceColor', tempColors(iClass,:),...
            'EdgeColor', [0.6 0.6 0.6],...
            'FaceAlpha', 1,...
            'LineWidth', 1,...
            'LineStyle','none');
    end
	text(0,0,{sprintf('n = %d',size(Group.Assem.TSsigLRSC{s},1));'Assemblies'},'HorizontalAlignment','center')

    title({P.names{s};'Assemblies';''})
    axis square
    axis off
end
%% scatter time significant vs. peak decoding
bins = 0:0.01:1;
sig = 0.05;
clear p
figure
temp_names = {'mPFC';'dCA1';{'Joint';'dCA1~mPFC'}};temp_col = {'b','r','g'};
for s= 1:3
    subplot(1,3,s); hold on
    u_BS=[];a_BS=[];
    
%     u_ = cell2mat(Group.Units.TSsigLRsum{s})'; %u_(isnan(u_))=[];
    u_time = cell2mat(Group.MemberUnits.TSsigLRsum{s})'; u_time(isnan(u_time))=[];
    u_time(u_time==0)=NaN;
    u_peak = cell2mat(cellfun(@max,Group.MemberUnits.LR_TS{s},'UniformOutput',false)');
    
    a_time = cell2mat(Group.Assem.TSsigLRsum{s})';
    a_peak = cell2mat(Group.Assem.TSPeak{s}');
    
     a_time(isnan(a_time))=[];
     a_peak(isnan(a_peak))=[];
    scatter(u_time,u_peak,15,[0.9 0.6 0.3],'filled')
    scatter(a_time,a_peak,15,'k','filled')
    axis([0 0.6 0 20])
end
%% Compare times as significant decoders - Ass vs Units
bins = 0:0.01:1;
sig = 0.05;
clear p
figure
temp_names = {'mPFC';'dCA1';{'Joint';'dCA1~mPFC'}};temp_col = {'b','r','g'};
for s= 1:3
    u_BS=[];a_BS=[];
    
%     u_ = cell2mat(Group.Units.TSsigLRsum{s})'; %u_(isnan(u_))=[];
    u_ = cell2mat(Group.MemberUnits.TSsigLRsum{s})';% u_(isnan(u_))=[];
    
    a_ = cell2mat(Group.Assem.TSsigLRsum{s})';
    no_Drawn = ceil(0.5*min([length(u_),length(a_)]));
    for i = 1:1000
        a_hBS(i,:) = cumsum(histc(randsample(a_, no_Drawn),bins));
        u_hBS(i,:) = cumsum(histc(randsample(u_, no_Drawn),bins));
    end
    a_hBSmean = nanmean(a_hBS);
    u_hBSmean = nanmean(u_hBS);
    
    for i=1:length(bins)
        a_sort = sort(a_hBS(:,i),'ascend');
        a_hBSsem(i,1) = a_sort(max(1,round(length(a_sort)*0.05)));
        a_hBSsem(i,2) = a_sort(max(1,round(length(a_sort)*(1-0.05))));

        u_sort = sort(u_hBS(:,i),'ascend');
        u_hBSsem(i,1) = u_sort(max(1,round(length(a_sort)*0.05)));
        u_hBSsem(i,2) = u_sort(max(1,round(length(a_sort)*(1-0.05))));
    end

    

    p(s) = kruskalwallis([u_;a_],...
                   [0*u_+1;0*a_+2],'off');
               
   subplot(3,1,s); hold on
%    u_h = cumsum(histc(u_,bins));
%    a_h = cumsum(histc(a_,bins));
   
%    plot(bins,u_h./max(u_h),'color','k','LineWidth',1.5)
%    plot(bins,a_h./max(a_h),'color',[0.9 0.6 0.3],'LineWidth',1.5)  

   plot(bins*30,a_hBSsem(:,1)./max(max(a_hBSsem)),'color','k','LineWidth',1,'LineStyle', '-')
   plot(bins*30,a_hBSsem(:,2)./max(max(a_hBSsem)),'color','k','LineWidth',1,'LineStyle', '-')

   plot(bins*30,u_hBSsem(:,1)./max(max(u_hBSsem)),'color',[0.9 0.6 0.3],'LineWidth',1,'LineStyle', '-')
   plot(bins*30,u_hBSsem(:,2)./max(max(u_hBSsem)),'color',[0.9 0.6 0.3],'LineWidth',1,'LineStyle', '-')
   
   
   plot(bins*30,a_hBSmean./max(a_hBSmean),'color','k','LineWidth',1.5)
   plot(bins*30,u_hBSmean./max(u_hBSmean),'color',[0.9 0.6 0.3],'LineWidth',1.5)
%    if p(s)<sig
%        text(0.6,0.5,{'***',sprintf('p=%03.4f',p(s))},'HorizontalAlignment','center')
%    else
%        text(0.6,0.5,{'n.s.',sprintf('p=%03.4f',p(s))},'HorizontalAlignment','center')
%    end
%    if s==2
%       ylabel('Distribution')
%    elseif s==3
%       xlabel({'Fraction of time';'Spent significant'})
%    end
%    text(0.95,0.2,temp_names{s},'color',temp_col{s},'HorizontalAlignment','right')
end
p_ = p<0.05;
clear s p p_ a_ a_h u_ u_h
%% Compare peak decoding
bins = 0:0.1:10;
sig = 0.05;
clear p
figure
temp_names = {'mPFC';'dCA1';{'Joint';'dCA1~mPFC'}};temp_col = {'b','r','g'};
for s= 1:3
    u_hBS=[];a_hBS=[];
    
    %     u_ = cell2mat(Group.Units.TSsigLRsum{s})'; %u_(isnan(u_))=[];

    u_peak = cell2mat(cellfun(@max,Group.MemberUnits.LR_TS{s},'UniformOutput',false)')';
    
    a_peak = cell2mat(Group.Assem.TSPeak{s}')';
    
        
    a_peak(isnan(a_peak))=[];
    u_ = u_peak;
    
    a_ = a_peak;
    no_Drawn = ceil(0.5*min([length(u_),length(a_)]));
    for i = 1:1000
        a_hBS(i,:) = cumsum(histc(randsample(a_, no_Drawn),bins));
        u_hBS(i,:) = cumsum(histc(randsample(u_, no_Drawn),bins));
    end
    a_hBSmean = nanmean(a_hBS);
    u_hBSmean = nanmean(u_hBS);
    
    for i=1:length(bins)
        a_sort = sort(a_hBS(:,i),'ascend');
        a_hBSsem(i,1) = a_sort(max(1,round(length(a_sort)*0.05)));
        a_hBSsem(i,2) = a_sort(max(1,round(length(a_sort)*(1-0.05))));

        u_sort = sort(u_hBS(:,i),'ascend');
        u_hBSsem(i,1) = u_sort(max(1,round(length(a_sort)*0.05)));
        u_hBSsem(i,2) = u_sort(max(1,round(length(a_sort)*(1-0.05))));
    end

    

    p(s) = kruskalwallis([u_;a_],...
                   [0*u_+1;0*a_+2],'off');
               
   subplot(3,1,s); hold on
%    u_h = cumsum(histc(u_,bins));
%    a_h = cumsum(histc(a_,bins));
   
%    plot(bins,u_h./max(u_h),'color','k','LineWidth',1.5)
%    plot(bins,a_h./max(a_h),'color',[0.9 0.6 0.3],'LineWidth',1.5)  

   plot(bins,a_hBSsem(:,1)./max(max(a_hBSsem)),'color','k','LineWidth',1,'LineStyle', '-')
   plot(bins,a_hBSsem(:,2)./max(max(a_hBSsem)),'color','k','LineWidth',1,'LineStyle', '-')

   plot(bins,u_hBSsem(:,1)./max(max(u_hBSsem)),'color',[0.9 0.6 0.3],'LineWidth',1,'LineStyle', '-')
   plot(bins,u_hBSsem(:,2)./max(max(u_hBSsem)),'color',[0.9 0.6 0.3],'LineWidth',1,'LineStyle', '-')
   
   
   plot(bins,a_hBSmean./max(a_hBSmean),'color','k','LineWidth',1.5)
   plot(bins,u_hBSmean./max(u_hBSmean),'color',[0.9 0.6 0.3],'LineWidth',1.5)
%    if p(s)<sig
%        text(0.3*max(bins),0.5,{'*',sprintf('p=%03.4f',p(s))},'HorizontalAlignment','center')
%    else
%        text(0.9*max(bins),0.5,{'n.s.',sprintf('p=%03.4f',p(s))},'HorizontalAlignment','center')
%    end
%    if s==2
%       ylabel('Distribution')
%    elseif s==3
%       xlabel({'Peak t-score'})
%    end
%    text(0.95*max(bins),0.2,temp_names{s},'color',temp_col{s},'HorizontalAlignment','right')
end
p_ = p<0.05;
clear s p p_ a_ a_h u_ u_h u_hBS a_hBS u_time u_time a_time a_peak 
clear a_hBSmean u_hBSmean a_hBSsem a_hBSsem u_sort u_hBSsem u_hBSsem
%% Compare decoding integral
bins = 0:60;
sig = 0.05;
clear p
figure
temp_names = {'mPFC';'dCA1';{'Joint';'dCA1~mPFC'}};temp_col = {'b','r','g'};
for s= 1:3
    u_hBS=[];a_hBS=[];
    
    %     u_ = cell2mat(Group.Units.TSsigLRsum{s})'; %u_(isnan(u_))=[];

    u_peak = cell2mat(cellfun(@sum,Group.MemberUnits.LR_TS{s},'UniformOutput',false)')';
    
%     a_peak = cell2mat(Group.Assem.TSIntegral{s}')';
        a_peak = cell2mat(cellfun(@sum,Group.Assem.LR_TS{s},'UniformOutput',false)')';

        
    a_peak(isnan(a_peak))=[];
    u_ = u_peak/30;
    
    a_ = a_peak/30;
    no_Drawn = ceil(0.5*min([length(u_),length(a_)]));
    for i = 1:1000
        a_hBS(i,:) = cumsum(histc(randsample(a_, no_Drawn),bins));
        u_hBS(i,:) = cumsum(histc(randsample(u_, no_Drawn),bins));
    end
    a_hBSmean = nanmean(a_hBS);
    u_hBSmean = nanmean(u_hBS);
    
    for i=1:length(bins)
        a_sort = sort(a_hBS(:,i),'ascend');
        a_hBSsem(i,1) = a_sort(max(1,round(length(a_sort)*0.05)));
        a_hBSsem(i,2) = a_sort(max(1,round(length(a_sort)*(1-0.05))));

        u_sort = sort(u_hBS(:,i),'ascend');
        u_hBSsem(i,1) = u_sort(max(1,round(length(a_sort)*0.05)));
        u_hBSsem(i,2) = u_sort(max(1,round(length(a_sort)*(1-0.05))));
    end

    

    p(s) = kruskalwallis([u_;a_],...
                   [0*u_+1;0*a_+2],'off');
               
   subplot(3,1,s); hold on
%    u_h = cumsum(histc(u_,bins));
%    a_h = cumsum(histc(a_,bins));
   
%    plot(bins,u_h./max(u_h),'color','k','LineWidth',1.5)
%    plot(bins,a_h./max(a_h),'color',[0.9 0.6 0.3],'LineWidth',1.5)  

   plot(bins,a_hBSsem(:,1)./max(max(a_hBSsem)),'color','k','LineWidth',1,'LineStyle', '-')
   plot(bins,a_hBSsem(:,2)./max(max(a_hBSsem)),'color','k','LineWidth',1,'LineStyle', '-')

   plot(bins,u_hBSsem(:,1)./max(max(u_hBSsem)),'color',[0.9 0.6 0.3],'LineWidth',1,'LineStyle', '-')
   plot(bins,u_hBSsem(:,2)./max(max(u_hBSsem)),'color',[0.9 0.6 0.3],'LineWidth',1,'LineStyle', '-')
   
   
   plot(bins,a_hBSmean./max(a_hBSmean),'color','k','LineWidth',1.5)
   plot(bins,u_hBSmean./max(u_hBSmean),'color',[0.9 0.6 0.3],'LineWidth',1.5)
   if p(s)<sig
       text(0.8*max(bins),0.5,{'*',sprintf('p=%03.3f',p(s))},'HorizontalAlignment','center')
   else
%        text(0.8*max(bins),0.5,{'n.s.',sprintf('p=%03.3f',p(s))},'HorizontalAlignment','center')
       text(0.8*max(bins),0.5,'n.s.','HorizontalAlignment','center')
   end
   if s==2
      ylabel('Distribution')
   elseif s==3
      xlabel({'Integral'})
   end
%    text(0.95*max(bins),0.2,temp_names{s},'color',temp_col{s},'HorizontalAlignment','right')
end
p_ = p<0.05;
clear s p p_ a_ a_h u_ u_h u_hBS a_hBS u_time u_time a_time a_peak 
clear a_hBSmean u_hBSmean a_hBSsem a_hBSsem u_sort u_hBSsem u_hBSsem

%% Compare peak decoding - areas overlaid
bins = 0:0.1:20;
sig = 0.05;
clear p
figure
temp_names = {'mPFC';'dCA1';{'Joint';'dCA1~mPFC'}};temp_col = {'b','r','g'};
subplot(2,1,1); hold on
for s= 1:3
    u_hBS=[];a_hBS=[];
    
    %     u_ = cell2mat(Group.Units.TSsigLRsum{s})'; %u_(isnan(u_))=[];

    u_peak = cell2mat(cellfun(@max,Group.MemberUnits.LR_TS{s},'UniformOutput',false)')';
    a_peak = cell2mat(Group.Assem.TSPeak{s}')'; a_peak(isnan(a_peak))=[];
    a_ = a_peak;
    u_ = u_peak;
    no_Drawn = ceil(0.5*min([length(u_),length(a_)]));
    for i = 1:1000
        a_hBS(i,:) = cumsum(histc(randsample(a_, no_Drawn),bins));
        u_hBS(i,:) = cumsum(histc(randsample(u_, no_Drawn),bins));
    end
    a_hBSmean = nanmean(a_hBS);
    u_hBSmean = nanmean(u_hBS);
    
    for i=1:length(bins)
        a_sort = sort(a_hBS(:,i),'ascend');
        a_hBSsem(i,1) = a_sort(max(1,round(length(a_sort)*0.05)));
        a_hBSsem(i,2) = a_sort(max(1,round(length(a_sort)*(1-0.05))));

        u_sort = sort(u_hBS(:,i),'ascend');
        u_hBSsem(i,1) = u_sort(max(1,round(length(a_sort)*0.05)));
        u_hBSsem(i,2) = u_sort(max(1,round(length(a_sort)*(1-0.05))));
    end

    


%    plot(bins,u_hBSsem(:,1)./max(max(u_hBSsem)),'color',temp_col{s},'LineWidth',1,'LineStyle', '-')
%    plot(bins,u_hBSsem(:,2)./max(max(u_hBSsem)),'color',temp_col{s},'LineWidth',1,'LineStyle', '-')
   
   
   plot(bins,u_hBSmean./max(u_hBSmean),'color',temp_col{s},'LineWidth',1.5)

   ylabel('Distribution')
   text(0.95*max(bins),0.2,temp_names{s},'color',temp_col{s},'HorizontalAlignment','right')
end

subplot(2,1,2); hold on
for s= 1:3
    u_hBS=[];a_hBS=[];
    
    %     u_ = cell2mat(Group.Units.TSsigLRsum{s})'; %u_(isnan(u_))=[];

    u_peak = cell2mat(cellfun(@max,Group.MemberUnits.LR_TS{s},'UniformOutput',false)')';
    a_peak = cell2mat(Group.Assem.TSPeak{s}')'; a_peak(isnan(a_peak))=[];
    
    u_ = u_peak;
    a_ = a_peak;
    no_Drawn = ceil(0.5*min([length(u_),length(a_)]));
    for i = 1:1000
        a_hBS(i,:) = cumsum(histc(randsample(a_, no_Drawn),bins));
        u_hBS(i,:) = cumsum(histc(randsample(u_, no_Drawn),bins));
    end
    a_hBSmean = nanmean(a_hBS);
    u_hBSmean = nanmean(u_hBS);
    
    for i=1:length(bins)
        a_sort = sort(a_hBS(:,i),'ascend');
        a_hBSsem(i,1) = a_sort(max(1,round(length(a_sort)*0.05)));
        a_hBSsem(i,2) = a_sort(max(1,round(length(a_sort)*(1-0.05))));

        u_sort = sort(u_hBS(:,i),'ascend');
        u_hBSsem(i,1) = u_sort(max(1,round(length(a_sort)*0.05)));
        u_hBSsem(i,2) = u_sort(max(1,round(length(a_sort)*(1-0.05))));
    end

%    plot(bins,a_hBSsem(:,1)./max(max(a_hBSsem)),'color',temp_col{s},'LineWidth',1,'LineStyle', '-')
%    plot(bins,a_hBSsem(:,2)./max(max(a_hBSsem)),'color',temp_col{s},'LineWidth',1,'LineStyle', '-')   
   
   plot(bins,a_hBSmean./max(a_hBSmean),'color',temp_col{s},'LineWidth',1.5)

   if s==2
      ylabel('Distribution')
   elseif s==3
      xlabel({'Peak t-score'})
   end
   text(0.95*max(bins),0.2,temp_names{s},'color',temp_col{s},'HorizontalAlignment','right')
end
xlabel({'Peak t-score'})

clear s p p_ a_ a_h u_ u_h u_hBS a_hBS u_time u_time a_time a_peak 
clear a_hBSmean u_hBSmean a_hBSsem a_hBSsem u_sort u_hBSsem u_hBSsem
%% Compare times as significant decoders x  peak decoding
bins = 0:0.1:5;
sig = 0.05;
clear p
figure
temp_names = {'mPFC';'dCA1';{'Joint';'dCA1~mPFC'}};temp_col = {'b','r','g'};
for s= 1:3
    u_hBS=[];a_hBS=[];
    
    %     u_ = cell2mat(Group.Units.TSsigLRsum{s})'; %u_(isnan(u_))=[];
    u_time = cell2mat(Group.MemberUnits.TSsigLRsum{s})'; u_time(isnan(u_time))=[];
    u_time(u_time==0)=NaN;
    u_peak = cell2mat(cellfun(@max,Group.MemberUnits.LR_TS{s},'UniformOutput',false)')';
    
    a_time = cell2mat(Group.Assem.TSsigLRsum{s})';
    a_peak = cell2mat(Group.Assem.TSPeak{s}')';
    
        
    a_time(isnan(a_time))=[];
    a_peak(isnan(a_peak))=[];
    u_ = u_peak .* u_time;
    
    a_ = a_peak .* a_time;
    no_Drawn = ceil(0.5*min([length(u_),length(a_)]));
    for i = 1:1000
        a_hBS(i,:) = cumsum(histc(randsample(a_, no_Drawn),bins));
        u_hBS(i,:) = cumsum(histc(randsample(u_, no_Drawn),bins));
    end
    a_hBSmean = nanmean(a_hBS);
    u_hBSmean = nanmean(u_hBS);
    
    for i=1:length(bins)
        a_sort = sort(a_hBS(:,i),'ascend');
        a_hBSsem(i,1) = a_sort(max(1,round(length(a_sort)*0.05)));
        a_hBSsem(i,2) = a_sort(max(1,round(length(a_sort)*(1-0.05))));

        u_sort = sort(u_hBS(:,i),'ascend');
        u_hBSsem(i,1) = u_sort(max(1,round(length(a_sort)*0.05)));
        u_hBSsem(i,2) = u_sort(max(1,round(length(a_sort)*(1-0.05))));
    end

    

    p(s) = kruskalwallis([u_;a_],...
                   [0*u_+1;0*a_+2],'off');
               
   subplot(3,1,s); hold on
%    u_h = cumsum(histc(u_,bins));
%    a_h = cumsum(histc(a_,bins));
   
%    plot(bins,u_h./max(u_h),'color','k','LineWidth',1.5)
%    plot(bins,a_h./max(a_h),'color',[0.9 0.6 0.3],'LineWidth',1.5)  

   plot(bins,a_hBSsem(:,1)./max(max(a_hBSsem)),'color','k','LineWidth',1,'LineStyle', '-')
   plot(bins,a_hBSsem(:,2)./max(max(a_hBSsem)),'color','k','LineWidth',1,'LineStyle', '-')

   plot(bins,u_hBSsem(:,1)./max(max(u_hBSsem)),'color',[0.9 0.6 0.3],'LineWidth',1,'LineStyle', '-')
   plot(bins,u_hBSsem(:,2)./max(max(u_hBSsem)),'color',[0.9 0.6 0.3],'LineWidth',1,'LineStyle', '-')
   
   
   plot(bins,a_hBSmean./max(a_hBSmean),'color','k','LineWidth',1.5)
   plot(bins,u_hBSmean./max(u_hBSmean),'color',[0.9 0.6 0.3],'LineWidth',1.5)
   if p(s)<sig
       text(0.9*max(bins),0.5,{'***',sprintf('p=%03.4f',p(s))},'HorizontalAlignment','center')
   else
       text(0.9*max(bins),0.5,{'n.s.',sprintf('p=%03.4f',p(s))},'HorizontalAlignment','center')
   end
   if s==2
      ylabel('Distribution')
   elseif s==3
      xlabel({'Fraction of time';'Spent significant'})
   end
   text(0.95*max(bins),0.2,temp_names{s},'color',temp_col{s},'HorizontalAlignment','right')
end
p_ = p<0.05;
clear s p p_ a_ a_h u_ u_h u_hBS a_hBS u_time u_time a_time a_peak 
clear a_hBSmean u_hBSmean a_hBSsem a_hBSsem u_sort u_hBSsem u_hBSsem
%% LR decoding SC decoding overlaid (units and assems separate)
col_ = {[0.5 0.5 0.9],[0.9 0.5 0.5],[0.5 0.9 0.5]};
col_2 = {'b','r','g'};
lims = [0 5];
figure('name','SC vs LR decoding','color','w')
subplot(1,2,1); hold on
    plot(lims,lims,':k')
    title('Units')
    for s=1:2
%         tempLR = cell2mat(Group.Units.TSPeak{s}');
%         tempSC = cell2mat(Group.Units.SC_TSPeak{s}');
        tempLR = mean(cell2mat(Group.MemberUnits.LR_TS{s}'));
        tempSC = mean(cell2mat(Group.MemberUnits.SC_TS{s}'));
%         tempLR = cell2mat(Group.Units.TSIntegral{s}')./P.Ltr*2;
%         tempSC = cell2mat(Group.Units.SC_TSIntegral{s}')./P.Ltr_SC;
        scatter(tempLR,tempSC,100,col_{s},'.')
    end
    for s=1:2
%         tempLR = cell2mat(Group.Units.TSPeak{s}');
%         tempSC = cell2mat(Group.Units.SC_TSPeak{s}');
        tempLR = mean(cell2mat(Group.MemberUnits.LR_TS{s}'));
        tempSC = mean(cell2mat(Group.MemberUnits.SC_TS{s}'));
%         tempLR = cell2mat(Group.Units.TSIntegral{s}')./P.Ltr*2;
%         tempSC = cell2mat(Group.Units.SC_TSIntegral{s}')./P.Ltr_SC;
        h= errorbarxy(mean(tempLR),mean(tempSC),nanstd(tempLR),nanstd(tempSC),{col_2{s}, col_2{s}, col_2{s}});
        for i=1:6
            set(h.hErrorbar(i),'LineWidth',2)
        end
        scatter(nanmean(tempLR),nanmean(tempSC),100,col_2{s},'filled')
    end
    xlabel('Spatial information (mean t-score)')    
    ylabel('Contextual information (mean t-score)')
    axis([lims lims])
    
subplot(1,2,2); hold on
    plot(lims,lims,':k')
    title('Assemblies')
    for s=1:3
%         tempLR = cell2mat(Group.Assem.TSPeak{s}');
%         tempSC = cell2mat(Group.Assem.SC_TSPeak{s}');
        tempLR = mean(cell2mat(Group.Assem.LR_TS{s}'));
        tempSC = mean(cell2mat(Group.Assem.SC_TS{s}'));
%         tempLR = cell2mat(Group.Assem.TSIntegral{s}')./P.Ltr*2;
%         tempSC = cell2mat(Group.Assem.SC_TSIntegral{s}')./P.Ltr_SC;
        scatter(tempLR,tempSC,100,col_{s},'.')
    end
    for s=1:3
%         tempLR = cell2mat(Group.Assem.TSPeak{s}');
%         tempSC = cell2mat(Group.Assem.SC_TSPeak{s}');
        tempLR = mean(cell2mat(Group.Assem.LR_TS{s}'));
        tempSC = mean(cell2mat(Group.Assem.SC_TS{s}'));
%         tempLR = cell2mat(Group.Assem.TSIntegral{s}')./P.Ltr*2;
%         tempSC = cell2mat(Group.Assem.SC_TSIntegral{s}')./P.Ltr_SC;
        h= errorbarxy(nanmean(tempLR),nanmean(tempSC),nanstd(tempLR),nanstd(tempSC),{col_2{s}, col_2{s}, col_2{s}});
        for i=1:6
            set(h.hErrorbar(i),'LineWidth',2)
        end
        scatter(nanmean(tempLR),nanmean(tempSC),100,col_2{s},'filled')
    end
    xlabel('Spatial information (mean t-score)')    
    ylabel('Contextual information (mean t-score)')
    axis([lims lims])
%% LR decoding SC decoding overlaid (units and assems compared)

col_point = {[0.5 0.5 0.9],[0.9 0.5 0.5],[0.5 0.9 0.5]};
col_point = {[0 0 0],[0 0 0],[0 0 0]};
col_mean = {'b','r','g'};
lims = [0 5];
figure('name','SC vs LR decoding','color','w')
for s=1:3
    subplot(1,3,s); hold on
    tempLR = mean(cell2mat(Group.MemberUnits.LR_TS{s}'));
    tempSC = mean(cell2mat(Group.MemberUnits.SC_TS{s}'));
%     scatter(tempLR,tempSC,10,col_point{s},'filled','Markerfacecolor',[0.6 0.6 0.6],'Markeredgecolor',[0 0 0])
    h= errorbarxy(nanmean(tempLR),nanmean(tempSC),nanstd(tempLR),nanstd(tempSC));
    for i=1:6
        set(h.hErrorbar(i),'LineWidth',2,'color',col_point{s})
    end
    scatter(nanmean(tempLR),nanmean(tempSC),100,col_point{s},'filled','Markeredgecolor',[0 0 0])

    tempLR = mean(cell2mat(Group.Assem.LR_TS{s}'));
    tempSC = mean(cell2mat(Group.Assem.SC_TS{s}'));
%     scatter(tempLR,tempSC,20,col_mean{s},'filled','Markeredgecolor',[0 0 0])
    h= errorbarxy(nanmean(tempLR),nanmean(tempSC),nanstd(tempLR),nanstd(tempSC));
    for i=1:6
        set(h.hErrorbar(i),'LineWidth',2,'color',col_mean{s})
    end
    scatter(nanmean(tempLR),nanmean(tempSC),100,col_mean{s},'filled')

    plot(lims,lims,':k')
    title(P.names{s})
    
    if s==1
        ylabel('Contextual information (mean t-score)')
    elseif s==2
        xlabel('Spatial information (mean t-score)')    
    end
        
    axis([lims lims]); axis square
    
end

clear h p
rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
for s=1:3
    tempAss   = mean(cell2mat(Group.Assem.SC_TS{s}'));
    tempUnits = mean(cell2mat(Group.MemberUnits.SC_TS{s}'));
    [p(s), h(s)]  = ttest2(tempAss,tempUnits)    
end
addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\nansuite')
%% Plot group LR and SC decoding


maxH = 40;              % yMax of decoding plots
LeverPressWidth = 2;    % Width of lever press bar (s)
MinMaxAdvantage = 10;   % yMin and yMax of assembly advantage plot
showVariance_ = false;
textAlign_ = 'middle';
SAMPLEln  = [0.6 0.95 0.6];
CHOICEln  = [0.95 0.6 0.6];
Press_line_col_=0.8*[1 1 1];
alpha = 0.8;

col_ = {'k','b','g','r'};
label_str_ = {'',... %All single units
    'Assemblies',...
    'Assembly member units',...
    'Assembly non-member units'};
label_pos_ = [-1,-2,-3,-4];

figure('name','Group decoding','color','w')

%%%% Plot Spatial information
x_temp = (1:P.Ltr*2)*P.bw;
bp=find(x_temp >=0 & x_temp <=7);% bp=1;
for s=1:3
    subplot(2,3,s); hold on;    
    title(P.names{s})
    plot([10 10],[0 0.86*maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
    plot([20 20],[0 0.86*maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
    text(10, maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
    text(20,maxH,'Choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
    data_={nan*Group.Units.Ft2x{s},...
           Group.Assem.Ft2x{s},...
           Group.MemberUnits.LR_Ft2x{s},...
           Group.NonMemberUnits.LR_Ft2x{s}};
    
    for a =1:length(data_)
        B = nanmean(data_{a}(:,bp),2);
        data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
        data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
        data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));
        
        plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
        if ~showVariance_
            h= area(x_temp,data_mean{a});
            set(h(1),'FaceColor', col_{a},'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
        else
               ciplot(data_mean{a}-data_SEM{a},...
                      data_mean{a}+data_SEM{a},...
                      x_temp,col_{a},alpha)
        end
    end
    
    plot([15 15],[0 5],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 5],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    % plot(x_temp,ones(size(x_temp)),'color',[0.6 0.6 0.6],'LineWidth',1)
    %     set(gca,'XTick',[])
    axis([min(x_temp) max(x_temp) 0 maxH ])
    if s==1
        ylabel({'Spatial Decoding (norm.)'})
        axis([0 Inf 0 maxH ])       
    elseif s==2
        %         xlabel('Time (s)')
        for a=1:length(label_str_)
            text(15,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','center' )
        end
        axis([0 Inf 0 maxH ])
    elseif s==3
        axis([0 Inf 0 maxH ])
    end
end

%%%% Plot Contextual information
x_temp = (1:P.Ltr_SC)*P.bw; x_temp = x_temp - 5;
bp=find(x_temp >=-5 & x_temp <=-4); bp=1;
for s=1:3
    subplot(2,3,s+3); hold on;    % title(P.names{s})
    plot([0 0],[0 0.58*maxH],'color',Press_line_col_,'LineWidth',LeverPressWidth)
    text(0,maxH,'Sample / Choice Press','Rotation',90,'Color',Press_line_col_,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')

    data_={nan*Group.Units.SC_Ft2x{s},...
        Group.Assem.SC_Ft2x{s},...
        Group.MemberUnits.SC_Ft2x{s},...
        Group.NonMemberUnits.SC_Ft2x{s}};
    
    for a =1:length(data_)
        B = nanmean(data_{a}(:,bp),2);
        data_norm{a} =data_{a}./(B*ones(1,length(x_temp)));
        data_mean{a} =nanmean(data_{a}./(B*ones(1,length(x_temp))));
        data_SEM{a}  =nansem(data_{a}./(B*ones(1,length(x_temp))));
        
        plot(x_temp,data_mean{a},'color', col_{a},'LineWidth',1.5)
        if ~showVariance_

            h= area(x_temp,data_mean{a});
            set(h(1),'FaceColor', col_{a},'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
        else
            ciplot(data_mean{a}-data_SEM{a},...
                   data_mean{a}+data_SEM{a},...
                   x_temp,col_{a},alpha)
        end
    end
    % plot(x_temp,ones(size(x_temp)),'color',[0.6 0.6 0.6],'LineWidth',1)

    %     set(gca,'XTick',[])
    axis([min(x_temp) max(x_temp) 0 maxH ])
    
    if s==1
        ylabel({'Contextual Decoding (norm.)'})
    elseif s==2
        xlabel('Time (s)')
        %         for a=1:length(label_str_)
        %             text(0,label_pos_(a),label_str_{a},'Color',col_{a},'HorizontalAlignment','center' )
        %         end
    elseif s==3
    end
end
%% Plot group Assembly LR and SC decoding overlaid - correct
maxH = 20;              % yMax of decoding plots
LeverPressWidth = 2;    % Width of lever press bar (s)
MinMaxAdvantage = 10;   % yMin and yMax of assembly advantage plot
showVariance_ = 0;
textAlign_ = 'middle';
SAMPLEln  = [0.6 0.95 0.6];
CHOICEln  = [0.95 0.6 0.6];
Press_line_col_=0.8*[1 1 1];
alpha = 0.8;
normDecoding = true;
col_ = {'b','r','g'};
label_pos_ = [-2,-3,-4];

figure('name','Correct decoding','color','w')

%%%% Plot Spatial information
x_temp = (1:P.Ltr*2)*P.bw;
bp=find(x_temp >=0 & x_temp <=7);% bp=1;
subplot(1,8,[1 6]); hold on
plot([10 10],[0 maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
plot([20 20],[0 maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
text(9,maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
text(19,maxH,'Correct choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
%  title({'Spatial';'Information'})
for s=[2 3 1]
    data_=Group.Assem.Ft2x{s};
    if normDecoding
        B = nanmean(data_(:,bp),2);
        data_norm =data_./(B*ones(1,length(x_temp)));
        data_mean =nanmean(data_./(B*ones(1,length(x_temp))));
        data_SEM  =nansem(data_./(B*ones(1,length(x_temp))));
    else
        data_norm =data_;
        data_mean =nanmean(data_);
        data_SEM  =nansem(data_);
    end
    plot(x_temp,data_mean,'color', col_{s},'LineWidth',1.5)
    if ~showVariance_
        h= area(x_temp,data_mean);
        set(h(1),'FaceColor', col_{s},'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
    else
        ciplot(data_mean-data_SEM,...
            data_mean+data_SEM,...
            x_temp,col_{s},alpha)
    end
end
plot([15 15],[0 5],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
plot([15 15],[0 5],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
axis([min(x_temp) max(x_temp) 0 maxH])
if normDecoding
    ylabel({'Assembly Decoding (norm.)'})
else
    ylabel({'Assembly Decoding'})
end
axis([0 Inf 0 maxH ])
xlabel('Time (s)')
% text(25, 27,'Assemblies','color','k','HorizontalAlignment','center','Fontsize',12)
text(25, 9,'PFC','color','b','HorizontalAlignment','left','Fontsize',12)
text(25, 7,'HP','color','r','HorizontalAlignment','left','Fontsize',12)
text(25, 5,{'Joint HP~PFC'},'color','g','HorizontalAlignment','left','Fontsize',12)


%%%% Plot Contextual information
x_temp = (1:P.Ltr_SC)*P.bw; x_temp = x_temp - 5;
bp=find(x_temp >=-5 & x_temp <=-4);
subplot(1,6,[6]); hold on;
plot([0 0],[0 0.6*maxH],'color',Press_line_col_,'LineWidth',LeverPressWidth)
text(0,0.62*maxH,'Sample / Choice','Rotation',90,'Color',Press_line_col_,'VerticalAlignment',textAlign_,'HorizontalAlignment','left')
    
% title({'Contextual';'Information'})
for s=1:3
    data_=Group.Assem.SC_Ft2x{s};
    if normDecoding
        B = nanmean(data_(:,bp),2);
        data_norm =data_./(B*ones(1,length(x_temp)));
        data_mean =nanmean(data_./(B*ones(1,length(x_temp))));
        data_SEM  =nansem(data_./(B*ones(1,length(x_temp))));
    else
        data_norm =data_;
        data_mean =nanmean(data_);
        data_SEM  =nansem(data_);
    end
    plot(x_temp,data_mean,'color', col_{s},'LineWidth',1.5)
    if ~showVariance_
        
        h= area(x_temp,data_mean);
        set(h(1),'FaceColor', col_{s},'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
    else
        ciplot(data_mean-data_SEM,...
            data_mean+data_SEM,...
            x_temp,col_{s},alpha)
    end
end
set(gca,'XTick',[-5 0 5],'XTickLabel',[-5 0 5])
axis([-5  5 0 maxH])

% ylabel({'Contextual Decoding (norm.)'})
xlabel('Time (s)')
%% Plot group Assembly LR and SC decoding overlaid - error
maxH = 20;              % yMax of decoding plots
LeverPressWidth = 2;    % Width of lever press bar (s)
MinMaxAdvantage = 10;   % yMin and yMax of assembly advantage plot
showVariance_ = 0;
textAlign_ = 'middle';
SAMPLEln  = [0.6 0.95 0.6];
CHOICEln  = [0.95 0.6 0.6];
Press_line_col_=0.8*[1 1 1];
alpha = 0.8;
normDecoding = true;
col_ = {'b','r','g'};
label_pos_ = [-2,-3,-4];

figure('name','Incorrect decoding','color','w')

%%%% Plot Spatial information
x_temp = (1:P.Ltr*2)*P.bw;
bp=find(x_temp >=0 & x_temp <=7);% bp=1;
subplot(1,8,[1 6]); hold on
plot([10 10],[0 maxH],'color',SAMPLEln,'LineWidth',LeverPressWidth)
plot([20 20],[0 maxH],'color',CHOICEln,'LineWidth',LeverPressWidth)
text(9,maxH,'Sample','Rotation',90,'Color',SAMPLEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
text(19,maxH,'Incorrect choice','Rotation',90,'Color',CHOICEln,'VerticalAlignment',textAlign_,'HorizontalAlignment','right')
%  title({'Spatial';'Information'})
for s=[2 3 1]
    data_=Group.Assem.Err.Ft2x{s};
    if normDecoding
        B = nanmean(data_(:,bp),2);
        data_norm =data_./(B*ones(1,length(x_temp)));
        data_mean =nanmean(data_./(B*ones(1,length(x_temp))));
        data_SEM  =nansem(data_./(B*ones(1,length(x_temp))));
    else
        data_norm =data_;
        data_mean =nanmean(data_);
        data_SEM  =nansem(data_);
    end
    plot(x_temp,data_mean,'color', col_{s},'LineWidth',1.5)
    if ~showVariance_
        h= area(x_temp,data_mean);
        set(h(1),'FaceColor', col_{s},'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
    else
        ciplot(data_mean-data_SEM,...
            data_mean+data_SEM,...
            x_temp,col_{s},alpha)
    end
end
plot([15 15],[0 5],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
plot([15 15],[0 5],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
axis([min(x_temp) max(x_temp) 0 maxH])
if normDecoding
    ylabel({'Assembly Decoding (norm.)'})
else
    ylabel({'Assembly Decoding'})
end
axis([0 Inf 0 maxH ])
xlabel('Time (s)')
% text(25, 27,'Assemblies','color','k','HorizontalAlignment','center','Fontsize',12)
text(25, 9,'PFC','color','b','HorizontalAlignment','left','Fontsize',12)
text(25, 7,'HP','color','r','HorizontalAlignment','left','Fontsize',12)
text(25, 5,{'Joint HP~PFC'},'color','g','HorizontalAlignment','left','Fontsize',12)


%%%% Plot Contextual information
x_temp = (1:P.Ltr_SC)*P.bw; x_temp = x_temp - 5;
bp=find(x_temp >=-5 & x_temp <=-4);
subplot(1,6,[6]); hold on;
plot([0 0],[0 0.6*maxH],'color',Press_line_col_,'LineWidth',LeverPressWidth)
text(0,0.62*maxH,'Sample / Choice','Rotation',90,'Color',Press_line_col_,'VerticalAlignment',textAlign_,'HorizontalAlignment','left')
    
% title({'Contextual';'Information'})
for s=1:3
    data_=Group.Assem.Err.SC_Ft2x{s};
    if normDecoding
        B = nanmean(data_(:,bp),2);
        data_norm =data_./(B*ones(1,length(x_temp)));
        data_mean =nanmean(data_./(B*ones(1,length(x_temp))));
        data_SEM  =nansem(data_./(B*ones(1,length(x_temp))));
    else
        data_norm =data_;
        data_mean =nanmean(data_);
        data_SEM  =nansem(data_);
    end
    plot(x_temp,data_mean,'color', col_{s},'LineWidth',1.5)
    if ~showVariance_
        
        h= area(x_temp,data_mean);
        set(h(1),'FaceColor', col_{s},'EdgeColor','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
    else
        ciplot(data_mean-data_SEM,...
            data_mean+data_SEM,...
            x_temp,col_{s},alpha)
    end
end
set(gca,'XTick',[-5 0 5],'XTickLabel',[-5 0 5])
axis([-5  5 0 maxH])

% ylabel({'Contextual Decoding (norm.)'})
xlabel('Time (s)')
%% Plot significant decoding time TS




temp_TS = [];
temp_TS_mask = [];

for s = 1:3
    
    figure('color','w');
    subplot(1,2,1)
        temp_TS  = cell2mat(Group.MemberUnits.LR_TS{s}');
        temp_TS_mask  = cell2mat(Group.MemberUnits.TSsigLR{s}');
%     temp_TS  = cell2mat(Group.Units.LR_TS{s}');
%     temp_TS_mask  = cell2mat(Group.Units.TSsigLR{s}');
    %     temp_TS=zscore(temp_TS);
    
    temp_TS(~temp_TS_mask)=NaN;
    [~,idx] = nanmax(temp_TS',2);
    [~,idx] = sort(idx);
    temp_TS = temp_TS(:,idx) ;
    temp_TS_{s} = temp_TS(:,idx) ;
    h= pcolor(temp_TS');
    set(h, 'EdgeColor', 'none');
    
    set(gca,'YDir','normal');
    colormap((hot))
    caxis([0 10])
    rectangle('Position',[0 0, size(temp_TS)],'LineWidth',2)
    axis off; axis square
    axis([-1 size(temp_TS,1) -1 size(temp_TS,2)])
    
    temp_TS = [];
    temp_TS_mask = [];
    subplot(1,2,2)
        tf= ~isempty_cell(Group.Assem.TSsigLR{s});

        temp_TS  = cell2mat(Group.Assem.LR_TS{s}(tf)');
        temp_TS_mask  = cell2mat(Group.Assem.TSsigLR{s}(tf)');
        

    %     temp_TS=zscore(temp_TS);
    
    temp_TS(~temp_TS_mask)=NaN;
    [~,idx] = nanmax(temp_TS',2);
    [~,idx] = sort(idx);
    temp_TS = temp_TS(:,idx) ;
    temp_TS_{s} = temp_TS(:,idx) ;
    h= pcolor(temp_TS');
    set(h, 'EdgeColor', 'none');
    
    set(gca,'YDir','normal');
    colormap((hot))
    caxis([0 10])
    rectangle('Position',[0 0, size(temp_TS)],'LineWidth',2)
    axis off; axis square
    axis([-1 size(temp_TS,1) -1 size(temp_TS,2)])
    
    
end
%% How many assemblies?
nos = [];
for iRat = 1:length(Rats)
    for s =1:3
        nos(iRat,s) = size(Group.Assem.SC_TS{s}{iRat},2);
    end
end
    
y = [mean(nos)'];                  %The data.
s = [nansem(nos)'];                %The standard deviation.
fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
colors = [0 0 1; 1 0 0;0 1 0];
for i = 1:numel(y)
    bar(i, y(i), 'parent', aHand, 'facecolor', 'w', 'edgecolor', colors(i,:),'LineWidth',1.5);
    errorbar(i,y(i),s(i),'color',colors(i,:),'LineWidth',1.5);

end
set(gca, 'XTick', 1:numel(y), 'XTickLabel', {'mPFC','CA1','Joint'})
ylabel('Assembly count')

% set(gca,'XTickLabelRotation',-45)
 
 

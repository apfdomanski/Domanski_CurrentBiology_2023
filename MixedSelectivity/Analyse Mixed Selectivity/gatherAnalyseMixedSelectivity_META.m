clear 
pat = 'C:\Analysis\AssemblyAnalysis\raw\';

% fileList_SHORT = dir([pat 'MixedSelectivity_SHORT' filesep '*SHORT*.mat']);
% fileList_MEDIUM = dir([pat 'MixedSelectivity_MEDIUM' filesep '*MEDIUM*.mat']);
% fileList_LONG = dir([pat 'MixedSelectivity' filesep '*LONG*','*PFC.mat']);

fileList_SHORT = dir([pat 'allTimestamps' filesep '*SHORT*.mat']);
fileList_MEDIUM = dir([pat 'allTimestamps' filesep '*MEDIUM*.mat']);
fileList_LONG = dir([pat 'allTimestamps' filesep '*LONG*.mat']);
Areas = {'PFC','HP','Joint'};

%% batch import unit data for meta-analysis and plotting - SHORT
% Import
fileList = fileList_SHORT;
Delays_ = {'Delay_0'};
for iArea = 1:2%length(Areas)
    D_.Delay_0.ANOVA.pThresh{iArea} = [];
    
    for iFile =1:length(fileList)
        try
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity_SHORT\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
        load(fnIn ,'D');
        
        for iDelay = 1:length(Delays_)
           
            % ANOVA: Fractions of tuned cells
            eval(sprintf('D_.%s.ANOVA.prcSig{iArea}(:,iFile)   = D.%s.ANOVA.prcSig;',Delays_{iDelay},Delays_{iDelay}))
            % ANOVA: Fractions of pure and mixed tuned cells
            eval(sprintf('D_.%s.ANOVA.prcTuned{iArea}(1:5,iFile)  = D.%s.ANOVA.prcTuned;',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA.prcTuned{iArea}(6,iFile)    = sum(D.%s.ANOVA.prcTuned(3:4));',Delays_{iDelay},Delays_{iDelay}))
                    
            eval(sprintf('D_.%s.ANOVA.pThresh{iArea}{iFile}     = D.%s.ANOVA.p<0.05;',Delays_{iDelay},Delays_{iDelay}))
            
            % Fano factor: Response selectivity for one condition
            eval(sprintf('D_.%s.Fano.FFR{iArea}{iFile,1}  = D.%s.Fano.FFR;',Delays_{iDelay},Delays_{iDelay}))
            % Fano factor: Trial-to-trial variability in firing rate for each trial type
            eval(sprintf('D_.%s.Fano.FFt{iArea}{iFile}  = D.%s.Fano.FFt;',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.Fano.FFtmean{iArea}(:,iFile)  = D.%s.Fano.FFtmean;',Delays_{iDelay},Delays_{iDelay}))
        end
        FFtNames = D.Delay_0.Fano.FFtNames;
        tuning   = [D.Delay_0.ANOVA.tuning,'Any Mixed'];
        Factors  = D.Delay_0.ANOVA.Factors;
        clear D
    end
    end
end

% Collapse population means/eerors
for iArea = 1:2%length(Areas)
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_.%s.ANOVA.prcSigMean{iArea}  = nanmean(D_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcSigSEM{iArea}   = nansem(D_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcTunedMean{iArea}  = nanmean(D_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcTunedSEM{iArea}  = nansem(D_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFRmean{iArea}  = cellfun(@nanmean,D_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFRsem{iArea}  = cellfun(@nansem,D_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFtmean_popmean{iArea}  = nanmean(D_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFtmean_popSEM{iArea}   = nansem(D_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
    end
end

clear D Stability D_temp D_tempL D_tempR B bp 
%% batch import unit data for meta-analysis and plotting - MEDIUM
% Import
fileList = fileList_MEDIUM;
Delays_ = {'Delay_2','Delay_4','Delay_6','Delay_8',};
for iArea = 1:2%length(Areas)
    
    for iFile =1:length(fileList)
        try
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity_MEDIUM\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
        load(fnIn ,'D');
        
        for iDelay = 1:length(Delays_)
           
            % ANOVA: Fractions of tuned cells
            eval(sprintf('D_.%s.ANOVA.prcSig{iArea}(:,iFile)   = D.%s.ANOVA.prcSig;',Delays_{iDelay},Delays_{iDelay}))
            % ANOVA: Fractions of pure and mixed tuned cells
            eval(sprintf('D_.%s.ANOVA.prcTuned{iArea}(1:5,iFile)  = D.%s.ANOVA.prcTuned;',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA.prcTuned{iArea}(6,iFile)    = sum(D.%s.ANOVA.prcTuned(3:4));',Delays_{iDelay},Delays_{iDelay}))

            eval(sprintf('D_.%s.ANOVA.pThresh{iArea}{iFile}     = D.%s.ANOVA.p<0.05;',Delays_{iDelay},Delays_{iDelay}))
            
            % Fano factor: Response selectivity for one condition
            eval(sprintf('D_.%s.Fano.FFR{iArea}{iFile,1}  = D.%s.Fano.FFR;',Delays_{iDelay},Delays_{iDelay}))
            % Fano factor: Trial-to-trial variability in firing rate for each trial type
            eval(sprintf('D_.%s.Fano.FFt{iArea}{iFile}  = D.%s.Fano.FFt;',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.Fano.FFtmean{iArea}(:,iFile)  = D.%s.Fano.FFtmean;',Delays_{iDelay},Delays_{iDelay}))
        end
      
        clear D
    end
    end
end

% Collapse population means/eerors
for iArea = 1:2%length(Areas)
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_.%s.ANOVA.prcSigMean{iArea}  = nanmean(D_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcSigSEM{iArea}   = nansem(D_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcTunedMean{iArea}  = nanmean(D_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcTunedSEM{iArea}  = nansem(D_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFRmean{iArea}  = cellfun(@nanmean,D_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFRsem{iArea}  = cellfun(@nansem,D_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFtmean_popmean{iArea}  = nanmean(D_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFtmean_popSEM{iArea}   = nansem(D_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
    end
end

clear D Stability D_temp D_tempL D_tempR B bp 
%% batch import unit data for meta-analysis and plotting - LONG
% Import
fileList = fileList_LONG;
Delays_ = {'Short','Medium','Long'};
for iArea = 1:2%length(Areas)
    
    for iFile =1:length(fileList)
        try
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
        load(fnIn ,'D');
        
        for iDelay = 1:length(Delays_)
           
            % ANOVA: Fractions of tuned cells
            eval(sprintf('D_.%s.ANOVA.prcSig{iArea}(:,iFile)   = D.%s.ANOVA.prcSig;',Delays_{iDelay},Delays_{iDelay}))
            % ANOVA: Fractions of pure and mixed tuned cells
            eval(sprintf('D_.%s.ANOVA.prcTuned{iArea}(1:5,iFile)  = D.%s.ANOVA.prcTuned;',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA.prcTuned{iArea}(6,iFile)    = sum(D.%s.ANOVA.prcTuned(3:4));',Delays_{iDelay},Delays_{iDelay}))

            eval(sprintf('D_.%s.ANOVA.pThresh{iArea}{iFile}     = D.%s.ANOVA.p<0.05;',Delays_{iDelay},Delays_{iDelay}))
            
            % Fano factor: Response selectivity for one condition
            eval(sprintf('D_.%s.Fano.FFR{iArea}{iFile,1}  = D.%s.Fano.FFR;',Delays_{iDelay},Delays_{iDelay}))
            % Fano factor: Trial-to-trial variability in firing rate for each trial type
            eval(sprintf('D_.%s.Fano.FFt{iArea}{iFile}  = D.%s.Fano.FFt;',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.Fano.FFtmean{iArea}(:,iFile)  = D.%s.Fano.FFtmean;',Delays_{iDelay},Delays_{iDelay}))
        end
        FFtNames = D.Short.Fano.FFtNames;
        tuning   = [D.Delay_0.ANOVA.tuning,'Any Mixed'];
        Factors  = D.Short.ANOVA.Factors;
        clear D
    end
    end
end

% Collapse population means/eerors
for iArea = 1:2%length(Areas)
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_.%s.ANOVA.prcSigMean{iArea}  = nanmean(D_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcSigSEM{iArea}   = nansem(D_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcTunedMean{iArea}  = nanmean(D_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcTunedSEM{iArea}  = nansem(D_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFRmean{iArea}  = cellfun(@nanmean,D_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFRsem{iArea}  = cellfun(@nansem,D_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFtmean_popmean{iArea}  = nanmean(D_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFtmean_popSEM{iArea}   = nansem(D_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
    end
end

clear D Stability D_temp D_tempL D_tempR B bp 
%% tidy up
clear Delays_ fileList fname fnIn iArea iDelay iFile 

%% Plot proportions of tuned/untuned Units - group by delay
Delays_ =  fieldnames(D_);
Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};

for iArea = 1:2
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    for iDelay = 1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('D = D_.%s;',Delays_{iDelay}))
        bar(D.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar(D.ANOVA.prcTunedMean{iArea},D.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
        title(Delays__{iDelay})
        axis([0 6 0 100])
        if iDelay==1
            ylabel('% of units')
        end    
                h=gca; h.Layer = 'top';

    end
end
clear D
%% Plot proportions of tuned/untuned Units - group by tuning
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};

for iArea = 1:2
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    for iTuning = 1:length(tuning)
        subplot(1,length(tuning),iTuning);hold on
        area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
        text(mean([2.2 6.8]),95,'Middle','HorizontalAlignment','center')
        text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D = D_.%s;',Delays_{iDelay}))
            m(iDelay) = D.ANOVA.prcTunedMean{iArea}(iTuning);
            e(iDelay) = D.ANOVA.prcTunedSEM{iArea}(iTuning);
        end
        bar([1,3,4,5,6,8,9,10],m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar([1,3,4,5,6,8,9,10],m,e,'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
        clear m e
        title(tuning{iTuning})
        axis([0 11 0 100])
        h=gca; h.Layer = 'top';

        if iTuning==1
            ylabel('% of units')
        end
    end
end
clear D
%% Plot proportions of tuned/untuned Units - group by factor


for iArea = 1:2
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    for iType = 1:length(Factors)
        subplot(1,length(Factors),iType);hold on
        area([0.2 1.8],0.5*90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],0.5*90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],0.5*90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        text(mean([0.2 1.8]),0.5*95,'Early','HorizontalAlignment','center')
        text(mean([2.2 6.8]),0.5*95,'Middle','HorizontalAlignment','center')
        text(mean([7.2 10.8]),0.5*95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D = D_.%s;',Delays_{iDelay}))
            m(iDelay) = D.ANOVA.prcSigMean{iArea}(iType);
            e(iDelay) = D.ANOVA.prcSigSEM{iArea}(iType);
        end
        bar([1,3,4,5,6,8,9,10],m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar([1,3,4,5,6,8,9,10],m,e,'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
        clear m e
        title(Factors{iType})
        axis([0 11 0 50])
        h=gca; h.Layer = 'top';

        if iType==1
            ylabel('% of units')
        end
    end
end
clear D        
   
%% Plot only medium and long delays, untuned, pure and mixed tuning
% tuning{2} = 
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
for iArea = 1:2
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]);
    tuning__ = [1 2 3 6 4 5];
    for iTuning = 1:length(tuning__)
        subplot(1,length(tuning__),iTuning);hold on
%         area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
%         text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
        text(mean([2.2 6.8]),95,'Early','HorizontalAlignment','center')
        text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D = D_.%s;',Delays_{iDelay}))
            m(iDelay) = D.ANOVA.prcTunedMean{iArea}(tuning__(iTuning));
            e(iDelay) = D.ANOVA.prcTunedSEM{iArea}(tuning__(iTuning));
        end
        
        bar([1,3,4,5,6,8,9,10],m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar([1,3,4,5,6,8,9,10],m,e,'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
        clear m e
        title(tuning{tuning__(iTuning)})
        axis([2 11 0 100])
        h=gca; h.Layer = 'top';

        if iTuning==1
            ylabel('% of units')
        end
        
        
    end
end

%% Run stats on pure and mixed tuning
h = [];
p=[];
tuning__ = [1 2 3 6 4 5];
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
for iArea = 1:2
    
    m_medium =[];
    for iDelay = 2:4
        eval(sprintf('D = D_.%s;',Delays_{iDelay}))
        m_medium = [m_medium,D.ANOVA.prcTuned{iArea}];
    end
    
    m_long =[];
    for iDelay = 5:8
        eval(sprintf('D = D_.%s;',Delays_{iDelay}))
        m_long = [m_long,D.ANOVA.prcTuned{iArea}];
    end
    
    
   
    
    path(rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite'))
    
   
    for iTuning = 1:length(tuning)
        [h(iTuning,iArea),p(iTuning,iArea)] =   ttest2(m_medium(tuning__(iTuning),:),m_long(tuning__(iTuning),:))
    end

    path(addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite'))
end

%% batch import Assembly data for meta-analysis and plotting - SHORT
% Import
fileList = fileList_SHORT;
Delays_ = {'Delay_0'};
%     D_Ass.Delay_0.ANOVA.pThresh{iArea} = [];
    
    for iFile =1:length(fileList)
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity_SHORT\\%s_MixedSelectivity_Ass.mat',pat,fname);
        load(fnIn ,'D_Ass');
        fprintf('Opening %s...\n',fname)

        for iDelay = 1:length(Delays_)
            for iArea = 1:length(Areas)
   
            % ANOVA: Fractions of tuned cells
            try
                eval(sprintf('D_Ass_.%s.ANOVA.prcSig{iArea}(:,iFile)   = D_Ass.%s.ANOVA.prcSig{iArea};',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.ANOVA.prcSig{iArea}(:,iFile)   = nan(1,7);',Delays_{iDelay}))
            end
            % ANOVA: Fractions of pure and mixed tuned cells
            try
                eval(sprintf('D_Ass_.%s.ANOVA.prcTuned{iArea}(1:5,iFile)  = D_Ass.%s.ANOVA.prcTuned{iArea};',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_Ass_.%s.ANOVA.prcTuned{iArea}(6,iFile)    = sum(D_Ass.%s.ANOVA.prcTuned{iArea}(3:4));',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.ANOVA.prcTuned{iArea}(1:6,iFile)  = nan(1,6);',Delays_{iDelay}))
            end
            try
                eval(sprintf('D_Ass_.%s.ANOVA.pThresh{iArea}{iFile}     = D_Ass.%s.ANOVA.p{iArea}<0.05;',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.ANOVA.pThresh{iArea}{iFile}     = nan(1,7);',Delays_{iDelay},Delays_{iDelay}))
            end
            
            % Fano factor: Response selectivity for one condition
            try
                eval(sprintf('D_Ass_.%s.Fano.FFR{iArea}{iFile,1}  = D_Ass.%s.Fano.FFR{iArea};',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.Fano.FFR{iArea}{iFile,1}  = nan;',Delays_{iDelay}))
            end
            
            % Fano factor: Trial-to-trial variability in firing rate for each trial type
            try
                eval(sprintf('D_Ass_.%s.Fano.FFt{iArea}{iFile}  = D_Ass.%s.Fano.FFt{iArea};',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_Ass_.%s.Fano.FFtmean{iArea}(:,iFile)  = D_Ass.%s.Fano.FFtmean{iArea};',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.Fano.FFt{iArea}{iFile}  = nan(1,4);',Delays_{iDelay}))
                eval(sprintf('D_Ass_.%s.Fano.FFtmean{iArea}(:,iFile)  = nan(1,4);',Delays_{iDelay}))
            end
        end
%         FFtNames = D_Ass.Delay_0.Fano.FFtNames;
%         tuning   = [D_Ass.Delay_0.ANOVA.tuning,'Any Mixed'];
%         Factors  = D_Ass.Delay_0.ANOVA.Factors;
        clear D_Ass
    end
end

% Collapse population means/eerors
for iArea = 1:length(Areas)
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_Ass_.%s.ANOVA.prcSigMean{iArea}  = nanmean(D_Ass_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.ANOVA.prcSigSEM{iArea}   = nansem(D_Ass_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.ANOVA.prcTunedMean{iArea}  = nanmean(D_Ass_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.ANOVA.prcTunedSEM{iArea}  = nansem(D_Ass_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFRmean{iArea}  = cellfun(@nanmean,D_Ass_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFRsem{iArea}  = cellfun(@nansem,D_Ass_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFtmean_popmean{iArea}  = nanmean(D_Ass_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFtmean_popSEM{iArea}   = nansem(D_Ass_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
    end
end

clear D_Ass Stability D_temp D_tempL D_tempR B bp 
%% batch import Assembly data for meta-analysis and plotting - MEDIUM
% Import
fileList = fileList_MEDIUM;
reject_list={'NorbertMEDIUM1_Events.mat','NorbertMEDIUM2_Events.mat','OnufryMEDIUM1_Events.mat','OnufryMEDIUM2_Events.mat'}; %'ALL_events.mat'
% reject_list={'MiroslawLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

Delays_ = {'Delay_2','Delay_4','Delay_6','Delay_8'};

for iFile =1:length(fileList)
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity_MEDIUM\\%s_MixedSelectivity_Ass.mat',pat,fname);
        load(fnIn ,'D_Ass');
        fprintf('Opening %s...\n',fname)

        for iDelay = 1:length(Delays_)
            for iArea = 1:length(Areas)
   
            % ANOVA: Fractions of tuned cells
            try
                eval(sprintf('D_Ass_.%s.ANOVA.prcSig{iArea}(1:7,iFile)   = D_Ass.%s.ANOVA.prcSig{iArea};',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.ANOVA.prcSig{iArea}(:,iFile)   = nan(1,7);',Delays_{iDelay}))
            end
            % ANOVA: Fractions of pure and mixed tuned cells
            try
                eval(sprintf('D_Ass_.%s.ANOVA.prcTuned{iArea}(1:5,iFile)  = D_Ass.%s.ANOVA.prcTuned{iArea};',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_Ass_.%s.ANOVA.prcTuned{iArea}(6,iFile)    = sum(D_Ass.%s.ANOVA.prcTuned{iArea}(3:4));',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.ANOVA.prcTuned{iArea}(1:6,iFile)  = nan(1,6);',Delays_{iDelay}))
            end
            try
                eval(sprintf('D_Ass_.%s.ANOVA.pThresh{iArea}{iFile}     = D_Ass.%s.ANOVA.p{iArea}<0.05;',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.ANOVA.pThresh{iArea}{iFile}     = nan(1,7);',Delays_{iDelay},Delays_{iDelay}))
            end
            
            % Fano factor: Response selectivity for one condition
            try
                eval(sprintf('D_Ass_.%s.Fano.FFR{iArea}{iFile,1}  = D_Ass.%s.Fano.FFR{iArea};',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.Fano.FFR{iArea}{iFile,1}  = nan;',Delays_{iDelay}))
            end
            
            % Fano factor: Trial-to-trial variability in firing rate for each trial type
            try
                eval(sprintf('D_Ass_.%s.Fano.FFt{iArea}{iFile}  = D_Ass.%s.Fano.FFt{iArea};',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_Ass_.%s.Fano.FFtmean{iArea}(:,iFile)  = D_Ass.%s.Fano.FFtmean{iArea};',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.Fano.FFt{iArea}{iFile}  = nan(1,4);',Delays_{iDelay}))
                eval(sprintf('D_Ass_.%s.Fano.FFtmean{iArea}(:,iFile)  = nan(1,4);',Delays_{iDelay}))
            end
        end
       
        clear D_Ass
    end
end

% Collapse population means/errors
for iArea = 1:length(Areas)
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_Ass_.%s.ANOVA.prcSigMean{iArea}  = nanmean(D_Ass_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.ANOVA.prcSigSEM{iArea}   = nansem(D_Ass_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.ANOVA.prcTunedMean{iArea}  = nanmean(D_Ass_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.ANOVA.prcTunedSEM{iArea}  = nansem(D_Ass_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFRmean{iArea}  = cellfun(@nanmean,D_Ass_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFRsem{iArea}  = cellfun(@nansem,D_Ass_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFtmean_popmean{iArea}  = nanmean(D_Ass_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFtmean_popSEM{iArea}   = nansem(D_Ass_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
    end
end

clear D_Ass Stability D_temp D_tempL D_tempR B bp 
%% batch import Assembly data for meta-analysis and plotting - LONG
% Import
fileList = fileList_LONG;
reject_list={'MiroslawLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

Delays_ = {'Short','Medium','Long'};

for iFile =1:length(fileList)
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity_LONG\\%s_MixedSelectivity_Ass.mat',pat,fname);
        load(fnIn ,'D_Ass');
        fprintf('Opening %s...\n',fname)

        for iDelay = 1:length(Delays_)
            for iArea = 1:length(Areas)
   
            % ANOVA: Fractions of tuned cells
            try
                eval(sprintf('D_Ass_.%s.ANOVA.prcSig{iArea}(:,iFile)   = D_Ass.%s.ANOVA.prcSig{iArea};',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.ANOVA.prcSig{iArea}(:,iFile)   = nan(1,7);',Delays_{iDelay}))
            end
            % ANOVA: Fractions of pure and mixed tuned cells
            try
                eval(sprintf('D_Ass_.%s.ANOVA.prcTuned{iArea}(1:5,iFile)  = D_Ass.%s.ANOVA.prcTuned{iArea};',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_Ass_.%s.ANOVA.prcTuned{iArea}(6,iFile)    = sum(D_Ass.%s.ANOVA.prcTuned{iArea}(3:4));',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.ANOVA.prcTuned{iArea}(1:6,iFile)  = nan(1,6);',Delays_{iDelay}))
            end
            try
                eval(sprintf('D_Ass_.%s.ANOVA.pThresh{iArea}{iFile}     = D_Ass.%s.ANOVA.p{iArea}<0.05;',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.ANOVA.pThresh{iArea}{iFile}     = nan(1,7);',Delays_{iDelay},Delays_{iDelay}))
            end
            
            % Fano factor: Response selectivity for one condition
            try
                eval(sprintf('D_Ass_.%s.Fano.FFR{iArea}{iFile,1}  = D_Ass.%s.Fano.FFR{iArea};',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.Fano.FFR{iArea}{iFile,1}  = nan;',Delays_{iDelay}))
            end
            
            % Fano factor: Trial-to-trial variability in firing rate for each trial type
            try
                eval(sprintf('D_Ass_.%s.Fano.FFt{iArea}{iFile}  = D_Ass.%s.Fano.FFt{iArea};',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_Ass_.%s.Fano.FFtmean{iArea}(:,iFile)  = D_Ass.%s.Fano.FFtmean{iArea};',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_Ass_.%s.Fano.FFt{iArea}{iFile}  = nan(1,4);',Delays_{iDelay}))
                eval(sprintf('D_Ass_.%s.Fano.FFtmean{iArea}(:,iFile)  = nan(1,4);',Delays_{iDelay}))
            end
        end
       
        clear D_Ass
    end
end

% Collapse population means/errors
for iArea = 1:length(Areas)
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_Ass_.%s.ANOVA.prcSigMean{iArea}  = nanmean(D_Ass_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.ANOVA.prcSigSEM{iArea}   = nansem(D_Ass_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.ANOVA.prcTunedMean{iArea}  = nanmean(D_Ass_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.ANOVA.prcTunedSEM{iArea}  = nansem(D_Ass_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFRmean{iArea}  = cellfun(@nanmean,D_Ass_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFRsem{iArea}  = cellfun(@nansem,D_Ass_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFtmean_popmean{iArea}  = nanmean(D_Ass_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_Ass_.%s.Fano.FFtmean_popSEM{iArea}   = nansem(D_Ass_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
    end
end

clear D_Ass Stability D_temp D_tempL D_tempR B bp 
%% tidy up
clear Delays_ fileList fname fnIn iArea iDelay iFile 
%% Plot proportions of tuned/untuned Units - group by delay
Delays_ =  fieldnames(D_);
Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};

for iArea = 1:3
    figure('name',['Percentage of tuned cells: ' Areas{iArea}],'Position',[334,572,1097,260]); 
    for iDelay = 1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('D = D_Ass_.%s;',Delays_{iDelay}))
        bar(D.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar(D.ANOVA.prcTunedMean{iArea},D.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
        title(Delays__{iDelay})
        axis([0 7 0 100])
        if iDelay==1
            ylabel('% of Assemblies')
        end    
        h=gca; h.Layer = 'top';
    end
end
clear D
%% Plot proportions of tuned/untuned Units - group by tuning
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};

for iArea = 1:3
    figure('name',['Percentage of tuned cells: ' Areas{iArea}],'Position',[334,572,1097,260]); 
    for iTuning = 1:length(tuning)
        subplot(1,length(tuning),iTuning);hold on
        area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
        text(mean([2.2 6.8]),95,'Middle','HorizontalAlignment','center')
        text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D = D_Ass_.%s;',Delays_{iDelay}))
            m(iDelay) = D.ANOVA.prcTunedMean{iArea}(iTuning);
            e(iDelay) = D.ANOVA.prcTunedSEM{iArea}(iTuning);
        end
        bar([1,3,4,5,6,8,9,10],m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar([1,3,4,5,6,8,9,10],m,e,'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
        clear m e
        title(tuning{iTuning})
        axis([0 11 0 100])
        h=gca; h.Layer = 'top';

        if iTuning==1
            ylabel('% of units')
        end
    end
end
clear D
%% Plot proportions of tuned/untuned Units - group by factor


for iArea = 1:3
    figure('name',['Percentage of tuned Assemblies: ' Areas{iArea}],'Position',[334,572,1097,260]); 
    for iType = 1:length(Factors)
        subplot(1,length(Factors),iType);hold on
        area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
        text(mean([2.2 6.8]),95,'Middle','HorizontalAlignment','center')
        text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D = D_Ass_.%s;',Delays_{iDelay}))
            m(iDelay) = D.ANOVA.prcSigMean{iArea}(iType);
            e(iDelay) = D.ANOVA.prcSigSEM{iArea}(iType);
        end
        bar([1,3,4,5,6,8,9,10],m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar([1,3,4,5,6,8,9,10],m,e,'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
        clear m e
        title(Factors{iType})
        axis([0 11 0 100])
        h=gca; h.Layer = 'top';

        if iType==1
            ylabel('% of units')
        end
    end
end
clear D  
%% Plot only medium and long delays, untuned, pure and mixed tuning
% tuning{2} = 
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
for iArea = 1:3
    figure('name',['Percentage of tuned assemblies: ' Areas{iArea}],'Position',[334,572,1097,260]);
    tuning__ = [1 2 3 6 4 5];
    for iTuning = 1:length(tuning__)
        subplot(1,length(tuning__),iTuning);hold on
%         area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
%         text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
        text(mean([2.2 6.8]),95,'Early','HorizontalAlignment','center')
        text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D = D_Ass_.%s;',Delays_{iDelay}))
            m(iDelay) = D.ANOVA.prcTunedMean{iArea}(tuning__(iTuning));
            e(iDelay) = D.ANOVA.prcTunedSEM{iArea}(tuning__(iTuning));
        end
        
        bar([1,3,4,5,6,8,9,10],m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar([1,3,4,5,6,8,9,10],m,e,'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
        clear m e
        title(tuning{tuning__(iTuning)})
        axis([2 11 0 100])
        h=gca; h.Layer = 'top';

        if iTuning==1
            ylabel('% of Assemblies')
        end
        
        
    end
end

clear 

pat = 'C:\Analysis\AssemblyAnalysis\raw\';
Areas = {'PFC','HP','Joint'};
Targets = {'SHORT','MEDIUM','LONG'};
Targets_ = {'Short','Medium','Long'};
DelaysList_ = {{'Delay_0'};{'Delay_2','Delay_4','Delay_6','Delay_8'};{'Short','Medium','Long'}};

DelaysList__ = {{'0s'};{'2s','4s','6s','8s'};{'4s','8s','16s'}};

Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};
useResampledTrials = true;


for iTarget = 1:length(Targets)
    if useResampledTrials
        eval(sprintf('fileList_%s = dir([pat ''MixedSelectivity\\TrialNosMatched\\*%s*PFC*Units*.mat'']);',Targets{iTarget},Targets{iTarget}) )
        eval(sprintf('AssfileList_%s = dir([pat ''MixedSelectivity\\TrialNosMatched\\*%s*PFC*Ass*.mat'']);',Targets{iTarget},Targets{iTarget}) )
    else
        eval(sprintf('fileList_%s = dir([pat ''MixedSelectivity\\*%s*PFC*Units*.mat'']);',Targets{iTarget},Targets{iTarget}) )
        eval(sprintf('AssfileList_%s = dir([pat ''MixedSelectivity\\*%s*PFC*Ass*.mat'']);',Targets{iTarget},Targets{iTarget}) )
    end
end

%% SINGLE UNITS batch import unit data for meta-analysis and plotting - 2-way
for iTarget = 1:length(Targets)
   
    eval(sprintf('fileList = fileList_%s;',Targets{iTarget}))
    Delays_ = DelaysList_{iTarget};
    
    % Import
    for iArea = 1:2
        for iFile =1:length(fileList)
            
                fname=strtok(fileList(iFile).name,'_');
                if useResampledTrials
                    fnIn = sprintf('%s\\MixedSelectivity\\TrialNosMatched\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
                else
                    fnIn = sprintf('%s\\MixedSelectivity\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
                end
                load(fnIn ,'D');
                
                for iDelay = 1:length(Delays_)
                    
                    % ANOVA:
                    %Fractions of tuned cells
                    eval(sprintf('D_.%s.ANOVA2.prcSig_cor{iArea}(:,iFile)   = D.%s.ANOVA2.prcSig_cor;',Delays_{iDelay},Delays_{iDelay}))
                    eval(sprintf('D_.%s.ANOVA2.prcSig_err{iArea}(:,iFile)   = D.%s.ANOVA2.prcSig_err;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: Fractions of tuned cells
                    eval(sprintf('D_.%s.ANOVA2.prcTuned_cor{iArea}(:,iFile)  = D.%s.ANOVA2.prcTuned_cor;',Delays_{iDelay},Delays_{iDelay}))
                    eval(sprintf('D_.%s.ANOVA2.prcTuned_err{iArea}(:,iFile)  = D.%s.ANOVA2.prcTuned_err;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: Thresholded p values
                    eval(sprintf('D_.%s.ANOVA2.pThresh_cor{iArea}{iFile}     = D.%s.ANOVA2.p_cor<0.05;',Delays_{iDelay},Delays_{iDelay}))
                    eval(sprintf('D_.%s.ANOVA2.pThresh_err{iArea}{iFile}     = D.%s.ANOVA2.p_err<0.05;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: F-values
                    eval(sprintf('D_.%s.ANOVA2.F_cor{iArea}{iFile}     = D.%s.ANOVA2.F_cor;',Delays_{iDelay},Delays_{iDelay}))
                    eval(sprintf('D_.%s.ANOVA2.F_err{iArea}{iFile}     = D.%s.ANOVA2.F_err;',Delays_{iDelay},Delays_{iDelay}))
                    
                    eval(sprintf('tuning   = D.%s.ANOVA2.tuning_cor;',Delays_{1}))
                    eval(sprintf('Factors   = D.%s.ANOVA2.Factors_cor;',Delays_{1}))
                   
                end
            clear D
        end
    end
    
   
    % Collapse population means/errors
    for iArea = 1:2
        for iDelay = 1:length(Delays_)
            
            eval(sprintf('D_.%s.ANOVA2.prcSig_corMean{iArea}  = nanmean(D_.%s.ANOVA2.prcSig_cor{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.prcSig_errMean{iArea}  = nanmean(D_.%s.ANOVA2.prcSig_err{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.prcSig_corSEM{iArea}   = nansem(D_.%s.ANOVA2.prcSig_cor{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.prcSig_errSEM{iArea}   = nansem(D_.%s.ANOVA2.prcSig_err{iArea},2);',Delays_{iDelay},Delays_{iDelay}))

            eval(sprintf('D_.%s.ANOVA2.prcTuned_corMean{iArea}  = nanmean(D_.%s.ANOVA2.prcTuned_cor{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.prcTuned_errMean{iArea}  = nanmean(D_.%s.ANOVA2.prcTuned_err{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.prcTuned_corSEM{iArea}   = nansem(D_.%s.ANOVA2.prcTuned_cor{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.prcTuned_errSEM{iArea}   = nansem(D_.%s.ANOVA2.prcTuned_err{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            
            
            
            eval(sprintf('D_.%s.ANOVA2.pThresh_corMean{iArea}  = cell2mat(cellfun(@(X) nanmean(X,1),D_.%s.ANOVA2.pThresh_cor{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.pThresh_errMean{iArea}  = cell2mat(cellfun(@(X) nanmean(X,1),D_.%s.ANOVA2.pThresh_err{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.pThresh_corSEM{iArea}   = cell2mat(cellfun(@(X) nansem(X,1),D_.%s.ANOVA2.pThresh_cor{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.pThresh_errSEM{iArea}   = cell2mat(cellfun(@(X) nansem(X,1),D_.%s.ANOVA2.pThresh_err{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            
            eval(sprintf('D_.%s.ANOVA2.F_corMean{iArea}  = cell2mat(cellfun(@(X) nanmean(X,1),D_.%s.ANOVA2.F_cor{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.F_errMean{iArea}  = cell2mat(cellfun(@(X) nanmean(X,1),D_.%s.ANOVA2.F_err{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.F_corSEM{iArea}   = cell2mat(cellfun(@(X) nansem(X,1),D_.%s.ANOVA2.F_cor{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA2.F_errSEM{iArea}   = cell2mat(cellfun(@(X) nansem(X,1),D_.%s.ANOVA2.F_err{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
        end
    end
    
end
%% pie chart of fractions of units
for s = 1:2
    figure('name',Areas{s})
    for iTarget = 1:3
        subplot(1,3,iTarget); hold on
        
        eval(sprintf('p_thresh = cell2mat(D_.%s.ANOVA2.pThresh_cor{s}'');',Targets_{iTarget}))
        %     p_thresh(p_thresh(:,3)==1,:)=1;
        
%         q_ = [sum(sum(p_thresh,2)==0),...              % No tuning
%             sum(sum(p_thresh == [1 0 0],2)==3),...     % CS for Context only
%             sum(sum(p_thresh == [0 1 0],2)==3),...     % CS for Location only
%             sum(sum(p_thresh == [1 1 0],2)==3),...     % LMS for Context x Location
%             sum(sum(p_thresh(:,3) == 1,2)),...   % NMS
%             ];
%         
%         pie(q_./size(p_thresh,1),{'Untuned','CS: Context','CS: Location','LMS: Context X Location','NMS: Context X Location'})
        
        q_ = [sum(sum(p_thresh,2)==0),...                                                 % No tuning
              sum(sum(p_thresh == [1 0 0],2)==3) + sum(sum(p_thresh == [0 1 0],2)==3),... % Any CS
              sum(sum(p_thresh == [1 1 0],2)==3),...                                      % LMS for Context x Location
              sum(sum(p_thresh(:,3) == 1,2)),...                                          % NMS
            ];
        
        pie(q_./size(p_thresh,1),{'Untuned','CS','LMS: Context X Location','NMS: Context X Location'})        
        axis off
    end
end
%% Plot proportions of 2-way tuned/untuned Units - group by delay
Delays_ =  fieldnames(D_);


for iArea = 1:2
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    for iDelay = 1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('D = D_.%s;',Delays_{iDelay}))
        bar(D.ANOVA2.prcTuned_corMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar(D.ANOVA2.prcTuned_corMean{iArea},D.ANOVA2.prcTuned_corSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
        title(Delays__{iDelay})
        axis([0 8 0 100])
        if iDelay==1
            ylabel('% of units')
        end    
                h=gca; h.Layer = 'top';

    end
end
clear D
%% Plot proportions of 2-way tuned/untuned Units - group by tuning
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
            m(iDelay) = D.ANOVA2.prcTuned_corMean{iArea}(iTuning);
            e(iDelay) = D.ANOVA2.prcTuned_corSEM{iArea}(iTuning);
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
%% Plot proportions of 2-way tuned/untuned Units - group by factor


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
            m(iDelay) = D.ANOVA2.prcSig_errMean{iArea}(iType);
            e(iDelay) = D.ANOVA2.prcSig_errSEM{iArea}(iType);
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
%% Plot only medium and long delays, 2-way untuned, pure and mixed tuning - correct
tuning = {'Untuned','Simple (Context)','Simple (Location)','Any Simple','Linear Mixed','Non-linear Mixed','Any Tuning'}
Delays_ =  fieldnames(D_);
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
for iArea = 1:2
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]);
    tuning__ = [1 2 3 5 6];
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
            m(iDelay) = D.ANOVA2.prcTuned_corMean{iArea}(tuning__(iTuning));
            e(iDelay) = D.ANOVA2.prcTuned_corSEM{iArea}(tuning__(iTuning));
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
%% Plot only medium and long delays, 2-way untuned, pure and mixed tuning - errors
tuning = {'Untuned','Simple (Context)','Simple (Location)','Any Simple','Linear Mixed','Non-linear Mixed','Any Tuning'};
Delays_ =  fieldnames(D_);
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
for iArea = 1:2
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]);
    tuning__ = [1 2 3 5 6];
    for iTuning = 1:length(tuning__)
        subplot(1,length(tuning__),iTuning);hold on
%         area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],100*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],100*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
%         text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
%         text(mean([2.2 6.8]),95,'Early','HorizontalAlignment','center')
%         text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D = D_.%s;',Delays_{iDelay}))
            m(iDelay) = D.ANOVA2.prcTuned_errMean{iArea}(tuning__(iTuning));
            e(iDelay) = D.ANOVA2.prcTuned_errSEM{iArea}(tuning__(iTuning));
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
%% Run stats on pure and mixed tuning, 2-way
h = [];
p=[];
    tuning__ = [1 2 3 5 6];
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
for iArea = 1:2
    
    m_medium =[];
    for iDelay = 2:4
        eval(sprintf('D = D_.%s;',Delays_{iDelay}))
        m_medium = [m_medium,D.ANOVA2.prcTuned_err{iArea}];
    end
    
    m_long =[];
    for iDelay = 5:8
        eval(sprintf('D = D_.%s;',Delays_{iDelay}))
        m_long = [m_long,D.ANOVA2.prcTuned_err{iArea}];
    end
    
    
   
    
    path(rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite'))
    
   
    for iTuning = 1:length(tuning__)
        [h(iTuning,iArea),p(iTuning,iArea)] = ttest2(m_medium(tuning__(iTuning),:),m_long(tuning__(iTuning),:))
%         [p(iTuning,iArea),h(iTuning,iArea)] = ranksum(m_medium(tuning__(iTuning),:),m_long(tuning__(iTuning),:))
    end

    path(addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite'))
end
%% batch import unit data for meta-analysis and plotting - 3-way
for iTarget = 1:length(Targets)
   
    eval(sprintf('fileList = fileList_%s;',Targets{iTarget}))
    Delays_ = DelaysList_{iTarget};
    
    % Import
    for iArea = 1:2
        for iFile =1:length(fileList)
            
                fname=strtok(fileList(iFile).name,'_');
                fnIn = sprintf('%s\\MixedSelectivity\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
                load(fnIn ,'D');
                
                for iDelay = 1:length(Delays_)
                    
                    
                    % ANOVA: Fractions of tuned cells
                    eval(sprintf('D_.%s.ANOVA3.prcTuned{iArea}(:,iFile)  = D.%s.ANOVA3.prcTuned2;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: Thresholded p values
                    eval(sprintf('D_.%s.ANOVA3.pThresh{iArea}{iFile}     = D.%s.ANOVA3.p<0.05;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: F-values
                    eval(sprintf('D_.%s.ANOVA3.F{iArea}{iFile}     = D.%s.ANOVA3.F;',Delays_{iDelay},Delays_{iDelay}))
                    
                    eval(sprintf('tuning3   = D.%s.ANOVA3.tuning2;',Delays_{1}))
                   
                end
            clear D
        end
    end
    
   
    % Collapse population means/errors
    for iArea = 1:2
        for iDelay = 1:length(Delays_)
            
            eval(sprintf('D_.%s.ANOVA3.prcTunedMean{iArea}  = nanmean(D_.%s.ANOVA3.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA3.prcTunedSEM{iArea}   = nansem(D_.%s.ANOVA3.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            
            eval(sprintf('D_.%s.ANOVA3.pThreshMean{iArea}  = cell2mat(cellfun(@nanmean,D_.%s.ANOVA3.pThresh{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA3.pThreshSEM{iArea}   = cell2mat(cellfun(@nansem,D_.%s.ANOVA3.pThresh{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            
            eval(sprintf('D_.%s.ANOVA3.FMean{iArea}  = cell2mat(cellfun(@nanmean,D_.%s.ANOVA3.F{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA3.FSEM{iArea}   = cell2mat(cellfun(@nansem,D_.%s.ANOVA3.F{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
        end
    end
    
end
%% Plot proportions of 3-way tuned/untuned Units - group by tuning
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
Delays_ =  fieldnames(D_);

for iArea = 1:2
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    for iTuning = 1:length(tuning3)
        subplot(1,length(tuning3),iTuning);hold on
        area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
        text(mean([2.2 6.8]),95,'Middle','HorizontalAlignment','center')
        text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D = D_.%s;',Delays_{iDelay}))
            m(iDelay) = D.ANOVA3.prcTunedMean{iArea}(iTuning);
            e(iDelay) = D.ANOVA3.prcTunedSEM{iArea}(iTuning);
        end
        bar([1,3,4,5,6,8,9,10],m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar([1,3,4,5,6,8,9,10],m,e,'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
        clear m e
        title(tuning3{iTuning})
        axis([0 11 0 100])
        h=gca; h.Layer = 'top';

        if iTuning==1
            ylabel('% of units')
        end
    end
end
clear D

%% ASSEMBLIES batch import Assembly data for meta-analysis and plotting - 2-way
for iTarget = 1:length(Targets)
   
    eval(sprintf('fileList = AssfileList_%s;',Targets{iTarget}))
    Delays_ = DelaysList_{iTarget};
    
    % Import
    for iArea = 1:3
        for iFile =1:length(fileList)
            
            fname=strtok(fileList(iFile).name,'_');
            if useResampledTrials
                fnIn = sprintf('%sMixedSelectivity\\TrialNosMatched\\%s_%s_MixedSelectivity_Ass.mat',pat,fname,Areas{iArea});
            else
                fnIn = sprintf('%sMixedSelectivity\\%s_%s_MixedSelectivity_Ass.mat',pat,fname,Areas{iArea});
            end
            if exist(fnIn)==2
                load(fnIn ,'D_Ass');
                
                for iDelay = 1:length(Delays_)
                    
                    % ANOVA:
                    %Fractions of tuned cells
                    eval(sprintf('D_Ass_.%s.ANOVA2.prcSig_cor{iArea}(:,iFile)   = D_Ass.%s.ANOVA2.prcSig_cor;',Delays_{iDelay},Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.ANOVA2.prcSig_err{iArea}(:,iFile)   = D_Ass.%s.ANOVA2.prcSig_err;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: Fractions of tuned cells
                    eval(sprintf('D_Ass_.%s.ANOVA2.prcTuned_cor{iArea}(:,iFile)  = D_Ass.%s.ANOVA2.prcTuned_cor;',Delays_{iDelay},Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.ANOVA2.prcTuned_err{iArea}(:,iFile)  = D_Ass.%s.ANOVA2.prcTuned_err;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: Thresholded p values
                    eval(sprintf('D_Ass_.%s.ANOVA2.pThresh_cor{iArea}{iFile}     = D_Ass.%s.ANOVA2.p_cor<0.05;',Delays_{iDelay},Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.ANOVA2.pThresh_err{iArea}{iFile}     = D_Ass.%s.ANOVA2.p_err<0.05;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: F-values
                    eval(sprintf('D_Ass_.%s.ANOVA2.F_cor{iArea}{iFile}     = D_Ass.%s.ANOVA2.F_cor;',Delays_{iDelay},Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.ANOVA2.F_err{iArea}{iFile}     = D_Ass.%s.ANOVA2.F_err;',Delays_{iDelay},Delays_{iDelay}))
                    
                    eval(sprintf('tuning   = D_Ass.%s.ANOVA2.tuning_cor;',Delays_{1}))
                    eval(sprintf('Factors   = D_Ass.%s.ANOVA2.Factors_cor;',Delays_{1}))
                    
                end
                clear D_Ass
                
            else
                for iDelay = 1:length(Delays_)
                    
                    % ANOVA:
                    %Fractions of tuned cells
                    eval(sprintf('D_Ass_.%s.ANOVA2.prcSig_cor{iArea}(:,iFile)   = nan(1,3);',Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.ANOVA2.prcSig_err{iArea}(:,iFile)   = nan(1,3);',Delays_{iDelay}))
                    % ANOVA: Fractions of tuned cells
                    eval(sprintf('D_Ass_.%s.ANOVA2.prcTuned_cor{iArea}(:,iFile)  = nan(1,7);',Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.ANOVA2.prcTuned_err{iArea}(:,iFile)  = nan(1,7);',Delays_{iDelay}))
                    % ANOVA: Thresholded p values
                    eval(sprintf('D_Ass_.%s.ANOVA2.pThresh_cor{iArea}{iFile}     = nan(1,3);',Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.ANOVA2.pThresh_err{iArea}{iFile}     = nan(1,3);',Delays_{iDelay}))
                    % ANOVA: F-values
                    eval(sprintf('D_Ass_.%s.ANOVA2.F_cor{iArea}{iFile}     = nan(1,3);',Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.ANOVA2.F_err{iArea}{iFile}     = nan(1,3);',Delays_{iDelay}))
                    
                    %                     eval(sprintf('tuning   = D_Ass.%s.ANOVA2.tuning_cor;',Delays_{1}))
                    %                     eval(sprintf('Factors   = D_Ass.%s.ANOVA2.Factors_cor;',Delays_{1}))
                    
                end
            end
        end
    end
    
    % Collapse population means/errors
    for iArea = 1:3
        for iDelay = 1:length(Delays_)
            
            eval(sprintf('D_Ass_.%s.ANOVA2.prcSig_corMean{iArea}  = nanmean(D_Ass_.%s.ANOVA2.prcSig_cor{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.prcSig_errMean{iArea}  = nanmean(D_Ass_.%s.ANOVA2.prcSig_err{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.prcSig_corSEM{iArea}   = nansem(D_Ass_.%s.ANOVA2.prcSig_cor{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.prcSig_errSEM{iArea}   = nansem(D_Ass_.%s.ANOVA2.prcSig_err{iArea},2);',Delays_{iDelay},Delays_{iDelay}))

            eval(sprintf('D_Ass_.%s.ANOVA2.prcTuned_corMean{iArea}  = nanmean(D_Ass_.%s.ANOVA2.prcTuned_cor{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.prcTuned_errMean{iArea}  = nanmean(D_Ass_.%s.ANOVA2.prcTuned_err{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.prcTuned_corSEM{iArea}   = nansem(D_Ass_.%s.ANOVA2.prcTuned_cor{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.prcTuned_errSEM{iArea}   = nansem(D_Ass_.%s.ANOVA2.prcTuned_err{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            
            
            eval(sprintf('D_Ass_.%s.ANOVA2.pThresh_corMean{iArea}  = cell2mat(cellfun(@(X) nanmean(X,1),D_Ass_.%s.ANOVA2.pThresh_cor{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.pThresh_errMean{iArea}  = cell2mat(cellfun(@(X) nanmean(X,1),D_Ass_.%s.ANOVA2.pThresh_err{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.pThresh_corSEM{iArea}   = cell2mat(cellfun(@(X) nansem(X,1),D_Ass_.%s.ANOVA2.pThresh_cor{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.pThresh_errSEM{iArea}   = cell2mat(cellfun(@(X) nansem(X,1),D_Ass_.%s.ANOVA2.pThresh_err{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            
            eval(sprintf('D_Ass_.%s.ANOVA2.F_corMean{iArea}  = cell2mat(cellfun(@(X) nanmean(X,1),D_Ass_.%s.ANOVA2.F_cor{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.F_errMean{iArea}  = cell2mat(cellfun(@(X) nanmean(X,1),D_Ass_.%s.ANOVA2.F_err{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.F_corSEM{iArea}   = cell2mat(cellfun(@(X) nansem(X,1),D_Ass_.%s.ANOVA2.F_cor{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA2.F_errSEM{iArea}   = cell2mat(cellfun(@(X) nansem(X,1),D_Ass_.%s.ANOVA2.F_err{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
        end
    end
    
end
%% pie chart of fractions of Assemblies
for s = 1:3
    figure('name',Areas{s})
    for iTarget = 1:3
        subplot(1,3,iTarget); hold on
        
        eval(sprintf('p_thresh = cell2mat(cellfun(@double,D_Ass_.%s.ANOVA2.pThresh_cor{s},''UniformOUtput'',false)'');',Targets_{iTarget}))
        %     p_thresh(p_thresh(:,3)==1,:)=1;
        p_thresh(isnan(p_thresh(:,1)),:)=[];
        q_ = [sum(sum(p_thresh,2)==0),...              % No tuning
            sum(sum(p_thresh == [1 0 0],2)==3),...     % CS for Context only
            sum(sum(p_thresh == [0 1 0],2)==3),...     % CS for Location only
            sum(sum(p_thresh == [1 1 0],2)==3),...     % LMS for Context x Location
            sum(sum(p_thresh(:,3) == 1,2)),...         % NMS
            ];
        pie(q_./size(p_thresh,1),{'Untuned','CS: Context','CS: Location','LMS: Context X Location','NMS: Context X Location'})
        
%         q_ = [sum(sum(p_thresh,2)==0),...                                                 % No tuning
%               sum(sum(p_thresh == [1 0 0],2)==3) + sum(sum(p_thresh == [0 1 0],2)==3),... % Any CS
%               sum(sum(p_thresh == [1 1 0],2)==3),...                                      % LMS for Context x Location
%               sum(sum(p_thresh(:,3) == 1,2)),...                                          % NMS
%             ];
%         pie(q_./size(p_thresh,1),{'Untuned','CS','LMS: Context X Location','NMS: Context X Location'})        
        axis off
    end
end
%% Plot proportions of 2-way tuned/untuned Units - group by tuning
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};

for iArea = 1:3
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
            eval(sprintf('D_Ass = D_Ass_.%s;',Delays_{iDelay}))
            m(iDelay) = D_Ass.ANOVA2.prcTuned_corMean{iArea}(iTuning);
            e(iDelay) = D_Ass.ANOVA2.prcTuned_corSEM{iArea}(iTuning);
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
%% Plot only medium and long delays, 2-way untuned, pure and mixed tuning - correct
tuning = {'Untuned','Simple (Context)','Simple (Location)','Any Simple','Linear Mixed','Non-linear Mixed','Any Tuning'}
Delays_ =  fieldnames(D_Ass_);
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
for iArea = 1:3
    figure('name',['Percentage of tuned Assemblies: ' Areas{iArea}]);
    tuning__ = [1 2 3 5 6];
    for iTuning = 1:length(tuning__)
        subplot(1,length(tuning__),iTuning);hold on
%         area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
%         text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
        text(mean([2.2 6.8]),95,'Early','HorizontalAlignment','center')
        text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D_Ass = D_Ass_.%s;',Delays_{iDelay}))
            m(iDelay) = D_Ass.ANOVA2.prcTuned_corMean{iArea}(tuning__(iTuning));
            e(iDelay) = D_Ass.ANOVA2.prcTuned_corSEM{iArea}(tuning__(iTuning));
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
%% Plot only medium and long delays, 2-way untuned, pure and mixed tuning - errors
tuning = {'Untuned','Simple (Context)','Simple (Location)','Any Simple','Linear Mixed','Non-linear Mixed','Any Tuning'};
Delays_ =  fieldnames(D_Ass_);
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
for iArea = 1:3
    figure('name',['Percentage of tuned units: ' Areas{iArea}]);
    tuning__ = [1 2 3 5 6];
    for iTuning = 1:length(tuning__)
        subplot(1,length(tuning__),iTuning);hold on
%         area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],100*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],100*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
%         text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
%         text(mean([2.2 6.8]),95,'Early','HorizontalAlignment','center')
%         text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D_Ass = D_Ass_.%s;',Delays_{iDelay}))
            m(iDelay) = D_Ass.ANOVA2.prcTuned_errMean{iArea}(tuning__(iTuning));
            e(iDelay) = D_Ass.ANOVA2.prcTuned_errSEM{iArea}(tuning__(iTuning));
        end
        
        bar([1,3,4,5,6,8,9,10],m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar([1,3,4,5,6,8,9,10],m,e,'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
        clear m e
        title(tuning{tuning__(iTuning)})
        axis([2 11 0 100])
        h=gca; h.Layer = 'top';

        if iTuning==1
            ylabel('% of assemblies')
        end
        
        
    end
end
%% batch import unit data for meta-analysis and plotting - 3-way
for iTarget = 1:length(Targets)
   
    eval(sprintf('fileList = AssfileList_%s;',Targets{iTarget}))
    Delays_ = DelaysList_{iTarget};
    
    % Import
    for iArea = 1:3
        for iFile =1:length(fileList)
            
            fname=strtok(fileList(iFile).name,'_');
            if useResampledTrials
                fnIn = sprintf('%sMixedSelectivity\\TrialNosMatched\\%s_%s_MixedSelectivity_Ass.mat',pat,fname,Areas{iArea});
            else
                fnIn = sprintf('%sMixedSelectivity\\%s_%s_MixedSelectivity_Ass.mat',pat,fname,Areas{iArea});
            end
             if exist(fnIn)==2
                load(fnIn ,'D_Ass');
                for iDelay = 1:length(Delays_)
                    % ANOVA: Fractions of tuned cells
                    eval(sprintf('D_Ass_.%s.ANOVA3.prcTuned{iArea}(:,iFile)  = D_Ass.%s.ANOVA3.prcTuned2;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: Thresholded p values
                    eval(sprintf('D_Ass_.%s.ANOVA3.pThresh{iArea}{iFile}     = D_Ass.%s.ANOVA3.p<0.05;',Delays_{iDelay},Delays_{iDelay}))
                    % ANOVA: F-values
                    eval(sprintf('D_Ass_.%s.ANOVA3.F{iArea}{iFile}           = D_Ass.%s.ANOVA3.F;',Delays_{iDelay},Delays_{iDelay}))
                    eval(sprintf('tuning3   = D_Ass.%s.ANOVA3.tuning2;',Delays_{1}))
                end
             else
                for iDelay = 1:length(Delays_)
                    % ANOVA: Fractions of tuned cells
                    eval(sprintf('D_Ass_.%s.ANOVA3.prcTuned{iArea}(:,iFile)  = nan(1,15);',Delays_{iDelay}))
                    % ANOVA: Thresholded p values
                    eval(sprintf('D_Ass_.%s.ANOVA3.pThresh{iArea}{iFile}     = nan(1,7);',Delays_{iDelay}))
                    % ANOVA: F-values
                    eval(sprintf('D_Ass_.%s.ANOVA3.F{iArea}{iFile}           = nan(1,7);,',Delays_{iDelay}))
                    
                end                
             end
            clear D
        end
    end
    
   
    % Collapse population means/errors
    for iArea = 1:3
        for iDelay = 1:length(Delays_)
            
            eval(sprintf('D_Ass_.%s.ANOVA3.prcTunedMean{iArea}  = nanmean(D_Ass_.%s.ANOVA3.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA3.prcTunedSEM{iArea}   = nansem(D_Ass_.%s.ANOVA3.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
            
            eval(sprintf('D_Ass_.%s.ANOVA3.pThreshMean{iArea}  = cell2mat(cellfun(@(x) nanmean(x,1),D_Ass_.%s.ANOVA3.pThresh{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA3.pThreshSEM{iArea}   = cell2mat(cellfun(@(x) nansem(x,1),D_Ass_.%s.ANOVA3.pThresh{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            
            eval(sprintf('D_Ass_.%s.ANOVA3.FMean{iArea}  = cell2mat(cellfun(@(x) nanmean(x,1),D_Ass_.%s.ANOVA3.F{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_Ass_.%s.ANOVA3.FSEM{iArea}   = cell2mat(cellfun(@(x) nansem(x,1),D_Ass_.%s.ANOVA3.F{iArea},''UniformOutput'',false)'');',Delays_{iDelay},Delays_{iDelay}))
        end
    end
    
end
%% Plot proportions of 3-way tuned/untuned Assemblies - group by tuning
Delays__ = {'0s';'2s';'4s';'6s';'8s';'4s';'8s';'16s'};
Delays_ =  fieldnames(D_Ass_);

for iArea = 1:3
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    for iTuning = 1:length(tuning3)
        subplot(1,length(tuning3),iTuning);hold on
        area([0.2 1.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([2.2 6.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        area([7.2 10.8],90*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
        text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
        text(mean([2.2 6.8]),95,'Middle','HorizontalAlignment','center')
        text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')
        for iDelay = 1:length(Delays_)
            eval(sprintf('D_Ass = D_Ass_.%s;',Delays_{iDelay}))
            m(iDelay) = D_Ass.ANOVA3.prcTunedMean{iArea}(iTuning);
            e(iDelay) = D_Ass.ANOVA3.prcTunedSEM{iArea}(iTuning);
        end
        bar([1,3,4,5,6,8,9,10],m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
        errorbar([1,3,4,5,6,8,9,10],m,e,'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
        set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
        clear m e
        title(tuning3{iTuning})
        axis([0 11 0 100])
        h=gca; h.Layer = 'top';

        if iTuning==1
            ylabel('% of units')
        end
    end
end
clear D

clear
pat = 'C:\Analysis\AssemblyAnalysis\raw\';

% fileList_SHORT = dir([pat 'MixedSelectivity_SHORT' filesep '*SHORT*.mat']);
% fileList_MEDIUM = dir([pat 'MixedSelectivity_MEDIUM' filesep '*MEDIUM*.mat']);
% fileList_LONG = dir([pat 'MixedSelectivity' filesep '*LONG*','*PFC.mat']);

fileList_SHORT = dir([pat 'allTimestamps' filesep '*SHORT*.mat']);
fileList_MEDIUM = dir([pat 'allTimestamps' filesep '*MEDIUM*.mat']);
fileList_LONG = dir([pat 'allTimestamps' filesep '*LONG*.mat']);
% fileList_LONG(end)=[]
Areas = {'PFC','HP'};
bw = 0.05;
sigSpanThreshold = 0.2;
%% batch import unit data for meta-analysis and plotting - SHORT
% Import
clear D_
fileList = fileList_SHORT;
Delays_ = {'Delay_0'};

for iArea = 1:2%length(Areas)
    
    
    for iFile =1:length(fileList)
        try
            fname=strtok(fileList(iFile).name,'_');
            fnIn = sprintf('%s\\MixedSelectivity_SHORT\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
            load(fnIn ,'D');
            
            for iDelay = 1:length(Delays_)
                % T-score for L vs R decoding
                eval(sprintf('D_.%s.LR.avgFR{iArea}{iFile}  = mean(D.%s.LR.avgFR{iArea},1);',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TS{iArea}{iFile}     = D.%s.LR.TS;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TSsig{iArea}{iFile}  = D.%s.LR.TSsig;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TS_drawnTrials{iArea}{iFile}  = D.%s.LR.TS_drawnTrials;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TSsig_drawnTrials{iArea}{iFile}  = D.%s.LR.TSsig_drawnTrials;',Delays_{iDelay},Delays_{iDelay}))

                eval(sprintf('D_.%s.LR.CVEindividual{iArea}{iFile}  = 1-D.%s.LR.CVEindividual;',Delays_{iDelay},Delays_{iDelay}))
            end
            
%             clear D
        end
    end
end

% Collapse population means/errors
for iArea = 1:2
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_.%s.LR.avgFR_collapse{iArea}      = cell2mat(D_.%s.LR.avgFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.TSpeak_collapse{iArea}     = cell2mat(cellfun(@max,D_.%s.LR.TS{iArea},''UniformOutput'',false));',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.TSsigSpan_collapse{iArea}  = cell2mat(cellfun(@sum,D_.%s.LR.TSsig{iArea},''UniformOutput'',false))*bw ;',Delays_{iDelay},Delays_{iDelay}))
        
       
        eval(sprintf('D_.%s.LR.maxCVE_collapse{iArea}     = cellfun(@(A) max(A,[],2),D_.%s.LR.CVEindividual{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.meanCVE_collapse{iArea}    = cellfun(@(A) mean(A,2),D_.%s.LR.CVEindividual{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.sigSpanCVE_collapse{iArea} = cellfun(@(A) sum(A>sigSpanThreshold,2)*bw,D_.%s.LR.CVEindividual{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))        
        
        eval(sprintf('D_.%s.LR.TSpeak_collapse_drawnTrials{iArea}     = cell2mat(cellfun(@max,D_.%s.LR.TS_drawnTrials{iArea},''UniformOutput'',false));',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.TSsigSpan_drawnTrials{iArea}  = cellfun(@(x) nansum(x>sigSpanThreshold)*bw,D_.%s.LR.TSsig_drawnTrials{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))

    end
    
end

% clear D Stability D_temp D_tempL D_tempR B bp
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
                % T-score for L vs R decoding
                eval(sprintf('D_.%s.LR.avgFR{iArea}{iFile}  = mean(D.%s.LR.avgFR{iArea},1);',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TS{iArea}{iFile}     = D.%s.LR.TS;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TSsig{iArea}{iFile}  = D.%s.LR.TSsig;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TS_drawnTrials{iArea}{iFile}  = D.%s.LR.TS_drawnTrials;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TSsig_drawnTrials{iArea}{iFile}  = D.%s.LR.TSsig_drawnTrials;',Delays_{iDelay},Delays_{iDelay}))
                
                eval(sprintf('D_.%s.LR.CVEindividual{iArea}{iFile}  = 1-D.%s.LR.CVEindividual;',Delays_{iDelay},Delays_{iDelay}))
            end
            
            clear D
        end
    end
end

% Collapse population means/errors
for iArea = 1:2
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_.%s.LR.avgFR_collapse{iArea}      = cell2mat(D_.%s.LR.avgFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.TSpeak_collapse{iArea}     = cell2mat(cellfun(@max,D_.%s.LR.TS{iArea},''UniformOutput'',false));',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.TSsigSpan_collapse{iArea}  = cell2mat(cellfun(@sum,D_.%s.LR.TSsig{iArea},''UniformOutput'',false))*bw ;',Delays_{iDelay},Delays_{iDelay}))
        
        eval(sprintf('D_.%s.LR.maxCVE_collapse{iArea}     = cellfun(@(A) max(A,[],2),D_.%s.LR.CVEindividual{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.meanCVE_collapse{iArea}    = cellfun(@(A) mean(A,2),D_.%s.LR.CVEindividual{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.sigSpanCVE_collapse{iArea} = cellfun(@(A) sum(A>sigSpanThreshold,2)*bw,D_.%s.LR.CVEindividual{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))        
        
        eval(sprintf('D_.%s.LR.TSpeak_collapse_drawnTrials{iArea}     = cell2mat(cellfun(@max,D_.%s.LR.TS_drawnTrials{iArea},''UniformOutput'',false));',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.TSsigSpan_drawnTrials{iArea}  = cellfun(@(x) nansum(x>sigSpanThreshold)*bw,D_.%s.LR.TSsig_drawnTrials{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))

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
            fnIn = sprintf('%s\\MixedSelectivity_LONG\\MixedSelectivity\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
            load(fnIn ,'D');
            
            for iDelay = 1:length(Delays_)
                % T-score for L vs R decoding
                eval(sprintf('D_.%s.LR.avgFR{iArea}{iFile}  = mean(D.%s.LR.avgFR{iArea},1);',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TS{iArea}{iFile}     = D.%s.LR.TS;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TSsig{iArea}{iFile}  = D.%s.LR.TSsig;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TS_drawnTrials{iArea}{iFile}  = D.%s.LR.TS_drawnTrials;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.LR.TSsig_drawnTrials{iArea}{iFile}  = D.%s.LR.TSsig_drawnTrials;',Delays_{iDelay},Delays_{iDelay}))
                
                eval(sprintf('D_.%s.LR.CVEindividual{iArea}{iFile}  = 1-D.%s.LR.CVEindividual;',Delays_{iDelay},Delays_{iDelay}))
            end
            
            clear D
        end
    end
end

% Collapse population means/errors
for iArea = 1:2
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_.%s.LR.avgFR_collapse{iArea}      = cell2mat(D_.%s.LR.avgFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.TSpeak_collapse{iArea}     = cell2mat(cellfun(@max,D_.%s.LR.TS{iArea},''UniformOutput'',false));',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.TSsigSpan_collapse{iArea}  = cell2mat(cellfun(@sum,D_.%s.LR.TSsig{iArea},''UniformOutput'',false))*bw ;',Delays_{iDelay},Delays_{iDelay}))
%         eval(sprintf('D_.%s.LR.maxCVE_collapse{iArea}     = cellfun(@(A) max(A,[],2),D_.%s.LR.CVEindividual{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))
%         eval(sprintf('D_.%s.LR.meanCVE_collapse{iArea}    = cellfun(@(A) mean(A,2),D_.%s.LR.CVEindividual{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))
%         eval(sprintf('D_.%s.LR.sigSpanCVE_collapse{iArea} = cellfun(@(A) sum(A>sigSpanThreshold,2)*bw,D_.%s.LR.CVEindividual{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))        
        eval(sprintf('D_.%s.LR.TSpeak_collapse_drawnTrials{iArea}     = cell2mat(cellfun(@max,D_.%s.LR.TS_drawnTrials{iArea},''UniformOutput'',false));',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.LR.TSsigSpan_drawnTrials{iArea}  = cellfun(@(x) nansum(x>sigSpanThreshold)*bw,D_.%s.LR.TSsig_drawnTrials{iArea},''UniformOutput'',false);',Delays_{iDelay},Delays_{iDelay}))

    end
end

clear D Stability D_temp D_tempL D_tempR B bp
%% Plot peak T score histograms
bins = 0:0.1:20;
conditions = fieldnames(D_);

Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
% col_ = [0 0 0;hot(4);cool(3)];
col_ = [0 0 0;jet(7)];
figure;
for iArea = 1:2
    
    subplot(1,2,iArea);hold on
    for iDelay = 1:length(conditions)
        eval(sprintf('avgFR = D_.%s.LR.avgFR{iArea};',conditions{iDelay}))
        eval(sprintf('TS    = cellfun(@max,D_.%s.LR.TS_drawnTrials{iArea},''UniformOutput'',false);',conditions{iDelay}))    
%         eval(sprintf('TS    = cellfun(@max,D_.%s.LR.TS{iArea},''UniformOutput'',false);',conditions{iDelay}))
        h_ = zeros(length(avgFR),length(bins));
        for iFile =1:length(avgFR)
            TS_ = TS{iFile};
            avgFR_ = avgFR{iFile};
            TS_(avgFR_<0.5)=[];
            try
             h_(iFile,:) = cumsum(histc(TS_,bins)./numel(TS_));
            catch
                h_(iFile,:) = nan(1,length(bins));
            end
        end
       
%         ciplot(nanmean(h_)+nansem(h_),nanmean(h_)-nansem(h_),bins,col_(iDelay,:),0.6)
               plot(bins,nanmean(h_),'Color',col_(iDelay,:),'LineWidth',1.5)

    end
    title([Areas{iArea},' units'])
    axis([min(bins) max(bins) 0 1])
%     xlabel('Peak T-score (L vs R cue discrimination)');
    if iArea==1
        
        ylabel('Fraction of units')
    end
end
legend(Delays__,'Location','southeast');legend boxoff
%% Plot peak T score histograms - with shaded bootstrap CIs
bins = 0:0.2:20;
clear h_ x i CI
conditions = fieldnames(D_);
nDraws = 5000;
Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
% col_ = [0 0 0;hot(4);cool(3)];
col_ = [0 0 0;jet(7)];
figure;
for iArea = 1:2
    
    for iDelay = 1:length(conditions)
        eval(sprintf(' n(iDelay) = length(D_.%s.LR.avgFR_collapse{iArea});',conditions{iDelay}))
    end
    N = ceil(0.8*min(n));
    
    
    subplot(1,2,iArea);hold on
    for iDelay = 1:length(conditions)
        eval(sprintf('avgFR = D_.%s.LR.avgFR_collapse{iArea};',conditions{iDelay}))
%         eval(sprintf('TS    = D_.%s.LR.TSpeak_collapse{iArea};',conditions{iDelay}))
        eval(sprintf('TS    = D_.%s.LR.TSpeak_collapse_drawnTrials{iArea};',conditions{iDelay}))
        %TS(avgFR<0.5)=NaN;
        h_=zeros(nDraws,length(bins));
        for iDraw = 1:nDraws
            i = randsample(n(iDelay),N);
            x = TS(i); x(avgFR(i)<0.5)=[];
            h_(iDraw,:) = cumsum(histc(x,bins)./numel(x));
        end
        for iBin = 1:length(bins)
            CI(:,iBin)=prctile(h_(:,iBin),[5,95]);
        end
%         ciplot(nanmean(h_)+nansem(h_),nanmean(h_)-nansem(h_),bins,col_(iDelay,:),0.6)
               ciplot(CI(1,:),CI(2,:),bins,col_(iDelay,:),0.6)
        %        plot(bins,h_,'Color',col_(iDelay,:),'LineWidth',1.5)
    end
    title([Areas{iArea},' units'])
    axis([min(bins) max(bins) 0 1])
    %     xlabel('Peak T-score (L vs R cue discrimination)');
    if iArea==1
        ylabel('Fraction of units')
    end
end
legend(Delays__,'Location','southeast');legend boxoff
%% Plot significant span of decoding histograms
bins = bw:0.5:20;
conditions = fieldnames(D_);
chooseLongest = true; %(if true, pick the longest significant block of decoding, else sum all)

Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
% col_ = [0 0 0;hot(4);cool(3)];
col_ = [0 0 0;jet(7)];
figure;
for iArea = 1:2
    
    subplot(1,2,iArea);hold on
    for iDelay = 1:length(conditions)
        eval(sprintf('avgFR = D_.%s.LR.avgFR{iArea};',conditions{iDelay}))
        eval(sprintf('TS    = D_.%s.LR.TSsig{iArea};',conditions{iDelay}))
%         eval(sprintf('TS    = D_.%s.LR.TSsig_drawnTrials{iArea};',conditions{iDelay}))
        h_ = zeros(length(avgFR),length(bins));
        for iFile =1:length(avgFR)
            %             TS_ = TS{iFile}*bw;
            TS_ = TS{iFile};
            avgFR_ = avgFR{iFile};
            TS_(:,avgFR_<0.5)=[];
            %             TS_(:,TS_==0)=[];
            
                for i=1:size(TS_,2)
                    if chooseLongest
                        [~,~,w,~] = findpeaks(double(TS_(:,i)));
                        if ~isempty(w)
                            a(i) = max(w)*bw;
                        else
                            a(i) = 0;
                        end
                    
                    else
                         a(i) = nansum(TS_(:,i),1).*bw
                    end
                end
                    try
                        h_(iFile,:) = cumsum((histc(a,bins))./length(a(a>0)));
                    catch
                        h_(iFile,:) = nan(1,length(bins));
                    end
                
        end
        
        %         ciplot(nanmean(h_)+nansem(h_),nanmean(h_)-nansem(h_),bins,col_(iDelay,:),0.6)
        plot(bins,nanmean(h_),'Color',col_(iDelay,:),'LineWidth',1.5)
        
    end
    title([Areas{iArea},' units'])
    axis([min(bins) max(bins) 0 1])
%      xlabel('Significant span of decoding (s)');
    if iArea==1
        
        ylabel('Fraction of units')
    end
end
legend(Delays__,'Location','southeast');legend boxoff
%% Plot significant span of decoding histograms - with shaded bootstrap CIs
bins = 0:0.1:30;
clear h_ x i CI
conditions = fieldnames(D_);
chooseLongest = true; %(if true, pick the longest significant block of decoding, else sum all)
nDraws = 5000;
Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
% col_ = [0 0 0;hot(4);cool(3)];
col_ = [0 0 0;jet(7)];
figure;
for iArea = 1:2
    
    for iDelay = 1:length(conditions)
        eval(sprintf(' n(iDelay) = length(D_.%s.LR.avgFR_collapse{iArea});',conditions{iDelay}))
    end
    N = ceil(0.8*min(n));
    
    
    subplot(1,2,iArea);hold on
    for iDelay = 1:length(conditions)
        eval(sprintf('avgFR = D_.%s.LR.avgFR_collapse{iArea};',conditions{iDelay}))
%         eval(sprintf('TS    = D_.%s.LR.TSsig{iArea};',conditions{iDelay}))
        eval(sprintf('TS    = cell2mat(D_.%s.LR.TSsigSpan_drawnTrials{iArea});',conditions{iDelay}))
        TS(avgFR<0.5)=NaN;
        h_=zeros(nDraws,length(bins));
        for iDraw = 1:nDraws
            i = randsample(n(iDelay),N);
            x = TS(i); x(avgFR(i)<0.5)=[];
%             x(x==0)=[];
            h_(iDraw,:) = cumsum(histc(x,bins)./numel(x));
        end
        for iBin = 1:length(bins)
            CI(:,iBin)=prctile(h_(:,iBin),[5,95]);
        end
        %        ciplot(nanmean(h_)+nansem(h_),nanmean(h_)-nansem(h_),bins,col_(iDelay,:),0.6)
        ciplot(CI(1,:),CI(2,:),bins,col_(iDelay,:),0.6)
%                plot(bins,nanmean(h_),'Color',col_(iDelay,:),'LineWidth',1.5)
    end
    title([Areas{iArea},' units'])
    axis([min(bins) max(bins) 0 1])
    xlabel('Significant span of decoding (s)');
    if iArea==1
        ylabel('Fraction of units')
    end
end
legend(Delays__,'Location','southeast');legend boxoff
%% Plot bars of significant span of decoding >1s
bins = 0:0.1:15;
conditions = fieldnames(D_);

Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
% col_ = [0 0 0;hot(4);cool(3)];
col_ = [0 0 0;jet(7)];
clear h_
figure;
for iArea = 1:2
    
   
    for iDelay = 1:length(conditions)
        eval(sprintf('avgFR = D_.%s.LR.avgFR{iArea};',conditions{iDelay}))
        eval(sprintf('TS    = cellfun(@sum,D_.%s.LR.TSsig{iArea},''UniformOutput'',false);',conditions{iDelay}))
        for iFile =1:length(avgFR)
            TS_ = TS{iFile}*bw;
            avgFR_ = avgFR{iFile};
            TS_(avgFR_<0.05)=[];
            
            try
                h_{iDelay}(iFile) = sum(TS_>1)./length(TS_);
            catch
                h_{iDelay}(iFile) = nan(1,length(bins));
            end
        end
    end
    
    subplot(1,2,iArea);hold on
    area([0.2 1.8],0.9*[1 1],'FaceColor',[0.95 0.95 0.95],'LineStyle','none')
    area([2.2 6.8],0.9*[1 1],'FaceColor',[0.95 0.95 0.95],'LineStyle','none')
    area([7.2 10.8],0.9*[1 1],'FaceColor',[0.95 0.95 0.95],'LineStyle','none')
    text(mean([0.2 1.8]),0.95,'Early','HorizontalAlignment','center')
    text(mean([2.2 6.8]),0.95,'Middle','HorizontalAlignment','center')
    text(mean([7.2 10.8]),0.95,'Late','HorizontalAlignment','center')
    m = cellfun(@nanmean,h_);
    e = cellfun(@nansem,h_);
    x = [1,3,4,5,6,8,9,10];
    for i=1:length(x)
        bar(x(i),m(i),'EdgeColor',col_(i,:),'FaceColor',col_(i,:),'FaceAlpha',0.6)
        errorbar(x(i),m(i),e(i),'LineStyle','none','LineWidth',1.5,'Color',col_(i,:))
    end
    set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTicklabel',Delays__,'XTickLabelRotation',50)
    
    title([Areas{iArea},' units'])
    clear m e
    axis([0 11 0 1])
    h=gca; h.Layer = 'top';
    
    if iArea==1
        ylabel('Fraction of units with siginificant L/R decoding')
    end
end
%% Plot peak single-cell CVE decoding
bins = 0:0.05:1;
conditions = fieldnames(D_);

Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
% col_ = [0 0 0;hot(4);cool(3)];
col_ = [0 0 0;jet(7)];
figure;
for iArea = 1:2
    
    subplot(1,2,iArea);hold on
    for iDelay = 1:length(conditions)
        eval(sprintf('avgFR = D_.%s.LR.avgFR{iArea};',conditions{iDelay}))
        eval(sprintf('TS    = D_.%s.LR.maxCVE_collapse{iArea};',conditions{iDelay}))
        h_ = zeros(length(avgFR),length(bins));
        for iFile =1:length(avgFR)
            TS_ = TS{iFile};
            avgFR_ = avgFR{iFile};
            TS_(avgFR_<0.5)=[];
            try
             h_(iFile,:) = cumsum(histc(TS_,bins)./numel(TS_));
            catch
                h_(iFile,:) = nan(1,length(bins));
            end
        end
       
%         ciplot(nanmean(h_)+nansem(h_),nanmean(h_)-nansem(h_),bins,col_(iDelay,:),0.6)
               plot(bins,nanmean(h_),'Color',col_(iDelay,:),'LineWidth',1.5)

    end
    title([Areas{iArea},' units'])
    axis([min(bins) max(bins) 0 1])
%     xlabel('Peak T-score (L vs R cue discrimination)');
    if iArea==1
        
        ylabel('Fraction of units')
    end
end
legend(Delays__,'Location','southeast');legend boxoff
%% Plot mean single-cell CVE decoding
bins = 0:0.2:20;
conditions = fieldnames(D_);

Delays__ = {'Early: No Delay';'Mid: 2s';'Mid: 4s';'Mid: 6s';'Mid: 8s';'Late: 4s';'Late: 8s';'Late: 16s'};
% col_ = [0 0 0;hot(4);cool(3)];
col_ = [0 0 0;jet(7)];
figure;
for iArea = 1:2
    
    subplot(1,2,iArea);hold on
    for iDelay = 1:length(conditions)
        eval(sprintf('avgFR = D_.%s.LR.avgFR{iArea};',conditions{iDelay}))
        eval(sprintf('TS    = D_.%s.LR.sigSpanCVE_collapse{iArea};',conditions{iDelay}))
        h_ = zeros(length(avgFR),length(bins));
        for iFile =1:length(avgFR)
            TS_ = TS{iFile};
            avgFR_ = avgFR{iFile};
%             TS_(avgFR_<0.5)=[];
            try
             h_(iFile,:) = cumsum(histc(TS_,bins)./numel(TS_));
            catch
                h_(iFile,:) = nan(1,length(bins));
            end
        end
       
        ciplot(nanmean(h_)+nansem(h_),nanmean(h_)-nansem(h_),bins,col_(iDelay,:),0.6)
%                plot(bins,nanmean(h_),'Color',col_(iDelay,:),'LineWidth',1.5)

    end
    title([Areas{iArea},' units'])
    axis([min(bins) max(bins) 0 1])
%     xlabel('Peak T-score (L vs R cue discrimination)');
    if iArea==1
        
        ylabel('Fraction of units')
    end
end
legend(Delays__,'Location','southeast');legend boxoff

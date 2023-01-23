%% %%%%%% PREAMBLE %%%%%%
targets = {'LONG2','AM251','CP55940','URB597','Physostigmine','Scopolamine','HU210','CPphyso','CPAM'};
targets2 = {'Vehicle','AM251','CP55940','URB597','Physostigmine','Scopolamine','HU210','CP55940 & Physostigmine','CP55940 & AM251'}; 
targets_ = {'Control','AM251','CP55940','URB597','Physostigmine','Scopolamine','HU210',{'CP55940 &';'Physostigmine'},{'CP55940';'& AM251'}};

% targets = {'LONG','AM251','CP55940','Scopolamine','Physostigmine','URB597'};
% targets2 = {'Vehicle','AM251','CP55940','Scopolamine','Physostigmine','URB597'}; 
% targets_ = {'Control','AM251','CP55940','Scopolamine','Physostigmine','URB597'};


for iTarget = 1:length(targets)
    Target = targets{iTarget};
    Delays_ = {'Short','Medium','Long'};
    Delays__ = {'4s','8s','16s'};
    
    
    
    plotOnline = false;
    
    warning ('off')
    if ispc
        pat = 'C:\Analysis\AssemblyAnalysis\raw\';
        cd(pat)
        fileList=dir(sprintf('allTimestamps\\*%s*.mat',Target));
    else ismac
        pat = '/Volumes/HDD2/DNMTP/raw/';
        cd(pat)
        fileList=dir(sprintf('allTimestamps/*%s*.mat',Target));
    end
    
    
    for iFile =1:length(fileList)
        
        fname=strtok(fileList(iFile).name,'_');
        disp(sprintf('Analysing %s: %s...',Target,fname))
        load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
        
        for iDelay=1:length(Delays_)

            Delay_ = Delays_{iDelay};
            
            % Sample Latency
            eval(sprintf('Slatency.%s.Correct{iFile} =  [(t.%s.SamplePress_LeftCorrect-t.%s.CueLight_LeftCorrect)/1e6;(t.%s.SamplePress_RightCorrect-t.%s.CueLight_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
            
            % Nosepoke latency
            eval(sprintf('NPlatency.%s.Correct{iFile} =  [(t.%s.NosePoke_LeftCorrect-t.%s.DelayEnd_LeftCorrect)/1e6;(t.%s.NosePoke_RightCorrect-t.%s.DelayEnd_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
            
            % Response latency
            % eval(sprintf('Rlatency.%s.Correct{iFile} =  [(t.%s.ChoicePress_LeftCorrect-t.%s.DelayEnd_LeftCorrect)/1e6;(t.%s.ChoicePress_RightCorrect-t.%s.DelayEnd_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
            eval(sprintf('Rlatency.%s.Correct{iFile} =  [(t.%s.ChoicePress_LeftCorrect-t.%s.NosePoke_LeftCorrect)/1e6;(t.%s.ChoicePress_RightCorrect-t.%s.NosePoke_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
            
            % Errors trials : This bugs out because some of the error trials are omissions,throwing out the alingmet
            try
                eval(sprintf('Slatency.%s.Error{iFile} =  [(t.%s.SamplePress_LeftError-t.%s.CueLight_LeftError)/1e6,(t.%s.SamplePress_RightError-t.%s.CueLight_RightError)/1e6]'';;',Delay_,Delay_,Delay_,Delay_,Delay_))
            end
            try
                eval(sprintf('NPlatency.%s.Error{iFile} =  [(t.%s.NosePoke_LeftError-t.%s.DelayEnd_LeftError)/1e6,(t.%s.NosePoke_RightError-t.%s.DelayEnd_RightError)/1e6]'';;',Delay_,Delay_,Delay_,Delay_,Delay_))
            end
            try
                eval(sprintf('Rlatency.%s.Error{iFile} =  [(t.%s.ChoicePress_LeftError-t.%s.NosePoke_LeftError)/1e6,(t.%s.ChoicePress_RightError-t.%s.NosePoke_RightError)/1e6]'';;',Delay_,Delay_,Delay_,Delay_,Delay_))
            end
        end
    end
    
    save(fullfile(pat,'allTimestamps','Behaviour analysis','Pharmacology',sprintf('Behav_%s.mat',Target)),'Slatency','NPlatency','Rlatency');
    clear fname NPlatency inputData t Rlatency Slatency iDelay iFile Target 
end

%% batch re-import
for iTarget = 1:length(targets)
    Target = targets{iTarget};
    fn = fullfile(pat,'allTimestamps','Behaviour analysis','Pharmacology',sprintf('Behav_%s.mat',Target));
    eval(sprintf('%s = load(fn,''Slatency'',''NPlatency'',''Rlatency'');',Target))
end
%% fractions correct/error/omitted
for iTarget = 1:length(targets)
    Target = targets{iTarget};
    fileList=dir(sprintf('allTimestamps%s*%s*.mat',filesep,Target));
    trialCounts{iTarget} = nan(length(fileList),3);
    for iFile =1:length(fileList)
        fname=strtok(fileList(iFile).name,'_');
        disp(sprintf('Analysing %s: %s...',Target,fname))
        fn = fullfile(pat,'allTimestamps',sprintf('%s_Events.mat',fname));
        load(fn,'inputData');

     
        trialCounts{iTarget}(iFile,:) = parseNLX_TTLsForRatios(inputData);
    end
end

%%
figure('color','w');
for iTarget = 1:length(targets)
    subplot(1,length(targets),iTarget); hold on
    x = sortrows(trialCounts{iTarget},1,'descend');
    b = bar(x,'stacked'); 
    b(1).FaceColor = 'b';
    b(2).FaceColor = 'r';
    b(3).FaceColor = [0.8 0.6 0.1];
    b(1).LineStyle	 = 'none';
    b(2).LineStyle	 = 'none';
    b(3).LineStyle	 = 'none';
    set(gca,'Xtick',[])
    if iTarget==1
        ylabel('Fraction of trials')
    else
        set(gca,'Ytick',[])
    end
    title(targets_{iTarget})
    xlabel(sprintf('(%d animals)',size(x,1)))
    axis([0 14 0 1])
end
legend({'Correct','Error','Omission'}); legend boxoff
%% Fraction of omisions
clear m e
for iTarget = 1:length(targets) 
    m(iTarget) = nanmean(trialCounts{iTarget}(:,3) ./ nanmean(trialCounts{1}(:,3)));
    e(iTarget) = nansem(trialCounts{iTarget}(:,3) ./ nanmean(trialCounts{1}(:,3)));        
end
% e=e./m(1);m=m./m(1);
figure; hold on
bar(m,'EdgeColor',[0.3 0.3 0.3],'FaceColor',[0.9 0.9 0.9])
errorbar(1:length(targets),m,e,'LineStyle','none','LineWidth',1.5,'Color',[0.3 0.3 0.3])
set(gca,'XTick',1:length(targets),'XTickLabel',targets2,'XTickLabelRotation',45)
plot([0,length(targets)+1],[1 1],':k')
ylabel('Norm. Fraction of omissions')

%% percentage correct respones 
ids = 1:length(targets) ;
% plot by delay duration
Delays = {'Short','Medium','Long'};
Delays_ = {'4s','8s','16s'};
prc_ = []; clear m e
for iTarget = ids%1:length(targets)
    Target = targets{iTarget};

    for iDelay = 1:3
        Delay_ = Delays{iDelay};
        eval(sprintf('C = %s.Slatency.%s.Correct;',Target,Delay_))
        eval(sprintf('E = %s.Slatency.%s.Error;',Target,Delay_))
        prc_{iDelay,iTarget} = cellfun(@numel,C)./(cellfun(@numel,C)+cellfun(@numel,E))*100;
    
        m(iDelay,iTarget) = nanmean(prc_{iDelay,iTarget});
        e(iDelay,iTarget) = nansem(prc_{iDelay,iTarget});
    end
end
% plot by delay duration
figure; 
for iTarget = ids%1:length(targets)
    subplot(1,length(targets),iTarget); hold on
    bar(m(:,iTarget),'EdgeColor',[0.3 0.3 0.3],'FaceColor',[0.9 0.9 0.9])
    errorbar(1:3,m(:,iTarget),e(:,iTarget),'LineStyle','none','LineWidth',1.5,'Color',[0.3 0.3 0.3])
    set(gca,'XTick',1:3,'XTickLabel',Delays_,'XTickLabelRotation',45)
    
     D = {prc_{:,iTarget}};
     [p,tbl,stats] =anova1(cell2mat(D)',char([ones(numel(D{1}),1);2*ones(numel(D{2}),1);3*ones(numel(D{3}),1)]),'off');
%      [p,tbl,stats] =kruskalwallis(cell2mat(D)',char([ones(numel(D{1}),1);2*ones(numel(D{2}),1);3*ones(numel(D{3}),1)]),'off');
     c = multcompare(stats,'Display','off');
     if c(1,6)<0.05 
         plot(2,m(2,iTarget)+e(2,iTarget)+5,'*r')
     end
     if c(2,6)<0.05 
         plot(3,m(3,iTarget)+e(3,iTarget)+5,'*r')
     end

    axis([0 4 0 100])
    if iTarget==1
        ylabel('% Correct responses')
    else
    end
    title(targets_{iTarget})

end

%% plot overlaid
figure; hold on
x = 1:3
subplot(1,2,1); hold on
errorbar(x,m(:,1),e(:,1),'k','lineWidth',3)     % control
errorbar(x+0.025,m(:,2),e(:,2),'r:','lineWidth',3) % AM
errorbar(x-0.025,m(:,3),e(:,3),'r','lineWidth',3)  % CP
errorbar(x,m(:,4),e(:,4),'-.r','lineWidth',3)  % URB
% errorbar(x+0.05,m(:,7),e(:,7),'-.r','lineWidth',3)  % HU
legend(targets2([1:4]),'Location','southeast');legend('boxoff')
plot([0.5 3.5],[50 50],':k','HandleVisibility','off')
ylabel('% Choices correct')
set(gca,'XTick',x,'XTickLabel',{'4s','8s','16s'})
xlabel('Delay')
axis([0 4 0 100])

subplot(1,2,2); hold on
errorbar(x,m(:,1),e(:,1),'k','lineWidth',3)     % contro
errorbar(x+0.025,m(:,5),e(:,5),'b:','lineWidth',3)   % Physo
errorbar(x-0.025,m(:,6),e(:,6),'b','lineWidth',3)  % Scop
legend(targets2([1,5,6]),'Location','southeast');legend('boxoff')
plot([0.5 3.5],[50 50],':k','HandleVisibility','off')
set(gca,'XTick',x,'XTickLabel',{'4s','8s','16s'})
xlabel('Delay')
axis([0 4 0 100])
%%
for iTarget = ids(2:end)
    errorbar(x,m(:,iTarget),e(:,iTarget))
end

% plot by drug condition
figure; 
for iDelay = 1:3
        subplot(1,3,iDelay); hold on
        bar(m(iDelay,:),'EdgeColor',[0.3 0.3 0.3],'FaceColor',[0.9 0.9 0.9])
        errorbar(ids,m(iDelay,ids),e(iDelay,ids),'LineStyle','none','LineWidth',1.5,'Color',[0.3 0.3 0.3])
        set(gca,'XTick',ids,'XTickLabel',targets2(ids),'XTickLabelRotation',45)
        
        D = {prc_{iDelay,ids}};
        names = [];
        for iTarget = ids
            if iTarget ==1
                names = [names;repmat({'Control'},length(D{iTarget}),1)];
            else
                names = [names;repmat({targets{iTarget}},length(D{iTarget}),1)];
            end                
        end
%         Fisherextest
        
        [p,tbl,stats] = anova1(cell2mat(D)',names,'off');
%         [p,tbl,stats] = kruskalwallis(cell2mat(D)',names,'off');
        c = multcompare(stats,'Display','off');
        for iTarget = ids
            if c(iTarget,6)<0.05
                plot(iTarget+1,m(iDelay,iTarget+1)+e(iDelay,iTarget+1)+5,'*r')
            end
        end
        
        axis([0 10 0 120])
        if iDelay==1
            ylabel('% Correct responses')
        end
        title(Delays_{iDelay})

    
end
%% percentage correct respones (normalised to vehicle)
% plot by delay duration
Delays = {'Short','Medium','Long'};
Delays_ = {'4s','8s','16s'};
prc_ = [];
for iTarget = 1:length(targets)
    Target = targets{iTarget};

    for iDelay = 1:3
        Delay_ = Delays{iDelay};
        eval(sprintf('C = %s.Slatency.%s.Correct;',Target,Delay_))
        eval(sprintf('E = %s.Slatency.%s.Error;',Target,Delay_))
        prc_{iDelay,iTarget} = cellfun(@numel,C)./(cellfun(@numel,C)+cellfun(@numel,E));
%         if iTarget>1
%             prc_{iDelay,iTarget} = prc_{iDelay,iTarget} ./nanmean(prc_{iDelay,1}); % per delay
%         end
        prc_{iDelay,iTarget} = prc_{iDelay,iTarget} ./nanmean(prc_{1,1}); % Norm. to shortest delay
        m(iDelay,iTarget) = nanmean(prc_{iDelay,iTarget});
        e(iDelay,iTarget) = nansem(prc_{iDelay,iTarget});
    end
end
iTarget = 1;
for iDelay = 1:3
    
   
        
        % prc_{iDelay,iTarget} = prc_{iDelay,iTarget} ./nanmean(prc_{iDelay,1});% per delay
        prc_{iDelay,iTarget} = prc_{iDelay,iTarget} ./nanmean(prc_{1,1});% Norm. to shortest delay
        m(iDelay,iTarget) = nanmean(prc_{iDelay,iTarget});
        e(iDelay,iTarget) = nansem(prc_{iDelay,iTarget});
        
end


% plot by delay duration
figure; 
for iTarget = 1:length(targets)
    subplot(1,length(targets),iTarget); hold on
    bar(m(:,iTarget),'EdgeColor',[0.3 0.3 0.3],'FaceColor',[0.9 0.9 0.9])
    errorbar(1:3,m(:,iTarget),e(:,iTarget),'LineStyle','none','LineWidth',1.5,'Color',[0.3 0.3 0.3])
    set(gca,'XTick',1:3,'XTickLabel',Delays_,'XTickLabelRotation',45)
    
     D = {prc_{:,iTarget}};
     [p,tbl,stats] =anova1(cell2mat(D)',char([ones(numel(D{1}),1);2*ones(numel(D{2}),1);3*ones(numel(D{3}),1)]),'off');
%      [p,tbl,stats] =kruskalwallis(cell2mat(D)',char([ones(numel(D{1}),1);2*ones(numel(D{2}),1);3*ones(numel(D{3}),1)]),'off');
     c = multcompare(stats,'Display','off');
     if c(1,6)<0.05 
         plot(2,m(2,iTarget)+e(2,iTarget)+5,'*r')
     end
     if c(2,6)<0.05 
         plot(3,m(3,iTarget)+e(3,iTarget)+5,'*r')
     end

    axis([0 4 0 2])
    if iTarget==1
        ylabel('% Correct responses')
    else
    end
    title(targets_{iTarget})

end

% plot by drug condition
ids = 1:length(targets) ;
figure; 
for iDelay = 1:3
        subplot(1,3,iDelay); hold on
        plot([0.5 max(ids)+0.5],[1 1],':k' )
        bar(m(iDelay,ids),'EdgeColor',[0.3 0.3 0.3],'FaceColor',[0.9 0.9 0.9])
        errorbar(ids,m(iDelay,ids),e(iDelay,ids),'LineStyle','none','LineWidth',1.5,'Color',[0.3 0.3 0.3])
        set(gca,'XTick',ids,'XTickLabel',targets_(ids),'XTickLabelRotation',45)
        
        D = {prc_{iDelay,ids}};
        names = [];
        for iTarget = ids
            if iTarget ==1
                names = [names;repmat({'Control'},length(D{iTarget}),1)];
            else
                names = [names;repmat({targets{iTarget}},length(D{iTarget}),1)];
            end                
        end
        [p,tbl,stats] = anova1(cell2mat(D(ids))',names,'off');
%         [p,tbl,stats] = kruskalwallis(cell2mat(D(ids))',names,'off');
        c = multcompare(stats,'Display','off');
        for iTarget = ids%1:length(targets)
        if c(iTarget,6)<0.05
            plot(iTarget+1,1.8,'*r')
        end
        end
         
%         axis([0 10 0 2])
        if iDelay==1
            ylabel('% Correct responses')
        end
        title(Delays_{iDelay})
        
        

    
end

% plot overlaid

figure; hold on

x = 1:3
errorbar(x,m(:,1),e(:,1),'k','lineWidth',1.5)
for iTarget = ids(2:end)
    errorbar(x,m(:,iTarget),e(:,iTarget))
end
legend(targets2);legend('boxoff')

%% Plot Sample latency distribution
figure('name','Sample latency');
Delays = {'Short','Medium','Long'};
Delays_ = {'4s delay','8s delay','16s delay'};
for iDelay = 1:3
    Delay_ = Delays{iDelay};
    subplot(1,3,iDelay); hold on
    D=[];
    names=[];
    for iTarget = 1:length(targets)
        D_=[];
        Target = targets{iTarget};
        eval(sprintf('D_ = cell2mat(%s.Slatency.%s.Correct'')',Target,Delay_))
        D = [D;D_];
        if strcmp(Target,'LONG')
            names = [names;repmat({'Control'},length(D_),1)];
        else
            names = [names;repmat({Target},length(D_),1)];
        end
    end
    
    boxplot(D,names,'Orientation','vertical','whisker',1000,'Colors','k')
    set(gca,'XTickLabelRotation',45)
    title(Delays_{iDelay})
    ylim([0 22])
    if iDelay==1
       ylabel('Sample latency (s)') 
    end
    clear p tbl stats c
        [p{iDelay},tbl{iDelay},stats{iDelay}] = kruskalwallis(D,names,'off');
%     [p{iDelay},tbl{iDelay},stats{iDelay}] = anova1(D,names,'off');
    c{iDelay} = multcompare(stats{iDelay},[],'off');
    
    for iTarget = 2:length(targets)
        if c{iDelay}(iTarget-1,6)<0.05
            scatter(iTarget,21,'*r')
        end
    end
end       
%% Plot Nose-poke latency distribution
figure('name','Nosepoke latency');
Delays = {'Short','Medium','Long'};
Delays_ = {'4s delay','8s delay','16s delay'};
for iDelay = 1:3
    Delay_ = Delays{iDelay};
    subplot(1,3,iDelay); hold on
    D=[];
    names=[];
    for iTarget = 1:length(targets)
        D_=[];
        Target = targets{iTarget};
        eval(sprintf('D_ = cell2mat(%s.NPlatency.%s.Correct'')',Target,Delay_))
        D = [D;D_];
        if strcmp(Target,'LONG')
            names = [names;repmat({'Control'},length(D_),1)];
        else
            names = [names;repmat({Target},length(D_),1)];
        end
    end
    
    boxplot(D,names,'Orientation','vertical','whisker',1000,'Colors','k')
    set(gca,'XTickLabelRotation',45)
    title(Delays_{iDelay})
    ylim([0 25])
    if iDelay==1
       ylabel('Nose-poke latency (s)') 
    end
    clear p tbl stats c
        [p{iDelay},tbl{iDelay},stats{iDelay}] = kruskalwallis(D,names,'off');
%     [p{iDelay},tbl{iDelay},stats{iDelay}] = anova1(D,names,'off');
    c{iDelay} = multcompare(stats{iDelay},[],'off');
    
    for iTarget = 2:length(targets)
        if c{iDelay}(iTarget-1,6)<0.05
            scatter(iTarget,22,'*r')
        end
    end
end
%% Plot Choice latency distribution
figure('name','Choice latency');
Delays = {'Short','Medium','Long'};
Delays_ = {'4s delay','8s delay','16s delay'};
for iDelay = 1:3
    Delay_ = Delays{iDelay};
    subplot(1,3,iDelay); hold on
    D=[];
    names=[];
    for iTarget = 1:length(targets)
        D_=[];
        Target = targets{iTarget};
        eval(sprintf('D_ = cell2mat(%s.Rlatency.%s.Correct'')',Target,Delay_))
        D = [D;D_];
        if strcmp(Target,'LONG')
            names = [names;repmat({'Control'},length(D_),1)];
        else
            names = [names;repmat({Target},length(D_),1)];
        end
    end
    
    boxplot(D,names,'Orientation','vertical','whisker',1000,'Colors','k')
    set(gca,'XTickLabelRotation',45)
    title(Delays_{iDelay})
    ylim([0 12])
    if iDelay==1
       ylabel('Choice latency (s)') 
    end
    clear p tbl stats c
        [p{iDelay},tbl{iDelay},stats{iDelay}] = kruskalwallis(D,names,'off');
%     [p{iDelay},tbl{iDelay},stats{iDelay}] = anova1(D,names,'off');
    c{iDelay} = multcompare(stats{iDelay},[],'off');
    
    for iTarget = 2:length(targets)
        if c{iDelay}(iTarget-1,6)<0.05
            scatter(iTarget,11,'*r')
        end
    end
end
%% Plot Sample latency distribution (collapsed across delays)
figure('name','Sample latency'); hold on
Delays = {'Short','Medium','Long'};
Delays_ = {'4s delay','8s delay','16s delay'};
    D=[];
    D__ = cell(1,length(targets));
    names=[];
    for iTarget = 1:length(targets)
        D_=[];
        Target = targets{iTarget};
        for iDelay = 1:3
            Delay_ = Delays{iDelay};
            eval(sprintf('D_ = cell2mat(%s.Slatency.%s.Correct'');',Target,Delay_))
            D = [D;D_];
            
            if strcmp(Target,'LONG')
                names = [names;repmat({'Control'},length(D_),1)];
            else
                names = [names;repmat({Target},length(D_),1)];
            end
            D__{iTarget}=D;
        end
    end
    
    boxplot(D,names,'Orientation','vertical','whisker',1000,'Colors','k')
    set(gca,'XTickLabelRotation',45)
    title(Delays_{iDelay})
    ylim([0 22])
    if iDelay==1
       ylabel('Sample latency (s)') 
    end
    clear p tbl stats c
%         [p,tbl,stats] = kruskalwallis(D,names,'off');
    [p,tbl,stats] = anova1(D,names,'off');
    c = multcompare(stats,[],'off');
    
    for iTarget = 2:length(targets)
        if c(iTarget-1,6)<0.05
            scatter(iTarget,21,'*r')
        end
    end
    
    
    
%     Nmax = max(cellfun(@length,D__));
%     Y = nan(Nmax,length(targets));
%     for iTarget = 1:length(targets)
%         Y(1:length(D__{iTarget}),iTarget) = D__{iTarget};
%     end
%     figure; hold on
%     violin (Y),'bw',0.05,'xlabel',targets_)
    
    figure; hold on
    
    for iTarget = 1:length(targets)
        subplot(1,length(targets),iTarget);
        violin(D__{iTarget},'bw',0.05);
        axis([ 0.5 1.5 0 20])
        legend off
        axis off
    end

bins = 0:0.05:20;
D__ = zeros(length(bins),length(targets));
for iTarget = 1:length(targets)
    D=[];
    Target = targets{iTarget};
    temp = [eval(sprintf('%s.Slatency.Short.Correct',Target));...
        eval(sprintf('%s.Slatency.Medium.Correct',Target));...
        eval(sprintf('%s.Slatency.Long.Correct',Target))];
    for iFile =1:size(temp,2)
        Y = cell2mat(temp(1:3,iFile));
        Y = histc(Y,bins); Y = Y./nansum(Y);
        D(:,iFile) = Y;
    end
    D__(:,iTarget) = nanmean(D,2);
end
    



figure; hold on
% cmap_ = jet(length(targets));
cmap_ = hsv(length(targets));cmap_(1,:)=0;
for iTarget = 1:length(targets)
    Y = D__(:,iTarget);
    if iTarget ==1
        stairs(bins,cumsum(Y),'LineWidth',5,'color',cmap_(iTarget,:))
    else
        stairs(bins,cumsum(Y),'LineWidth',3,'color',cmap_(iTarget,:))
    end
end
legend(targets2,'location','southeast');legend('boxoff')
axis([0 max(bins) 0 1])

Y = D__(:,1);
        stairs(bins,cumsum(Y),'LineWidth',5,'color',cmap_(1,:))
%% Plot Nose-poke latency distribution (collapsed across delays)
figure('name','Nose-poke latency'); hold on
Delays = {'Short','Medium','Long'};
Delays_ = {'4s delay','8s delay','16s delay'};
    D=[];
    D__ = cell(1,length(targets));
    names=[];
    for iTarget = 1:length(targets)
        D_=[];
        Target = targets{iTarget};
        for iDelay = 1:3
            Delay_ = Delays{iDelay};
            eval(sprintf('D_ = cell2mat(%s.NPlatency.%s.Correct'');',Target,Delay_))
            D = [D;D_];
            
            if strcmp(Target,'LONG')
                names = [names;repmat({'Control'},length(D_),1)];
            else
                names = [names;repmat({Target},length(D_),1)];
            end
            D__{iTarget}=D;
        end
    end
    
    boxplot(D,names,'Orientation','vertical','whisker',1000,'Colors','k')
    set(gca,'XTickLabelRotation',45)
    title(Delays_{iDelay})
    ylim([0 22])
    if iDelay==1
       ylabel('Sample latency (s)') 
    end
    clear p tbl stats c
%         [p,tbl,stats] = kruskalwallis(D,names,'off');
    [p,tbl,stats] = anova1(D,names,'off');
    c = multcompare(stats,[],'off');
    
    for iTarget = 2:length(targets)
        if c(iTarget-1,6)<0.05
            scatter(iTarget,21,'*r')
        end
    end
    
    
    
%     Nmax = max(cellfun(@length,D__));
%     Y = nan(Nmax,length(targets));
%     for iTarget = 1:length(targets)
%         Y(1:length(D__{iTarget}),iTarget) = D__{iTarget};
%     end
%     figure; hold on
%     violin (Y),'bw',0.05,'xlabel',targets_)
    
    figure; hold on
    
    for iTarget = 1:length(targets)
        subplot(1,length(targets),iTarget);
        violin(D__{iTarget},'bw',0.05);
        axis([ 0.5 1.5 0 20])
        legend off
        axis off
    end

bins = 0:0.05:20;
D__ = zeros(length(bins),length(targets));
for iTarget = 1:length(targets)
    D=[];
    Target = targets{iTarget};
    temp = [eval(sprintf('%s.NPlatency.Short.Correct',Target));...
        eval(sprintf('%s.NPlatency.Medium.Correct',Target));...
        eval(sprintf('%s.NPlatency.Long.Correct',Target))];
    for iFile =1:size(temp,2)
        Y = cell2mat(temp(1:3,iFile));
        Y = histc(Y,bins); Y = Y./nansum(Y);
        D(:,iFile) = Y;
    end
    D__(:,iTarget) = nanmean(D,2);
end
    



figure; hold on
% cmap_ = jet(length(targets));
cmap_ = hsv(length(targets));cmap_(1,:)=0;
for iTarget = 1:length(targets)
    Y = D__(:,iTarget);
    if iTarget ==1
        stairs(bins,cumsum(Y),'LineWidth',5,'color',cmap_(iTarget,:))
    else
        stairs(bins,cumsum(Y),'LineWidth',3,'color',cmap_(iTarget,:))
    end
end
legend(targets2,'location','southeast');legend('boxoff')
axis([0 max(bins) 0 1])

Y = D__(:,1);
        stairs(bins,cumsum(Y),'LineWidth',5,'color',cmap_(1,:))
%% Plot Choice latency distribution (collapsed across delays)
figure('name','Choice latency'); hold on
Delays = {'Short','Medium','Long'};
Delays_ = {'4s delay','8s delay','16s delay'};
    D=[];
    D__ = cell(1,length(targets));
    names=[];
    for iTarget = 1:length(targets)
        D_=[];
        Target = targets{iTarget};
        for iDelay = 1:3
            Delay_ = Delays{iDelay};
            eval(sprintf('D_ = cell2mat(%s.Rlatency.%s.Correct'');',Target,Delay_))
            D = [D;D_];
                
            if strcmp(Target,'LONG')
                names = [names;repmat({'Control'},length(D_),1)];
            else
                names = [names;repmat({Target},length(D_),1)];
            end
            D__{iTarget}=D;
        end
    end
    
    boxplot(D,names,'Orientation','vertical','whisker',1000,'Colors','k')
    set(gca,'XTickLabelRotation',45)
    title(Delays_{iDelay})
    ylim([0 12])
    if iDelay==1
       ylabel('Sample latency (s)') 
    end
    clear p tbl stats c
        [p,tbl,stats] = kruskalwallis(D,names,'off');
%     [p,tbl,stats] = anova1(D,names,'off');
    c = multcompare(stats,[],'off');
    scatter(1:length(targets),cellfun(@nanmean,D__),'*r')
    for iTarget = 2:length(targets)
        if c(iTarget-1,6)<(0.05/(length(targets)-1))
            scatter(iTarget,11,'*r')
        end
    end
    
    
    
%     Nmax = max(cellfun(@length,D__));
%     Y = nan(Nmax,length(targets));
%     for iTarget = 1:length(targets)
%         Y(1:length(D__{iTarget}),iTarget) = D__{iTarget};
%     end
%     figure; hold on
%     violin (Y),'bw',0.05,'xlabel',targets_)
    
    figure; hold on
    
    for iTarget = 1:length(targets)
        subplot(1,length(targets),iTarget);
        violin(D__{iTarget},'bw',0.05);
        axis([ 0.5 1.5 0 20])
        legend off
        axis off
    end

bins = 0:0.05:10;
D__ = zeros(length(bins),length(targets));
for iTarget = 1:length(targets)
    D=[];
    Target = targets{iTarget};
    temp = [eval(sprintf('%s.Rlatency.Short.Correct',Target));...
        eval(sprintf('%s.Rlatency.Medium.Correct',Target));...
        eval(sprintf('%s.Rlatency.Long.Correct',Target))];
    for iFile =1:size(temp,2)
        Y = cell2mat(temp(1:3,iFile));
        Y = histc(Y,bins); Y = Y./nansum(Y);
        D(:,iFile) = Y;
    end
    D__(:,iTarget) = nanmean(D,2);
end
    



figure; hold on
% cmap_ = jet(length(targets));
cmap_ = hsv(length(targets));cmap_(1,:)=0;
for iTarget = 1:length(targets)
    Y = D__(:,iTarget);
    if iTarget ==1
        stairs(bins,cumsum(Y),'LineWidth',5,'color',cmap_(iTarget,:))
    else
        stairs(bins,cumsum(Y),'LineWidth',3,'color',cmap_(iTarget,:))
    end
end
legend(targets2,'location','southeast');legend('boxoff')
axis([0 max(bins) 0 1])

Y = D__(:,1);
        stairs(bins,cumsum(Y),'LineWidth',5,'color',cmap_(1,:))
 
%% Plot Sample latency distribution - errors
figure('name','Sample latency - errors');
Delays = {'Short','Medium','Long'};
Delays_ = {'4s delay','8s delay','16s delay'};
for iDelay = 1:3
    Delay_ = Delays{iDelay};
    subplot(1,3,iDelay); hold on
    D=[];
    names=[];
    for iTarget = 1:length(targets)
        D_=[];
        Target = targets{iTarget};
        eval(sprintf('D_ = cell2mat(%s.Slatency.%s.Error'')',Target,Delay_))
        D = [D;D_];
        if strcmp(Target,'LONG')
            names = [names;repmat({'Control'},length(D_),1)];
        else
            names = [names;repmat({Target},length(D_),1)];
        end
    end
    
    boxplot(D,names,'Orientation','vertical','whisker',1000,'Colors','k')
    set(gca,'XTickLabelRotation',45)
    title(Delays_{iDelay})
    ylim([0 20])
    if iDelay==1
       ylabel('Sample latency (s)') 
    end
    
% %     [p{iDelay},tbl{iDelay},stats{iDelay}] = kruskalwallis(D,names,'off');
%     [p{iDelay},tbl{iDelay},stats{iDelay}] = anova1(D,names,'off');
%     c{iDelay} = multcompare(stats{iDelay});

end
%% Plot Nose-poke latency distribution - errors
figure('name','Nosepoke latency- errors');
Delays = {'Short','Medium','Long'};
Delays_ = {'4s delay','8s delay','16s delay'};
for iDelay = 1:3
    Delay_ = Delays{iDelay};
    subplot(1,3,iDelay); hold on
    D=[];
    names=[];
    for iTarget = 1:length(targets)
        D_=[];
        Target = targets{iTarget};
        eval(sprintf('D_ = cell2mat(%s.NPlatency.%s.Error'')',Target,Delay_))
        D = [D;D_];
        if strcmp(Target,'LONG')
            names = [names;repmat({'Control'},length(D_),1)];
        else
            names = [names;repmat({Target},length(D_),1)];
        end
    end
    
    boxplot(D,names,'Orientation','vertical','whisker',1000,'Colors','k')
    set(gca,'XTickLabelRotation',45)
    title(Delays_{iDelay})
    ylim([0 25])
    if iDelay==1
       ylabel('Nose-poke latency (s)') 
    end
% %     [p{iDelay},tbl{iDelay},stats{iDelay}] = kruskalwallis(D,names,'off');
%     [p{iDelay},tbl{iDelay},stats{iDelay}] = anova1(D,names,'off');
%     c{iDelay} = multcompare(stats{iDelay});
end
%% Plot Choice latency distribution - errors
figure('name','Choice latency - errors');
Delays = {'Short','Medium','Long'};
Delays_ = {'4s delay','8s delay','16s delay'};
for iDelay = 1:3
    Delay_ = Delays{iDelay};
    subplot(1,3,iDelay); hold on
    D=[];
    names=[];
    for iTarget = 1:length(targets)
        D_=[];
        Target = targets{iTarget};
        eval(sprintf('D_ = cell2mat(%s.Rlatency.%s.Error'')',Target,Delay_))
        D = [D;D_];
        if strcmp(Target,'LONG')
            names = [names;repmat({'Control'},length(D_),1)];
        else
            names = [names;repmat({Target},length(D_),1)];
        end
    end
    
    boxplot(D,names,'Orientation','vertical','whisker',1000,'Colors','k')
    set(gca,'XTickLabelRotation',45)
    title(Delays_{iDelay})
    ylim([0 10])
    if iDelay==1
       ylabel('Choice latency (s)') 
    end
% %     [p{iDelay},tbl{iDelay},stats{iDelay}] = kruskalwallis(D,names,'off');
%     [p{iDelay},tbl{iDelay},stats{iDelay}] = anova1(D,names,'off');
%     c{iDelay} = multcompare(stats{iDelay});
end

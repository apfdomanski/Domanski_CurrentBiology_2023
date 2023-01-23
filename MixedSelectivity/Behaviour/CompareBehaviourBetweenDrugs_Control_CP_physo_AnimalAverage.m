%% %%%%%% PREAMBLE %%%%%%
% targets = {'LONG','AM251','CP55940','URB597','Physostigmine','Scopolamine','HU210','CPphyso','CPAM'};
% targets2 = {'Vehicle','AM251','CP55940','URB597','Physostigmine','Scopolamine','HU210','CP55940 & Physostigmine','CP55940 & AM251'}; 
% targets_ = {'Control','AM251','CP55940','URB597','Physostigmine','Scopolamine','HU210',{'CP55940 &';'Physostigmine'},{'CP55940';'& AM251'}};

% targets = {'LONG','AM251','CP55940','Scopolamine','Physostigmine','URB597'};
% targets2 = {'Vehicle','AM251','CP55940','Scopolamine','Physostigmine','URB597'}; 
% targets_ = {'Control','AM251','CP55940','Scopolamine','Physostigmine','URB597'};


targets = {'LONG','CP55940','Physostigmine'};
targets2 = {'Vehicle','CP55940','Physostigmine'}; 
targets_ = {'Control','CP55940','Physostigmine'};

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
%         pat = '/Volumes/HDD2/DNMTP/raw/';
        pat = '/Volumes/Akasa/DNMTP/raw/'
        cd(pat)
        fileList=dir(sprintf('allTimestamps/*%s*.mat',Target));
    end
    
    
    
    reject_list={'MiroslawLONG2_Events.mat'}; % These only have 5s delays
    name_flag=zeros(numel(fileList),1);
    for idx=1:numel(fileList)
        fnames_{idx,1}=fileList(idx).name;
        name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
    end
    try
        fileList(find(name_flag))=[];% fnames_(name_flag)=[];
    end
    clear  reject_list idx fnames_ name_flag
    
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
%Control,CP:  jar Krzes Kryz Miro Norb Onu
%Physo:  jar Krzes Kryz Miro Norb

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

figure('color','w');
for iTarget = 1:length(targets)
    subplot(1,length(targets),iTarget); hold on
    tc = trialCounts{iTarget};
    if iTarget ==3
        tc = [tc(1:4,:);nan(1,3);tc(5,:)];
    end
    if iTarget ==1
        [x,idx] = sortrows(tc,1,'descend');
    else
        x = tc(idx,:);
    end
    b = bar(x,'stacked'); 
    b(1).FaceColor = 'k';
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
%     title(targets_{iTarget})
%     xlabel(sprintf('(%d sessions)',size(x,1)))
    axis([0 14 0 1])
end
plots=get(gca, 'Children');
legend(plots([1,2,3]),fliplr({'Correct','Error','Omission'}));legend boxoff

% legend({'Correct','Error','Omission'}); legend boxoff
%% Fraction of omisions
clear m e
tempOut = [];
for iTarget = 1:length(targets) 
    temp = trialCounts{iTarget}(:,3) ./ nanmean(trialCounts{1}(:,3));
    m(iTarget) = nanmean(temp);
    e(iTarget) = nansem(temp);        
    
    tempOut = [tempOut;[temp,iTarget*ones(size(temp))]];
end
% e=e./m(1);m=m./m(1);
figure; hold on
bar(m,'EdgeColor',[0.3 0.3 0.3],'FaceColor',[0.9 0.9 0.9])
errorbar(1:length(targets),m,e,'LineStyle','none','LineWidth',1.5,'Color',[0.3 0.3 0.3])
set(gca,'XTick',1:length(targets),'XTickLabel',targets2,'XTickLabelRotation',45)
plot([0,length(targets)+1],[1 1],':k')
ylabel('Norm. Fraction of omissions')

anova1(tempOut(:,1),tempOut(:,2))

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

%% plot overlaid - CP and Physo only
figure; hold on
x = 1:3
hold on
errorbar(x,m(:,1),e(:,1),'k','lineWidth',3)                     % control
errorbar(x-0.025,m(:,2),e(:,2),'r','lineWidth',3)               % CP
errorbar(x+0.025,m(:,3),e(:,3),'b','lineWidth',3)               % Physo
% errorbar(x+0.025,m(:,8),e(:,8),'g','lineWidth',3)             % CP+ Physo
legend(targets2,'Location','southeast');legend('boxoff')
plot([0.5 3.5],[50 50],':k','HandleVisibility','off')
ylabel('% Choices correct')
set(gca,'XTick',x,'XTickLabel',{'4s','8s','16s'})
xlabel('Delay')
axis([0 4 30 100])

xlabel('Delay')
axis([0 4 30 100])
legend(targets2(),'Location','southwest');legend('boxoff')
%% Test performance - CP and Physo only
ids = [1:3];
% prc_{iDelay,iTarget} 
% AnIDs = {[7,1,1,2,2,3,3,4,4,4,5,5,6,6],[1,2,3,4,5,6],[7,1,2,3,4,6]}
figure; hold on
for iTarget = ids(1:end)
    errorbar(x,m(:,iTarget),e(:,iTarget))
end

% plot by drug condition
% figure; 
for iDelay = 1:3
%         subplot(1,3,iDelay); hold on
%         bar(m(iDelay,:),'EdgeColor',[0.3 0.3 0.3],'FaceColor',[0.9 0.9 0.9])
%         errorbar(ids,m(iDelay,ids),e(iDelay,ids),'LineStyle','none','LineWidth',1.5,'Color',[0.3 0.3 0.3])
%         set(gca,'XTick',ids,'XTickLabel',targets2(ids),'XTickLabelRotation',45)
        
        D = {prc_{iDelay,ids}};
        names = [];
        for iTarget = ids
            if iTarget ==1
                names = [names;repmat({'Control'},length(D{iTarget}),1)];
            else
                names = [names;repmat({targets{iTarget}},length(D{iTarget}),1)];
            end                
        end
        
        
%         [p,tbl,stats] = anova1(cell2mat(D)',names,'off');
[p,tbl,stats] = kruskalwallis(cell2mat(D)',names,'off');
        p_(iDelay)=p;
        stats_{iDelay}=stats;
%         Fisherextest
        

        
       
        axis([0.5 3.5 0 120])
        if iDelay==1
            ylabel('% Correct responses')
        end
        title(Delays_{iDelay})
end

for iDelay=1:3
        
        figure
        c = multcompare(stats_{iDelay},'Display','on');
        c_{iDelay}=c;
%         for iTarget = ids
%             if c(iTarget,6)<0.05
%                 plot(iTarget,m(iDelay,iTarget)+e(iDelay,iTarget)+5,'*r')
%             end
%         end
        
end
% Analysing LONG1: JaroslawLONG1... 1
% Analysing LONG1: KrzesimirLONG1... 2
% Analysing LONG1: KrzysztofLONG1... 3
% Analysing LONG1: MiroslawLONG1... 4
% Analysing LONG1: NorbertLONG1... 5
% Analysing LONG1: OnufryLONG1... 6
% Analysing CP55940: JaroslawCP55940... 1
% Analysing CP55940: KrzesimirCP55940... 2
% Analysing CP55940: KrzysztofCP55940... 3
% Analysing CP55940: MiroslawCP55940... 4
% Analysing CP55940: NorbertCP55940... 5
% Analysing CP55940: OnufryCP55940... 6
% Analysing Physostigmine: Jaroslawphysostigmine... 1
% Analysing Physostigmine: Krzesimirphysostigmine... 2
% Analysing Physostigmine: Krzysztofphysostigmine... 3
% Analysing Physostigmine: Miroslawphysostigmine... 4
% Analysing Physostigmine: Onufryphysostigmine... 6

% AnIDs = {[7,1,1,2,2,3,3,4,4,4,5,5,6,6],[1,2,3,4,5,6],[7,1,2,3,4,6]}
AnIDs = {[1:6],[1:6],[1,2,3,4,6]};
ids= 1:3;
clear h p 
for iDelay = 1:3      
    D = {prc_{iDelay,ids}};
%     [h(iDelay,1),p(iDelay,1)]=ttest(D{1},D{2},'Alpha',0.05/2)
%     [h(iDelay,2),p(iDelay,2)]=ttest(D{1},[D{3}(1:4),NaN,D{3}(5)],'Alpha',0.05/2)
    [h(iDelay,1),p(iDelay,1)]=signrank(D{1},D{2},'Alpha',0.05/2)
    [h(iDelay,2),p(iDelay,2)]=signrank(D{1},[D{3}(1:4),NaN,D{3}(5)],'Alpha',0.05/2)
    
end

%% RM-ANOVA - CP and Physo only
AnIDs = {[1:6],[1:6],[1,2,3,4,6]};
cellfun(@length,prc_,'UniformOutput',false)
% prc_{iDelay,ids}
% control_ = [prc_{1,1}';prc_{2,1}';prc_{3,1}'];
% CP_      = [prc_{1,2}';prc_{2,2}';prc_{3,2}'];
% physo_   = [[prc_{1,3}(1:4),NaN,prc_{1,3}(5)]';[prc_{2,3}(1:4),NaN,prc_{2,3}(5)]';[prc_{3,3}(1:4),NaN,prc_{3,3}(5)]'];
control_4_1 = prc_{1,1}([1,3,5,7,9,11])';
control_4_2 = prc_{1,1}([2,4,6,8,10,12])';
control_8_1 = prc_{2,1}([1,3,5,7,9,11])';
control_8_2 = prc_{2,1}([2,4,6,8,10,12])';
control_16_1 = prc_{3,1}([1,3,5,7,9,11])';
control_16_2 = prc_{3,1}([2,4,6,8,10,12])';

CP_4  = prc_{1,2}';
CP_8  = prc_{2,2}';
CP_16 = prc_{3,2}';
physo_4 = [prc_{1,3}(1:4),NaN,prc_{1,3}(5)]';
physo_8 = [prc_{2,3}(1:4),NaN,prc_{2,3}(5)]';
physo_16 = [prc_{3,3}(1:4),NaN,prc_{3,3}(5)]';


performanceTbl_ = table(control_4_1,control_4_2,control_8_1,control_8_2,control_16_1,control_16_2,CP_4,CP_8,CP_16,physo_4,physo_8,physo_16) % create table
% within = table([1 1 1 2 2 2 3 3 3]',[1 2 3 1 2 3 1 2 3]','VariableNames',{'Drug' 'Delay'}) % within model
within = table([1 2 3 1 2 3 1 2 3 1 2 3]',[1 1 1 1 1 1 2 2 2 3 3 3]','VariableNames',{'Delay' 'Drug'}) % within model

performance_rm = fitrm(performanceTbl_,'control_4_1-physo_16~1','WithinDesign',within,'WithinModel','Drug*Delay')

[ranovatbl,A,C,D] = ranova(performance_rm,'WithinModel','Drug*Delay')
multcompare(performance_rm,'Drug','By','Delay')
% multcompare(performance_rm,'Delay','By','Drug')

% A = load("repeatedmeas.mat")

% 'JaroslawLONG1_Events.mat'; 1/1   1
% 'JaroslawLONG2_Events.mat'; 1/2       2
% 'KrzesimirLONG1_Events.mat';2/1   3
% 'KrzesimirLONG2_Events.mat';2/2       4
% 'KrzysztofLONG1_Events.mat';3/1   5
% 'KrzysztofLONG2_Events.mat';3/2       6
% 'MiroslawLONG0_Events.mat';4/1    7
% 'MiroslawLONG1_Events.mat';4/2        8
% 'NorbertLONG1_Events.mat';5/1     9
% 'NorbertLONG2_Events.mat';5/2         10
% 'OnufryLONG1_Events.mat';6/1      11
% 'OnufryLONG2_Events.mat';6/2          12

% Analysing CP55940: JaroslawCP55940... 1
% Analysing CP55940: KrzesimirCP55940... 2
% Analysing CP55940: KrzysztofCP55940... 3
% Analysing CP55940: MiroslawCP55940... 4
% Analysing CP55940: NorbertCP55940... 5
% Analysing CP55940: OnufryCP55940... 6
% Analysing Physostigmine: Jaroslawphysostigmine... 1
% Analysing Physostigmine: Krzesimirphysostigmine... 2
% Analysing Physostigmine: Krzysztofphysostigmine... 3
% Analysing Physostigmine: Miroslawphysostigmine... 4
% Analysing Physostigmine: Onufryphysostigmine... 6



%% percentage correct responses (normalised to vehicle - group)
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
%% percentage correct responses (normalised to vehicle - individual animals)
% plot by delay duration
Delays = {'Short','Medium','Long'};
Delays_ = {'4s','8s','16s'};
AnIDs = {[1:6],[1:6],[1,2,3,4,6]};
    
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
%         prc_{iDelay,iTarget} = prc_{iDelay,iTarget} ./nanmean(prc_{1,1}); % Norm. to shortest delay
        if iTarget ==3
            prc_{iDelay,iTarget}  = [prc_{iDelay,iTarget}(1:4),NaN,prc_{iDelay,iTarget}(5)];
        end
        m(iDelay,iTarget) = nanmean(prc_{iDelay,iTarget});
        e(iDelay,iTarget) = nansem(prc_{iDelay,iTarget});
    end
end
prcNorm_ = [];
for iTarget = 1:length(targets)
    for iDelay = 1:3
        prcNorm__{iDelay,iTarget} = prc_{iDelay,iTarget} ./prc_{iDelay,1};
        m(iDelay,iTarget) = nanmean(prcNorm__{iDelay,iTarget});
        e(iDelay,iTarget) = nansem(prcNorm__{iDelay,iTarget});
    end
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
%%
rmpath('/Users/domansa/Dropbox/MATLAB/PC MATLAB Path/functions/nansuite');
        

for iDelay = 1:3
    x = prcNorm__{iDelay,2};
    x(isnan(x))=[];
    [~,p_] = ttest(x,1)
%     p_ = signrank(x,1)
    p(iDelay,1) = p_<(0.05)
    
    x = prcNorm__{iDelay,3};
    x(isnan(x))=[];
    [~,p_] = ttest(x,1)
%     p_ = signrank(x,1)
    p(iDelay,2) = p_<(0.05)
end

addpath('/Users/domansa/Dropbox/MATLAB/PC MATLAB Path/functions/nansuite');
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

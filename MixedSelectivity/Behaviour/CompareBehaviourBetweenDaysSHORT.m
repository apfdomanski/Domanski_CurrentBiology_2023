%% %%%%%% PREAMBLE %%%%%%
clear
Target = 'SHORT';
switch Target
    case 'SHORT'
        Delays_ = {'Delay_0'};
        Delays__ = {'0s'};
    case 'MEDIUM'
        Delays_ = {'Delay_2','Delay_4','Delay_6','Delay_8'};
        Delays__ = {'2s','4s','6s','8s'};
    case 'LONG'
        Delays_ = {'Short','Medium','Long'};
        Delays__ = {'4s','8s','16s'};
end


AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time

tlimsAll = [-5 5];
tlimsShort=[-4 10];
tlimsMedium=[-8 10];
tlimsLong=[-16 10];
tlimsANOVA = [-4 4];

shift = 0;
plotOnline = false;
bw=0.05;
Nbs = 500;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbANOVA=tlimsANOVA(1):bw:tlimsANOVA(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
pat = 'C:\Analysis\AssemblyAnalysis\raw';
cd(pat)
fileList=dir(sprintf('allTimestamps\\*%s*.mat',Target));
fileListAss = fileList;

% reject_list={'NorbertMEDIUM1_Events.mat','NorbertMEDIUM2_Events.mat','OnufryMEDIUM1_Events.mat','OnufryMEDIUM2_Events.mat'}; % These only have 5s delays
reject_list=[];
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileListAss(find(name_flag))=[];% fnames_(name_flag)=[];
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

normaliseFscores = false;
normWin = [-5 -3];
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};

%% %%%%%% ANALYSE BEHAVIOUR %%%%%%
%% Behaviour - Response Latency (time from Sample press to lever press)
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
    for iDelay=1:length(Delays_)
        Delay_ = Delays_{iDelay};
        
        % Sample latency
        eval(sprintf('Slatency.%s.Correct{iFile} =  [(t.%s.SamplePress_LeftCorrect-t.%s.CueLight_LeftCorrect)/1e6;(t.%s.SamplePress_RightCorrect-t.%s.CueLight_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))

        % Nosepoke latency
        eval(sprintf('NPlatency.%s.Correct{iFile} =  [(t.%s.NosePoke_LeftCorrect-t.%s.SamplePress_LeftCorrect)/1e6;(t.%s.NosePoke_RightCorrect-t.%s.SamplePress_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
        
        % Response latency
%         eval(sprintf('Rlatency.%s.Correct{iFile} =  [(t.%s.ChoicePress_LeftCorrect-t.%s.DelayEnd_LeftCorrect)/1e6;(t.%s.ChoicePress_RightCorrect-t.%s.DelayEnd_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
        eval(sprintf('Rlatency.%s.Correct{iFile} =  [(t.%s.ChoicePress_LeftCorrect-t.%s.NosePoke_LeftCorrect)/1e6;(t.%s.ChoicePress_RightCorrect-t.%s.NosePoke_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
        
        % Errors trials : This bugs out because some of the error trials are omissions,throwing out the alingmet
        try 
            eval(sprintf('Slatency.%s.Error{iFile} =  [(t.%s.SamplePress_LeftError-t.%s.CueLight_LeftError)/1e6,(t.%s.SamplePress_RightError-t.%s.CueLight_RightError)/1e6]'';;',Delay_,Delay_,Delay_,Delay_,Delay_))
        end
        try
            eval(sprintf('NPlatency.%s.Error{iFile} =  [(t.%s.NosePoke_LeftError-t.%s.SamplePress_LeftError)/1e6,(t.%s.NosePoke_RightError-t.%s.SamplePress_RightError)/1e6]'';;',Delay_,Delay_,Delay_,Delay_,Delay_))
        end
        try
            eval(sprintf('Rlatency.%s.Error{iFile} =  [(t.%s.ChoicePress_LeftError-t.%s.NosePoke_LeftError)/1e6,(t.%s.ChoicePress_RightError-t.%s.NosePoke_RightError)/1e6]'';;',Delay_,Delay_,Delay_,Delay_,Delay_))
        end

    end
end
    save(sprintf('%s\\allTimestamps\\Behaviour analysis\\DayByDay_SHORT.mat',pat),'Slatency','NPlatency','Rlatency');

%% Plot response latency distribution
bins = 0:0.5:30;
figure; hold on
x=histc(cell2mat(Rlatency.Delay_0.Correct'),bins);
area(bins,x./sum(x)*100,'Facecolor',[0.1 0.1 0.1],'FaceAlpha',0.3,'LineStyle','none')

bins = 0:0.1:30;
x=histc(cell2mat(NPlatency.Delay_0.Correct'),bins);
plot(bins,x./sum(x)*100,'color',[0.2 0.7 0.2 0.4],'LineStyle','-')

plot([0 0],[0 6],':k')
plot([10 10],[0 6],':k')
text(0,6.5,'Tone','Rotation',90,'HorizontalAlignment','left')
text(10,6.5,'Time out','Rotation',90,'HorizontalAlignment','left')
axis([-1 25 0 50])
% [bins(x==max(x)),bins(y==max(y)),bins(z==max(z))]
% [median(cell2mat(Rlatency.Short.Correct')),median(cell2mat(Rlatency.Medium.Correct')),median(cell2mat(Rlatency.Long.Correct'));]

legend([cellfun(@(x) [x ' Choices'], Delays__,'UniformOutput',false),cellfun(@(x) [x ' Nose-pokes'], Delays__,'UniformOutput',false)]); legend boxoff

ylabel('Percentage of trials')
xlabel('Choice latency (s)')
% kruskalwallis([cell2mat(NPlatency.Short.Correct');cell2mat(NPlatency.Medium.Correct');cell2mat(NPlatency.Long.Correct')],...
%               [ones(length(cell2mat(NPlatency.Short.Correct')),1);2*ones(length(cell2mat(NPlatency.Medium.Correct')),1);3*ones(length(cell2mat(NPlatency.Long.Correct')),1)])
% 
% [p,tbl,stats] = kruskalwallis([cell2mat(Rlatency.Short.Correct');cell2mat(Rlatency.Medium.Correct');cell2mat(Rlatency.Long.Correct')],...
%               [ones(length(cell2mat(Rlatency.Short.Correct')),1);2*ones(length(cell2mat(Rlatency.Medium.Correct')),1);3*ones(length(cell2mat(Rlatency.Long.Correct')),1)])
% c = multcompare(stats);
%% Plot sample latency distribution

bins = 0:0.5:20;
figure; hold on
x=histc(cell2mat(Slatency.Delay_0.Correct'),bins);
area(bins,x./sum(x)*100,'Facecolor',[0.1 0.1 0.7],'FaceAlpha',0.3,'LineStyle','none')

plot([0 0],[0 6],':k')
plot([20 20],[0 6],':k')
text(0,6.5,'Cue light','Rotation',90,'HorizontalAlignment','left')
text(20,6.5,'Time out','Rotation',90,'HorizontalAlignment','left')
axis([-1 21 0 20])
% [bins(x==max(x)),bins(y==max(y)),bins(z==max(z))];
% [median(cell2mat(Rlatency.Short.Correct')),median(cell2mat(Rlatency.Medium.Correct')),median(cell2mat(Rlatency.Long.Correct'));]
legend(cellfun(@(x) [x ' Trials'], Delays__,'UniformOutput',false)); legend boxoff

ylabel('Percentage of trials')
xlabel('Sample latency (s)')
% kruskalwallis([cell2mat(Slatency.Short.Correct');cell2mat(Slatency.Medium.Correct');cell2mat(Slatency.Long.Correct')],...
%               [ones(length(cell2mat(Slatency.Short.Correct')),1);2*ones(length(cell2mat(Slatency.Medium.Correct')),1);3*ones(length(cell2mat(Slatency.Long.Correct')),1)])
%% Behaviour - Plot all sorted correct trials
Trialraster.Delay_0 = [];
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
    
    for iDelay=1:length(Delays_)
        Delay_ = Delays_{iDelay};
        
        eval(sprintf('t_ = [t.%s.SamplePress_LeftCorrect,t.%s.NosePoke_LeftCorrect,t.%s.ChoicePress_LeftCorrect;t.%s.SamplePress_RightCorrect,t.%s.NosePoke_RightCorrect,t.%s.ChoicePress_RightCorrect;]*1e-6;',Delay_,Delay_,Delay_,Delay_,Delay_,Delay_));

%         eval(sprintf('t_ = [t.%s.NosePoke_LeftCorrect,t.%s.ChoicePress_LeftCorrect;t.%s.NosePoke_RightCorrect,t.%s.ChoicePress_RightCorrect;]*1e-6;',Delay_,Delay_,Delay_,Delay_));
        t_ = t_-t_(:,1);
        eval(sprintf('Trialraster.%s = [Trialraster.%s ;t_];',Delay_,Delay_))
    end
end

for iDelay=1:length(Delays_)
    Delay_ = Delays_{iDelay};
    eval(sprintf('Trialraster.%s  = sortrows(Trialraster.%s,2);',Delay_,Delay_))
end

Ns = size(Trialraster.Delay_0,1);
X =  Trialraster.Delay_0;
Y = repmat((1:size(X,1))',1,3);
figure('color','w'); hold on

scatter(0,-200,[],[0.9 0.4 0.2],'x')
scatter(0,-200,[],[1 0 0],'x')

scatter(X(:,1),Y(:,1),[],0.6*[0 0 0],'x')
scatter(X(:,2),Y(:,2),[],[0.9 0.4 0.2],'x')
scatter(X(:,3),Y(:,3),[],[1 0 0],'.')
% %%
% 
plot([-3 0],[-40 -40],'k','LineWidth',2)
plot([0 10],[-40 -40],':k','LineWidth',2)
plot([10 15],[-40 -40],'k','LineWidth',2)

plot([12 14],[150 150],'k','LineWidth',2)
plot([12 12],[150 250],'k','LineWidth',2)
text(13,120,'2s','HorizontalAlignment','center')
text(11,200,'100 trials','HorizontalAlignment','center','Rotation',90)
text(-1,10,'Sample Press','Rotation',90)

axis off

axis([-5 16 -50 800])

% plot([28 28],[0 Ns(1)],'color',[0.1 0.1 0.7 0.3],'LineWidth',2)
% plot([28 28],[Ns(1) Ns(1)+Ns(2) ],'color',[0.2 0.7 0.2 0.4],'LineWidth',2)
% plot([28 28],[Ns(1)+Ns(2) Ns(1)+Ns(2)+Ns(3)],'color',[0.7 0.2 0.2 0.5],'LineWidth',2)
% plot([28 28],[Ns(1)+Ns(2)+Ns(3) sum(Ns)],'color',[0.7 0.5 0.6 0.5],'LineWidth',2)
% 
% text(29,Ns(1)/2,'2s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.1 0.1 0.7 0.3])
% text(29,Ns(1) + Ns(2)/2,'4s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.2 0.7 0.2 0.4])
% text(29,Ns(1)+Ns(2) + Ns(3)/2,'6s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.7 0.2 0.2 0.5])
% text(29,Ns(1)+Ns(2) + Ns(3) + Ns(4)/2,'8s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.7 0.5 0.6 0.5])
% legend('Nose-poke','Choice Press','Location','southoutside','Orientation','horizontal'); legend boxoff
%% Behaviour - Plot all sorted error trials
TrialrasterError.Delay_0 = [];
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
    
    for iDelay=1:length(Delays_)
        Delay_ = Delays_{iDelay};
        
        eval(sprintf('t_ = [t.%s.SamplePress_LeftError'',t.%s.NosePoke_LeftError'',t.%s.ChoicePress_LeftError'';t.%s.SamplePress_RightError'',t.%s.NosePoke_RightError'',t.%s.ChoicePress_RightError'';]*1e-6;',Delay_,Delay_,Delay_,Delay_,Delay_,Delay_));
        t_ = t_-t_(:,1);
        eval(sprintf('TrialrasterError.%s = [TrialrasterError.%s ;t_];',Delay_,Delay_))
    end
end

for iDelay=1:length(Delays_)
    Delay_ = Delays_{iDelay};
    eval(sprintf('TrialrasterError.%s  = sortrows(TrialrasterError.%s,2);',Delay_,Delay_))
end

Ns = size(TrialrasterError.Delay_0,1);
X =  TrialrasterError.Delay_0;
Y = repmat((1:size(X,1))',1,3);
figure('color','w'); hold on

scatter(0,-200,[],[0.9 0.4 0.2],'x')
scatter(0,-200,[],[1 0 0],'x')

scatter(X(:,1),Y(:,1),[],0.6*[0 0 0],'x')
scatter(X(:,2),Y(:,2),[],[0.9 0.4 0.2],'x')
scatter(X(:,3),Y(:,3),[],[1 0 0],'.')

% %%
% 
plot([-3 0],[-40 -40],'k','LineWidth',2)
plot([0 25],[-40 -40],':k','LineWidth',2)
plot([25 30],[-40 -40],'k','LineWidth',2)

plot([12 14],[0 0],'k','LineWidth',2)
plot([12 12],[0 20],'k','LineWidth',2)
text(13,-10,'2s','HorizontalAlignment','center')
text(11,10,'100 trials','HorizontalAlignment','center','Rotation',90)
text(-1,10,'Sample Press','Rotation',90)


axis off

axis([-5 30 -50 200])

% plot([28 28],[0 Ns(1)],'color',[0.1 0.1 0.7 0.3],'LineWidth',2)
% plot([28 28],[Ns(1) Ns(1)+Ns(2) ],'color',[0.2 0.7 0.2 0.4],'LineWidth',2)
% plot([28 28],[Ns(1)+Ns(2) Ns(1)+Ns(2)+Ns(3)],'color',[0.7 0.2 0.2 0.5],'LineWidth',2)
% plot([28 28],[Ns(1)+Ns(2)+Ns(3) sum(Ns)],'color',[0.7 0.5 0.6 0.5],'LineWidth',2)
% 
% text(29,Ns(1)/2,'2s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.1 0.1 0.7 0.3])
% text(29,Ns(1) + Ns(2)/2,'4s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.2 0.7 0.2 0.4])
% text(29,Ns(1)+Ns(2) + Ns(3)/2,'6s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.7 0.2 0.2 0.5])
% text(29,Ns(1)+Ns(2) + Ns(3) + Ns(4)/2,'8s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.7 0.5 0.6 0.5])
legend('Nose-poke','Choice Press','Location','southoutside','Orientation','horizontal'); legend boxoff
%% Cumulative distributions of latencies - correct trials
figure
x1= diff(Trialraster.Delay_0,[],2);
subplot(1,2,1); hold on
    bins =0:0.2:15;
    x = hist(x1(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.1 0.1 0.7],'LineWidth',1.5)
    x = hist(x2(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.2 0.7 0.2],'LineWidth',1.5)
    x = hist(x3(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.7 0.2 0.2],'LineWidth',1.5)
    x = hist(x4(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.7 0.5 0.6],'LineWidth',1.5)
    axis square
    ylabel('Cumulative fraction of trials')
    xlabel('Tone-Nosepoke latency (s)')
subplot(1,2,2); hold on
    bins =0:0.2:10;
    x = hist(x1(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.1 0.1 0.7],'LineWidth',1.5)
    x = hist(x2(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.2 0.7 0.2],'LineWidth',1.5)
    x = hist(x3(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.7 0.2 0.2],'LineWidth',1.5)
    x = hist(x4(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.7 0.5 0.6],'LineWidth',1.5)
    axis square
    xlabel('Nosepoke-Choice latency (s)')

%     hold on
%     bins =-20:20;
%     x = hist(diff(x1,[],2),bins);x=(x)./sum(x);
%     plot(bins,x,'Color',[0.1 0.1 0.7],'LineWidth',1.5)
%     x = hist(diff(x2,[],2),bins);x=(x)./sum(x);
%     plot(bins,x,'Color',[0.2 0.7 0.2],'LineWidth',1.5)
%     x = hist(diff(x3,[],2),bins);x=(x)./sum(x);
%     plot(bins,x,'Color',[0.7 0.2 0.2],'LineWidth',1.5)
%     axis square
%     xlabel('Nosepoke-Choice latency (s)')
    
    
legend([cellfun(@(x) [x ' Trials'], Delays__,'UniformOutput',false),cellfun(@(x) [x ' Nose-pokes'], Delays__,'UniformOutput',false)],'Location','southeast'); legend boxoff
%% Cumulative distributions of latencies - error trials
% figure
% x1= diff(TrialrasterError.Delay_0,[],2);
% subplot(1,2,1); hold on
%     bins =0:0.2:25;
%     x = hist(x1(:,1),bins);x=cumsum(x)./sum(x);
%     plot(bins,x,'Color',[0.1 0.1 0.7],'LineWidth',1.5)
%     axis square
%     ylabel('Cumulative fraction of trials')
%     xlabel('Tone-Nosepoke latency (s)')
% subplot(1,2,2); hold on
%     bins =0:0.2:10;
%     x = hist(x1(:,2),bins);x=cumsum(x)./sum(x);
%     plot(bins,x,'Color',[0.1 0.1 0.7],'LineWidth',1.5)
%     axis square
%     xlabel('Nosepoke-Choice latency (s)')

%     hold on
%     bins =-20:20;
%     x = hist(diff(x1,[],2),bins);x=(x)./sum(x);
%     plot(bins,x,'Color',[0.1 0.1 0.7],'LineWidth',1.5)
%     x = hist(diff(x2,[],2),bins);x=(x)./sum(x);
%     plot(bins,x,'Color',[0.2 0.7 0.2],'LineWidth',1.5)
%     x = hist(diff(x3,[],2),bins);x=(x)./sum(x);
%     plot(bins,x,'Color',[0.7 0.2 0.2],'LineWidth',1.5)
%     axis square
%     xlabel('Nosepoke-Choice latency (s)')
    
    
% legend([cellfun(@(x) [x ' Trials'], Delays__,'UniformOutput',false),cellfun(@(x) [x ' Nose-pokes'], Delays__,'UniformOutput',false)],'Location','southeast'); legend boxoff

x1=  diff(Trialraster.Delay_0,[],2);
x1e= diff(TrialrasterError.Delay_0,[],2);
[h(1,1),p(1,1)] = kstest2(x1(:,1),x1e(:,1));
[h(2,1),p(2,1)] = kstest2(x1(:,2),x1e(:,2));
figure('color','w','Name','Tone to nosepoke latency');

subplot(1,2,1); hold on
    bins =0:0.2:20;
    x = hist(x1(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x1e(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
%     title('Tone - Nosepoke Delay')
    axis square
    ylabel('Cumulative fraction of trials')
    xlabel('Latency (s)')
    title('Nosepoke Latency')
    
    
subplot(1,2,2); hold on
    bins =0:0.2:20;
    x = hist(x1(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x1e(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
%     title('Tone - Nosepoke Delay')
    text(max(bins),0.1,'No delay trials','HorizontalAlignment','right','color','k')
    axis square
    xlim([0 10])
    xlabel('Latency (s)')
    title('Choice Latency')
    legend('Correct','Error','Location','southeast'); legend boxoff
%% Sample vs Choice Latency - points
figure;
subplot(1,2,1);hold on
for iFile =1:length(fileList)
    scatter(Slatency.Delay_0.Correct{iFile},Rlatency.Delay_0.Correct{iFile},10,'k','filled')
    try scatter(Slatency.Delay_0.Error{iFile},Rlatency.Delay_0.Error{iFile},10,'r','filled'); end
   
end
axis square
axis([0 20 0 10])
ylabel('Choice press latency (s)')
xlabel('Sample press latency (s)')

subplot(1,2,2);hold on
for iFile =1:length(fileList)
    plot([mean(Slatency.Delay_0.Correct{iFile}),mean(Slatency.Delay_0.Error{iFile})],[mean(Rlatency.Delay_0.Correct{iFile}), mean(Rlatency.Delay_0.Error{iFile})],'k')

    H = errorbarxy(mean(Slatency.Delay_0.Correct{iFile}),mean(Rlatency.Delay_0.Correct{iFile}),...
               nansem(Slatency.Delay_0.Correct{iFile}),nansem(Rlatency.Delay_0.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = 'k';
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.1 0.1 0.1]);
    end
    
    H = errorbarxy(mean(Slatency.Delay_0.Error{iFile}),mean(Rlatency.Delay_0.Error{iFile}),...
               nansem(Slatency.Delay_0.Error{iFile}),nansem(Rlatency.Delay_0.Error{iFile}),{'r','r','r'});
    H.hMain.Color = 'r';    
    for i=1:6
        set(H.hErrorbar(i),'Color','r');
    end
end  
xlabel('Sample latency (s)')
ylabel('Choice Latency (s)')
axis square
axis([0 11 0 5])
%     legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','south'); legend boxoff
  
    %% bars of mean latencies - SHORT
figure; 
subplot(1,3,1);hold on
   
    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=nanmean(Slatency.Delay_0.Correct{iFile});

    end
    bar(nanmean(x,2),'FaceColor','w','EdgeColor','k','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.k','Marker','none','LineWidth',1.5)
 
    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=nanmean(Slatency.Delay_0.Error{iFile});
    end
    bar(mean(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar(mean(x,2),nansem(x,2),'.r','Marker','none','LineWidth',1.5)

    set(gca,'XTick',[1 ],'XTickLabel',Delays__)
    xlabel('Delay duration')
    ylabel('Sample Latency (s)')
    ylim([0 10])

subplot(1,3,2);hold on
    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=nanmean(NPlatency.Delay_0.Error{iFile});
    end
    bar(nanmean(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.r','Marker','none','LineWidth',1.5)

    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=nanmean(NPlatency.Delay_0.Correct{iFile});
    end
    bar(nanmean(x,2),'FaceColor','w','EdgeColor','k','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.k','Marker','none','LineWidth',1.5)
    
    set(gca,'XTIck',[1],'XTickLabel',Delays__)
    xlabel('Delay duration')
    ylabel('Nosepoke Latency (s)')
    ylim([0 6.5])

subplot(1,3,3);hold on
    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=nanmean(Rlatency.Delay_0.Error{iFile});
    end
    bar(mean(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.r','Marker','none','LineWidth',1.5)
    
    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=nanmean(Rlatency.Delay_0.Correct{iFile});
    end
    bar(nanmean(x,2),'FaceColor','none','EdgeColor','k','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.k','Marker','none','LineWidth',1.5)
    
    set(gca,'XTIck',[1],'XTickLabel',Delays__)
    xlabel('Delay duration')
    ylabel('Lever-press Latency (s)')
    ylim([0 6.5])
    
    
    
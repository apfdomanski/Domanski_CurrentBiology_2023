%% %%%%%% PREAMBLE %%%%%%
clear
Target = 'MEDIUM';
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

reject_list={'NorbertMEDIUM1_Events.mat','NorbertMEDIUM2_Events.mat','OnufryMEDIUM1_Events.mat','OnufryMEDIUM2_Events.mat'}; % These only have 5s delays
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
%% Behaviour - Response Latency (time from end-of-delay signal to lever press)
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
    for iDelay=1:length(Delays_)
        Delay_ = Delays_{iDelay};
        
        % Sample latency
        eval(sprintf('Slatency.%s.Correct{iFile} =  [(t.%s.SamplePress_LeftCorrect-t.%s.CueLight_LeftCorrect)/1e6;(t.%s.SamplePress_RightCorrect-t.%s.CueLight_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))

        % Nosepoke latency
        eval(sprintf('NPlatency.%s.Correct{iFile} =  [(t.%s.NosePoke_LeftCorrect-t.%s.DelayEnd_LeftCorrect)/1e6;(t.%s.NosePoke_RightCorrect-t.%s.DelayEnd_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
        
        % Response latency
%         eval(sprintf('Rlatency.%s.Correct{iFile} =  [(t.%s.ChoicePress_LeftCorrect-t.%s.DelayEnd_LeftCorrect)/1e6;(t.%s.ChoicePress_RightCorrect-t.%s.DelayEnd_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
        eval(sprintf('Rlatency.%s.Correct{iFile} =  [(t.%s.ChoicePress_LeftCorrect-t.%s.NosePoke_LeftCorrect)/1e6;(t.%s.ChoicePress_RightCorrect-t.%s.NosePoke_RightCorrect)/1e6];',Delay_,Delay_,Delay_,Delay_,Delay_))
        
        % Errors trials : This bugs out because some of the error trials are omissions,throwing out the alingmet
        try
            eval(sprintf('NPlatency.%s.Error{iFile} =  [(t.%s.NosePoke_LeftError-t.%s.DelayEnd_LeftError)/1e6,(t.%s.NosePoke_RightError-t.%s.DelayEnd_RightError)/1e6]'';;',Delay_,Delay_,Delay_,Delay_,Delay_))
        end
        try
            eval(sprintf('Rlatency.%s.Error{iFile} =  [(t.%s.ChoicePress_LeftError-t.%s.NosePoke_LeftError)/1e6,(t.%s.ChoicePress_RightError-t.%s.NosePoke_RightError)/1e6]'';;',Delay_,Delay_,Delay_,Delay_,Delay_))
        end
        try 
            eval(sprintf('Slatency.%s.Error{iFile} =  [(t.%s.SamplePress_LeftError-t.%s.CueLight_LeftError)/1e6,(t.%s.SamplePress_RightError-t.%s.CueLight_RightError)/1e6]'';;',Delay_,Delay_,Delay_,Delay_,Delay_))
        end
    end
end
    save(sprintf('%s\\allTimestamps\\Behaviour analysis\\DayByDay_MEDIUM.mat',pat),'Slatency','NPlatency','Rlatency');

%%  Plot Nosepoke and Choice Response latency distribution
bins = 0:0.5:30;
figure; hold on
x=histc(cell2mat(Rlatency.Delay_4.Correct'),bins);
area(bins,x./sum(x)*100,'Facecolor',[0.1 0.1 0.7],'FaceAlpha',0.3,'LineStyle','none')
y=histc(cell2mat(Rlatency.Delay_8.Correct'),bins);
area(bins,y./sum(y)*100,'Facecolor',[0.2 0.7 0.2],'FaceAlpha',0.4,'LineStyle','none')

bins = 0:0.1:30;
x=histc(cell2mat(NPlatency.Delay_4.Correct'),bins);
plot(bins,x./sum(x)*100,'color',[0.1 0.1 0.7 0.3],'LineStyle','-')

x=histc(cell2mat(NPlatency.Delay_8.Correct'),bins);
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
%%  Plot Sample latency distribution


bins = 0:0.5:20;
figure; hold on
x=histc(cell2mat(Slatency.Delay_2.Correct'),bins);
area(bins,x./sum(x)*100,'Facecolor',[0.1 0.1 0.7],'FaceAlpha',0.3,'LineStyle','none')
y=histc(cell2mat(Slatency.Delay_4.Correct'),bins);
area(bins,y./sum(y)*100,'Facecolor',[0.2 0.7 0.2],'FaceAlpha',0.4,'LineStyle','none')

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
Trialraster.Delay_2 = [];
Trialraster.Delay_4 = [];
Trialraster.Delay_6 = [];
Trialraster.Delay_8 = [];

for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
    
    for iDelay=1:length(Delays_)
        Delay_ = Delays_{iDelay};
        
        eval(sprintf('t_ = [t.%s.DelayEnd_LeftCorrect,t.%s.NosePoke_LeftCorrect,t.%s.ChoicePress_LeftCorrect;t.%s.DelayEnd_RightCorrect,t.%s.NosePoke_RightCorrect,t.%s.ChoicePress_RightCorrect;]*1e-6;',Delay_,Delay_,Delay_,Delay_,Delay_,Delay_));
        t_ = t_-t_(:,1);
        eval(sprintf('Trialraster.%s = [Trialraster.%s ;t_];',Delay_,Delay_))
    end
    

    
  
end

for iDelay=1:length(Delays_)
    Delay_ = Delays_{iDelay};
    eval(sprintf('Trialraster.%s  = sortrows(Trialraster.%s,2);',Delay_,Delay_))
end

Ns = [size(Trialraster.Delay_2,1),...
      size(Trialraster.Delay_4,1),...
      size(Trialraster.Delay_6,1),...
      size(Trialraster.Delay_8,1)];
X = [Trialraster.Delay_2;...
     Trialraster.Delay_4;...
     Trialraster.Delay_6;...
     Trialraster.Delay_8];
Y = repmat((1:size(X,1))',1,3);
figure('color','w'); hold on

scatter(0,-200,[],[0.9 0.4 0.2],'x')
scatter(0,-200,[],[1 0 0],'.')

scatter(X(:,1),Y(:,1),[],0.6*[0 0 0],'.')
scatter(X(:,2),Y(:,2),[],[0.9 0.4 0.2],'x')
scatter(X(:,3),Y(:,3),[],[1 0 0],'.')

plot([-3 0],[-40 -40],'k','LineWidth',2)
plot([0 25],[-40 -40],':k','LineWidth',2)
plot([25 28],[-40 -40],'k','LineWidth',2)


plot([20 25],[150 150],'k','LineWidth',2)
plot([20 20],[150 250],'k','LineWidth',2)
text(22.5,120,'5s','HorizontalAlignment','center')
text(19,200,'100 trials','HorizontalAlignment','center','Rotation',90)
text(-1,10,'End-of-Delay tone','Rotation',90)

axis off

axis([-5 30 -50 400])

plot([28 28],[0 Ns(1)],'color',[0.1 0.1 0.7 0.3],'LineWidth',2)
plot([28 28],[Ns(1) Ns(1)+Ns(2) ],'color',[0.2 0.7 0.2 0.4],'LineWidth',2)
plot([28 28],[Ns(1)+Ns(2) Ns(1)+Ns(2)+Ns(3)],'color',[0.7 0.2 0.2 0.5],'LineWidth',2)
plot([28 28],[Ns(1)+Ns(2)+Ns(3) sum(Ns)],'color',[0.7 0.5 0.6 0.5],'LineWidth',2)

text(29,Ns(1)/2,'2s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.1 0.1 0.7 0.3])
text(29,Ns(1) + Ns(2)/2,'4s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.2 0.7 0.2 0.4])
text(29,Ns(1)+Ns(2) + Ns(3)/2,'6s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.7 0.2 0.2 0.5])
text(29,Ns(1)+Ns(2) + Ns(3) + Ns(4)/2,'8s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.7 0.5 0.6 0.5])
legend('Nose-poke','Choice Press','Location','southoutside','Orientation','horizontal'); legend boxoff
%% Behaviour - Plot all sorted error trials
TrialrasterError.Delay_2 = [];
TrialrasterError.Delay_4 = [];
TrialrasterError.Delay_6 = [];
TrialrasterError.Delay_8 = [];

for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
    for iDelay=1:length(Delays_)
        Delay_ = Delays_{iDelay};
        
        eval(sprintf('t_ = [t.%s.DelayEnd_LeftError,t.%s.DelayEnd_RightError;t.%s.NosePoke_LeftError,t.%s.NosePoke_RightError;t.%s.ChoicePress_LeftError,t.%s.ChoicePress_RightError]''*1e-6;',Delay_,Delay_,Delay_,Delay_,Delay_,Delay_))
        
        t_ = t_-t_(:,1);
        eval(sprintf('TrialrasterError.%s = [TrialrasterError.%s ;t_];',Delay_,Delay_))
        
    end
    
   
end
for iDelay=1:length(Delays_)
    Delay_ = Delays_{iDelay};
    eval(sprintf('TrialrasterError.%s  = sortrows(TrialrasterError.%s,2);',Delay_,Delay_))
end

Ns = [size(TrialrasterError.Delay_2,1),...
      size(TrialrasterError.Delay_4,1),...
      size(TrialrasterError.Delay_6,1),...
      size(TrialrasterError.Delay_8,1)];
X = [TrialrasterError.Delay_2;...
     TrialrasterError.Delay_4;...
     TrialrasterError.Delay_6;...
     TrialrasterError.Delay_8];
Y = repmat((1:size(X,1))',1,3);
figure('color','w'); hold on

scatter(0,-200,[],[0.9 0.4 0.2],'x')
scatter(0,-200,[],[1 0 0],'.')

scatter(X(:,1),Y(:,1),[],0.6*[0 0 0],'.')
scatter(X(:,2),Y(:,2),[],[0.9 0.4 0.2],'x')
scatter(X(:,3),Y(:,3),[],[1 0 0],'.')

plot([-3 0],[-40 -40],'k','LineWidth',2)
plot([0 25],[-40 -40],':k','LineWidth',2)
plot([25 28],[-40 -40],'k','LineWidth',2)


plot([20 25],[150 150],'k','LineWidth',2)
plot([20 20],[150 250],'k','LineWidth',2)
text(22.5,120,'5s','HorizontalAlignment','center')
text(19,200,'100 trials','HorizontalAlignment','center','Rotation',90)
text(-1,10,'End-of-Delay tone','Rotation',90)

axis off

axis([-5 30 -50 300])

plot([28 28],[0 Ns(1)],'color',[0.1 0.1 0.7 0.3],'LineWidth',2)
plot([28 28],[Ns(1) Ns(1)+Ns(2) ],'color',[0.2 0.7 0.2 0.4],'LineWidth',2)
plot([28 28],[Ns(1)+Ns(2) sum(Ns)],'color',[0.7 0.2 0.2 0.5],'LineWidth',2)
text(29,Ns(1)/2,'4s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.1 0.1 0.7 0.3])
text(29,Ns(1) + Ns(2)/2,'8s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.2 0.7 0.2 0.4])
text(29,Ns(1)+Ns(2) + Ns(3)/2,'16s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.7 0.2 0.2 0.5])
legend('Nose-poke','Choice Press','Location','southoutside','Orientation','horizontal'); legend boxoff
%% Cumulative distributions of latencies - correct trials
figure
x1= diff(Trialraster.Delay_2,[],2);
x2= diff(Trialraster.Delay_4,[],2);
x3= diff(Trialraster.Delay_6,[],2);
x4= diff(Trialraster.Delay_8,[],2);
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
figure
x1= diff(TrialrasterError.Delay_2,[],2);
x2= diff(TrialrasterError.Delay_4,[],2);
x3= diff(TrialrasterError.Delay_6,[],2);
x4= diff(TrialrasterError.Delay_8,[],2);
subplot(1,2,1); hold on
    bins =0:0.2:25;
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

x1= diff(Trialraster.Delay_2,[],2);
x2= diff(Trialraster.Delay_4,[],2);
x3= diff(Trialraster.Delay_6,[],2);
x4= diff(Trialraster.Delay_8,[],2);
x1e= diff(TrialrasterError.Delay_2,[],2);
x2e= diff(TrialrasterError.Delay_4,[],2);
x3e= diff(TrialrasterError.Delay_6,[],2);
x4e= diff(TrialrasterError.Delay_8,[],2);
    [h(1,1),p(1,1)] = kstest2(x1(:,1),x1e(:,1));
    [h(1,2),p(1,2)] = kstest2(x2(:,1),x2e(:,1));
    [h(1,3),p(1,3)] = kstest2(x3(:,1),x3e(:,1));
    [h(1,4),p(1,4)] = kstest2(x4(:,1),x4e(:,1));

    [h(2,1),p(2,1)] = kstest2(x1(:,2),x1e(:,2));
    [h(2,2),p(2,2)] = kstest2(x2(:,2),x2e(:,2));
    [h(2,3),p(2,3)] = kstest2(x3(:,2),x3e(:,2));
    [h(2,4),p(2,4)] = kstest2(x4(:,2),x4e(:,2));

figure('color','w','Name','Tone to nosepoke latency');

subplot(1,4,1); hold on
    bins =0:0.2:20;
    x = hist(x1(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x1e(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
%     title('Tone - Nosepoke Delay')
    text(max(bins),0.1,'2s delay trials','HorizontalAlignment','right','color',[0.1 0.1 0.7])
    axis square
    ylabel('Cumulative fraction of trials')
    xlabel('Latency (s)')
    
subplot(1,4,2); hold on
    x = hist(x2(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x2e(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
    text(max(bins),0.1,'4s delay trials','HorizontalAlignment','right','color',[0.2 0.7 0.2])
    xlabel('Latency (s)')
    axis square

subplot(1,4,3); hold on
    x = hist(x3(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x3e(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)      
    text(max(bins),0.1,'6s delay trials','HorizontalAlignment','right','color',[0.7 0.2 0.2])
    xlabel('Latency (s)')
    axis square

subplot(1,4,4); hold on
    x = hist(x4(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x4e(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)      
    text(max(bins),0.1,'8s delay trials','HorizontalAlignment','right','color',[0.7 0.5 0.6])
    xlabel('Latency (s)')
    axis square
    
    
figure('color','w','Name','Nosepoke to choice latency');

subplot(1,4,1); hold on
    bins =0:0.2:20;
    x = hist(x1(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x1e(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
%     title('Tone - Nosepoke Delay')
    text(10,0.1,'2s delay trials','HorizontalAlignment','right','color',[0.1 0.1 0.7])
    axis square
    xlim([0 10])
    ylabel('Cumulative fraction of trials')
    xlabel('Latency (s)')
    
subplot(1,4,2); hold on
    x = hist(x2(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x2e(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
    text(10,0.1,'4s delay trials','HorizontalAlignment','right','color',[0.2 0.7 0.2])
    xlabel('Latency (s)')
    axis square
    xlim([0 10])

subplot(1,4,3); hold on
    x = hist(x3(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x3e(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)      
    text(10,0.1,'6s delay trials','HorizontalAlignment','right','color',[0.7 0.2 0.2])
    xlabel('Latency (s)')
    axis square
    xlim([0 10])

subplot(1,4,4); hold on
    x = hist(x4(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x4e(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)      
    text(10,0.1,'8s delay trials','HorizontalAlignment','right','color',[0.7 0.5 0.6])
    xlabel('Latency (s)')
    axis square
    xlim([0 10])
%% Sample vs Choice Latency - points
figure;
for iFile =1:length(fileList)
    subplot(1,4,1);hold on
    scatter(Slatency.Delay_2.Correct{iFile},Rlatency.Delay_2.Correct{iFile},10,[0.1 0.1 0.7],'filled')
    try scatter(Slatency.Delay_2.Error{iFile},Rlatency.Delay_2.Error{iFile},10,[0.1 0.1 0.7]); end
    axis square
    axis([0 30 0 30])
    ylabel('Choice press latency (s)')
    
    subplot(1,4,2);hold on
    scatter(Slatency.Delay_4.Correct{iFile},Rlatency.Delay_4.Correct{iFile},10,[0.2 0.7 0.2],'filled')
    try scatter(Slatency.Delay_4.Error{iFile},Rlatency.Delay_4.Error{iFile},10,[0.2 0.7 0.2]); end
    axis square
    axis([0 30 0 30])
    xlabel('Sample press latency (s)')
    
    subplot(1,4,3);hold on
    scatter(Slatency.Delay_6.Correct{iFile},Rlatency.Delay_6.Correct{iFile},10,[0.7 0.2 0.2],'filled')
    try scatter(Slatency.Delay_6.Error{iFile},Rlatency.Delay_6.Error{iFile},10,[0.7 0.2 0.2]); end
    axis square
    axis([0 30 0 30])
    xlabel('Sample press latency (s)')
    
    subplot(1,4,4);hold on
    scatter(Slatency.Delay_8.Correct{iFile},Rlatency.Delay_8.Correct{iFile},10,[0.7 0.5 0.6],'filled')
    try scatter(Slatency.Delay_8.Error{iFile},Rlatency.Delay_8.Error{iFile},10,[0.7 0.5 0.6]); end
    axis square
    axis([0 30 0 30])
    xlabel('Sample press latency (s)')
end
%% Sample vs Choice Latency - errorbars 
figure; hold on
for iFile =1:length(fileList)
    subplot(1,4,1);hold on
    H = errorbarxy(mean(Slatency.Delay_2.Correct{iFile}),mean(Rlatency.Delay_2.Correct{iFile}),...
               nansem(Slatency.Delay_2.Correct{iFile}),nansem(Rlatency.Delay_2.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.1 0.1 0.7];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.1 0.1 0.7]);
    end
    xlabel('Sample latency (s)')
    ylabel('Choice Latency (s)')
    axis square
    axis([0 15 0 5])
%     legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','south'); legend boxoff
    
    subplot(1,4,2);hold on
    H = errorbarxy(mean(Slatency.Delay_4.Correct{iFile}),mean(Rlatency.Delay_4.Correct{iFile}),...
               nansem(Slatency.Delay_4.Correct{iFile}),nansem(Rlatency.Delay_4.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.2 0.7 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.2 0.7 0.2]);
    end    
    xlabel('Sample latency (s)')
    ylabel('Choice Latency (s)')
    axis square
    axis([0 15 0 5])
%     legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','south'); legend boxoff

    subplot(1,4,3);hold on    
    H = errorbarxy(mean(Slatency.Delay_6.Correct{iFile}),mean(Rlatency.Delay_6.Correct{iFile}),...
               nansem(Slatency.Delay_6.Correct{iFile}),nansem(Rlatency.Delay_6.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.7 0.2 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.7 0.2 0.2]);
    end
    xlabel('Sample latency (s)')
    ylabel('Choice Latency (s)')
    axis square
    axis([0 15 0 5])
%     legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','south'); legend boxoff
    
    subplot(1,4,4);hold on
    H = errorbarxy(mean(Slatency.Delay_8.Correct{iFile}),mean(Rlatency.Delay_8.Correct{iFile}),...
               nansem(Slatency.Delay_8.Correct{iFile}),nansem(Rlatency.Delay_8.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.7 0.5 0.6];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.7 0.5 0.6]);
    end    
    
    xlabel('Sample latency (s)')
    ylabel('Choice Latency (s)')
    axis square
    axis([0 15 0 5])
%     legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','south'); legend boxoff
end
%% Sample vs Choice Latency - linked errorbars by animal
figure; hold on
for iFile =1:length(fileList)
    subplot(1,length(fileList),iFile);hold on
    H = errorbarxy(mean(Slatency.Delay_2.Correct{iFile}),mean(Rlatency.Delay_2.Correct{iFile}),...
               nansem(Slatency.Delay_2.Correct{iFile}),nansem(Rlatency.Delay_2.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.1 0.1 0.7];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.1 0.1 0.7]);
    end
   
    
    H = errorbarxy(mean(Slatency.Delay_4.Correct{iFile}),mean(Rlatency.Delay_4.Correct{iFile}),...
               nansem(Slatency.Delay_4.Correct{iFile}),nansem(Rlatency.Delay_4.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.2 0.7 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.2 0.7 0.2]);
    end    
  

    H = errorbarxy(mean(Slatency.Delay_6.Correct{iFile}),mean(Rlatency.Delay_6.Correct{iFile}),...
               nansem(Slatency.Delay_6.Correct{iFile}),nansem(Rlatency.Delay_6.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.7 0.2 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.7 0.2 0.2]);
    end
  
    
    H = errorbarxy(mean(Slatency.Delay_8.Correct{iFile}),mean(Rlatency.Delay_8.Correct{iFile}),...
               nansem(Slatency.Delay_8.Correct{iFile}),nansem(Rlatency.Delay_8.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.7 0.5 0.6];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.7 0.5 0.6]);
    end    
    
    
    plot([mean(Slatency.Delay_2.Correct{iFile}),mean(Slatency.Delay_4.Correct{iFile}),mean(Slatency.Delay_6.Correct{iFile}),mean(Slatency.Delay_8.Correct{iFile})],...
         [mean(Rlatency.Delay_2.Correct{iFile}),mean(Rlatency.Delay_4.Correct{iFile}),mean(Rlatency.Delay_6.Correct{iFile}),mean(Rlatency.Delay_8.Correct{iFile})],':k')
    
    xlabel('Sample latency (s)')
    if iFile==1
        ylabel('Choice Latency (s)')
    end
    axis square
    axis([0 10 0 5])
        fname=strtok(fileList(iFile).name,'_');
        title(fname)

end
    legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','northoutside','Orientation','horizontal'); legend boxoff

    
    %% bars of mean latencies - MEDIUM
figure; 
subplot(1,3,1);hold on
    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=[nanmean(Slatency.Delay_2.Error{iFile}),...
            nanmean(Slatency.Delay_4.Error{iFile}),...
            nanmean(Slatency.Delay_6.Error{iFile}),...
            nanmean(Slatency.Delay_8.Error{iFile})];
    end
    bar(mean(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar(mean(x,2),nansem(x,2),'.r','Marker','none','LineWidth',1.5)

    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=[nanmean(Slatency.Delay_2.Correct{iFile}),...
                    nanmean(Slatency.Delay_4.Correct{iFile}),...
                    nanmean(Slatency.Delay_6.Correct{iFile}),...
                    nanmean(Slatency.Delay_8.Correct{iFile})];
    end
    bar(nanmean(x,2),'FaceColor','w','EdgeColor','k','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.k','Marker','none','LineWidth',1.5)

    set(gca,'XTick',[1 2 3 4],'XTickLabel',Delays__)
    xlabel('Delay duration')
    ylabel('Sample Latency (s)')
    ylim([0 10])

subplot(1,3,2);hold on
    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=[nanmean(NPlatency.Delay_2.Error{iFile}),...
                    nanmean(NPlatency.Delay_4.Error{iFile}),...
                    nanmean(NPlatency.Delay_6.Error{iFile}),...
                    nanmean(NPlatency.Delay_8.Error{iFile})];
    end
    bar(nanmean(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.r','Marker','none','LineWidth',1.5)

    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=[nanmean(NPlatency.Delay_2.Correct{iFile}),...
                    nanmean(NPlatency.Delay_4.Correct{iFile}),...
                    nanmean(NPlatency.Delay_6.Correct{iFile}),...
                    nanmean(NPlatency.Delay_8.Correct{iFile})];
    end
    bar(nanmean(x,2),'FaceColor','w','EdgeColor','k','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.k','Marker','none','LineWidth',1.5)
    
    set(gca,'XTIck',[1 2 3 4],'XTickLabel',Delays__)
    xlabel('Delay duration')
    ylabel('Nosepoke Latency (s)')
    ylim([0 10])

subplot(1,3,3);hold on
    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=[nanmean(Rlatency.Delay_2.Error{iFile}),...
                    nanmean(Rlatency.Delay_4.Error{iFile}),...
                    nanmean(Rlatency.Delay_6.Error{iFile}),...
                    nanmean(Rlatency.Delay_8.Error{iFile})];
    end
    bar(mean(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.r','Marker','none','LineWidth',1.5)
    
    x=[];
    for iFile =1:length(fileList)
        x(:,iFile)=[nanmean(Rlatency.Delay_2.Correct{iFile}),...
                    nanmean(Rlatency.Delay_4.Correct{iFile}),...
                    nanmean(Rlatency.Delay_6.Correct{iFile}),...
                    nanmean(Rlatency.Delay_8.Correct{iFile})];
    end
    bar(nanmean(x,2),'FaceColor','none','EdgeColor','k','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.k','Marker','none','LineWidth',1.5)
    
    set(gca,'XTIck',[1 2 3 4],'XTickLabel',Delays__)
    xlabel('Delay duration')
    ylabel('Lever-press Latency (s)')
    ylim([0 10])    
    
    
    
%% %%%%%% PREAMBLE %%%%%%
clear
Target = 'LONG';
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
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
else
    pat = '/Volumes/HDD2/DNMTP/raw';
end    
cd(pat)
fileList=dir(sprintf('allTimestamps%s*%s*Events.mat',filesep,Target));
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
%% Behaviour - Sample, Nosepoke and Response Latencies (time from end-of-delay signal to lever press)
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s%sallTimestamps%s%s_Events.mat',pat,filesep,filesep,fname));
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

    save(fullfile(pat, 'allTimestamps','Behaviour analysis','DayByDay_LONG.mat'),'Slatency','NPlatency','Rlatency');
%% Plot latency distribution - correct Samples

bins = 0:0.5:20;
figure; hold on
x=histc(cell2mat(Slatency.Short.Correct'),bins);
area(bins,x./sum(x)*100,'Facecolor',[0.1 0.1 0.7],'FaceAlpha',0.3,'LineStyle','none')
y=histc(cell2mat(Slatency.Medium.Correct'),bins);
area(bins,y./sum(y)*100,'Facecolor',[0.2 0.7 0.2],'FaceAlpha',0.4,'LineStyle','none')
z=histc(cell2mat(Slatency.Long.Correct'),bins);
area(bins,z./sum(z)*100,'Facecolor',[0.7 0.2 0.2],'FaceAlpha',0.5,'LineStyle','none')

plot([0 0],[0 6],':k')
plot([20 20],[0 6],':k')
text(0,6.5,'Cue light','Rotation',90,'HorizontalAlignment','left')
text(20,6.5,'Time out','Rotation',90,'HorizontalAlignment','left')
axis([-1 25 0 20])
[bins(x==max(x)),bins(y==max(y)),bins(z==max(z))];
[median(cell2mat(Rlatency.Short.Correct')),median(cell2mat(Rlatency.Medium.Correct')),median(cell2mat(Rlatency.Long.Correct'));]
legend(cellfun(@(x) [x ' Trials'], Delays__,'UniformOutput',false)); legend boxoff

ylabel('Percentage of trials')
xlabel('Sample latency (s)')
[kstest(cell2mat(Slatency.Short.Correct'));kstest(cell2mat(Slatency.Medium.Correct'));kstest(cell2mat(Slatency.Long.Correct'))]

[p,tbl,stats] = kruskalwallis([cell2mat(Slatency.Short.Correct');cell2mat(Slatency.Medium.Correct');cell2mat(Slatency.Long.Correct')],...
              [ones(length(cell2mat(Slatency.Short.Correct')),1);2*ones(length(cell2mat(Slatency.Medium.Correct')),1);3*ones(length(cell2mat(Slatency.Long.Correct')),1)])
[c,~,~,gnames] = multcompare(stats);
%% Plot latency distribution - correct Nosepokes

bins = 0:0.5:20;
figure; hold on
x=histc(cell2mat(NPlatency.Short.Correct'),bins);
area(bins,x./sum(x)*100,'Facecolor',[0.1 0.1 0.7],'FaceAlpha',0.3,'LineStyle','none')
y=histc(cell2mat(NPlatency.Medium.Correct'),bins);
area(bins,y./sum(y)*100,'Facecolor',[0.2 0.7 0.2],'FaceAlpha',0.4,'LineStyle','none')
z=histc(cell2mat(NPlatency.Long.Correct'),bins);
area(bins,z./sum(z)*100,'Facecolor',[0.7 0.2 0.2],'FaceAlpha',0.5,'LineStyle','none')

plot([0 0],[0 6],':k')
plot([20 20],[0 6],':k')
text(0,6.5,'Cue light','Rotation',90,'HorizontalAlignment','left')
text(20,6.5,'Time out','Rotation',90,'HorizontalAlignment','left')
axis([-1 25 0 20])
[bins(x==max(x)),bins(y==max(y)),bins(z==max(z))];
[median(cell2mat(NPlatency.Short.Correct')),median(cell2mat(NPlatency.Medium.Correct')),median(cell2mat(NPlatency.Long.Correct'));]
legend(cellfun(@(x) [x ' Trials'], Delays__,'UniformOutput',false)); legend boxoff

ylabel('Percentage of trials')
xlabel('Nose-poke latency (s)')

[kstest(cell2mat(NPlatency.Short.Correct'));kstest(cell2mat(NPlatency.Medium.Correct'));kstest(cell2mat(NPlatency.Long.Correct'))]
[p,tbl,stats] =kruskalwallis([cell2mat(NPlatency.Short.Correct');cell2mat(NPlatency.Medium.Correct');cell2mat(NPlatency.Long.Correct')],...
                             [ones(length(cell2mat(NPlatency.Short.Correct')),1);2*ones(length(cell2mat(NPlatency.Medium.Correct')),1);3*ones(length(cell2mat(NPlatency.Long.Correct')),1)])
                         
                         
[c,~,~,gnames] = multcompare(stats);
%% Plot latency distribution - correct Choices

bins = 0:0.5:20;
figure; hold on
x=histc(cell2mat(Rlatency.Short.Correct'),bins);
area(bins,x./sum(x)*100,'Facecolor',[0.1 0.1 0.7],'FaceAlpha',0.3,'LineStyle','none')
y=histc(cell2mat(Rlatency.Medium.Correct'),bins);
area(bins,y./sum(y)*100,'Facecolor',[0.2 0.7 0.2],'FaceAlpha',0.4,'LineStyle','none')
z=histc(cell2mat(Rlatency.Long.Correct'),bins);
area(bins,z./sum(z)*100,'Facecolor',[0.7 0.2 0.2],'FaceAlpha',0.5,'LineStyle','none')

plot([0 0],[0 6],':k')
plot([20 20],[0 6],':k')
text(0,6.5,'Cue light','Rotation',90,'HorizontalAlignment','left')
text(20,6.5,'Time out','Rotation',90,'HorizontalAlignment','left')
axis([-1 25 0 20])
[bins(x==max(x)),bins(y==max(y)),bins(z==max(z))];
[median(cell2mat(Rlatency.Short.Correct')),median(cell2mat(Rlatency.Medium.Correct')),median(cell2mat(Rlatency.Long.Correct'));]
legend(cellfun(@(x) [x ' Trials'], Delays__,'UniformOutput',false)); legend boxoff

ylabel('Percentage of trials')
xlabel('Nose-poke latency (s)')
[kstest(cell2mat(Rlatency.Short.Correct'));kstest(cell2mat(Rlatency.Medium.Correct'));kstest(cell2mat(Rlatency.Long.Correct'))]

[p,tbl,stats] =kruskalwallis([cell2mat(Rlatency.Short.Correct');cell2mat(Rlatency.Medium.Correct');cell2mat(Rlatency.Long.Correct')],...
                             [ones(length(cell2mat(Rlatency.Short.Correct')),1);2*ones(length(cell2mat(Rlatency.Medium.Correct')),1);3*ones(length(cell2mat(Rlatency.Long.Correct')),1)])
[c,~,~,gnames] = multcompare(stats);

%% Plot Nosepoke and Choice Response latency distribution

bins = 0:0.2:10;
figure; hold on
x=histc(cell2mat(Rlatency.Short.Correct'),bins);
area(bins,x./sum(x)*100,'Facecolor',[0.1 0.1 0.7],'FaceAlpha',0.3,'LineStyle','none')
y=histc(cell2mat(Rlatency.Medium.Correct'),bins);
area(bins,y./sum(y)*100,'Facecolor',[0.2 0.7 0.2],'FaceAlpha',0.4,'LineStyle','none')
z=histc(cell2mat(Rlatency.Long.Correct'),bins);
area(bins,z./sum(z)*100,'Facecolor',[0.7 0.2 0.2],'FaceAlpha',0.5,'LineStyle','none')

bins = 0:0.2:30;
x=histc(cell2mat(NPlatency.Short.Correct'),bins);
plot(bins,-x./sum(x)*100,'color',[0.1 0.1 0.7 0.3],'LineStyle','-')

x=histc(cell2mat(NPlatency.Medium.Correct'),bins);
plot(bins,-x./sum(x)*100,'color',[0.2 0.7 0.2 0.4],'LineStyle','-')

x=histc(cell2mat(NPlatency.Long.Correct'),bins);
plot(bins,-x./sum(x)*100,'color',[0.7 0.2 0.2 0.5],'LineStyle','-')


plot([0 0],[0 6],':k')
plot([10 10],[0 6],':k')
text(0,6.5,'Tone','Rotation',90,'HorizontalAlignment','left')
text(10,6.5,'Time out','Rotation',90,'HorizontalAlignment','left')
axis([-1 25 -15 20])
[bins(x==max(x)),bins(y==max(y)),bins(z==max(z))]
[median(cell2mat(Rlatency.Short.Correct')),median(cell2mat(Rlatency.Medium.Correct')),median(cell2mat(Rlatency.Long.Correct'));]

legend([cellfun(@(x) [x ' Choices'], Delays__,'UniformOutput',false),cellfun(@(x) [x ' Nose-pokes'], Delays__,'UniformOutput',false)]); legend boxoff

ylabel('Percentage of trials')
xlabel('Choice latency (s)')
% kruskalwallis([cell2mat(NPlatency.Short.Correct');cell2mat(NPlatency.Medium.Correct');cell2mat(NPlatency.Long.Correct')],...
%               [ones(length(cell2mat(NPlatency.Short.Correct')),1);2*ones(length(cell2mat(NPlatency.Medium.Correct')),1);3*ones(length(cell2mat(NPlatency.Long.Correct')),1)])
% 
% [p,tbl,stats] = kruskalwallis([cell2mat(Rlatency.Short.Correct');cell2mat(Rlatency.Medium.Correct');cell2mat(Rlatency.Long.Correct')],...
%               [ones(length(cell2mat(Rlatency.Short.Correct')),1);2*ones(length(cell2mat(Rlatency.Medium.Correct')),1);3*ones(length(cell2mat(Rlatency.Long.Correct')),1)])
% c = multcompare(stats);
%% Behaviour - Plot all sorted correct trials
Trialraster.Short = [];
Trialraster.Medium= [];
Trialraster.Long = [];

for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s%sallTimestamps%s%s_Events.mat',pat,filesep,filesep,fname));
    
    t_ = [t.Short.DelayEnd_LeftCorrect,t.Short.NosePoke_LeftCorrect,t.Short.ChoicePress_LeftCorrect;...
    t.Short.DelayEnd_RightCorrect,t.Short.NosePoke_RightCorrect,t.Short.ChoicePress_RightCorrect;]*1e-6;
    t_ = bsxfun(@minus,t_,t_(:,1));
    Trialraster.Short = [Trialraster.Short ;t_];
    
    t_ = [t.Medium.DelayEnd_LeftCorrect,t.Medium.NosePoke_LeftCorrect,t.Medium.ChoicePress_LeftCorrect;...
    t.Medium.DelayEnd_RightCorrect,t.Medium.NosePoke_RightCorrect,t.Medium.ChoicePress_RightCorrect;]*1e-6;
    t_ = bsxfun(@minus,t_,t_(:,1));
    Trialraster.Medium = [Trialraster.Medium ;t_];
    
    t_ = [t.Long.DelayEnd_LeftCorrect,t.Long.NosePoke_LeftCorrect,t.Long.ChoicePress_LeftCorrect;...
    t.Long.DelayEnd_RightCorrect,t.Long.NosePoke_RightCorrect,t.Long.ChoicePress_RightCorrect;]*1e-6;
    t_ = bsxfun(@minus,t_,t_(:,1));
    Trialraster.Long = [Trialraster.Long ;t_];    
end
Trialraster.Short  = sortrows(Trialraster.Short,2);
Trialraster.Medium = sortrows(Trialraster.Medium,2);
Trialraster.Long   = sortrows(Trialraster.Long,2);

Ns = [size(Trialraster.Short,1),size(Trialraster.Medium,1),size(Trialraster.Long,1)];
X = [Trialraster.Short;Trialraster.Medium;Trialraster.Long];
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

axis([-5 30 -50 1300])

plot([28 28],[0 Ns(1)],'color',[0.1 0.1 0.7 0.3],'LineWidth',2)
plot([28 28],[Ns(1) Ns(1)+Ns(2) ],'color',[0.2 0.7 0.2 0.4],'LineWidth',2)
plot([28 28],[Ns(1)+Ns(2) sum(Ns)],'color',[0.7 0.2 0.2 0.5],'LineWidth',2)
text(29,Ns(1)/2,'4s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.1 0.1 0.7 0.3])
text(29,Ns(1) + Ns(2)/2,'8s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.2 0.7 0.2 0.4])
text(29,Ns(1)+Ns(2) + Ns(3)/2,'16s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.7 0.2 0.2 0.5])
legend('Nose-poke','Choice Press','Location','southoutside','Orientation','horizontal'); legend boxoff
%% Behaviour - Plot all sorted error trials
TrialrasterError.Short = [];
TrialrasterError.Medium= [];
TrialrasterError.Long = [];

for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s%sallTimestamps%s%s_Events.mat',pat,filesep,filesep,fname));
    
    t_ =     [t.Short.DelayEnd_LeftError,t.Short.DelayEnd_RightError; ... 
              t.Short.NosePoke_LeftError,t.Short.NosePoke_RightError; ...
              t.Short.ChoicePress_LeftError,t.Short.ChoicePress_RightError]'*1e-6;
    t_ = bsxfun(@minus,t_,t_(:,1));
    TrialrasterError.Short = [TrialrasterError.Short ;t_];    
    
    t_ =     [t.Medium.DelayEnd_LeftError,t.Medium.DelayEnd_RightError; ... 
              t.Medium.NosePoke_LeftError,t.Medium.NosePoke_RightError; ...
              t.Medium.ChoicePress_LeftError,t.Medium.ChoicePress_RightError]'*1e-6;
    t_ = bsxfun(@minus,t_,t_(:,1));
    TrialrasterError.Medium = [TrialrasterError.Medium ;t_];
    
    t_ =     [t.Long.DelayEnd_LeftError,t.Long.DelayEnd_RightError; ... 
              t.Long.NosePoke_LeftError,t.Long.NosePoke_RightError; ...
              t.Long.ChoicePress_LeftError,t.Long.ChoicePress_RightError]'*1e-6;
    t_ = bsxfun(@minus,t_,t_(:,1));
    TrialrasterError.Long = [TrialrasterError.Long ;t_];    
end
TrialrasterError.Short  = sortrows(TrialrasterError.Short,3);
TrialrasterError.Medium = sortrows(TrialrasterError.Medium,3);
TrialrasterError.Long   = sortrows(TrialrasterError.Long,3);

Ns = [size(TrialrasterError.Short,1),size(TrialrasterError.Medium,1),size(TrialrasterError.Long,1)];
X = [TrialrasterError.Short;TrialrasterError.Medium;TrialrasterError.Long];
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

axis([-5 30 -50 1300])

plot([28 28],[0 Ns(1)],'color',[0.1 0.1 0.7 0.3],'LineWidth',2)
plot([28 28],[Ns(1) Ns(1)+Ns(2) ],'color',[0.2 0.7 0.2 0.4],'LineWidth',2)
plot([28 28],[Ns(1)+Ns(2) sum(Ns)],'color',[0.7 0.2 0.2 0.5],'LineWidth',2)
text(29,Ns(1)/2,'4s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.1 0.1 0.7 0.3])
text(29,Ns(1) + Ns(2)/2,'8s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.2 0.7 0.2 0.4])
text(29,Ns(1)+Ns(2) + Ns(3)/2,'16s Delay trials','Rotation',90,'HorizontalAlignment','center','Color',[0.7 0.2 0.2 0.5])
legend('Nose-poke','Choice Press','Location','southoutside','Orientation','horizontal'); legend boxoff
%% Cumulative distributions of latencies - correct trials
figure
x1= diff(Trialraster.Short,[],2);
x2= diff(Trialraster.Medium,[],2);
x3= diff(Trialraster.Long,[],2);
subplot(1,2,1); hold on
    bins =0:0.2:15;
    x = hist(x1(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.1 0.1 0.7],'LineWidth',1.5)
    x = hist(x2(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.2 0.7 0.2],'LineWidth',1.5)
    x = hist(x3(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.7 0.2 0.2],'LineWidth',1.5)
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
x1= diff(TrialrasterError.Short,[],2);
x2= diff(TrialrasterError.Medium,[],2);
x3= diff(TrialrasterError.Long,[],2);
subplot(1,2,1); hold on
    bins =0:0.2:15;
    x = hist(x1(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.1 0.1 0.7],'LineWidth',1.5)
    x = hist(x2(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.2 0.7 0.2],'LineWidth',1.5)
    x = hist(x3(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.7 0.2 0.2],'LineWidth',1.5)
    axis square
    ylabel('Cumulative fraction of trials')
    xlabel('Tone-Nosepoke latency (s)')
    legend([cellfun(@(x) [x ' Trials'], Delays__,'UniformOutput',false),cellfun(@(x) [x ' Nose-pokes'], Delays__,'UniformOutput',false)],'Location','southeast'); legend boxoff

subplot(1,2,2); hold on
    bins =0:0.2:10;
    x = hist(x1(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.1 0.1 0.7],'LineWidth',1.5)
    x = hist(x2(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.2 0.7 0.2],'LineWidth',1.5)
    x = hist(x3(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color',[0.7 0.2 0.2],'LineWidth',1.5)
    axis square
    xlabel('Nosepoke-Choice latency (s)')
legend([cellfun(@(x) [x ' Trials'], Delays__,'UniformOutput',false),cellfun(@(x) [x ' Nose-pokes'], Delays__,'UniformOutput',false)],'Location','southeast'); legend boxoff

x1= diff(Trialraster.Short,[],2);
x2= diff(Trialraster.Medium,[],2);
x3= diff(Trialraster.Long,[],2);
x1e= diff(TrialrasterError.Short,[],2);
x2e= diff(TrialrasterError.Medium,[],2);
x3e= diff(TrialrasterError.Long,[],2);
    [h(1,1),p(1,1)] = kstest2(x1(:,1),x1e(:,1));
    [h(1,2),p(1,2)] = kstest2(x2(:,1),x2e(:,1));
    [h(1,3),p(1,3)] = kstest2(x3(:,1),x3e(:,1));

    [h(2,1),p(2,1)] = kstest2(x1(:,2),x1e(:,2));
    [h(2,2),p(2,2)] = kstest2(x2(:,2),x2e(:,2));
    [h(2,3),p(2,3)] = kstest2(x3(:,2),x3e(:,2));

figure('color','w','name','Tone-Nosepoke latency');

subplot(1,3,1); hold on
    bins =0:0.2:20;
    x = hist(x1(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x1e(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
%     title('Tone - Nosepoke Delay')
    text(max(bins),0.1,'4s delay trials','HorizontalAlignment','right','color',[0.1 0.1 0.7])
    axis square
    xlabel('Latency (s)')
    ylabel('Cumulative fraction of trials')
    
subplot(1,3,2); hold on
    x = hist(x2(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x2e(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
    text(max(bins),0.1,'8s delay trials','HorizontalAlignment','right','color',[0.2 0.7 0.2])
    xlabel('Latency (s)')
    axis square

subplot(1,3,3); hold on
    x = hist(x3(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x3e(:,1),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)      
    text(max(bins),0.1,'16s delay trials','HorizontalAlignment','right','color',[0.7 0.2 0.2])
    xlabel('Latency (s)')
    axis square


    
figure('color','w','name','Nosepoke-Choice latency');

subplot(1,3,1); hold on
    bins =0:0.2:10;
    x = hist(x1(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x1e(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
%     title('Nosepoke-Choice Delay')
    text(max(bins),0.1,'4s delay trials','HorizontalAlignment','right','color',[0.1 0.1 0.7])
    axis square
    ylabel('Cumulative fraction of trials')
    xlabel('Latency (s)')
    [h,p] = kstest2(x1(:,2),x1e(:,2))
    
subplot(1,3,2); hold on
    x = hist(x2(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x2e(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)
    text(max(bins),0.1,'8s delay trials','HorizontalAlignment','right','color',[0.2 0.7 0.2])
    axis square
    xlabel('Latency (s)')

subplot(1,3,3); hold on
    x = hist(x3(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','k','LineWidth',1.5)
    x = hist(x3e(:,2),bins);x=cumsum(x)./sum(x);
    plot(bins,x,'Color','r','LineWidth',1.5)      
    text(max(bins),0.1,'16s delay trials','HorizontalAlignment','right','color',[0.7 0.2 0.2])
    xlabel('Latency (s)')
    axis square  
%% Sample vs Choice Latency - points
figure;
for iFile =1:length(fileList)
    subplot(1,3,1);hold on
    scatter(Slatency.Short.Correct{iFile},Rlatency.Short.Correct{iFile},10,[0.1 0.1 0.7],'filled')
    try scatter(Slatency.Short.Error{iFile},Rlatency.Short.Error{iFile},10,[0.1 0.1 0.7]); end
    axis square
    axis([0 30 0 30])
    ylabel('Choice press latency (s)')
    
    subplot(1,3,2);hold on
    scatter(Slatency.Medium.Correct{iFile},Rlatency.Medium.Correct{iFile},10,[0.2 0.7 0.2],'filled')
    try scatter(Slatency.Medium.Error{iFile},Rlatency.Medium.Error{iFile},10,[0.2 0.7 0.2]); end
    axis square
    axis([0 30 0 30])
    xlabel('Sample press latency (s)')
    
    subplot(1,3,3);hold on
    scatter(Slatency.Long.Correct{iFile},Rlatency.Long.Correct{iFile},10,[0.7 0.2 0.2],'filled')
    try scatter(Slatency.Long.Error{iFile},Rlatency.Long.Error{iFile},10,[0.7 0.2 0.2]); end
    axis square
    axis([0 30 0 30])
    
end
%% Sample vs Choice Latency - errorbars 
figure; hold on
for iFile =1:length(fileList)
    H = errorbarxy(mean(Slatency.Short.Correct{iFile}),mean(Rlatency.Short.Correct{iFile}),...
               nansem(Slatency.Short.Correct{iFile}),nansem(Rlatency.Short.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.1 0.1 0.7];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.1 0.1 0.7]);
    end
    
    H = errorbarxy(mean(Slatency.Medium.Correct{iFile}),mean(Rlatency.Medium.Correct{iFile}),...
               nansem(Slatency.Medium.Correct{iFile}),nansem(Rlatency.Medium.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.2 0.7 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.2 0.7 0.2]);
    end    

    H = errorbarxy(mean(Slatency.Long.Correct{iFile}),mean(Rlatency.Long.Correct{iFile}),...
               nansem(Slatency.Long.Correct{iFile}),nansem(Rlatency.Long.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.7 0.2 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.7 0.2 0.2]);
    end    
    
    xlabel('Sample latency (s)')
    ylabel('Choice Latency (s)')
    axis square
    axis([0 15 0 10])
    legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','south'); legend boxoff
end
%% Sample vs Choice Latency - linked errorbars by animal
files = setdiff(1:length(fileList),7);
figure('color','w'); hold on
for iFile = 1:length(files)
%     subplot(ceil(sqrt(length(fileList))),ceil(sqrt(length(fileList))),iFile);hold on
%     subaxis(ceil(sqrt(length(fileList))),ceil(sqrt(length(fileList))),iFile,'Spacing',0.09);hold on
    subaxis(6,2,iFile,'SpacingVert',0.05,'MarginBottom',0.1,'Padding',0.001);hold on
%     subplot(1,length(fileList),iFile);hold on
    H = errorbarxy(mean(Slatency.Short.Correct{files(iFile)}),mean(Rlatency.Short.Correct{files(iFile)}),...
               nansem(Slatency.Short.Correct{files(iFile)}),nansem(Rlatency.Short.Correct{files(iFile)}),{'k','k','k'});
    H.hMain.Color = [0.1 0.1 0.7];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.1 0.1 0.7]);
    end
   
    
    H = errorbarxy(mean(Slatency.Medium.Correct{files(iFile)}),mean(Rlatency.Medium.Correct{files(iFile)}),...
               nansem(Slatency.Medium.Correct{files(iFile)}),nansem(Rlatency.Medium.Correct{files(iFile)}),{'k','k','k'});
    H.hMain.Color = [0.2 0.7 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.2 0.7 0.2]);
    end    
  

    H = errorbarxy(mean(Slatency.Long.Correct{files(iFile)}),mean(Rlatency.Long.Correct{files(iFile)}),...
               nansem(Slatency.Long.Correct{files(iFile)}),nansem(Rlatency.Long.Correct{files(iFile)}),{'k','k','k'});
    H.hMain.Color = [0.7 0.2 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.7 0.2 0.2]);
    end
  
    plot([mean(Slatency.Short.Correct{files(iFile)}),mean(Slatency.Medium.Correct{files(iFile)}),mean(Slatency.Long.Correct{files(iFile)})],...
         [mean(Rlatency.Short.Correct{files(iFile)}),mean(Rlatency.Medium.Correct{files(iFile)}),mean(Rlatency.Long.Correct{files(iFile)})],':k')
  
    axis square
    axis([0 12 0 5])
    fname=strtok(fileList(files(iFile)).name,'_');
    title(fname,'FontSize',8)
%     text(0.5,0.5,fname,'FontSize',8)
end

% xlabel('Sample latency (s)')
% ylabel('Choice Latency (s)')
% legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','eastoutside','Orientation','vertical'); legend boxoff
%% Sample vs Choice Latency - linked errorbars by animal - median
files = setdiff(1:length(fileList),7);
figure('color','w'); hold on
for iFile = 1:length(files)
%     subplot(ceil(sqrt(length(fileList))),ceil(sqrt(length(fileList))),iFile);hold on
%     subaxis(ceil(sqrt(length(fileList))),ceil(sqrt(length(fileList))),iFile,'Spacing',0.09);hold on
    subaxis(6,2,iFile,'SpacingVert',0.05,'MarginBottom',0.1,'Padding',0.001);hold on
%     subplot(1,length(fileList),iFile);hold on
    H = errorbarxy(median(Slatency.Short.Correct{files(iFile)}),median(Rlatency.Short.Correct{files(iFile)}),...
               iqr(Slatency.Short.Correct{files(iFile)}),iqr(Rlatency.Short.Correct{files(iFile)}),{'k','k','k'});
    H.hMain.Color = [0.1 0.1 0.7];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.1 0.1 0.7]);
    end
   
    
    H = errorbarxy(mean(Slatency.Medium.Correct{files(iFile)}),mean(Rlatency.Medium.Correct{files(iFile)}),...
               nansem(Slatency.Medium.Correct{files(iFile)}),nansem(Rlatency.Medium.Correct{files(iFile)}),{'k','k','k'});
    H.hMain.Color = [0.2 0.7 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.2 0.7 0.2]);
    end    
  

    H = errorbarxy(mean(Slatency.Long.Correct{files(iFile)}),mean(Rlatency.Long.Correct{files(iFile)}),...
               nansem(Slatency.Long.Correct{files(iFile)}),nansem(Rlatency.Long.Correct{files(iFile)}),{'k','k','k'});
    H.hMain.Color = [0.7 0.2 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.7 0.2 0.2]);
    end
  
    plot([mean(Slatency.Short.Correct{files(iFile)}),mean(Slatency.Medium.Correct{files(iFile)}),mean(Slatency.Long.Correct{files(iFile)})],...
         [mean(Rlatency.Short.Correct{files(iFile)}),mean(Rlatency.Medium.Correct{files(iFile)}),mean(Rlatency.Long.Correct{files(iFile)})],':k')
  
    axis square
    axis([0 12 0 5])
    fname=strtok(fileList(files(iFile)).name,'_');
    title(fname,'FontSize',8)
%     text(0.5,0.5,fname,'FontSize',8)
end

% xlabel('Sample latency (s)')
% ylabel('Choice Latency (s)')
% legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','eastoutside','Orientation','vertical'); legend boxoff
%% Nosepoke vs Choice Latency - linked errorbars by animal
figure('color','w'); hold on
for iFile =1:length(fileList)
    subplot(ceil(sqrt(length(fileList))),ceil(sqrt(length(fileList))),iFile);hold on
    subaxis(ceil(sqrt(length(fileList))),ceil(sqrt(length(fileList))),iFile,'Spacing',0.09);hold on
%     subplot(1,length(fileList),iFile);hold on
    H = errorbarxy(mean(NPlatency.Short.Correct{iFile}),mean(Rlatency.Short.Correct{iFile}),...
               nansem(NPlatency.Short.Correct{iFile}),nansem(Rlatency.Short.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.1 0.1 0.7];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.1 0.1 0.7]);
    end
   
    
    H = errorbarxy(mean(NPlatency.Medium.Correct{iFile}),mean(Rlatency.Medium.Correct{iFile}),...
               nansem(NPlatency.Medium.Correct{iFile}),nansem(Rlatency.Medium.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.2 0.7 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.2 0.7 0.2]);
    end    
  

    H = errorbarxy(mean(NPlatency.Long.Correct{iFile}),mean(Rlatency.Long.Correct{iFile}),...
               nansem(NPlatency.Long.Correct{iFile}),nansem(Rlatency.Long.Correct{iFile}),{'k','k','k'});
    H.hMain.Color = [0.7 0.2 0.2];
    for i=1:6
        set(H.hErrorbar(i),'Color',[0.7 0.2 0.2]);
    end
  
    plot([mean(NPlatency.Short.Correct{iFile}),mean(NPlatency.Medium.Correct{iFile}),mean(NPlatency.Long.Correct{iFile})],...
         [mean(Rlatency.Short.Correct{iFile}),mean(Rlatency.Medium.Correct{iFile}),mean(Rlatency.Long.Correct{iFile})],':k')
  
    axis square
    axis([0 5 0 5])
    fname=strtok(fileList(iFile).name,'_');
    title(fname)
end

xlabel('Nosepoke latency (s)')
ylabel('Choice Latency (s)')
legend(cellfun(@(x) [x ' trials'], Delays__,'UniformOutput',false),'Location','eastoutside','Orientation','vertical'); legend boxoff    %% bars of mean latencies - LONG

%% bar plots 
figure; 
subplot(3,1,1);hold on
    x=[];
    for iFile =1:length(fileList)
%         x(:,iFile)=[nanmean(Slatency.Short.Error{iFile}),...
%             nanmean(Slatency.Medium.Error{iFile}),...
%             nanmean(Slatency.Long.Error{iFile})];
        x(:,iFile)=[nanmedian(Slatency.Short.Error{iFile}),...
            nanmedian(Slatency.Medium.Error{iFile}),...
            nanmedian(Slatency.Long.Error{iFile})];
    end
    bar(mean(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar(mean(x,2),nansem(x,2),'.r','Marker','none','LineWidth',1.5)
%     bar(median(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
%     errorbar(median(x,2),iqr(x,2),'.r','Marker','none','LineWidth',1.5)

    x=[];
    for iFile =1:length(fileList)
%         x(:,iFile)=[nanmean(Slatency.Short.Correct{iFile}),...
%                     nanmean(Slatency.Medium.Correct{iFile}),...
%                     nanmean(Slatency.Long.Correct{iFile})];
        x(:,iFile)=[nanmedian(Slatency.Short.Correct{iFile}),...
                    nanmedian(Slatency.Medium.Correct{iFile}),...
                    nanmedian(Slatency.Long.Correct{iFile})];
    end
    bar(nanmean(x,2),'FaceColor','w','EdgeColor','k','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.k','Marker','none','LineWidth',1.5) 
%     bar(median(x,2),'FaceColor','w','EdgeColor','k','LineWidth',1.5)
%     errorbar(median(x,2),iqr(x,2),'.k','Marker','none','LineWidth',1.5)

    set(gca,'XTick',[1 2 3],'XTickLabel',Delays__)
%     xlabel('Delay duration')
%     ylabel('Latency (s)')
    title('Sample')
    ylim([0 10])

subplot(3,1,2);hold on
    x=[];
    for iFile =1:length(fileList)
%         x(:,iFile)=[nanmean(NPlatency.Short.Error{iFile}),...
%                     nanmean(NPlatency.Medium.Error{iFile}),...
%                     nanmean(NPlatency.Long.Error{iFile})];        
        x(:,iFile)=[nanmedian(NPlatency.Short.Error{iFile}),...
                    nanmedian(NPlatency.Medium.Error{iFile}),...
                    nanmedian(NPlatency.Long.Error{iFile})];
    end
    bar(nanmean(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.r','Marker','none','LineWidth',1.5)
%     bar(median(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
%     errorbar(median(x,2),iqr(x,2),'.r','Marker','none','LineWidth',1.5)

    x=[];
    for iFile =1:length(fileList)
%         x(:,iFile)=[nanmean(NPlatency.Short.Correct{iFile}),...
%                     nanmean(NPlatency.Medium.Correct{iFile}),...
%                     nanmean(NPlatency.Long.Correct{iFile})];        
        x(:,iFile)=[nanmedian(NPlatency.Short.Correct{iFile}),...
                    nanmedian(NPlatency.Medium.Correct{iFile}),...
                    nanmedian(NPlatency.Long.Correct{iFile})];
    end
    bar(nanmean(x,2),'FaceColor','w','EdgeColor','k','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.k','Marker','none','LineWidth',1.5)    
%     bar(median(x,2),'FaceColor','w','EdgeColor','k','LineWidth',1.5)
%     errorbar(median(x,2),iqr(x,2),'.k','Marker','none','LineWidth',1.5)
    
    set(gca,'XTick',[1 2 3],'XTickLabel',Delays__)
%     xlabel('Delay duration')
%     ylabel('Nosepoke Latency (s)')
    title('Nosepoke')
    ylim([0 8])
    ylabel(' Latency (Median +/- IQR, s)')


subplot(3,1,3);hold on
    x=[];
    for iFile =1:length(fileList)
%         x(:,iFile)=[nanmean(Rlatency.Short.Error{iFile}),...
%                     nanmean(Rlatency.Medium.Error{iFile}),...
%                     nanmean(Rlatency.Long.Error{iFile})];        
                x(:,iFile)=[nanmedian(Rlatency.Short.Error{iFile}),...
                    nanmedian(Rlatency.Medium.Error{iFile}),...
                    nanmedian(Rlatency.Long.Error{iFile})];
    end
    bar(mean(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.r','Marker','none','LineWidth',1.5)
%     bar(median(x,2),'FaceColor','w','EdgeColor','r','LineWidth',1.5)
%     errorbar(median(x,2),iqr(x,2),'.r','Marker','none','LineWidth',1.5)    
    x=[];
    for iFile =1:length(fileList)
%         x(:,iFile)=[nanmean(Rlatency.Short.Correct{iFile}),...
%                     nanmean(Rlatency.Medium.Correct{iFile}),...
%                     nanmean(Rlatency.Long.Correct{iFile})];       
        x(:,iFile)=[nanmedian(Rlatency.Short.Correct{iFile}),...
                    nanmedian(Rlatency.Medium.Correct{iFile}),...
                    nanmedian(Rlatency.Long.Correct{iFile})];
    end
    bar(nanmean(x,2),'FaceColor','none','EdgeColor','k','LineWidth',1.5)
    errorbar(nanmean(x,2),nansem(x,2),'.k','Marker','none','LineWidth',1.5)   
%     bar(median(x,2),'FaceColor','none','EdgeColor','k','LineWidth',1.5)
%     errorbar(median(x,2),iqr(x,2),'.k','Marker','none','LineWidth',1.5)
    
    set(gca,'XTIck',[1 2 3],'XTickLabel',Delays__)
    xlabel('Delay duration')
%     ylabel('Lever-press Latency (s)')
    title('Choice')
    ylim([0 8])
%% box plots
    NOTCH = true;
    SYM   = [];
    VERT  = true;
    WHIS  = 0.5;
    FILLIT = true;
    LINEWIDTH = 0.5;
    KOUTLINE  = false;
    FATTER    = true;
    LESSOVERLAP = 0.01;
    MEANSYM   = [];
   grouping = [repmat({'4s'},sum(cellfun(@length,Slatency.Short.Correct)),1);...
               repmat({'8s'},sum(cellfun(@length,Slatency.Medium.Correct)),1);...
               repmat({'16s'},sum(cellfun(@length,Slatency.Long.Correct)),1)];
 
    groupingE = [repmat({'4s'},sum(cellfun(@length,Slatency.Short.Error)),1);...
                 repmat({'8s'},sum(cellfun(@length,Slatency.Medium.Error)),1);...
                 repmat({'16s'},sum(cellfun(@length,Slatency.Long.Error)),1)];
             
    x = cell(3,1); xE = cell(3,1);
    for iFile =1:length(fileList)

        x{1} = [x{1};(Slatency.Short.Correct{iFile})];
        x{2} = [x{2};(Slatency.Medium.Correct{iFile})];
        x{3} = [x{3};(Slatency.Long.Correct{iFile})];
        
        xE{1} = [xE{1};(Slatency.Short.Error{iFile})];
        xE{2} = [xE{2};(Slatency.Medium.Error{iFile})];
        xE{3} = [xE{3};(Slatency.Long.Error{iFile})];

    end
    

    figure; 
    subplot(3,1,1);hold on
    x = cell(3,1); xE = cell(3,1);
    for iFile =1:length(fileList)

        x{1} = [x{1};(Slatency.Short.Correct{iFile})];
        x{2} = [x{2};(Slatency.Medium.Correct{iFile})];
        x{3} = [x{3};(Slatency.Long.Correct{iFile})];
        
        xE{1} = [xE{1};(Slatency.Short.Error{iFile})];
        xE{2} = [xE{2};(Slatency.Medium.Error{iFile})];
        xE{3} = [xE{3};(Slatency.Long.Error{iFile})];

    end
	boxplotCsub(cell2mat(x), grouping ,NOTCH,SYM,VERT,WHIS,'k',FILLIT,LINEWIDTH,KOUTLINE,[1 2],FATTER,LESSOVERLAP,MEANSYM)
	boxplotCsub(cell2mat(xE),groupingE,NOTCH,SYM,VERT,WHIS,'r',FILLIT,LINEWIDTH,KOUTLINE,[2 2],FATTER,LESSOVERLAP,MEANSYM)
    title('Cue - Sample')
    
    if VERT
        ylim([0 20])
        ylabel('')

    else
        xlabel('')
        xlim([0 20])
    end
    
    
    subplot(3,1,2);hold on
    x = cell(3,1); xE = cell(3,1);
    for iFile =1:length(fileList)

        x{1} = [x{1};(NPlatency.Short.Correct{iFile})];
        x{2} = [x{2};(NPlatency.Medium.Correct{iFile})];
        x{3} = [x{3};(NPlatency.Long.Correct{iFile})];
        
        xE{1} = [xE{1};(NPlatency.Short.Error{iFile})];
        xE{2} = [xE{2};(NPlatency.Medium.Error{iFile})];
        xE{3} = [xE{3};(NPlatency.Long.Error{iFile})];

    end
	boxplotCsub(cell2mat(x), grouping ,NOTCH,SYM,VERT,WHIS,'k',FILLIT,LINEWIDTH,KOUTLINE,[1 2],FATTER,LESSOVERLAP,MEANSYM)
	boxplotCsub(cell2mat(xE),groupingE,NOTCH,SYM,VERT,WHIS,'r',FILLIT,LINEWIDTH,KOUTLINE,[2 2],FATTER,LESSOVERLAP,MEANSYM)
    title('Tone - Nosepoke')
    ylabel('Delay Length')
     if VERT
        ylim([0 10])
    else
        xlim([0 10])
        xlabel('')
        ylabel('Latency (s)')
    end
    
    
    subplot(3,1,3);hold on
    x = cell(3,1); xE = cell(3,1);
    for iFile =1:length(fileList)

        x{1} = [x{1};(Rlatency.Short.Correct{iFile})];
        x{2} = [x{2};(Rlatency.Medium.Correct{iFile})];
        x{3} = [x{3};(Rlatency.Long.Correct{iFile})];
        
        xE{1} = [xE{1};(Rlatency.Short.Error{iFile})];
        xE{2} = [xE{2};(Rlatency.Medium.Error{iFile})];
        xE{3} = [xE{3};(Rlatency.Long.Error{iFile})];

    end
	boxplotCsub(cell2mat(x), grouping ,NOTCH,SYM,VERT,WHIS,'k',FILLIT,LINEWIDTH,KOUTLINE,[1 2],FATTER,LESSOVERLAP,MEANSYM)
	boxplotCsub(cell2mat(xE),groupingE,NOTCH,SYM,VERT,WHIS,'r',FILLIT,LINEWIDTH,KOUTLINE,[2 2],FATTER,LESSOVERLAP,MEANSYM)
    title('Nosepoke - Choice')
 if VERT
        ylim([0 6])
        ylabel('')
        xlabel('Delay length')
 else
        xlabel('Latency (s)')
        xlim([0 6])
 end
%% big ANOVA
    [p,tbl,stats] = kruskalwallis([cell2mat(Slatency.Short.Correct');...
                                   cell2mat(Slatency.Short.Error');...
                                   cell2mat(Slatency.Medium.Correct');...
                                   cell2mat(Slatency.Medium.Error');...
                                   cell2mat(Slatency.Long.Correct');...
                                   cell2mat(Slatency.Long.Error')],...
                                  [  ones(length(cell2mat(Slatency.Short.Correct')),1);...
                                   2*ones(length(cell2mat(Slatency.Short.Error')),1);...
                                   3*ones(length(cell2mat(Slatency.Medium.Correct')),1);...
                                   4*ones(length(cell2mat(Slatency.Medium.Error')),1);...
                                   5*ones(length(cell2mat(Slatency.Long.Correct')),1);...
                                   6*ones(length(cell2mat(Slatency.Long.Error')),1)]);
    c = multcompare(stats);
    
    [p,tbl,stats] = kruskalwallis([cell2mat(NPlatency.Short.Correct');...
                                   cell2mat(NPlatency.Short.Error');...
                                   cell2mat(NPlatency.Medium.Correct');...
                                   cell2mat(NPlatency.Medium.Error');...
                                   cell2mat(NPlatency.Long.Correct');...
                                   cell2mat(NPlatency.Long.Error')],...
                                  [  ones(length(cell2mat(NPlatency.Short.Correct')),1);...
                                   2*ones(length(cell2mat(NPlatency.Short.Error')),1);...
                                   3*ones(length(cell2mat(NPlatency.Medium.Correct')),1);...
                                   4*ones(length(cell2mat(NPlatency.Medium.Error')),1);...
                                   5*ones(length(cell2mat(NPlatency.Long.Correct')),1);...
                                   6*ones(length(cell2mat(NPlatency.Long.Error')),1)]);
    c = multcompare(stats);
    
    [p,tbl,stats] = kruskalwallis([cell2mat(Rlatency.Short.Correct');...
                                   cell2mat(Rlatency.Short.Error');...
                                   cell2mat(Rlatency.Medium.Correct');...
                                   cell2mat(Rlatency.Medium.Error');...
                                   cell2mat(Rlatency.Long.Correct');...
                                   cell2mat(Rlatency.Long.Error')],...
                                  [  ones(length(cell2mat(Rlatency.Short.Correct')),1);...
                                   2*ones(length(cell2mat(Rlatency.Short.Error')),1);...
                                   3*ones(length(cell2mat(Rlatency.Medium.Correct')),1);...
                                   4*ones(length(cell2mat(Rlatency.Medium.Error')),1);...
                                   5*ones(length(cell2mat(Rlatency.Long.Correct')),1);...
                                   6*ones(length(cell2mat(Rlatency.Long.Error')),1)]);
    c = multcompare(stats);

%% Sample stats
     [p,tbl,stats] = kruskalwallis([cell2mat(Slatency.Short.Correct');...
                                   cell2mat(Slatency.Medium.Correct');...
                                   cell2mat(Slatency.Long.Correct')],...
                                  [ones(length(cell2mat(Slatency.Short.Correct')),1);...
                                   2*ones(length(cell2mat(Slatency.Medium.Correct')),1);...
                                   3*ones(length(cell2mat(Slatency.Long.Correct')),1)]);
    c = multcompare(stats);
    
     [p,tbl,stats] = kruskalwallis([cell2mat(Slatency.Short.Error');...
                                   cell2mat(Slatency.Medium.Error');...
                                   cell2mat(Slatency.Long.Error')],...
                                  [ones(length(cell2mat(Slatency.Short.Error')),1);...
                                   2*ones(length(cell2mat(Slatency.Medium.Error')),1);...
                                   3*ones(length(cell2mat(Slatency.Long.Error')),1)]);
    c = multcompare(stats);
    
    p=[ranksum(cell2mat(Slatency.Short.Correct'),cell2mat(Slatency.Short.Error')),...
        ranksum(cell2mat(Slatency.Medium.Correct'),cell2mat(Slatency.Medium.Error')),...
        ranksum(cell2mat(Slatency.Long.Correct'),cell2mat(Slatency.Long.Error'))]<0.05
%% NP stats
     [p,tbl,stats] = kruskalwallis([cell2mat(NPlatency.Short.Correct');...
                                   cell2mat(NPlatency.Medium.Correct');...
                                   cell2mat(NPlatency.Long.Correct')],...
                                  [ones(length(cell2mat(NPlatency.Short.Correct')),1);...
                                   2*ones(length(cell2mat(NPlatency.Medium.Correct')),1);...
                                   3*ones(length(cell2mat(NPlatency.Long.Correct')),1)]);
    c = multcompare(stats);
    
     [p,tbl,stats] = kruskalwallis([cell2mat(NPlatency.Short.Error');...
                                   cell2mat(NPlatency.Medium.Error');...
                                   cell2mat(NPlatency.Long.Error')],...
                                  [ones(length(cell2mat(NPlatency.Short.Error')),1);...
                                   2*ones(length(cell2mat(NPlatency.Medium.Error')),1);...
                                   3*ones(length(cell2mat(NPlatency.Long.Error')),1)]);
    c = multcompare(stats);    
    
    p=[ranksum(cell2mat(NPlatency.Short.Correct'),cell2mat(NPlatency.Short.Error')),...
        ranksum(cell2mat(NPlatency.Medium.Correct'),cell2mat(NPlatency.Medium.Error')),...
        ranksum(cell2mat(NPlatency.Long.Correct'),cell2mat(NPlatency.Long.Error'))]<0.05
%% Choice stats
     [p,tbl,stats] = kruskalwallis([cell2mat(Rlatency.Short.Correct');...
                                   cell2mat(Rlatency.Medium.Correct');...
                                   cell2mat(Rlatency.Long.Correct')],...
                                  [ones(length(cell2mat(Rlatency.Short.Correct')),1);...
                                   2*ones(length(cell2mat(Rlatency.Medium.Correct')),1);...
                                   3*ones(length(cell2mat(Rlatency.Long.Correct')),1)]);
    c = multcompare(stats);
    
     [p,tbl,stats] = kruskalwallis([cell2mat(Rlatency.Short.Error');...
                                   cell2mat(Rlatency.Medium.Error');...
                                   cell2mat(Rlatency.Long.Error')],...
                                  [ones(length(cell2mat(Rlatency.Short.Error')),1);...
                                   2*ones(length(cell2mat(Rlatency.Medium.Error')),1);...
                                   3*ones(length(cell2mat(Rlatency.Long.Error')),1)]);
    c = multcompare(stats);    
    
    p=[ranksum(cell2mat(Rlatency.Short.Correct'),cell2mat(Rlatency.Short.Error')),...
        ranksum(cell2mat(Rlatency.Medium.Correct'),cell2mat(Rlatency.Medium.Error')),...
        ranksum(cell2mat(Rlatency.Long.Correct'),cell2mat(Rlatency.Long.Error'))]<0.05    
%% Compare first and second trial
clear p h stats chi2stat
comparisons_ = [1,2;3,4;5,6;8,9;10,11;12,13];
for iFile = 1:size(comparisons_,1)
    n1 = length([Slatency.Short.Correct{comparisons_(iFile,1)};...
                 Slatency.Medium.Correct{comparisons_(iFile,1)};...
                 ]);
    N1 = n1 + length([Slatency.Short.Error{comparisons_(iFile,1)};...
                      Slatency.Medium.Error{comparisons_(iFile,1)};...
                      ]);

    n2 = length([Slatency.Short.Correct{comparisons_(iFile,2)};...
                 Slatency.Medium.Correct{comparisons_(iFile,2)};...
                 ]);
    N2 = n2 + length([Slatency.Short.Error{comparisons_(iFile,2)};...
                      Slatency.Medium.Error{comparisons_(iFile,2)};...
                      ]);
    
    x1 = [repmat('a',N1,1); repmat('b',N2,1)];
    x2 = [repmat(1,n1,1); 
          repmat(2,N1-n1,1); 
          repmat(1,n2,1); 
          repmat(2,N2-n2,1)];
       [~,chi2stat(iFile),p(iFile)] = crosstab(x1,x2);

   

       
end
stouffer(p)
Fisher(p)
%%
comparisons_ = [1,2;3,4;5,6;8,9;10,11;12,13];
for iFile = 1:size(comparisons_,1)
    C1 = [length(Slatency.Short.Correct{comparisons_(iFile,1)});...
          length(Slatency.Medium.Correct{comparisons_(iFile,1)});...
          length(Slatency.Long.Correct{comparisons_(iFile,1)})];
    E1 = [length(Slatency.Short.Error{comparisons_(iFile,1)});...
          length(Slatency.Medium.Error{comparisons_(iFile,1)});...
          length(Slatency.Long.Error{comparisons_(iFile,1)})];
      
    C2 = [length(Slatency.Short.Correct{comparisons_(iFile,2)});...
          length(Slatency.Medium.Correct{comparisons_(iFile,2)});...
          length(Slatency.Long.Correct{comparisons_(iFile,2)})];
    E2 = [length(Slatency.Short.Error{comparisons_(iFile,2)});...
          length(Slatency.Medium.Error{comparisons_(iFile,2)});...
          length(Slatency.Long.Error{comparisons_(iFile,2)})];
    
    fracCorr(iFile,:) =[mean(C1./(C1+E1)),mean(C2./(C2+E2))];
      
end

h = ttest(fracCorr(:,1),fracCorr(:,2))

figure; plot(fracCorr'); axis([1 2 0 1])
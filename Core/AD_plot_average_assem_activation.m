% Load and plot trial-averaged activation of all detected assemblies during different phases of behavioural tasks. Includes time-resolved and time-averaged subplots
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%


%% Example plot 
clear 
cd ('C:\Analysis\AssemblyAnalysis\raw\KDE_bins\LONG')    
load('KrzysztofLONG2_iFR50_FSC.mat')

Tinv=[2 8; 9 11; 13 17; 19 21; 24 30]; % Time intervals for averaging
bw=0.05;
cmp=colormap('jet');
clr1={cmp(18,:),cmp(33,:),cmp(50,:)};
clr2={'b','g','r'};

% Ya holds mean values for each time period: Ya{time period}{region}
Ya=cell(1,5); 
for t=1:5
    Ya{t}=cell(1,3); 
end;    

figure%, hold off cla
subplot(2,1,1)
ntr=length(EvtTs)/2;              % no trials
Ltr=round(size(FSCsel{1},1)/ntr); % length of trials
evt0=EvtLs(((1:ntr)-1)*2+1)-2;    % Indices of sample events
Tv=(1:Ltr)*bw-bw/2;               % time axis for plotting

for area=1:3
    na=size(FSCsel{area},2);      % no. assemblies
    Xavg=zeros(na,Ltr);           % constructor for average assembly activation
    Xsem=Xavg;
    for i=1:na
%             subplot(3,1,area)
%         try
            A=reshape(FSCsel{area}(:,i),Ltr,ntr)';
            Xavg(i,:) = mean(A);
            Xsem(i,:) = nansem(A);
            plot(Tv,Xavg(i,:),'LineStyle','-',...
                'LineWidth',2.5,...
                'Color',clr2{area}); hold on
            ciplot(Xavg(i,:)+Xsem(i,:),...
                Xavg(i,:)-Xsem(i,:),Tv,clr1{area})
            
            % calculate mean values in each time period
            for t=1:size(Tinv,1)
                k=find(Tv>Tinv(t,1) & Tv<=Tinv(t,2));
                Ya{t}{area}=[Ya{t}{area},mean(Xavg(i,k))];
            end;
            
%         end
    end;
end
set(gca, 'XTick',5:5:25,'XTickLabel',{'pre-S','S','Delay','C','post-C'});
ylabel({'Assembly activation';'(avg. Factor score)'})


%     for s=1:3%unique(AssemType{f})
%         k=s%find(AssemType{f}==s);
%        %%%%%%%%%%%%%%%%%%%% AD replaced k with s+1
%        plot(Tv,mean(Xavg(k,:),1),'Color',clr2{s+1},'LineWidth',2.5); hold on
%        plot(Tv,mean(Xavg(s+1,:),1),'Color',clr2{s+1},'LineWidth',2.5); hold on
%         plot(Tv,mean(Xavg(s+1,:),1),'LineWidth',2.5); hold on
%     end;
    %a=input([num2str(f) '...']);


subplot(2,1,2); hold on
for t=1:5,
    for s=1:3
        YaAvg(t,s)=mean(Ya{t}{s});
        YaSE(t,s)=std(Ya{t}{s})/sqrt(length(Ya{t}{s}));
    end;
end;
% hold off cla,
errorbar(YaAvg,YaSE,'o-','LineWidth',2)
% xlim([0 6])
set(gca,'FontSize',20)
legend('PFC','HP','PFC-HP','Location','SouthWest'), legend('boxoff'),set(gca,'FontSize',10)

set(gca, 'XTick',1:5,'XTickLabel',{'pre-S','S','Delay','C','post-C'});
ylabel({'Assembly activation';'(avg. Factor score)'})
%ylabel('read-out (avg. {\itt}-score)')
text(1.85,1.9,'*','FontSize',42)
text(2.85,1.9,'*','FontSize',42)
text(3.85,1.9,'*','FontSize',42)
%
% for t=1:5
%     [~,p1]=ttest2(Ya{t}{1},Ya{t}{3});
%     [~,p2]=ttest2(Ya{t}{2},Ya{t}{3});
%     [p1 p2]
% end;
%% separate left and right choices and trials
load('KrzysztofLONG2_iFR50_FSC.mat')

Tinv=[2 8; 9 11; 13 17; 19 21; 24 30]; % Time intervals for averaging
bw=0.05;
cmp=colormap('jet');
clr1={cmp(18,:),cmp(33,:),cmp(50,:)};
clr2={'b','g','r'};

% Ya holds mean values for each time period: Ya{time period}{region}
Ya=cell(1,5); 
for t=1:5
    Ya{t}=cell(1,3); 
end;    

figure(2)%, hold off cla
subplot(2,1,1)
ntr=length(EvtTs)/2;              % no trials
Ltr=round(size(FSCsel{1},1)/ntr); % length of trials
evt0=EvtLs(((1:ntr)-1)*2+1)-2;    % Indices of sample events
Tv=(1:Ltr)*bw-bw/2;               % time axis for plotting

for area=1:3
    na=size(FSCsel{area},2);      % no. assemblies
    Xavg=zeros(na,Ltr);           % constructor for average assembly activation
    Xsem=Xavg;
    for i=1:na
%             subplot(3,1,area)
%         try
            A=reshape(FSCsel{area}(:,i),Ltr,ntr)';
            Xavg(i,:) = mean(A);
            Xsem(i,:) = nansem(A);
            plot(Tv,Xavg(i,:),'LineStyle','-',...
                'LineWidth',2.5,...
                'Color',clr2{area}); hold on
            ciplot(Xavg(i,:)+Xsem(i,:),...
                Xavg(i,:)-Xsem(i,:),Tv,clr1{area})
            
            % calculate mean values in each time period
            for t=1:size(Tinv,1)
                k=find(Tv>Tinv(t,1) & Tv<=Tinv(t,2));
                Ya{t}{area}=[Ya{t}{area},mean(Xavg(i,k))];
            end;
            
%         end
    end;
end
set(gca, 'XTick',5:5:25,'XTickLabel',{'pre-S','S','Delay','C','post-C'});
ylabel({'Assembly activation';'(avg. Factor score)'})


%     for s=1:3%unique(AssemType{f})
%         k=s%find(AssemType{f}==s);
%        %%%%%%%%%%%%%%%%%%%% AD replaced k with s+1
%        plot(Tv,mean(Xavg(k,:),1),'Color',clr2{s+1},'LineWidth',2.5); hold on
%        plot(Tv,mean(Xavg(s+1,:),1),'Color',clr2{s+1},'LineWidth',2.5); hold on
%         plot(Tv,mean(Xavg(s+1,:),1),'LineWidth',2.5); hold on
%     end;
    %a=input([num2str(f) '...']);


subplot(2,1,2); hold on
for t=1:5,
    for s=1:3
        YaAvg(t,s)=mean(Ya{t}{s});
        YaSE(t,s)=std(Ya{t}{s})/sqrt(length(Ya{t}{s}));
    end;
end;
% hold off cla,
errorbar(YaAvg,YaSE,'o-','LineWidth',2)
% xlim([0 6])
set(gca,'FontSize',20)
legend('PFC','HP','PFC-HP','Location','SouthWest'), legend('boxoff'),set(gca,'FontSize',10)

set(gca, 'XTick',1:5,'XTickLabel',{'pre-S','S','Delay','C','post-C'});
ylabel({'Assembly activation';'(avg. Factor score)'})
%ylabel('read-out (avg. {\itt}-score)')
text(1.85,1.9,'*','FontSize',42)
text(2.85,1.9,'*','FontSize',42)
text(3.85,1.9,'*','FontSize',42)
%
% for t=1:5
%     [~,p1]=ttest2(Ya{t}{1},Ya{t}{3});
%     [~,p2]=ttest2(Ya{t}{2},Ya{t}{3});
%     [p1 p2]
% end;


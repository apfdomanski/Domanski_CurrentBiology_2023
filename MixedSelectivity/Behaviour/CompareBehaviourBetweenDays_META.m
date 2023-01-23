clear 
pat = 'C:\Analysis\AssemblyAnalysis\raw';
SHORT = load(sprintf('%s\\allTimestamps\\Behaviour analysis\\DayByDay_SHORT.mat',pat),'Slatency','NPlatency','Rlatency');
MEDIUM = load(sprintf('%s\\allTimestamps\\Behaviour analysis\\DayByDay_MEDIUM.mat',pat),'Slatency','NPlatency','Rlatency');
LONG = load(sprintf('%s\\allTimestamps\\Behaviour analysis\\DayByDay_LONG.mat',pat),'Slatency','NPlatency','Rlatency');

Delays__ = {'0s','2s','4s','6s','8s','4s','8s','16s'}    ;
Slatency.Delay_0 = SHORT.Slatency.Delay_0;
Slatency.Delay_0 = SHORT.Slatency.Delay_0;
Slatency.Delay_2 = MEDIUM.Slatency.Delay_2;
Slatency.Delay_4 = MEDIUM.Slatency.Delay_4;
Slatency.Delay_6 = MEDIUM.Slatency.Delay_6;
Slatency.Delay_8 = MEDIUM.Slatency.Delay_8;
Slatency.Short	 = LONG.Slatency.Short;
Slatency.Medium	 = LONG.Slatency.Medium;
Slatency.Long	 = LONG.Slatency.Long;

Rlatency.Delay_0 = SHORT.Rlatency.Delay_0;
Rlatency.Delay_0 = SHORT.Rlatency.Delay_0;  
Rlatency.Delay_2 = MEDIUM.Rlatency.Delay_2;
Rlatency.Delay_4 = MEDIUM.Rlatency.Delay_4;
Rlatency.Delay_6 = MEDIUM.Rlatency.Delay_6;
Rlatency.Delay_8 = MEDIUM.Rlatency.Delay_8;
Rlatency.Short	 = LONG.Rlatency.Short;
Rlatency.Medium	 = LONG.Rlatency.Medium;
Rlatency.Long	 = LONG.Rlatency.Long;

NPlatency.Delay_0 = SHORT.NPlatency.Delay_0;
NPlatency.Delay_0 = SHORT.NPlatency.Delay_0;
NPlatency.Delay_2 = MEDIUM.NPlatency.Delay_2;
NPlatency.Delay_4 = MEDIUM.NPlatency.Delay_4;
NPlatency.Delay_6 = MEDIUM.NPlatency.Delay_6;
NPlatency.Delay_8 = MEDIUM.NPlatency.Delay_8;
NPlatency.Short	  = LONG.NPlatency.Short;
NPlatency.Medium  = LONG.NPlatency.Medium;
NPlatency.Long	  = LONG.NPlatency.Long;

clear SHORT MEDIUM LONG pat
%% bars of correct percentages
clear x

x{1} = cellfun(@length,Slatency.Delay_0.Correct)./(cellfun(@length,Slatency.Delay_0.Correct) + cellfun(@length,Slatency.Delay_0.Error))*100;
x{2} = cellfun(@length,Slatency.Delay_2.Correct)./(cellfun(@length,Slatency.Delay_2.Correct) + cellfun(@length,Slatency.Delay_2.Error))*100;
x{3} = cellfun(@length,Slatency.Delay_4.Correct)./(cellfun(@length,Slatency.Delay_4.Correct) + cellfun(@length,Slatency.Delay_4.Error))*100;
x{4} = cellfun(@length,Slatency.Delay_6.Correct)./(cellfun(@length,Slatency.Delay_6.Correct) + cellfun(@length,Slatency.Delay_6.Error))*100;
x{5} = cellfun(@length,Slatency.Delay_8.Correct)./(cellfun(@length,Slatency.Delay_8.Correct) + cellfun(@length,Slatency.Delay_8.Error))*100;
x{6} = cellfun(@length,Slatency.Short.Correct)./(cellfun(@length,Slatency.Short.Correct) + cellfun(@length,Slatency.Short.Error))*100;
x{7} = cellfun(@length,Slatency.Medium.Correct)./(cellfun(@length,Slatency.Medium.Correct) + cellfun(@length,Slatency.Medium.Error))*100;
x{8} = cellfun(@length,Slatency.Long.Correct)./(cellfun(@length,Slatency.Long.Correct) + cellfun(@length,Slatency.Long.Error))*100;

figure('color','w'); hold on
h = gca;
% plot([0.5 1.5],[90 90],'k','LineWidth',1.5)
% plot([2.5 6.5],[90 90],'k','LineWidth',1.5)
% plot([7.5 10.5],[90 90],'k','LineWidth',1.5)
plot([0 11],[50 50],'r')
area([0.2 1.8],[90 90],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
area([2.2 6.8],[90 90],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
area([7.2 10.8],[90 90],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
text(mean([0.2 1.8]),95,'Early','HorizontalAlignment','center')
text(mean([2.2 6.8]),95,'Middle','HorizontalAlignment','center')
text(mean([7.2 10.8]),95,'Late','HorizontalAlignment','center')

bar([1,3,4,5,6,8,9,10],cellfun(@nanmean,x),'FaceColor','w','EdgeColor','k','LineWidth',1.5)
% errorbar([1,3,4,5,6,8,9,10],cellfun(@nanmean,x),cellfun(@nansem,x),'.k','Marker','none','LineWidth',1.5)

errorbar([1],nanmean(x{1}),nansem(x{1}),'.k','Marker','none','LineWidth',1.5)
errorbar([3,4,5,6],cellfun(@nanmean,x(2:5)),cellfun(@nansem,x(2:5)),'-k','Marker','none','LineWidth',1.5)
errorbar([8,9,10],cellfun(@nanmean,x(6:8)),cellfun(@nansem,x(6:8)),'-k','Marker','none','LineWidth',1.5)
for iSession = 1:length(x{2})
    plot([3,4,5,6],[x{2}(iSession),x{3}(iSession),x{4}(iSession),x{5}(iSession)],':k')
end
for iSession = 1:length(x{6})
    plot([8,9,10],[x{6}(iSession),x{7}(iSession),x{8}(iSession)],':k')
end    

set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTickLabel',Delays__,'XTickLabelRotation',45)
    xlabel('Delay duration')
    ylabel('% Correct choices')
    ylim([0 100])
%     uistack(h, 'top')
    h.Layer = 'top';
%%
clear y_ y__
y_ = [x{6};x{7}]; 
% strip Miroslaw LONG0
% x{6}(7)=NaN;x{7}(7)=NaN;
% y_(:,7)=NaN;

for i = 1:size(y_,2)
    m = diff(y_(:,i))./4;
    y__(i) = y_(2,i)+8*m;
end
figure
subplot(1,2,1); hold on
errorbar(0,nanmean(x{1}),nansem(x{1}),'Color',[0 0 0.6],'Marker','none','LineWidth',1.5)
errorbar([2,4,6,8],cellfun(@nanmean,x(2:5)),cellfun(@nansem,x(2:5)),'Color',[0.1 0.6 0.9],'Marker','none','LineWidth',1.5)
errorbar([4,8,16],cellfun(@nanmean,x(6:8)),cellfun(@nansem,x(6:8)),'Color',[0.6 0.6 0.9],'Marker','none','LineWidth',1.5)
plot([0 16],[50 50],':r')
xlabel('Delay duration (s)')
ylabel('% Correct choices')
title('Performance vs. Training stage')
axis([-1 17 0 100])
legend('Early','Middle','Late','Location','south','Orientation','vertical'); legend boxoff

subplot(1,2,2); hold on


axis([-1 18 0 100])
path(rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite'))
[h,p] = ttest2(y__',x{8}');
% [h,p] = ttest2(y__',x{8}');
[p,h] = ranksum(y__',x{8}');
path(addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite'))



errorbar(16,nanmean(y__),nansem(y__),'Color',[0.6 0.6 0.6 0.9],'Marker','o','LineWidth',1.5)
errorbar([4,8,16],cellfun(@nanmean,x(6:8)),cellfun(@nansem,x(6:8)),'Color',[0.6 0.6 0.9],'Marker','o','LineWidth',1.5)
plot([8,16],[nanmean(x{7}),nanmean(y__)],'color',[0.6 0.6 0.6],'LineWidth',1.5)
plot([0 16],[50 50],':r')


% scatter(16*ones(length(y__),1),y__)
plot(repmat([4,8],length(x{6}),1)',[x{6};x{7}],'Color',[0.6 0.6 0.9 0.9])
% plot(repmat([4,8,16],length(x{6}),1)',[x{6};x{7};x{8}],'Color',[0.6 0.6 0.6 0.9])
plot(repmat([8,16],size(y_,2),1)',[x{7};y__],'Color',[0.6 0.6 0.6 0.9],'LineStyle',':')

if ~h
    txt = sprintf('Data vs linear extrapolation: n.s.(p=%0.3f)',p);
else
    txt = sprintf('Data vs linear extrapolation: * (p=%0.3f)',p);
end
text(0,20,txt)
xlabel('Delay duration (s)')
ylabel('% Correct choices')
title({'Linear decline with delay duration ';'or performance collapse?'})
legend('[4, 8, 16s] data','[4, 8s] linear extrapolation to 16s','Location','south','Orientation','vertical'); legend boxoff

%% bars of latencies
figure; 
subplot(1,3,1);hold on
x=[]; y= [];

x=[nanmean(cell2mat(Slatency.Delay_0.Error')),...
   nanmean(cell2mat(Slatency.Delay_2.Error')),...
   nanmean(cell2mat(Slatency.Delay_4.Error')),...
   nanmean(cell2mat(Slatency.Delay_6.Error')),...
   nanmean(cell2mat(Slatency.Delay_8.Error')),...
   nanmean(cell2mat(Slatency.Short.Error')),...
   nanmean(cell2mat(Slatency.Medium.Error')),...
   nanmean(cell2mat(Slatency.Long.Error'))];

y=[nansem(cell2mat(Slatency.Delay_0.Error')),...
   nansem(cell2mat(Slatency.Delay_2.Error')),...
   nansem(cell2mat(Slatency.Delay_4.Error')),...
   nansem(cell2mat(Slatency.Delay_6.Error')),...
   nansem(cell2mat(Slatency.Delay_8.Error')),...
   nansem(cell2mat(Slatency.Short.Error')),...
   nansem(cell2mat(Slatency.Medium.Error')),...
   nansem(cell2mat(Slatency.Long.Error'))];
area([0.2 1.8],[9 9],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
area([2.2 6.8],[9 9],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
area([7.2 10.8],[9 9],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
text(mean([0.2 1.8]),9.5,'Early','HorizontalAlignment','center')
text(mean([2.2 6.8]),9.5,'Middle','HorizontalAlignment','center')
text(mean([7.2 10.8]),9.5,'Late','HorizontalAlignment','center')

bar([1,3,4,5,6,8,9,10],x,'FaceColor','w','EdgeColor','r','LineWidth',1.5)
errorbar([1,3,4,5,6,8,9,10],x,y,'.r','Marker','none','LineWidth',1.5)

x=[]; y= [];

x=[nanmean(cell2mat(Slatency.Delay_0.Correct')),...
   nanmean(cell2mat(Slatency.Delay_2.Correct')),...
   nanmean(cell2mat(Slatency.Delay_4.Correct')),...
   nanmean(cell2mat(Slatency.Delay_6.Correct')),...
   nanmean(cell2mat(Slatency.Delay_8.Correct')),...
   nanmean(cell2mat(Slatency.Short.Correct')),...
   nanmean(cell2mat(Slatency.Medium.Correct')),...
   nanmean(cell2mat(Slatency.Long.Correct'))];

y=[nansem(cell2mat(Slatency.Delay_0.Correct')),...
   nansem(cell2mat(Slatency.Delay_2.Correct')),...
   nansem(cell2mat(Slatency.Delay_4.Correct')),...
   nansem(cell2mat(Slatency.Delay_6.Correct')),...
   nansem(cell2mat(Slatency.Delay_8.Correct')),...
   nansem(cell2mat(Slatency.Short.Correct')),...
   nansem(cell2mat(Slatency.Medium.Correct')),...
   nansem(cell2mat(Slatency.Long.Correct'))];


bar([1,3,4,5,6,8,9,10],x,'FaceColor','w','EdgeColor','k','LineWidth',1.5)
errorbar([1,3,4,5,6,8,9,10],x,y,'.k','Marker','none','LineWidth',1.5)

    set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTickLabel',Delays__,'XTickLabelRotation',45)
%     xlabel('Delay duration')
    ylabel('Latency (s)')
    title('Sample press latency')
    ylim([0 10])
    h=gca; h.Layer = 'top';

subplot(1,3,2);hold on

x=[]; y= [];

x=[nanmean(cell2mat(NPlatency.Delay_0.Error')),...
   nanmean(cell2mat(NPlatency.Delay_2.Error')),...
   nanmean(cell2mat(NPlatency.Delay_4.Error')),...
   nanmean(cell2mat(NPlatency.Delay_6.Error')),...
   nanmean(cell2mat(NPlatency.Delay_8.Error')),...
   nanmean(cell2mat(NPlatency.Short.Error')),...
   nanmean(cell2mat(NPlatency.Medium.Error')),...
   nanmean(cell2mat(NPlatency.Long.Error'))];

y=[nansem(cell2mat(NPlatency.Delay_0.Error')),...
   nansem(cell2mat(NPlatency.Delay_2.Error')),...
   nansem(cell2mat(NPlatency.Delay_4.Error')),...
   nansem(cell2mat(NPlatency.Delay_6.Error')),...
   nansem(cell2mat(NPlatency.Delay_8.Error')),...
   nansem(cell2mat(NPlatency.Short.Error')),...
   nansem(cell2mat(NPlatency.Medium.Error')),...
   nansem(cell2mat(NPlatency.Long.Error'))];

area([0.2 1.8],7.2*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
area([2.2 6.8],7.2*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
area([7.2 10.8],7.2*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
text(mean([0.2 1.8]),7.6,'Early','HorizontalAlignment','center')
text(mean([2.2 6.8]),7.6,'Middle','HorizontalAlignment','center')
text(mean([7.2 10.8]),7.6,'Late','HorizontalAlignment','center')

bar([1,3,4,5,6,8,9,10],x,'FaceColor','w','EdgeColor','r','LineWidth',1.5)
errorbar([1,3,4,5,6,8,9,10],x,y,'.r','Marker','none','LineWidth',1.5)

x=[]; y= [];

x=[nanmean(cell2mat(NPlatency.Delay_0.Correct')),...
   nanmean(cell2mat(NPlatency.Delay_2.Correct')),...
   nanmean(cell2mat(NPlatency.Delay_4.Correct')),...
   nanmean(cell2mat(NPlatency.Delay_6.Correct')),...
   nanmean(cell2mat(NPlatency.Delay_8.Correct')),...
   nanmean(cell2mat(NPlatency.Short.Correct')),...
   nanmean(cell2mat(NPlatency.Medium.Correct')),...
   nanmean(cell2mat(NPlatency.Long.Correct'))];

y=[nansem(cell2mat(NPlatency.Delay_0.Correct')),...
   nansem(cell2mat(NPlatency.Delay_2.Correct')),...
   nansem(cell2mat(NPlatency.Delay_4.Correct')),...
   nansem(cell2mat(NPlatency.Delay_6.Correct')),...
   nansem(cell2mat(NPlatency.Delay_8.Correct')),...
   nansem(cell2mat(NPlatency.Short.Correct')),...
   nansem(cell2mat(NPlatency.Medium.Correct')),...
   nansem(cell2mat(NPlatency.Long.Correct'))];

bar([1,3,4,5,6,8,9,10],x,'FaceColor','w','EdgeColor','k','LineWidth',1.5)
errorbar([1,3,4,5,6,8,9,10],x,y,'.k','Marker','none','LineWidth',1.5)

set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTickLabel',Delays__,'XTickLabelRotation',45)
xlabel('Delay duration')
ylabel('Latency (s)')
title('Nosepoke latency')

ylim([0 8])
h=gca; h.Layer = 'top';

subplot(1,3,3);hold on
area([0.2 1.8],3.6*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
area([2.2 6.8],3.6*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
area([7.2 10.8],3.6*[1 1],'FaceColor',[0.9 0.9 0.9],'LineStyle','none')
text(mean([0.2 1.8]),3.8,'Early','HorizontalAlignment','center')
text(mean([2.2 6.8]),3.8,'Middle','HorizontalAlignment','center')
text(mean([7.2 10.8]),3.8,'Late','HorizontalAlignment','center')
x=[]; y= [];

x=[nanmean(cell2mat(Rlatency.Delay_0.Error')),...
   nanmean(cell2mat(Rlatency.Delay_2.Error')),...
   nanmean(cell2mat(Rlatency.Delay_4.Error')),...
   nanmean(cell2mat(Rlatency.Delay_6.Error')),...
   nanmean(cell2mat(Rlatency.Delay_8.Error')),...
   nanmean(cell2mat(Rlatency.Short.Error')),...
   nanmean(cell2mat(Rlatency.Medium.Error')),...
   nanmean(cell2mat(Rlatency.Long.Error'))];

y=[nansem(cell2mat(Rlatency.Delay_0.Error')),...
   nansem(cell2mat(Rlatency.Delay_2.Error')),...
   nansem(cell2mat(Rlatency.Delay_4.Error')),...
   nansem(cell2mat(Rlatency.Delay_6.Error')),...
   nansem(cell2mat(Rlatency.Delay_8.Error')),...
   nansem(cell2mat(Rlatency.Short.Error')),...
   nansem(cell2mat(Rlatency.Medium.Error')),...
   nansem(cell2mat(Rlatency.Long.Error'))];
    bar([1,3,4,5,6,8,9,10],x,'FaceColor','w','EdgeColor','r','LineWidth',1.5)
    errorbar([1,3,4,5,6,8,9,10],x,y,'.r','Marker','none','LineWidth',1.5)

x=[]; y= [];

x=[nanmean(cell2mat(Rlatency.Delay_0.Correct')),...
   nanmean(cell2mat(Rlatency.Delay_2.Correct')),...
   nanmean(cell2mat(Rlatency.Delay_4.Correct')),...
   nanmean(cell2mat(Rlatency.Delay_6.Correct')),...
   nanmean(cell2mat(Rlatency.Delay_8.Correct')),...
   nanmean(cell2mat(Rlatency.Short.Correct')),...
   nanmean(cell2mat(Rlatency.Medium.Correct')),...
   nanmean(cell2mat(Rlatency.Long.Correct'))];

y=[nansem(cell2mat(Rlatency.Delay_0.Correct')),...
   nansem(cell2mat(Rlatency.Delay_2.Correct')),...
   nansem(cell2mat(Rlatency.Delay_4.Correct')),...
   nansem(cell2mat(Rlatency.Delay_6.Correct')),...
   nansem(cell2mat(Rlatency.Delay_8.Correct')),...
   nansem(cell2mat(Rlatency.Short.Correct')),...
   nansem(cell2mat(Rlatency.Medium.Correct')),...
   nansem(cell2mat(Rlatency.Long.Correct'))];
    bar([1,3,4,5,6,8,9,10],x,'FaceColor','w','EdgeColor','k','LineWidth',1.5)
    errorbar([1,3,4,5,6,8,9,10],x,y,'.k','Marker','none','LineWidth',1.5)

set(gca,'XTick',[1,3,4,5,6,8,9,10],'XTickLabel',Delays__,'XTickLabelRotation',45)
%     xlabel('Delay duration')
    ylabel('Latency (s)')
    title('Choice press latency')
    ylim([0 4])
h=gca; h.Layer = 'top';

% Plot the L/R preference of each assembly.
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Need data from runFAassem first.

clear 
cd('C:\Analysis\AssemblyAnalysis\raw\KDE_bins\LONG')
fn='JaroslawLONG2_iFR50_FSC.mat';
% fn='KrzesimirLONG2_iFR50_FSC.mat';

twin=10;       % Time window
bw=0.05;
ko=round((twin+5)/bw);
% [TmtxS,iFRs,EvtLs,EvtTs,usel_out] = SelTrialsCells(fn,twin,minFR,critCvBW,'iFR'); %'iFRsc'
[TmtxS,FSCs,EvtLs,EvtTs]          = SelTrialsFSCs(fn,twin);

Vstr={'left_choice','right_choice','left_sample','right_sample'};

no_assem_types=length(FSCs);
ntr=length(EvtTs)/2;
Ltr=round(size(FSCs{1},1)/ntr);
evt0=EvtLs(((1:ntr)-1)*2+1)-2;
for s=1:no_assem_types
    for i=1:length(TmtxS{s})
        TmtxS{s}{i}=TmtxS{s}{i}([1:ko end-ko+1:end]);
        FSCs{s}{i}=FSCs{s}{i}([1:ko end-ko+1:end],:)';
    end;
    FSCsel{s}=cell2mat(FSCs{s})'; n_assem(s)=size(FSCs{s},2);
end;


% plot(EvtTs,EvtLs)
for s=1:no_assem_types
    for i=1:length(FSCs{s})
        for u=1:size(FSCs{s}{i},1)
            [s i u] %[area trial_no assem_no] 
            FSCs_{s}{u}(i,:)=FSCs{s}{i}(u,:);
        end
    end
end

for s=1:no_assem_types
    for u=1:length(FSCs_{s})
            FSCs_sort.mean1{s}(u,:) = nanmean(FSCs_{s}{u}(evt0==1,:),1);
            FSCs_sort.mean2{s}(u,:) = nanmean(FSCs_{s}{u}(evt0==2,:),1);
            FSCs_sort.SEM1{s}(u,:)  = nansem(FSCs_{s}{u}(evt0==1,:),1);
            FSCs_sort.SEM2{s}(u,:)  = nansem(FSCs_{s}{u}(evt0==2,:),1);

    end
end
%%
tb=((1:ko*2))*bw*2/3;
area={'PFC','HP','HP-PFC'};
for s=1:no_assem_types
    for u=1:length(FSCs_{s})
        figure('name',[area{s} ' Assem no. ' num2str(u)],'color','w'); hold on

        %     subplot(ceil(sqrt(length(FSCs_{s}))),ceil(sqrt(length(FSCs_{s}))),u); hold on
%         subaxis(floor(sqrt(length(FSCs_{s}))),ceil(sqrt(length(FSCs_{s}))), u ,...
%             'Spacing', 0.05, 'Padding', 0.05, 'Margin', 0.05); hold on
        plot(tb,FSCs_sort.mean1{s}(u,:),'color',[0.3 0.3 0.3],'LineWidth',1.5)
        plot(tb,FSCs_sort.mean2{s}(u,:),'color',[0.6 0.6 0.6],'LineWidth',1.5)
        ciplot(FSCs_sort.mean1{s}(u,:)+FSCs_sort.SEM1{s}(u,:),...
               FSCs_sort.mean1{s}(u,:)-FSCs_sort.SEM1{s}(u,:),tb,[0.3 0.3 0.3])
        ciplot(FSCs_sort.mean2{s}(u,:)+FSCs_sort.SEM2{s}(u,:),...
               FSCs_sort.mean2{s}(u,:)-FSCs_sort.SEM2{s}(u,:),tb,[0.6 0.6 0.6])
        plot(tb,FSCs_sort.mean1{s}(u,:),'color',[0.3 0.3 0.3],'LineWidth',1.5)
        plot(tb,FSCs_sort.mean2{s}(u,:),'color',[0.6 0.6 0.6],'LineWidth',1.5)
        y_range=[min([FSCs_sort.mean1{s}(u,:) FSCs_sort.mean2{s}(u,:)]) max([FSCs_sort.mean1{s}(u,:) FSCs_sort.mean2{s}(u,:)]) ];
        text(4.5,y_range(2),'Sample','color','g','Rotation',90);
        text(14.5,y_range(2),'Choice','color','r','Rotation',90);

        plot([twin twin],1.5*y_range,':k','LineWidth', 1.5)
        plot([5 5],1.5*y_range,'g','LineWidth', 1.5)
        plot([twin*1.5 twin*1.5],1.5*y_range,'r','LineWidth', 1.5)
         %     text(5, 1600,num2str(u))
        	axis off
        set(gca,'Xtick',[]) 
        set(gca,'Ytick',[])
        legend({'Left Trials' 'Right Trials'},'location','SouthEast'); legend boxoff
%         plot([1 1],[y_range(2)-4 y_range(2)-2],'k','LineWidth',2)
%         plot([1 3],[y_range(2)-4 y_range(2)-4],'k','LineWidth',2)
%         text(1.5,y_range(2)-3,{'2s'})
%         axis([0 20 2*y_range])
        plot([1 3],[5 5],'k','LineWidth',2)
        text(1.8,5.5,{'2s'})

        axis([0 20 -5 10])


    end
    
end




       
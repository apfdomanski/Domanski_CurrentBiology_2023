% Plot the L/R preference of each single unit.
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
clear 
fn='JaroslawLONG2_PFC_iFR50.mat';
% fn = 'JaroslawLONG2_PFC_iFR50_behavOnly.mat';
critCvBW=1e6;  % critical max bound on kernel width drift over time (variance)
minFR=0.01;     % minimal acceptable firing rate
twin=10;       % Time window
bw=0.05;
ko=round((twin+5)/bw);
[TmtxS,iFRs,EvtLs,EvtTs,usel_out]=SelTrialsCells(fn,twin,minFR,critCvBW,'iFR'); %'iFRsc'
Vstr={'left_choice','right_choice','left_sample','right_sample'};

ntr=length(EvtTs)/2;
Ltr=round(size(iFRs{1},1)/ntr);
evt0=EvtLs(((1:ntr)-1)*2+1)-2;

for s=1:2
    for i=1:length(TmtxS{s})
        TmtxS{s}{i}=TmtxS{s}{i}([1:ko end-ko+1:end]);
        iFRs{s}{i}=iFRs{s}{i}([1:ko end-ko+1:end],:)';
    end;
    iFRsel{s}=cell2mat(iFRs{s})'; n_assem(s)=size(iFRs{s},2);
end;

    
% plot(EvtTs,EvtLs)
for s=1:2
    for i=1:length(iFRs{s})
        for u=1:size(iFRs{s}{i},1)
            iFRs_{s}{u}(i,:)=iFRs{s}{i}(u,:);
        end
    end
end

for s=1:2
    for u=1:length(iFRs_{s})
            iFRs_sort.mean1{s}(u,:) = nanmean(iFRs_{s}{u}(evt0==1,:),1);
            iFRs_sort.mean2{s}(u,:) = nanmean(iFRs_{s}{u}(evt0==2,:),1);
            iFRs_sort.SEM1{s}(u,:)  = nansem(iFRs_{s}{u}(evt0==1,:),1);
            iFRs_sort.SEM2{s}(u,:)  = nansem(iFRs_{s}{u}(evt0==2,:),1);

    end
end
%%
tb=(1:ko*2)/(ko*2)*2*twin;
area={'PFC','HP'};
for s=1:2
    figure('name',area{s})
    for u=1:length(iFRs_{s})
        %     subplot(ceil(sqrt(length(iFRs_{s}))),ceil(sqrt(length(iFRs_{s}))),u); hold on
        subaxis(ceil(sqrt(length(iFRs_{s}))),ceil(sqrt(length(iFRs_{s}))),u,...
                'Spacing', 0.05, 'Padding', 0.005, 'Margin', 0.005); hold on
        ciplot(iFRs_sort.mean1{s}(u,:)+iFRs_sort.SEM1{s}(u,:),...
            iFRs_sort.mean1{s}(u,:)-iFRs_sort.SEM1{s}(u,:),tb,'r')
        ciplot(iFRs_sort.mean2{s}(u,:)+iFRs_sort.SEM2{s}(u,:),...
            iFRs_sort.mean2{s}(u,:)-iFRs_sort.SEM2{s}(u,:),tb,'g')
        plot(tb,iFRs_sort.mean1{s}(u,:),'r','LineWidth',1.5)
        plot(tb,iFRs_sort.mean2{s}(u,:),'g','LineWidth',1.5)
        axis([0 Inf 0 Inf])
        %     text(5, 1600,num2str(u))
        %	axis off
        set(gca,'Xtick',[])
        set(gca,'Ytick',[])
        
        plot([twin twin],1.5*[ 0 max([iFRs_sort.mean1{s}(u,:) iFRs_sort.mean2{s}(u,:)])],':k','LineWidth', 1.5)
        plot([twin/2 twin/2],1.5*[ 0 max([iFRs_sort.mean1{s}(u,:) iFRs_sort.mean2{s}(u,:)])],':k','LineWidth', 1.5)
        plot([twin*1.5 twin*1.5],1.5*[ 0 max([iFRs_sort.mean1{s}(u,:) iFRs_sort.mean2{s}(u,:)])],':k','LineWidth', 1.5)
    end   
end




       
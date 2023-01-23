% Explore single unit vs. assembly performance in decoding behaviour
clear all; close all
target = 'LONG';
pat  = 'C:\Analysis\AssemblyAnalysis\raw\KDE_bins';

% A=dir([pat '*PFC_iFR50*.mat']);
A=dir([pat filesep '*' target '*', '*PFC_iFR50*.mat']);
for i=1:length(A)
    k=findstr(A(i).name,'_'); r=findstr(A(i).name,'.');
    Name{i}=[A(i).name(1:k(1)) A(i).name(k(2)+1:r-1)];
end;
reg  = 0.05;
bw   = 0.05;
Tinv = [0 6;    ...
        9 10.5; ...
        13 17;  ...
        19 20.5;...
        24 30];
twin = 10; 
minFR = 0.1;
critCvBW = 1e6;

%% 

% average assembly decoding
prFt2=cell(1,2); Rt2all=cell(1,2);
for f=1:length(A)
    fnOut=[A(f).name]; load(fnOut);
    fnOut=[pat filesep target filesep Name{f} '_FSC']; load(fnOut);
%     EvtTs=EvtT; %AD Added
%     EvtLs=EvtL; %AD Added
    ntr=length(EvtTs)/2; Ltr=round(size(FSCsel{1},1)/ntr);
    evt0=EvtLs(((1:ntr)-1)*2+1)-2;
%     evt0=EvtLs(((1:ntr)-1)*2+1);
    Tv=(1:Ltr)*bw-bw/2;
    disp(['... Working on file ' Name{f} '.'])
    fprintf('... %G / %G / %G assemblies detected.\n',cellfun(@length,units))
    for s=1:3
        if s==3
            uL=false(1,length(units{s}));
            AssemType=zeros(1,length(units{s}));
            for u=1:length(units{s})
                uL(u)=1-(isempty(find(units{s}{u}>nu(1))) | ...
                         isempty(find(units{s}{u}<=nu(1))));
                AssemType(u)=1-isempty(find(units{s}{u}>nu(1)))+uL(u);
            end;
            FSCsel{s}=FSCsel{s}(:,uL);
        end;
       
        if ~isempty(FSCsel{s})
            [~,~,Ft2,Rt2,~,~,dfeff]=DecodeStats(FSCsel{s},evt0,reg);
            q=size(FSCsel{s},2); 
            n1=length(find(evt0==1)); % no trials for outcome=1 
            n2=length(find(evt0==2)); % no trials for outcome=2
            prFt2{s}(f,:)=1-fcdf(Ft2,q,n1+n2-dfeff-1);
            Rt2all{s}(f,:)=Rt2;
        else
            prFt2{s}(f,:)=zeros(size(Ft2))+nan; Rt2all{s}(f,:)=zeros(size(Ft2))+nan;
        end;
    end;
end;

subplot(1,2,2), hold off cla
ms='Rt2all';
P=cell(1,2); Px=P;
Clr={'b','r','g'};
ksel=[]; for s=1:3, v=eval([ms '{s}']); ksel{s}=find(~isnan(v(:,1))); end;
ks=intersect(ksel{1},ksel{2});
for s=1:3
    v=eval([ms '{s}']);
    k=find(~isnan(v(:,1))); v=v(k,:); N=size(v,1);
    %v=v(ks,:); N=size(v,1);
    %errorbar(Tv,mean(v),std(v)./sqrt(N),'o-','Color',Clr{s},'LineWidth',2);
    %plot(Tv,mean(v),'o-','Color',Clr{s},'LineWidth',2);
    y0=LocLinReg(Tv',mean(v)',Tv',0.5);
    plot(Tv,y0,'Color',Clr{s},'LineWidth',2);
    hold on
    for i=1:size(Tinv,1)
        k=find(Tv>=Tinv(i,1) & Tv<Tinv(i,2));
        P{s}(:,i)=mean(v(:,k)')';
        Px{s}(:,i)=max(v(:,k)')';
        %plot(mean(Tv(k)),mean(Px{s}(:,i)),'o','Color',Clr{s},'LineWidth',2);
    end;
end;
%plot(Tv,zeros(size(Tv))+0.5,'k','LineWidth',2)
box off, set(gca,'FontSize',16); xlabel('time (s)');
%ylabel('pr_F'); ylim([0 0.6])
ylabel('<R^2>'); ylim([0 0.2]);
legend('PFC','HP','PFC-HP'); legend boxoff
title('Assem-Decoding');


% average rate decoding

prFt2=cell(1,2); Rt2all=cell(1,2);
for f=1:length(A)
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fn=[pat A(f).name];
    [TmtxS,iFR0,EvtLs,EvtTs]=SelTrialsCells(fn,twin,minFR,critCvBW,'iFR');
    ko=round((twin+5)/bw);
    iFRsel=cell(1,2);
    for s=1:2
        for i=1:length(TmtxS{s})
            TmtxS{s}{i}=TmtxS{s}{i}([1:ko end-ko+1:end]);
            iFR0{s}{i}=iFR0{s}{i}([1:ko end-ko+1:end],:)';
        end;
        iFRsel{s}=cell2mat(iFR0{s})';
    end;
    Tmtx=cell2mat(TmtxS{1})';
    
    ntr=length(EvtTs)/2; Ltr=round(size(iFRsel{1},1)/ntr);
    evt0=EvtLs(((1:ntr)-1)*2+1)-2;
    Tv=(1:Ltr)*bw-bw/2;
    n1=length(find(evt0==1)); n2=length(find(evt0==2));
    for s=1:2
        [~,~,Ft2,Rt2,~,~,dfeff]=DecodeStats(iFRsel{s},evt0,reg);
        q=size(iFRsel{s},2); 
        prFt2{s}(f,:)=1-fcdf(Ft2,q,n1+n2-dfeff-1);
        Rt2all{s}(f,:)=Rt2;
    end;
end;


subplot(1,2,1), hold off cla
ms='Rt2all';
Q=cell(1,2); Qx=Q;
Clr={'b','r','g'};
for s=1:2
    v=eval([ms '{s}']);
    k=find(~isnan(v(:,1))); v=v(k,:); N=size(v,1);
    y0=LocLinReg(Tv',mean(v)',Tv',0.5);
    plot(Tv,y0,'Color',Clr{s},'LineWidth',2);
    hold on
    for i=1:size(Tinv,1)
        k=find(Tv>=Tinv(i,1) & Tv<Tinv(i,2));
        Q{s}(:,i)=mean(v(:,k)')';
        Qx{s}(:,i)=max(v(:,k)')';
        %plot(mean(Tv(k)),mean(Qx{s}(:,i)),'o','Color',Clr{s},'LineWidth',2);
    end;
end;
%plot(Tv,zeros(size(Tv))+0.5,'k','LineWidth',2)
box off, set(gca,'FontSize',16); xlabel('time (s)'); 
%ylabel('pr_F'); ylim([0 0.6]);
ylabel('<R^2>'); ylim([0 0.2]);
legend('PFC','HP'); legend boxoff
title('Pop-Decoding');


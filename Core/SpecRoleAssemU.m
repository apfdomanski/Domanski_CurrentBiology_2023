% Plot the specific role of Assembly-member units in decoding behavioural out come
%
% (c) 2014, Daniel Durstewitz, Dept. Theoretical Neuroscience, CIMH
% Mannheim/ Heidelberg University

clear all
close all


pat='data/'; 
% A=dir([pat '*TAUpairs10.mat']);
A=dir('*.mat')
avgIFin=zeros(length(A),4); avgIFout=avgIFin; avgTSin=avgIFin; avgTSout=avgIFin;
rIF_TS=zeros(length(A),2); pIF_TS=rIF_TS; TSall=cell(1,2); IFall=cell(1,2);
for f=1:length(A)
%     fn=[pat A(f).name]; load(fn);
    k=min(strfind(fn,'_'));
    if str2num(fn(k+2))==2, ext='L2'; else ext=''; end;
    fnF=[fn(1:k) 'iFR50' ext '_FSC']; load(fnF);
    fnFR=[fn(1:k) 'PFC_iFR50' ext];
    
    % Decoding scores for all sgl-units
    twin=10; minFR=0.1; critCvBW=1e6; bw=0.05; reg=0.05;
    [TmtxS,iFR0,EvtLs,EvtTs]=SelTrialsCells(fnFR,twin,minFR,critCvBW,'iFR');
    ko=round((twin+5)/bw);
    iFRsel=cell(1,2);
    for s=1:2
        for i=1:length(TmtxS{s})
            TmtxS{s}{i}=TmtxS{s}{i}([1:ko end-ko+1:end]);
            iFR0{s}{i}=iFR0{s}{i}([1:ko end-ko+1:end],:)';
        end;
        iFRsel{s}=cell2mat(iFR0{s})';
    end;
    ntr=length(EvtTs)/2;
    evt0=EvtLs(((1:ntr)-1)*2+1)-2;
    n1=length(find(evt0==1)); n2=length(find(evt0==2)); TS=cell(1,2);
    for s=1:2, [~,~,~,~,~,~,TS{s}]=DecodeStats(iFRsel{s},evt0,reg); end;
    
    % IF's & D's for assembly vs. non-assembly units
    rs=cell(1,2);
    w=repmat(1:nu(1),nu(2),1); rs{1}=w(1:end)';
    rs{2}=repmat((1:nu(2))',nu(1),1);
    for s=1:2
        U=unique(cell2mat(units{s}));
        k=find(ismember(rs{s},U)); q=IF{s}(k,:); avgIFin(f,s)=mean(q(~isnan(q)));
        k=find(~ismember(rs{s},U)); q=IF{s}(k,:); avgIFout(f,s)=mean(q(~isnan(q)));
        q=TS{s}(:,U); avgTSin(f,s)=mean(q(~isnan(q)));
        q=TS{s}(:,setdiff(1:nu(s),U)); avgTSout(f,s)=mean(q(~isnan(q)));
    end;
    y=cellfun(@(x)(x>nu(1)),units{3},'UniformOutput',0);
    w=cellfun(@sum,y); L=cellfun(@length,y);
    R=find(w./L>0 & w./L<1);
    ks=[]; ku={[],[]};
    for r=1:length(R)
        m1=units{3}{R(r)}<=nu(1); m2=units{3}{R(r)}>nu(1);
        ks=union(ks,find(ismember(rs{1},units{3}{R(r)}(m1)) & ...
            ismember(rs{2},units{3}{R(r)}(m2)-nu(1))));
        ku{1}=union(ku{1},units{3}{R(r)}(m1));
        ku{2}=union(ku{2},units{3}{R(r)}(m2)-nu(1));
    end;
    kn=setdiff(1:nu(1)*nu(2),ks);
    for s=1:2
        q=IF{s}(ks,:); avgIFin(f,s+2)=mean(q(~isnan(q)));
        q=IF{s}(kn,:); avgIFout(f,s+2)=mean(q(~isnan(q)));
        q=TS{s}(:,ku{s}); avgTSin(f,s+2)=mean(q(~isnan(q)));
        kk=setdiff(1:nu(s),ku{s});
        q=TS{s}(:,kk); avgTSout(f,s+2)=mean(q(~isnan(q)));
    end;
    
end;
save([pat 'SpecAssemRole'],'avgIFin','avgIFout','avgTSin','avgTSout','rIF_TS','pIF_TS');


load([pat 'SpecAssemRole']);
xIN=[]; xOUT=[];
for s=1:4
    kk=find(~isnan(avgIFin(:,s)) & ~isnan(avgIFout(:,s)));
    disp((avgIFin(kk,s)>avgIFout(kk,s))')
    [h,p]=ttest(avgIFin(kk,s),avgIFout(kk,s)); disp(p)
    if s>2, xIN=[xIN avgIFin(kk,s)']; xOUT=[xOUT avgIFout(kk,s)']; end;
end;
yIN=[]; yOUT=[];
for s=1:4
    kk=find(~isnan(avgTSin(:,s)) & ~isnan(avgTSout(:,s)));
    disp((avgTSin(kk,s)>avgTSout(kk,s))')
    [h,p]=ttest(avgTSin(kk,s),avgTSout(kk,s)); disp(p)
    if s>2, yIN=[yIN avgTSin(kk,s)']; yOUT=[yOUT avgTSout(kk,s)']; end;
end;
subplot(1,2,1), hold off cla
plot([1 2]'*ones(size(xIN)),[xIN;xOUT],'bo-','LineWidth',2); xlim([0 3]);
[h,p]=ttest(xIN',xOUT'); xlabel(['p <= ' num2str(p)]); box off
title('GC'); set(gca,'FontSize',20,'XTick',[1 2],'XTickLabel',{'assem','non-assem'});
subplot(1,2,2), hold off cla
plot([1 2]'*ones(size(yIN)),[yIN;yOUT],'bo-','LineWidth',2); xlim([0 3]);
[h,p]=ttest(yIN',yOUT'); xlabel(['p <= ' num2str(p)]); box off
title('Decoding'); set(gca,'FontSize',20,'XTick',[1 2],'XTickLabel',{'assem','non-assem'});


% (c) 2014, Daniel Durstewitz, Dept. Theoretical Neuroscience, CIMH
% Mannheim/ Heidelberg University

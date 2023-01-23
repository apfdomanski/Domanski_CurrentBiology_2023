% Plot individual/population assembly decoding performance
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk

clear all; close all
% start in directory holding decoding results
bw=0.05; iFRtyp='iFR';
load(['ASSEMPopDecodeStats_' iFRtyp num2str(bw*1e3)]);
% load([pat 'ErrTrDecodeStats_' iFRtyp num2str(bw*1e3)]);
Font_Size = 12;

%which files are LONG2 recordings?
A=dir('*LONG*');
for f_id=1:length(A)
       B= strfind(A(f_id).name,'LONG2');
       if ~isempty(B), C(f_id)=f_id; else C(f_id)=NaN;end
end
C(isnan(C))=[];
%... or, use all recordings 
C=1:length(A);

%%
f=6;
CL=zeros(1,length(Tv));
k=find(Tv>=9 & Tv<=11);  CL(k)=1;
k=find(Tv>=19 & Tv<=21); CL(k)=2;
X=[TS{1}{f} TS{2}{f} TS{3}{f}]; mxT=max(X(1:end));
mxF=max([Ft2{1}(f,:) Ft2{2}(f,:) Ft2{3}(f,:)]);
r=round(length(Tv)/2);
X1=1:r; X2=r+1:2*r;
tt={'mPFC','Hipp', 'Joint mPFC-Hipp'};
figure
for s=1:3
    % Assembly decoding
    subplot(3,3,s)%, hold off cla
    %figure(1), hold off cla
    [~,i1]=find(TS{s}{f}>1.5); tsel=unique(i1);
    TS0=TS{s}{f}(:,tsel);
    [~,k]=max(TS0); [~,r]=sort(k);
    TS0(1,1)=mxT;   % just to ensure visualization with common max
    [evtB,Clr]=ShowMASSEM(Tv,TS0(:,r),CL);
    xlim([5 25]); title(tt{s},'FontSize',Font_Size);
end
for s=1:3
    % pop. decoding
    subplot(3,3,s+3)
    %figure(2), hold off cla
    Y1=[Ft2ciL{s}(f,X1);Ft2ciH{s}(f,X1)-Ft2ciL{s}(f,X1)];
    Y2=[Ft2ciL{s}(f,X2);Ft2ciH{s}(f,X2)-Ft2ciL{s}(f,X2)];
    h=area(Tv(X1),Y1'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','y'); hold on
    h=area(Tv(X2),Y2'); set(h(1),'FaceColor','w'); set(h(2),'FaceColor','y'); hold on
    %plot(Tv(X1),Ft2{s}(f,X1),'ko-',Tv(X2),Ft2{s}(f,X2),'ko-','LineWidth',2);
    plot(Tv(X1),Ft2{s}(f,X1),'k',Tv(X2),Ft2{s}(f,X2),'k','LineWidth',2);
    axis([5 25 0 1.1*mxF]); box off
%     axis([5 25 0 plot_max]); box off

    for i=1:2:3
        k0=evtB(i)+1; k1=evtB(i+1); c=CL(k0);
        plot([Tv(k0)-bw/2 Tv(k0)-bw/2],[0 1.1*mxF],'Color',Clr{c},'LineWidth',2);
        plot([Tv(k1)+bw/2 Tv(k1)+bw/2],[0 1.1*mxF],'Color',Clr{c},'LineWidth',2);
    end;
    plot([mean(Tv) mean(Tv)],[0 1.1*mxF],'k--','LineWidth',2);
    set(gca,'FontSize',Font_Size);
    ylabel('Decode score (F-value)'); xlabel('Time (sec)');
end;
%% Plot correct decoding pop average
subplot(3,2,[5 6]), hold off cla
plot(Tv,ones(size(Tv)),'Color',[0.7 0.7 0.7],'LineWidth',3); hold on

% normalize F-value to baseline
bp=find(Tv>=5 & Tv<=7);  
Ft2nrm=cell(1,3);
for s=1:3
    B=mean(Ft2{s}(:,bp)');
    Ft2nrm{s}=Ft2{s}./(B'*ones(1,length(Tv)));  
    Ft2nrm{s}=Ft2nrm{s}(C,:);
end;

U1=nanmean(Ft2nrm{1}); U2=nanmean(Ft2nrm{2}); U3=nanmean(Ft2nrm{3});
% N=sum(~isnan(sum(Ft2nrm{1},2))); SE1=std(Ft2nrm{1})./sqrt(N); 
SE1=nansem(Ft2nrm{1}); 
SE2=nansem(Ft2nrm{2});
SE3=nansem(Ft2nrm{3});
%miZ=-1.1*abs(min([U1 U2])); mxZ=1.1*max([U1 U2]);
%plot(Tv,U1,'r',Tv,U2,'b','LineWidth',2); hold on
miZ=-1.5*abs(min([U1 U2 U3])); mxZ=1.5*max([U1 U2 U3]);
% errorbar(Tv,U1,SE1,'ro-','LineWidth',2); hold on
% errorbar(Tv,U2,SE2,'bo-','LineWidth',2);
cmp=colormap('jet');
hold on
h=area([Tv;Tv]',[U1-SE1;2*SE1]');
set(h(1),'FaceColor','none','EdgeColor','none');
set(h(2),'FaceColor',[1 0.6 0.6],'EdgeColor','none');%set(h(2),'FaceColor',cmp(51,:),'EdgeColor','none');

h=area([Tv;Tv]',[U2-SE2;2*SE2]');
set(h(1),'FaceColor','none','EdgeColor','none');
set(h(2),'FaceColor',[0.6 0.6 1],'EdgeColor','none');%set(h(2),'FaceColor',cmp(22,:),'EdgeColor','none');

h=area([Tv;Tv]',[U3-SE3;2*SE3]');
set(h(1),'FaceColor','none','EdgeColor','none');
set(h(2),'FaceColor',[0.6 1 0.6],'EdgeColor','none');%set(h(2),'FaceColor',cmp(22,:),'EdgeColor','none');


plot(Tv,U1,'r',Tv,U2,'b',Tv,U3,'g','LineWidth',2);
axis([5 25 miZ mxZ]); box off

% axis([5 25 0 plot_max]); box off
for i=1:2:3
    k0=evtB(i)+1; k1=evtB(i+1); c=CL(k0);
    plot([Tv(k0)-bw/2 Tv(k0)-bw/2],[miZ mxZ],'Color',Clr{c},'LineWidth',2);
    plot([Tv(k1)+bw/2 Tv(k1)+bw/2],[miZ mxZ],'Color',Clr{c},'LineWidth',2);
end;
plot([mean(Tv) mean(Tv)],[miZ mxZ],'k--','LineWidth',2);
set(gca,'FontSize',Font_Size);
ylabel('Norm. decode score'); xlabel('Time (sec)');
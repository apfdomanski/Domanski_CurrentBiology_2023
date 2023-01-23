clear all
close all
target= 'LONG';
if ispc
    home = 'C:\'; % [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    pat = [home 'Analysis\AssemblyAnalysis\raw\KDE_bins\'];
elseif ismac
    path(path,'/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/MichalDataAna')
    home = getenv('HOME');
    cd ([home '/Documents/Bristol_work/assembly_analysis/MichalDataAna'])
    pat = [home '/Documents/Bristol_work/assembly_analysis/Michal_data/raw/KDE_bins/'];
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/5KDS_bins/'];    
end

if ~isdir([pat '/' target])
    mkdir([pat '/' target])
end
% run on 1) iFRsc, 2) KDE iFR

%bw=0.01; iFRtyp='iFRsc';
%bw=0.05; iFRtyp='iFRsc';

bw=0.05; iFRtyp='iFR';
A=dir([pat strcat('*', target) '*PFC_iFR' num2str(bw*1e3) '*.mat']);

for event_id=1:length(A)
    k=findstr(A(event_id).name,'_'); r=findstr(A(event_id).name,'.');
    Name{event_id}=[A(event_id).name(1:k(1)) A(event_id).name(k(2)+1:r-1)];
end;
twin     = 10;   % Time window
minFR    = 0.1;  % minimal acceptable firing rate
critCvBW = 1e6;  % critical max bound on kernel width drift over time (variance)
alpha    = 0.01; % significance threshold
kmax     = 30;   % number of factors to check up to
Nbs      = 500;  % number of bootstrap draws

%% run Factor analysis

for f=1:length(A)
    
    %%% Prepare output vectors
    nassem=cell(1,3); 
    FL=nassem; 
    FSC=nassem; LL=FL; AIC=LL; BIC=LL; psixBS=LL; perm=LL;
    Pr=LL; Chi2=LL; prBS=LL; LLbs=LL; Chi2bs=LL; nuBS=LL; psix=LL; FLbs=LL; FSCbs=LL;
    %%%%
    
    %%%% Prepare data
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fn=[pat A(f).name];
    [TmtxS,SCM,~,~,usel_out]=SelTrialsCells(fn,twin,minFR,critCvBW,'iFR'); %'iFRsc'
    bw=mean(diff(TmtxS{1}{1})); % get bin-width
    ko=round((twin+5)/bw);      % no points surrounding events
    for area_id=1:2  % ...loop across all brain areas
        for event_id=1:length(TmtxS{area_id}) % ...loop across all events
            TmtxS{area_id}{event_id} = TmtxS{area_id}{event_id}([1:ko end-ko+1:end]);
            SCM{area_id}{event_id}=SCM{area_id}{event_id}([1:ko end-ko+1:end],:)';
        end;
    end;
    %%%%
    
    %%%% Run Assembly detection
    % (Intra area): Run FA-based assembly detection algo separately for each area
    for s=1:2
        [nassem{s},FL{s},~,LL{s},AIC{s},BIC{s},Pr{s},Chi2{s}, ...
            prBS{s},LLbs{s},Chi2bs{s},psix{s},FLbs{s},~,psixBS{s},perm{s}]  =   FAassem_(SCM{s},kmax,Nbs,alpha);
    end;
    % (Inter area): Combine PFC and HP, run Factor analysis again
        s=3; SCMcomb=cell(1,length(SCM{1}));
    for event_id=1:length(SCM{1}), SCMcomb{event_id}=[SCM{1}{event_id};SCM{2}{event_id}]; end;
        [nassem{s},FL{s},~,LL{s},AIC{s},BIC{s},Pr{s},Chi2{s}, ...
         prBS{s},LLbs{s},Chi2bs{s},psix{s},FLbs{s},~,psixBS{s},perm{s}]      =   FAassem_(SCMcomb,kmax,Nbs,alpha);
    %%%% 
    
    %%%% Save this recording
    fnOut=[pat target '/' Name{f} '_AssemRes2'];
    save(fnOut,'TmtxS','nassem','FL','LL','AIC','BIC', ...
        'Pr','prBS','LLbs','FLbs','Chi2bs','Chi2','psix','psixBS','perm','usel_out');
    %%%%
    
end;

%% compute factor scores and BS confidence limits

for f=1:length(A)
    
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fn=[pat A(f).name];
    [TmtxS,SCM,~,~,usel_out]=SelTrialsCells(fn,twin,minFR,critCvBW,'iFR'); %'iFRsc'
    
    fnOut=[pat target '/' Name{f} '_AssemRes2'];
    [FSC,FSCbs,ciLld,ciHld,ciLsc,ciHsc]=BSci(TmtxS,SCM,fnOut,twin);
    
    k=findstr(fnOut,'_');
    save([fnOut(1:k(end)) '_FSCtemp']','FSC','FSCbs','ciLld','ciHld','ciLsc','ciHsc','usel_out');

end;
%% extract assemblies (choose BS as reasonable but relatively conservative)
% extract assembly activations as f(t) through factor scores
for f=1:length(A)
    fnOut=[pat target '/' Name{f} '_AssemRes2']; load(fnOut);
    fprintf([Name{f}, ' AssemRes2: %G, %G, %G, factors.\n'],nassem{1}(3),nassem{2}(3),nassem{3}(3))
    fnOut=[pat target '/' Name{f} '__FSCtemp']; load(fnOut);    
    units=cell(1,length(nassem)); FSCsel=units; Tmtx=FSCsel;
    for s=1:length(nassem)
        
        numf=nassem{s}(3);
        for n=1:numf, 
            units{s}{n}=find(abs(FL{s}{numf}(:,n))>ciHld(s,3))'; 
        end;
        % Weed out factors with only one unit
        k=find(cellfun(@length,units{s})>1);
        units{s}=units{s}(k);
        FSCsel{s}=FSC{s}(:,k);
        
    end;
    
    fn=[pat A(f).name];
    [TmtxS,SCM,~,~,usel_out]=SelTrialsCells(fn,twin,minFR,critCvBW,'iFR'); %'iFRsc'
    
    bw=mean(diff(TmtxS{1}{1})); ko=round((twin+5)/bw);
    
    nu=[]; 
    for s=1:2, 
        nu(s)=size(SCM{s}{1},2); 
    end;
    for event_id=1:length(TmtxS{1}),
        TmtxS{1}{event_id} = TmtxS{1}{event_id}([1:ko end-ko+1:end]);
    end;
    Tmtx=cell2mat(TmtxS{1})';
    k=findstr(fnOut,'_');
    save([fnOut(1:k(end-1)) 'FSC'],'units','FSCsel','Tmtx','Name','nassem', ...
        'ciLld','ciHld','ciLsc','ciHsc','nu');
end;
for f=1:length(A)
    
    fnOut=[pat target '/' Name{f} '_FSC']; load(fnOut);    fn=[pat A(f).name];
    [TmtxS,SCM,EvtLs,EvtTs]=SelTrialsCells(fn,twin,minFR,critCvBW,'iFR'); %'iFRsc'
    save(fnOut,'units','FSCsel','Tmtx','Name','nassem', ...
        'ciLld','ciHld','ciLsc','ciHsc','nu','EvtTs','EvtLs');
end;
%% Check factor no.s line up

for f=1:length(A)
    clear FSC
    fnOut=[pat target '/' Name{f} '_AssemRes2']; load(fnOut,'nassem');
    fprintf([Name{f}, ' AssemRes2: %G, %G, %G, factors.\n'],nassem{1}(3),nassem{2}(3),nassem{3}(3))
    fnOut=[pat target '/' Name{f} '_FSC']; load(fnOut,'FSCsel');    
	fprintf([Name{f}, ' _FSC: %:G, %G, %G, factors.\n'],nassem{1}(3),nassem{2}(3),nassem{3}(3))
    fprintf([Name{f}, ' _FSCCsel %G %G %G\n\n'],size(FSCsel{1},2),size(FSCsel{2},2),size(FSCsel{3},2))
end

%% visualize results
for f=1:length(A)
    figure
    fnOut=[pat target '/' Name{f} '_AssemRes2']; load(fnOut);
    fnOut=[pat target '/' Name{f} '__FSCtemp']; load(fnOut);
    for s=1:length(FSC)
        z=diff(LL{s});
        Zbs=diff(LLbs{s}')';
        x=2:length(z)+1;
        subplot(1,2,1), hold off cla, plot(x,z,'b','LineWidth',2);
        %hold on, plot(Zbs','r');
        r=round(length(Zbs)*alpha); zz=sort(Zbs,'descend');
        hold on, plot(x,zz(r)*ones(1,length(z)),'r--','LineWidth',2);
        set(gca,'FontSize',22), box off
        title(['(' num2str(f) ',' num2str(s) '): ' num2str(nassem{s})]);
        xlabel('Factor (# assemblies)'), ylabel('log-likelihood-ratio')
        legend('original','trial-permutation'), legend('boxoff')
        subplot(1,2,2), hold off cla
        vs=cell2mat(FLbs{s}{2});
        x=-0.5:0.005:0.5; h=histc(vs(1:end),x)./length(vs(1:end));
        plot(x,h,'b',[ciHld(s,2) ciHld(s,2)],[0 0.05],'g','LineWidth',2);
        hold on, plot(FL{s}{2}(1:end),0.01,'r.'); xlim([min(x) max(x)]);
        set(gca,'FontSize',22), box off
        
%         subplot(1,3,3), hold off cla
%         ws=sort(FSCbs{s}(1:end),'ascend');
%         ciV=[0.1 0.05 0.01 5e-3 1e-3];
%         for i=1:length(ciV)
%             ciLsc(s,i)=ws(round(length(ws)*ciV(i)));
%             ciHsc(s,i)=ws(round(length(ws)*(1-ciV(i))));
%         end;
%         x=-3:0.1:3; h=histc(ws,x)./length(ws);
%         plot(x,h,'b',[ciHsc(s,2) ciHsc(s,2)],[0 0.05],'g','LineWidth',2);
%         hold on, plot(FSC{s}(1:end),0.01,'r.'); xlim([min(x) max(x)]);
%         set(gca,'FontSize',22), box off
%         
%         a=input('...');
    end;
end;
%% visualize factor scores/ activations in time & distribution
for f=1:length(A)
    fnOut=[pat target '/' Name{f} '_AssemRes2']; load(fnOut);
    fnOut=[pat target '/' Name{f} '__FSCtemp']; load(fnOut);
    for s=1:length(FSC)
        [ntp,nu]=size(FSC{s});
        for u=1:nu
            [s u]
            figure
            n=find(abs(FL{s}{nassem{s}(3)}(:,u))>ciHld(s,3));
            subplot(1,2,1), hold off cla
            t=(1:ntp)*0.05+0.05/2;
            plot(t,FSC{s}(:,u),'b',t,ciHsc(s,2),'r-','LineWidth',2);
            set(gca,'FontSize',22), box off, xlim([680 700]);
            title(['#units: ' num2str(n')]);
            subplot(1,2,2), hold off cla
            x=-3:0.01:3; h=histc(FSC{s}(:,u),x)./length(FSC{s}(:,u));
            plot(x,h,'b',[ciHsc(s,2) ciHsc(s,2)],[0 0.05],'g','LineWidth',2);
            xlim([-0.5 1]);
            set(gca,'FontSize',22), box off
            title(['#units: ' num2str(n')]);
%             a=input([num2str([f s u]) '...']);
        end;
    end;
end;
%% visualize average activity of PFC, HP, & joint assemblies as func. of task phase
% Tinv=[2 8; 9 11; 13 17; 19 21; 24 30];
% load([pat target '/' 'CVEres'])
% bw=0.05;
% cmp=colormap('jet');
% clr1={cmp(18,:),cmp(33,:),cmp(50,:)};
% clr2={'b','g','r'};
% Ya=cell(1,5); for t=1:5, Ya{t}=cell(1,3); end;
% for f=[1]%1:length(A)
%     figure(1), hold off cla
%     fnOut=['data/' Name{f} '_FSC']; load(fnOut);
%     ntr=length(EvtTs)/2; Ltr=round(size(FSCsel{1},1)/ntr);
%     evt0=EvtLs(((1:ntr)-1)*2+1)-2;
%     Tv=(1:Ltr)*bw-bw/2;
%     na=size(FSCsel{3},2);
%     Xavg=zeros(na,Ltr);
%     for i=1:na
%         A=reshape(FSCsel{3}(:,i),Ltr,ntr)';
%         Xavg(i,:)= mean(A);
%         plot(Tv,Xavg(i,:),'Color',clr1{AssemType{f}(i)+1},'LineStyle','--','LineWidth',0.1); hold on
% %         for t=1:size(Tinv,1)
% %             k=find(Tv>Tinv(t,1) & Tv<=Tinv(t,2));
% %             Ya{t}{AssemType{f}(i)+1}=[Ya{t}{AssemType{f}(i)+1} mean(Xavg(i,k))];
% %         end;
%     end;
%     for s=unique(AssemType{f})
%         k=find(AssemType{f}==s);
%        %%%%%%%%%%%%%%%%%%%% AD replaced k with s+1
%        % plot(Tv,mean(Xavg(k,:),1),'Color',clr2{s+1},'LineWidth',2.5); hold on
%         plot(Tv,mean(Xavg(s+1,:),1),'Color',clr2{s+1},'LineWidth',2.5); hold on 
% 
%     end;
%     %a=input([num2str(f) '...']);
% end;
% for t=1:5, for s=1:3
%         YaAvg(t,s)=mean(Ya{t}{s}); YaSE(t,s)=std(Ya{t}{s})/sqrt(length(Ya{t}{s}));
%     end; end;
% % hold off cla, 
% errorbar(YaAvg,YaSE,'o-','LineWidth',2), xlim([0 6]),set(gca,'FontSize',20)
%     
%     legend('PFC','HP','PFC-HP'), legend('boxoff'),set(gca,'FontSize',20)
%     set(gca, 'XTick',1:5,'XTickLabel',{'pre-S','S','Delay','C','post-C'});
%     ylabel('read-out (avg. {\itt}-score)')
%     text(1.85,1.9,'*','FontSize',42), text(2.85,1.9,'*','FontSize',42), text(3.85,1.9,'*','FontSize',42)
% 
% for t=1:5
%     [~,p1]=ttest2(Ya{t}{1},Ya{t}{3});
%     [~,p2]=ttest2(Ya{t}{2},Ya{t}{3});
%     [p1 p2]
% end;


% % *** overlap among assemblies
% load('data/FAassem');
% Ovlp=cell(size(units));
% Ov=cell(2,1);
% for i=1:size(units,1)
%     for j=1:size(units,2)
%         k=length(units{i,j});
%         Ovlp{i,j}=ones(k);
%         for n1=1:k, for n2=1:k,
%                 Ovlp{i,j}(n1,n2)=length(intersect(units{i,j}{n1},units{i,j}{n2}))/length(units{i,j}{n1});
%         end; end;
%         if k>1
%             y=squareform(triu(Ovlp{i,j},1)+triu(Ovlp{i,j},1)');
%             Ov{j}(end+1:end+length(y))=y;
%         end;
%     end;
% end;


% (c) 2014, Daniel Durstewitz, Dept. Theoretical Neuroscience, CIMH
% Mannheim/ Heidelberg University

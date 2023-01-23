clear all
close all
target= 'physostigmine';
if ispc
    home = 'C:\'; % [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    pat = [home 'Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\'];
    pat_ = [home 'Analysis\AssemblyAnalysis\raw\'];
elseif ismac
%     path(path,'/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/MichalDataAna')
    
%     pat_ = '/Volumes/HDD2/DNMTP/raw/';
%     pat = '/Volumes/HDD2/DNMTP/raw/KDE_binsTaskonly/';
    
%     pat = '/Volumes/Aleks Data Portable/AssemblyAnalysis/raw/';
%     pat_ = '/Volumes/Aleks Data Portable/AssemblyAnalysis/raw/KDE_binsTaskonly/';

    pat_ = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw/';
    pat = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw/KDE_binsTaskonly/';

else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/KDE_binsTaskonly/'];   
    pat_ = [home '/raw/'];   
end

if ~isdir([pat target  'TaskOnly'])
    mkdir([pat target 'TaskOnly'])
end

bw=0.05; iFRtyp='iFR';
A=dir([pat strcat('*', target) '*PFC_iFR' num2str(bw*1e3) '*.mat']);

 reject_list={'._'}; 
    name_flag=zeros(numel(A),1);
    for idx=1:numel(A)
        fnames_{idx,1}=A(idx).name;
        name_flag(idx,1) = cellfun(@(s) ~isempty(strfind(A(idx).name, s)), reject_list);
%         name_flag(idx,1)= logical(sum(ismember(reject_list,A(idx).name))); %AD inverted logical flag for testing
    end
    try
        A(find(name_flag))=[];% fnames_(name_flag)=[];
    end
    clear  reject_list idx fnames_ name_flag


for event_id=1:length(A)
    k=findstr(A(event_id).name,'_'); r=findstr(A(event_id).name,'.');
    Name{event_id}=[A(event_id).name(1:k(1)) A(event_id).name(k(2)+1:r-1)];
end;
twin     = 10;   % Time window
minFR    = 0.1;  % minimal acceptable firing rate
critCvBW = 1e6;  % critical max bound on kernel width drift over time (variance)
alpha    = 0.01; % significance threshold
kmax     = 30;   % number of factors to check up to
Nbs      = 10;  % number of bootstrap draws
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
    fname=strtok(A(f).name,'_');

    load(sprintf('%sallTimestamps%s%s_Events.mat',pat_,filesep,fname));
    tlims_ = [...
        t.Short.CueLight_LeftCorrect,t.Short.ChoicePress_LeftCorrect;...
        t.Short.CueLight_RightCorrect,t.Short.ChoicePress_RightCorrect;...
        t.Short.CueLight_LeftError',t.Short.ChoicePress_LeftError';...
        t.Short.CueLight_RightError',t.Short.ChoicePress_RightError';...
        t.Medium.CueLight_LeftCorrect,t.Medium.ChoicePress_LeftCorrect;...
        t.Medium.CueLight_RightCorrect,t.Medium.ChoicePress_RightCorrect;...
        t.Medium.CueLight_LeftError',t.Medium.ChoicePress_LeftError';...
        t.Medium.CueLight_RightError',t.Medium.ChoicePress_RightError';...
        t.Long.CueLight_LeftCorrect,t.Long.ChoicePress_LeftCorrect;...
        t.Long.CueLight_RightCorrect,t.Long.ChoicePress_RightCorrect;...
        t.Long.CueLight_LeftError',t.Long.ChoicePress_LeftError';...
        t.Long.CueLight_RightError',t.Long.ChoicePress_RightError']*1e-6;
    
    tlims_ = bsxfun(@plus,tlims_,[-twin twin]);
    

    fn=[pat A(f).name];
    [TmtxS,SCM,usel_out]=SelTrialsCellsWholeTrial(fn,tlims_,minFR,critCvBW); %'iFRsc'
    bw=mean(diff(TmtxS{1}{1})); % get bin-width
%     for s=1:2  % ...loop across all brain areas
%        TmtxS{s}=cell2mat(TmtxS{s});
%        SCM{s}=cell2mat(SCM{s}');
%     end
    %%%%
    
    %%%% Run Assembly detection
    % (Intra area): Run FA-based assembly detection algo separately for each area
    for s=1:2
        try
            [nassem{s},FL{s},~,LL{s},AIC{s},BIC{s},Pr{s},Chi2{s}, ...
                prBS{s},LLbs{s},Chi2bs{s},psix{s},FLbs{s},~,psixBS{s},perm{s}]  =   FAassem_(SCM{s},kmax,Nbs,alpha);
        end
    end
    % (Inter area): Combine PFC and HP, run Factor analysis again
    try
        s=3; SCMcomb=cell(1,length(SCM{1}));
        for event_id=1:length(SCM{1}), SCMcomb{event_id}=[SCM{1}{event_id};SCM{2}{event_id}]; end;
        [nassem{s},FL{s},~,LL{s},AIC{s},BIC{s},Pr{s},Chi2{s}, ...
            prBS{s},LLbs{s},Chi2bs{s},psix{s},FLbs{s},~,psixBS{s},perm{s}]      =   FAassem_(SCMcomb,kmax,Nbs,alpha);
    end
    %%%% 
    
    %%%% Save this recording
    fnOut=[pat target  'Taskonly' filesep Name{f} '_AssemRes2'];
    save(fnOut,'TmtxS','nassem','FL','LL','AIC','BIC', ...
        'Pr','prBS','LLbs','FLbs','Chi2bs','Chi2','psix','psixBS','perm','usel_out');
    %%%%
    
end
%% compute factor scores and BS confidence limits
 
for f=1:length(A)
    
     %%%% Prepare data
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fname=strtok(A(f).name,'_');

    load(sprintf('%sallTimestamps%s%s_Events.mat',pat_,filesep,fname));
    tlims_ = [...
        t.Short.CueLight_LeftCorrect,t.Short.ChoicePress_LeftCorrect;...
        t.Short.CueLight_RightCorrect,t.Short.ChoicePress_RightCorrect;...
        t.Short.CueLight_LeftError',t.Short.ChoicePress_LeftError';...
        t.Short.CueLight_RightError',t.Short.ChoicePress_RightError';...
        t.Medium.CueLight_LeftCorrect,t.Medium.ChoicePress_LeftCorrect;...
        t.Medium.CueLight_RightCorrect,t.Medium.ChoicePress_RightCorrect;...
        t.Medium.CueLight_LeftError',t.Medium.ChoicePress_LeftError';...
        t.Medium.CueLight_RightError',t.Medium.ChoicePress_RightError';...
        t.Long.CueLight_LeftCorrect,t.Long.ChoicePress_LeftCorrect;...
        t.Long.CueLight_RightCorrect,t.Long.ChoicePress_RightCorrect;...
        t.Long.CueLight_LeftError',t.Long.ChoicePress_LeftError';...
        t.Long.CueLight_RightError',t.Long.ChoicePress_RightError']*1e-6;
    tlims_ = tlims_ + repmat([-twin twin],size(tlims_,1),1);
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fn=[pat A(f).name];
    [TmtxS,SCM,usel_out]=SelTrialsCellsWholeTrial(fn,tlims_,minFR,critCvBW); %'iFRsc'
    
    fnOut=[pat target  'Taskonly' filesep Name{f} '_AssemRes2'];
    
        [FSC,FSCbs,ciLld,ciHld,ciLsc,ciHsc]=BSciTaskOnly(TmtxS,SCM,fnOut,twin);   
    
    k=findstr(fnOut,'_');
    save([fnOut(1:k(end)) '_FSCtemp']','FSC','FSCbs','ciLld','ciHld','ciLsc','ciHsc','usel_out');

end;
%% extract assemblies (choose BS as reasonable but relatively conservative)
%extract assembly activations as f(t) through factor scores
for f=1:length(A)
     %%%% Prepare data
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fname=strtok(A(f).name,'_');

    load(sprintf('%sallTimestamps%s%s_Events.mat',pat_,filesep,fname));
    tlims_ = [...
        t.Short.CueLight_LeftCorrect,t.Short.ChoicePress_LeftCorrect;...
        t.Short.CueLight_RightCorrect,t.Short.ChoicePress_RightCorrect;...
        t.Short.CueLight_LeftError',t.Short.ChoicePress_LeftError';...
        t.Short.CueLight_RightError',t.Short.ChoicePress_RightError';...
        t.Medium.CueLight_LeftCorrect,t.Medium.ChoicePress_LeftCorrect;...
        t.Medium.CueLight_RightCorrect,t.Medium.ChoicePress_RightCorrect;...
        t.Medium.CueLight_LeftError',t.Medium.ChoicePress_LeftError';...
        t.Medium.CueLight_RightError',t.Medium.ChoicePress_RightError';...
        t.Long.CueLight_LeftCorrect,t.Long.ChoicePress_LeftCorrect;...
        t.Long.CueLight_RightCorrect,t.Long.ChoicePress_RightCorrect;...
        t.Long.CueLight_LeftError',t.Long.ChoicePress_LeftError';...
        t.Long.CueLight_RightError',t.Long.ChoicePress_RightError']*1e-6;
    tlims_ = tlims_ + repmat([-twin twin],size(tlims_,1),1);
    
    nu=[];
    ciLld=[];
    ciHld=[];
    ciLsc=[];
    ciHsc=[];
    nassem = cell(1,3);
    units  = cell(1,length(nassem)); 
    FSCsel = cell(1,length(nassem)); 
    Tmtx   = cell(1,length(nassem)); 
    
    
    fnOut=[pat target 'Taskonly' filesep Name{f} ]; 
    disp(sprintf('Processing: %s [...] \n',fnOut));
    
    load([fnOut '_AssemRes2']);
    load([fnOut '__FSCtemp']);
    
thresh_ = 3;
    for s=1:length(nassem)
            if ~isempty(nassem{s})
                numf=nassem{s}(thresh_);
                for n=1:numf
                    units{s}{n}=find(abs(FL{s}{numf}(:,n))>ciHld(s,thresh_))';
%                     units{s}{n}=find(abs(FL{s}{numf}(:,n))>0.01)';
                end
                % Weed out factors with only one unit
                k=find(cellfun(@length,units{s})>1);
                units{s}=units{s}(k);
                FSCsel{s}=FSC{s}(:,k);
            else
                FSCsel{s} = []; units{s};
            end
    end
   
    % how many units went in?
    fn=[pat A(f).name];
    [TmtxS,SCM,usel_out]=SelTrialsCellsWholeTrial(fn,tlims_,minFR,critCvBW); %'iFRsc'
    nu=[]; 
    for s=1:2
        nu(s)=size(SCM{s}{1},1); 
    end
    
    % Make sure that inter-area assems truely span the two areas
    if ~isempty(FSCsel{3})
        %numf=nassem{length(nassem)}(3);
        numf=size(FSCsel{3},2);
        for i=1:numf
            units_ = units{3}{i};
            k(i) = sum(ismember(units_,1:nu(1)))>0 & sum(ismember(units_,(nu(1)+1):sum(nu)))>0;
        end
        units{3}(~k)=[];
        FSCsel{3}(:,~k)=[];
    end
    try
        Tmtx=cell2mat(TmtxS{1}');
    catch
        Tmtx=cell2mat(TmtxS{2}');
    end

    save([fnOut '_FSC'],'units','FSCsel','Tmtx','Name','nassem', ...
        'ciLld','ciHld','ciLsc','ciHsc','nu');
end
%% visualize results
% for f=1:length(A)
%     figure
%     fnOut=[pat target 'Taskonly' filesep Name{f} '_AssemRes2']; load(fnOut);    
%     fnOut=[pat target 'Taskonly' filesep Name{f} '__FSCtemp']; load(fnOut);    
% 
%   
%     for s=1:length(FSC)
%         z=diff(LL{s});
%         Zbs=diff(LLbs{s}')';
%         x=2:length(z)+1;
%         subplot(1,2,1), hold off cla, plot(x,z,'b','LineWidth',2);
%         %hold on, plot(Zbs','r');
%         r=round(length(Zbs)*alpha); zz=sort(Zbs,'descend');
%         hold on, plot(x,zz(r)*ones(1,length(z)),'r--','LineWidth',2);
%         set(gca,'FontSize',22), box off
%         title(['(' num2str(f) ',' num2str(s) '): ' num2str(nassem{s})]);
%         xlabel('Factor (# assemblies)'), ylabel('log-likelihood-ratio')
%         legend('original','trial-permutation'), legend('boxoff')
%         subplot(1,2,2), hold off cla
%         vs=cell2mat(FLbs{s}{2});
%         x=-0.5:0.005:0.5; h=histc(vs(1:end),x)./length(vs(1:end));
%         plot(x,h,'b',[ciHld(s,2) ciHld(s,2)],[0 0.05],'g','LineWidth',2);
%         hold on, plot(FL{s}{2}(1:end),0.01,'r.'); xlim([min(x) max(x)]);
%         set(gca,'FontSize',22), box off
%         
% %         subplot(1,3,3), hold off cla
% %         ws=sort(FSCbs{s}(1:end),'ascend');
% %         ciV=[0.1 0.05 0.01 5e-3 1e-3];
% %         for i=1:length(ciV)
% %             ciLsc(s,i)=ws(round(length(ws)*ciV(i)));
% %             ciHsc(s,i)=ws(round(length(ws)*(1-ciV(i))));
% %         end;
% %         x=-3:0.1:3; h=histc(ws,x)./length(ws);
% %         plot(x,h,'b',[ciHsc(s,2) ciHsc(s,2)],[0 0.05],'g','LineWidth',2);
% %         hold on, plot(FSC{s}(1:end),0.01,'r.'); xlim([min(x) max(x)]);
% %         set(gca,'FontSize',22), box off
% %         
% %         a=input('...');
%     end;
% end;
%% visualize factor scores/ activations in time & distribution
% for f=1:length(A)
%     fnOut=[pat target '/' Name{f} '_AssemRes2']; load(fnOut);
%     fnOut=[pat target '/' Name{f} '__FSCtemp']; load(fnOut);
%     for s=1:length(FSC)
%         [ntp,nu]=size(FSC{s});
%         for u=1:nu
%             [s u]
%             figure
%             n=find(abs(FL{s}{nassem{s}(3)}(:,u))>ciHld(s,3));
%             subplot(1,2,1), hold off cla
%             t=(1:ntp)*0.05+0.05/2;
%             plot(t,FSC{s}(:,u),'b',t,ciHsc(s,2),'r-','LineWidth',2);
%             set(gca,'FontSize',22), box off, xlim([680 700]);
%             title(['#units: ' num2str(n')]);
%             subplot(1,2,2), hold off cla
%             x=-3:0.01:3; h=histc(FSC{s}(:,u),x)./length(FSC{s}(:,u));
%             plot(x,h,'b',[ciHsc(s,2) ciHsc(s,2)],[0 0.05],'g','LineWidth',2);
%             xlim([-0.5 1]);
%             set(gca,'FontSize',22), box off
%             title(['#units: ' num2str(n')]);
% %             a=input([num2str([f s u]) '...']);
%         end;
%     end;
% end;

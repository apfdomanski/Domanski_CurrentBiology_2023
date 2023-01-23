%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Target = 'LONG';
Areas = {'PFC','HP','Joint'};
TimeSpan = 4;

bw=0.05;
tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

maxRange=20;    % Upper limit of assembly size to explore
noTrials = 5;   % How many trials to consider simultanrously
nReps = 50;     % how many repeated draws among trials to perform
        
        
plotOnline = false;
useWholeTaskPeriod = true;
ResampleTrials = true;
InfoCriteria = 'max'; % 'mean','max'


clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);


% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
if sum(name_flag)>0
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag
AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time
switch AssemblyChoice
    case 1
        pat2 = [pat 'KDE_bins' filesep Target filesep];
    case 2
        pat2 = [pat 'KDE_binsTaskonly' filesep 'LONGTaskonly' filesep];
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end
%% Batch Process
for iFile = 1:length(fileList)
    %% Get spike times and rates, behaviour
    fn = fullfile(pat,strtok(fileList(iFile).name,'_'));
    spikes = load(fn);
    if useWholeTaskPeriod
        fn = fullfile(pat,'KDE_binsTaskOnly',[strtok(fileList(iFile).name,'_'),'_PFC_iFR50_behavOnly.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{1} = temp.iFR;clear temp
        fn = fullfile(pat,'KDE_binsTaskOnly',[strtok(fileList(iFile).name,'_'),'_HP_iFR50_behavOnly.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{2} = temp.iFR;
        Tmtx=temp.Tmtx;clear temp
    else
        fn = fullfile(pat,'KDE_bins',[strtok(fileList(iFile).name,'_'),'_PFC_iFR50.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{1} = temp.iFR;clear temp
        fn = fullfile(pat,'KDE_bins',[strtok(fileList(iFile).name,'_'),'_HP_iFR50.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{2} = temp.iFR;
        Tmtx=temp.Tmtx; clear temp
    end
    iFR{3} = [iFR{1},iFR{2}];
    
    fname=strtok(fileList(iFile).name,'_')
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    %% Aggregate all trials regardless of delay length
    LeftTrials=[];RightTrials=[];
    for iDelay =1:length(Delays_)
        eval(sprintf('LeftTrials = [LeftTrials;[t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect]];',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf('RightTrials = [RightTrials;[t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect]];',Delays_{iDelay},Delays_{iDelay}));
    end
    % Run decoder once on all trials to rank order the neurons for later selection
  
    for s=1:2
        Ltrials_{s} = {}; Rtrials_{s} = {};
        Ltrials = [];nL = 0;
        for iTrial =1:size(LeftTrials,1)
            try
                tlims_  = [LeftTrials(iTrial,1)/1e6;LeftTrials(iTrial,2)/1e6]+tlimsAll(1);
                tlims_  = closest(Tmtx,tlims_);
                tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                Ltrials = [Ltrials;iFR{s}(tlims_,:)];
                Ltrials_{s}{1,nL+1} = iFR{s}(tlims_,:);
                nL=nL+1;
            end
        end
        
        Rtrials = [];nR = 0; 
        for iTrial =1:size(RightTrials,1)
            try
                tlims_  = [RightTrials(iTrial,1)/1e6;RightTrials(iTrial,2)/1e6]+tlimsAll(1);
                tlims_  = closest(Tmtx,tlims_);
                tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                Rtrials = [Rtrials;iFR{s}(tlims_,:)];
                Rtrials_{s}{1,nR+1} = iFR{s}(tlims_,:);
                nR=nR+1;
            end
        end
        
        FR{s} = [Ltrials;Rtrials];
        evt0 = [ones(nL,1);2*ones(nR,1)];
        Ltr = 2*length(tbAll);
    end
   
    Ltrials_{3}=cellfun(@(x,y)[x,y],Ltrials_{1},Ltrials_{2},'UniformOutput',false);
    Rtrials_{3}=cellfun(@(x,y)[x,y],Rtrials_{1},Rtrials_{2},'UniformOutput',false);
    clear Ltrials Rtrials iTrial
    %% Get assembly information
    
    if useWholeTaskPeriod
        fn = fullfile(pat,'KDE_binsTaskOnly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_FSC.mat',strtok(fileList(iFile).name,'_')));
        A  = load(fn);
        fn = fullfile(pat,'KDE_binsTaskOnly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_AssemRes2.mat',strtok(fileList(iFile).name,'_')));
        B  = load(fn);
        % recalculate included units
        fn = fullfile(pat,'KDE_binsTaskOnly',sprintf('%s_PFC_iFR50_behavOnly.mat',strtok(fileList(iFile).name,'_')));
        usel_out=SelCells(fn,0.1,1e6);
    else
        fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',strtok(fileList(iFile).name,'_')));
        A  = load(fn);
        fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',strtok(fileList(iFile).name,'_')));
        B  = load(fn);
        % recalculate included units
        fn = fullfile(pat,'KDE_bins',sprintf('%s_PFC_iFR50.mat',strtok(fileList(iFile).name,'_')));
        [~,~,~,~,B.usel_out]=SelTrialsCellsWholeTrial(fn,10,0.1,1e6,'iFR');
    end
    
    spikes.jointcells=[spikes.PFCcells;spikes.HPcells]';
    A.nu(3) = sum(A.nu(1:2));
    %     B.usel_out{3}= [B.usel_out{1},B.usel_out{2}+max(B.usel_out{1})];
    B.usel_out{3}= [B.usel_out{1},B.usel_out{2}+size(iFR{1},2)];
    %% Sort units by assembly membership types
        types_={'Nonmember','Localmember','Jointmember'};
        for s = 1:2
            
            % Global = joint + local
            % Nomember = ~Global
            switch s
                case 1
                    joint_ = unique(cell2mat(A.units{3})); joint_(joint_>A.nu(1))=[];
                case 2
                    joint_ = unique(cell2mat(A.units{3})); joint_ = joint_(joint_>A.nu(1))-A.nu(1);
            end
            local_ = unique(cell2mat(A.units{s}));
            all_   = 1:A.nu(s);
            LocalmemberUnits    = setdiff(local_,joint_);                      % local only members
            JointmemberUnits    = joint_;                                      % inter-area only members
            GlobalmemberUnits   = unique([JointmemberUnits,LocalmemberUnits]); % membership of any assembly
            NonmemberUnits      = setdiff(all_,GlobalmemberUnits);             % global nonmembers
            
            LocalmemberUnits    = B.usel_out{s}(LocalmemberUnits);
            NonmemberUnits      = B.usel_out{s}(NonmemberUnits);
            JointmemberUnits    = B.usel_out{s}(JointmemberUnits);
            GlobalmemberUnits   = B.usel_out{s}(GlobalmemberUnits);
            clear joint_ local_ all_
            
            D.fname{iFile} = fname;
            D.nonmemberUnits{iFile}{s} = NonmemberUnits;
            D.localmemberUnits{iFile}{s} = LocalmemberUnits;
            D.jointmemberUnits{iFile}{s} = JointmemberUnits;
            
        end %/s
  
	%% Draw random groups of 1:N units and work out average decoding from each type
    fprintf('Calculating average single units: Assembly member neurons\n')
    D_.AvDecoding.Nonmember{iFile} = MakeAverageAssems(D.nonmemberUnits{iFile}, Ltrials_,Rtrials_,20);
    D_.AvDecoding.Localmember{iFile} = MakeAverageAssems(D.localmemberUnits{iFile}, Ltrials_,Rtrials_,20);
    D_.AvDecoding.Jointmember{iFile} = MakeAverageAssems(D.jointmemberUnits{iFile}, Ltrials_,Rtrials_,20);
    %%
    col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
    tb = (1:length(tbAll)*2)*bw;
    
    figure
    for s=1:3
        subplot(1,3,s); hold on
        plot(1:20,D_.AvDecoding.Nonmember{iFile}.Score{s},'color',col_{1});
        plot(1:20,D_.AvDecoding.Localmember{iFile}.Score{s},'color',col_{2});
        plot(1:20,D_.AvDecoding.Jointmember{iFile}.Score{s},'color',col_{3});
        axis([2 20 0.4 1])
    end
    %%
    figure
    for s=1:3
        subplot(1,3,s); hold on
%         plot(tb,D_.AvDecoding.Nonmember{iFile}.CVE{1},'color',col_{1});
%         plot(tb,D_.AvDecoding.Localmember{iFile}.CVE{2},'color',col_{2});
%         plot(tb,D_.AvDecoding.Jointmember{iFile}.CVE{3},'color',col_{3});
        x = D_.AvDecoding.Nonmember{iFile}.CVE{1}';
        ciplot(nanmean(x)+nansem(x),nanmean(x)-nansem(x),tb,col_{1},0.5)
        x = D_.AvDecoding.Nonmember{iFile}.CVE{2}';
        ciplot(nanmean(x)+nansem(x),nanmean(x)-nansem(x),tb,col_{2},0.5)
        x = D_.AvDecoding.Nonmember{iFile}.CVE{3}';
        ciplot(nanmean(x)+nansem(x),nanmean(x)-nansem(x),tb,col_{3},0.5)
        
    end
end

%% Batch plot time axis
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
tb = (1:length(tbAll)*2)*bw;
figure
N = [2 10];
for s=1:3
    subplot(1,3,s); hold on
    Nonmember{s} = [];
    Member{s} = [];
    Jointmember{s} = [];
    
    Nonmember_H{s} = [];
    Member_H{s} = [];
    Jointmember_H{s} = [];
    
    Nonmember_L{s} = [];
    Member_L{s} = [];
    Jointmember_L{s} = [];
    
    Nonmember_HCVE{s} = [];
    Member_HCVE{s} = [];
    Jointmember_HCVE{s} = [];
    
    Nonmember_LCVE{s} = [];
    Member_LCVE{s} = [];
    Jointmember_LCVE{s} = [];
    
    Nonmember_sig{s} = [];
    Member_sig{s} = [];
    Jointmember_sig{s} = [];
    
    
    for iFile = 1:length(fileList)
        Nonmember{s}           = [Nonmember{s};D_.AvDecoding.Nonmember{iFile}.CVE{s}(:,N(1):N(2))'];
        Member{s}              = [Member{s};D_.AvDecoding.Localmember{iFile}.CVE{s}(:,N(1):N(2))'];
        Jointmember{s}         = [Jointmember{s};D_.AvDecoding.Jointmember{iFile}.CVE{s}(:,N(1):N(2))'];
        
        Nonmember_H{s}         = [Nonmember_H{s},D_.AvDecoding.Nonmember{iFile}.Score_shuffled_95pc{s}(N(1):N(2))];
        Member_H{s}            = [Member_H{s},D_.AvDecoding.Localmember{iFile}.Score_shuffled_95pc{s}(N(1):N(2))];
        Jointmember_H{s}       = [Jointmember_H{s},D_.AvDecoding.Jointmember{iFile}.Score_shuffled_95pc{s}(N(1):N(2))];
        
        Nonmember_L{s}         = [Nonmember_L{s},D_.AvDecoding.Nonmember{iFile}.Score_shuffled_5pc{s}(N(1):N(2))];
        Member_L{s}            = [Member_L{s},D_.AvDecoding.Localmember{iFile}.Score_shuffled_5pc{s}(N(1):N(2))];
        Jointmember_L{s}       = [Jointmember_L{s},D_.AvDecoding.Jointmember{iFile}.Score_shuffled_5pc{s}(N(1):N(2))];
        
        
        x = cell2mat(D_.AvDecoding.Nonmember{iFile}.CVE_shuffled{s}); y=x; x(:,1:2:(end-1))=[]; y(:,2:2:(end))=[];
        Nonmember_HCVE{s}      = [Nonmember_HCVE{s}; x(:,(N(1):N(2))-1)'];
        Nonmember_LCVE{s}      = [Nonmember_LCVE{s}; y(:,(N(1):N(2))-1)'];
        
        x = cell2mat(D_.AvDecoding.Localmember{iFile}.CVE_shuffled{s}); y=x; x(:,1:2:(end-1))=[]; y(:,2:2:(end))=[];
        Member_HCVE{s}         = [Member_HCVE{s}; x(:,(N(1):N(2))-1)'];
        Member_LCVE{s}         = [Member_LCVE{s}; y(:,(N(1):N(2))-1)'];
        
        x = cell2mat(D_.AvDecoding.Jointmember{iFile}.CVE_shuffled{s}); y=x; x(:,1:2:(end-1))=[]; y(:,2:2:(end))=[];
        Jointmember_HCVE{s}    = [Jointmember_HCVE{s}; x(:,(N(1):N(2))-1)'];
        Jointmember_LCVE{s}    = [Jointmember_LCVE{s}; y(:,(N(1):N(2))-1)'];
        
        
        
        %          Nonmember{s}=[Nonmember{s};nanmean(D_.AvDecoding.Nonmember{iFile}.CVE{s}(:,N(1):N(2))')];
        %          Member{s}=[Member{s};nanmean(D_.AvDecoding.Localmember{iFile}.CVE{s}(:,N(1):N(2))')];
        %          Jointmember{s}=[Jointmember{s};nanmean(D_.AvDecoding.Jointmember{iFile}.CVE{s}(:,N(1):N(2))')];
        
        
        
    end
    
    FracSig=0.3;
    Nonmember_sig{s}    = Nonmember{s}>Nonmember_HCVE{s};
    Nonmember_sig{s}(sum(isnan(Nonmember{s}),2)>0,:)=[];
    Nonmember_sig{s}    = nansum(Nonmember_sig{s})./size(Nonmember{s},1);
    
    Member_sig{s}    = Member{s}>Member_HCVE{s};
    Member_sig{s}(sum(isnan(Member{s}),2)>0,:)=[];
    Member_sig{s}    = nansum(Member_sig{s})./size(Member_sig{s},1);
    
    Jointmember_sig{s}    = Jointmember{s}>Jointmember_HCVE{s};
    Jointmember_sig{s}(sum(isnan(Jointmember{s}),2)>0,:)=[];
    Jointmember_sig{s}    = nansum(Jointmember_sig{s})./size(Jointmember{s},1);
    
    
    %     ciplot(nanmean(Nonmember_LCVE{s}),nanmean(Nonmember_HCVE{s}),tb,col_{1})
    %     ciplot(nanmean(Member_LCVE{s}),nanmean(Member_HCVE{s}),tb,col_{2})
    %     ciplot(nanmean(Jointmember_LCVE{s}),nanmean(Jointmember_HCVE{s}),tb,col_{3})
    
    H = nanmean(nanmean(Nonmember_H{s}));
    L = nanmean(nanmean(Nonmember_L{s}));
    Y1 = [L,H-L].*ones(length(tb),2);
    h=area([tb;tb]',Y1); set(h(1),'FaceAlpha',0,'LineStyle','none'); set(h(2),'FaceColor',col_{1},'FaceAlpha',0.2,'LineStyle','none')
    
    H = nanmean(nanmean(Member_H{s}));
    L = nanmean(nanmean(Member_L{s}));
    Y1 = [L,H-L].*ones(length(tb),2);
    h=area([tb;tb]',Y1); set(h(1),'FaceAlpha',0,'LineStyle','none'); set(h(2),'FaceColor',col_{2},'FaceAlpha',0.2,'LineStyle','none')
    
    H = nanmean(nanmean(Jointmember_H{s}));
    L = nanmean(nanmean(Jointmember_L{s}));
    Y1 = [L,H-L].*ones(length(tb),2);
    h=area([tb;tb]',Y1); set(h(1),'FaceAlpha',0,'LineStyle','none'); set(h(2),'FaceColor',col_{3},'FaceAlpha',0.2,'LineStyle','none')
    
    ciplot(nanmean(Nonmember{s})+nansem(Nonmember{s}),nanmean(Nonmember{s})-nansem(Nonmember{s}),tb,col_{1},0.8)
    ciplot(nanmean(Member{s})+nansem(Member{s}),nanmean(Member{s})-nansem(Member{s}),tb,col_{2},0.8)
    ciplot(nanmean(Jointmember{s})+nansem(Jointmember{s}),nanmean(Jointmember{s})-nansem(Jointmember{s}),tb,col_{3},0.8)
    plot([min(tb) max(tb)],[0.5 0.5],':k','LineWidth',1.5)
    
    
    z = [Nonmember_sig{s};Member_sig{s};Jointmember_sig{s}];
    imagesc(tb,[1 1.05 1.1],z); set(gca,'ydir','normal');
    colormap(flipud(gray))
    
    rectangle('Position',[0 0.975 max(tb) 0.05],'LineWidth',1.5,'Edgecolor',col_{1})
    rectangle('Position',[0 1.025 max(tb) 0.05],'LineWidth',1.5,'Edgecolor',col_{2})
    rectangle('Position',[0 1.075 max(tb) 0.05],'LineWidth',1.5,'Edgecolor',col_{3})
    axis([-1 21 0.3 1.2])
    
    title([Areas{s} ' Neurons'])
    if s==1
        ylabel('Fraction correct decoding')
        text(20,0.46,'Chance','HorizontalAlignment','Right')

    elseif s==2
        xlabel('Time (s)')
    end
    caxis([0 1])
    set(gca,'YTick',[0.3:0.1:1])
end
% legend
text(20,1.16,'Frac.(real>shuffled)','HorizontalAlignment','Right')
text(20,0.31,'Non-members','Color',col_{1},'HorizontalAlignment','Right','VerticalAlignment','baseline')
text(20,0.37,'Local Members','Color',col_{2},'HorizontalAlignment','Right','VerticalAlignment','baseline')
text(20,0.44,'Joint Members','Color',col_{3},'HorizontalAlignment','Right','VerticalAlignment','baseline')

%% Histogram of peak decoding

col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
bins = 0:0.01:1;
figure
N = [2 20];
for s=1:3
    subplot(1,3,s); hold on
    Nonmember{s} = [];
    Member{s} = [];
    Jointmember{s} = [];
    for iFile = 1:length(fileList)
         Nonmember{s}=[Nonmember{s};D_.AvDecoding.Nonmember{iFile}.Score{s}(N(1):N(2))];
         Member{s}=[Member{s};D_.AvDecoding.Localmember{iFile}.Score{s}(N(1):N(2))];
         Jointmember{s}=[Jointmember{s};D_.AvDecoding.Jointmember{iFile}.Score{s}(N(1):N(2))];
    end
    x = cumsum(histc(Nonmember{s},bins));
    x=x./max(x);
    plot(bins,x,'Color',col_{1},'LineWidth',1.5)
    
    x = cumsum(histc(Member{s},bins));
    x=x./max(x);
    plot(bins,x,'Color',col_{2},'LineWidth',1.5)
    
    x = cumsum(histc(Jointmember{s},bins));
    x=x./max(x);
    plot(bins,x,'Color',col_{3},'LineWidth',1.5)
    axis([0.5 1 0 1])
end
    
 %% Decoding wrt ensemble size

col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
figure
for s=1:3
    subplot(1,3,s); hold on
    Nonmember{s} = [];
    Member{s} = [];
    Jointmember{s} = [];
    for iFile = 1:length(fileList)
         Nonmember{s}=[Nonmember{s},D_.AvDecoding.Nonmember{iFile}.Score{s}(2:end)];
         Member{s}=[Member{s},D_.AvDecoding.Localmember{iFile}.Score{s}(2:end)];
         Jointmember{s}=[Jointmember{s},D_.AvDecoding.Jointmember{iFile}.Score{s}(2:end)];
    end
  
    x =  Nonmember{s}';
    ciplot(nanmean(x)+nansem(x),nanmean(x)-nansem(x),2:20,col_{1})
    
    x =  Member{s}';
    ciplot(nanmean(x)+nansem(x),nanmean(x)-nansem(x),2:20,col_{2})
    
    x =  Jointmember{s}';
    ciplot(nanmean(x)+nansem(x),nanmean(x)-nansem(x),2:20,col_{3})
    axis([2 20 0.6 1])
end   
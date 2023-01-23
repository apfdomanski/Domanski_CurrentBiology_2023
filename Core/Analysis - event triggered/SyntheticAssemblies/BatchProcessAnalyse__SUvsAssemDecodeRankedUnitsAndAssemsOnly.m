%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Target = 'LONG';
Areas = {'PFC','HP','Joint'};

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
    
    fname=strtok(fileList(iFile).name,'_');
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
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
    %% Mark up units as members/nonmembers
    
    % First strip any Joint assembly that is exclusively single area
    % contributed - 
    % THIS HAS ALREADY BEEN PERFORMED ON OUTPUT FROM F.A. ASSEMBLY DETECTION
%     for iAss =1:length(A.units{3})
%         if nansum(A.units{3}{iAss}<=A.nu(1))==0 | nansum(A.units{3}{iAss}>A.nu(1))==0
%             A.units{3}{iAss}=[];
%         end
%     end
        
    % Indices of included units
    % PFC
    x = cell2mat(A.units{3});
    x(x>length(B.usel_out{1}))=[];
    x = unique([cell2mat(A.units{1}), x]);
    globalmemberUnits{1} = x;
    nonmemberUnits{1}    = setdiff(1:A.nu(1),x);
    
    % HP
    x = cell2mat(A.units{3});
    x = x(x>length(B.usel_out{1}))-length(B.usel_out{1});
    x = unique([cell2mat(A.units{2}), x]);
    globalmemberUnits{2} = x;
    nonmemberUnits{2}    = setdiff(1:A.nu(2),x);
    
    % Joint
    x = unique([cell2mat(A.units{1}), cell2mat(A.units{2})+length(B.usel_out{1}),cell2mat(A.units{3})]);
    globalmemberUnits{3} = x;
    nonmemberUnits{3}   = setdiff(1:A.nu(3),x);
    JointMemberUnits{3} = unique(cell2mat(A.units{3}));
    
    for s=1:3
        globalmemberUnitsReal{s}    = B.usel_out{s}(globalmemberUnits{s});
        globalnonmemberUnitsReal{s} = B.usel_out{s}(nonmemberUnits{s});
    end
    
    JointMemberUnitsReal{3} = B.usel_out{3}(JointMemberUnits{3});
    
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
    %% (1) Decoding from REAL ASSEMBLIES
    for s=1:3
        nAss = length(A.units{s});
        D_.AssReal.CVE{s}                   = nan(Ltr,nAss);
        D_.AssReal.score{s}                 = nan(nAss,2);
        D_.AssReal.score_shuffled{s}        = nan(nAss,2);    
        D_.AssReal.score_shuffled_5pc{s}    = nan(nAss,2);    
        D_.AssReal.score_shuffled_95pc{s}   = nan(nAss,2);   
        for iAss = 1:nAss
            if ~isempty(A.units{s}{iAss})
                fprintf('Decoding Real %s Assemblies : Assembly %d (of %d)...\n',Areas{s},iAss,nAss)
                
                unitIds = B.usel_out{s}(A.units{s}{iAss});
                noUnits = length(unitIds);
                if ResampleTrials
                    [D_.AssReal.CVE{s}(:,iAss),...
                        score,...
                        score_shuffled,...
                        score_shuffled_5pc,...
                        score_shuffled_95pc] = RunDecodeGroups(unitIds,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
                    
                    D_.AssReal.score{s}(iAss,:)                    =  [noUnits,score];
                    D_.AssReal.score_shuffled{s}(iAss,:)           =  [noUnits,score_shuffled];
                    D_.AssReal.score_shuffled_5pc{s}(iAss,:)       =  [noUnits,score_shuffled_5pc];
                    D_.AssReal.score_shuffled_95pc{s}(iAss,:)      =  [noUnits,score_shuffled_95pc];
                else
                    CVE_ = 1-DecodeCVE_mex(FR{s}(:,unitIds),evt0,0.05);
                    CVEshuf_ = 1-DecodeCVE_mex(FR{s}(:,unitIds),evt0(randperm(length(evt0))),0.05);
                    D_.AssReal.CVE{s}(:,iAss) = CVE_;
                    switch InfoCriteria
                        case 'mean'
                            D_.AssReal.score{s}(iAss,:)          = [noUnits,mean(CVE_,2)];
                            D_.AssReal.score_shuffled{s}(iAss,:) = [noUnits,mean(CVEshuf_,2)];
                        case 'max'
                            D_.AssReal.score{s}(iAss,:)          = [noUnits,max(CVE_,[],2)];
                            D_.AssReal.score_shuffled{s}(iAss,:) = [noUnits,max(CVEshuf_,[],2)];
                    end
                end
            end
        end
    end
    %% (2) Rank all single units by decoding power 
    for s=1:2
        noUnits             = size(FR{s},2);
        CVE                 = zeros(noUnits,Ltr);
        CVEscore            = zeros(noUnits,1);
        CVEscore_Shuffled   = zeros(noUnits,1);
        CVErank             = zeros(noUnits,1);
        
        
        parfor iUnit =1:noUnits
            fprintf('Decoding all units: %s unit %d of %d...\n',Areas{s},iUnit,noUnits)
            [CVE(iUnit,:),CVEscore(iUnit),CVEscore_Shuffled(iUnit)]= RunDecodeGroups(iUnit,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
            CVErank(iUnit) = iUnit;
        end
        switch InfoCriteria
            case 'max'
                CVEscore = max(CVE,[],2);
            case 'mean'
                CVEscore = mean(CVE,2);
        end
        
        % ALL neurons
        unitIDs = 1:noUnits;
        [CVEscore_, rank_] = sort(CVEscore(unitIDs),'descend');
        CVE_ = CVE(rank_,:);
        D_.Units.All.CVE{s}                     = CVE_;
        D_.Units.All.SortedScore{s}             = CVEscore_;
        D_.Units.All.SortedScore_shuffled{s}    = CVEscore_Shuffled(rank_);
        D_.Units.All.Rank{s}                    = unitIDs(rank_);
        
        % MEMBER neurons
        unitIDs = globalmemberUnitsReal{s};
        CVE_ = CVE(unitIDs,:);
        CVEscore_ = CVEscore(unitIDs);
        [CVEscore_, rank_] = sort(CVEscore_,'descend');
        CVE_ = CVE_(rank_,:);
        D_.Units.Members.CVE{s}                  = CVE_;
        D_.Units.Members.SortedScore{s}          = CVEscore_;
        D_.Units.Members.SortedScore_shuffled{s} = CVEscore_Shuffled(rank_);
        D_.Units.Members.Rank{s}                 = unitIDs(rank_);
        
        % NON-MEMBER neurons
        unitIDs = globalnonmemberUnitsReal{s};
        CVE_ = CVE(unitIDs,:);
        CVEscore_ = CVEscore(unitIDs);
        [CVEscore_, rank_] = sort(CVEscore_,'descend');
        CVE_ = CVE_(rank_,:);
        D_.Units.Nonmembers.CVE{s}                  = CVE_;
        D_.Units.Nonmembers.SortedScore{s}          = CVEscore_;
        D_.Units.Nonmembers.SortedScore_shuffled{s} = CVEscore_Shuffled(rank_);
        D_.Units.Nonmembers.Rank{s}                 = unitIDs(rank_);
        
        % Joint HP-PFC member neurons
        if s==1
            unitIDs = JointMemberUnitsReal{3}( JointMemberUnitsReal{3}<=size(FR{1},2));
        elseif s==2
            unitIDs = JointMemberUnitsReal{3}( JointMemberUnitsReal{3}>size(FR{1},2))-size(FR{1},2);
        end
        CVE_ = CVE(unitIDs,:);
        CVEscore_ = CVEscore(unitIDs);
        [CVEscore_, rank_] = sort(CVEscore_,'descend');
        CVE_ = CVE_(rank_,:);
        D_.Units.JointMembers.CVE{s}                  = CVE_;
        D_.Units.JointMembers.SortedScore{s}          = CVEscore_;
        D_.Units.JointMembers.SortedScore_shuffled{s} = CVEscore_Shuffled(rank_);
        D_.Units.JointMembers.Rank{s}                 = unitIDs(rank_);

        clear CVE CVEscore CVEscore_ rank_ CVE_
       
    end
    %% (4) Save results
    fnOut = sprintf('%sSyntheticAssemblies%s%s_AssemblyComparisonMvsNMvsAss.mat',pat,filesep,fname);
    save(fnOut,'D_','A','B','Areas','maxRange','InfoCriteria','tlimsAll','tbAll','Delays_','Ltrials_','Rtrials_','-v7.3')
    
    %% clear up
    clear A B CVE CVE_ CVErank CVEscore CVEscore CVEshuf_ evt0 FR iFR globalmemberUnits globalmemberUnitsReal globalnonmemberUnitsReal JointMemberUnits nonmemberUnits x
    clear inputData LeftTrials Ltrials_ nL nAss noUnits  nR Rtrials_  RightTrials score iUnit score_shuffled score_shuffled_5pc score_shuffled_95pc spikes t unitIds Tmtx unitIDs usel_out
end
%% Reimport for analysis
for s=1:2
    D.Units.All.SortedScore{s} = [];
    D.Units.Members.SortedScore{s} = [];
    D.Units.Nonmembers.SortedScore{s} = [];
    D.Units.JointMembers.SortedScore{s} = [];
    
    D.Units.All.SortedScore_shuffled{s} = [];
    D.Units.Members.SortedScore_shuffled{s} = [];
    D.Units.Nonmembers.SortedScore_shuffled{s} = [];
    D.Units.JointMembers.SortedScore_shuffled{s} = [];
end
for s=1:3
   D.AssReal.score{s}=[]; 
   D.AssReal.score_shuffled_5pc{s}=[]; 
   D.AssReal.score_shuffled_95pc{s}=[]; 
end
for iFile = 1:length(fileList)
    fname=strtok(fileList(iFile).name,'_');

    fnIn = sprintf('%sSyntheticAssemblies%s%s_AssemblyComparisonMvsNMvsAss.mat',pat,filesep,fname);
    load(fnIn,'D_');
 
    for s=1:2
        D.Units.All.SortedScore{s} = [D.Units.All.SortedScore{s}; D_.Units.All.SortedScore{s}];
        D.Units.Members.SortedScore{s} = [D.Units.Members.SortedScore{s}; D_.Units.Members.SortedScore{s}];
        D.Units.Nonmembers.SortedScore{s} = [D.Units.Nonmembers.SortedScore{s}; D_.Units.Nonmembers.SortedScore{s}];
        D.Units.JointMembers.SortedScore{s} = [D.Units.JointMembers.SortedScore{s}; D_.Units.JointMembers.SortedScore{s}];
        
        D.Units.All.SortedScore_shuffled{s} = [D.Units.All.SortedScore_shuffled{s}; D_.Units.All.SortedScore_shuffled{s}];
        D.Units.Members.SortedScore_shuffled{s} = [D.Units.Members.SortedScore_shuffled{s}; D_.Units.Members.SortedScore_shuffled{s}];
        D.Units.Nonmembers.SortedScore_shuffled{s} = [D.Units.Nonmembers.SortedScore_shuffled{s}; D_.Units.Nonmembers.SortedScore_shuffled{s}];
        D.Units.JointMembers.SortedScore_shuffled{s} = [D.Units.JointMembers.SortedScore_shuffled{s}; D_.Units.JointMembers.SortedScore_shuffled{s}];
        
        D.Units.All.SortedScore_all{s}{iFile}           = D_.Units.All.SortedScore{s};
        D.Units.Members.SortedScore_all{s}{iFile}       = D_.Units.Members.SortedScore{s};
        D.Units.Nonmembers.SortedScore_all{s}{iFile}    = D_.Units.Nonmembers.SortedScore{s};
        D.Units.JointMembers.SortedScore_all{s}{iFile}  = D_.Units.JointMembers.SortedScore{s}; 
    end
    
    for s=1:3
        D.AssReal.score{s}=[D.AssReal.score{s}; D_.AssReal.score{s}]; 
        D.AssReal.score_shuffled_5pc{s}=[D.AssReal.score_shuffled_5pc{s}; D_.AssReal.score_shuffled_5pc{s}]; 
        D.AssReal.score_shuffled_95pc{s}=[D.AssReal.score_shuffled_95pc{s}; D_.AssReal.score_shuffled_95pc{s}]; 
    end
end
%% Plot distribution of decoding scores for single unit
bins = 0.4:0.02:1;
figure;
for s=1:2
    
    subplot(2,1,s); hold on
    
    x= histc(D.Units.All.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.All.SortedScore_shuffled{s});
    area(bins,x,'HandleVisibility','off','FaceColor','k','EdgeColor','none','FaceAlpha',0.3)
    
    x= histc(D.Units.Nonmembers.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.Nonmembers.SortedScore_shuffled{s});
    area(bins,x,'HandleVisibility','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.3)

    x= histc(D.Units.Members.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.Members.SortedScore_shuffled{s});
    area(bins,x,'HandleVisibility','off','FaceColor','b','EdgeColor','none','FaceAlpha',0.3)

    x= histc(D.Units.JointMembers.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.JointMembers.SortedScore_shuffled{s});
    area(bins,x,'HandleVisibility','off','FaceColor','g','EdgeColor','none','FaceAlpha',0.3)    
    
    x= histc(D.Units.All.SortedScore{s},bins);
    x= x./numel(D.Units.All.SortedScore{s});
    stairs(bins,cumsum((x)),'k','LineWidth',1.5)
    
    x= histc(D.Units.Nonmembers.SortedScore{s},bins);
    x= x./numel(D.Units.Nonmembers.SortedScore{s});
    stairs(bins,cumsum((x)),'r','LineWidth',1.5)

    x= histc(D.Units.Members.SortedScore{s},bins);
    x= x./numel(D.Units.Members.SortedScore{s});
    stairs(bins,cumsum((x)),'b','LineWidth',1.5)
    
    x= histc(D.Units.JointMembers.SortedScore{s},bins);
    x= x./numel(D.Units.JointMembers.SortedScore{s});
    stairs(bins,cumsum((x)),'g','LineWidth',1.5)    
    
    text(0.49,1,'Chance/Shuffled','HorizontalAlignment','right','Rotation',90)
    plot([0.5 0.5],[0 1],':k')
    axis([min(bins) max(bins) 0 1])
    title([Areas{s} ' units'])
    ylabel('Fraction of units')
    if s==2
        xlabel({'Peak decoding';'(% Correct classification)'})    
    end
    legend({'All Units','Non-members','Members','Joint members','Chance'},'Location','southeast');legend boxoff
end
%% Plot distribution of decoding scores for single unit
BinWidth = 0.02;
bins = 0.4:BinWidth:1;
figure;
for s=1:2
    
    subplot(2,1,s); hold on
    title([Areas{s} ' units'])
    x= histc(D.Units.All.SortedScore{s},bins);
    x= x./numel(D.Units.All.SortedScore{s});
%     stairs(bins,smooth_hist(x),'k')
    stairs(bins,x,'k','LineWidth',1.5)
%     bar(bins,(x),'FaceColor','k','FaceAlpha',0.5,'EdgeColor','none')
    
    x= histc(D.Units.All.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.All.SortedScore_shuffled{s});
%     stairs(bins,smooth_hist(x),'k')
    stairs(bins,((x)),'k','HandleVisibility','off')
    
    x= histc(D.Units.Nonmembers.SortedScore{s},bins);
    x= x./numel(D.Units.Nonmembers.SortedScore{s});
%     stairs(bins,smooth_hist(x),'r')
    stairs(bins,((x)),'r','LineWidth',1.5)
%     bar(bins,(x),'FaceColor','r','FaceAlpha',0.5,'EdgeColor','none')
    
     x= histc(D.Units.Nonmembers.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.Nonmembers.SortedScore_shuffled{s});
%     stairs(bins,smooth_hist(x),'r')
    stairs(bins,((x)),'r','HandleVisibility','off')

    x= histc(D.Units.Members.SortedScore{s},bins);
    x= x./numel(D.Units.Members.SortedScore{s});
%     stairs(bins,smooth_hist(x),'b')
    stairs(bins,((x)),'b','LineWidth',1.5)
%         bar(bins,(x),'FaceColor','b','FaceAlpha',0.5,'EdgeColor','none')

    x= histc(D.Units.Members.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.Members.SortedScore_shuffled{s});
%     stairs(bins,smooth_hist(x),'b')
    stairs(bins,((x)),'b','HandleVisibility','off')
    
    x= histc(D.Units.JointMembers.SortedScore{s},bins);
    x= x./numel(D.Units.JointMembers.SortedScore{s});
%     stairs(bins,smooth_hist(x),'g')
    stairs(bins,((x)),'g','LineWidth',1.5)    
%     bar(bins,(x),'FaceColor','g','FaceAlpha',0.5,'EdgeColor','none')

    x= histc(D.Units.JointMembers.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.JointMembers.SortedScore_shuffled{s});
%     stairs(bins,smooth_hist(x),'g')
    stairs(bins,((x)),'g','HandleVisibility','off')    
    
    plot([0.5 0.5],[0 1],':k')
    axis([min(bins) max(bins) 0 0.6])
    ylabel('Fraction of units')
    if s==2
        xlabel({'Peak decoding';'(% Correct classification)'})    
    end
    
end
legend({'All Units','Non-members','Members','Joint members','Chance'},'Location','southeast');legend boxoff

%% Plot distribution of decoding scores for single units - shaded regions
BinWidth = 0.05;
bins = 0.4:BinWidth:1;
nDrawn = 5;

figure;
for s = 1:2
    subplot(2,1,s); hold on
    title([Areas{s} ' units'])

%     x=[];
%     for iFile = 1:length(fileList)
%         x(:,iFile)= histc(D.Units.All.SortedScore_all{s}{iFile},bins);
%         x(:,iFile)= x(:,iFile)./numel(D.Units.All.SortedScore_all{s}{iFile});
%     end
%     x(:,isnan(sum(x)))=[];
%     ciplot(nanmean(cumsum(x),2)+nansem(cumsum(x),2),nanmean(cumsum(x),2)-nansem(cumsum(x),2),bins,'k')
    
    x=[];x_ = nan(length(bins),50);
    for iFile = 1:length(fileList)

            
        x(:,iFile)= histc(D.Units.Members.SortedScore_all{s}{iFile},bins);
        x(:,iFile)= x(:,iFile)./numel(D.Units.Members.SortedScore_all{s}{iFile});
        plot(bins,cumsum(x(:,iFile)),'color',[0 0 1 0.1])
    end
    x(:,isnan(sum(x)))=[];
    
    ciplot(nanmean(cumsum(x),2)+nansem(cumsum(x),2),nanmean(cumsum(x),2)-nansem(cumsum(x),2),bins,'b')

    x=[];
    for iFile = 1:length(fileList)
        x(:,iFile)= histc(D.Units.Nonmembers.SortedScore_all{s}{iFile},bins);
        x(:,iFile)= x(:,iFile)./numel(D.Units.Nonmembers.SortedScore_all{s}{iFile});
        plot(bins,cumsum(x(:,iFile)),'color',[1 0 0 0.1])
    end
    ciplot(nanmean(cumsum(x),2)+nansem(cumsum(x),2),nanmean(cumsum(x),2)-nansem(cumsum(x),2),bins,'r')
    
    x=[];
    for iFile = 1:length(fileList)
        x(:,iFile)= histc(D.Units.JointMembers.SortedScore_all{s}{iFile},bins);
        x(:,iFile)= x(:,iFile)./numel(D.Units.JointMembers.SortedScore_all{s}{iFile});
        plot(bins,cumsum(x(:,iFile)),'color',[0 1 0 0.1])
    end
    ciplot(nanmean(cumsum(x),2)+nansem(cumsum(x),2),nanmean(cumsum(x),2)-nansem(cumsum(x),2),bins,'g')
    
    ylabel('Fraction of units')
    if s==2
       xlabel({'Peak decoding score';'(% Correct Decoding)'})
    end        
        
end
    legend({'Non-members','Members','Joint members','Chance'},'Location','southeast');legend boxoff
%% Plot distribution of decoding scores for single units - matched numbers
BinWidth = 0.02;
bins = 0.4:BinWidth:1;
nDrawn = 1;
nBS = 500;
opacity  = 0.6;

col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};


figure;
for s = 1:2
    subplot(2,1,s); hold on
    title([Areas{s} ' units'])

%     x= histc(D.Units.All.SortedScore_shuffled{s},bins);
%     x= x./numel(D.Units.All.SortedScore_shuffled{s});
%     area(bins,x,'HandleVisibility','off','FaceColor','k','EdgeColor','none','FaceAlpha',0.2)
    
    x= histc(D.Units.Nonmembers.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.Nonmembers.SortedScore_shuffled{s});
    area(bins,x,'HandleVisibility','off','FaceColor',col_{1},'EdgeColor','none','FaceAlpha',opacity)

    x= histc(D.Units.Members.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.Members.SortedScore_shuffled{s});
    area(bins,x,'HandleVisibility','off','FaceColor',col_{2},'EdgeColor','none','FaceAlpha',opacity)

    x= histc(D.Units.JointMembers.SortedScore_shuffled{s},bins);
    x= x./numel(D.Units.JointMembers.SortedScore_shuffled{s});
    area(bins,x,'HandleVisibility','off','FaceColor',col_{3},'EdgeColor','none','FaceAlpha',opacity)    
    
%  x=[];x_ = nan(length(bins),nBS);
%     for iFile = 1:length(fileList)
%         for iDraw =1:nBS
%             try
%                 ids = randsample(1:length(D.Units.All.SortedScore_all{s}{iFile}),nDrawn);
%             catch 
%                 ids = 1:length(D.Units.All.SortedScore_all{s}{iFile});
%             end
%             x_(:,iDraw)= (histc(D.Units.All.SortedScore_all{s}{iFile}(ids),bins));
%         end
%         x_= x_./length(ids);
%         x(:,iFile)= nanmean(x_,2);
%         plot(bins,cumsum(x(:,iFile)),'color',[0 0 0 0.1])
%     end
%     ciplot(nanmean(cumsum(x),2)+nansem(cumsum(x),2),nanmean(cumsum(x),2)-nansem(cumsum(x),2),bins,'k')

    x=[];x_ = nan(length(bins),nBS); 
    for iFile = 1:length(fileList)
        for iDraw =1:nBS
            try
                ids = randsample(1:length(D.Units.Nonmembers.SortedScore_all{s}{iFile}),nDrawn);
            catch 
                ids = 1:length(D.Units.Nonmembers.SortedScore_all{s}{iFile});
            end
            x_(:,iDraw)= histc(D.Units.Nonmembers.SortedScore_all{s}{iFile}(ids),bins);
        end
        x_= x_./length(ids);
        x(:,iFile)= nanmean(x_,2);
%         plot(bins,cumsum(x(:,iFile)),'color',[1 0 0 0.1])
    end
    ciplot(nanmean(cumsum(x),2)+nansem(cumsum(x),2),nanmean(cumsum(x),2)-nansem(cumsum(x),2),bins,col_{1},opacity)

    x=[];x_ = nan(length(bins),nBS);
    for iFile = 1:length(fileList)
        for iDraw =1:nBS
            try
                ids = randsample(1:length(D.Units.Members.SortedScore_all{s}{iFile}),nDrawn);
            catch 
                ids = 1:length(D.Units.Members.SortedScore_all{s}{iFile});
            end
            x_(:,iDraw)= histc(D.Units.Members.SortedScore_all{s}{iFile}(ids),bins);
        end
        x_= x_./length(ids);
        x(:,iFile)= nanmean(x_,2);
%         plot(bins,cumsum(x(:,iFile)),'color',[0 0 1 0.1])
    end
    ciplot(nanmean(cumsum(x),2)+nansem(cumsum(x),2),nanmean(cumsum(x),2)-nansem(cumsum(x),2),bins,col_{2},opacity)
    
    x=[];x_ = nan(length(bins),nBS);
    for iFile = 1:length(fileList)
        for iDraw =1:nBS
            try
                ids = randsample(1:length(D.Units.JointMembers.SortedScore_all{s}{iFile}),nDrawn);
            catch 
                ids = 1:length(D.Units.JointMembers.SortedScore_all{s}{iFile});
            end
            x_(:,iDraw)= histc(D.Units.JointMembers.SortedScore_all{s}{iFile}(ids),bins);
        end
        x_= x_./length(ids);
        x(:,iFile)= nanmean(x_,2);
%         plot(bins,cumsum(x(:,iFile)),'color',[0 1 0 0.1])
    end
    ciplot(nanmean(cumsum(x),2)+nansem(cumsum(x),2),nanmean(cumsum(x),2)-nansem(cumsum(x),2),bins,col_{3},opacity)
    
    text(0.48,1,'Chance/Shuffled','HorizontalAlignment','right','Rotation',90)
    plot([0.5 0.5],[0 1],':k')

    
    ylabel('Fraction of units')
    if s==2
       xlabel({'Peak decoding score';'(% Correct Decoding)'})
    end        
    axis([min(bins) max(bins) 0 1])

end
    legend({'Non-members','Members','Joint members','Chance'},'Location','southeast');legend boxoff

%% plot decoding for assemblies
figure
for s=1:3
   subplot(1,3,s); hold on
   scatter(D.AssReal.score{s}(:,1),D.AssReal.score{s}(:,2),'og','filled')
   scatter(D.AssReal.score_shuffled_5pc{s}(:,1),mean([D.AssReal.score_shuffled_5pc{s}(:,2),D.AssReal.score_shuffled_95pc{s}(:,2)],2),'og')
   plot([2 11],[0.5 0.5],':k')
   axis([1 11 0 1])
   if s==1
       ylabel({'Peak decoding score';'(% Correct Decoding)'})
   elseif s==2
       xlabel('Assembly size (no. member units)')
   end
   title([Areas{s} ' assemblies'])
end
legend({'Assemblies','Scrambled trial labels'},'Location','southeast'); legend boxoff
text(11,0.53,'Chance','HorizontalAlignment','right')

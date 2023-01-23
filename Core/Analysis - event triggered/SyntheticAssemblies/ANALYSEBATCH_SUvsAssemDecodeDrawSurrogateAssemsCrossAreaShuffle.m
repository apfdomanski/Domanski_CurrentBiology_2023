%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Target = 'LONG';
Areas = {'PFC','HP','Joint'};
TimeSpan = 4;

bw=0.05;
tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
maxRange = 20;    % Upper limit of assembly size to explore
noTrials = 5;   % How many trials to consider simultanrously
nReps = 50;     % how many repeated draws among trials to perform        
plotOnline = false;
useWholeTaskPeriod = true;
ResampleTrials = true;
InfoCriteria = 'max'; % 'mean','max'
RunBestSU = true;
RunBestEnsemble = true;

clear Av
warning ('off')
if ispc
    
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
%     pat = '/Volumes/HDD2/DNMTP/raw/';
    pat = '/Volumes/Aleks Data Portable/AssemblyAnalysis/raw/';
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
for iFile =6:length(fileList)
    %% Get spike times and rates, behaviour
    fname=strtok(fileList(iFile).name,'_');
    fprintf('Loading file %d/%d: %s...\n',iFile,length(fileList),fname)
    fn = fullfile(pat,fname);
    spikes = load(fn);
    if useWholeTaskPeriod
        fn = fullfile(pat,'KDE_binsTaskOnly',[fname,'_PFC_iFR50_behavOnly.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{1} = temp.iFR;clear temp
        fn = fullfile(pat,'KDE_binsTaskOnly',[fname,'_HP_iFR50_behavOnly.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{2} = temp.iFR;
        Tmtx=temp.Tmtx;clear temp
    else
        fn = fullfile(pat,'KDE_bins',[fname,'_PFC_iFR50.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{1} = temp.iFR;clear temp
        fn = fullfile(pat,'KDE_bins',[fname,'_HP_iFR50.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{2} = temp.iFR;
        Tmtx=temp.Tmtx; clear temp
    end
    iFR{3} = [iFR{1},iFR{2}];
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
    A.nu(3) =sum(A.nu(1:2));
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
    
    for s=1:3
        JointMemberUnits{s} = unique(cell2mat(A.units{s}));
    end
    
    for s=1:3
        globalmemberUnitsReal{s}    = B.usel_out{s}(globalmemberUnits{s});
        globalnonmemberUnitsReal{s} = B.usel_out{s}(nonmemberUnits{s});
        JointMemberUnitsReal{s}     = B.usel_out{s}(JointMemberUnits{s});
    end
    
    %% Aggregate all trials regardless of delay length ** EDIT: use only medium or long trials
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
    %% 1/ Decoding from REAL ASSEMBLIES
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
    %% 2/ Single units ranked by decoding power 
    for s=1:2 
        % NB for s==3, ranked units are drawn from the two individual areas in a first-past-the-post manner
        noUnits     = size(FR{s},2);
        CVE         = zeros(noUnits,Ltr);
        CVEscore    = zeros(noUnits,1);
        CVErank     = zeros(noUnits,1);
        
        parfor iUnit =1:noUnits
            fprintf('Decoding all units: %s unit %d of %d...\n',Areas{s},iUnit,noUnits)
            [CVE(iUnit,:),CVEscore(iUnit)]= RunDecodeGroups(iUnit,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
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
        D_.Units.All.CVE{s}             = CVE_;
        D_.Units.All.SortedScore{s}     = CVEscore_;
        D_.Units.All.Rank{s}            = unitIDs(rank_);
        
        % MEMBER neurons
        unitIDs = globalmemberUnitsReal{s};
        CVE_ = CVE(unitIDs,:);
        CVEscore_ = CVEscore(unitIDs);
        [CVEscore_, rank_] = sort(CVEscore_,'descend');
        CVE_ = CVE_(rank_,:);
        D_.Units.Members.CVE{s}             = CVE_;
        D_.Units.Members.SortedScore{s}     = CVEscore_;
        D_.Units.Members.Rank{s}            = unitIDs(rank_);
        
        % NON-MEMBER neurons
        unitIDs = globalnonmemberUnitsReal{s};
        CVE_ = CVE(unitIDs,:);
        CVEscore_ = CVEscore(unitIDs);
        [CVEscore_, rank_] = sort(CVEscore_,'descend');
        CVE_ = CVE_(rank_,:);
        D_.Units.Nonmembers.CVE{s}             = CVE_;
        D_.Units.Nonmembers.SortedScore{s}     = CVEscore_;
        D_.Units.Nonmembers.Rank{s}            = unitIDs(rank_);
        
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
        D_.Units.JointMembers.CVE{s}             = CVE_;
        D_.Units.JointMembers.SortedScore{s}     = CVEscore_;
        D_.Units.JointMembers.Rank{s}            = unitIDs(rank_);

        clear CVE CVEscore CVEscore_ rank_ CVE_
       
    end
    %% 3/ Decoding from BEST SINGLE UNITS
    if RunBestSU
        % 1) with noise correlations intact
        
        % All neurons
        fprintf('Calculating best single units: All units\n')
        D_.AssBestSUs.All = MakebestSUAssems(D_.Units.All.Rank,D_.Units.All.SortedScore, Ltrials_,Rtrials_,20,false);
        % Cell assembly member neurons
        fprintf('Calculating best single units: Assembly member neurons\n')
        D_.AssBestSUs.Members = MakebestSUAssems(D_.Units.Members.Rank,D_.Units.Members.SortedScore, Ltrials_,Rtrials_,20,false);
        % Cell assembly nonmember neurons
        fprintf('Calculating best single units: Assembly non-member neurons\n')
        D_.AssBestSUs.Nonmembers = MakebestSUAssems(D_.Units.Nonmembers.Rank,D_.Units.Nonmembers.SortedScore, Ltrials_,Rtrials_,20,false);
        % Members of joint-area assemblies
        fprintf('Calculating best single units: Joint area assembly member neurons\n')
        D_.AssBestSUs.JointMembers = MakebestSUAssems(D_.Units.JointMembers.Rank,D_.Units.JointMembers.SortedScore, Ltrials_,Rtrials_,20,false);
        
        % 2) with noise correlations smashed
        
        % All neurons
        fprintf('Calculating best single units: All units without noise correlations\n')
        D_.AssBestSUs.All_noNoiseCorrs = MakebestSUAssems(D_.Units.All.Rank,D_.Units.All.SortedScore, Ltrials_,Rtrials_,20,true);
        % Cell assembly member neurons
        fprintf('Calculating best single units: Assembly member neurons without noise correlations\n')
        D_.AssBestSUs.Members_noNoiseCorrs = MakebestSUAssems(D_.Units.Members.Rank,D_.Units.Members.SortedScore, Ltrials_,Rtrials_,20,true);
        % Cell assembly nonmember neurons
        fprintf('Calculating best single units: Assembly non-member neurons without noise correlations\n')
        D_.AssBestSUs.Nonmembers_noNoiseCorrs = MakebestSUAssems(D_.Units.Nonmembers.Rank,D_.Units.Nonmembers.SortedScore, Ltrials_,Rtrials_,20,true);
        % Members of joint-area assemblies
        %fprintf('Calculating best single units: Joint area assembly member neurons\n')
        %D_.AssBestSUs.JointMembers_noNoiseCorrs = MakebestSUAssems(D_.Units.JointMembers.Rank,D_.Units.JointMembers.SortedScore, Ltrials_,Rtrials_,20,true);
        fprintf('Calculating best single units: Joint area assembly member neurons without noise correlations\n')
        D_.AssBestSUs.JointMembers_noNoiseCorrs = MakebestSUAssemsCrossShuffle(D_.Units.JointMembers.Rank,        ...
            D_.Units.JointMembers.SortedScore, ...
            Ltrials_,Rtrials_,                 ...
            20,                                ...
            [1,1,1]);
        
        if plotOnline
            for s=1:3
                figure('color','w','name',Areas{s}); hold on
                plot((2:length(D_.AssBestSUs.All.Score{s})),D_.AssBestSUs.All.Score{s}(2:end),'k','LineWidth',1.5);
                plot((2:length(D_.AssBestSUs.Members.Score{s})),D_.AssBestSUs.Members.Score{s}(2:end),'b','LineWidth',1.5);
                plot((2:length(D_.AssBestSUs.Nonmembers.Score{s})),D_.AssBestSUs.Nonmembers.Score{s}(2:end),'r','LineWidth',1.5);
                plot((2:length(D_.AssBestSUs.JointMembers.Score{s})),D_.AssBestSUs.JointMembers.Score{s}(2:end),'g','LineWidth',1.5);
                
                ciplot(D_.AssBestSUs.All.Score_shuffled_5pc{s}(2:length(D_.AssBestSUs.All.Score_shuffled_5pc{s})),...
                    D_.AssBestSUs.All.Score_shuffled_95pc{s}(2:length(D_.AssBestSUs.All.Score_shuffled_95pc{s})),...
                    2:length(D_.AssBestSUs.All.Score_shuffled_5pc{s}),'k')
                ciplot(D_.AssBestSUs.Members.Score_shuffled_5pc{s}(2:length(D_.AssBestSUs.Members.Score_shuffled_5pc{s})),...
                    D_.AssBestSUs.Members.Score_shuffled_95pc{s}(2:length(D_.AssBestSUs.Members.Score_shuffled_95pc{s})),...
                    2:length(D_.AssBestSUs.Members.Score_shuffled_5pc{s}),'b')
                ciplot(D_.AssBestSUs.Nonmembers.Score_shuffled_5pc{s}(2:length(D_.AssBestSUs.Nonmembers.Score_shuffled_5pc{s})),...
                    D_.AssBestSUs.Nonmembers.Score_shuffled_95pc{s}(2:length(D_.AssBestSUs.Nonmembers.Score_shuffled_95pc{s})),...
                    2:length(D_.AssBestSUs.Nonmembers.Score_shuffled_5pc{s}),'r')
                if s<3
                    ciplot(D_.AssBestSUs.JointMembers.Score_shuffled_5pc{s}(2:length(D_.AssBestSUs.JointMembers.Score_shuffled_5pc{s})),...
                        D_.AssBestSUs.JointMembers.Score_shuffled_95pc{s}(2:length(D_.AssBestSUs.JointMembers.Score_shuffled_95pc{s})),...
                        2:length(D_.AssBestSUs.JointMembers.Score_shuffled_5pc{s}),'g')
                end
                if s<3
                    %     scatter(1:maxRange,D_.CVEindividualSorted(1:maxRange),'ok','filled');
                    scatter(1:length(D_.Units.Members.SortedScore{s}),D_.Units.Members.SortedScore{s},'ob','filled');
                    scatter(1:length(D_.Units.Nonmembers.SortedScore{s}),D_.Units.Nonmembers.SortedScore{s},'or','filled');
                else
                    x=sort(cell2mat(D_.Units.Members.SortedScore'),'Descend');
                    scatter(1:length(x),x,'ob','filled');
                    x=sort(cell2mat(D_.Units.Nonmembers.SortedScore'),'Descend');
                    scatter(1:length(x),x,'or','filled');
                end
                
                scatter(D_.AssReal.score{s}(:,1),D_.AssReal.score{s}(:,2),'og','filled')
                scatter(D_.AssReal.score_shuffled{s}(:,1),D_.AssReal.score_shuffled{s}(:,2),'og')
                plot([1 maxRange],[0.5 0.5],':k')
                axis([0 maxRange 0.35 1])
                xlabel('Assembly Size (no. Units)')
                ylabel('% Correct decoding')
                legend({'Synthetic Assemblies (Best member units)','Synthetic Assemblies (Best nonmember units)','Shuffled Synthetic Assemblies (members)','Shuffled Synthetic Assemblies (nonmembers)',...
                    'Ranked member units','Ranked nonmember units','Real Assemblies','Real Assemblies (shuffled)','chance'},'Location','eastoutside'); legend boxoff
            end
            
            col_ = gray(10);
            for s=1:3
                figure('color','w','name',Areas{s})
                
                subplot(1,3,1); hold on
                plot((2:length(D_.AssBestSUs.All.Score{s})),D_.AssBestSUs.All.Score{s}(2:end),'k','LineWidth',2);
                for rank_=2:10
                    plot(2:size(D_.AssBestSUs.All.ScoreRanks{s},1),D_.AssBestSUs.All.ScoreRanks{s}(2:end,rank_),'color',col_(rank_,:),'LineWidth',1.5);
                end
                if s<3
                    scatter(1:length(D_.Units.All.SortedScore{s}),D_.Units.All.SortedScore{s},'ob','filled');
                end
                scatter(D_.AssReal.score{s}(:,1),D_.AssReal.score{s}(:,2),'og','filled')
                scatter(D_.AssReal.score_shuffled{s}(:,1),D_.AssReal.score_shuffled{s}(:,2),'og')
                title('All neurons')
                axis([1 20 0 1])
                
                
                subplot(1,3,2); hold on
                plot((2:length(D_.AssBestSUs.Members.Score{s})),D_.AssBestSUs.Members.Score{s}(2:end),'k','LineWidth',2);
                for rank_=2:10
                    plot(2:size(D_.AssBestSUs.Members.ScoreRanks{s},1),D_.AssBestSUs.Members.ScoreRanks{s}(2:end,rank_),'color',col_(rank_,:),'LineWidth',1.5);
                end
                if s<3
                    scatter(1:length(D_.Units.Members.SortedScore{s}),D_.Units.Members.SortedScore{s},'ob','filled');
                end
                scatter(D_.AssReal.score{s}(:,1),D_.AssReal.score{s}(:,2),'og','filled')
                scatter(D_.AssReal.score_shuffled{s}(:,1),D_.AssReal.score_shuffled{s}(:,2),'og')
                title('Assembly members')
                axis([1 20 0 1])
                
                subplot(1,3,3); hold on
                plot((2:length(D_.AssBestSUs.Nonmembers.Score{s})),D_.AssBestSUs.Nonmembers.Score{s}(2:end),'k','LineWidth',2);
                for rank_=2:10
                    plot(2:size(D_.AssBestSUs.Nonmembers.ScoreRanks{s},1),D_.AssBestSUs.Nonmembers.ScoreRanks{s}(2:end,rank_),'color',col_(rank_,:),'LineWidth',1.5);
                end
                if s<3
                    scatter(1:length(D_.Units.Nonmembers.SortedScore{s}),D_.Units.Nonmembers.SortedScore{s},'ob','filled');
                end
                scatter(D_.AssReal.score{s}(:,1),D_.AssReal.score{s}(:,2),'og','filled')
                scatter(D_.AssReal.score_shuffled{s}(:,1),D_.AssReal.score_shuffled{s}(:,2),'og')
                title('Assembly non-members')
                axis([1 20 0 1])
            end
            figure;
            for s=1:3
                subplot(1,3,s); hold on
                %             surf(smooth2a(D_.AssBestSUs.All.ScoreRanks{s},2,2),'EdgeColor','k','FaceColor','k','FaceAlpha',0.3)
                surf(smooth2a(D_.AssBestSUs.Members.ScoreRanks{s},2,2),'EdgeColor','b','FaceColor','b','FaceAlpha',0.5)
                surf(smooth2a(D_.AssBestSUs.Nonmembers.ScoreRanks{s},2,2),'EdgeColor','r','FaceColor','r','FaceAlpha',0.5)
                surf(smooth2a(D_.AssBestSUs.JointMembers.ScoreRanks{s},2,2),'EdgeColor','g','FaceColor','g','FaceAlpha',0.5)
                view(3)
            end
        end
    end
    %% 4/ Decoding from BEST ASSEMBLIES
    if RunBestEnsemble
        % 1) with noise correlations intact

        % All neurons
        fprintf('Calculating best assemblies: All units\n')
        D_.AssBestAssem.All = MakebestEnsembleAssems(B.usel_out,Ltrials_,Rtrials_,20,false,false);
        % Cell assembly member neurons
        fprintf('Calculating best assemblies: Assembly member neurons\n')
        D_.AssBestAssem.Members = MakebestEnsembleAssems(globalmemberUnitsReal,Ltrials_,Rtrials_,20,false,false);
        % Cell assembly nonmember neurons
        fprintf('Calculating best assemblies: Assembly non-member neurons\n')
        D_.AssBestAssem.Nonmembers = MakebestEnsembleAssems(globalnonmemberUnitsReal,Ltrials_,Rtrials_,20,false,false);
        % Members of joint-area assemblies
        fprintf('Calculating best assemblies: Joint area assembly member neurons\n')
        JointOnly = false;
        D_.AssBestAssem.JointMembers = MakebestEnsembleAssems(JointMemberUnitsReal,Ltrials_,Rtrials_,20,JointOnly,false);
        
        % ... NB Last arguement skips the check for units in PFC and HP individually as these are ignored for this
        
        % 2) with noise correlations smashed
        
        % All neurons
        fprintf('Calculating best assemblies: All units without noise correlations\n')
        D_.AssBestAssem.All_noNoiseCorrs = MakebestEnsembleAssems(B.usel_out,Ltrials_,Rtrials_,20,false,true);
        % Cell assembly member neurons
        fprintf('Calculating best assemblies: Assembly member neurons without noise correlations\n')
        D_.AssBestAssem.Members_noNoiseCorrs = MakebestEnsembleAssems(globalmemberUnitsReal,Ltrials_,Rtrials_,20,false,true);
        % Cell assembly nonmember neurons
        fprintf('Calculating best assemblies: Assembly non-member neurons without noise correlations\n')
        D_.AssBestAssem.Nonmembers_noNoiseCorrs = MakebestEnsembleAssems(globalnonmemberUnitsReal,Ltrials_,Rtrials_,20,false,true);
        % Members of joint-area assemblies
        fprintf('Calculating best assemblies: Joint area assembly member neurons without noise correlations\n')
        %     JointOnly = false;
        %     D_.AssBestAssem.JointMembers_noNoiseCorrs = MakebestEnsembleAssems(JointMemberUnitsReal,Ltrials_,Rtrials_,20,JointOnly,true);
        
        % ... NB Last arguement skips the check for units in PFC and HP individually as these are ignored for this
        fprintf('Calculating best assemblies: Joint area assembly member neurons without noise correlations\n')
        D_.AssBestAssem.JointMembers_noNoiseCorrs = MakebestEnsembleAssemsCrossShuffle(JointMemberUnitsReal,        ...
            Ltrials_,Rtrials_,                 ...
            20,                                ...
            [1,1,1]);
    end
    %% 3/Plot results
    if plotOnline
        for s=1:3
            
            figure('color','w','name',Areas{s}); hold on
            %     plot(1:maxRange,D_.CVEBestSU(1:maxRange),'k');
%             plot((2:length(D_.AssBestSUs.Members.Score{s})),D_.AssBestSUs.Members.Score{s}(2:end),'b','LineWidth',1.5);
%             plot((2:length(D_.AssBestSUs.Nonmembers.Score{s})),D_.AssBestSUs.Nonmembers.Score{s}(2:end),'r','LineWidth',1.5);
            
%             plot((2:length(D_.AssBestAssem.Members.Score{s})),D_.AssBestAssem.Members.Score{s}(2:end),':b','LineWidth',1.5);
            plot((2:length(D_.AssBestAssem.Nonmembers.Score{s})),D_.AssBestAssem.Nonmembers.Score{s}(2:end),':r','LineWidth',1.5);
            
            %     plot(1:maxRange,D_.CVEBestSU_shuffled(1:maxRange),'r');
            
            %     ciplot(D_.AssBestSUs.Members.Score_shuffled_5pc{s}(2:length(D_.AssBestSUs.Members.Score_shuffled_5pc{s})),...
            %            D_.AssBestSUs.Members.Score_shuffled_95pc{s}(2:length(D_.AssBestSUs.Members.Score_shuffled_95pc{s})),...
            %            2:length(D_.AssBestSUs.Members.Score_shuffled_5pc{s}),'b')
            %     ciplot(D_.AssBestSUs.Nonmembers.Score_shuffled_5pc{s}(2:length(D_.AssBestSUs.Nonmembers.Score_shuffled_5pc{s})),...
            %            D_.AssBestSUs.Nonmembers.Score_shuffled_95pc{s}(2:length(D_.AssBestSUs.Nonmembers.Score_shuffled_95pc{s})),...
            %            2:length(D_.AssBestSUs.Nonmembers.Score_shuffled_5pc{s}),'r')
            %        if s<3
            %            %     scatter(1:maxRange,D_.CVEindividualSorted(1:maxRange),'ok','filled');
            %            scatter(1:length(D_.Units.Members.SortedScore{s}),D_.Units.Members.SortedScore{s},'ob','filled');
            %            scatter(1:length(D_.Units.Nonmembers.SortedScore{s}),D_.Units.Nonmembers.SortedScore{s},'or','filled');
            %        else
            %            x=sort(cell2mat(D_.Units.Members.SortedScore'),'Descend');
            %            scatter(1:length(x),x,'ob','filled');
            %            x=sort(cell2mat(D_.Units.Nonmembers.SortedScore'),'Descend');
            %            scatter(1:length(x),x,'or','filled');
            %        end
            %
            scatter(D_.AssReal.score{s}(:,1),D_.AssReal.score{s}(:,2),'og','filled')
            scatter(D_.AssReal.score_shuffled{s}(:,1),D_.AssReal.score_shuffled{s}(:,2),'og')
            %     plot([1 maxRange],[0.5 0.5],':k')
            %     axis([0 maxRange 0.35 1])
            %     xlabel('Assembly Size (no. Units)')
            %     ylabel('% Correct decoding')
            %     legend({'Synthetic Assemblies (Best member units)','Synthetic Assemblies (Best nonmember units)','Shuffled Synthetic Assemblies (members)','Shuffled Synthetic Assemblies (nonmembers)',...
            %             'Ranked member units','Ranked nonmember units','Real Assemblies','Real Assemblies (shuffled)','chance'},'Location','eastoutside'); legend boxoff
        end
    end
    %% (4) Save results
    fnOut = sprintf('%sSyntheticAssemblies%s%s_AssemblyComparison_redux4.mat',pat,filesep,fname);
    save(fnOut,'D_','A','B','Areas','maxRange','InfoCriteria','tlimsAll','tbAll','Delays_','Ltrials_','Rtrials_','-v7.3')
    clear D_
end










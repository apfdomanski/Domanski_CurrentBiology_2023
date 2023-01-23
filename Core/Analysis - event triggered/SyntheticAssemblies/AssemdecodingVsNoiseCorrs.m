%% %%%%%% PREAMBLE %%%%%%
Delays_ = {'Short','Medium','Long'};
Target = 'LONG';
Areas = {'PFC','HP','Joint'};
TimeSpan = 4;

bw=0.05;
tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
maxRange = 20;       % Upper limit of assembly size to explore
noTrials = Inf;      % How many trials to consider simultaneously
nReps    = 1;        % how many repeated draws among trials to perform
nRepsNoise = 10;
plotOnline = false;
useWholeTaskPeriod = true;
if isinf(noTrials)
    ResampleTrials = false;
else
    ResampleTrials = true;
end
InfoCriteria    = 'max'; % 'mean','max'
RunBestSU       = true;
RunBestEnsemble = true;
TrialsToUse     = 'Error'; %'Correct','Error'

clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
    %     cd (['~/MATLAB/MichalDataAna'])
    %     pat = ['~/raw/'];
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
for iFile =1:length(fileList)
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
        switch TrialsToUse
            case 'Correct'
                eval(sprintf('LeftTrials = [LeftTrials;[t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect]];',Delays_{iDelay},Delays_{iDelay}));
                eval(sprintf('RightTrials = [RightTrials;[t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect]];',Delays_{iDelay},Delays_{iDelay}));
            case 'Error'
                eval(sprintf('LeftTrials = [LeftTrials;[t.%s.SamplePress_LeftError'',t.%s.ChoicePress_LeftError'']];',Delays_{iDelay},Delays_{iDelay}));
                eval(sprintf('RightTrials = [RightTrials;[t.%s.SamplePress_RightError'',t.%s.ChoicePress_RightError'']];',Delays_{iDelay},Delays_{iDelay}));
        end
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
                
                [D_.AssReal.CVE{s}(:,iAss),...
                    score,...
                    score_shuffled,...
                    score_shuffled_5pc,...
                    score_shuffled_95pc] = RunDecodeGroups(unitIds,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
                
                D_.AssReal.score{s}(iAss,:)                    =  [noUnits,score];
                D_.AssReal.score_shuffled{s}(iAss,:)           =  [noUnits,score_shuffled];
                D_.AssReal.score_shuffled_5pc{s}(iAss,:)       =  [noUnits,score_shuffled_5pc];
                D_.AssReal.score_shuffled_95pc{s}(iAss,:)      =  [noUnits,score_shuffled_95pc];
                clear score score_shuffled score_shuffled_5pc score_shuffled_95pc
                 
                [D_.AssReal_noNoiseCorrs.CVE{s}(:,iAss),...
                    score,...
                    score_shuffled,...
                    score_shuffled_5pc,...
                    score_shuffled_95pc] = RunDecodeGroupsDestroyNoiseCorrs(unitIds,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
                D_.AssReal_noNoiseCorrs.score{s}(iAss,:)                    =  [noUnits,score];
                D_.AssReal_noNoiseCorrs.score_shuffled{s}(iAss,:)           =  [noUnits,score_shuffled];
                D_.AssReal_noNoiseCorrs.score_shuffled_5pc{s}(iAss,:)       =  [noUnits,score_shuffled_5pc];
                D_.AssReal_noNoiseCorrs.score_shuffled_95pc{s}(iAss,:)      =  [noUnits,score_shuffled_95pc];
            
            clear score score_shuffled score_shuffled_5pc score_shuffled_95pc
            if s==3
                for s_=1:3
                    [D_.AssReal_noNoiseCorrs.CVE_crossShuffle{s}{s_}(:,iAss),...
                     score] = RunDecodeGroupsDestroyNoiseCorrsCrossArea(unitIds,10,Ltrials_,Rtrials_,noTrials,'max',s_);
                     D_.AssReal_noNoiseCorrs.score_crossShuffle{s}{s_}(iAss,:)                    =  [noUnits,score];
                     clear score 
                end
            end
            end
        end
    end

    %% (4) Save results
    fnOut = sprintf('%sSyntheticAssemblies%s%s_AssDecodingVsNoiseCorrs.mat',pat,filesep,fname);
    save(fnOut,'D_','A','B','Areas','maxRange','InfoCriteria','tlimsAll','tbAll','Delays_','Ltrials_','Rtrials_','-v7.3')
    clear D_
end










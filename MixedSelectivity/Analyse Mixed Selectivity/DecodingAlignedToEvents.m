%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Target = 'LONG';


runCVEdecoder = false;
ResampleTS    = false;
ResampleCVE   = false;
minTrialCount = 5;

bw=0.05;
Nbs = 10;
Nbs_TScore = 10;

tlimsAll = [-2 5];
tlimsNormWin = [-2 -0.2];

tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
tlims_  = closest(tbAll,tlimsNormWin);
IdxNorm = tlims_(1):tlims_(2); clear tlims_
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
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);
fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileListAss(find(name_flag))=[];% fnames_(name_flag)=[];
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag
pat2 = 'C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\LONGTaskonly';
normaliseFscores = false;
normWin = [-5 -3];
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};

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

mkdir([pat 'MixedSelectivity_LONG'])
MemberClasses ={'LocalMembers','JointMembers','NonMembers'};
Outcome = {'Correct','Error'};
Events_ = {'CueLight','SamplePress','DelayEnd','NosePoke','ChoicePress','RewardConsume'};
Events__ = {'Cue Light','Sample Press','Delay End','Nose Poke','Choice Press','Reward Collection'};
no_Events_ = [6,5];
%%%%%
%% Batch process units
for iFile = 1:length(fileList)
    %% Get the single unit files
    fname=strtok(fileList(iFile).name,'_');
    
    for iArea = 1:2%length(Areas)
        fprintf('Getting units: %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
        load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
        
        iFR_{iArea} = iFR;
%         iFR_{iArea} = zscore(iFR);
    end
    fprintf('Getting units: %d/%d %s (Joint)...\n',iFile,length(fileList),fname)
    iFR_{3} = [iFR_{1},iFR_{2}];
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    noPFC = length(PFCcells);
    %% Get the assembly membership info
    fname=strtok(fileListAss(iFile).name,'_');
    fprintf('Loading run %d/%d %s ...\n',iFile,length(fileList),fname)
    % for assemblies describing whole task period (cue  - reward)
    switch AssemblyChoice
        case 1
            load(sprintf('%s%s%s_iFR50_FSC.mat',pat2,filesep,fname),'units','nu');
            %             load(sprintf('%s%s%s_iFR50_AssemRes2.mat',pat2,filesep,fname),'usel_out');
            usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
        case 2
            load(sprintf('%s%s%s_iFR50_BehavOnly_FSC.mat',pat2,filesep,fname),'units','nu');
            load(sprintf('%s%s%s_iFR50_BehavOnly_AssemRes2.mat',pat2,filesep,fname),'usel_out');
        case 3
            load(sprintf('%s%s%s_iFR50_Task_FSC.mat',pat2,filesep,fname),'units','nu');
            %             load(sprintf('%s%s%s_iFR50_Task_AssemRes2.mat',pat2,filesep,fname),'usel_out');
            usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
            
    end
    Ass.usel_out = usel_out;
    Ass.usel_out{3} = [Ass.usel_out{1},Ass.usel_out{2}+noPFC];
    Ass.nu = nu; Ass.nu(3) = sum(Ass.nu(1:2));
    Ass.units = units;
    for iArea = 1:3
        Ass.LocalMembers{iArea}    = unique(cell2mat(Ass.units{iArea}));
    end
    
    
    Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3}<=Ass.nu(1)));
    Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3}>Ass.nu(1))-Ass.nu(1));
    for iArea = 1:2
        Ass.JointMembers{iArea} = setdiff(unique(cell2mat(Ass.units{iArea})),Ass.LocalMembers{iArea});
        Ass.NonMembers{iArea} = setdiff(1:Ass.nu(iArea),[Ass.LocalMembers{iArea},Ass.JointMembers{iArea}]);
    end
    for iArea = 1:2
        Ass.NonMembers{iArea}   = Ass.usel_out{iArea}(Ass.NonMembers{iArea});
        Ass.LocalMembers{iArea} = Ass.usel_out{iArea}(Ass.LocalMembers{iArea});
        Ass.JointMembers{iArea} = Ass.usel_out{iArea}(Ass.JointMembers{iArea});
    end
    
    clear units usel_out nu
    %% Get event times
    clear LeftTrials RightTrials
    for iOutcome =1:2
        for iDelay =1:length(Delays_)
            % Get event times
            noEvents = no_Events_(iOutcome);
            for iEvent = 1:noEvents
                eval(sprintf('LeftTrials{iOutcome}{iEvent}{iDelay}  =  t.%s.%s_Left%s(:);',Delays_{iDelay},Events_{iEvent},Outcome{iOutcome}));
                eval(sprintf('RightTrials{iOutcome}{iEvent}{iDelay} =  t.%s.%s_Right%s(:);',Delays_{iDelay},Events_{iEvent},Outcome{iOutcome}));
            end
        end
    end
    %% Batch Process
    nu = [];
    for iArea =1:2
        for iClass =1:length(MemberClasses)
            eval(sprintf('idx = Ass.%s{iArea};',MemberClasses{iClass}));
            nu = [nu;length(idx)];
        end
    end
    nu(nu==0)=[];
    
    for iOutcome =1:2
        for iDelay =1:length(Delays_)
            for iEvent = 1:no_Events_(iOutcome)
                for iArea = 1:2
                    for iClass =1:length(MemberClasses)
                        eval(sprintf('idx = Ass.%s{iArea};',MemberClasses{iClass}))
                        if ~isempty(idx)
                            
                            LeftTrials_ = LeftTrials{iOutcome}{iEvent}{iDelay};   LeftTrials_(isnan(LeftTrials_))=[];
                            RightTrials_ = RightTrials{iOutcome}{iEvent}{iDelay}; RightTrials_(isnan(RightTrials_))=[];
                            
                            Ltrials = [];nL = 0;Ltrials_ = {};
                            for iUnit = 1:length(idx)
                                Ltrials__{iUnit}=[];
                            end
                            for iTrial =1:size(LeftTrials_,1)
                                try
                                    tlims_  = LeftTrials_(iTrial)/1e6+tlimsAll(1);
                                    tlims_  = closest(Tmtx,tlims_);
                                    tlims_  = [tlims_:(tlims_+length(tbAll)-1)];
                                    
                                    Ltrials = [Ltrials;iFR_{iArea}(tlims_,idx)];
                                    Ltrials_{1,nL+1} = iFR_{iArea}(tlims_,idx);
                                    for iUnit = 1:length(idx)
                                        Ltrials__{iUnit} = [Ltrials__{iUnit}, iFR_{iArea}(tlims_,idx(iUnit))];
                                    end
                                    nL=nL+1;
                                end
                            end
                            
                            Rtrials = [];nR = 0; Rtrials_ = {};
                            for iUnit = 1:length(idx)
                                Rtrials__{iUnit}=[];
                            end
                            for iTrial =1:size(RightTrials_,1)
                                try
                                    tlims_  = RightTrials_(iTrial)/1e6+tlimsAll(1);
                                    tlims_  = closest(Tmtx,tlims_);
                                    tlims_  = [tlims_:(tlims_+length(tbAll)-1)];
                                    Rtrials = [Rtrials;iFR_{iArea}(tlims_,idx)];
                                    Rtrials_{1,nR+1} = iFR_{iArea}(tlims_,idx);
                                    for iUnit = 1:length(idx)
                                        Rtrials__{iUnit} = [Rtrials__{iUnit}, iFR_{iArea}(tlims_,idx(iUnit))];
                                    end
                                    nR=nR+1;
                                end
                            end
                            
                            if nL >= minTrialCount && nR >= minTrialCount
                                %% Decode L vs. R
                                FR = [Ltrials;Rtrials];
                                evt0 = [ones(nL,1);2*ones(nR,1)];
                                Ltr = length(tbAll);
                                
                                % (1) Multivariate F-score decoder
                                fprintf('Decoding file %d/%d (%s) %s %s , %s %s trials: %s, F-decoder...\n',iFile,length(fileList),fname,Areas{iArea},MemberClasses{iClass},Outcome{iOutcome}, Delays_{iDelay},Events__{iClass})
                                
                                [D_.avgFR,D_.seFR,...
                                    D_.Ft2,D_.Rt2,...
                                    D_.Ft2ciL,D_.Ft2ciH,...
                                    D_.TS,...
                                    D_.dfnum,D_.dfden] = DecodeStats(FR,evt0,0.05);
                                D_.TSsig = tpdf(D_.TS,nL+nR-2)<0.05;
                                
                                if ResampleTS
                                    for nrep = 1:Nbs_TScore
                                        fprintf('Decoding file %d/%d (%s) %s %s , %s %s trials, F-decoder... bs %d/%d \n',iFile,length(fileList),fname,Areas{iArea},MemberClasses{iClass},Outcome{iOutcome}, Delays_{iDelay},nrep,Nbs_TScore)
                                        
                                        
                                        rsL = randsample(length(Ltrials_),minTrialCount);
                                        rsR = randsample(length(Rtrials_),minTrialCount);
                                        
                                        FR   = [cell2mat(Ltrials_(rsL)');cell2mat(Rtrials_(rsR)')];
                                        evt0 = [ones(minTrialCount,1);2*ones(minTrialCount,1)];
                                        
                                        % Multivariate F-score decoder
                                        [E.avgFR(:,:,nrep),E.seFR(:,:,nrep),...
                                            E.Ft2(:,nrep),E.Rt2(:,nrep),...
                                            E.Ft2ciL(:,nrep),E.Ft2ciH(:,nrep),...
                                            E.TS(:,:,nrep),...
                                            ~,~] = DecodeStats(FR,evt0,0.05);
                                        E.TSsig(:,:,nrep) = tpdf(E.TS(:,:,nrep),2*minTrialCount-2)<0.05;
                                    end
                                    D_.TS_drawnTrials     =  mean(E.TS,3);
                                    D_.TSsig_drawnTrials  =  sum(E.TSsig,3)./Nbs_TScore;
                                    D_.Ft2_drawnTrials    =  sum(E.Ft2,2)./Nbs_TScore;
                                    D_.Rt2_drawnTrials    =  sum(E.Rt2,2)./Nbs_TScore;
                                    D_.Rt2_drawnTrials    =  sum(E.Rt2,2)./Nbs_TScore;
                                    D_.Ft2ciH_drawnTrials =  sum(E.Ft2ciH,2)./Nbs_TScore;
                                    clear E
                                end
                                
                                % (2) Cross-validation error decoder
                                if runCVEdecoder
                                    fprintf('Decoding file %d/%d (%s) %s %s , %s %s trials: %s, F-decoder... CVE \n',iFile,length(fileList),fname,Areas{iArea},MemberClasses{iClass},Outcome{iOutcome}, Delays_{iDelay},Events__{iClass})
                                    
                                    for iUnit =1:size(FR,2)
                                        D_.CVEindividual(iUnit,:) = DecodeCVE_mex(FR(:,iUnit),evt0,0.05);
                                        %                                     D_.CVEindividual(iUnit,:) = DecodeCVE(FR(:,iUnit),evt0,0.05);
                                    end
                                    if size(Ltrials,2)==min(nu), nrep=1; else nrep = 10; end;
                                    %                  nrep=1;
                                    for i=1:nrep
                                        try
                                            rs{i} = randsample(size(Ltrials,2),3);
                                            %                                             rs{i} = randsample(size(Ltrials,2),min(nu));
                                        catch
                                            rs{i} = 1:size(Ltrials,2);
                                        end
                                        %                     rs{i} = randsample(size(Ltrials,2),nu(iArea));
                                        %                     D_.CVE(i,:) = DecodeCVE(FR(:,rs{i}),evt0,0.05);
                                        %                         D_.CVE(i,:) = DecodeCVE_mex(FR(:,rs{i}),evt0,0.05);
                                        D_.CVE(i,:) = DecodeCVE(FR(:,rs{i}),evt0,0.05);
                                    end
                                    D_.CVE = nanmean(D_.CVE,1);
                                    if ResampleCVE
                                        D_.cveBSciH=zeros(1,Ltr);
                                        D_.cveBSciL=zeros(1,Ltr);
                                        D_.CVEbs=[];
                                        CVEbs = zeros(Nbs,Ltr);
                                        parfor b=1:Nbs
                                            t1=tic;
                                            k=randperm(length(evt0)); evt1=evt0(k);
                                            cve0=zeros(nrep,Ltr);
                                            %                                     for i=1:nrep, cve0(i,:)=DecodeCVE(FR,evt1,0.05); end
                                            for i=1:nrep, cve0(i,:)=DecodeCVE_mex(FR(:,rs{i}),evt1,0.05); end
                                            %                                     for i=1:nrep, cve0(i,:)=DecodeCVE    (FR(:,rs{i}),evt1,0.05); end
                                            CVEbs(b,:) = nanmean(cve0,1);
                                            fprintf('Decoding file %d/%d (%s) %s %s , %s %s trials: %s, F-decoder... CVE draw %d of %d (%s elapsed)... \n',iFile,length(fileList),fname,Areas{iArea},MemberClasses{iClass},Outcome{iOutcome}, Delays_{iDelay},Events__{iClass},b,Nbs,seconds2human(toc(t1)))
                                            
                                        end
                                        for tr=1:Ltr
                                            cves=sort(CVEbs(:,tr),'ascend');
                                            D_.cveBSciH(tr)=cves(round(0.95*Nbs));
                                            D_.cveBSciL(tr)=cves(round(0.05*Nbs));
                                        end
                                        D_.CVEbs = nanmean(CVEbs,1);
                                    end
                                end
                                clear rs evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
                                %% Mean Firing rate and Fano
                                fprintf('Decoding file %d/%d (%s) %s %s , %s %s trials: %s, Fano Factor...\n',iFile,length(fileList),fname,Areas{iArea},MemberClasses{iClass},Outcome{iOutcome}, Delays_{iDelay},Events__{iClass})

                                Lmean  = cellfun(@(x) nanmean(x,2),Ltrials__,'UniformOutput',false);
                                Rmean  = cellfun(@(x) nanmean(x,2),Rtrials__,'UniformOutput',false);
                                Lvar   = cellfun(@(x) nanvar(x,[],2),Ltrials__,'UniformOutput',false);
                                Rvar   = cellfun(@(x) nanvar(x,[],2),Rtrials__,'UniformOutput',false);
                                
                                LFano     = cellfun(@(x,y) x./y,Lvar,Lmean,'UniformOutput',false);
                                RFano     = cellfun(@(x,y) x./y,Rvar,Rmean,'UniformOutput',false);
                                Fano_     = cellfun(@(x,y)(x+y)./2,LFano,RFano,'UniformOutput',false);
                                
%                                 Lbaseline = cellfun(@(x) repmat(nanmean(x(IdxNorm,:)),length(tbAll),1),Ltrials__,'UniformOutput',false);
%                                 Lmean     = cellfun(@(x,y) ((x-y)./y), Ltrials__,Lbaseline,'UniformOutput',false);
%                                 Rbaseline = cellfun(@(x)   repmat(nanmean(x(IdxNorm,:)),length(tbAll),1),Rtrials__,'UniformOutput',false);
%                                 Rmean     = cellfun(@(x,y) ((x-y)./y), Rtrials__,Rbaseline,'UniformOutput',false);
%                                 Mean_     = cellfun(@(x,y) (nanmean(x,2) + nanmean(y,2))./2, Lmean,Rmean,'UniformOutput',false);
                                
                                Mean_     = cellfun(@(x,y) (nanmean(x,2) + nanmean(y,2))./2, Ltrials__,Rtrials__,'UniformOutput',false);


                                clear Lmean Rmean Lvar Rvar LFano RFano Lbaseline Rbaseline
                                D_.Mean = cell2mat(Mean_);
                                D_.Fano = cell2mat(Fano_);
                                %% bag 'em
                                D{iOutcome}{iEvent}{iDelay}{iArea}{iClass} = D_;
                                
                                
                            end
                        end
                        clear D_ Ltrials nL Rtrials nR iTrial
                    end
                end
            end
        end
    end
    %% Save results
    if exist('D')
        fnOut = sprintf('%sMixedSelectivity_LONG%s%s_MixedSelectivity_MembershipsortedUnitsEventAnalysis_Fano.mat',pat,filesep,fname);
        save(fnOut,'D','tbAll','-v7.3');
        fprintf('Done.\n')
        clear D
    end
end
clearvars -except Nbs AssemblyChoice color_ Delays_ bw tbShort tbAll tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin

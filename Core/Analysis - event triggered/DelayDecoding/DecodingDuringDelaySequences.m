% function DecodingDuringDelayIndividualUnits(iFile)
% if isstr(iFile)
%     iFile = str2num(iFile);
% end
%% %%%%%% PREAMBLE %%%%%%

Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};

Target = 'LONG';
RunErrors = false;
tlimsAll = [0 5];
tlimsShort=[0 4];
tlimsMedium=[0 8];
tlimsLong=[0 16];

shift = 0;
plotOnline = false;
bw    = 0.05;
Nbs = 500;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

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

normaliseFscores = false;
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
        pat2 = [pat 'KDE_binsTaskOnly' filesep 'LONGTaskonly' filesep];
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end
mkdir([pat 'MixedSelectivity'])

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

TrialsSubsample = Inf;
if isinf(TrialsSubsample)
    nReps = 1;
else
    nReps =500;
end
%% batch process...
for iFile =1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
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
    noPFC=nu(1);
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
    
    clear units usel_out
    %% Batch process behavioural performance
    L = [length(t.Short.ChoicePress_LeftCorrect)./(length(t.Short.ChoicePress_LeftCorrect)+length(t.Short.ChoicePress_LeftError)) , ...
        length(t.Medium.ChoicePress_LeftCorrect)./(length(t.Medium.ChoicePress_LeftCorrect)+length(t.Medium.ChoicePress_LeftError)) , ...
        length(t.Long.ChoicePress_LeftCorrect)./(length(t.Long.ChoicePress_LeftCorrect)+length(t.Long.ChoicePress_LeftError))];
    
    
    R = [length(t.Short.ChoicePress_RightCorrect)./(length(t.Short.ChoicePress_RightCorrect)+length(t.Short.ChoicePress_RightError)) , ...
        length(t.Medium.ChoicePress_RightCorrect)./(length(t.Medium.ChoicePress_RightCorrect)+length(t.Medium.ChoicePress_RightError)) , ...
        length(t.Long.ChoicePress_RightCorrect)./(length(t.Long.ChoicePress_RightCorrect)+length(t.Long.ChoicePress_RightError))];
    
    D.Accuracy = (L+R)./2;
    
    C_ =  [length(t.Short.ChoicePress_LeftCorrect)  + length(t.Short.ChoicePress_RightCorrect),...
        length(t.Medium.ChoicePress_LeftCorrect) + length(t.Medium.ChoicePress_RightCorrect),...
        length(t.Long.ChoicePress_LeftCorrect)   + length(t.Long.ChoicePress_RightCorrect)];
    
    E_ =  [length(t.Short.ChoicePress_LeftError)    + length(t.Short.ChoicePress_RightError),...
        length(t.Medium.ChoicePress_LeftError)   + length(t.Medium.ChoicePress_RightError),...
        length(t.Long.ChoicePress_LeftError)     + length(t.Long.ChoicePress_RightError)];
    
    D.pCorr = C_./(C_+E_);
    D.AboveChance = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;
    clear L R C_ E_
    
end

%% batch process units...
for iFile =1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
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
    noPFC=nu(1);
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
    
    clear units usel_out
    %% Batch process units
    for iArea = 1%:length(Areas)
        %% Get the files
        fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
        if iArea < 3
            % load(sprintf('%sKDE_bins%s%s_%s_iFR50.mat',pat,filesep,fname,Areas{iArea}));
            load(sprintf('%sKDE_binsTaskOnly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
            iFR_ = iFR;
        else
            U_{1} = load(sprintf('%sKDE_binsTaskOnly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{1}));
            U_{2} = load(sprintf('%sKDE_binsTaskOnly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{2}));
            iFR_ = [U_{1}.iFR,U_{2}.iFR];
            Tmtx = U_{1}.Tmtx;
            clear U_
        end
        
        % iFR_ = zscore(iFR_);
        %% Delay period decoding for units
        shuffleCorrect = true;
        for iDelay =1:length(Delays_)
            
            eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftCorrect];',Delays_{iDelay}));
            eval(sprintf('RightTrials = [t.%s.SamplePress_RightCorrect];',Delays_{iDelay}));
            eval(sprintf('LeftTrialsE = [t.%s.SamplePress_LeftError''];',Delays_{iDelay}));
            eval(sprintf('RightTrialsE = [t.%s.SamplePress_RightError''];',Delays_{iDelay}));
            eval(sprintf('tlims_X = tlims%s;',Delays_{iDelay}));
            length_ = sum(abs(tlims_X))/bw;
            eval(sprintf('tb_ = tb%s;',Delays_{iDelay}))
            
            Ltrials = [];nL = 0; Ltrials_={};
            for iTrial =1:size(LeftTrials,1)
                try
                    tlims_  = LeftTrials(iTrial,1)/1e6+tlims_X(1);
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    Ltrials = [Ltrials;iFR_(tlims_,:)];
                    Ltrials_{iTrial} = iFR_(tlims_,:);
                    nL=nL+1;
                end
            end
            Rtrials = [];nR = 0;  Rtrials_={};
            for iTrial =1:size(RightTrials,1)
                try
                    tlims_  = RightTrials(iTrial,1)/1e6+tlims_X(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    Rtrials = [Rtrials;iFR_(tlims_,:)];
                    Rtrials_{iTrial} = iFR_(tlims_,:);
                    
                    nR=nR+1;
                end
            end
            LtrialsE = [];nLe = 0; LtrialsE_={};
            for iTrial =1:size(LeftTrialsE,1)
                try
                    tlims_  = LeftTrialsE(iTrial,1)/1e6+tlims_X(1) + shift;
                    tlims_ = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    LtrialsE = [LtrialsE;iFR_(tlims_,:)];
                    LtrialsE_{iTrial} = iFR_(tlims_,:);
                    
                    nLe=nLe+1;
                end
            end
            RtrialsE = [];nRe = 0; RtrialsE_={};
            for iTrial =1:size(RightTrialsE,1)
                try
                    tlims_  = RightTrialsE(iTrial,1)/1e6+tlims_X(1) + shift;
                    tlims_ = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    RtrialsE = [RtrialsE;iFR_(tlims_,:)];
                    RtrialsE_{iTrial} = iFR_(tlims_,:);
                    nRe=nRe+1;
                end
            end
            flag_ = 1;
            if ~isinf(TrialsSubsample)
                if sum([nL,nR]>=TrialsSubsample)==2;
                    flag_ = 0;
                end
            end
            if flag_
                %%%
                % Sort by peak mean firing rate time
                %                 FR_   = [cell2mat(Ltrials_');cell2mat(Rtrials_')];
                % %                 FR_ = zscore(FR_);
                %
                %                 evt0  = [ones(length(Ltrials_),1);2*ones(length(Rtrials_),1)];
                %                 [~,meanFR_]  = DecodeStatsTSonly(FR_,evt0,0.05);
                %                 meanFR_ = (meanFR_{1}+meanFR_{2})./2;
                %                 [~,idx] = max(meanFR_,[],1);
                %                 [~,idx] = sort(idx);
                %                 tmax = size(meanFR_,1);
                nu_ = length(Ass.NonMembers{iArea});
                idx = 1:nu_;
                tmax = size(Ltrials_{1},1);
                zscoreYN = true;
                nBS=500;
                maxLag = 16/bw;
                xcL  = nan(nu_,nu_,nL);
                xcLE = nan(nu_,nu_,nLe);
                xcR  = nan(nu_,nu_,nR);
                xcRE = nan(nu_,nu_,nRe);
                if shuffleCorrect
                    xcLshufH = xcL;
                    xcLshufL = xcL;
                    xcLsig   = xcL;
                    
                    xcLEshufH = xcLE;
                    xcLEshufL = xcLE;
                    xcLEsig   = xcLE;
                    
                    xcRshufH = xcR;
                    xcRshufL = xcR;
                    xcRsig   = xcR;
                    
                    xcREshufH = xcRE;
                    xcREshufL = xcRE;
                    xcREsig   = xcRE;
                end
                
                for iTrial =1:nL
                    FR_ = Ltrials_{iTrial}(:,Ass.NonMembers{iArea});
                    if zscoreYN
                        FR_ = zscore(FR_);
                    end
                    for ii=1:nu_
                        [iTrial ii]
                        parfor jj=ii+1:nu_
                            % Process Correct
                            xcL(ii,jj,iTrial) = max(xcorr(FR_(1:tmax,idx(ii)),FR_(1:tmax,idx(jj)),maxLag,'unbiased'));
                            if shuffleCorrect
                                xcshuf  = nan(1,nBS);
                                for iDraw=1:nBS
                                    FR__ = FR_(randperm(size(FR_,1)),idx([ii,jj]));
                                    xcshuf(iDraw) = max(xcorr(FR__(:,1),FR__(:,2),maxLag,'unbiased'));
                                end
                                xcLshufL(ii,jj,iTrial)  = prctile(xcshuf,5);
                                xcLshufH(ii,jj,iTrial)  = prctile(xcshuf,95);
                            end
                            
                        end
                    end
                end
                for iTrial =1:nLe
                    FR_ = LtrialsE_{iTrial}(:,Ass.NonMembers{iArea});
                    if zscoreYN
                        FR_ = zscore(FR_);
                    end
                    for ii=1:nu_
                        [iTrial ii]
                        parfor jj=ii+1:nu_
                            % Process Correct
                            xcLE(ii,jj,iTrial) = max(xcorr(FR_(1:tmax,idx(ii)),FR_(1:tmax,idx(jj)),maxLag,'unbiased'));
                            if shuffleCorrect
                                xcshuf  = nan(1,nBS);
                                for iDraw=1:nBS
                                    FR__ = FR_(randperm(size(FR_,1)),idx([ii,jj]));
                                    xcshuf(iDraw) = max(xcorr(FR__(:,1),FR__(:,2),maxLag,'unbiased'));
                                end
                                xcLEshufL(ii,jj,iTrial)  = prctile(xcshuf,5);
                                xcLEshufH(ii,jj,iTrial)  = prctile(xcshuf,95);
                            end
                            
                        end
                    end
                end
                
                for iTrial =1:nR
                    FR_ = Rtrials_{iTrial}(:,Ass.NonMembers{iArea});
                    if zscoreYN
                        FR_ = zscore(FR_);
                    end
                    for ii=1:nu_
                        [iTrial ii]
                        parfor jj=ii+1:nu_
                            % Process Correct
                            xcR(ii,jj,iTrial) = max(xcorr(FR_(1:tmax,idx(ii)),FR_(1:tmax,idx(jj)),maxLag,'unbiased'));
                            if shuffleCorrect
                                xcshuf  = nan(1,nBS);
                                for iDraw=1:nBS
                                    FR__ = FR_(randperm(size(FR_,1)),idx([ii,jj]));
                                    xcshuf(iDraw) = max(xcorr(FR__(:,1),FR__(:,2),maxLag,'unbiased'));
                                end
                                xcRshufL(ii,jj,iTrial)  = prctile(xcshuf,5);
                                xcRshufH(ii,jj,iTrial)  = prctile(xcshuf,95);
                            end
                            
                        end
                    end
                end
                for iTrial =1:nRe
                    FR_ = RtrialsE_{iTrial}(:,Ass.NonMembers{iArea});
                    if zscoreYN
                        FR_ = zscore(FR_);
                    end
                    for ii=1:nu_
                        [iTrial ii]
                        parfor jj=ii+1:nu_
                            % Process Correct
                            xcRE(ii,jj,iTrial) = max(xcorr(FR_(1:tmax,idx(ii)),FR_(1:tmax,idx(jj)),maxLag,'unbiased'));
                            if shuffleCorrect
                                xcshuf  = nan(1,nBS);
                                for iDraw=1:nBS
                                    FR__ = FR_(randperm(size(FR_,1)),idx([ii,jj]));
                                    xcshuf(iDraw) = max(xcorr(FR__(:,1),FR__(:,2),maxLag,'unbiased'));
                                end
                                xcREshufL(ii,jj,iTrial)  = prctile(xcshuf,5);
                                xcREshufH(ii,jj,iTrial)  = prctile(xcshuf,95);
                            end
                            
                        end
                    end
                end
                
                if shuffleCorrect
                    xcLsig = xcL<xcLshufL | xcL>xcLshufH;
                    xcRsig = xcR<xcRshufL | xcR>xcRshufH;
                    xcLEsig = xcLE<xcLEshufL | xcLE>xcLEshufH;
                    xcREsig = xcRE<xcREshufL | xcRE>xcREshufH;            
                    
                    xcL_  = xcL;
                    xcLE_ = xcLE; 
                    xcR_  = xcR;
                    xcRE_ = xcRE;
                    
                    xcL_(~xcLsig)=nan;
                    xcR_(~xcRsig)=nan;
                    xcLE_(~xcLEsig)=nan;
                    xcRE_(~xcREsig)=nan;
                    
                    xcL_  = nanmean(xcL_,3);
                    xcLE_ = nanmean(xcLE_,3);
                    xcR_  = nanmean(xcR_,3);
                    xcRE_ = nanmean(xcRE_,3);
                else
                    xcL_  = nanmean(xcL,3);
                    xcLE_ = nanmean(xcLE,3);
                    xcR_  = nanmean(xcR,3);
                    xcRE_ = nanmean(xcRE,3);
                end
                %%
                if plotOnline
                    diffL = (xcLE_-xcL_);diffL=diffL(:);diffL=sort(diffL(~isnan(diffL)));
                    diffR = (xcRE_-xcR_);diffR=diffR(:);diffR=sort(diffR(~isnan(diffR)));
                    
                    %                 % Sort by peak mean firing rate time
                    %                 FR_   = zscore([cell2mat(Ltrials_');cell2mat(Rtrials_')]);
                    %                 evt0  = [ones(length(Ltrials_),1);2*ones(length(Rtrials_),1)];
                    %                 [~,meanFR_]  = DecodeStatsTSonly(FR_,evt0,0.05);
                    %                 meanFR_ = (meanFR_{1}+meanFR_{2})./2;
                    %                 meanFR_=meanFR_(:,Ass.NonMembers{iArea});
                    %                 [~,idx] = max(meanFR_,[],1);
                    %                 [~,idx] = sort(idx);
                    
                    subplot(2,3,1)
                    imagesc(xcL_(idx,idx))
                    caxis([0 100])
                    subplot(2,3,2)
                    imagesc(xcLE_(idx,idx))
                    caxis([0 100])
                    subplot(2,3,3); hold on
                    plot(diffL)
                    
                    axis([0 Inf -10 10])
                    subplot(2,3,4)
                    imagesc(xcR_(idx,idx))
                    caxis([0 100])
                    subplot(2,3,5)
                    imagesc(xcRE_(idx,idx))
                    caxis([0 100])
                    subplot(2,3,6); hold on
                    plot(diffR)
                    axis([0 Inf -10 10])
                end
                %%
                
                D.xcCorrMean{iArea}{iDelay}{1}{iFile,1} = xcL_;
                D.xcCorrMean{iArea}{iDelay}{2}{iFile,1} = xcR_;
                D.xcErrMean{iArea}{iDelay}{1}{iFile,1} = xcLE_;
                D.xcErrMean{iArea}{iDelay}{2}{iFile,1} = xcRE_;
                
                D.xcCorrRaw{iArea}{iDelay}{1}{iFile,1} = xcL;
                D.xcCorrRaw{iArea}{iDelay}{2}{iFile,1} = xcR;
                D.xcErrRaw{iArea}{iDelay}{1}{iFile,1} = xcLE;
                D.xcErrRaw{iArea}{iDelay}{2}{iFile,1} = xcRE;
                
                if shuffleCorrect
                    D.xcCorrSig{iArea}{iDelay}{1}{iFile,1} = xcLsig;
                    D.xcCorrSig{iArea}{iDelay}{2}{iFile,1} = xcRsig;
                    D.xcErrSig{iArea}{iDelay}{1}{iFile,1} = xcLEsig;
                    D.xcErrSig{iArea}{iDelay}{2}{iFile,1} = xcREsig;
                end
                
                %             clear Ltrials nL Rtrials nR LeftTrials RightTrials iTrial clear TS TSsig TSshufL TSshufH
            end
        end
        %% Save results
        %         if exist('D')
        %
        %             fnOut = sprintf('%sMixedSelectivity%s%s_%s_DelayDecoding_UnitSpan_%s.mat',pat,filesep,fname,Areas{iArea},UnitSelection);
        %             save(fnOut,'D','tbShort','tbMedium','tbLong','-v7.3');
        %             fprintf('Done.\n')
        %             clear D
        %         end
        %
        % end
    end
end
%% batch process performance
for iFile =1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
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
    noPFC=nu(1);
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
    
    clear units usel_out
    %% Batch process behavioural performance
    L = [length(t.Short.ChoicePress_LeftCorrect)./(length(t.Short.ChoicePress_LeftCorrect)+length(t.Short.ChoicePress_LeftError)) , ...
        length(t.Medium.ChoicePress_LeftCorrect)./(length(t.Medium.ChoicePress_LeftCorrect)+length(t.Medium.ChoicePress_LeftError)) , ...
        length(t.Long.ChoicePress_LeftCorrect)./(length(t.Long.ChoicePress_LeftCorrect)+length(t.Long.ChoicePress_LeftError))];
    
    
    R = [length(t.Short.ChoicePress_RightCorrect)./(length(t.Short.ChoicePress_RightCorrect)+length(t.Short.ChoicePress_RightError)) , ...
        length(t.Medium.ChoicePress_RightCorrect)./(length(t.Medium.ChoicePress_RightCorrect)+length(t.Medium.ChoicePress_RightError)) , ...
        length(t.Long.ChoicePress_RightCorrect)./(length(t.Long.ChoicePress_RightCorrect)+length(t.Long.ChoicePress_RightError))];
    
    D.Accuracy(iFile,:) = (L+R)./2;
    D.AccuracyL(iFile,:) = L;
    D.AccuracyR(iFile,:) = R;
    
    C_ =  [length(t.Short.ChoicePress_LeftCorrect)  + length(t.Short.ChoicePress_RightCorrect),...
        length(t.Medium.ChoicePress_LeftCorrect) + length(t.Medium.ChoicePress_RightCorrect),...
        length(t.Long.ChoicePress_LeftCorrect)   + length(t.Long.ChoicePress_RightCorrect)];
    
    E_ =  [length(t.Short.ChoicePress_LeftError)    + length(t.Short.ChoicePress_RightError),...
        length(t.Medium.ChoicePress_LeftError)   + length(t.Medium.ChoicePress_RightError),...
        length(t.Long.ChoicePress_LeftError)     + length(t.Long.ChoicePress_RightError)];
    
    D.pCorr(iFile,:) = C_./(C_+E_);
    D.AboveChance(iFile,:) = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;

    C_ =  [length(t.Short.ChoicePress_LeftCorrect),...
        length(t.Medium.ChoicePress_LeftCorrect),...
        length(t.Long.ChoicePress_LeftCorrect)];
    
    E_ =  [length(t.Short.ChoicePress_LeftError),...
        length(t.Medium.ChoicePress_LeftError),...
        length(t.Long.ChoicePress_LeftError)];
       
    D.pCorrL(iFile,:) = C_./(C_+E_);
    D.AboveChanceL(iFile,:) = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;


    C_ =  [length(t.Short.ChoicePress_RightCorrect),...
        length(t.Medium.ChoicePress_RightCorrect),...
        length(t.Long.ChoicePress_RightCorrect)];
    
    E_ =  [length(t.Short.ChoicePress_RightError),...
        length(t.Medium.ChoicePress_RightError),...
        length(t.Long.ChoicePress_RightError)];
       
    D.pCorrR(iFile,:) = C_./(C_+E_);
    D.AboveChanceR(iFile,:) = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;
    clear L R C_ E_
    
end

%% Bulk analyse - correct vs error difference
bins = 0:0.1:10;
for iArea =1
    for iDelay =1:length(Delays_)
        for dir_=1:2
            Corr_ = cellfun(@(x) x(:),D.xcCorrMean{iArea}{iDelay}{dir_},'UniformOutput',false);
            Err_  = cellfun(@(x) x(:),D.xcErrMean{iArea}{iDelay}{dir_},'UniformOutput',false);
            Corr_ = cell2mat(cellfun(@(x) histc(x(~isnan(x)),bins),Corr_,'UniformOutput',false)');
            Err_  = cell2mat(cellfun(@(x) histc(x(~isnan(x)),bins),Err_,'UniformOutput',false)');
            Corr{dir_}=cumsum(Corr_)./repmat(sum(Corr_),length(bins),1);
            Err{dir_}=cumsum(Err_)./repmat(sum(Err_),length(bins),1);
        end
        CorrHist{iArea}{iDelay} = (Corr{1}+Corr{2})./2;
        if size(Err{1},2)==size(Err{2},2)
            ErrHist{iArea}{iDelay} = (Err{1}+Err{2})./2;
        else
            if size(Err{1},2)>size(Err{2},2)
                ErrHist{iArea}{iDelay} = Err{1};
            else
                ErrHist{iArea}{iDelay} = Err{2};
            end
        end
    end
    
    for iDelay=1:3
            figure; hold on 

        ciplot(nanmean(CorrHist{iArea}{iDelay},2)+nansem(CorrHist{iArea}{iDelay},2),...
               nanmean(CorrHist{iArea}{iDelay},2)-nansem(CorrHist{iArea}{iDelay},2),...
               bins,color_{iDelay},1)
           
        ciplot(nanmean(ErrHist{iArea}{iDelay},2)+nansem(ErrHist{iArea}{iDelay},2),...
               nanmean(ErrHist{iArea}{iDelay},2)-nansem(ErrHist{iArea}{iDelay},2),...
               bins,color_{iDelay},0.2)           
    end
    plot([0 0],[0 1],':k')
    axis([0 10 0 1])
    
end
%% Bulk analyse - correct vs error difference
bins = -20:0.1:20;
for iArea =1
    for iDelay =1:length(Delays_)
        for dir_=1:2
            diff_ = cellfun(@(x,y) x(:)-y(:), D.xcErrMean{iArea}{iDelay}{dir_},D.xcCorrMean{iArea}{iDelay}{dir_},'UniformOutput',false);
            diff_ = cell2mat(cellfun(@(x) histc(x(~isnan(x)),bins),diff_,'UniformOutput',false)');
            diffHist{dir_}=diff_./size(diff_,2);
        end
        diffHist_{iArea}{iDelay} = (diffHist{1}+diffHist{1})./2;
    end
    
    figure; hold on 
    for iDelay=1:3
        ciplot(nanmean(diffHist_{iArea}{iDelay},2)+nansem(diffHist_{iArea}{iDelay},2),...
               nanmean(diffHist_{iArea}{iDelay},2)-nansem(diffHist_{iArea}{iDelay},2),...
               bins,color_{iDelay},0.8)
    end
    plot([0 0],[0 10],':k')
    axis([-10 10 0 2])
    
end

%% correlation structure vs. performance
for iArea = 1
    figure
    for iDelay =1:length(Delays_)
        for dir_=1:2
            Corr_ = cellfun(@(x) x(:),D.xcCorrMean{iArea}{iDelay}{dir_},'UniformOutput',false);
            CorrOut(:,dir_) = cell2mat(cellfun(@(x) nanmedian(x(~isnan(x))),Corr_,'UniformOutput',false));
%             PerfOut(:,1) = D.AccuracyL(:,iDelay);
%             PerfOut(:,2) = D.AccuracyR(:,iDelay);
            PerfOut(:,1) = D.pCorrL(:,iDelay);
            PerfOut(:,2) = D.pCorrR(:,iDelay);
            AboveChance(:,1) = D.AboveChanceL(:,iDelay);
            AboveChance(:,2) = D.AboveChanceR(:,iDelay);
        end
        subplot(1,3,iDelay);hold on
        x_ = CorrOut(:);
        y_ = PerfOut(:);
        idx = AboveChance(:);
        scatter(x_(~idx),y_(~idx),20,'b','o')
        scatter(x_(idx),y_(idx),20,'b','o','filled')
        
    end
end
%%
figure; hold on 

scatter(D.pCorrL(:,1),D.pCorrL(:,2))
scatter(D.pCorrR(:,1),D.pCorrR(:,2))
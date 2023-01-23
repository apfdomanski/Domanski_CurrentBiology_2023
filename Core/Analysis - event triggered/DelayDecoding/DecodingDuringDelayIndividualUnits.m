% function DecodingDuringDelayIndividualUnits(iFile)
% if isstr(iFile)
%     iFile = str2num(iFile);
% end
%% %%%%%% PREAMBLE %%%%%%
clear 
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
%     pat = '/Volumes/HDD2/DNMTP/raw/';
    pat = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw/';

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
    %% Batch process units
    for iArea = 1:2%length(Areas)
        %% Get the files
        fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
        if iArea < 3
            % load(sprintf('%sKDE_bins%s%s_%s_iFR50.mat',pat,filesep,fname,Areas{iArea}));
            load(sprintf('%sKDE_binsTaskOnly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
            iFR_ = iFR;
            D_.avgFR{iArea}{iFile} = avgFR;
            D_.KDEsigma{iArea}{iFile} = hucv.^0.5;
        else
            U_{1} = load(sprintf('%sKDE_binsTaskOnly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{1}));
            U_{2} = load(sprintf('%sKDE_binsTaskOnly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{2}));
            iFR_ = [U_{1}.iFR,U_{2}.iFR];
            Tmtx = U_{1}.Tmtx;
            clear U_
        end
        
        % iFR_ = zscore(iFR_);
        
    %% Delay period decoding for units
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
                TS_ = cell(nReps,1);meanFR_ = cell(nReps,1);
                TSshuf_ = cell(nReps,1);
                for nrep = 1:nReps
                    if isinf(TrialsSubsample)
                        rsL = 1:length(Ltrials_);
                        rsR = 1:length(Rtrials_);
                        evt0  = [ones(length(Ltrials_),1);2*ones(length(Rtrials_),1)];
                    else
                        rsL = randsample(length(Ltrials_),TrialsSubsample);
                        rsR = randsample(length(Rtrials_),TrialsSubsample);
                        evt0  = [ones(TrialsSubsample,1);2*ones(TrialsSubsample,1)];
                    end
                    
                    FR_   = [cell2mat(Ltrials_(rsL)');cell2mat(Rtrials_(rsR)')];
                    
                    [TS_{nrep},meanFR_{nrep}]  = DecodeStatsTSonly(FR_,evt0,0.05);
                    TSshuf_{nrep}  = DecodeStatsTSonly(FR_,evt0(randperm(length(evt0))),0.05);
                end
                TS__=[]; TSshuf__=[]; meanFR__{1}=[];meanFR__{2}=[];
                for nrep = 1:nReps
                    meanFR__{1} = cat(3,meanFR__{1},meanFR_{nrep}{1});
                    meanFR__{2} = cat(3,meanFR__{2},meanFR_{nrep}{2});
                    TS__ = cat(3,TS__,TS_{nrep});
                    TSshuf__ = cat(3,TSshuf__,TSshuf_{nrep});
                end
                meanFR{1}     = nanmean(meanFR__{1},3);
                meanFR{2}     = nanmean(meanFR__{2},3);
                TS     = nanmean(TS__,3);
                TSshufL = prctile(TSshuf__,5,3);
                TSshufH = prctile(TSshuf__,95,3);
                if isinf(TrialsSubsample)
                    TSsig = tpdf(TS,length(evt0)-2)<0.05;
                else
                    TSsig = TS>TSshufH;
                end
                clear nrep TS__ TSshuf__ TS_ TSshuf_ FR_ rsL rsR evt0
                
                D.meanFR{iArea}{iDelay} = meanFR;
                D.TS{iArea}{iDelay}     = TS;
                D.TSsig{iArea}{iDelay}  = TSsig;
                D.TSshufL{iArea}{iDelay}  = TSshufL;
                D.TSshufH{iArea}{iDelay}  = TSshufH;
                
                D_.meanFR{iArea}{iDelay}{1}{iFile,1}  = meanFR{1};
                D_.meanFR{iArea}{iDelay}{2}{iFile,1}  = meanFR{2};
                D_.TS{iArea}{iDelay}{iFile,1}     = TS;
                D_.TSsig{iArea}{iDelay}{iFile,1}  = TSsig;
                D_.TSshufL{iArea}{iDelay}{iFile,1}  = TSshufL;
                D_.TSshufH{iArea}{iDelay}{iFile,1}  = TSshufH;
                
            else
                D.meanFR{iArea}{iDelay}{1} = nan(length(tb_),nu(iArea));
                D.meanFR{iArea}{iDelay}{2} = nan(length(tb_),nu(iArea));
                D.TS{iArea}{iDelay}       = nan(length(tb_),nu(iArea));
                D.TSsig{iArea}{iDelay}    = false(length(tb_),nu(iArea));
                D.TSshufL{iArea}{iDelay}  = nan(length(tb_),nu(iArea));
                D.TSshufH{iArea}{iDelay}  = nan(length(tb_),nu(iArea));
                
                D_.meanFR{iArea}{iDelay}{1}{iFile,1} = nan(length(tb_),nu(iArea));
                D_.meanFR{iArea}{iDelay}{2}{iFile,1} = nan(length(tb_),nu(iArea));
                D_.TS{iArea}{iDelay}{iFile,1}       = nan(length(tb_),nu(iArea));
                D_.TSsig{iArea}{iDelay}{iFile,1}    = false(length(tb_),nu(iArea));
                D_.TSshufL{iArea}{iDelay}{iFile,1}  = nan(length(tb_),nu(iArea));
                D_.TSshufH{iArea}{iDelay}{iFile,1}  = nan(length(tb_),nu(iArea));
            end
            %%%%
            
            %             FR = [Ltrials;Rtrials];
            %             evt0 = [ones(nL,1);2*ones(nR,1)];
            %             Ltr = length(tb_);
            %             reg = 0.05;
            %             TS  = DecodeStatsTSonly(FR,evt0,reg);
            %             TSsig = tpdf(TS,length(evt0)-2)<0.05;
            %             D.TS{iArea}{iDelay}     = TS;
            %             D.TSsig{iArea}{iDelay}  = TSsig;
            %
            %             D_.TS{iArea}{iDelay}{iFile,1}     = TS;
            %             D_.TSsig{iArea}{iDelay}{iFile,1}  = TSsig;
            clear Ltrials nL Rtrials nR LeftTrials RightTrials iTrial clear TS TSsig TSshufL TSshufH
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
%% batch process performance
for iFile =1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    

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

%% Get the assembly membership info
 for iFile =1:length(fileList)
      for iArea = 1:2
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
    Ass_{iFile} = Ass;
    clear units usel_out
      end
 end
%% Collapse
iArea = 1;
temp =[];
% remove any dead entries
for iDelay=1:length(Delays_)
    temp=[temp,cell2mat(cellfun(@(x) nansum(nansum(isnan(x)))>0, D_.TS{iArea}{iDelay},'UniformOutput', false ))];
end
temp = sum(temp,2)>0;

for iArea = 1:2%length(Areas)
    
    for iDelay=1:length(Delays_)
        
        D_.TS_Collapsed{iArea}{iDelay}    = cell2mat(D_.TS{iArea}{iDelay}(~temp)');
        D_.TSsig_Collapsed{iArea}{iDelay} = cell2mat(D_.TSsig{iArea}{iDelay}(~temp)');
        D_.meanFR_Collapsed{iArea}{iDelay}{1} = cell2mat(D_.meanFR{iArea}{iDelay}{1}(~temp)');
        D_.meanFR_Collapsed{iArea}{iDelay}{2} = cell2mat(D_.meanFR{iArea}{iDelay}{2}(~temp)');
        
        
        
    end
    D_.UnitsAboveChance{iArea} = [];
    nUnits = cellfun(@(x) size(x,2), D_.TS{iArea}{1});

    for iFile = 1:length(fileList)
        if temp(iFile)~=1
            D_.UnitsAboveChance{iArea} = [D_.UnitsAboveChance{iArea}; repmat(double(D.AboveChance(iFile,:)),nUnits(iFile),1)];
        end
    end

    
    D_.avgFR_Collapsed{iArea}= cell2mat(D_.avgFR{iArea}(~temp));
    D_.KDEsigma_Collapsed{iArea}= cell2mat(D_.KDEsigma{iArea}(~temp));
    
    
end

% Collapse assembly annotation
    
for iArea =1:2
    for iFile = 1:length(fileList)
        
        AssTemp = zeros(1,size(D_.meanFR{iArea}{iDelay}{1}{iFile},2));
        AssTempGlobal = -ones(1,size(D_.meanFR{iArea}{iDelay}{1}{iFile},2));
        AssTempGlobal(Ass_{iFile}.NonMembers{iArea})=0;
        AssTempGlobal(Ass_{iFile}.LocalMembers{iArea})=1;
        AssTempGlobal(Ass_{iFile}.JointMembers{iArea})=2;
        
        AssTemp(unique([Ass_{iFile}.LocalMembers{iArea},Ass_{iFile}.JointMembers{iArea}]))=1;
        
        AssCollapsed.class{iArea}{iFile} = AssTemp;
        AssCollapsed.classGlobal{iArea}{iFile} = AssTempGlobal;
    end
     AssCollapsed.classCollapsed{iArea} = cell2mat( AssCollapsed.class{iArea}(~temp))';
     AssCollapsed.classGlobalCollapsed{iArea} = cell2mat( AssCollapsed.classGlobal{iArea}(~temp))';
end
clear temp
%% plot decoding staggered - all files collapsed 

sortCriterion = 'weighted';  % 'first','peak','weighted'
shiftNS = false;
maskNS = true;
for iArea = 1%:2%length(Areas)
    clear idx
    figure
    iDelay = 1;
    
    TS_= D_.TS_Collapsed{iArea}{iDelay}';
    TSsig_= D_.TSsig_Collapsed{iArea}{iDelay}';
    switch sortCriterion
        case 'peak'
            
            if shiftNS % Sort by peak time, separate insignificant units                
                idxNS = find(nansum(TSsig_,2)==0);
                [~,idx] = max(TS_,[],2);
            else % Sort by peak time only
                idxNS= zeros(size(TSsig_,1),1);
                [~,idx] = max(TS_,[],2);
            end
            
        case 'first' %%%% Time of first signifncant encoding
            idx= zeros(size(TSsig_,1),1);
            for i=1:size(TSsig_,1)
                pos_ = find(double(TSsig_(i,:)),1,'first');
                if isempty(pos_)
                    pos_ = size(TSsig_,2);
                end
                idx(i) = pos_;
            end
            clear i pos_
            
        case 'weighted'
                            idxNS= zeros(size(TSsig_,1),1);

            idx = zeros(size(TSsig_,1),1);
            for i=1:size(TSsig_,1)
                [~,pos_,w,~] = findpeaks([0,double(TSsig_(i,:)),0]');
                
                
                if isempty(pos_)
                    pos_ = size(TSsig_,2);
                else
                    pos_ = pos_(find(w==max(w),1,'first')) + round(w(find(w==max(w),1,'first'))/2) - 1;
                end
                idx(i) = pos_;
            end
            clear i pos_ w
    end
    
    [~,idx] = sort(idx);
    idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];
    idxOut = idx;
    
    for iDelay=1:length(Delays_)
        TS_= D_.TS_Collapsed{iArea}{iDelay}';
        TSsig_= D_.TSsig_Collapsed{iArea}{iDelay}';
        
        
        TS_     = TS_(idx,:);
        TSsig_  = TSsig_ (idx,:);
        if maskNS
            TS_(~TSsig_ )=0;
        end
        subplot(3,1,iDelay); hold on
        x_  = repmat((1:size(TS_,2)).*bw,size(TS_,1),1);
        y_  = repmat(1:size(TS_,1),size(TS_,2),1);
        imagesc(x_(:),y_(:),TS_);
        set(gca,'YDir','normal')
            cmap =([1 1 1;(hot)]);
%             cmap =([1 1 1;(jet)]);
%         cmap =(jet);
        colormap (cmap)
        caxis([0 5])
        
        plot([4 4],[1 size(TS_,1)],'color',0.6*[1 1 1],'LineWidth',1.2)
        plot([8 8],[1 size(TS_,1)],'color',0.6*[1 1 1],'LineWidth',1.2)
        plot([16 16],[1 size(TS_,1)],'color',0.6*[1 1 1],'LineWidth',1.2)
        axis([0 20 0 350])
        if iArea==1 && iDelay ==3
            plot([18 20],[20 20],'k','LineWidth',1.5)
            plot([18 18],[20 70],'k','LineWidth',1.5)
        end
        % axis off
    end
end
%% plot decoding staggered - each recording
delaytime = [4,8,16];
sortCriterion = 'weighted';  % 'first','peak','weighted'
shiftNS = false;
maskNS = true;
for iFile = 5%:length(fileList)
for iArea = 1%:length(Areas)
    clear idx
    figure
    iDelay = 1;
    
    TS_= D_.TS{iArea}{iDelay}{iFile}';
    TSsig_= D_.TSsig{iArea}{iDelay}{iFile}';
    switch sortCriterion
        case 'peak'
            
            if shiftNS % Sort by peak time, separate insignificant units                
                idxNS = find(nansum(TSsig_,2)==0);
                [~,idx] = max(TS_,[],2);
            else % Sort by peak time only
                idxNS= zeros(size(TSsig_,1),1);
                [~,idx] = max(TS_,[],2);
            end
            
        case 'first' %%%% Time of first signifncant encoding
            idx= zeros(size(TSsig_,1),1);
            for i=1:size(TSsig_,1)
                pos_ = find(double(TSsig_(i,:)),1,'first');
                if isempty(pos_)
                    pos_ = size(TSsig_,2);
                end
                idx(i) = pos_;
            end
            clear i pos_
            
        case 'weighted'
            idx = zeros(size(TSsig_,1),1);
            for i=1:size(TSsig_,1)
                [~,pos_,w,~] = findpeaks([0,double(TSsig_(i,:)),0]');
                
                
                if isempty(pos_)
                    pos_ = size(TSsig_,2);
                else
                    pos_ = pos_(find(w==max(w),1,'first')) + round(w(find(w==max(w),1,'first'))/2) - 1;
                end
                idx(i) = pos_;
            end
            clear i pos_ w
    end
    
    [~,idx] = sort(idx);
    idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];
    
    
    for iDelay=1:length(Delays_)
        TS_= D_.TS{iArea}{iDelay}{iFile}';
        TSsig_= D_.TSsig{iArea}{iDelay}{iFile}';
        
        
        TS_     = TS_(idx,:);
        TSsig_  = TSsig_ (idx,:);
        if maskNS
            TS_(~TSsig_ )=0;
        end
        subplot(3,1,iDelay); hold on
        x_  = repmat((1:size(TS_,2)).*bw,size(TS_,1),1);
        y_  = repmat(1:size(TS_,1),size(TS_,2),1);
        imagesc(x_(:),y_(:),TS_);
                plot(delaytime(iDelay)*[1 1],[1 size(TS_,1)],'color',0.6*[1 1 1],'LineWidth',1.5,'LineStyle',':')

        set(gca,'YDir','normal')
        %     cmap =([1 1 1;(hot)]);
        if maskNS
            cmap =([1 1 1;(jet)]);
        else
            cmap =(jet);
        end
        colormap (cmap)
        caxis([0 5])
       
        axis([0 20 0 size(TS_,1)])
        if iArea==1 && iDelay ==1
            plot([6 8],[10 10],'k','LineWidth',1.5)
            plot([6 6],[10 25],'k','LineWidth',1.5)
        end
        set(gca,'XTick',[0 delaytime(iDelay)])
        % axis off
    end
end
end
%% plot mean FR staggered - all files collapsed 

delaytime = [4,8,16];
dir_ = 1;
for iArea = 1%length(Areas)
    clear idx
    figure
    iDelay = 3;
    
    meanFR_= D_.meanFR_Collapsed{iArea}{iDelay}{dir_};
%     meanFR_ = zscore(meanFR_);
%     xc=pdist(meanFR_,'correlation');
%     xc = -(squareform(xc)-1);
%     xc=xcorr(meanFR_,10,'unbiased');
%     squareform(xc(101,:))
%     idx = 1:size(meanFR_,2)
    [~,idx] = max(meanFR_,[],1);
    [~,idx] = sort(idx);
    
    AssCollapsed_ = AssCollapsed.classGlobalCollapsed{iArea}(idx');
%     AssCollapsed_ = AssCollapsed.classCollapsed{iArea}(idx');
    avgFRCollapsed_ = D_.avgFR_Collapsed{iArea}(idx');
    for iDelay=1:length(Delays_)
        meanFR_ = D_.meanFR_Collapsed{iArea}{iDelay}{dir_};
        meanFR_ = zscore(meanFR_);
%         [~,idx] = max(meanFR_,[],1);
%         [~,idx] = sort(idx);
    
        subplot(3,1,iDelay); hold on
        x_  = repmat((1:size(meanFR_,1)).*bw,size(meanFR_,1),1);
        y_  = repmat(1:size(meanFR_,2),size(meanFR_,2),1);
        imagesc(x_(:),y_(:),meanFR_(:,idx)');
        set(gca,'YDir','normal')
        
%         Plot nonmembers
        idxClass = find(AssCollapsed_ == 0 );
%         scatter(zeros(size(idxClass)),idxClass,5,col_{1},'filled')
        plot(repmat([-0.8 -0.2],size(idxClass,1),1)',[idxClass,idxClass]','color',[col_{1},0.6],'LineWidth',1.2)        
%         Plot joint members
        idxClass = find(AssCollapsed_ == 2);
%         scatter(zeros(size(idxClass)),idxClass,5,col_{3},'filled')
        plot(repmat([-0.8 -0.2],size(idxClass,1),1)',[idxClass,idxClass]','color',[col_{3},0.6],'LineWidth',1.2)
%         Plot local members
        idxClass = find(AssCollapsed_ == 1);
%         scatter(zeros(size(idxClass)),idxClass,5,col_{2},'filled')
        plot(repmat([-0.8 -0.2],size(idxClass,1),1)',[idxClass,idxClass]','color',[col_{2},0.6],'LineWidth',1.2)


%         % Plot nonmembers
%         idxClass = find(AssCollapsed_ == 0);
%         plot(repmat([-1.1 0],size(idxClass,1),1)',[idxClass,idxClass]','color',0.6*[1 1 1],'LineWidth',1.2)
% 
%         % Plot members
%         idxClass = find(AssCollapsed_ ==1 );
%         plot(repmat([-1.1 0],size(idxClass,1),1)',[idxClass,idxClass]','color','b','LineWidth',1.5)
        
        cmap =(jet);
        colormap (cmap)
        caxis([0 2])
        plot(delaytime(iDelay)*[1 1],[1 size(TS_,1)],'color',0.6*[1 1 1],'LineWidth',1.5,'LineStyle',':')
        set(gca,'XTick',[0 delaytime(iDelay)])

                        
        axis([-1 20 0 300])
        if iArea==1 && iDelay ==1
%             plot([18 20],[20 20],'k','LineWidth',1.5)
%             plot([18 18],[20 70],'k','LineWidth',1.5)
             plot([6 8],[20 20],'k','LineWidth',1.5)
            plot([6 6],[20 120],'k','LineWidth',1.5)
        end
        % axis off
    end
end
%% plot peak info time vs. mean FR/intrinsic properties


delaytime = [4,8,16];
dir_ = 1;
for iArea = 1%length(Areas)
    clear idx
    figure
    for iDelay = 1:3
        subplot(3,1,iDelay); hold on
        
        meanFR_= D_.meanFR_Collapsed{iArea}{iDelay}{dir_};
        %     TS_= D_.TS_Collapsed{iArea}{iDelay};
        %     meanFR_ = zscore(meanFR_);
        %     xc=pdist(meanFR_,'correlation');
        %     xc = -(squareform(xc)-1);
        %     xc=xcorr(meanFR_,10,'unbiased');
        %     squareform(xc(101,:))
        %     idx = 1:size(meanFR_,2)
        [peak_,idx] = max(meanFR_,[],1);
        [~,idx] = sort(idx);
        avgFRCollapsed_ = D_.avgFR_Collapsed{iArea}(idx');
        %         scatter(avgFRCollapsed_,peak_)
        
        KDEsigmaCollapsed_ = D_.KDEsigma_Collapsed{iArea}(idx');
        scatter(avgFRCollapsed_,idx)
        %         scatter(avgFRCollapsed_,KDEsigmaCollapsed_)
%         set(gca,'XScale','log','YScale','linear')
    end
end
%% plot peak info time vs. mean FR/intrinsic properties (1)

delaytime = [4,8,16];
dir_ = 1;
for iArea = 1%length(Areas)
    clear idx
    figure
    for iDelay = 1:3
        subplot(3,1,iDelay); hold on
        
%             includedUnits = find(D_.UnitsAboveChance{iArea}(:,iDelay));
            includedUnits = 1:length(D_.UnitsAboveChance{iArea}(:,iDelay));
            TS_= D_.TS_Collapsed{iArea}{iDelay}(:,includedUnits)';
            TSsig_= D_.TSsig_Collapsed{iArea}{iDelay}(:,includedUnits)';
            KDE_= D_.KDEsigma_Collapsed{iArea}(includedUnits)';
            meanFR_= D_.avgFR_Collapsed{iArea}(includedUnits)';
            clear t_max t_min
            
            
            duration      = sum(TSsig_,2).*bw;
            
            [~,t_max_peak]  = max(TS_,[],2);
            t_max_peak=t_max_peak*bw ;
            t_min = nan(size(TSsig_,1),1);
            t_max = nan(size(TSsig_,1),1);
            t_average = nan(size(TSsig_,1),1);
            for iUnit = 1:size(TSsig_,1)
                a = TSsig_(iUnit ,:);
                first_idx = find(a,1,'first'); if isempty(first_idx), first_idx=NaN; end
                last_idx = find(a,1,'last'); if isempty(last_idx), last_idx=NaN; end
                t_av = find(a); if isempty(t_av), t_av=NaN; else t_av = mean(t_av); end
                t_min(iUnit)=first_idx*bw;
                t_max(iUnit)=last_idx*bw;
                t_average(iUnit)=t_av*bw;
             
%                 % OR: Time of longest contiguous block
%                 [~,locs,w] = findpeaks(double([0,a,0]));
%                 if ~isempty(locs)
%                     if sum(w==max(w))>1
%                         locs = max(locs);
%                         w = max(w);
%                     end
%                     t_min(iUnit) = (locs(w==max(w))-1)*bw;
%                     t_max(iUnit) = (locs(w==max(w))-1+max(w))*bw;
%                 end
                
            end
            t_mid = t_min+(t_max-t_min)./2;
            
            x = t_max; %1./KDE_ 
            y = duration ;
            idx = isnan(x) | y==0;
            x(idx)=[];y(idx)=[];
            scatter(x,y,'k')
            
            
            mdl   = fitlm(x,y,'RobustOpts','on' );
            pFit  = mdl.Coefficients{2,4};
            Rsq   = mdl.Rsquared.Adjusted;
            slope = mdl.Coefficients{2,1};
            intercept = mdl.Coefficients{1,1};
            
            plot(mdl.VariableInfo.Range{1},(mdl.VariableInfo.Range{1}) * slope + intercept,'k','LineWidth',1.5)
            [pFit, Rsq]
            
            axis([0 20 0 10])
%             set(gca,'Xscale','log')
            switch iDelay
                case  1
                    title('4s delay')
                case  2
                    title('8s delay')                    
                    ylabel('Span of significant decoding (s)')
                case  3
                    title('16s delay')                                        
                    xlabel('last significant timepoint (s)')
            end
    end
            
    
end
% @TODO: sort by above/at chance
% @TODO: decoding time Vs. time constant
%% plot peak info time vs. mean FR/intrinsic properties (2)

delaytime = [4,8,16];
dir_ = 1;
for iArea = 1%length(Areas)
    clear idx
    figure
    for iDelay = 1:3
        subplot(3,1,iDelay); hold on
        
%             includedUnits = find(D_.UnitsAboveChance{iArea}(:,iDelay));
            includedUnits = 1:length(D_.UnitsAboveChance{iArea}(:,iDelay));
            TS_= D_.TS_Collapsed{iArea}{iDelay}(:,includedUnits)';
            TSsig_= D_.TSsig_Collapsed{iArea}{iDelay}(:,includedUnits)';
            KDE_= D_.KDEsigma_Collapsed{iArea}(includedUnits)';
            meanFR_= D_.avgFR_Collapsed{iArea}(includedUnits)';
            clear t_max t_min
            
            
            duration      = sum(TSsig_,2).*bw;
            
            [~,t_max_peak]  = max(TS_,[],2);
            t_max_peak=t_max_peak*bw ;
            t_min = nan(size(TSsig_,1),1);
            t_max = nan(size(TSsig_,1),1);
            t_average = nan(size(TSsig_,1),1);
            for iUnit = 1:size(TSsig_,1)
                a = TSsig_(iUnit ,:);
                first_idx = find(a,1,'first'); if isempty(first_idx), first_idx=NaN; end
                last_idx = find(a,1,'last'); if isempty(last_idx), last_idx=NaN; end
                t_av = find(a); if isempty(t_av), t_av=NaN; else t_av = mean(t_av); end
                t_min(iUnit)=first_idx*bw;
                t_max(iUnit)=last_idx*bw;
                t_average(iUnit)=t_av*bw;
             
%                 % OR: Time of longest contiguous block
%                 [~,locs,w] = findpeaks(double([0,a,0]));
%                 if ~isempty(locs)
%                     if sum(w==max(w))>1
%                         locs = max(locs);
%                         w = max(w);
%                     end
%                     t_min(iUnit) = (locs(w==max(w))-1)*bw;
%                     t_max(iUnit) = (locs(w==max(w))-1+max(w))*bw;
%                 end
                
            end
            t_mid = t_min+(t_max-t_min)./2;
            
            x = 1./sqrt(KDE_); 
            y = duration ;
            idx = isnan(x) | y==0;
            x(idx)=[];y(idx)=[];
            scatter(x,y,'k')
            
            
            mdl   = fitlm(x,y,'RobustOpts','on' );
            pFit  = mdl.Coefficients{2,4};
            Rsq   = mdl.Rsquared.Adjusted;
            slope = mdl.Coefficients{2,1};
            intercept = mdl.Coefficients{1,1};
            
            plot(mdl.VariableInfo.Range{1},(mdl.VariableInfo.Range{1}) * slope + intercept,'k','LineWidth',1.5)
            [pFit, Rsq]
            
            axis([0 3 0 10])
%             set(gca,'Xscale','log')
            switch iDelay
                case  1
                    title('4s delay')
                case  2
                    title('8s delay')                    
                    ylabel('Span of significant decoding (s)')
                case  3
                    title('16s delay')                                        
                    xlabel('KDE time-constant (1/\surd(\sigma) )')
            end
    end
end
% @TODO: sort by above/at chance
% @TODO: decoding time Vs. time constant
%%
figure; 
scatter(meanFR_,sqrt(KDE_))
%% plot mean FR  sequence - all files collapsed 
dir_ = 1;
for iArea = 1%length(Areas)
    clear idx
   
    iDelay = 1;
    maxBin = 1;
    bins =(-maxBin:bw:maxBin);
    meanFR_ = D_.meanFR_Collapsed{iArea}{iDelay}{dir_};
    meanFR_ = zscore(meanFR_);
    [~,idx] = max(meanFR_,[],1);
    [~,idx] = sort(idx);
    
    for iDelay=1:length(Delays_)
        meanFR_= D_.meanFR_Collapsed{iArea}{iDelay}{dir_};
%         meanFR_ = zscore(meanFR_);
        
%         [~,idx] = max(meanFR_,[],1);
%         [~,idx] = sort(idx);
        meanFR_ = meanFR_(:,idx);
%         meanFR_ = meanFR_(:,idxOut); import previous sort 

        xc = zeros(size(meanFR_,2));        xcTime  = zeros(size(meanFR_,2));

        for ii=1:size(meanFR_,2)
            ii
            parfor jj=ii+1:size(meanFR_,2)
%                 [ii jj]
                temp = [meanFR_(:,ii);meanFR_(:,jj)];
                temp = zscore(temp);
                unit1 = temp(1:size(meanFR_,1));
                unit2 = temp(size(meanFR_,1)+1:2*size(meanFR_,1));
                temp = xcorr(unit2,unit1,maxBin/bw,'unbiased');
                temp(maxBin/bw+1)=NaN;
                [xc(ii,jj), xcTime(ii,jj)]  = nanmax(temp);
                xcTime(ii,jj) = bins(xcTime(ii,jj));
            end
        end
        xc_{iDelay}     = xc;
        xcTime_{iDelay} = xcTime;
    end
    %%
     figure
        for iDelay=1:length(Delays_)

        xc__ = smooth2a(xc_{iDelay},0,0);
        xc__(xc__<0)=0;
        xc__ = smooth2a(xc__,2,2);
        subplot(3,1,iDelay); hold on
%         x_  = repmat((1:size(meanFR_,1)).*bw,size(meanFR_,1),1);
%         y_  = repmat(1:size(meanFR_,2),size(meanFR_,2),1);
%         imagesc(x_(:),y_(:),meanFR_(:,idx)');
        imagesc(xc__)
        set(gca,'YDir','normal')
            
%         cmap =(jet);
        cmap = load('blue_white_red.mat');cmap =cmap.cmap;
        cmap(1:size(cmap,1)/2,:)=[];
%         cmap =  ([summer;[1 1 1];flipud(autumn)])

        colormap (cmap)
        caxis([0 0.2])
        axis([0 350 0 350])
%         caxis([-Inf Inf])
        axis off
        axis square
        end
       %%
     figure
        for iDelay=1:length(Delays_)

        xcTime__ = smooth2a(xcTime_{iDelay},0,0);
        xcTime__(xc_{iDelay}<0)=0;
        xcTime__ = smooth2a(xcTime__,2,2);
        subplot(3,1,iDelay); hold on
%         x_  = repmat((1:size(meanFR_,1)).*bw,size(meanFR_,1),1);
%         y_  = repmat(1:size(meanFR_,2),size(meanFR_,2),1);
%         imagesc(x_(:),y_(:),meanFR_(:,idx)');
        imagesc(xcTime__)
        set(gca,'YDir','normal')
            
        cmap =(jet);
%         cmap = load('blue_white_red.mat');cmap =cmap.cmap;
        colormap ([summer;[1 1 1];flipud(autumn)])
        colormap ([[1 1 1];flipud(autumn)])
        caxis([0 0.5])
        
%         plot([4 4],[1 size(meanFR_,1)],'color',0.6*[1 1 1],'LineWidth',1.2)
%         plot([8 8],[1 size(meanFR_,1)],'color',0.6*[1 1 1],'LineWidth',1.2)
%         plot([16 16],[1 size(meanFR_,1)],'color',0.6*[1 1 1],'LineWidth',1.2)
        axis([0 350 0 350])
%         caxis([-Inf Inf])
%         if iArea==1 && iDelay ==3
%             plot([18 20],[20 20],'k','LineWidth',1.5)
%             plot([18 18],[20 70],'k','LineWidth',1.5)
%         end
        axis off
        axis square

        
        
    end
end
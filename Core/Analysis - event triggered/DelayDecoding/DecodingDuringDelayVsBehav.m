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
    pat = '/Volumes/B001/Bristol/AssemblyAnalysis/raw/';

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
UnitSelection = 'all'; %{'pairs','groups','all'}
groupSize = 5;
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
mkdir([pat 'MixedSelectivity'])

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
%% Batch process 
for iFile = 1:length(fileList)
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
    %% Batch process behaviour
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
    for iArea = 2%1:length(Areas)
        %% Get the files
        fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
        if iArea < 3
            % load(sprintf('%sKDE_bins%s%s_%s_iFR50.mat',pat,filesep,fname,Areas{iArea}));
            load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
            iFR_ = iFR;
        else
            U_{1} = load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{1}));
            U_{2} = load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{2}));
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
            
            Ltrials = [];nL = 0;
            for iTrial =1:size(LeftTrials,1)
                try
                    tlims_  = LeftTrials(iTrial,1)/1e6+tlims_X(1);
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    Ltrials = [Ltrials;iFR_(tlims_,:)];
                    nL=nL+1;
                end
            end
            Rtrials = [];nR = 0;
            for iTrial =1:size(RightTrials,1)
                try
                    tlims_  = RightTrials(iTrial,1)/1e6+tlims_X(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    Rtrials = [Rtrials;iFR_(tlims_,:)];
                    nR=nR+1;
                end
            end
            LtrialsE = [];nLe = 0;
            for iTrial =1:size(LeftTrialsE,1)
                try
                    tlims_  = LeftTrialsE(iTrial,1)/1e6+tlims_X(1) + shift;
                    tlims_ = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    LtrialsE = [LtrialsE;iFR_(tlims_,:)];
                    nLe=nLe+1;
                end
            end
            RtrialsE = [];nRe = 0;
            for iTrial =1:size(RightTrialsE,1)
                try
                    tlims_  = RightTrialsE(iTrial,1)/1e6+tlims_X(1) + shift;
                    tlims_ = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    RtrialsE = [RtrialsE;iFR_(tlims_,:)];
                    nRe=nRe+1;
                end
            end
            
            %%%% Cross-validation error decoder (correct trials)
            if nL>=2 && nR>=2
                fprintf('File [%d/%d], %s Units, %s delay: Decoding correct trials...\n',iFile,length(fileList),Areas{iArea},Delays_{iDelay})
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = length(tb_);
                switch UnitSelection
                    case {'groups','all'} % decode from drawn groupings of units
                        %if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
                        if strcmp(UnitSelection,'all')
                            groupSize_ = nu(iArea);
                            nrep=1;
                        else
                            nrep=10; 
                            groupSize_ = min([groupSize,nu(iArea)]);
                        end
                       
                        CVE = cell(nrep,1);
                        SigDecodeTime= zeros(nrep,1);
                        SigDecodeFracDraws = cell(nrep,1);
                        if strcmp(UnitSelection,'groups')
                            parfor i=1:nrep
                                if iArea<3
                                    rs{i} = randsample(size(Ltrials,2),groupSize_);
                                else
                                    rs{i} = [randsample(nu(1),floor(groupSize_/2));...
                                        randsample(nu(2),ceil(groupSize_/2))+nu(1)];
                                end
%                                 CVE{i,1} = DecodeCVE(FR(:,rs{i}),evt0,0.05);
                                
                                L_ = cell(1,nL);
                                R_ = cell(1,nR);
                                for iL = 1:nL
                                    idx=(1:Ltr)+(iL-1)*Ltr;
                                    L_{iL} = Ltrials(idx,rs{i});
                                end
                                for iR = 1:nR
                                    idx=(1:Ltr)+(iR-1)*Ltr;
                                    R_{iR} = Rtrials(idx,rs{i});
                                end
                                
                                [CVE_,~,~,~,~,CVE_shuffledHL] = RunDecodeGroups(1:groupSize_,50,L_,R_,inf,'max');
                                CVE{i,1}=CVE_';
                                
                                SigDecodeTime(i)                  = sum(CVE_>CVE_shuffledHL(:,2));
                                SigDecodeFracDraws{i,1}(1:Ltr)    = CVE_>CVE_shuffledHL(:,2);
                                CVE_shuffledCIL{i,1}(1:Ltr)       = CVE_shuffledHL(:,1);
                                CVE_shuffledCIH{i,1}(1:Ltr)       = CVE_shuffledHL(:,2);
                            end
                        elseif strcmp(UnitSelection,'all')
%                             CVE{1} = DecodeCVE(FR,evt0,0.05);
                            L_ = cell(1,nL);
                            R_ = cell(1,nR);
                            for iL = 1:nL
                                idx=(1:Ltr)+(iL-1)*Ltr;
                                L_{iL} = Ltrials(idx,1:nu(iArea));
                            end
                            for iR = 1:nR
                                idx=(1:Ltr)+(iR-1)*Ltr;
                                R_{iR} = Rtrials(idx,1:nu(iArea));
                            end
                            
                            [CVE_,~,~,~,~,CVE_shuffledHL] = RunDecodeGroups(1:nu(iArea),1,L_,R_,inf,'max');
                            CVE{1}=CVE_';
                            SigDecodeTime(1)                = sum(CVE_>CVE_shuffledHL(:,2));
                            SigDecodeFracDraws{1}(1:Ltr)    = CVE_>CVE_shuffledHL(:,2);
                            CVE_shuffledCIL{1}(1:Ltr)       = CVE_shuffledHL(:,1);
                            CVE_shuffledCIH{1}(1:Ltr)       = CVE_shuffledHL(:,2);
                        end
                        
                    case 'pairs' % Average decoding across all pairs of units
                        % Prepare data arrays
                        nrep=1;
                        CVE = cell(nu(iArea),nu(iArea)); 
                        SigDecodeTime = nan(nu(iArea),nu(iArea));
                        SigDecodeFracDraws = cell(nu(iArea),nu(iArea));
                        CVE_shuffledCIH = cell(nu(iArea),nu(iArea));
                        CVE_shuffledCIL = cell(nu(iArea),nu(iArea));
                        fprintf('Decoding pair: ')
                        L_ = cell(1,nL);
                        R_ = cell(1,nR);
                        for iL = 1:nL
                            idx=(1:Ltr)+(iL-1)*Ltr;
                            L_{iL} = Ltrials(idx,1:nu(iArea));
                        end
                        for iR = 1:nR
                            idx=(1:Ltr)+(iR-1)*Ltr;
                            R_{iR} = Rtrials(idx,1:nu(iArea));
                        end
                        % Run decoding matrices
                        if iArea < 3
                            
                            %CVE{i,j} = DecodeCVE(FR(:,[i,j]),evt0,0.05);
                            
                            for i =1:nu(iArea)
                                fprintf('%d,',i)
                                parfor j =i+1:nu(iArea)
                                    
                                    [CVE_,~,~,~,~,CVE_shuffledHL] = RunDecodeGroups([i,j],1,L_,R_,inf,'max');
                                    CVE{i,j}=CVE_';
                                    
                                    % Number of significant decoding timesteps for this draw
                                    SigDecodeTime(i,j)                = sum(CVE_>CVE_shuffledHL(:,2)); 
                                    SigDecodeFracDraws{i,j}(1:Ltr)    = CVE_>CVE_shuffledHL(:,2);
                                    CVE_shuffledCIL{i,j}(1:Ltr)       = CVE_shuffledHL(:,1)
                                    CVE_shuffledCIH{i,j}(1:Ltr)       = CVE_shuffledHL(:,2)
                                    
                                end
                            end
                            
                        else
                            for i =1:nu(1)
                                fprintf('%d,',i)
                                parfor j = nu(1)+1:nu(3)
                                    %CVE{i,j} = DecodeCVE(FR(:,[i,j]),evt0,0.05);
                                    [CVE_,~,~,~,~,CVE_shuffledHL] = RunDecodeGroups([i,j],1,L_,R_,inf,'max');
                                    CVE{i,j}=CVE_';
                                    
                                    SigDecodeTime(i,j)                = sum(CVE_>CVE_shuffledHL(:,2)); % Number of significant decoding timesteps for this draw
                                    SigDecodeFracDraws{i,j}(1:Ltr)    = CVE_>CVE_shuffledHL(:,2);
                                end
                            end
                        end
                        fprintf('\n,',i)
                end
                
                D_.CVE_raw = CVE;
                D_.CVE_shuffledL_raw = CVE_shuffledCIL;
                D_.CVE_shuffledH_raw = CVE_shuffledCIH;
                
                
                CVE=CVE(~isempty_cell(CVE));
                SigDecodeFracDraws=SigDecodeFracDraws(~isempty_cell(SigDecodeFracDraws));

                D_.CVE_shuffledCILmean=nanmean(cell2mat(CVE_shuffledCIL(~isempty_cell(CVE_shuffledCIL))));
                D_.CVE_shuffledCIHmean=nanmean(cell2mat(CVE_shuffledCIH(~isempty_cell(CVE_shuffledCIH))));
                
                D_.CVE      = nanmean(cell2mat(CVE),1);
                D_.CVE_SEM  = nansem(cell2mat(CVE),1);
                D_.SigDecodeTime = nanmean(nanmean(SigDecodeTime,1))*bw;
                D_.SigDecodeFracDraws = nansum(cell2mat(SigDecodeFracDraws),1)./numel(SigDecodeFracDraws);
                figure; hold on 
                plot(D_.CVE_shuffledCILmean,'r');plot(D_.CVE_shuffledCIHmean,'r');
                plot(D_.CVE+D_.CVE_SEM,'k');plot(D_.CVE-D_.CVE_SEM,'k');
                plot(D_.SigDecodeFracDraws,'g')
                clear  evt0 FR cve0 b k cve0 i CVEbs cves nrep trCVE

            end
            if exist('D_')
                eval(sprintf('D.%s.DelaySpan =  D_;',Delays_{iDelay}))
            end
            clear D_
            if RunErrors
            %%%% Cross-validation error decoder (error trials)
            if nLe>=2 && nRe>=2 && nL>=2 && nR>=2
                fprintf('File [%d/%d], %s Units, %s delay: Decoding error trials...\n',iFile,length(fileList),Areas{iArea},Delays_{iDelay})
                FR = [Ltrials;Rtrials;LtrialsE;RtrialsE];
                evt0 = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
                Ltr = length(tb_);
                switch UnitSelection
                    case {'groups','all'} % decode from drawn groupings of units                       
                        %if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
                        if strcmp(UnitSelection,'all')
                            groupSize = nu(iArea);
                            nrep=1;
                        else
                            nrep=10;
                        end
                         
                        CVE = cell(nrep,1);
                        PredErr = cell(nrep,1);
                        FR = [LtrialsE;RtrialsE];
                        evt0 = [ones(nLe,1);2*ones(nRe,1)];
                        FR_ = [Ltrials;Rtrials;LtrialsE;RtrialsE];
                        evt0_ = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
                        
                        if strcmp(UnitSelection,'groups')
                            parfor i=1:nrep
                                CVE{i,1} = DecodeCVE(FR(:,rs{i}),evt0,0.05);
                                PredErr{i,1} = DecodePredErr(FR_(:,rs{i}),evt0_,0.05);
                            end
                        elseif strcmp(UnitSelection,'all')
                            CVE{1} = DecodeCVE(FR,evt0,0.05);
                            PredErr{1} = DecodePredErr(FR_,evt0_,0.05);
                        end
                    case 'pairs' % Average decoding across all pairs of units
                        % CVE on error trials 
                        FR = [LtrialsE;RtrialsE];
                        evt0 = [ones(nLe,1);2*ones(nRe,1)];
                        CVE = cell(nu(iArea),nu(iArea));
                        fprintf('Decoding pair: ')
                        if iArea<3
                            for i =1:nu(iArea)
                                fprintf('%d,',i)
                                parfor j =i+1:nu(iArea)
                                    CVE{i,j} = DecodeCVE(FR(:,[i,j]),evt0,0.05);
                                end
                            end
                        else
                            for i =1:nu(1)
                                fprintf('%d,',i)
                                parfor j = nu(1)+1:nu(3)
                                    CVE{i,j} = DecodeCVE(FR(:,[i,j]),evt0,0.05);
                                end
                            end
                        end
                        fprintf('\n,',i)                          
                        % Prediction error on error trials 
                        FR = [Ltrials;Rtrials;LtrialsE;RtrialsE];
                        evt0 = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
                        PredErr = cell(nu(iArea),nu(iArea));
                        fprintf('Decoding pair: ')
                        if iArea<3
                            for i =1:nu(iArea)
                                fprintf('%d,',i)
                                parfor j =i+1:nu(iArea)
                                    PredErr{i,j} = DecodePredErr(FR(:,[i,j]),evt0,0.05);
                                end
                            end
                        else
                            for i =1:nu(1)
                                fprintf('%d,',i)
                                parfor j =nu(1)+1:nu(3)
                                    PredErr{i,j} = DecodePredErr(FR(:,[i,j]),evt0,0.05);
                                end
                            end  
                        end
                        fprintf('\n,',i)
                end
                
                D_.CVE_raw = CVE;
                D_.CVE = nanmean(cell2mat(CVE(:)),1);
                D_.CVE_SEM = nansem(cell2mat(CVE(:)),1);
                D_.PredErr_raw = PredErr;
                D_.PredErr = nanmean(cell2mat(PredErr(:)),1);
                D_.PredErr_SEM = nansem(cell2mat(PredErr(:)),1);
                clear rs evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
            end
            
            if exist('D_')
                eval(sprintf('D.%s.DelaySpan_err =  D_;',Delays_{iDelay}))
            end
            clear D_
        end
            clear Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
        end
        try 
            figure; hold on
            plot(tbShort,smooth2a(1-D.Short.DelaySpan.CVE,1,10))
            plot(tbMedium,smooth2a(1-D.Medium.DelaySpan.CVE,1,10))
            plot(tbLong,smooth2a(1-D.Long.DelaySpan.CVE,1,10))
        end
        %% Save results
        if exist('D')
            
            fnOut = sprintf('%sMixedSelectivity%s%s_%s_DelayDecoding_UnitSpan_%s.mat',pat,filesep,fname,Areas{iArea},UnitSelection);
            save(fnOut,'D','tbShort','tbMedium','tbLong','-v7.3');
            fprintf('Done.\n')
            clear D
        end
    end
end
 clearvars -except pat2 Nbs AssemblyChoice color_ Delays__ Delays_ bw tlimsShort tlimsMedium UnitSelection tlimsLong tbShort tbAll tbShort tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin %% Batch Process assemblies

%% batch import unit data for meta-analysis and plotting
clear D_
for iArea = 1:3
    try
        for iFile = 1:length(fileList)
            fname=strtok(fileList(iFile).name,'_');
            fnIn = sprintf('%sMixedSelectivity%sDelayDecodingVsBehav%s%s_%s_DelayDecoding_UnitSpan_%s.mat',pat,filesep,filesep,fname,Areas{iArea},UnitSelection);
            D_{iArea}{iFile} = load(fnIn);
            fprintf('Done.\n')
        end
    end
end

%% Collapse - Short delays
clear D_Collapsed
iArea=2

for iFile = 1:length(fileList)
    for iDelay =1:3
        D_Collapsed.meanDecoding.Short(iFile,:) = D_{iArea}{iFile}.D.Short.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Medium(iFile,:) = D_{iArea}{iFile}.D.Medium.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Long(iFile,:) = D_{iArea}{iFile}.D.Long.DelaySpan.CVE;
        
        D_Collapsed.meanDecoding.ShortCIL(iFile,:)  = nanmean(cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledL_raw(:)),1);
        D_Collapsed.meanDecoding.MediumCIL(iFile,:) = nanmean(cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledL_raw(:)),1);
        D_Collapsed.meanDecoding.LongCIL(iFile,:)   = nanmean(cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledL_raw(:)),1);
        
        D_Collapsed.meanDecoding.ShortCIH(iFile,:)  = nanmean(cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledH_raw(:)),1);
        D_Collapsed.meanDecoding.MediumCIH(iFile,:) = nanmean(cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledH_raw(:)),1);
        D_Collapsed.meanDecoding.LongCIH(iFile,:)   = nanmean(cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledH_raw(:)),1);
        
        D_Collapsed.meanDecoding.tbShort = D_{iArea}{iFile}.tbShort;
        D_Collapsed.meanDecoding.tbMedium = D_{iArea}{iFile}.tbMedium;
        D_Collapsed.meanDecoding.tbLong = D_{iArea}{iFile}.tbLong;
    end
    
        D_Collapsed.SigDecodeTimeDelayMeans(iFile) = nanmean([D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime ./diff(tlimsShort), ...
                                                   D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeTime ./diff(tlimsMedium), ...
                                                   D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeTime   ./diff(tlimsLong)]);
        D_Collapsed.SigDecodeFracDrawsDelayMeans(iFile) = nanmean([sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsShort) , ...
                                                          sum(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)./diff(tlimsMedium), ...
                                                          sum(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong)]);
                                      
        D_Collapsed.AccuracyDelayMeans(iFile) = mean(D_{1}{iFile}.D.Accuracy);
        D_Collapsed.pCorrDelayMeans(iFile) = mean(D_{1}{iFile}.D.pCorr);

        D_Collapsed.SigDecodeTime(iFile)        = D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime./diff(tlimsShort);
        D_Collapsed.SigDecodeFracDraws(iFile)   = sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsShort);
        D_Collapsed.LastSigTime(iFile)          = max([0,find(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)])*bw;
                                      
        D_Collapsed.Accuracy(iFile) = D_{1}{iFile}.D.Accuracy(1);
        D_Collapsed.pCorr(iFile)    = D_{1}{iFile}.D.pCorr(1);
        
   D_Collapsed.AboveChance(iFile,:) = D_{1}{iFile}.D.AboveChance;
end
%% Plot performance vs decoding - Short delays
figure
    subplot(1,2,1); hold on

    x = D_Collapsed.LastSigTime;
    y = D_Collapsed.pCorr;
    % y = D_Collapsed.Accuracy;
    scatter(x(~D_Collapsed.AboveChance(:,1)),y(~D_Collapsed.AboveChance(:,1)),'r')
    scatter(x(D_Collapsed.AboveChance(:,1)),y(D_Collapsed.AboveChance(:,1)),'r','filled')
    [R,P]= corrcoef([x;y]')

    mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
    pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
    slope(iArea) = mdl{iArea}.Coefficients{2,1};
    intercept(iArea) = mdl{iArea}.Coefficients{1,1};

    plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
    [pFit(iArea), Rsq(iArea)]

    h = plot(mdl{iArea});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'k';
    h(3).Color = 'k';
    h(4).Color = 'k';
    axis([0 5 0.4 1])
    legend off
%     xlabel({'Last significant timepoint (s)'})
%     ylabel('Fraction of choices correct')
%     title(Areas{iArea})
    xlabel('')
    ylabel('')
    title('')
    set(gca,'YTick',[0.4 0.6 0.8 1])
    axis square

subplot(1,2,2); hold on
    x = D_Collapsed.SigDecodeTime*4;
    % x = D_Collapsed.SigDecodeFracDraws;
    y = D_Collapsed.pCorr;
    % y = D_Collapsed.Accuracy;
    scatter(x(~D_Collapsed.AboveChance(:,2)),y(~D_Collapsed.AboveChance(:,2)),'r')
    scatter(x(D_Collapsed.AboveChance(:,2)),y(D_Collapsed.AboveChance(:,2)),'r','filled')
    [R,P]= corrcoef([x;y]');

    mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
    pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
    slope(iArea) = mdl{iArea}.Coefficients{2,1};
    intercept(iArea) = mdl{iArea}.Coefficients{1,1};

    plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
    [pFit(iArea), Rsq(iArea)]

    h = plot(mdl{iArea});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'k';
    h(3).Color = 'k';
    h(4).Color = 'k';
    axis([0 5 0.4 1])
    legend off
%     xlabel({'Significant cue decoding span (s)'})
%     ylabel('Fraction of choices correct')
%     title(Areas{iArea})
    xlabel('')
    ylabel('')
    title('')
    set(gca,'YTick',[0.4 0.6 0.8 1])
    axis square
%% Collapse - Medium delays
clear D_Collapsed
iArea=2

for iFile = 1:length(fileList)
    for iDelay =1:3
        D_Collapsed.meanDecoding.Short(iFile,:) = D_{iArea}{iFile}.D.Short.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Medium(iFile,:) = D_{iArea}{iFile}.D.Medium.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Long(iFile,:) = D_{iArea}{iFile}.D.Long.DelaySpan.CVE;
        
        D_Collapsed.meanDecoding.ShortCIL(iFile,:)  = nanmean(cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledL_raw(:)),1);
        D_Collapsed.meanDecoding.MediumCIL(iFile,:) = nanmean(cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledL_raw(:)),1);
        D_Collapsed.meanDecoding.LongCIL(iFile,:)   = nanmean(cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledL_raw(:)),1);
        
        D_Collapsed.meanDecoding.ShortCIH(iFile,:)  = nanmean(cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledH_raw(:)),1);
        D_Collapsed.meanDecoding.MediumCIH(iFile,:) = nanmean(cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledH_raw(:)),1);
        D_Collapsed.meanDecoding.LongCIH(iFile,:)   = nanmean(cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledH_raw(:)),1);
        
        D_Collapsed.meanDecoding.tbShort = D_{iArea}{iFile}.tbShort;
        D_Collapsed.meanDecoding.tbMedium = D_{iArea}{iFile}.tbMedium;
        D_Collapsed.meanDecoding.tbLong = D_{iArea}{iFile}.tbLong;
    end
    
        D_Collapsed.SigDecodeTimeDelayMeans(iFile) = nanmean([D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime ./diff(tlimsShort), ...
                                                              D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeTime ./diff(tlimsMedium), ...
                                                              D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeTime   ./diff(tlimsLong)]);
        D_Collapsed.SigDecodeFracDrawsDelayMeans(iFile) = nanmean([sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsShort) , ...
                                                                   sum(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)./diff(tlimsMedium), ...
                                                                   sum(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong)]);
                                      
        D_Collapsed.AccuracyDelayMeans(iFile) = mean(D_{1}{iFile}.D.Accuracy);
        D_Collapsed.pCorrDelayMeans(iFile) = mean(D_{1}{iFile}.D.pCorr);

        D_Collapsed.SigDecodeTime(iFile)        = D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeTime./diff(tlimsMedium);
        D_Collapsed.SigDecodeFracDraws(iFile)   = sum(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)./diff(tlimsMedium);
        D_Collapsed.LastSigTime(iFile)          = max([0,find(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)])*bw;
                                      
        D_Collapsed.Accuracy(iFile) = D_{1}{iFile}.D.Accuracy(2);
        D_Collapsed.pCorr(iFile)    = D_{1}{iFile}.D.pCorr(2);
        
   D_Collapsed.AboveChance(iFile,:) = D_{1}{iFile}.D.AboveChance;
end
%% Plot performance vs decoding - Medium delays
figure
    subplot(1,2,1); hold on

    x = D_Collapsed.LastSigTime;
    y = D_Collapsed.pCorr;
    % y = D_Collapsed.Accuracy;
    scatter(x(~D_Collapsed.AboveChance(:,2)),y(~D_Collapsed.AboveChance(:,2)),'r')
    scatter(x(D_Collapsed.AboveChance(:,2)),y(D_Collapsed.AboveChance(:,2)),'r','filled')
    [R,P]= corrcoef([x;y]')

    mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
    pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
    slope(iArea) = mdl{iArea}.Coefficients{2,1};
    intercept(iArea) = mdl{iArea}.Coefficients{1,1};

    plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
    [pFit(iArea), Rsq(iArea)]

    h = plot(mdl{iArea});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'k';
    h(3).Color = 'k';
    h(4).Color = 'k';
    axis([0 9 0.4 1])
    legend off
%     xlabel({'Last significant timepoint (s)'})
%     ylabel('Fraction of choices correct')
%     title(Areas{iArea})
    xlabel('')
    ylabel('')
    title('')
    set(gca,'YTick',[0.4 0.6 0.8 1])
    axis square

subplot(1,2,2); hold on
    x = D_Collapsed.SigDecodeTime*8;
    % x = D_Collapsed.SigDecodeFracDraws;
    y = D_Collapsed.pCorr;
    % y = D_Collapsed.Accuracy;
    scatter(x(~D_Collapsed.AboveChance(:,2)),y(~D_Collapsed.AboveChance(:,2)),'r')
    scatter(x(D_Collapsed.AboveChance(:,2)),y(D_Collapsed.AboveChance(:,2)),'r','filled')
    [R,P]= corrcoef([x;y]');

    mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
    pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
    slope(iArea) = mdl{iArea}.Coefficients{2,1};
    intercept(iArea) = mdl{iArea}.Coefficients{1,1};

    plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
    [pFit(iArea), Rsq(iArea)]

    h = plot(mdl{iArea});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'k';
    h(3).Color = 'k';
    h(4).Color = 'k';
   axis([0 9 0.4 1])
    legend off
%     xlabel({'Significant cue decoding span (s)'})
%     ylabel('Fraction of choices correct')
%     title(Areas{iArea})
    xlabel('')
    ylabel('')
    title('')
    set(gca,'YTick',[0.4 0.6 0.8 1])
    axis square  
%% Collapse - Long delays
clear D_Collapsed
iArea=2

for iFile = 1:length(fileList)
    for iDelay =1:3
        D_Collapsed.meanDecoding.Short(iFile,:) = D_{iArea}{iFile}.D.Short.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Medium(iFile,:) = D_{iArea}{iFile}.D.Medium.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Long(iFile,:) = D_{iArea}{iFile}.D.Long.DelaySpan.CVE;
        
        D_Collapsed.meanDecoding.ShortCIL(iFile,:)  = nanmean(cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledL_raw(:)),1);
        D_Collapsed.meanDecoding.MediumCIL(iFile,:) = nanmean(cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledL_raw(:)),1);
        D_Collapsed.meanDecoding.LongCIL(iFile,:)   = nanmean(cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledL_raw(:)),1);
        
        D_Collapsed.meanDecoding.ShortCIH(iFile,:)  = nanmean(cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledH_raw(:)),1);
        D_Collapsed.meanDecoding.MediumCIH(iFile,:) = nanmean(cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledH_raw(:)),1);
        D_Collapsed.meanDecoding.LongCIH(iFile,:)   = nanmean(cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledH_raw(:)),1);
        
        D_Collapsed.meanDecoding.tbShort = D_{iArea}{iFile}.tbShort;
        D_Collapsed.meanDecoding.tbMedium = D_{iArea}{iFile}.tbMedium;
        D_Collapsed.meanDecoding.tbLong = D_{iArea}{iFile}.tbLong;
    end
    
        D_Collapsed.SigDecodeTimeDelayMeans(iFile) = nanmean([D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime ./diff(tlimsShort), ...
                                                   D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeTime ./diff(tlimsMedium), ...
                                                   D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeTime   ./diff(tlimsLong)]);
        D_Collapsed.SigDecodeFracDrawsDelayMeans(iFile) = nanmean([sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsShort) , ...
                                                          sum(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)./diff(tlimsMedium), ...
                                                          sum(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong)]);
                                      
        D_Collapsed.AccuracyDelayMeans(iFile) = mean(D_{1}{iFile}.D.Accuracy);
        D_Collapsed.pCorrDelayMeans(iFile) = mean(D_{1}{iFile}.D.pCorr);

        D_Collapsed.SigDecodeTime(iFile)        = D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeTime./diff(tlimsLong);
        D_Collapsed.SigDecodeFracDraws(iFile)   = sum(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong);
        D_Collapsed.LastSigTime(iFile)          = max([0,find(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)])*bw;
                                      
        D_Collapsed.Accuracy(iFile) = D_{1}{iFile}.D.Accuracy(3);
        D_Collapsed.pCorr(iFile)    = D_{1}{iFile}.D.pCorr(3);
        
   D_Collapsed.AboveChance(iFile,:) = D_{1}{iFile}.D.AboveChance;
end
%% Plot performance vs decoding - Long delays
figure
    subplot(1,2,1); hold on

    x = D_Collapsed.LastSigTime;
    y = D_Collapsed.pCorr;
    % y = D_Collapsed.Accuracy;
    scatter(x(~D_Collapsed.AboveChance(:,3)),y(~D_Collapsed.AboveChance(:,3)),'r')
    scatter(x(D_Collapsed.AboveChance(:,3)),y(D_Collapsed.AboveChance(:,3)),'r','filled')
    [R,P]= corrcoef([x;y]')

    mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
    pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
    slope(iArea) = mdl{iArea}.Coefficients{2,1};
    intercept(iArea) = mdl{iArea}.Coefficients{1,1};

    plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
    [pFit(iArea), Rsq(iArea)]

    h = plot(mdl{iArea});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'k';
    h(3).Color = 'k';
    h(4).Color = 'k';
    axis([0 17 0.4 1])
    legend off
%     xlabel({'Last significant timepoint (s)'})
%     ylabel('Fraction of choices correct')
%     title(Areas{iArea})
    xlabel('')
    ylabel('')
    title('')
    set(gca,'YTick',[0.4 0.6 0.8 1])
    axis square

subplot(1,2,2); hold on
    x = D_Collapsed.SigDecodeTime*16;
    % x = D_Collapsed.SigDecodeFracDraws;
    y = D_Collapsed.pCorr;
    % y = D_Collapsed.Accuracy;
    scatter(x(~D_Collapsed.AboveChance(:,3)),y(~D_Collapsed.AboveChance(:,3)),'r')
    scatter(x(D_Collapsed.AboveChance(:,3)),y(D_Collapsed.AboveChance(:,3)),'r','filled')
    [R,P]= corrcoef([x;y]');

    mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
    pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
    slope(iArea) = mdl{iArea}.Coefficients{2,1};
    intercept(iArea) = mdl{iArea}.Coefficients{1,1};

    plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
    [pFit(iArea), Rsq(iArea)]

    h = plot(mdl{iArea});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'k';
    h(3).Color = 'k';
    h(4).Color = 'k';
    axis([0 9 0.4 1])
    legend off
%     xlabel({'Significant cue decoding span (s)'})
%     ylabel('Fraction of choices correct')
%     title(Areas{iArea})
    xlabel('')
    ylabel('')
    title('')
    set(gca,'YTick',[0.4 0.6 0.8 1])
    axis square  
     
    
%% Plot mean performance - all recordings
figure
subplot(3,1,1); hold on
Ltr = length(D_Collapsed.meanDecoding.tbShort);

ciplot(nanmean(D_Collapsed.meanDecoding.Short)+nansem(D_Collapsed.meanDecoding.Short),...
       nanmean(D_Collapsed.meanDecoding.Short)-nansem(D_Collapsed.meanDecoding.Short),...
       D_Collapsed.meanDecoding.tbShort,'k')
 
ciplot(nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),...
       nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),...
       D_Collapsed.meanDecoding.tbShort,'r')
axis([0 16 0.4 0.8]) 

subplot(3,1,2); hold on
Ltr = length(D_Collapsed.meanDecoding.tbMedium);

ciplot(nanmean(D_Collapsed.meanDecoding.Medium)+nansem(D_Collapsed.meanDecoding.Medium),...
       nanmean(D_Collapsed.meanDecoding.Medium)-nansem(D_Collapsed.meanDecoding.Medium),...
       D_Collapsed.meanDecoding.tbMedium,'k')
ciplot(nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),...
       nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),...
       D_Collapsed.meanDecoding.tbMedium,'r')
axis([0 16 0.4 0.8]) 

subplot(3,1,3); hold on
Ltr = length(D_Collapsed.meanDecoding.tbLong);

ciplot(nanmean(D_Collapsed.meanDecoding.Long)+nansem(D_Collapsed.meanDecoding.Long),...
        nanmean(D_Collapsed.meanDecoding.Long)-nansem(D_Collapsed.meanDecoding.Long),...
     D_Collapsed.meanDecoding.tbLong,'k')
  
 ciplot(nanmean((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2),...
       nanmean((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2),...
       D_Collapsed.meanDecoding.tbLong,'r')
axis([0 16 0.4 0.8]) 
%% Plot mean performance - all recordings overlaid
figure
hold on
Ltr = length(D_Collapsed.meanDecoding.tbLong);
ciplot(nanmean(D_Collapsed.meanDecoding.Long)+nansem(D_Collapsed.meanDecoding.Long),...
        nanmean(D_Collapsed.meanDecoding.Long)-nansem(D_Collapsed.meanDecoding.Long),...
     D_Collapsed.meanDecoding.tbLong,'k',0.7)
ciplot(nanmean((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2),...
       nanmean((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2),...
       D_Collapsed.meanDecoding.tbLong,'k')

Ltr = length(D_Collapsed.meanDecoding.tbMedium);
ciplot(nanmean(D_Collapsed.meanDecoding.Medium)+nansem(D_Collapsed.meanDecoding.Medium),...
       nanmean(D_Collapsed.meanDecoding.Medium)-nansem(D_Collapsed.meanDecoding.Medium),...
       D_Collapsed.meanDecoding.tbMedium,'g',0.9)
ciplot(nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),...
       nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),...
       D_Collapsed.meanDecoding.tbMedium,'g')
   
Ltr = length(D_Collapsed.meanDecoding.tbShort);
ciplot(nanmean(D_Collapsed.meanDecoding.Short)+nansem(D_Collapsed.meanDecoding.Short),...
       nanmean(D_Collapsed.meanDecoding.Short)-nansem(D_Collapsed.meanDecoding.Short),...
       D_Collapsed.meanDecoding.tbShort,'r',0.9)
ciplot(nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),...
       nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),...
       D_Collapsed.meanDecoding.tbShort,'r')

plot([0 16],[0.5 0.5],':k','LineWidth',1.5)
axis([0 16 0.4 0.8]) 
%% Plot mean performance - all recordings overlaid (each recording)
figure
hold on

plot(D_Collapsed.meanDecoding.tbLong,smooth2a(D_Collapsed.meanDecoding.Long,0,2),'color',[0 0 0 0.5],'Linewidth',1)
plot(D_Collapsed.meanDecoding.tbShort,smooth2a(D_Collapsed.meanDecoding.Short,0,2),'color',[1 0 0 0.5],'Linewidth',1)
plot(D_Collapsed.meanDecoding.tbMedium,smooth2a(D_Collapsed.meanDecoding.Medium,0,2),'color',[0 1 0 0.5],'Linewidth',1)
plot([0 16],[0.5 0.5],':k','LineWidth',1.5)
axis([0 16 0 1]) 

%%
twin = [3,4];
idx = find(D_Collapsed.meanDecoding.tbShort>=twin(1) & D_Collapsed.meanDecoding.tbShort<twin(2));
short_ = nanmean(D_Collapsed.meanDecoding.Short(:,idx),2);

idx = find(D_Collapsed.meanDecoding.tbMedium>=twin(1) & D_Collapsed.meanDecoding.tbMedium<twin(2));
medium_ = nanmean(D_Collapsed.meanDecoding.Medium(:,idx),2);

idx = find(D_Collapsed.meanDecoding.tbLong>=twin(1) & D_Collapsed.meanDecoding.tbLong<twin(2));
long_ = nanmean(D_Collapsed.meanDecoding.Long(:,idx),2);

figure; hold on
% scatter(1*ones(size(short_)),short_,'k');
% scatter(2*ones(size(medium_)),medium_,'k');
% scatter(3*ones(size(long_)),long_,'k');

plot([short_,medium_,long_]',':ok')

e = errorbar([1,2,3],nanmean([short_,medium_,long_]),nansem([short_,medium_,long_]),'o-k','LineWidth',1.5)
plot([0 16],[0.5 0.5],':k','LineWidth',1.5)
set(gca,'Xtick',[1:3],'Ytick',[0:0.25:1],'XtickLabel',{'4s','8s','16'})
axis([0.5 3.5 0 1]) 

friedman([short_,medium_,long_],1)
%%
twin = [7.5,8];
idx = find(D_Collapsed.meanDecoding.tbLong>=twin(1) & D_Collapsed.meanDecoding.tbLong<twin(2));
long_ = nanmean(D_Collapsed.meanDecoding.Long(:,idx),2);

idx = find(D_Collapsed.meanDecoding.tbMedium>=twin(1) & D_Collapsed.meanDecoding.tbMedium<twin(2));
medium_ = nanmean(D_Collapsed.meanDecoding.Medium(:,idx),2);
figure; hold on
scatter(medium_,long_,'k','filled');
plot([0,1],[0,1],':k','LineWidth',1.5)
% e = errorbar(mean(medium_),mean(long_),nansem(medium_),nansem(medium_),nansem(long_),nansem(long_),'o');
% e.Color='k';
% e.LineWidth=1.5
set(gca,'Xtick',[0:0.25:1],'Ytick',[0:0.25:1])
% [h,p,ci,stats] = ttest(long_,medium_)
[p,h,stats] = signrank(long_,medium_,'method','approximate')

% ciplot(nanmean(D_Collapsed.meanDecoding.Long)+nansem(D_Collapsed.meanDecoding.Long),...
%         nanmean(D_Collapsed.meanDecoding.Long)-nansem(D_Collapsed.meanDecoding.Long),...
%      D_Collapsed.meanDecoding.tbLong,'k',0.7)
% ciplot(nanmean((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2),...
%        nanmean((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.LongCIL(:,1:Ltr)+D_Collapsed.meanDecoding.LongCIH(:,1:Ltr))./2),...
%        D_Collapsed.meanDecoding.tbLong,'k')
% 
% Ltr = length(D_Collapsed.meanDecoding.tbMedium);
% ciplot(nanmean(D_Collapsed.meanDecoding.Medium)+nansem(D_Collapsed.meanDecoding.Medium),...
%        nanmean(D_Collapsed.meanDecoding.Medium)-nansem(D_Collapsed.meanDecoding.Medium),...
%        D_Collapsed.meanDecoding.tbMedium,'g',0.9)
% ciplot(nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),...
%        nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),...
%        D_Collapsed.meanDecoding.tbMedium,'g')
%    
% Ltr = length(D_Collapsed.meanDecoding.tbShort);
% ciplot(nanmean(D_Collapsed.meanDecoding.Short)+nansem(D_Collapsed.meanDecoding.Short),...
%        nanmean(D_Collapsed.meanDecoding.Short)-nansem(D_Collapsed.meanDecoding.Short),...
%        D_Collapsed.meanDecoding.tbShort,'r',0.9)
% ciplot(nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),...
%        nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),...
%        D_Collapsed.meanDecoding.tbShort,'r')

%% Plot mean performance - above chance separate
figure
subplot(3,1,1); hold on
Ltr = length(D_Collapsed.meanDecoding.tbShort);

ciplot(nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),...
       nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),...
       D_Collapsed.meanDecoding.tbShort,'k')
plot(D_Collapsed.meanDecoding.tbShort,nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) + ...
                                      nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),'color',0.6*[1 1 1],'LineWidth',1.2,'LineStyle','--')
plot(D_Collapsed.meanDecoding.tbShort,nanmean((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2) - ...
                                      nansem((D_Collapsed.meanDecoding.ShortCIL(:,1:Ltr)+D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./2),'color',0.6*[1 1 1],'LineWidth',1.2,'LineStyle','--')

                                  
ciplot(nanmean(D_Collapsed.meanDecoding.Short)+nansem(D_Collapsed.meanDecoding.Short),...
       nanmean(D_Collapsed.meanDecoding.Short)-nansem(D_Collapsed.meanDecoding.Short),...
       D_Collapsed.meanDecoding.tbShort,'r',1)
plot(D_Collapsed.meanDecoding.tbShort,nanmean(D_Collapsed.meanDecoding.Short)+nansem(D_Collapsed.meanDecoding.Short),'color','r','LineWidth',1.2)
plot(D_Collapsed.meanDecoding.tbShort,nanmean(D_Collapsed.meanDecoding.Short)-nansem(D_Collapsed.meanDecoding.Short),'color','r','LineWidth',1.2)

z = sum(D_Collapsed.meanDecoding.Short(:,1:Ltr)>D_Collapsed.meanDecoding.ShortCIH(:,1:Ltr))./size(D_Collapsed.meanDecoding.Short(:,1:Ltr),1); 
x = D_Collapsed.meanDecoding.tbShort;
y = [0.95*ones(size(x));1*ones(size(x))]; 
x=[x;x];
z=[z;z];

pc = pcolor(x,y,z);
pc.EdgeColor = 'none';
caxis([0 1]) 
colormap(flipud(gray))

plot([4 4],[0.4 0.9],'LineWidth',1.5,'LineStyle',':','color',0.6*[1 1 1])
set(gca,'XTick',[0 4],'XTickLabel',{'0' '4s'})
axis([0 17 0.3 1]) 

subplot(3,1,2); hold on
Ltr = length(D_Collapsed.meanDecoding.tbMedium);

ciplot(nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) + nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),...
       nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) - nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),...
       D_Collapsed.meanDecoding.tbMedium,'k')
plot(D_Collapsed.meanDecoding.tbMedium,nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) + ...
                                      nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),'color',0.6*[1 1 1],'LineWidth',1.2,'LineStyle','--')
plot(D_Collapsed.meanDecoding.tbMedium ,nanmean((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2) - ...
                                      nansem((D_Collapsed.meanDecoding.MediumCIL(:,1:Ltr)+D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./2),'color',0.6*[1 1 1],'LineWidth',1.2,'LineStyle','--')
                                  
aboveChance_ = D_Collapsed.AboveChance(:,2);
ciplot(nanmean(D_Collapsed.meanDecoding.Medium(aboveChance_,:))+nansem(D_Collapsed.meanDecoding.Medium(aboveChance_,:)),...
       nanmean(D_Collapsed.meanDecoding.Medium(aboveChance_,:))-nansem(D_Collapsed.meanDecoding.Medium(aboveChance_,:)),...
       D_Collapsed.meanDecoding.tbMedium,'r',1)
plot(D_Collapsed.meanDecoding.tbMedium,nanmean(D_Collapsed.meanDecoding.Medium)+nansem(D_Collapsed.meanDecoding.Medium),'color','r','LineWidth',1.2)
plot(D_Collapsed.meanDecoding.tbMedium,nanmean(D_Collapsed.meanDecoding.Medium)-nansem(D_Collapsed.meanDecoding.Medium),'color','r','LineWidth',1.2)   

z = sum(D_Collapsed.meanDecoding.Medium(:,1:Ltr)>D_Collapsed.meanDecoding.MediumCIH(:,1:Ltr))./size(D_Collapsed.meanDecoding.Medium(:,1:Ltr),1); 
x = D_Collapsed.meanDecoding.tbMedium;
y = [0.95*ones(size(x));1*ones(size(x))]; 
x=[x;x];
z=[z;z];

pc = pcolor(x,y,z);
pc.EdgeColor = 'none';
caxis([0 1]) 
colormap(flipud(gray))

plot([8 8],[0.4 0.9],'LineWidth',1.5,'LineStyle',':','color',0.6*[1 1 1])   
set(gca,'XTick',[0 8],'XTickLabel',{'0' '8s'})
axis([0 17 0.3 1]) 

subplot(3,1,3); hold on
Ltr = length(D_Collapsed.meanDecoding.tbLong);

ciplot(nanmean((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2) + nansem((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2),...
       nanmean((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2) - nansem((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2),...
       D_Collapsed.meanDecoding.tbLong,'k')
   
plot(D_Collapsed.meanDecoding.tbLong ,nanmean((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2) + ...
                                      nansem((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2),'color',0.6*[1 1 1],'LineWidth',1.2,'LineStyle','--')
plot(D_Collapsed.meanDecoding.tbLong ,nanmean((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2) - ...
                                      nansem((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2),'color',0.6*[1 1 1],'LineWidth',1.2,'LineStyle','--')

aboveChance_ = D_Collapsed.AboveChance(:,3);
ciplot(nanmean(D_Collapsed.meanDecoding.Long(aboveChance_,:))+nansem(D_Collapsed.meanDecoding.Long(aboveChance_,:)),...
       nanmean(D_Collapsed.meanDecoding.Long(aboveChance_,:))-nansem(D_Collapsed.meanDecoding.Long(aboveChance_,:)),...
       D_Collapsed.meanDecoding.tbLong,'r',1)
plot(D_Collapsed.meanDecoding.tbLong,nanmean(D_Collapsed.meanDecoding.Long(aboveChance_,:))+nansem(D_Collapsed.meanDecoding.Long(aboveChance_,:)),'color','r','LineWidth',1.2)
plot(D_Collapsed.meanDecoding.tbLong,nanmean(D_Collapsed.meanDecoding.Long(aboveChance_,:))-nansem(D_Collapsed.meanDecoding.Long(aboveChance_,:)),'color','r','LineWidth',1.2) 
   
ciplot(nanmean(D_Collapsed.meanDecoding.Long(~aboveChance_,:))+nansem(D_Collapsed.meanDecoding.Long(~aboveChance_,:)),...
       nanmean(D_Collapsed.meanDecoding.Long(~aboveChance_,:))-nansem(D_Collapsed.meanDecoding.Long(~aboveChance_,:)),...
       D_Collapsed.meanDecoding.tbLong,'r',0.1)
   
%    plot( D_Collapsed.meanDecoding.tbLong,D_Collapsed.meanDecoding.Long(aboveChance_,:),'b')
%    plot( D_Collapsed.meanDecoding.tbLong,D_Collapsed.meanDecoding.Long(~aboveChance_,:),'r')
plot(D_Collapsed.meanDecoding.tbLong,nanmean(D_Collapsed.meanDecoding.Long(~aboveChance_,:))+nansem(D_Collapsed.meanDecoding.Long(~aboveChance_,:)),'color',[1 0 0 0.4],'LineWidth',1.2)
plot(D_Collapsed.meanDecoding.tbLong,nanmean(D_Collapsed.meanDecoding.Long(~aboveChance_,:))-nansem(D_Collapsed.meanDecoding.Long(~aboveChance_,:)),'color',[1 0 0 0.4],'LineWidth',1.2) 
z = sum(D_Collapsed.meanDecoding.Long>D_Collapsed.meanDecoding.LongCIH)./size(D_Collapsed.meanDecoding.Long,1); 
x = D_Collapsed.meanDecoding.tbLong;
y = [0.95*ones(size(x));1*ones(size(x))]; 
x=[x;x];
z=[z;z];

pc = pcolor(x,y,z);
pc.EdgeColor = 'none';
caxis([0 1]) 
colormap(flipud(gray))

plot([16 16],[0.4 0.9],'LineWidth',1.5,'LineStyle',':','color',0.6*[1 1 1])
set(gca,'XTick',[0 16],'XTickLabel',{'0' '16s'})

axis([0 17 0.3 1]) 
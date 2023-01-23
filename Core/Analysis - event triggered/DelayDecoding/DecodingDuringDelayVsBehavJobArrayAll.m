function DecodingDuringDelayVsBehavJobArrayAll(iFile)
if isstr(iFile)
    iFile = str2num(iFile);
end
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
        pat2 = [pat 'KDE_binsTaskOnly' filesep 'LONGTaskonly' filesep];
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end
mkdir([pat 'MixedSelectivity'])

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
%% Batch process 
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
    for iArea = 1:length(Areas)
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
                                    CVE_shuffledCIL{i,j}(1:Ltr)       = CVE_shuffledHL(:,1);
                                    CVE_shuffledCIH{i,j}(1:Ltr)       = CVE_shuffledHL(:,2);
                                    
                                end
                            end
                            
                        else
                            for i =1:nu(1)
                                fprintf('%d,',i)
                                parfor j = nu(1)+1:nu(3)
                                    %CVE{i,j} = DecodeCVE(FR(:,[i,j]),evt0,0.05);
                                    [CVE_,~,~,~,~,CVE_shuffledHL] = RunDecodeGroups([i,j],1,L_,R_,inf,'max');
                                    CVE{i,j}=CVE_';
                                    % Number of significant decoding timesteps for this draw
                                    SigDecodeTime(i,j)                = sum(CVE_>CVE_shuffledHL(:,2)); 
                                    SigDecodeFracDraws{i,j}(1:Ltr)    = CVE_>CVE_shuffledHL(:,2);
                                    CVE_shuffledCIL{i,j}(1:Ltr)       = CVE_shuffledHL(:,1);
                                    CVE_shuffledCIH{i,j}(1:Ltr)       = CVE_shuffledHL(:,2);
                                    
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
%         try 
%             figure; hold on
%             plot(tbShort,smooth2a(1-D.Short.DelaySpan.CVE,1,10))
%             plot(tbMedium,smooth2a(1-D.Medium.DelaySpan.CVE,1,10))
%             plot(tbLong,smooth2a(1-D.Long.DelaySpan.CVE,1,10))
%         end
        %% Save results
        if exist('D')
            
            fnOut = sprintf('%sMixedSelectivity%s%s_%s_DelayDecoding_UnitSpan_%s.mat',pat,filesep,fname,Areas{iArea},UnitSelection);
            save(fnOut,'D','tbShort','tbMedium','tbLong','-v7.3');
            fprintf('Done.\n')
            clear D
        end
    end

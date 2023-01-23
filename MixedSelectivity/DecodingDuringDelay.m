%%
% Compares continuous L/R discrimination during the delay period,
% calcluated on short/medium/long delay trials independently
%% %%%%%% PREAMBLE %%%%%%

clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};

Target = 'LONG';

tlimsAll = [-5 5];
tlimsShort=[-5 6];
tlimsMedium=[-5 10];
tlimsLong=[-5 20];

shift = 0;
plotOnline = false;
bw=0.05;
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
UnitSelection = 'pairs'; %{'pairs','groups'}
groupSize = 8;
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
mkdir([pat 'MixedSelectivity'])

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
%% Batch process units
for iFile =1:length(fileList)
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
    for iArea = 1:length(Areas)
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
            if nL>2 && nR>2
                fprintf('File [%d/%d], %s Units, %s delay: Decoding correct trials...\n',iFile,length(fileList),Areas{iArea},Delays_{iDelay})
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = length(tb_);
                switch UnitSelection
                    case 'groups' % decode from drawn groupings of units
                        %if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
                        nrep=10;
                        CVE = cell(nrep,1);
                        parfor i=1:nrep
                            if iArea<3
                                rs{i} = randsample(size(Ltrials,2),groupSize);
                            else
                                 rs{i} = [randsample(nu(1),floor(groupSize/2));...
                                          randsample(nu(2),ceil(groupSize/2))+nu(1)];
                            end
                            CVE{i,1} = DecodeCVE(FR(:,rs{i}),evt0,0.05);
                        end
                        
                    case 'pairs' % Average decoding across all pairs of units
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
                end
                D_.CVE_raw = CVE;
                D_.CVE = nanmean(cell2mat(CVE(:)),1);
                D_.CVE_SEM = nansem(cell2mat(CVE(:)),1);
                clear CVE
                %hold on; plot(1-D_.CVE+D_.CVE_SEM);plot(1-D_.CVE-D_.CVE_SEM);
                clear  evt0 FR cve0 b k cve0 i CVEbs cves nrep tr
            end
            if exist('D_')
                eval(sprintf('D.%s.DelaySpan =  D_;',Delays_{iDelay}))
            end
            clear D_
            
            %%%% Cross-validation error decoder (error trials)
            if nLe>2 && nRe>2 && nL>2 && nR>2
                fprintf('File [%d/%d], %s Units, %s delay: Decoding error trials...\n',iFile,length(fileList),Areas{iArea},Delays_{iDelay})
                FR = [Ltrials;Rtrials;LtrialsE;RtrialsE];
                evt0 = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
                Ltr = length(tb_);
                switch UnitSelection
                    case 'groups' % decode from drawn groupings of units
                        fprintf('File [%d/%d], %s Units, %s delay: Decoding error trials...\n',iFile,length(fileList),Areas{iArea},Delays_{iDelay})

                        %if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
                        nrep=10;
                        FR = [LtrialsE;RtrialsE];
                        evt0 = [ones(nLe,1);2*ones(nRe,1)];
                        CVE = cell(nrep,1);
                        parfor i=1:nrep
                            CVE{i,1} = DecodeCVE(FR(:,rs{i}),evt0,0.05);
                        end
                        
                        FR = [Ltrials;Rtrials;LtrialsE;RtrialsE];
                        evt0 = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
                        PredErr = cell(nrep,1);
                        parfor i=1:nrep
                            PredErr{i,1} = DecodePredErr(FR(:,rs{i}),evt0,0.05);
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
            
            clear Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
        end
        %% Save results
        if exist('D')
            
            fnOut = sprintf('%sMixedSelectivity%s%s_%s_DelayDecoding_Units_%s.mat',pat,filesep,fname,Areas{iArea},UnitSelection);
            save(fnOut,'D','tbShort','tbMedium','tbLong','-v7.3');
            fprintf('Done.\n')
            clear D
        end
    end
end
% clearvars -except pat2 Nbs AssemblyChoice color_ Delays__ Delays_ bw tlimsShort tlimsMedium tlimsLong tbShort tbAll tbShort tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin %% Batch Process assemblies
%% Batch process assemblies
for iFile =1:length(fileListAss)
    %% Get the files
    fname=strtok(fileListAss(iFile).name,'_');
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    fprintf('Analysing run %d/%d %s ...',iFile,length(fileList),fname)
    switch AssemblyChoice
        case 1
            Ass = load(sprintf('%s%s_iFR50_FSC.mat',pat2,fname));
        case 2
            Ass = load(sprintf('%s%s_iFR50_behavOnly_FSC.mat',pat2,fname));
        case 3
            Ass = load(sprintf('%s%s_iFR50_Task_FSC.mat',pat2,fname));
    end
    
    %%  Delay decoding for assemblies
    for iDelay = 1:length(Delays_)
        eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftCorrect];',Delays_{iDelay}));
        eval(sprintf('RightTrials = [t.%s.SamplePress_RightCorrect];',Delays_{iDelay}));
        eval(sprintf('LeftTrialsE = [t.%s.SamplePress_LeftError''];',Delays_{iDelay}));
        eval(sprintf('RightTrialsE = [t.%s.SamplePress_RightError''];',Delays_{iDelay}));
        eval(sprintf('tlims_X = tlims%s;',Delays_{iDelay}));
        length_ = sum(abs(tlims_X))/bw;
        eval(sprintf('tb_ = tb%s;',Delays_{iDelay}))
        for iArea = 1:length(Ass.FSCsel)
            
            FSC = Ass.FSCsel{iArea};
            if ~isempty(FSC)
                
                Ltrials = [];nL = 0;
                for iTrial =1:size(LeftTrials,1)
                    try
                        tlims_  = LeftTrials(iTrial,1)/1e6+tlims_X(1);
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                        Ltrials = [Ltrials;FSC(tlims_,:)];
                        nL=nL+1;
                    end
                end
                
                Rtrials = [];nR = 0;
                for iTrial =1:size(RightTrials,1)
                    try
                        tlims_  = RightTrials(iTrial,1)/1e6+tlims_X(1);
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                        
                        Rtrials = [Rtrials;FSC(tlims_,:)];
                        nR=nR+1;
                    end
                end
                
                LtrialsE = [];nLe = 0;
                for iTrial =1:size(LeftTrialsE,1)
                    try
                        tlims_  = LeftTrialsE(iTrial,1)/1e6+tlims_X(1);
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                        
                        LtrialsE = [LtrialsE;FSC(tlims_,:)];
                        nLe=nLe+1;
                    end
                end
                
                RtrialsE = [];nRe = 0;
                for iTrial =1:size(RightTrialsE,1)
                    try
                        tlims_  = RightTrialsE(iTrial,1)/1e6+tlims_X(1);
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                        
                        RtrialsE = [RtrialsE;FSC(tlims_,:)];
                        nRe=nRe+1;
                    end
                end
                
                
                
                % Cross-validation error decoder (correct trials)
                %                 if nL>2 && nR>2
                
                FSC_ = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = length(tb_);
                nAss_ = cellfun(@(x)size(x,2),Ass.FSCsel);
                % remove zero entries
                idx = nAss_>0;
                nAss_(~idx) = min(nAss_(idx));
                
                
                % Prediction error for cross-validation decoder
                %                     if size(Ltrials,2)==min(nAss_), nrep=1; else nrep=10; end
                %                     %                  nrep=1;
                %                     for i=1:nrep
                %                         rs{i} = randsample(size(Ltrials,2),min(nAss_));
                %                         D_.DelaySpan{iArea}.CVE(i,:) = DecodeCVE_mex(FSC_(:,rs{i}),evt0,0.05);
                %                     end
                %                     D_.DelaySpan{iArea}.CVE = nanmean(D_.DelaySpan{iArea}.CVE,1);
                
                CVE = cell(nAss_(iArea),nAss_(iArea));
                for i =1:nAss_(iArea)
                    parfor j =1:nAss_(iArea)
                        [i j]
                        if i~=j
                            CVE{i,j} = DecodeCVE_mex(FSC_(:,[i,j]),evt0,0.05);
                        end
                    end
                end
                D_.DelaySpan{iArea}.CVE_raw = CVE;
                D_.DelaySpan{iArea}.CVE = nanmean(cell2mat(CVE(:)),1);
                D_.DelaySpan{iArea}.CVE_SEM = nansem(cell2mat(CVE(:)),1);
                clear CVE
                clear rs nAss idx evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
                %                 end
                
                
                
                
                
                % Cross-validation error decoder (error trials)
                %                 if nLe>2 && nRe>2
                
                FSC_ = [Ltrials;Rtrials;LtrialsE;RtrialsE];
                evt0 = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
                Ltr = length(tb_);
                nAss_ = cellfun(@(x)size(x,2),Ass.FSCsel);
                % remove zero entries
                idx = nAss_>0;
                nAss_(~idx) = min(nAss_(idx));
                
                
                % Prediction error for cross-validation decoder
                %                     if size(Ltrials,2)==min(nAss_), nrep=1; else nrep=10; end
                %                     %                  nrep=1;
                %                     for i=1:nrep
                %                         rs{i} = randsample(size(Ltrials,2),min(nAss_));
                %                         D_.DelaySpan{iArea}.CVE_err(i,:)=DecodePredErr(FSC_(:,rs{i}),evt0,0.05);
                %                     end
                %                     D_.DelaySpan{iArea}.CVE_err = nanmean(D_.DelaySpan{iArea}.CVE_err,1);
                
                
                CVE = cell(nAss_(iArea),nAss_(iArea));
                for i =1:nAss_(iArea)
                    parfor j =1:nAss_(iArea)
                        [i j]
                        if i~=j
                            CVE{i,j} = DecodeCVE_mex(FSC_(:,[i,j]),evt0,0.05);
                        end
                    end
                end
                D_.DelaySpan{iArea}.CVE_err_raw = CVE;
                D_.DelaySpan{iArea}.CVE_err = nanmean(cell2mat(CVE(:)),1);
                D_.DelaySpan{iArea}.CVE_err_SEM = nansem(cell2mat(CVE(:)),1);
                clear CVE
                clear rs nAss idx evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
                %                 end
                
            end
        end
        if exist('D_')
            eval(sprintf('D_Ass.%s.DelaySpan =  D_.DelaySpan;',Delays_{iDelay}))
        end
        clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
    end
    
    %% Save
    if exist('D_Ass')
        
        fnOut = sprintf('%sMixedSelectivity%s%s_DelayDecoding_Ass.mat',pat,filesep,fname);
        save(fnOut,'D_Ass','tbShort','tbMedium','tbLong','-v7.3');
        fprintf('Done.\n')
        clear D_Ass
    end
end
% clearvars -except pat2 Nbs AssemblyChoice color_ Delays_ bw tbShort tbAll tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin
%% batch import unit data for meta-analysis and plotting
% Import
clear D_
for iArea = 1:length(Areas)
    
    for iFile =1:length(fileList)
        
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%sMixedSelectivity\\%s_%s_DelayDecoding_Units_%s.mat',pat,fname,Areas{iArea},UnitSelection);
        load(fnIn ,'D');
        
        if strcmp(UnitSelection,'pairs')  %% Get the assembly membership info
           
            fname=strtok(fileListAss(iFile).name,'_');
            fprintf('Loading Assemblies %d/%d %s ...\n',iFile,length(fileList),fname)
            % for assemblies describing whole task period (cue  - reward)
            switch AssemblyChoice
                case 1
                    load(sprintf('%s%s%s_iFR50_FSC.mat',pat2,filesep,fname),'units','nu');
                    %             load(sprintf('%s%s%s_iFR50_AssemRes2.mat',pat2,filesep,fname),'usel_out');
                    usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
                case 2
                    load(sprintf('%s%s_iFR50_BehavOnly_FSC.mat',pat2,fname),'units','nu','FSCsel','Tmtx');
                    load(sprintf('%s%s_iFR50_BehavOnly_AssemRes2.mat',pat2,fname),'usel_out');
                case 3
                    load(sprintf('%s%s%s_iFR50_Task_FSC.mat',pat2,filesep,fname),'units','nu');
                    %             load(sprintf('%s%s%s_iFR50_Task_AssemRes2.mat',pat2,filesep,fname),'usel_out');
                    usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
            end
            Ass.usel_out    = usel_out;
            Ass.FSC         = FSCsel;
            Ass.Tmtx        = Tmtx;
            Ass.usel_out{3} = [Ass.usel_out{1},Ass.usel_out{2}+nu(1)];
            Ass.units       = units;
            Ass.nu = nu; Ass.nu(3) = sum(Ass.nu(1:2));
            % Check that inter-area assemblies actually span the two areas
            if ~isempty(Ass.FSC{3})
                nAss = length(Ass.units{3});
                idx = false(nAss,1);
                for iAss = 1:nAss
                    U_ = Ass.usel_out{3}(Ass.units{3}{iAss});
                    
                    if  min(U_)>max(Ass.usel_out{1})
                        idx(iAss) = true;
                    end
                    Ass.units{3}(idx)=[];
                    Ass.FSC{3}(:,idx)=[];
                    
                end
                
            end
            for iArea_ = 1:3
                Ass.LocalMembers{iArea_}    = unique(cell2mat(Ass.units{iArea_}));
            end
            
            
            Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3}<=Ass.nu(1)));
            Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3}>Ass.nu(1))-Ass.nu(1));
            for iArea_ = 1:2
                Ass.JointMembers{iArea_} = setdiff(unique(cell2mat(Ass.units{iArea_})),Ass.LocalMembers{iArea_});
                Ass.NonMembers{iArea_}   = setdiff(1:Ass.nu(iArea_),[Ass.LocalMembers{iArea_},Ass.JointMembers{iArea_}]);
            end
            for iArea_ = 1:2
                Ass.NonMembers{iArea_}   = Ass.usel_out{iArea_}(Ass.NonMembers{iArea_});
                Ass.LocalMembers{iArea_} = Ass.usel_out{iArea_}(Ass.LocalMembers{iArea_});
                Ass.JointMembers{iArea_} = Ass.usel_out{iArea_}(Ass.JointMembers{iArea_});
            end
            
            for iArea_ = 1:2
                %Membership_{iArea_} = -ones( Ass.nu(iArea_),1);
                Membership_{iArea_} = -ones(nu(iArea_),1);
                Membership_{iArea_}(Ass.NonMembers{iArea_})=0;
                Membership_{iArea_}(Ass.LocalMembers{iArea_})=1;
                Membership_{iArea_}(Ass.JointMembers{iArea_})=2;
                Membership_{iArea_}(setdiff(1:nu(iArea_),Ass.usel_out{iArea_}))=[];
            end

            for iArea_ = 1:2
                % Pad this out to ensure that there's an entry even if there's no detected assembly
                LocalMembersMatrix_{iArea_} = nan(Ass.nu(iArea_),length(Ass.units{iArea_})+1);
                if ~isempty(Ass.units{iArea_})
                    for iAss=1:length(Ass.units{iArea_})
                        LocalMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_),Ass.units{iArea_}{iAss});
                    end
                end
                
                JointMembersMatrix_{iArea_} = nan(Ass.nu(iArea_),length(Ass.units{3})+1);
                for iAss=1:length(Ass.units{3})
                    if ~isempty(Ass.units{3})
                        units_ = Ass.units{3}{iAss};
                        if iArea_==1
                            JointMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_), units_(units_<=Ass.nu(1)));
                        else
                            JointMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_), units_(units_>Ass.nu(1))-Ass.nu(1));
                        end
                    end
                    
                end
            end
            iArea_ = 3;
            JointMembersMatrix_{iArea_} = nan(Ass.nu(iArea_),length(Ass.units{3})+1);
            for iAss=1:length(Ass.units{3})
                if ~isempty(Ass.units{3})
                    units_ = Ass.units{3}{iAss};
                    JointMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_), units_);
                end
            end
            clear units usel_out nu FSCsel Tmtx  iFR noPFC idx nAss U_ iAss
            nu = cellfun(@length,Membership_);nu(3) = sum(nu);
        end
        
        for iDelay = 1:length(Delays_)
            
            try
                eval(sprintf('D_.%s.DelaySpan.CVE{iArea}(:,iFile)      = D.%s.DelaySpan.CVE;',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_.%s.DelaySpan.CVE{iArea}(:,iFile)      = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
            end
            try
                eval(sprintf('D_.%s.DelaySpan_err.CVE{iArea}(:,iFile)      = D.%s.DelaySpan_err.CVE;',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.DelaySpan_err.PredErr{iArea}(:,iFile)  = D.%s.DelaySpan_err.PredErr;',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_.%s.DelaySpan_err.CVE{iArea}(:,iFile)      = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                eval(sprintf('D_.%s.DelaySpan_err.PredErr{iArea}(:,iFile)  = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
            end
            
            if strcmp(UnitSelection,'pairs')
                clear temp
                eval(sprintf('temp.CVE_Raw     = D.%s.DelaySpan.CVE_raw;',Delays_{iDelay}))
                try
                    eval(sprintf('temp.CVE_RawErr  = D.%s.DelaySpan_err.CVE_raw;',Delays_{iDelay}))
                    eval(sprintf('temp.CVE_PredErr = D.%s.DelaySpan_err.PredErr_raw;',Delays_{iDelay}))
                end

                if iArea<3
                    
                    % stash results specifically of NON assembly members
                    % ...get indices of nonmembers
                    lookup_ = [];
                    members_ = find(Membership_{iArea}==0);
                    for i = 1:length(members_)
                        for j = 1:length(members_)
                            if i~=j 
                                lookup_= [lookup_;[members_(i),members_(j)]];
                            end
                        end
                    end
                    if ~isempty(lookup_)
                        lookup_ = sub2ind(size(temp.CVE_Raw),lookup_(:,1),lookup_(:,2));
                        eval(sprintf('D_.%s.DelaySpan.CVEnonMembers{iArea}(:,iFile)= nanmean(cell2mat({temp.CVE_Raw{lookup_}}''),1);',Delays_{iDelay}))
                        try
                            eval(sprintf('D_.%s.DelaySpan.CVEErrnonMembers{iArea}(:,iFile)  = nanmean(cell2mat({temp.CVE_RawErr{lookup_}}''),1);',Delays_{iDelay}))
                            eval(sprintf('D_.%s.DelaySpan.PredErrnonMembers{iArea}(:,iFile) = nanmean(cell2mat({temp.CVE_PredErr{lookup_}}''),1);',Delays_{iDelay}))
                        end
                    else
                        eval(sprintf('D_.%s.DelaySpan.CVEnonMembers{iArea}(:,iFile)= nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.CVEErrnonMembers{iArea}(:,iFile)= nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.PredErrnonMembers{iArea}(:,iFile)= nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                    end
                    
                    % stash results specifically of LOCAL assembly members
                    % ...get indices of local members
                    lookup_ = [];
                    for iAss=1:size(LocalMembersMatrix_{iArea},2)-1
                       members_ = find(LocalMembersMatrix_{iArea}(:,iAss+1));
                       for i = 1:length(members_)
                           for j = 1:length(members_)
                               if i~=j
                                    lookup_= [lookup_;[members_(i),members_(j)]];
                               end
                           end
                       end
                    end
                    if ~isempty(lookup_)
                        lookup_ = sub2ind(size(temp.CVE_Raw),lookup_(:,1),lookup_(:,2));
                        eval(sprintf('D_.%s.DelaySpan.CVElocalMembers{iArea}(:,iFile)= nanmean(cell2mat({temp.CVE_Raw{lookup_}}''),1);',Delays_{iDelay}))
                        try
                            eval(sprintf('D_.%s.DelaySpan.CVEErrlocalMembers{iArea}(:,iFile)  = nanmean(cell2mat({temp.CVE_RawErr{lookup_}}''),1);',Delays_{iDelay}))
                            eval(sprintf('D_.%s.DelaySpan.PredErrlocalMembers{iArea}(:,iFile) = nanmean(cell2mat({temp.CVE_PredErr{lookup_}}''),1);',Delays_{iDelay}))
                        end
                    else
                        eval(sprintf('D_.%s.DelaySpan.CVElocalMembers{iArea}(:,iFile)= nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.CVEErrlocalMembers{iArea}(:,iFile)= nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.PredErrlocalMembers{iArea}(:,iFile)= nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                    end
                    
                    % stash results specifically of JOINT assembly members
                    % ...get indices of joint members
                    lookup_ = [];
                    for iAss=1:size(JointMembersMatrix_{iArea},2)-1
                       members_ = find(JointMembersMatrix_{iArea}(:,iAss+1));
                       for i = 1:length(members_)
                           for j = 1:length(members_)
                               if i~=j
                                    lookup_= [lookup_;[members_(i),members_(j)]];
                               end
                           end
                       end
                    end
                    if ~isempty(lookup_)
                        lookup_ = sub2ind(size(temp.CVE_Raw),lookup_(:,1),lookup_(:,2));
                        eval(sprintf('D_.%s.DelaySpan.CVEjointMembers{iArea}(:,iFile)= nanmean(cell2mat({temp.CVE_Raw{lookup_}}''),1);',Delays_{iDelay}))
                        try
                            eval(sprintf('D_.%s.DelaySpan.CVEErrjointMembers{iArea}(:,iFile) = nanmean(cell2mat({temp.CVE_RawErr{lookup_}}''),1);',Delays_{iDelay}))
                            eval(sprintf('D_.%s.DelaySpan.PredErrjointMembers{iArea}(:,iFile)= nanmean(cell2mat({temp.CVE_PredErr{lookup_}}''),1);',Delays_{iDelay}))
                        end
                    else
                        eval(sprintf('D_.%s.DelaySpan.CVEjointMembers{iArea}(:,iFile)     = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.CVEErrjointMembers{iArea}(:,iFile)  = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.PredErrjointMembers{iArea}(:,iFile) = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                    end
                
                else
                    % stash results specifically of NON assembly members
                    % ...get indices of nonmembers
                    lookup_ = [];
                    members_ = find([Membership_{1};Membership_{2}]==0); 
                    for i = 1:length(members_)
                        for j = 1:length(members_)
                            if i~=j && members_(i)<=nu(1) && members_(j)>nu(1)
                                lookup_= [lookup_;[members_(i),members_(j)]];
                            end
                        end
                    end
                    if ~isempty(lookup_)
                        lookup_ = sub2ind(size(temp.CVE_Raw),lookup_(:,1),lookup_(:,2));
                        eval(sprintf('D_.%s.DelaySpan.CVEnonMembers{iArea}(:,iFile)= nanmean(cell2mat({temp.CVE_Raw{lookup_}}''),1);',Delays_{iDelay}))
                        try
                            eval(sprintf('D_.%s.DelaySpan.CVEErrnonMembers{iArea}(:,iFile)  = nanmean(cell2mat({temp.CVE_RawErr{lookup_}}''),1);',Delays_{iDelay}))
                            eval(sprintf('D_.%s.DelaySpan.PredErrnonMembers{iArea}(:,iFile) = nanmean(cell2mat({temp.CVE_PredErr{lookup_}}''),1);',Delays_{iDelay}))
                        end
                    else
                        eval(sprintf('D_.%s.DelaySpan.CVEnonMembers{iArea}(:,iFile)= nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.CVEErrnonMembers{iArea}(:,iFile)= nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.PredErrnonMembers{iArea}(:,iFile)= nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                    end
                    
                    % stash results specifically of JOINT assembly members
                    % ...get indices of joint members
                    lookup_ = [];
                    for iAss=1:size(JointMembersMatrix_{iArea},2)-1
                       members_ = find(JointMembersMatrix_{iArea}(:,iAss+1));
                       for i = 1:length(members_)
                           for j = 1:length(members_)
                               if i~=j && members_(i)<=nu(1) && members_(j)>nu(1)
                                    lookup_= [lookup_;[members_(i),members_(j)]];
                               end
                           end
                       end
                    end
                    if ~isempty(lookup_)
                        lookup_ = sub2ind(size(temp.CVE_Raw),lookup_(:,1),lookup_(:,2));
                        eval(sprintf('D_.%s.DelaySpan.CVEjointMembers{iArea}(:,iFile)= nanmean(cell2mat({temp.CVE_Raw{lookup_}}''),1);',Delays_{iDelay}))
                        try
                            eval(sprintf('D_.%s.DelaySpan.CVEErrjointMembers{iArea}(:,iFile) = nanmean(cell2mat({temp.CVE_RawErr{lookup_}}''),1);',Delays_{iDelay}))
                            eval(sprintf('D_.%s.DelaySpan.PredErrjointMembers{iArea}(:,iFile)= nanmean(cell2mat({temp.CVE_PredErr{lookup_}}''),1);',Delays_{iDelay}))
                        end
                    else
                        eval(sprintf('D_.%s.DelaySpan.CVEjointMembers{iArea}(:,iFile)     = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.CVEErrjointMembers{iArea}(:,iFile)  = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        eval(sprintf('D_.%s.DelaySpan.PredErrjointMembers{iArea}(:,iFile) = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                    end
                    
                end
                
            end
            
            
            clear D_temp
            
            
        end
        
        clear D
    end
    
end
clear D Stability D_temp D_tempL D_tempR B bp
%% batch import Assem data for meta-analysis and plotting
for iFile = 1:length(fileListAss)
    try
        fname=strtok(fileListAss(iFile).name,'_');
        %         fnIn = sprintf('%s\\MixedSelectivity\\%s_DelayDecoding_Ass.mat',pat,fname);
        fnIn = sprintf('%s\\MixedSelectivity_LONG\\DelayDecoding\\%s_%s_DelayDecoding_Ass.mat',pat,fname,Areas{iArea});
        load(fnIn ,'D_Ass');
        
        for iArea = 1:3
            for iDelay = 1:length(Delays_)
                
                % L/R: Continuous decoding of Positional information - correct trials
                flag = 0;
                if eval(sprintf('isfield(D_Ass.%s,''DelaySpan'')',Delays_{iDelay}))
                    if eval(sprintf('~isempty(D_Ass.%s.DelaySpan{iArea})',Delays_{iDelay}))
                        flag = 1;
                        eval(sprintf('D_temp = D_Ass.%s.DelaySpan{iArea};',Delays_{iDelay}))
                        
                        try
                            eval(sprintf('D_Ass_.%s.DelaySpan.CVE{iArea}(:,iFile)      = D_temp.CVE;',Delays_{iDelay}))
                        catch
                            eval(sprintf('D_Ass_.%s.DelaySpan.CVE{iArea}(:,iFile)     = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        end
                        try
                            eval(sprintf('D_Ass_.%s.DelaySpan.CVE_err{iArea}(:,iFile)      = D_temp.CVE_err;',Delays_{iDelay}))
                        catch
                            eval(sprintf('D_Ass_.%s.DelaySpan.CVE_err{iArea}(:,iFile)     = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                        end
                        
                        clear D_temp
                    end
                end
                
            end
        end
    catch
    end
end

%% Strip non-responders
for iDelay = 1:length(Delays_)
    for iArea = 1:length(Areas)
        try
            eval(sprintf('D_temp = D_.%s.DelaySpan.CVE{iArea};',Delays_{iDelay}))
            D_temp(:,sum(D_temp)==0)=NaN;
            D_temp(:,isnan(D_temp(1,:)))=[];
            eval(sprintf('D_.%s.DelaySpan.CVE{iArea} = D_temp;',Delays_{iDelay}))
        end
    end
    for iArea = 1:length(Areas)
        try
            eval(sprintf('D_temp = D_.%s.DelaySpan_err.CVE{iArea};',Delays_{iDelay}))
            D_temp(:,sum(D_temp)==0)=NaN;
            D_temp(:,isnan(D_temp(1,:)))=[];
            eval(sprintf('D_.%s.DelaySpan_err.CVE{iArea} = D_temp;',Delays_{iDelay}))
            eval(sprintf('D_temp = D_.%s.DelaySpan_err.PredErr{iArea};',Delays_{iDelay}))
            D_temp(:,sum(D_temp)==0)=NaN;
            D_temp(:,isnan(D_temp(1,:)))=[];
            eval(sprintf('D_.%s.DelaySpan_err.PredErr{iArea} = D_temp;',Delays_{iDelay}))
        end
    end
    
    %     for iArea = 1:3
    %         eval(sprintf('D_temp = D_Ass_.%s.DelaySpan.CVE{iArea};',Delays_{iDelay}))
    %         D_temp(:,sum(D_temp)==0)=NaN;
    %         D_temp(:,isnan(D_temp(1,:)))=[];
    %         eval(sprintf('D_Ass_.%s.DelaySpan.CVE{iArea} = D_temp;',Delays_{iDelay}))
    %
    %         eval(sprintf('D_temp = D_Ass_.%s.DelaySpan.CVE_err{iArea};',Delays_{iDelay}))
    %         D_temp(:,sum(D_temp)==0)=NaN;
    %         D_temp(:,isnan(D_temp(1,:)))=[];
    %         eval(sprintf('D_Ass_.%s.DelaySpan.CVE_err{iArea} = D_temp;',Delays_{iDelay}))
    %     end
end

%% Plot L/R decoding on correct trials for different delays overlaid - CVE vs CVEerr - units
alpha_ = 0.3;
figure;
col_delay = flipud(hsv(3));
for iArea =1:length(Areas)
    subplot(1,2,iArea);hold on
    for iDelay = 1:length(Delays_)
        eval(sprintf('tb = tb%s;',Delays_{iDelay}))

        plot([0 0],[10 90],':k','LineWidth',1.5,'HandleVisibility','off');%'g'
        plot(str2double(strtok(Delays__{iDelay},'s'))*[1 1],[10 90],'color',[col_delay(iDelay,:), 0.5],'LineWidth',1.5,'HandleVisibility','off')
        plot([min(tb) max(tb)],[50 50],':k','HandleVisibility','off')
        
        
       
        eval(sprintf('corr = (1-D_.%s.DelaySpan.CVE{iArea})*100;',Delays_{iDelay}));
        eval(sprintf('err = (1-D_.%s.DelaySpan_err.PredErr{iArea})*100;',Delays_{iDelay}));
%                 corr = smooth2a(corr,10,0);
%                 err = smooth2a(err,10,0);
        
%         plot(tb,corr,'color',[0 0 0 0.1])
%         plot(tb,err,'color',[1 0 0 0.1])
        
        mC = nanmean(corr,2);
        eC = nansem(corr,2);
        mE = nanmean(err,2);
        eE = nansem(err,2);
        
        ciplot(mC+eC,mC-eC,tb,col_delay(iDelay,:),alpha_)
%         ciplot(mE+eE,mE-eE,tb,col_delay(iDelay,:),alpha_)
        

        title(Areas{iArea})
        axis([-5 20 0 100])
        if iArea==1
            ylabel('% Correct decoding')
        else
            legend(Delays__,'Orientation','horizontal'); legend boxoff
        end
        xlabel('Time (s)')
    end
end

%% Plot L/R decoding on correct and errors - CVE vs CVEerr - units
alpha_ = 0.6;
for iArea =1:length(Areas)
    figure('name',['L/R decoding ',Areas{iArea}]);
    for iDelay = 1:length(Delays_)
        
        eval(sprintf('tb = tb%s;',Delays_{iDelay}))
        
        subplot(length(Delays_),1,iDelay);hold on
        eval(sprintf('corr = (1-D_.%s.DelaySpan.CVE{iArea})*100;',Delays_{iDelay}));
        eval(sprintf('err = (1-D_.%s.DelaySpan_err.CVE{iArea})*100;',Delays_{iDelay}));
        %         corr = smooth2a(corr,10,0);
        %         err = smooth2a(err,10,0);
        
        plot(tb,corr,'color',[0 0 0 0.1])
        plot(tb,err,'color',[1 0 0 0.1])
        
        mC = nanmean(corr,2);
        eC = nansem(corr,2);
        mE = nanmean(err,2);
        eE = nansem(err,2);
        
        ciplot(mC+eC,mC-eC,tb,'k',alpha_)
        ciplot(mE+eE,mE-eE,tb,'r',alpha_)
        
        [sig,crit] = permtest2vec(corr,err,500,0.05);
        
        a = nan(size(sig));a(sig)=1;
        plot(tb,repmat(105,length(tb),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)
        plot(tb,a*105,'color',[0 0 0],'LineWidth',3)
        
        
        plot([0 0],[0 100],'g','LineWidth',1.5)
        plot(str2double(strtok(Delays__{iDelay},'s'))*[1 1],[0 100],'color',[0 0 0 0.1],'LineWidth',1.5)
        plot([min(tb) max(tb)],[50 50],':k')
        title(Delays_{iDelay})
        axis([-5 20 0 110])
        if iDelay==2
            ylabel('% Correct decoding')
        end
        xlabel('Time (s)')
    end
    legend({'Correct','Errors'})
end
%% Plot L/R decoding on correct and errors - CVE vs PredErr - units
alpha_ = 0.6;
for iArea =1:length(Areas)
    figure('name',['L/R decoding ',Areas{iArea}]);
    for iDelay = 1:length(Delays_)
        
        eval(sprintf('tb = tb%s;',Delays_{iDelay}))
        
        subplot(length(Delays_),1,iDelay);hold on
        eval(sprintf('corr = (1-D_.%s.DelaySpan.CVE{iArea})*100;',Delays_{iDelay}));
        eval(sprintf('err = (1-D_.%s.DelaySpan_err.PredErr{iArea})*100;',Delays_{iDelay}));
%                  corr = smooth2a(corr,10,0);
%                  err = smooth2a(err,10,0);
        
        plot(tb,corr,'color',[0 0 0 0.1])
        plot(tb,err,'color',[1 0 0 0.1])
        
        mC = nanmean(corr,2);
        eC = nansem(corr,2);
        mE = nanmean(err,2);
        eE = nansem(err,2);
        
        ciplot(mC+eC,mC-eC,tb,'k',alpha_)
        ciplot(mE+eE,mE-eE,tb,'r',alpha_)
        
        [sig,crit] = permtest2vec(corr,err,100,0.05);
        
        a = nan(size(sig));a(sig)=1;
        plot(tb,repmat(105,length(tb),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)
        plot(tb,a*105,'color',[0 0 0],'LineWidth',3)
        
        
        plot([0 0],[0 100],'g','LineWidth',1.5)
        plot(str2double(strtok(Delays__{iDelay},'s'))*[1 1],[0 100],'color',[0 0 0 0.1],'LineWidth',1.5)
        plot([min(tb) max(tb)],[50 50],':k')
        title(Delays_{iDelay})
        axis([-5 20 0 110])
        if iDelay==2
            ylabel('% Correct decoding')
        end
        xlabel('Time (s)')
    end
    legend({'Correct','Errors'})
end
%% Plot delay-dependent decoding: separate assembly classes

alpha_ = 0.6;
for iArea =1:length(Areas)
    figure('name',['L/R decoding ',Areas{iArea}]);
    for iDelay = 1:length(Delays_)
        
        eval(sprintf('tb = tb%s;',Delays_{iDelay}))
        
        subplot(length(Delays_),1,iDelay);hold on
        eval(sprintf('corrNon   = (1-D_.%s.DelaySpan.CVEnonMembers{iArea})*100;',Delays_{iDelay}));
        if iArea<3
            eval(sprintf('corrLocal = (1-D_.%s.DelaySpan.CVElocalMembers{iArea})*100;',Delays_{iDelay}));
        end
        eval(sprintf('corrJoint = (1-D_.%s.DelaySpan.CVEjointMembers{iArea})*100;',Delays_{iDelay}));
%         eval(sprintf('err = (1-D_.%s.DelaySpan_err.CVE{iArea})*100;',Delays_{iDelay}));
        %         corr = smooth2a(corr,10,0);
        %         err = smooth2a(err,10,0);
        
%         plot(tb,corr,'color',[0 0 0 0.3])
%         plot(tb,err,'color',[1 0 0 0.3])
        
        mC_non = nanmean(corrNon,2);
        eC_non = nansem(corrNon,2);
         if iArea<3
             mC_local = nanmean(corrLocal,2);
             eC_local = nansem(corrLocal,2);
         end
        mC_joint = nanmean(corrJoint,2);
        eC_joint = nansem(corrJoint,2);
        
%         mE = nanmean(err,2);
%         eE = nansem(err,2);
        
        ciplot(mC_non+eC_non,mC_non-eC_non,tb,col_{1},alpha_)
        if iArea<3
            ciplot(mC_local+eC_local,mC_local-eC_local,tb,col_{2},alpha_)
        end
        ciplot(mC_joint+eC_joint,mC_joint-eC_joint,tb,col_{3},alpha_)
        
%         [sig,crit] = permtest2vec(corr,err,500,0.05);       
%         a = nan(size(sig));a(sig)=1;
%         plot(tb,repmat(105,length(tb),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)
%         plot(tb,a*105,'color',[0 0 0],'LineWidth',3)
        
        
        plot([0 0],[0 100],'g','LineWidth',1.5)
        plot(str2double(strtok(Delays__{iDelay},'s'))*[1 1],[0 100],'color',[0 0 0 0.1],'LineWidth',1.5)
        plot([min(tb) max(tb)],[50 50],':k')
        title(Delays_{iDelay})
        axis([-5 20 0 110])
        if iDelay==2
            ylabel('% Correct decoding')
        end
        xlabel('Time (s)')
    end
    legend({'Correct','Errors'})
end

%% Plot L/R decoding on correct and errors - CVE - Assemblies
alpha_ = 0.6;
for iArea =1:length(Areas)
    figure('name',['L/R decoding ',Areas{iArea}]);
    for iDelay = 1:length(Delays_)
        eval(sprintf('tb = tb%s;',Delays_{iDelay}))
        
        subplot(length(Delays_),1,iDelay);hold on
        eval(sprintf('corr = (1-D_Ass_.%s.DelaySpan.CVE{iArea})*100;',Delays_{iDelay}));
        eval(sprintf('err = (1-D_Ass_.%s.DelaySpan.CVE_err{iArea})*100;',Delays_{iDelay}));
        %         corr = smooth2a(corr,10,0);
        %         err = smooth2a(err,10,0);
        
        m = nanmean(corr,2);
        e = nansem(corr,2);
        %         plot(tb,corr,'k')
        ciplot(m+e,m-e,tb,'k',alpha_)
        m = nanmean(err,2);
        e = nansem(err,2);
        %         plot(tb,err,'r')
        ciplot(m+e,m-e,tb,'r',alpha_)
        
        [sig,crit] = permtest2vec(corr,err,500,0.05);
        
        a = nan(size(sig));a(sig)=1;
        plot(tb,repmat(105,length(tb),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)
        plot(tb,a*105,'color',[0 0 0],'LineWidth',3)
        
        
        plot([0 0],[0 100],'g','LineWidth',1.5)
        plot(str2double(strtok(Delays__{iDelay},'s'))*[1 1],[0 100],'color',[0 0 0 0.1],'LineWidth',1.5)
        plot([min(tb) max(tb)],[50 50],':k')
        title(Delays_{iDelay})
        axis([-5 20 0 110])
        if iDelay==2
            ylabel('% Correct decoding')
        end
        xlabel('Time (s)')
    end
    legend({'Correct','Errors'})
end
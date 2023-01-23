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
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];    
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);

fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat'}; %'ALL_events.mat'
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
%% Batch process units
for iFile =1:length(fileList)
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    
    for iArea = 1:2%length(Areas)
        %% Get the files
        fprintf('Analysing run %d/%d %s (%s)...',iFile,length(fileList),fname,Areas{iArea})
        % load(sprintf('%sKDE_bins%s%s_%s_iFR50.mat',pat,filesep,fname,Areas{iArea}));
        load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
        
        iFR_ = iFR;
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
                
            % Cross-validation error decoder (correct trials)    
%             if nL>2 && nR>2
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = length(tb_);
                nu = [length(PFCcells),length(HPcells)];
%                 
%                  if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
%                  
%                  for i=1:nrep
%                      rs{i} = randsample(size(Ltrials,2),min(nu));
%                      D_.CVE(i,:) = DecodeCVE_mex(FR(:,rs{i}),evt0,0.05);
%                  end
%                  D_.CVE = nanmean(D_.CVE,1);
%                 

                 CVE = cell(nu(iArea),nu(iArea));
                 for i =1:nu(iArea)
                    parfor j =1:nu(iArea)
                        [i j]
                        if i~=j
                            CVE{i,j} = DecodeCVE_mex(FR(:,[i,j]),evt0,0.05);
                        end
                    end
                 end
                 D_.CVE_raw = CVE;
                 D_.CVE = nanmean(cell2mat(CVE(:)),1);
                 D_.CVE_SEM = nansem(cell2mat(CVE(:)),1);
                 clear CVE
                 %hold on; plot(1-D_.CVE+D_.CVE_SEM);plot(1-D_.CVE-D_.CVE_SEM);
                 clear rs nu evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr

%             end
            if exist('D_')
                eval(sprintf('D.%s.DelaySpan =  D_;',Delays_{iDelay}))
            end
            clear D_
        
             % Cross-validation error decoder (error trials)    
%              if nLe>2 && nRe>2
                FR = [Ltrials;Rtrials;LtrialsE;RtrialsE];
                evt0 = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
                Ltr = length(tb_);
                nu = [length(PFCcells),length(HPcells)];

%                  if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
%                  
%                  for i=1:nrep
%                      rs{i} = randsample(size(Ltrials,2),min(nu));
%                      D_.CVE(i,:) = DecodePredErr(FR(:,rs{i}),evt0,0.05);
%                  end
%                  D_.CVE = nanmean(D_.CVE,1);
                 
                 CVE = cell(nu(iArea),nu(iArea));
                 for i =1:nu(iArea)
                    parfor j =1:nu(iArea)
                        [i j]
                        if i~=j
                            CVE{i,j} = DecodeCVE_mex(FR(:,[i,j]),evt0,0.05);
                        end
                    end
                 end
                 D_.CVE_raw = CVE;
                 D_.CVE = nanmean(cell2mat(CVE(:)),1);
                 D_.CVE_SEM = nansem(cell2mat(CVE(:)),1);
                 clear CVE
                 clear rs nu evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
                 
%              end
             if exist('D_')
                 eval(sprintf('D.%s.DelaySpan_err =  D_;',Delays_{iDelay}))
             end
             clear D_
             clear Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
        end
        %% Save results
        if exist('D')
            
            fnOut = sprintf('%sMixedSelectivity%s%s_%s_DelayDecoding_Units.mat',pat,filesep,fname,Areas{iArea});
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


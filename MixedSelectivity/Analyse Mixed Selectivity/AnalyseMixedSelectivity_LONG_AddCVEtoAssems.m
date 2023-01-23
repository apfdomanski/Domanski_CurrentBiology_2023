%% %%%%%% PREAMBLE %%%%%%
clear 
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Target = 'LONG';



tlimsAll = [-5 5];
tlimsShort=[-4 10];
tlimsMedium=[-8 10];
tlimsLong=[-16 10];
tlimsANOVA = [-2 2];

shift = 0;
plotOnline = false;
runCVEdecoder = true;
bw=0.05;
Nbs = 100;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbANOVA=tlimsANOVA(1):bw:tlimsANOVA(2);
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
%% Batch Process assemblies 

for iFile =4:length(fileListAss)
    %% Get the files
    fname=strtok(fileListAss(iFile).name,'_');
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    fprintf('Analysing run %d/%d %s ...',iFile,length(fileList),fname)
    switch AssemblyChoice
        case 1
            Ass = load(sprintf('%s%s_iFR50_FSC.mat',pat2,fname));    
        case 2
            Ass = load(sprintf('%s%s_iFR50_BehavOnly_FSC.mat',pat2,fname));    
        case 3
            Ass = load(sprintf('%s%s_iFR50_Task_FSC.mat',pat2,fname));
    end
     fnOut = sprintf('%sMixedSelectivity_LONG%sLeftvsRight%s%s_MixedSelectivity_Ass.mat',pat,filesep,filesep,fname);
    load(fnOut,'D_Ass','tbShort','tbMedium','tbLong');
    %% Decode Left/Right on correct trials
    for iDelay = 1:length(Delays_)
        eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect];',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf('RightTrials = [t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}));
        for iArea = 1:length(Ass.FSCsel)
            
            FSC = Ass.FSCsel{iArea};
            
            if ~isempty(FSC)
                
                Ltrials = [];nL = 0; Ltrials_ = {};

                for iTrial =1:size(LeftTrials,1)
                    try
                        tlims_  = [LeftTrials(iTrial,1)/1e6;LeftTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                        Ltrials = [Ltrials;FSC(tlims_,:)];
                        Ltrials_{1,nL+1} = FSC(tlims_,:);
                        nL=nL+1;
                    end
                end
                
                Rtrials = [];nR = 0; Rtrials_ = {};
                            
                for iTrial =1:size(RightTrials,1)
                    try
                        tlims_  = [RightTrials(iTrial,1)/1e6;RightTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                        Rtrials = [Rtrials;FSC(tlims_,:)];
                        Rtrials_{1,nR+1} = FSC(tlims_,:);
                        nR=nR+1;
                    end
                end
                
                if nL>5 && nR>5
                   
                    FSC_ = [Ltrials;Rtrials];
                    evt0 = [ones(nL,1);2*ones(nR,1)];
                    Ltr = 2*length(tbAll);
                    nAss_ = cellfun(@(x)size(x,2),Ass.FSCsel);
                    % remove zero entries
                    idx = nAss_>0;
                    nAss_(~idx) = min(nAss_(idx));
                    
                     
                     % Cross-validation error decoder
                     FSC_ = [Ltrials;Rtrials];
                     evt0 = [ones(nL,1);2*ones(nR,1)];
                     for iAss =1:size(FSC,2)
                         iAss
%                          D_.LR{iArea}.CVEindividual(iAss,:) = DecodeCVE_mex(FSC_(:,iAss),evt0,0.05);
                         D_.LR{iArea}.CVEindividual(iAss,:) = DecodeCVE(FSC_(:,iAss),evt0,0.05);
                     end
                     if runCVEdecoder
                         if size(Ltrials,2)==min(nAss_), nrep=1; else nrep = 10; end
                         % nrep=1;
                         for i=1:nrep
                             rs{i} = randsample(size(Ltrials,2),min(nAss_));
                             % rs{i} = randsample(size(Ltrials,2),nu(iArea));
                             %                         D_.LR{iArea}.CVE(i,:) = DecodeCVE(FSC_,evt0,0.05);
%                              D_.LR{iArea}.CVE(i,:) = DecodeCVE_mex(FSC_(:,rs{i}),evt0,0.05);
                             D_.LR{iArea}.CVE(i,:) = DecodeCVE(FSC_(:,rs{i}),evt0,0.05);
                         end
                         D_.LR{iArea}.CVE = nanmean(D_.LR{iArea}.CVE,1);
                         
                         D_.LR{iArea}.cveBSciH=zeros(1,Ltr);  D_.LR{iArea}.cveBSciL=zeros(1,Ltr);
                         D_.LR{iArea}.CVEbs=[];
                         CVEbs = zeros(Nbs,Ltr);
                         
                         parfor b=1:Nbs
                             t1=tic;
                             k=randperm(length(evt0)); evt1=evt0(k);
                             cve0=zeros(nrep,Ltr);
                             %for i=1:nrep, cve0(i,:)=DecodeCVE(FSC_,evt1,0.05); end
%                              for i=1:nrep, cve0(i,:)=DecodeCVE_mex(FSC_(:,rs{i}),evt1,0.05); end
                             for i=1:nrep, cve0(i,:)=DecodeCVE(FSC_(:,rs{i}),evt1,0.05); end
                             CVEbs(b,:) = nanmean(cve0,1);
                             fprintf('%s: Decoded %s assemblies, correct trials, %s delay, draw %d of %d (%s elapsed)...\n',fname, Areas{iArea},Delays_{iDelay},b,Nbs,seconds2human(toc(t1)))
                         end
                         for tr=1:Ltr
                             cves=sort(CVEbs(:,tr),'ascend');
                             D_.LR{iArea}.cveBSciH(tr)=cves(round(0.95*Nbs));
                             D_.LR{iArea}.cveBSciL(tr)=cves(round(0.05*Nbs));
                         end
                         D_.LR{iArea}.CVEbs = nanmean(CVEbs,1);
                     end
                     clear rs nAss idx evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
                end
            end
        end
        
        if exist('D_')
            for iArea=1:3
                if(~isempty(D_.LR{iArea}))
                eval(sprintf('D_Ass.%s.LR{iArea}.CVEindividual =  D_.LR{iArea}.CVEindividual;',Delays_{iDelay}))
                eval(sprintf('D_Ass.%s.LR{iArea}.CVE           =  D_.LR{iArea}.CVE;',Delays_{iDelay}))
                eval(sprintf('D_Ass.%s.LR{iArea}.cveBSciH      =  D_.LR{iArea}.cveBSciH;',Delays_{iDelay}))
                eval(sprintf('D_Ass.%s.LR{iArea}.cveBSciL      =  D_.LR{iArea}.cveBSciL;',Delays_{iDelay}))
                eval(sprintf('D_Ass.%s.LR{iArea}.CVEbs         =  D_.LR{iArea}.CVEbs;',Delays_{iDelay}))
                end
            end
        end
        clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
    end
    %% Decode Left/Right on error trials
    for iDelay = 1:length(Delays_)
        eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftError'',t.%s.ChoicePress_LeftError''];',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf('RightTrials = [t.%s.SamplePress_RightError'',t.%s.ChoicePress_RightError''];',Delays_{iDelay},Delays_{iDelay}));
        for iArea = 1:length(Ass.FSCsel)
            
            FSC = Ass.FSCsel{iArea};
            if ~isempty(FSC)
                
                Ltrials = [];nL = 0;

                for iTrial =1:size(LeftTrials,1)
                    try
                        tlims_  = [LeftTrials(iTrial,1)/1e6;LeftTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];

                        Ltrials = [Ltrials;FSC(tlims_,:)];
                        nL=nL+1;
                    end
                end
                
                Rtrials = [];nR = 0;
                            
                for iTrial =1:size(RightTrials,1)
                    try
                        tlims_  = [RightTrials(iTrial,1)/1e6;RightTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];

                        Rtrials = [Rtrials;FSC(tlims_,:)];
                        nR=nR+1;
                    end
                end
                
                if nL>2 && nR>2
                    
                    FSC_ = [Ltrials;Rtrials];
                    evt0 = [ones(nL,1);2*ones(nR,1)];
                    Ltr = 2*length(tbAll);
                    nAss_ = cellfun(@(x)size(x,2),Ass.FSCsel);
                    % remove zero entries
                    idx = nAss_>0;
                    nAss_(~idx) = min(nAss_(idx));
                    
                                     
                    % Cross-validation error decoder
                    FSC_ = [Ltrials;Rtrials];
                    evt0 = [ones(nL,1);2*ones(nR,1)];
                    for iAss =1:size(FSC,2)
                        iAss
                        %                          D_.LR{iArea}.CVEindividual(iAss,:) = DecodeCVE_mex(FSC_(:,iAss),evt0,0.05);
                        D_.LR{iArea}.CVEindividual(iAss,:) = DecodeCVE(FSC_(:,iAss),evt0,0.05);
                    end
                    if runCVEdecoder
                        if size(Ltrials,2)==min(nAss_), nrep=1; else nrep=10; end
                        %                  nrep=1;
                        for i=1:nrep
                            rs{i} = randsample(size(Ltrials,2),min(nAss_));
                            %                         D_.LR{iArea}.CVE(i,:) = DecodeCVE(FSC_,evt0,0.05);
%                             D_.LR{iArea}.CVE(i,:) = DecodeCVE_mex(FSC_(:,rs{i}),evt0,0.05);
                            D_.LR{iArea}.CVE(i,:) = DecodeCVE(FSC_(:,rs{i}),evt0,0.05);
                        end
                        D_.LR{iArea}.CVE = nanmean(D_.LR{iArea}.CVE,1);
                        
                        D_.LR{iArea}.cveBSciH=zeros(1,Ltr);  D_.LR{iArea}.cveBSciL=zeros(1,Ltr);
                        D_.LR{iArea}.CVEbs=[];
                        CVEbs = zeros(Nbs,Ltr);
                        
                        parfor b=1:Nbs
                            t1=tic;
                            k=randperm(length(evt0)); evt1=evt0(k);
                            cve0=zeros(nrep,Ltr);
                            %                         for i=1:nrep, cve0(i,:)=DecodeCVE(FSC_,evt1,0.05); end
%                             for i=1:nrep, cve0(i,:)=DecodeCVE_mex(FSC_,evt1,0.05); end
                            for i=1:nrep, cve0(i,:)=DecodeCVE(FSC_,evt1,0.05); end
                            CVEbs(b,:) = nanmean(cve0,1);
                            fprintf('%s: Decoded %s assemblies, error trials, %s delay, draw %d of %d (%s elapsed)...\n',fname,Areas{iArea},Delays_{iDelay},b,Nbs,seconds2human(toc(t1)))
                        end
                        for tr=1:Ltr
                            cves=sort(CVEbs(:,tr),'ascend');
                            D_.LR{iArea}.cveBSciH(tr)=cves(round(0.95*Nbs));
                            D_.LR{iArea}.cveBSciL(tr)=cves(round(0.05*Nbs));
                        end
                        D_.LR{iArea}.CVEbs = nanmean(CVEbs,1);
                    end
                 clear rs nAss idx evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
                 
                end
            end
        end
        if exist('D_')
%             eval(sprintf('D_Ass.%s.LR_err =  D_.LR;',Delays_{iDelay}))
            for iArea=1:3
                if(~isempty(D_.LR{iArea}))
                eval(sprintf('D_Ass.%s.LR_err{iArea}.CVEindividual =  D_.LR{iArea}.CVEindividual;',Delays_{iDelay}))
                eval(sprintf('D_Ass.%s.LR_err{iArea}.CVE           =  D_.LR{iArea}.CVE;',Delays_{iDelay}))
                eval(sprintf('D_Ass.%s.LR_err{iArea}.cveBSciH      =  D_.LR{iArea}.cveBSciH;',Delays_{iDelay}))
                eval(sprintf('D_Ass.%s.LR_err{iArea}.cveBSciL      =  D_.LR{iArea}.cveBSciL;',Delays_{iDelay}))
                eval(sprintf('D_Ass.%s.LR_err{iArea}.CVEbs         =  D_.LR{iArea}.CVEbs;',Delays_{iDelay}))
                end
            end
        end
        clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
    end
    %% Save
    fnOut = sprintf('%sMixedSelectivity_LONG%sLeftvsRight%s%s_MixedSelectivity_Ass.mat',pat,filesep,filesep,fname);
    save(fnOut,'D_Ass','tbShort','tbMedium','tbLong','-v7.3');
    fprintf([fname, ': Done.\n'])
    clear D_Ass
end
clearvars -except Nbs AssemblyChoice color_ Delays_ bw tbShort tbAll tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin 

%% %%%%%% PREAMBLE %%%%%%
clear 
Delays_ = {'Delay_0'};
Delays__ = {'0s'};
Target = 'SHORT';



tlimsAll = [-5 5];
tlimsDelay_0=[-2 2];
tlimsShort=[-4 10];

tlimsANOVA = [-2 2];

shift = 0;
plotOnline = false;
runCVEdecoder = false;
bw=0.05;
Nbs = 100;

tbDelay_0=tlimsDelay_0(1):bw:tlimsDelay_0(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
tbANOVA = tlimsANOVA(1):bw:tlimsANOVA(2);%0:bw:2*sum(abs(tlimsAll));

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

reject_list={''}; %'ALL_events.mat'
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
        pat2 = [pat 'KDE_binsTaskonly' filesep 'SHORTTaskonly' filesep];
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end

mkdir([pat 'MixedSelectivity_SHORT'])

%% Batch process units
for iFile =1:length(fileList)
    fname=strtok(fileList(iFile).name,'_');
    try
        
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    
    for iArea = 1:2%length(Areas)
        %% Get the files
        fprintf('Analysing run %d/%d %s (%s)...',iFile,length(fileList),fname,Areas{iArea})
%         load(sprintf('%sKDE_bins%s%s_%s_iFR50.mat',pat,filesep,fname,Areas{iArea}));
        load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
        
        iFR_ = iFR;
        % iFR_ = zscore(iFR_);
        %% Decode Left/Right on correct trials
        for iDelay =1:length(Delays_)
            
            eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect];',Delays_{iDelay},Delays_{iDelay}));
            eval(sprintf('RightTrials = [t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}));
            
            Ltrials = [];nL = 0;Ltrials_ = {};
            
            for iTrial =1:size(LeftTrials,1)
                try
                    tlims_  = [LeftTrials(iTrial,1)/1e6;LeftTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                    Ltrials = [Ltrials;iFR_(tlims_,:)];
                    Ltrials_{1,nL+1} = iFR_(tlims_,:);
                    nL=nL+1;
                end
            end
            Rtrials = [];nR = 0; Rtrials_ = {};
            
            for iTrial =1:size(RightTrials,1)
                try
                    tlims_  = [RightTrials(iTrial,1)/1e6;RightTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                    Rtrials = [Rtrials;iFR_(tlims_,:)];
                    Rtrials_{1,nR+1} = iFR_(tlims_,:);
                    nR=nR+1;
                end
            end
            
            % Run the decoders
            if nL>5 && nR>5
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = 2*length(tbAll);
                nu = [length(PFCcells),length(HPcells)];
                
                % Multivariate F-score decoder
                [D_.avgFR,D_.seFR,...
                 D_.Ft2,D_.Rt2,...
                 D_.Ft2ciL,D_.Ft2ciH,...
                 D_.TS,...
                 D_.dfnum,D_.dfden] = DecodeStats(FR,evt0,0.05);
                 D_.TSsig = tpdf(D_.TS,nL+nR-2)<0.05;
                 
                 % make trial counts consistent for comparisons
                 clear E
                 for nrep = 1:100
                     nrep
                       rsL = randsample(length(Ltrials_),5);
                       rsR = randsample(length(Rtrials_),5);
                       
                       FR   = [cell2mat(Ltrials_(rsL)');cell2mat(Rtrials_(rsR)')];
                       evt0 = [ones(5,1);2*ones(5,1)];
                       
                       % Multivariate F-score decoder
                       [E.avgFR(:,:,nrep),E.seFR(:,:,nrep),...
                        E.Ft2(:,nrep),E.Rt2(:,nrep),...
                        E.Ft2ciL(:,nrep),E.Ft2ciH(:,nrep),...
                        E.TS(:,:,nrep),...
                        ~,~] = DecodeStats(FR,evt0,0.05);
                        E.TSsig(:,:,nrep) = tpdf(E.TS(:,:,nrep),8)<0.05;                
                 end
                 D_.TS_drawnTrials =  mean(E.TS,3);
                 D_.TSsig_drawnTrials =  sum(E.TSsig,3)./100;
                 
                 
                 % Cross-validation error decoder
                 for iUnit =1:size(FR,2)
                     iUnit
                      D_.CVEindividual(iUnit,:) = DecodeCVE_mex(FR(:,iUnit),evt0,0.05);
                 end
                 if runCVEdecoder
                 if size(Ltrials,2)==min(nu), nrep=1; else nrep = 10; end;
%                  nrep=1;
                 for i=1:nrep
                    rs{i} = randsample(size(Ltrials,2),min(nu));
                    % rs{i} = randsample(size(Ltrials,2),nu(iArea));
                    % D_.CVE(i,:) = DecodeCVE(FR(:,rs{i}),evt0,0.05);
                    D_.CVE(i,:) = DecodeCVE_mex(FR(:,rs{i}),evt0,0.05);
                 end
                 D_.CVE = nanmean(D_.CVE,1);
                 
                 D_.cveBSciH=zeros(1,Ltr);  D_.cveBSciL=zeros(1,Ltr);
                 D_.CVEbs=[];
                 CVEbs = zeros(Nbs,Ltr);
                 parfor b=1:Nbs
                     t1=tic;
                     k=randperm(length(evt0)); evt1=evt0(k);
                     cve0=zeros(nrep,Ltr);
%                      for i=1:nrep, cve0(i,:)=DecodeCVE(FR,evt1,0.05); end
                     for i=1:nrep, cve0(i,:)=DecodeCVE_mex(FR(:,rs{i}),evt1,0.05); end
                     CVEbs(b,:) = nanmean(cve0,1);
                     fprintf('Decoded %s units, correct trials, %s delay, draw %d of %d (%s elapsed)...\n',Areas{iArea},Delays_{iDelay},b,Nbs,seconds2human(toc(t1)))
                 end
                 for tr=1:Ltr
                      cves=sort(CVEbs(:,tr),'ascend');
                      D_.cveBSciH(tr)=cves(round(0.95*Nbs));
                      D_.cveBSciL(tr)=cves(round(0.05*Nbs));
                 end
                 D_.CVEbs = nanmean(CVEbs,1);
                 end
                 clear rs nu evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
                 % figure; hold on
                 %  plot(D_.CVEbs);plot(D_.cveBSciL);plot(D_.cveBSciH);
                 %  plot(D_.CVE)                 
                  
                  eval(sprintf('D.%s.LR = D_;',Delays_{iDelay}));
            end
            
            clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
        end
        %% Decode Left/Right on error trials
        for iDelay =1:length(Delays_)
            
            eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftError'',t.%s.ChoicePress_LeftError''];',Delays_{iDelay},Delays_{iDelay}));
            eval(sprintf('RightTrials = [t.%s.SamplePress_RightError'',t.%s.ChoicePress_RightError''];',Delays_{iDelay},Delays_{iDelay}));
            
            Ltrials = [];nL = 0;            
            for iTrial =1:size(LeftTrials,1)
                 try
                    tlims_  = [LeftTrials(iTrial,1)/1e6;LeftTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                    Ltrials = [Ltrials;iFR_(tlims_,:)];
                    nL=nL+1;
                 end
            end
            
            Rtrials = [];nR = 0;            
            for iTrial =1:size(RightTrials,1)
                try
                    tlims_  = [RightTrials(iTrial,1)/1e6;RightTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                    Rtrials = [Rtrials;iFR_(tlims_,:)];
                    nR=nR+1;
                end
            end
            
            if nL>2 && nR>2
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = 2*length(tbAll);
                nu = [length(PFCcells),length(HPcells)];
                
                % Multivariate F-score decoder
                [D_.avgFR,D_.seFR,...
                 D_.Ft2,D_.Rt2,...
                 D_.Ft2ciL,D_.Ft2ciH,...
                 D_.TS,...
                 D_.dfnum,D_.dfden] = DecodeStats(FR,evt0,0.05);
                 D_.TSsig = tpdf(D_.TS,nL+nR-2)<0.05;
                 
                 % Cross-validation error decoder
                 if runCVEdecoder
                 if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
%                  nrep=1;
                 for i=1:nrep
                    rs{i} = randsample(size(Ltrials,2),min(nu));
%                     D_.CVE(i,:) = DecodeCVE(FR,[ones(nL,1);2*ones(nR,1)],0.05);
                    D_.CVE(i,:) = DecodeCVE_mex(FR(:,rs{i}),evt0,0.05);
                 end
                 D_.CVE(i,:) = nanmean(D_.CVE,1);
                 
                 D_.cveBSciH=zeros(1,Ltr);  D_.cveBSciL=zeros(1,Ltr);
                 D_.CVEbs=[];
                 CVEbs = zeros(Nbs,Ltr);
                 parfor b=1:Nbs
                     t1= tic;
                     k=randperm(length(evt0)); evt1=evt0(k);
                     cve0=zeros(nrep,Ltr);
%                      for i=1:nrep, cve0(i,:)=DecodeCVE(FR,evt1,0.05); end;
                     for i=1:nrep, cve0(i,:)=DecodeCVE_mex(FR(:,rs{i}),evt1,0.05); end;
                     CVEbs(b,:) = nanmean(cve0,1);
                     fprintf('Decoded %s units, error trials, %s delay, draw %d of %d (%s elapsed)...\n',Areas{iArea},Delays_{iDelay},b,Nbs,seconds2human(toc(t1)))
                 end
                 for tr=1:Ltr

                      cves=sort(CVEbs(:,tr),'ascend');
                      D_.cveBSciH(tr)=cves(round(0.95*Nbs));
                      D_.cveBSciL(tr)=cves(round(0.05*Nbs));
                 end
                 D_.CVEbs = nanmean(CVEbs,1);
                 end
                 clear rs nu evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
%                  figure; hold on
%                  plot(D_.CVEbs);plot(D_.cveBSciL);plot(D_.cveBSciH);
%                  plot(D_.CVE)
                 
                 eval(sprintf('D.%s.LR_err= D_;',Delays_{iDelay}));
            end
            clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
        end
        %% Decode Correct/Error
        for iDelay =1:length(Delays_)
            
            % Run once for L trials
            eval(sprintf('CorrectTrials = [t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect];',Delays_{iDelay},Delays_{iDelay}));
            eval(sprintf('ErrorTrials =   [t.%s.SamplePress_LeftError'',t.%s.ChoicePress_LeftError''];',Delays_{iDelay},Delays_{iDelay}));
            
            
            Ctrials = [];nC = 0;
            for iTrial =1:size(CorrectTrials,1)
                try
                    tlims_  = [CorrectTrials(iTrial,1)/1e6;CorrectTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                    Ctrials = [Ctrials;iFR_(tlims_,:)];
                    nC=nC+1;
                end
            end
            Etrials = [];nE = 0;
            for iTrial =1:size(ErrorTrials,1)
                try
                    tlims_  = [ErrorTrials(iTrial,1)/1e6;ErrorTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                    Etrials = [Etrials;iFR_(tlims_,:)];
                    nE=nE+1;
                end
            end
            if nC>1 && nE>1
                [D_.avgFR,D_.seFR,...
                    D_.Ft2,D_.Rt2,...
                    D_.Ft2ciL,D_.Ft2ciH,...
                    D_.TS,...
                    D_.dfnum,D_.dfden] = DecodeStats([Ctrials;Etrials],[ones(nC,1);2*ones(nE,1)],0.05);
                D_.TSsig = tpdf(D_.TS,nC+nE-2)<0.05;
                
             eval(sprintf('D.%s.CE_L= D_;',Delays_{iDelay}));
            end
            clear D_ Ctrials nC Etrials nE CorrectTrials ErrorTrials iTrial
 
            % Run once for R trials
            eval(sprintf('CorrectTrials = [t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}));
            eval(sprintf('ErrorTrials =   [t.%s.SamplePress_RightError'',t.%s.ChoicePress_RightError''];',Delays_{iDelay},Delays_{iDelay}));
            
            
            Ctrials = [];nC = 0;
            for iTrial =1:size(CorrectTrials,1)
                try
                    tlims_  = [CorrectTrials(iTrial,1)/1e6;CorrectTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                    Ctrials = [Ctrials;iFR_(tlims_,:)];
                    nC=nC+1;
                end
            end
            Etrials = [];nE = 0;
            for iTrial =1:size(ErrorTrials,1)
                try
                    tlims_  = [ErrorTrials(iTrial,1)/1e6;ErrorTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                    Etrials = [Etrials;iFR_(tlims_,:)];
                    nE=nE+1;
                end
            end
            
            if nC>1 && nE>1
                [D_.avgFR,D_.seFR,...
                    D_.Ft2,D_.Rt2,...
                    D_.Ft2ciL,D_.Ft2ciH,...
                    D_.TS,...
                    D_.dfnum,D_.dfden] = DecodeStats([Ctrials;Etrials],[ones(nC,1);2*ones(nE,1)],0.05);
                D_.TSsig = tpdf(D_.TS,nC+nE-2)<0.05;
                eval(sprintf('D.%s.CE_R= D_;',Delays_{iDelay}));
            end
            clear D_ Ctrials nC Etrials nE CorrectTrials ErrorTrials iTrial
         
        end        
        %% Decode Sample/Choice on correct trials
        for iDelay =1:length(Delays_)
            % Run once for L trials
            eval(sprintf('SampleTrials = t.%s.SamplePress_LeftCorrect;',Delays_{iDelay}));
            eval(sprintf('ChoiceTrials = t.%s.ChoicePress_LeftCorrect;',Delays_{iDelay}));
           
            Strials = [];nS = 0;
            for iTrial =1:length(SampleTrials)
                try
                    eval(sprintf('tlims_ = SampleTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}));
                    tlims_ = closest(Tmtx,tlims_);
                    eval(sprintf('Strials = [Strials;iFR_(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}));
                    nS=nS+1;
                end
            end
            
            Ctrials = []; nC=0;
            for iTrial =1:length(ChoiceTrials)
                try
                    eval(sprintf('tlims_ = ChoiceTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}));
                    tlims_ = closest(Tmtx,tlims_);
                    eval(sprintf('Ctrials = [Ctrials;iFR_(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}));
                    nC=nC+1;
                end
            end
            if nS>1 && nC>1
                [D_.avgFR,D_.seFR,...
                    D_.Ft2,D_.Rt2,...
                    D_.Ft2ciL,D_.Ft2ciH,...
                    D_.TS,...
                    D_.dfnum,D_.dfden] = DecodeStats([Strials;Ctrials],[ones(nS,1);2*ones(nC,1)],0.05);
                D_.TSsig = tpdf(D_.TS,nS+nC-2)<0.05;
                eval(sprintf('D.%s.SC_L= D_;',Delays_{iDelay}));
            end
             clear D_ Strials nS Ctrials nC SampleTrials ChoiceTrials iTrial
    
            % Run once for R trials
            eval(sprintf('SampleTrials = t.%s.SamplePress_RightCorrect;',Delays_{iDelay}));
            eval(sprintf('ChoiceTrials = t.%s.ChoicePress_RightCorrect;',Delays_{iDelay}));
            %eval(sprintf(';',Delays_{iDelay},Delays_{iDelay}));
           
            Strials = [];nS = 0;
            for iTrial =1:length(SampleTrials)
                try
                    eval(sprintf('tlims_ = SampleTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}));
                    tlims_ = closest(Tmtx,tlims_);
                    eval(sprintf('Strials = [Strials;iFR_(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}));
                    nS=nS+1;
                end
            end
            
            Ctrials = []; nC=0;
            for iTrial =1:length(ChoiceTrials)
                try
                    eval(sprintf('tlims_ = ChoiceTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}));
                    tlims_ = closest(Tmtx,tlims_);
                    eval(sprintf('Ctrials = [Ctrials;iFR_(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}));
                    nC=nC+1;
                end
            end
            if nS>1 && nC>1
            [D_.avgFR,D_.seFR,...
             D_.Ft2,D_.Rt2,...
             D_.Ft2ciL,D_.Ft2ciH,...
             D_.TS,...
             D_.dfnum,D_.dfden] = DecodeStats([Strials;Ctrials],[ones(nS,1);2*ones(nC,1)],0.05);
             D_.TSsig = tpdf(D_.TS,nS+nC-2)<0.05;
             eval(sprintf('D.%s.SC_R= D_;',Delays_{iDelay}));
            end
             clear D_ Strials nS Ctrials nC SampleTrials ChoiceTrials iTrial

        end
        %% Decode Sample/Choice on error trials
        for iDelay =1:length(Delays_)
            % Run once for L trials
            eval(sprintf('SampleTrials = t.%s.SamplePress_LeftError;',Delays_{iDelay}));
            eval(sprintf('ChoiceTrials = t.%s.ChoicePress_LeftError;',Delays_{iDelay}));
            %eval(sprintf(';',Delays_{iDelay},Delays_{iDelay}));
           
            Strials = [];nS = 0;
            for iTrial =1:length(SampleTrials)
                try
                    eval(sprintf('tlims_ = SampleTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}));
                    tlims_ = closest(Tmtx,tlims_);
                    eval(sprintf('Strials = [Strials;iFR_(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}));
                    nS=nS+1;
                end
            end
            
            Ctrials = []; nC=0;
            for iTrial =1:length(ChoiceTrials)
                try
                    eval(sprintf('tlims_ = ChoiceTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}));
                    tlims_ = closest(Tmtx,tlims_);
                    eval(sprintf('Ctrials = [Ctrials;iFR_(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}));
                    nC=nC+1;
                end
            end
            if nS>2 && nC>2
                [D_.avgFR,D_.seFR,...
                    D_.Ft2,D_.Rt2,...
                    D_.Ft2ciL,D_.Ft2ciH,...
                    D_.TS,...
                    D_.dfnum,D_.dfden] = DecodeStats([Strials;Ctrials],[ones(nS,1);2*ones(nC,1)],0.05);
                D_.TSsig = tpdf(D_.TS,nS+nC-2)<0.05;
                eval(sprintf('D.%s.SCerr_L= D_;',Delays_{iDelay}));
            end
             clear D_ Strials nS Ctrials nC SampleTrials ChoiceTrials iTrial
    
            % Run once for R trials
            eval(sprintf('SampleTrials = t.%s.SamplePress_RightError;',Delays_{iDelay}));
            eval(sprintf('ChoiceTrials = t.%s.ChoicePress_RightError;',Delays_{iDelay}));
            %eval(sprintf(';',Delays_{iDelay},Delays_{iDelay}));
           
            Strials = [];nS = 0;
            for iTrial =1:length(SampleTrials)
                try
                    eval(sprintf('tlims_ = SampleTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}));
                    tlims_ = closest(Tmtx,tlims_);
                    eval(sprintf('Strials = [Strials;iFR_(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}));
                    nS=nS+1;
                end
            end
            
            Ctrials = []; nC=0;
            for iTrial =1:length(ChoiceTrials)
                try
                    eval(sprintf('tlims_ = ChoiceTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}));
                    tlims_ = closest(Tmtx,tlims_);
                    eval(sprintf('Ctrials = [Ctrials;iFR_(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}));
                    nC=nC+1;
                end
            end
            if nS>2 && nC>2
            [D_.avgFR,D_.seFR,...
             D_.Ft2,D_.Rt2,...
             D_.Ft2ciL,D_.Ft2ciH,...
             D_.TS,...
             D_.dfnum,D_.dfden] = DecodeStats([Strials;Ctrials],[ones(nS,1);2*ones(nC,1)],0.05);
             D_.TSsig = tpdf(D_.TS,nS+nC-2)<0.05;
             eval(sprintf('D.%s.SCerr_R= D_;',Delays_{iDelay}));
            end
             clear D_ Strials nS Ctrials nC SampleTrials ChoiceTrials iTrial

        end        
        %% Plot L/R decoding
        if plotOnline
            figure; hold on
            plot([0:(length(tbAll)*2-1)]*bw,D.Delay_2.LR.Ft2,'r')
            plot([0:(length(tbAll)*2-1)]*bw,D.Delay_4.LR.Ft2,'g')
            plot([0:(length(tbAll)*2-1)]*bw,D.Delay_6.LR.Ft2,'b')
            legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
            % plot(tbShort,D.ShortLR.Ft2ciH,':r')
            % plot(tbMedium,D.Medium.Ft2ciH,':g')
            % plot(tbLong,D.Long.Ft2ciH,':b')
        end
        %% Plot C/E decoding
        if plotOnline
            figure; hold on
            plot([0:(length(tbAll)*2-1)]*bw,mean([D.Delay_2.CE_L.Ft2;D.Delay_2.CE_R.Ft2]),'r')
            plot([0:(length(tbAll)*2-1)]*bw,mean([D.Delay_4.CE_L.Ft2;D.Delay_4.CE_R.Ft2]),'g')
            plot([0:(length(tbAll)*2-1)]*bw,mean([D.Delay_6.CE_L.Ft2;D.Delay_6.CE_R.Ft2]),'b')
            legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
            % plot(tbShort,D.ShortLR.Ft2ciH,':r')
            % plot(tbMedium,D.Medium.Ft2ciH,':g')
            % plot(tbLong,D.Long.Ft2ciH,':b')
        end
        %% Plot S/C decoding
        if plotOnline
            figure; hold on
            plot(tbShort,   mean([D.Delay_2.SC_L.Ft2;D.Delay_2.SC_R.Ft2]),'r')
            plot(tbMedium, mean([D.Medium.Delay_4.Ft2;D.Delay_4.SC_R.Ft2]),'g')
            plot(tbLong,    mean([D.Delay_6.SC_L.Ft2;D.Delay_6.SC_R.Ft2]),'b')
            legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
        end
        %% Trial Fano factor
         for iDelay =1:length(Delays_)
            eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
            fields = fieldnames(t_);
            fields_ = fields(cellfun(@length,strfind(fieldnames(t_),'Sample')) | cellfun(@length,strfind(fieldnames(t_),'Choice')));
            
                for iType = 1:length(fields_)
                
                eval(sprintf('Trials = t_.%s;',fields_{iType}));
                FR=[];
                for iTrial =1:length(Trials)
                    try
                        tlims_ =Trials(iTrial)/1e6+tlimsANOVA + shift;
                        tlims_ = closest(Tmtx,tlims_);
                        FR(iTrial,1:size(iFR_,2))= nanmean(iFR_(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
                    catch
                        FR(iTrial,1:size(iFR_,2))= nan(1,size(iFR_,2));
                    end
                end
                FFR(iType,:) = nanmean(FR);
                FFt(iType,:) = nanvar(FR)./nanmean(FR); % firing rate variability across all trials of each type for each neuron
            end
            A = [];
            idx = cellfun(@length,strfind(fields_,'Sample')) & cellfun(@length,strfind(fields_,'Correct'));
            A = [A,mean(FFt(idx,:))'];
            idx = cellfun(@length,strfind(fields_,'Choice')) & cellfun(@length,strfind(fields_,'Correct'));
            A = [A,mean(FFt(idx,:))'];
            idx = cellfun(@length,strfind(fields_,'Sample')) & cellfun(@length,strfind(fields_,'Error'));
            A = [A,mean(FFt(idx,:))'];
            idx = cellfun(@length,strfind(fields_,'Choice')) & cellfun(@length,strfind(fields_,'Error'));
            A = [A,mean(FFt(idx,:))'];            
            
            eval(sprintf('D.%s.Fano.FFt = A;',Delays_{iDelay}));
            eval(sprintf('D.%s.Fano.FFtmean = nanmean(A);',Delays_{iDelay}));
            eval(sprintf('D.%s.Fano.FFtSEM = nansem(A);',Delays_{iDelay}));
            eval(sprintf('D.%s.Fano.FFtNames = {''Sample (Correct)'',''Choice (Correct)'',''Sample (Error)'',''Choice (Error)''};',Delays_{iDelay}));
            eval(sprintf('D.%s.Fano.FFR = nanvar(FR)./nanmean(FR);',Delays_{iDelay}));
            
            clear A idx FFR FFt
            
         end
         
        %% ANOVA on L/R S/C C/E 
        q_names = {'Untuned','Pure only','Mixed only','Mixed and Pure','Any tuning'}; 
        q2_names = {'Untuned','CS','LMS','NMS',...
                    'LMS CxO','LMS CxP','LMS OxP','LMS CxOxP',...
                    'NMS CxO','NMS CxP','NMS OxP','NMS CxOxP'}; 

        for iDelay =1:length(Delays_)
            eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
            Trials = [t_.SamplePress_LeftCorrect;...
                t_.SamplePress_LeftError';...
                t_.SamplePress_RightCorrect;...
                t_.SamplePress_RightError';...
                t_.ChoicePress_LeftCorrect;...
                t_.ChoicePress_LeftError';...
                t_.ChoicePress_RightCorrect;...
                t_.ChoicePress_RightError'];
            SC = [repmat({'Sample'},length(t_.SamplePress_LeftCorrect),1);...
                repmat({'Sample'},length(t_.SamplePress_LeftError),1);...
                repmat({'Sample'},length(t_.SamplePress_RightCorrect),1);...
                repmat({'Sample'},length(t_.SamplePress_RightError),1);...
                repmat({'Choice'},length(t_.ChoicePress_LeftCorrect),1);...
                repmat({'Choice'},length(t_.ChoicePress_LeftError),1);...
                repmat({'Choice'},length(t_.ChoicePress_RightCorrect),1);...
                repmat({'Choice'},length(t_.ChoicePress_RightError),1)];
            
            CE = [repmat({'Correct'},length(t_.SamplePress_LeftCorrect),1);...
                repmat({'Error'},length(t_.SamplePress_LeftError),1);...
                repmat({'Correct'},length(t_.SamplePress_RightCorrect),1);...
                repmat({'Error'},length(t_.SamplePress_RightError),1);...
                repmat({'Correct'},length(t_.ChoicePress_LeftCorrect),1);...
                repmat({'Error'},length(t_.ChoicePress_LeftError),1);...
                repmat({'Correct'},length(t_.ChoicePress_RightCorrect),1);...
                repmat({'Error'},length(t_.ChoicePress_RightError),1)];
            
            LR = [repmat({'Left'},length(t_.SamplePress_LeftCorrect),1);...
                repmat({'Left'},length(t_.SamplePress_LeftError),1);...
                repmat({'Right'},length(t_.SamplePress_RightCorrect),1);...
                repmat({'Right'},length(t_.SamplePress_RightError),1);...
                repmat({'Left'},length(t_.ChoicePress_LeftCorrect),1);...
                repmat({'Left'},length(t_.ChoicePress_LeftError),1);...
                repmat({'Right'},length(t_.ChoicePress_RightCorrect),1);...
                repmat({'Right'},length(t_.ChoicePress_RightError),1)];
            FR =[];
            for iTrial =1:length(Trials)
                try
                    tlims_ =Trials(iTrial)/1e6+tlimsANOVA + shift;
                    tlims_ = closest(Tmtx,tlims_);
                    FR(iTrial,1:size(iFR_,2))= mean(iFR_(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
                catch
                    FR(iTrial,1:size(iFR_,2))= nan(1,size(iFR_,2));
                end
            end
            SC(isnan(sum(FR,2)))=[];
            CE(isnan(sum(FR,2)))=[];
            LR(isnan(sum(FR,2)))=[];
            FR(isnan(sum(FR,2)),:)=[];
            p_ = zeros(size(iFR_,2),7); F_=p_;
            

            for iUnit=1:size(iFR_,2)
                [p_(iUnit,:),tbl_] = anovan(FR(:,iUnit),{SC CE LR},'model','full',...
                    'varnames',{'Context','Outcome','Position'},'display','off');
                F_(iUnit,:)=[tbl_{2:8,6}];
                % [results,~,~,gnames] = multcompare(stats,'Dimension',[1 2 3])
            end
            p_thresh = p_<0.05;
            
            %Fractions of units with different tunings:
            q_ = [sum(sum(p_thresh,2)==0),...                                      % No tuning
                  sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)==0),...   % Pure tuning only
                  sum(sum(p_thresh(:,1:3),2)==0 & sum(p_thresh(:,4:7),2)>0),...    % Mixed tuning only
                  sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)>0),...    % Mixed and pure tuning
                  sum(sum(p_thresh,2)>0)];                                          % Any tuning
              
            q2_= [sum(sum(p_thresh,2)==0),...                                       % No tuning
                  sum(sum(p_thresh(:,1:3),2)==1 & sum(p_thresh(:,4:7),2)==0),...    % CS:  Classical selectivity for one variable
                  sum(sum(p_thresh(:,1:3),2)>1  & sum(p_thresh(:,4:7),2)==0),...    % LMS: Linear mixed selectivity for 2+ variables
                  sum(sum(p_thresh(:,4:7),2)>1),...                                      % NMS: Non-linear mixed selectivity
                  sum(sum(p_thresh == [1 0 0 0 0 0 0],2)==7),...     CS for Context
                  sum(sum(p_thresh == [0 1 0 0 0 0 0],2)==7),...     CS for Outcome
                  sum(sum(p_thresh == [0 0 1 0 0 0 0],2)==7),...     CS for Position
                  sum(sum(p_thresh == [1 1 0 0 0 0 0],2)==7),...     LMS for Context X Outcome
                  sum(sum(p_thresh == [1 0 1 0 0 0 0],2)==7),...     LMS for Context X Position
                  sum(sum(p_thresh == [0 1 1 0 0 0 0],2)==7),...     LMS for Outcome X Position
                  sum(sum(p_thresh == [1 1 1 0 0 0 0],2)==7),...     LMS for Context X Outcome X Position
                  sum(sum(p_thresh(:,[4 7])==[1 0],2)==2),...        NMS for Context X Outcome
                  sum(sum(p_thresh(:,[5 7])==[1 0],2)==2),...        NMS for Context X Position
                  sum(sum(p_thresh(:,[6 7])==[1 0],2)==2),...        NMS for Outcome X Position
                  sum(sum(p_thresh(:,7) == 1 ,2)==2),...             NMS for Outcome X Position X Position
                  ]; 
%             q_ = [sum((sum(p_'<0.05)<1)), ...        % No tuning
%                 sum((sum(p_(:,1:3)'<0.05)>1)), ... % Pure tuning only
%                 sum((sum(p_(:,4:7)'<0.05)>1))];    % Mixed tuning
            eval(sprintf('D.%s.ANOVA.p=p_;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA.F=F_;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA.Factors=transpose({tbl_{2:8,1}});',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA.prcSig=sum(p_<0.05)./size(iFR_,2)*100;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA.prcTuned=q_./size(iFR_,2)*100;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA.prcTuned2=q2_./size(iFR_,2)*100;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA.tuning=q_names;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA.tuning2=q2_names;',Delays_{iDelay}));
                        
        end
        
        clear p_ q_ q2_ q_names q2_names FR SC LR CE p_thresh
        %% ANOVA on L/R S/C correct
        q_names = {'Untuned','CS','LMS','NMS'}; 

        for iDelay =1:length(Delays_)
            
            eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
            
            Trials = [t_.SamplePress_LeftCorrect;...
                      t_.SamplePress_RightCorrect;...
                      t_.ChoicePress_LeftCorrect;...
                      t_.ChoicePress_RightCorrect];
                  
            SC = [repmat({'Sample'},length(t_.SamplePress_LeftCorrect),1);...
                  repmat({'Sample'},length(t_.SamplePress_RightCorrect),1);...
                  repmat({'Choice'},length(t_.ChoicePress_LeftCorrect),1);...
                  repmat({'Choice'},length(t_.ChoicePress_RightCorrect),1)];
            
            
            LR = [repmat({'Left'},length(t_.SamplePress_LeftCorrect),1);...
                  repmat({'Right'},length(t_.SamplePress_RightCorrect),1);...
                  repmat({'Left'},length(t_.ChoicePress_LeftCorrect),1);...
                  repmat({'Right'},length(t_.ChoicePress_RightCorrect),1)];
              
            FR =[];
            for iTrial =1:length(Trials)
                try
                    tlims_ =Trials(iTrial)/1e6+tlimsANOVA + shift;
                    tlims_ = closest(Tmtx,tlims_);
                    FR(iTrial,1:size(iFR_,2))= mean(iFR_(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
                catch
                    FR(iTrial,1:size(iFR_,2))= nan(1,size(iFR_,2));
                end
            end
            
            SC(isnan(sum(FR,2)))=[];
            LR(isnan(sum(FR,2)))=[];
            FR(isnan(sum(FR,2)),:)=[];
            p_ = zeros(size(iFR_,2),3); F_=p_;
            

            for iUnit=1:size(iFR_,2)
                [p_(iUnit,:),tbl_] = anovan(FR(:,iUnit),{SC LR},'model','full',...
                    'varnames',{'Context','Position'},'display','off');
                F_(iUnit,:)=[tbl_{2:4,6}];
                % [results,~,~,gnames] = multcompare(stats,'Dimension',[1 2 3])
            end
            p_thresh = p_<0.05;
            
            %Fractions of units with different tunings:
            q_ = [sum(sum(p_thresh,2)==0),...                % No tuning
                  sum(sum(p_thresh == [1 0 0],2)==3),...     % CS for Context
                  sum(sum(p_thresh == [0 1 0],2)==3),...     % CS for Outcome
                  sum(sum(p_thresh == [1 1 0],2)==3),...     % LMS for Context x Outcome 
                  sum(sum(p_thresh(:,3) == 1,2)),...         % NMS for Context x Outcome 
                  sum(sum(p_thresh,2)>0)];                   % Any tuning
              
            
            eval(sprintf('D.%s.ANOVA2.p_cor=p_;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.F_cor=F_;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.Factors_cor=transpose({tbl_{2:4,1}});',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.prcSig_cor=sum(p_<0.05)./size(iFR_,2)*100;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.prcTuned_cor=q_./size(iFR_,2)*100;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.tuning_cor=q_names;',Delays_{iDelay}));
                        
        end
        clear p_ q_ q_names  FR SC LR CE p_thresh
        %% ANOVA on L/R S/C error
        q_names = {'Untuned','CS','LMS','NMS'}; 

        for iDelay =1:length(Delays_)
            
            eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
            
            Trials = [t_.SamplePress_LeftError';...
                      t_.SamplePress_RightError';...
                      t_.ChoicePress_LeftError';...
                      t_.ChoicePress_RightError'];
            SC = [repmat({'Sample'},length(t_.SamplePress_LeftError),1);...
                repmat({'Sample'},length(t_.SamplePress_RightError),1);...
                repmat({'Choice'},length(t_.ChoicePress_LeftError),1);...
                repmat({'Choice'},length(t_.ChoicePress_RightError),1)];
            
            LR = [repmat({'Left'},length(t_.SamplePress_LeftError),1);...
                repmat({'Right'},length(t_.SamplePress_RightError),1);...
                repmat({'Left'},length(t_.ChoicePress_LeftError),1);...
                repmat({'Right'},length(t_.ChoicePress_RightError),1)];
              
            FR =[];
            for iTrial =1:length(Trials)
                try
                    tlims_ =Trials(iTrial)/1e6+tlimsANOVA + shift;
                    tlims_ = closest(Tmtx,tlims_);
                    FR(iTrial,1:size(iFR_,2))= mean(iFR_(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
                catch
                    FR(iTrial,1:size(iFR_,2))= nan(1,size(iFR_,2));
                end
            end
            
            SC(isnan(sum(FR,2)))=[];
            LR(isnan(sum(FR,2)))=[];
            FR(isnan(sum(FR,2)),:)=[];
            p_ = zeros(size(iFR_,2),3); F_=p_;
            

            for iUnit=1:size(iFR_,2)
                [p_(iUnit,:),tbl_] = anovan(FR(:,iUnit),{SC LR},'model','full',...
                    'varnames',{'Context','Position'},'display','off');
                F_(iUnit,:)=[tbl_{2:4,6}];
                % [results,~,~,gnames] = multcompare(stats,'Dimension',[1 2 3])
            end
            p_thresh = p_<0.05;
            
            %Fractions of units with different tunings:
            q_ = [sum(sum(p_thresh,2)==0),...                % No tuning
                  sum(sum(p_thresh == [1 0 0],2)==3),...     % CS for Context
                  sum(sum(p_thresh == [0 1 0],2)==3),...     % CS for Outcome
                  sum(sum(p_thresh == [1 1 0],2)==3),...     % LMS for Context x Outcome 
                  sum(sum(p_thresh(:,3) == 1,2)),...         % NMS for Context x Outcome 
                  sum(sum(p_thresh,2)>0)];                   % Any tuning
              
            
            eval(sprintf('D.%s.ANOVA2.p_err=p_;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.F_err=F_;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.Factors_err=transpose({tbl_{2:4,1}});',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.prcSig_err=sum(p_<0.05)./size(iFR_,2)*100;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.prcTuned_err=q_./size(iFR_,2)*100;',Delays_{iDelay}));
            eval(sprintf('D.%s.ANOVA2.tuning_err=q_names;',Delays_{iDelay}));
                        
        end
        clear p_ q_ q_names  FR SC LR CE p_thresh
        %% Save results
        fnOut = sprintf('%sMixedSelectivity_SHORT%s%s_%s_MixedSelectivity_Units.mat',pat,filesep,fname,Areas{iArea});
        save(fnOut,'D','tbDelay_0','tbAll','-v7.3');
        
        fprintf('Done.\n')
    end
    end
end
clearvars -except Nbs AssemblyChoice color_ Delays_ bw tbShort tbAll tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin 

%% Batch Process assemblies 

for iFile =1:length(fileListAss)
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
                    
                    % Multivariate F-score decoder  
                    [D_.LR{iArea}.avgFR,...
                        D_.LR{iArea}.seFR,...
                        D_.LR{iArea}.Ft2,...
                        D_.LR{iArea}.Rt2,...
                        D_.LR{iArea}.Ft2ciL,...
                        D_.LR{iArea}.Ft2ciH,...
                        D_.LR{iArea}.TS,...
                        D_.LR{iArea}.dfnum,...
                        D_.LR{iArea}.dfden] = DecodeStats(FSC_,evt0,0.05);
                    D_.LR{iArea}.TSsig = tpdf(D_.LR{iArea}.TS,nL+nR-2)<0.05;
                    
                     % make trial counts consistent for comparisons
                     clear E
                     for nrep = 1:100
                         nrep
                         rsL = randsample(length(Ltrials_),5);
                         rsR = randsample(length(Rtrials_),5);
                         
                         FR   = [cell2mat(Ltrials_(rsL)');cell2mat(Rtrials_(rsR)')];
                         evt0 = [ones(5,1);2*ones(5,1)];
                         
                         % Multivariate F-score decoder
                         [E.avgFR(:,:,nrep),E.seFR(:,:,nrep),...
                             E.Ft2(:,nrep),E.Rt2(:,nrep),...
                             E.Ft2ciL(:,nrep),E.Ft2ciH(:,nrep),...
                             E.TS(:,:,nrep),...
                             ~,~] = DecodeStats(FR,evt0,0.05);
                         E.TSsig(:,:,nrep) = tpdf(E.TS(:,:,nrep),8)<0.05;
                     end
                     D_.TS_drawnTrials     =  mean(E.TS,3);
                     D_.TSsig_drawnTrials  =  sum(E.TSsig,3)./100;
                 
                 
                    % Cross-validation error decoder
                    FSC_ = [Ltrials;Rtrials];
                    evt0 = [ones(nL,1);2*ones(nR,1)];
                    for iAss =1:size(FSC,2)
                        iAss
                        D_.CVEindividual(iAss,:) = DecodeCVE_mex(FSC_(:,iAss),evt0,0.05);
                    end
                    if runCVEdecoder
                        if size(Ltrials,2)==min(nAss_), nrep=1; else nrep = 10; end;
                    % nrep=1;
                    for i=1:nrep
                        i
                        rs{i} = randsample(size(Ltrials,2),min(nAss_));
                        % rs{i} = randsample(size(Ltrials,2),nu(iArea));
%                         D_.LR{iArea}.CVE(i,:) = DecodeCVE(FSC_,evt0,0.05);
                        D_.LR{iArea}.CVE(i,:) = DecodeCVE_mex(FSC_(:,rs{i}),evt0,0.05);
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
                        for i=1:nrep, cve0(i,:)=DecodeCVE_mex(FSC_(:,rs{i}),evt1,0.05); end
                        CVEbs(b,:) = nanmean(cve0,1);
                        fprintf('Decoded %s assemblies, correct trials, %s delay, draw %d of %d (%s elapsed)...\n',Areas{iArea},Delays_{iDelay},b,Nbs,seconds2human(toc(t1)))                    
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
        
%         figure; hold on
%         for s=1:3
%             plot(1-D_.LR{s}.CVE)
%         end
        if exist('D_')
            eval(sprintf('D_Ass.%s.LR =  D_.LR;',Delays_{iDelay}))
        end
        clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
    end
%     figure
%     color_ = {[0 0 0],0.3*[1 1 1],0.6*[1 1 1]};
%     for iArea = 1:3
%         subplot(3,1,iArea); hold on
%         for iDelay =1:3
%             eval(sprintf('y = D_Ass.%s.LR{iArea}.CVE;',Delays_{iDelay}))
%             plot(1-y,'color',color_{iDelay})
%         end
%         axis([ 0 Inf 0 1])
%     end
%         
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
                
                if nL>5 && nR>5
                    
                    FSC_ = [Ltrials;Rtrials];
                    evt0 = [ones(nL,1);2*ones(nR,1)];
                    Ltr = 2*length(tbAll);
                    nAss_ = cellfun(@(x)size(x,2),Ass.FSCsel);
                    % remove zero entries
                    idx = nAss_>0;
                    nAss_(~idx) = min(nAss_(idx));
                    
                    % Multivariate F-score decoder
                    [D_.LR{iArea}.avgFR,...
                        D_.LR{iArea}.seFR,...
                        D_.LR{iArea}.Ft2,...
                        D_.LR{iArea}.Rt2,...
                        D_.LR{iArea}.Ft2ciL,...
                        D_.LR{iArea}.Ft2ciH,...
                        D_.LR{iArea}.TS,...
                        D_.LR{iArea}.dfnum,...
                        D_.LR{iArea}.dfden] = DecodeStats(FSC_,evt0,0.05);
                    D_.LR{iArea}.TSsig = tpdf(D_.LR{iArea}.TS,nL+nR-2)<0.05;
                    
                    % Cross-validation error decoder
                    if runCVEdecoder
                        if size(Ltrials,2)==min(nAss_), nrep=1; else nrep=10; end
                        %                  nrep=1;
                        for i=1:nrep
                            rs{i} = randsample(size(Ltrials,2),min(nAss_));
                            %                         D_.LR{iArea}.CVE(i,:) = DecodeCVE(FSC_,evt0,0.05);
                            D_.LR{iArea}.CVE(i,:) = DecodeCVE_mex(FSC_(:,rs{i}),evt0,0.05);
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
                            for i=1:nrep, cve0(i,:)=DecodeCVE_mex(FSC_,evt1,0.05); end
                            CVEbs(b,:) = nanmean(cve0,1);
                            fprintf('Decoded %s assemblies, error trials, %s delay, draw %d of %d (%s elapsed)...\n',Areas{iArea},Delays_{iDelay},b,Nbs,seconds2human(toc(t1)))
                            
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
            eval(sprintf('D_Ass.%s.LR_err =  D_.LR;',Delays_{iDelay}))
        end
        clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
    end
    %% Decode Correct/Error
     for iDelay =1:length(Delays_)
         
         % Run once for L trials
         eval(sprintf('CorrectTrials = [t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect];',Delays_{iDelay},Delays_{iDelay}));
         eval(sprintf('ErrorTrials =   [t.%s.SamplePress_LeftError'',t.%s.ChoicePress_LeftError''];',Delays_{iDelay},Delays_{iDelay}));
            
         for iArea = 1:length(Ass.FSCsel)
             
             
             FSC = Ass.FSCsel{iArea};
             if ~isempty(FSC)
                 
                 Ctrials = [];nC = 0;
                 for iTrial =1:size(CorrectTrials,1)
                     try
                         %eval(sprintf('tlims_ = CorrectTrials(iTrial,1)/1e6+tlims%s + shift;',Delays_{iDelay}))
                         tlims_  = [CorrectTrials(iTrial,1)/1e6;CorrectTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                         tlims_ = closest(Ass.Tmtx,tlims_);
                         tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                         %eval(sprintf('Ctrials = [Ctrials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                         Ctrials = [Ctrials;FSC(tlims_,:)];
                         nC=nC+1;
                     end
                 end
                 
                 Etrials = []; nE=0;           
                 for iTrial =1:size(ErrorTrials,1)
                     try
                         %eval(sprintf('tlims_ = ErrorTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                         tlims_  = [ErrorTrials(iTrial,1)/1e6;ErrorTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                         tlims_ = closest(Ass.Tmtx,tlims_);
                         tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                         %eval(sprintf('Etrials = [Etrials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                         Etrials = [Etrials;FSC(tlims_,:)];
                         nE=nE+1;
                     end
                 end
                 
                 if nC>1 && nE>1
                    [D_.CE{iArea}.avgFR,...
                     D_.CE{iArea}.seFR,...
                     D_.CE{iArea}.Ft2,...
                     D_.CE{iArea}.Rt2,...
                     D_.CE{iArea}.Ft2ciL,...
                     D_.CE{iArea}.Ft2ciH,...
                     D_.CE{iArea}.TS,...
                     D_.CE{iArea}.dfnum,...
                     D_.CE{iArea}.dfden] = DecodeStats([Ctrials;Etrials],[ones(nC,1);2*ones(nE,1)],0.05);
                     D_.CE{iArea}.TSsig = tpdf(D_.CE{iArea}.TS,nC+nE-2)<0.05;
                 end
             end
         end
         if exist('D_')
             eval(sprintf('D_Ass.%s.CE_L =  D_.CE;',Delays_{iDelay}))
         end
         clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
         
         % Run once for L trials
         eval(sprintf('CorrectTrials = [t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}));
         eval(sprintf('ErrorTrials =   [t.%s.SamplePress_RightError'',t.%s.ChoicePress_RightError''];',Delays_{iDelay},Delays_{iDelay}));
            
         for iArea = 1:length(Ass.FSCsel)
             
             FSC = Ass.FSCsel{iArea};
             if ~isempty(FSC)
                 
                 Ctrials = [];nC = 0;
                 for iTrial =1:size(CorrectTrials,1)
                     try
                         %eval(sprintf('tlims_ = CorrectTrials(iTrial,1)/1e6+tlims%s + shift;',Delays_{iDelay}))
                         tlims_  = [CorrectTrials(iTrial,1)/1e6;CorrectTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                         tlims_ = closest(Ass.Tmtx,tlims_);
                         tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                         %eval(sprintf('Ctrials = [Ctrials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                         Ctrials = [Ctrials;FSC(tlims_,:)];
                         nC=nC+1;
                     end
                 end
                 
                 Etrials = []; nE=0;           
                 for iTrial =1:size(ErrorTrials,1)
                     try
                         %eval(sprintf('tlims_ = ErrorTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                         tlims_  = [ErrorTrials(iTrial,1)/1e6;ErrorTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                         tlims_ = closest(Ass.Tmtx,tlims_);
                         tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                         %eval(sprintf('Etrials = [Etrials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                         Etrials = [Etrials;FSC(tlims_,:)];
                         nE=nE+1;
                     end
                 end
                 
                 if nC>1 && nE>1
                    [D_.CE{iArea}.avgFR,...
                     D_.CE{iArea}.seFR,...
                     D_.CE{iArea}.Ft2,...
                     D_.CE{iArea}.Rt2,...
                     D_.CE{iArea}.Ft2ciL,...
                     D_.CE{iArea}.Ft2ciH,...
                     D_.CE{iArea}.TS,...
                     D_.CE{iArea}.dfnum,...
                     D_.CE{iArea}.dfden] = DecodeStats([Ctrials;Etrials],[ones(nC,1);2*ones(nE,1)],0.05);
                     D_.CE{iArea}.TSsig = tpdf(D_.CE{iArea}.TS,nC+nE-2)<0.05;
                 end
             end
         end
         if exist('D_')
             eval(sprintf('D_Ass.%s.CE_R =  D_.CE;',Delays_{iDelay}))
         end
         clear D_ Ctrials nC Etrials nE CorrectTrials ErrorTrials iTrial

    end
    %% Decode Sample/Choice on correct trials
    for iDelay =1:length(Delays_)
        % Run once for L trials 
        eval(sprintf('SampleTrials = t.%s.SamplePress_LeftCorrect;',Delays_{iDelay}));
        eval(sprintf('ChoiceTrials = t.%s.ChoicePress_LeftCorrect;',Delays_{iDelay}));
        for iArea = 1:length(Ass.FSCsel)

            FSC = Ass.FSCsel{iArea};
            if ~isempty(FSC)
                
                Strials = [];nS = 0;
                for iTrial =1:length(SampleTrials)
                     try
                        eval(sprintf('tlims_ = SampleTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        eval(sprintf('Strials = [Strials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                        nS=nS+1;
                     end
                end
                Ctrials = [];nC = 0;
                for iTrial =1:length(ChoiceTrials)
                     try
                        eval(sprintf('tlims_ = ChoiceTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        eval(sprintf('Ctrials = [Ctrials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                        nC=nC+1;
                     end
                 end

                 if nS>1 && nC>1
                    [D_.SC{iArea}.avgFR,...
                     D_.SC{iArea}.seFR,...
                     D_.SC{iArea}.Ft2,...
                     D_.SC{iArea}.Rt2,...
                     D_.SC{iArea}.Ft2ciL,...
                     D_.SC{iArea}.Ft2ciH,...
                     D_.SC{iArea}.TS,...
                     D_.SC{iArea}.dfnum,...
                     D_.SC{iArea}.dfden] = DecodeStats([Strials;Ctrials],[ones(nS,1);2*ones(nC,1)],0.05);
                     D_.SC{iArea}.TSsig = tpdf(D_.SC{iArea}.TS,nS+nC-2)<0.05;
                 end
            end
        end
        if exist('D_')
        	eval(sprintf('D_Ass.%s.SC_L =  D_.SC;',Delays_{iDelay}))
            clear D_
        end
        clear D_ Strials nS Ctrials nC SampleTrials ChoiceTrials iTrial

        % Run once for R trials 
        eval(sprintf('SampleTrials = t.%s.SamplePress_RightCorrect;',Delays_{iDelay}));
        eval(sprintf('ChoiceTrials = t.%s.ChoicePress_RightCorrect;',Delays_{iDelay}));
        for iArea = 1:length(Ass.FSCsel)

            FSC = Ass.FSCsel{iArea};
            if ~isempty(FSC)
                
                Strials = [];nS = 0;
                for iTrial =1:length(SampleTrials)
                     try
                        eval(sprintf('tlims_ = SampleTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        eval(sprintf('Strials = [Strials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                        nS=nS+1;
                     end
                end
                Ctrials = [];nC = 0;
                for iTrial =1:length(ChoiceTrials)
                     try
                        eval(sprintf('tlims_ = ChoiceTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        eval(sprintf('Ctrials = [Ctrials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                        nC=nC+1;
                     end
                 end

                 if nS>1 && nC>1
                    [D_.SC{iArea}.avgFR,...
                     D_.SC{iArea}.seFR,...
                     D_.SC{iArea}.Ft2,...
                     D_.SC{iArea}.Rt2,...
                     D_.SC{iArea}.Ft2ciL,...
                     D_.SC{iArea}.Ft2ciH,...
                     D_.SC{iArea}.TS,...
                     D_.SC{iArea}.dfnum,...
                     D_.SC{iArea}.dfden] = DecodeStats([Strials;Ctrials],[ones(nS,1);2*ones(nC,1)],0.05);
                     D_.SC{iArea}.TSsig = tpdf(D_.SC{iArea}.TS,nS+nC-2)<0.05;
                 end
            end
        end
        if exist('D_')
        	eval(sprintf('D_Ass.%s.SC_R =  D_.SC;',Delays_{iDelay}))
            clear D_
        end
        clear D_ Strials nS Ctrials nC SampleTrials ChoiceTrials iTrial
    end
    %% Decode Sample/Choice on error trials
    for iDelay =1:length(Delays_)
        % Run once for L trials 
        eval(sprintf('SampleTrials = t.%s.SamplePress_LeftError;',Delays_{iDelay}));
        eval(sprintf('ChoiceTrials = t.%s.ChoicePress_LeftError;',Delays_{iDelay}));
        for iArea = 1:length(Ass.FSCsel)

            FSC = Ass.FSCsel{iArea};
            if ~isempty(FSC)
                
                Strials = [];nS = 0;
                for iTrial =1:length(SampleTrials)
                     try
                        eval(sprintf('tlims_ = SampleTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        eval(sprintf('Strials = [Strials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                        nS=nS+1;
                     end
                end
                Ctrials = [];nC = 0;
                for iTrial =1:length(ChoiceTrials)
                     try
                        eval(sprintf('tlims_ = ChoiceTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        eval(sprintf('Ctrials = [Ctrials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                        nC=nC+1;
                     end
                 end

                 if nS>1 && nC>1
                    [D_.SC{iArea}.avgFR,...
                     D_.SC{iArea}.seFR,...
                     D_.SC{iArea}.Ft2,...
                     D_.SC{iArea}.Rt2,...
                     D_.SC{iArea}.Ft2ciL,...
                     D_.SC{iArea}.Ft2ciH,...
                     D_.SC{iArea}.TS,...
                     D_.SC{iArea}.dfnum,...
                     D_.SC{iArea}.dfden] = DecodeStats([Strials;Ctrials],[ones(nS,1);2*ones(nC,1)],0.05);
                     D_.SC{iArea}.TSsig = tpdf(D_.SC{iArea}.TS,nS+nC-2)<0.05;
                 end
            end
        end
        if exist('D_')
        	eval(sprintf('D_Ass.%s.SCerr_L =  D_.SC;',Delays_{iDelay}))
            clear D_
        end
        clear D_ Strials nS Ctrials nC SampleTrials ChoiceTrials iTrial

        % Run once for R trials 
        eval(sprintf('SampleTrials = t.%s.SamplePress_RightError;',Delays_{iDelay}));
        eval(sprintf('ChoiceTrials = t.%s.ChoicePress_RightError;',Delays_{iDelay}));
        for iArea = 1:length(Ass.FSCsel)

            FSC = Ass.FSCsel{iArea};
            if ~isempty(FSC)
                
                Strials = [];nS = 0;
                for iTrial =1:length(SampleTrials)
                     try
                        eval(sprintf('tlims_ = SampleTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        eval(sprintf('Strials = [Strials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                        nS=nS+1;
                     end
                end
                Ctrials = [];nC = 0;
                for iTrial =1:length(ChoiceTrials)
                     try
                        eval(sprintf('tlims_ = ChoiceTrials(iTrial)/1e6+tlims%s + shift;',Delays_{iDelay}))
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        eval(sprintf('Ctrials = [Ctrials;FSC(tlims_(1):tlims_(1)+length(tb%s)-1,:)];',Delays_{iDelay}))
                        nC=nC+1;
                     end
                 end

                 if nS>1 && nC>1
                    [D_.SC{iArea}.avgFR,...
                     D_.SC{iArea}.seFR,...
                     D_.SC{iArea}.Ft2,...
                     D_.SC{iArea}.Rt2,...
                     D_.SC{iArea}.Ft2ciL,...
                     D_.SC{iArea}.Ft2ciH,...
                     D_.SC{iArea}.TS,...
                     D_.SC{iArea}.dfnum,...
                     D_.SC{iArea}.dfden] = DecodeStats([Strials;Ctrials],[ones(nS,1);2*ones(nC,1)],0.05);
                     D_.SC{iArea}.TSsig = tpdf(D_.SC{iArea}.TS,nS+nC-2)<0.05;
                 end
            end
        end
        if exist('D_')
        	eval(sprintf('D_Ass.%s.SCerr_R =  D_.SC;',Delays_{iDelay}))
            clear D_
        end
        clear D_ Strials nS Ctrials nC SampleTrials ChoiceTrials iTrial
    end      
    %% Trial Fano factor
    for iDelay =1:length(Delays_)
        eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
        fields = fieldnames(t_);
        fields_ = fields(cellfun(@length,strfind(fieldnames(t_),'Sample')) | cellfun(@length,strfind(fieldnames(t_),'Choice')));
        for iArea = 1:3
            for iType = 1:length(fields_)
                eval(sprintf('Trials = t_.%s;',fields_{iType}));
                nAss = size(Ass.FSCsel{iArea},2);
                FSC=[];
                for iTrial =1:length(Trials)
                    try
                        tlims_ = Trials(iTrial)/1e6+tlimsANOVA + shift;
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        FSC(iTrial,1:nAss) = nanmean(Ass.FSCsel{iArea}(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
                    catch
                        FSC(iTrial,1:nAss)= nan(1,nAss);
                    end
                end
                
                FFR{iArea}(iType,:)= abs(nanmean(FSC));
                FFt{iArea}(iType,:)= abs(nanvar(FSC)./nanmean(FSC)); % firing rate variability across all trials of each type for each neuron
                
            end
            
            A{iArea}=[];
            idx = cellfun(@length,strfind(fields_,'Sample')) & cellfun(@length,strfind(fields_,'Correct'));
            A{iArea} = [A{iArea},(nanmean(FFt{iArea}(idx,:)))'];
            idx = cellfun(@length,strfind(fields_,'Choice')) & cellfun(@length,strfind(fields_,'Correct'));
            A{iArea} = [A{iArea},(nanmean(FFt{iArea}(idx,:)))'];
            idx = cellfun(@length,strfind(fields_,'Sample')) & cellfun(@length,strfind(fields_,'Error'));
            A{iArea} = [A{iArea},(nanmean(FFt{iArea}(idx,:)))'];
            idx = cellfun(@length,strfind(fields_,'Choice')) & cellfun(@length,strfind(fields_,'Error'));
            A{iArea} = [A{iArea},(nanmean(FFt{iArea}(idx,:)))'];
            eval(sprintf('D_Ass.%s.Fano.FFt{iArea} = A{iArea};',Delays_{iDelay}));
            eval(sprintf('D_Ass.%s.Fano.FFtmean{iArea} = nanmean(A{iArea},1);',Delays_{iDelay}));
            eval(sprintf('D_Ass.%s.Fano.FFtSEM{iArea}  = nansem(A{iArea},1);',Delays_{iDelay}));
            eval(sprintf('D_Ass.%s.Fano.FFtNames = {''Sample (Correct)'',''Choice (Correct)'',''Sample (Error)'',''Choice (Error)''};',Delays_{iDelay}));
            eval(sprintf('D_Ass.%s.Fano.FFR{iArea} = nanvar(FSC,1)./nanmean(FSC,1);',Delays_{iDelay}));
        end
        
        clear A idx FFR FFt
        
    end
    %% ANOVA on L/R S/C C/E
    q_names = {'Untuned','Pure only','Mixed only','Mixed and Pure','Any tuning'};
    for iDelay =1:length(Delays_)
        
        eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
        
        Trials = [t_.SamplePress_LeftCorrect;...
            t_.SamplePress_LeftError';...
            t_.SamplePress_RightCorrect;...
            t_.SamplePress_RightError';...
            t_.ChoicePress_LeftCorrect;...
            t_.ChoicePress_LeftError';...
            t_.ChoicePress_RightCorrect;...
            t_.ChoicePress_RightError'];
        
        for iArea = 1:length(Ass.FSCsel)
            SC = [repmat({'Sample'},length(t_.SamplePress_LeftCorrect),1);...
                repmat({'Sample'},length(t_.SamplePress_LeftError),1);...
                repmat({'Sample'},length(t_.SamplePress_RightCorrect),1);...
                repmat({'Sample'},length(t_.SamplePress_RightError),1);...
                repmat({'Choice'},length(t_.ChoicePress_LeftCorrect),1);...
                repmat({'Choice'},length(t_.ChoicePress_LeftError),1);...
                repmat({'Choice'},length(t_.ChoicePress_RightCorrect),1);...
                repmat({'Choice'},length(t_.ChoicePress_RightError),1)];
            
            CE = [repmat({'Correct'},length(t_.SamplePress_LeftCorrect),1);...
                repmat({'Error'},length(t_.SamplePress_LeftError),1);...
                repmat({'Correct'},length(t_.SamplePress_RightCorrect),1);...
                repmat({'Error'},length(t_.SamplePress_RightError),1);...
                repmat({'Correct'},length(t_.ChoicePress_LeftCorrect),1);...
                repmat({'Error'},length(t_.ChoicePress_LeftError),1);...
                repmat({'Correct'},length(t_.ChoicePress_RightCorrect),1);...
                repmat({'Error'},length(t_.ChoicePress_RightError),1)];
            
            LR = [repmat({'Left'},length(t_.SamplePress_LeftCorrect),1);...
                repmat({'Left'},length(t_.SamplePress_LeftError),1);...
                repmat({'Right'},length(t_.SamplePress_RightCorrect),1);...
                repmat({'Right'},length(t_.SamplePress_RightError),1);...
                repmat({'Left'},length(t_.ChoicePress_LeftCorrect),1);...
                repmat({'Left'},length(t_.ChoicePress_LeftError),1);...
                repmat({'Right'},length(t_.ChoicePress_RightCorrect),1);...
                repmat({'Right'},length(t_.ChoicePress_RightError),1)];
            
            
            FSC = [];
           
            nAss = size(Ass.FSCsel{iArea},2);
            if nAss>0
                 for iTrial =1:length(Trials)
                try
                    tlims_ = Trials(iTrial)/1e6+tlimsANOVA + shift;
                    tlims_ = closest(Ass.Tmtx,tlims_);
                    FSC(iTrial,1:nAss) = mean(Ass.FSCsel{iArea}(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
                catch
                    FSC(iTrial,1:nAss)= nan(1,nAss);
                end
            end
            SC(isnan(sum(FSC,2)))=[];
            CE(isnan(sum(FSC,2)))=[];
            LR(isnan(sum(FSC,2)))=[];
            FSC(isnan(sum(FSC,2)),:)=[];
            p_ = zeros(nAss,7); F_=p_;
            
            for iAss=1:nAss
                [p_(iAss,:),tbl_] = anovan(FSC(:,iAss),{SC CE LR},'model','full',...
                    'varnames',{'Context','Outcome','Position'},'display','off');
                F_(iAss,:)=[tbl_{2:8,6}];
                % [results,~,~,gnames] = multcompare(stats,'Dimension',[1 2 3])
            end
            p_thresh = p_<0.05;
            q_ = [sum(sum(p_thresh,2)==0),...                                    % No tuning
                sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)==0),...   % Pure tuning only
                sum(sum(p_thresh(:,1:3),2)==0 & sum(p_thresh(:,4:7),2)>0),...    % Mixed tuning only
                sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)>0),...    % Mixed and pure tuning
                sum(sum(p_thresh,2)>0)];                                         % Any tuning
            %             q_ = [sum((sum(p_'<0.05)<1)), ...        % No tuning
            %                 sum((sum(p_(:,1:3)'<0.05)>1)), ... % Pure tuning only
            %                 sum((sum(p_(:,4:7)'<0.05)>1))];    % Mixed tuning
            eval(sprintf('D_Ass.%s.ANOVA.p{iArea}=p_;',Delays_{iDelay}));
            eval(sprintf('D_Ass.%s.ANOVA.F{iArea}=F_;',Delays_{iDelay}));
            eval(sprintf('D_Ass.%s.ANOVA.Factors{iArea}=transpose({tbl_{2:8,1}});',Delays_{iDelay}));
            eval(sprintf('D_Ass.%s.ANOVA.prcSig{iArea}=sum(p_<0.05)./nAss*100;',Delays_{iDelay}));
            eval(sprintf('D_Ass.%s.ANOVA.prcTuned{iArea}=q_./nAss*100;',Delays_{iDelay}));
            eval(sprintf('D_Ass.%s.ANOVA.tuning{iArea}=q_names;',Delays_{iDelay}));
            clear D_
            end
        end
    end
    
    clear p_ q_ q_names FR SC LR CE p_thresh
    %% Save
    fnOut = sprintf('%sMixedSelectivity_SHORT%s%s_MixedSelectivity_Ass.mat',pat,filesep,fname);
    save(fnOut,'D_Ass','tbDelay_0','tbAll','-v7.3');
    fprintf('Done.\n')
    clear D_Ass
end
clearvars -except Nbs AssemblyChoice color_ Delays_ bw tbShort tbAll tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin 

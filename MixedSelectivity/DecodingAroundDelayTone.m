%% %%%%%% PREAMBLE %%%%%%
clear 
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Target = 'LONG';



tlimsAll = [-5 5];
tlimsShort=[-4 10];
tlimsMedium=[-8 10];
tlimsLong=[-16 10];
tlimsANOVA = [-4 4];

shift = 0;
plotOnline = false;
bw=0.05;
Nbs = 500;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbANOVA=tlimsANOVA(1):bw:tlimsANOVA(2);
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
        %% DelayEnd Decoding for units
        for iDelay =1:length(Delays_)
            
            eval(sprintf('LeftTrials = [t.%s.DelayEnd_LeftCorrect];',Delays_{iDelay}));
            eval(sprintf('RightTrials = [t.%s.DelayEnd_RightCorrect];',Delays_{iDelay}));
            eval(sprintf('LeftTrialsE = [t.%s.DelayEnd_LeftError''];',Delays_{iDelay}));
            eval(sprintf('RightTrialsE = [t.%s.DelayEnd_RightError''];',Delays_{iDelay}));

            Ltrials = [];nL = 0;            
            for iTrial =1:size(LeftTrials,1)
                 try
                    tlims_  = LeftTrials(iTrial,1)/1e6+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];
                    Ltrials = [Ltrials;iFR_(tlims_,:)];
                    nL=nL+1;
                 end
            end
            
            Rtrials = [];nR = 0;            
            for iTrial =1:size(RightTrials,1)
                try
                    tlims_  = RightTrials(iTrial,1)/1e6+tlimsAll(1) + shift;
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];
                    Rtrials = [Rtrials;iFR_(tlims_,:)];
                    nR=nR+1;
                end
            end
            
            LtrialsE = [];nLe = 0;
            for iTrial =1:size(LeftTrialsE,1)
                try
                    tlims_  = LeftTrialsE(iTrial,1)/1e6+tlimsAll(1) + shift;
                    tlims_ = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];
                    
                    LtrialsE = [LtrialsE;iFR_(tlims_,:)];
                    nLe=nLe+1;
                end
            end
            
            RtrialsE = [];nRe = 0;
            for iTrial =1:size(RightTrialsE,1)
                    try
                        tlims_  = RightTrialsE(iTrial,1)/1e6+tlimsAll(1) + shift;
                        tlims_ = closest(Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];

                        RtrialsE = [RtrialsE;iFR_(tlims_,:)];
                        nRe=nRe+1;
                    end
            end
                
            % Cross-validation error decoder (correct trials)    
%             if nL>2 && nR>2
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = length(tbAll);
                nu = [length(PFCcells),length(HPcells)];
                
                 if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
                 
                 for i=1:nrep
                     rs{i} = randsample(size(Ltrials,2),min(nu));
                     D_.CVE(i,:) = DecodeCVE_mex(FR(:,rs{i}),evt0,0.05);
                 end
                 D_.CVE = nanmean(D_.CVE,1);
                
                 clear rs nu evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr

%             end
            if exist('D_')
                eval(sprintf('D.%s.Tone =  D_;',Delays_{iDelay}))
            end
            clear D_
        
             % Cross-validation error decoder (error trials)    
%              if nLe>2 && nRe>2
                  FR = [Ltrials;Rtrials;LtrialsE;RtrialsE];
                evt0 = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
                Ltr = 2*length(tbAll);
                nu = [length(PFCcells),length(HPcells)];

                 if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
                 
                 for i=1:nrep
                     rs{i} = randsample(size(Ltrials,2),min(nu));
                     D_.CVE(i,:) = DecodePredErr(FR(:,rs{i}),evt0,0.05);
                 end
                 D_.CVE = nanmean(D_.CVE,1);
                 
                 clear rs nu evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
                 
%              end
             if exist('D_')
                 eval(sprintf('D.%s.Tone_err =  D_;',Delays_{iDelay}))
             end
             clear D_
             clear Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
        end
        %% Save results
        if exist('D')
            
            fnOut = sprintf('%sMixedSelectivity%s%s_%s_DelayEndDecoding_Units.mat',pat,filesep,fname,Areas{iArea});
            save(fnOut,'D','tbShort','tbMedium','tbLong','-v7.3');
            fprintf('Done.\n')
            clear D
        end
    end
end
clearvars -except pat2 Nbs AssemblyChoice color_ Delays_ bw tbShort tbAll tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin %% Batch Process assemblies 
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
    
    %%  DelayEnd Decoding for units
     for iDelay = 1:length(Delays_)
         eval(sprintf('LeftTrials = [t.%s.DelayEnd_LeftCorrect];',Delays_{iDelay}));
            eval(sprintf('RightTrials = [t.%s.DelayEnd_RightCorrect];',Delays_{iDelay}));
            eval(sprintf('LeftTrialsE = [t.%s.DelayEnd_LeftError''];',Delays_{iDelay}));
            eval(sprintf('RightTrialsE = [t.%s.DelayEnd_RightError''];',Delays_{iDelay}));
        for iArea = 1:length(Ass.FSCsel)
            
            FSC = Ass.FSCsel{iArea};
            if ~isempty(FSC)
                
                Ltrials = [];nL = 0;
                for iTrial =1:size(LeftTrials,1)
                    try
                        tlims_  = LeftTrials(iTrial,1)/1e6+tlimsAll(1) + shift;
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];
                        Ltrials = [Ltrials;FSC(tlims_,:)];
                        nL=nL+1;
                    end
                end
                
                Rtrials = [];nR = 0;
                for iTrial =1:size(RightTrials,1)
                    try
                        tlims_  = RightTrials(iTrial,1)/1e6+tlimsAll(1) + shift;
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];

                        Rtrials = [Rtrials;FSC(tlims_,:)];
                        nR=nR+1;
                    end
                end
                
                LtrialsE = [];nLe = 0;
                for iTrial =1:size(LeftTrialsE,1)
                    try
                        tlims_  = LeftTrialsE(iTrial,1)/1e6+tlimsAll(1) + shift;
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];

                        LtrialsE = [LtrialsE;FSC(tlims_,:)];
                        nLe=nLe+1;
                    end
                end
                
                RtrialsE = [];nRe = 0;
                for iTrial =1:size(RightTrialsE,1)
                    try
                        tlims_  = RightTrialsE(iTrial,1)/1e6+tlimsAll(1) + shift;
                        tlims_ = closest(Ass.Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];

                        RtrialsE = [RtrialsE;FSC(tlims_,:)];
                        nRe=nRe+1;
                    end
                end
                
                
                
                % Cross-validation error decoder (correct trials)
%                 if nL>2 && nR>2
                    
                    FSC_ = [Ltrials;Rtrials];
                    evt0 = [ones(nL,1);2*ones(nR,1)];
                    Ltr = 2*length(tbAll);
                    nAss_ = cellfun(@(x)size(x,2),Ass.FSCsel);
                    % remove zero entries
                    idx = nAss_>0;
                    nAss_(~idx) = min(nAss_(idx));
                    
                    
                    % Prediction error for cross-validation decoder
                    if size(Ltrials,2)==min(nAss_), nrep=1; else nrep=10; end
                    %                  nrep=1;
                    for i=1:nrep
                        rs{i} = randsample(size(Ltrials,2),min(nAss_));
                        D_.Tone{iArea}.CVE(i,:) = DecodeCVE_mex(FSC_(:,rs{i}),evt0,0.05);
                    end
                    D_.Tone{iArea}.CVE = nanmean(D_.Tone{iArea}.CVE,1);
                   
                    clear rs nAss idx evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
%                 end
                
                
                
                
                
                % Cross-validation error decoder (error trials)
%                 if nLe>2 && nRe>2
                    
                    FSC_ = [Ltrials;Rtrials;LtrialsE;RtrialsE];
                    evt0 = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
                    Ltr = 2*length(tbAll);
                    nAss_ = cellfun(@(x)size(x,2),Ass.FSCsel);
                    % remove zero entries
                    idx = nAss_>0;
                    nAss_(~idx) = min(nAss_(idx));
                    
                    
                    % Prediction error for cross-validation decoder
                    if size(Ltrials,2)==min(nAss_), nrep=1; else nrep=10; end
                    %                  nrep=1;
                    for i=1:nrep
                        rs{i} = randsample(size(Ltrials,2),min(nAss_));
                        D_.Tone{iArea}.CVE_err(i,:)=DecodePredErr(FSC_(:,rs{i}),evt0,0.05);
                    end
                    D_.Tone{iArea}.CVE_err = nanmean(D_.Tone{iArea}.CVE_err,1);
                   
                    clear rs nAss idx evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
%                 end
                
            end
        end
        if exist('D_')
            eval(sprintf('D_Ass.%s.Tone =  D_.Tone;',Delays_{iDelay}))
        end
        clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
     end
      
    %% Save
    if exist('D_Ass')

        fnOut = sprintf('%sMixedSelectivity%s%s_DelayEndDecoding_Units.mat',pat,filesep,fname);
        save(fnOut,'D_Ass','tbShort','tbMedium','tbLong','-v7.3');
        fprintf('Done.\n')
        clear D_Ass
    end
end
clearvars -except pat2 Nbs AssemblyChoice color_ Delays_ bw tbShort tbAll tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin 
%% batch import unit data for meta-analysis and plotting
% Import
clear D_
for iArea = 1:2%length(Areas)
 
    for iFile =1:length(fileList)
        
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity\\%s_%s_DelayEndDecoding_Units.mat',pat,fname,Areas{iArea});
        load(fnIn ,'D');
        
        for iDelay = 1:length(Delays_)
            
            % L/R: Continuous decoding of Positional information - correct trials
            try
                eval(sprintf('D_.%s.Tone.CVE{iArea}(:,iFile)      = D.%s.Tone.CVE;',Delays_{iDelay},Delays_{iDelay}))
            catch
                 eval(sprintf('D_.%s.Tone.CVE{iArea}(:,iFile)     = nan(length(tbAll),1);',Delays_{iDelay}))
            end
            try
                eval(sprintf('D_.%s.Tone_err.CVE{iArea}(:,iFile)  = D.%s.Tone_err.CVE;',Delays_{iDelay},Delays_{iDelay}))
            catch
                eval(sprintf('D_.%s.Tone_err.CVE{iArea}(:,iFile)  = nan(length(tbAll),1);',Delays_{iDelay}))
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
        fnIn = sprintf('%s\\MixedSelectivity\\%s_DelayEndDecoding_Units.mat',pat,fname);
        load(fnIn ,'D_Ass');
        
        for iArea = 1:3
            for iDelay = 1:length(Delays_)
                
                % L/R: Continuous decoding of Positional information - correct trials
                flag = 0;
                if eval(sprintf('isfield(D_Ass.%s,''Tone'')',Delays_{iDelay})) 
                    if eval(sprintf('~isempty(D_Ass.%s.Tone{iArea})',Delays_{iDelay}))
                        flag = 1;
                        eval(sprintf('D_temp = D_Ass.%s.Tone{iArea};',Delays_{iDelay}))
                     
                        
                         try
                            eval(sprintf('D_Ass_.%s.Tone.CVE{iArea}(:,iFile)      = D_temp.CVE;',Delays_{iDelay}))
                         catch                         
                             eval(sprintf('D_Ass_.%s.Tone.CVE{iArea}(:,iFile)     = nan(length(tbAll),1);',Delays_{iDelay}))
                         end
                         try
                            eval(sprintf('D_Ass_.%s.Tone.CVE_err{iArea}(:,iFile)      = D_temp.CVE_err;',Delays_{iDelay}))
                         catch                         
                             eval(sprintf('D_Ass_.%s.Tone.CVE_err{iArea}(:,iFile)     = nan(length(tbAll),1);',Delays_{iDelay}))
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
    for iArea = 1:2
        eval(sprintf('D_temp = D_.%s.Tone.CVE{iArea};',Delays_{iDelay}))
        D_temp(:,sum(D_temp)==0)=NaN;
        D_temp(:,isnan(D_temp(1,:)))=[];
        eval(sprintf('D_.%s.Tone.CVE{iArea} = D_temp;',Delays_{iDelay}))
    end
    for iArea = 1:2
        eval(sprintf('D_temp = D_.%s.Tone_err.CVE{iArea};',Delays_{iDelay}))
        D_temp(:,sum(D_temp)==0)=NaN;
        D_temp(:,isnan(D_temp(1,:)))=[];
        eval(sprintf('D_.%s.Tone_err.CVE{iArea} = D_temp;',Delays_{iDelay}))
     end
    
    for iArea = 1:3
        eval(sprintf('D_temp = D_Ass_.%s.Tone.CVE{iArea};',Delays_{iDelay}))
        D_temp(:,sum(D_temp)==0)=NaN;
        D_temp(:,isnan(D_temp(1,:)))=[];
        eval(sprintf('D_Ass_.%s.Tone.CVE{iArea} = D_temp;',Delays_{iDelay}))
        
        eval(sprintf('D_temp = D_Ass_.%s.Tone.CVE_err{iArea};',Delays_{iDelay}))
        D_temp(:,sum(D_temp)==0)=NaN;
        D_temp(:,isnan(D_temp(1,:)))=[];
        eval(sprintf('D_Ass_.%s.Tone.CVE_err{iArea} = D_temp;',Delays_{iDelay}))
    end
end
 
%% Plot L/R decoding on correct and errors - CVE - units
alpha_ = 0.6;
tb = [0:(length(tbAll)-1)]*bw + tlimsAll(1);
for iArea =1:2%length(Areas)
    figure('name',['L/R decoding ',Areas{iArea}]);
    for iDelay = 1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('corr = (1-D_.%s.Tone.CVE{iArea})*100;',Delays_{iDelay}));
        eval(sprintf('err = (1-D_.%s.Tone_err.CVE{iArea})*100;',Delays_{iDelay}));
%         corr = smooth2a(corr,10,0);
%         err = smooth2a(err,10,0);
        
        m = nanmean(corr,2);
        e = nansem(corr,2);
        
        ciplot(m+e,m-e,tb,'k',alpha_)
        m = nanmean(err,2);
        e = nansem(err,2);
        ciplot(m+e,m-e,tb,'r',alpha_)
        
        [sig,crit] = permtest2vec(corr,err,500,0.05);
        
        a = nan(size(sig));a(sig)=1;
        plot(tb,repmat(105,length(tb),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)
        plot(tb,a*105,'color',[0 0 0],'LineWidth',3)   
        
      
        plot([0 0],[0 100],':k')
        title(Delays_{iDelay})
        axis([-5 5 0 110])
        if iDelay==1
            ylabel('% Correct decoding')
        end
        xlabel('Time (s)')
    end
    legend({'Correct','Errors'})
end
%% Plot L/R decoding on correct and errors - CVE - Ass
alpha_ = 0.6;
tb = [0:(length(tbAll)-1)]*bw + tlimsAll(1);
for iArea =1:length(Areas)
    figure('name',['N/P decoding ',Areas{iArea}]);
    for iDelay = 1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('corr = (1-D_Ass_.%s.Tone.CVE{iArea})*100;',Delays_{iDelay}));
        eval(sprintf('err = (1-D_Ass_.%s.Tone.CVE_err{iArea})*100;',Delays_{iDelay}));
%         corr = smooth2a(corr,10,0);
%         err = smooth2a(err,10,0);
        
        m = nanmean(corr,2);
        e = nansem(corr,2);
        
        ciplot(m+e,m-e,tb,'k',alpha_)
        m = nanmean(err,2);
        e = nansem(err,2);
        ciplot(m+e,m-e,tb,'r',alpha_)
        
        [sig,crit] = permtest2vec(corr,err,500,0.05);
        
        a = nan(size(sig));a(sig)=1;
        plot(tb,repmat(105,length(tb),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)
        plot(tb,a*105,'color',[0 0 0],'LineWidth',3)   
        
      
        plot([0 0],[0 100],':k')
        title(Delays_{iDelay})
        axis([-5 5 0 110])
        if iDelay==1
            ylabel('% Correct decoding')
        end
        xlabel('Time (s)')
    end
    legend({'Correct','Errors'})
end



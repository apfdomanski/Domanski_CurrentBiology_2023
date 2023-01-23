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
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/'
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
        %% Prediction error for units
        LeftTrials = [];
        RightTrials = [];
        LeftTrialsE = [];
        RightTrialsE = [];
        for iDelay = 1:length(Delays_)
            eval(sprintf('LeftTrials = [LeftTrials;[t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect]];',Delays_{iDelay},Delays_{iDelay}));
            eval(sprintf('RightTrials = [RightTrials;[t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect]];',Delays_{iDelay},Delays_{iDelay}));
            eval(sprintf('LeftTrialsE = [LeftTrialsE;[t.%s.SamplePress_LeftError'',t.%s.ChoicePress_LeftError'']];',Delays_{iDelay},Delays_{iDelay}));
            eval(sprintf('RightTrialsE = [RightTrialsE;[t.%s.SamplePress_RightError'',t.%s.ChoicePress_RightError'']];',Delays_{iDelay},Delays_{iDelay}));
        end
        
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
        
        LtrialsE = [];nLe = 0;
        for iTrial =1:size(LeftTrialsE,1)
            try
                tlims_  = [LeftTrialsE(iTrial,1)/1e6;LeftTrialsE(iTrial,2)/1e6]+tlimsAll(1) + shift;
                tlims_ = closest(Tmtx,tlims_);
                tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                
                LtrialsE = [LtrialsE;iFR_(tlims_,:)];
                nLe=nLe+1;
            end
        end
        
        RtrialsE = [];nRe = 0;
        for iTrial =1:size(RightTrialsE,1)
            try
                tlims_  = [RightTrialsE(iTrial,1)/1e6;RightTrialsE(iTrial,2)/1e6]+tlimsAll(1) + shift;
                tlims_ = closest(Tmtx,tlims_);
                tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                
                RtrialsE = [RtrialsE;iFR_(tlims_,:)];
                nRe=nRe+1;
            end
        end
        
        % (1) Correct trials first
        FR = [Ltrials;Rtrials];
        evt0 = [ones(nL,1);2*ones(nR,1)];
        Ltr = 2*length(tbAll);
        nu = [length(PFCcells),length(HPcells)];
        
        % Cross-validation error decoder
        if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
        for i=1:nrep
            rs{i} = randsample(size(Ltrials,2),min(nu));
            D_.LR{iArea}.CVE(i,:)=DecodeCVE(FR(:,rs{i}),evt0,0.05);
        end
        D_.LR{iArea}.CVE = nanmean( D_.LR{iArea}.CVE,1);
        
        % (2) Now errors
        FR = [Ltrials;Rtrials;LtrialsE;RtrialsE];
        evt0 = [ones(nL,1);2*ones(nR,1);10*ones(nLe,1);20*ones(nRe,1)];
        Ltr = 2*length(tbAll);
        nu = [length(PFCcells),length(HPcells)];
        % Cross-validation prediction decoder
        if size(Ltrials,2)==min(nu), nrep=1; else nrep=10; end
        for i=1:nrep
            rs{i} = randsample(size(Ltrials,2),min(nu));
            D_.LR_err{iArea}.pe0(i,:)=DecodePredErr(FR(:,rs{i}),evt0,0.05);
        end
        D_.LR_err{iArea}.pe0 = nanmean(D_.LR_err{iArea}.pe0,1);
        
        
        clear rs nu evt0 Ltr FR cve0 b k cve0 i CVEbs cves nrep tr
        
        if exist('D_')
            D.All.LR =  D_.LR;
            try
                D.All.LR_err =  D_.LR_err;
            catch
                disp('no errors!')
            end
        end
        clear Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
        
        %% Save results
        if exist('D')
            
            fnOut = sprintf('%sMixedSelectivity%s%s_%s_PredictionError_UnitsRedux.mat',pat,filesep,fname,Areas{iArea});
            save(fnOut,'D','tbShort','tbMedium','tbLong','-v7.3');
            fprintf('Done.\n')
            clear D
        end
        
    end
end
clearvars -except Nbs AssemblyChoice color_ Delays_ bw tbShort tbAll tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA tlimsAll tbShort tbMedium tbLong  shift normaliseFscores normWin %% Batch Process assemblies 
%% batch import unit data 

clear D_
for iArea = 1:2%length(Areas)
 
    for iFile =1:length(fileList)
        
        fname=strtok(fileList(iFile).name,'_');
        fnIn =fullfile(pat,'MixedSelectivity',sprintf('%s_%s_PredictionError_UnitsRedux.mat',fname,Areas{iArea}));
        load(fnIn ,'D');
        
        
            % L/R: Continuous decoding of Positional information - correct trials
            
            D_.LR.CVE{iArea}(:,iFile)      = D.All.LR{iArea}.CVE;
            D_.LR.pe0{iArea}(:,iFile)      = D.All.LR_err{iArea}.pe0;
            
%         clear D
    end
    
end
% clear D Stability D_temp D_tempL D_tempR B bp 

%% Plot L/R decoding on correct and errors - CVE - units
alpha_ = 0.2;
tb = [0:(length(tbAll)*2-1)]*bw;
col_ = {[0 0 1],[1 0 0]};
figure;
for iArea =1:2%length(Areas)
    subplot(1,2,iArea); 
    hold on 
    
        plot([5 5],[0 100],'color',[0 1 0 0.6],'LineWidth',1.5)
        plot([15 15],[0 100],'color',[1 0 0 0.6],'LineWidth',1.5)
        plot([0 20],[50 50],':k')
        corr = 100*(1-D_.LR.CVE{iArea});
        err  = 100*(1-D_.LR.pe0{iArea});
%         corr = full(smooth2a(corr,5,0));
%         err  = full(smooth2a(err,5,0));
        
        m = nanmean(corr,2);
        e = nansem(corr,2);
        
        ciplot(m+e,m-e,tb,col_{iArea},0.6)
        plot(tb,m+e,'LineWidth',1.5,'LineStyle','-','color',col_{iArea}) 
        plot(tb,m-e,'LineWidth',1.5,'LineStyle','-','color',col_{iArea}) 
%         plot(tb,m,'LineWidth',1.5,'LineStyle','-','color',col_{iArea}) 
        
        m = nanmean(err,2);
        e = nansem(err,2);
        ciplot(m+e,m-e,tb,col_{iArea},alpha_)
        plot(tb,m+e,'LineWidth',1.5,'LineStyle','-','color',col_{iArea}) 
        plot(tb,m-e,'LineWidth',1.5,'LineStyle','-','color',col_{iArea}) 
%         plot(tb,m,'LineWidth',1.5,'LineStyle','-','color',col_{iArea}) 
        
        [sig,crit] = permtest2vec(full(smooth2a(corr,0,0)),full(smooth2a(err,0,0)),1000,0.05);
        
        a = nan(size(sig));a(sig)=1;
        plot(tb,repmat(105,length(tb),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)
        plot(tb,a*105,'color',[0 0 0],'LineWidth',3)   
        
        plot([10 10],[30 106],'color','w','LineWidth',2)
        
        
        axis([0 20 0 106])
        if iArea==1
            ylabel('% Correct decoding')
        end
        xlabel('Time (s)')
end
%     legend({'Correct','Errors'})

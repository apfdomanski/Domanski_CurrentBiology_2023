% function DecodingDuringDelayIndividualUnits(iFile)
% if isstr(iFile)
%     iFile = str2num(iFile);
% end
%% %%%%%% PREAMBLE %%%%%%
clear 

Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Del = [4 8 16];
Target = 'LONG';
RunErrors = false;
tlimsShort=[0 4];
tlimsMedium=[0 8];
tlimsLong=[0 16];
tlimsShort=[-5 9];
tlimsMedium=[-5 13];
tlimsLong=[-5 21];
plotOnline = false;
bw    = 0.05;
Nbs = 500;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);

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

%% batch process performance

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
%% batch process units...
DelPer=[];
for iFile =1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
    %% Batch process units
    for iArea = 1:2%length(Areas)
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
        %% Delay period firing rates for units
      
        for iDelay =1:length(Delays_)
            
            eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftCorrect];',Delays_{iDelay}));
            eval(sprintf('RightTrials = [t.%s.SamplePress_RightCorrect];',Delays_{iDelay}));
            eval(sprintf('LeftTrialsE = [t.%s.SamplePress_LeftError''];',Delays_{iDelay}));
            eval(sprintf('RightTrialsE = [t.%s.SamplePress_RightError''];',Delays_{iDelay}));
            eval(sprintf('tlims_X = tlims%s;',Delays_{iDelay}));
            length_ = sum(abs(tlims_X))/bw;
            eval(sprintf('tb_ = tb%s;',Delays_{iDelay}))
            
            Ltrials = [];nL = 0; Ltrials_={}; 
            LtrialsMean =zeros(length(tb_),size(iFR_,2));
            LtrialsVar  =nan(length(tb_),size(iFR_,2));
            for iTrial =1:size(LeftTrials,1)
                try
                    tlims_  = LeftTrials(iTrial,1)/1e6+tlims_X(1);
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    Ltrials = [Ltrials;iFR_(tlims_,:)];
                    Ltrials_{iTrial} = iFR_(tlims_,:);
                    LtrialsVar(:,:,nL+1) = iFR_(tlims_,:);
                    LtrialsMean = LtrialsMean + Ltrials_{iTrial};
                    nL=nL+1;
                end
            end
            LtrialsMean=LtrialsMean./nL;
            LtrialsVar = nanvar(LtrialsVar,[],3);
            
            Rtrials = [];nR = 0;  Rtrials_={}; 
            RtrialsMean =zeros(length(tb_),size(iFR_,2));
            RtrialsVar  =nan(length(tb_),size(iFR_,2));
            for iTrial =1:size(RightTrials,1)
                try
                    tlims_  = RightTrials(iTrial,1)/1e6+tlims_X(1);
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    Rtrials = [Rtrials;iFR_(tlims_,:)];
                    Rtrials_{iTrial} = iFR_(tlims_,:);
                    RtrialsVar(:,:,nR+1) = iFR_(tlims_,:);
                    RtrialsMean = RtrialsMean + Rtrials_{iTrial};
                    nR=nR+1;
                end
            end
            RtrialsMean=RtrialsMean./nR;
            RtrialsVar = nanvar(RtrialsVar,[],3);                        
            
            LtrialsE = [];nLe = 0; LtrialsE_={}; 
            LtrialsMeanE =zeros(length(tb_),size(iFR_,2));
            LtrialsVarE  =nan(length(tb_),size(iFR_,2));
            for iTrial =1:size(LeftTrialsE,1)
                try
                    tlims_  = LeftTrialsE(iTrial,1)/1e6+tlims_X(1);
                    tlims_ = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    LtrialsE = [LtrialsE;iFR_(tlims_,:)];
                    LtrialsE_{iTrial} = iFR_(tlims_,:);
                    LtrialsVarE(:,:,nLe+1) = iFR_(tlims_,:);
                    LtrialsMeanE = LtrialsMeanE + LtrialsE_{iTrial};
                    nLe=nLe+1;
                end
            end
            LtrialsMeanE=LtrialsMeanE./nLe;
            LtrialsVarE = nanvar(LtrialsVarE,[],3);
            
            RtrialsE = [];nRe = 0; RtrialsE_={}; 
            RtrialsMeanE = zeros(length(tb_),size(iFR_,2));
            RtrialsVarE  = nan(length(tb_),size(iFR_,2));
            for iTrial =1:size(RightTrialsE,1)
                try
                    tlims_  = RightTrialsE(iTrial,1)/1e6+tlims_X(1) ;
                    tlims_ = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    RtrialsE = [RtrialsE;iFR_(tlims_,:)];
                    RtrialsE_{iTrial} = iFR_(tlims_,:);
                    RtrialsVarE(:,:,nRe+1) = iFR_(tlims_,:);
                    RtrialsMeanE = RtrialsMeanE + RtrialsE_{iTrial};
                    nRe=nRe+1;
                end
            end
            RtrialsMeanE=RtrialsMeanE./nRe;
            RtrialsVarE = nanvar(RtrialsVarE,[],3);
            
            LtrialsFano  = LtrialsVar./LtrialsMean;
            RtrialsFano  = RtrialsVar./RtrialsMean;
            LtrialsFanoE = LtrialsVarE./LtrialsMeanE;
            RtrialsFanoE = RtrialsVarE./RtrialsMeanE;
            
            DelPer.LtrialsMean{iArea}{iDelay}{iFile}  = LtrialsMean;
            DelPer.RtrialsMean{iArea}{iDelay}{iFile}  = RtrialsMean;
            DelPer.LtrialsMeanE{iArea}{iDelay}{iFile} = LtrialsMeanE;
            DelPer.RtrialsMeanE{iArea}{iDelay}{iFile} = RtrialsMeanE;
            
            DelPer.LtrialsFano{iArea}{iDelay}{iFile}  = LtrialsFano;
            DelPer.RtrialsFano{iArea}{iDelay}{iFile}  = RtrialsFano;
            DelPer.LtrialsFanoE{iArea}{iDelay}{iFile} = LtrialsFanoE;
            DelPer.RtrialsFanoE{iArea}{iDelay}{iFile} = RtrialsFanoE;
            
            DelPer.AlltrialsMean{iArea}{iDelay}{iFile}  = nansum(cat(3,LtrialsMean, RtrialsMean),3)./(~isnan(nansum(nansum(LtrialsMean))) + ~isnan(nansum(nansum(RtrialsMean))));
            DelPer.AlltrialsMeanE{iArea}{iDelay}{iFile} = nansum(cat(3,LtrialsMeanE, RtrialsMeanE),3)./(~isnan(nansum(nansum(LtrialsMeanE))) + ~isnan(nansum(nansum(RtrialsMeanE))));
            
            DelPer.AlltrialsFano{iArea}{iDelay}{iFile}  = nansum(cat(3,LtrialsFano, RtrialsFano),3)./(~isnan(nansum(nansum(LtrialsFano))) + ~isnan(nansum(nansum(RtrialsFano))));
            DelPer.AlltrialsFanoE{iArea}{iDelay}{iFile} = nansum(cat(3,LtrialsFanoE, RtrialsFanoE),3)./(~isnan(nansum(nansum(LtrialsFanoE))) + ~isnan(nansum(nansum(RtrialsFanoE))));
            
            best_ = (nanmax(LtrialsMean,2) < nanmax(RtrialsMean,2))+1;
            X     = cat(3,LtrialsMean,RtrialsMean);
            XErr  = cat(3,LtrialsMeanE,RtrialsMeanE);
            
            Y     = cat(3,LtrialsFano,RtrialsFano);
            YErr  = cat(3,LtrialsFanoE,RtrialsFanoE);
            
            for iUnit =1:length(best_) 
                DelPer.CorrTrialsPrefMean{iArea}{iDelay}{iFile} = X(:,:,best_(iUnit));
                DelPer.ErrTrialsPrefMean{iArea}{iDelay}{iFile}  = XErr(:,:,best_(iUnit));
                
                DelPer.CorrTrialsPrefFano{iArea}{iDelay}{iFile} = Y(:,:,best_(iUnit));
                DelPer.ErrTrialsPrefFano{iArea}{iDelay}{iFile}  = YErr(:,:,best_(iUnit));
            end
            
        end
        
    end
end
%% Plot Mean
for iArea=1:2
    figure
    for iDelay=1:3
        subplot(3,1,iDelay); hold on
        eval(sprintf('tb_ = tb%s;',Delays_{iDelay}))
%         eval(sprintf('tlims_X = tlims%s;',Delays_{iDelay}));
%             X = [cell2mat(DelPer.AlltrialsMean{iArea}{iDelay});cell2mat(DelPer.AlltrialsMeanE{iArea}{iDelay});];
            X = [cell2mat(DelPer.CorrTrialsPrefMean{iArea}{iDelay});cell2mat(DelPer.ErrTrialsPrefMean{iArea}{iDelay});];
%             X(:,sum(isnan(X))>0)=[];
            X = zscore(X);
            XErr = X(length(tb_)+1:end,:);
            X = X(1:length(tb_),:);
            ciplot(nanmean(X,2)+nansem(X,2),nanmean(X,2)-nansem(X,2),tb_',color_{iArea},1)
            ciplot(nanmean(XErr,2)+nansem(XErr,2),nanmean(XErr,2)-nansem(XErr,2),tb_',color_{iArea},0.6)
            
            axis([tlimsLong, -2 2])
            if iDelay==1
                plot([15 17],[0 0],'Linewidth',1.5,'color','k')
                plot([15 15],[0 1],'Linewidth',1.5,'color','k')
            end
            plot([0,0],[-2 -1],'Linewidth',1.5,'color','g')
            plot(Del(iDelay).*[1,1],[-2 -1],'Linewidth',1.5,'color','k')
                
    end
end
%% Plot Fano
for iArea=1:2
    figure
    for iDelay=1:3
        subplot(3,1,iDelay); hold on
        eval(sprintf('tb_ = tb%s;',Delays_{iDelay}))
%         eval(sprintf('tlims_X = tlims%s;',Delays_{iDelay}));
            X = [cell2mat(DelPer.AlltrialsFano{iArea}{iDelay});cell2mat(DelPer.AlltrialsFanoE{iArea}{iDelay});];
%             X = [cell2mat(DelPer.CorrTrialsPrefFano{iArea}{iDelay});cell2mat(DelPer.ErrTrialsPrefFano{iArea}{iDelay});];
%             X(:,sum(isnan(X))>0)=[];
%             X = zscore(X);
            XErr = X(length(tb_)+1:end,:);
            X = X(1:length(tb_),:);
            ciplot(nanmean(X,2)+nansem(X,2),nanmean(X,2)-nansem(X,2),tb_',color_{iArea},1)
            ciplot(nanmean(XErr,2)+nansem(XErr,2),nanmean(XErr,2)-nansem(XErr,2),tb_',color_{iArea},0.6)
            
            axis([tlimsLong, 0 5])
            if iDelay==1
                plot([15 17],[0 0],'Linewidth',1.5,'color','k')
                plot([15 15],[0 1],'Linewidth',1.5,'color','k')
            end
            plot([0,0],[-2 -1],'Linewidth',1.5,'color','g')
            plot(Del(iDelay).*[1,1],[-2 -1],'Linewidth',1.5,'color','k')
                
    end
end
            
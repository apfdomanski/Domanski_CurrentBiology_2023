%% Preamble
clear
warning('off')

Target = 'Long';
UseResampledResults = true;
UseTaskOnlyFiringRates = true;
ProcessUnits = true;
ProcessAssemblies = true;
useWholeTaskPeriod = true;
tlimsAll = [-5 5];

shift = 0;
bw=0.05;
tbAll = tlimsAll(1):bw:tlimsAll(2);

Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};


%%%%%%%%%
%Parameters for random draws
noDrawnTrials = 5;
noTrialsDraws = 100;
nDrawnUnits = 2;
noUnitsDraws = 100;

nDrawnAss = 2;
noAssDraws = 25;
%%%%%%%%%

if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
    plotOnline = true;
else
    plotOnline = false;
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB/')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'allTimestamps' filesep '*' upper(Target) '*.mat']);
% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end

if UseResampledResults
    eval(sprintf('fileList_%s = dir([pat ''MixedSelectivity%sTrialNosMatched%s*%s*PFC*Units*.mat'']);',upper(Target),filesep,filesep,upper(Target)))
    eval(sprintf('AssfileList_%s = dir([pat ''MixedSelectivity%sTrialNosMatched%s*%s*PFC*Ass*.mat'']);',upper(Target),filesep,filesep,upper(Target)))
else
    eval(sprintf('fileList_%s = dir([pat ''MixedSelectivity%s*%s*PFC*Units*.mat'']);',upper(Target),filesep,upper(Target)))
    eval(sprintf('AssfileList_%s = dir([pat ''MixedSelectivity%s*%s*PFC*Ass*.mat'']);',upper(Target),filesep,upper(Target)))
end
AssfileList = intersect(strtok({AssfileList_LONG.name},'_')',strtok({fileList.name},'_')');

fileList = intersect(strtok({fileList_LONG.name},'_')',strtok({fileList.name},'_')');
clear fileList_LONG reject_list idx fnames_ name_flag
tuning = {'Untuned','Simple (Context)','Simple (Location)','Simple','Linear Mixed','Non-linear Mixed','Any Tuning'};
%% Process Single units
clear Batch
if ProcessUnits
    %% Batch import unit results
    
    for iFile =1:length(fileList)
        fname=fileList{iFile};
        
        load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
        load(sprintf('%s%s.mat',pat,fname));
        
        for iArea = 1:2
            %% Get the files
            fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
            % Get firing rates
            if UseTaskOnlyFiringRates
                load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
            else
                load(sprintf('%sKDE_bins%s%s_%s_iFR50.mat',pat,filesep,fname,Areas{iArea}));
            end
            
            iFR_ = iFR;
            % Get Mixed selectivity results
            if UseResampledResults
                fnIn = sprintf('%sMixedSelectivity%sTrialNosMatched%s%s_%s_MixedSelectivity_Units.mat',pat,filesep,filesep,fname,Areas{iArea});
            else
                fnIn = sprintf('%sMixedSelectivity%s%s_%s_MixedSelectivity_Units.mat',pat,filesep,fname,Areas{iArea});
            end
            MixedSelectivity = load(fnIn ,'D');
            Batch.p{iArea}{iFile} =  eval(sprintf('MixedSelectivity.D.%s.ANOVA2.p_cor',Target));
            Batch.F{iArea}{iFile} =  eval(sprintf('MixedSelectivity.D.%s.ANOVA2.F_cor',Target));
            %% Extract the firing rate cutouts
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
                
                Batch.Ltrials{iArea}{iFile} = Ltrials_;
                Batch.Rtrials{iArea}{iFile} = Rtrials_;
                
                Batch.nR(iFile)=nR;
                Batch.nL(iFile)=nL;
                
                clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
            end
        end
    end
    clear s ERRORtrangeleft_choice ERRORtrangeleft_sample ERRORtrangeright_choice ERRORtrangeright_sample
    clear CvBW avgFR E evt0 FR fnIn HPcells HPinter HPtonic iFR iFR_ inputData iUnit Ltr Ltrials_ MISE Nbs
    clear normWin normaliseFscores nu PFCcells PFCinter PFCtonic rsL rsR Rtrials_ t Tmtx
    clear trangeleft_choice trangeleft_sample trangeright_choice trangeright_sample MixedSelectivity
    %% Analyse decoding contributions from different subsets of encoding units - Drawing from both trials and neurons
    MaxTrials = min([Batch.nL Batch.nR]);
    % Draw trials for all neurons first, then iteratively draw sub-pools of
    % neurons based on their classification type
    for iArea = 1:2
        [TrLen,noUnits] = size(Batch.Rtrials{iArea}{1}{1});
        
        for iClass = 1:length(tuning)
            %     D.Ft2{iClass} = zeros(noTrialsDraws,TrLen);
            %     D.Rt2{iClass} = zeros(noTrialsDraws,TrLen);
            D.CVE{iClass} = zeros(noTrialsDraws,TrLen);
        end
        
        for iDrawTrials = 1:noTrialsDraws
            clear FR p_
            p_=nan(1,3);
            for iFile =1:length(fileList)
                
                noTrialsL = length(Batch.Ltrials{iArea}{iFile});
                noTrialsR = length(Batch.Rtrials{iArea}{iFile});
                
                [TrLen,noUnits] = size(Batch.Rtrials{iArea}{iFile}{1});
                
                if min([noTrialsL noTrialsR])>=noDrawnTrials
                    
                    IdxL = randsample(noTrialsL,noDrawnTrials);
                    IdxR = randsample(noTrialsR,noDrawnTrials);
                    
                    FR{iFile,1} = cell2mat(cellfun(@transpose,[{Batch.Ltrials{iArea}{iFile}{IdxL}},...
                        {Batch.Rtrials{iArea}{iFile}{IdxR}}],'UniformOutput',false));
                    p_ = [p_;Batch.p{iArea}{iFile}<0.05];
                    
                else
                    
                    FR{iFile,1} = nan(noUnits,2*TrLen*noDrawnTrials);
                    p_ = [p_;nan(noUnits,3)];
                    
                end
                
            end
            p_(1,:)=[];
            
            FR   = cell2mat(FR);
            evt0 = [ones(1,noDrawnTrials),2*ones(1,noDrawnTrials)];
            
            idx = isnan(FR(:,1));FR(idx,:)=[]; p_(idx,:)=[];
            
            %Fractions of units with different tunings:
            UnitIDs = ...
                {find(sum(p_,2)==0),...               % No tuning
                find(sum(p_ == [1 0 0],2)==3),...     % CS for Context
                find(sum(p_ == [0 1 0],2)==3),...     % CS for Location
                sort([find(sum(p_ == [1 0 0],2)==3);find(sum(p_ == [0 1 0],2)==3)]),... % Any CS
                find(sum(p_ == [1 1 0],2)==3),...     % LMS for Context x Location
                find(sum(p_(:,3) == 1,2)),...         % NMS for Context x Location
                find(sum(p_,2)>0)};                   % Any tuning
            
            
            for iClass = 1:length(UnitIDs)
                [iDrawTrials,iClass]
                id  = zeros(nDrawnUnits,1);
                
                Ft2 = zeros(nDrawnUnits,TrLen);
                Rt2 = zeros(nDrawnUnits,TrLen);
                CVE = zeros(nDrawnUnits,TrLen);
                
                parfor iDrawUnits = 1:noUnitsDraws
                    id = UnitIDs{iClass}(randsample(1:length(UnitIDs{iClass}),nDrawnUnits))
                    [~,~,Ft2(iDrawUnits,:),Rt2(iDrawUnits,:)] = DecodeStats(FR(id,:)',evt0,0.05);
                    [CVE(iDrawUnits,:)] = DecodeCVE(FR(id,:)',evt0,0.05);
                end
                
                D.Ft2{iClass}(iDrawTrials,:) = mean(Ft2);
                D.Rt2{iClass}(iDrawTrials,:) = mean(Rt2);
                D.CVE{iClass}(iDrawTrials,:) = 1-mean(CVE);
                
            end
        end
        if UseResampledResults
            eval(sprintf('fnOut = ''%sMixedSelectivity%sDecodingByUnitClass%sUnitDecodingByType_%s_%s_ResampledMS.mat'';',pat,filesep,filesep,Target,Areas{iArea}));
        else
            eval(sprintf('fnOut = ''%sMixedSelectivity%sDecodingByUnitClass%sUnitDecodingByType_%s_%s.mat'';',pat,filesep,filesep,Target,Areas{iArea}));
        end
        save(fnOut,'D','tbAll')
    end
    %% optional plotting  - units
    if plotOnline
        %% Reimport for plotting
        TrLen = size(D.CVE{1},2);
        %% Plot all classes overlaid
        tb = (1:TrLen)*bw;
        col_ = flipud(copper(length(tuning)));
        figure; hold on
        
        for iClass = 1:length(tuning)
            ciplot(nanmean(D.CVE{iClass})+nansem(D.CVE{iClass}),...
                nanmean(D.CVE{iClass})-nansem(D.CVE{iClass}),...)
                tb,col_(iClass,:),0.8)
            %         plot(mean(D.CVE{iClass}))
        end
        xlabel('Time (s)')
        ylabel('Fraction Correct decoders')
        plot([min(tb) max(tb)], [0.5 0.5],':k')
        legend([tuning,'Chance'])
        axis([min(tb) 30 0 1])
        
        %% Plot all classes separately
        tb = (1:TrLen)*bw;
        col_ = jet(length(tuning));
        figure; hold on
        
        for iClass = 1:length(tuning)
            subplot(1,length(tuning),iClass); hold on
            ciplot(nanmean(D.CVE{iClass})+nansem(D.CVE{iClass}),...
                nanmean(D.CVE{iClass})-nansem(D.CVE{iClass}),...)
                tb,'b',0.8)
            plot([min(tb) max(tb)], [0.5 0.5],':k')
            xlabel('Time (s)')
            if iClass==1
                ylabel('Fraction Correct Decoding')
            end
            title(tuning{iClass})
            axis([min(tb) 30 0 1])
        end
        
        %% Plot sorted by tuning type - CVE
        figure
        subplot(1,2,1); hold on
        col_ = [0.3 0.3 0.3;0.6 0.6 0.6];
        ClasstoPlot = [7,1];
        for iClass = 1:length(ClasstoPlot)
            ciplot(nanmean(D.CVE{ClasstoPlot(iClass)})+nansem(D.CVE{ClasstoPlot(iClass)}),...
                nanmean(D.CVE{ClasstoPlot(iClass)})-nansem(D.CVE{ClasstoPlot(iClass)}),...)
                tb,col_(iClass,:),0.8)
        end
        plot([min(tb) max(tb)], [0.5 0.5],':k','LineWidth',1.5)
        xlabel('Time (s)')
        ylabel('Fraction Correct Decoding')
        axis([min(tb) max(tb) 0 1.25])
        plot([5 5],[0.95 1],'Color',[0 1 0 0.3],'LineWidth',4)
        plot([15 15],[0.95 1],'Color',[1 0 0 0.3],'LineWidth',4)
        text(5,1.05,'Sample','HorizontalAlignment','center')
        text(15,1.05,'Choice','HorizontalAlignment','center')
        legend([tuning(ClasstoPlot),'Chance','' ,''],'Location','South'); legend boxoff

        subplot(1,2,2); hold on
        col_ = flipud(copper(3));
        ClasstoPlot = [6 5 4];
        for iClass = 1:length(ClasstoPlot)
            ciplot(nanmean(D.CVE{ClasstoPlot(iClass)})+nansem(D.CVE{ClasstoPlot(iClass)}),...
                nanmean(D.CVE{ClasstoPlot(iClass)})-nansem(D.CVE{ClasstoPlot(iClass)}),...)
                tb,col_(iClass,:),0.8)
        end
        plot([min(tb) max(tb)], [0.5 0.5],':k','LineWidth',1.5)
        axis([min(tb) max(tb) 0 1.25])
        plot([5 5],[0.95 1],'Color',[0 1 0 0.3],'LineWidth',4)
        plot([15 15],[0.95 1],'Color',[1 0 0 0.3],'LineWidth',4)
        text(5,1.05,'Sample','HorizontalAlignment','center')
        text(15,1.05,'Choice','HorizontalAlignment','center')
        legend([tuning(ClasstoPlot),'Chance','' ,''],'Location','South'); legend boxoff
        xlabel('Time (s)')
        %% Plot sorted by tuning type - F-score
        figure
        subplot(1,2,1); hold on
        col_ = [0.3 0.3 0.3;0.6 0.6 0.6];
        ClasstoPlot = [7,1];
        for iClass = 1:length(ClasstoPlot)
            ciplot(nanmean(D.Rt2{ClasstoPlot(iClass)})+nansem(D.Rt2{ClasstoPlot(iClass)}),...
                nanmean(D.Rt2{ClasstoPlot(iClass)})-nansem(D.Rt2{ClasstoPlot(iClass)}),...)
                tb,col_(iClass,:),0.8)
        end
        xlabel('Time (s)')
        ylabel('Fraction Correct Decoding')
        axis([min(tb) max(tb) 0 0.5])
        plot([5 5],[0  0.05],'Color',[0 1 0 0.3],'LineWidth',4)
        plot([15 15],[0  0.05],'Color',[1 0 0 0.3],'LineWidth',4)
        legend([tuning(ClasstoPlot),'' ,''],'Location','North'); legend boxoff
        text(4,0.005,'Sample','Rotation',90)
        text(14,0.005,'Choice','Rotation',90)
        
        subplot(1,2,2); hold on
        col_ = flipud(copper(3));
        ClasstoPlot = [6 5 4];
        for iClass = 1:length(ClasstoPlot)
            ciplot(nanmean(D.Rt2{ClasstoPlot(iClass)})+nansem(D.Rt2{ClasstoPlot(iClass)}),...
                nanmean(D.Rt2{ClasstoPlot(iClass)})-nansem(D.Rt2{ClasstoPlot(iClass)}),...)
                tb,col_(iClass,:),0.8)
        end
        plot([5 5],[0  0.05],'Color',[0 1 0 0.3],'LineWidth',4)
        plot([15 15],[0  0.05],'Color',[1 0 0 0.3],'LineWidth',4)
        legend([tuning(ClasstoPlot),'' ,''],'Location','North'); legend boxoff
        text(4,0.005,'Sample','Rotation',90)
        text(14,0.005,'Choice','Rotation',90)
        xlabel('Time (s)')
        axis([min(tb) max(tb) 0 0.5])
        
    end
end
%% Process cell assemblies
clear Batch
if ProcessAssemblies
    %% Batch import assembly results
    for iFile =1:length(AssfileList)
        fname=AssfileList{iFile};
        
        load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
        load(sprintf('%s%s.mat',pat,fname));
        
        % Load Assemblies
        if useWholeTaskPeriod
            fn = fullfile(pat,'KDE_binsTaskonly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_FSC.mat',fname));
            A  = load(fn);
            fn = fullfile(pat,'KDE_binsTaskonly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_AssemRes2.mat',fname));
            B  = load(fn);
            % recalculate included units
            fn = fullfile(pat,'KDE_binsTaskonly',sprintf('%s_PFC_iFR50_behavOnly.mat',fname));
            usel_out=SelCells(fn,0.1,1e6);
        else
            fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',fname));
            A  = load(fn);
            fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',fname));
            B  = load(fn);
            % recalculate included units
            fn = fullfile(pat,'KDE_bins',sprintf('%s_PFC_iFR50.mat',fname));
            [~,~,B.usel_out]=SelTrialsCellsWholeTrial(fn,10,0.1,1e8);
        end
        
        A.nu(3) =sum(A.nu(1:2));
        B.usel_out{3}= [B.usel_out{1},B.usel_out{2}+max(B.usel_out{1})];
        for iArea = 1:3
            %% Get the files
            fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(AssfileList),fname,Areas{iArea})
            % Get firing rates
            
            FSC =   A.FSCsel{iArea};
            % Get Mixed selectivity results
            if UseResampledResults
                fnIn = sprintf('%sMixedSelectivity%sTrialNosMatched%s%s_%s_MixedSelectivity_Ass.mat',pat,filesep,filesep,fname,Areas{iArea});
            else
                fnIn = sprintf('%sMixedSelectivity%s%s_%s_MixedSelectivity_Ass.mat',pat,filesep,fname,Areas{iArea});
            end
            if exist(fnIn)==2
                MixedSelectivity = load(fnIn ,'D_Ass');
                Batch.p{iArea}{iFile} =  eval(sprintf('MixedSelectivity.D_Ass.%s.ANOVA2.p_cor',Target));
                Batch.F{iArea}{iFile} =  eval(sprintf('MixedSelectivity.D_Ass.%s.ANOVA2.F_cor',Target));
                %% Extract the firing rate cutouts
                for iDelay =1:length(Delays_)
                    
                    eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('RightTrials = [t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}));
                    
                    Ltrials = [];nL = 0;Ltrials_ = {};
                    for iTrial =1:size(LeftTrials,1)
                        try
                            tlims_  = [LeftTrials(iTrial,1)/1e6;LeftTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                            tlims_  = closest(A.Tmtx,tlims_);
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
                            tlims_  = closest(A.Tmtx,tlims_);
                            tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                            Rtrials = [Rtrials;FSC(tlims_,:)];
                            Rtrials_{1,nR+1} = FSC(tlims_,:);
                            nR=nR+1;
                        end
                    end
                    
                    Batch.Ltrials{iArea}{iFile} = Ltrials_;
                    Batch.Rtrials{iArea}{iFile} = Rtrials_;
                    
                    Batch.nR(iFile)=nR;
                    Batch.nL(iFile)=nL;
                    
                    clear D_ Ltrials nL Rtrials nR LeftTrials RightTrials iTrial
                end
                
            end
        end
    end
    clear s ERRORtrangeleft_choice ERRORtrangeleft_sample ERRORtrangeright_choice ERRORtrangeright_sample
    clear CvBW avgFR E evt0 FR fnIn HPcells HPinter HPtonic FSC A  B iFR_ inputData iUnit Ltr Ltrials_ MISE Nbs
    clear normWin normaliseFscores nu PFCcells PFCinter PFCtonic rsL rsR Rtrials_ t Tmtx
    clear trangeleft_choice trangeleft_sample trangeright_choice trangeright_sample MixedSelectivity
    %% Analyse decoding contributions from different subsets of encoding assemblies - Drawing from both trials and assemblies
    MaxTrials = min([Batch.nL Batch.nR]);
    % Draw trials for all neurons first, then iteratively draw sub-pools of
    % neurons based on their classification type
    for iArea = 2:3
        [TrLen,~] = size(Batch.Rtrials{iArea}{1}{1});
        
        for iClass = 1:length(tuning)
            D_Ass.Ft2{iClass} = zeros(noTrialsDraws,TrLen);
            D_Ass.Rt2{iClass} = zeros(noTrialsDraws,TrLen);
            D_Ass.CVE{iClass} = zeros(noTrialsDraws,TrLen);
        end
        
        for iDrawTrials = 1:noTrialsDraws
            clear FSC p_
            p_=nan(1,3);
%             for iFile =1:length(AssfileList)
            for iFile =1:length(Batch.Ltrials{iArea})
                
                noTrialsL = length(Batch.Ltrials{iArea}{iFile});
                noTrialsR = length(Batch.Rtrials{iArea}{iFile});
                
                [TrLen,noUnits] = size(Batch.Rtrials{iArea}{iFile}{1});
                
                if min([noTrialsL noTrialsR])>=noDrawnTrials
                    
                    IdxL = randsample(noTrialsL,noDrawnTrials);
                    IdxR = randsample(noTrialsR,noDrawnTrials);
                    
                    FSC{iFile,1} = cell2mat(cellfun(@transpose,[{Batch.Ltrials{iArea}{iFile}{IdxL}},...
                        {Batch.Rtrials{iArea}{iFile}{IdxR}}],'UniformOutput',false));
                    p_ = [p_;Batch.p{iArea}{iFile}<0.05];
                    
                else
                    FSC{iFile,1} = nan(noUnits,2*TrLen*noDrawnTrials);
                    p_ = [p_;nan(noUnits,3)];
                    
                end
                
            end
            p_(1,:)=[];
            
            FSC   = cell2mat(FSC);
            evt0 = [ones(1,noDrawnTrials),2*ones(1,noDrawnTrials)];
            
            idx = isnan(FSC(:,1));FSC(idx,:)=[]; p_(idx,:)=[];
            
            %Fractions of units with different tunings:
            AssIDs = ...
                {find(sum(p_,2)==0),...               % No tuning
                find(sum(p_ == [1 0 0],2)==3),...     % CS for Context
                find(sum(p_ == [0 1 0],2)==3),...     % CS for Location
                sort([find(sum(p_ == [1 0 0],2)==3);find(sum(p_ == [0 1 0],2)==3)]),... % Any CS
                find(sum(p_ == [1 1 0],2)==3),...     % LMS for Context x Location
                find(sum(p_(:,3) == 1,2)),...         % NMS for Context x Location
                find(sum(p_,2)>0)};                   % Any tuning
            
            for iClass = 1:length(AssIDs)
                [iDrawTrials,iClass]
                id  = zeros(nDrawnAss,1);
                
                Ft2 = nan(nDrawnAss,TrLen);
                Rt2 = nan(nDrawnAss,TrLen);
                CVE = nan(nDrawnAss,TrLen);
                if length(AssIDs{iClass})>=nDrawnAss
                    parfor iDrawAss= 1:noAssDraws
                        id = AssIDs{iClass}(randsample(1:length(AssIDs{iClass}),nDrawnAss))
                        [~,~,Ft2(iDrawAss,:),Rt2(iDrawAss,:)] = DecodeStats(FSC(id,:)',evt0,0.05);
                        [CVE(iDrawAss,:)] = DecodeCVE(FSC(id,:)',evt0,0.05);
                    end
                end
                D_Ass.Ft2{iClass}(iDrawTrials,:) = mean(Ft2);
                D_Ass.Rt2{iClass}(iDrawTrials,:) = mean(Rt2);
                D_Ass.CVE{iClass}(iDrawTrials,:) = 1-mean(CVE);
            end
        end
        if UseResampledResults
            eval(sprintf('fnOut = ''%sMixedSelectivity%sDecodingByUnitClass%sAssemblyDecodingByType_%s_%s_ResampledMS.mat'';',pat,filesep,filesep,Target,Areas{iArea}));
        else
            eval(sprintf('fnOut = ''%sMixedSelectivity%sDecodingByUnitClass%sAssemblyDecodingByType_%s_%s.mat'';',pat,filesep,filesep,Target,Areas{iArea}));
        end
        save(fnOut,'D_Ass','tbAll')
    end
    %% optional plotting  - assemblies
    if plotOnline
        %% Reimport for plotting
        TrLen = size(D_Ass.CVE{1},2);
        %% Plot all classes overlaid
        tb = (1:TrLen)*bw;
        col_ = flipud(copper(length(tuning)));
        figure; hold on
        
        for iClass = 1:length(tuning)
            ciplot(nanmean(D_Ass.CVE{iClass})+nansem(D_Ass.CVE{iClass}),...
                nanmean(D_Ass.CVE{iClass})-nansem(D_Ass.CVE{iClass}),...)
                tb,col_(iClass,:),0.8)
            %         plot(mean(D.CVE{iClass}))
        end
        xlabel('Time (s)')
        ylabel('Fraction Correct decoders')
        plot([min(tb) max(tb)], [0.5 0.5],':k')
        legend([tuning,'Chance'])
        axis([min(tb) 30 0 1])
        
        %% Plot all classes separately
        tb = (1:TrLen)*bw;
        col_ = jet(length(tuning));
        figure; hold on
        
        for iClass = 1:length(tuning)
            subplot(1,length(tuning),iClass); hold on
            ciplot(nanmean(D_Ass.CVE{iClass})+nansem(D_Ass.CVE{iClass}),...
                nanmean(D_Ass.CVE{iClass})-nansem(D_Ass.CVE{iClass}),...)
                tb,'b',0.8)
            plot([min(tb) max(tb)], [0.5 0.5],':k')
            xlabel('Time (s)')
            if iClass==1
                ylabel('Fraction Correct Decoding')
            end
            title(tuning{iClass})
%             axis([min(tb) 30 0.3 0.8])
        end
        
        %% Plot sorted by tuning type - CVE
        figure
        subplot(1,2,1); hold on
        col_ = [0.3 0.3 0.3;0.6 0.6 0.6];
        ClasstoPlot = [7,1];
        for iClass = 1:length(ClasstoPlot)
            ciplot(nanmean(D_Ass.CVE{ClasstoPlot(iClass)})+nansem(D_Ass.CVE{ClasstoPlot(iClass)}),...
                nanmean(D_Ass.CVE{ClasstoPlot(iClass)})-nansem(D_Ass.CVE{ClasstoPlot(iClass)}),...)
                tb,col_(iClass,:),0.8)
        end
        plot([min(tb) max(tb)], [0.5 0.5],':k','LineWidth',1.5)
        xlabel('Time (s)')
        ylabel('Fraction Correct Decoding')
        axis([min(tb) max(tb) 0 1.25])
        plot([5 5],[0.95 1],'Color',[0 1 0 0.3],'LineWidth',4)
        plot([15 15],[0.95 1],'Color',[1 0 0 0.3],'LineWidth',4)
        text(5,1.05,'Sample','HorizontalAlignment','center')
        text(15,1.05,'Choice','HorizontalAlignment','center')
        legend([tuning(ClasstoPlot),'Chance','' ,''],'Location','South'); legend boxoff

        
        subplot(1,2,2); hold on
        col_ = (copper(3));
        ClasstoPlot = [6 5 4];
        for iClass = 1:length(ClasstoPlot)
            ciplot(nanmean(D_Ass.CVE{ClasstoPlot(iClass)})+nansem(D_Ass.CVE{ClasstoPlot(iClass)}),...
                nanmean(D_Ass.CVE{ClasstoPlot(iClass)})-nansem(D_Ass.CVE{ClasstoPlot(iClass)}),...)
                tb,col_(iClass,:),0.8)
        end
        plot([min(tb) max(tb)], [0.5 0.5],':k','LineWidth',1.5)
        plot([5 5],[0.95 1],'Color',[0 1 0 0.3],'LineWidth',4)
        plot([15 15],[0.95 1],'Color',[1 0 0 0.3],'LineWidth',4)
        text(5,1.05,'Sample','HorizontalAlignment','center')
        text(15,1.05,'Choice','HorizontalAlignment','center')
        legend([tuning(ClasstoPlot),'Chance','' ,''],'Location','South'); legend boxoff
        xlabel('Time (s)')
        axis([min(tb) max(tb) 0 1.25])
        %% Plot sorted by tuning type - F-score
        figure
        subplot(1,2,1); hold on
        col_ = [0.3 0.3 0.3;0.6 0.6 0.6];
        ClasstoPlot = [7,1];
        for iClass = 1:length(ClasstoPlot)
            ciplot(nanmean(D_Ass.Rt2{ClasstoPlot(iClass)})+nansem(D_Ass.Rt2{ClasstoPlot(iClass)}),...
                nanmean(D_Ass.Rt2{ClasstoPlot(iClass)})-nansem(D_Ass.Rt2{ClasstoPlot(iClass)}),...)
                tb,col_(iClass,:),0.8)
        end
        xlabel('Time (s)')
        ylabel('Population decoding (F-Score)')
        axis([min(tb) max(tb) 0 1])
        plot([5 5],[0  0.05],'Color',[0 1 0 0.3],'LineWidth',4)
        plot([15 15],[0  0.05],'Color',[1 0 0 0.3],'LineWidth',4)
        legend([tuning(ClasstoPlot),'' ,''],'Location','North'); legend boxoff
        text(4,0.005,'Sample','Rotation',90)
        text(14,0.005,'Choice','Rotation',90)
        
        subplot(1,2,2); hold on
        col_ = flipud(copper(3));
        ClasstoPlot = [6 5 4];
        for iClass = 1:length(ClasstoPlot)
            ciplot(nanmean(D_Ass.Rt2{ClasstoPlot(iClass)})+nansem(D_Ass.Rt2{ClasstoPlot(iClass)}),...
                nanmean(D_Ass.Rt2{ClasstoPlot(iClass)})-nansem(D_Ass.Rt2{ClasstoPlot(iClass)}),...)
                tb,col_(iClass,:),0.9)
        end
        plot([5 5],[0  0.05],'Color',[0 1 0 0.3],'LineWidth',4)
        plot([15 15],[0  0.05],'Color',[1 0 0 0.3],'LineWidth',4)
        legend([tuning(ClasstoPlot),'' ,''],'Location','North'); legend boxoff
        text(4,0.005,'Sample','Rotation',90)
        text(14,0.005,'Choice','Rotation',90)
        xlabel('Time (s)')
        axis([min(tb) max(tb) 0 1])
        
    end
    
end
% function DecodingDuringDelayIndividualUnits(iFile)
% if isstr(iFile)
%     iFile = str2num(iFile);
% end
%% %%%%%% PREAMBLE %%%%%%
clear 

Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};

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

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

TrialsSubsample = 8;
if isinf(TrialsSubsample)
    nReps = 1;
else
    nReps =500;
end
%% batch process performance
for iFile =1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
    %% Get the assembly membership info
    fname=strtok(fileListAss(iFile).name,'_');
    fprintf('Loading run %d/%d %s ...\n',iFile,length(fileList),fname)
    % for assemblies describing whole task period (cue  - reward)
    switch AssemblyChoice
        case 1
            load(sprintf('%s%s%s_iFR50_FSC.mat',pat2,filesep,fname),'units','nu');
            %             load(sprintf('%s%s%s_iFR50_AssemRes2.mat',pat2,filesep,fname),'usel_out');
            usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
        case 2
            load(sprintf('%s%s%s_iFR50_BehavOnly_FSC.mat',pat2,filesep,fname),'units','nu');
            load(sprintf('%s%s%s_iFR50_BehavOnly_AssemRes2.mat',pat2,filesep,fname),'usel_out');
        case 3
            load(sprintf('%s%s%s_iFR50_Task_FSC.mat',pat2,filesep,fname),'units','nu');
            %             load(sprintf('%s%s%s_iFR50_Task_AssemRes2.mat',pat2,filesep,fname),'usel_out');
            usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
            
    end
    noPFC=nu(1);
    Ass.usel_out = usel_out;
    Ass.usel_out{3} = [Ass.usel_out{1},Ass.usel_out{2}+noPFC];
    Ass.nu = nu; Ass.nu(3) = sum(Ass.nu(1:2));
    Ass.units = units;
    for iArea = 1:3
        Ass.LocalMembers{iArea}    = unique(cell2mat(Ass.units{iArea}));
    end
    
    
    Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3}<=Ass.nu(1)));
    Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3}>Ass.nu(1))-Ass.nu(1));
    for iArea = 1:2
        Ass.JointMembers{iArea} = setdiff(unique(cell2mat(Ass.units{iArea})),Ass.LocalMembers{iArea});
        Ass.NonMembers{iArea} = setdiff(1:Ass.nu(iArea),[Ass.LocalMembers{iArea},Ass.JointMembers{iArea}]);
    end
    for iArea = 1:2
        Ass.NonMembers{iArea}   = Ass.usel_out{iArea}(Ass.NonMembers{iArea});
        Ass.LocalMembers{iArea} = Ass.usel_out{iArea}(Ass.LocalMembers{iArea});
        Ass.JointMembers{iArea} = Ass.usel_out{iArea}(Ass.JointMembers{iArea});
    end
    
    clear units usel_out
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
for iFile =1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
    %% Get the assembly membership info
    fname=strtok(fileListAss(iFile).name,'_');
    fprintf('Loading run %d/%d %s ...\n',iFile,length(fileList),fname)
    % for assemblies describing whole task period (cue  - reward)
    switch AssemblyChoice
        case 1
            load(sprintf('%s%s%s_iFR50_FSC.mat',pat2,filesep,fname),'units','nu');
            %             load(sprintf('%s%s%s_iFR50_AssemRes2.mat',pat2,filesep,fname),'usel_out');
            usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
        case 2
            load(sprintf('%s%s%s_iFR50_BehavOnly_FSC.mat',pat2,filesep,fname),'units','nu');
            load(sprintf('%s%s%s_iFR50_BehavOnly_AssemRes2.mat',pat2,filesep,fname),'usel_out');
        case 3
            load(sprintf('%s%s%s_iFR50_Task_FSC.mat',pat2,filesep,fname),'units','nu');
            %             load(sprintf('%s%s%s_iFR50_Task_AssemRes2.mat',pat2,filesep,fname),'usel_out');
            usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
            
    end
    noPFC=nu(1);
    Ass.usel_out = usel_out;
    Ass.usel_out{3} = [Ass.usel_out{1},Ass.usel_out{2}+noPFC];
    Ass.nu = nu; Ass.nu(3) = sum(Ass.nu(1:2));
    Ass.units = units;
    for iArea = 1:3
        Ass.LocalMembers{iArea}    = unique(cell2mat(Ass.units{iArea}));
    end
    
    
    Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3}<=Ass.nu(1)));
    Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3}>Ass.nu(1))-Ass.nu(1));
    for iArea = 1:2
        Ass.JointMembers{iArea} = setdiff(unique(cell2mat(Ass.units{iArea})),Ass.LocalMembers{iArea});
        Ass.NonMembers{iArea} = setdiff(1:Ass.nu(iArea),[Ass.LocalMembers{iArea},Ass.JointMembers{iArea}]);
    end
    for iArea = 1:2
        Ass.NonMembers{iArea}   = Ass.usel_out{iArea}(Ass.NonMembers{iArea});
        Ass.LocalMembers{iArea} = Ass.usel_out{iArea}(Ass.LocalMembers{iArea});
        Ass.JointMembers{iArea} = Ass.usel_out{iArea}(Ass.JointMembers{iArea});
    end
    
    clear units usel_out
    %% Batch process units
    for iArea = 1%:length(Areas)
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
        shuffleCorrect = true;
        for iDelay =1:length(Delays_)
            
            eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftCorrect];',Delays_{iDelay}));
            eval(sprintf('RightTrials = [t.%s.SamplePress_RightCorrect];',Delays_{iDelay}));
            eval(sprintf('LeftTrialsE = [t.%s.SamplePress_LeftError''];',Delays_{iDelay}));
            eval(sprintf('RightTrialsE = [t.%s.SamplePress_RightError''];',Delays_{iDelay}));
            eval(sprintf('tlims_X = tlims%s;',Delays_{iDelay}));
            length_ = sum(abs(tlims_X))/bw;
            eval(sprintf('tb_ = tb%s;',Delays_{iDelay}))
            
            Ltrials = [];nL = 0; Ltrials_={};
            for iTrial =1:size(LeftTrials,1)
                try
                    tlims_  = LeftTrials(iTrial,1)/1e6+tlims_X(1);
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    Ltrials = [Ltrials;iFR_(tlims_,:)];
                    Ltrials_{iTrial} = iFR_(tlims_,:);
                    nL=nL+1;
                end
            end
            Rtrials = [];nR = 0;  Rtrials_={};
            for iTrial =1:size(RightTrials,1)
                try
                    tlims_  = RightTrials(iTrial,1)/1e6+tlims_X(1);
                    tlims_  = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    Rtrials = [Rtrials;iFR_(tlims_,:)];
                    Rtrials_{iTrial} = iFR_(tlims_,:);
                    
                    nR=nR+1;
                end
            end
            LtrialsE = [];nLe = 0; LtrialsE_={};
            for iTrial =1:size(LeftTrialsE,1)
                try
                    tlims_  = LeftTrialsE(iTrial,1)/1e6+tlims_X(1);
                    tlims_ = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    LtrialsE = [LtrialsE;iFR_(tlims_,:)];
                    LtrialsE_{iTrial} = iFR_(tlims_,:);
                    
                    nLe=nLe+1;
                end
            end
            RtrialsE = [];nRe = 0; RtrialsE_={};
            for iTrial =1:size(RightTrialsE,1)
                try
                    tlims_  = RightTrialsE(iTrial,1)/1e6+tlims_X(1) ;
                    tlims_ = closest(Tmtx,tlims_);
                    tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                    RtrialsE = [RtrialsE;iFR_(tlims_,:)];
                    RtrialsE_{iTrial} = iFR_(tlims_,:);
                    nRe=nRe+1;
                end
            end
            
            % run correct trials
            flag_ = true;
            if ~isinf(TrialsSubsample)
                if sum([nL,nR]<=TrialsSubsample)==2
                    flag_ = false;
                end
            end
            if flag_
                FR   = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                reg  = 0.05;
                [PE,PE_CIl,PE_CIh] = DecodePredErrCrossTemporal_CV_GPU(FR,evt0,reg,true);%(:,Ass.JointMembers{iArea})
                D_{iArea}{iDelay}(:,:,iFile)    = PE;
                D_CIh{iArea}{iDelay}(:,:,iFile) = PE_CIh;
                D_CIl{iArea}{iDelay}(:,:,iFile) = PE_CIl;
            else
                D_{iArea}{iDelay}(:,:,iFile)    = nan(length_+1,length_+1);
                D_CIh{iArea}{iDelay}(:,:,iFile) = nan(length_+1,length_+1);
                D_CIl{iArea}{iDelay}(:,:,iFile) = nan(length_+1,length_+1);
            end
            
            % Run error trials
%             flag_ = 1;
%             if ~isinf(TrialsSubsample)
%                 if sum([nLe,nRe]<=TrialsSubsample)==2
%                     flag_ = 0;
%                 end
%             end
%             if flag_
%                 FR = [LtrialsE;RtrialsE];
%                 evt0 = [ones(nLe,1);2*ones(nRe,1)];
%                 reg = 0.05;
%                 PE=DecodePredErrCrossTemporal_CV_GPU(FR,evt0,reg);
%                 D_Err{iArea}{iDelay}(:,:,iFile) = PE;
%             else
%                 D_Err{iArea}{iDelay}(:,:,iFile) = nan(length_+1,length_+1);
%             end
            
        end
        %% Save results
                if exist('D_')
                    fnOut = sprintf('%sMixedSelectivity%s%s_%s_DelayDecoding_UnitSpan_ErrorBounds.mat',pat,filesep,fname,Areas{iArea});
                    save(fnOut,'D_','D_CIh','D_CIl','tbShort','tbMedium','tbLong','-v7.3');
                    fprintf('Done.\n')
                    clear D_ PE_CIl PE_CIh
                end
        
    end
end

%% plot Correct - mean across sessions
iArea = 2;
figure
% Ticks_ = {[-5 0 4 9],[-5 0 4 8 13],[-5 0 4 8 16 21]}';
Ticks_ = {[0:4:4],[0:4:8],[0:4:16]}';
for iDelay =1:length(Delays_)
    subplot(1,3,iDelay); hold on 
    eval(sprintf('tlims_X = tb%s;',Delays_{iDelay}));
    aboveChance_ = D.AboveChance(:,iDelay);
    decoding_ = D_{iArea}{iDelay}(:,:,aboveChance_);
    for i=1:size(decoding_,3)
        decoding_(:,:,i) = smooth2a(decoding_(:,:,i) ,5,5);
    end
    decoding_ = nanmean(decoding_,3);
%     decoding_ = smooth2a(decoding_ ,5,5);
    imagesc(tlims_X,tlims_X,decoding_);
    plot([0 0],[-5 0],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w','LineWidth',1.5)  
    plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w','LineWidth',1.5);
    plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w') 
%     plot([-5 max(tlims_X)+5],[-5 max(tlims_X)+5],'w','LineWidth',1.5) 

    rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w','LineWidth',1.5)
    rectangle('Position',[-5 -5  max(tlims_X)+5 max(tlims_X)+5],'EdgeColor','w','LineWidth',1.5)
    
%     scatter(0,0,50,'k')
%     scatter(max(tlims_X)-5,max(tlims_X)-5,50,'k')
    
    set(gca,'XTick',Ticks_{iDelay},'YTick',Ticks_{iDelay})
    axis([-5 21 -5 21]); 
    set(gca,'Colormap',jet); caxis([0.4 0.6])
    axis square
    ylabel('Training time (s)')
    xlabel('Testing time (s)')
end
%% difference between animals above and below chance at 16s
iDelay = 3;
% Ticks_ = {[-5 0 4 9],[-5 0 4 8 13],[-5 0 4 8 16 21]}';
Ticks_ = {[0:4:4],[0:4:8],[0:4:16]}';
figure
eval(sprintf('tlims_X = tb%s;',Delays_{iDelay}));
aboveChance_ = D.AboveChance(:,iDelay);

subplot(1,3,1); hold on
    decoding_ = D_{iArea}{iDelay}(:,:,aboveChance_);
    for i=1:size(decoding_,3)
        decoding_(:,:,i) = smooth2a(decoding_(:,:,i) ,5,5);
    end
    decoding_ = nanmean(decoding_,3);
%     decoding_ = smooth2a(decoding_ ,5,5);
    imagesc(tlims_X,tlims_X,decoding_);
    plot([0 0],[-5 0],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w','LineWidth',1.5)  
    plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w','LineWidth',1.5);
    plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w') 
%     plot([-5 max(tlims_X)+5],[-5 max(tlims_X)+5],':k','LineWidth',1.5) 

    rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w','LineWidth',1.5)
    rectangle('Position',[-5 -5  max(tlims_X)+5 max(tlims_X)+5],'EdgeColor','w','LineWidth',1.5)
    axis([-5 21 -5 21]);
    set(gca,'XTick',Ticks_{iDelay},'YTick',Ticks_{iDelay})
    set(gca,'Colormap',jet); caxis([0.4 0.6])
    axis square
    ylabel('Training time (s)')
    xlabel('Testing time (s)')

subplot(1,3,2); hold on
    decoding_ = D_{iArea}{iDelay}(:,:,~aboveChance_);
    for i=1:size(decoding_,3)
        decoding_(:,:,i) = smooth2a(decoding_(:,:,i) ,5,5);
    end
    decoding_ = nanmean(decoding_,3);
%     decoding_ = smooth2a(decoding_ ,5,5);
    imagesc(tlims_X,tlims_X,decoding_);
    plot([0 0],[-5 0],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w','LineWidth',1.5)  
    plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w','LineWidth',1.5);
    plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w') 
%     plot([-5 max(tlims_X)+5],[-5 max(tlims_X)+5],':k','LineWidth',1.5) 

    rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w','LineWidth',1.5)
    rectangle('Position',[-5 -5  max(tlims_X)+5 max(tlims_X)+5],'EdgeColor','w','LineWidth',1.5)
    axis([-5 21 -5 21])
    set(gca,'Colormap',jet); caxis([0.4 0.6])
  set(gca,'XTick',Ticks_{iDelay},'YTick',Ticks_{iDelay})
    axis square
    ylabel('Training time (s)')
    xlabel('Testing time (s)')
subplot(1,3,3); hold on
    corr_ = D_{iArea}{iDelay}(:,:,aboveChance_);
    for i=1:size(corr_,3)
        corr_(:,:,i) = smooth2a(corr_(:,:,i) ,5,5);
    end
    corr_ = nanmean(corr_,3);
%     corr_ = smooth2a(corr_ ,5,5);

    err_ = D_{iArea}{iDelay}(:,:,~aboveChance_);
    for i=1:size(err_,3)
        err_(:,:,i) = smooth2a(err_(:,:,i) ,5,5);
    end
    err_ = nanmean(err_,3);
%     err_ = smooth2a(decoding_ ,5,5);

    decoding_ = (corr_ - err_);%./(corr_ + err_);
%     decoding_ = smooth2a(decoding_ ,20,20);
    imagesc(tlims_X,tlims_X,decoding_);
    plot([0 0],[-5 0],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w','LineWidth',1.5)  
    plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w','LineWidth',1.5);
    plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w') 
%     plot([-5 max(tlims_X)+5],[-5 max(tlims_X)+5],':k','LineWidth',1.5) 

    rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w','LineWidth',1.5)
    rectangle('Position',[-5 -5  max(tlims_X)+5 max(tlims_X)+5],'EdgeColor','w','LineWidth',1.5)
    axis([-5 21 -5 21]); caxis([-0.2 0.2])
    cmap = load('blue_white_red.mat');cmap = flipud(cmap.cmap);
    set(gca,'Colormap',cmap)
  set(gca,'XTick',Ticks_{iDelay},'YTick',Ticks_{iDelay})
    axis square
    ylabel('Training time (s)')
    xlabel('Testing time (s)')
    
%% Plot diagonal
Delays___ = [4 8 16];
figure
for iArea =1:2
    for iDelay =1:length(Delays_)
        subplot(3,1,iDelay);hold on
        plot([-5 22],[0.5 0.5],'--k','linewidth',1.5)
        eval(sprintf('tlims_X = tb%s;',Delays_{iDelay}));
%         idx = true(length(fileList),1);
        idx = D.AboveChance(:,iDelay);
        if sum(~idx)>0
            m = smooth2a(nanmean(D_{iArea}{iDelay}(:,:,~idx),3),2,2);
            e = smooth2a(nansem( D_{iArea}{iDelay}(:,:,~idx),3),2,2);
            m = m(find(eye(length(m))));
            e = e(find(eye(length(m))));
            ciplot(m+e,m-e,tlims_X,color_{iArea},0.3);
        end
        m = smooth2a(nanmean(D_{iArea}{iDelay}(:,:,idx),3),2,2);
        e = smooth2a(nansem( D_{iArea}{iDelay}(:,:,idx),3),2,2);
        m = m(find(eye(length(m))));
        e = e(find(eye(length(m))));
        ciplot(m+e,m-e,tlims_X,color_{iArea},0.6);
        
        plot([Delays___(iDelay),Delays___(iDelay)],[0.2 0.8],':k','LineWidth',1.5);
        
        axis([-5 22 0 1])
    end
end
%% plot Correct - each session

iArea = 1;
for iFile =1:length(fileList)
    figure
for iDelay =1:length(Delays_)
    subplot(1,3,iDelay); hold on 
      eval(sprintf('tlims_X = tb%s;',Delays_{iDelay}));
    imagesc(tlims_X,tlims_X,D_{iArea}{iDelay}(:,:,iFile));
    plot([0 0],[-5 0],'w');
    plot([-5 0],[0 0],'w')  
    plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w');
    plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w');
    plot([-5 0],[0 0],'w') 
    rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w')
    
    axis([-5 21 -5 21]); caxis([0.5 0.6])
    axis square
    
end
end
%% plot Error
iArea = 1;
figure
for iDelay =1:length(Delays_)
    subplot(1,3,iDelay); hold on 
      eval(sprintf('tlims_X = tb%s;',Delays_{iDelay}));
    imagesc(tlims_X,tlims_X,nanmean(D_Err{iArea}{iDelay},3));
    plot([0 0],[-5 0],'w');
    plot([-5 0],[0 0],'w')  
    plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w');
    plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w');
    plot([-5 0],[0 0],'w') 
    rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w')
    
    axis([-5 21 -5 21]); %caxis([0 1])
    axis square
    
end
%% plot Error - each session

iArea = 1;
for iFile =1:length(fileList)
    figure
for iDelay =1:length(Delays_)
    subplot(1,3,iDelay); hold on 
      eval(sprintf('tlims_X = tb%s;',Delays_{iDelay}));
    imagesc(tlims_X,tlims_X,D_Err{iArea}{iDelay}(:,:,iFile));
    plot([0 0],[-5 0],'w');
    plot([-5 0],[0 0],'w')  
    plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w');
    plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w');
    plot([-5 0],[0 0],'w') 
    rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w')
    
    axis([-5 21 -5 21]); caxis([0.5 1])
    axis square
    
end
end
%% plot change on Errors - each session

iArea = 1;
for iFile =1:length(fileList)
    figure
for iDelay =1:length(Delays_)
    subplot(1,3,iDelay); hold on 
      eval(sprintf('tlims_X = tb%s;',Delays_{iDelay}));
    imagesc(tlims_X,tlims_X,D_Err{iArea}{iDelay}(:,:,iFile) - D_{iArea}{iDelay}(:,:,iFile));
    plot([0 0],[-5 0],'w');
    plot([-5 0],[0 0],'w')  
    plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w');
    plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w');
    plot([-5 0],[0 0],'w') 
    rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w')
    
    axis([-5 21 -5 21]); caxis([-0.5 0.5])
    axis square
    
end
end
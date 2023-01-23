%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Target = 'LONG';


shift = 0;
minTrialCount = 3;
subsampleTrials = false;
BootstrapCI_TS = true;
nDrawsTrials = 50;
Nbs = 500;
bw=0.05;
tlimsAll = [-4 4];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
    pat2 = 'C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\LONGTaskonly';
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
    pat2 = '/Volumes/HDD2/DNMTP/raw/KDE_binsTaskonly/LONGTaskonly/';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);
% AnimalID  = {[1,1,2,2,3,3,4,4,5,5,6];[1,1,2,2,3,3,4,4,5,5,6,6];[1,1,2,2,3,3,4,4,5,5,6,6]};
% SessionID = {[1,2,1,2,1,2,1,2,1,2,1];[1,2,1,2,1,2,1,2,1,2,1,2];[1,2,1,2,1,2,1,2,1,2,1,2]};

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

Areas = {'PFC','HP','Joint'};
Areas_ = {'mPFC','dCA1','dCA1-mPFC'};
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

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6],};
Outcome = {'Correct','Errors'};
DelayRange = 1:3;

% load('C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\ClusteringAssemsVSUnits\BSProcessing.mat')
%% Batch process units and Assems
for iFile = 1:length(fileList)
    %% Get the single unit files
    fname = strtok(fileList(iFile).name,'_');
    
    for iArea = 1:2%length(Areas)
        fprintf('Getting units: %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
        load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}),'iFR','Tmtx');
        
        iFR_{iArea} = iFR;
        Tmtx_Units = Tmtx;
    end
    
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname),'t');
    load(sprintf('%s%s.mat',pat,fname),'PFCcells','HPcells');
    noPFC = length(PFCcells);
    clear Tmtx
    %% Get the assembly membership info
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
    Ass.usel_out{3} = [Ass.usel_out{1},Ass.usel_out{2}+noPFC];
    Ass.units       = units;
    Ass.nu = nu; Ass.nu(3) = sum(Ass.nu(1:2));
    % Check that inter-area assems actually span the two areas
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
    for iArea = 1:3
        Ass.LocalMembers{iArea}    = unique(cell2mat(Ass.units{iArea}));
    end
    
    
    Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3}<=Ass.nu(1)));
    Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3}>Ass.nu(1))-Ass.nu(1));
    for iArea = 1:2
        Ass.JointMembers{iArea} = setdiff(unique(cell2mat(Ass.units{iArea})),Ass.LocalMembers{iArea});
        Ass.NonMembers{iArea}   = setdiff(1:Ass.nu(iArea),[Ass.LocalMembers{iArea},Ass.JointMembers{iArea}]);
    end
    for iArea = 1:2
        Ass.NonMembers{iArea}   = Ass.usel_out{iArea}(Ass.NonMembers{iArea});
        Ass.LocalMembers{iArea} = Ass.usel_out{iArea}(Ass.LocalMembers{iArea});
        Ass.JointMembers{iArea} = Ass.usel_out{iArea}(Ass.JointMembers{iArea});
    end
    
    for iArea = 1:2
        %Membership_{iArea} = -ones( Ass.nu(iArea),1);
        Membership_{iArea} = -ones(size(iFR_{iArea},2),1);
        Membership_{iArea}(Ass.NonMembers{iArea})=0;
        Membership_{iArea}(Ass.LocalMembers{iArea})=1;
        Membership_{iArea}(Ass.JointMembers{iArea})=2;
        Membership_{iArea}(setdiff(1:size(iFR_{iArea},2),Ass.usel_out{iArea}))=[];
    end
    for iArea = 1:2
        iFR_{iArea} = iFR_{iArea}(:,Ass.usel_out{iArea});
    end
    
    
    for iArea = 1:2
        % Pad this out to ensure that there's an entry even if there's no detected assembly
        LocalMembersMatrix_{iArea} = nan(Ass.nu(iArea),length(Ass.units{iArea})+1);
        if ~isempty(Ass.units{iArea})
            for iAss=1:length(Ass.units{iArea})
                LocalMembersMatrix_{iArea}(:,iAss+1)=ismember(1:Ass.nu(iArea),Ass.units{iArea}{iAss});
            end
        end
        
        JointMembersMatrix_{iArea} = nan(Ass.nu(iArea),length(Ass.units{3})+1);
        for iAss=1:length(Ass.units{3})
            if ~isempty(Ass.units{3})
                units_ = Ass.units{3}{iAss};
                if iArea==1
                    JointMembersMatrix_{iArea}(:,iAss+1)=ismember(1:Ass.nu(iArea), units_(units_<=Ass.nu(1)));
                else
                    JointMembersMatrix_{iArea}(:,iAss+1)=ismember(1:Ass.nu(iArea), units_(units_>Ass.nu(1))-Ass.nu(1));
                end
            end
            
        end
    end
    iArea = 3;
    JointMembersMatrix_{iArea} = nan(Ass.nu(iArea),length(Ass.units{3})+1);
    for iAss=1:length(Ass.units{3})
        if ~isempty(Ass.units{3})
            units_ = Ass.units{3}{iAss};
            JointMembersMatrix_{iArea}(:,iAss+1)=ismember(1:Ass.nu(iArea), units_);
        end
    end   
    clear units usel_out nu FSCsel Tmtx  iFR noPFC idx nAss U_
                nu = cellfun(@length,Membership_);nu(3) = sum(nu);
%%%%%%%
 nonMemberUnits = cell(3,1);
    for iArea =1:2
        nonMemberUnits{iArea} = find(Membership_{iArea}==0);
    end
    nonMemberUnits{2} = nonMemberUnits{2}+nu(1);
    nonMemberUnits{3} = find([Membership_{1};Membership_{2}]==0);
     
    % Remove local assemblies that are subsets of joint assemblies
    for iArea=1:2
        toremove = zeros(1,size(LocalMembersMatrix_{iArea},2));
        for iAss = 2:size(LocalMembersMatrix_{iArea},2)
            members_=find(LocalMembersMatrix_{iArea}(:,iAss)>0);
            toremove_= zeros(1,size(JointMembersMatrix_{iArea},2));
            for jAss=2:size(JointMembersMatrix_{iArea},2)
                toremove_(jAss) = sum( sum(ismember(members_,find(JointMembersMatrix_{iArea}(:,jAss)))) == length(members_));
            end
            toremove(iAss) = sum(toremove_)>0;
        end
       LocalMembersMatrix_{iArea}(:,toremove==1) =[];
       
       Membership_{iArea} = zeros(size(LocalMembersMatrix_{iArea},1),1);
       Membership2_{iArea} = zeros(size(LocalMembersMatrix_{iArea},1),1);
       
       % 0=nonmember, 1=localmember, 2 = jointmember | local&jointmember
       Membership_{iArea} = max([ nansum(LocalMembersMatrix_{iArea},2),...
                                  2*(nansum(JointMembersMatrix_{iArea},2)>0)...
                                   ],[],2);
       
       % 0=nonmember, 1=localmember, 2 = jointmember, 3=local&jointmember
       Membership2_{iArea} = max([1*(nansum(LocalMembersMatrix_{iArea},2)==1 & nansum(JointMembersMatrix_{iArea},2)==0),... 
                             2*(nansum(LocalMembersMatrix_{iArea},2)<1 & nansum(JointMembersMatrix_{iArea},2)>0),... 
                             3*(nansum(LocalMembersMatrix_{iArea},2)>0 & nansum(JointMembersMatrix_{iArea},2)>0)],[],2); 
    end
    
    clear iAss jAss toremove toremove_ 
    %%%%%%%
    %% Decode Left/Right
    for iOutcome = 1:2
        % get trials
        SampleTrials=[];ChoiceTrials=[];
        for iDelay = DelayRange
            switch iOutcome
                case 1
                    eval(sprintf('SampleTrials = [SampleTrials;t.%s.SamplePress_LeftCorrect;t.%s.SamplePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('ChoiceTrials = [ChoiceTrials;t.%s.ChoicePress_LeftCorrect;t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}));
                case 2
                    eval(sprintf('SampleTrials = [SampleTrials;t.%s.SamplePress_LeftError'';t.%s.SamplePress_LeftError''];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('ChoiceTrials = [ChoiceTrials;t.%s.ChoicePress_RightError'';t.%s.ChoicePress_RightError''];',Delays_{iDelay},Delays_{iDelay}));
            end
        end
        SampleTrials(sum(isnan(SampleTrials),2)>0,:)=[];
        ChoiceTrials(sum(isnan(ChoiceTrials),2)>0,:)=[];
        %% Decode units,  aggregate
        for iArea = 1:2
            clear D_
            if subsampleTrials
                % Draw counterbalanced trial numbers to standardise degrees of freedom
               
                for iDraw = 1:nDrawsTrials
                    
                    Strials = [];nS = 0;Strials_ = {};
                    trials_ = randsample(size(SampleTrials,1),minTrialCount);
                    for iTrial = 1:minTrialCount
                        try
                            tlims_  = SampleTrials(trials_(iTrial),1)/1e6+tlimsAll(1) + shift;
                            tlims_  = closest(Tmtx_Units,tlims_);
                            tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];
                            Strials = [Strials;iFR_{iArea}(tlims_,:)];
                            Strials_{1,nS+1} = iFR_{iArea}(tlims_,:);
                            nS = nS+1;
                        end
                    end
                    
                    Ctrials = [];nC = 0; Ctrials_ = {};
                    trials_ = randsample(size(ChoiceTrials,1),minTrialCount);
                    for iTrial = 1:minTrialCount
                        try
                            tlims_  = ChoiceTrials(trials_(iTrial),1)/1e6+tlimsAll(1) + shift;
                            tlims_  = closest(Tmtx_Units,tlims_);
                            tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1)];
                            Ctrials = [Ctrials;iFR_{iArea}(tlims_,:)];
                            Ctrials_{1,nC+1} = iFR_{iArea}(tlims_,:);
                            nC = nC+1;
                        end
                    end
                    
                    FR   = [Strials;Ctrials];
                    evt0 = [ones(nS,1);2*ones(nC,1)];
                    Ltr  = length(tbAll);
                    %             nu = [length(PFCcells),length(HPcells)];
                    if min([nS,nC])==minTrialCount
                        % (1) Multivariate F-score decoder
                        fprintf('Decoding file %d/%d (%s) %s units, %s %s trials, F-decoder: draw #%d/%d  \n',iFile,length(fileList),fname,Areas{iArea},Outcome{iOutcome}, Delays_{iDelay},iDraw,nDrawsTrials)
                        
                        [D_.avgFR{iDraw},D_.seFR{iDraw},...
                            D_.Ft2(iDraw,:),D_.Rt2(iDraw,:),...
                            D_.Ft2ciL(iDraw,:),D_.Ft2ciH(iDraw,:),...
                            D_.TS{iDraw},...
                            D_.dfnum(iDraw),D_.dfden(iDraw,:)] = DecodeStats(FR,evt0,0.05);
                        
                        
                        
                    else
                        D_.avgFR{iDraw} = nan(Ltr,size(iFR_{iArea},2));
                        D_.seFR{iDraw}  = nan(Ltr,size(iFR_{iArea},2));
                        D_.Ft2(iDraw,:) = nan(1,Ltr);
                        D_.Rt2(iDraw,:) = nan(1,Ltr);
                        D_.Ft2ciL(iDraw,:) = nan(1,Ltr);
                        D_.Ft2ciH(iDraw,:) = nan(1,Ltr);
                        D_.TS{iDraw} = nan(Ltr,size(iFR_{iArea},2));
                        D_.dfnum(iDraw) = NaN;
                        D_.dfden(iDraw,:) = nan(1,Ltr);
                        
                    end
                end
                for Dir = 1:2
                    avgFR_{Dir} = zeros(Ltr,size(iFR_{iArea},2));
                    seFR_{Dir}  = zeros(Ltr,size(iFR_{iArea},2));
                    for iDraw = 1:nDrawsTrials
                        avgFR_{Dir} =  avgFR_{Dir} + D_.avgFR{iDraw}{Dir};
                        seFR_{Dir}  =  seFR_{Dir}  + D_.seFR{iDraw}{Dir};
                    end
                     avgFR_{Dir} =  avgFR_{Dir}./nDrawsTrials;
                     seFR_{Dir}  =  seFR_{Dir}./nDrawsTrials;
                end
                D_.avgFR = avgFR_;
                D_.seFR  = seFR_;
                D_.Ft2(iDraw,:) = nanmean(D_.Ft2,1);
                D_.Rt2(iDraw,:) = nanmean(D_.Rt2,1);
                D_.Ft2ciL(iDraw,:) = nanmean(D_.Ft2ciL,1);
                D_.Ft2ciH(iDraw,:) = nanmean(D_.Ft2ciH,1);
                TS_    = zeros(Ltr,size(iFR_{iArea},2));
                TSsig_ = zeros(Ltr,size(iFR_{iArea},2));
                for iDraw = 1:nDrawsTrials
                    TS_ = TS_ + D_.TS{iDraw};
                    TSsig_  = TSsig_ + tpdf(D_.TS{iDraw},2*minTrialCount-2);
                end
                D_.TS = TS_./nDrawsTrials;
                D_.TSsig = (TSsig_./nDrawsTrials)<0.05;
%                 D_.TSsig = tpdf(D_.TS,2*minTrialCount-2)<0.05;
                
            else % Use all available trials
                    
                Strials = [];nS = 0;Strials_ = {};
                for iTrial =1:size(SampleTrials,1)
                    try
                        tlims_  = SampleTrials(iTrial,1)/1e6 + tlimsAll(1) + shift;
                        tlims_  = closest(Tmtx_Units,tlims_);
                        tlims_  = tlims_(1):(tlims_(1)+length(tbAll)-1);
                        Strials = [Strials;iFR_{iArea}(tlims_,:)];
                        Strials_{1,nS+1} = iFR_{iArea}(tlims_,:);
                        nS = nS+1;
                    end
                end
                
                Ctrials = [];nC = 0;Ctrials_ = {};
                for iTrial =1:size(ChoiceTrials,1)
                    try
                        tlims_  = ChoiceTrials(iTrial,1)/1e6 + tlimsAll(1) + shift;
                        tlims_  = closest(Tmtx_Units,tlims_);
                        tlims_  = tlims_(1):(tlims_(1)+length(tbAll)-1);
                        Ctrials = [Ctrials;iFR_{iArea}(tlims_,:)];
                        Ctrials_{1,nC+1} = iFR_{iArea}(tlims_,:);
                        nC = nC+1;
                    end
                end
                
                FR = [Strials;Ctrials];
                evt0 = [ones(nS,1);2*ones(nC,1)];
                Ltr = length(tbAll);
                if min([nS,nC])>=minTrialCount
                    % (1) Multivariate F-score decoder
                    fprintf('Decoding file %d/%d (%s) %s units, %s %s trials, F-decoder...\n',iFile,length(fileList),fname,Areas{iArea},Outcome{iOutcome}, Delays_{iDelay})
                    
                    [D_.avgFR,D_.seFR,...
                        D_.Ft2,D_.Rt2,...
                        D_.Ft2ciL,D_.Ft2ciH,...
                        D_.TS,...
                        D_.dfnum,D_.dfden] = DecodeStats(FR,evt0,0.05);
                    D_.TSsig = tpdf(D_.TS,nS+nC-2)<0.05;
                    
                    
                    %%%%%%%%%%%%%%%
                    %%% (2) bootstraps (scramble assignments of trajectories to conditions)
                    % shuffle trial outcomes and run DecodeStats for each draw 
                    if BootstrapCI_TS
                        TSbs=zeros(Ltr,size(FR,2),Nbs);
                        parfor b=1:Nbs
%                             disp(b);
                            evt1=evt0;
                            while sum(evt1==evt0)==length(evt0)
                                k=randperm(length(evt0)); evt1=evt0(k);
                            end
                            [~,~,~,~,~,~,TSbs(:,:,b)]=DecodeStats(FR,evt1,0.05);
                        end
                        D_.TSBSciL=zeros(Ltr,size(FR,2));D_.TSBSciH=zeros(Ltr,size(FR,2));
                        for iUnit = 1:size(FR,2)
                            for t_=1:Ltr
                                TSs=sort(squeeze(TSbs(t_,iUnit,:)),'ascend');
                                D_.TSBSciH(t_,iUnit)=TSs(round(0.95*Nbs));
                                D_.TSBSciL(t_,iUnit)=TSs(round(0.05*Nbs));
                            end
                        end
                        D_.TSsigBS = D_.TS>D_.TSBSciH;
                    else
                        D_.TSsigBS = nan(Ltr,size(iFR_{iArea},2));
                    end
                    %%%%%%%%%%%%%%%%%%%
                        
                        
                else
                    D_.TS = nan(Ltr,size(iFR_{iArea},2));
                    D_.TSsig = false(size(D_.TS));
                    D_.TSBSciH = nan(Ltr,size(iFR_{iArea},2));
                    D_.TSBSciL = nan(Ltr,size(iFR_{iArea},2));
                    D_.TSsigBS = nan(Ltr,size(iFR_{iArea},2));
                    D_.Ft2ciL = nan(1,Ltr);
                    D_.Ft2ciH = nan(1,Ltr);
                    D_.Ft2 = nan(1,Ltr);
                    D_.Rt2 = nan(1,Ltr);
                    
                end
            end
            switch iOutcome
                case 1
                    D_SC.TS{iArea}{iFile}          = D_.TS;
                    D_SC.TSsig{iArea}{iFile}       = D_.TSsig;
                    D_SC.TSsigBS{iArea}{iFile}     = D_.TSsigBS;
                    D_SC.avgFR{iArea}{iFile}       = D_.avgFR;
                    D_SC.seFR{iArea}{iFile}        = D_.seFR;
                    
                    D_SC.Membership{iArea}{iFile}  = Membership_{iArea};
                    D_SC.LocalMembersMatrix{iArea}{iFile} = LocalMembersMatrix_{iArea};
                    D_SC.JointMembersMatrix{iArea}{iFile} = JointMembersMatrix_{iArea};
                    D_SC.JointMembersMatrix{3}{iFile} = JointMembersMatrix_{3};
                    
                case 2
                    D_SC.TS_err{iArea}{iFile}      = D_.TS;
                    D_SC.TSsig_err{iArea}{iFile}   = D_.TSsig;
                    D_SC.TSsigBS_err{iArea}{iFile} = D_.TSsigBS;
            end
            
        end
        %% Decode assems, aggregate

        for iArea = 1:3
            FSC = (Ass.FSC{iArea});
            if ~isempty(FSC)
                
                Strials = [];nS = 0;Strials_ = {};
                for iTrial =1:size(SampleTrials,1)
                    try
                        tlims_  = SampleTrials(iTrial,1)/1e6 + tlimsAll(1) + shift;
                        tlims_  = closest(Tmtx_Units,tlims_);
                        tlims_  = tlims_(1):(tlims_(1)+length(tbAll)-1);
                        Strials = [Strials;FSC(tlims_,:)];
                        Strials_{1,nS+1} = FSC(tlims_,:);
                        nS = nS+1;
                    end
                end
                
                Ctrials = [];nC = 0;Ctrials_ = {};
                for iTrial =1:size(ChoiceTrials,1)
                    try
                        tlims_  = ChoiceTrials(iTrial,1)/1e6 + tlimsAll(1) + shift;
                        tlims_  = closest(Tmtx_Units,tlims_);
                        tlims_  = tlims_(1):(tlims_(1)+length(tbAll)-1);
                        Ctrials = [Ctrials;FSC(tlims_,:)];
                        Ctrials_{1,nC+1} = FSC(tlims_,:);
                        nC = nC+1;
                    end
                end
               
                
                
                if min([nS,nC])>=minTrialCount
                    % (1) Multivariate F-score decoder
                    FSC_ = [Strials;Ctrials];
                    evt0 = [ones(nS,1);2*ones(nC,1)];
                    Ltr = length(tbAll);
                    nAss_ = cellfun(@(x)size(x,2),Ass.FSC);
                    % remove zero entries
                    idx = nAss_>0;
                    nAss_(~idx) = min(nAss_(idx));
                    % Multivariate F-score decoder
                    [D_.SC{iArea}.avgFR,...
                        D_.SC{iArea}.seFR,...
                        D_.SC{iArea}.Ft2,...
                        D_.SC{iArea}.Rt2,...
                        D_.SC{iArea}.Ft2ciL,...
                        D_.SC{iArea}.Ft2ciH,...
                        D_.SC{iArea}.TS,...
                        D_.SC{iArea}.dfnum,...
                        D_.SC{iArea}.dfden] = DecodeStats(FSC_,evt0,0.05);
                    D_.SC{iArea}.TSsig = tpdf(D_.SC{iArea}.TS,nS+nC-2)<0.05;
                    
                    %%%%%%%%%%%%%%%
                    %%% (2) bootstraps (scramble assignments of trajectories to conditions)
                    % shuffle trial outcomes and run DecodeStats for each draw 
                    if BootstrapCI_TS
                        TSbs=zeros(Ltr,nAss_(iArea),Nbs);
                        parfor b=1:Nbs
%                             disp(b);
                            evt1=evt0;
                            while sum(evt1==evt0)==length(evt0)
                                k=randperm(length(evt0)); evt1=evt0(k);
                            end
                            [~,~,~,~,~,~,TSbs(:,:,b)]=DecodeStats(FSC_,evt1,0.05);
                        end
                        D_.SC{iArea}.TSBSciL=zeros(Ltr,nAss_(iArea));D_.SC{iArea}.TSBSciH=zeros(Ltr,nAss_(iArea));
                        for iAss = 1:nAss_(iArea)
                            for t_=1:Ltr
                                TSs=sort(squeeze(TSbs(t_,iAss,:)),'ascend');
                                D_.SC{iArea}.TSBSciH(t_,iAss)=TSs(round(0.95*Nbs));
                                D_.SC{iArea}.TSBSciL(t_,iAss)=TSs(round(0.05*Nbs));
                            end
                        end
                        D_.SC{iArea}.TSsigBS = D_.SC{iArea}.TS>D_.SC{iArea}.TSBSciH;
                    else
                        D_.SC{iArea}.TSsigBS = nan(Ltr,nAss_(iArea));
                    end
                    %%%%%%%%%%%%%%%%%%%
                    
                else
                    D_.SC{iArea}.TS = nan(Ltr,size(FSC_,2));
                    D_.SC{iArea}.TSsig = nan(Ltr,size(FSC_,2));
                    D_.SC{iArea}.TSsigBS = nan(Ltr,size(FSC_,2));
                    D_.SC{iArea}.TSBSciH = nan(Ltr,size(FSC_,2));
                    D_.SC{iArea}.TSBSciL = nan(Ltr,size(FSC_,2));
                    D_.SC{iArea}.Ft2ciL = nan(1,Ltr);
                    D_.SC{iArea}.Ft2ciH = nan(1,Ltr);
                    D_.SC{iArea}.Ft2 = nan(1,Ltr);
                    D_.SC{iArea}.Rt2 = nan(1,Ltr);
                    D_.SC{iArea}.avgFR{1} = nan(Ltr,size(FSC_,2));
                    D_.SC{iArea}.avgFR{2} = nan(Ltr,size(FSC_,2));
                end
                
                switch iOutcome
                    case 1
                        D_Ass_SC.TS{iArea}{iFile}          = D_.SC{iArea}.TS;
                        D_Ass_SC.TSsig{iArea}{iFile}       = D_.SC{iArea}.TSsig;
                        D_Ass_SC.TSsigBS{iArea}{iFile}     = D_.SC{iArea}.TSsigBS;
                        D_Ass_SC.ActMean{iArea}{iFile}     = D_.SC{iArea}.avgFR;
                        D_Ass_SC.avgFR{iArea}{iFile}       = D_.SC{iArea}.avgFR;
                        D_Ass_SC.seFR{iArea}{iFile}        = D_.SC{iArea}.seFR;
                        
                    case 2
                        D_Ass_SC.TS_err{iArea}{iFile}      = D_.SC{iArea}.TS;
                        D_Ass_SC.TSsig_err{iArea}{iFile}   = D_.SC{iArea}.TSsig;
                        D_Ass_SC.TSsigBS_err{iArea}{iFile} = D_.SC{iArea}.TSsigBS;
                        D_Ass_SC.ActMean_err{iArea}{iFile} = D_.SC{iArea}.avgFR;
                end
            end
            D_Ass.units{iArea}{iFile} = Ass.units{iArea};
            D_Ass.usel_out{iArea}{iFile} = Ass.usel_out{iArea};
            clear D_
        end      
    end  
end
%% Collapse
% units
for iArea = 1:2
    D_SC.TScollapsed{iArea}             = cell2mat(D_SC.TS{iArea});
    D_SC.TS_errcollapsed{iArea}         = cell2mat(D_SC.TS_err{iArea});
    if BootstrapCI_TS
        D_SC.TSsigcollapsed{iArea}          = cell2mat(D_SC.TSsigBS{iArea});
        D_SC.TS_errsigcollapsed{iArea}      = cell2mat(D_SC.TSsigBS_err{iArea});
    else
        D_SC.TSsigcollapsed{iArea}          = cell2mat(D_SC.TSsig{iArea});
        D_SC.TS_errsigcollapsed{iArea}      = cell2mat(D_SC.TSsig_err{iArea});
    end
    D_SC.membershipsigcollapsed{iArea}  = cell2mat(D_SC.Membership{iArea}');
    
    % Collapse local members
    nAss     = cellfun(@(x) size(x,2),D_SC.LocalMembersMatrix{iArea});
    nAssSum  = cumsum(nAss,2) - nAss;
    nAssSum_ = cumsum(nAss,2,'reverse') - nAss;
    
    nUnits = cellfun(@(x) size(x,1),D_SC.LocalMembersMatrix{iArea});
    D_SC.nUnits{iArea} = nUnits;
    
    for iFile = 1:length(nAss)
        D_SC.LocalMembersMatrixPad{iArea}{iFile} =[nan(nUnits(iFile),nAssSum(iFile)),D_SC.LocalMembersMatrix{iArea}{iFile},nan(nUnits(iFile),nAssSum_(iFile))];
    end
    D_SC.LocalMembersMatrixPadCollapse{iArea} = cell2mat( D_SC.LocalMembersMatrixPad{iArea}');
    D_SC.LocalMembersMatrixPadCollapse{iArea}(:,isnan(nansum(D_SC.LocalMembersMatrixPadCollapse{iArea}))) = [];
    D_SC.LocalMembersMatrixPadCollapse{iArea}(isnan(D_SC.LocalMembersMatrixPadCollapse{iArea})) = 0;
    
    
    % Collapse joint members
    nAss     = cellfun(@(x) size(x,2),D_SC.JointMembersMatrix{iArea});
    nAssSum  = cumsum(nAss,2) - nAss;
    nAssSum_ = cumsum(nAss,2,'reverse') - nAss;
    
    nUnits = cellfun(@(x) size(x,1),D_SC.JointMembersMatrix{iArea});
    
    for iFile = 1:length(nAss)
        D_SC.JointMembersMatrixPad{iArea}{iFile} =[nan(nUnits(iFile),nAssSum(iFile)),D_SC.JointMembersMatrix{iArea}{iFile},nan(nUnits(iFile),nAssSum_(iFile))];
    end
    D_SC.JointMembersMatrixPadCollapse{iArea} = cell2mat( D_SC.JointMembersMatrixPad{iArea}');
    D_SC.JointMembersMatrixPadCollapse{iArea}(:,isnan(nansum(D_SC.JointMembersMatrixPadCollapse{iArea}))) = [];
    D_SC.JointMembersMatrixPadCollapse{iArea}(isnan(D_SC.JointMembersMatrixPadCollapse{iArea})) = 0;
    
    % Finally remove any local assemblies that are also joint
    idxUnits = find(sum(D_SC.LocalMembersMatrixPadCollapse{iArea},2) + sum(D_SC.JointMembersMatrixPadCollapse{iArea},2)==2);
    idxAss = find((sum(D_SC.LocalMembersMatrixPadCollapse{iArea}(idxUnits,:)))>0);
    D_SC.LocalMembersMatrixPadCollapse{iArea}(:,idxAss)=[];
end

iArea = 3;
% Collapse joint members
nAss     = cellfun(@(x) size(x,2),D_SC.JointMembersMatrix{iArea});
nAssSum  = cumsum(nAss,2) - nAss;
nAssSum_ = cumsum(nAss,2,'reverse') - nAss;

nUnits = cellfun(@(x) size(x,1),D_SC.JointMembersMatrix{iArea});

for iFile = 1:length(nAss)
    D_SC.JointMembersMatrixPad{iArea}{iFile} =[nan(nUnits(iFile),nAssSum(iFile)),D_SC.JointMembersMatrix{iArea}{iFile},nan(nUnits(iFile),nAssSum_(iFile))];
end
D_SC.JointMembersMatrixPadCollapse{iArea} = cell2mat( D_SC.JointMembersMatrixPad{iArea}');
D_SC.JointMembersMatrixPadCollapse{iArea}(:,isnan(nansum(D_SC.JointMembersMatrixPadCollapse{iArea}))) = [];
D_SC.JointMembersMatrixPadCollapse{iArea}(isnan(D_SC.JointMembersMatrixPadCollapse{iArea})) = 0;

% assems

for iArea = 1:3
    D_Ass_SC.TScollapsed{iArea}    = cell2mat(D_Ass_SC.TS{iArea}(~isempty_cell(D_Ass_SC.TS{iArea})));
    D_Ass_SC.TS_errcollapsed{iArea}    = cell2mat(D_Ass_SC.TS_err{iArea}(~isempty_cell(D_Ass_SC.TS_err{iArea})));
    if BootstrapCI_TS
        D_Ass_SC.TSsigcollapsed{iArea} = cell2mat(D_Ass_SC.TSsigBS{iArea}(~isempty_cell(D_Ass_SC.TSsigBS{iArea})));
        D_Ass_SC.TS_errsigcollapsed{iArea} = cell2mat(D_Ass_SC.TSsigBS_err{iArea}(~isempty_cell(D_Ass_SC.TSsigBS_err{iArea}) & cellfun(@islogical ,D_Ass_SC.TSsigBS_err{iArea})));
    else
        D_Ass_SC.TSsigcollapsed{iArea} = cell2mat(D_Ass_SC.TSsig{iArea}(~isempty_cell(D_Ass_SC.TSsig{iArea})));
        D_Ass_SC.TS_errsigcollapsed{iArea} = cell2mat(D_Ass_SC.TSsig_err{iArea}(~isempty_cell(D_Ass_SC.TSsig_err{iArea}) & cellfun(@islogical ,D_Ass_SC.TSsig_err{iArea})));
    end
    temp = D_Ass_SC.ActMean{iArea}(~isempty_cell(D_Ass_SC.ActMean{iArea})); temp = temp(:);
    temp_=[];
    for i = 1:length(temp)
        if iscell(temp{i})
            %(1) Mean across all conditions
            temp_ = [temp_; cellfun(@(x) (x{1}+x{2})./2,temp(i,:),'UniformOutput',false)'];
            %(2) Preferred condition only
%             idx = (max(temp{i}{1})<max(temp{i}{2}))+1;
%             temp__ = [];
%             for j = 1:size(temp{i}{1},2)
%                 temp__= [temp__,temp{i}{idx(j)}(:,j)];
%             end
%             temp_ = [temp_;{temp__}];
        else
            temp_ = [temp_; temp{i}];
        end
    end
%     D_Ass_SC.ActMeancollapsed{iArea} = cellfun(@zscore,temp_,'UniformOutput',false);
%     D_Ass_SC.ActMeancollapsed{iArea} = cellfun(@(x) x-nanmean(x(1:10,:)),temp_,'UniformOutput',false);
    D_Ass_SC.ActMeancollapsed{iArea} = temp_;
    
    temp = D_Ass_SC.ActMean_err{iArea}(~isempty_cell(D_Ass_SC.ActMean_err{iArea}))'; temp = temp(:);
    temp_=[];
    for i = 1:length(temp)
        if iscell(temp{i})
            temp_ = [temp_; cellfun(@(x) (x{1}+x{2})./2,temp(i,:),'UniformOutput',false)'];
%             idx = (max(temp{i}{1})<max(temp{i}{2}))+1;
%             temp__ = [];
%             for j = 1:size(temp{i}{1},2)
%                 temp__= [temp__,temp{i}{idx(j)}(:,j)];
%             end
%             temp_ = [temp_;{temp__}];        
        else
            temp_ = [temp_; temp{i}];
        end
    end
%     D_Ass_SC.ActMeanErrcollapsed{iArea} = cellfun(@zscore,temp_,'UniformOutput',false);
%     D_Ass_SC.ActMeanErrcollapsed{iArea} = cellfun(@(x) x-nanmean(x(1:10,:)),temp_,'UniformOutput',false);
    D_Ass_SC.ActMeanErrcollapsed{iArea} = temp_;
clear temp temp_ temp__
end
%% plot some stats
bins = 0:0.5:20;
tlimsAllSC = [-4 4];
tbAllSC = tlimsAllSC(1):bw:tlimsAllSC(2);%0:bw:2*sum(abs(tlimsAll));
Ltr=length(tbAllSC);
t0 = FindClosestIndex(tbAllSC,0); 
% TRange = 1:Ltr; % Whole trial
TRange = 1:round(t0/2); % Pre-press only
% TRange = (round(t0/2)+1):Ltr; % Post-press only
figure
temp_Peak = cell(1,3);
for iArea = 1:2
        subplot(1,2,iArea) ; hold on
        D_   =  cell2mat(D_SC.TS{iArea});
        %     D_   =  double(cell2mat(D_SC.TSsigBS{iArea}));
        %     D_ = D_(range_,:);
        Didx =  cell2mat(D_SC.Membership{iArea}');
        temp_Peak{1} = nanmax(D_(TRange,Didx ==0),1);
        temp_Peak{2} = nanmax(D_(TRange,Didx ==1),1);
        temp_Peak{3} = nanmax(D_(TRange,Didx ==2),1);
        
        temp_y = cellfun(@(x) cumsum(histc(x,bins))./numel(x), temp_Peak, 'UniformOutput', false);
       
        stairs(bins,temp_y{1},'Color',col_{1},'LineWidth',1.5)
        stairs(bins,temp_y{2},'Color',col_{2},'LineWidth',1.5)
        stairs(bins,temp_y{3},'Color',col_{3},'LineWidth',1.5)
        
       

end
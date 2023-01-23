%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Target = 'LONG';


shift = 0;
minTrialCount = 3;
subsampleTrials = false;
nDrawsTrials = 50;

bw=0.05;
tlimsAll = [-5 5];
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
pat2 = 'C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\LONGTaskonly';
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

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6],};
Outcome = {'Correct','Errors'};
DelayRange = 1:3;
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
        LeftTrials=[];RightTrials=[];
        for iDelay = DelayRange
            switch iOutcome
                case 1
                    eval(sprintf('LeftTrials =  [LeftTrials;t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('RightTrials = [RightTrials;t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}));
                case 2
                    eval(sprintf('LeftTrials = [LeftTrials;t.%s.SamplePress_LeftError'',t.%s.ChoicePress_LeftError''];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('RightTrials = [RightTrials;t.%s.SamplePress_RightError'',t.%s.ChoicePress_RightError''];',Delays_{iDelay},Delays_{iDelay}));
            end
            
        end
        LeftTrials(sum(isnan(LeftTrials),2)>0,:)=[];
        RightTrials(sum(isnan(RightTrials),2)>0,:)=[];
        %% Decode units, aggregate
        for iArea = 1:2
            clear D_
            if subsampleTrials
                % Draw counterbalanced trial numbers to standardise degrees of freedom
               
                for iDraw = 1:nDrawsTrials
                    
                    Ltrials = [];nL = 0;Ltrials_ = {};
                    trials_ = randsample(size(LeftTrials,1),minTrialCount);
                    for iTrial = 1:minTrialCount
                        try
                            tlims_  = [LeftTrials(trials_(iTrial),1)/1e6;LeftTrials(trials_(iTrial),2)/1e6]+tlimsAll(1) + shift;
                            tlims_  = closest(Tmtx_Units,tlims_);
                            tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                            Ltrials = [Ltrials;iFR_{iArea}(tlims_,:)];
                            Ltrials_{1,nL+1} = iFR_{iArea}(tlims_,:);
                            nL = nL+1;
                        end
                    end
                    
                    Rtrials = [];nR = 0; Rtrials_ = {};
                    trials_ = randsample(size(RightTrials,1),minTrialCount);
                    for iTrial = 1:minTrialCount
                        try
                            tlims_  = [RightTrials(trials_(iTrial),1)/1e6;RightTrials(trials_(iTrial),2)/1e6]+tlimsAll(1) + shift;
                            tlims_  = closest(Tmtx_Units,tlims_);
                            tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                            Rtrials = [Rtrials;iFR_{iArea}(tlims_,:)];
                            Rtrials_{1,nR+1} = iFR_{iArea}(tlims_,:);
                            nR = nR+1;
                        end
                    end
                    
                    FR   = [Ltrials;Rtrials];
                    evt0 = [ones(nL,1);2*ones(nR,1)];
                    Ltr  = 2*length(tbAll);
                    %             nu = [length(PFCcells),length(HPcells)];
                    if min([nL,nR])==minTrialCount
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
                
                Ltrials = [];nL = 0;Ltrials_ = {};
                for iTrial =1:size(LeftTrials,1)
                    try
                        tlims_  = [LeftTrials(iTrial,1)/1e6;LeftTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                        tlims_  = closest(Tmtx_Units,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                        Ltrials = [Ltrials;iFR_{iArea}(tlims_,:)];
                        Ltrials_{1,nL+1} = iFR_{iArea}(tlims_,:);
                        nL = nL+1;
                    end
                end
                
                Rtrials = [];nR = 0; Rtrials_ = {};
                for iTrial =1:size(RightTrials,1)
                    try
                        tlims_  = [RightTrials(iTrial,1)/1e6;RightTrials(iTrial,2)/1e6]+tlimsAll(1) + shift;
                        tlims_  = closest(Tmtx_Units,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                        Rtrials = [Rtrials;iFR_{iArea}(tlims_,:)];
                        Rtrials_{1,nR+1} = iFR_{iArea}(tlims_,:);
                        nR = nR+1;
                    end
                end
                
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = 2*length(tbAll);
                if min([nL,nR])>=minTrialCount
                    % (1) Multivariate F-score decoder
                    fprintf('Decoding file %d/%d (%s) %s units, %s %s trials, F-decoder...\n',iFile,length(fileList),fname,Areas{iArea},Outcome{iOutcome}, Delays_{iDelay})
                    
                    [D_.avgFR,D_.seFR,...
                        D_.Ft2,D_.Rt2,...
                        D_.Ft2ciL,D_.Ft2ciH,...
                        D_.TS,...
                        D_.dfnum,D_.dfden] = DecodeStats(FR,evt0,0.05);
                    D_.TSsig = tpdf(D_.TS,nL+nR-2)<0.05;
                else
                    D_.TS = nan(Ltr,size(iFR_{iArea},2));
                    D_.TSsig = false(size(D_.TS));
                    D_.Ft2ciL = nan(1,Ltr);
                    D_.Ft2ciH = nan(1,Ltr);
                    D_.Ft2 = nan(1,Ltr);
                    D_.Rt2 = nan(1,Ltr);
                    
                end
            end
            switch iOutcome
                case 1
                    D.TS{iArea}{iFile}          = D_.TS;
                    D.TSsig{iArea}{iFile}       = D_.TSsig;
                    D.Membership{iArea}{iFile}  = Membership_{iArea};
                    D.LocalMembersMatrix{iArea}{iFile} = LocalMembersMatrix_{iArea};
                    D.JointMembersMatrix{iArea}{iFile} = JointMembersMatrix_{iArea};
                    D.JointMembersMatrix{3}{iFile} = JointMembersMatrix_{3};
                    
                case 2
                    D.TS_err{iArea}{iFile}      = D_.TS;
                    D.TSsig_err{iArea}{iFile}   = D_.TSsig;
            end
            
        end
        %% Decode assems, aggregate

        for iArea = 1:3
            FSC = (Ass.FSC{iArea});
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
                if min([nL,nR])>=minTrialCount
                    % (1) Multivariate F-score decoder
                    FSC_ = [Ltrials;Rtrials];
                    evt0 = [ones(nL,1);2*ones(nR,1)];
                    Ltr = 2*length(tbAll);
                    nAss_ = cellfun(@(x)size(x,2),Ass.FSC);
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
                    
                else
                    D_.LR{iArea}.TS = nan(Ltr,size(FSC_,2));
                    D_.LR{iArea}.TSsig = nan(Ltr,size(FSC_,2));
                    D_.LR{iArea}.Ft2ciL = nan(1,Ltr);
                    D_.LR{iArea}.Ft2ciH = nan(1,Ltr);
                    D_.LR{iArea}.Ft2 = nan(1,Ltr);
                    D_.LR{iArea}.Rt2 = nan(1,Ltr);
                    D_.LR{iArea}.avgFR{1} = nan(Ltr,size(FSC_,2));
                    D_.LR{iArea}.avgFR{2} = nan(Ltr,size(FSC_,2));
                end
                
                
                switch iOutcome
                    case 1
                        D_Ass.TS{iArea}{iFile}          = D_.LR{iArea}.TS;
                        D_Ass.TSsig{iArea}{iFile}       = D_.LR{iArea}.TSsig;
                        D_Ass.ActMean{iArea}{iFile}     = D_.LR{iArea}.avgFR;
                        
                    case 2
                        D_Ass.TS_err{iArea}{iFile}      = D_.LR{iArea}.TS;
                        D_Ass.TSsig_err{iArea}{iFile}   = D_.LR{iArea}.TSsig;
                        D_Ass.ActMean_err{iArea}{iFile} =  D_.LR{iArea}.avgFR;
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
    D.TScollapsed{iArea}             = cell2mat(D.TS{iArea});
    D.TSsigcollapsed{iArea}          = cell2mat(D.TSsig{iArea});
    D.TS_errcollapsed{iArea}         = cell2mat(D.TS_err{iArea});
    D.TS_errsigcollapsed{iArea}      = cell2mat(D.TSsig_err{iArea});
    D.membershipsigcollapsed{iArea}  = cell2mat(D.Membership{iArea}');
    
    % Collapse local members
    nAss     = cellfun(@(x) size(x,2),D.LocalMembersMatrix{iArea});
    nAssSum  = cumsum(nAss,2) - nAss;
    nAssSum_ = cumsum(nAss,2,'reverse') - nAss;
    
    nUnits = cellfun(@(x) size(x,1),D.LocalMembersMatrix{iArea});
    D.nUnits{iArea} = nUnits;
    
    for iFile = 1:length(nAss)
        D.LocalMembersMatrixPad{iArea}{iFile} =[nan(nUnits(iFile),nAssSum(iFile)),D.LocalMembersMatrix{iArea}{iFile},nan(nUnits(iFile),nAssSum_(iFile))];
    end
    D.LocalMembersMatrixPadCollapse{iArea} = cell2mat( D.LocalMembersMatrixPad{iArea}');
    D.LocalMembersMatrixPadCollapse{iArea}(:,isnan(nansum(D.LocalMembersMatrixPadCollapse{iArea}))) = [];
    D.LocalMembersMatrixPadCollapse{iArea}(isnan(D.LocalMembersMatrixPadCollapse{iArea})) = 0;
    
    
    % Collapse joint members
    nAss     = cellfun(@(x) size(x,2),D.JointMembersMatrix{iArea});
    nAssSum  = cumsum(nAss,2) - nAss;
    nAssSum_ = cumsum(nAss,2,'reverse') - nAss;
    
    nUnits = cellfun(@(x) size(x,1),D.JointMembersMatrix{iArea});
    
    for iFile = 1:length(nAss)
        D.JointMembersMatrixPad{iArea}{iFile} =[nan(nUnits(iFile),nAssSum(iFile)),D.JointMembersMatrix{iArea}{iFile},nan(nUnits(iFile),nAssSum_(iFile))];
    end
    D.JointMembersMatrixPadCollapse{iArea} = cell2mat( D.JointMembersMatrixPad{iArea}');
    D.JointMembersMatrixPadCollapse{iArea}(:,isnan(nansum(D.JointMembersMatrixPadCollapse{iArea}))) = [];
    D.JointMembersMatrixPadCollapse{iArea}(isnan(D.JointMembersMatrixPadCollapse{iArea})) = 0;
    
    % Finally remove any local assemblies that are also joint
    idxUnits = find(sum(D.LocalMembersMatrixPadCollapse{iArea},2) + sum(D.JointMembersMatrixPadCollapse{iArea},2)==2);
    idxAss = find((sum(D.LocalMembersMatrixPadCollapse{iArea}(idxUnits,:)))>0);
    D.LocalMembersMatrixPadCollapse{iArea}(:,idxAss)=[];
end

iArea = 3;
% Collapse joint members
nAss     = cellfun(@(x) size(x,2),D.JointMembersMatrix{iArea});
nAssSum  = cumsum(nAss,2) - nAss;
nAssSum_ = cumsum(nAss,2,'reverse') - nAss;

nUnits = cellfun(@(x) size(x,1),D.JointMembersMatrix{iArea});

for iFile = 1:length(nAss)
    D.JointMembersMatrixPad{iArea}{iFile} =[nan(nUnits(iFile),nAssSum(iFile)),D.JointMembersMatrix{iArea}{iFile},nan(nUnits(iFile),nAssSum_(iFile))];
end
D.JointMembersMatrixPadCollapse{iArea} = cell2mat( D.JointMembersMatrixPad{iArea}');
D.JointMembersMatrixPadCollapse{iArea}(:,isnan(nansum(D.JointMembersMatrixPadCollapse{iArea}))) = [];
D.JointMembersMatrixPadCollapse{iArea}(isnan(D.JointMembersMatrixPadCollapse{iArea})) = 0;


% assems

for iArea = 1:3
    D_Ass.TScollapsed{iArea}    = cell2mat(D_Ass.TS{iArea}(~isempty_cell(D_Ass.TS{iArea})));
    D_Ass.TSsigcollapsed{iArea} = cell2mat(D_Ass.TSsig{iArea}(~isempty_cell(D_Ass.TSsig{iArea})));
    
    D_Ass.TS_errcollapsed{iArea}    = cell2mat(D_Ass.TS_err{iArea}(~isempty_cell(D_Ass.TS_err{iArea})));
    D_Ass.TS_errsigcollapsed{iArea} = cell2mat(D_Ass.TSsig_err{iArea}(~isempty_cell(D_Ass.TSsig_err{iArea}) & cellfun(@islogical ,D_Ass.TSsig_err{iArea})));
    
    temp = D_Ass.ActMean{iArea}(~isempty_cell(D_Ass.ActMean{iArea})); temp = temp(:);
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
    D_Ass.ActMeancollapsed{iArea} = cellfun(@zscore,temp_,'UniformOutput',false);
%     D_Ass.ActMeancollapsed{iArea} = temp_;
    
    temp = D_Ass.ActMean_err{iArea}(~isempty_cell(D_Ass.ActMean_err{iArea}))'; temp = temp(:);
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
    D_Ass.ActMeanErrcollapsed{iArea} = cellfun(@zscore,temp_,'UniformOutput',false);
%     D_Ass.ActMeanErrcollapsed{iArea} = temp_;
clear temp temp_ temp__
end
%% Plot average assembly activation
figure;
for iArea = 1:3
    subplot(1,3,iArea); hold on
    plot([5 5],[-1.5 1],'color',[0 1 0 0.5],'LineWidth',1.5,'HandleVisibility','off')
    plot([15 15],[-1.5 1],'color',[1 0 0 0.5],'LineWidth',1.5,'HandleVisibility','off')
    plot([10 10],[-1.5 1],'color',[0 0 0 0.5],'LineStyle',':','LineWidth',1.5,'HandleVisibility','off')
    
    y = cell2mat(cellfun(@(x) mean(x,2),D_Ass.ActMeancollapsed{iArea},'UniformOutput',false)')';
    yE = cell2mat(cellfun(@(x) mean(x,2),D_Ass.ActMeanErrcollapsed{iArea},'UniformOutput',false)')';
    x = (1:size(y,2)).*bw;

    ciplot(nanmean(yE)+nansem(yE),nanmean(yE)-nansem(yE),x,[1 0 0],0.6) 
    ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),x,0*[1 1 1],1)        

    title([Areas{iArea}, ' assemblies'])
    if iArea==1 
       ylabel({'Assembly activity';'(Mean Factor score)'}) 
    elseif iArea==2
%        xlabel('Time (s)') 

    end
    axis([0 20 -2 2])
    
end
legend('Error','Correct','Location','NorthEast'); legend boxoff

figure; hold on
plot([5 5],[-1 1],'color',[0 1 0 0.5],'LineWidth',1.5,'HandleVisibility','off')
plot([15 15],[-1 1],'color',[1 0 0 0.5],'LineWidth',1.5,'HandleVisibility','off')
plot([10 10],[-1 1],'color',[0 0 0 0.5],'LineStyle',':','LineWidth',1.5,'HandleVisibility','off')
for iArea = 1:3
    
    y = cell2mat(cellfun(@(x) mean(x,2),D_Ass.ActMeancollapsed{iArea},'UniformOutput',false)')';
    x = (1:size(y,2)).*bw;

    ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),x,color_{iArea},0.8)        
    
end
title('Average assembly activation')
    ylabel({'Assembly activity';'(Mean Factor score)'}) 
    xlabel('Time (s)') 
    axis([0 20 -2 2])
    legend(Areas,'Location','SouthWest'); legend boxoff
%%
%% Peak decoding - units
x = 0:0.1:20;
% x2 = -20:1:20;
x2 = -1:0.2:1;

compfunction = @(a,b) (b-a)./(a+b);
clear histPeak histPeak_Err histPeak_Delta shiftPeak PeakCollapse PeakCollapseErr PeakCollapseDelta
PeakCollapse = cell(2,3);
PeakCollapseErr = cell(2,3);
PeakCollapseDelta = cell(2,3);
for iArea = 1:2
    for iFile = 1:length(fileList)
        TS_ = D.TS{iArea}{iFile};
        TS_Err = D.TS_err{iArea}{iFile};
        for iClass = 1:3
            idx    =  D.Membership{iArea}{iFile} == (iClass-1);
            temp   = TS_(:,idx);
            a = nanmax(temp);
            y = cumsum(histc(a,x));
            y = y./nanmax(y);
            histPeak{iArea,iClass}(iFile,:) = y;
            PeakCollapse{iArea,iClass} = [PeakCollapse{iArea,iClass}; a'];
            
            tempErr = TS_Err(:,idx);
            b = nanmax(tempErr);
            y = cumsum(histc(b,x));
            y = y./nanmax(y);
          
            histPeak_Err{iArea,iClass}(iFile,:) = y;
            PeakCollapseErr{iArea,iClass} = [PeakCollapseErr{iArea,iClass}; b'];

            
            if ~isempty(a)
                y = compfunction(a,b);
                PeakCollapseDelta{iArea,iClass}=[PeakCollapseDelta{iArea,iClass}; y'];
                shiftPeak{iArea,iClass}{iFile} = y';
                y = histc(y,x2);
                y = y./nansum(y);
                histPeak_Delta{iArea,iClass}(iFile,:) = y;
            else
                histPeak_Delta{iArea,iClass}(iFile,:) = nan(1,length(x2));
            end
        end
        
        
    end
end
histPeakMean = cellfun(@nanmean,histPeak,'UniformOutput',false);
histPeakSEM  = cellfun(@nansem,histPeak,'UniformOutput',false);
histPeakMean_Err = cellfun(@nanmean,histPeak_Err,'UniformOutput',false);
histPeakSEM_Err  = cellfun(@nansem,histPeak_Err,'UniformOutput',false);

histPeakMean_Delta = cellfun(@nanmean,histPeak_Delta,'UniformOutput',false);
histPeakSEM_Delta  = cellfun(@nansem,histPeak_Delta,'UniformOutput',false);
%% Peak decoding - correct hists
figure('name','Peak decoding')
for iArea = 1:2
    subplot(1,2,iArea);hold on
    title([Areas{iArea} ' units'])
    for iClass = 1:3
        ciplot(histPeakMean{iArea,iClass}+histPeakSEM{iArea,iClass},...
            histPeakMean{iArea,iClass}-histPeakSEM{iArea,iClass},...
            x,col_{iClass},0.6)
        
%                 ciplot(histPeakMean_Err{iArea,iClass}+histPeakSEM_Err{iArea,iClass},...
%                     histPeakMean_Err{iArea,iClass}-histPeakSEM_Err{iArea,iClass},...
%                     x,col_{iClass},0.5)
        xlabel('Peak t-score')
        if iArea==1
           ylabel('Fraction of Units') 
        end
    end
    axis tight
end
% legend(MemberClasses_,'Location','SouthEast');legend boxoff
%% Peak decoding - change on errors
clear h p 
rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')
figure;
for iArea = 1:2
    subplot(1,2,iArea);hold on
    title([Areas{iArea} ' units'])
    for iClass = [3 2 1]
        ciplot(histPeakMean_Delta{iArea,iClass}+histPeakSEM_Delta{iArea,iClass},...
            histPeakMean_Delta{iArea,iClass}-histPeakSEM_Delta{iArea,iClass},...
            x2,col_{iClass},0.9)
%         [h(iArea,iClass),p(iArea,iClass)] =ttest(PeakCollapse{iArea,iClass},PeakCollapseErr{iArea,iClass});
%         [h(iArea,iClass),p(iArea,iClass)] =ttest(PeakCollapseDelta{iArea,iClass});
        [h(iArea,iClass),p(iArea,iClass)] =signrank(PeakCollapseDelta{iArea,iClass});
%         [h(iArea,iClass),p(iArea,iClass)] =signrank(PeakCollapse{iArea,iClass},PeakCollapseErr{iArea,iClass});

    end
    plot([0 0],[0 0.7],'k','LineWidth',1.5)

    axis([min(x2) max(x2) 0 1])
    if iArea ==1
       ylabel('Fraction of units')        
    end
    xlabel('Shift on errors')
end
addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')

%% Time sig - units
x = bw:bw:20;
x2 = -20:2:20;
% x2 = -3:0.25:3;
compfunction = @(a,b) (b-a);%./(a+b);
chooseLongest = false; %(if true, pick the longest significant block of decoding, else sum all)
clear histTime histTime_Err histTime_Delta shiftTime TimeCollapse TimeCollapseErr TimeCollapseDelta
TimeCollapse = cell(2,3);
TimeCollapseErr = cell(2,3);
TimeCollapseDelta = cell(2,3);
for iArea = 1:2
    for iFile = 1:length(fileList)
        TSsig_ = D.TSsig{iArea}{iFile};
        TSsig_Err = D.TSsig_err{iArea}{iFile};
        for iClass = 1:3
            idx =  D.Membership{iArea}{iFile} == (iClass-1);
            if sum(idx>0)
                temp = TSsig_(:,idx);
                tempErr = TSsig_Err(:,idx);
                a=[];b=[];
                if chooseLongest
                for i=1:size(temp,2)
                    [~,~,w,~] = findpeaks(single(temp(:,i)));
                    if ~isempty(w)
                        a(i) = max(w)*bw;
                    else
                        a(i) = 0;
                    end
                end
                for i=1:size(tempErr,2)
                    [~,~,w,~] = findpeaks(single(tempErr(:,i)));
                    if ~isempty(w)
                        b(i) = max(w)*bw;
                    else
                        b(i) = 0;
                    end
                end
                else
                    a = nansum(temp,1).*bw;
                    b = nansum(tempErr,1).*bw;
                end
                    y = cumsum(histc(a,x));
                    y = y./nanmax(y);
                    histTime{iArea,iClass}(iFile,:) = y;
                    TimeCollapse{iArea,iClass} = [TimeCollapse{iArea,iClass}; a'];

                    y = cumsum(histc(b,x));
                    y = y./nanmax(y);
                    histTime_Err{iArea,iClass}(iFile,:) = y;
                    TimeCollapseErr{iArea,iClass} = [TimeCollapseErr{iArea,iClass}; b'];

                if ~isempty(a) && ~isempty(b)
                    y = compfunction(a,b);
                    TimeCollapseDelta{iArea,iClass}=[TimeCollapseDelta{iArea,iClass}; y'];

                    shiftTime{iArea,iClass}{iFile} = y';
                    y = histc(y,x2);
                    y = y./nansum(y);
                    histTime_Delta{iArea,iClass}(iFile,:) = y;
                else
                    histTime_Delta{iArea,iClass}(iFile,:) = nan(1,length(x2));
                end
                
            else
                histTime{iArea,iClass}(iFile,:) = nan(size(x));
                histTime_Err{iArea,iClass}(iFile,:) = nan(size(x));
                histTime_Delta{iArea,iClass}(iFile,:) = nan(size(x2));
            end
            %         plot(x,y,'Color',col_{iClass});
        end
    end
end
histTimeMean = cellfun(@nanmean,histTime,'UniformOutput',false);
histTimeSEM = cellfun(@nansem,histTime,'UniformOutput',false);
histTimeMean_Err = cellfun(@nanmean,histTime_Err,'UniformOutput',false);
histTimeSEM_Err = cellfun(@nansem,histTime_Err,'UniformOutput',false);

histTimeMean_Delta = cellfun(@nanmean,histTime_Delta,'UniformOutput',false);
histTimeSEM_Delta = cellfun(@nansem,histTime_Delta,'UniformOutput',false);
%% Time sig - correct hists
figure('name','Time significant');
for iArea = 1:2
    subplot(1,2,iArea);hold on
    title([Areas{iArea} ' units'])
    for iClass = 1:3
        ciplot(histTimeMean{iArea,iClass}+histTimeSEM{iArea,iClass},...
            histTimeMean{iArea,iClass}-histTimeSEM{iArea,iClass},...
            x,col_{iClass},1)
        %         ciplot(histTimeMean_Err{iArea,iClass}+histTimeSEM_Err{iArea,iClass},...
        %             histTimeMean_Err{iArea,iClass}-histTimeSEM_Err{iArea,iClass},...
        %             x,col_{iClass},0.5)
        
         xlabel('Significant coding span (s)')
        if iArea==1
           ylabel('Fraction of Units') 
        end
    end
end
%% Time sig - change on errors
rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')

figure;
for iArea = 1:2
    subplot(1,2,iArea);hold on
    title([Areas{iArea} ' units'])
    plot([0 0],[0 0.7],'k','LineWidth',1.5)
    for iClass = [1 3 2]%1:3
        ciplot(histTimeMean_Delta{iArea,iClass}+histTimeSEM_Delta{iArea,iClass},...
            histTimeMean_Delta{iArea,iClass}-histTimeSEM_Delta{iArea,iClass},...
            x2,col_{iClass},0.9)
        
        %         [h(iArea,iClass),p(iArea,iClass)] =ttest(TimeCollapse{iArea,iClass},TimeCollapseErr{iArea,iClass});
%         [h(iArea,iClass),p(iArea,iClass)] =ttest(TimeCollapseDelta{iArea,iClass});
        [h(iArea,iClass),p(iArea,iClass)] =signrank(TimeCollapseDelta{iArea,iClass});
%         [h(iArea,iClass),p(iArea,iClass)] =signrank(TimeCollapse{iArea,iClass},TimeCollapseErr{iArea,iClass});

if iArea ==1
    ylabel('Fraction of units')
end
xlabel('Shift on errors')

    end
    axis([min(x2) max(x2) 0 1])
end
addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')

%% Units: time vs peak shift
clear Mx My Sx Sy
figure; 
for iArea = 1:2
    subplot(1,2,iArea);hold on
    title([Areas{iArea} ' units'])

    plot([-20 20],[0 0],':k')
    plot([0 0],[-20 20],':k')
    for iClass = [1 3 2]
        y = cell2mat( shiftPeak{iArea,iClass}');
        x = cell2mat( shiftTime{iArea,iClass}');
        Mx{iArea,iClass} = nanmean(x);Sx{iArea,iClass} = nanstd(x);
        My{iArea,iClass} = nanmean(y);Sy{iArea,iClass} = nanstd(y);
        scatter(x,y,40,col_{iClass},'filled','MarkerFaceAlpha',0.5)
        
    end
%     for iClass = [1 3 2]
%         
%         errorbar(Mx{iArea,iClass},My{iArea,iClass},...
%                  Sy{iArea,iClass},Sy{iArea,iClass},...
%                  Sx{iArea,iClass},Sx{iArea,iClass},...
%             '-o','MarkerSize',5,...
%             'MarkerEdgeColor','k','MarkerFaceColor',col_{iClass},'Color','k')
%     end
    
%     axis([-2 2 -2 2])
    axis([-20 20 -2 2])
    if iArea==1
        ylabel({'Change in';'Peak decoding'})
    end
    xlabel({'Change in span of';' significant decoding'})

end

%% Peak decoding - assems
x = 0:0.1:20;
% x2 = -20:1:20;
x2 = -1:0.2:1;
compfunction = @(a,b) (b-a)./(a+b);

clear histPeak_Ass histPeak_AssErr histPeak_AssDelta shiftPeak_Ass PeakCollapseAss PeakCollapseErrAss PeakCollapseDeltaAss
PeakCollapseAss = cell(3,1);
PeakCollapseErrAss = cell(3,1);
PeakCollapseDeltaAss = cell(3,1);
for iArea = 1:3
    temp = D_Ass.TS{iArea}(~isempty_cell(D_Ass.TS{iArea}));
    tempErr = D_Ass.TS_err{iArea}(~isempty_cell(D_Ass.TS_err{iArea}));
    for iFile = 1:length(temp)
        TS_    = temp{iFile};
        TS_Err = tempErr{iFile};
        
        a = nanmax(TS_);
        PeakCollapseAss{iArea} = [PeakCollapseAss{iArea}; a'];

        b = nanmax(TS_Err);
        PeakCollapseErrAss{iArea} = [PeakCollapseErrAss{iArea}; b'];

        y = cumsum(histc(a,x));
        y = y./nanmax(y);
        histPeak_Ass{iArea}(iFile,:) = y;
        
        
        y = cumsum(histc(b,x));
        y = y./nanmax(y);
        histPeak_AssErr{iArea}(iFile,:) = y;
        
        y = compfunction(a,b);
        PeakCollapseDeltaAss{iArea} = [PeakCollapseDeltaAss{iArea}; y'];

        shiftPeak_Ass{iArea}{iFile} = y';
        y = histc(y,x2);
        y = y./nansum(y);
        histPeak_AssDelta{iArea}(iFile,:) = y;
    end
    
end
histPeakMean_Ass    = cellfun(@nanmean,histPeak_Ass,'UniformOutput',false);
histPeakSEM_Ass     = cellfun(@nansem,histPeak_Ass,'UniformOutput',false);
histPeakMean_AssErr = cellfun(@nanmean,histPeak_AssErr,'UniformOutput',false);
histPeakSEM_AssErr  = cellfun(@nansem,histPeak_AssErr,'UniformOutput',false);
histPeakMean_AssDelta = cellfun(@nanmean,histPeak_AssDelta,'UniformOutput',false);
histPeakSEM_AssDelta  = cellfun(@nansem,histPeak_AssDelta,'UniformOutput',false);
%% Peak decoding - correct hists
figure('name','Peak decoding');hold on
for iArea = 1:3
    %     subplot(1,3,iArea);hold on
    %     title(Areas{iArea})
    ciplot(histPeakMean_Ass{iArea}+histPeakSEM_Ass{iArea},...
        histPeakMean_Ass{iArea}-histPeakSEM_Ass{iArea},...
        x,color_{iArea},0.8)
    
    
%              ciplot(histPeakMean_AssErr{iArea}+histPeakSEM_AssErr{iArea},...
%                 histPeakMean_AssErr{iArea}-histPeakSEM_AssErr{iArea},...
%                 x,color_{iArea},0.3)
end
legend(Areas,'Location','SouthEast'); legend boxoff
xlabel('Peak t-score')
ylabel('Fraction of assemblies')
%% Peak decoding - change on errors
clear h p 
rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')
figure;hold on
plot([0 0],[0 0.7],'k','LineWidth',1.5)
for iArea = 1:3
    ciplot(histPeakMean_AssDelta{iArea}+histPeakSEM_AssDelta{iArea},...
        histPeakMean_AssDelta{iArea}-histPeakSEM_AssDelta{iArea},...
        x2,color_{iArea},0.9)
    %         [h(iArea),p(iArea)] =ttest(PeakCollapseAss{iArea},PeakCollapseErrAss{iArea});
    %         [h(iArea),p(iArea)] =ttest(PeakCollapseDeltaAss{iArea});
    [h(iArea),p(iArea)] =signrank(PeakCollapseDeltaAss{iArea});
    %         [h(iArea),p(iArea)] =signrank(PeakCollapseAss{iArea},PeakCollapseErrAss{iArea});

    

end
axis([min(x2) max(x2) 0 1])
ylabel('Fraction of assemblies')
xlabel('Shift on errors')
  
addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')

%% time sig - assems
x = bw:bw:20;
x2 = -20:2:20;
% x2 = -3:0.25:3;
compfunction = @(a,b) (b-a);%./(a+b);
chooseLongest = false; %(if true, pick the longest significant block of decoding, else sum all)

clear histTime_Ass histTime_AssErr histTime_AssDelta shiftTime_Ass TimeCollapseAss TimeCollapseErrAss TimeCollapseDeltaAss
TimeCollapseAss = cell(3,1);
TimeCollapseErrAss = cell(3,1);
TimeCollapseDeltaAss = cell(3,1);
for iArea = 1:3
    temp = D_Ass.TSsig{iArea}(~isempty_cell(D_Ass.TSsig{iArea}));
    tempErr = D_Ass.TSsig_err{iArea}(~isempty_cell(D_Ass.TSsig_err{iArea}));
    for iFile = 1:length(temp)
        TSsig_    = temp{iFile};
        TSsig_Err = tempErr{iFile};
        
        a=[];b=[];
        if chooseLongest
            for i=1:size(TSsig_,2)
                [~,~,w,~] = findpeaks(single(TSsig_(:,i)));
                if ~isempty(w)
                    a(i) = max(w)*bw;
                else
                    a(i) = 0;
                end
            end
            for i=1:size(TSsig_Err,2)
                [~,~,w,~] = findpeaks(single(TSsig_Err(:,i)));
                if ~isempty(w)
                    b(i) = max(w)*bw;
                else
                    b(i) = 0;
                end
            end
        else
            a = nansum(TSsig_,1).*bw;
            b = nansum(TSsig_Err,1).*bw;
        end
        TimeCollapseAss{iArea} = [TimeCollapseAss{iArea}; a'];
        TimeCollapseErrAss{iArea} = [TimeCollapseErrAss{iArea}; b'];

        y = cumsum(histc(a,x));
        y = y./max(y);
        histTime_Ass{iArea}(iFile,:) = y;
        
        y = cumsum(histc(b,x));
        y = y./max(y);
        histTime_AssErr{iArea}(iFile,:) = y;
        
        y = compfunction(a,b);
        TimeCollapseDeltaAss{iArea} = [TimeCollapseDeltaAss{iArea}; y'];

        
        shiftTime_Ass{iArea}{iFile} = y';

        y = histc(y,x2);
        y = y./nansum(y);
        histTime_AssDelta{iArea}(iFile,:) = y;
    end  
end
histTimeMean_Ass = cellfun(@nanmean,histTime_Ass,'UniformOutput',false);
histTimeSEM_Ass = cellfun(@nansem,histTime_Ass,'UniformOutput',false);
histTimeMean_AssErr = cellfun(@nanmean,histTime_AssErr,'UniformOutput',false);
histTimeSEM_AssErr = cellfun(@nansem,histTime_AssErr,'UniformOutput',false);
histTimeMean_AssDelta = cellfun(@nanmean,histTime_AssDelta,'UniformOutput',false);
histTimeSEM_AssDelta = cellfun(@nansem,histTime_AssDelta,'UniformOutput',false);
%% Time sig - correct hists
figure('name','Time significant');hold on
for iArea = 1:3
    ciplot(histTimeMean_Ass{iArea}+histTimeSEM_Ass{iArea},...
        histTimeMean_Ass{iArea}-histTimeSEM_Ass{iArea},...
        x,color_{iArea},1)
%             ciplot(histTimeMean_AssErr{iArea}+histTimeSEM_AssErr{iArea},...
%                 histTimeMean_AssErr{iArea}-histTimeSEM_AssErr{iArea},...
%                 x,color_{iArea},0.3)
%     plot(x,histTime_Ass{iArea},'color',color_{iArea})
end
axis tight
legend(Areas,'Location','SouthEast'); legend boxoff
xlabel('Significant coding span (s)')
ylabel('Fraction of assemblies')
%% Time sig - change on errors
clear h p 
rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')
figure;hold on
plot([0 0],[0 0.7],'k','LineWidth',1.5)
for iArea = 1:3
    
    
    ciplot(histTimeMean_AssDelta{iArea}+histTimeSEM_AssDelta{iArea},...
        histTimeMean_AssDelta{iArea}-histTimeSEM_AssDelta{iArea},...
        x2,color_{iArea},0.9)
     %         [h(iArea),p(iArea)] =ttest(TimeCollapseAss{iArea},TimeCollapseErrAss{iArea});
    %         [h(iArea),p(iArea)] =ttest(TimeCollapseDeltaAss{iArea});
    [h(iArea),p(iArea)] =signrank(TimeCollapseDeltaAss{iArea});
    %         [h(iArea),p(iArea)] =signrank(TimeCollapseAss{iArea},TimeCollapseErrAss{iArea});
end
axis([min(x2) max(x2) 0 1])
ylabel('Fraction of assemblies')
xlabel('Shift on errors')
addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')

%% Assems: time vs peak shift
clear Mx My Sx Sy
figure; hold on
plot([-20 20],[0 0],':k')
plot([0 0],[-20 20],':k')
for iArea = 1:3
    y = cell2mat( shiftPeak_Ass{iArea}');
    x = cell2mat( shiftTime_Ass{iArea}');
    Mx{iArea,iClass} = nanmean(x);Sx{iArea,iClass} = nanstd(x);
    My{iArea,iClass} = nanmean(y);Sy{iArea,iClass} = nanstd(y);
    scatter(x,y,40,color_{iArea},'filled','MarkerFaceAlpha',0.5)
end
% for iArea = 1:3
%     errorbar(Mx{iArea},My{iArea},Sy{iArea},Sy{iArea},Sx{iArea},Sx{iArea},...
%         '-o','MarkerSize',5,...
%         'MarkerEdgeColor','k','MarkerFaceColor',color_{iArea},'Color','k')
% end
axis([-20 20 -2 2])
ylabel({'Change in';'Peak decoding'})
xlabel({'Change in span of';' significant decoding'})
%% Ass vs units shift
clear Mx My Sx Sy xAss yAss MxAss MyAss SxAss SyAss
iClass = 2;
figure
for iArea = 1:2
    subplot(1,3,iArea); hold on 
    plot([-20 20],[0 0],':k')
    plot([0 0],[-20 20],':k')
    title([Areas{iArea} ' Units/Assemblies'])
    y = cell2mat( shiftPeak{iArea,iClass}');
    x = cell2mat( shiftTime{iArea,iClass}');
    Mx{iArea} = nanmean(x);Sx{iArea} = nanstd(x);
    My{iArea} = nanmean(y);Sy{iArea} = nanstd(y);
    scatter(x,y,40,color_{iArea},'filled','MarkerFaceAlpha',0.3)
    
    yAss = cell2mat( shiftPeak_Ass{iArea}');
    xAss = cell2mat( shiftTime_Ass{iArea}');
    MxAss{iArea,iClass} = nanmean(x);SxAss{iArea} = nanstd(x);
    MyAss{iArea,iClass} = nanmean(y);SyAss{iArea} = nanstd(y);
    scatter(xAss,yAss,40,color_{iArea},'filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k')
    axis([-20 20 -2 2])
    switch iArea 
        case 1
            ylabel({'Change in';'Peak decoding'})
        case 2
            xlabel({'Change in span of';' significant decoding'})
    end
end

iArea = 3;
subplot(1,3,iArea); hold on
    title([Areas{iArea} ' Units/Assemblies'])

plot([-20 20],[0 0],':k')
plot([0 0],[-20 20],':k')

y = cell2mat( shiftPeak{1,3}');
x = cell2mat( shiftTime{1,3}');
Mx{iArea} = nanmean(x);Sx{iArea} = nanstd(x);
My{iArea} = nanmean(y);Sy{iArea} = nanstd(y);
scatter(x,y,40,'b','filled','MarkerFaceAlpha',0.3)

y = cell2mat( shiftPeak{2,3}');
x = cell2mat( shiftTime{2,3}');
Mx{iArea} = nanmean(x);Sx{iArea} = nanstd(x);
My{iArea} = nanmean(y);Sy{iArea} = nanstd(y);
scatter(x,y,40,'r','filled','MarkerFaceAlpha',0.3)

    yAss = cell2mat( shiftPeak_Ass{iArea}');
    xAss = cell2mat( shiftTime_Ass{iArea}');
    MxAss{iArea,iClass} = nanmean(x);SxAss{iArea} = nanstd(x);
    MyAss{iArea,iClass} = nanmean(y);SyAss{iArea} = nanstd(y);
    scatter(xAss,yAss,40,color_{iArea},'filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k')
    axis([-20 20 -2 2])
%% Peak, time sig histograms - collapsed
% Peak decoding
%
% x = 0:0.1:20;
% figure('name','Peak decoding');
% for iArea = 1:2
%     subplot(1,2,iArea); hold on
%     title(Areas{iArea})
%     TS_ = D.TScollapsed{iArea};
%     TSsig_ = D.TSsigcollapsed{iArea};
%     %     idx = nansum(TSsig_,2).*bw<=1;
%     %     TS_(idx,:)=[];
%     %     TSsig_(idx,:)=[];
%     for iClass = 1:3
%         idx =  D.membershipsigcollapsed{iArea} == (iClass-1);
%         temp = TS_(:,idx);
%         y = cumsum(histc(nanmax(temp),x));
%         y=y./max(y);
%         plot(x,y,'Color',col_{iClass});
%     end
% end
%
% %Time sig - collapsed
% x = 0:0.01:20;
% figure('name','Time significant');
% for iArea = 1:2
%     subplot(1,2,iArea); hold on
%     title(Areas{iArea})
%     TS_ = D.TScollapsed{iArea};
%     TSsig_ = double(D.TSsigcollapsed{iArea});
%     %     idx = nansum(TSsig_,2).*bw<=1;
%     %     TS_(idx,:)=[];
%     %     TSsig_(idx,:)=[];
%     for iClass = 1:3
%         idx =  D.membershipsigcollapsed{iArea} == (iClass-1);
%         temp = TSsig_(:,idx);
%         y = cumsum(histc(nansum(temp,1).*bw,x));
%         y=y./max(y);
%         plot(x,y,'Color',col_{iClass});
%     end
% end

%% PCA on joint error and correct data - T scores
clear inputC inputE
normalizePCA = false;
distanceMetric = 'euclidean';
softenNorm =  10;
for iArea = 1:3
    if iArea<3
        Membership_   =  D.membershipsigcollapsed{iArea};
        inputC{iArea} = (D.TScollapsed{iArea});
        inputE{iArea} = (D.TS_errcollapsed{iArea});
    else
        Membership_   =  [D.membershipsigcollapsed{1};D.membershipsigcollapsed{2}];
        inputC{iArea} = [D.TScollapsed{1},D.TScollapsed{2}];
        inputE{iArea} = [D.TS_errcollapsed{1},D.TS_errcollapsed{2}];
    end
    %         for i=1:size(inputC{iArea},2)
    %             inputC{iArea}(:,i)= mat2gray(inputC{iArea}(:,i));
    %             inputE{iArea}(:,i)= mat2gray(inputE{iArea}(:,i));
    %         end
    %          inputC{iArea}=zscore( inputC{iArea});
    %          inputE{iArea}=zscore( inputE{iArea});
    
    
    %     inputC{iArea}(isnan(inputC{iArea}) | isinf(inputC{iArea}))=0;
    %     inputE{iArea}(isnan(inputE{iArea}) | isinf(inputE{iArea}))=0;
    
    input = [inputC{iArea}; inputE{iArea}];
    input(isnan(input) | isinf(input))=0;
    N =length(Membership_);
    
    idx = []; 
    
%     if iArea<3
%         idx = find(nansum(D.TSsigcollapsed{iArea}).*bw<=1);
%     else
%         idx = find(nansum([D.TSsigcollapsed{1},D.TSsigcollapsed{2}]).*bw<=1);
%     end
    
    Membership_(idx) = [];
    if normalizePCA
        
        % Shared mean between correct and error
        
        %         ranges = range(input);
        %         normFactors = (ranges+softenNorm);
        %         input = bsxfun(@times, input, 1./normFactors);  % normalize
        %
        % %         suminput = 0;
        % %         suminput = suminput + bsxfun(@times, D.TScollapsed{iArea}, 1./normFactors);  % using the same normalization as above
        % %         suminput = suminput + bsxfun(@times, D.TS_errcollapsed{iArea}, 1./normFactors);  % using the same normalization as above
        % %         meaninput = suminput./2;
        % %         input = input - repmat(meaninput,2,1);
        %
        %         inputC{iArea} = input(1:N,:);
        %         inputE{iArea} = input(N+1:end,:);
        
        
        % Independent mean for correct and error
        ranges = range(inputC{iArea});
        normFactors = (ranges+softenNorm);
        inputC{iArea} = bsxfun(@times, inputC{iArea}, 1./normFactors);  % normalize
        
        ranges = range(inputE{iArea});
        normFactors = (ranges+softenNorm);
        inputE{iArea} = bsxfun(@times, inputE{iArea}, 1./normFactors);  % normalize
        %         suminput = 0;
        %         suminput = suminput + bsxfun(@times, D.TScollapsed{iArea}, 1./normFactors);  % using the same normalization as above
        %         suminput = suminput + bsxfun(@times, D.TS_errcollapsed{iArea}, 1./normFactors);  % using the same normalization as above
        %         meaninput = suminput./2;
        %         input = input - repmat(meaninput,2,1);
    else
        
%         inputC{iArea}=zscore( inputC{iArea});
%         inputE{iArea}=zscore( inputE{iArea})
        inputC{iArea}=inputC{iArea} - nanmean(inputC{iArea});
        inputE{iArea}=inputE{iArea} - nanmean(inputE{iArea});
        
%         inputC{iArea}(isnan(inputC{iArea}) | isinf(inputC{iArea}))=0;
%         inputE{iArea}(isnan(inputE{iArea}) | isinf(inputE{iArea}))=0;
    end
    
end
figure;

subplot(2,3,1);plot(inputC{1})
subplot(2,3,2);plot(inputC{2})
subplot(2,3,3);plot(inputE{3})
subplot(2,3,4);plot(inputE{1})
subplot(2,3,5);plot(inputE{2})
subplot(2,3,6);plot(inputE{3})
figure;
subplot(2,3,1);imagesc(inputC{1}')
subplot(2,3,2);imagesc(inputC{2}')
subplot(2,3,3);imagesc(inputC{3}')
subplot(2,3,4);imagesc(inputE{1}')
subplot(2,3,5);imagesc(inputE{2}')
subplot(2,3,6);imagesc(inputE{3}')
%% PCA on joint error and correct data - Time sig
clear inputC inputE
normalizePCA = true;
distanceMetric = 'hamming';
softenNorm = 10;
for iArea =1:3
    
     if iArea<3
        Membership_   =  D.membershipsigcollapsed{iArea};
        inputC{iArea} = double(D.TSsigcollapsed{iArea});
        inputE{iArea} = double(D.TS_errsigcollapsed{iArea});
    else
        Membership_   =  [D.membershipsigcollapsed{1};D.membershipsigcollapsed{2}];
        inputC{iArea} = double([D.TSsigcollapsed{1},D.TSsigcollapsed{2}]);
        inputE{iArea} = double([D.TS_errsigcollapsed{1},D.TS_errsigcollapsed{2}]);
     end

    
    %     for i=1:size(inputC{iArea},2)
    %         inputC{iArea}(:,i)= mat2gray(inputC{iArea}(:,i));
    %         inputE{iArea}(:,i)= mat2gray(inputE{iArea}(:,i));
    %     end
%     inputC{iArea}=zscore( inputC{iArea});
%     inputE{iArea}=zscore( inputE{iArea});
    
    
    inputC{iArea}(isnan(inputC{iArea}) | isinf(inputC{iArea}))=0;
    inputE{iArea}(isnan(inputE{iArea}) | isinf(inputE{iArea}))=0;
    
    input = [inputC{iArea}; inputE{iArea}];
    input(isnan(input) | isinf(input))=0;
    N =length(Membership_);
    idx = []; %idx = find(nansum(D.TSsigcollapsed{iArea}).*bw<=1);
    Membership_(idx) = [];
end
figure;

subplot(2,3,1);imagesc(inputC{1}')
subplot(2,3,2);imagesc(inputC{2}')
subplot(2,3,3);imagesc(inputC{3}')
subplot(2,3,4);imagesc(inputE{1}')
subplot(2,3,5);imagesc(inputE{2}')
subplot(2,3,6);imagesc(inputE{3}')
%% Pool and PCA  - correct dots and Assems
clear coeff score PCA_correct explained
for iArea = 1:3
    if iArea<3
        Membership_ = D.membershipsigcollapsed{iArea};
    else
        Membership_ = [D.membershipsigcollapsed{1};D.membershipsigcollapsed{2}];
    end
    input = inputC{iArea}; %  input(isnan(input) | isinf(input))=0;
    idx = nansum(input)==0;
    
    %     idx = find(nansum(D.TSsigcollapsed{iArea}).*bw<=1);
    Membership_(idx) = [];
    input(:,idx)=[];
    %     if normalizePCA
    %         input =zscore(input);
    %     end
%         for i=1:size(input,2)
%             input(:,i)= mat2gray(input(:,i));
%         end
        input(isnan(input) | isinf(input)) = 0;
    s_ =size(input);
    if s_(1)<s_(2)
        input = [input;input];
        [coeff{iArea},score{iArea},~,~,explained{iArea}] = pca(input);
        score{iArea}=score{iArea}(1:s_,:);
    else
        [coeff{iArea},score{iArea},~,~,explained{iArea}] = pca(input);
    end
    
    %     latent{iArea} =100*(latent{iArea}/sum(latent{iArea}));
    PCA_correct.score = score;
    PCA_correct.coeff = coeff;
    PCA_correct.explained = explained;
    PCA_correct.Membership_{iArea} = Membership_;
    PCA_correct.PD2{iArea} = squareform(pdist(mat2gray(input)',distanceMetric));
    idx_{iArea} = idx;    
end
%% Pool and PCA  - error dots and Assems
clear coeff score PCA_error explained
for iArea = 1:3
    if iArea<3
        Membership_ = D.membershipsigcollapsed{iArea};
    else
        Membership_ = [D.membershipsigcollapsed{1};D.membershipsigcollapsed{2}];
    end
    input = inputE{iArea}; %  input(isnan(input) | isinf(input))=0;
    idx = nansum(input)==0;
    %     idx = find(nansum(D.TSsigcollapsed{iArea}).*bw<=1);
    Membership_(idx) = [];
    input(:,idx)=[];
    %     if normalizePCA
    %         input =zscore(input);
    %     end
%         for i=1:size(input,2)
%             input(:,i)= mat2gray(input(:,i));
%         end
        input(isnan(input) | isinf(input)) = 0;
    
    s_ =size(input);
    if s_(1)<s_(2)
        input = [input;input];
        [coeff{iArea},score{iArea},~,~,explained{iArea}] = pca(input);
        score{iArea}=score{iArea}(1:s_,:);
    else
        [coeff{iArea},score{iArea},~,~,explained{iArea}] = pca(input);
    end
    %     latent{iArea} =100*(latent{iArea}/sum(latent{iArea}));
    coeff{iArea} = inputE{iArea}'/PCA_correct.score{iArea}';

    
    PCA_error.score = score;
    PCA_error.coeff = coeff;
    PCA_error.explained = explained;
    PCA_error.Membership_{iArea} = Membership_;
    PCA_error.PD2{iArea} = squareform(pdist(mat2gray(input)',distanceMetric));
    idx_{iArea} = idx;    
end
%% Pool and PCA  - Project errors into correct trials
clear coeff score PCA_correct PCA_error explained
for iArea = 1:3
    if iArea<3
        Membership_ = D.membershipsigcollapsed{iArea};
    else
        Membership_ = [D.membershipsigcollapsed{1};D.membershipsigcollapsed{2}];
    end
    inputE_ = inputE{iArea}; 
    inputC_ = inputC{iArea}; 
    idx = nansum([inputC_;inputE_])==0;
    Membership_(idx) = [];
    inputC_(:,idx)=[];inputC_(isnan(inputC_))=0;
    inputE_(:,idx)=[];inputE_(isnan(inputE_))=0;
    for i=1:size(inputC_,2)
        inputC_(:,i)= mat2gray(inputC_(:,i));
        inputE_(:,i)= mat2gray(inputE_(:,i));
    end
     
    [PCA_correct.score{iArea},PCA_correct.coeff{iArea},...
     PCA_error.score{iArea},PCA_error.coeff{iArea},PCA_correct.explained{iArea},PCA_error.explained{iArea}] =PCAeval(inputC_,inputE_);

    PCA_correct.Membership_{iArea} = Membership_;
    PCA_error.Membership_{iArea} = Membership_;
    PCA_correct.PD2{iArea} = squareform(pdist(mat2gray(inputC_)',distanceMetric));
    PCA_error.PD2{iArea} = squareform(pdist(mat2gray(inputE_)',distanceMetric));
    idx_{iArea} = idx;    
end
%%
no_PCs = 6;
A = [0,1,1;0,0,1;0,1,0;0,0,0;0,1,0;0,0,0];
B = [0,0,0;0,1,1;1,1,1;0,1,1;0,1,0;1,0,0];
figure
for iArea = 1:3
    subplot(1,3,iArea); hold on
    scoreC = PCA_correct.score{iArea}(:,1:no_PCs);
    scoreE = PCA_error.score{iArea}(:,1:no_PCs);
    for i=1:no_PCs
%        if max(zscore(scoreC(:,i)))<abs(min(zscore(scoreC(:,i))))
        if A(i,iArea)==1
           scoreC(:,i)=-scoreC(:,i);
%            scoreE(:,i)=scoreE(:,i);
        end
       
        if B(i,iArea)==1
            B(i,iArea)
            scoreE(:,i)=-scoreE(:,i);
        end
    end
    plot(staggerplot(scoreC,0,40),'k','LineWidth',1.5)
    plot(staggerplot(scoreE,0,40),':k','LineWidth',1.5)
    axis([0 402 4 Inf]); axis off
end
%% plot PCA bases
figure
for iArea = 1:3
    subplot(1,3,iArea); hold on
    plot(staggerplot(PCA_correct.score{iArea}(:,1:4),0,30),'b')
    plot(staggerplot(PCA_error.score{iArea}(:,1:4),0,30),'r')
end

no_PCs = 10;
cm = (hsv(no_PCs));
cm = flipud(jet(no_PCs));
figure
for iArea = 1:3
    subplot(1,3,iArea); hold on
    for iPC = 1:no_PCs
        PC = no_PCs +1 - iPC;
        PC = iPC;
        y = mat2gray( PCA_correct.score{iArea}(:,PC));
        x = (1:length(y)).*bw;
        PCcol_ = cm(PC,:);
        %         plot(x,y,'color',[cm(PC,:),0.1+0.01*PCA_correct.explained{iArea}(PC)],'LineWidth',2);
        %        h(iPC) = area(x,1.5*iPC+y,'FaceColor',PCcol_,'FaceAlpha',0.6);
        patch([min(x),x, max(x)],[1.2*iPC;1.2*iPC+y;1.2*iPC],PCcol_,'FaceAlpha',0.05+0.01*PCA_correct.explained{iArea}(PC));
    end
    
    
    axis([0 20 -0.5 Inf])
end
%% Var explained
nCs = 1:5;
figure
for iArea = 1:3
    subplot(1,3,iArea); hold on
    plot(nCs,cumsum(PCA_correct.explained{iArea}(nCs)),'k','LineWidth',2);
    plot(nCs,cumsum(PCA_error.explained{iArea}(nCs)),':k','LineWidth',2);
    axis([min(nCs) max(nCs) 0 101])
end
%% plot correct PCA
figure('name','Correct Trials')
for iArea = 1:3
    subplot(1,3,iArea); hold on
    for iClass = 1%:3
        x = PCA_correct.coeff{iArea}(1,PCA_correct.Membership_{iArea}==iClass-1);
        y = PCA_correct.coeff{iArea}(2,PCA_correct.Membership_{iArea}==iClass-1);
        scatter(x,y,20,col_{iClass},'filled','MarkerFaceAlpha',0.6);
    end
    
    
    if iArea < 3
        % Plot joint
        temp = D.JointMembersMatrixPadCollapse{iArea};
        temp(idx_{iArea},:)=[];
        for iAss = 1:size(temp,2)
            x_ = PCA_correct.coeff{iArea}(1,find(temp(:,iAss)));
            y_ = PCA_correct.coeff{iArea}(2,find(temp(:,iAss)));
            for i=1:length(x_)
                for j=1:length(x_)
                    if i~=j
                        plot(x_([i,j]),y_([i,j]),'Color',col_{3})
                        scatter(x_([i,j]),y_([i,j]),20,col_{3},'filled','MarkerFaceAlpha',0.6);
                    end
                end
            end
        end
        
        % Plot local
        temp = D.LocalMembersMatrixPadCollapse{iArea};
        temp(idx_{iArea},:)=[];
        for iAss = 1:size(temp,2)
            x_ = PCA_correct.coeff{iArea}(1,find(temp(:,iAss)));
            y_ = PCA_correct.coeff{iArea}(2,find(temp(:,iAss)));
            for i=1:length(x_)
                for j=1:length(x_)
                    if i~=j
                        plot(x_([i,j]),y_([i,j]),'Color',col_{2})
                        scatter(x_([i,j]),y_([i,j]),20,col_{2},'filled','MarkerFaceAlpha',0.6);
                    end
                end
            end
        end
    else
        % Plot joint
        temp = D.JointMembersMatrixPadCollapse{iArea};
        temp(idx_{iArea},:)=[];
        for iAss = 1:size(temp,2)
            x_ = PCA_correct.coeff{iArea}(1,find(temp(:,iAss)));
            y_ = PCA_correct.coeff{iArea}(2,find(temp(:,iAss)));
            for i=1:length(x_)
                for j=1:length(x_)
                    if i~=j
                        plot(x_([i,j]),y_([i,j]),'Color',col_{3})
                        scatter(x_([i,j]),y_([i,j]),20,col_{3},'filled','MarkerFaceAlpha',0.6);
                    end
                end
            end
        end
        
    end
    
%         axis([-1 1 -1 1]);
        axis([-0.4 0.4 -0.4 0.5]);
    axis square; grid on
    %     clear idx x y x_ y_ iClass
    title(Areas{iArea})
    xlabel('PC1')
    ylabel('PC2')
end
%% plot error PCA
figure('name','Error Trials')
for iArea = 1:3
    subplot(1,3,iArea); hold on
    for iClass = 1%:3
        x = PCA_error.coeff{iArea}(1,PCA_error.Membership_{iArea}==iClass-1);
        y = PCA_error.coeff{iArea}(2,PCA_error.Membership_{iArea}==iClass-1);
        scatter(x,y,20,col_{iClass},'filled','MarkerFaceAlpha',0.6);
    end
    
     if iArea < 3
        % Plot joint
        temp = D.JointMembersMatrixPadCollapse{iArea};
        temp(idx_{iArea},:)=[];
        for iAss = 1:size(temp,2)
            x_ = PCA_error.coeff{iArea}(1,find(temp(:,iAss)));
            y_ = PCA_error.coeff{iArea}(2,find(temp(:,iAss)));
            for i=1:length(x_)
                for j=1:length(x_)
                    if i~=j
                        plot(x_([i,j]),y_([i,j]),'Color',col_{3})
                        scatter(x_([i,j]),y_([i,j]),20,col_{3},'filled','MarkerFaceAlpha',0.6);
                    end
                end
            end
        end
        
        % Plot local
        temp = D.LocalMembersMatrixPadCollapse{iArea};
        temp(idx_{iArea},:)=[];
        for iAss = 1:size(temp,2)
            x_ = PCA_error.coeff{iArea}(1,find(temp(:,iAss)));
            y_ = PCA_error.coeff{iArea}(2,find(temp(:,iAss)));
            for i=1:length(x_)
                for j=1:length(x_)
                    if i~=j
                        plot(x_([i,j]),y_([i,j]),'Color',col_{2})
                        scatter(x_([i,j]),y_([i,j]),20,col_{2},'filled','MarkerFaceAlpha',0.6);
                    end
                end
            end
        end
        
    else
        % Plot joint
        temp = D.JointMembersMatrixPadCollapse{iArea};
        temp(idx_{iArea},:)=[];
        for iAss = 1:size(temp,2)
            x_ = PCA_error.coeff{iArea}(1,find(temp(:,iAss)));
            y_ = PCA_error.coeff{iArea}(2,find(temp(:,iAss)));
            for i=1:length(x_)
                for j=1:length(x_)
                    if i~=j
                        plot(x_([i,j]),y_([i,j]),'Color',col_{3})
                        scatter(x_([i,j]),y_([i,j]),20,col_{3},'filled','MarkerFaceAlpha',0.6);
                    end
                end
            end
        end
        
    end
    
        axis([-1 1 -1 1]);
%         axis([-0.4 0.4 -0.4 0.5]);
    rectangle('Position',[-0.4 -0.4 0.8 0.9])
    axis square; grid on
    %     clear idx x y x_ y_ iClass
    title(Areas{iArea})
    xlabel('PC1')
    ylabel('PC2')
end
%% t-SNE - correct dots and Assems
figure
for iArea = 1:2
    Membership_ = D.membershipsigcollapsed{iArea};
    input = (inputC{iArea});   input(isnan(input) | isinf(input))=0;
    %     for i=1:size(input,2)
    %         input(:,i)= mat2gray(input(:,i));
    %     end
    %     input = zscore(input);
    %     ydata = tsne(input',[],2,2,5);%size(input,2)
    ydata = tsne(input','Algorithm','exact','Distance','hamming');
    
    
    subplot(1,2,iArea); hold on
    for iClass = 1%:3
        x = ydata(Membership_==iClass-1,1);
        y = ydata(Membership_==iClass-1,2);
        scatter(x,y,20,col_{iClass},'filled','MarkerFaceAlpha',0.6);
    end
    
    %     Plot local
    for iAss = 1:size(D.LocalMembersMatrixPadCollapse{iArea},2)
        x_ =  ydata(find(D.LocalMembersMatrixPadCollapse{iArea}(:,iAss)),1);
        y_ =  ydata(find(D.LocalMembersMatrixPadCollapse{iArea}(:,iAss)),2);
        for i=1:length(x_)
            for j=1:length(x_)
                if i~=j
                    plot(x_([i,j]),y_([i,j]),'Color',[col_{2},0.1])
                    scatter(x_([i,j]),y_([i,j]),20,col_{2},'filled','MarkerFaceAlpha',0.6);
                end
            end
        end
    end
    % Plot joint
    for iAss = 1:size(D.JointMembersMatrixPadCollapse{iArea},2)
        x_ =  ydata(find(D.JointMembersMatrixPadCollapse{iArea}(:,iAss)),1);
        y_ =  ydata(find(D.JointMembersMatrixPadCollapse{iArea}(:,iAss)),2);
        for i=1:length(x_)
            for j=1:length(x_)
                if i~=j
                    plot(x_([i,j]),y_([i,j]),'Color',[col_{3},0.1])
                    scatter(x_([i,j]),y_([i,j]),20,col_{3},'filled','MarkerFaceAlpha',0.6);
                end
            end
        end
    end
    
    
    %             axis([-40 40 -40 40]);
    %     axis([-0.4 0.4 -0.4 0.5]);
    axis square; grid on
    %     clear idx x y x_ y_ iClass
    title(Areas{iArea})
    xlabel('PC1')
    ylabel('PC2')
    
end
%% t-SNE - error dots and Assems
figure
for iArea = 1:2
    Membership_ = D.membershipsigcollapsed{iArea};
    input = (inputE{iArea});   input(isnan(input) | isinf(input))=0;
    %     for i=1:size(input,2)
    %         input(:,i)= mat2gray(input(:,i));
    %     end
    %     input = zscore(input);
    %     ydata = tsne(input',[],2,2,5);%size(input,2)
    ydata = tsne(input','Algorithm','exact','Distance','hamming');
    
    
    subplot(1,2,iArea); hold on
    for iClass = 1%:3
        x = ydata(Membership_==iClass-1,1);
        y = ydata(Membership_==iClass-1,2);
        scatter(x,y,20,col_{iClass},'filled','MarkerFaceAlpha',0.6);
    end
    
    %     Plot local
    for iAss = 1:size(D.LocalMembersMatrixPadCollapse{iArea},2)
        x_ =  ydata(find(D.LocalMembersMatrixPadCollapse{iArea}(:,iAss)),1);
        y_ =  ydata(find(D.LocalMembersMatrixPadCollapse{iArea}(:,iAss)),2);
        for i=1:length(x_)
            for j=1:length(x_)
                if i~=j
                    plot(x_([i,j]),y_([i,j]),'Color',[col_{2},0.1])
                    scatter(x_([i,j]),y_([i,j]),20,col_{2},'filled','MarkerFaceAlpha',0.6);
                end
            end
        end
    end
    % Plot joint
    for iAss = 1:size(D.JointMembersMatrixPadCollapse{iArea},2)
        x_ =  ydata(find(D.JointMembersMatrixPadCollapse{iArea}(:,iAss)),1);
        y_ =  ydata(find(D.JointMembersMatrixPadCollapse{iArea}(:,iAss)),2);
        for i=1:length(x_)
            for j=1:length(x_)
                if i~=j
                    plot(x_([i,j]),y_([i,j]),'Color',[col_{3},0.1])
                    scatter(x_([i,j]),y_([i,j]),20,col_{3},'filled','MarkerFaceAlpha',0.6);
                end
            end
        end
    end
    
    
    %             axis([-40 40 -40 40]);
    %     axis([-0.4 0.4 -0.4 0.5]);
    axis square; grid on
    %     clear idx x y x_ y_ iClass
    title(Areas{iArea})
    xlabel('PC1')
    ylabel('PC2')
    
end

%% Distribution of pairwise PC loadings between member units - correct
clear Pdist_
try
    D = rmfield(D,'PDist');
end
bins = 0:0.01:1;
for iArea = 1:2
%     coeff = triu(squareform(pdist(inputC{iArea}')));
    coeff      = PCA_correct.coeff{iArea};
    Membership = PCA_correct.Membership_{iArea};
    find(idx_{iArea})
    coeff(idx_{iArea},:)=[];coeff(:,idx_{iArea})=[];
    Membership(idx_{iArea})=[];
    for iFile = 1:length(D.nUnits{iArea})
        
        nUnits = D.nUnits{iArea}(iFile);
     
        
        coeff_ = coeff(1:nUnits,1:nUnits); %[eigenvector/component, eigenvalue]
        Membership_ = Membership(1:nUnits);
        JointMembers_ = D.JointMembersMatrix{iArea}{iFile};
        LocalMembers_ = D.LocalMembersMatrix{iArea}{iFile};
        
        % Strip from remaining archive
        coeff =  coeff(nUnits+1:end,nUnits+1:end);
        Membership(1:nUnits) = [];
        
        % non-members
        try
            x = coeff_(1,Membership_==0);
            y = coeff_(2,Membership_==0);
            dist_ = triu(distmat([x;y]'));
            Pdist_ = dist_(dist_>0);
            y_ = histc(Pdist_,bins);
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Corr_NonMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Corr_NonMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
        % local members
        try
            x = coeff_(1,Membership_==1);
            y = coeff_(2,Membership_==1);
            dist_ = triu(distmat([x;y]'));
            Pdist_ = dist_(dist_>0);
            y_ = histc(Pdist_,bins);
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Corr_LocalMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Corr_LocalMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
        % joint members
        try
            x = coeff_(1,Membership_==2);
            y = coeff_(2,Membership_==2);
            dist_ = triu(distmat([x;y]'));
            Pdist_ = dist_(dist_>0);
            y_ = histc(Pdist_,bins);
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Corr_JointMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Corr_JointMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
    end
    
    D.PDist.Corr_NonMembersMEAN{iArea}   = nanmean(D.PDist.Corr_NonMembers{iArea},2);
    D.PDist.Corr_NonMembersSEM{iArea}    = nansem(D.PDist.Corr_NonMembers{iArea},2);
    D.PDist.Corr_LocalMembersMEAN{iArea} = nanmean(D.PDist.Corr_LocalMembers{iArea},2);
    D.PDist.Corr_LocalMembersSEM{iArea}  = nansem(D.PDist.Corr_LocalMembers{iArea},2);
    D.PDist.Corr_JointMembersMEAN{iArea} = nanmean(D.PDist.Corr_JointMembers{iArea},2);
    D.PDist.Corr_JointMembersSEM{iArea}  = nansem(D.PDist.Corr_JointMembers{iArea},2);
    
end

figure;
for iArea = 1:2
    subplot(1,2,iArea); hold on
    M = D.PDist.Corr_NonMembersMEAN{iArea};
    E = D.PDist.Corr_NonMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{1})
    
    M = D.PDist.Corr_LocalMembersMEAN{iArea};
    E = D.PDist.Corr_LocalMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{2})
    
    M = D.PDist.Corr_JointMembersMEAN{iArea};
    E = D.PDist.Corr_JointMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{3})
    axis([min(bins) max(bins) 0 1])
    
end
% plot as hists from all cross-animals
clear Pdist_
figure
for iArea = 1:2
    subplot(1,2,iArea); hold on
    for iClass = 1:3
        coeff = PCA_correct.coeff{iArea};
        Membership_ = PCA_correct.Membership_{iArea};
        PD2 = PCA_correct.PD2{iArea};
        x = coeff(1,Membership_==iClass-1);
        y = coeff(2,Membership_==iClass-1);
        %         x = PD2(1,Membership_==iClass-1);
        %         y = PD2(2,Membership_==iClass-1);
        dist_ = triu(distmat([x;y]'));
        Pdist_{iArea}{iClass} = dist_(dist_>0);
        y_ = histc(Pdist_{iArea}{iClass},bins);
        y_ = y_ ./nansum(y_);
        plot(bins,y_,'color',col_{iClass})
        
        
    end
end
%% Distribution of pairwise PC loadings between member units - error

clear Pdist_
try
    D = rmfield(D,'PDist');
end
bins = 0:0.01:1;
for iArea = 1:2
%     coeff = triu(squareform(pdist(inputE{iArea}')));
    coeff      = PCA_error.coeff{iArea};
    Membership = PCA_error.Membership_{iArea};
    
    for iFile = 1:length(D.nUnits{iArea})
        
        nUnits = D.nUnits{iArea}(iFile);
        coeff_ = coeff(1:nUnits,1:nUnits); %[eigenvector/component, eigenvalue]
        Membership_ = Membership(1:nUnits);
        JointMembers_ = D.JointMembersMatrix{iArea}{iFile};
        LocalMembers_ = D.LocalMembersMatrix{iArea}{iFile};
        
        % Strip from remaining archive
        coeff =  coeff(nUnits+1:end,nUnits+1:end);
        Membership(1:nUnits) = [];
        
        % non-members
        try
            x = coeff_(1,Membership_==0);
            y = coeff_(2,Membership_==0);
            dist_ = triu(distmat([x;y]'));
            Pdist_ = dist_(dist_>0);
            y_ = histc(Pdist_,bins);
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Err_NonMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Err_NonMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
        % local members
        try
            x = coeff_(1,Membership_==1);
            y = coeff_(2,Membership_==1);
            dist_ = triu(distmat([x;y]'));
            Pdist_ = dist_(dist_>0);
            y_ = histc(Pdist_,bins);
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Err_LocalMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Err_LocalMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
        % joint members
        try
            x = coeff_(1,Membership_==2);
            y = coeff_(2,Membership_==2);
            dist_ = triu(distmat([x;y]'));
            Pdist_ = dist_(dist_>0);
            y_ = histc(Pdist_,bins);
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Err_JointMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Err_JointMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
    end
    
    D.PDist.Err_NonMembersMEAN{iArea}   = nanmean(D.PDist.Err_NonMembers{iArea},2);
    D.PDist.Err_NonMembersSEM{iArea}    = nansem(D.PDist.Err_NonMembers{iArea},2);
    D.PDist.Err_LocalMembersMEAN{iArea} = nanmean(D.PDist.Err_LocalMembers{iArea},2);
    D.PDist.Err_LocalMembersSEM{iArea}  = nansem(D.PDist.Err_LocalMembers{iArea},2);
    D.PDist.Err_JointMembersMEAN{iArea} = nanmean(D.PDist.Err_JointMembers{iArea},2);
    D.PDist.Err_JointMembersSEM{iArea}  = nansem(D.PDist.Err_JointMembers{iArea},2);
    
end

figure;
for iArea = 1:2
    subplot(1,2,iArea); hold on
    M = D.PDist.Err_NonMembersMEAN{iArea};
    E = D.PDist.Err_NonMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{1})
    
    M = D.PDist.Err_LocalMembersMEAN{iArea};
    E = D.PDist.Err_LocalMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{2})
    
    M = D.PDist.Err_JointMembersMEAN{iArea};
    E = D.PDist.Err_JointMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{3})
    axis([min(bins) max(bins) 0 1])
end
% plot as hists from all cross-animals
clear Pdist_
figure
for iArea = 1:2
    subplot(1,2,iArea); hold on
    for iClass = 1:3
        coeff = PCA_error.coeff{iArea};
        Membership_ = PCA_error.Membership_{iArea};
        PD2 = PCA_error.PD2{iArea};
        x = coeff(1,Membership_==iClass-1);
        y = coeff(2,Membership_==iClass-1);
        %         x = PD2(1,Membership_==iClass-1);
        %         y = PD2(2,Membership_==iClass-1);
        dist_ = triu(distmat([x;y]'));
        Pdist_{iArea}{iClass} = dist_(dist_>0);
        y_ = histc(Pdist_{iArea}{iClass},bins);
        y_ = y_ ./nansum(y_);
        plot(bins,y_,'color',col_{iClass})
        
        
    end
end
%% Distribution of pairwise distances PC loadings between member units - correct vs error
clear Pdist_
bins = 0:0.001:0.25;
for iArea = 1:2
    figure('Name',Areas{iArea}) ;
    
    
    
    for iClass = 1:3
        subplot(1,3,iClass); hold on
        title(MemberClasses_{iClass})
        coeff       = PCA_correct.coeff{iArea};
        Membership_ = PCA_correct.Membership_{iArea};
        x = coeff(1,Membership_==iClass-1);
        y = coeff(2,Membership_==iClass-1);
        dist_ = triu(distmat([x;y]'));
        Pdist_ = dist_(dist_>0);
        y_ = cumsum(histc(Pdist_,bins));
        y_ = y_ ./max(y_);
        plot(bins,y_,'color',col_{iClass},'LineWidth',1.5)
        
        coeff       = PCA_error.coeff{iArea};
        Membership_ = PCA_error.Membership_{iArea};
        x = coeff(1,Membership_==iClass-1);
        y = coeff(2,Membership_==iClass-1);
        dist_ = triu(distmat([x;y]'));
        Pdist_ = dist_(dist_>0);
        y_ = cumsum(histc(Pdist_,bins));
        y_ = y_ ./max(y_);
        plot(bins,y_,'color',col_{iClass},'LineWidth',1.5,'LineStyle',':')
        axis tight
    end
end
%% Distribution of pairwise distances between D' scores - correct vs error
noDraws = 5;
noBS = 1000;
clear Pdist_
try
    D = rmfield(D,'PDist');
end
bins = 0:1:40;
for iArea = 1:2
    coeff = triu(squareform(pdist(inputC{iArea}',distanceMetric)));
    Membership = PCA_correct.Membership_{iArea};
    
    for iFile = 1:length(D.nUnits{iArea})
        
        nUnits = D.nUnits{iArea}(iFile);
        coeff_ = coeff(1:nUnits,1:nUnits); %[eigenvector/component, eigenvalue]
        Membership_ = Membership(1:nUnits);
        JointMembers_ = D.JointMembersMatrix{iArea}{iFile};
        LocalMembers_ = D.LocalMembersMatrix{iArea}{iFile};
        
        % Strip from remaining archive
        coeff =  coeff(nUnits+1:end,nUnits+1:end);
        Membership(1:nUnits) = [];
        
        % non-members
        try
            idx = find(Membership_==0);
            PDist_ = [];
            for i = 1:length(idx)
                for j =1:length(idx)
                    if i~=j & i<j
                        PDist_=[PDist_;coeff_(idx(i),idx(j))];
                    end
                end
            end
            
            y_ = zeros(size(bins));
            for bs=1:noBS
                y_ = y_ + histc(PDist_(randsample(1:length(PDist_),noDraws)),bins)';
            end
            y_ = y_./noBS;
            
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Corr_NonMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Corr_NonMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
        % local members
        try
            idx = find(Membership_==1);
            PDist_ = [];
            for i = 1:length(idx)
                for j =1:length(idx)
                    if i~=j & i<j
                        PDist_=[PDist_;coeff_(idx(i),idx(j))];
                    end
                end
            end
            y_ = zeros(size(bins));
            for bs=1:noBS
                y_ = y_ + histc(PDist_(randsample(1:length(PDist_),noDraws)),bins)';
            end
            y_ = y_./noBS;
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Corr_LocalMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Corr_LocalMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
        % joint members
        try
            idx = find(Membership_==2);
            PDist_ = [];
            for i = 1:length(idx)
                for j =1:length(idx)
                    if i~=j & i<j
                        PDist_=[PDist_;coeff_(idx(i),idx(j))];
                    end
                end
            end
            y_ = zeros(size(bins));
            for bs=1:noBS
                y_ = y_ + histc(PDist_(randsample(1:length(PDist_),noDraws)),bins)';
            end
            y_ = y_./noBS;
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Corr_JointMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Corr_JointMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
    end
    
    D.PDist.Corr_NonMembersMEAN{iArea}   = nanmean(D.PDist.Corr_NonMembers{iArea},2);
    D.PDist.Corr_NonMembersSEM{iArea}    = nansem(D.PDist.Corr_NonMembers{iArea},2);
    D.PDist.Corr_LocalMembersMEAN{iArea} = nanmean(D.PDist.Corr_LocalMembers{iArea},2);
    D.PDist.Corr_LocalMembersSEM{iArea}  = nansem(D.PDist.Corr_LocalMembers{iArea},2);
    D.PDist.Corr_JointMembersMEAN{iArea} = nanmean(D.PDist.Corr_JointMembers{iArea},2);
    D.PDist.Corr_JointMembersSEM{iArea}  = nansem(D.PDist.Corr_JointMembers{iArea},2);
    
end
for iArea = 1:2
    coeff = triu(squareform(pdist(inputE{iArea}',distanceMetric)));
    Membership = PCA_error.Membership_{iArea};
    
    for iFile = 1:length(D.nUnits{iArea})
        
        nUnits = D.nUnits{iArea}(iFile);
        coeff_ = coeff(1:nUnits,1:nUnits); %[eigenvector/component, eigenvalue]
        Membership_ = Membership(1:nUnits);
        JointMembers_ = D.JointMembersMatrix{iArea}{iFile};
        LocalMembers_ = D.LocalMembersMatrix{iArea}{iFile};
        
        % Strip from remaining archive
        coeff =  coeff(nUnits+1:end,nUnits+1:end);
        Membership(1:nUnits) = [];
        
        % non-members
        try
            idx = find(Membership_==0);
            PDist_ = [];
            for i = 1:length(idx)
                for j =1:length(idx)
                    if i~=j & i<j
                        PDist_=[PDist_;coeff_(idx(i),idx(j))];
                    end
                end
            end
            
            y_ = zeros(size(bins));
            for bs=1:noBS
                y_ = y_ + histc(PDist_(randsample(1:length(PDist_),noDraws)),bins)';
            end
            y_ = y_./noBS;
            
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Err_NonMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Err_NonMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
        % local members
        try
            idx = find(Membership_==1);
            PDist_ = [];
            for i = 1:length(idx)
                for j =1:length(idx)
                    if i~=j & i<j
                        PDist_=[PDist_;coeff_(idx(i),idx(j))];
                    end
                end
            end
            y_ = zeros(size(bins));
            for bs=1:noBS
                y_ = y_ + histc(PDist_(randsample(1:length(PDist_),noDraws)),bins)';
            end
            y_ = y_./noBS;
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Err_LocalMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Err_LocalMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
        % joint members
        try
            idx = find(Membership_==2);
            PDist_ = [];
            for i = 1:length(idx)
                for j =1:length(idx)
                    if i~=j & i<j
                        PDist_=[PDist_;coeff_(idx(i),idx(j))];
                    end
                end
            end
            y_ = zeros(size(bins));
            for bs=1:noBS
                y_ = y_ + histc(PDist_(randsample(1:length(PDist_),noDraws)),bins)';
            end
            y_ = y_./noBS;
            y_ = cumsum(y_ ./nansum(y_));
            D.PDist.Err_JointMembers{iArea}(:,iFile) = y_;
        catch
            D.PDist.Err_JointMembers{iArea}(:,iFile) = nan(size(bins));
        end
        
    end
    
    D.PDist.Err_NonMembersMEAN{iArea}   = nanmean(D.PDist.Err_NonMembers{iArea},2);
    D.PDist.Err_NonMembersSEM{iArea}    = nansem(D.PDist.Err_NonMembers{iArea},2);
    D.PDist.Err_LocalMembersMEAN{iArea} = nanmean(D.PDist.Err_LocalMembers{iArea},2);
    D.PDist.Err_LocalMembersSEM{iArea}  = nansem(D.PDist.Err_LocalMembers{iArea},2);
    D.PDist.Err_JointMembersMEAN{iArea} = nanmean(D.PDist.Err_JointMembers{iArea},2);
    D.PDist.Err_JointMembersSEM{iArea}  = nansem(D.PDist.Err_JointMembers{iArea},2);
    
end
for iArea = 1:2
    figure;

    subplot(1,3,1); hold on
    M = D.PDist.Corr_NonMembersMEAN{iArea};
    E = D.PDist.Corr_NonMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{1},1)

    M = D.PDist.Err_NonMembersMEAN{iArea};
    E = D.PDist.Err_NonMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{1})
    
    subplot(1,3,2); hold on
    M = D.PDist.Corr_LocalMembersMEAN{iArea};
    E = D.PDist.Corr_LocalMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{2},1)
    M = D.PDist.Err_LocalMembersMEAN{iArea};
    E = D.PDist.Err_LocalMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{2})
    
    subplot(1,3,3); hold on
    M = D.PDist.Corr_JointMembersMEAN{iArea};
    E = D.PDist.Corr_JointMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{3},1)
    
    M = D.PDist.Err_JointMembersMEAN{iArea};
    E = D.PDist.Err_JointMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{3})
    axis([min(bins) max(bins) 0 1])
end
figure;

for iArea = 1:2
    subplot(1,2,iArea); hold on
    M = D.PDist.Err_NonMembersMEAN{iArea};
    E = D.PDist.Err_NonMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{1})
    
    M = D.PDist.Err_LocalMembersMEAN{iArea};
    E = D.PDist.Err_LocalMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{2})
    
    M = D.PDist.Err_JointMembersMEAN{iArea};
    E = D.PDist.Err_JointMembersSEM{iArea};
    ciplot(M+E,M-E,bins,col_{3})
    axis([min(bins) max(bins) 0 1])
end

%% plot sorted independently - correct
clear idx_
figure;
for iArea = 1:2
    subplot(1,2,iArea); hold on
    sorted_T = [];
    sorted_Tsig = [];
    for iClass=1:3
        
        TS_ =    D.TScollapsed{iArea}(:,D.membershipsigcollapsed{iArea}==iClass-1)';
        TSsig_ = D.TSsigcollapsed{iArea}(:,D.membershipsigcollapsed{iArea}==iClass-1)';
        FracSig{iArea}(iClass,:) = sum(TSsig_)./size(TSsig_,1);
        [~,idx] = max(TS_,[],2);
        [~,idx] = sort(idx);
        idx_{iArea}{iClass} = idx; % to maintain sort order across correct/error conditions
        TS_ = TS_(idx,:);
        TSsig_ = TSsig_ (idx,:);
        %
        %         idx = nansum(TSsig_,2).*bw<=1;
        %         TS_(idx,:)=[];
        %         TSsig_(idx,:)=[];
        %noMembers(iClass) = sum(D.membershipsigcollapsed{iArea}==iClass-1);
        noMembers(iClass) = size(TSsig_,1);
        
        sorted_T = [sorted_T;TS_];
        sorted_Tsig = [sorted_Tsig;TSsig_];
        
    end
    sorted_T(~sorted_Tsig )= 0;
    sorted_TCorrect{iArea} = sorted_T;
    x = (1:size(sorted_T,2)).*bw;
    y = 1:size(sorted_T,1);
    x_  = repmat(x,size(sorted_T,1),1);
    y_  = repmat(y,size(sorted_T,2),1);
    imagesc(x_(:),y_(:),sorted_T);
    set(gca,'YDir','normal')
    %         cmap =([1 1 1;flipud(hot)]);
    cmap =([1 1 1;flipud(gray)]);
    colormap (cmap)
    caxis([0 1])
    plot([10 10],[1 size(sorted_T,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(sorted_T,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(sorted_T,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    rectangle('Position',[min(x), 0, max(x), noMembers(1)],'EdgeColor',col_{1}','LineWidth',2)
    rectangle('Position',[min(x), noMembers(1)+1, max(x), noMembers(2)],'EdgeColor',col_{2}','LineWidth',2)
    rectangle('Position',[min(x), sum(noMembers(1:2))+2, max(x), noMembers(3)],'EdgeColor',col_{3}','LineWidth',2)
    %                   axis([0 20 0 350])
    axis tight
    
end

% figure;
% for iArea = 1:2
%     subplot(1,2,iArea); hold on
%     for iClass = 1:3
%         plot(x,FracSig{iArea}(iClass,:),'color',col_{iClass},'LineWidth',1.5)
%     end
%     plot([10 10],[0 1],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[0 1],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[0 1],'color',[1 0 0 0.6],'LineWidth',1.5)
%     axis tight
% end
%% plot sorted independently - error

figure;
for iArea = 1:2
    subplot(1,2,iArea); hold on
    sorted_T = [];
    sorted_Tsig = [];
    for iClass=1:3
        
        TS_ =    D.TS_errcollapsed{iArea}(:,D.membershipsigcollapsed{iArea}==iClass-1)';
        TSsig_ = D.TS_errsigcollapsed {iArea}(:,D.membershipsigcollapsed{iArea}==iClass-1)';
        FracSig_err{iArea}(iClass,:) = sum(TSsig_)./size(TSsig_,1);
        
        %         [~,idx] = max(TS_,[],2);
        %         [~,idx] = sort(idx);
        idx = idx_{iArea}{iClass};
        TS_ = TS_(idx,:);
        TSsig_ = TSsig_ (idx,:);
        
        %         idx = nansum(TSsig_,2).*bw<=0;
        %         TS_(idx,:)=[];
        %         TSsig_(idx,:)=[];
        %noMembers(iClass) = sum(D.membershipsigcollapsed{iArea}==iClass-1);
        noMembers(iClass) = size(TSsig_,1);
        
        sorted_T = [sorted_T;TS_];
        sorted_Tsig = [sorted_Tsig;TSsig_];
        
    end
    sorted_T(~sorted_Tsig )=0;
    sorted_TError{iArea} = sorted_T;
    x = (1:size(sorted_T,2)).*bw;
    y = 1:size(sorted_T,1);
    x_  = repmat(x,size(sorted_T,1),1);
    y_  = repmat(y,size(sorted_T,2),1);
    imagesc(x_(:),y_(:),sorted_T);
    set(gca,'YDir','normal')
    %         cmap =([1 1 1;flipud(hot)]);
    cmap =([1 1 1;flipud(gray)]);
    colormap (cmap)
    caxis([0 1])
    plot([10 10],[1 size(sorted_T,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(sorted_T,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(sorted_T,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    rectangle('Position',[min(x), 0, max(x), noMembers(1)],'EdgeColor',col_{1}','LineWidth',2)
    rectangle('Position',[min(x), noMembers(1)+1, max(x), noMembers(2)],'EdgeColor',col_{2}','LineWidth',2)
    rectangle('Position',[min(x), sum(noMembers(1:2))+2, max(x), noMembers(3)],'EdgeColor',col_{3}','LineWidth',2)
    %                   axis([0 20 0 350])
    axis tight
end
% figure;
% for iArea = 1:2
%     subplot(1,2,iArea); hold on
%     for iClass = 1:3
%         plot(x,FracSig_err{iArea}(iClass,:),'color',col_{iClass},'LineWidth',1.5)
%     end
%     plot([10 10],[0 1],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[0 1],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[0 1],'color',[1 0 0 0.6],'LineWidth',1.5)
%     axis tight
% end
%% plot sorted independently - correct vs error
figure;
for iArea = 1:2
    subplot(1,2,iArea); hold on
    sorted_T = [];
    sorted_Tmean{iArea} = [];
    
    for iClass=1:3
        %         tempC = (D.TS_errcollapsed{iArea}-D.TScollapsed{iArea});
        %         tempC = (D.TS_errcollapsed{iArea}-D.TScollapsed{iArea}) ./ (D.TScollapsed{iArea}+D.TS_errcollapsed{iArea});
        
        C =  (D.TSsigcollapsed{iArea});
        E =  (D.TS_errsigcollapsed{iArea});
        tempC = zeros(size(C));
        tempC(~C & E)=1;
        tempC(C & ~E)=-1;
        tempC(C & E) = 0;
        
        tempCsig = D.TSsigcollapsed{iArea} & D.TS_errsigcollapsed{iArea};
        %         tempC(~tempCsig) = 0;
        TS_  =    tempC(:,D.membershipsigcollapsed{iArea}==iClass-1)';
        TS__ =    D.TScollapsed{iArea}(:,D.membershipsigcollapsed{iArea}==iClass-1)';
        
        %         [~,idx] = max(TS__,[],2);
        %         [~,idx] = sort(idx);
        idx = idx_{iArea}{iClass};
        
        TS_ = TS_(idx,:);
        noMembers(iClass) = size(TS__,1);
        sorted_T = [sorted_T;TS_];
        sorted_Tmean{iArea} = [sorted_Tmean{iArea}; nanmean(TS_)];
        
    end
    sorted_T(isnan(sorted_T))=0;
    x = (1:size(sorted_T,2)).*bw;
    y = 1:size(sorted_T,1);
    x_  = repmat(x,size(sorted_T,1),1);
    y_  = repmat(y,size(sorted_T,2),1);
    imagesc(x_(:),y_(:),sorted_T);
    set(gca,'YDir','normal')
    %     colormap(redbluecmap);
    colormap([1 0 0; 1 1 1 ;0 0 1]);
    caxis([-1 1])
    plot([10 10],[1 size(sorted_T,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(sorted_T,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(sorted_T,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    rectangle('Position',[min(x), 0, max(x), noMembers(1)],'EdgeColor',col_{1}','LineWidth',2)
    rectangle('Position',[min(x), noMembers(1)+1, max(x), noMembers(2)],'EdgeColor',col_{2}','LineWidth',2)
    rectangle('Position',[min(x), sum(noMembers(1:2))+2, max(x), noMembers(3)],'EdgeColor',col_{3}','LineWidth',2)
    %                   axis([0 20 0 350])
    axis tight
    
end

figure;
for iArea=1:2
    subplot(1,2,iArea); hold on
    for iClass=1:3
        plot(x,sorted_Tmean{iArea}(iClass,:),'color',col_{iClass},'LineWidth',1.5)
        %                 area(x,sorted_Tmean{iArea}(iClass,:),'Facecolor',col_{iClass},'FaceAlpha',0.5)
    end
    axis([0 20 -1 1])
end

%% Pool specific effects of errors - collapsed
x = (1:Ltr).*bw;

for iClass = 1:3
    figure;
    for iArea = 1:2
        subplot(1,2,iArea); hold on
        y =  FracSig{iArea}(iClass,:);
        plot(x,y,'color',col_{iClass},'LineWidth',3)% plot(x,smooth2a(y,1,5),'color',col_{iClass},'LineWidth',3)
        y =  FracSig_err{iArea}(iClass,:);
        plot(x,y,'color',col_{iClass},'LineWidth',1.5,'LineStyle','-') %plot(x,smooth2a(y,1,5),'color',col_{iClass},'LineWidth',1.5,'LineStyle','-')
        plot([10 10],[0 1],'color',[0 0 0 0.6],'LineWidth',1.5)
        plot([5 5],[0 1],'color',[0 1 0 0.6],'LineWidth',1.5)
        plot([15 15],[0 1],'color',[1 0 0 0.6],'LineWidth',1.5)
        axis([0 20 0 0.5])
    end
end
%% Pool specific effects of errors
x = (1:size(sorted_T,2)).*bw;
clear TSsig_ FracSig_Corr FracSig_Err FracSig_Corr_M FracSig_Corr_E FracSig_Err_M FracSig_Err_E

for iArea = 1:2
    for iClass=1:3
        for iFile = 1:length(fileList)
            TSsig_ = D.TSsig{iArea}{iFile}(:,D.Membership{iArea}{iFile}==iClass-1)';
            if ~isempty(TSsig_)
                FracSig_Corr{iArea}{iClass}(:,iFile) = nansum(TSsig_)./size(TSsig_,1);
            else
                FracSig_Corr{iArea}{iClass}(:,iFile) = nan(1,length(x));
            end
            
            TSsig_ = D.TSsig_err{iArea}{iFile}(:,D.Membership{iArea}{iFile}==iClass-1)';
            if ~isempty(TSsig_)
                FracSig_Err{iArea}{iClass}(:,iFile) = nansum(TSsig_)./size(TSsig_,1);
            else
                FracSig_Err{iArea}{iClass}(:,iFile) = nan(1,length(x));
            end
        end
        FracSig_Corr_M{iArea}(:,iClass) = nanmean(smooth2a(FracSig_Corr{iArea}{iClass},5,0),2);
        FracSig_Corr_E{iArea}(:,iClass) = nansem(smooth2a(FracSig_Corr{iArea}{iClass},5,0),2);
        
        FracSig_Err_M{iArea}(:,iClass)  = nanmean(smooth2a(FracSig_Err{iArea}{iClass},5,0),2);
        FracSig_Err_E{iArea}(:,iClass)  = nansem(smooth2a(FracSig_Err{iArea}{iClass},5,0),2);
    end
end

for iClass = 1:3
    figure;
    for iArea = 1:2
        subplot(1,2,iArea); hold on
        
        Y1 = smooth2a(FracSig_Corr{iArea}{iClass},0,0);Y1(:,isnan(nansum(Y1)))=[];
        Y2 = smooth2a(FracSig_Err{iArea}{iClass},0,0);Y2(:,isnan(nansum(Y2)))=[];
        [sig,~] = permtest2vec(Y1,Y2,100,0.05);
        
        a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
        plot(x_,a*1,'color','k','LineWidth',5)

        y  =  FracSig_Corr_M{iArea}(:,iClass);
        yE =  FracSig_Corr_E{iArea}(:,iClass);
%         plot(x,y,'color',col_{iClass},'LineWidth',3)% plot(x,smooth2a(y,1,5),'color',col_{iClass},'LineWidth',3)
        ciplot(y+yE,y-yE,x,col_{iClass},0.8)
        y  =  FracSig_Err_M{iArea}(:,iClass);
        yE =  FracSig_Err_E{iArea}(:,iClass);
%         plot(x,y,'color',col_{iClass},'LineWidth',1.5,'LineStyle',':') %plot(x,smooth2a(y,1,5),'color',col_{iClass},'LineWidth',1.5,'LineStyle','-')
        ciplot(y+yE,y-yE,x,col_{iClass},0.5)
        
        plot([10 10],[0 1],'color',[0 0 0 0.6],'LineWidth',1.5)
        plot([5 5],[0 1],'color',[0 1 0 0.6],'LineWidth',1.5)
        plot([15 15],[0 1],'color',[1 0 0 0.6],'LineWidth',1.5)
        axis([0 20 0 1.1])
    end
end
%% Pool specific effects of errors
scaling = 0.5;
figure;
for iArea = 1:2
    
    subplot(1,2,iArea); hold on
    plot([10 10],[0 10],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 10],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 10],'color',[1 0 0 0.3],'LineWidth',1.5)
    for iClass = 1:3
        y = FracSig{iArea}(iClass,:);
        y =  smooth2a(y,0,5)';
        patch([min(x),x, max(x)],[scaling*iClass;scaling*iClass+y;scaling*iClass],col_{iClass},'FaceAlpha',1);
        
        
        %         y =  FracSig{iArea}(iClass,:);
        %         plot(x,smooth2a(y,0,5),'color',col_{iClass},'LineWidth',3)
        %         y =  FracSig_err{iArea}(iClass,:);
        %         plot(x,smooth2a(y,1,5),'color',col_{iClass},'LineWidth',1.5,'LineStyle',':')
        
        
    end
    axis([0 20 0.4 2]); axis off
end
%% Plot example Ass vs member units - params
% iFile = 1;
% iArea = 1;
% iAss = 5;

iFile = 5;
iArea = 3;
iAss = 1;
offset = 10;
step = 2;
%% Plot example Ass vs member units - correct
clear temp
temp.Ass = D_Ass.TS{iArea}{iFile}(:,iAss);
if iArea==3
    temp.units = [D.TS{1}{iFile},D.TS{2}{iFile}];
else
    temp.units = D.TS{iArea}{iFile};
end
temp.members  = D_Ass.units{iArea}{iFile}{iAss};
temp.usel_out = D_Ass.usel_out{iArea}{iFile};

temp.nPFC = length(D_Ass.usel_out{1}{iFile});
x = (1:length(temp.Ass)).*bw;

figure; hold on
area(x,(temp.Ass),'FaceColor','g','EdgeColor','g')
clear y_
if iArea ==3
    y_{1} = temp.units(:,temp.members(temp.members<=temp.nPFC));
    y_{2} = temp.units(:,temp.members(temp.members>temp.nPFC));
    for i =1:size(y_{1},2)
        plot(x,(y_{1}(:,i))+i*step+offset,'b','LineWidth',1.5)
%         plot(x,(y_{1}(:,i)),'b','LineWidth',1.5)
    end
    for i =1:size(y_{2},2)
        plot(x,(y_{2}(:,i))+i*step+size(y_{1},2)+1.5*offset,'r','LineWidth',1.5)       
%         plot(x,(y_{2}(:,i)),'r','LineWidth',1.5)
    end
else
    y_{iArea} = temp.units(:,temp.members);
    for i =1:size(y_{iArea},2)
        plot(x,(y_{iArea}(:,i))+i*step+offset,color_{iArea})
    end
end
axis([0 20 0 40])
axis off
%% Plot example Ass vs member units - error

clear temp
temp.Ass = D_Ass.TS_err{iArea}{iFile}(:,iAss);
if iArea==3
    temp.units = [D.TS_err{1}{iFile},D.TS_err{2}{iFile}];
else
    temp.units = D.TS_err{iArea}{iFile};
end
temp.members  = D_Ass.units{iArea}{iFile}{iAss};
temp.usel_out = D_Ass.usel_out{iArea}{iFile};

temp.nPFC = length(D_Ass.usel_out{1}{iFile});
x = (1:length(temp.Ass)).*bw;

figure; hold on
area(x,(temp.Ass),'FaceColor','g','EdgeColor','g')
clear y_
if iArea ==3
    y_{1} = temp.units(:,temp.members(temp.members<=temp.nPFC));
    y_{2} = temp.units(:,temp.members(temp.members>temp.nPFC));
    for i =1:size(y_{1},2)
        plot(x,(y_{1}(:,i))+i*step+offset,'b','LineWidth',1.5)
    end
    for i = 1:size(y_{2},2)
        plot(x,(y_{2}(:,i))+i*step+size(y_{1},2)+1.5*offset,'r','LineWidth',1.5)
    end
else
    y_{iArea} = temp.units(:,temp.members);
    for i = 1:size(y_{iArea},2)
        plot(x,(y_{iArea}(:,i))+i*step+offset,color_{iArea},'LineWidth',1.5)
    end
end
axis([0 20 0 40])
axis off

%% Distribution of pairwise distances between D' scores of units/units and units/assemblies - collapsed
clear PDist_  temp
sigtimeThreshSec = 0.05;
compfunction = @(a,b) (b-a)./(a+b); % a=correct b=error
aggregatefun = @(x) (x);%mean
normalizefun = @(x) zscore(x);% zscore mat2gray2D
% normalizefun = @(x) x-nanmean(x);% zscore mat2gray2D
% normalizefun = @(x) x;% zscore mat2gray2D
% distanceMetric = 'zscore';
% bins = 0:0.005:1;
% bins2 = -5:.25:5;

distanceMetric = 'cosine';
bins = 0:0.1:2;
bins2 = -5:.2:5;
% distanceMetric = 'euclidean';
% bins = 0:0.05:50;
% bins2 = -50:5:50;
% % bins2 = -1:0.1:1;

for iArea = 1:3
    for iFile = 1:length(D_Ass.TS{iArea})
        
        nPFC = size(D.TS{1}{iFile},2);
        
        % Get data correct
        temp.TS_Ass   = D_Ass.TS{iArea}{iFile};
        temp.TS_AssSig   = D_Ass.TSsig {iArea}{iFile};
        if iArea < 3
            temp.TS_Units = D.TS{iArea}{iFile};
            temp.TS_UnitsSig = D.TSsig{iArea}{iFile};
        else
            temp.TS_Units = [D.TS{1}{iFile}, D.TS{2}{iFile}];
            temp.TS_UnitsSig = [D.TSsig{1}{iFile}, D.TSsig{2}{iFile}];
        end
        
        % Get data error
        temp.TS_AssErr    = D_Ass.TS_err{iArea}{iFile};
        temp.TS_AssSigErr = D_Ass.TSsig_err{iArea}{iFile};
        if iArea < 3
            temp.TS_UnitsErr = D.TS_err{iArea}{iFile};
            temp.TS_UnitsSigErr = D.TSsig_err{iArea}{iFile};
        else
            temp.TS_UnitsErr = [D.TS_err{1}{iFile},D.TS_err{2}{iFile}];
            temp.TS_UnitsSigErr = [D.TSsig_err{1}{iFile},D.TSsig_err{2}{iFile}];
        end
        
        % Preallocate correct
        if iArea < 3
            PDist_.A2Umembers{iArea}{iFile}    = []; % mean Assembly to unit distance (members)
            PDist_.A2Uothermembers{iArea}{iFile} = [];
            PDist_.A2Unonmembers{iArea}{iFile} = [];
            PDist_.U2U_M2M{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-member)
            PDist_.U2U_M2N{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-nonmember)
            PDist_.U2U_N2N{iArea}{iFile} = [];       % mean inter-unit coherence in information (nonmember-to-nonmember)
        else
            PDist_.A2Umembers{iArea}{1}{iFile}    = [];  PDist_.A2Umembers{iArea}{2}{iFile}    = [];
            PDist_.A2Uothermembers{iArea}{1}{iFile} = [];  PDist_.A2Uothermembers{iArea}{2}{iFile} = [];
            PDist_.A2Unonmembers{iArea}{1}{iFile} = [];  PDist_.A2Unonmembers{iArea}{2}{iFile} = [];
            PDist_.U2U_M2M{iArea}{1}{iFile} = []; %PFC only portion
            PDist_.U2U_M2M{iArea}{2}{iFile} = []; %HP only portion
            PDist_.U2U_M2M{iArea}{3}{iFile} = []; %HP-PFC only portion
            PDist_.U2U_N2N{iArea}{1}{iFile} = [];
            PDist_.U2U_M2N{iArea}{1}{iFile} = [];    PDist_.U2U_M2N{iArea}{2}{iFile} = [];    %member to nonmember {Mpfc2Nhp}{Mhp2Npfc}
        end
        
        % Preallocate error
        if iArea < 3
            PDist_.A2UmembersError{iArea}{iFile}    = []; % mean Assembly to unit distance (members)
            PDist_.A2UothermembersError{iArea}{iFile} = [];
            PDist_.A2UnonmembersError{iArea}{iFile} = [];
            PDist_.U2U_M2MError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-member)
            PDist_.U2U_M2NError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-nonmember)
            PDist_.U2U_N2NError{iArea}{iFile} = [];       % mean inter-unit coherence in information (nonmember-to-nonmember)
        else
            PDist_.A2UmembersError{iArea}{1}{iFile}    = [];            PDist_.A2UmembersError{iArea}{2}{iFile}    = [];
            PDist_.A2UothermembersError{iArea}{1}{iFile} = [];          PDist_.A2UothermembersError{iArea}{2}{iFile} = [];
            PDist_.A2UnonmembersError{iArea}{1}{iFile} = [];            PDist_.A2UnonmembersError{iArea}{2}{iFile} = [];
            PDist_.U2U_M2MError{iArea}{1}{iFile} = [];
            PDist_.U2U_M2MError{iArea}{2}{iFile} = [];
            PDist_.U2U_M2MError{iArea}{3}{iFile} = [];
            PDist_.U2U_N2NError{iArea}{1}{iFile} = [];
            PDist_.U2U_M2NError{iArea}{1}{iFile} = [];    PDist_.U2U_M2NError{iArea}{2}{iFile} = [];    %member to nonmember {Mpfc2Nhp}{Mhp2Npfc}
        end
        
        % Process correct and error
        for iAss = 1:size(temp.TS_Ass,2)
            if sum(temp.TS_AssSig(:,iAss),1)*bw>=sigtimeThreshSec %& sum(temp.TS_AssSigErr(:,iAss),1)*bw>=sigtimeThreshSec
                
                temp.members      = D_Ass.units{iArea}{iFile}{iAss};
                temp.othermembers = setdiff(unique(cell2mat(D_Ass.units{iArea}{iFile})),temp.members);
                temp.allmembers   = unique(cell2mat(D_Ass.units{iArea}{iFile}));
                temp.nonmembers   = setdiff(1:size(temp.TS_Units,2),temp.allmembers);
                
                input     = normalizefun([temp.TS_Ass(:,iAss),temp.TS_Units]);
                inputsig  = [sum(temp.TS_AssSig(:,iAss),1),sum(temp.TS_UnitsSig,1)];
                
                inputErr     = normalizefun([temp.TS_AssErr(:,iAss),temp.TS_UnitsErr]);
                inputsigErr  = [sum(temp.TS_AssSigErr(:,iAss),1),sum(temp.TS_UnitsSigErr,1)];
                
                %idx = find(inputsig(2:end)*bw<=sigtimeThreshSec & inputsigErr(2:end)*bw<=sigtimeThreshSec);
                idx = find(inputsig(2:end)*bw<=sigtimeThreshSec);

                temp.members      = setdiff(temp.members,idx);
                temp.othermembers = setdiff(temp.othermembers,idx);
                temp.allmembers   = setdiff(temp.allmembers,idx);
                temp.nonmembers   = setdiff(temp.nonmembers,idx);
                
                % Process correct
%                 input = input((size(input,1))/2+1:end,:);
%                 input = input(1:(size(input,1))/2,:);
                coeff = triu(squareform(pdist((input)',distanceMetric)));
                
                A2U = coeff(1,2:end);
                U2U = coeff(2:end,2:end);
                if iArea < 3
                    % Assembly to unit distances
                    PDist_.A2Umembers{iArea}{iFile}      = [PDist_.A2Umembers{iArea}{iFile},     aggregatefun(A2U(temp.members))];
                    PDist_.A2Uothermembers{iArea}{iFile} = [PDist_.A2Uothermembers{iArea}{iFile},  aggregatefun(A2U(temp.othermembers))];
                    PDist_.A2Unonmembers{iArea}{iFile}   = [PDist_.A2Unonmembers{iArea}{iFile},  aggregatefun(A2U(temp.nonmembers))];
                    
                    % Unit to unit distances
                    
                    D_ = []; % member to member
                    for i = temp.members
                        for j = temp.members
                            if i<j
                                D_=[D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_M2M{iArea}{iFile} = [PDist_.U2U_M2M{iArea}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % member to nonmember
                    for i = temp.members
                        for j = temp.nonmembers
                            if i<j
                                D_=[D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_M2N{iArea}{iFile} = [PDist_.U2U_M2N{iArea}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % nonmember to nonmember
                    for i =  temp.nonmembers
                        for j = temp.nonmembers
                            if i<j
                                D_=[D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_N2N{iArea}{iFile} = [PDist_.U2U_N2N{iArea}{iFile}, aggregatefun(D_(:)')];
                    
                    
                else
                    % Assembly to unit distances
                    PDist_.A2Umembers{iArea}{1}{iFile} = [PDist_.A2Umembers{iArea}{1}{iFile},     aggregatefun(A2U(temp.members(temp.members<=nPFC)))];
                    PDist_.A2Umembers{iArea}{2}{iFile} = [PDist_.A2Umembers{iArea}{2}{iFile},     aggregatefun(A2U(temp.members(temp.members>nPFC)))];
                    
                    PDist_.A2Uothermembers{iArea}{1}{iFile} = [PDist_.A2Uothermembers{iArea}{1}{iFile},aggregatefun(A2U(temp.othermembers(temp.othermembers<=nPFC)))];
                    PDist_.A2Uothermembers{iArea}{2}{iFile} = [PDist_.A2Uothermembers{iArea}{2}{iFile},aggregatefun(A2U(temp.othermembers(temp.othermembers>nPFC)))];
                    
                    PDist_.A2Unonmembers{iArea}{1}{iFile} = [PDist_.A2Unonmembers{iArea}{1}{iFile},aggregatefun(A2U(temp.nonmembers(temp.nonmembers<=nPFC)))];
                    PDist_.A2Unonmembers{iArea}{2}{iFile} = [PDist_.A2Unonmembers{iArea}{2}{iFile},aggregatefun(A2U(temp.nonmembers(temp.nonmembers>nPFC)))];
                    
                    % Unit to unit distances
                    D_ = []; % member to member (PFC to PFC exclusively)
                    for i = temp.members(temp.members<=nPFC)
                        for j = temp.members(temp.members<=nPFC)
                            if i<j
                                D_ = [D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_M2M{iArea}{1}{iFile} = [PDist_.U2U_M2M{iArea}{1}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % member to member (HP to HP exclusively)
                    for i = temp.members(temp.members>nPFC)
                        for j = temp.members(temp.members>nPFC)
                            if i<j
                                D_ = [D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_M2M{iArea}{2}{iFile} = [PDist_.U2U_M2M{iArea}{2}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % member to member (HP to PFC exclusively)
                    for i = temp.members(temp.members<=nPFC)
                        for j = temp.members(temp.members>nPFC)
                            D_ = [D_;U2U(i,j)];
                        end
                    end
                    PDist_.U2U_M2M{iArea}{3}{iFile} = [PDist_.U2U_M2M{iArea}{3}{iFile}, aggregatefun(D_(:)')];
                    
                    
                    D_ = []; % PFC member to HP nonmember
                    for i = temp.members(temp.members<=nPFC)
                        for j = temp.nonmembers(temp.nonmembers>nPFC)
                            D_ = [D_;U2U(i,j)];
                        end
                    end
                    PDist_.U2U_M2N{iArea}{1}{iFile} = [PDist_.U2U_M2N{iArea}{1}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % PFC nonmember to HP member
                    for i = temp.nonmembers(temp.nonmembers<=nPFC)
                        for j = temp.members(temp.members>nPFC)
                            D_ = [D_;U2U(i,j)];
                        end
                    end
                    PDist_.U2U_M2N{iArea}{2}{iFile} = [PDist_.U2U_M2N{iArea}{2}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % PFC nonmember to HP nonmember
                    for i = temp.nonmembers(temp.nonmembers<=nPFC)
                        for j = temp.nonmembers(temp.nonmembers>nPFC)
                            D_ = [D_;U2U(i,j)];
                        end
                    end
                    PDist_.U2U_N2N{iArea}{1}{iFile} = [PDist_.U2U_N2N{iArea}{1}{iFile}, aggregatefun(D_(:)')];
                    
                end
                
                % Process error
%                 input = input((size(input,1))/2+1:end,:);
%                 input = input(1:(size(input,1))/2,:);
                coeff = triu(squareform(pdist((inputErr)',distanceMetric)));
                
                A2U = coeff(1,2:end);
                U2U = coeff(2:end,2:end);
                if iArea < 3
                    PDist_.A2UmembersError{iArea}{iFile}    = [PDist_.A2UmembersError{iArea}{iFile},     aggregatefun(A2U(temp.members))];
                    PDist_.A2UothermembersError{iArea}{iFile} = [PDist_.A2UothermembersError{iArea}{iFile},  aggregatefun(A2U(temp.othermembers))];
                    PDist_.A2UnonmembersError{iArea}{iFile} = [PDist_.A2UnonmembersError{iArea}{iFile},  aggregatefun(A2U(temp.nonmembers))];
                    
                    % Unit to unit distances
                    
                    D_ = []; % member to member
                    for i = temp.members
                        for j = temp.members
                            if i<j
                                D_=[D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_M2MError{iArea}{iFile} = [PDist_.U2U_M2MError{iArea}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % member to nonmember
                    for i = temp.members
                        for j = temp.nonmembers
                            if i<j
                                D_=[D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_M2NError{iArea}{iFile} = [PDist_.U2U_M2NError{iArea}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % nonmember to nonmember
                    for i =  temp.nonmembers
                        for j = temp.nonmembers
                            if i<j
                                D_=[D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_N2NError{iArea}{iFile} = [PDist_.U2U_N2NError{iArea}{iFile}, aggregatefun(D_(:)')];
                    
                    
                else
                    PDist_.A2UmembersError{iArea}{1}{iFile} = [PDist_.A2UmembersError{iArea}{1}{iFile},     aggregatefun(A2U(temp.members(temp.members<=nPFC)))];
                    PDist_.A2UmembersError{iArea}{2}{iFile} = [PDist_.A2UmembersError{iArea}{2}{iFile},     aggregatefun(A2U(temp.members(temp.members>nPFC)))];
                    
                    PDist_.A2UothermembersError{iArea}{1}{iFile} = [PDist_.A2UothermembersError{iArea}{1}{iFile},aggregatefun(A2U(temp.othermembers(temp.othermembers<=nPFC)))];
                    PDist_.A2UothermembersError{iArea}{2}{iFile} = [PDist_.A2UothermembersError{iArea}{2}{iFile},aggregatefun(A2U(temp.othermembers(temp.othermembers>nPFC)))];
                    
                    PDist_.A2UnonmembersError{iArea}{1}{iFile} = [PDist_.A2UnonmembersError{iArea}{1}{iFile},aggregatefun(A2U(temp.nonmembers(temp.nonmembers<=nPFC)))];
                    PDist_.A2UnonmembersError{iArea}{2}{iFile} = [PDist_.A2UnonmembersError{iArea}{2}{iFile},aggregatefun(A2U(temp.nonmembers(temp.nonmembers>nPFC)))];
                    
                    
                    % Unit to unit distances
                    D_ = []; % member to member (PFC to PFC exclusively)
                    for i = temp.members(temp.members<=nPFC)
                        for j = temp.members(temp.members<=nPFC)
                            if i<j
                                D_ = [D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_M2MError{iArea}{1}{iFile} = [PDist_.U2U_M2MError{iArea}{1}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % member to member (HP to HP exclusively)
                    for i = temp.members(temp.members>nPFC)
                        for j = temp.members(temp.members>nPFC)
                            if i<j
                                D_ = [D_;U2U(i,j)];
                            end
                        end
                    end
                    PDist_.U2U_M2MError{iArea}{2}{iFile} = [PDist_.U2U_M2MError{iArea}{2}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % member to member (HP to PFC exclusively)
                    for i = temp.members(temp.members<=nPFC)
                        for j = temp.members(temp.members>nPFC)
                            D_ = [D_;U2U(i,j)];
                        end
                    end
                    PDist_.U2U_M2MError{iArea}{3}{iFile} = [PDist_.U2U_M2MError{iArea}{3}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % PFC member to HP nonmember
                    for i = temp.members(temp.members<=nPFC)
                        for j = temp.nonmembers(temp.nonmembers>nPFC)
                            D_ = [D_;U2U(i,j)];
                        end
                    end
                    PDist_.U2U_M2NError{iArea}{1}{iFile} = [PDist_.U2U_M2NError{iArea}{1}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % PFC nonmember to HP member
                    for i = temp.nonmembers(temp.nonmembers<=nPFC)
                        for j = temp.members(temp.members>nPFC)
                            D_ = [D_;U2U(i,j)];
                        end
                    end
                    PDist_.U2U_M2NError{iArea}{2}{iFile} = [PDist_.U2U_M2NError{iArea}{2}{iFile}, aggregatefun(D_(:)')];
                    
                    D_ = []; % PFC nonmember to HP nonmember
                    for i =  temp.nonmembers(temp.nonmembers<=nPFC)
                        for j =  temp.nonmembers(temp.nonmembers>nPFC)
                            D_ = [D_;U2U(i,j)];
                        end
                    end
                    PDist_.U2U_N2NError{iArea}{1}{iFile} = [PDist_.U2U_N2NError{iArea}{1}{iFile}, aggregatefun(D_(:)')];
                    
                end
            end
        end
        
        % Compare correct and error trial results
        if iArea < 3
            a = PDist_.A2Umembers{iArea}{iFile};
            b = PDist_.A2UmembersError{iArea}{iFile};
            PDist_.A2UmembersDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDist_.A2Uothermembers{iArea}{iFile};
            b = PDist_.A2UothermembersError{iArea}{iFile};
            PDist_.A2UothermembersDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDist_.A2Unonmembers{iArea}{iFile};
            b = PDist_.A2UnonmembersError{iArea}{iFile};
            PDist_.A2UnonmembersDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDist_.U2U_M2M{iArea}{iFile};
            b = PDist_.U2U_M2MError{iArea}{iFile};
            PDist_.U2U_M2MDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDist_.U2U_M2N{iArea}{iFile};
            b = PDist_.U2U_M2NError{iArea}{iFile};
            PDist_.U2U_M2NDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDist_.U2U_N2N{iArea}{iFile};
            b = PDist_.U2U_N2NError{iArea}{iFile};
            PDist_.U2U_N2NDelta{iArea}{iFile} =  compfunction(a,b);
            
            
        else
            for i=1:3
                a = PDist_.U2U_M2M{iArea}{i}{iFile};
                b = PDist_.U2U_M2MError{iArea}{i}{iFile};
                PDist_.U2U_M2MDelta{iArea}{i}{iFile} =  compfunction(a,b);
            end
            a = PDist_.U2U_N2N{iArea}{1}{iFile};
            b = PDist_.U2U_N2NError{iArea}{1}{iFile};
            PDist_.U2U_N2NDelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            
            for i=1:2
                a = PDist_.A2Umembers{iArea}{i}{iFile};
                b = PDist_.A2UmembersError{iArea}{i}{iFile};
                PDist_.A2UmembersDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
                a = PDist_.A2Uothermembers{iArea}{i}{iFile};
                b = PDist_.A2UothermembersError{iArea}{i}{iFile};
                PDist_.A2UothermembersDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
                a = PDist_.A2Unonmembers{iArea}{i}{iFile};
                b = PDist_.A2UnonmembersError{iArea}{i}{iFile};
                PDist_.A2UnonmembersDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
                a = PDist_.U2U_M2N{iArea}{i}{iFile};
                b = PDist_.U2U_M2NError{iArea}{i}{iFile};
                PDist_.U2U_M2NDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
            end
        end
    end
end   
%% Plot distribution of unit-assembly distances - correct
figure('name','Unit to Assembly distance - correct');
for iArea = 1:3
    
    
    if iArea<3
        subplot(2,2,iArea); hold on
        title(Areas{iArea})
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2Umembers{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,[0 0 0],1)
        
%         hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UmembersError{iArea},'UniformOutput',false);
%         hist_ = cell2mat(hist_(~isempty_cell(hist_))');
%         ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,[1 0 0],1)
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2Uothermembers{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.3*[1 1 1],1)
        
%         hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UothermembersError{iArea},'UniformOutput',false);
%         hist_ = cell2mat(hist_(~isempty_cell(hist_))');
%         ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.3*[1 0 0],1)

        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2Unonmembers{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.6*[1 1 1],1)
        
%         hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UnonmembersError{iArea},'UniformOutput',false);
%         hist_ = cell2mat(hist_(~isempty_cell(hist_))');
%         ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.6*[1 0 0],1)
         
%         legend('Member units to parent assembly ','Member units to non-parent assembly ','Nonmember units to parent assembly','Location','SouthEast'); legend boxoff

    else
        for iArea_ =1:2
            subplot(2,2,2+iArea_); hold on
            hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2Umembers{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,color_{iArea_},1)
            
            hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2Uothermembers{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,color_{iArea_},0.6)
            
            hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2Unonmembers{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,color_{iArea_},0.3)
            title([Areas{iArea} ' (' Areas{iArea_} ')'])

                ylabel('Fraction of units')

        end
    end
%     xlabel('Distance between Assembly and Unit decoding scores')
    ylabel('Fraction of units')
    axis([min(bins) max(bins) 0 1])

    %     hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UmembersError{iArea},'UniformOutput',false);
    %     hist_ = cell2mat(hist_(~isempty_cell(hist_))');
    %     ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'k',0.5)
    %
    %     hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UnonmembersError{iArea},'UniformOutput',false);
    %     hist_ = cell2mat(hist_(~isempty_cell(hist_))');
    %     ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'r',0.5)
end
%% Plot distribution of unit-assembly distances - errors
figure('name','Unit to Assembly distance - error');
for iArea = 1:3
    
    
    if iArea<3
        subplot(2,2,iArea); hold on
        title(Areas{iArea})
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UmembersError{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,[0 0 0],1)
        
%         hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UmembersError{iArea},'UniformOutput',false);
%         hist_ = cell2mat(hist_(~isempty_cell(hist_))');
%         ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,[1 0 0],1)
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UothermembersError{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.3*[1 1 1],1)
        
%         hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UothermembersError{iArea},'UniformOutput',false);
%         hist_ = cell2mat(hist_(~isempty_cell(hist_))');
%         ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.3*[1 0 0],1)

        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UnonmembersError{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.6*[1 1 1],1)
        
%         hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UnonmembersError{iArea},'UniformOutput',false);
%         hist_ = cell2mat(hist_(~isempty_cell(hist_))');
%         ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.6*[1 0 0],1)
         
%         legend('Member units to parent assembly ','Member units to non-parent assembly ','Nonmember units to parent assembly','Location','SouthEast'); legend boxoff

    else
        for iArea_ =1:2
            subplot(2,2,2+iArea_); hold on
            hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UmembersError{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,color_{iArea_},1)
            
            hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UothermembersError{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,color_{iArea_},0.6)
            
            hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UnonmembersError{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,color_{iArea_},0.3)
            title([Areas{iArea} ' (' Areas{iArea_} ')'])

                ylabel('Fraction of units')

        end
    end
%     xlabel('Distance between Assembly and Unit decoding scores')
    ylabel('Fraction of units')
    axis([min(bins) max(bins) 0 1])

    %     hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UmembersError{iArea},'UniformOutput',false);
    %     hist_ = cell2mat(hist_(~isempty_cell(hist_))');
    %     ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'k',0.5)
    %
    %     hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.A2UnonmembersError{iArea},'UniformOutput',false);
    %     hist_ = cell2mat(hist_(~isempty_cell(hist_))');
    %     ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'r',0.5)
end
%% Plot change in unit-assembly distance on errors
figure
for iArea = 1:3
    
    subplot(3,1,iArea); hold on
    if iArea<3
        hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.A2UmembersDelta{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,'k',1)
        
        hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.A2UothermembersDelta{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,'k',0.6)
        
        hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.A2UnonmembersDelta{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,'k',0.3)
    else
        for iArea_ =1:2
            hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.A2UmembersDelta{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,color_{iArea_},1)
            
            hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.A2UnonmembersDelta{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,color_{iArea_},0.6)
        end
    end
    plot([0 0],[0 1],':k','LineWidth',1.5)
    title(Areas{iArea})
    axis([min(bins2) max(bins2) 0 0.6])
end
%% plot average pairwise distances between units - correct
figure; %correct
for iArea = 1:3
%     figure; hold on
    subplot(3,1,iArea); hold on
    title(['Unit-Unit distance: ', Areas{iArea}])
    if iArea<3
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2M{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0*[1 1 1],1)
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2N{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.4*[1 1 1],1)
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_N2N{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.6*[1 1 1],1)
        legend('Between members','Member to Non-member', 'Between non-members','Location','SouthEast'); legend boxoff
        
    else
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2M{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0*[1 1 1],1)
        
        for iArea_ =1:2
            hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2N{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,color_{iArea_},0.6)
        end
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_N2N{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.6*[1 1 1],1)
        legend('Between members','PFC Member to HP Non-member','HP Member to PFC Non-member', 'Between non-members','Location','SouthEast'); legend boxoff
        
    end
    axis([min(bins) max(bins) 0 1])
end
%% plot average pairwise distances between units - errors
figure; %error
for iArea = 1:3
    
    subplot(3,1,iArea); hold on
    title(['Unit-Unit distance: ', Areas{iArea}])
    if iArea<3
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2MError{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0*[1 1 1],1)
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2NError{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.4*[1 1 1],1)
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_N2NError{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.6*[1 1 1],1)
        legend('Between members','Member to Non-member', 'Between non-members','Location','SouthEast'); legend boxoff
        
    else
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2MError{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0*[1 1 1],1)
        
        for iArea_ =1:2
            hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2NError{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,color_{iArea_},0.6)
        end
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_N2NError{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,0.6*[1 1 1],1)
        legend('Between members','PFC Member to HP Non-member','HP Member to PFC Non-member', 'Between non-members','Location','SouthEast'); legend boxoff
        
    end
    axis([min(bins) max(bins) 0 1])
end
%% Plot histograms of unit-unit distances
figure('name','Unit-Unit distances- Member/Member');
for iArea = 1:3
    subplot(1,3,iArea); hold on
    title(Areas{iArea})
    if iArea<3
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2MError{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'r',0.8)
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2M{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'k',0.8)
        plot(bins,nanmean(hist_),'Color','k','LineWidth',1.5)
%         x = cell2mat(cellfun(@(x) cell2mat(x')),PDist_.U2U_M2M{iArea}{1})
    else
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2MError{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'r',0.8)

        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2M{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'k',0.8)
%         plot(bins,nanmean(hist_),'Color',color_{iArea},'LineWidth',1.5)
        
    end
    axis square

end

figure('name','Unit-Unit distances- Member/NonMember');
for iArea = 1:3
    subplot(1,3,iArea); hold on
    title(Areas{iArea})
    if iArea<3
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2NError{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'r',0.8)
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2N{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'k',0.8)
        plot(bins,nanmean(hist_),'Color','k','LineWidth',1.5)
%         x = cell2mat(cellfun(@(x) cell2mat(x')),PDist_.U2U_M2M{iArea}{1})
    else
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2NError{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'r',0.8)

        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_M2N{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'k',0.8)
%         plot(bins,nanmean(hist_),'Color',color_{iArea},'LineWidth',1.5)
        
    end
    axis square

end

figure('name','Unit-Unit distances- NonMember/NonMember');
for iArea = 1:3
    subplot(1,3,iArea); hold on
    title(Areas{iArea})
    if iArea<3
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_N2NError{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'r',0.8)
        
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_N2N{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'k',0.8)
        plot(bins,nanmean(hist_),'Color','k','LineWidth',1.5)
%         x = cell2mat(cellfun(@(x) cell2mat(x')),PDist_.U2U_M2M{iArea}{1})
    else
        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_N2NError{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'r',0.8)

        hist_ = cellfun(@(x) cumsum(histc(x,bins))./numel(x),PDist_.U2U_N2N{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins,'k',0.8)
%         plot(bins,nanmean(hist_),'Color',color_{iArea},'LineWidth',1.5)
        
    end
    axis square
end
%% Plot change in unit/unit distance on errors
figure
for iArea = 1:3
    subplot(3,1,iArea); hold on
    if iArea<3
        hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.U2U_M2MDelta{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,0*[1 1 1],1)
        
        hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.U2U_M2NDelta{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,0.4*[1 1 1],1)
        
        hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.U2U_N2NDelta{iArea},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,0.6*[1 1 1],1)
    else
        
        hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.U2U_M2MDelta{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,0*[1 1 1],1)
        
        hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.U2U_M2NDelta{iArea}{1},'UniformOutput',false);
        hist_ = cell2mat(hist_(~isempty_cell(hist_))');
        ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,0.6*[1 1 1],1)
        for iArea_ =1:2
            hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.A2UnonmembersDelta{iArea}{iArea_},'UniformOutput',false);
            hist_ = cell2mat(hist_(~isempty_cell(hist_))');
            ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,color_{iArea_},0.6)
        end
    end
    plot([0 0],[0 1],':k','LineWidth',1.5)
    title(Areas{iArea})
    axis([min(bins2) max(bins2) 0 0.6])

end
%% Plot change in unit/unit distance on errors - within Joint assems only
iArea = 3;
figure; hold on

for iArea_=1:3
    hist_ = cellfun(@(x) (histc(x,bins2))./numel(x),PDist_.U2U_M2MDelta{iArea}{iArea_},'UniformOutput',false);
    hist_ = cell2mat(hist_(~isempty_cell(hist_))');
    ciplot(nanmean(hist_)+nansem(hist_),nanmean(hist_)-nansem(hist_),bins2,color_{iArea_},1)
  
end
    plot([0 0],[0 1],':k','LineWidth',1.5)

title(Areas{iArea})
axis([min(bins2) max(bins2) 0 0.6])
%% Distribution of pairwise distances between D' scores of units/units and units/assemblies - continuous 
clear PDistCont_ temp
compfunction = @(a,b) (b-a)./(a+b); % a=correct b=error
aggregatefun = @(x) (x);

% normfun = @(x) (x);
normfun = @(x) zscore(x);
% normfun = @(x) x - repmat(nanmean(x),size(x,1),1);
% aggregatefun = @(x) mean(x,1);

% distanceMetric = 'Hamming';
% bins = 0:0.005:1;
% bins2 = -5:.25:5;

% distanceMetric = 'cosine';
% bins = 0:0.1:1.5;
% bins2 = -5:.25:5;
distanceMetric = 'euclidean';
bins = 0:0.05:50;
bins2 = -50:2:50;
% bins2 = -1:0.1:1;

for iArea = 1:3
    for iFile = 1:length(D_Ass.TS{iArea})
        [iFile,iArea]
        nPFC = size(D.TS{1}{iFile},2);
        
        % Get data correct
        temp.TS_Ass   = D_Ass.TS{iArea}{iFile};
        if iArea < 3
            temp.TS_Units = D.TS{iArea}{iFile};
        else
            temp.TS_Units = [D.TS{1}{iFile}, D.TS{2}{iFile}];
        end
        % Preallocate correct
        if iArea < 3
            PDistCont_.A2Umembers{iArea}{iFile}    = []; % mean Assembly to unit distance (members)
            PDistCont_.A2Unonmembers{iArea}{iFile} = [];
            PDistCont_.U2U_M2M{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-member)
            PDistCont_.U2U_M2N{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-nonmember)
            PDistCont_.U2U_N2N{iArea}{iFile} = [];       % mean inter-unit coherence in information (nonmember-to-nonmember)
            PDistCont_.U2U_M2O{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-other ass member)
        else
            PDistCont_.A2Umembers{iArea}{1}{iFile}    = [];  PDistCont_.A2Umembers{iArea}{2}{iFile}    = [];
            PDistCont_.A2Unonmembers{iArea}{1}{iFile} = [];  PDistCont_.A2Unonmembers{iArea}{2}{iFile} = [];
            PDistCont_.U2U_M2M{iArea}{1}{iFile} = [];
            PDistCont_.U2U_N2N{iArea}{1}{iFile} = [];
            PDistCont_.U2U_M2N{iArea}{1}{iFile} = [];    PDistCont_.U2U_M2N{iArea}{2}{iFile} = [];    %member to nonmember {Mpfc2Nhp}{Mhp2Npfc}
            PDistCont_.U2U_M2O{iArea}{1}{iFile} = [];   % mean inter-unit coherence in information (member-to-other ass member)
        end
        % Process correct
        for iAss = 1:size(temp.TS_Ass,2)
            clear coeff
            temp.members = D_Ass.units{iArea}{iFile}{iAss};
            temp.Othermembers = setdiff(cell2mat(D_Ass.units{iArea}{iFile}),temp.members);
            
            input = ([temp.TS_Ass(:,iAss),temp.TS_Units]);
            input = normfun(input);
            
            for t =1:size(input,1)
                coeff(:,:,t) = triu(squareform(pdist(input(t,:)',distanceMetric)));
            end
            
%             squareform(pdist(input(:,[1:2])',distanceMetric))
%             x = input(:,1);
%             y = input(:,2);
%             xy   = dot(x,y);
%             nx   = norm(x);
%             ny   = norm(y);
%             nxny = nx*ny;
%             Cs   = xy/nxny;
            
            A2U = squeeze(coeff(1,2:end,:));
            U2U = coeff(2:end,2:end,:);
            if iArea < 3
                % Assembly to unit distances
                PDistCont_.A2Umembers{iArea}{iFile}    = [PDistCont_.A2Umembers{iArea}{iFile};     aggregatefun(A2U(temp.members,:))];
                PDistCont_.A2Unonmembers{iArea}{iFile} = [PDistCont_.A2Unonmembers{iArea}{iFile};  aggregatefun(A2U(setdiff(1:size(A2U,1),temp.members),:))];
                
                % Unit to unit distances
                
                D_ = []; % member to member
                for i = temp.members
                    for j = temp.members
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2M{iArea}{iFile} = [PDistCont_.U2U_M2M{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to nonmember
                for i = temp.members
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2N{iArea}{iFile} = [PDistCont_.U2U_M2N{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % nonmember to nonmember
                for i =  setdiff(1:size(A2U,1),temp.members)
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_N2N{iArea}{iFile} = [PDistCont_.U2U_N2N{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2O{iArea}{iFile} = [PDistCont_.U2U_M2O{iArea}{iFile}; aggregatefun(D_)];
                
            else
                PDistCont_.A2Umembers{iArea}{1}{iFile} = [PDistCont_.A2Umembers{iArea}{1}{iFile};     aggregatefun(A2U(temp.members(temp.members<=nPFC),:))];
                PDistCont_.A2Umembers{iArea}{2}{iFile} = [PDistCont_.A2Umembers{iArea}{2}{iFile};     aggregatefun(A2U(temp.members(temp.members>nPFC),:))];
                
                PDistCont_.A2Unonmembers{iArea}{1}{iFile} = [PDistCont_.A2Unonmembers{iArea}{1}{iFile};...
                    aggregatefun(A2U(setdiff(1:nPFC,temp.members(temp.members<=nPFC)),:))];
                
                PDistCont_.A2Unonmembers{iArea}{2}{iFile} = [PDistCont_.A2Unonmembers{iArea}{2}{iFile};...
                    aggregatefun(A2U(setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC)),:))];
                
                D_ = []; % member to member (HP to PFC exclusively)
                for i = temp.members(temp.members<=nPFC)
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2M{iArea}{1}{iFile} = [PDistCont_.U2U_M2M{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC member to HP nonmember
                for i = temp.members(temp.members<=nPFC)
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2N{iArea}{1}{iFile} = [PDistCont_.U2U_M2N{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP member
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2N{iArea}{2}{iFile} = [PDistCont_.U2U_M2N{iArea}{2}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP nonmember
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_N2N{iArea}{1}{iFile} = [PDistCont_.U2U_N2N{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2O{iArea}{1}{iFile} = [PDistCont_.U2U_M2O{iArea}{1}{iFile}; aggregatefun(D_)];
                
            end
        end
        
        % Get data error
        temp.TS_Ass   = D_Ass.TS_err{iArea}{iFile};
        if iArea < 3
            temp.TS_Units = D.TS_err{iArea}{iFile};
        else
            temp.TS_Units = [D.TS_err{1}{iFile},D.TS_err{2}{iFile}];
        end
        % Preallocate error
        if iArea < 3
            PDistCont_.A2UmembersError{iArea}{iFile}    = []; % mean Assembly to unit distance (members)
            PDistCont_.A2UnonmembersError{iArea}{iFile} = [];
            PDistCont_.U2U_M2MError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-member)
            PDistCont_.U2U_M2NError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-nonmember)
            PDistCont_.U2U_N2NError{iArea}{iFile} = [];       % mean inter-unit coherence in information (nonmember-to-nonmember)
            PDistCont_.U2U_M2OError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-other ass member)
        else
            PDistCont_.A2UmembersError{iArea}{1}{iFile}    = [];            PDistCont_.A2UmembersError{iArea}{2}{iFile}    = [];
            PDistCont_.A2UnonmembersError{iArea}{1}{iFile} = [];            PDistCont_.A2UnonmembersError{iArea}{2}{iFile} = [];
            PDistCont_.U2U_M2MError{iArea}{1}{iFile} = [];
            PDistCont_.U2U_N2NError{iArea}{1}{iFile} = [];
            PDistCont_.U2U_M2NError{iArea}{1}{iFile} = [];    PDistCont_.U2U_M2NError{iArea}{2}{iFile} = [];    %member to nonmember {Mpfc2Nhp}{Mhp2Npfc}
            PDistCont_.U2U_M2OError{iArea}{1}{iFile} = [];   % mean inter-unit coherence in information (member-to-other ass member)

        end
        % Process error
        for iAss = 1:size(temp.TS_Ass,2)
            temp.members = D_Ass.units{iArea}{iFile}{iAss};
            temp.Othermembers = setdiff(cell2mat(D_Ass.units{iArea}{iFile}),temp.members);
            input = ([temp.TS_Ass(:,iAss),temp.TS_Units]);
            input = normfun(input);
            for t =1:size(input,1)
                coeff(:,:,t) = triu(squareform(pdist(input(t,:)',distanceMetric)));
            end
            
            A2U = squeeze(coeff(1,2:end,:));
            U2U = coeff(2:end,2:end,:);
            if iArea < 3
                PDistCont_.A2UmembersError{iArea}{iFile}    = [PDistCont_.A2UmembersError{iArea}{iFile};     aggregatefun(A2U(temp.members,:))];
                PDistCont_.A2UnonmembersError{iArea}{iFile} = [PDistCont_.A2UnonmembersError{iArea}{iFile};  aggregatefun(A2U(setdiff(1:size(A2U,1),temp.members),:))];
                
                % Unit to unit distances
                
                D_ = []; % member to member
                for i = temp.members
                    for j = temp.members
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2MError{iArea}{iFile} = [PDistCont_.U2U_M2MError{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to nonmember
                for i = temp.members
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2NError{iArea}{iFile} = [PDistCont_.U2U_M2NError{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % nonmember to nonmember
                for i =  setdiff(1:size(A2U,1),temp.members)
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_N2NError{iArea}{iFile} = [PDistCont_.U2U_N2NError{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2OError{iArea}{iFile} = [PDistCont_.U2U_M2OError{iArea}{iFile}; aggregatefun(D_)];
                
                
            else
                PDistCont_.A2UmembersError{iArea}{1}{iFile} = [PDistCont_.A2UmembersError{iArea}{1}{iFile};     aggregatefun(A2U(temp.members(temp.members<=nPFC),:))];
                PDistCont_.A2UmembersError{iArea}{2}{iFile} = [PDistCont_.A2UmembersError{iArea}{2}{iFile};     aggregatefun(A2U(temp.members(temp.members>nPFC),:))];
                
                PDistCont_.A2UnonmembersError{iArea}{1}{iFile} = [PDistCont_.A2UnonmembersError{iArea}{1}{iFile};...
                    aggregatefun(A2U(setdiff(1:nPFC,temp.members(temp.members<=nPFC)),:))];
                
                PDistCont_.A2UnonmembersError{iArea}{2}{iFile} = [PDistCont_.A2UnonmembersError{iArea}{2}{iFile};...
                    aggregatefun(A2U(setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC)),:))];
                
                
                D_ = []; % member to member (HP to PFC exclusively)
                for i = temp.members(temp.members<=nPFC)
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2MError{iArea}{1}{iFile} = [PDistCont_.U2U_M2MError{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC member to HP nonmember
                for i = temp.members(temp.members<=nPFC)
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2NError{iArea}{1}{iFile} = [PDistCont_.U2U_M2NError{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP member
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2NError{iArea}{2}{iFile} = [PDistCont_.U2U_M2NError{iArea}{2}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP nonmember
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_N2NError{iArea}{1}{iFile} = [PDistCont_.U2U_N2NError{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2OError{iArea}{1}{iFile} = [PDistCont_.U2U_M2OError{iArea}{1}{iFile}; aggregatefun(D_)];
                
            end
        end
        
        % Compare correct and error trial results
        if iArea < 3
            a = PDistCont_.A2Umembers{iArea}{iFile};
            b = PDistCont_.A2UmembersError{iArea}{iFile};
            PDistCont_.A2UmembersDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.A2Unonmembers{iArea}{iFile};
            b = PDistCont_.A2UnonmembersError{iArea}{iFile};
            PDistCont_.A2UnonmembersDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2M{iArea}{iFile};
            b = PDistCont_.U2U_M2MError{iArea}{iFile};
            PDistCont_.U2U_M2MDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2N{iArea}{iFile};
            b = PDistCont_.U2U_M2NError{iArea}{iFile};
            PDistCont_.U2U_M2NDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_N2N{iArea}{iFile};
            b = PDistCont_.U2U_N2NError{iArea}{iFile};
            PDistCont_.U2U_N2NDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2O{iArea}{iFile};
            b = PDistCont_.U2U_M2OError{iArea}{iFile};
            PDistCont_.U2U_M2ODelta{iArea}{iFile} =  compfunction(a,b);
            
            
        else
            a = PDistCont_.U2U_M2M{iArea}{1}{iFile};
            b = PDistCont_.U2U_M2MError{iArea}{1}{iFile};
            PDistCont_.U2U_M2MDelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_N2N{iArea}{1}{iFile};
            b = PDistCont_.U2U_N2NError{iArea}{1}{iFile};
            PDistCont_.U2U_N2NDelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2O{iArea}{1}{iFile};
            b = PDistCont_.U2U_M2OError{iArea}{1}{iFile};
            PDistCont_.U2U_M2ODelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            for i=1:2
                a = PDistCont_.A2Umembers{iArea}{i}{iFile};
                b = PDistCont_.A2UmembersError{iArea}{i}{iFile};
                PDistCont_.A2UmembersDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
                a = PDistCont_.A2Unonmembers{iArea}{i}{iFile};
                b = PDistCont_.A2UnonmembersError{iArea}{i}{iFile};
                PDistCont_.A2UnonmembersDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
                a = PDistCont_.U2U_M2N{iArea}{i}{iFile};
                b = PDistCont_.U2U_M2NError{iArea}{i}{iFile};
                PDistCont_.U2U_M2NDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
            end
        end
    end
    
end
%% Distribution of pairwise distances between D' scores of units/units and units/assemblies - time blocks
clear PDistCont_ temp
compfunction = @(a,b) (b-a)%./(a+b); % a=correct b=error
aggregatefun = @(x) (x);
% aggregatefun = @(x) mean(x,1);

% normfun = @(x) (x);
normfun = @(x) zscore(x);
% normfun = @(x) x - repmat(nanmean(x),size(x,1),1);

blocksize = 10;
nData = size(D.TS{1}{1},1);
blockOffsets =(blocksize:blocksize:nData)-blocksize+1;

% distanceMetric = 'Hamming';
% bins = 0:0.005:1;
% bins2 = -5:.25:5;

% distanceMetric = 'cosine';
% bins = 0:0.1:1.5;
% bins2 = -5:.25:5;
distanceMetric = 'euclidean';
bins = 0:0.05:50;
bins2 = -50:2:50;
% bins2 = -1:0.1:1;

for iArea = 1:3
    for iFile = 1:length(D_Ass.TS{iArea})
        [iFile,iArea]
        nPFC = size(D.TS{1}{iFile},2);
        
        % Get data correct
        temp.TS_Ass   = D_Ass.TS{iArea}{iFile};
        if iArea < 3
            temp.TS_Units = D.TS{iArea}{iFile};
        else
            temp.TS_Units = [D.TS{1}{iFile}, D.TS{2}{iFile}];
        end
        % Get data error
        temp.TS_Ass   = D_Ass.TS_err{iArea}{iFile};
        if iArea < 3
            temp.TS_Units = D.TS_err{iArea}{iFile};
        else
            temp.TS_Units = [D.TS_err{1}{iFile},D.TS_err{2}{iFile}];
        end
        
        % Preallocate correct
        if iArea < 3
            PDistCont_.A2Umembers{iArea}{iFile}    = []; % mean Assembly to unit distance (members)
            PDistCont_.A2Unonmembers{iArea}{iFile} = [];
            PDistCont_.U2U_M2M{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-member)
            PDistCont_.U2U_M2N{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-nonmember)
            PDistCont_.U2U_N2N{iArea}{iFile} = [];       % mean inter-unit coherence in information (nonmember-to-nonmember)
            PDistCont_.U2U_M2O{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-other ass member)
        else
            PDistCont_.A2Umembers{iArea}{1}{iFile}    = [];  PDistCont_.A2Umembers{iArea}{2}{iFile}    = [];
            PDistCont_.A2Unonmembers{iArea}{1}{iFile} = [];  PDistCont_.A2Unonmembers{iArea}{2}{iFile} = [];
            PDistCont_.U2U_M2M{iArea}{1}{iFile} = [];
            PDistCont_.U2U_N2N{iArea}{1}{iFile} = [];
            PDistCont_.U2U_M2N{iArea}{1}{iFile} = [];    PDistCont_.U2U_M2N{iArea}{2}{iFile} = [];    %member to nonmember {Mpfc2Nhp}{Mhp2Npfc}
            PDistCont_.U2U_M2O{iArea}{1}{iFile} = [];   % mean inter-unit coherence in information (member-to-other ass member)
        end
        % Preallocate error
        if iArea < 3
            PDistCont_.A2UmembersError{iArea}{iFile}    = []; % mean Assembly to unit distance (members)
            PDistCont_.A2UnonmembersError{iArea}{iFile} = [];
            PDistCont_.U2U_M2MError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-member)
            PDistCont_.U2U_M2NError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-nonmember)
            PDistCont_.U2U_N2NError{iArea}{iFile} = [];       % mean inter-unit coherence in information (nonmember-to-nonmember)
            PDistCont_.U2U_M2OError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-other ass member)
        else
            PDistCont_.A2UmembersError{iArea}{1}{iFile}    = [];            PDistCont_.A2UmembersError{iArea}{2}{iFile}    = [];
            PDistCont_.A2UnonmembersError{iArea}{1}{iFile} = [];            PDistCont_.A2UnonmembersError{iArea}{2}{iFile} = [];
            PDistCont_.U2U_M2MError{iArea}{1}{iFile} = [];
            PDistCont_.U2U_N2NError{iArea}{1}{iFile} = [];
            PDistCont_.U2U_M2NError{iArea}{1}{iFile} = [];    PDistCont_.U2U_M2NError{iArea}{2}{iFile} = [];    %member to nonmember {Mpfc2Nhp}{Mhp2Npfc}
            PDistCont_.U2U_M2OError{iArea}{1}{iFile} = [];   % mean inter-unit coherence in information (member-to-other ass member)

        end
        
        
        % Process correct
        for iAss = 1:size(temp.TS_Ass,2)
            clear coeff
            temp.members = D_Ass.units{iArea}{iFile}{iAss};
            temp.Othermembers = setdiff(cell2mat(D_Ass.units{iArea}{iFile}),temp.members);
            
            input = ([temp.TS_Ass(:,iAss),temp.TS_Units]);
            input = normfun(input);
            
            for iBlock  = 1:length(blockOffsets)
                t= blockOffsets(iBlock):(blockOffsets(iBlock)+blocksize-1);
                coeff(:,:,iBlock) = triu(squareform(pdist(input(t,:)',distanceMetric)));
            end
            
            
            A2U = squeeze(coeff(1,2:end,:));
            U2U = coeff(2:end,2:end,:);
            if iArea < 3
                % Assembly to unit distances
                PDistCont_.A2Umembers{iArea}{iFile}    = [PDistCont_.A2Umembers{iArea}{iFile};     aggregatefun(A2U(temp.members,:))];
                PDistCont_.A2Unonmembers{iArea}{iFile} = [PDistCont_.A2Unonmembers{iArea}{iFile};  aggregatefun(A2U(setdiff(1:size(A2U,1),temp.members),:))];
                
                % Unit to unit distances
                
                D_ = []; % member to member
                for i = temp.members
                    for j = temp.members
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2M{iArea}{iFile} = [PDistCont_.U2U_M2M{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to nonmember
                for i = temp.members
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2N{iArea}{iFile} = [PDistCont_.U2U_M2N{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % nonmember to nonmember
                for i =  setdiff(1:size(A2U,1),temp.members)
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_N2N{iArea}{iFile} = [PDistCont_.U2U_N2N{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2O{iArea}{iFile} = [PDistCont_.U2U_M2O{iArea}{iFile}; aggregatefun(D_)];
                
            else
                PDistCont_.A2Umembers{iArea}{1}{iFile} = [PDistCont_.A2Umembers{iArea}{1}{iFile};     aggregatefun(A2U(temp.members(temp.members<=nPFC),:))];
                PDistCont_.A2Umembers{iArea}{2}{iFile} = [PDistCont_.A2Umembers{iArea}{2}{iFile};     aggregatefun(A2U(temp.members(temp.members>nPFC),:))];
                
                PDistCont_.A2Unonmembers{iArea}{1}{iFile} = [PDistCont_.A2Unonmembers{iArea}{1}{iFile};...
                    aggregatefun(A2U(setdiff(1:nPFC,temp.members(temp.members<=nPFC)),:))];
                
                PDistCont_.A2Unonmembers{iArea}{2}{iFile} = [PDistCont_.A2Unonmembers{iArea}{2}{iFile};...
                    aggregatefun(A2U(setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC)),:))];
                
                D_ = []; % member to member (HP to PFC exclusively)
                for i = temp.members(temp.members<=nPFC)
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2M{iArea}{1}{iFile} = [PDistCont_.U2U_M2M{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC member to HP nonmember
                for i = temp.members(temp.members<=nPFC)
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2N{iArea}{1}{iFile} = [PDistCont_.U2U_M2N{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP member
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2N{iArea}{2}{iFile} = [PDistCont_.U2U_M2N{iArea}{2}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP nonmember
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_N2N{iArea}{1}{iFile} = [PDistCont_.U2U_N2N{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2O{iArea}{1}{iFile} = [PDistCont_.U2U_M2O{iArea}{1}{iFile}; aggregatefun(D_)];
                
            end
        end
        
        % Process error
        for iAss = 1:size(temp.TS_Ass,2)
            temp.members = D_Ass.units{iArea}{iFile}{iAss};
            temp.Othermembers = setdiff(cell2mat(D_Ass.units{iArea}{iFile}),temp.members);
            input = ([temp.TS_Ass(:,iAss),temp.TS_Units]);
            input = normfun(input);
            for iBlock  = 1:length(blockOffsets)
                t= blockOffsets(iBlock):(blockOffsets(iBlock)+blocksize-1);
                coeff(:,:,iBlock) = triu(squareform(pdist(input(t,:)',distanceMetric)));
            end
            
            A2U = squeeze(coeff(1,2:end,:));
            U2U = coeff(2:end,2:end,:);
            if iArea < 3
                PDistCont_.A2UmembersError{iArea}{iFile}    = [PDistCont_.A2UmembersError{iArea}{iFile};     aggregatefun(A2U(temp.members,:))];
                PDistCont_.A2UnonmembersError{iArea}{iFile} = [PDistCont_.A2UnonmembersError{iArea}{iFile};  aggregatefun(A2U(setdiff(1:size(A2U,1),temp.members),:))];
                
                % Unit to unit distances
                
                D_ = []; % member to member
                for i = temp.members
                    for j = temp.members
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2MError{iArea}{iFile} = [PDistCont_.U2U_M2MError{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to nonmember
                for i = temp.members
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2NError{iArea}{iFile} = [PDistCont_.U2U_M2NError{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % nonmember to nonmember
                for i =  setdiff(1:size(A2U,1),temp.members)
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_N2NError{iArea}{iFile} = [PDistCont_.U2U_N2NError{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2OError{iArea}{iFile} = [PDistCont_.U2U_M2OError{iArea}{iFile}; aggregatefun(D_)];
                
                
            else
                PDistCont_.A2UmembersError{iArea}{1}{iFile} = [PDistCont_.A2UmembersError{iArea}{1}{iFile};     aggregatefun(A2U(temp.members(temp.members<=nPFC),:))];
                PDistCont_.A2UmembersError{iArea}{2}{iFile} = [PDistCont_.A2UmembersError{iArea}{2}{iFile};     aggregatefun(A2U(temp.members(temp.members>nPFC),:))];
                
                PDistCont_.A2UnonmembersError{iArea}{1}{iFile} = [PDistCont_.A2UnonmembersError{iArea}{1}{iFile};...
                    aggregatefun(A2U(setdiff(1:nPFC,temp.members(temp.members<=nPFC)),:))];
                
                PDistCont_.A2UnonmembersError{iArea}{2}{iFile} = [PDistCont_.A2UnonmembersError{iArea}{2}{iFile};...
                    aggregatefun(A2U(setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC)),:))];
                
                
                D_ = []; % member to member (HP to PFC exclusively)
                for i = temp.members(temp.members<=nPFC)
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2MError{iArea}{1}{iFile} = [PDistCont_.U2U_M2MError{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC member to HP nonmember
                for i = temp.members(temp.members<=nPFC)
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2NError{iArea}{1}{iFile} = [PDistCont_.U2U_M2NError{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP member
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2NError{iArea}{2}{iFile} = [PDistCont_.U2U_M2NError{iArea}{2}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP nonmember
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_N2NError{iArea}{1}{iFile} = [PDistCont_.U2U_N2NError{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2OError{iArea}{1}{iFile} = [PDistCont_.U2U_M2OError{iArea}{1}{iFile}; aggregatefun(D_)];
                
            end
        end
        
        % Compare correct and error trial results
        if iArea < 3
            a = PDistCont_.A2Umembers{iArea}{iFile};
            b = PDistCont_.A2UmembersError{iArea}{iFile};
            PDistCont_.A2UmembersDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.A2Unonmembers{iArea}{iFile};
            b = PDistCont_.A2UnonmembersError{iArea}{iFile};
            PDistCont_.A2UnonmembersDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2M{iArea}{iFile};
            b = PDistCont_.U2U_M2MError{iArea}{iFile};
            PDistCont_.U2U_M2MDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2N{iArea}{iFile};
            b = PDistCont_.U2U_M2NError{iArea}{iFile};
            PDistCont_.U2U_M2NDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_N2N{iArea}{iFile};
            b = PDistCont_.U2U_N2NError{iArea}{iFile};
            PDistCont_.U2U_N2NDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2O{iArea}{iFile};
            b = PDistCont_.U2U_M2OError{iArea}{iFile};
            PDistCont_.U2U_M2ODelta{iArea}{iFile} =  compfunction(a,b);
            
            
        else
            a = PDistCont_.U2U_M2M{iArea}{1}{iFile};
            b = PDistCont_.U2U_M2MError{iArea}{1}{iFile};
            PDistCont_.U2U_M2MDelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_N2N{iArea}{1}{iFile};
            b = PDistCont_.U2U_N2NError{iArea}{1}{iFile};
            PDistCont_.U2U_N2NDelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2O{iArea}{1}{iFile};
            b = PDistCont_.U2U_M2OError{iArea}{1}{iFile};
            PDistCont_.U2U_M2ODelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            for i=1:2
                a = PDistCont_.A2Umembers{iArea}{i}{iFile};
                b = PDistCont_.A2UmembersError{iArea}{i}{iFile};
                PDistCont_.A2UmembersDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
                a = PDistCont_.A2Unonmembers{iArea}{i}{iFile};
                b = PDistCont_.A2UnonmembersError{iArea}{i}{iFile};
                PDistCont_.A2UnonmembersDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
                a = PDistCont_.U2U_M2N{iArea}{i}{iFile};
                b = PDistCont_.U2U_M2NError{iArea}{i}{iFile};
                PDistCont_.U2U_M2NDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
            end
        end
    end
    
end
%% Distribution of pairwise distances between D' scores of units/units and units/assemblies - continuous/time blocks shared variance
clear PDistCont_ temp
compfunction = @(a,b) (b-a)%./(a+b); % a=correct b=error
aggregatefun = @(x) (x);
% aggregatefun = @(x) mean(x,1);

% normfun = @(x) (x);
normfun = @(x) zscore(x);
% normfun = @(x) x - repmat(nanmean(x),size(x,1),1);

blocksize = 10;
nData = size(D.TS{1}{1},1);
blockOffsets =(blocksize:blocksize:nData)-blocksize+1;
runBlocks=true;
if runBlocks
    tbBlocks = (1:length(blockOffsets))./bw/40
else
    tbBlocks = (1:Ltr)*bw;
end
    
% distanceMetric = 'Hamming';
% bins = 0:0.005:1;
% bins2 = -5:.25:5;

% distanceMetric = 'cosine';
% bins = 0:0.1:1.5;
% bins2 = -5:.25:5;
distanceMetric = 'euclidean';
bins = 0:0.05:50;
bins2 = -50:2:50;
% bins2 = -1:0.1:1;

    for iArea = 1:3
    for iFile = 1:length(D_Ass.TS{iArea})
        [iFile,iArea]
        nPFC = size(D.TS{1}{iFile},2);
        
        % Get data correct
        temp.TS_Ass   = D_Ass.TS{iArea}{iFile};
        if iArea < 3
            temp.TS_Units = D.TS{iArea}{iFile};
        else
            temp.TS_Units = [D.TS{1}{iFile}, D.TS{2}{iFile}];
        end
        % Get data error
        temp.TS_AssErr   = D_Ass.TS_err{iArea}{iFile};
        if iArea < 3
            temp.TS_UnitsErr = D.TS_err{iArea}{iFile};
        else
            temp.TS_UnitsErr = [D.TS_err{1}{iFile},D.TS_err{2}{iFile}];
        end
         % Preallocate correct
        if iArea < 3
            PDistCont_.A2Umembers{iArea}{iFile}    = []; % mean Assembly to unit distance (members)
            PDistCont_.A2Unonmembers{iArea}{iFile} = [];
            PDistCont_.U2U_M2M{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-member)
            PDistCont_.U2U_M2N{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-nonmember)
            PDistCont_.U2U_N2N{iArea}{iFile} = [];       % mean inter-unit coherence in information (nonmember-to-nonmember)
            PDistCont_.U2U_M2O{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-other ass member)
        else
            PDistCont_.A2Umembers{iArea}{1}{iFile}    = [];  PDistCont_.A2Umembers{iArea}{2}{iFile}    = [];
            PDistCont_.A2Unonmembers{iArea}{1}{iFile} = [];  PDistCont_.A2Unonmembers{iArea}{2}{iFile} = [];
            PDistCont_.U2U_M2M{iArea}{1}{iFile} = [];
            PDistCont_.U2U_N2N{iArea}{1}{iFile} = [];
            PDistCont_.U2U_M2N{iArea}{1}{iFile} = [];    PDistCont_.U2U_M2N{iArea}{2}{iFile} = [];    %member to nonmember {Mpfc2Nhp}{Mhp2Npfc}
            PDistCont_.U2U_M2O{iArea}{1}{iFile} = [];   % mean inter-unit coherence in information (member-to-other ass member)
        end
        % Preallocate error
        if iArea < 3
            PDistCont_.A2UmembersError{iArea}{iFile}    = []; % mean Assembly to unit distance (members)
            PDistCont_.A2UnonmembersError{iArea}{iFile} = [];
            PDistCont_.U2U_M2MError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-member)
            PDistCont_.U2U_M2NError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-nonmember)
            PDistCont_.U2U_N2NError{iArea}{iFile} = [];       % mean inter-unit coherence in information (nonmember-to-nonmember)
            PDistCont_.U2U_M2OError{iArea}{iFile} = [];       % mean inter-unit coherence in information (member-to-other ass member)
        else
            PDistCont_.A2UmembersError{iArea}{1}{iFile}    = [];            PDistCont_.A2UmembersError{iArea}{2}{iFile}    = [];
            PDistCont_.A2UnonmembersError{iArea}{1}{iFile} = [];            PDistCont_.A2UnonmembersError{iArea}{2}{iFile} = [];
            PDistCont_.U2U_M2MError{iArea}{1}{iFile} = [];
            PDistCont_.U2U_N2NError{iArea}{1}{iFile} = [];
            PDistCont_.U2U_M2NError{iArea}{1}{iFile} = [];    PDistCont_.U2U_M2NError{iArea}{2}{iFile} = [];    %member to nonmember {Mpfc2Nhp}{Mhp2Npfc}
            PDistCont_.U2U_M2OError{iArea}{1}{iFile} = [];   % mean inter-unit coherence in information (member-to-other ass member)

        end
        
        
         % Process 
        for iAss = 1:size(temp.TS_Ass,2)
            
            temp.members = D_Ass.units{iArea}{iFile}{iAss};
            temp.Othermembers = setdiff(cell2mat(D_Ass.units{iArea}{iFile}),temp.members);
            
            % process correct
            clear coeff

            input = ([temp.TS_Ass(:,iAss),temp.TS_Units;...
                      temp.TS_AssErr(:,iAss),temp.TS_UnitsErr]);
            input = normfun(input);
            input = input(1:Ltr,:);
            
            if runBlocks
                for iBlock  = 1:length(blockOffsets)
                    t= blockOffsets(iBlock):(blockOffsets(iBlock)+blocksize-1);
                    coeff(:,:,iBlock) = triu(squareform(pdist(input(t,:)',distanceMetric)));
                end
            else
                for t =1:size(input,1)
                    coeff(:,:,t) = triu(squareform(pdist(input(t,:)',distanceMetric)));
                end
            end
            
            A2U = squeeze(coeff(1,2:end,:));
            U2U = coeff(2:end,2:end,:);
            if iArea < 3
                % Assembly to unit distances
                PDistCont_.A2Umembers{iArea}{iFile}    = [PDistCont_.A2Umembers{iArea}{iFile};     aggregatefun(A2U(temp.members,:))];
                PDistCont_.A2Unonmembers{iArea}{iFile} = [PDistCont_.A2Unonmembers{iArea}{iFile};  aggregatefun(A2U(setdiff(1:size(A2U,1),temp.members),:))];
                 % Unit to unit distances
                
                D_ = []; % member to member
                for i = temp.members
                    for j = temp.members
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2M{iArea}{iFile} = [PDistCont_.U2U_M2M{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to nonmember
                for i = temp.members
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2N{iArea}{iFile} = [PDistCont_.U2U_M2N{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % nonmember to nonmember
                for i =  setdiff(1:size(A2U,1),temp.members)
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_N2N{iArea}{iFile} = [PDistCont_.U2U_N2N{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2O{iArea}{iFile} = [PDistCont_.U2U_M2O{iArea}{iFile}; aggregatefun(D_)];
                
            else
                PDistCont_.A2Umembers{iArea}{1}{iFile} = [PDistCont_.A2Umembers{iArea}{1}{iFile};     aggregatefun(A2U(temp.members(temp.members<=nPFC),:))];
                PDistCont_.A2Umembers{iArea}{2}{iFile} = [PDistCont_.A2Umembers{iArea}{2}{iFile};     aggregatefun(A2U(temp.members(temp.members>nPFC),:))];
                
                PDistCont_.A2Unonmembers{iArea}{1}{iFile} = [PDistCont_.A2Unonmembers{iArea}{1}{iFile};...
                    aggregatefun(A2U(setdiff(1:nPFC,temp.members(temp.members<=nPFC)),:))];
                
                PDistCont_.A2Unonmembers{iArea}{2}{iFile} = [PDistCont_.A2Unonmembers{iArea}{2}{iFile};...
                    aggregatefun(A2U(setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC)),:))];
                
                D_ = []; % member to member (HP to PFC exclusively)
                for i = temp.members(temp.members<=nPFC)
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2M{iArea}{1}{iFile} = [PDistCont_.U2U_M2M{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC member to HP nonmember
                for i = temp.members(temp.members<=nPFC)
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2N{iArea}{1}{iFile} = [PDistCont_.U2U_M2N{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP member
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2N{iArea}{2}{iFile} = [PDistCont_.U2U_M2N{iArea}{2}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP nonmember
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_N2N{iArea}{1}{iFile} = [PDistCont_.U2U_N2N{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2O{iArea}{1}{iFile} = [PDistCont_.U2U_M2O{iArea}{1}{iFile}; aggregatefun(D_)];
                
            end
        
            % Process error
            clear coeff

            input = ([temp.TS_Ass(:,iAss),temp.TS_Units;...
                      temp.TS_AssErr(:,iAss),temp.TS_UnitsErr]);
            input = normfun(input);
            input = input((Ltr+1):2*Ltr,:);
            if runBlocks
                for iBlock  = 1:length(blockOffsets)
                    t= blockOffsets(iBlock):(blockOffsets(iBlock)+blocksize-1);
                    coeff(:,:,iBlock) = triu(squareform(pdist(input(t,:)',distanceMetric)));
                end
            else
                for t =1:size(input,1)
                    coeff(:,:,t) = triu(squareform(pdist(input(t,:)',distanceMetric)));
                end
            end
            A2U = squeeze(coeff(1,2:end,:));
            U2U = coeff(2:end,2:end,:);
            if iArea < 3
                PDistCont_.A2UmembersError{iArea}{iFile}    = [PDistCont_.A2UmembersError{iArea}{iFile};     aggregatefun(A2U(temp.members,:))];
                PDistCont_.A2UnonmembersError{iArea}{iFile} = [PDistCont_.A2UnonmembersError{iArea}{iFile};  aggregatefun(A2U(setdiff(1:size(A2U,1),temp.members),:))];
                
                % Unit to unit distances
                
                D_ = []; % member to member
                for i = temp.members
                    for j = temp.members
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2MError{iArea}{iFile} = [PDistCont_.U2U_M2MError{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to nonmember
                for i = temp.members
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2NError{iArea}{iFile} = [PDistCont_.U2U_M2NError{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % nonmember to nonmember
                for i =  setdiff(1:size(A2U,1),temp.members)
                    for j = setdiff(1:size(A2U,1),temp.members)
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_N2NError{iArea}{iFile} = [PDistCont_.U2U_N2NError{iArea}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2OError{iArea}{iFile} = [PDistCont_.U2U_M2OError{iArea}{iFile}; aggregatefun(D_)];
                
                
            else
                PDistCont_.A2UmembersError{iArea}{1}{iFile} = [PDistCont_.A2UmembersError{iArea}{1}{iFile};     aggregatefun(A2U(temp.members(temp.members<=nPFC),:))];
                PDistCont_.A2UmembersError{iArea}{2}{iFile} = [PDistCont_.A2UmembersError{iArea}{2}{iFile};     aggregatefun(A2U(temp.members(temp.members>nPFC),:))];
                
                PDistCont_.A2UnonmembersError{iArea}{1}{iFile} = [PDistCont_.A2UnonmembersError{iArea}{1}{iFile};...
                    aggregatefun(A2U(setdiff(1:nPFC,temp.members(temp.members<=nPFC)),:))];
                
                PDistCont_.A2UnonmembersError{iArea}{2}{iFile} = [PDistCont_.A2UnonmembersError{iArea}{2}{iFile};...
                    aggregatefun(A2U(setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC)),:))];
                
                
                D_ = []; % member to member (HP to PFC exclusively)
                for i = temp.members(temp.members<=nPFC)
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2MError{iArea}{1}{iFile} = [PDistCont_.U2U_M2MError{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC member to HP nonmember
                for i = temp.members(temp.members<=nPFC)
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2NError{iArea}{1}{iFile} = [PDistCont_.U2U_M2NError{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP member
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = temp.members(temp.members>nPFC)
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_M2NError{iArea}{2}{iFile} = [PDistCont_.U2U_M2NError{iArea}{2}{iFile}; aggregatefun(D_)];
                
                D_ = []; % PFC nonmember to HP nonmember
                for i = setdiff(1:nPFC,temp.members(temp.members<=nPFC))
                    for j = setdiff((nPFC+1):size(A2U,1),temp.members(temp.members>nPFC))
                        D_ = [D_;squeeze(U2U(i,j,:))'];
                    end
                end
                PDistCont_.U2U_N2NError{iArea}{1}{iFile} = [PDistCont_.U2U_N2NError{iArea}{1}{iFile}; aggregatefun(D_)];
                
                D_ = []; % member to other assembly member
                for i =  temp.members
                    for j = temp.Othermembers
                        if i<j
                            D_=[D_;squeeze(U2U(i,j,:))'];
                        end
                    end
                end
                PDistCont_.U2U_M2OError{iArea}{1}{iFile} = [PDistCont_.U2U_M2OError{iArea}{1}{iFile}; aggregatefun(D_)];
                
            end
        end
        % Compare correct and error trial results
        if iArea < 3
            a = PDistCont_.A2Umembers{iArea}{iFile};
            b = PDistCont_.A2UmembersError{iArea}{iFile};
            PDistCont_.A2UmembersDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.A2Unonmembers{iArea}{iFile};
            b = PDistCont_.A2UnonmembersError{iArea}{iFile};
            PDistCont_.A2UnonmembersDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2M{iArea}{iFile};
            b = PDistCont_.U2U_M2MError{iArea}{iFile};
            PDistCont_.U2U_M2MDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2N{iArea}{iFile};
            b = PDistCont_.U2U_M2NError{iArea}{iFile};
            PDistCont_.U2U_M2NDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_N2N{iArea}{iFile};
            b = PDistCont_.U2U_N2NError{iArea}{iFile};
            PDistCont_.U2U_N2NDelta{iArea}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2O{iArea}{iFile};
            b = PDistCont_.U2U_M2OError{iArea}{iFile};
            PDistCont_.U2U_M2ODelta{iArea}{iFile} =  compfunction(a,b);
            
            
        else
            a = PDistCont_.U2U_M2M{iArea}{1}{iFile};
            b = PDistCont_.U2U_M2MError{iArea}{1}{iFile};
            PDistCont_.U2U_M2MDelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_N2N{iArea}{1}{iFile};
            b = PDistCont_.U2U_N2NError{iArea}{1}{iFile};
            PDistCont_.U2U_N2NDelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            a = PDistCont_.U2U_M2O{iArea}{1}{iFile};
            b = PDistCont_.U2U_M2OError{iArea}{1}{iFile};
            PDistCont_.U2U_M2ODelta{iArea}{1}{iFile} =  compfunction(a,b);
            
            for i=1:2
                a = PDistCont_.A2Umembers{iArea}{i}{iFile};
                b = PDistCont_.A2UmembersError{iArea}{i}{iFile};
                PDistCont_.A2UmembersDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
                a = PDistCont_.A2Unonmembers{iArea}{i}{iFile};
                b = PDistCont_.A2UnonmembersError{iArea}{i}{iFile};
                PDistCont_.A2UnonmembersDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
                a = PDistCont_.U2U_M2N{iArea}{i}{iFile};
                b = PDistCont_.U2U_M2NError{iArea}{i}{iFile};
                PDistCont_.U2U_M2NDelta{iArea}{i}{iFile} =  compfunction(a,b);
                
            end
        end
    end
    
end
%% Plot average pairwise distances between units/assemblies - members vs. non-members
figure('name','Distance-bewteen Units and Assemblies (correct)');
for iArea = 1:3
    
    subplot(1,3,iArea); hold on
    title(Areas{iArea})
    if iArea<3
        mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2Umembers{iArea}(~isempty_cell(PDistCont_.A2Umembers{iArea}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        x = (1:size(mean_,2)).*bw;
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,'k',1)
        
        mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2Unonmembers{iArea}(~isempty_cell(PDistCont_.A2Unonmembers{iArea}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,'k',0.3)
    else
        for iArea_ =1:2
            mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2Umembers{iArea}{iArea_}(~isempty_cell(PDistCont_.A2Umembers{iArea}{iArea_}))','UniformOutput',false);
            mean_ = cell2mat(mean_);
            ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,color_{iArea_},1)
            
            mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2Unonmembers{iArea}{iArea_}(~isempty_cell(PDistCont_.A2Unonmembers{iArea}{iArea_}))','UniformOutput',false);
            mean_ = cell2mat(mean_);
            ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,color_{iArea_},0.6)
           
        end
    end
    
    axis([0 Inf 0 20])
end
%% Plot average pairwise distances between units/assemblies - correct vs. error
figure('name','Distance-bewteen Units and Assemblies (correct vs error)');
for iArea = 1:3
    
    subplot(1,3,iArea); hold on
    title(Areas{iArea})
    if iArea<3
        mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2Umembers{iArea}(~isempty_cell(PDistCont_.A2Umembers{iArea}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,'k',1)
        
        mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2UmembersError{iArea}(~isempty_cell(PDistCont_.A2UmembersError{iArea}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,'k',0.3)
    
       axis([0 Inf 0 10]) 
    end
    iArea = 3;
    for iArea_ =1:2
            subplot(2,3,iArea_*3); hold on
            
            mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2Umembers{iArea}{iArea_}(~isempty_cell(PDistCont_.A2Umembers{iArea}{iArea_}))','UniformOutput',false);
            mean_ = cell2mat(mean_);
            ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,color_{iArea_},1)
            
            mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2UmembersError{iArea}{iArea_}(~isempty_cell(PDistCont_.A2UmembersError{iArea}{iArea_}))','UniformOutput',false);
            mean_ = cell2mat(mean_);
            ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,color_{iArea_},0.3)
               axis([0 Inf 0 10])
        end

end
%% Plot average pairwise distances between units/assemblies - correct vs. error - comparison
figure('name','Distance-bewteen Units and Assemblies (correct vs error)');
for iArea = 1:3
    
    subplot(1,3,iArea); hold on
    title(Areas{iArea})
    if iArea<3
        mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2UmembersDelta{iArea}(~isempty_cell(PDistCont_.A2UmembersDelta{iArea}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,'k',1)
        
         mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2UnonmembersDelta{iArea}(~isempty_cell(PDistCont_.A2UnonmembersDelta{iArea}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,0.6*[1 1 1],1)
        
        
      
    else
        for iArea_ =1:2
            mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2UmembersDelta{iArea}{iArea_}(~isempty_cell(PDistCont_.A2UmembersDelta{iArea}{iArea_}))','UniformOutput',false);
            mean_ = cell2mat(mean_);
            ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,color_{iArea_},1)
            
              mean_ = cellfun(@(x) mean(x,1),PDistCont_.A2UnonmembersDelta{iArea}{iArea_}(~isempty_cell(PDistCont_.A2UnonmembersDelta{iArea}{iArea_}))','UniformOutput',false);
            mean_ = cell2mat(mean_);
            ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,color_{iArea_},0.6)
            

           
        end
    end
    
    axis([0 20 -5 5])
end
%% plot average pairwise distances between units - members vs. non-members
figure('name', 'Unit-Unit distance (correct)');
for iArea = 1:3
    subplot(1,3,iArea); hold on
    title([Areas{iArea}])
    if iArea<3
        mean_ = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2M{iArea}(~isempty_cell(PDistCont_.U2U_M2M{iArea}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,0*[1 1 1],1)
        
%         mean_ = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2N{iArea}(~isempty_cell(PDistCont_.U2U_M2N{iArea}))','UniformOutput',false);
%         mean_ = cell2mat(mean_);
%         ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),x,0.4*[1 1 1],1)        
       
        mean_ = cellfun(@(x) mean(x,1),PDistCont_.U2U_N2N{iArea}(~isempty_cell(PDistCont_.U2U_N2N{iArea}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,0.6*[1 1 1],1)
      
%         legend('Between members','Member to Non-member', 'Between non-members','Location','SouthEast')
        
    else
        mean_ = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2M{iArea}{1}(~isempty_cell(PDistCont_.U2U_M2M{iArea}{1}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,0*[1 1 1],1)
%         
        for iArea_ =1:2
            mean_ = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2N{iArea}{iArea_}(~isempty_cell(PDistCont_.U2U_M2N{iArea}{iArea_}))','UniformOutput',false);
            mean_ = cell2mat(mean_);
            ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,color_{iArea_},0.6)
        end
        mean_ = cellfun(@(x) mean(x,1),PDistCont_.U2U_N2N{iArea}{1}(~isempty_cell(PDistCont_.U2U_N2N{iArea}{1}))','UniformOutput',false);
        mean_ = cell2mat(mean_);
        ciplot(nanmean(mean_)+nansem(mean_),nanmean(mean_)-nansem(mean_),tbBlocks,0.6*[1 1 1],1)
        legend('Between members','PFC Member to HP Non-member','HP Member to PFC Non-member', 'Between non-members','Location','SouthEast')
        
    end
    axis([0 20 0 20])
end
%% plot average pairwise distances between units - correct vs. error
figure('name','inter-unit distance: members')
for iArea = 1:3
    
    subplot(1,3,iArea); hold on
    plot([5 5],[0.8 1.6],'g','LineWidth',1.5)
    plot([15 15],[0.8 1.6],'r','LineWidth',1.5)
    title([Areas{iArea}])
    if iArea<3
        mean_C = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2M{iArea}(~isempty_cell(PDistCont_.U2U_M2M{iArea}))','UniformOutput',false);
        mean_C = cell2mat(mean_C);
        ciplot(nanmean(mean_C)+nansem(mean_C),nanmean(mean_C)-nansem(mean_C),tbBlocks,0*[1 1 1],1)
        
        mean_E = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2MError{iArea}(~isempty_cell(PDistCont_.U2U_M2MError{iArea}))','UniformOutput',false);
        mean_E = cell2mat(mean_E);
        ciplot(nanmean(mean_E)+nansem(mean_E),nanmean(mean_E)-nansem(mean_E),tbBlocks,'r',0.6)
        Y1 = smooth2a(mean_C,0,10);Y1(:,isnan(nansum(Y1)))=[];
        Y2 = smooth2a(mean_E,0,10);Y2(:,isnan(nansum(Y2)))=[];
        [sig,~] = permtest2vec(Y1',Y2',100,0.05);
        
        a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
        plot(tbBlocks,a*1.8,'color','k','LineWidth',5)
    else
        mean_C = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2M{iArea}{1}(~isempty_cell(PDistCont_.U2U_M2M{iArea}{1}))','UniformOutput',false);
        mean_C = cell2mat(mean_C);
        ciplot(nanmean(mean_C)+nansem(mean_C),nanmean(mean_C)-nansem(mean_C),tbBlocks,0*[1 1 1],1)
        
        mean_E = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2MError{iArea}{1}(~isempty_cell(PDistCont_.U2U_M2MError{iArea}{1}))','UniformOutput',false);
        mean_E = cell2mat(mean_E);
        ciplot(nanmean(mean_E)+nansem(mean_E),nanmean(mean_E)-nansem(mean_E),tbBlocks,'r',0.6)
        Y1 = smooth2a(mean_C,0,10);Y1(:,isnan(nansum(Y1)))=[];
        Y2 = smooth2a(mean_E,0,10);Y2(:,isnan(nansum(Y2)))=[];
        [sig,~] = permtest2vec(Y1',Y2',100,0.05);
        a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
        plot(tbBlocks,a*1.8,'color','k','LineWidth',5)
    end
    axis([0 20 0 20])
end

figure('name','inter-unit distance: non-members to members')
for iArea = 1:3
    subplot(1,3,iArea); hold on
    plot([5 5],[0.8 1.6],'g','LineWidth',1.5)
    plot([15 15],[0.8 1.6],'r','LineWidth',1.5)
    title([Areas{iArea}])
    if iArea<3
        mean_C = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2N{iArea}(~isempty_cell(PDistCont_.U2U_M2N{iArea}))','UniformOutput',false);
        mean_C = cell2mat(mean_C);
        ciplot(nanmean(mean_C)+nansem(mean_C),nanmean(mean_C)-nansem(mean_C),tbBlocks,0*[1 1 1],1)
        
        mean_E = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2NError{iArea}(~isempty_cell(PDistCont_.U2U_M2NError{iArea}))','UniformOutput',false);
        mean_E = cell2mat(mean_E);
        ciplot(nanmean(mean_E)+nansem(mean_E),nanmean(mean_E)-nansem(mean_E),tbBlocks,'r',0.6)
        Y1 = smooth2a(mean_C,0,10);Y1(:,isnan(nansum(Y1)))=[];
        Y2 = smooth2a(mean_E,0,10);Y2(:,isnan(nansum(Y2)))=[];
        [sig,~] = permtest2vec(Y1',Y2',100,0.05);
        
        a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
        plot(tbBlocks,a*1.8,'color','k','LineWidth',5)
        
    else
        for iArea_ =1:2
            mean_C = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2N{iArea}{iArea_}(~isempty_cell(PDistCont_.U2U_M2N{iArea}{iArea_}))','UniformOutput',false);
            mean_C = cell2mat(mean_C);
            ciplot(nanmean(mean_C)+nansem(mean_C),nanmean(mean_C)-nansem(mean_C),tbBlocks,color_{iArea_},1)
            
            mean_E = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2NError{iArea}{iArea_}(~isempty_cell(PDistCont_.U2U_M2NError{iArea}{iArea_}))','UniformOutput',false);
            mean_E = cell2mat(mean_E);
            ciplot(nanmean(mean_E)+nansem(mean_E),nanmean(mean_E)-nansem(mean_E),tbBlocks,color_{iArea_},0.6)
             Y1 = smooth2a(mean_C,0,10);Y1(:,isnan(nansum(Y1)))=[];
        Y2 = smooth2a(mean_E,0,10);Y2(:,isnan(nansum(Y2)))=[];
        [sig,~] = permtest2vec(Y1',Y2',100,0.05);
        
        a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
        plot(tbBlocks,a*1.6+0.1*iArea_,'color',color_{iArea_},'LineWidth',5)
        end

    end
    axis([0 20 0 20])
end

iArea =3;
figure
for iArea_ =1:2
    subplot(1,3,iArea_); hold on
    plot([5 5],[0.5 1.6],'g','LineWidth',1.5)
    plot([15 15],[0.5 1.6],'r','LineWidth',1.5)
    title([Areas{iArea}])
    mean_C = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2N{iArea}{iArea_}(~isempty_cell(PDistCont_.U2U_M2N{iArea}{iArea_}))','UniformOutput',false);
    mean_C = cell2mat(mean_C);
    ciplot(nanmean(mean_C)+nansem(mean_C),nanmean(mean_C)-nansem(mean_C),tbBlocks,'k',1)
    
    mean_E = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2NError{iArea}{iArea_}(~isempty_cell(PDistCont_.U2U_M2NError{iArea}{iArea_}))','UniformOutput',false);
    mean_E = cell2mat(mean_E);
    ciplot(nanmean(mean_E)+nansem(mean_E),nanmean(mean_E)-nansem(mean_E),tbBlocks,'r',0.6)
    Y1 = smooth2a(mean_C,0,5);Y1(:,isnan(nansum(Y1)))=[];
    Y2 = smooth2a(mean_E,0,5);Y2(:,isnan(nansum(Y2)))=[];
    [sig,~] = permtest2vec(Y1',Y2',100,0.05);
    
    a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
    plot(tbBlocks,a*1.9,'color',color_{iArea_},'LineWidth',5)
    axis([0 20 0 20])
end



figure('name','inter-unit distance: non-members')
for iArea = 1:3
    
    subplot(1,3,iArea);hold on
    plot([5 5],[0.8 1.6],'g','LineWidth',1.5)
    plot([15 15],[0.8 1.6],'r','LineWidth',1.5)
    title([Areas{iArea}])
    if iArea<3
        mean_C = cellfun(@(x) mean(x,1),PDistCont_.U2U_N2N{iArea}(~isempty_cell(PDistCont_.U2U_N2N{iArea}))','UniformOutput',false);
        mean_C = cell2mat(mean_C);
        ciplot(nanmean(mean_C)+nansem(mean_C),nanmean(mean_C)-nansem(mean_C),tbBlocks,0*[1 1 1],1)
        
        mean_E = cellfun(@(x) mean(x,1),PDistCont_.U2U_N2NError{iArea}(~isempty_cell(PDistCont_.U2U_N2NError{iArea}))','UniformOutput',false);
        mean_E = cell2mat(mean_E);
        ciplot(nanmean(mean_E)+nansem(mean_E),nanmean(mean_E)-nansem(mean_E),tbBlocks,'r',0.6)
        Y1 = smooth2a(mean_C,0,10);Y1(:,isnan(nansum(Y1)))=[];
        Y2 = smooth2a(mean_E,0,10);Y2(:,isnan(nansum(Y2)))=[];
        [sig,~] = permtest2vec(Y1',Y2',100,0.05);
        a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
        plot(tbBlocks,a*1.8,'color','k','LineWidth',5)
    else
        mean_C = cellfun(@(x) mean(x,1),PDistCont_.U2U_N2N{iArea}{1}(~isempty_cell(PDistCont_.U2U_N2N{iArea}{1}))','UniformOutput',false);
        mean_C = cell2mat(mean_C);
        ciplot(nanmean(mean_C)+nansem(mean_C),nanmean(mean_C)-nansem(mean_C),tbBlocks,0*[1 1 1],1)
        
        mean_E = cellfun(@(x) mean(x,1),PDistCont_.U2U_N2NError{iArea}{1}(~isempty_cell(PDistCont_.U2U_N2NError{iArea}{1}))','UniformOutput',false);
        mean_E = cell2mat(mean_E);
        ciplot(nanmean(mean_E)+nansem(mean_E),nanmean(mean_E)-nansem(mean_E),tbBlocks,'r',0.6)
        Y1 = smooth2a(mean_C,0,10);Y1(:,isnan(nansum(Y1)))=[];
        Y2 = smooth2a(mean_E,0,10);Y2(:,isnan(nansum(Y2)))=[];
        [sig,~] = permtest2vec(Y1',Y2',100,0.05);
        a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
        plot(tbBlocks,a*1.8,'color','k','LineWidth',5)
    end
    axis([0 20 0 20])

end
%%
figure('name','inter-unit distance: members to other assembly members')
for iArea = 1:3
    subplot(1,3,iArea); hold on
    plot([5 5],[0.8 1.6],'g','LineWidth',1.5)
    plot([15 15],[0.8 1.6],'r','LineWidth',1.5)
    title([Areas{iArea}])
    if iArea<3
        mean_C = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2O{iArea}(~isempty_cell(PDistCont_.U2U_M2O{iArea}))','UniformOutput',false);
        mean_C = cell2mat(mean_C);
        ciplot(nanmean(mean_C)+nansem(mean_C),nanmean(mean_C)-nansem(mean_C),tbBlocks,0*[1 1 1],1)
        
        mean_E = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2OError{iArea}(~isempty_cell(PDistCont_.U2U_M2OError{iArea}))','UniformOutput',false);
        mean_E = cell2mat(mean_E);
        ciplot(nanmean(mean_E)+nansem(mean_E),nanmean(mean_E)-nansem(mean_E),tbBlocks,'r',0.6)
        Y1 = smooth2a(mean_C,0,10);Y1(:,isnan(nansum(Y1)))=[];
        Y2 = smooth2a(mean_E,0,10);Y2(:,isnan(nansum(Y2)))=[];
        [sig,p] = permtest2vec(Y1',Y2',100,0.05);
        a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
        plot(tbBlocks,a*1.8,'color','k','LineWidth',5)        
    else
        mean_C = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2O{iArea}{1}(~isempty_cell(PDistCont_.U2U_M2O{iArea}{1}))','UniformOutput',false);
        mean_C = cell2mat(mean_C);
        ciplot(nanmean(mean_C)+nansem(mean_C),nanmean(mean_C)-nansem(mean_C),tbBlocks,0*[1 1 1],1)
        
        mean_E = cellfun(@(x) mean(x,1),PDistCont_.U2U_M2OError{iArea}{1}(~isempty_cell(PDistCont_.U2U_M2OError{iArea}{1}))','UniformOutput',false);
        mean_E = cell2mat(mean_E);
        ciplot(nanmean(mean_E)+nansem(mean_E),nanmean(mean_E)-nansem(mean_E),tbBlocks,'r',0.6)
        Y1 = smooth2a(mean_C,0,10);Y1(:,isnan(nansum(Y1)))=[];
        Y2 = smooth2a(mean_E,0,10);Y2(:,isnan(nansum(Y2)))=[];
        [sig,~] = permtest2vec(Y1',Y2',100,0.05);
        a = nan(size(sig));a(sig)=1;% a(YZ)=NaN;
        plot(tbBlocks,a*1.8,'color','k','LineWidth',5)
    end
    axis([0 20 0 20])
end
%% Clustering - hiearchical
h_.distance_ = 'euclidean'
h_.method_   = 'ward'

for iArea = 1:2
    Membership_ = D.membershipsigcollapsed{iArea};
    input = D.TScollapsed{iArea};    input(isnan(input))=0;
    idx = [];
    %     idx = find(nansum(D.TSsigcollapsed{iArea}).*bw<=1);
    Membership_(idx) = [];
    input(:,idx)=[];
    
    input = zscore(input)';
    
    
    
    h_.eucD = pdist(input,h_.distance_);
    h_.clustTreeEuc = linkage(h_.eucD,h_.method_);
    
    h_.diagnostic.sqEucD=squareform(h_.eucD);
    
    [h_.diagnostic.coph, h_.diagnostic.sqCoph]=cophenet(h_.clustTreeEuc,h_.eucD);
    h_.diagnostic.sqCoph=squareform(h_.diagnostic.sqCoph);
    
    h_.clustTreeEuc(:,3)=h_.clustTreeEuc(:,3)./(max(h_.clustTreeEuc(:,3))); % link/link_max
    
    % figure('name','Decoding shape clustering','Color','w'); hold on
    [h,nodes] = dendrogram(h_.clustTreeEuc,0,'colorthreshold',0.5);
    set(gca,'TickDir','out',...
        'TickLength',[.002 0],...
        'XTickLabel',[]);
    ylabel('linkage / max. linkage','interpreter','none');
    set(gcf,'name','Histogram shape clustering','Color','w');
    view (-90,90)
    axis off
    set(gcf,'color','w')
end
% hierarchical plot
% T = clusterdata(clustTreeEuc,...
%                 'maxclust',noClusts,...
%                 'distance',distance_,...
%                 'linkage',method_); % agglomerative
noClusts=3
h_.T = cluster(h_.clustTreeEuc,noClusts); %



% figure; silhouette(cluster_input,h_.T,h_.distance_)

%%% remap data by clusters peak time (list)
cluster_input_=[];
% list = [4 3 2 5 1];
for clID=1:noClusts
    temp = find(h_.T==(clID));
    cluster_input_=[cluster_input_ ;input(temp,:)];
end
h_.ReorderedInput = cluster_input_;

figure('color','w')
for clID=1:noClusts
    clID
    subplot(ceil(noClusts^0.5),ceil(noClusts^0.5),clID); hold on
    temp = find(h_.T==(clID));
    %     plot(temp_x,cluster_input(temp,:),'color',[0.5 0.5 0.5]);
    imagesc(input(temp,:))
    plot(5*(1+mean(input(temp,:))),'k','LineWidth',2);
    plot([0 0],[0 1],':k')
    % 	axis([0 30 -1 10])
    axis off
end



% show pairwise similarity plot
m = mat2gray(squareform(pdist(cluster_input_,h_.distance_)));



% eucD = pdist(cluster_input_,h_.distance_);
% clustTreeEuc_ = linkage(h_.eucD,h_.method_);
% [~, m]=cophenet(clustTreeEuc_,eucD_);
% m=squareform(m);



ml=triu(logical(ones(size(m))),1);
m(~ml)=nan;
load('BlueWhite.mat')

figure('color','w');
h=pcolor(m);
set(h,'edgealpha',0)
colormap(cmap)
axis off;
axis ij;
% box on


s=3
input = zscore(Group.UnitFeatures{s}',[],2); % All units from all rats in area s

clob=clustergram(input,...
    'RowPDist',h_.distance_,...
    'Cluster','column',...
    'Linkage',h_.method_,...
    'Dendrogram',1000,... %100
    'DisplayRatio',0.2,...
    'Colormap',parula,...
    'Symmetric',true,...
    'OptimalLeafOrder',true)

% group_markers = struct('GroupNumber', {369,363,287,393,  391,},...
%                       'Annotation', ClassNames,...
%                       'Color', ClassColors);
%
%                   clob.RowGroupMarker = group_markers;

clear h nodes s cluster_input iRat temp ml m group_markers cmap cluster_input_ clID clob

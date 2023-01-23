%% %%%%%% PREAMBLE %%%%%%
clear
PlotOnline = false;
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Target = 'LONG';

trialTimeoutSecs = 5;
bw=0.05;

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
DelayRange = 2:3;

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
            load(sprintf('%s%s_%s_iFR50_Task.mat',pat2,fname,Areas{iArea}),'units','nu');
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
    %% Assembly power spectra
    for iOutcome = 1:2
        % get trials
        LeftTrials=[];RightTrials=[];
        SlatencyL = []; NPlatencyL = []; ClatencyL = [];
        SlatencyR = []; NPlatencyR = []; ClatencyR = [];
        for iDelay = DelayRange
            switch iOutcome
                case 1
                    eval(sprintf('LeftTrials =  [LeftTrials;t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('RightTrials = [RightTrials;t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}));
                    
                    eval(sprintf('SlatencyL =  [SlatencyL;t.%s.SamplePress_LeftCorrect*1e-6-t.%s.CueLight_LeftCorrect*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('SlatencyR =  [SlatencyR;t.%s.SamplePress_RightCorrect*1e-6-t.%s.CueLight_RightCorrect*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    
                    eval(sprintf('NPlatencyL =  [NPlatencyL;t.%s.NosePoke_LeftCorrect*1e-6-t.%s.DelayEnd_LeftCorrect*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('NPlatencyR =  [NPlatencyR;t.%s.NosePoke_RightCorrect*1e-6-t.%s.DelayEnd_RightCorrect*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    
                    eval(sprintf('ClatencyL =  [ClatencyL;t.%s.ChoicePress_LeftCorrect*1e-6-t.%s.NosePoke_LeftCorrect*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('ClatencyR =  [ClatencyR;t.%s.ChoicePress_RightCorrect*1e-6-t.%s.NosePoke_RightCorrect*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    
                case 2
                    eval(sprintf('LeftTrials = [LeftTrials;t.%s.SamplePress_LeftError'',t.%s.ChoicePress_LeftError''];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('RightTrials = [RightTrials;t.%s.SamplePress_RightError'',t.%s.ChoicePress_RightError''];',Delays_{iDelay},Delays_{iDelay}));
                    
                    eval(sprintf('SlatencyL =  [SlatencyL;t.%s.SamplePress_LeftError''*1e-6-t.%s.CueLight_LeftError''*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('SlatencyR =  [SlatencyR;t.%s.SamplePress_RightError''*1e-6-t.%s.CueLight_RightError''*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    
                    eval(sprintf('NPlatencyL =  [NPlatencyL;t.%s.NosePoke_LeftError''*1e-6-t.%s.DelayEnd_LeftError''*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('NPlatencyR =  [NPlatencyR;t.%s.NosePoke_RightError''*1e-6-t.%s.DelayEnd_RightError''*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    
                    eval(sprintf('ClatencyL =  [ClatencyL;t.%s.ChoicePress_LeftError''*1e-6-t.%s.NosePoke_LeftError''*1e-6];',Delays_{iDelay},Delays_{iDelay}));
                    eval(sprintf('ClatencyR =  [ClatencyR;t.%s.ChoicePress_RightError''*1e-6-t.%s.NosePoke_RightError''*1e-6];',Delays_{iDelay},Delays_{iDelay}));
            end
        end
        idx = sum(isnan(LeftTrials),2)>0 | sum([SlatencyL,NPlatencyL,ClatencyL]>trialTimeoutSecs,2)>0;
        LeftTrials(idx,:)=[]; SlatencyL(idx,:)=[]; NPlatencyL(idx,:)=[]; ClatencyL(idx,:)=[];
        idx = sum(isnan(RightTrials),2)>0 | sum([SlatencyR,NPlatencyR,ClatencyR]>trialTimeoutSecs,2)>0;
        RightTrials(idx,:)=[]; SlatencyR(idx,:)=[]; NPlatencyR(idx,:)=[]; ClatencyR(idx,:)=[];
        clear idx

        % Get Assem cutouts
        FSC_ = SelTimesAssem2([5],Ass.Tmtx,Ass.FSC,[LeftTrials;RightTrials]/1e6,bw);
        nTrials = size([LeftTrials;RightTrials],1);
        
        
        for iArea=1:3
            if ~isempty(FSC_{iArea})
                nAss = length(FSC_{iArea});
                for iAss=1:nAss
                    clear A period f P_ P__
                    for  iTrial =1:nTrials
                        LTr =size(FSC_{iArea}{iAss},1);
                        S1 = (FSC_{iArea}{iAss}(1:LTr/2,iTrial));
                        S2 = (FSC_{iArea}{iAss}(LTr/2+1:LTr,iTrial));
                        %% Wavelet spectrogram
                        %                         [A_1,f] = cwt(S1,'morse',1/bw,'FrequencyLimits',[1 10],'VoicesPerOctave',24,'WaveletParameters',[3,120]);
                        %                         [A_2,f] = cwt(S2,'morse',1/bw,'FrequencyLimits',[1 10],'VoicesPerOctave',24,'WaveletParameters',[3,120]);
                        %                         A = [A_1,A_2];
                        %                         P__ = [abs(A_1),abs(A_2)].^2;
                        %                         P__ = [10*log10(abs(A_1)),10*log10(abs(A_2))];
                        S1 = zscore(S1);
                        S2 = zscore(S2);
                        S2 = S2-(S2(1)-S1(end));
                        S_ = zscore([S1;S2]);
                        S_ = detrend(S_);
                        %                         plot(S_)
                        [A,f] = cwt(S_,'morse',1/bw,'FrequencyLimits',[1 10],'VoicesPerOctave',48,'WaveletParameters',[3,120]);
%                         [A,f] = cwt(S_,'morse',milliseconds(bw*1000),'PeriodLimits',[2 200]*milliseconds(bw*1000),'VoicesPerOctave',48,'WaveletParameters',[3,120]);
%                         f = 1./seconds(f);
                        P__ = abs(A).^2;
%                         P__ = 10*log10(abs(A));
                        t_ = (0:LTr-1)*bw;
                        
                        P_(:,:,iTrial) = P__;
                        
                        % Online plot
                        if PlotOnline
                            figure;
                            subplot(2,1,1); hold on
                            %                             plot(t_,FSC_{iArea}{iAss}(:,iTrial),'color',0.6*[1 1 1])
                            plot(t_,S_,'color',0.6*[1 1 1])
                            axis tight
                            subplot(2,1,2)
                            %imagesc(t,f,nanmean(10*log10(MP.rEnergy{1}),3))
                            imagesc(t_,f,P__)
                            set(gca,'YDir','normal')
                            caxis([0 0.005])
                        end
                    end
                    if iOutcome == 1
                        P.Corr{iArea}{iFile}{iAss} = P_;
                        P.SlatencyCorr{iFile}   = [SlatencyL;SlatencyR];
                        P.NPlatencyCorr{iFile}  = [NPlatencyL;NPlatencyR];
                        P.ClatencyCorr{iFile}   = [ClatencyL;ClatencyR];
                    else
                        P.Err{iArea}{iFile}{iAss} = P_;
                        P.SlatencyErr{iFile}   = [SlatencyL;SlatencyR];
                        P.NPlatencyErr{iFile}  = [NPlatencyL;NPlatencyR];
                        P.ClatencyErr{iFile}   = [ClatencyL;ClatencyR];
                    end
                    P.f = f;
                    P.t = t_;
                end
            end
        end
    end
end
%% Collapse
% iArea =3;
% P_ = P.Corr{iArea}(~isempty_cell(P.Corr{iArea}));
% for iFile =1:length(P_)
%     for iAss =1:length(P_{iFile})
%         figure
%         imagesc(P.t,P.f,nanmean(P_{iFile}{iAss},3));
%         set(gca,'YDir','normal')
%         %     caxis([-40 -10])
%             caxis([0 0.01])
%     end
% end
figure;
for iArea =1:3
    P_mean5Hz{iArea}   = [];
    P_meanTheta{iArea} = [];
    
    P_ = P.Corr{iArea}(~isempty_cell(P.Corr{iArea}));
    
    Pave_ = [];
    for iFile =1:length(P_)
        for iAss =1:length(P_{iFile})
            Pave_ = cat(3,Pave_,nanmean(P_{iFile}{iAss},3));

        end
    end
    Pave.Corr{iArea} = Pave_;
    subplot(2,3,iArea); hold on
    imagesc(P.t,P.f,nanmean(Pave_,3));
    set(gca,'YDir','normal')
%         caxis([-40 -10])
%     caxis([0 0.01])
    patch([0 0 20 20],[4 5 5 4],'k','facealpha',0.2)
    patch([0 0 20 20],[6 8 8 6],'k','facealpha',0.2)
    axis([0 20 1 10])
    plot([5 5],[0 20],'g')
    plot([15 15],[0 20],'r')
    plot([10 10],[0 20],':k')
    
    P_ = P.Err{iArea}(~isempty_cell(P.Err{iArea}));
    
    Pave_ = [];
    for iFile =1:length(P_)
        for iAss =1:length(P_{iFile})
            Pave_ = cat(3,Pave_,nanmean(P_{iFile}{iAss},3));
        end
    end
    Pave.Err{iArea} = Pave_;
    subplot(2,3,iArea+3);hold on
    imagesc(P.t,P.f,nanmean(Pave_,3));
    set(gca,'YDir','normal')
%         caxis([-40 -10])
%     caxis([0 0.01])
    patch([0 0 20 20],[4 5 5 4],'k','facealpha',0.2)
    patch([0 0 20 20],[6 8 8 6],'k','facealpha',0.2)
    axis([0 20 1 10])
    plot([5 5],[0 20],'g')
    plot([15 15],[0 20],'r')
    plot([10 10],[0 20],':k')
end
colormap (jet)

%% Difference spectra
figure;
for iArea=1:3
%     y_ = ((Pave.Err{iArea})-Pave.Corr{iArea})./(Pave.Err{iArea}+Pave.Corr{iArea});
%     y_ = Pave.Err{iArea}./Pave.Corr{iArea};
    y_ = Pave.Err{iArea}-Pave.Corr{iArea};
    y_ = squeeze(nanmean(y_,3));
    subplot(1,3,iArea)
    imagesc(P.t,P.f,y_);
    set(gca,'YDir','normal')
%         caxis([-1 1])
%     caxis([0.1 1.9])
    caxis([-0.01 0.01])
end
% colormap([1 0 0; 1 1 1; 0 0 1])
load blue_white_red.mat
colormap(cmap)

%% Spectral power
fRange_5Hz = (P.f>=4 & P.f<5);
fRange_Theta = (P.f>=6 & P.f<8);
Tblank = false(size(P.t));
Tblank([P.t>=0 & P.t<1.5] | ...
    [P.t>=8.5 & P.t<11.5] | ...
    [P.t>=18.5 & P.t<21]) = true;
blankWidth = 2.5;

% plot 5Hz
figure('color','w','name','5Hz power')
for iArea =1:3
%     P_Corr = squeeze(nanmean(Pave.Corr{iArea}(fRange_5Hz,:,:),1));
%     P_Err = squeeze(nanmean(Pave.Err{iArea}(fRange_5Hz,:,:),1));
    P_Corr = squeeze(nanmax(Pave.Corr{iArea}(fRange_5Hz,:,:),1));
    P_Err = squeeze(nanmax(Pave.Err{iArea}(fRange_5Hz,:,:),1));
    %P_Corr(Tblank,:)=NaN; P_Err(Tblank,:)=NaN;
    subplot(1,3,iArea); hold on

%     plot(P.t,P_Corr,'k')
%     plot(P.t,P_Err,'r')
    ciplot(nanmean(P_Corr,2)+nansem(P_Corr,2),nanmean(P_Corr,2)-nansem(P_Corr,2),P.t,'k')
    ciplot(nanmean(P_Err,2)+nansem(P_Err,2),nanmean(P_Err,2)-nansem(P_Err,2),P.t,'r')
    if iArea ==2
        plot([14 19], [0.002 0.002],'k','LineWidth',1.5)
        plot([14 14], [0.002 0.003],'k','LineWidth',1.5)
    end
        axis([min(P.t) 20 0  0.01])
    axis off
end

% plot Theta
figure('color','w','name','Theta power')
for iArea =1:3
%     P_Corr = squeeze(nanmean(Pave.Corr{iArea}(fRange_Theta,:,:),1));
%     P_Err = squeeze(nanmean(Pave.Err{iArea}(fRange_Theta,:,:),1));    
    P_Corr = squeeze(nanmax(Pave.Corr{iArea}(fRange_Theta,:,:),1));
    P_Err = squeeze(nanmax(Pave.Err{iArea}(fRange_Theta,:,:),1));    
    
    subplot(1,3,iArea); hold on
    %     plot([5 5],[-100 10],'g','LineWidth',1.5)
    %     plot([15 15],[-100 10],'r','LineWidth',1.5)
%         plot(P.t,P_Corr,'k')
%         plot(P.t,P_Err,'r')
    ciplot(nanmean(P_Corr,2)+nansem(P_Corr,2),nanmean(P_Corr,2)-nansem(P_Corr,2),P.t,'k')
    ciplot(nanmean(P_Err,2)+nansem(P_Err,2),nanmean(P_Err,2)-nansem(P_Err,2),P.t,'r')
    
    if iArea ==2
        plot([14 19], [0.002 0.002],'k','LineWidth',1.5)
        plot([14 14], [0.002 0.003],'k','LineWidth',1.5)
    end

        axis([min(P.t) 20 0  0.01])
    axis off
end
%% Spectral differences
fRange_5Hz = (P.f>=3.5 & P.f<6.5);
fRange_Theta = (P.f>=7.5 & P.f<12);
figure
for iArea =1:3
    P_Corr = squeeze(nanmean(Pave.Corr{iArea}(fRange_5Hz,:,:),1));
    P_Err = squeeze(nanmean(Pave.Err{iArea}(fRange_5Hz,:,:),1));
    subplot(1,3,iArea); hold on
    %     plot(P.t,P_Corr,'k')
    %     plot(P.t,P_Err,'r')
    P_Diff = (P_Err - P_Corr)./(P_Err + P_Corr);
    %     P_Diff = P_Err./P_Corr;
    %     plot(P.t,P_Diff,'g')
    ciplot(nanmean(P_Diff,2)+nansem(P_Diff,2),nanmean(P_Diff,2)-nansem(P_Diff,2),P.t,'g')
    
    P_Corr = squeeze(nanmean(Pave.Corr{iArea}(fRange_Theta,:,:),1));
    P_Err = squeeze(nanmean(Pave.Err{iArea}(fRange_Theta,:,:),1));
    
    %     plot(P.t,P_Corr,'k')
    %     plot(P.t,P_Err,'r')
    P_Diff = (P_Err - P_Corr)./(P_Err + P_Corr);
    %     P_Diff = P_Err./P_Corr;
    
    %     plot(P.t,P_Diff,'g')
    ciplot(nanmean(P_Diff,2)+nansem(P_Diff,2),nanmean(P_Diff,2)-nansem(P_Diff,2),P.t,[1 0.5 0.2])
    
    rectangle('Position',[0-blankWidth/2 -5 blankWidth  50],'LineStyle','none','FaceColor','w')
    rectangle('Position',[10-blankWidth/2 -5 blankWidth  50],'LineStyle','none','FaceColor','w')
    rectangle('Position',[20-blankWidth/2 -5 blankWidth  50],'LineStyle','none','FaceColor','w')
    % axis([min(P.t) 20 0.8 1.2])
    axis([min(P.t) 20 -1 1])
    % axis([min(P.t) 20 -1 1])
end
%% correlation between behavioural latencies and band power - correct trials 
% fRange_ = fRange_5Hz;% fRange_5Hz;%fRange_Theta
fRange_ = fRange_Theta;% fRange_5Hz;%fRange_Theta
tRange_Sample = (P.t>=0 &P.t<5);
tRange_Choice = (P.t>10 &P.t<15);
% tRange_Sample = (P.t>=2.5 & P.t<4);
% tRange_Choice = (P.t>12.5 & P.t<14);
% Sample/correct    Choice/correct     Sample/error        Choice/error
Smdl_.pFitC  = [];  Cmdl_.pFitC  = []; Smdl_.pFitE  = [];  Cmdl_.pFitE  = [];
Smdl_.RsqC   = [];  Cmdl_.RsqC   = []; Smdl_.RsqE   = [];  Cmdl_.RsqE   = [];
Smdl_.slopeC = [];  Cmdl_.slopeC = []; Smdl_.slopeE = [];  Cmdl_.slopeE = [];

            
for iArea =3%1:3
    files = find(~isempty_cell(P.Corr{iArea}));
    tempList = P.Corr{iArea}(files);
    figure('name','Correct trials');hold on
    % Scatter points first 
    for iFile = 1:length(files)
        Slatency = P.SlatencyCorr{files(iFile)};
        Clatency = P.ClatencyCorr{files(iFile)};
        NPlatency = P.NPlatencyCorr{files(iFile)};
        idx = sum([Slatency ,Clatency ]>trialTimeoutSecs,2)>0;
        Slatency(idx)=[];NPlatency(idx)=[];Clatency(idx)=[];
        for iAss = 1:length(tempList{iFile})
            %colorAss = rand(1,3);
            colorAss = [1 1 1];
            tempP = tempList{iFile}{iAss};
            
            subplot(1,2,1);hold on
            powCutoutS = nanmean(squeeze(nanmean(tempP(fRange_,tRange_Sample,:),1)),1);
            powCutoutS = log10(powCutoutS);
            powCutoutS(idx)=[];
            scatter(powCutoutS,Slatency,5,0.6*colorAss,'filled');
            
            subplot(1,2,2);hold on
            
            powCutoutC = nanmean(squeeze(nanmean(tempP(fRange_,tRange_Choice,:),1)),1);
            powCutoutC = log10(powCutoutC);
            powCutoutC(idx)=[];
            scatter(powCutoutC,Clatency,5,0.6*colorAss,'filled');
            
        end
    end
    % Plot regressions on top
    for iFile = 1:length(files)
        Slatency = P.SlatencyCorr{files(iFile)};
        Clatency = P.ClatencyCorr{files(iFile)};
        NPlatency = P.NPlatencyCorr{files(iFile)};
        idx = sum([Slatency ,Clatency ]>trialTimeoutSecs,2)>0;
        Slatency(idx)=[];NPlatency(idx)=[];Clatency(idx)=[];
        for iAss = 1:length(tempList{iFile})
            %colorAss = rand(1,3);
            colorAss = [0 0 0];
            tempP = tempList{iFile}{iAss};
            
            subplot(1,2,1);hold on
            powCutoutS = nanmean(squeeze(nanmean(tempP(fRange_,tRange_Sample,:),1)),1);
            powCutoutS = log10(powCutoutS);
            powCutoutS(idx)=[];
            mdl = fitlm(powCutoutS,Slatency);
%             Smdl_.pFitC   = [Smdl_.pFitC;mdl.Coefficients{2,4}];
            Smdl_.pFitC   = [Smdl_.pFitC;coefTest(mdl)];
            Smdl_.RsqC    = [Smdl_.RsqC;mdl.Rsquared.Adjusted];
%             Smdl_.RsqC    = [Smdl_.RsqC;mdl.Rsquared.Ordinary];
            Smdl_.slopeC  = [Smdl_.slopeC;mdl.Coefficients{2,1}];
            
            h = plot(mdl);
            if Smdl_.pFitC(end)<0.05
                h(1).Marker = 'none';
                h(2).LineWidth=1.5;
                h(2).LineStyle='-';
                h(2).Color = colorAss;
                h(3).LineStyle='none';
                h(4).LineStyle='none';
            else
                 h(1).Marker = 'none';
                h(1).Marker = 'none';
                h(2).LineWidth=1.5;
                h(2).LineStyle='-';
                h(2).Color = colorAss+0.6;
                h(3).LineStyle='none';
                h(4).LineStyle='none';
            end
            legend off
            
            subplot(1,2,2);hold on
            
            powCutoutC = nanmean(squeeze(nanmean(tempP(fRange_,tRange_Choice,:),1)),1);
            powCutoutC = log10(powCutoutC);
            powCutoutC(idx)=[];
            mdl = fitlm(powCutoutC,Clatency);
%             Cmdl_.pFitC   = [Cmdl_.pFitC;mdl.Coefficients{2,4}];
            Cmdl_.pFitC   = [Cmdl_.pFitC;coefTest(mdl)];
            Cmdl_.RsqC    = [Cmdl_.RsqC;mdl.Rsquared.Adjusted];
%             Cmdl_.RsqC    = [Cmdl_.RsqC;mdl.Rsquared.Ordinary];
            Cmdl_.slopeC  = [Cmdl_.slopeC;mdl.Coefficients{2,1}];
            
            h = plot(mdl);
            if Cmdl_.pFitC(end)<0.05
                h(1).Marker = 'none';
                h(2).LineWidth=1.5;
                h(2).LineStyle='-';
                h(2).Color =colorAss;
                h(3).LineStyle='none';
                h(4).LineStyle='none';
            else
             h(1).Marker = 'none';
                h(1).Marker = 'none';
                h(2).LineWidth=1.5;
                h(2).LineStyle='-';
                h(2).Color = colorAss+0.6;
                h(3).LineStyle='none';
                h(4).LineStyle='none';
            end
            legend off
        end
    end
end
subplot(1,2,1)
title('Sample')
xlabel('Log(Mean band power)')
ylabel('Sample latency (s)')
ylim([0 5])

subplot(1,2,2)
title('Choice')
xlabel('Log(Mean band power)')
ylabel('Choice latency (s)')
ylim([0 5])
%% correlation between behavioural latencies and band power - error trials 
for iArea =3%1:3
    files = find(~isempty_cell(P.Err{iArea}));
    tempList = P.Err{iArea}(files);
    figure('name','Error trials');hold on
    % Scatter points first 
    for iFile = 1:length(files)
        Slatency = P.SlatencyErr{files(iFile)};
        Clatency = P.ClatencyErr{files(iFile)};
        NPlatency = P.NPlatencyErr{files(iFile)};
        idx = sum([Slatency ,Clatency ]>trialTimeoutSecs,2)>0;
        Slatency(idx)=[];NPlatency(idx)=[];Clatency(idx)=[];
        for iAss = 1:length(tempList{iFile})
            %colorAss = rand(1,3);
            colorAss = [1 1 1];
            tempP = tempList{iFile}{iAss};
            
            subplot(1,2,1);hold on
            powCutoutS = nanmean(squeeze(nanmean(tempP(fRange_,tRange_Sample,:),1)),1);
            powCutoutS = log10(powCutoutS);
            powCutoutS(idx)=[];
            scatter(powCutoutS,Slatency,5,0.6*colorAss,'filled');
            
            subplot(1,2,2);hold on
            
            powCutoutC = nanmean(squeeze(nanmean(tempP(fRange_,tRange_Choice,:),1)),1);
            powCutoutC = log10(powCutoutC);
            powCutoutC(idx)=[];
            scatter(powCutoutC,Clatency,5,0.6*colorAss,'filled');
            
        end
    end
    
    % Plot regressions on top
    for iFile = 1:length(files)
        Slatency = P.SlatencyErr{files(iFile)};
        Clatency = P.ClatencyErr{files(iFile)};
        NPlatency = P.NPlatencyErr{files(iFile)};
        idx = sum([Slatency ,Clatency ]>trialTimeoutSecs,2)>0;
        Slatency(idx)=[];NPlatency(idx)=[];Clatency(idx)=[];
        for iAss = 1:length(tempList{iFile})
            %colorAss = rand(1,3);
            colorAss = [0 0 0];
            tempP = tempList{iFile}{iAss};
            
            subplot(1,2,1);hold on
            powCutoutS = nanmean(squeeze(nanmean(tempP(fRange_,tRange_Sample,:),1)),1);
            powCutoutS = log10(powCutoutS);
            powCutoutS(idx) = [];
            mdl = fitlm(powCutoutS,Slatency);
%             Smdl_.pFitE   = [Smdl_.pFitE;mdl.Coefficients{2,4}];
            try    
                Smdl_.pFitE   = [Smdl_.pFitE;coefTest(mdl)];
            catch
                Smdl_.pFitE   = [Smdl_.pFitE;1];
            end
            Smdl_.RsqE    = [Smdl_.RsqE;mdl.Rsquared.Adjusted];
%             Smdl_.RsqE    = [Smdl_.RsqE;mdl.Rsquared.Ordinary];
            Smdl_.slopeE  = [Smdl_.slopeE;mdl.Coefficients{2,1}];
            
            h = plot(mdl);
            if Smdl_.pFitE(end)<0.05
                h(1).Marker = 'none';
                h(2).LineWidth=1.5;
                h(2).LineStyle='-';
                h(2).Color = colorAss;
                h(3).LineStyle='none';
                h(4).LineStyle='none';
            else
                 h(1).Marker = 'none';
                h(1).Marker = 'none';
                h(2).LineWidth=1.5;
                h(2).LineStyle='-';
                h(2).Color = colorAss+0.6;
                h(3).LineStyle='none';
                h(4).LineStyle='none';
            end
            legend off
            
            subplot(1,2,2);hold on
            
            powCutoutC = nanmean(squeeze(nanmean(tempP(fRange_,tRange_Choice,:),1)),1);
            powCutoutC = log10(powCutoutC);
            powCutoutC(idx) = [];
            mdl = fitlm(powCutoutC,Clatency);
%             Cmdl_.pFitE   = [Cmdl_.pFitE;mdl.Coefficients{2,4}];
            try
                Cmdl_.pFitE   = [Cmdl_.pFitE;coefTest(mdl)];
            catch
                Cmdl_.pFitE   = [Cmdl_.pFitE;1];
            end
            Cmdl_.RsqE    = [Cmdl_.RsqE;mdl.Rsquared.Adjusted];
%             Cmdl_.RsqE    = [Cmdl_.RsqE;mdl.Rsquared.Ordinary];
            Cmdl_.slopeE  = [Cmdl_.slopeE;mdl.Coefficients{2,1}];
            
            h = plot(mdl);
            if Cmdl_.pFitE(end)<0.05
                h(1).Marker = 'none';
                h(2).LineWidth=1.5;
                h(2).LineStyle='-';
                h(2).Color =colorAss;
                h(3).LineStyle='none';
                h(4).LineStyle='none';
            else
             h(1).Marker = 'none';
                h(1).Marker = 'none';
                h(2).LineWidth=1.5;
                h(2).LineStyle='-';
                h(2).Color = colorAss+0.6;
                h(3).LineStyle='none';
                h(4).LineStyle='none';
            end
            legend off
        end
    end       
end
subplot(1,2,1)
title('Sample')
xlabel('Band power')
ylabel('Sample latency (s)')
ylim([0 5])

subplot(1,2,2)
title('Choice')
xlabel('Band power')
ylabel('Choice latency (s)')
ylim([0 5])
clear tempP tempList files
%% plot distributions of model fits - slope

bins = -2.01:0.1:2.01;
figure;
subplot(1,2,1); hold on
temp  = Smdl_.slopeC;
tempP = Smdl_.pFitC;
HC = histc(temp(tempP<0.05),bins);HC=HC./numel(temp);
stairs(bins,HC,'k','LineWidth',1.5)
HC = histc(temp(tempP>0.05),bins);HC=HC./numel(temp); 
stairs(bins,HC,':k','LineWidth',1)
temp  = Smdl_.slopeE;
tempP = Smdl_.pFitE;
HC = histc(temp(tempP<0.05),bins);HC=HC./numel(temp);
stairs(bins,HC,'r','LineWidth',1.5)
HC = histc(temp(tempP>0.05),bins);HC=HC./numel(temp); 
stairs(bins,HC,':r','LineWidth',1)
title('Samples')

subplot(1,2,2); hold on
temp  = Cmdl_.slopeC;
tempP = Cmdl_.pFitC;
HC = histc(temp(tempP<0.05),bins);HC=HC./numel(temp);
stairs(bins,HC,'k','LineWidth',1.5)
HC = histc(temp(tempP>0.05),bins);HC=HC./numel(temp); 
stairs(bins,HC,':k','LineWidth',1)
temp  = Cmdl_.slopeE;
tempP = Cmdl_.pFitE;
HC = histc(temp(tempP<0.05),bins);HC=HC./numel(temp);
stairs(bins,HC,'r','LineWidth',1.5)
HC = histc(temp(tempP>0.05),bins);HC=HC./numel(temp); 
stairs(bins,HC,':r','LineWidth',1)
title('Choices')

%% plot distributions of model fits - R^2

bins = -1.01:0.05:1.01;
figure;
subplot(1,2,1); hold on
temp  = Smdl_.RsqC;
tempP = Smdl_.pFitC ;
HC = histc(temp(tempP<0.05),bins);HC=HC./numel(temp);
stairs(bins,HC,'k','LineWidth',1.5)
HC = histc(temp(tempP>0.05),bins);HC=HC./numel(temp); 
stairs(bins,HC,':k','LineWidth',1)
temp  = Smdl_.RsqE;
tempP = Smdl_.pFitE;
HC = histc(temp(tempP<0.05),bins);HC=HC./numel(temp);
stairs(bins,HC,'r','LineWidth',1.5)
HC = histc(temp(tempP>0.05),bins);HC=HC./numel(temp); 
stairs(bins,HC,':r','LineWidth',1)
title('Samples')

subplot(1,2,2); hold on
temp  = Cmdl_.RsqC;
tempP = Cmdl_.pFitC;
HC = histc(temp(tempP<0.05),bins);HC=HC./numel(temp);
stairs(bins,HC,'k','LineWidth',1.5)
HC = histc(temp(tempP>0.05),bins);HC=HC./numel(temp); 
stairs(bins,HC,':k','LineWidth',1)
temp  = Cmdl_.RsqE;
tempP = Cmdl_.pFitE;
HC = histc(temp(tempP<0.05),bins);HC=HC./numel(temp);
stairs(bins,HC,'r','LineWidth',1.5)
HC = histc(temp(tempP>0.05),bins);HC=HC./numel(temp); 
stairs(bins,HC,':r','LineWidth',1)
title('Choices')




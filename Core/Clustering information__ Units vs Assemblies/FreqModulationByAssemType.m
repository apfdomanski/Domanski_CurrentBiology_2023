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
bw = 0.05;
tlimsAll = [-5 5];
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
    %% Limit spike times by behavioural epoch
    STranges=[];
    tRanges=[]; % Cue-Sample
    for iDelay = DelayRange
        tRanges = [tRanges;eval(sprintf('[t.%s.CueLight_LeftCorrect,t.%s.SamplePress_LeftCorrect;t.%s.CueLight_RightCorrect,t.%s.SamplePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))];        
    end    
    tRanges=tRanges*1e-6;
    [~,STranges.PFC.SampleCorr] = restrictTranges(PFCcells,tRanges);
    [~,STranges.HP.SampleCorr]  = restrictTranges(HPcells,tRanges);
    
    tRanges=[]; % Cue-Sample (Error)
    for iDelay = DelayRange
        tRanges = [tRanges;eval(sprintf('[t.%s.CueLight_LeftError'',t.%s.SamplePress_LeftError'';t.%s.CueLight_RightError'',t.%s.SamplePress_RightError''];',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))];        
    end    
    tRanges=tRanges*1e-6;
    [~,STranges.PFC.SampleErr] = restrictTranges(PFCcells,tRanges);
    [~,STranges.HP.SampleErr]  = restrictTranges(HPcells,tRanges);
    
    tRanges=[]; % Delay
    for iDelay = DelayRange
        tRanges = [tRanges;eval(sprintf('[t.%s.SamplePress_LeftCorrect,t.%s.DelayEnd_LeftCorrect;t.%s.SamplePress_RightCorrect,t.%s.DelayEnd_RightCorrect];',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))];        
    end    
    tRanges=tRanges*1e-6;
    [~,STranges.PFC.DelayCorr] = restrictTranges(PFCcells,tRanges);
    [~,STranges.HP.DelayCorr] = restrictTranges(HPcells,tRanges);
    
    tRanges=[]; % Delay (Error)
    for iDelay = DelayRange
        tRanges = [tRanges;eval(sprintf('[t.%s.SamplePress_LeftError'',t.%s.DelayEnd_LeftError'';t.%s.SamplePress_RightError'',t.%s.DelayEnd_RightError''];',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))];        
    end    
    tRanges=tRanges*1e-6;
    [~,STranges.PFC.DelayErr] = restrictTranges(PFCcells,tRanges);
    [~,STranges.HP.DelayErr] = restrictTranges(HPcells,tRanges);
    
    tRanges=[]; % Nosepoke
    for iDelay = DelayRange
        tRanges = [tRanges;eval(sprintf('[t.%s.DelayEnd_LeftCorrect,t.%s.NosePoke_LeftCorrect;t.%s.DelayEnd_RightCorrect,t.%s.NosePoke_RightCorrect];',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))];        
    end    
    tRanges=tRanges*1e-6;
    [~,STranges.PFC.NosePokeCorr] = restrictTranges(PFCcells,tRanges);
    [~,STranges.HP.NosePokeCorr] = restrictTranges(HPcells,tRanges);
    
    tRanges=[]; % Nosepoke (Error)
    for iDelay = DelayRange
        tRanges = [tRanges;eval(sprintf('[t.%s.DelayEnd_LeftError'',t.%s.NosePoke_LeftError'';t.%s.DelayEnd_RightError'',t.%s.NosePoke_RightError''];',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))];        
    end    
    tRanges=tRanges*1e-6;
    [~,STranges.PFC.NosePokeErr] = restrictTranges(PFCcells,tRanges);
    [~,STranges.HP.NosePokeErr] = restrictTranges(HPcells,tRanges);    
    
    tRanges=[]; % Choice
    for iDelay = DelayRange
        tRanges = [tRanges;eval(sprintf('[t.%s.NosePoke_LeftCorrect,t.%s.ChoicePress_LeftCorrect;t.%s.NosePoke_RightCorrect,t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))];        
    end    
    tRanges=tRanges*1e-6;
    [~,STranges.PFC.ChoiceCorr] = restrictTranges(PFCcells,tRanges);
    [~,STranges.HP.ChoiceCorr] = restrictTranges(HPcells,tRanges);
    
    tRanges=[]; % Choice (Error)
    for iDelay = DelayRange
        tRanges = [tRanges;eval(sprintf('[t.%s.NosePoke_LeftError'',t.%s.ChoicePress_LeftError'';t.%s.NosePoke_RightError'',t.%s.ChoicePress_RightError''];',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))];        
    end    
    tRanges=tRanges*1e-6;
    [~,STranges.PFC.ChoiceErr] = restrictTranges(PFCcells,tRanges);
    [~,STranges.HP.ChoiceErr] = restrictTranges(HPcells,tRanges);    
    
    % Reward
    Reward_rangeCorr = [];
    Reward_rangeErr=[];
    Cue_range = [];
    for iDelay = DelayRange
        Reward_rangeCorr = [Reward_rangeCorr; eval(sprintf('[t.%s.ChoicePress_LeftCorrect;t.%s.ChoicePress_RightCorrect];',Delays_{iDelay},Delays_{iDelay}))];  
        Reward_rangeErr  = [Reward_rangeErr;  eval(sprintf('[t.%s.ChoicePress_LeftError'';t.%s.ChoicePress_RightError''];',Delays_{iDelay},Delays_{iDelay}))];  
        Cue_range        = [Cue_range;        eval(sprintf('[t.%s.CueLight_LeftCorrect;t.%s.CueLight_LeftError'';t.%s.CueLight_RightCorrect;t.%s.CueLight_RightError''];',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))];          
    end 
    Reward_rangeCorr = sort(Reward_rangeCorr*1e-6); Reward_rangeCorr(isnan(Reward_rangeCorr))=[];
    Reward_rangeErr  = sort(Reward_rangeErr*1e-6); Reward_rangeErr(isnan(Reward_rangeErr))=[];
%     Cue_range        = sort(Cue_range);

tRangeCorr = [Reward_rangeCorr,Reward_rangeCorr+5];
tRangeErr = [Reward_rangeErr,Reward_rangeErr+5];

[~,STranges.PFC.RewardCorr] = restrictTranges(PFCcells,tRangeCorr);
[~,STranges.HP.RewardCorr] = restrictTranges(HPcells,tRangeCorr);

[~,STranges.PFC.RewardErr] = restrictTranges(PFCcells,tRangeErr);
[~,STranges.HP.RewardErr] = restrictTranges(HPcells,tRangeErr);
%     %Find the next cue time after each choice is made
%     tRangeCorr = nan(length(Choice_rangeCorr),2); 
%     for i=1:length(Choice_rangeCorr)
%         t_ = Cue_range-Choice_rangeCorr(i);
%         tRangeCorr(i,:) = [Choice_rangeCorr(i),Choice_rangeCorr(i)+ min(t_(t_>0))];
%     end
%     
%     tRangeErr = nan(length(Choice_rangeErr),2); 
%     for i=1:length(Choice_rangeErr)
%         t_ = Cue_range-Choice_rangeErr(i);
%         tRangeErr(i,:) = [Choice_rangeErr(i),Choice_rangeErr(i)+ min(t_(t_>0))];
%     end
%     diff(tRangeCorr,[],2)*1e-6
%     diff(tRangeErr,[],2)*1e-6
%     figure; hold on
%     scatter(PFCcells{1}.t*1e-4,ones(size(PFCcells{1}.t)))

    %% Spike Autocorr - All spikes
%     clear AC tlims_
%     t_ = [collapseStructure(t.Short);...
%           collapseStructure(t.Medium);...
%           collapseStructure(t.Long)]*1e-6;
%     t_(isnan(t_))=[];
%     tlims_ = [min(t_), max(t_)]+[-5 5];
%     [AC{1},lags] = GetSpikeAutocorr(PFCcells,0.01,1,1e-4,tlims_);
%     [AC{2},lags] = GetSpikeAutocorr(HPcells,0.01,1,1e-4,tlims_);
%     [P{1},f] = GetSpectrumfromAutocorr(AC{1},lags);
%     [P{2},f] = GetSpectrumfromAutocorr(AC{2},lags);
%     
%     AC{3} = [AC{1};AC{2}];
%     P{3} = [P{1};P{2}];
%     %%
%     figure;hold on
%     for iArea=1:3
% %         subplot(1,3,iArea);hold on 
%         AC_ = AC{iArea}./sum(AC{iArea},2);
%         AC_(sum(isnan(AC_),2)>0,:)=[];
% %         AC_ = staggerplot(AC_,0,0.1);
%         lags_ = repmat(lags,size(AC_,1),1);
% %         stairs(lags_',AC_')
%         stairs(lags,mean(AC_)+0.01*iArea,'color',color_{iArea})
%     end
% %     set(gca,'XTick'
% smoothing_ = 0;
% figure; hold on
% for iArea=1:3   
%     x = P{iArea};
%     x = x - mean(x(:,1:2),2);
%     x = smooth2a(x,0,smoothing_);
%     x(isinf(x)) = nan;
%     plot(f,nanmean(x),'color',color_{iArea},'LineWidth',1.5)
% %     plot(f,x,'color',color_{iArea},'LineWidth',1.5)
% end
% xlabel('Freq (Hz)')
% ylabel('power (dB)')
% % axis([0.6 100 -30 20])
% % set(gca,'XScale','log')
    %% Spike frequency following - all spikes
    clear ModIndex powWelch tlims_
    t_ = [collapseStructure(t.Short);collapseStructure(t.Medium);collapseStructure(t.Long)]*1e-6;
    t_(isnan(t_))=[];
    tlims_ = [min(t_), max(t_)]+[-5 5];
    [ModIndex{1},f,~,powWelch{1},r0{1}] = CalculateFrequencyModulation2(PFCcells,1e-4,2,tlims_);
    [ModIndex{2},f,~,powWelch{2},r0{2}] = CalculateFrequencyModulation2(HPcells,1e-4,2,tlims_);
%     figure;
%     subplot(2,1,1);hold on
%     plot(f,ModIndex{1},'b')
%     plot(f,ModIndex{2},'r')
%     set(gca,'Xscale','log')
%     subplot(2,1,2);hold on
%     plot(f, powWelch{1},'b')
%     plot(f, powWelch{2},'r')
%     set(gca,'Xscale','log')

    % bins = 0:0.5:10;
    % figure; hold on
    % stairs(plot(bins,(histc(r0{1},bins)),'b'))
    % stairs(plot(bins,(histc(r0{2},bins)),'r'))
    FreqMod.ModIndex{iFile}{1} = ModIndex{1};
    FreqMod.ModIndex{iFile}{2} = ModIndex{2};
    FreqMod.f = f;
    clear ModIndex powWelch r0 t_ tlims_ f
    %% Assembly power spectra
%     for iArea =1:3
%         if ~isempty(Ass.FSC{iArea})
%             [~,f,powWelch]=getAssemblyPowerSpectrum(Ass.FSC{iArea},Ass.Tmtx,600);
%         end
%     end
    %% Spike frequency following - By event
% %     clear ModIndex powWelch
% %     events = fieldnames(STranges.PFC);
% %     for iEvent=1:length(events)
% %         ST_PFC = eval(sprintf('STranges.PFC.%s',events{iEvent}));
% %         ST_HP  = eval(sprintf('STranges.HP.%s',events{iEvent}));
% %         ModIndex=cell(2,1);
% %         powWelch=cell(2,1);
% %         r0=cell(2,1);
% %         [ModIndex{1},f,~,powWelch{1},r0{1}] = CalculateFrequencyModulation2(ST_PFC,1e-4,2);
% %         [ModIndex{2},f,~,powWelch{2},r0{2}] = CalculateFrequencyModulation2(ST_HP,1e-4,2);
% %         
% %         eval(sprintf('FreqModEvt.ModIndex.%s{iFile}{1} = ModIndex{1};',events{iEvent}))
% %         eval(sprintf('FreqModEvt.ModIndex.%s{iFile}{2} = ModIndex{2};',events{iEvent}))
% %         eval(sprintf('FreqModEvt.powWelch.%s{iFile}{1} = powWelch{1};',events{iEvent}))
% %         eval(sprintf('FreqModEvt.powWelch.%s{iFile}{2} = powWelch{2};',events{iEvent}))
% %         eval(sprintf('FreqModEvt.r0.%s{iFile}{1} = r0{1};',events{iEvent}))
% %         eval(sprintf('FreqModEvt.r0.%s{iFile}{2} = r0{2};',events{iEvent}))
% %         
% %         FreqModEvt.f = f;
% %         clear ModIndex powWelch r0  f
% %     end
    %% Assembly power spectra
%     for iArea =1:3
%         if ~isempty(Ass.FSC{iArea})
%             [~,f,powWelch]=getAssemblyPowerSpectrum(Ass.FSC{iArea},Ass.Tmtx,600);
%         end
%     end
    %% Spike-Spike coherence - all spikes
    clear ModIndex powWelch tlims_
    t_ = [collapseStructure(t.Short);collapseStructure(t.Medium);collapseStructure(t.Long)]*1e-6;
    t_(isnan(t_))=[];
    tlims_ = [min(t_), max(t_)]+[-5 5];
    [coh_,f] = CalculateSTcoherence([PFCcells;HPcells],1e-4,5,tlims_);
    FreqMod.SpkSpkcoh{iFile} = coh_;
    FreqMod.SpkSpkcoh_f = f;
    %% Spike-Spike coherence - event windows
%     clear ModIndex powWelch
%     events = fieldnames(STranges.PFC);
%     for iEvent=1:length(events)
%         ST_PFC = eval(sprintf('STranges.PFC.%s',events{iEvent}));
%         ST_HP  = eval(sprintf('STranges.HP.%s',events{iEvent}));
%         [coh_,f] = CalculateSTcoherence([ST_PFC;ST_HP],1e-4,5);
%         
%         eval(sprintf('FreqModEvt.SpkSpkcoh.%s{iFile} = coh_;',events{iEvent}))
%         FreqModEvt.SpkSpkcoh_f = f;
%         
%     end
end

%% Collapse - all spikes

% units - firing rate modulation
for iArea = 1:2
    for iFile = 1:length(fileList)            
        idx = D_Ass.usel_out{iArea}{iFile};
        FreqMod_ = FreqMod.ModIndex{iFile}{iArea}(:,idx);
        FreqMod.ModIndexNon{iArea}{iFile}   = FreqMod_(:,D.Membership{iArea}{iFile}==0);
        FreqMod.ModIndexLocal{iArea}{iFile} = FreqMod_(:,D.Membership{iArea}{iFile}==1);
        FreqMod.ModIndexJoint{iArea}{iFile} = FreqMod_(:,D.Membership{iArea}{iFile}==2);
    end        
end


% units  - spike time coherence
for iFile = 1:length(fileList)
    for iArea = 1:2
        
        % Collate non-member pairs
        idx = find(D.Membership{iArea}{iFile}==0);
%         idx = find(nansum([D.LocalMembersMatrix{iArea}{iFile},D.JointMembersMatrix{iArea}{iFile}],2)==0);
        if iArea == 2
            idx=idx+D.nUnits{1}(iFile);
        end
        temp = [];
        IDs  = [];
        for i = 1:length(idx)
            for j = 1:length(idx)
                if i < j
                    temp = [temp;FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)}];
                    if ~isempty(FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)})
                        IDs  = [IDs;idx(i),idx(j)];
                    end
                end
            end
        end
        FreqMod.SpkSpkcoh_Collapsed.Non{iArea}{iFile}     = temp;
        FreqMod.SpkSpkcoh_Collapsed.Non_IDs{iArea}{iFile} = IDs;
        
        % Collate local member pairs
        membership_ = D.LocalMembersMatrix{iArea}{iFile};
        if size(membership_,2)>1
            temp = [];
            IDs  = [];
            for iAss = 2:size(membership_,2)
                idx = find(membership_(:,iAss));
                if iArea ==2
                    idx=idx+D.nUnits{1}(iFile);
                end
                for i = 1:length(idx)
                    for j = 1:length(idx)
                        if i < j
                            temp = [temp;FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)}];
                            if ~isempty(FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)})
                                IDs  = [IDs;idx(i),idx(j)];
                            end
                        end
                    end
                end
            end
        end
        if ~isempty(temp)
            FreqMod.SpkSpkcoh_Collapsed.Local{iArea}{iFile}     = temp;
            FreqMod.SpkSpkcoh_Collapsed.Local_IDs{iArea}{iFile} = IDs;
        else
%             FreqMod.SpkSpkcoh_Collapsed.Local{iArea}{iFile} = nan(length(FreqMod.SpkSpkcoh_f),1);
        end
    end
    
    for iArea = 1:3
        % Collate inter-area member pairs
        membership_ = D.JointMembersMatrix{iArea}{iFile};
        if size(membership_,2)>1
            temp = [];
            IDs = [];
            for iAss = 2:size(membership_,2)
                idx = find(membership_(:,iAss));
                if iArea == 2
                    idx=idx+D.nUnits{1}(iFile);
                end
                for i = 1:length(idx)
                    for j = 1:length(idx)
                        if iArea < 3
                            if i < j
                                temp = [temp;FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)}];
                                if ~isempty(FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)})
                                    IDs  = [IDs;idx(i),idx(j)];
                                end
                            end
                        else

                            if idx(i) <= D.nUnits{1}(iFile) && idx(j) > D.nUnits{1}(iFile)      
                                temp = [temp;FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)}];
                                if ~isempty(FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)})
                                    IDs  = [IDs;idx(i),idx(j)];
                                end
                            end
                        end
                        
                    end
                end
            end
        end
        if ~isempty(temp)
            FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}{iFile}     = temp;
            FreqMod.SpkSpkcoh_Collapsed.Joint_IDs{iArea}{iFile} = IDs;
        else
            % FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}{iFile} = nan(length(FreqMod.SpkSpkcoh_f),1);
        end
    end
    
    iArea = 3;
    % Collate inter-area non-member pairs
    membership_ = D.JointMembersMatrix{iArea}{iFile};
    if size(membership_,2)>1
        idx = find(sum(membership_(:,2:end),2)==0);
        temp = [];
        IDs  = [];
        for i =1:length(idx)
            for j =1:length(idx)
                if idx(i) <= D.nUnits{1}(iFile) && idx(j)>D.nUnits{1}(iFile)
                    temp = [temp;FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)}];
                    if ~isempty(FreqMod.SpkSpkcoh{iFile}{idx(i),idx(j)})
                        IDs  = [IDs;idx(i),idx(j)];
                    end
                end
            end
        end
    end
    if ~isempty(temp)
        FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers{iArea}{iFile}     = temp;
        FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers_IDs{iArea}{iFile} = IDs;
    else
        %  FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers{iArea}{iFile} = nan(length(FreqMod.SpkSpkcoh_f),1);
    end
    
    % Collate inter-area cross assembly-member pairs
    membership_ = D.JointMembersMatrix{iArea}{iFile};
    if size(membership_,2)>1
        temp = [];
        IDs  = [];
        for iAss = 2:size(membership_,2)
            members = find(membership_(:,iAss));
            othermembers = setdiff(find(sum(membership_(:,2:end),2)),members);
            for i =1:length(members)
                for j =1:length(othermembers)
                    if members(i) <= D.nUnits{1}(iFile) && othermembers(j)>D.nUnits{1}(iFile)
                        
                        temp = [temp;FreqMod.SpkSpkcoh{iFile}{members(i),othermembers(j)}];
                        if ~isempty(FreqMod.SpkSpkcoh{iFile}{members(i),othermembers(j)})
                            IDs  = [IDs;members(i),othermembers(j)];
                        end
                    end
                end
            end
        end
    end
    if ~isempty(temp)
        FreqMod.SpkSpkcoh_Collapsed.Joint_CrossMembers{iArea}{iFile} = temp;
        FreqMod.SpkSpkcoh_Collapsed.Joint_CrossMembers_IDs{iArea}{iFile} = IDs;
    else
        %  FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers{iArea}{iFile} = nan(length(FreqMod.SpkSpkcoh_f),1);
    end
end
%% Plot SC vs LR decoding
%% Plot freq modulation (3) [Band]Hz vs. Peak SC decoding
tlimsAllSC = [-4 4];
tbAllSC = tlimsAllSC(1):bw:tlimsAllSC(2);%0:bw:2*sum(abs(tlimsAll));
Ltr=length(tbAllSC);
t0 = FindClosestIndex(tbAllSC,0);
% TRangeSC = 1:Ltr; % Whole trial
TRangeSC = 1:round(t0/2); % Pre-press only
% TRangeSC = (round(t0/2)+1):Ltr; % Post-press only

tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
Ltr = length(tbAll)*2;
TRangeLR = 1:Ltr; % Whole trial
% TRangeLR = 1:Ltr/2; % Sample only
% TRangeLR = (Ltr/2+1):Ltr; % Choice only
% TRangeLR = round(Ltr/4):round(Ltr*3/4); % Delay only

figure; 
for iArea =1:2
    D_LR   =  cell2mat(D.TS{iArea});
    D_LRidx =  cell2mat(D.Membership{iArea}');
    
    D_SC_   =  cell2mat(D_SC.TS{iArea});
    D_SCidx =  cell2mat(D_SC.Membership{iArea}');
    subplot(1,2,iArea); hold on

    temp_x = nanmax(D_LR(TRangeLR,D_LRidx ==0),1);
    temp_y = nanmax(D_SC_(TRangeSC,D_SCidx ==0),1);
    scatter(temp_x,temp_y,10,col_{1},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{1} = fitlm(temp_x,temp_y);
    
    temp_x = nanmax(D_LR(TRangeLR,D_LRidx ==1),1);
    temp_y = nanmax(D_SC_(TRangeSC,D_SCidx ==1),1);
    scatter(temp_x,temp_y,10,col_{2},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{2} = fitlm(temp_x,temp_y);

    temp_x = nanmax(D_LR(TRangeLR,D_LRidx ==2),1);
    temp_y = nanmax(D_SC_(TRangeSC,D_SCidx ==2),1);
    scatter(temp_x,temp_y,10,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{3} = fitlm(temp_x,temp_y);

    for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='none';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='none';
        h(4).Color = col_{i};
    end
    legend off
    title(Areas{iArea})
    
    ylabel('Peak Contextual decoding')    
    xlabel('Peak Spatial decoding')
    axis([0 20 0 20])
end

%% Plot freq modulation (1) spectra
f = FreqMod.f;
CollapseFun=@(x) cell2mat(x)';
% CollapseFun=@(x) cell2mat(cellfun(@(x) nanmean(x,2),x,'UniformOutput',false))';
figure
for iArea =1:2
    subplot(1,2,iArea); hold on
    y = CollapseFun(FreqMod.ModIndexNon{iArea});
%     plot(f,y,'Color',col_{1})
    plot(f,nanmean(y),'Color',col_{1},'Linewidth',1.5,'HandleVisibility','off')
    ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{1})
%     y = cell2mat(cellfun(@(x) mean(x,2),FreqMod.ModIndexNon{iArea},'UniformOutput',false))';
%     plot(f,y,'Color',col_{1},'Linewidth',1,'HandleVisibility','off')

    y = CollapseFun(FreqMod.ModIndexLocal{iArea}(~isempty_cell(FreqMod.ModIndexLocal{iArea})));
%     plot(f,y,'Color',col_{2}) 
    plot(FreqMod.f,nanmean(y),'Color',col_{2},'Linewidth',1.5,'HandleVisibility','off')
    ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{2})
%     y = cell2mat(cellfun(@(x) mean(x,2),FreqMod.ModIndexLocal{iArea},'UniformOutput',false))';
%     plot(f,y,'Color',col_{2},'Linewidth',1,'HandleVisibility','off')
    
    y = CollapseFun(FreqMod.ModIndexJoint{iArea}(~isempty_cell(FreqMod.ModIndexJoint{iArea})));
%     plot(f,y,'Color',col_{3})
    plot(FreqMod.f,nanmean(y),'Color',col_{3},'Linewidth',1.5,'HandleVisibility','off')
    ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{3})
%     y = cell2mat(cellfun(@(x) mean(x,2),FreqMod.ModIndexJoint{iArea},'UniformOutput',false))';
%     plot(f,y,'Color',col_{3},'Linewidth',1,'HandleVisibility','off')
    
    axis([1 20 0 2.5])
%     set(gca,'Xscale','log')
%     ylabel('Spike train modulation index')
%     xlabel('Frequency (Hz)')
%     title(Areas{iArea})
end
% legend(MemberClasses_,'Location','south'); legend boxoff
%% Plot freq modulation (2) 5Hz vs. theta

f = FreqMod.f;
fRange_5Hz = (f>=3.5 & f<6.5);
fRange_Theta = (f>=7.5 & f<12);
% fRange_5Hz = find(f>=4 & f<5.5);
% fRange_Theta = find(f>=8 & f<9);

figure; 
for iArea =1:2
    subplot(1,2,iArea); hold on
    plot([1 1],[0 6],':k')
    plot([0 6],[1 1],':k')
    plot([0 6],[0 6],':k')
    temp_ = cell2mat(FreqMod.ModIndexNon{iArea})';
    x = nanmax(temp_(:,fRange_5Hz),2);
    y = nanmax(temp_(:,fRange_Theta),2);
    mdl{1} = fitlm(x,y);
    scatter(x,y,10,col_{1},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    
    temp_ = cell2mat(FreqMod.ModIndexLocal{iArea})';
    x = nanmax(temp_(:,fRange_5Hz),2);
    y = nanmax(temp_(:,fRange_Theta),2);
    mdl{2} = fitlm(x,y);
    scatter(x,y,10,col_{2},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    
    temp_ = cell2mat(FreqMod.ModIndexJoint{iArea})';
    x = nanmax(temp_(:,fRange_5Hz),2);
    y = nanmax(temp_(:,fRange_Theta),2);
    mdl{3} = fitlm(x,y);
    scatter(x,y,10,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    
     for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='-';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='-';
        h(4).Color = col_{i};
    end
    legend off
    title(Areas{iArea})
  
    xlabel('4-6Hz modulation')    
    ylabel('theta modulation')
    if iArea==1
        axis([0 6 0 6])
    else
        axis([0 6 0 6])
    end
        
end

%% Plot freq modulation (3) [Band]Hz vs. Peak LR decoding
f = FreqMod.f;
fRange_ = (f>=3.5 & f<6.5);   %5Hz
% fRange_ = (f>=7.5 & f<9);  %Theta

tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
Ltr = length(tbAll)*2;
% TRange = 1:Ltr; % Whole trial
% TRange = 1:Ltr/2; % Sample only
% TRange = (Ltr/2+1):Ltr; % Choice only
TRange = round(Ltr/4):round(Ltr*3/4); % Delay only

figure; 
for iArea =1:2
    D_   =  cell2mat(D.TS{iArea});
    Didx =  cell2mat(D.Membership{iArea}');
    subplot(1,2,iArea); hold on
    
    temp_x = cell2mat(FreqMod.ModIndexNon{iArea})';
    temp_x = nanmax(temp_x(:,fRange_),2);
    temp_y = nanmax(D_(TRange,Didx ==0),1);
    scatter(temp_x,temp_y,20,col_{1},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{1} = fitlm(temp_x,temp_y);
    
    temp_x = cell2mat(FreqMod.ModIndexLocal{iArea})';
    temp_x = nanmax(temp_x(:,fRange_),2);
    temp_y = nanmax(D_(TRange,Didx ==1),1);
    
    scatter(temp_x,temp_y,20,col_{2},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{2} = fitlm(temp_x,temp_y);
    
    temp_x = cell2mat(FreqMod.ModIndexJoint{iArea})';
    temp_x = nanmax(temp_x(:,fRange_),2);
    temp_y = nanmax(D_(TRange,Didx ==2),1);
    scatter(temp_x,temp_y,20,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{3} = fitlm(temp_x,temp_y);
    
    
    pFit{iArea} = [];RFit{iArea} = [];
    for i=1:3
       pFit{iArea} = [pFit{iArea};coefTest(mdl{i})];
       RFit{iArea} = [RFit{iArea};mdl{i}.Rsquared.Adjusted];
    end
    for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='-';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='-';
        h(4).Color = col_{i};
    end
    
    legend off
    title(Areas{iArea})
    
    xlabel('Freq. modulation')    
    ylabel('Peak spatial decoding')
%     axis([0 2 0 20])
end
%% Plot freq modulation (4) ratio vs  Peak LR decoding
CompFun = @(x,y) (x-y)./(x+y);
f = FreqMod.f;
fRange_5Hz   = (f>=3.5 & f<6.5);
fRange_Theta = (f>=7.5 & f<12);


tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
Ltr = length(tbAll)*2;
% TRange = 1:Ltr; % Whole trial
% TRange = 1:Ltr/2; % Sample only
% TRange = (Ltr/2+1):Ltr; % Choice only
TRange = round(Ltr/4):round(Ltr*3/4); % Delay only

figure; 
for iArea =1:2
    D_   =  cell2mat(D.TS{iArea});
    Didx =  cell2mat(D.Membership{iArea}');
    subplot(1,2,iArea); hold on
    temp_x = cell2mat(FreqMod.ModIndexNon{iArea})';
    temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));
    temp_y = max(D_(TRange,Didx ==0));
    scatter(temp_x,temp_y,20,col_{1},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{1} = fitlm(temp_x,temp_y);
    
    temp_x = cell2mat(FreqMod.ModIndexLocal{iArea})';
    temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));
    temp_y = max(D_(TRange,Didx ==1));
    scatter(temp_x,temp_y,20,col_{2},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{2} = fitlm(temp_x,temp_y);
    
    temp_x = cell2mat(FreqMod.ModIndexJoint{iArea})';
    temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));
    temp_y = max(D_(TRange,Didx ==2));
    scatter(temp_x,temp_y,20,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{3} = fitlm(temp_x,temp_y);
    
    for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='none';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='none';
        h(4).Color = col_{i};
    end
    legend off
    title(Areas{iArea})
end
%% Plot freq modulation (3) [Band]Hz vs. Peak SC decoding
f = FreqMod.f;
% fRange_ = (f>=3.5 & f<6.5);   %5Hz
fRange_ = (f>=7.5 & f<9);  %Theta


tlimsAllSC = [-4 4];
tbAllSC = tlimsAllSC(1):bw:tlimsAllSC(2);%0:bw:2*sum(abs(tlimsAll));
Ltr=length(tbAllSC);
t0 = FindClosestIndex(tbAllSC,0);
% TRange = 1:Ltr; % Whole trial
TRange = 1:round(t0/2); % Pre-press only
% TRange = (round(t0/2)+1):Ltr; % Post-press only

figure; 
for iArea =1:2
    D_   =  cell2mat(D_SC.TS{iArea});

    Didx =  cell2mat(D_SC.Membership{iArea}');
    subplot(1,2,iArea); hold on

    temp_x = cell2mat(FreqMod.ModIndexNon{iArea})';
    temp_x = nanmax(temp_x(:,fRange_),2);
    temp_y = nanmax(D_(TRange,Didx ==0),1);
    scatter(temp_x,temp_y,20,col_{1},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{1} = fitlm(temp_x,temp_y);
    
    temp_x = cell2mat(FreqMod.ModIndexLocal{iArea})';
    temp_x = nanmax(temp_x(:,fRange_),2);
    temp_y = nanmax(D_(TRange,Didx ==1),1);
    scatter(temp_x,temp_y,20,col_{2},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{2} = fitlm(temp_x,temp_y);
    
    temp_x = cell2mat(FreqMod.ModIndexJoint{iArea})';
    temp_x = nanmax(temp_x(:,fRange_),2);
    temp_y = nanmax(D_(TRange,Didx ==2),1);
    scatter(temp_x,temp_y,20,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{3} = fitlm(temp_x,temp_y);
    
    pFit{iArea} = [];RFit{iArea} = [];
    for i=1:3
       pFit{iArea} = [pFit{iArea};coefTest(mdl{i})];
       RFit{iArea} = [RFit{iArea};mdl{i}.Rsquared.Adjusted];
    end
    for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='-';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='-';
        h(4).Color = col_{i};
    end
    legend off
    title(Areas{iArea})
    
    xlabel('Freq. modulation')    
    ylabel('Peak contextual decoding')
%     axis([0 6 0 40])
end
%% Plot freq modulation (4) ratio vs  Peak SC decoding
CompFun = @(x,y) (x-y)./(x+y);
f = FreqMod.f;
fRange_5Hz   = (f>=3.5 & f<6.5);
fRange_Theta = (f>=7.5 & f<12);
% fRange_5Hz = find(f>=4 & f<5.5);
% fRange_Theta = find(f>=8 & f<9);

tlimsAllSC = [-4 4];
tbAllSC = tlimsAllSC(1):bw:tlimsAllSC(2);%0:bw:2*sum(abs(tlimsAll));
Ltr=length(tbAllSC);
t0 = FindClosestIndex(tbAllSC,0);
% TRange = 1:Ltr; % Whole trial
TRange = 1:round(t0/2); % Pre-press only
% TRange = (round(t0/2)+1):Ltr; % Post-press only

figure; 
for iArea =1:2
    D_   =  cell2mat(D_SC.TS{iArea});
    Didx =  cell2mat(D_SC.Membership{iArea}');
    subplot(1,2,iArea); hold on
    temp_x = cell2mat(FreqMod.ModIndexNon{iArea})';
    temp_x = CompFun(nanmax(temp_x(:,fRange_5Hz),2),nanmax(temp_x(:,fRange_Theta),2));
    temp_y = max(D_(TRange,Didx ==0));
    scatter(temp_x,temp_y,20,col_{1},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{1} = fitlm(temp_x,temp_y);
    
    temp_x = cell2mat(FreqMod.ModIndexLocal{iArea})';
    temp_x = CompFun(nanmax(temp_x(:,fRange_5Hz),2),nanmax(temp_x(:,fRange_Theta),2));
    temp_y = max(D_(TRange,Didx ==1));
    scatter(temp_x,temp_y,20,col_{2},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{2} = fitlm(temp_x,temp_y);
    
    temp_x = cell2mat(FreqMod.ModIndexJoint{iArea})';
    temp_x = CompFun(nanmax(temp_x(:,fRange_5Hz),2),nanmax(temp_x(:,fRange_Theta),2));
    temp_y = max(D_(TRange,Didx ==2));
    scatter(temp_x,temp_y,20,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{3} = fitlm(temp_x,temp_y);
    
    for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='none';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='none';
        h(4).Color = col_{i};
    end
    legend off
    title(Areas{iArea})
end

%% Plot spike train coherence   (1) spectra
smooth_ = 0;
yMax = 0.03;
CollapseFun=@(x) cell2mat(x')';
% CollapseFun=@(x) cell2mat(cellfun(@(x) nanmean(x,1),x,'UniformOutput',false)')';
f = FreqMod.SpkSpkcoh_f;
figure; 
for iArea = 1:2
    subplot(1,3,iArea); hold on
%     title
    y = CollapseFun(FreqMod.SpkSpkcoh_Collapsed.Non{iArea});
    y = smooth2a(y,smooth_,0);
    ciplot(nanmean(y,2)+nansem(y,2),nanmean(y,2)-nansem(y,2),f,col_{1},0.6)
	plot(f,nanmean(y,2),'color',col_{1},'LineWidth',1.5)
    y = CollapseFun(FreqMod.SpkSpkcoh_Collapsed.Local{iArea}(~isempty_cell(FreqMod.SpkSpkcoh_Collapsed.Local{iArea})));
    y = smooth2a(y,smooth_,0);
    plot(f,nanmean(y,2),'color',col_{2},'LineWidth',1.5)
    ciplot(nanmean(y,2)+nansem(y,2),nanmean(y,2)-nansem(y,2),f,col_{2},0.6)
    y = CollapseFun(FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}(~isempty_cell(FreqMod.SpkSpkcoh_Collapsed.Joint{iArea})));
    y = smooth2a(y,smooth_,0);
    ciplot(nanmean(y,2)+nansem(y,2),nanmean(y,2)-nansem(y,2),f,col_{3},0.6)
    plot(f,nanmean(y,2),'color',col_{3},'LineWidth',1.5)
%     plot(f,smooth2a(y,smooth_,0),'color',col_{3},'LineWidth',1)
    axis([min(FreqMod.SpkSpkcoh_f) max(FreqMod.SpkSpkcoh_f) 0 yMax])
end
subplot(1,3,3); hold on
y = CollapseFun(FreqMod.SpkSpkcoh_Collapsed.Joint{3});
y = smooth2a(y,smooth_,0);
ciplot(nanmean(y,2)+nansem(y,2),nanmean(y,2)-nansem(y,2),f,col_{3},0.6)
plot(f,nanmean(y,2),'color',col_{3},'LineWidth',1.5)

y =  CollapseFun(FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers{3});
y = smooth2a(y,smooth_,0);
ciplot(nanmean(y,2)+nansem(y,2),nanmean(y,2)-nansem(y,2),f,col_{1},0.6)
plot(f,nanmean(y,2),'color',col_{1},'LineWidth',1.5)

y =  CollapseFun(FreqMod.SpkSpkcoh_Collapsed.Joint_CrossMembers{3}(~isempty_cell(FreqMod.SpkSpkcoh_Collapsed.Joint_CrossMembers{3})));
y = smooth2a(y,smooth_,0);
ciplot(nanmean(y,2)+nansem(y,2),nanmean(y,2)-nansem(y,2),f,col_{3},0.2)
plot(f,nanmean(y,2),'color',col_{3},'LineWidth',1.5,'LineStyle',':')

axis([min(FreqMod.SpkSpkcoh_f) max(FreqMod.SpkSpkcoh_f) 0 yMax])
%% Plot spike train coherence   (2) 5Hz vs. theta
xMax = 0.6;
yMax = 0.3;
f = FreqMod.SpkSpkcoh_f;
fRange_5Hz = (f>=3.5 & f<6.5);
fRange_Theta = (f>=7.5 & f<12);
% fRange_5Hz = find(f>=4 & f<5.5);
% fRange_Theta = find(f>=8 & f<9);

figure; 
for iArea =1:2
    subplot(1,3,iArea); hold on
    plot([0 1],[0 1],':k')
    
    temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Non{iArea}')';
    x = nanmax(temp_(fRange_5Hz,:),1);
    y = nanmax(temp_(fRange_Theta,:),1);
    idx = y>0.5;x(idx)=[];y(idx)=[];
    scatter(x,y,20,col_{1},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{1} = fitlm(x,y);
    
    temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Local{iArea}')';
    x = nanmax(temp_(fRange_5Hz,:),1);
    y = nanmax(temp_(fRange_Theta,:),1);
    idx = y>0.5;x(idx)=[];y(idx)=[];
    scatter(x,y,20,col_{2},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{2} = fitlm(x,y);
    
    temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}')';
    x = nanmax(temp_(fRange_5Hz,:),1);
    y = nanmax(temp_(fRange_Theta,:),1);
    scatter(x,y,20,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{3} = fitlm(x,y);
    axis([0 xMax 0 yMax])    
    
     for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='-';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='-';
        h(4).Color = col_{i};
    end
    legend off
    title(Areas{iArea})
    if iArea==1
        ylabel('Theta spike train coherence')
        xlabel('')
    elseif iArea==2
        xlabel('5Hz spike train coherence')    
        ylabel('')
    end

%     set(gca,'XScale','log','YScale','log')
%     if iArea ==1
%         axis([0 0.02 0 0.02])
%     else
%         axis([0 0.02 0 0.02])        
%     end
end

iArea = 3;
subplot(1,3,iArea); hold on
plot([0 1],[0 1],':k')

temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers{iArea}')';
x = nanmax(temp_(fRange_5Hz,:),1);
y = nanmax(temp_(fRange_Theta,:),1);
scatter(x,y,20,col_{1},'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
    mdl{2} = fitlm(x,y);

temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_CrossMembers{iArea}')';
x = nanmax(temp_(fRange_5Hz,:),1);
y = nanmax(temp_(fRange_Theta,:),1);
scatter(x,y,20,col_{3},'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
    mdl{2} = fitlm(x,y);
        
        
temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}')';
x = nanmax(temp_(fRange_5Hz,:),1);
y = nanmax(temp_(fRange_Theta,:),1);
scatter(x,y,20,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{3} = fitlm(x,y);
    
     for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='-';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='-';
        h(4).Color = col_{i};
    end
    legend off
    title(Areas{iArea})
%     xlabel('5Hz Coherence')
%     ylabel('Theta Coherence')
    xlabel('');ylabel('')
axis([0 xMax 0 yMax])    
% set(gca,'XScale','log','YScale','log')
%% Plot spike train coherence   (3) [Band]Hz vs Peak LR decoding
AggregateFun =@(x,y) sqrt(x.^2 + y.^2);
% AggregateFun =@(x,y) (x + y)./2;
% AggregateFun =@(x,y) sum([x; y])';
% AggregateFun =@(x,y) max([x; y])';
xMin = -5;
xMax = 0;
yMax = 7;
f = FreqMod.SpkSpkcoh_f;

% fRange_ = (f>=3.5 & f<6.5); % 5Hz
fRange_ = (f>=7.5 & f<12); %Theta

TRange = 1:Ltr; % Whole trial
% TRange = 1:Ltr/2; % Sample only
% TRange = (Ltr/2+1):Ltr; % Choice only
figure; 
for iArea =1:2
    subplot(1,3,iArea); hold on
    
    temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Non{iArea}');
    x = nanmean(temp_(:,fRange_),2);x=log10(x);
    D_ = cell2mat(D.TS{iArea});
    D_   =  nanmax(D_(TRange,:),1);
    Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Non_IDs{iArea}');
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
    scatter(x,y,20,col_{1},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{1} = fitlm(x,y);
    
    temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Local{iArea}');
    x = nanmean(temp_(:,fRange_),2);x=log10(x);
    D_ = cell2mat(D.TS{iArea});
    D_   =  nanmax(D_(TRange,:),1);
    Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Local_IDs{iArea}');
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
    scatter(x,y,20,col_{2},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{2} = fitlm(x,y);
    
    temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}');
    x = nanmean(temp_(:,fRange_),2);x=log10(x);
    D_ = cell2mat(D.TS{iArea});
    D_   =  nanmax(D_(TRange,:),1);
    Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_IDs{iArea}');
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
    scatter(x,y,20,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{3} = fitlm(x,y);
    
    pFit{iArea} = [];RFit{iArea} = [];
    for i=1:3
       pFit{iArea} = [pFit{iArea};coefTest(mdl{i})];
       RFit{iArea} = [RFit{iArea};mdl{i}.Rsquared.Adjusted];
    end
    for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='-';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='-';
        h(4).Color = col_{i};
    end
    legend off
    title(Areas{iArea})
    
    axis([xMin xMax 0 yMax]) 
    if iArea ==2
        xlabel('log(Coherence Magnitude)')
        ylabel('')
    else
        ylabel('Joint max. spatial decoding')
        xlabel('')
    end
%     set(gca,'Xscale','log')

end
iArea = 3;
subplot(1,3,iArea); hold on

temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers{iArea}');
x = nanmax(temp_(:,fRange_),2);x=log10(x);
D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D.TS{1},D.TS{2},'UniformOutput',false));
D_   =  nanmax(D_(TRange,:),1);
Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers_IDs{iArea}');
y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
scatter(x,y,20,col_{1},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{1} = fitlm(x,y);

temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_CrossMembers{iArea}');
x = nanmax(temp_(:,fRange_),2);x=log10(x);
D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D.TS{1},D.TS{2},'UniformOutput',false));
D_   =  nanmax(D_(TRange,:));
Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_CrossMembers_IDs{iArea}');
y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
scatter(x,y,20,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{2} = fitlm(x,y);
    
temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}');
x = nanmax(temp_(:,fRange_),2);x=log10(x);
D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D.TS{1},D.TS{2},'UniformOutput',false));
D_   =  nanmax(D_(TRange,:));
Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_IDs{iArea}');
y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
scatter(x,y,20,col_{3},'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0)
    mdl{3} = fitlm(x,y);
    pFit{iArea} = [];RFit{iArea} = [];
    for i=1:3
       pFit{iArea} = [pFit{iArea};coefTest(mdl{i})];
       RFit{iArea} = [RFit{iArea};mdl{i}.Rsquared.Adjusted];
    end
    for i=1:3
        if i~=2
            h = plot(mdl{i});
            h(1).Marker = 'none';
            h(2).LineWidth=2;
            h(2).LineStyle='-';
            h(2).Color = col_{i};
            h(3).LineWidth=1;
            h(3).LineStyle='-';
            h(3).Color = col_{i};
            h(4).LineWidth=1;
            h(4).LineStyle='-';
            h(4).Color = col_{i};
        else
            h = plot(mdl{i});
            h(1).Marker = 'none';
            h(2).LineWidth=2;
            h(2).LineStyle=':';
            h(2).Color = col_{3};
            h(3).LineWidth=1;
            h(3).LineStyle='-';
            h(3).Color = col_{3};
            h(4).LineWidth=1;
            h(4).LineStyle='-';
            h(4).Color = col_{3};
        end
    end
    legend off
    title(Areas{iArea})
    
axis([xMin xMax 0 yMax])    
% set(gca,'Xscale','log')
xlabel('');ylabel('')
%% Plot spike train coherence   (3) [Band]Hz vs Peak SC decoding
AggregateFun =@(x,y) sqrt(x.^2 + y.^2);
% AggregateFun =@(x,y) (x + y)./2;
% AggregateFun =@(x,y) max([x; y])';
xMin = -5;
xMax = 0;
yMax = 12;
f = FreqMod.SpkSpkcoh_f;
% fRange_ = (f>=3.5 & f<6.5); % 5Hz
fRange_ = (f>=7.5 & f<12); %Theta

tlimsAllSC = [-4 4];
tbAllSC = tlimsAllSC(1):bw:tlimsAllSC(2);%0:bw:2*sum(abs(tlimsAll));
Ltr=length(tbAllSC);
t0 = FindClosestIndex(tbAllSC,0);
% TRange = 1:Ltr; % Whole trial
TRange = 1:round(t0/2); % Pre-press only
% TRange = (round(t0/2)+1):Ltr; % Post-press only

figure; 
for iArea =1:2
    subplot(1,3,iArea); hold on
    
    temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Non{iArea}');
    x = nanmax(temp_(:,fRange_),2);x = log10(x);
    D_ = cell2mat(D_SC.TS{iArea});
    D_   =  nanmax(D_(TRange,:),1);
    Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Non_IDs{iArea}');
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
    scatter(x,y,20,col_{1},'filled')
	mdl{1} = fitlm(x,y);
    
    temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Local{iArea}');
    x = nanmax(temp_(:,fRange_),2);x = log10(x);
    D_ = cell2mat(D_SC.TS{iArea});
    D_   =  nanmax(D_(TRange,:),1);
    Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Local_IDs{iArea}');
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
    scatter(x,y,20,col_{2},'filled')
    mdl{2} = fitlm(x,y);
    
    temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}');
    x = nanmax(temp_(:,fRange_),2);x = log10(x);
    D_ = cell2mat(D_SC.TS{iArea});
    D_   =  nanmax(D_(TRange,:),1);
    Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_IDs{iArea}');
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
    scatter(x,y,20,col_{3},'filled')
    axis([xMin xMax 0 yMax])        
%     set(gca,'Xscale','log')
    mdl{3} = fitlm(x,y);

    pFit{iArea} = [];RFit{iArea} = [];
    for i=1:3
       pFit{iArea} = [pFit{iArea};coefTest(mdl{i})];
       RFit{iArea} = [RFit{iArea};mdl{i}.Rsquared.Adjusted];
    end
    for i=1:3
        h = plot(mdl{i});
        h(1).Marker = 'none';
        h(2).LineWidth=2;
        h(2).LineStyle='-';
        h(2).Color = col_{i};
        h(3).LineWidth=1;
        h(3).LineStyle='-';
        h(3).Color = col_{i};
        h(4).LineWidth=1;
        h(4).LineStyle='-';
        h(4).Color = col_{i};
    end
    legend off
    title(Areas{iArea})
    
    axis([xMin xMax 0 yMax]) 
    if iArea ==2
        xlabel('log(Coherence Magnitude)')
        ylabel('')
    else
        ylabel('Joint max. contextual decoding')
        xlabel('')
    end
%     set(gca,'Xscale','log')

end
iArea = 3;
subplot(1,3,iArea); hold on

temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers{iArea}');
x = nanmax(temp_(:,fRange_),2);x = log10(x);
D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D_SC.TS{1},D_SC.TS{2},'UniformOutput',false));
D_   =  nanmax(D_(TRange,:),1);
Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers_IDs{iArea}');
y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
scatter(x,y,20,col_{1})

temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_CrossMembers{iArea}');
x = nanmax(temp_(:,fRange_),2);x = log10(x);
D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D_SC.TS{1},D_SC.TS{2},'UniformOutput',false));
D_   =  nanmax(D_(TRange,:));
Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_CrossMembers_IDs{iArea}');
y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
scatter(x,y,20,col_{3})

temp_ = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}');
x = nanmax(temp_(:,fRange_),2);x = log10(x);
D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D_SC.TS{1},D_SC.TS{2},'UniformOutput',false));
D_   =  nanmax(D_(TRange,:));
Didx = cell2mat(FreqMod.SpkSpkcoh_Collapsed.Joint_IDs{iArea}');
y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
idx = isinf(x) |isinf(y');x(idx)=[]; y(idx)=[];
scatter(x,y,20,col_{3},'filled')
    mdl{3} = fitlm(x,y);
    pFit{iArea} = [];RFit{iArea} = [];
    for i=1:3
       pFit{iArea} = [pFit{iArea};coefTest(mdl{i})];
       RFit{iArea} = [RFit{iArea};mdl{i}.Rsquared.Adjusted];
    end
    for i=1:3
        if i~=2
            h = plot(mdl{i});
            h(1).Marker = 'none';
            h(2).LineWidth=2;
            h(2).LineStyle='-';
            h(2).Color = col_{i};
            h(3).LineWidth=1;
            h(3).LineStyle='-';
            h(3).Color = col_{i};
            h(4).LineWidth=1;
            h(4).LineStyle='-';
            h(4).Color = col_{i};
        else
            h = plot(mdl{i});
            h(1).Marker = 'none';
            h(2).LineWidth=2;
            h(2).LineStyle=':';
            h(2).Color = col_{3};
            h(3).LineWidth=1;
            h(3).LineStyle='-';
            h(3).Color = col_{3};
            h(4).LineWidth=1;
            h(4).LineStyle='-';
            h(4).Color = col_{3};
        end
    end
    legend off
    title(Areas{iArea})
    
axis([xMin xMax 0 yMax])    
% set(gca,'Xscale','log')
xlabel('');ylabel('')
%% Collapse - event-limited spikes
events = fieldnames(FreqModEvt.ModIndex);
for iEvent = 1:length(events)
    % units - firing rate modulation
    for iArea = 1:2
        for iFile = 1:length(fileList)
            idx = D_Ass.usel_out{iArea}{iFile};
            eval(sprintf('FreqMod_ = FreqModEvt.ModIndex.%s{iFile}{iArea}(:,idx);',events{iEvent}))
            eval(sprintf('FreqModEvt.ModIndexcollapsed.%s.ModIndexNon{iArea}{iFile}   = FreqMod_(:,D.Membership{iArea}{iFile}==0);',events{iEvent}))
            eval(sprintf('FreqModEvt.ModIndexcollapsed.%s.ModIndexLocal{iArea}{iFile} = FreqMod_(:,D.Membership{iArea}{iFile}==1);',events{iEvent}))
            eval(sprintf('FreqModEvt.ModIndexcollapsed.%s.ModIndexJoint{iArea}{iFile} = FreqMod_(:,D.Membership{iArea}{iFile}==2);',events{iEvent}))
            
            eval(sprintf('FreqMod_ = FreqModEvt.powWelch.%s{iFile}{iArea}(:,idx);',events{iEvent}))
            eval(sprintf('FreqModEvt.powWelchcollapsed.%s.powWelchNon{iArea}{iFile}   = FreqMod_(:,D.Membership{iArea}{iFile}==0);',events{iEvent}))
            eval(sprintf('FreqModEvt.powWelchcollapsed.%s.powWelchLocal{iArea}{iFile} = FreqMod_(:,D.Membership{iArea}{iFile}==1);',events{iEvent}))
            eval(sprintf('FreqModEvt.powWelchcollapsed.%s.powWelchJoint{iArea}{iFile} = FreqMod_(:,D.Membership{iArea}{iFile}==2);',events{iEvent}))
            
            
        end
    end
    
    
% units  - spike time coherence
for iFile = 1:length(fileList)
    for iArea = 1:2
        
        % Collate non-member pairs
        idx = find(D.Membership{iArea}{iFile}==0);

        if iArea == 2
            idx=idx+D.nUnits{1}(iFile);
        end
        temp = [];
        IDs  = [];
        for i = 1:length(idx)
            for j = 1:length(idx)
                if i < j
                    eval(sprintf('temp_ = FreqModEvt.SpkSpkcoh.%s{iFile}{idx(i),idx(j)};',events{iEvent}))
                    temp = [temp;temp_];
                    if ~isempty(temp_)
                        IDs  = [IDs;idx(i),idx(j)];
                    end
                end
            end
        end
        eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Non{iArea}{iFile}     = temp;',events{iEvent}));
        eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Non_IDs{iArea}{iFile} = IDs;',events{iEvent}));
        
        % Collate local member pairs
        membership_ = D.LocalMembersMatrix{iArea}{iFile};
        if size(membership_,2)>1
            temp = [];
            IDs  = [];
            for iAss = 2:size(membership_,2)
                idx = find(membership_(:,iAss));
                if iArea ==2
                    idx=idx+D.nUnits{1}(iFile);
                end
                for i = 1:length(idx)
                    for j = 1:length(idx)
                        if i < j
                            eval(sprintf('temp_ = FreqModEvt.SpkSpkcoh.%s{iFile}{idx(i),idx(j)};',events{iEvent}))
                            
                            temp = [temp;temp_];
                            if ~isempty(temp_)
                                IDs  = [IDs;idx(i),idx(j)];
                            end
                        end
                    end
                end
            end
        end
        if ~isempty(temp)
            eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Local{iArea}{iFile}     = temp;',events{iEvent}));
            eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Local_IDs{iArea}{iFile} = IDs;',events{iEvent}));
            
        else
%             FreqMod.SpkSpkcoh_Collapsed.Local{iArea}{iFile} = nan(length(FreqMod.SpkSpkcoh_f),1);
        end
    end
    
    for iArea = 1:3
        % Collate inter-area member pairs
        membership_ = D.JointMembersMatrix{iArea}{iFile};
        if size(membership_,2)>1
            temp = [];
            IDs = [];
            for iAss = 2:size(membership_,2)
                idx = find(membership_(:,iAss));
                if iArea == 2
                    idx=idx+D.nUnits{1}(iFile);
                end
                for i = 1:length(idx)
                    for j = 1:length(idx)
                        if iArea < 3
                            if i < j
                                eval(sprintf('temp_ = FreqModEvt.SpkSpkcoh.%s{iFile}{idx(i),idx(j)};',events{iEvent}))
                                temp = [temp;temp_];
                                if ~isempty(temp_)
                                    IDs  = [IDs;idx(i),idx(j)];
                                end
                            end
                        else

                            if idx(i) <= D.nUnits{1}(iFile) && idx(j) > D.nUnits{1}(iFile)      
                                eval(sprintf('temp_ = FreqModEvt.SpkSpkcoh.%s{iFile}{idx(i),idx(j)};',events{iEvent}))
                                temp = [temp;temp_];
                                if ~isempty(temp_)
                                    IDs  = [IDs;idx(i),idx(j)];
                                end
                            end
                        end
                        
                    end
                end
            end
        end
        if ~isempty(temp)
            eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint{iArea}{iFile}     = temp;',events{iEvent}));
            eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_IDs{iArea}{iFile} = IDs;',events{iEvent}));
        else
            % FreqMod.SpkSpkcoh_Collapsed.Joint{iArea}{iFile} = nan(length(FreqMod.SpkSpkcoh_f),1);
        end
    end
    
    iArea = 3;
    % Collate inter-area non-member pairs
    membership_ = D.JointMembersMatrix{iArea}{iFile};
    if size(membership_,2)>1
        idx = find(sum(membership_(:,2:end),2)==0);
        temp = [];
        IDs  = [];
        for i =1:length(idx)
            for j =1:length(idx)
                if idx(i) <= D.nUnits{1}(iFile) && idx(j)>D.nUnits{1}(iFile)
                    eval(sprintf('temp_ = FreqModEvt.SpkSpkcoh.%s{iFile}{idx(i),idx(j)};',events{iEvent}))
                    temp = [temp;temp_];
                    if ~isempty(temp_)
                        IDs  = [IDs;idx(i),idx(j)];
                    end
                end
            end
        end
    end
    if ~isempty(temp)
        eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_nonMembers{iArea}{iFile}     = temp;',events{iEvent}));
        eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_nonMembers_IDs{iArea}{iFile} = IDs;',events{iEvent}));

    else
        %  FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers{iArea}{iFile} = nan(length(FreqMod.SpkSpkcoh_f),1);
    end
    
    % Collate inter-area cross assembly-member pairs
    membership_ = D.JointMembersMatrix{iArea}{iFile};
    if size(membership_,2)>1
        temp = [];
        IDs  = [];
        for iAss = 2:size(membership_,2)
            members = find(membership_(:,iAss));
            othermembers = setdiff(find(sum(membership_(:,2:end),2)),members);
            for i =1:length(members)
                for j =1:length(othermembers)
                    if members(i) <= D.nUnits{1}(iFile) && othermembers(j)>D.nUnits{1}(iFile)
                        eval(sprintf('temp_ = FreqModEvt.SpkSpkcoh.%s{iFile}{members(i),othermembers(j)};',events{iEvent}))
                        temp = [temp;temp_];
                        if ~isempty(temp_)
                            IDs  = [IDs;members(i),othermembers(j)];
                        end
                    end
                end
            end
        end
    end
    if ~isempty(temp)
        eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_CrossMembers{iArea}{iFile}     = temp;',events{iEvent}));
        eval(sprintf('FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_CrossMembers_IDs{iArea}{iFile} = IDs;',events{iEvent}));
    else
        %  FreqMod.SpkSpkcoh_Collapsed.Joint_nonMembers{iArea}{iFile} = nan(length(FreqMod.SpkSpkcoh_f),1);
    end
end


end
%% Plot average freq modulation (1) spectra (event-limited)
Events_ ={'Sample','Delay','NosePoke','Choice','Reward'};
Events_Names ={'Cue-Sample','Delay Period','Tone-NosePoke','NosePoke-Choice','Reward'};
f = FreqModEvt.f;
CollapseFun=@(x) cell2mat(x)';
% CollapseFun=@(x) cell2mat(cellfun(@(x) nanmean(x,2),x,'UniformOutput',false))';
Types_ = {'Non','Local','Joint'};
maxY = 2;
for iType =1:length(Types_)
    figure
    for iArea =1:2
        for iEvent = 1:length(Events_)
            if iArea ==1
                subplot(2,length(Events_),iEvent); hold on
            else
                subplot(2,length(Events_),iEvent+length(Events_)); hold on
            end
            
            eval(sprintf('y = FreqModEvt.ModIndexcollapsed.%sCorr.ModIndex%s{iArea};',Events_{iEvent},Types_{iType}))
%                         eval(sprintf('y = FreqModEvt.powWelchcollapsed.%sCorr.powWelch%s{iArea};',Events_{iEvent},Types_{iType}))
            y = CollapseFun(y(~isempty_cell(y)));
%             y = 10*log10(y);
            plot(f,nanmean(y),'Color',col_{iType},'Linewidth',1.5,'HandleVisibility','off')
            ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{iType},0.9)
            
            eval(sprintf('y_ = FreqModEvt.ModIndexcollapsed.%sErr.ModIndex%s{iArea};',Events_{iEvent},Types_{iType}))
%             eval(sprintf('y_ = FreqModEvt.powWelchcollapsed.%sErr.powWelch%s{iArea};',Events_{iEvent},Types_{iType}))
            y_ = CollapseFun(y_(~isempty_cell(y_)));
%             y = 10*log10(y_);
            plot(f,nanmean(y_),'Color',col_{iType},'Linewidth',1.5,'LineStyle',':','HandleVisibility','off')
            ciplot(nanmean(y_)+nansem(y_),nanmean(y_)-nansem(y_),f,col_{iType},0.3)
            
            sig1 = permtest2vec(y',y_',100,0.05);
                    a = nan(size(sig1));a(sig1)=1;
%     plot(f,maxY*a-0.1,'-k','LineWidth',2.5)
    a(isnan(a))=0;
    [~,Xpk,Wpk] = findpeaks(a)
    for iPeak=1:length(Xpk)
       plot(f([Xpk(iPeak) Xpk(iPeak)+Wpk(iPeak)]),ones(1,2)*maxY-0.1,'-k','LineWidth',2.5) 
    end
%             title(Events_{iEvent})
            axis([1 20 0.5 maxY])
%             axis([1 20 -35 -30])
            %     set(gca,'Xscale','log')
            %     ylabel('Spike train modulation index')
            %     xlabel('Frequency (Hz)')
            %     title(Areas{iArea})
        end
    end
    % legend(MemberClasses_,'Location','south'); legend boxoff
end
%% Plot spike train coherence   (1) spectra (event-limited)
smooth_ = 0;
yMax = 0.03;
CollapseFun=@(x) cell2mat(x')';
% CollapseFun=@(x) cell2mat(cellfun(@(x) nanmean(x,1),x,'UniformOutput',false)')';
f = FreqModEvt.SpkSpkcoh_f;

for iType =1:length(Types_)
     figure
% hold on
    for iArea = 1:3
        for iEvent = 1:length(Events_)
            if iArea ==1
                subplot(3,length(Events_),iEvent); hold on
            else
                subplot(3,length(Events_),iEvent+(iArea-1)*length(Events_)); hold on
            end
            
            if iArea<3
            eval(sprintf('y = FreqModEvt.SpkSpkcoh_Collapsed.%sCorr.%s{iArea};',Events_{iEvent},Types_{iType}))
            y = CollapseFun(y(~isempty_cell(y)));
            y = smooth2a(y,smooth_,0)';    
            plot(f,nanmean(y),'Color',col_{iType},'Linewidth',1.5,'HandleVisibility','off')
            ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{iType},0.9)
            
%             eval(sprintf('y = FreqModEvt.SpkSpkcoh_Collapsed.%sErr.%s{iArea};',Events_{iEvent},Types_{iType}))
%             y = CollapseFun(y(~isempty_cell(y)));
%             y = smooth2a(y,smooth_,0)';    
%             plot(f,nanmean(y),'Color',col_{iType},'Linewidth',1.5,'LineStyle',':','HandleVisibility','off')
%             ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{iType},0.3)
            
            else
               if iType==3
                   eval(sprintf('y = FreqModEvt.SpkSpkcoh_Collapsed.%sCorr.%s{iArea};',Events_{iEvent},Types_{iType}))
                   y = CollapseFun(y(~isempty_cell(y)));
                   y = smooth2a(y,smooth_,0)';
                   plot(f,nanmean(y),'Color',col_{iType},'Linewidth',1.5,'HandleVisibility','off')
                   ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{iType},0.9)
                   
%                    eval(sprintf('y = FreqModEvt.SpkSpkcoh_Collapsed.%sErr.%s{iArea};',Events_{iEvent},Types_{iType}))
%                    y = CollapseFun(y(~isempty_cell(y)));
%                    y = smooth2a(y,smooth_,0)';
%                    plot(f,nanmean(y),'Color',col_{iType},'Linewidth',1.5,'LineStyle',':','HandleVisibility','off')
%                    ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{iType},0.3)

                   eval(sprintf('y = FreqModEvt.SpkSpkcoh_Collapsed.%sCorr.%s_nonMembers{iArea};',Events_{iEvent},Types_{iType}))
                   y = CollapseFun(y(~isempty_cell(y)));
                   y = smooth2a(y,smooth_,0)';
                   plot(f,nanmean(y),'Color',col_{1},'Linewidth',1.5,'HandleVisibility','off')
                   ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{1},0.9)
                   
                   eval(sprintf('y = FreqModEvt.SpkSpkcoh_Collapsed.%sCorr.%s_CrossMembers{iArea};',Events_{iEvent},Types_{iType}))
                   y = CollapseFun(y(~isempty_cell(y)));
                   y = smooth2a(y,smooth_,0)';
                   plot(f,nanmean(y),'Color',col_{3},'Linewidth',1.5,'LineStyle',':','HandleVisibility','off')
                   ciplot(nanmean(y)+nansem(y),nanmean(y)-nansem(y),f,col_{3},0.3)
               end
            end
            title(Events_{iEvent})
            axis([min(FreqModEvt.SpkSpkcoh_f) max(FreqModEvt.SpkSpkcoh_f) 0 yMax])

        end
    end
end

%% Plot freq modulation (3) [Band]Hz vs. Peak LR decoding
f = FreqMod.f;
fRange_ = (f>=3.5 & f<6.5); %5Hz
% fRange_ = (f>=7.5 & f<12); % Theta

tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
Ltr = length(tbAll)*2;
TRange = 1:Ltr; % Whole trial
% TRange = 1:Ltr/2; % Sample only
% TRange = (Ltr/2+1):Ltr; % Choice only
% TRange = round(Ltr/4):round(Ltr*3/4); % Delay only
for iEvent =1:length(events)
    figure('name',events{iEvent});
    for iArea =1:2
        D_   =  cell2mat(D.TS{iArea});
        %     D_   =  double(cell2mat(D.TSsigBS{iArea}));
        %     D_ = D_(range_,:);
        Didx =  cell2mat(D.Membership{iArea}');
        subplot(1,2,iArea); hold on
        %     plot([1 1],[0 6],':k')
        %     plot([0 6],[1 1],':k')
        %     temp_x = cell2mat(FreqMod.ModIndexNon{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexNon{iArea})'';',events{iEvent}));
        temp_x = nanmean(temp_x(:,fRange_),2);
        temp_y = nanmax(D_(TRange,Didx ==0),1);
        %     temp_y = sum(double(D_(:,Didx ==0)))*bw;
        scatter(temp_x,temp_y,20,col_{1},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{1},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        %     temp_x = cell2mat(FreqMod.ModIndexLocal{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexLocal{iArea})'';',events{iEvent}));
        
        temp_x = nanmax(temp_x(:,fRange_),2);
        temp_y = nanmax(D_(TRange,Didx ==1),1);
        %     temp_y = sum(double(D_(:,Didx ==1)))*bw;
        scatter(temp_x,temp_y,20,col_{2},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{2},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        %     temp_x = cell2mat(FreqMod.ModIndexJoint{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexJoint{iArea})'';',events{iEvent}));
        
        temp_x = nanmax(temp_x(:,fRange_),2);
        temp_y = nanmax(D_(TRange,Didx ==2),1);
        %     temp_y = sum(double(D_(:,Didx ==2)))*bw;
        scatter(temp_x,temp_y,20,col_{3},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        xlabel('Freq. modulation')
        ylabel('Peak Spatial decoding')
        axis([0 6 0 40])
    end
end
%% Plot freq modulation (3) Theta/5Hz mod ratio vs. Peak LR decoding
f = FreqMod.f;
fRange_5Hz   = (f>=3.5 & f<6.5); %5Hz
fRange_Theta = (f>=7.5 & f<12); % Theta

CompFun = @(x,y) (x-y)./(x+y);

tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
Ltr = length(tbAll)*2;
TRange = 1:Ltr; % Whole trial
% TRange = 1:Ltr/2; % Sample only
% TRange = (Ltr/2+1):Ltr; % Choice only
% TRange = round(Ltr/4):round(Ltr*3/4); % Delay only
for iEvent =1:length(events)
    figure('name',events{iEvent});
    for iArea =1:2
        D_   =  cell2mat(D.TS{iArea});
        %     D_   =  double(cell2mat(D.TSsigBS{iArea}));
        %     D_ = D_(range_,:);
        Didx =  cell2mat(D.Membership{iArea}');
        subplot(1,2,iArea); hold on
        %     plot([1 1],[0 6],':k')
        %     plot([0 6],[1 1],':k')
        %     temp_x = cell2mat(FreqMod.ModIndexNon{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexNon{iArea})'';',events{iEvent}));
        
        temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));

        
        temp_y = nanmax(D_(TRange,Didx ==0),1);
        %     temp_y = sum(double(D_(:,Didx ==0)))*bw;
        scatter(temp_x,temp_y,20,col_{1},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{1},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        %     temp_x = cell2mat(FreqMod.ModIndexLocal{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexLocal{iArea})'';',events{iEvent}));
        
        temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));
        temp_y = nanmax(D_(TRange,Didx ==1),1);
        %     temp_y = sum(double(D_(:,Didx ==1)))*bw;
        scatter(temp_x,temp_y,20,col_{2},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{2},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        %     temp_x = cell2mat(FreqMod.ModIndexJoint{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexJoint{iArea})'';',events{iEvent}));
        
        temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));
        temp_y = nanmax(D_(TRange,Didx ==2),1);
        %     temp_y = sum(double(D_(:,Didx ==2)))*bw;
        scatter(temp_x,temp_y,20,col_{3},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        xlabel('Freq. modulation')
        ylabel('Peak Spatial decoding')
        axis([-0.5 0.5 0 40])
    end
end

%% Plot freq modulation (3) [Band]Hz vs. Peak SC decoding
f = FreqMod.f;
fRange_ = (f>=3.5 & f<6.5); %5Hz
% fRange_ = (f>=7.5 & f<12); % Theta

tlimsAllSC = [-4 4];
tbAllSC = tlimsAllSC(1):bw:tlimsAllSC(2);%0:bw:2*sum(abs(tlimsAll));
Ltr=length(tbAllSC);
t0 = FindClosestIndex(tbAllSC,0); 
% TRange = 1:Ltr; % Whole trial
TRange = 1:round(t0/2); % Pre-press only
% TRange = (round(t0/2)+1):Ltr; % Post-press only

for iEvent =1:length(events)
    
    figure('name',events{iEvent});
    for iArea =1:2
        D_   =  cell2mat(D_SC.TS{iArea});
        Didx =  cell2mat(D_SC.Membership{iArea}');
        subplot(1,2,iArea); hold on
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexNon{iArea})'';',events{iEvent}));
        temp_x = nanmean(temp_x(:,fRange_),2);
        temp_y = nanmax(D_(TRange,Didx ==0),1);
        scatter(temp_x,temp_y,20,col_{1},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{1},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexLocal{iArea})'';',events{iEvent}));
        
        temp_x = nanmean(temp_x(:,fRange_),2);
        temp_y = nanmax(D_(TRange,Didx ==1),1);
        scatter(temp_x,temp_y,20,col_{2},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{2},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexJoint{iArea})'';',events{iEvent}));
        
        temp_x = nanmean(temp_x(:,fRange_),2);
        temp_y = nanmax(D_(TRange,Didx ==2),1);
        scatter(temp_x,temp_y,20,col_{3},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        xlabel('Freq. modulation')
        ylabel('Peak Spatial decoding')
        axis([0.5 2.5 0 40])
    end
end
%% Plot freq modulation (3) Theta/5Hz mod ratio vs. Peak LR decoding
f = FreqMod.f;
fRange_5Hz   = (f>=3.5 & f<6.5); %5Hz
fRange_Theta = (f>=7.5 & f<12); % Theta

CompFun = @(x,y) (x-y)./(x+y);

tlimsAllSC = [-4 4];
tbAllSC = tlimsAllSC(1):bw:tlimsAllSC(2);%0:bw:2*sum(abs(tlimsAll));
Ltr=length(tbAllSC);
t0 = FindClosestIndex(tbAllSC,0); 
% TRange = 1:Ltr; % Whole trial
TRange = 1:round(t0/2); % Pre-press only
% TRange = (round(t0/2)+1):Ltr; % Post-press only

for iEvent =1:length(events)
    figure('name',events{iEvent});
    for iArea = 1:2
        D_   =  cell2mat(D_SC.TS{iArea});
        %     D_   =  double(cell2mat(D_SC.TSsigBS{iArea}));
        %     D_ = D_(range_,:);
        Didx =  cell2mat(D_SC.Membership{iArea}');
        subplot(1,2,iArea); hold on
        %     plot([1 1],[0 6],':k')
        %     plot([0 6],[1 1],':k')
        %     temp_x = cell2mat(FreqMod.ModIndexNon{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexNon{iArea})'';',events{iEvent}));
        
        temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));

        
        temp_y = nanmax(D_(TRange,Didx ==0),1);
        %     temp_y = sum(double(D_(:,Didx ==0)))*bw;
        scatter(temp_x,temp_y,20,col_{1},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{1},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        %     temp_x = cell2mat(FreqMod.ModIndexLocal{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexLocal{iArea})'';',events{iEvent}));
        
        temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));
        temp_y = nanmax(D_(TRange,Didx ==1),1);
        %     temp_y = sum(double(D_(:,Didx ==1)))*bw;
        scatter(temp_x,temp_y,20,col_{2},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{2},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        %     temp_x = cell2mat(FreqMod.ModIndexJoint{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexJoint{iArea})'';',events{iEvent}));
        
        temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));
        temp_y = nanmax(D_(TRange,Didx ==2),1);
        %     temp_y = sum(double(D_(:,Didx ==2)))*bw;
        scatter(temp_x,temp_y,20,col_{3},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        xlabel('Freq. modulation')
        ylabel('Peak Spatial decoding')
        axis([-0.5 0.5 0 40])
    end
end

f = FreqMod.f;
fRange_5Hz   = (f>=3.5 & f<6.5); %5Hz
fRange_Theta = (f>=7.5 & f<12); % Theta

CompFun = @(x,y) (x-y)./(x+y);

tlimsAllSC = [-4 4];
tbAllSC = tlimsAllSC(1):bw:tlimsAllSC(2);%0:bw:2*sum(abs(tlimsAll));
Ltr=length(tbAllSC);
t0 = FindClosestIndex(tbAllSC,0); 
% TRange = 1:Ltr; % Whole trial
TRange = 1:round(t0/2); % Pre-press only
% TRange = (round(t0/2)+1):Ltr; % Post-press only

for iEvent =1:length(events)
    figure('name',events{iEvent});
    for iArea = 1:2
        D_   =  cell2mat(D_SC.TS{iArea});
        %     D_   =  double(cell2mat(D_SC.TSsigBS{iArea}));
        %     D_ = D_(range_,:);
        Didx =  cell2mat(D_SC.Membership{iArea}');
        subplot(1,2,iArea); hold on
        %     plot([1 1],[0 6],':k')
        %     plot([0 6],[1 1],':k')
        %     temp_x = cell2mat(FreqMod.ModIndexNon{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexNon{iArea})'';',events{iEvent}));
        
        temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));

        
        temp_y = nanmax(D_(TRange,Didx ==0),1);
        %     temp_y = sum(double(D_(:,Didx ==0)))*bw;
        scatter(temp_x,temp_y,20,col_{1},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{1},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        %     temp_x = cell2mat(FreqMod.ModIndexLocal{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexLocal{iArea})'';',events{iEvent}));
        
        temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));
        temp_y = nanmax(D_(TRange,Didx ==1),1);
        %     temp_y = sum(double(D_(:,Didx ==1)))*bw;
        scatter(temp_x,temp_y,20,col_{2},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{2},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        %     temp_x = cell2mat(FreqMod.ModIndexJoint{iArea})';
        eval(sprintf('temp_x = cell2mat(FreqModEvt.ModIndexcollapsed.%s.ModIndexJoint{iArea})'';',events{iEvent}));
        
        temp_x = CompFun(nanmean(temp_x(:,fRange_5Hz),2),nanmean(temp_x(:,fRange_Theta),2));
        temp_y = nanmax(D_(TRange,Didx ==2),1);
        %     temp_y = sum(double(D_(:,Didx ==2)))*bw;
        scatter(temp_x,temp_y,20,col_{3},'filled')
        fit_ = ezfit(temp_x,temp_y,'affine');
        showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        xlabel('Freq. modulation')
        ylabel('Peak Spatial decoding')
        axis([-0.5 0.5 0 40])
    end
end

%% Plot spike train coherence   (3) [Band]Hz vs Peak LR decoding
AggregateFun =@(x,y) sqrt(x.^2 + y.^2);
% AggregateFun =@(x,y) (x + y)./2;
% AggregateFun =@(x,y) max([x; y])';

xMax = 0.2;
yMax = 20;

% fRange_ = (f>=3.5 & f<6.5); %5Hz
fRange_ = (f>=7.5 & f<12); % Theta

tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
Ltr = length(tbAll)*2;
TRange = 1:Ltr; % Whole trial
% TRange = 1:Ltr/2; % Sample only
% TRange = (Ltr/2+1):Ltr; % Choice only
% TRange = round(Ltr/4):round(Ltr*3/4); % Delay only
for iEvent =1:length(events)
    
    figure('name',events{iEvent});
    for iArea =1:2
        subplot(1,3,iArea); hold on
        
        eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Non{iArea}'');',events{iEvent}));
        x = nanmean(temp_(:,fRange_),2);
        D_ = cell2mat(D.TS{iArea});
        D_   =  nanmax(D_(TRange,:),1);
        eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Non_IDs{iArea}'');',events{iEvent}));
        y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
        scatter(x,y,20,col_{1},'filled')
        fit_ = ezfit(x,y,'affine');
        %     showfit(fit_,'fitcolor',col_{1},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Local{iArea}'');',events{iEvent}));
        x = nanmean(temp_(:,fRange_),2);
        D_ = cell2mat(D.TS{iArea});
        D_   =  nanmax(D_(TRange,:),1);
        eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Local_IDs{iArea}'');',events{iEvent}));
        y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
        scatter(x,y,20,col_{2},'filled')
        fit_ = ezfit(x,y,'affine');
        %     showfit(fit_,'fitcolor',col_{2},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint{iArea}'');',events{iEvent}));
        x = nanmean(temp_(:,fRange_),2);
        D_ = cell2mat(D.TS{iArea});
        D_   =  nanmax(D_(TRange,:),1);
        eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_IDs{iArea}'');',events{iEvent}));
        y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
        scatter(x,y,20,col_{3},'filled')
        fit_ = ezfit(x,y,'affine');
        %     showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        axis([0 xMax 0 yMax])
%         set(gca,'Xscale','log')
        
    end
    iArea = 3;
    subplot(1,3,iArea); hold on
    
    eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_nonMembers{iArea}'');',events{iEvent}));
    x = nanmean(temp_(:,fRange_),2);
    D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D.TS{1},D.TS{2},'UniformOutput',false));
    D_   =  nanmax(D_(TRange,:),1);
    eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_nonMembers_IDs{iArea}'');',events{iEvent}));
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    scatter(x,y,20,col_{1})
    fit_ = ezfit(x,y,'affine');
    % showfit(fit_,'fitcolor',col_{1},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
    
    eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_CrossMembers{iArea}'');',events{iEvent}));    
    x = nanmean(temp_(:,fRange_),2);
    D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D.TS{1},D.TS{2},'UniformOutput',false));
    D_   =  nanmax(D_(TRange,:));
    eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_CrossMembers_IDs{iArea}'');',events{iEvent}));
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    scatter(x,y,20,col_{3})
    fit_ = ezfit(x,y,'affine');
    % showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
    
    eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint{iArea}'');',events{iEvent}));    
    x = nanmean(temp_(:,fRange_),2);
    D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D.TS{1},D.TS{2},'UniformOutput',false));
    D_   =  nanmax(D_(TRange,:));
    eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_IDs{iArea}'');',events{iEvent}));
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    scatter(x,y,20,col_{3},'filled')
    fit_ = ezfit(x,y,'affine');
    % showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
    axis([0 xMax 0 yMax])
%     set(gca,'Xscale','log')
end

%% Plot spike train coherence   (3) [Band]Hz vs Peak SC decoding
AggregateFun =@(x,y) sqrt(x.^2 + y.^2);
% AggregateFun =@(x,y) (x + y)./2;
% AggregateFun =@(x,y) max([x; y])';

xMax = 0.2;
yMax = 20;

fRange_ = (f>=3.5 & f<6.5); %5Hz
% fRange_ = (f>=7.5 & f<12); % Theta

tlimsAllSC = [-4 4];
tbAllSC = tlimsAllSC(1):bw:tlimsAllSC(2);%0:bw:2*sum(abs(tlimsAll));
Ltr=length(tbAllSC);
t0 = FindClosestIndex(tbAllSC,0); 
TRange = 1:Ltr; % Whole trial
% TRange = 1:round(t0/2); % Pre-press only
% TRange = (round(t0/2)+1):Ltr; % Post-press only

for iEvent =1:length(events)
    
    figure('name',events{iEvent});
    for iArea =1:2
        subplot(1,3,iArea); hold on
        
        eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Non{iArea}'');',events{iEvent}));
        x = nanmean(temp_(:,fRange_),2);
        D_ = cell2mat(D_SC.TS{iArea});
        D_   =  nanmax(D_(TRange,:),1);
        eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Non_IDs{iArea}'');',events{iEvent}));
        y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
        scatter(x,y,20,col_{1},'filled')
        fit_ = ezfit(x,y,'affine');
        %     showfit(fit_,'fitcolor',col_{1},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Local{iArea}'');',events{iEvent}));
        x = nanmean(temp_(:,fRange_),2);
        D_ = cell2mat(D_SC.TS{iArea});
        D_   =  nanmax(D_(TRange,:),1);
        eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Local_IDs{iArea}'');',events{iEvent}));
        y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
        scatter(x,y,20,col_{2},'filled')
        fit_ = ezfit(x,y,'affine');
        %     showfit(fit_,'fitcolor',col_{2},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        
        eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint{iArea}'');',events{iEvent}));
        x = nanmean(temp_(:,fRange_),2);
        D_ = cell2mat(D_SC.TS{iArea});
        D_   =  nanmax(D_(TRange,:),1);
        eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_IDs{iArea}'');',events{iEvent}));
        y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
        scatter(x,y,20,col_{3},'filled')
        fit_ = ezfit(x,y,'affine');
        %     showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
        axis([0 xMax 0 yMax])
%         set(gca,'Xscale','log')
        
    end
    iArea = 3;
    subplot(1,3,iArea); hold on
    
    eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_nonMembers{iArea}'');',events{iEvent}));
    x = nanmean(temp_(:,fRange_),2);
    D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D_SC.TS{1},D_SC.TS{2},'UniformOutput',false));
    D_   =  nanmax(D_(TRange,:),1);
    eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_nonMembers_IDs{iArea}'');',events{iEvent}));
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    scatter(x,y,20,col_{1})
    fit_ = ezfit(x,y,'affine');
    % showfit(fit_,'fitcolor',col_{1},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
    
    eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_CrossMembers{iArea}'');',events{iEvent}));    
    x = nanmean(temp_(:,fRange_),2);
    D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D_SC.TS{1},D_SC.TS{2},'UniformOutput',false));
    D_   =  nanmax(D_(TRange,:));
    eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_CrossMembers_IDs{iArea}'');',events{iEvent}));
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    scatter(x,y,20,col_{3})
    fit_ = ezfit(x,y,'affine');
    % showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
    
    eval(sprintf('temp_ = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint{iArea}'');',events{iEvent}));    
    x = nanmean(temp_(:,fRange_),2);
    D_   = cell2mat(cellfun(@(x,y) cat(2,x,y), D_SC.TS{1},D_SC.TS{2},'UniformOutput',false));
    D_   =  nanmax(D_(TRange,:));
    eval(sprintf('Didx = cell2mat(FreqModEvt.SpkSpkcoh_Collapsed.%s.Joint_IDs{iArea}'');',events{iEvent}));
    y = AggregateFun( D_(Didx(:,1)) , D_(Didx(:,2)) );
    scatter(x,y,20,col_{3},'filled')
    fit_ = ezfit(x,y,'affine');
    % showfit(fit_,'fitcolor',col_{3},'fitlinewidth',2,'dispfitlegend','off','dispeqboxmode','off');
    axis([0 xMax 0 yMax])
%     set(gca,'Xscale','log')
end



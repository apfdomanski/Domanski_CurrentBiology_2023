%% %%%%%% PREAMBLE %%%%%%
clear

TaskOnly = true;
plotOnline = false;
useDiscreteActTimes = true;
Epochs = {'PreSleep','Task','PostSleep'};

Epoch_ = 'Task';
Target = 'LONG';

warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
    addpath(genpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\AssemblyCode\Eleonora\Mic\programs'))
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    addpath(genpath('/panfs/panasas01/phph/ad15419/MATLAB/Mic'))
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat '*' Target '*.mat']);


% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={};%{'MiroslawLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileListAss(find(name_flag))=[];% fnames_(name_flag)=[];
    fileList(find(name_flag))=[];   % fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};

Targets = {'SHORT','MEDIUM','LONG'};
Delays_ = {'Short','Medium','Long'};
Delays = [4 8 16];

tlimsSample = [-5 5];
tlimsChoice = [-5 5];
bw = 0.05;
if ~exist([pat 'MICmetaAnalysis'])
    mkdir([pat 'MICmetaAnalysis'])
end

BinSizes=[0.005 0.01 0.0150    0.0200    0.0300    0.0500    0.0800    0.1200    0.2000    0.3500    0.5000   0.7000    1.0000];
MaxLags=[10   10   10   10   10    10    10    10    10     10     10     10     10];
% alph=0.05; %Sig level
% Dc=100;    % Bin length for variance #abba calculation
% No_th=10;  % Minimal no. activations
% O_th = Inf;   % Max assembly order
% bytelimit = 200e9; %Memory limit
% display='ordunits';%'clustered''raw'
% criteria = 'distance'; % 'biggest' or 'distance'
% act_count = 'full'; % 'full' 'partial' 'combined'
% lagChoice = 'duration'; % 'duration' 'beginning'

%% Test save
testSave = true;
if useDiscreteActTimes
    save(sprintf('%sMICmetaAnalysis%s%sAssembliesMetaTest.mat',pat,filesep,Epoch_),'testSave','-v7.3')
else
    save(sprintf('%sMICresults%sMICmetaAnalysis%s%sAssembliesMetaPartialTest.mat',pat,filesep,filesep,Epoch_),'testSave','-v7.3')
end
%% Batch process
for iFile =1:length(fileList)
     try
        fname=fileList(iFile).name;
        %% Get the spike times and events
        % Get the spike times
        fname=fileList(iFile).name;
        load(sprintf('%s%s',pat,fname),'PFCcells','HPcells');
        load(sprintf('%s%s%s%s_Events.mat',pat,'allTimestamps',filesep,strtok(fname,'.')),'t');
        
        % Wrangle data
        for i =1:length(PFCcells)
            ST{1}{i}  =  PFCcells{i}.t*1e-4;
            N_{1}(i)  = length(ST{1}{i}); % No spikes
            M_{1}(i)  = max(ST{1}{i});    % last spike time
            M__{1}(i) = min(ST{1}{i});    % first spike time
        end
        for i =1:length(HPcells)
            ST{2}{i}   =  HPcells{i}.t*1e-4;
            N_{2}(i)   = length(ST{2}{i}); % No spikes
            M_{2}(i)   = max(ST{2}{i});    % last spike time
            M__{2}(i)  = min(ST{2}{i});    % first spike time
        end
        
        ST{3} = [ST{1},ST{2}];    % Spike times
        N_{3} = [N_{1},N_{2}];    % No spikes
        M_{3} = [M_{1},M_{2}];    % Last spike
        M__{3} = [M__{1},M__{2}]; % First spike
        
        % epoch = [0, max(cell2mat(M_))];
        epochTotal = [min(cell2mat(M__)), max(cell2mat(M_))];
        
        % Restrict assembly analysis to range of task timestamps +/-30s
        times_ = [];
        
        for iDelay = 1:length(Delays_)
            eval(sprintf('temp = fieldnames(t.%s);',Delays_{iDelay}));
            for i=1:length(temp)
                eval(sprintf('times_ = [times_; t.%s.%s(:)];',Delays_{iDelay},temp{i}));
            end
        end
        
        times_ = sort(times_)*1e-6;
        epoch = [min(times_)-30, max(times_)+30];
        
        %     figure; hold on
        %     plot(epoch,[-1 -1],'*r')
        %     plot(times_,0.*times_,'.k')
        %
        %     x = sort(cell2mat(ST{3}'));
        %     plot(x,3+0.*x,'.k')
        
        clear temp times_ i
        
        Tmtx = epoch(1):0.05:epoch(2);
        
        nUnits = [length(PFCcells) length(HPcells) length(PFCcells)+length(HPcells)];
        clear i ST N_ M_ PFCcells HPcells  
        %% Batch process Assems to get activation times
        for s =1:3
            if strcmp(Epoch_,'Task')
                if TaskOnly
%                     fname_ = [pat  'MICresults' filesep 'MICresultsTask' filesep strtok(fname,'.') '_' Areas{s} '_AssembliesTask.mat'];
                    if useDiscreteActTimes
                        fname_ = [pat  'MICresults' filesep 'MICresultsTask' filesep 'UnorderedAss' filesep strtok(fname,'.') '_' Areas{s} '_AssembliesTask_Unordered.mat'];
                    else
                        fname_ = [pat  'MICresults' filesep 'MICresultsTask' filesep 'UnorderedAssPartialActivation' filesep strtok(fname,'.') '_' Areas{s} '_AssembliesTask_UnorderedPartial.mat'];
                    end
                else
                    fname_ = [pat 'MICresults' filesep 'MICresults' strtok(fname,'.') '_' Areas{s} '_Assemblies.mat'];
                end
            else
                fname_ = [pat 'MICresults' filesep 'MICresults' Epoch_ filesep strtok(fname,'.') '_' Areas{s} '_Assemblies' Epoch_ '.mat'];
            end
            Ass{s} = load(fname_);
            Ass{s}.nAss = length(Ass{s}.assembly_activity);
            
            %% Process assembly members and timing
            ClusterID = 1:nUnits(s);
            for iAss = 1:Ass{s}.nAss
                
                thisAssem =  Ass{s}.As_order(iAss);
                idx = find(~isnan(Ass{s}.Amatrix(: ,thisAssem)));
                Ass{s}.members{iAss} = Ass{s}.Unit_order(idx);
                Ass{s}.timing{iAss}  = Ass{s}.Amatrix(idx,thisAssem);
                
                Ass{s}.memberstring{iAss} = '';
                units_ = sort(Ass{s}.members{iAss});
                for imember =1:length(units_)-1
                    Ass{s}.memberstring{iAss} = strcat(Ass{s}.memberstring{iAss},[num2str(ClusterID(units_(imember))) ,', ']);
                end
                Ass{s}.memberstring{iAss} = strcat(Ass{s}.memberstring{iAss}, num2str(ClusterID(units_(end))));
                
            end
            clear idx iAss units_ thisAssem imember iAss
            
            % Correct timing variable to multiply by binwidth to indicate total duration of assembly activation
            Ass{s}.timingAdjusted = cellfun(@(timing,binwidth) (timing+1) * binwidth,...
                                        Ass{s}.timing,...
                                        num2cell(Ass{s}.Binvector(Ass{s}.As_order)),'UniformOutput',false);
                                    
            Ass{s}.Duration = (cellfun(@max,Ass{s}.timing)+1) .* Ass{s}.Binvector(Ass{s}.As_order);
            
            % Prune single-area assemblies from joint assemblies
            if s==3
                
                idx_act = [];
                idx_Matrix = [];
                for iAss =1:Ass{s}.nAss
                    thisAssem = Ass{s}.As_order(iAss);
                    members_  = find(~isnan(Ass{s}.Amatrix(:,thisAssem)));
                    Ass{s}.Unit_order(members_)
                    if ~(sum(Ass{s}.Unit_order(members_)>nUnits(1))>0 && sum(Ass{s}.Unit_order(members_)<=nUnits(1))>0)
                        idx_act    = [idx_act; iAss] ;
                        idx_Matrix = [idx_Matrix ;thisAssem] ;
                    end
                end
                
                
                Ass{s}.Amatrix(:,idx_Matrix)=[];
                Ass{s}.As_across_bins_index_pr(idx_act)=[];
                Ass{s}.As_across_bins_pr(idx_act)=[];
                            
                Ass{s}.Binvector(idx_Matrix)=[];
                Ass{s}.assembly_activity(idx_act)=[];
                Ass{s}.nAss =length(Ass{s}.Binvector);
                Ass{s}.As_order = 1:Ass{s}.nAss;
                Ass{s}.members(idx_act)=[];
                Ass{s}.memberstring(idx_act) = [];
                Ass{s}.timing(idx_act) = [];
                Ass{s}.timingAdjusted(idx_act) = [];
                Ass{s}.Duration(idx_act) = [];

               
                
%                 for i = 1:length(idx_act)
%                     thisAss = idx_act(i);
%                     Ass{s}.assembly_activity{thisAss} = nan(1,2);
%                     Ass{s}.members{thisAss} = [];
%                     Ass{s}.memberstring{thisAss} = [];
%                     Ass{s}.timing{thisAss} = NaN;
%                 end
%                 Ass{s}.Amatrix(:,idx_Matrix) = NaN;
%                 Ass{s}.Binvector(idx_Matrix)=NaN;
                
            end
%             if s==3
                
%                 idx = ~ cell2mat(cellfun(@(members,nUnits1) sum(members>nUnits1)>0 & sum(members<=nUnits1)>0 , ...
%                     Ass{s}.members,num2cell(repmat(nUnits(1),1,Ass{s}.nAss)),'UniformOutput',false));
%                 
%                 Ass{s}.members(idx) = NaN;%[];
%                 Ass{s}.memberstring(idx) = NaN% [];
%                 
%                 Ass{s}.timing(idx)  = [];
%                 Ass{s}.As_order(ismember(find(idx),Ass{s}.As_order)) = [];
%                 Ass{s}.Binvector((idx)) = [];
%                 Ass{s}.assembly_activity(idx)=[];
%                 Ass{s}.nAss = length(Ass{s}.assembly_activity);
%                 Ass{s}.Amatrix(:,idx) =[];
%             end
            
            %% Process activation times
            
            for iAss = 1:Ass{s}.nAss
                fprintf('Analysing Assem no %d/%d\n',iAss,Ass{s}.nAss)
                %% Get activation times

                    temp = Ass{s}.assembly_activity{iAss};
                    times_  = temp((find(temp(:,2)>0)),1);
                    counts_ = temp((find(temp(:,2)>0)),2);
                    counts_(times_<epoch(1) | times_>epoch(2) )=[];
                    times_ (times_<epoch(1) | times_>epoch(2) )=[];
                    if useDiscreteActTimes
                        Ass{s}.ActTimes{iAss} = repelem(times_,counts_);
                    else
                        counts_ = round(counts_*1000);
                        Ass{s}.ActTimes{iAss} = repelem(times_,counts_);
                    end
                                        
                
%                 end
           %% prune unused variables
                
                %          Ass{s}.As_across_bins{iAss}.Time=[];
                %             Ass{s}.As_across_bins_pr{iAss}.Time=[];
                %             for i = 1:length(Ass{s}.assemblies.bin)
                %                 try
                %                     Ass{s}.assemblies.bin{i}.bin_edges=[];
                %                 end
                %                 try
                %                     for j = 1:length( Ass{s}.assemblies.bin{i}.n)
                %
                %                         Ass{s}.assemblies.bin{i}.n{j}.Time= [];
                %                     end
                %                 end
                %             end
            end
            %% prune unused variables
            Ass{s} = rmfield(Ass{s},{'As_across_bins','As_across_bins_index','As_across_bins_pr','As_across_bins_index_pr','assemblies','assembly_activity'});
            clear iAss  times_ counts_ temp
            
        end
        %% calculate average activation on L and R Trials (correct and errors)
        
        
        for s =1:3
            Ass{s}.trialHists = [];
            for iDelay = 1:length(Delays_)
                Delay_ = Delays_{iDelay};
                DelayDur = Delays(iDelay);
                Cue_L= eval(sprintf('t.%s.CueLight_LeftCorrect*1e-6',Delay_));
                Cue_R= eval(sprintf('t.%s.CueLight_RightCorrect*1e-6',Delay_));
                SP_L = eval(sprintf('t.%s.SamplePress_LeftCorrect*1e-6',Delay_));
                SP_R = eval(sprintf('t.%s.SamplePress_RightCorrect*1e-6',Delay_));
                NP_L = eval(sprintf('t.%s.NosePoke_LeftCorrect*1e-6',Delay_));
                NP_R = eval(sprintf('t.%s.NosePoke_RightCorrect*1e-6',Delay_));
                CP_L = eval(sprintf('t.%s.ChoicePress_LeftCorrect*1e-6',Delay_));
                CP_R = eval(sprintf('t.%s.ChoicePress_RightCorrect*1e-6',Delay_));
                Rew_L = eval(sprintf('t.%s.RewardConsume_LeftCorrect*1e-6',Delay_));
                Rew_R = eval(sprintf('t.%s.RewardConsume_RightCorrect*1e-6',Delay_));
                
                Cue_Le= eval(sprintf('t.%s.CueLight_LeftError*1e-6',Delay_))';
                Cue_Re= eval(sprintf('t.%s.CueLight_RightError*1e-6',Delay_))';
                SP_Le = eval(sprintf('t.%s.SamplePress_LeftError*1e-6',Delay_))';
                SP_Re = eval(sprintf('t.%s.SamplePress_RightError*1e-6',Delay_))';
                NP_Le = eval(sprintf('t.%s.NosePoke_LeftError*1e-6',Delay_))';
                NP_Re = eval(sprintf('t.%s.NosePoke_RightError*1e-6',Delay_))';
                CP_Le = eval(sprintf('t.%s.ChoicePress_LeftError*1e-6',Delay_))';
                CP_Re = eval(sprintf('t.%s.ChoicePress_RightError*1e-6',Delay_))';
%                 if useDiscreteActTimes
                    for iAss = 1:Ass{s}.nAss
                        [s iAss iDelay]
                        
                        LC=[];RC=[];LCe=[];RCe=[];
                        for iTrial = 1:length(SP_L)
                            %@TODO: should the binwidth change to match the assembly's own
                            tbS = [tlimsSample(1)+(SP_L(iTrial)):bw:tlimsSample(2)+(SP_L(iTrial))+bw];
                            tbC = [tlimsChoice(1)+(CP_L(iTrial)):bw:tlimsChoice(2)+(CP_L(iTrial))+bw];
                            tempS = histc(Ass{s}.ActTimes{iAss},tbS);
                            tempC = histc(Ass{s}.ActTimes{iAss},tbC);
                            if ~useDiscreteActTimes
                                tempS = tempS/1000;
                                tempC = tempC/1000;
                            end
                            LC(:,iTrial) = [tempS(1:end-1);tempC(1:end-1)];
                        end
                        tempS = [];tempC = [];
                        
                        for iTrial = 1:length(SP_R)
                            %@TODO: should the binwidth change to match the assembly's own
                            tbS = [tlimsSample(1)+(SP_R(iTrial)):bw:tlimsSample(2)+(SP_R(iTrial))+bw];
                            tbC = [tlimsChoice(1)+(CP_R(iTrial)):bw:tlimsChoice(2)+(CP_R(iTrial))+bw];
                            tempS = histc(Ass{s}.ActTimes{iAss},tbS);
                            tempC = histc(Ass{s}.ActTimes{iAss},tbC);
                            if ~useDiscreteActTimes
                                tempS = tempS/1000;
                                tempC = tempC/1000;
                            end
                            RC(:,iTrial) = [tempS(1:end-1);tempC(1:end-1)];
                        end
                        tempS = [];tempC = [];
                        
                        for iTrial = 1:length(SP_Le)
                            %@TODO: should the binwidth change to match the assembly's own
                            try
                                tbS = [tlimsSample(1)+(SP_Le(iTrial)):bw:tlimsSample(2)+(SP_Le(iTrial))+bw];
                                tbC = [tlimsChoice(1)+(CP_Le(iTrial)):bw:tlimsChoice(2)+(CP_Le(iTrial))+bw];
                                tempS = histc(Ass{s}.ActTimes{iAss},tbS);
                                tempC = histc(Ass{s}.ActTimes{iAss},tbC);
                                if ~useDiscreteActTimes
                                    tempS = tempS/1000;
                                    tempC = tempC/1000;
                                end
                                LCe(:,iTrial) = [tempS(1:end-1);tempC(1:end-1)];
                            catch
                                LCe(:,iTrial) = nan(size(LC,1),1);
                            end
                        end
                        tempS = [];tempC = [];
                        
                        for iTrial = 1:length(SP_Re)
                            %@TODO: should the binwidth change to match the assembly's own
                            try
                                tbS = [tlimsSample(1)+(SP_Re(iTrial)):bw:tlimsSample(2)+(SP_Re(iTrial))+bw];
                                tbC = [tlimsChoice(1)+(CP_Re(iTrial)):bw:tlimsChoice(2)+(CP_Re(iTrial))+bw];
                                tempS = histc(Ass{s}.ActTimes{iAss},tbS);
                                tempC = histc(Ass{s}.ActTimes{iAss},tbC);
                                if ~useDiscreteActTimes
                                    tempS = tempS/1000;
                                    tempC = tempC/1000;
                                end
                                RCe(:,iTrial) = [tempS(1:end-1);tempC(1:end-1)];
                            catch
                                RCe(:,iTrial) = nan(size(LC,1),1);
                            end
                        end
                        tempS = [];tempC = [];
                        
                        Ass{s}.trialHists.LC{iDelay}{iAss} = LC;
                        Ass{s}.trialHists.RC{iDelay}{iAss} = RC;
                        Ass{s}.trialHists.LCe{iDelay}{iAss} = LCe;
                        Ass{s}.trialHists.RCe{iDelay}{iAss} = RCe;
                        
                        Ass{s}.trialHists.LC_mean{iDelay}(iAss,:)  =  nanmean(LC,2);
                        Ass{s}.trialHists.RC_mean{iDelay}(iAss,:)  =  nanmean(RC,2);
                        Ass{s}.trialHists.LCe_mean{iDelay}(iAss,:) =  nanmean(LCe,2);
                        Ass{s}.trialHists.RCe_mean{iDelay}(iAss,:) =  nanmean(RCe,2);
                        
                        Ass{s}.trialHists.LC_sem{iDelay}(iAss,:)  =  nansem(LC,2);
                        Ass{s}.trialHists.RC_sem{iDelay}(iAss,:)  =  nansem(RC,2);
                        Ass{s}.trialHists.LCe_sem{iDelay}(iAss,:) =  nansem(LCe,2);
                        Ass{s}.trialHists.RCe_sem{iDelay}(iAss,:) =  nansem(RCe,2);                        
                        
                        %Delay activations
                        LD=[];RD=[]; LDe=[];RDe=[]; 
                        for iTrial = 1:length(SP_L)
                            %@TODO: should the binwidth change to match the assembly's own
                            try
                                tbD = tlimsSample(1)+SP_L(iTrial):bw:SP_L(iTrial)+DelayDur+bw;
                                tempD = histc(Ass{s}.ActTimes{iAss},tbD);
                                if ~useDiscreteActTimes
                                    tempD = tempD/1000;
                                end
                                LD(:,iTrial) = tempD(1:end-1);
                            catch
                                LD(:,iTrial) = nan(size(tbD,1),1);
                            end
                        end
                        tempD = [];
                        for iTrial = 1:length(SP_R)
                            %@TODO: should the binwidth change to match the assembly's own
                            try
                                tbD = tlimsSample(1)+SP_R(iTrial):bw:SP_R(iTrial)+DelayDur+bw;
                                tempD = histc(Ass{s}.ActTimes{iAss},tbD);
                                if ~useDiscreteActTimes
                                    tempD = tempD/1000;
                                end
                                RD(:,iTrial) = tempD(1:end-1);
                            catch
                                RD(:,iTrial) = nan(size(tbD,1),1);
                            end
                        end
                        tempD = [];
                        for iTrial = 1:length(SP_Le)
                            %@TODO: should the binwidth change to match the assembly's own
                            try
                                tbD = tlimsSample(1)+SP_Le(iTrial):bw:SP_Le(iTrial)+DelayDur+bw;
                                tempD = histc(Ass{s}.ActTimes{iAss},tbD);
                                if ~useDiscreteActTimes
                                    tempD = tempD/1000;
                                end
                                LDe(:,iTrial) = tempD(1:end-1);
                            catch
                                LDe(:,iTrial) = nan(size(tbD,1),1);
                            end
                        end
                        tempD = [];
                        for iTrial = 1:length(SP_Re)
                            %@TODO: should the binwidth change to match the assembly's own
                            try
                                tbD = tlimsSample(1)+SP_Re(iTrial):bw:SP_Re(iTrial)+DelayDur+bw;
                                tempD = histc(Ass{s}.ActTimes{iAss},tbD);
                                if ~useDiscreteActTimes
                                    tempD = tempD/1000;
                                end
                                RDe(:,iTrial) = tempD(1:end-1);
                            catch
                                RDe(:,iTrial) = nan(size(tbD,1),1);
                            end
                        end
                        tempD = [];

                        Ass{s}.trialHists.LD{iDelay}{iAss} = LD;
                        Ass{s}.trialHists.RD{iDelay}{iAss} = RD;
                        Ass{s}.trialHists.LDe{iDelay}{iAss} = LDe;
                        Ass{s}.trialHists.RDe{iDelay}{iAss} = RDe;
                        
                        Ass{s}.trialHists.LD_mean{iDelay}(iAss,:)  =  nanmean(LD,2);
                        Ass{s}.trialHists.RD_mean{iDelay}(iAss,:)  =  nanmean(RD,2);
                        Ass{s}.trialHists.LDe_mean{iDelay}(iAss,:) =  nanmean(LDe,2);
                        Ass{s}.trialHists.RDe_mean{iDelay}(iAss,:) =  nanmean(RDe,2);
                        
                        Ass{s}.trialHists.LD_sem{iDelay}(iAss,:)  =  nansem(LD,2);
                        Ass{s}.trialHists.RD_sem{iDelay}(iAss,:)  =  nansem(RD,2);
                        Ass{s}.trialHists.LDe_sem{iDelay}(iAss,:) =  nansem(LDe,2);
                        Ass{s}.trialHists.RDe_sem{iDelay}(iAss,:) =  nansem(RDe,2);                          

                        %               ciplot(smooth_hist(nanmean(LC,2)+nansem(LC,2)),smooth_hist(nanmean(LC,2)-nansem(LC,2)),(1:size(LC,1))*bw,[0.5 0 0.5])
                        %               ciplot(smooth_hist(nanmean(RC,2)+nansem(RC,2)),smooth_hist(nanmean(RC,2)-nansem(RC,2)),(1:size(RC,1))*bw,'b')
                        % %               plot((1:size(LCe,1))*bw,nanmean(LCe,2),'color',[0.5 0 0.5],'LineStyle',':')
                        % %               plot((1:size(RCe,1))*bw,nanmean(RCe,2),'color',[0.5 0 0.5],'LineStyle',':')
                        %
                    end
%                 else
%                     bins_ = unique(Ass{s}.Binvector);
%                     for iBin=1:length(bins_)
%                         AssList = find(Ass{s}.Binvector==bins_(iBin));
%                         times_ = Ass{s}.assembly_activity{AssList(1)}(:,1);
%                         upsamplefactor = bins_(iBin)/0.05;
%                         %%%
%                         iAss = 1;
%                         actn_ = Ass{s}.assembly_activity{AssList(iAss)}(:,2);
%                         %%%
%                         actnResampled = interp1(times_,actn_,min(times_):0.05:max(times_),'cubic');
% %                         actnResampled = resample(actn_,bins_(iBin)*1000,0.05*1000,1000);
% %                         timesResampled = resample(times_,bins_(iBin)*1000,0.05*1000,1000);
%                         timesResampled = min(times_):0.05:max(times_);
%                         figure; hold on;
%                         plot(times_,actn_);
%                         plot(timesResampled,actnResampled);
%                         actnResampled  = repmat(actn_,[upsamplefactor,1]);%y=y(:)';
%                         
%                         N = [length(SP_L);length(SP_R);length(SP_Le);length(SP_Re)];
%                         timeLookup=[];
%                         for iTrial = 1:length(SP_L)
%                             tbS = [tlimsSample(1)+(SP_L(iTrial)):bw:tlimsSample(2)+(SP_L(iTrial))+bw];
%                             tbC = [tlimsChoice(1)+(CP_L(iTrial)):bw:tlimsChoice(2)+(CP_L(iTrial))+bw];
%                             timeLookup=[timeLookup,tbS,tbC];
%                             
%                         end
%                         [~,Idx_LC] = min(abs(bsxfun(@minus,timeLookup, times_)));
% interp1
%                         for iTrial = 1:length(SP_R)
%                             tbS = [tlimsSample(1)+(SP_R(iTrial)):bw:tlimsSample(2)+(SP_R(iTrial))+bw];
%                             tbC = [tlimsChoice(1)+(CP_R(iTrial)):bw:tlimsChoice(2)+(CP_R(iTrial))+bw];
%                             timeLookup=[timeLookup,tbS,tbC];
%                         end
%                         
%                         for iTrial = 1:length(SP_Le)
%                             tbS = [tlimsSample(1)+(SP_Le(iTrial)):bw:tlimsSample(2)+(SP_Le(iTrial))+bw];
%                             tbC = [tlimsChoice(1)+(CP_Le(iTrial)):bw:tlimsChoice(2)+(CP_Le(iTrial))+bw];
%                             timeLookup=[timeLookup,tbS,tbC];
%                         end
%                         
%                         for iTrial = 1:length(SP_Re)
%                             tbS = [tlimsSample(1)+(SP_Re(iTrial)):bw:tlimsSample(2)+(SP_Re(iTrial))+bw];
%                             tbC = [tlimsChoice(1)+(CP_Re(iTrial)):bw:tlimsChoice(2)+(CP_Re(iTrial))+bw];
%                             timeLookup=[timeLookup,tbS,tbC];
%                         end
% %                     times_ = Ass{s}.assembly_activity{iAss}(:,1);
% %                         actn_  = Ass{s}.assembly_activity{iAss}(:,2);
%                     [~,closestIndex] = min(abs(bsxfun(@minus,timeLookup, times_)));
%                     LC(:,iTrial) = actn_(closestIndex);
%                 end
%             end
            
            end
        end
        tb = (1:size(LC,1))*bw;
        for s = 1:3
            Ass{s}.trialHists.tb = tb;
        end
        if plotOnline
            for s =1:3
                figure('name',Areas{s})
                for iDelay = 1:length(Delays_)
                    subplot(1,3,iDelay); hold on
                    z =(Ass{s}.trialHists.LC_mean{iDelay})';
                    %            smooth2a(z,0,10)
                    z=zscore(z);
                    %            if iDelay ==1
                    [~,I] = max(z);
                    [~,I] = sort(I,'descend');
                    N= max(I);
                    %            end
                    imagesc(tb,1:N,z(:,I)')
                    %            imagesc(tb,1:N,zscore(z(:,I))')
                    colormap hot
                    %         colormap(cmap)
                    %         caxis([-1 1])
                    title(Delays_{iDelay})
                    
                    % Vertical bars
                    plot(sum(abs(tlimsSample))/2*[1,1],[0 N+1],'color',[0 1 0 0.3],'LineWidth',2)
                    plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0 N+1],'color',[1 0 0 0.3],'LineWidth',2)
                    plot(sum(abs(tlimsSample))*[1,1],[0 N+1],'k')
                    axis([min(tb) max(tb) 0 max(I)])
                end
            end
        end
        %% Decode L vs R whole trial
        for s = 1:3
            for iDelay =1:length(Delays_)
                x = cellfun(@(a) smooth2a(a,10,0),Ass{s}.trialHists.LC{iDelay},'UniformOutput',false);
                Ltrials = cell2mat(cellfun(@(a) a(:),x,'UniformOutput',false));
                x = cellfun(@(a) smooth2a(a,10,0),Ass{s}.trialHists.RC{iDelay},'UniformOutput',false);
                Rtrials = cell2mat(cellfun(@(a) a(:),x,'UniformOutput',false));
                
                %         Ltrials = cell2mat(cellfun(@(a) a(:),Ass{s}.trialHists.LC{iDelay},'UniformOutput',false));
                %         Rtrials = cell2mat(cellfun(@(a) a(:),Ass{s}.trialHists.RC{iDelay},'UniformOutput',false));
                %                     Ltrials = cell2mat(cellfun(@(a) a(:),cellfun(@transpose,Ass{s}.trialHists.LC{iDelay},'UniformOutput',false),'UniformOutput',false));
                %                     Rtrials = cell2mat(cellfun(@(a) a(:),cellfun(@transpose,Ass{s}.trialHists.RC{iDelay},'UniformOutput',false),'UniformOutput',false));
                nL = size(Ass{s}.trialHists.LC{iDelay}{1},2);
                nR = size(Ass{s}.trialHists.RC{iDelay}{1},2);
                
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = size(Ass{s}.trialHists.LC{iDelay}{1},1);
                
                % Multivariate F-score decoder
                [D_.avgFR,D_.seFR,...
                    D_.Ft2,D_.Rt2,...
                    D_.Ft2ciL,D_.Ft2ciH,...
                    D_.TS,...
                    D_.dfnum,D_.dfden] = DecodeStats(FR,evt0,0.05);
                D_.TSsig = tpdf(D_.TS,nL+nR-2)<0.05;
                D_.TS(isnan(D_.TS))=0;
                
                %D_.CVE = DecodeCVE_mex(FR,evt0,0.05);
                Ass{s}.trialHists.D_LR{iDelay}=D_;
                clear D_
            end
        end
        
        if plotOnline
            col_ = gray(4);
            for s = 1:3
                figure('name',Areas{s})
                for iDelay =1:length(Delays_)
                    D_ = Ass{s}.trialHists.D_LR{iDelay};
                    subplot(1,4,iDelay);hold on
                    title(Delays_{iDelay})
                    % Vertical bars
                    plot(sum(abs(tlimsSample))/2*[1,1],[0 N+1],'color',[0 1 0 0.3],'LineWidth',2,'HandleVisibility','off')
                    plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0 N+1],'color',[1 0 0 0.3],'LineWidth',2,'HandleVisibility','off')
                    plot(sum(abs(tlimsSample))*[1,1],[0 N+1],'k','HandleVisibility','off')
                    axis([min(tb) max(tb) 0 max(I)])
                    
                    
                    x = (1:Ltr).*bw;
                    y = 1:Ass{s}.nAss;
                    z = D_.TS;
                    z_= D_.TSsig;
                    
                    z = zscore(z);
                    
                    %         if iDelay ==1
                    [~,I] = max(z);
                    [~,I] = sort(I,'descend');
                    %         end
                    z=z(:,I)';
                    z_=z_(:,I)';
                    z(~z_)=0;
                    imagesc(x,y,z)
                    % Vertical bars
                    plot(sum(abs(tlimsSample))/2*[1,1],[0 N+1],'color',[0 1 0 0.3],'LineWidth',2,'HandleVisibility','off')
                    plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0 N+1],'color',[1 0 0 0.3],'LineWidth',2,'HandleVisibility','off')
                    plot(sum(abs(tlimsSample))*[1,1],[0 N+1],'k','HandleVisibility','off')
                    axis([min(tb) max(tb) 0 max(I)])
                    colormap([1 1 1; jet])
                    caxis([0 5])
                    subplot(1,4,4);hold on
                    plot((1:Ltr).*bw,D_.Ft2,'color',col_(iDelay,:),'LineWidth',1.5)
                    plot((1:Ltr).*bw,D_.Ft2ciH,'color',col_(iDelay,:),'lineStyle',':','HandleVisibility','off')
                    %         plot((1:Ltr).*bw,1-D_.CVE,'color',col_(iDelay,:))
                    
                    % Vertical bars
                    plot(sum(abs(tlimsSample))/2*[1,1],[0 N+1],'color',[0 1 0 0.3],'LineWidth',2,'HandleVisibility','off')
                    plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0 N+1],'color',[1 0 0 0.3],'LineWidth',2,'HandleVisibility','off')
                    plot(sum(abs(tlimsSample))*[1,1],[0 N+1],'k','HandleVisibility','off')
                    %         axis([min(tb) max(tb) 0 1])
                    axis([min(tb) max(tb) 0 10])
                    
                    
                end
                subplot(1,4,4)
                legend(Delays_)
            end
        end
        %% Decode L vs R delay only
          for s = 1:3
            for iDelay =1:length(Delays_)
                x = cellfun(@(a) smooth2a(a,10,0),Ass{s}.trialHists.LD{iDelay},'UniformOutput',false);
                Ltrials = cell2mat(cellfun(@(a) a(:),x,'UniformOutput',false));
                x = cellfun(@(a) smooth2a(a,10,0),Ass{s}.trialHists.RD{iDelay},'UniformOutput',false);
                Rtrials = cell2mat(cellfun(@(a) a(:),x,'UniformOutput',false));
                
                %         Ltrials = cell2mat(cellfun(@(a) a(:),Ass{s}.trialHists.LC{iDelay},'UniformOutput',false));
                %         Rtrials = cell2mat(cellfun(@(a) a(:),Ass{s}.trialHists.RC{iDelay},'UniformOutput',false));
                %                     Ltrials = cell2mat(cellfun(@(a) a(:),cellfun(@transpose,Ass{s}.trialHists.LC{iDelay},'UniformOutput',false),'UniformOutput',false));
                %                     Rtrials = cell2mat(cellfun(@(a) a(:),cellfun(@transpose,Ass{s}.trialHists.RC{iDelay},'UniformOutput',false),'UniformOutput',false));
                nL = size(Ass{s}.trialHists.LD{iDelay}{1},2);
                nR = size(Ass{s}.trialHists.RD{iDelay}{1},2);
                
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = size(Ass{s}.trialHists.LC{iDelay}{1},1);
                
                % Multivariate F-score decoder
                [D_.avgFR,D_.seFR,...
                    D_.Ft2,D_.Rt2,...
                    D_.Ft2ciL,D_.Ft2ciH,...
                    D_.TS,...
                    D_.dfnum,D_.dfden] = DecodeStats(FR,evt0,0.05);
                D_.TSsig = tpdf(D_.TS,nL+nR-2)<0.05;
                D_.TS(isnan(D_.TS))=0;
                
                %D_.CVE = DecodeCVE_mex(FR,evt0,0.05);
                Ass{s}.trialHists.D_Delay{iDelay}=D_;
                clear D_
            end
          end
         if plotOnline
            col_ = gray(4);
            for s = 1:3
                figure('name',Areas{s})
                for iDelay =1:length(Delays_)
                    D_ = Ass{s}.trialHists.D_Delay{iDelay};
                    N = size(D_.TS,2);
                    subplot(1,4,iDelay);hold on
                    title(Delays_{iDelay})
                   
                    
                    Ltr_ = size(D_.TS,1);
                    
                    x = (1:Ltr_).*bw;
                    y = 1:Ass{s}.nAss;
                    z = D_.TS;
                    z_= D_.TSsig;
                    
                    z = zscore(z);
                    
                    %         if iDelay ==1
                    [~,I] = max(z);
                    [~,I] = sort(I,'descend');
                    %         end
                    z=z(:,I)';
                    z_=z_(:,I)';
                    z(~z_)=0;
                    imagesc(x,y,z)
                    % Vertical bars
                    axis([min(tb) max(tb) 0 max(I)])
                    colormap([1 1 1; jet])
                    caxis([0 5])
                    subplot(1,4,4);hold on
                    plot((1:Ltr_).*bw,D_.Ft2,'color',col_(iDelay,:),'LineWidth',1.5)
                    plot((1:Ltr_).*bw,D_.Ft2ciH,'color',col_(iDelay,:),'lineStyle',':','HandleVisibility','off')
                    %         plot((1:Ltr).*bw,1-D_.CVE,'color',col_(iDelay,:))
                    plot(sum(abs(tlimsSample))/2*[1,1],[0 N+1],'k','HandleVisibility','off')


                    axis([min(tb) 20 0 10])
                    
                    
                end
                subplot(1,4,4)
                legend(Delays_)
            end
        end
       
        %% Nosepoke-triggered assem rates
        for s =1:3
            %           figure
            for iDelay = 1:length(Delays_)
                Delay_ = Delays_{iDelay};
                
                NP_L = eval(sprintf('t.%s.NosePoke_LeftCorrect*1e-6',Delay_));
                NP_R = eval(sprintf('t.%s.NosePoke_RightCorrect*1e-6',Delay_));
                
                NP_Le = eval(sprintf('t.%s.NosePoke_LeftError*1e-6',Delay_))';
                NP_Re = eval(sprintf('t.%s.NosePoke_RightError*1e-6',Delay_))';
                
                for iAss = 1:Ass{s}.nAss
                    [s iAss iDelay]
                    
                    
                    LP=[];RP=[];LPe=[];RPe=[];
                    
                    for iTrial = 1:length(NP_L)
                        %@TODO: should the binwidth change to match the assembly's own
                        tbNP = [tlimsSample(1)+(NP_L(iTrial)):bw:tlimsSample(2)+(NP_L(iTrial))+bw];
                        tempNP = histc(Ass{s}.ActTimes{iAss},tbNP);
                        if ~useDiscreteActTimes
                            tempNP = tempNP/1000;
                        end
                    LP(:,iTrial) = [tempNP(1:end-1)];
                    end
                    tempNP = [];
                    
                    for iTrial = 1:length(NP_R)
                        %@TODO: should the binwidth change to match the assembly's own
                        tbNP = [tlimsSample(1)+(NP_R(iTrial)):bw:tlimsSample(2)+(NP_R(iTrial))+bw];
                        tempNP = histc(Ass{s}.ActTimes{iAss},tbNP);
                        if ~useDiscreteActTimes
                            tempNP = tempNP/1000;
                        end
                        RP(:,iTrial) = [tempNP(1:end-1)];
                    end
                    tempNP = [];
                    
                    for iTrial = 1:length(NP_Le)
                        try
                            %@TODO: should the binwidth change to match the assembly's own
                            tbNP = [tlimsSample(1)+(NP_Le(iTrial)):bw:tlimsSample(2)+(NP_Le(iTrial))+bw];
                            tempNP = histc(Ass{s}.ActTimes{iAss},tbNP);
                            if ~useDiscreteActTimes
                                tempNP = tempNP/1000;
                            end
                            LPe(:,iTrial) = [tempNP(1:end-1)];
                        catch
                            LPe(:,iTrial) = nan(size(LP,1),1);
                        end
                    end
                    tempNP = [];
                    
                    for iTrial = 1:length(NP_Re)
                        try
                            %@TODO: should the binwidth change to match the assembly's own
                            tbNP = [tlimsSample(1)+(NP_Re(iTrial)):bw:tlimsSample(2)+(NP_Re(iTrial))+bw];
                            tempNP = histc(Ass{s}.ActTimes{iAss},tbNP);
                            if ~useDiscreteActTimes
                                tempNP = tempNP/1000;
                            end
                            RPe(:,iTrial) = [tempNP(1:end-1)];
                        catch
                            RPe(:,iTrial) = nan(size(RP,1),1);
                        end
                    end
                    tempNP = [];
                    
                    Ass{s}.trialHists.LP{iDelay}{iAss} = LP;
                    Ass{s}.trialHists.RP{iDelay}{iAss} = RP;
                    Ass{s}.trialHists.LPe{iDelay}{iAss} = LPe;
                    Ass{s}.trialHists.RPe{iDelay}{iAss} = RPe;
                    
                    Ass{s}.trialHists.LP_mean{iDelay}(iAss,:)  =  nanmean(LP,2);
                    Ass{s}.trialHists.RP_mean{iDelay}(iAss,:)  =  nanmean(RP,2);
                    Ass{s}.trialHists.LPe_mean{iDelay}(iAss,:) =  nanmean(LPe,2);
                    Ass{s}.trialHists.RPe_mean{iDelay}(iAss,:) =  nanmean(RPe,2);
                    
                    Ass{s}.trialHists.LP_sem{iDelay}(iAss,:)  =  nansem(LP,2);
                    Ass{s}.trialHists.RP_sem{iDelay}(iAss,:)  =  nansem(RP,2);
                    Ass{s}.trialHists.LPe_sem{iDelay}(iAss,:) =  nansem(LPe,2);
                    Ass{s}.trialHists.RPe_sem{iDelay}(iAss,:) =  nansem(RPe,2);
                    %               ciplot(smooth_hist(nanmean(LC,2)+nansem(LC,2)),smooth_hist(nanmean(LC,2)-nansem(LC,2)),(1:size(LC,1))*bw,[0.5 0 0.5])
                    %               ciplot(smooth_hist(nanmean(RC,2)+nansem(RC,2)),smooth_hist(nanmean(RC,2)-nansem(RC,2)),(1:size(RC,1))*bw,'b')
                    % %               plot((1:size(LCe,1))*bw,nanmean(LCe,2),'color',[0.5 0 0.5],'LineStyle',':')
                    % %               plot((1:size(RCe,1))*bw,nanmean(RCe,2),'color',[0.5 0 0.5],'LineStyle',':')
                    %
                end
                
            end
        end
        tb = (1:size(LP,1))*bw;
        
        if plotOnline
            for s =1:3
                figure
                for iDelay = 1:length(Delays_)
                    subplot(1,3,iDelay); hold on
                    z =((Ass{s}.trialHists.LP_mean{iDelay})' + (Ass{s}.trialHists.RP_mean{iDelay})')./2;
                    %            smooth2a(z,0,10)
                    %            z=zscore(z);
                    %            if iDelay ==1
                    [~,I] = max(z);
                    [~,I] = sort(I,'descend');
                    N= max(I);
                    %            end
                    imagesc(tb,1:N,z(:,I)')
                    %            imagesc(tb,1:N,zscore(z(:,I))')
                    colormap jet
                    %         colormap(cmap)
                    %         caxis([-1 1])
                    
                    % Vertical bars
                    plot(sum(abs(tlimsSample))/2*[1,1],[0 N+1],'color',[0 1 0 0.3],'LineWidth',2)
                    plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0 N+1],'color',[1 0 0 0.3],'LineWidth',2)
                    plot(sum(abs(tlimsSample))*[1,1],[0 N+1],'k')
                    axis([min(tb) max(tb) 0 max(I)])
                end
            end
        end
        %% Decode L vs R  - Nosepokes
        for s = 1:3
            for iDelay =1:length(Delays_)
                x = cellfun(@(a) smooth2a(a,10,0),Ass{s}.trialHists.LP{iDelay},'UniformOutput',false);
                Ltrials = cell2mat(cellfun(@(a) a(:),x,'UniformOutput',false));
                x = cellfun(@(a) smooth2a(a,10,0),Ass{s}.trialHists.RP{iDelay},'UniformOutput',false);
                Rtrials = cell2mat(cellfun(@(a) a(:),x,'UniformOutput',false));
                
                %         Ltrials = cell2mat(cellfun(@(a) a(:),Ass{s}.trialHists.LC{iDelay},'UniformOutput',false));
                %         Rtrials = cell2mat(cellfun(@(a) a(:),Ass{s}.trialHists.RC{iDelay},'UniformOutput',false));
                %                     Ltrials = cell2mat(cellfun(@(a) a(:),cellfun(@transpose,Ass{s}.trialHists.LC{iDelay},'UniformOutput',false),'UniformOutput',false));
                %                     Rtrials = cell2mat(cellfun(@(a) a(:),cellfun(@transpose,Ass{s}.trialHists.RC{iDelay},'UniformOutput',false),'UniformOutput',false));
                nL = size(Ass{s}.trialHists.LP{iDelay}{1},2);
                nR = size(Ass{s}.trialHists.RP{iDelay}{1},2);
                
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = size(Ass{s}.trialHists.LP{iDelay}{1},1);
                
                % Multivariate F-score decoder
                [D_.avgFR,D_.seFR,...
                    D_.Ft2,D_.Rt2,...
                    D_.Ft2ciL,D_.Ft2ciH,...
                    D_.TS,...
                    D_.dfnum,D_.dfden] = DecodeStats(FR,evt0,0.05);
                D_.TSsig = tpdf(D_.TS,nL+nR-2)<0.05;
                D_.TS(isnan(D_.TS))=0;
                
                %         D_.CVE = DecodeCVE_mex(FR,evt0,0.05);
                Ass{s}.trialHists.D_NP{iDelay} = D_;
            end
        end
        
        if plotOnline
            col_ = gray(4);
            for s = 1:3
                figure('name',Areas{s})
                for iDelay =1:length(Delays_)
                    D_ = Ass{s}.trialHists.D_NP{iDelay};

                    subplot(1,4,iDelay);hold on
                    % Vertical bars
                    plot([0 0],[0 N+1],'k','HandleVisibility','off')
                    axis([min(tb) max(tb) 0 max(I)])
                    
                    
                    x = (1:Ltr).*bw + tlimsSample(1);
                    y = 1:Ass{s}.nAss;
                    z = D_.TS;
                    z_= D_.TSsig;
                    %         z = zscore(z);
                    
                            if iDelay ==1
                    [~,I] = max(z);
                    [~,I] = sort(I,'descend');
                            end
                    z=z(:,I)';
                    z_=z_(:,I)';
                    z(~z_)=0;
                    imagesc(x,y,z)
                    plot([0 0],[0 N+1],'k','HandleVisibility','off')
                    axis([min(x) max(x) 0 max(I)])
                    colormap([1 1 1; jet])
                    title(Delays_{iDelay})
                    
                    subplot(1,4,4);hold on
                    plot(x,D_.Ft2,'color',col_(iDelay,:),'Linewidth',1.5)
                    plot(x,D_.Ft2ciH,'color',col_(iDelay,:),'lineStyle',':','HandleVisibility','off')
                    %  plot(x,1-D_.CVE,'color',color_{iDelay})
                    plot([0 0],[0 1],'k','HandleVisibility','off')
                    axis([min(x) max(x) 0 10])
                    legend(Delays_)
                    
                end
            end
        end        
        %% Assembly diagnostics
         if plotOnline
                figure;   
         end
        for s=1:3
            Ass{s}.membercountHist = histc(cellfun(@length,Ass{s}.members),2:20)./Ass{s}.nAss;                                          % How many member units @TODO: normalize by no. unit
            Ass{s}.membercountHistAjusted  = histc(cellfun(@length,Ass{s}.members)/nUnits(s),0:0.05:0.5)./Ass{s}.nAss;                  % What fraction of recorded units counts make up assemblies
            Ass{s}.membersharingHist = histc(histc(cell2mat(Ass{s}.members'),1:nUnits(s)),1:20)./Ass{s}.nAss;                           % How many assemblies are units shared by
            Ass{s}.membersharingHistAjusted = histc(histc(cell2mat(Ass{s}.members'),1:nUnits(s))./Ass{s}.nAss,0:0.05:0.5)./Ass{s}.nAss; % How many assemblies are units shared by
            
            Ass{s}.TimeBinHist = histc(Ass{s}.Binvector,BinSizes)./Ass{s}.nAss;
            Ass{s}.DurationHist = histc(Ass{s}.Duration,logspace(-1.3,1.3,10))./Ass{s}.nAss; % 10 log bins from 0.05 to 20s
            if plotOnline
                subplot(1,3,1); hold on
                stairs(BinSizes,Ass{s}.TimeBinHist,'Color',color_{s},'LineWidth',1.5)
                set(gca,'Xscale','log','Xlim',[min(BinSizes) max(BinSizes)])
                xlabel('Interaction timescale (s)')
                ylabel('Fraction of Assemblies')
                subplot(1,3,2); hold on
                stairs(2:20,Ass{s}.membercountHist,'Color',color_{s},'LineWidth',1.5)
                xlabel('No. interacting units')
                subplot(1,3,3); hold on
                stairs(1:20,Ass{s}.membersharingHist,'Color',color_{s},'LineWidth',1.5)
                xlabel({'Promiscuity of units';'between assemblies'})
                legend(Areas)
            end
            
            Ass_{s}.membercountHist(:,iFile)            = Ass{s}.membercountHist;
            Ass_{s}.membercountHistAjusted(:,iFile)     = Ass{s}.membercountHistAjusted;
            Ass_{s}.membersharingHist(:,iFile)          = Ass{s}.membersharingHist;
            Ass_{s}.membersharingHistAjusted(:,iFile)   = Ass{s}.membersharingHistAjusted;
            Ass_{s}.TimeBinHist(:,iFile)                = Ass{s}.TimeBinHist;
            Ass_{s}.DurationHist(:,iFile)               = Ass{s}.DurationHist;
            Ass_{s}.Duration{iFile}                     = Ass{s}.Duration;
            Ass_{s}.noMembers{iFile}                    = cellfun(@length,Ass{s}.members);
            
            [~,I] = sort(Ass{s}.As_order);
            Ass_{s}.BinVector{iFile}                    = Ass{s}.Binvector(I);
            Ass_{s}.BinVector_{iFile}                   = Ass{s}.Binvector(Ass{s}.As_order);
            
            Ass_{s}.actTimes{iFile}     = Ass{s}.ActTimes;
            Ass_{s}.trialHists{iFile}   = Ass{s}.trialHists;
            %         catch
            %             Ass_{s}.membercountHist(:,iFile)   = nan(11,1);
            %             Ass_{s}.membersharingHist(:,iFile) = nan(11,1);
            %             Ass_{s}.TimeBinHist(:,iFile)       = nan(length(BinSizes),1);
            %
        end
        %%
        fname
        if useDiscreteActTimes
            save(sprintf('%sMICmetaAnalysis%s%s_%sAssemblies.mat',pat,filesep,strtok(fname,'.'),Epoch_),'Ass','-v7.3')
        else
            save(sprintf('%sMICresults%sMICmetaAnalysis%s%s_%sAssembliesPartial.mat',pat,filesep,filesep,strtok(fname,'.'),Epoch_),'Ass','-v7.3')
        end
        %%
        clear Ass
     end
end
%% Save metadata
if useDiscreteActTimes
    save(sprintf('%sMICmetaAnalysis%s%sAssembliesMeta.mat',pat,filesep,Epoch_),'Ass_','-v7.3')
else
    save(sprintf('%sMICresults%sMICmetaAnalysis%s%sAssembliesMetaPartial.mat',pat,filesep,filesep,Epoch_),'Ass_','-v7.3')
end
%% Posthoc decoding - delay
for iFile =1:length(fileList)
     try
          for s = 1:3
            for iDelay =1:length(Delays_)
                [iFile s iDelay]
                x = cellfun(@(a) smooth2a(a,10,0),Ass_{s}.trialHists{iFile}.LD{iDelay},'UniformOutput',false);
                Ltrials = cell2mat(cellfun(@(a) a(:),x,'UniformOutput',false));
                x = cellfun(@(a) smooth2a(a,10,0),Ass_{s}.trialHists{iFile}.RD{iDelay},'UniformOutput',false);
                Rtrials = cell2mat(cellfun(@(a) a(:),x,'UniformOutput',false));
                
                %         Ltrials = cell2mat(cellfun(@(a) a(:),Ass{s}.trialHists.LC{iDelay},'UniformOutput',false));
                %         Rtrials = cell2mat(cellfun(@(a) a(:),Ass{s}.trialHists.RC{iDelay},'UniformOutput',false));
                %                     Ltrials = cell2mat(cellfun(@(a) a(:),cellfun(@transpose,Ass{s}.trialHists.LC{iDelay},'UniformOutput',false),'UniformOutput',false));
                %                     Rtrials = cell2mat(cellfun(@(a) a(:),cellfun(@transpose,Ass{s}.trialHists.RC{iDelay},'UniformOutput',false),'UniformOutput',false));
                nL = size(Ass_{s}.trialHists{iFile}.LD{iDelay}{1},2);
                nR = size(Ass_{s}.trialHists{iFile}.RD{iDelay}{1},2);
                
                FR = [Ltrials;Rtrials];
                evt0 = [ones(nL,1);2*ones(nR,1)];
                Ltr = size(Ass_{s}.trialHists{iFile}.RD{iDelay}{1},1);
                
                % Multivariate F-score decoder
                [D_.avgFR,D_.seFR,...
                    D_.Ft2,D_.Rt2,...
                    D_.Ft2ciL,D_.Ft2ciH,...
                    D_.TS,...
                    D_.dfnum,D_.dfden] = DecodeStats(FR,evt0,0.05);
                D_.TSsig = tpdf(D_.TS,nL+nR-2)<0.05;
                D_.TS(isnan(D_.TS))=0;
                
                %D_.CVE = DecodeCVE_mex(FR,evt0,0.05);
                [D_.CVE, D_.CVEindividual] = DecodeCVEunivariate(FR,evt0,0.05) ;
                Ass_{s}.trialHists{iFile}.D_Delay{iDelay}=D_;
                clear D_
            end
          end
          
     catch
         
     end
end

%% plot group stats for assembly diagnostics
if plotOnline
%% Plot errorbars
figure
for s=1:3
    
    Ass_{s}.membercountHist(:,sum(Ass_{s}.membercountHist)==0) = nan;
    
    subplot(2,3,1); hold on % raw bins
    ciplot(nanmean(Ass_{s}.TimeBinHist,2)+nansem(Ass_{s}.TimeBinHist,2),nanmean(Ass_{s}.TimeBinHist,2)-nansem(Ass_{s}.TimeBinHist,2),BinSizes,color_{s},'0.5')
    % stairs(BinSizes,Ass_{s}.TimeBinHist,'Color',color_{s},'LineWidth',1.5)
    set(gca,'Xscale','log','Xlim',[min(BinSizes) max(BinSizes)],'Ylim',[0 Inf])
    xlabel('Binwidth (s)')
    ylabel('Fraction of Assemblies')
    subplot(2,3,4); hold on % bin-adjusted total pattern duration
    ciplot(nanmean(Ass_{s}.DurationHist,2)+nansem(Ass_{s}.DurationHist,2),nanmean(Ass_{s}.DurationHist,2)-nansem(Ass_{s}.DurationHist,2),logspace(-1.3,1.3,10),color_{s},'0.5')
    %stairs(logspace(-1.3,1.3,10),Ass_{s}.DurationHist,'Color',color_{s},'LineWidth',1.5)
    set(gca,'Xscale','log','Ylim',[0 Inf])
    xlabel('Pattern duration (s)')
    ylabel('Fraction of Assemblies')
    
    subplot(2,3,2); hold on
    %stairs(0:10,Ass_{s}.membercountHist,'Color',color_{s},'LineWidth',1.5)
    ciplot(nanmean(Ass_{s}.membercountHist,2)+nansem(Ass_{s}.membercountHist,2),...
           nanmean(Ass_{s}.membercountHist,2)-nansem(Ass_{s}.membercountHist,2),2:20,color_{s},'0.5')
    set(gca,'Xlim',[2 20],'Ylim',[0 Inf])
    xlabel('No. interacting units')
    subplot(2,3,5); hold on
    %stairs(0:10,Ass_{s}.membercountHistAdjusted,'Color',color_{s},'LineWidth',1.5)
    ciplot(nanmean(Ass_{s}.membercountHistAjusted,2)+nansem(Ass_{s}.membercountHistAjusted,2),nanmean(Ass_{s}.membercountHistAjusted,2)-nansem(Ass_{s}.membercountHistAjusted,2),0:0.05:0.5,color_{s},'0.5')
    set(gca,'Xlim',[0 0.5],'Ylim',[0 Inf])
    xlabel('Fraction of recorded units interacting ')
    
    subplot(2,3,3); hold on
    %stairs(2:20,Ass_{s}.membersharingHist,'Color',color_{s},'LineWidth',1.5)
    ciplot(nanmean(Ass_{s}.membersharingHist,2)+nansem(Ass_{s}.membersharingHist,2),nanmean(Ass_{s}.membersharingHist,2)-nansem(Ass_{s}.membersharingHist,2),1:20,color_{s},'0.5')
    xlabel({'Promiscuity of units between assemblies'})
    set(gca,'Xlim',[0 20],'Ylim',[0 Inf])
    subplot(2,3,6); hold on
    %stairs(1:20,Ass_{s}.membersharingHist,'Color',color_{s},'LineWidth',1.5)
    ciplot(nanmean(Ass_{s}.membersharingHistAjusted,2)+nansem(Ass_{s}.membersharingHistAjusted,2),nanmean(Ass_{s}.membersharingHistAjusted,2)-nansem(Ass_{s}.membersharingHistAjusted,2),0:0.05:0.5,color_{s},'0.5')
    xlabel({'Promiscuity of units between assemblies';'(As Fraction of total Assemblies)'})
    set(gca,'Xlim',[0 0.5],'Ylim',[0 Inf])
    
end
legend(Areas)
%% Plot errorbars - mPFC only
figure
for s=1
    
    Ass_{s}.membercountHist(:,sum(Ass_{s}.membercountHist)==0) = nan;
    
    subplot(3,1,1); hold on % raw bins
    ciplot(100*(nanmean(Ass_{s}.TimeBinHist,2)+nansem(Ass_{s}.TimeBinHist,2)),...
           100*(nanmean(Ass_{s}.TimeBinHist,2)-nansem(Ass_{s}.TimeBinHist,2)),BinSizes,color_{s},'0.5')
    % stairs(BinSizes,Ass_{s}.TimeBinHist,'Color',color_{s},'LineWidth',1.5)
    set(gca,'Xscale','log','Xlim',[min(BinSizes) max(BinSizes)],'Ylim',[0 25])
%     xlabel('Bin-width (s)')
%     ylabel('Fraction of Assemblies')
    subplot(3,1,2); hold on % bin-adjusted total pattern duration
    ciplot(100*(nanmean(Ass_{s}.DurationHist,2)+nansem(Ass_{s}.DurationHist,2)),...
           100*(nanmean(Ass_{s}.DurationHist,2)-nansem(Ass_{s}.DurationHist,2)),logspace(-1.3,1.3,10),color_{s},'0.5')
    %stairs(logspace(-1.3,1.3,10),Ass_{s}.DurationHist,'Color',color_{s},'LineWidth',1.5)
    set(gca,'Xscale','log','Xlim',[0.05 11],'Ylim',[0 25])
%     xlabel('Pattern duration (s)')
%     ylabel('Fraction of Assemblies')
    
    subplot(3,1,3); hold on
    %stairs(0:10,Ass_{s}.membercountHist,'Color',color_{s},'LineWidth',1.5)
    ciplot(100*(nanmean(Ass_{s}.membercountHist,2)+nansem(Ass_{s}.membercountHist,2)),...
           100*(nanmean(Ass_{s}.membercountHist,2)-nansem(Ass_{s}.membercountHist,2)),2:20,color_{s},'0.5')
    set(gca,'Xlim',[2 10],'Ylim',[0 Inf])
%     xlabel('No. interacting units')
   
    
end

%% plot spectrum of assembly binwidth/no.Units
noMembers = 0:20;
figure
for s=1:3
    %     subplot(1,3,s)
    %     scatter(cell2mat(cellfun(@transpose,Ass_{s}.BinVector','UniformOutput',false)),cell2mat(cellfun(@transpose,Ass_{s}.noMembers','UniformOutput',false)))
    z{s} = zeros(length(BinSizes),max(noMembers));
    x = cell2mat(cellfun(@transpose,Ass_{s}.BinVector(~isempty_cell(Ass_{s}.BinVector))','UniformOutput',false));
    y = cell2mat(cellfun(@transpose,Ass_{s}.noMembers(~isempty_cell(Ass_{s}.noMembers))','UniformOutput',false));
    
    for i = 1:length(x)
        z{s}(x(i) == BinSizes,y(i) == noMembers) = z{s}(x(i) == BinSizes,y(i) == noMembers)+1;
    end
    z{s}= z{s}./length(x);
    
    subplot(1,3,s) ; hold on
    %     imagesc(noMembers,BinSizes,z{s}');
    p = pcolor(repmat(BinSizes,length(noMembers)-1,1),repmat(noMembers(2:end),length(BinSizes),1)',z{s}');
    p.LineStyle = 'none';
    set(gca,'YDir','normal','XScale','log','XLim',[min(BinSizes) max(BinSizes)],'YLim',[2 max(noMembers)])
    if s==1
        ylabel('No. member units')
    elseif s==2
        xlabel('Assembly binwidth (s)')
    end
    
    title(Areas{s})
    axis square
    
    colormap([1 1 1; parula])
    caxis([0 0.05])
end
%% plot spectrum of assembly duration/no.Units
noMembers = 0:20;
Duration_ = [0 0.1 0.5 1 2 5 10 20];

figure
for s=1:3
    %     subplot(1,3,s)
    %     scatter(cell2mat(cellfun(@transpose,Ass_{s}.BinVector','UniformOutput',false)),cell2mat(cellfun(@transpose,Ass_{s}.noMembers','UniformOutput',false)))
    z{s} = zeros(length(Duration_),max(noMembers));
    x = cell2mat(cellfun(@transpose,Ass_{s}.BinVector(~isempty_cell(Ass_{s}.Duration))','UniformOutput',false));
    y = cell2mat(cellfun(@transpose,Ass_{s}.noMembers(~isempty_cell(Ass_{s}.noMembers))','UniformOutput',false));
    for i = 1:length(x)
%         FindClosestIndex(Duration_,x(i))
        xIDX = FindClosestIndex(Duration_,x(i));
        yIDX = y(i) == noMembers;
        z{s}(xIDX ,yIDX) = z{s}(xIDX,yIDX)+1;
    end
    z{s}= z{s}./length(x);
    
    subplot(1,3,s) ; hold on
    %     imagesc(noMembers,BinSizes,z{s}');
    p = pcolor(repmat(Duration_,length(noMembers)-1,1),repmat(noMembers(2:end),length(Duration_),1)',z{s}');
    p.LineStyle = 'none';
    set(gca,'YDir','normal','XScale','log','XLim',[min(Duration_) max(Duration_)],'YLim',[2 max(noMembers)])
    if s==1
        ylabel('No. member units')
    elseif s==2
        xlabel('Assembly duration (s)')
    end
    
    title(Areas{s})
    axis square
    
    colormap([1 1 1; parula])
    caxis([0 0.05])
end
%% plot spectrum of assembly timing/binwidth
y=[];
x = [];
z={};
% for s=1:3
%     y= [ y; cell2mat(cellfun(@transpose,Ass_{s}.Duration','UniformOutput',false))];
% end
% unique(y)
timing = logspace(-1.4,2,30);
% timing = round(timing*10)/10;
% timing = min(BinSizes):min(BinSizes):max(BinSizes)*10;
% timing = [];
% for i = 1:length(BinSizes)
%     timing = [timing; BinSizes(i)*(1:10)];
% end
% timing = unique(timing(:))';
% figure; hold on
% plot((1:length(timing))/length(timing),timing)
% plot((1:length(timing_))/length(timing_),timing_)

% timing = logspace(-1.3,1,30);


figure;
for s=1:3
    %     subplot(1,3,s)
    %     scatter(cell2mat(cellfun(@transpose,Ass_{s}.BinVector','UniformOutput',false)),...
    %             cell2mat(cellfun(@transpose,Ass_{s}.noMembers','UniformOutput',false)))
    z{s} = zeros(length(BinSizes),length(timing));
    x = cell2mat(cellfun(@transpose,Ass_{s}.BinVector','UniformOutput',false));
    y = cell2mat(cellfun(@transpose,Ass_{s}.Duration','UniformOutput',false));
    %     scatter(x,y)
    for i = 1:length(x)
        idx = closest(timing,y(i));
        z{s}(x(i)==BinSizes,idx) = z{s}(x(i)==BinSizes,idx)+1;
    end
    z{s}= z{s}./length(x);
    
    subplot(1,3,s); hold on
    %     imagesc(noMembers,BinSizes,z{s}')
    p = pcolor(repmat(BinSizes,length(timing),1)',repmat(timing,length(BinSizes),1),100*z{s});
    p.LineStyle = 'none';
    set(gca,'YDir','normal',...
        'XScale','log','XLim',[min(BinSizes) max(BinSizes)],...
        'YScale','log','YLim',[min(timing) 20],...
        'XTick',[1e-2 1e-1 1],'XTickLabel',[1e-2 1e-1 1],...
        'YTick',[1e-2 1e-1 1 10 20], 'YTickLabel',[1e-2 1e-1 1 10 20])
    plot([min(BinSizes) max(BinSizes)],[min(BinSizes) max(BinSizes)],'color',0*[1 1 1],'LineWidth',3)
    for iLag = 2:10
%         plot(BinSizes,iLag*(BinSizes),'-','color',0*[1 1 1], 'LineWidth',1+0.1*(10-iLag))
        plot(BinSizes,iLag*(BinSizes),'-','color',0*[1 1 1], 'LineWidth',1)
    end
    
    if s==1
        ylabel('Pattern Duration (s)')
    elseif s==2
        xlabel('Assembly Binwidth (s)')
    end
    title(Areas{s})
    axis square
        colormap([1 1 1 ; (parula)])
    caxis([0 5])
    text(0.14, 0.11,'Synchronous','Rotation',42,'Color',0*[1 1 1])
    text(0.008, 0.11,'lagged by 10 bins','Rotation',41,'Color',0*[1 1 1])
    
end

%% timing vs duration

figure;
for s=1:3
    subplot(1,3,s); hold on
    %     scatter(cell2mat(cellfun(@transpose,Ass_{s}.BinVector','UniformOutput',false)),...
    %             cell2mat(cellfun(@transpose,Ass_{s}.noMembers','UniformOutput',false)))
    z{s} = zeros(length(BinSizes),length(timing));
    x = cell2mat(cellfun(@transpose,Ass_{s}.BinVector','UniformOutput',false));
    y = cell2mat(cellfun(@transpose,Ass_{s}.Duration','UniformOutput',false));
    
   
    for iLag = 1:10
        plot(BinSizes,iLag*(BinSizes),'-','color',0*[1 1 1], 'LineWidth',1);% +0.1*(10-iLag)
    end
    
    for i = 1:length(BinSizes)
        idx = x == BinSizes(i);
        scatter(BinSizes(i),nanmean(y(idx)),10,color_{s},'o','filled')
        errorbar(BinSizes(i),nanmean(y(idx)),nansem(y(idx)),'LineWidth',1.5,'color',color_{s})
    end
    set(gca,'XScale','log','YScale','log')
    if s==1
        ylabel('Pattern Duration (s)')
    elseif s==2
        xlabel('Pattern Timing (s)')
    end
end

figure; hold on
for iLag = 1:10
    plot(BinSizes,iLag*(BinSizes),'-','color',0.8*[1 1 1], 'LineWidth',1,'Handlevisibility','off');% +0.1*(10-iLag)
end
for s=1:3
    
    %     scatter(cell2mat(cellfun(@transpose,Ass_{s}.BinVector','UniformOutput',false)),...
    %             cell2mat(cellfun(@transpose,Ass_{s}.noMembers','UniformOutput',false)))
    z{s} = zeros(length(BinSizes),length(timing));
    x = cell2mat(cellfun(@transpose,Ass_{s}.BinVector','UniformOutput',false));
    y = cell2mat(cellfun(@transpose,Ass_{s}.Duration','UniformOutput',false));
    
    
    
    y_ = [];
    for i = 1:length(BinSizes)
        idx = x == BinSizes(i);
        y_(i) = nanmean(y(idx));
        scatter(BinSizes(i),nanmean(y(idx)),10,color_{s},'filled','Handlevisibility','off')
        errorbar(BinSizes(i),nanmean(y(idx)),nansem(y(idx)),'LineWidth',1.5,'color',color_{s},'Handlevisibility','off')
    end
    plot(BinSizes,y_,'color',color_{s},'Handlevisibility','on','LineWidth',1.5)
        set(gca,'XScale','log','YScale','log')
    ylabel('Pattern Duration (s)')
    xlabel('Pattern Binwidth (s)')
end
legend(Areas)
%% duration violin plots

figure;
for s=1:3
    subplot(1,3,s); hold on
    z{s} = zeros(length(BinSizes),length(timing));
    x = cell2mat(cellfun(@transpose,Ass_{s}.BinVector','UniformOutput',false));
    y = cell2mat(cellfun(@transpose,Ass_{s}.Duration','UniformOutput',false));
    
    for i = 1:length(BinSizes)
        idx = x == BinSizes(i);
        Ns(i) = nansum(idx);
    end
    Y=[];
    for i = 1:length(BinSizes)
        idx = x == BinSizes(i);
        Y = [Y,[y(idx);nan(max(Ns) - Ns(i),1)]];
    end
    violin (Y,'bw',0.5,'xlabel',string(BinSizes))
    set(gca,'Xtick',1:length(BinSizes),'xticklabel',string(BinSizes),'XtickLabelRotation',-45 )
    axis()
end

%% Group decoding - L/R
tb = Ass_{1, 1}.trialHists{1, 1}.tb;
for s =1:3
    for iDelay = 1:3
        for iFile =1:length(fileList)
            if ~isempty(Ass_{s}.trialHists{iFile})
                AssCat.D_LR.Ft2{s}{iDelay}(iFile,:) = Ass_{s}.trialHists{iFile}.D_LR{iDelay}.Ft2;
            else
                AssCat.D_LR.Ft2{s}{iDelay}(iFile,:) = nan(1,length(tb));
            end
        end
        AssCat.D_LR.Ft2mean{s}(iDelay,:)=nanmean(AssCat.D_LR.Ft2{s}{iDelay});
        AssCat.D_LR.Ft2SEM{s}(iDelay,:)=nansem(AssCat.D_LR.Ft2{s}{iDelay});
    end
end
col_ = (gray(4));

figure;
for s =1:3
    subplot(1,3,s); hold on
    for iDelay = 1:3
        
%         plot(x,D_.Ft2,'color',col_(iDelay,:),'Linewidth',1.5)
%         plot(x,D_.Ft2ciH,'color',col_(iDelay,:),'lineStyle',':','HandleVisibility','off')
        %  plot(x,1-D_.CVE,'color',color_{iDelay})
        plot([5 5],[0 5],'color',[0 1 0 0.1],'LineWidth',2,'HandleVisibility','off')
        plot([15 15],[0 5],'color',[1 0 0 0.1],'LineWidth',2,'HandleVisibility','off')
        plot([10 10],[0 5],'color',[0 0 0 0.1],'HandleVisibility','off')
                   
                    
        M =  AssCat.D_LR.Ft2mean{s}(iDelay,:);
        E =  AssCat.D_LR.Ft2SEM{s}(iDelay,:);
        ciplot(M+E,M-E,tb,col_(iDelay,:),1)
        plot([0 0],[0 1],'k','HandleVisibility','off')
        axis([min(tb) max(tb) 0 5])
        title(Areas{s})
    end
end
        legend(Delays_)
%% Group decoding - NP
tb = (1:length(Ass_{1, 1}.trialHists{1, 1}.D_NP{1, 1}.Ft2))*0.1;tb = tb-mean(tb);
for s =1:3
    for iDelay = 1:3
        for iFile =1:length(fileList)
            if ~isempty(Ass_{s}.trialHists{iFile})
                AssCat.D_NP.Ft2{s}{iDelay}(iFile,:) = Ass_{s}.trialHists{iFile}.D_NP{iDelay}.Ft2;
            else
                AssCat.D_NP.Ft2{s}{iDelay}(iFile,:) = nan(1,length(tb));
            end
        end
        AssCat.D_NP.Ft2mean{s}(iDelay,:)=nanmean(AssCat.D_NP.Ft2{s}{iDelay});
        AssCat.D_NP.Ft2SEM{s}(iDelay,:)=nansem(AssCat.D_NP.Ft2{s}{iDelay});
    end
end
col_ = (gray(4));

figure;
for s =1:3
    subplot(1,3,s); hold on
    for iDelay = 1:3
        
%         plot(x,D_.Ft2,'color',col_(iDelay,:),'Linewidth',1.5)
%         plot(x,D_.Ft2ciH,'color',col_(iDelay,:),'lineStyle',':','HandleVisibility','off')
        %  plot(x,1-D_.CVE,'color',color_{iDelay})
%         plot([5 5],[0 5],'color',[0 1 0 0.1],'LineWidth',2,'HandleVisibility','off')
%         plot([15 15],[0 5],'color',[1 0 0 0.1],'LineWidth',2,'HandleVisibility','off')
        plot([0 0],[0 5],'color',[0 0 0 0.1],'HandleVisibility','off')
                   
                    
        M =  AssCat.D_NP.Ft2mean{s}(iDelay,:);
        E =  AssCat.D_NP.Ft2SEM{s}(iDelay,:);
        ciplot(M+E,M-E,tb,col_(iDelay,:),0.8)
        plot([0 0],[0 1],'k','HandleVisibility','off')
        axis([min(tb) max(tb) 0 5])
        title(Areas{s})
    end
end
        legend(Delays_)
%% Peak cue information vs. binwidth
clear M E
figure; hold on
for s = 1:3
    temp=[];
    for iFile =1:length(fileList)
        if ~isempty(Ass_{s}.BinVector{iFile})
            A =  [nanmax([Ass_{s}.trialHists{iFile}.D_LR{1}.TS;...
                Ass_{s}.trialHists{iFile}.D_LR{2}.TS;...
                Ass_{s}.trialHists{iFile}.D_LR{3}.TS]);Ass_{s}.BinVector{iFile}]';
            temp = [temp;A];
        end
    end
    BinSizes_ = unique(A(:,2));
    for i = 1:length(BinSizes)
        M(i) = nanmean(temp(A==BinSizes(i),1));
        E(i) = nansem(temp(A==BinSizes(i),1));
        scatter(BinSizes(i),M(i),10,color_{s},'filled','Handlevisibility','off')
        errorbar(BinSizes(i),M(i),E(i),'LineWidth',1.5,'color',color_{s})
    end
    plot(BinSizes,M,'LineWidth',1.5,'color',color_{s})

end
set(gca,'XScale','log')
xlabel('Pattern binwidth (s)')
ylabel('Peak cue information (max. t-score)')

%% Peak cue information vs. duration
clear M E
figure; hold on
for s = 1:3
    temp=[];
    for iFile =1:length(fileList)
        if ~isempty(Ass_{s}.BinVector{iFile})
            A =  [nanmax([Ass_{s}.trialHists{iFile}.D_LR{1}.TS;...
                Ass_{s}.trialHists{iFile}.D_LR{2}.TS;...
                Ass_{s}.trialHists{iFile}.D_LR{3}.TS]);Ass_{s}.Duration{iFile}]';
            temp = [temp;A];
        end
    end
    BinSizes_ = unique(A(:,2));
    clear M E
    for i = 1:length(BinSizes_)
        M(i) = nanmean(temp(A==BinSizes_(i),1));
        E(i) = nansem(temp(A==BinSizes_(i),1));
        scatter(BinSizes_(i),M(i),10,color_{s},'filled','Handlevisibility','off')
        errorbar(BinSizes_(i),M(i),E(i),'LineWidth',1.5,'color',color_{s})
    end
    plot(BinSizes_,M,'LineWidth',1.5,'color',color_{s})

end
set(gca,'XScale','log')
xlabel('Assembly duration (s)')
ylabel('Peak cue information (max. t-score)')
%% Peak cue information vs. no members
clear M E
figure; hold on
for s = 1:3
    temp=[];
    for iFile =1:length(fileList)
        if ~isempty(Ass_{s}.BinVector{iFile})
            A =  [nanmax([Ass_{s}.trialHists{iFile}.D_LR{1}.TS;...
                Ass_{s}.trialHists{iFile}.D_LR{2}.TS;...
                Ass_{s}.trialHists{iFile}.D_LR{3}.TS]);Ass_{s}.noMembers{iFile}]';
            temp = [temp;A];
        end
    end
    noMembers = unique(A(:,2));
    clear M E
    for i = 1:length(noMembers)
        M(i) = nanmean(temp(A==noMembers(i),1));
        E(i) = nansem( temp(A==noMembers(i),1));
        scatter(noMembers(i),M(i),10,color_{s},'filled','Handlevisibility','off')
        errorbar(noMembers(i),M(i),E(i),'LineWidth',1.5,'color',color_{s})
    end
    plot(noMembers,M,'LineWidth',1.5,'color',color_{s})

end
xlabel('Assembly size (No. member units)')
ylabel('Peak cue information (max. t-score)')
%% Assembly activation time vs. bin width
clear temp
tb = Ass_{1}.trialHists{1}.tb;
col_ = jet(length(BinSizes));

for s=1:3
    temp{s} = cell(3,length(BinSizes));
    for i =1:length(BinSizes)
        for iFile =1:length(fileList)
            if ~isempty(Ass_{s}.BinVector{iFile})
                idx = find(Ass_{s}.BinVector{iFile}==BinSizes(i));
              
                for iDelay = 1:3
                    a = (Ass_{s}.trialHists{iFile}.LC_mean{iDelay}(idx,:)+Ass_{s}.trialHists{iFile}.RC_mean{iDelay}(idx,:))./2;
                    a = smooth2a(a,0,10);
%                     a = full(a);
%                     for iUnit = 1:size(a,1)
%                         a(iUnit,:) = mat2gray(a(iUnit,:));
%                     end
%                     a = zscore(a,[],2);
                    a = nanmean(a);
%                     a = zscore(a,[],2);
                    temp{s}{iDelay,i}(iFile,1:length(tb)) = a;
                end
            else
                temp{s}{iDelay,i}(iFile,:) = nan(size(tb));
            end
        end
    end
    b{s} = cell(1,length(BinSizes));
    for i =1:length(BinSizes)
%         b{s}{iBin} = (temp{s}{1,iBin} + temp{s}{2,iBin} + temp{s}{3,iBin})./3
        b{s}{i} = (temp{s}{2,i} + temp{s}{3,i})./3
    end
end

figure
for s=1:3
    % Plot specific delay
    iDelay = 3
    subplot(1,3,s); hold on
    for i =1:length(BinSizes)
        a = temp{s}{iDelay,i};
        a(isnan(a(:,1)),:)= [];
        a=zscore(a')';
        plot(tb,nanmean(a),'color',[col_(i,:) 0.5],'LineWidth',1.5)
%         ciplot(nanmean(a)+nansem(a)+iBin,nanmean(a)-nansem(a)+iBin,tb,col_(iBin,:))
    end
end

figure;
for s = 1:3
    subplot(1,3,s); hold on
    for i =1:length(BinSizes)
        a = b{s}{i};
        a(isnan(a(:,1)),:)= [];
        a=zscore(a')';
%         plot(tb,nanmean(a),'color',[col_(iBin,:) 0.5],'LineWidth',1.5)
        ciplot(nanmean(a)+nansem(a)+i,nanmean(a)-nansem(a)+i,tb,col_(i,:))
    end
end

%% Assembly decoding vs. bin width
clear temp
tb = Ass_{1}.trialHists{1}.tb;
col_ = jet(length(BinSizes));

for s=1:3
    temp{s} = cell(3,length(BinSizes));
    for i =1:length(BinSizes)
        for iFile =1:length(fileList)
            if ~isempty(Ass_{s}.BinVector{iFile})
                idx = find(Ass_{s}.BinVector{iFile}==BinSizes(i));
              
                for iDelay = 1:3
%                    a = (Ass_{s}.trialHists{iFile}.LC_mean{iDelay}(idx,:)+Ass_{s}.trialHists{iFile}.RC_mean{iDelay}(idx,:))./2;
                    a =  Ass_{s}.trialHists{iFile}.D_LR{iDelay}.TS(:,idx)';
                    a = smooth2a(a,0,2);
%                     a = full(a);
%                     for iUnit = 1:size(a,1)
%                         a(iUnit,:) = mat2gray(a(iUnit,:));
%                     end
%                     a = zscore(a,[],2);
                    a = nanmean(a);
%                     a = zscore(a,[],2);
                    temp{s}{iDelay,i}(iFile,1:length(tb)) = a;
                end
            else
                temp{s}{iDelay,i}(iFile,:) = nan(size(tb));
            end
        end
    end
    b{s} = cell(1,length(BinSizes));
    for i =1:length(BinSizes)
%         b{s}{iBin} = (temp{s}{1,iBin} + temp{s}{2,iBin} + temp{s}{3,iBin})./3
        b{s}{i} = (temp{s}{2,i} + temp{s}{3,i})./3
    end
end

figure
for s=1:3
    % Plot specific delay
    iDelay = 2
    subplot(1,3,s); hold on
    for i =1:length(BinSizes)
        a = temp{s}{iDelay,i};
        a(isnan(a(:,1)),:)= [];
%         a=zscore(a')';
%         plot(tb,nanmean(a)+iBin,'color',[col_(iBin,:) 0.5],'LineWidth',1.5)
        ciplot(nanmean(a)+nansem(a)+i,nanmean(a)-nansem(a)+i,tb,col_(i,:),0.6)
    end
end

figure;
for s = 1:3
    subplot(1,3,s); hold on
    plot([10 10],[0 15],':k')
    plot([5 5],[0 15],'g')
    plot([15 15],[0 15],'r')
    for i =1:length(BinSizes)
        a = b{s}{i};
        a(isnan(a(:,1)),:)= [];
        a=zscore(a')';
        plot(tb,nanmean(a)+2*i,'color',[col_(i,:) 0.5],'LineWidth',1.5)
%         ciplot(nanmean(a)+nansem(a)+iBin,nanmean(a)-nansem(a)+iBin,tb,col_(iBin,:),0.6)
    end
end

%% Assembly activation vs. [synchrony, bin width] 
clear temp
tb = Ass_{1}.trialHists{1}.tb;
col_ = jet(length(BinSizes));

for s=1:3
    for i =1:length(BinSizes)
        for iFile =1:length(fileList)
            for iLag = 1:max(MaxLags)
                for iDelay = 1:3
                    temp{s}{iDelay}{i,iLag} = nan(length(fileList),length(tb));
                end
            end
        end
    end
end

for s=1:3
    for i =1:length(BinSizes)
        for iFile =1:length(fileList)
            if ~isempty(Ass_{s}.BinVector{iFile})
                % first match binwidth
                idx = find(Ass_{s}.BinVector{iFile}==BinSizes(i));
                lags = Ass_{s}.Duration{iFile}(idx)./Ass_{s}.BinVector{iFile}(idx);
                
                for iLag = 1:max(MaxLags)
%                     [s i iLag iDelay iFile]
                    % next match synchrony
                    idx2 = find(lags==iLag);
                    AssList =idx(idx2);
                    for iDelay = 1:3
                        % mean activity:
                        a = (Ass_{s}.trialHists{iFile}.LC_mean{iDelay}(AssList,:) +...
                             Ass_{s}.trialHists{iFile}.RC_mean{iDelay}(AssList,:) )./2;
                         
                        
%                         a = smooth2a(a,0,10);
                        %                     a = full(a);
                        %                     for iUnit = 1:size(a,1)
                        %                         a(iUnit,:) = mat2gray(a(iUnit,:));
                        %                     end
                        %                     a = zscore(a,[],2);
                        a = nanmean(a,1);
                        %                     a = zscore(a,[],2);
                        temp{s}{iDelay}{i,iLag}(iFile,1:length(tb)) = a;
                    end
                end
                
            end
        end
    end
end
%% Plot activation (1)
iDelay = 3;
for s=1:3
    figure;hold on
    
    for i =1:length(BinSizes)
        for iLag = 1:max(MaxLags)
            a = temp{s}{iDelay}{i,iLag};
            a(isnan(a(:,1)),:)= [];
            a=zscore(a')';
%             plot(tb+(iLag-1)*25,nanmean(a)+2*i,'color',[col_(i,:) 0.5],'LineWidth',1.5)
            plot(tb+(iLag-1)*25,nanmean(a)+2*i,'color',[0 0 0 1],'LineWidth',1.5)
        end
    end
    
    set(gca,'ytick',2*(1:length(BinSizes)),'YTickLabel',BinSizes,...
        'xtick',mean(tb)+((1:max(MaxLags))-1)*25,'XTickLabel',1:10)
    
    ylabel('Assembly binwidth (s)')
    xlabel('Assembly duration (no. bins)')
    title([Areas{s} ' Assemblies'])
end
%% Plot activation (2)
a=[];
for s=1:3
    for iDelay =1:3
        a = [a,max(cellfun(@(X) nanmax(nanmax(X)), temp{s}{iDelay}))];
    end
end
a = nanmax(a);
maxCol = 0.5;%ceil(a);
ColRange = (0.1:0.1:1)*maxCol;
cmap_ = jet(10);
iDelay = 2;
clear a
for s=1:3
    
      a = max(cellfun(@(X) nanmax(nanmean(X)), temp{s}{iDelay}));
a = nanmax(a);
maxCol = a %ceil(a);
ColRange = (0.1:0.1:1)*maxCol;
cmap_ = jet(10);
clear a


    figure;hold on
    
    for i =1:length(BinSizes)
        for iLag = 1:max(MaxLags)
            a = temp{s}{iDelay}{i,iLag};
            a(isnan(a(:,1)),:) = [];
            if ~isempty(a)
                colIdx = FindClosestIndex(ColRange,nanmax(nanmean(a)));
                a = 10*(a-nanmin(nanmean(a)));
                %plot(tb+(iLag-1)*25,nanmean(a)+2*i,'color',[col_(i,:) 0.5],'LineWidth',1.5)
%                 ciplot(nanmean(a,1)+nansem(a,1)+2*i,...
%                        nanmean(a,1)-nansem(a,1)+2*i,...
%                        tb+(iLag-1)*25,cmap_(colIdx,:))
                plot(tb+(iLag-1)*25,nanmean(a,1)+2*i,'color',[cmap_(colIdx,:) 0.5],'LineWidth',1.5)

            end
        end
    end
    
    set(gca,'ytick',2*(1:length(BinSizes)),'YTickLabel',BinSizes,...
        'xtick',mean(tb)+((1:max(MaxLags))-1)*25,'XTickLabel',1:10)
    
    ylabel('Assembly binwidth (s)')
    xlabel('Assembly duration (no. bins)')
    title([Areas{s} ' Assemblies'])
end

%% Assembly activation vs. [synchrony, bin width] 
clear temp
tb = Ass_{1}.trialHists{1}.tb;
col_ = jet(length(BinSizes));

for s=1:3
    for i =1:length(BinSizes)
        for iFile =1:length(fileList)
            for iLag = 1:max(MaxLags)
                for iDelay = 1:3
                    temp{s}{iDelay}{i,iLag} = nan(length(fileList),length(tb));
                end
            end
        end
    end
end

for s=1:3
    for i =1:length(BinSizes)
        for iFile =1:length(fileList)
                                [s i]

            if ~isempty(Ass_{s}.BinVector{iFile})
                % first match binwidth
                idx = find(Ass_{s}.BinVector{iFile} == BinSizes(i));
                lags = Ass_{s}.Duration{iFile}(idx)./Ass_{s}.BinVector{iFile}(idx);
                
                for iLag = 1:max(MaxLags)
%                     [s i iLag iDelay iFile]
                    % next match synchrony
                    idx2 = find(lags==iLag);
                    AssList =idx(idx2);
                    for iDelay = 1:3

                        % L/R T-score:
                        a = Ass_{s}.trialHists{iFile}.D_LR{iDelay}.TS(:,AssList)';
                        
                        
                        a = smooth2a(a,0,10);
                        %                     a = full(a);
                        %                     for iUnit = 1:size(a,1)
                        %                         a(iUnit,:) = mat2gray(a(iUnit,:));
                        %                     end
                        %                     a = zscore(a,[],2);
                        a = nanmean(a,1);
                        %                     a = zscore(a,[],2);
                        temp{s}{iDelay}{i,iLag}(iFile,1:length(tb)) = a;
                    end
                end
                
            end
        end
    end
end
%% Plot activation (1)
iDelay = 3;
for s=1:3
    figure;hold on
    
    for i =1:length(BinSizes)
        for iLag = 1:max(MaxLags)
            a = temp{s}{iDelay}{i,iLag};
            a(isnan(a(:,1)),:)= [];
            a=zscore(a')';
%             plot(tb+(iLag-1)*25,nanmean(a)+2*i,'color',[col_(i,:) 0.5],'LineWidth',1.5)
            plot(tb+(iLag-1)*25,nanmean(a)+2*i,'color',[0 0 0 1],'LineWidth',1.5)
        end
    end
    
    set(gca,'ytick',2*(1:length(BinSizes)),'YTickLabel',BinSizes,...
        'xtick',mean(tb)+((1:max(MaxLags))-1)*25,'XTickLabel',1:10)
    
    ylabel('Assembly binwidth (s)')
    xlabel('Assembly duration (no. bins)')
    title([Areas{s} ' Assemblies'])
end
%% Plot activation (2)
a=[];
for s=1:3
    for iDelay =1:3
        a = [a,max(cellfun(@(X) nanmax(nanmax(X)), temp{s}{iDelay}))];
    end
end
a = nanmax(a);
maxCol = 0.5;%ceil(a);
ColRange = (0.1:0.1:1)*maxCol;
cmap_ = jet(10);
iDelay = 2;
clear a
for s=1:3
    
      a = max(cellfun(@(X) nanmax(nanmean(X)), temp{s}{iDelay}));
a = nanmax(a);
maxCol = a %ceil(a);
ColRange = (0.1:0.1:1)*maxCol;
cmap_ = jet(10);
clear a


    figure;hold on
    
    for i =1:length(BinSizes)
        for iLag = 1:max(MaxLags)
            a = temp{s}{iDelay}{i,iLag};
            a(isnan(a(:,1)),:) = [];
            if ~isempty(a)
                colIdx = FindClosestIndex(ColRange,nanmax(nanmean(a)));
                a = 1*(a-nanmin(nanmean(a)));
                plot(tb+(iLag-1)*25,nanmean(a)+2*i,'color',[col_(i,:) 0.5],'LineWidth',1.5)
%                 ciplot(nanmean(a,1)+nansem(a,1)+2*i,...
%                        nanmean(a,1)-nansem(a,1)+2*i,...
%                        tb+(iLag-1)*25,cmap_(colIdx,:))
%                 plot(tb+(iLag-1)*25,nanmean(a,1)+2*i,'color',[cmap_(colIdx,:) 0.5],'LineWidth',1.5)

            end
        end
    end
    
    set(gca,'ytick',2*(1:length(BinSizes)),'YTickLabel',BinSizes,...
        'xtick',mean(tb)+((1:max(MaxLags))-1)*25,'XTickLabel',1:10)
    
    ylabel('Assembly binwidth (s)')
    xlabel('Assembly duration (no. bins)')
    title([Areas{s} ' Assemblies'])
end

%% Delay period decoding  - group decoding
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    for s =1:3
        
        for iFile =1:length(fileList)
            if ~isempty(Ass_{s}.trialHists{iFile})
%                 AssCat.D_Delay.Ft2{s}{iDelay}(iFile,:) =  Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.Ft2;
                AssCat.D_Delay.CVE{s}{iDelay}(iFile,:) =  1-Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.CVE;
            else
%                 AssCat.D_Delay.Ft2{s}{iDelay}(iFile,:) = nan(1,length(tb));
                AssCat.D_Delay.CVE{s}{iDelay}(iFile,:) = nan(1,length(tb));
            end
        end
        temp = AssCat.D_Delay.CVE{s}{iDelay};
%         temp = AssCat.D_Delay.Ft2{s}{iDelay};
        temp = smooth2a(temp,0,5);
%         AssCat.D_Delay.Ft2mean{s}{iDelay} = nanmean(temp);
%         AssCat.D_Delay.Ft2SEM{s}{iDelay}  = nansem(temp);       
        AssCat.D_Delay.CVEmean{s}{iDelay} = nanmean(temp);
        AssCat.D_Delay.CVESEM{s}{iDelay}  = nansem(temp);
    end
end

figure;
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    subplot(length(Delays),1,iDelay); hold on
    for s =1:3
%         m = AssCat.D_Delay.Ft2mean{s}{iDelay};
%         e = AssCat.D_Delay.Ft2SEM{s}{iDelay};
        m = AssCat.D_Delay.CVEmean{s}{iDelay};
        e = AssCat.D_Delay.CVESEM{s}{iDelay};
        ciplot(m+e,m-e,tb,color_{s},0.6);
    end
    plot([0 0],[0 5],'g','LineWidth',1.5)
    for iDelay =1:3
       plot([Delays(iDelay) Delays(iDelay)],[0 5],':k') 
    end
%     axis([-5 16 0 5])
    axis([-5 16 0 1])
end
%% Delay period decoding  - % sig decoding assemblies
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    for s =1:3
        
        for iFile =1:length(fileList)
            if ~isempty(Ass_{s}.trialHists{iFile})
                temp = Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.TSsig;
                temp = nansum(temp,2)./size(temp,2);
                AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) =  temp;
            else
                AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) = nan(1,length(tb));
            end
        end
        temp = AssCat.D_Delay.fracSig{s}{iDelay};
%         temp = smooth2a(temp,0,10);
        AssCat.D_Delay.fracSigmean{s}{iDelay} = nanmean(temp);
        AssCat.D_Delay.fracSigSEM{s}{iDelay}  = nansem(temp);
    end
end

figure;
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    subplot(length(Delays),1,iDelay); hold on
    for s =1:3
        m = AssCat.D_Delay.fracSigmean{s}{iDelay};
        e = AssCat.D_Delay.fracSigSEM{s}{iDelay};
        ciplot(m+e,m-e,tb,color_{s},0.6);
    end
    plot([0 0],[0 5],'g','LineWidth',1.5)
    for iDelay =1:3
       plot([Delays(iDelay) Delays(iDelay)],[0 5],':k') 
    end
    axis([-5 16 0 0.2])
end
%% Delay period decoding  - % sig decoding assemblies by bin
for iBin = 1:length(BinSizes)
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    for s =1:3
        
        for iFile =1:length(fileList) 
            AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) = nan(1,length(tb));
            if ~isempty(Ass_{s}.trialHists{iFile})
                temp = Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.TSsig;
                idx =  Ass_{s}.BinVector_{iFile} == BinSizes(iBin);
                temp(:,~idx)=[];
                if ~isempty(temp)
                    temp = nansum(temp,2)./size(temp,2);
                    AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) =  temp; 
                end
                
            else
               
            end
        end
        temp = AssCat.D_Delay.fracSig{s}{iDelay};
%         temp = smooth2a(temp,0,10);
        AssCat.D_Delay.fracSigmean{s}{iDelay} = nanmean(temp);
        AssCat.D_Delay.fracSigSEM{s}{iDelay}  = nansem(temp);
    end
end

figure('name', sprintf('%fs bins',BinSizes(iBin)));
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    subplot(length(Delays),1,iDelay); hold on
    for s =1:3
        m = AssCat.D_Delay.fracSigmean{s}{iDelay};
        e = AssCat.D_Delay.fracSigSEM{s}{iDelay};
        ciplot(m+e,m-e,tb,color_{s},0.6);
    end
    plot([0 0],[0 5],'g','LineWidth',1.5)
    for iDelay =1:3
       plot([Delays(iDelay) Delays(iDelay)],[0 5],':k') 
    end
    axis([-5 16 0 0.5])
end

end
%% Delay period decoding  - % group decoding assemblies by bin (CVE)
for iBin = 1:length(BinSizes)
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    for s =1:3
        
        for iFile =1:length(fileList) 
            AssCat.D_Delay.CVEindividualMean{s}{iDelay}(iFile,:) = nan(1,length(tb));
            if ~isempty(Ass_{s}.trialHists{iFile})
                temp = Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.CVEindividual;
                idx =  Ass_{s}.BinVector_{iFile} == BinSizes(iBin);
                temp(:,~idx)=[];
                if ~isempty(temp)
                    temp = nansum(temp,2)./size(temp,2);
                    AssCat.D_Delay.CVEindividualMean{s}{iDelay}(iFile,:) =  temp; 
                end
                
            else
               
            end
        end
        temp = AssCat.D_Delay.CVEindividualMean{s}{iDelay};
%         temp = smooth2a(temp,0,10);
        AssCat.D_Delay.CVEindividualMeanmean{s}{iDelay} = nanmean(temp);
        AssCat.D_Delay.CVEindividualMeanSEM{s}{iDelay}  = nansem(temp);
    end
end

figure('name', sprintf('%fs bins',BinSizes(iBin)));
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    subplot(length(Delays),1,iDelay); hold on
    for s =1:3
        m = AssCat.D_Delay.fracSigmean{s}{iDelay};
        e = AssCat.D_Delay.fracSigSEM{s}{iDelay};
        ciplot(m+e,m-e,tb,color_{s},0.6);
    end
    plot([0 0],[0 5],'g','LineWidth',1.5)
    for iDelay =1:3
       plot([Delays(iDelay) Delays(iDelay)],[0 5],':k') 
    end
    axis([-5 16 0 0.5])
end

end
%% Delay period decoding  - % sig decoding assemblies by bin lumped
BinSizes_ = [0, 0.2, 0.5, 1]
for iBin = 1:length(BinSizes_)-1
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        for s =1:3
            
            for iFile =1:length(fileList)
                AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) = nan(1,length(tb));
                if ~isempty(Ass_{s}.trialHists{iFile})
                    temp = Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.TSsig;
                    idx =  find(Ass_{s}.BinVector_{iFile} >= BinSizes_(iBin) & ...
                        Ass_{s}.BinVector_{iFile} < BinSizes_(iBin+1));
                    temp = temp(:,idx);
                    if ~isempty(temp)
                        temp = nansum(temp,2)./size(temp,2);
                        AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) =  temp;
                    end
                end
            end
            temp = AssCat.D_Delay.fracSig{s}{iDelay};
            %         temp = smooth2a(temp,0,10);
            AssCat.D_Delay.fracSigmean{s}{iDelay} = nanmean(temp);
            AssCat.D_Delay.fracSigSEM{s}{iDelay}  = nansem(temp);
        end
    end
    
    figure('name', sprintf('%gs~%gs bins',BinSizes_(iBin),BinSizes_(iBin+1)));
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        subplot(length(Delays),1,iDelay); hold on
        for s =1:3
            m = AssCat.D_Delay.fracSigmean{s}{iDelay};
            e = AssCat.D_Delay.fracSigSEM{s}{iDelay};
            ciplot(m+e,m-e,tb,color_{s},0.6);
        end
        plot([0 0],[0 5],'g','LineWidth',1.5)
        for iDelay =1:3
            plot([Delays(iDelay) Delays(iDelay)],[0 5],':k')
        end
        axis([-5 16 0 0.5])
    end
    
end
%% Delay period decoding  - % active assemblies by duration lumped
BinSizes_ = [0:10]
for iBin = 1:length(BinSizes_)-1
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        for s =1:3
            
            for iFile =1:length(fileList)
                AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) = nan(1,length(tb));
                if ~isempty(Ass_{s}.trialHists{iFile})
                    temp = Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.TSsig;
                    idx =  find(Ass_{s}.Duration{iFile} >= BinSizes_(iBin) & ...
                                Ass_{s}.Duration{iFile} < BinSizes_(iBin+1));
                    temp = temp(:,idx);
                    if ~isempty(temp)
                        temp = nansum(temp,2)./size(temp,2);
                        AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) =  temp;
                    end
                end
            end
            temp = AssCat.D_Delay.fracSig{s}{iDelay};
                    temp = smooth2a(temp,0,2);
            AssCat.D_Delay.fracSigmean{s}{iDelay} = nanmean(temp);
            AssCat.D_Delay.fracSigSEM{s}{iDelay}  = nansem(temp);
        end
    end
    
    figure('name', sprintf('%gs~%gs duration',BinSizes_(iBin),BinSizes_(iBin+1)));
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        subplot(length(Delays),1,iDelay); hold on
        for s =1:3
            m = AssCat.D_Delay.fracSigmean{s}{iDelay};
            e = AssCat.D_Delay.fracSigSEM{s}{iDelay};
            ciplot(m+e,m-e,tb,color_{s},0.6);
        end
        plot([0 0],[0 5],'g','LineWidth',1.5)
        for iDelay =1:3
            plot([Delays(iDelay) Delays(iDelay)],[0 5],':k')
        end
        axis([-5 16 0 1])
    end
    
end
%% Delay period decoding  - % active assemblies
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    for s =1:3
        
        for iFile =1:length(fileList)
            if ~isempty(Ass_{s}.trialHists{iFile})
                temp = (Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.avgFR{1}+...
                       Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.avgFR{1})./2;
%                 AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) =  (nanmean(zscore(temp,[],1),2));
                AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) =  (nanmean(zscore(temp,[],1),2));
            else
                AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) = nan(1,length(tb));
            end
        end
        temp = AssCat.D_Delay.fracSig{s}{iDelay};
%         temp = smooth2a(temp,0,10);
        AssCat.D_Delay.fracSigmean{s}{iDelay} = nanmean(temp);
        AssCat.D_Delay.fracSigSEM{s}{iDelay}  = nansem(temp);
    end
end

figure;
for iDelay = 1:3
    tb = tlimsSample(1):bw:Delays(iDelay);
    subplot(length(Delays),1,iDelay); hold on
    for s =1:3
        m = AssCat.D_Delay.fracSigmean{s}{iDelay};
        e = AssCat.D_Delay.fracSigSEM{s}{iDelay};
        ciplot(m+e,m-e,tb,color_{s},0.6);
    end
    plot([0 0],[0 5],'g','LineWidth',1.5)
    for iDelay =1:3
       plot([Delays(iDelay) Delays(iDelay)],[0 5],':k') 
    end
    axis([-5 16 -2 2])
end
%% Delay period decoding  - % active assemblies by bin
for iBin = 1:length(BinSizes)
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        for s =1:3
            
            for iFile =1:length(fileList)
                AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) = nan(1,length(tb));
                if ~isempty(Ass_{s}.trialHists{iFile})
                    temp = (Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.avgFR{1}+...
                            Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.avgFR{1})./2;
                    
                    idx =  Ass_{s}.BinVector_{iFile} == BinSizes(iBin);
                    temp(:,~idx)=[];
                    if ~isempty(temp)
                        temp = zscore(temp,[],1);
                        temp = nanmean(temp,2);
                        AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) =  temp;
                    end
                    
                else
                    
                end
            end
            temp = AssCat.D_Delay.fracSig{s}{iDelay};
            %         temp = smooth2a(temp,0,10);
            AssCat.D_Delay.fracSigmean{s}{iDelay} = nanmean(temp);
            AssCat.D_Delay.fracSigSEM{s}{iDelay}  = nansem(temp);
        end
    end
    
    figure('name', sprintf('%fs bins',BinSizes(iBin)));
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        subplot(length(Delays),1,iDelay); hold on
        for s =1:3
            m = AssCat.D_Delay.fracSigmean{s}{iDelay};
            e = AssCat.D_Delay.fracSigSEM{s}{iDelay};
            ciplot(m+e,m-e,tb,color_{s},0.6);
        end
        plot([0 0],[0 5],'g','LineWidth',1.5)
        for iDelay =1:3
            plot([Delays(iDelay) Delays(iDelay)],[0 5],':k')
        end
    axis([-5 16 -2 2])
%         axis([-5 16 0 0.1])
    end
    
end
%% Delay period decoding  - % active assemblies by bin lumped
BinSizes_ = [0, 0.2, 0.5, 1]
for iBin = 1:length(BinSizes_)-1
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        for s =1:3
            
            for iFile =1:length(fileList)
                AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) = nan(1,length(tb));
                if ~isempty(Ass_{s}.trialHists{iFile})
                    temp = (Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.avgFR{1}+...
                            Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.avgFR{1})./2;
                    idx =  find(Ass_{s}.BinVector_{iFile} >= BinSizes_(iBin) & ...
                        Ass_{s}.BinVector_{iFile} < BinSizes_(iBin+1));
                    temp = temp(:,idx);
                    if ~isempty(temp)
                        temp = zscore(temp,[],1);
                        temp = nanmean(temp,2);
                        AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) =  temp;
                    end
                end
            end
            temp = AssCat.D_Delay.fracSig{s}{iDelay};
            %         temp = smooth2a(temp,0,10);
            AssCat.D_Delay.fracSigmean{s}{iDelay} = nanmean(temp);
            AssCat.D_Delay.fracSigSEM{s}{iDelay}  = nansem(temp);
        end
    end
    
    figure('name', sprintf('%gs~%gs bins',BinSizes_(iBin),BinSizes_(iBin+1)));
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        subplot(length(Delays),1,iDelay); hold on
        for s =1:3
            m = AssCat.D_Delay.fracSigmean{s}{iDelay};
            e = AssCat.D_Delay.fracSigSEM{s}{iDelay};
            ciplot(m+e,m-e,tb,color_{s},0.6);
        end
        plot([0 0],[0 5],'g','LineWidth',1.5)
        for iDelay =1:3
            plot([Delays(iDelay) Delays(iDelay)],[0 5],':k')
        end
        axis([-5 16 -2 2])
    end
    
end
%% Delay period decoding  - % active assemblies by duration lumped
BinSizes_ = [0:20]
for iBin = 1:length(BinSizes_)-1
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        for s =1:3
            
            for iFile =1:length(fileList)
                AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) = nan(1,length(tb));
                if ~isempty(Ass_{s}.trialHists{iFile})
                    temp = (Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.avgFR{1}+...
                            Ass_{s}.trialHists{iFile}.D_Delay{iDelay}.avgFR{1})./2;
                    idx =  find(Ass_{s}.Duration{iFile} >= BinSizes_(iBin) & ...
                        Ass_{s}.Duration{iFile} < BinSizes_(iBin+1));
                    temp = temp(:,idx);
                    if ~isempty(temp)
                        temp = zscore(temp,[],1);
                        temp = nanmean(temp,2);
                        AssCat.D_Delay.fracSig{s}{iDelay}(iFile,:) =  temp;
                    end
                end
            end
            temp = AssCat.D_Delay.fracSig{s}{iDelay};
            %         temp = smooth2a(temp,0,10);
            AssCat.D_Delay.fracSigmean{s}{iDelay} = nanmean(temp);
            AssCat.D_Delay.fracSigSEM{s}{iDelay}  = nansem(temp);
        end
    end
    
    figure('name', sprintf('%gs~%gs duration',BinSizes_(iBin),BinSizes_(iBin+1)));
    for iDelay = 1:3
        tb = tlimsSample(1):bw:Delays(iDelay);
        subplot(length(Delays),1,iDelay); hold on
        for s =1:3
            m = AssCat.D_Delay.fracSigmean{s}{iDelay};
            e = AssCat.D_Delay.fracSigSEM{s}{iDelay};
            ciplot(m+e,m-e,tb,color_{s},0.6);
        end
        plot([0 0],[0 5],'g','LineWidth',1.5)
        for iDelay =1:3
            plot([Delays(iDelay) Delays(iDelay)],[0 5],':k')
        end
        axis([-5 16 -2 2])
    end
    
end
end
%%
% Compares continuous L/R discrimination during the delay period,
% calcluated on short/medium/long delay trials independently
%% %%%%%% PREAMBLE %%%%%%

clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
DelaysDur = [4,8,16];

Target = 'LONG';

tlimsAll = [-30 30];
tlimsTrials = [-5 5];
tlimsShort=[-5 6];
tlimsMedium=[-5 10];
tlimsLong=[-5 20];

plotOnline = false;
bw=0.05;
Nbs = 500;

minFR = 0.1;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
%     pat = '/Volumes/HDD2/DNMTP/raw/';
    pat = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw/';
    home = getenv('HOME');
else
    
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    addpath(genpath('/panfs/panasas01/phph/ad15419/MATLAB'))
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'SpikeTrainCoherence' filesep '*' Target '*coherenceMean.mat']);

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

UnitSelection = 'pairs'; %{'pairs','groups'}
groupSize = 8;
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

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
EventList = {'CueLight','SamplePress','DelayEnd','NosePoke','ChoicePress','RewardConsume'};
Direction = {'Left','Right'};
Outcome   = {'Correct','Error'};
%% Batch process units
for iFile = 1:length(fileList)
   
    %% Get the files
    fname=strtok(fileList(iFile).name,'_');
    fprintf('Analysing run %d/%d %s...\n',iFile,length(fileList),fname)

    load(fullfile(pat,'SpikeTrainCoherence',sprintf('%s_coherenceMean.mat',fname)));
    load(sprintf('%s%s%s.mat',pat,filesep,fname),'PFCcells','HPcells');
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu); clear PFCcells HPcells
    %% Get the assembly membership info
           
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
            Ass.usel_out{3} = [Ass.usel_out{1},Ass.usel_out{2}+nu(1)];
            Ass.units       = units;
            Ass.nu = nu; Ass.nu(3) = sum(Ass.nu(1:2));
            % Check that inter-area assemblies actually span the two areas
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
            for iArea_ = 1:3
                Ass.LocalMembers{iArea_}    = unique(cell2mat(Ass.units{iArea_}));
            end
            
            
            Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3}<=Ass.nu(1)));
            Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3}>Ass.nu(1))-Ass.nu(1));
            for iArea_ = 1:2
                Ass.JointMembers{iArea_} = setdiff(unique(cell2mat(Ass.units{iArea_})),Ass.LocalMembers{iArea_});
                Ass.NonMembers{iArea_}   = setdiff(1:Ass.nu(iArea_),[Ass.LocalMembers{iArea_},Ass.JointMembers{iArea_}]);
            end
            for iArea_ = 1:2
                Ass.NonMembers{iArea_}   = Ass.usel_out{iArea_}(Ass.NonMembers{iArea_});
                Ass.LocalMembers{iArea_} = Ass.usel_out{iArea_}(Ass.LocalMembers{iArea_});
                Ass.JointMembers{iArea_} = Ass.usel_out{iArea_}(Ass.JointMembers{iArea_});
            end
            
            for iArea_ = 1:2
                %Membership_{iArea_} = -ones( Ass.nu(iArea_),1);
                Membership_{iArea_} = -ones(nu(iArea_),1);
                Membership_{iArea_}(Ass.NonMembers{iArea_})=0;
                Membership_{iArea_}(Ass.LocalMembers{iArea_})=1;
                Membership_{iArea_}(Ass.JointMembers{iArea_})=2;
                Membership_{iArea_}(setdiff(1:nu(iArea_),Ass.usel_out{iArea_}))=[];
            end

            for iArea_ = 1:2
                % Pad this out to ensure that there's an entry even if there's no detected assembly
                LocalMembersMatrix_{iArea_} = nan(Ass.nu(iArea_),length(Ass.units{iArea_})+1);
                if ~isempty(Ass.units{iArea_})
                    for iAss=1:length(Ass.units{iArea_})
                        LocalMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_),Ass.units{iArea_}{iAss});
                    end
                end
                
                JointMembersMatrix_{iArea_} = nan(Ass.nu(iArea_),length(Ass.units{3})+1);
                for iAss=1:length(Ass.units{3})
                    if ~isempty(Ass.units{3})
                        units_ = Ass.units{3}{iAss};
                        if iArea_==1
                            JointMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_), units_(units_<=Ass.nu(1)));
                        else
                            JointMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_), units_(units_>Ass.nu(1))-Ass.nu(1));
                        end
                    end
                end
            end
            iArea_ = 3;
            JointMembersMatrix_{iArea_} = nan(Ass.nu(iArea_),length(Ass.units{3})+1);
            for iAss=1:length(Ass.units{3})
                if ~isempty(Ass.units{3})
                    units_ = Ass.units{3}{iAss};
                    JointMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_), units_);
                end
            end
            clear units usel_out nu FSCsel Tmtx  iFR noPFC idx nAss U_ iAss
            nu = cellfun(@length,Membership_);nu(3) = sum(nu);
    %% (1) Collect member classes
    nonMemberUnits = cell(3,1);
    for iArea =1:2
        nonMemberUnits{iArea} = find(Membership_{iArea}==0);
    end
    nonMemberUnits{2} = nonMemberUnits{2}+nu(1);
    nonMemberUnits{3} = find([Membership_{1};Membership_{2}]==0);
     
    %non-member pairs
    nonMemberUnitsLookup = cell(3,1);
    for iArea =1:3
        units_ = nonMemberUnits{iArea};
        if iArea<3
            for i = 1:length(units_)
                for j = i+1:length(units_)
                    nonMemberUnitsLookup{iArea}=[nonMemberUnitsLookup{iArea}; units_(i),units_(j)];
                end
            end
        else
            unitsJoint{1} = units_(units_<=nu(1));
            unitsJoint{2} = units_(units_>nu(1));
            for i = 1:length(unitsJoint{1})
                for j = 1:length(unitsJoint{2})
                    nonMemberUnitsLookup{iArea}=[nonMemberUnitsLookup{iArea}; unitsJoint{1}(i),unitsJoint{2}(j)];
                end
            end
        end
    end
    
    % Within-assembly pairs
    localMemberUnitsLookup = cell(2,1);
    for iArea =1:2
        for iAss = 2:size(LocalMembersMatrix_{iArea},2)
            units_ = find(LocalMembersMatrix_{iArea}(:,iAss));
            if iArea==2 units_=units_ + nu(1);end
            for i = 1:length(units_)
                for j=i+1:length(units_)
                    localMemberUnitsLookup{iArea}=[localMemberUnitsLookup{iArea}; units_(i),units_(j)];
                end
            end
            
        end
    end
    
    jointMemberUnitsLookup = cell(3,1);
    for iArea =1:3
        for iAss = 2:size(JointMembersMatrix_{iArea},2)
            units_ = find(JointMembersMatrix_{iArea}(:,iAss));
            if iArea==2 units_=units_ + nu(1);end
            if iArea<3
                for i = 1:length(units_)
                    for j=i+1:length(units_)
                        jointMemberUnitsLookup{iArea}=[jointMemberUnitsLookup{iArea}; units_(i),units_(j)];
                    end
                end
            else
                unitsJoint{1} = units_(units_<=nu(1));
                unitsJoint{2} = units_(units_>nu(1));
                 for i = 1:length(unitsJoint{1})
                    for j=1:length(unitsJoint{2})
                        jointMemberUnitsLookup{iArea}=[jointMemberUnitsLookup{iArea}; unitsJoint{1}(i),unitsJoint{2}(j)];
                    end
                end
            end
        end
    end
    
    % Cross-assembly pairs
    localMemberUnitsOtherLookup = cell(2,1);
    for iArea =1:2
        Allmembers_ = [];
        for iAss = 2:size(LocalMembersMatrix_{iArea},2)
            Allmembers_ = [ Allmembers_;  find(LocalMembersMatrix_{iArea}(:,iAss))] ;
        end
        Allmembers_ = unique(Allmembers_);
        
        for iAss = 2:size(LocalMembersMatrix_{iArea},2)
            units_ = find(LocalMembersMatrix_{iArea}(:,iAss));
            otherunits_ = setdiff(Allmembers_,units_);
            for i = 1:length(units_)
                for j=1:length(otherunits_)
                    localMemberUnitsOtherLookup{iArea}=[localMemberUnitsOtherLookup{iArea}; units_(i),otherunits_(j)];
                end
            end
        end
    end
    
    jointMemberUnitsOtherLookup = cell(3,1);
    for iArea =1:3
        
        Allmembers_ = [];
        for iAss = 2:size(JointMembersMatrix_{iArea},2)
            Allmembers_ = [ Allmembers_;  find(JointMembersMatrix_{iArea}(:,iAss))] ;
        end
        Allmembers_ = unique(Allmembers_);
        
        for iAss = 2:size(JointMembersMatrix_{iArea},2)
            
            units_ = find(JointMembersMatrix_{iArea}(:,iAss));
            otherunits_ = setdiff(Allmembers_,units_);
            
            if iArea==2 units_=units_ + nu(1);end
            if iArea <3 
            for i = 1:length(units_)
                for j=1:length(otherunits_)
                    jointMemberUnitsOtherLookup{iArea}=[jointMemberUnitsOtherLookup{iArea}; units_(i),otherunits_(j)];
                end
            end
            else
                for i = 1:length(units_)
                    for j=1:length(otherunits_)
                        if (units_(i)<=nu(1) && otherunits_(j)>nu(1)) || (otherunits_(j)<=nu(1) && units_(i)>nu(1))
                            jointMemberUnitsOtherLookup{iArea}=[jointMemberUnitsOtherLookup{iArea}; units_(i),otherunits_(j)];
                        end
                    end
                end
                
            end
            
        end
    end
    clear units_ unitsJoint Allmembers_ otherunits_
    %% (2) Loop over member classes to extract correct entries from coherence matrix
    %     clear Archive
    MemberClasses = {'nonMemberUnitsLookup','localMemberUnitsLookup','jointMemberUnitsLookup','localMemberUnitsOtherLookup','jointMemberUnitsOtherLookup';...
        'nonMembers','localMembers','jointMembers','localMembersCross','jointMembersCross';...
        3,2,3,2,3};
    for iClass =1:size(MemberClasses,2)
        lookup_ = eval(MemberClasses{1,iClass});
        
        for iArea = 1:MemberClasses{3,iClass}
            
            if ~isempty(lookup_{iArea})
                
                idx = sub2ind(size(coherence_mean),lookup_{iArea}(:,1),lookup_{iArea}(:,2));
                temp = MeanRate(idx);
                idx(cell2mat(cellfun(@min, temp,'UniformOutput',false))<minFR)=[];
                %idx(min([],[],1)<minFR)=[];
                collapsed_ = coherence_mean(idx);
                collapsed_ = collapsed_(~isempty_cell(collapsed_));
                
                for i=1:length(collapsed_) %Loop across unit pairs
                    
                    
                    for iDelay = 1:length(Delays_)
                        
                        for iOutcome = 1:2 %Correct | Error
                            if iOutcome==1
                                nEvents = length(EventList);
                            else
                                nEvents = length(EventList)-1;
                            end
                            for iEvent = 1:nEvents  % Loop across events
                                
                                evt_ = sprintf('%s.%s.%s',Delays_{iDelay},EventList{iEvent},Outcome{iOutcome});
                                temp = eval(sprintf('collapsed_{i}.%s',evt_));
                                if isempty(temp{1})
                                    temp = temp{2};
                                elseif isempty(temp{2})
                                    temp = temp{1};
                                else
                                    try
                                        temp = (temp{1}+temp{2})./2;
                                    catch 
                                        warning('Cutout size mismatch!')
                                        sz = cellfun(@(x) size(x,1),temp);
                                        if sz(1)>sz(2)
                                            temp = temp{1};
                                        else
                                            temp = temp{2};
                                        end
                                    end
                                end
                                
                                eval(sprintf('Archive.%s.%s{iArea}{iFile,1}(1:size(temp,1),1:size(temp,2),i) = temp;',MemberClasses{2,iClass},evt_))
                                
                            end
                            evt_ = sprintf('%s.%s.%s',Delays_{iDelay},'DelayPeriod',Outcome{iOutcome});
                            temp = eval(sprintf('collapsed_{i}.%s',evt_));
                            if isempty(temp{1})
                                temp = temp{2};
                            elseif isempty(temp{2})
                                temp = temp{1};
                            else
                                try
                                    temp = (temp{1}+temp{2})./2;
                                catch
                                    warning('Cutout size mismatch!')
                                    sz = cellfun(@(x) size(x,1),temp);
                                    if sz(1)>sz(2)
                                        temp = temp{1};
                                    else
                                        temp = temp{2};
                                    end
                                end
                            end
                            eval(sprintf('Archive.%s.%s{iArea}{iFile,1}(1:size(temp,1),1:size(temp,2),i) = temp;',MemberClasses{2,iClass},evt_))
                            
                        end
                    end
                end
            end
            %         temp = nanmean(Archive.nonmembers.Long.DelayPeriod.Correct{1},3);
            %         imagesc(1:size(temp,1),f_,temp');set(gca, 'YDir','normal')
        end
    end
end

%% Collapse and plot - [ sample | choice ]
for iClass =1:3%size(MemberClasses,2)
    figure('name',MemberClasses{2,iClass}); hold on
    for iArea = 1:MemberClasses{3,iClass}
        iDelay = 3;
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}',MemberClasses{2,iClass},Delays_{iDelay},'SamplePress'))
        eval(sprintf('B = Archive.%s.%s.%s.Correct{iArea}',MemberClasses{2,iClass},Delays_{iDelay},'ChoicePress'))
        idx = isempty_cell(A) | isempty_cell(B);
        A(idx)=[]; B(idx)=[];
%         A = (cellfun(@(x) nanmean(x,3),A,'UniformOutput',false));
%         B = (cellfun(@(x) nanmean(x,3),B,'UniformOutput',false));
        

        temp = [];
        for i=1:length(A)
            try
                temp = cat(3,temp,[A{i};B{i}]);
            end
        end

        temp_mean = nanmean(temp,3);
        temp_mean = smooth2a(temp_mean,0,0);
        subplot(1,MemberClasses{3,iClass},iArea); title(Areas{iArea})
        imagesc(1:size(temp_mean,1),f_,temp_mean');set(gca, 'YDir','normal')
        ylim([0.5 100]);
%         caxis([0.4 0.5])
        
        
    end
end
%% Collapse and plot - loop over events
for iClass =1:3%size(MemberClasses,2)
    for iArea = 1:MemberClasses{3,iClass}
        
        figure('name',sprintf('%s: %s',Areas{iArea},MemberClasses{2,iClass}));
          iDelay = 2;
            for iOutcome = 1:2 %Correct | Error
                if iOutcome==1
                    nEvents = length(EventList);
                else
                    nEvents = length(EventList)-1;
                end
                for iEvent = 1:nEvents  % Loop across events
                    evt_ = sprintf('%s.%s.%s',Delays_{iDelay},EventList{iEvent},Outcome{iOutcome});
                    eval(sprintf('A = Archive.%s.%s{iArea}',MemberClasses{2,iClass},evt_))
                    idx = isempty_cell(A);
                    A(idx)=[]; 
                    
                    
                    
                    % Across animal average
%                     temp = cellfun(@(x) nanmean(x,3),A,'UniformOutput',0);
%                     temp2 = [];
%                     for i=1:length(temp)
%                         if numel(temp{i})>1 && size(temp{i},1)>1
%                             temp2 = cat(3,temp2,temp{i});
%                         end
%                     end
%                     temp_mean = nanmean(temp2,3);
%                     
                    % Across all pairs average
                    temp = [];
                    for i=1:length(A)
                        if numel(A{i})>1 && size(A{i},1)>1
                            temp = cat(3,temp,A{i});
                        end
                    end
                    temp_mean = nanmean(temp,3);
                    
                    
                    temp_mean = smooth2a(temp_mean,0,4);                    
                    if iOutcome==1
                        subplot(2,length(EventList),iEvent);
                    else
                        subplot(2,length(EventList),iEvent+length(EventList));

                    end
                    
                    imagesc(1:size(temp_mean,1),f_,temp_mean');set(gca, 'YDir','normal')
                    ylim([1 90]);
%                     caxis([0.2 0.4])
                    title(sprintf('%s(%s)',EventList{iEvent},Outcome{iOutcome}))
                end
            end
    end
end
%% Collapse and plot - Delay
for iClass =1:3%size(MemberClasses,2)
    figure('name',MemberClasses{2,iClass}); hold on
    for iArea = 1:MemberClasses{3,iClass}
        iDelay = 3;
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}',MemberClasses{2,iClass},Delays_{iDelay},'DelayPeriod'))
        idx = isempty_cell(A) ;
        A(idx)=[];
        
        % Across animal average
        temp = cellfun(@(x) nanmean(x,3),A,'UniformOutput',0);
        temp2 = [];
        for i=1:length(temp)
            temp2 = cat(3,temp2,temp{i});
        end
        temp_mean = nanmean(temp2,3);
        
        % Across all pairs average
%         temp = [];
%         for i=1:length(A)
%             temp = cat(3,temp,A{i});
%         end
%         temp_mean = nanmean(temp,3);


        %         temp_mean = smooth2a(temp_mean,1,1);
        subplot(1,MemberClasses{3,iClass},iArea); title(Areas{iArea})
        imagesc(1:size(temp_mean,1),f_,temp_mean');set(gca, 'YDir','normal')
        ylim([0.5 100]);
        caxis([0.4 0.5])
        
        
    end
end

%% Collapse and plot - [ sample | choice ]

figure('name',MemberClasses{2,iClass}); hold on
iArea =3;
temp = [];

for iDelay = 1:3
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}','nonMembers',Delays_{iDelay},'SamplePress'))
        eval(sprintf('B = Archive.%s.%s.%s.Correct{iArea}','nonMembers',Delays_{iDelay},'ChoicePress'))
        idx = isempty_cell(A) | isempty_cell(B);
        A(idx)=[];
%         A = (cellfun(@(x) nanmean(x,3),A,'UniformOutput',false));
       

        for i=1:length(A)
            try
%                 temp = cat(3,temp,[nanmean(A{i})]);
                temp = cat(3,temp,[(A{i})]);
            end
        end
end
temp= squeeze(temp);

%         temp_mean = nanmean(temp,3);
        temp_mean = nanmean(temp,1);

        temp_mean = smooth2a(temp_mean,0,0);
        m_ = nanmean(temp_mean,2);
        e_ = nansem(temp_mean,2);
        ciplot(m_+e_,m_-e_,f_,col_{1})
        
temp = [];

for iDelay = 1:3;
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}','jointMembers',Delays_{iDelay},'SamplePress'))
        idx = isempty_cell(A);
        A(idx)=[];
%         A = (cellfun(@(x) nanmean(x,3),A,'UniformOutput',false));
       

        for i=1:length(A)
            try
%                 temp = cat(3,temp,[nanmean(A{i})]);
                temp = cat(3,temp,[(A{i})]);
            end
        end
        
end
temp= squeeze(temp);
        temp_mean = nanmean(temp,3);
        temp_mean = smooth2a(temp_mean,0,0);
        m_ = nanmean(temp_mean,2);
        e_ = nansem(temp_mean,2);
        ciplot(m_+e_,m_-e_,f_,col_{3},0.9)        
        
temp = [];

for iDelay = 1:3;
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}','jointMembersCross',Delays_{iDelay},'SamplePress'))
        idx = isempty_cell(A);
         A(idx)=[];
%         A = (cellfun(@(x) nanmean(x,3),A,'UniformOutput',false));
       

        for i=1:length(A)
            try
%                 temp = cat(3,temp,[nanmean(A{i})]);
                   temp = cat(3,temp,[(A{i})]);
            end
        end
   
end
temp= squeeze(temp);
        temp_mean = nanmean(temp,3);
        temp_mean = smooth2a(temp_mean,0,0);
        m_ = nanmean(temp_mean,2);
        e_ = nansem(temp_mean,2);
        ciplot(m_+e_,m_-e_,f_,col_{3},0.2)        
        
        
axis([0.5 30 -Inf Inf]);

%% Collapse and plot - Sample collapse all
figure('name',MemberClasses{2,iClass}); hold on
iArea = 3;
temp = [];
iF = [0,f_]>3.5 & [0,f_]<5.5;

for iDelay = 1:3
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}','nonMembers',Delays_{iDelay},'SamplePress'))
        idx = isempty_cell(A);
        A(idx)=[];
        A = (cellfun(@(x) squeeze(nanmean(x,1)),A,'UniformOutput',false));
        idx =cell2mat(cellfun(@(x) size(x,2),A,'UniformOutput',false))==1;
        A(idx)=[];
        for i=1:length(A)
            try
%                 temp = cat(3,temp,nanmean(A{i},1));
                temp = cat(2,temp,[(A{i})]);
            end
        end
end


%         temp_mean = nanmean(temp,3);

%         temp_mean = smooth2a(temp_mean,0,0);
        coh_{1} = nanmean(temp(iF,:),2);

        m_ = nanmean(temp,2);
        e_ = nansem(temp,2);
        ciplot(m_+e_,m_-e_,[0,f_],col_{1})
        
temp = [];

for iDelay = 1:3
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}','jointMembers',Delays_{iDelay},'SamplePress'))
        idx = isempty_cell(A);
        A(idx)=[];
        A = (cellfun(@(x) squeeze(nanmean(x,1)),A,'UniformOutput',false));
        idx =cell2mat(cellfun(@(x) size(x,2),A,'UniformOutput',false))==1;
        A(idx)=[];
        for i=1:length(A)
            try
%                 temp = cat(3,temp,nanmean(A{i},1));
                temp = cat(2,temp,[(A{i})]);
            end
        end
end

temp= squeeze(temp);
        coh_{2} = nanmean(temp(iF,:),2);

        temp_mean = nanmean(temp,3);
        temp_mean = smooth2a(temp_mean,0,0);
        m_ = nanmean(temp_mean,2);
        e_ = nansem(temp_mean,2);
        ciplot(m_+e_,m_-e_,[0,f_],col_{3},0.9)        
        
temp = [];

for iDelay = 1:3
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}','jointMembersCross',Delays_{iDelay},'SamplePress'))
        idx = isempty_cell(A);
        A(idx)=[];
        A = (cellfun(@(x) squeeze(nanmean(x,1)),A,'UniformOutput',false));
        idx =cell2mat(cellfun(@(x) size(x,2),A,'UniformOutput',false))==1;
        A(idx)=[];
        for i=1:length(A)
            try
%                 temp = cat(3,temp,nanmean(A{i},1));
                temp = cat(2,temp,[(A{i})]);
            end
        end
end


temp= squeeze(temp);
        coh_{3} = nanmean(temp(iF,:),2);

        temp_mean = nanmean(temp,3);
        temp_mean = smooth2a(temp_mean,0,0);
        m_ = nanmean(temp_mean,2);
        e_ = nansem(temp_mean,2);
        ciplot(m_+e_,m_-e_,[0,f_],col_{3},0.2)        
        
        
axis([0.5 30 0.29 0.31]);
[p,tbl,stats]=anova1(cell2mat(coh_'),[ones(length(coh_{1}),1);2*ones(length(coh_{2}),1);3*ones(length(coh_{3}),1)])
multcompare(stats)

%% Collapse and plot - Sample animal averages
figure('name',MemberClasses{2,iClass}); hold on
iArea = 3;
iF = [0,f_]>3.5 & [0,f_]<5.5;
temp = [];
for iDelay = 1:3
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}','nonMembers',Delays_{iDelay},'SamplePress'))
        idx = isempty_cell(A);
        A(idx)=[];
        idx =cell2mat(cellfun(@(x) size(x,1),A,'UniformOutput',false))==1;
        A(idx)=[];
        
        A = (cellfun(@(x) nanmean(squeeze(nanmean(x,1)),2),A,'UniformOutput',false));

        for i=1:length(A)
            try
%                 temp = cat(3,temp,nanmean(A{i},1));
                temp = cat(2,temp,[(A{i})]);
            end
        end
end


%         temp_mean = nanmean(temp,3);

%         temp_mean = smooth2a(temp_mean,0,0);
        coh_{1} = nanmean(temp(iF,:),2);
        m_ = nanmean(temp,2);
        e_ = nansem(temp,2);
        ciplot(m_+e_,m_-e_,[0,f_],col_{1},0.9)
        
temp = [];

temp = [];
for iDelay = 1:3
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}','jointMembers',Delays_{iDelay},'SamplePress'))
        idx = isempty_cell(A);
        A(idx)=[];
        idx =cell2mat(cellfun(@(x) size(x,1),A,'UniformOutput',false))==1;
        A(idx)=[];
        
        A = (cellfun(@(x) nanmean(squeeze(nanmean(x,1)),2),A,'UniformOutput',false));

        for i=1:length(A)
            try
%                 temp = cat(3,temp,nanmean(A{i},1));
                temp = cat(2,temp,[(A{i})]);
            end
        end
end


%         temp_mean = nanmean(temp,3);

%         temp_mean = smooth2a(temp_mean,0,0);
        coh_{2} = nanmean(temp(iF,:),2)
        m_ = nanmean(temp,2);
        e_ = nansem(temp,2);
        ciplot(m_+e_,m_-e_,[0,f_],col_{3},0.9)
        
temp = [];

temp = [];
for iDelay = 1:3
        
        eval(sprintf('A = Archive.%s.%s.%s.Correct{iArea}','jointMembersCross',Delays_{iDelay},'SamplePress'))
        idx = isempty_cell(A);
        A(idx)=[];
        idx =cell2mat(cellfun(@(x) size(x,1),A,'UniformOutput',false))==1;
        A(idx)=[];
        
        A = (cellfun(@(x) nanmean(squeeze(nanmean(x,1)),2),A,'UniformOutput',false));

        for i=1:length(A)
            try
%                 temp = cat(3,temp,nanmean(A{i},1));
                temp = cat(2,temp,[(A{i})]);
            end
        end
end


%         temp_mean = nanmean(temp,3);

%         temp_mean = smooth2a(temp_mean,0,0);
        coh_{3} = nanmean(temp(iF,:),2)
        m_ = nanmean(temp,2);
        e_ = nansem(temp,2);
        ciplot(m_+e_,m_-e_,[0,f_],col_{3},0.3)

        
axis([0.5 30 0.29 0.31]);


[p,tbl,stats]=anova1(cell2mat(coh_'),[ones(length(coh_{1}),1);2*ones(length(coh_{2}),1);3*ones(length(coh_{3}),1)])
multcompare(stats)

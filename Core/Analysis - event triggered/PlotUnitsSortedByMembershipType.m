%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Target = 'LONG';
Areas = {'PFC','HP','Joint'};
TimeSpan = 4;

bw=0.05;
tlimsAll = [-5 5];
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

maxRange=20;    % Upper limit of assembly size to explore
noTrials = 5;   % How many trials to consider simultanrously
nReps = 50;     % how many repeated draws among trials to perform
        
        
plotOnline = false;
useWholeTaskPeriod = true;
ResampleTrials = true;
InfoCriteria = 'max'; % 'mean','max'


clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);


% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); % AD inverted logical flag for testing
end
if sum(name_flag)>0
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag
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
%%
for iFile =1:length(fileList)
    %% Get spike times and rates, behaviour
    fn = fullfile(pat,strtok(fileList(iFile).name,'_'));
    spikes = load(fn);
    if useWholeTaskPeriod
        fn = fullfile(pat,'KDE_binsTaskOnly',[strtok(fileList(iFile).name,'_'),'_PFC_iFR50_behavOnly.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{1} = temp.iFR;clear temp
        fn = fullfile(pat,'KDE_binsTaskOnly',[strtok(fileList(iFile).name,'_'),'_HP_iFR50_behavOnly.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{2} = (temp.iFR);
        Tmtx=temp.Tmtx;clear temp
    else
        fn = fullfile(pat,'KDE_bins',[strtok(fileList(iFile).name,'_'),'_PFC_iFR50.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{1} = temp.iFR;clear temp
        fn = fullfile(pat,'KDE_bins',[strtok(fileList(iFile).name,'_'),'_HP_iFR50.mat']);
        temp = load(fn,'iFR','Tmtx');
        iFR{2} = (temp.iFR);
        Tmtx=temp.Tmtx; clear temp
    end
    iFR{3} = [iFR{1},iFR{2}];
    
    fname=strtok(fileList(iFile).name,'_')
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    %% Aggregate all trials regardless of delay length
    LeftTrials=[];RightTrials=[];
    for iDelay =2:length(Delays_)
%         eval(sprintf('LeftTrials = [LeftTrials;[t.%s.SamplePress_LeftCorrect,t.%s.ChoicePress_LeftCorrect]];',Delays_{iDelay},Delays_{iDelay}));
%         eval(sprintf('RightTrials = [RightTrials;[t.%s.SamplePress_RightCorrect,t.%s.ChoicePress_RightCorrect]];',Delays_{iDelay},Delays_{iDelay}));
        
        eval(sprintf('LeftTrials = [LeftTrials;[t.%s.SamplePress_LeftCorrect,t.%s.NosePoke_LeftCorrect]];',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf('RightTrials = [RightTrials;[t.%s.SamplePress_RightCorrect,t.%s.NosePoke_RightCorrect]];',Delays_{iDelay},Delays_{iDelay}));
        
    end
    % Run decoder once on all trials to rank order the neurons for later selection
  
    for s=1:2
        Ltrials_{s} = []; Rtrials_{s} = [];
        Ltrials = [];nL = 0;
        for iTrial =1:size(LeftTrials,1)
            try
                tlims_  = [LeftTrials(iTrial,1)/1e6;LeftTrials(iTrial,2)/1e6]+tlimsAll(1);
                tlims_  = closest(Tmtx,tlims_);
                tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                Ltrials = [Ltrials;iFR{s}(tlims_,:)];
                Ltrials_{s}(:,:,nL+1) = iFR{s}(tlims_,:);
                nL=nL+1;
            end
        end
        
        Rtrials = [];nR = 0; 
        for iTrial =1:size(RightTrials,1)
            try
                tlims_  = [RightTrials(iTrial,1)/1e6;RightTrials(iTrial,2)/1e6]+tlimsAll(1);
                tlims_  = closest(Tmtx,tlims_);
                tlims_  = [tlims_(1):(tlims_(1)+length(tbAll)-1), tlims_(2):(tlims_(2)+length(tbAll)-1)];
                Rtrials = [Rtrials;iFR{s}(tlims_,:)];
                Rtrials_{s}(:,:,nR+1) = iFR{s}(tlims_,:);
                nR=nR+1;
            end
        end
        
        FR{s} = [Ltrials;Rtrials];
        evt0 = [ones(nL,1);2*ones(nR,1)];
        Ltr = 2*length(tbAll);
    end
    Ltrials_{3}=cat(2,Ltrials_{1},Ltrials_{2});
    Rtrials_{3}=cat(2,Rtrials_{1},Rtrials_{2});
   
    clear Ltrials Rtrials iTrial
    %% Get assembly information
    
    if useWholeTaskPeriod
        fn = fullfile(pat,'KDE_binsTaskOnly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_FSC.mat',strtok(fileList(iFile).name,'_')));
        A  = load(fn);
        fn = fullfile(pat,'KDE_binsTaskOnly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_AssemRes2.mat',strtok(fileList(iFile).name,'_')));
        B  = load(fn);
        % recalculate included units
        fn = fullfile(pat,'KDE_binsTaskOnly',sprintf('%s_PFC_iFR50_behavOnly.mat',strtok(fileList(iFile).name,'_')));
        usel_out=SelCells(fn,0.1,1e6);
    else
        fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',strtok(fileList(iFile).name,'_')));
        A  = load(fn);
        fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',strtok(fileList(iFile).name,'_')));
        B  = load(fn);
        % recalculate included units
        fn = fullfile(pat,'KDE_bins',sprintf('%s_PFC_iFR50.mat',strtok(fileList(iFile).name,'_')));
        [~,~,~,~,B.usel_out]=SelTrialsCellsWholeTrial(fn,10,0.1,1e6,'iFR');
    end
    
    spikes.jointcells=[spikes.PFCcells;spikes.HPcells]';
    A.nu(3) = sum(A.nu(1:2));
    %     B.usel_out{3}= [B.usel_out{1},B.usel_out{2}+max(B.usel_out{1})];
    B.usel_out{3}= [B.usel_out{1},B.usel_out{2}+size(iFR{1},2)];
    %% Sort units by assembly membership types
        types_={'Nonmember','Localmember','Jointmember'};
        for s = 1:2
            
            % Global = joint + local
            % Nomember = ~Global
            switch s
                case 1
                    joint_ = unique(cell2mat(A.units{3})); joint_(joint_>A.nu(1))=[];
                case 2
                    joint_ = unique(cell2mat(A.units{3})); joint_ = joint_(joint_>A.nu(1))-A.nu(1);
            end
            local_ = unique(cell2mat(A.units{s}));
            all_   = 1:A.nu(s);
            LocalmemberUnits    = setdiff(local_,joint_);                      % local only members
            JointmemberUnits    = joint_;                                      % inter-area only members
            GlobalmemberUnits   = unique([JointmemberUnits,LocalmemberUnits]); % membership of any assembly
            NonmemberUnits      = setdiff(all_,GlobalmemberUnits);             % global nonmembers
            
            LocalmemberUnits    = B.usel_out{s}(LocalmemberUnits);
            NonmemberUnits      = B.usel_out{s}(NonmemberUnits);
            JointmemberUnits    = B.usel_out{s}(JointmemberUnits);
            GlobalmemberUnits   = B.usel_out{s}(GlobalmemberUnits);
            clear joint_ local_ all_
            
            D.fname{iFile} = fname;
            D.nonmemberUnits{iFile}{s} = NonmemberUnits;
            D.localmemberUnits{iFile}{s} = LocalmemberUnits;
            D.jointmemberUnits{iFile}{s} = JointmemberUnits;
            
        end %/s
    %% plot trial-averaged firing rates separated by member class
    figure;
    
    col_ = {[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
    tb   = (1:length(tbAll)*2)*bw;

    
    
    for s=1:2
        
        subplot(1,2,s); hold on

        NonmemberUnits   = D.nonmemberUnits{iFile}{s};
        LocalmemberUnits = D.localmemberUnits{iFile}{s};
        JointmemberUnits = D.jointmemberUnits{iFile}{s};
        L_ = [1, length(NonmemberUnits)];
        L_ = [L_; [max(max(L_))+1, length(LocalmemberUnits)+max(max(L_))]];
        L_ = [L_; [max(max(L_))+1, length(JointmemberUnits)+max(max(L_))]];
        
        
        NonmemberUnits_m   = cat(1,nanmean(zscore(Ltrials_{s}(:,NonmemberUnits,:),[],1),3),nanmean(zscore(Rtrials_{s}(:,NonmemberUnits,:),[],1),3));
        LocalmemberUnits_m = cat(1,nanmean(zscore(Ltrials_{s}(:,LocalmemberUnits,:),[],1),3),nanmean(zscore(Rtrials_{s}(:,LocalmemberUnits,:),[],1),3));
        JointmemberUnits_m = cat(1,nanmean(zscore(Ltrials_{s}(:,JointmemberUnits,:),[],3),3),nanmean(zscore(Rtrials_{s}(:,JointmemberUnits,:),[],1),3));
        
        NonmemberUnits_s   = cat(1,nansem(zscore(Ltrials_{s}(:,NonmemberUnits,:),[],1),3),nansem(zscore(Rtrials_{s}(:,NonmemberUnits,:),[],1),3));
        LocalmemberUnits_s = cat(1,nansem(zscore(Ltrials_{s}(:,LocalmemberUnits,:),[],1),3),nansem(zscore(Rtrials_{s}(:,LocalmemberUnits,:),[],1),3));
        JointmemberUnits_s= cat(1,nansem(zscore(Ltrials_{s}(:,JointmemberUnits,:),[],1),3),nansem(zscore(Rtrials_{s}(:,JointmemberUnits,:),[],1),3));
        
       
        
        x = [NonmemberUnits_m,LocalmemberUnits_m,JointmemberUnits_m];
        x = staggerplot(x,0,2);
        
        y = [NonmemberUnits_s,LocalmemberUnits_s,JointmemberUnits_s];
        y = staggerplot(y,0,0);
        
        for i=1:3
%             plot(tb,x(1:length(tb),L_(i,1):L_(i,2)),'color',col_{i})
%             plot(tb,x(length(tb)+(1:length(tb)),L_(i,1):L_(i,2)),'color',col_{i})
            for j=L_(i,1) : L_(i,2) 
                ciplot(x(1:length(tb),j) + y(1:length(tb),j),...
                       x(1:length(tb),j) - y(1:length(tb),j),...
                       tb,col_{i})
                   
                ciplot(x(length(tb)+(1:length(tb)),j) + y(length(tb)+(1:length(tb)),j),...
                       x(length(tb)+(1:length(tb)),j) - y(length(tb)+(1:length(tb)),j),...
                       tb,col_{i})                   
            end
           
        end
        
        
        
    end

end
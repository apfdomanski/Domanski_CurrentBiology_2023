% %%%%%% PREAMBLE %%%%%%
clear                              

Delays_ = {'Short','Medium','Long'};
Target = 'LONG';
Areas = {'PFC','HP','Joint'};

bw=0.05; tlimsAll = [-5 5]; tbAll = tlimsAll(1):bw:tlimsAll(2);

maxRange=20;    % Upper limit of assembly size to explore
noTrials = 5;   % How many trials to consider simultanrously
nReps = 50;     % how many repeated draws among trials to perform


InfoCriteria = 'max'; % 'mean','max'

if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
    pat = '/Volumes/Aleks Data Portable/AssemblyAnalysis/raw/';
elseif isunix
    pat ='/Volumes/Data/DNMTP/raw';
end

fileList = dir([pat 'SyntheticAssemblies' filesep '*' Target '*AssemblyComparison_redux2.mat']);


%% Batch import
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
for iFile =1:length(fileList)
    
    fnIn = fullfile(pat, 'SyntheticAssemblies', fileList(iFile).name);
    fnIn_ = strsplit(char(fileList(iFile).name),'2.mat');fnIn_=fnIn_{1};
    try
        D{iFile} = load([fnIn_,'3.mat'],'D_');
    catch
        D{iFile} = load(fnIn,'D_');
    end
   
end
clear temp 
%% Recalculate score
for iFile =1:length(D)
    for s=1:3
        D{iFile}.D_.AssBestSUs.All.Score{s} = nanmax(D{iFile}.D_.AssBestSUs.All.CVE{s});
        D{iFile}.D_.AssBestSUs.Members.Score{s} = nanmax(D{iFile}.D_.AssBestSUs.Members.CVE{s});
        D{iFile}.D_.AssBestSUs.Nonmembers.Score{s} = nanmax(D{iFile}.D_.AssBestSUs.Nonmembers.CVE{s});
        D{iFile}.D_.AssBestSUs.JointMembers.Score{s} = nanmax(D{iFile}.D_.AssBestSUs.JointMembers.CVE{s});
        
        D{iFile}.D_.AssBestSUs.All_noNoiseCorrs.Score{s} = nanmax(D{iFile}.D_.AssBestSUs.All_noNoiseCorrs.CVE{s});
        D{iFile}.D_.AssBestSUs.Members_noNoiseCorrs.Score{s} = nanmax(D{iFile}.D_.AssBestSUs.Members_noNoiseCorrs.CVE{s});
        D{iFile}.D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score{s} = nanmax(D{iFile}.D_.AssBestSUs.Nonmembers_noNoiseCorrs.CVE{s});
        try
            D{iFile}.D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s} = nanmax(D{iFile}.D_.AssBestSUs.JointMembers_noNoiseCorrs.CVE{s});
        catch
            D{iFile}.D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s} = nan(maxRange,1);
        end
    end
end
%% Summary stats - unshuffled trials
for iFile =1:length(D)
    for s=1:3
        
        % Scores
        D_.AssBestSUs.All.Score{s}(iFile,:)             = D{iFile}.D_.AssBestSUs.All.Score{s};
        D_.AssBestSUs.Members.Score{s}(iFile,:)         = D{iFile}.D_.AssBestSUs.Members.Score{s};
        D_.AssBestSUs.Nonmembers.Score{s}(iFile,:)      = D{iFile}.D_.AssBestSUs.Nonmembers.Score{s};
        D_.AssBestSUs.JointMembers.Score{s}(iFile,:)    = D{iFile}.D_.AssBestSUs.JointMembers.Score{s};
         
        % Score_shuffle
        D_.AssBestSUs.All.Score_shuffled{s}(iFile,:)             = D{iFile}.D_.AssBestSUs.All.Score_shuffled{s};
        D_.AssBestSUs.Members.Score_shuffled{s}(iFile,:)         = D{iFile}.D_.AssBestSUs.Members.Score_shuffled{s};
        D_.AssBestSUs.Nonmembers.Score_shuffled{s}(iFile,:)      = D{iFile}.D_.AssBestSUs.Nonmembers.Score_shuffled{s};
        D_.AssBestSUs.JointMembers.Score_shuffled{s}(iFile,:)    = D{iFile}.D_.AssBestSUs.JointMembers.Score_shuffled{s};
        
        % Ranks
        D_.AssBestSUs.All.ScoreRanks{s}(:,:,iFile)              = D{iFile}.D_.AssBestSUs.All.ScoreRanks{s};
        D_.AssBestSUs.Members.ScoreRanks{s}(:,:,iFile)          = D{iFile}.D_.AssBestSUs.Members.ScoreRanks{s};
        D_.AssBestSUs.Nonmembers.ScoreRanks{s}(:,:,iFile)       = D{iFile}.D_.AssBestSUs.Nonmembers.ScoreRanks{s};
        D_.AssBestSUs.JointMembers.ScoreRanks{s}(:,:,iFile)     = D{iFile}.D_.AssBestSUs.JointMembers.ScoreRanks{s};
           
    end
end

for iFile =1:length(D)
    for s=1:3
        % Best single unit assembly decoding
        D_.AssBestSUs.All.Score_mean{s} = nanmean(D_.AssBestSUs.All.Score{s});
        D_.AssBestSUs.All.Score_SEM{s} = nansem(D_.AssBestSUs.All.Score{s});
        
        D_.AssBestSUs.Members.Score_mean{s} = nanmean(D_.AssBestSUs.Members.Score{s});
        D_.AssBestSUs.Members.Score_SEM{s} = nansem(D_.AssBestSUs.Members.Score{s});
        
        D_.AssBestSUs.Nonmembers.Score_mean{s} = nanmean(D_.AssBestSUs.Nonmembers.Score{s});
        D_.AssBestSUs.Nonmembers.Score_SEM{s} = nansem(D_.AssBestSUs.Nonmembers.Score{s});
        
        D_.AssBestSUs.JointMembers.Score_mean{s} = nanmean(D_.AssBestSUs.JointMembers.Score{s});
        D_.AssBestSUs.JointMembers.Score_SEM{s} = nansem(D_.AssBestSUs.JointMembers.Score{s});
        
        D_.AssBestSUs.All.Score_shuffled_mean{s} = nanmean(D_.AssBestSUs.All.Score_shuffled{s});
        D_.AssBestSUs.All.Score_shuffled_SEM{s} = nansem(D_.AssBestSUs.All.Score_shuffled{s});
        
        D_.AssBestSUs.Members.Score_shuffled_mean{s} = nanmean(D_.AssBestSUs.Members.Score_shuffled{s});
        D_.AssBestSUs.Members.Score_shuffled_SEM{s} = nansem(D_.AssBestSUs.Members.Score_shuffled{s});
        
        D_.AssBestSUs.Nonmembers.Score_shuffled_mean{s} = nanmean(D_.AssBestSUs.Nonmembers.Score_shuffled{s});
        D_.AssBestSUs.Nonmembers.Score_shuffled_SEM{s} = nansem(D_.AssBestSUs.Nonmembers.Score_shuffled{s});
        
        D_.AssBestSUs.JointMembers.Score_shuffled_mean{s} = nanmean(D_.AssBestSUs.JointMembers.Score_shuffled{s});
        D_.AssBestSUs.JointMembers.Score_shuffled_SEM{s} = nansem(D_.AssBestSUs.JointMembers.Score_shuffled{s});
        
        D_.AssBestSUs.All.ScoreRanks_mean{s} = nanmean(D_.AssBestSUs.All.ScoreRanks{s},3);
        D_.AssBestSUs.Members.ScoreRanks_mean{s} = nanmean(D_.AssBestSUs.Members.ScoreRanks{s},3);
        D_.AssBestSUs.Nonmembers.ScoreRanks_mean{s} = nanmean(D_.AssBestSUs.Nonmembers.ScoreRanks{s},3);
        D_.AssBestSUs.JointMembers.ScoreRanks_mean{s} = nanmean(D_.AssBestSUs.JointMembers.ScoreRanks{s},3);
        
        D_.AssBestSUs.All.ScoreRanks_sem{s} = nansem(D_.AssBestSUs.All.ScoreRanks{s},3);
        D_.AssBestSUs.Members.ScoreRanks_sem{s} = nansem(D_.AssBestSUs.Members.ScoreRanks{s},3);
        D_.AssBestSUs.Nonmembers.ScoreRanks_sem{s} = nansem(D_.AssBestSUs.Nonmembers.ScoreRanks{s},3);
        D_.AssBestSUs.JointMembers.ScoreRanks_sem{s} = nansem(D_.AssBestSUs.JointMembers.ScoreRanks{s},3);
        
        
    end
end
%% Summary stats - shuffled trials
for iFile =1:length(D)
    for s=1:3
        % Scores
        D_.AssBestSUs.All_noNoiseCorrs.Score{s}(iFile,:)             = D{iFile}.D_.AssBestSUs.All_noNoiseCorrs.Score{s};
        D_.AssBestSUs.Members_noNoiseCorrs.Score{s}(iFile,:)         = D{iFile}.D_.AssBestSUs.Members_noNoiseCorrs.Score{s};
        D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score{s}(iFile,:)      = D{iFile}.D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score{s};
        try
            D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s}(iFile,:)    = D{iFile}.D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s};
        catch
            D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s}(iFile,:)    = nan(1,maxRange);
        end
    end
end

    for s=1:3
        % Best single unit assembly decoding
        D_.AssBestSUs.All_noNoiseCorrs.Score_mean{s} = nanmean(D_.AssBestSUs.All_noNoiseCorrs.Score{s});
        D_.AssBestSUs.All_noNoiseCorrs.Score_SEM{s} = nansem(D_.AssBestSUs.All_noNoiseCorrs.Score{s});
        
        D_.AssBestSUs.Members_noNoiseCorrs.Score_mean{s} = nanmean(D_.AssBestSUs.Members_noNoiseCorrs.Score{s});
        D_.AssBestSUs.Members_noNoiseCorrs.Score_SEM{s} = nansem(D_.AssBestSUs.Members_noNoiseCorrs.Score{s});
        
        D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score_mean{s} = nanmean(D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score{s});
        D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score_SEM{s} = nansem(D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score{s});
        
        D_.AssBestSUs.JointMembers_noNoiseCorrs.Score_mean{s} = nanmean(D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s});
        D_.AssBestSUs.JointMembers_noNoiseCorrs.Score_SEM{s} = nansem(D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s});
%                 
        
       
    end
%% Compare benefits of noise correlations
filt_ = 0;
for s=1:3
    % Best single unit assembly decoding
    D_.AssBestSUs.AllSubtracted.Score{s}          = smooth2a(D_.AssBestSUs.All_noNoiseCorrs.Score{s},0,filt_) - smooth2a(D_.AssBestSUs.All.Score{s},0,filt_);
    D_.AssBestSUs.MemberSubtracted.Score{s}       = smooth2a(D_.AssBestSUs.Members_noNoiseCorrs.Score{s},0,filt_) - smooth2a(D_.AssBestSUs.Members.Score{s},0,filt_);
    D_.AssBestSUs.NonmembersSubtracted.Score{s}   = smooth2a(D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score{s},0,filt_) - smooth2a(D_.AssBestSUs.Nonmembers.Score{s},0,filt_);
    D_.AssBestSUs.JointMembersSubtracted.Score{s} = smooth2a(D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s},0,filt_) - smooth2a(D_.AssBestSUs.JointMembers.Score{3},0,filt_);
    
      
    % Best single unit assembly decoding (collapse)
    D_.AssBestSUs.AllSubtracted.ScoreMean{s} = nanmean(D_.AssBestSUs.AllSubtracted.Score{s});
    D_.AssBestSUs.AllSubtracted.ScoreSEM{s}  = nansem(D_.AssBestSUs.AllSubtracted.Score{s});
    
    D_.AssBestSUs.NonmembersSubtracted.ScoreMean{s} = nanmean(D_.AssBestSUs.NonmembersSubtracted.Score{s});
    D_.AssBestSUs.NonmembersSubtracted.ScoreSEM{s}  = nansem(D_.AssBestSUs.NonmembersSubtracted.Score{s});
    
    D_.AssBestSUs.MemberSubtracted.ScoreMean{s} = nanmean(D_.AssBestSUs.MemberSubtracted.Score{s});
    D_.AssBestSUs.MemberSubtracted.ScoreSEM{s}  = nansem(D_.AssBestSUs.MemberSubtracted.Score{s});
    
    D_.AssBestSUs.JointMembersSubtracted.ScoreMean{s} = nanmean(D_.AssBestSUs.JointMembersSubtracted.Score{s});
    D_.AssBestSUs.JointMembersSubtracted.ScoreSEM{s}  = nansem(D_.AssBestSUs.JointMembersSubtracted.Score{s});
    
   
end

%% Find closest performing rank curve for each assembly
for iFile =1:length(D)
    for s=1:3
        x = D{iFile}.D_.AssReal.score{s};
        for iAss =1:size(x,1)
            nU    = x(iAss,1);
            Score = x(iAss,2);
            
            ranks_  = D{iFile}.D_.AssBestSUs.Members.ScoreRanks{s}(nU,:);
            try
                D_.AssReal.ClosestRank_BestSU{s}{iFile}(iAss,1) =  FindClosestIndex(ranks_,Score);
            catch
                D_.AssReal.ClosestRank_BestSU{s}{iFile}(iAss,1) =  NaN;
            end
            
           
        end
    end
end

%% Plot real/synthetic assemblies continuous decoding - Members
figure
subplot(1,3,1); hold on
col_ ={'b','r','g'};
tb = (1:length(tbAll)*2)*bw;
AssSize = 2:4;

for s=1:3
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssReal.CVE{s}];
    end
    
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{s})
    axis([min(tb) max(tb) 0.4 1])
end
title('Real assemblies')
legend(Areas); legend boxoff
plot([min(tb) max(tb)],[0.5 0.5],':k','HandleVisibility','off')
xlabel('Time (s)')
ylabel('Fraction correct decoding')

% Plot best ensemble assemblies
col_ ={'b','r','g'};
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,2); hold on
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestSUs.Members.CVE{s}(:,AssSize)];
    end
    
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{s})
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
    title('Best Single Unit Assemblies')
    xlabel('Time (s)')

    subplot(1,3,3); hold on
    
   

    
end
%% Plot real/synthetic assemblies continuous decoding - Non-members
figure
subplot(1,3,1); hold on
col_ ={'b','r','g'};
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssReal.CVE{s}];
    end
    
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{s})
    axis([min(tb) max(tb) 0.4 1])
end
title('Real assemblies')
legend(Areas); legend boxoff
plot([min(tb) max(tb)],[0.5 0.5],':k','HandleVisibility','off')
xlabel('Time (s)')
ylabel('Fraction correct decoding')

% Plot best ensemble assemblies
col_ ={'b','r','g'};
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,2); hold on
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestSUs.Nonmembers.CVE{s}(:,2)];
    end
    
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{s})
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
    title('Assemblies made from best single units')
    xlabel('Time (s)')

    subplot(1,3,3); hold on
    
   

end
%% Plot real/synthetic assemblies continuous decoding - Joint members
figure
subplot(1,3,1); hold on
col_ ={'b','r','g'};
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssReal.CVE{s}];
    end
    
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{s})
    axis([min(tb) max(tb) 0.4 1])
end
title('Real assemblies')
xlabel('Time (s)')
ylabel('Fraction correct decoding')
legend(Areas); legend boxoff
plot([min(tb) max(tb)],[0.5 0.5],':k','HandleVisibility','off')

% Plot best ensemble assemblies
col_ ={'b','r','g'};
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,2); hold on
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestSUs.JointMembers.CVE{s}(:,2)];
    end
    
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{s})
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
    title('Assemblies made from best single units')
    xlabel('Time (s)')

    subplot(1,3,3); hold on
    
   

end

%% *** Plot synthetic assemblies continuous decoding - overlaid
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
tb = (1:length(tbAll)*2)*bw;
AssSize = 2:4;
for s=1:3
    figure('name',Areas{s})
    
    % Plot best SU assemblies
    subplot(1,2,1); hold on
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestSUs.Nonmembers.CVE{s}(:,AssSize)];
    end
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{1},0.6)
%     plot(tb,nanmean(x,2),'color',col_{1},'LineWidth',3)

    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestSUs.Members.CVE{s}(:,AssSize)];
    end
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{2},0.6)
%     plot(tb,nanmean(x,2),'color',col_{2},'LineWidth',3)

    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestSUs.JointMembers.CVE{s}(:,AssSize)];
    end
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{3},0.6)
%     plot(tb,nanmean(x,2),'color',col_{3},'LineWidth',3)

    plot([min(tb) max(tb)],[0.5 0.5],':k','HandleVisibility','off')
    axis([min(tb) max(tb) 0.4 1])
    title('Best Single-Unit Assemblies')
    legend({'Non-members','Assembly Members','Inter-area members'}); legend boxoff
    
end
%% Plot best single unit assemblies by size - overlaid
figure('Name','Best Single units (including trial noise)'); hold on
col_ = flipud(jet(20));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for size_=1:20
        i=21-size_;
        x=[];
        for iFile = 1:length(D)
            x=[x,D{iFile}.D_.AssBestSUs.Nonmembers.CVE{s}(:,i)];
        end
        
        %         ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(size_,:))
        plot(tb,nanmean(x,2),'color',[col_(size_,:),0.6],'LineWidth',1.5)
        
    end
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
end

figure('Name','Best Single units (trial noise shuffled)'); hold on
col_ = flipud(jet(20));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for size_=1:20
        i=21-size_;
        x=[];
        for iFile = 1:length(D)
            x=[x,D{iFile}.D_.AssBestSUs.Nonmembers_noNoiseCorrs.CVE{s}(:,i)];
        end
        
        %         ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(size_,:))
        plot(tb,nanmean(x,2),'color',[col_(size_,:),0.6],'LineWidth',1.5)
        
    end
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
end

figure('Name','Best Single units (effect of removing trial correlations)'); hold on
col_ = flipud(jet(20));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for size_=1:20
        i=21-size_;
        x=[];
        for iFile = 1:length(D)
%             x=[x,D{iFile}.D_.AssBestSUs.All_noNoiseCorrs.CVE{s}(:,i) - D{iFile}.D_.AssBestSUs.All.CVE{s}(:,i)];
%             x=[x,D{iFile}.D_.AssBestSUs.Members_noNoiseCorrs.CVE{s}(:,i) - D{iFile}.D_.AssBestSUs.Members.CVE{s}(:,i)];
            try
                x=[x,D{iFile}.D_.AssBestSUs.JointMembers_noNoiseCorrs.CVE{s}(:,i) - D{iFile}.D_.AssBestSUs.JointMembers.CVE{s}(:,i)];
            catch
                x=[x,nan(402,1)];
            end
        end
        
        %         ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(size_,:))
        plot(tb,nanmean(x,2),'color',[col_(size_,:),0.6],'LineWidth',1.5)
        
    end
    plot([min(tb) max(tb)],[0 0],':k')
    axis([min(tb) max(tb) -0.2 0.2])
end

%% Plot best single unit assemblies by ranks - overlaid
figure('Name','Best Single units (including trial noise)'); hold on
col_ = flipud(jet(9));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for rank_=2:10
        
        x=[];
        for iFile = 1:length(D)
            x=[x,D{iFile}.D_.AssBestSUs.All.CVE{s}(:,rank_)];
        end
        
        ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(rank_-1,:))
        plot([min(tb) max(tb)],[0.5 0.5],':k')
        axis([min(tb) max(tb) 0.4 1])
    end
end

figure('Name','Best Single units (removing trial noise)'); hold on
col_ = flipud(jet(9));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for rank_=2:10
        
        x=[];
        for iFile = 1:length(D)
            x=[x,D{iFile}.D_.AssBestSUs.All_noNoiseCorrs.CVE{s}(:,rank_)];
        end
        
        ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(rank_-1,:))
        plot([min(tb) max(tb)],[0.5 0.5],':k')
        axis([min(tb) max(tb) 0.4 1])
    end
end
figure('Name','Best Single units (removing trial noise)'); hold on
col_ = flipud(jet(9));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for rank_=2:10
        
        x=[];
        for iFile = 1:length(D)
            x=[x,D{iFile}.D_.AssBestSUs.All_noNoiseCorrs.CVE{s}(:,rank_)-D{iFile}.D_.AssBestSUs.All.CVE{s}(:,rank_)];
        end
        
        ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(rank_-1,:))
        plot([min(tb) max(tb)],[0 0],':k')
        axis([min(tb) max(tb) -0.4 0.4])
    end
end
%% Plot closest rank histograms
figure; hold on
bins = 1:11;
col_ ={'b','r','g'};
for s=1:3
    subplot(1,2,1); hold on
    x = cell2mat( D_.AssReal.ClosestRank_BestSU{s}');
    x = histc(x,bins); x=x./sum(x);
    plot([0 10],[s s]-1,':k','HandleVisibility','off','LineWidth',1.5)
    stairs(bins-0.5,s-1+x,'Color',col_{s},'LineWidth',1.5)
    title('Best Single Unit aggregation')
    ylabel('Fraction of Assemblies')
    xlabel('Closest ranked synthetic assembly')
    axis([0.5 10.5 0 3])
    set(gca,'YTick',[])
    
   
end
legend(Areas,'Orientation','horizontal');legend boxoff
%% *** Plot real and shuffled assemblies
x=[]; x_=[];

figure
for s=1:3
    for iFile =1:length(D)
        subplot(1,3,s); hold on
        x = [x; D{iFile}.D_.AssReal.score{s}];
        x_= [x_;D{iFile}.D_.AssReal.score_shuffled{s}];
    end
    
    % plot real assemblies
    ranks = 1:10;
    for rank_ = ranks
        x_mean(rank_)=nanmean(x(x(:,1)==rank_,2));
        x_sem(rank_) =nansem(x(x(:,1)==rank_,2));
    end
    x_sem(isnan(x_mean))=[];
    ranks(isnan(x_mean))=[];
    x_mean(isnan(x_mean))=[];

    scatter(x(:,1),x(:,2),10,[0.6 0.6 0.6],'o')
    scatter(x_(:,1),x_(:,2),'.k')
    errorbar(ranks,x_mean,x_sem,'k','LineWidth',1.5)
    scatter(ranks,x_mean,'ok','filled')
    % plot shuffled assemblies
    ranks = 1:10;
    for rank_ = ranks
        x_mean(rank_)=nanmean(x_(x_(:,1)==rank_,2));
        x_sem(rank_) =nansem(x_(x_(:,1)==rank_,2));
    end
    x_sem(isnan(x_mean))=[];
    ranks(isnan(x_mean))=[];
    x_mean(isnan(x_mean))=[];
    errorbar(ranks,x_mean,x_sem,':k','LineWidth',1.5)
    scatter(ranks,x_mean,'.k','filled')
    %     scatter(x_(:,1),x_(:,2),'.k')
    
    if s==1
        ylabel({'Peak decoding score';'(% Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0.5 0.5],':k')

    if s==3
        axis([2 16 0.4 1])     
    else
        axis([2 6 0.4 1])     
    end
end
clear x x_
text(16,0.53,'Chance','HorizontalAlignment','right')
legend({'Assemblies','Scrambled trial labels'},'Location','southeast'); legend boxoff
%% Plot best SU decoding from all units
x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
  
    m = D_.AssBestSUs.All.Score_mean{s}(2:end);
    e = D_.AssBestSUs.All.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
    
    m = D_.AssBestSUs.All_noNoiseCorrs.Score_mean{s}(2:end);
    e = D_.AssBestSUs.All_noNoiseCorrs.Score_SEM{s}(2:end);
    plot(x,m+e,'color',[0 0 0],'LineWidth',1.5)
    plot(x,m-e,'color',[0 0 0],'LineWidth',1.5)
    
    if s==1
        ylabel({'Peak decoding score';'(% Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0.5 0.5],':k')
    axis([2 20 0 1])     
end
legend({'Best Single Unit Aggregation'},'Location','southeast');legend boxoff
clear x x_
%% Plot best SU decoding from non-member units
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
  
    m = D_.AssBestSUs.Nonmembers.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Nonmembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
  
    m = D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score_SEM{s}(2:end);
    plot(x,m+e,'color',[0 0 0],'LineWidth',1.5)
    plot(x,m-e,'color',[0 0 0],'LineWidth',1.5)
    
    
    if s==1
        ylabel({'Peak decoding score';'(% Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0.5 0.5],':k')
    axis([2 20 0 1])     
end
legend({'Best Single Unit Aggregation'},'Location','southeast');legend boxoff
clear x x_
%% Plot best SU decoding from member units
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
       
    m = D_.AssBestSUs.Members.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Members.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
        
    m = D_.AssBestSUs.Members_noNoiseCorrs.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Members_noNoiseCorrs.Score_SEM{s}(2:end);
    plot(x,m+e,'color',[0 0 0],'LineWidth',1.5)
    plot(x,m-e,'color',[0 0 0],'LineWidth',1.5)
    
    if s==1
        ylabel({'Peak decoding score';'(% Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0.5 0.5],':k')
    axis([2 20 0 1])     
end
legend({'Best Single Unit Aggregation'},'Location','southeast');legend boxoff
clear x x_
%% Plot best SU decoding from joint HP-PFC units

x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    
    m = D_.AssBestSUs.JointMembers.Score_mean{s}(2:end);
    e = D_.AssBestSUs.JointMembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
       
    m = D_.AssBestSUs.JointMembers_noNoiseCorrs.Score_mean{s}(2:end);
    e = D_.AssBestSUs.JointMembers_noNoiseCorrs.Score_SEM{s}(2:end);
    plot(x,m+e,'color',[0 0 0],'LineWidth',1.5)
    plot(x,m-e,'color',[0 0 0],'LineWidth',1.5)
    
    
    if s==1
        ylabel({'Peak decoding score';'(% Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0.5 0.5],':k')
    axis([2 20 0 1])     
end
legend({'Location','southeast'});legend boxoff
clear x x_
%% Plot best SU decoding from all classes of units
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    m = D_.AssBestSUs.Nonmembers.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Nonmembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = D_.AssBestSUs.Members.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Members.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.9 0.6 0],0.8)
    
    m = D_.AssBestSUs.JointMembers.Score_mean{s}(2:end);
    e = D_.AssBestSUs.JointMembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0 0.6],0.8)
    
    if s==1
        ylabel({'Peak decoding score';'(% Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0.5 0.5],':k')
    axis([2 20 0 1])     
end
clear x x_
legend({'Non-member units','Member units','Joint area members'},'Location','southeast');legend boxoff

%% Plot real and shuffled assemblies with best synthetics

col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

x=[]; x_=[];
nU = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    % Plot synthetic assemblies
    m = D_.AssBestSUs.Members.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Members.Score_SEM{s}(2:end);
    %     ciplot(m+e,m-e,nU,'b',0.8)
    plot(nU,m,'color',[0.9 0.6 0 0.8])
    
    m = D_.AssBestSUs.Nonmembers.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Nonmembers.Score_SEM{s}(2:end);
    %     ciplot(m+e,m-e,nU,'r',0.8)
    plot(nU,m,'color',[0.6 0.6 0.6 0.8])
    
   
    
    m = D_.AssBestSUs.Members.Score_shuffled_mean{s}(2:end);
    e = D_.AssBestSUs.Members.Score_shuffled_SEM{s}(2:end);
    plot(nU,m+e,'color',[0.6 0.6 0.6 0.8],'LineStyle',':')
    plot(nU,m-e,'color',[0.6 0.6 0.6 0.8],'LineStyle',':')
    m = D_.AssBestSUs.Nonmembers.Score_shuffled_mean{s}(2:end);
    e = D_.AssBestSUs.Nonmembers.Score_shuffled_SEM{s}(2:end);
    plot(nU,m+e,'color',[0.6 0.6 0.6],'LineStyle',':')
    plot(nU,m-e,'color',[0.6 0.6 0.6],'LineStyle',':')
    
    for iFile =1:length(D)
        x = [x; D{iFile}.D_.AssReal.score{s}];
        x_= [x_;D{iFile}.D_.AssReal.score_shuffled{s}];
    end
    
    % plot real assemblies
    ranks = 1:10;
    for rank_ = ranks
        x_mean(rank_)=nanmean(x(x(:,1)==rank_,2));
        x_sem(rank_) =nansem(x(x(:,1)==rank_,2));
    end
    x_sem(isnan(x_mean))=[];
    ranks(isnan(x_mean))=[];
    x_mean(isnan(x_mean))=[];
    errorbar(ranks,x_mean,x_sem,'g','LineWidth',1.5)
    scatter(ranks,x_mean,'og','filled')
    scatter(x(:,1),x(:,2),10,[0 1 0],'o')
    scatter(x_(:,1),x_(:,2),'.g')
    
    % plot shuffled assemblies
    ranks = 1:10;
    for rank_ = ranks
        x_mean(rank_)=nanmean(x_(x_(:,1)==rank_,2));
        x_sem(rank_) =nansem(x_(x_(:,1)==rank_,2));
    end
    x_sem(isnan(x_mean))=[];
    ranks(isnan(x_mean))=[];
    x_mean(isnan(x_mean))=[];
    errorbar(ranks,x_mean,x_sem,':g','LineWidth',1.5)
    scatter(ranks,x_mean,'og','filled')
    %     scatter(x_(:,1),x_(:,2),'.k')
%     plot([1 20],[0.5 0.5],':k','HandleVisibility','off')
    
    if s==3
        axis([1 15 0.4 1])
        text(15,0.52,'Chance','HorizontalAlignment','right')
    else
        axis([1 10 0.4 1])
        text(10,0.52,'Chance','HorizontalAlignment','right')
    end
    title([Areas{s} ' Assemblies'])
end
clear x x_
%% Plot best single unit synthetic assemblies
x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    %     m = D_.AssBestSUs.All.Score_mean{s}(2:end);
    %     e = D_.AssBestSUs.All.Score_SEM{s}(2:end);
    %     ciplot(m+e,m-e,x,'k',0.8)
    
    m = D_.AssBestSUs.Members.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Members.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'b',0.8)
    
    m = D_.AssBestSUs.Nonmembers.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Nonmembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'r',0.8)
    
        m = D_.AssBestSUs.JointMembers.Score_mean{s}(2:end);
        e = D_.AssBestSUs.JointMembers.Score_SEM{s}(2:end);
        ciplot(m+e,m-e,x,'g',0.8)
end
for iFile =1:length(D)
    for s=1:3
        subplot(1,3,s); hold on
        x = D{iFile}.D_.AssReal.score{s};
        x_= D{iFile}.D_.AssReal.score_shuffled{s};
        scatter(x(:,1),x(:,2),'og','filled')
        scatter(x_(:,1),x_(:,2),'og')
    end
end
clear x x_

%% Plot rank surfaces - best SU
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

figure
for s=1:3
    subplot(1,3,s); hold on
    
    %         mesh(smooth2a(D_.AssBestAssem.All.ScoreRanks_mean{s},2,2),'EdgeColor','k')
    x = smooth2a(D_.AssBestSUs.Nonmembers.ScoreRanks_mean{s},1,1);
    mesh(x,'EdgeColor',col_{1},'FaceColor',col_{1},'FaceAlpha',0.5)
    plot3(ones(1,20),1:20,x(:,1)','Color',col_{1},'LineWidth',1.5,'HandleVisibility','off')
    x = smooth2a(D_.AssBestSUs.Members.ScoreRanks_mean{s},1,1);
    mesh(x,'EdgeColor',col_{2},'FaceColor',col_{2},'FaceAlpha',0.5)
    plot3(ones(1,20),1:20,x(:,1)','Color',col_{2},'LineWidth',1.5,'HandleVisibility','off')
    x = smooth2a(D_.AssBestSUs.JointMembers.ScoreRanks_mean{s},1,1);
    mesh(x,'EdgeColor',col_{3},'FaceColor',col_{3},'FaceAlpha',0.5)
    plot3(ones(1,20),1:20,x(:,1)','Color',col_{3},'LineWidth',1.5,'HandleVisibility','off')
    
    set(gca,'XTick',[1,2:2:10],'XTickLabel',{'Best','2^n^d','4^t^h','6^t^h','8^t^h','10^t^h'},'XTickLabelRotation',-60)

    view(3)
    axis([ 1 4 2 20 0.9 1])
    if s==1
        ylabel({'Assembly size';'No. neurons'})
        zlabel({'Peak decoding';'(% Correct classification)'})
    elseif s==3
        xlabel('Assembly rank')
    end
    title({Areas{s} 'Assemblies'})
    %     box on
end
legend({'Non-member units','Member units','Joint area members'},'Location','southeast');legend boxoff
%% *** Plot rank surfaces - best SU with assemblies
col_ = flipud(jet(10));
x=[];x_=[];
figure
for s=1:3
    subplot(1,3,s); hold on
    
%         m = D_.AssBestSUs.Members.ScoreRanks_mean{s};
%         e = D_.AssBestSUs.Members.ScoreRanks_sem{s};
%         m(:,1) = D_.AssBestSUs.Members.Score_mean{s};
%         e(:,1) = D_.AssBestSUs.Members.Score_SEM{s};

        m = D_.AssBestSUs.All.ScoreRanks_mean{s};
        e = D_.AssBestSUs.All.ScoreRanks_sem{s};
        m(:,1) = D_.AssBestSUs.All.Score_mean{s};
        e(:,1) = D_.AssBestSUs.All.Score_SEM{s};
        
        m = smooth2a(m,1,0); e = smooth2a(e,1,0);
        
    for i=1:10
        if ismember(i,[1,10])
            plot(2:20,m(2:end,i),'color',col_(i,:),'LineWidth',1.5)
        else
            plot(2:20,m(2:end,i),'color',col_(i,:),'LineWidth',1.5,'HandleVisibility','off')
        end
%         errorbar(2:20,m(2:end,i),e(2:end,i),'color',col_(i-1,:),'LineWidth',1.5)
    end
    
    
    for iFile =1:length(D)
        x = [x; D{iFile}.D_.AssReal.score{s}];
        x_= [x_;D{iFile}.D_.AssReal.score_shuffled{s}];
    end
    % plot real assemblies
    ranks = 1:10;
    for rank_ = ranks
        x_mean(rank_)=nanmean(x(x(:,1)==rank_,2));
        x_sem(rank_) =nansem(x(x(:,1)==rank_,2));
    end
    x_sem(isnan(x_mean))=[];
    ranks(isnan(x_mean))=[];
    x_mean(isnan(x_mean))=[];
    errorbar(ranks,x_mean,x_sem,'k','LineWidth',1.5)
    scatter(ranks,x_mean,'ok','filled','HandleVisibility','off')
    scatter(x(:,1),x(:,2),10,[0.6 0.6 0.6],'o','HandleVisibility','off')
    scatter(x_(:,1),x_(:,2),5,'.k','HandleVisibility','off')
    
    % plot shuffled assemblies
    ranks = 1:10;
    for rank_ = ranks
        x_mean(rank_)= nanmean(x_(x_(:,1)==rank_,2));
        x_sem(rank_) = nansem(x_(x_(:,1)==rank_,2));
    end
    x_sem(isnan(x_mean))=[];
    ranks(isnan(x_mean))=[];
    x_mean(isnan(x_mean))=[];
    errorbar(ranks,x_mean,x_sem,':k','LineWidth',1)
    scatter(ranks,x_mean,5,'ok','filled','HandleVisibility','off')
    %     scatter(x_(:,1),x_(:,2),'.k')
    plot([1 20],[0.5 0.5],':k','HandleVisibility','off')
    axis([0 10 0.4 1])
    
        title({Areas{s}, 'Assemblies'})

      if s==1
        ylabel({'Peak decoding';'(% Correct classification)'})
      elseif s==2
        xlabel('Assembly size (No. Units)')
      end

end
legend({'Best Synthetic Assembly','10^t^h best','Real Assemblies','Shuffled real assemblies'})
    text(10,0.52,'Chance','HorizontalAlignment','right')

legend boxoff
clear x x_

%% *** Plot average benefit of noise correlations
x = 2:20;
ymax = 15;
figure('Name','All units')
for s=1:3
    subplot(1,3,s); hold on
    
   
    
    m = 100*D_.AssBestSUs.AllSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestSUs.AllSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
  
    
    if s==1
        ylabel({'Benefit over shuffled noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -ymax ymax])     
end
legend({'Best Single Unit Aggregation'},'Location','southeast');legend boxoff

figure('Name','Nonmember units')
for s=1:3
    subplot(1,3,s); hold on
    
    m = 100*D_.AssBestSUs.NonmembersSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestSUs.NonmembersSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
  
    
    if s==1
        ylabel({'Benefit over shuffled noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -ymax ymax])     
end
legend({'Best Single Unit Aggregation'},'Location','southeast');legend boxoff

figure('Name','Member units')
for s=1:3
    subplot(1,3,s); hold on
    m = 100*D_.AssBestSUs.MemberSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestSUs.MemberSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
  
    
    if s==1
        ylabel({'Benefit over shuffled noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -ymax ymax])     
end
legend({'Best Single Unit Aggregation'},'Location','southeast');legend boxoff

figure('Name','Joint-area units')
for s=1:3
    subplot(1,3,s); hold on
    
        m = 100*D_.AssBestSUs.JointMembersSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestSUs.JointMembersSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
  
    
    if s==1
        ylabel({'Benefit over shuffled noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -ymax ymax])     
end
legend({'Best Single Unit Aggregation'},'Location','southeast');legend boxoff


clear x x_


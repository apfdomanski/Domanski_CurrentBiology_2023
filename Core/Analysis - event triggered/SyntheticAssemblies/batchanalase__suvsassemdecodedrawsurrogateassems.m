% %%%%%% PREAMBLE %%%%%%
clear                              

Delays_ = {'Short','Medium','Long'};
Target = 'LONG';
Areas = {'PFC','HP','Joint'};

bw=0.05; tlimsAll = [-5 5]; tbAll = tlimsAll(1):bw:tlimsAll(2);

maxRange = 20;    % Upper limit of assembly size to explore
noTrials = 5;   % How many trials to consider simultanrously
nReps = 50;     % how many repeated draws among trials to perform


InfoCriteria = 'max'; % 'mean','max'

if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
    %pat = '/Volumes/HDD2/DNMTP/raw/'
    pat = '/Volumes/Aleks Data Portable/AssemblyAnalysis/raw/';
elseif isunix
    pat ='/Volumes/Data/DNMTP/raw';
end

fileList = dir([pat 'SyntheticAssemblies' filesep '*' Target '*AssemblyComparison_redux8.mat']);
% fileList = dir([pat 'SyntheticAssemblies' filesep '*' Target '*AssemblyComparison_redux7PosthocShuffle.mat']);
%% Batch import
for iFile =1:length(fileList)
    
    %fnIn = fullfile(pat, 'SyntheticAssemblies', fileList(iFile).name);
    fnIn = fullfile(pat, 'SyntheticAssemblies',fileList(iFile).name); 
    D{iFile} = load(fnIn,'D_');

end
clear temp 
%% Batch import (separate noise runs)
mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
for iFile =1:length(fileList)
    
    fnIn = fullfile(pat, 'SyntheticAssemblies', fileList(iFile).name);
    D{iFile} = load(fnIn,'D_');
    fnIn = fullfile(pat, 'SyntheticAssemblies', [strtok(fileList(iFile).name,'.'),'_NoNoiseCorrs.mat']);
    temp = load(fnIn,'D_');
    
    D{iFile}.D_.AssBestSUs   = mergestructs(D{iFile}.D_.AssBestSUs,temp.D_.AssBestSUs);
    D{iFile}.D_.AssBestAssem = mergestructs(D{iFile}.D_.AssBestAssem,temp.D_.AssBestAssem);
end
clear temp 
%% Batch import (separate SU and ensemble runs)
for iFile =1:length(fileList)
    
    %fnIn = fullfile(pat, 'SyntheticAssemblies', fileList(iFile).name);
    fnIn = fullfile(pat, 'SyntheticAssemblies', [strtok(fileList(iFile).name,'.'),'_redux4.mat']); % Ass, SU, best assem aggrn
    D{iFile} = load(fnIn,'D_');
    fnIn = fullfile(pat, 'SyntheticAssemblies', [strtok(fileList(iFile).name,'.'),'_redux3.mat']); % best SU decoding aggrn
    temp = load(fnIn,'D_');
    D{iFile}.D_.AssBestSUs =  temp.D_.AssBestSUs;
end
clear temp 
%% Batch import (separate shuffling on ensemble runs)
for iFile =1:length(fileList)
    
    %fnIn = fullfile(pat, 'SyntheticAssemblies', fileList(iFile).name);
    fnIn = fullfile(pat, 'SyntheticAssemblies', [strtok(fileList(iFile).name,'8'),'8.mat']); % Ass, SU, best assem aggrn
    D{iFile} = load(fnIn,'D_');
    fnIn = fullfile(pat, 'SyntheticAssemblies', fileList(iFile).name); % best SU decoding aggrn
    temp = load(fnIn,'D_');
    D{iFile}.D_.AssBestAssem.All_noNoiseCorrs       =  temp.D_.AssBestAssem.All_noNoiseCorrs;
    D{iFile}.D_.AssBestAssem.Members_noNoiseCorrs =  temp.D_.AssBestAssem.Members_noNoiseCorrs;
    D{iFile}.D_.AssBestAssem.Nonmembers_noNoiseCorrs =  temp.D_.AssBestAssem.Nonmembers_noNoiseCorrs;
    D{iFile}.D_.AssBestAssem.JointMembers_noNoiseCorrs =  temp.D_.AssBestAssem.JointMembers_noNoiseCorrs;
end
clear temp 
%% Summary stats - unshuffled trials
for iFile =1:length(D)
    for s=1:3
        % Scores
        D_.AssBestSUs.All.Score{s}(iFile,:)             = D{iFile}.D_.AssBestSUs.All.Score{s};
        D_.AssBestSUs.Members.Score{s}(iFile,:)         = D{iFile}.D_.AssBestSUs.Members.Score{s};
        D_.AssBestSUs.Nonmembers.Score{s}(iFile,:)      = D{iFile}.D_.AssBestSUs.Nonmembers.Score{s};
        D_.AssBestSUs.JointMembers.Score{s}(iFile,:)    = D{iFile}.D_.AssBestSUs.JointMembers.Score{s};
        
        D_.AssBestAssem.All.Score{s}(iFile,:)             = D{iFile}.D_.AssBestAssem.All.Score{s};
        D_.AssBestAssem.Members.Score{s}(iFile,:)         = D{iFile}.D_.AssBestAssem.Members.Score{s};
        D_.AssBestAssem.Nonmembers.Score{s}(iFile,:)      = D{iFile}.D_.AssBestAssem.Nonmembers.Score{s};
        D_.AssBestAssem.JointMembers.Score{s}(iFile,:)    = D{iFile}.D_.AssBestAssem.JointMembers.Score{s};
        
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
        
        D_.AssBestAssem.All.ScoreRanks{s}(:,:,iFile)        = D{iFile}.D_.AssBestAssem.All.ScoreRanks{s};
        D_.AssBestAssem.Members.ScoreRanks{s}(:,:,iFile)        = D{iFile}.D_.AssBestAssem.Members.ScoreRanks{s};
        D_.AssBestAssem.Nonmembers.ScoreRanks{s}(:,:,iFile)        = D{iFile}.D_.AssBestAssem.Nonmembers.ScoreRanks{s};
        D_.AssBestAssem.JointMembers.ScoreRanks{s}(:,:,iFile)        = D{iFile}.D_.AssBestAssem.JointMembers.ScoreRanks{s};
        
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
        
        % Best ensemble Assembly decoding
        D_.AssBestAssem.All.Score_mean{s} = nanmean(D_.AssBestAssem.All.Score{s});
        D_.AssBestAssem.All.Score_SEM{s} = nansem(D_.AssBestAssem.All.Score{s});
        
        D_.AssBestAssem.Members.Score_mean{s} = nanmean(D_.AssBestAssem.Members.Score{s});
        D_.AssBestAssem.Members.Score_SEM{s} = nansem(D_.AssBestAssem.Members.Score{s});
        
        D_.AssBestAssem.Nonmembers.Score_mean{s} = nanmean(D_.AssBestAssem.Nonmembers.Score{s});
        D_.AssBestAssem.Nonmembers.Score_SEM{s} = nansem(D_.AssBestAssem.Nonmembers.Score{s});
        
        D_.AssBestAssem.JointMembers.Score_mean{s} = nanmean(D_.AssBestAssem.JointMembers.Score{s});
        D_.AssBestAssem.JointMembers.Score_SEM{s} = nansem(D_.AssBestAssem.JointMembers.Score{s});
        
        D_.AssBestAssem.All.ScoreRanks_mean{s} = nanmean(D_.AssBestAssem.All.ScoreRanks{s},3);
        D_.AssBestAssem.Members.ScoreRanks_mean{s} = nanmean(D_.AssBestAssem.Members.ScoreRanks{s},3);
        D_.AssBestAssem.Nonmembers.ScoreRanks_mean{s} = nanmean(D_.AssBestAssem.Nonmembers.ScoreRanks{s},3);
        D_.AssBestAssem.JointMembers.ScoreRanks_mean{s} = nanmean(D_.AssBestAssem.JointMembers.ScoreRanks{s},3);
        
        D_.AssBestAssem.All.ScoreRanks_sem{s} = nansem(D_.AssBestAssem.All.ScoreRanks{s},3);
        D_.AssBestAssem.Members.ScoreRanks_sem{s} = nansem(D_.AssBestAssem.Members.ScoreRanks{s},3);
        D_.AssBestAssem.Nonmembers.ScoreRanks_sem{s} = nansem(D_.AssBestAssem.Nonmembers.ScoreRanks{s},3);
        D_.AssBestAssem.JointMembers.ScoreRanks_sem{s} = nansem(D_.AssBestAssem.JointMembers.ScoreRanks{s},3);
        
    end
end
%% Summary stats - shuffled trials
for iFile =1:length(D)
    for s=1:3
        % Scores
        D_.AssBestSUs.All_noNoiseCorrs.Score{s}(iFile,:)             = D{iFile}.D_.AssBestSUs.All_noNoiseCorrs.Score{s};
        D_.AssBestSUs.Members_noNoiseCorrs.Score{s}(iFile,:)         = D{iFile}.D_.AssBestSUs.Members_noNoiseCorrs.Score{s};
        D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score{s}(iFile,:)      = D{iFile}.D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score{s};
        D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s}(iFile,:)    = D{iFile}.D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s};
        
        try
            D_.AssBestAssem.All_noNoiseCorrs.Score{s}(iFile,:)             = D{iFile}.D_.AssBestAssem.All_noNoiseCorrs.Score{s};
        catch
             D_.AssBestAssem.All_noNoiseCorrs.Score{s}(iFile,:)             = nan(20,1);
        end
        try
            D_.AssBestAssem.Members_noNoiseCorrs.Score{s}(iFile,:)         = D{iFile}.D_.AssBestAssem.Members_noNoiseCorrs.Score{s};
        catch
            D_.AssBestAssem.Members_noNoiseCorrs.Score{s}(iFile,:)         = nan(20,1);
        end
        try
            D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score{s}(iFile,:)      = D{iFile}.D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score{s};
        catch
            D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score{s}(iFile,:)      = nan(20,1);
        end
        % Ranks
        D_.AssBestSUs.All_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)              = D{iFile}.D_.AssBestSUs.All_noNoiseCorrs.ScoreRanks{s};
        D_.AssBestSUs.Members_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)          = D{iFile}.D_.AssBestSUs.Members_noNoiseCorrs.ScoreRanks{s};
        D_.AssBestSUs.Nonmembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)       = D{iFile}.D_.AssBestSUs.Nonmembers_noNoiseCorrs.ScoreRanks{s};
        D_.AssBestSUs.JointMembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)     = D{iFile}.D_.AssBestSUs.JointMembers_noNoiseCorrs.ScoreRanks{s};
        
        try
            D_.AssBestAssem.All_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)              = D{iFile}.D_.AssBestAssem.All_noNoiseCorrs.ScoreRanks{s};
        catch
            D_.AssBestAssem.All_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)              = nan(20,10); 
        end
        try
            D_.AssBestAssem.Members_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)          = D{iFile}.D_.AssBestAssem.Members_noNoiseCorrs.ScoreRanks{s};
        catch
            D_.AssBestAssem.Members_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)          = nan(20,10);
        end
        try
            D_.AssBestAssem.Nonmembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)       = D{iFile}.D_.AssBestAssem.Nonmembers_noNoiseCorrs.ScoreRanks{s};
        catch
            D_.AssBestAssem.Nonmembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)       = nan(20,10);
        end
        try
            D_.AssBestAssem.JointMembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)     = D{iFile}.D_.AssBestAssem.JointMembers_noNoiseCorrs.ScoreRanks{s};
        catch
            D_.AssBestAssem.JointMembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile)     = nan(20,10);
        end
        
        try
            D_.AssBestAssem.JointMembers_noNoiseCorrs.Score{s}(iFile,:)    = D{iFile}.D_.AssBestAssem.JointMembers_noNoiseCorrs.Score{s};
        catch
            D_.AssBestAssem.JointMembers_noNoiseCorrs.Score{s}(iFile,:)    = nan(1,maxRange);

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
                
        
        % Best ensemble Assembly decoding
        D_.AssBestAssem.All_noNoiseCorrs.Score_mean{s} = nanmean(D_.AssBestAssem.All_noNoiseCorrs.Score{s});
        D_.AssBestAssem.All_noNoiseCorrs.Score_SEM{s} = nansem(D_.AssBestAssem.All_noNoiseCorrs.Score{s});
        
        D_.AssBestAssem.Members_noNoiseCorrs.Score_mean{s} = nanmean(D_.AssBestAssem.Members_noNoiseCorrs.Score{s});
        D_.AssBestAssem.Members_noNoiseCorrs.Score_SEM{s} = nansem(D_.AssBestAssem.Members_noNoiseCorrs.Score{s});
        
        D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score_mean{s} = nanmean(D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score{s});
        D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score_SEM{s} = nansem(D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score{s});
        
        D_.AssBestAssem.JointMembers_noNoiseCorrs.Score_mean{s} = nanmean(D_.AssBestAssem.JointMembers_noNoiseCorrs.Score{s});
        D_.AssBestAssem.JointMembers_noNoiseCorrs.Score_SEM{s} = nansem(D_.AssBestAssem.JointMembers_noNoiseCorrs.Score{s});
        
    end
%% Calculate performance
for iFile =1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
    %% Batch process behavioural performance
    L = [length(t.Short.ChoicePress_LeftCorrect)./(length(t.Short.ChoicePress_LeftCorrect)+length(t.Short.ChoicePress_LeftError)) , ...
        length(t.Medium.ChoicePress_LeftCorrect)./(length(t.Medium.ChoicePress_LeftCorrect)+length(t.Medium.ChoicePress_LeftError)) , ...
        length(t.Long.ChoicePress_LeftCorrect)./(length(t.Long.ChoicePress_LeftCorrect)+length(t.Long.ChoicePress_LeftError))];
    
    
    R = [length(t.Short.ChoicePress_RightCorrect)./(length(t.Short.ChoicePress_RightCorrect)+length(t.Short.ChoicePress_RightError)) , ...
        length(t.Medium.ChoicePress_RightCorrect)./(length(t.Medium.ChoicePress_RightCorrect)+length(t.Medium.ChoicePress_RightError)) , ...
        length(t.Long.ChoicePress_RightCorrect)./(length(t.Long.ChoicePress_RightCorrect)+length(t.Long.ChoicePress_RightError))];
    
    D_.Accuracy = (L+R)./2;
    
    C_ =  [length(t.Short.ChoicePress_LeftCorrect)  + length(t.Short.ChoicePress_RightCorrect),...
        length(t.Medium.ChoicePress_LeftCorrect) + length(t.Medium.ChoicePress_RightCorrect),...
        length(t.Long.ChoicePress_LeftCorrect)   + length(t.Long.ChoicePress_RightCorrect)];
    
    E_ =  [length(t.Short.ChoicePress_LeftError)    + length(t.Short.ChoicePress_RightError),...
        length(t.Medium.ChoicePress_LeftError)   + length(t.Medium.ChoicePress_RightError),...
        length(t.Long.ChoicePress_LeftError)     + length(t.Long.ChoicePress_RightError)];
    
    D_.pCorr = C_./(C_+E_);
    D_.AboveChance = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;
    clear L R C_ E_
    
end
clear PFCcells PFCinter PFCtonic trangeleft_choice trangeleft_sample trangeright_sample trangeright_choice HPtonic HPinter HPcells ERRORtrangeright_sample ERRORtrangeright_choice ERRORtrangeleft_sample ERRORtrangeleft_choice
    %% Compare benefits of noise correlations
    smooth_ = 0;
    smooth_2 = 0;
    smooth_3 = 1;
    maxRank = 5;
    
            % Best single unit assembly decoding (collapse)
            
            for s=1:3
                for iFile =1:length(D)
                    temp = smooth2a(D_.AssBestSUs.All.ScoreRanks{s}(:,:,iFile),0,smooth_) - ...
                           smooth2a(D_.AssBestSUs.All_noNoiseCorrs.ScoreRanks{s}(:,:,iFile),0,smooth_);
                    D_.AssBestSUs.AllSubtracted.Score{s}(iFile,:) = nanmean(smooth2a(temp(:,1:maxRank),smooth_2,0),2);
                    
                    temp = smooth2a(D_.AssBestSUs.Members.ScoreRanks{s}(:,:,iFile),0,smooth_) - ...
                           smooth2a(D_.AssBestSUs.Members_noNoiseCorrs.ScoreRanks{s}(:,:,iFile),0,smooth_);
                    D_.AssBestSUs.MemberSubtracted.Score{s}(iFile,:) =  nanmean(smooth2a(temp(:,1:maxRank),smooth_2,0),2);
                    
                    temp = smooth2a(D_.AssBestSUs.Nonmembers.ScoreRanks{s}(:,:,iFile),0,smooth_) - ...
                           smooth2a(D_.AssBestSUs.Nonmembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile),0,smooth_);
                    D_.AssBestSUs.NonmembersSubtracted.Score{s}(iFile,:) =  nanmean(smooth2a(temp(:,1:maxRank),smooth_2,0),2);
                    
                    temp = smooth2a(D_.AssBestSUs.JointMembers.ScoreRanks{s}(:,:,iFile),0,smooth_) - ...
                           smooth2a(D_.AssBestSUs.JointMembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile),0,smooth_);
                    D_.AssBestSUs.JointMembersSubtracted.Score{s}(iFile,:) =  nanmean(smooth2a(temp(:,1:maxRank),smooth_2,0),2);
                end
                D_.AssBestSUs.AllSubtracted.ScoreMean{s} = nanmean(smooth2a(D_.AssBestSUs.AllSubtracted.Score{s},0,smooth_3));
                D_.AssBestSUs.AllSubtracted.ScoreSEM{s}  = nansem(smooth2a(D_.AssBestSUs.AllSubtracted.Score{s},0,smooth_3));
                
                D_.AssBestSUs.NonmembersSubtracted.ScoreMean{s} = nanmean(smooth2a(D_.AssBestSUs.NonmembersSubtracted.Score{s},0,smooth_3));
                D_.AssBestSUs.NonmembersSubtracted.ScoreSEM{s}  = nansem(smooth2a(D_.AssBestSUs.NonmembersSubtracted.Score{s},0,smooth_3));
                
                D_.AssBestSUs.MemberSubtracted.ScoreMean{s} = nanmean(smooth2a(D_.AssBestSUs.MemberSubtracted.Score{s},0,smooth_3));
                D_.AssBestSUs.MemberSubtracted.ScoreSEM{s}  = nansem(smooth2a(D_.AssBestSUs.MemberSubtracted.Score{s},0,smooth_3));
                
                D_.AssBestSUs.JointMembersSubtracted.ScoreMean{s} = nanmean(smooth2a(D_.AssBestSUs.JointMembersSubtracted.Score{s},0,smooth_3));
                D_.AssBestSUs.JointMembersSubtracted.ScoreSEM{s}  = nansem(smooth2a(D_.AssBestSUs.JointMembersSubtracted.Score{s},0,smooth_3));
            end
            
            % Best ensemble assembly decoding (collapse)

    for s=1:3
       for iFile =1:length(D)
                    temp = smooth2a(D_.AssBestAssem.All.ScoreRanks{s}(:,:,iFile),0,smooth_) - ...
                           smooth2a(D_.AssBestAssem.All_noNoiseCorrs.ScoreRanks{s}(:,:,iFile),0,smooth_);
                    D_.AssBestAssem.AllSubtracted.Score{s}(iFile,:) =  nanmean(smooth2a(temp(:,1:maxRank),smooth_2,0),2);
                    
                    temp = smooth2a(D_.AssBestAssem.Members.ScoreRanks{s}(:,:,iFile),0,smooth_) - ...
                           smooth2a(D_.AssBestAssem.Members_noNoiseCorrs.ScoreRanks{s}(:,:,iFile),0,smooth_);
                    D_.AssBestAssem.MemberSubtracted.Score{s}(iFile,:) =  nanmean(smooth2a(temp(:,1:maxRank),smooth_2,0),2);
                    
                    temp = smooth2a(D_.AssBestAssem.Nonmembers.ScoreRanks{s}(:,:,iFile),0,smooth_) - ...
                           smooth2a(D_.AssBestAssem.Nonmembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile),0,smooth_);
                    D_.AssBestAssem.NonmembersSubtracted.Score{s}(iFile,:) =  nanmean(smooth2a(temp(:,1:maxRank),smooth_2,0),2);
                    
                    temp = smooth2a(D_.AssBestAssem.JointMembers.ScoreRanks{s}(:,:,iFile),0,smooth_) - ...
                           smooth2a(D_.AssBestAssem.JointMembers_noNoiseCorrs.ScoreRanks{s}(:,:,iFile),0,smooth_);
                    D_.AssBestAssem.JointMembersSubtracted.Score{s}(iFile,:) =  nanmean(smooth2a(temp(:,1:maxRank),smooth_2,0),2);
                end
        D_.AssBestAssem.AllSubtracted.ScoreMean{s} = nanmean(smooth2a(D_.AssBestAssem.AllSubtracted.Score{s},0,smooth_3));
        D_.AssBestAssem.AllSubtracted.ScoreSEM{s}  = nansem(smooth2a(D_.AssBestAssem.AllSubtracted.Score{s},0,smooth_3));
        
        D_.AssBestAssem.NonmembersSubtracted.ScoreMean{s} = nanmean(smooth2a(D_.AssBestAssem.NonmembersSubtracted.Score{s},0,smooth_3));
        D_.AssBestAssem.NonmembersSubtracted.ScoreSEM{s}  = nansem(smooth2a(D_.AssBestAssem.NonmembersSubtracted.Score{s},0,smooth_3));
        
        D_.AssBestAssem.MemberSubtracted.ScoreMean{s} = nanmean(smooth2a(D_.AssBestAssem.MemberSubtracted.Score{s},0,smooth_3));
        D_.AssBestAssem.MemberSubtracted.ScoreSEM{s}  = nansem(smooth2a(D_.AssBestAssem.MemberSubtracted.Score{s},0,smooth_3));
        
        D_.AssBestAssem.JointMembersSubtracted.ScoreMean{s} = nanmean(smooth2a(D_.AssBestAssem.JointMembersSubtracted.Score{s},0,smooth_3));
        D_.AssBestAssem.JointMembersSubtracted.ScoreSEM{s}  = nansem(smooth2a(D_.AssBestAssem.JointMembersSubtracted.Score{s},0,smooth_3));
end
 
%% Compare benefits of noise correlations - best only 
smooth_ = 1;
for s=1:3
    % Best single unit assembly decoding
    D_.AssBestSUs.AllSubtracted.Score{s}          = smooth2a(D_.AssBestSUs.All.Score{s},0,smooth_) - smooth2a(D_.AssBestSUs.All_noNoiseCorrs.Score{s},0,smooth_);
    D_.AssBestSUs.MemberSubtracted.Score{s}       = smooth2a(D_.AssBestSUs.Members.Score{s},0,smooth_) - smooth2a(D_.AssBestSUs.Members_noNoiseCorrs.Score{s},0,smooth_);
    D_.AssBestSUs.NonmembersSubtracted.Score{s}   = smooth2a(D_.AssBestSUs.Nonmembers.Score{s},0,smooth_) - smooth2a(D_.AssBestSUs.Nonmembers_noNoiseCorrs.Score{s},0,smooth_);
    D_.AssBestSUs.JointMembersSubtracted.Score{s} = smooth2a(D_.AssBestSUs.JointMembers.Score{s},0,smooth_) - smooth2a(D_.AssBestSUs.JointMembers_noNoiseCorrs.Score{s},0,smooth_);
    
    
    % Best assem assembly decoding
    D_.AssBestAssem.AllSubtracted.Score{s}          = smooth2a(D_.AssBestAssem.All.Score{s},0,smooth_) - smooth2a(D_.AssBestAssem.All_noNoiseCorrs.Score{s},0,smooth_);
    D_.AssBestAssem.MemberSubtracted.Score{s}       = smooth2a(D_.AssBestAssem.Members.Score{s},0,smooth_) - smooth2a(D_.AssBestAssem.Members_noNoiseCorrs.Score{s},0,smooth_);
    D_.AssBestAssem.NonmembersSubtracted.Score{s}   = smooth2a(D_.AssBestAssem.Nonmembers.Score{s},0,smooth_) - smooth2a(D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score{s},0,smooth_);
    D_.AssBestAssem.JointMembersSubtracted.Score{s} = smooth2a(D_.AssBestAssem.JointMembers.Score{s},0,smooth_) - smooth2a(D_.AssBestAssem.JointMembers_noNoiseCorrs.Score{s},0,smooth_);
    
    % Best single unit assembly decoding (collapse)
    D_.AssBestSUs.AllSubtracted.ScoreMean{s} = nanmean(D_.AssBestSUs.AllSubtracted.Score{s});
    D_.AssBestSUs.AllSubtracted.ScoreSEM{s}  = nansem(D_.AssBestSUs.AllSubtracted.Score{s});
    
    D_.AssBestSUs.NonmembersSubtracted.ScoreMean{s} = nanmean(D_.AssBestSUs.NonmembersSubtracted.Score{s});
    D_.AssBestSUs.NonmembersSubtracted.ScoreSEM{s}  = nansem(D_.AssBestSUs.NonmembersSubtracted.Score{s});
    
    D_.AssBestSUs.MemberSubtracted.ScoreMean{s} = nanmean(D_.AssBestSUs.MemberSubtracted.Score{s});
    D_.AssBestSUs.MemberSubtracted.ScoreSEM{s}  = nansem(D_.AssBestSUs.MemberSubtracted.Score{s});
    
    D_.AssBestSUs.JointMembersSubtracted.ScoreMean{s} = nanmean(D_.AssBestSUs.JointMembersSubtracted.Score{s});
    D_.AssBestSUs.JointMembersSubtracted.ScoreSEM{s}  = nansem(D_.AssBestSUs.JointMembersSubtracted.Score{s});
    
    % Best assem assembly decoding (collapse)
    D_.AssBestAssem.AllSubtracted.ScoreMean{s} = nanmean(D_.AssBestAssem.AllSubtracted.Score{s});
    D_.AssBestAssem.AllSubtracted.ScoreSEM{s}  = nansem(D_.AssBestAssem.AllSubtracted.Score{s});
    
    D_.AssBestAssem.NonmembersSubtracted.ScoreMean{s} = nanmean(D_.AssBestAssem.NonmembersSubtracted.Score{s});
    D_.AssBestAssem.NonmembersSubtracted.ScoreSEM{s}  = nansem(D_.AssBestAssem.NonmembersSubtracted.Score{s});
    
    D_.AssBestAssem.MemberSubtracted.ScoreMean{s} = nanmean(D_.AssBestAssem.MemberSubtracted.Score{s});
    D_.AssBestAssem.MemberSubtracted.ScoreSEM{s}  = nansem(D_.AssBestAssem.MemberSubtracted.Score{s});
    
    D_.AssBestAssem.JointMembersSubtracted.ScoreMean{s} = nanmean(D_.AssBestAssem.JointMembersSubtracted.Score{s});
    D_.AssBestAssem.JointMembersSubtracted.ScoreSEM{s}  = nansem(D_.AssBestAssem.JointMembersSubtracted.Score{s});
end

%% Find closest performing rank curve for each assembly
for iFile =1:length(D)
    for s=1:3
        D_.AssReal.ClosestRank_BestSU{s}{iFile}=[];
        D_.AssReal.ClosestRank_BestAssem{s}{iFile}=[];
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
            
            ranks_  = D{iFile}.D_.AssBestAssem.Members.ScoreRanks{s}(nU,:);
            try
                D_.AssReal.ClosestRank_BestAssem{s}{iFile}(iAss,1) =  min(FindClosestIndex(ranks_,Score));
            catch
                D_.AssReal.ClosestRank_BestAssem{s}{iFile}(iAss,1) =  NaN;
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
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestAssem.Members.CVE{s}(:,AssSize)];
    end
    
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{s})
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
    title('Best Ensemble Assemblies')
    xlabel('Time (s)')

    
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
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestAssem.Nonmembers.CVE{s}(:,2)];
    end
    
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{s})
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
    title('Optimal assemblies')
    xlabel('Time (s)')

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
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestAssem.JointMembers.CVE{s}(:,2)];
    end
    
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{s})
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
    title('Optimal assemblies')
    xlabel('Time (s)')

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
    
    % Plot best ensemble assemblies
    subplot(1,2,2); hold on
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestAssem.Nonmembers.CVE{s}(:,AssSize)];
    end
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{1},0.6)
%     plot(tb,nanmean(x,2),'color',col_{1},'LineWidth',3)

    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestAssem.Members.CVE{s}(:,AssSize)];
    end
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{2},0.6)
%     plot(tb,nanmean(x,2),'color',col_{2},'LineWidth',3)
    
    x=[];
    for iFile = 1:length(D)
        x=[x,D{iFile}.D_.AssBestAssem.JointMembers.CVE{s}(:,AssSize)];
    end
    ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_{3},0.6)
%     plot(tb,nanmean(x,2),'color',col_{3},'LineWidth',3)

    plot([min(tb) max(tb)],[0.5 0.5],':k','HandleVisibility','off')
    axis([min(tb) max(tb) 0.4 1])
    title('Best Ensemble Assemblies')
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
            x=[x,D{iFile}.D_.AssBestSUs.JointMembers_noNoiseCorrs.CVE{s}(:,i) - D{iFile}.D_.AssBestSUs.JointMembers.CVE{s}(:,i)];
        end
        
        %         ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(size_,:))
        plot(tb,nanmean(x,2),'color',[col_(size_,:),0.6],'LineWidth',1.5)
        
    end
    plot([min(tb) max(tb)],[0 0],':k')
    axis([min(tb) max(tb) -0.2 0.2])
end
%% Plot best ensemble assemblies by size - overlaid
figure('Name','Best Assem (trial noise intact)'); hold on
col_ = (jet(20));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for size_=1:20
        i=21-size_;
        x=[];
        for iFile = 1:length(D)
            try
                x=[x,D{iFile}.D_.AssBestAssem.JointMembers.CVE{s}(:,i)];
            end
        end
        
        %         ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(size_,:))
        plot(tb,nanmean(x,2),'color',[col_(size_,:),0.6],'LineWidth',1.5)
        
    end
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
end

figure('Name','Best Assem (trial noise shuffled)'); hold on
col_ = flipud(jet(20));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for size_=1:20
        i=21-size_;
        x=[];
        for iFile = 1:length(D)
            try
                x=[x,D{iFile}.D_.AssBestAssem.JointMembers_noNoiseCorrs.CVE{s}(:,i)];
            end
        end
        
        %         ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(size_,:))
        plot(tb,nanmean(x,2),'color',[col_(size_,:),0.6],'LineWidth',1.5)
        
    end
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
end

figure('Name','Best Assems (effect of removing trial correlations)'); hold on
col_ = flipud(jet(20));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for size_=1:20
        i=21-size_;
        x=[];
        for iFile = 1:length(D)
            try
%             x=[x,D{iFile}.D_.AssBestAssem.All_noNoiseCorrs.CVE{s}(:,i) - D{iFile}.D_.AssBestAssem.All.CVE{s}(:,i)];
%             x=[x,D{iFile}.D_.AssBestAssem.Members_noNoiseCorrs.CVE{s}(:,i) - D{iFile}.D_.AssBestAssem.Members.CVE{s}(:,i)];
            x=[x,D{iFile}.D_.AssBestAssem.Nonmembers_noNoiseCorrs.CVE{s}(:,i) - D{iFile}.D_.AssBestAssem.Nonmembers.CVE{s}(:,i)];
%             x=[x,D{iFile}.D_.AssBestAssem.JointMembers_noNoiseCorrs.CVE{s}(:,i) - D{iFile}.D_.AssBestAssem.JointMembers.CVE{s}(:,i)];
            end
        end
        
        %         ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(size_,:))
        plot(tb,nanmean(x,2),'color',[col_(size_,:),0.6],'LineWidth',1.5)
        
    end
    plot([min(tb) max(tb)],[0 0],':k')
    axis([min(tb) max(tb) -0.2 0.2])
end
%% *** Plot best single unit assemblies by ranks - overlaid
figure('Name','Best Single units (including trial noise)'); hold on
col_ = flipud(copper(20));
tb = (1:length(tbAll)*2)*bw;
for s=1%:3
%     subplot(1,3,s);
    hold on
   
    for rank_=2:20
        
        x=[];
        for iFile = 1:length(D)
            try
                x=[x,D{iFile}.D_.AssBestSUs.All.CVE{s}(:,rank_)];
%             x=[x,D{iFile}.D_.AssBestAssem.All.CVE{s}(:,rank_)];
            end
        end
        
%         ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(rank_-1,:))
%         ciplot(nanmean(smooth2a(x,5,0),2)+nansem(smooth2a(x,5,0),2), ...
%                nanmean(smooth2a(x,5,0),2)-nansem(smooth2a(x,5,0),2),tb,...
%                col_(rank_-1,:),0.5)
        plot(tb,nanmean(smooth2a(x,5,0),2),'color',col_(rank_-1,:),'LineWidth',2)
%         area(tb,nanmean(smooth2a(x,5,0),2),'FaceColor',col_(rank_-1,:),'BaseValue',0.5,'LineStyle','none','FaceAlpha',0.25
        axis([min(tb) max(tb) 0.4 0.9])
        plot([1 3],[0.7 0.7],'k','LineWidth',3)
        plot([1 1],[0.7 0.8],'k','LineWidth',3)
         plot([5 5],[0.6 0.8],'g','LineWidth',3)
         plot([10 10],[0.52 0.9],'w','LineWidth',3)
        plot([15 15],[0.6 0.8],'r','LineWidth',3)
    end
end
 plot([min(tb) max(tb)],[0.5 0.5],':k')
 axis off
%% *** Ranked best performing single units
figure; hold on
col_ = copper(20);

for s = 1%:2
x_ = [];
for iFile =1:length(fileList)
   x = D{iFile}.D_.Units.All.SortedScore{s};
   x =[x; nan(20-length(x),1)];
%     x=cumsum(x);   x = x./max(x);
   x_ = [x_,flipud(x(1:20))];
end

% plot(x_)
% ciplot(nanmean(x_,2)+nansem(x_,2),nanmean(x_,2)-nansem(x_,2),1:20,'k',0.2)
plot(1:20,nanmean(x_,2),'-.k','LineWidth',1.5)
scatter([2,10,20],nanmean(x_([2,10,20],:),2),30,col_([2,10,20],:),'filled','o')
errorbar(2,nanmean(x_(2,:)),nansem(x_(2,:)),'Color',col_(2,:),'LineWidth',1.5)
errorbar(10,nanmean(x_(10,:)),nansem(x_(10,:)),'Color',col_(10,:),'LineWidth',1.5)
errorbar(20,nanmean(x_(20,:)),nansem(x_(20,:)),'Color',col_(20,:),'LineWidth',1.5)

plot([1 20],[0.5 0.5],':k')
axis([0 21 0.45 1])
set(gca,'xtick',[2 10 20],'XTickLabel',{'20^t^h','10^t^h','2^n^d'},'ytick',[0.5 0.75 1])
end
 %%
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
            x=[x,full(smooth2a(D{iFile}.D_.AssBestSUs.All_noNoiseCorrs.CVE{s}(:,rank_),1,0)-smooth2a(D{iFile}.D_.AssBestSUs.All.CVE{s}(:,rank_),1,0))];
        end
        
        ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(rank_-1,:))
        plot([min(tb) max(tb)],[0 0],':k')
        axis([min(tb) max(tb) -0.4 0.4])
    end
end
%% Plot best ensemble assemblies by ranks - overlaid
figure; hold on
col_ = flipud(jet(9));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for rank_=2:10
        
        x=[];
        for iFile = 1:length(D)
            try
                x=[x,D{iFile}.D_.AssBestAssem.All.CVE{s}(:,rank_)];
            catch
                x=[x,nan(402,1)];
            end
        end
        
        ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(rank_-1,:))
        %             plot(tb,nanmean(x,2),'color',col_(rank_-1,:))
        
    end
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
end

figure; hold on
col_ = flipud(jet(9));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for rank_=2:10
        
        x=[];
        for iFile = 1:length(D)
            try
                x=[x,D{iFile}.D_.AssBestAssem.All_noNoiseCorrs.CVE{s}(:,rank_)];
            catch
                x=[x,nan(402,1)]
            end
        end
        
        ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(rank_-1,:))
        %             plot(tb,nanmean(x,2),'color',col_(rank_-1,:))
        
    end
    plot([min(tb) max(tb)],[0.5 0.5],':k')
    axis([min(tb) max(tb) 0.4 1])
end

figure; hold on
col_ = flipud(jet(9));
tb = (1:length(tbAll)*2)*bw;
for s=1:3
    subplot(1,3,s); hold on
    for rank_=2:10
        
        x=[];
        for iFile = 1:length(D)
            try
%                 x=[x,D{iFile}.D_.AssBestAssem.All_noNoiseCorrs.CVE{s}(:,rank_) - D{iFile}.D_.AssBestAssem.All.CVE{s}(:,rank_)];
%                 x=[x,D{iFile}.D_.AssBestAssem.Members_noNoiseCorrs.CVE{s}(:,rank_) - D{iFile}.D_.AssBestAssem.Members.CVE{s}(:,rank_)];
%                 x=[x,D{iFile}.D_.AssBestAssem.Nonmembers_noNoiseCorrs.CVE{s}(:,rank_) - D{iFile}.D_.AssBestAssem.Nonmembers.CVE{s}(:,rank_)];
                x=[x,D{iFile}.D_.AssBestAssem.JointMembers_noNoiseCorrs.CVE{s}(:,rank_) - D{iFile}.D_.AssBestAssem.JointMembers.CVE{s}(:,rank_)];
            catch
                x=[x,nan(402,1)];
            end
        end
        
        ciplot(nanmean(x,2)+nansem(x,2),nanmean(x,2)-nansem(x,2),tb,col_(rank_-1,:))
        %             plot(tb,nanmean(x,2),'color',col_(rank_-1,:))
        
    end
    plot([min(tb) max(tb)],[0 0],':k')
    axis([min(tb) max(tb) -0.4 0.5])
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
    
    subplot(1,2,2); hold on
    x = cell2mat(D_.AssReal.ClosestRank_BestAssem{s}');
    x = histc(x,bins); x=x./sum(x);
    plot([0 10],[s s]-1,':k','HandleVisibility','off','LineWidth',1.5)
    stairs(bins-0.5,s-1+x,'Color',col_{s},'LineWidth',1.5)
    title('Best Ensemble aggregation')
    xlabel('Closest ranked synthetic assembly')
    axis([0.5 10.5 0 3])
    set(gca,'YTick',[])
end
legend(Areas,'Orientation','horizontal');legend boxoff
%% *** Plot closest rank histograms
col_ = flipud(jet(10));
x_= [];
grp_ = [];
figure; hold on
bins = 1:10;
% col_ ={'b','r','g'};
for s=1:3
    subplot(3,1,s);hold on
    x = cell2mat( D_.AssReal.ClosestRank_BestSU{s}');
%     x = cell2mat( D_.AssReal.ClosestRank_BestAssem{s}');
    xm = nanmean(x);
%     xe=iqr(x);
    xe=nansem(x);
    x_ = [x_;x(:)];
    grp_ = [grp_;repmat(s,numel(x),1)];
    x = histc(x,bins); x=x./sum(x);
%     plot([0 10],[s s]-1,':k','HandleVisibility','off','LineWidth',1.5)
%     stairs(bins-0.5,s-1+x,'Color',col_{s},'LineWidth',1.5)

    for iBin = 1:length(bins)
        bar(bins(iBin),x(iBin),'FaceColor',col_(iBin,:),'LineStyle','none','LineWidth',1.5)
    end
    
% % % % %     bar(bins,x,'FaceColor',col_{s},'LineStyle','none','LineWidth',1.5)
errorbar_x(xm,0.35,xe,'k')

%     title('Best Single Unit aggregation')
%     ylabel('Fraction of Assemblies')
%     xlabel('Closest ranked synthetic assembly')
    axis([0.5 max(bins)+0.5 0 0.4])
    if s==3
        set(gca,'YTick',[0 0.4],'Xtick', [1 5 10],'XtickLabel',{'Best','5^t^h','10^t^h'},'XTickLabelRotation',90)
    else
        set(gca,'YTick',[0 0.4],'Xtick', [])
    end
   
end

% legend(Areas,'Orientation','horizontal');legend boxoff
[T,~,stats] = kruskalwallis(x_,grp_);
multcompare(stats)
%% Get performance
for iFile = 1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
%% Batch process behavioural performance
    L = [length(t.Short.ChoicePress_LeftCorrect)./(length(t.Short.ChoicePress_LeftCorrect)+length(t.Short.ChoicePress_LeftError)) , ...
        length(t.Medium.ChoicePress_LeftCorrect)./(length(t.Medium.ChoicePress_LeftCorrect)+length(t.Medium.ChoicePress_LeftError)) , ...
        length(t.Long.ChoicePress_LeftCorrect)./(length(t.Long.ChoicePress_LeftCorrect)+length(t.Long.ChoicePress_LeftError))];
    
    
    R = [length(t.Short.ChoicePress_RightCorrect)./(length(t.Short.ChoicePress_RightCorrect)+length(t.Short.ChoicePress_RightError)) , ...
        length(t.Medium.ChoicePress_RightCorrect)./(length(t.Medium.ChoicePress_RightCorrect)+length(t.Medium.ChoicePress_RightError)) , ...
        length(t.Long.ChoicePress_RightCorrect)./(length(t.Long.ChoicePress_RightCorrect)+length(t.Long.ChoicePress_RightError))];
    
    
    
    C_ =  [length(t.Short.ChoicePress_LeftCorrect)  + length(t.Short.ChoicePress_RightCorrect),...
        length(t.Medium.ChoicePress_LeftCorrect) + length(t.Medium.ChoicePress_RightCorrect),...
        length(t.Long.ChoicePress_LeftCorrect)   + length(t.Long.ChoicePress_RightCorrect)];
    
    E_ =  [length(t.Short.ChoicePress_LeftError)    + length(t.Short.ChoicePress_RightError),...
        length(t.Medium.ChoicePress_LeftError)   + length(t.Medium.ChoicePress_RightError),...
        length(t.Long.ChoicePress_LeftError)     + length(t.Long.ChoicePress_RightError)];
    
    D_.Accuracy(iFile,:) = (L+R)./2;
    D_.pCorr(iFile,:) = C_./(C_+E_);
    D_.AboveChance(iFile,:) = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;
    clear L R C_ E_
    
end
clear ERRORtrangeleft_choice ERRORtrangeleft_sample ERRORtrangeright_choice ERRORtrangeright_sample HPcells HPtonic PFCtonic trangeleft_choice trangeleft_sample trangeright_choice trangeright_sample)

%% closest rank assems vs. performance
figure;
for iArea=1:3
    temp = D_.AssReal.ClosestRank_BestSU{iArea};
    tempidx = isempty_cell(temp);
    
%     y = mean(D_.Accuracy(~tempidx,:),2);
    y = D_.Accuracy(~tempidx,2);
    x = cellfun(@sum,temp(~tempidx));
    subplot(1,3,iArea); hold on
    
scatter(x,y)
end

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


for iArea =1:3
    
end

    %% ************** Noise correlations on top 
%% Plot best ensemble and SU decoding from all units
x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    m = D_.AssBestAssem.All.Score_mean{s}(2:end);
    e = D_.AssBestAssem.All.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = D_.AssBestSUs.All.Score_mean{s}(2:end);
    e = D_.AssBestSUs.All.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
    m = D_.AssBestAssem.All_noNoiseCorrs.Score_mean{s}(2:end);
    e = D_.AssBestAssem.All_noNoiseCorrs.Score_SEM{s}(2:end);
    plot(x,m+e,'color',[0.6 0.6 0.6],'LineWidth',1.5)
    plot(x,m-e,'color',[0.6 0.6 0.6],'LineWidth',1.5)
    
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
    axis([2 20 0.5 1])     
end
legend({'Best Ensemble Aggregation','Best Single Unit Aggregation'},'Location','southeast');legend boxoff
clear x x_
%% Plot best ensemble and SU decoding from non-member units
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    m = D_.AssBestAssem.Nonmembers.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Nonmembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = D_.AssBestSUs.Nonmembers.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Nonmembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
    m = D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Nonmembers_noNoiseCorrs.Score_SEM{s}(2:end);
    plot(x,m+e,'color',[0.6 0.6 0.6],'LineWidth',1.5)
    plot(x,m-e,'color',[0.6 0.6 0.6],'LineWidth',1.5)
    
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
    axis([2 20 0.5 1])     
end
legend({'Best Ensemble Aggregation','Best Single Unit Aggregation'},'Location','southeast');legend boxoff
clear x x_
%% Plot best ensemble and SU decoding from member units
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    m = D_.AssBestAssem.Members.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Members.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = D_.AssBestSUs.Members.Score_mean{s}(2:end);
    e = D_.AssBestSUs.Members.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
    m = D_.AssBestAssem.Members_noNoiseCorrs.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Members_noNoiseCorrs.Score_SEM{s}(2:end);
    plot(x,m+e,'color',[0.6 0.6 0.6],'LineWidth',1.5)
    plot(x,m-e,'color',[0.6 0.6 0.6],'LineWidth',1.5)
    
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
    axis([2 20 0.5 1])     
end
legend({'Best Ensemble Aggregation','Best Single Unit Aggregation'},'Location','southeast');legend boxoff
clear x x_
%% Plot best ensemble and SU decoding from joint HP-PFC units

x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    m = D_.AssBestAssem.JointMembers.Score_mean{s}(2:end);
    e = D_.AssBestAssem.JointMembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = D_.AssBestSUs.JointMembers.Score_mean{s}(2:end);
    e = D_.AssBestSUs.JointMembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
     m = D_.AssBestAssem.JointMembers_noNoiseCorrs.Score_mean{s}(2:end);
    e = D_.AssBestAssem.JointMembers_noNoiseCorrs.Score_SEM{s}(2:end);
    plot(x,m+e,'color',[0.6 0.6 0.6],'LineWidth',1.5)
    plot(x,m-e,'color',[0.6 0.6 0.6],'LineWidth',1.5)
    
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
    axis([2 20 0.5 1])     
end
legend({'Best Ensemble Aggregation','Best Single Unit Aggregation'},'Location','southeast');legend boxoff
clear x x_

%%************** 
%% *** Plot best SU decoding from non-member units
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

x = 2:20;
figure
for s=1:3
%     subplot(1,3,s); 
    figure('Name',Areas{s})
hold on
    
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
%         ylabel({'Peak decoding score';'(% Correct Decoding)'})
    elseif s==2
%         xlabel('Assembly size (no. member units)')
        legend({'Non-member units','Member units','Inter-area members'},'Location','southeast');legend boxoff

    end
%     title([Areas{s} ' assemblies'])
    plot([2 20],[0.5 0.5],':k','HandleVisibility','off')
    axis([2 20 0.4 1])     
end
clear x x_

%% Plot best ensemble decoding from all classes of units
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    m = D_.AssBestAssem.Nonmembers.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Nonmembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = D_.AssBestAssem.Members.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Members.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.9 0.6 0],0.8)
    
    m = D_.AssBestAssem.JointMembers.Score_mean{s}(2:end);
    e = D_.AssBestAssem.JointMembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0 0.6],0.8)
    
    if s==1
        ylabel({'Peak decoding score';'(% Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
        legend({'Non-member units','Member units','Joint area members'},'Location','southeast');legend boxoff

    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0.5 0.5],':k')
    axis([2 20 0 1])     
end
clear x x_
% legend({'Non-member units','Member units','Joint area members'},'Location','southeast');legend boxoff

%% *** Plot real and shuffled assemblies with best synthetics

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
    
    m = D_.AssBestAssem.Members.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Members.Score_SEM{s}(2:end);
    %     ciplot(m+e,m-e,nU,'b',0.8)
    plot(nU,m,'color',[0.9 0.6 0 0.8],'LineWidth',1.5)
    
    m = D_.AssBestAssem.Nonmembers.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Nonmembers.Score_SEM{s}(2:end);
    %     ciplot(m+e,m-e,nU,'r',0.8)
    plot(nU,m,'color',[0.6 0.6 0.6 0.8],'LineWidth',1.5)
    
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
    
    %     m = D_.AssBestSUs.JointMembers.Score_mean{s}(2:end);
    %     e = D_.AssBestSUs.JointMembers.Score_SEM{s}(2:end);
    %     ciplot(m+e,m-e,x,'g',0.8)
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
%% Plot best ensemble synthetic assemblies
x = 2:20;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    %         m = D_.AssBestAssem.All.Score_mean{s}(2:end);
    %         e = D_.AssBestAssem.All.Score_SEM{s}(2:end);
    %         ciplot(m+e,m-e,x,'k',0.8)
    
    m = D_.AssBestAssem.Members.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Members.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'b',0.8)
    
    m = D_.AssBestAssem.Nonmembers.Score_mean{s}(2:end);
    e = D_.AssBestAssem.Nonmembers.Score_SEM{s}(2:end);
    ciplot(m+e,m-e,x,'r',0.8)
    
    %         m = D_.AssBestAssem.JointMembers.Score_mean{s}(2:end);
    %         e = D_.AssBestAssem.JointMembers.Score_SEM{s}(2:end);
    %         ciplot(m+e,m-e,x,'g',0.8)
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
%% *** Plot rank surfaces - best SU
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
smooth_ = 1;
smooth_2 = 1;
smooth_3 = 1;
smooth_4 = 1;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    %     x = D_.AssBestSUs.Nonmembers.ScoreRanks_mean{s};
%     x = nanmean(D_.AssBestSUs.Nonmembers.ScoreRanks{s},3);
%     x = smooth2a(x,smooth_,smooth_2);

    x = D_.AssBestSUs.Nonmembers.ScoreRanks{s};
    for i=1:size(x,3)
        x(:,:,i) =  smooth2a(x(:,:,i),smooth_,smooth_2);
    end
    x = nanmean(x,3);
    x =  smooth2a(x,smooth_3,smooth_4);

    mesh(x,'EdgeColor',col_{1},'FaceColor',col_{1},'FaceAlpha',0.5)
    plot3(ones(1,20),1:20,x(:,1)','Color',col_{1},'LineWidth',1.5,'HandleVisibility','off')
    
    %     x = D_.AssBestSUs.Members.ScoreRanks_mean{s};
%     x = nanmean(D_.AssBestSUs.Members.ScoreRanks{s},3);
%     x = smooth2a(x,smooth_,smooth_2);
    x = D_.AssBestSUs.Members.ScoreRanks{s};
    for i=1:size(x,3)
        x(:,:,i) =  smooth2a(x(:,:,i),smooth_,smooth_2);
    end
    x =  nanmean(x,3);
    x =  smooth2a(x,smooth_3,smooth_4);

    mesh(x,'EdgeColor',col_{2},'FaceColor',col_{2},'FaceAlpha',0.5)
    plot3(ones(1,20),1:20,x(:,1)','Color',col_{2},'LineWidth',1.5,'HandleVisibility','off')
    
    %     x = D_.AssBestSUs.JointMembers.ScoreRanks_mean{s};
%     x = nanmean(D_.AssBestSUs.JointMembers.ScoreRanks{s},3);
%     x = smooth2a(x,smooth_,smooth_2);
    x = D_.AssBestSUs.JointMembers.ScoreRanks{s};
    for i=1:size(x,3)
        x(:,:,i) =  smooth2a(x(:,:,i),smooth_,smooth_2);
    end
    x = nanmean(x,3);
    x =  smooth2a(x,smooth_3,smooth_4);

    mesh(x,'EdgeColor',col_{3},'FaceColor',col_{3},'FaceAlpha',0.5)
    plot3(ones(1,20),1:20,x(:,1)','Color',col_{3},'LineWidth',1.5,'HandleVisibility','off')
    
    set(gca,'XTick',[1,2:2:10],'XTickLabel',{'Best','2^n^d','4^t^h','6^t^h','8^t^h','10^t^h'},'XTickLabelRotation',-60)

    view(3)
    grid on
    axis([ 1 10 2 20 0.5 1])
    if s==1
%         ylabel({'Assembly size';'No. neurons'})
%         zlabel({'Peak decoding';'(% Correct classification)'})
    elseif s==3
%         xlabel('Assembly rank')
    end
%     title({Areas{s} 'Assemblies'})
        box on
end
% legend({'Non-member units','Member units','Joint area members'},'Location','southeast');legend boxoff
%% *** Plot rank surfaces - best Ensemble
col_ ={0.3*[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
[x y] = meshgrid(1:20,1:10); 

opacity = 0.5;
opacity2 = 1;
figure
for s=1:3
    subplot(1,3,s); hold on
    
    %         mesh(smooth2a(D_.AssBestAssem.All.ScoreRanks_mean{s},2,2),'EdgeColor','k')
    z = smooth2a(D_.AssBestAssem.Nonmembers.ScoreRanks_mean{s},1,1);
    m = mesh(z,'EdgeColor',col_{1},'FaceColor',col_{1},'FaceAlpha',opacity,'EdgeAlpha',opacity2);
    plot3(ones(1,20),1:20,z(:,1)','Color',col_{1},'LineWidth',1.5,'HandleVisibility','off')
    plot3(1:10,2*ones(1,10),z(2,:)','Color',col_{1},'LineWidth',1.5,'HandleVisibility','off')
    
    z = smooth2a(D_.AssBestAssem.Members.ScoreRanks_mean{s},1,1);
    mesh(z,'EdgeColor',col_{2},'FaceColor',col_{2},'FaceAlpha',opacity,'EdgeAlpha',opacity2)
    plot3(ones(1,20),1:20,z(:,1)','Color',col_{2},'LineWidth',1.5,'HandleVisibility','off')
    plot3(1:10,2*ones(1,10),z(2,:)','Color',col_{2},'LineWidth',1.5,'HandleVisibility','off')

    z = smooth2a(D_.AssBestAssem.JointMembers.ScoreRanks_mean{s},1,1);
    mesh(z,'EdgeColor',col_{3},'FaceColor',col_{3},'FaceAlpha',opacity,'EdgeAlpha',opacity2)
    plot3(ones(1,20),1:20,z(:,1)','Color',col_{3},'LineWidth',1.5,'HandleVisibility','off')
    plot3(1:10,2*ones(1,10),z(2,:)','Color',col_{3},'LineWidth',1.5,'HandleVisibility','off')

    view(3); 
    axis([ 1 10 2 20 0.5 1])
    set(gca,'XTick',[1,2:2:10],'XTickLabel',{'Best','2^n^d','4^t^h','6^t^h','8^t^h','10^t^h'},'XTickLabelRotation',-60)
    if s==1
        ylabel({'Assembly size';'No. neurons'})
        zlabel({'Peak decoding';'(% Correct classification)'})
    elseif s==3
        xlabel('Assembly rank')
    end
    title({Areas{s} 'Assemblies'})
    
%      light('Position',[0.5 0.5 2],'Style','local')
%     lighting gouraud
    
end
% legend({'Non-member units','Member units','Joint area members'},'Location','southeast');legend boxoff
%% *** Plot rank surfaces - best Ensemble with assembles
  
    
col_ = flipud(jet(10));
x=[];x_=[];
figure
for s=1:3
    subplot(1,3,s); hold on
    
        m = D_.AssBestAssem.Members.ScoreRanks_mean{s};
        e = D_.AssBestAssem.Members.ScoreRanks_sem{s};
        m(:,1) = D_.AssBestAssem.Members.Score_mean{s};
        e(:,1) = D_.AssBestAssem.Members.Score_SEM{s};
        m = smooth2a(m,1,0); e = smooth2a(e,1,0);
    for i=[1,2,5,10]
        plot(2:20,m(2:end,i),'color',col_(i,:),'LineWidth',1.5)
%         errorbar(2:20,m(2:end,i),e(2:end,i),'color',col_(i,:),'LineWidth',1.5)
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
        x_mean(rank_)=nanmean(x_(x_(:,1)==rank_,2));
        x_sem(rank_) =nansem(x_(x_(:,1)==rank_,2));
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
legend({'Best Synthetic Assembly','2^n^d best','5^t^h best','10^t^h best','Real Assemblies','Shuffled real assemblies'})
    text(10,0.52,'Chance','HorizontalAlignment','right')

legend boxoff
clear x x_
%% *** Plot rank surfaces - best SU with assemblies
col_ = flipud(jet(10));
figure
for s=1:3
    subplot(1,3,s); hold on
    
        m = D_.AssBestAssem.Members.ScoreRanks_mean{s};
        e = D_.AssBestAssem.Members.ScoreRanks_sem{s};
        m(:,1) = D_.AssBestAssem.Members.Score_mean{s};
        e(:,1) = D_.AssBestAssem.Members.Score_SEM{s};
        m = smooth2a(m,1,0); e = smooth2a(e,1,0);
        


%     for i=1:10
%         x__ = 2:20;
%         x__=[x__,fliplr(x__)];
%         y_ = [m(2:end,i)',fliplr(m(2:end,i+1)')];
%         idx = isnan(y_);
%         x__(idx)=[];y_(idx)=[];
%         patch( x__,y_,col_(i,:),'LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off')
%     end
%     clear x__ y_ idx
    for i=1:10
        if ismember(i,[1,5])
            plot(2:20,m(2:end,i),'color',col_(i,:),'LineWidth',1.5)
        else
            plot(2:20,m(2:end,i),'color',col_(i,:),'LineWidth',1.5,'HandleVisibility','off')
        end
%         errorbar(2:20,m(2:end,i),e(2:end,i),'color',col_(i-1,:),'LineWidth',1.5)
    end
    
    x=[];x_=[];

    for iFile =1:length(D)
        x = [x; D{iFile}.D_.AssReal.score{s}];
        x_= [x_;D{iFile}.D_.AssReal.score_shuffled{s}];
    end
    
    % plot real assemblies
    ranks = 1:20;
    for rank_ = ranks
        x_mean(rank_) = nanmean(x(x(:,1)==rank_,2));
        x_sem(rank_)  = nansem( x(x(:,1)==rank_,2));
    end
    x_sem(isnan(x_mean))=[];
    ranks(isnan(x_mean))=[];
    x_mean(isnan(x_mean))=[];
    errorbar(ranks,x_mean,x_sem,'.k','LineWidth',1.5)
%     scatter(ranks,x_mean,'.k','filled','HandleVisibility','off')
    scatter(x(:,1),x(:,2),12,[0 0 0],'o','HandleVisibility','off')
%     scatter(x_(:,1),x_(:,2),5,'.k','HandleVisibility','off')
    
    mdl{s} = fitlm(x(:,1),x(:,2));
    pFit(s)   = mdl{s}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s)   =  mdl{s}.Rsquared.Adjusted;
    slope(s) = mdl{s}.Coefficients{2,1};
    intercept(s) = mdl{s}.Coefficients{1,1};
    plot(mdl{s}.VariableInfo.Range{1},(mdl{s}.VariableInfo.Range{1})* slope(s)+ intercept(s),'k','LineWidth',1.5)
%     h = plot(mdl{s});
%     h(1).Marker = 'none';
%     h(1).Marker = 'none';
%     h(2).LineWidth=1.5;
%     h(2).LineStyle='-';
%     h(2).Color = 'k';
%     h(3).Color = 'k';
%     h(4).Color = 'k';
    
    % plot shuffled assemblies
    ranks = 1:20;
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
    if s<3
        axis([1 20 0.4 1])
    else
        axis([1 20 0.4 1])
    end
    
%         title({Areas{s}, 'Assemblies'})

      if s==1
%         ylabel({'Peak decoding';'(% Correct classification)'})
      elseif s==2
%         xlabel('Assembly size (No. Units)')
      end
    set(gca,'YTick',[0.5,0.75,1])
end
% legend({'Best Synthetic Assembly','5^t^h best','Real Assemblies','Shuffled real assemblies'})
%     text(10,0.52,'Chance','HorizontalAlignment','right')

% legend boxoff
clear x x_

%% *** Plot rank surfaces - best SU with assemblies - with shaded regions
col_ = flipud(jet(10));
x=[];x_=[];
figure
for s=1:3
    subplot(1,3,s); hold on
    
        m = D_.AssBestSUs.Members.ScoreRanks_mean{s};
        e = D_.AssBestSUs.Members.ScoreRanks_sem{s};
        m(:,1) = D_.AssBestSUs.Members.Score_mean{s};
        e(:,1) = D_.AssBestSUs.Members.Score_SEM{s};
        m = smooth2a(m,1,0); e = smooth2a(e,1,0);
        


    for i=1:6
        x__ = 2:20;
        x__=[x__,fliplr(x__)];
        y_ = [m(2:end,i)',fliplr(m(2:end,i+1)')];
        idx = isnan(y_);
        x__(idx)=[];y_(idx)=[];
        patch( x__,y_,col_(i,:),'LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off')
    end
    clear x__ y_ idx
    for i=1:7
        if ismember(i,[1,5])
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
    ranks = 1:15;
    for rank_ = ranks
        x_mean(rank_)=nanmean(x(x(:,1)==rank_,2));
        x_sem(rank_) =nansem(x(x(:,1)==rank_,2));
    end
    x_sem(isnan(x_mean))=[];
    ranks(isnan(x_mean))=[];
    x_mean(isnan(x_mean))=[];
    errorbar(ranks,x_mean,x_sem,'.k','LineWidth',1.5)
%     scatter(ranks,x_mean,'.k','filled','HandleVisibility','off')
    scatter(x(:,1),x(:,2),12,[0 0 0],'o','HandleVisibility','off')
%     scatter(x_(:,1),x_(:,2),5,'.k','HandleVisibility','off')
    
    mdl{s} = fitlm(x(:,1),x(:,2));
    pFit(s)   = mdl{s}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s)   =  mdl{s}.Rsquared.Adjusted;
    slope(s) = mdl{s}.Coefficients{2,1};
    intercept(s) = mdl{s}.Coefficients{1,1};
    plot(mdl{s}.VariableInfo.Range{1},(mdl{s}.VariableInfo.Range{1})* slope(s)+ intercept(s),'k','LineWidth',1.5)
%     h = plot(mdl{s});
%     h(1).Marker = 'none';
%     h(1).Marker = 'none';
%     h(2).LineWidth=1.5;
%     h(2).LineStyle='-';
%     h(2).Color = 'k';
%     h(3).Color = 'k';
%     h(4).Color = 'k';
    
    % plot shuffled assemblies
    ranks = 1:15;
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
    if s<3
        axis([1 7 0.4 1])
    else
        axis([1 16 0.4 1])
    end
    
%         title({Areas{s}, 'Assemblies'})

      if s==1
%         ylabel({'Peak decoding';'(% Correct classification)'})
      elseif s==2
%         xlabel('Assembly size (No. Units)')
      end
    set(gca,'YTick',[0.5,0.75,1])
end
% legend({'Best Synthetic Assembly','5^t^h best','Real Assemblies','Shuffled real assemblies'})
%     text(10,0.52,'Chance','HorizontalAlignment','right')

% legend boxoff
clear x x_
%% *** Plot rank surfaces - best SU with assemblies
col_ = flipud(jet(10));
figure
for s=1:3
    subplot(1,3,s); hold on
    
        m = D_.AssBestSUs.Members.ScoreRanks_mean{s};
        e = D_.AssBestSUs.Members.ScoreRanks_sem{s};
        m(:,1) = D_.AssBestSUs.Members.Score_mean{s};
        e(:,1) = D_.AssBestSUs.Members.Score_SEM{s};
        m = smooth2a(m,1,0); e = smooth2a(e,1,0);
        


%     for i=1:10
%         x__ = 2:20;
%         x__=[x__,fliplr(x__)];
%         y_ = [m(2:end,i)',fliplr(m(2:end,i+1)')];
%         idx = isnan(y_);
%         x__(idx)=[];y_(idx)=[];
%         patch( x__,y_,col_(i,:),'LineStyle','none','FaceAlpha',0.2,'HandleVisibility','off')
%     end
%     clear x__ y_ idx
    for i=1:10
        if ismember(i,[1,5])
            plot(2:20,m(2:end,i),'color',col_(i,:),'LineWidth',1.5)
        else
            plot(2:20,m(2:end,i),'color',col_(i,:),'LineWidth',1.5,'HandleVisibility','off')
        end
%         errorbar(2:20,m(2:end,i),e(2:end,i),'color',col_(i-1,:),'LineWidth',1.5)
    end
    
    x=[];x_=[];

    for iFile =1:length(D)
        x = [x; D{iFile}.D_.AssReal.score{s}];
        x_= [x_;D{iFile}.D_.AssReal.score_shuffled{s}];
    end
    
    % plot real assemblies
    ranks = 1:20;
    for rank_ = ranks
        x_mean(rank_) = nanmean(x(x(:,1)==rank_,2));
        x_sem(rank_)  = nansem( x(x(:,1)==rank_,2));
    end
    x_sem(isnan(x_mean))=[];
    ranks(isnan(x_mean))=[];
    x_mean(isnan(x_mean))=[];
    errorbar(ranks,x_mean,x_sem,'.k','LineWidth',1.5)
%     scatter(ranks,x_mean,'.k','filled','HandleVisibility','off')
    scatter(x(:,1),x(:,2),12,[0 0 0],'o','HandleVisibility','off')
%     scatter(x_(:,1),x_(:,2),5,'.k','HandleVisibility','off')
    
    mdl{s} = fitlm(x(:,1),x(:,2));
    pFit(s)   = mdl{s}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s)   =  mdl{s}.Rsquared.Adjusted;
    slope(s) = mdl{s}.Coefficients{2,1};
    intercept(s) = mdl{s}.Coefficients{1,1};
    plot(mdl{s}.VariableInfo.Range{1},(mdl{s}.VariableInfo.Range{1})* slope(s)+ intercept(s),'k','LineWidth',1.5)
%     h = plot(mdl{s});
%     h(1).Marker = 'none';
%     h(1).Marker = 'none';
%     h(2).LineWidth=1.5;
%     h(2).LineStyle='-';
%     h(2).Color = 'k';
%     h(3).Color = 'k';
%     h(4).Color = 'k';
    
    % plot shuffled assemblies
    ranks = 1:20;
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
    if s<3
        axis([1 20 0.4 1])
    else
        axis([1 20 0.4 1])
    end
    
%         title({Areas{s}, 'Assemblies'})

      if s==1
%         ylabel({'Peak decoding';'(% Correct classification)'})
      elseif s==2
%         xlabel('Assembly size (No. Units)')
      end
    set(gca,'YTick',[0.5,0.75,1])
end
% legend({'Best Synthetic Assembly','5^t^h best','Real Assemblies','Shuffled real assemblies'})
%     text(10,0.52,'Chance','HorizontalAlignment','right')

% legend boxoff
clear x x_

%% Plot average benefit of noise correlations
x = 2:20;
ymax = 20;
figure('Name','All units')
for s=1:3
    subplot(1,3,s); hold on
    
    m = -100*D_.AssBestAssem.AllSubtracted.ScoreMean{s}(2:end);
    e = -100*D_.AssBestAssem.AllSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = -100*D_.AssBestSUs.AllSubtracted.ScoreMean{s}(2:end);
    e = -100*D_.AssBestSUs.AllSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
  
    
    if s==1
        ylabel({'Benefit of shuffled noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -ymax ymax])     
end
legend({'Best Ensemble Aggregation','Best Single Unit Aggregation'},'Location','southeast');legend boxoff

figure('Name','Nonmember units')
for s=1:3
    subplot(1,3,s); hold on
    
    m = -100*D_.AssBestAssem.NonmembersSubtracted.ScoreMean{s}(2:end);
    e = -100*D_.AssBestAssem.NonmembersSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = -100*D_.AssBestSUs.NonmembersSubtracted.ScoreMean{s}(2:end);
    e = -100*D_.AssBestSUs.NonmembersSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
  
    
    if s==1
        ylabel({'Benefit of shuffled noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -ymax ymax])     
end
legend({'Best Ensemble Aggregation','Best Single Unit Aggregation'},'Location','southeast');legend boxoff

figure('Name','Member units')
for s=1:3
    subplot(1,3,s); hold on
    
    m = -100*D_.AssBestAssem.MemberSubtracted.ScoreMean{s}(2:end);
    e = -100*D_.AssBestAssem.MemberSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = -100*D_.AssBestSUs.MemberSubtracted.ScoreMean{s}(2:end);
    e = -100*D_.AssBestSUs.MemberSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
  
    
    if s==1
        ylabel({'Benefit of shuffled noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -ymax ymax])     
end
legend({'Best Ensemble Aggregation','Best Single Unit Aggregation'},'Location','southeast');legend boxoff
figure('Name','Joint-area units')
for s=1:3
    subplot(1,3,s); hold on
    
    m = -100*D_.AssBestAssem.JointMembersSubtracted.ScoreMean{s}(2:end);
    e = -100*D_.AssBestAssem.JointMembersSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,[0.6 0.6 0.6],0.8)
    
    m = -100*D_.AssBestSUs.JointMembersSubtracted.ScoreMean{s}(2:end);
    e = -100*D_.AssBestSUs.JointMembersSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,'k',0.8)
    
  
    
    if s==1
        ylabel({'Benefit of noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        xlabel('Assembly size (no. member units)')
    end
    title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -ymax ymax])     
end
legend({'Best Ensemble Aggregation','Best Single Unit Aggregation'},'Location','southeast');legend boxoff


% clear x x_
%% *** Benefit of noise 2: Best SUs local/nonmembers
x = 2:20;
ymax = 15;
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
rmpath('/Users/aleksanderdomanski/Dropbox/MATLAB/PC MATLAB Path/functions/nansuite')
clear p stats
for s=1:2
%     subplot(1,3,s); hold on
        figure('Name',Areas{s}); hold on

    m = 100*D_.AssBestSUs.NonmembersSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestSUs.NonmembersSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,col_{1},0.8)
    
    m = 100*D_.AssBestSUs.MemberSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestSUs.MemberSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,col_{2},0.8)
    
    [~,p(s,1),~,stats(s,1)] = ttest(nanmean(D_.AssBestSUs.NonmembersSubtracted.Score{s}(:,2:end)'));
	[~,p(s,2),~,stats(s,2)] = ttest(nanmean(D_.AssBestSUs.MemberSubtracted.Score{s}(:,2:end)'));
%     m = 100*D_.AssBestSUs.JointMembersSubtracted.ScoreMean{s}(2:end);
%     e = 100*D_.AssBestSUs.JointMembersSubtracted.ScoreSEM{s}(2:end);
%     ciplot(m+e,m-e,x,col_{3},0.8)
    
    if s==1
%         ylabel({'Benefit of shuffled noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        %         xlabel('Assembly size (no. member units)')
    end
%     title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -ymax 5])
%     axis([1 7 -5 ymax])

end
%% *** Benefit of noise 2: Best Assems local/nonmembers
x = 2:20;
ymax = 10;
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

for s=1:2
%     subplot(1,3,s); hold on
        figure('Name',Areas{s}); hold on

    m = 100*D_.AssBestAssem.NonmembersSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestAssem.NonmembersSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,col_{1},0.8)
    
    m = 100*D_.AssBestAssem.MemberSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestAssem.MemberSubtracted.ScoreSEM{s}(2:end);
    ciplot(m+e,m-e,x,col_{2},0.8)
    
%     m = 100*D_.AssBestAssem.JointMembersSubtracted.ScoreMean{s}(2:end);
%     e = 100*D_.AssBestAssem.JointMembersSubtracted.ScoreSEM{s}(2:end);
%     ciplot(m+e,m-e,x,col_{3},0.8)
    
    if s==1
%         ylabel({'Benefit of shuffled noise correlations';'(\Delta% Peak Correct Decoding)'})
    elseif s==2
        %         xlabel('Assembly size (no. member units)')
    end
%     title([Areas{s} ' assemblies'])
    plot([2 20],[0 0],':k')
    axis([2 20 -15 10])
%     axis([1 7 -5 ymax])

end

%% inter-area noise shuffling - best SUs
s=3;
col_ ={'b','r','g'};

% subplot(1,3,s); hold on
figure('Name',Areas{s}); hold on
for s=1:3
    plot(x,100*D_.AssBestSUs.JointMembersSubtracted.Score{s}(:,2:end),':','color',col_{s});
end

for s=1:3
    
    m = 100*D_.AssBestSUs.JointMembersSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestSUs.JointMembersSubtracted.ScoreSEM{s}(2:end);
    idx=e==0;
    ciplot(m(~idx)+e(~idx),m(~idx)-e(~idx),x(~idx), col_{s},0.8)
end
 plot([2 20],[0 0],':k')
%     axis([1 16 -5 ymax])
    axis([2 20 -ymax 5])

    x_ = []; T= [];
    for s=1:3
        y = nanmean(100*smooth2a(D_.AssBestSUs.JointMembersSubtracted.Score{s}(:,2:end),0,0),2);
%         y(isnan(y))=[];
        x_ = [x_;y];
        T = [T;s*ones(length(y),1)];
    end
     [T,~,stats] = anova1(x_,T);
multcompare(stats)
%%
rmpath('/Users/aleksanderdomanski/Dropbox/MATLAB/PC MATLAB Path/functions/nansuite')

clear p stats
for s=1:3
    
%     m(:,s) = 100*D_.AssBestSUs.JointMembersSubtracted.ScoreMean{s}(2:end);
    [~,p(s),~,stats(s)]= ttest(nanmean(D_.AssBestSUs.JointMembersSubtracted.Score{s}(:,2:end)'))
end

%% inter-area noise shuffling - best assems
s=3;
col_ ={'b','r','g'};

% subplot(1,3,s); hold on
figure('Name',Areas{s}); hold on
for s=1:3
    plot(x,100*D_.AssBestAssem.JointMembersSubtracted.Score{s}(:,2:end),':','color',col_{s});
end

for s=1:1
    
    m = 100*D_.AssBestAssem.JointMembersSubtracted.ScoreMean{s}(2:end);
    e = 100*D_.AssBestAssem.JointMembersSubtracted.ScoreSEM{s}(2:end);
    idx=e==0;
    ciplot(m(~idx)+e(~idx),m(~idx)-e(~idx),x(~idx), col_{s},0.8)
end
 plot([2 20],[0 0],':k')
%     axis([1 16 -5 ymax])
    axis([2 20 -ymax 5])

    x_ = []; T= [];
    for s=1:3
        y = nanmean(100*smooth2a(D_.AssBestAssem.JointMembersSubtracted.Score{s}(:,2:end),0,0),2);
        y(isnan(y))=[];
        x_ = [x_;y];
        T = [T;s*ones(length(y),1)];
    end
     [T,~,stats] = anova1(x_,T);
multcompare(stats)
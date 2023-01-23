function Assemblies = MakebestSUAssems(RankIn_,ScoreIn_, Ltrials_,Rtrials_,MaxSize,DestroyNoiseCorrs,noTrials,nReps)
%@TODO: add shuffling of trial labels between cells to non- trial-resampled version
Areas = {'PFC','HP','Joint'};

InfoCriteria = 'max';
verbose = false;
Ltr = size(Ltrials_{1}{1},1);
noUnitsReal = [size(Ltrials_{1}{1},2),size(Ltrials_{2}{1},2),size(Ltrials_{3}{1},2)];
if nargin<6
    DestroyNoiseCorrs = false; % keep L or R trial identity but shuffle trial labels to remove within-trial correlations
end
if nargin<7
    noTrials = Inf; % Use all trials if noTrials is Inf
end
if nargin<8
    nReps = 10;
end
if isinf(noTrials)
    nReps = 1;
end
if DestroyNoiseCorrs
    nReps = 10;
end
%%

for s=1:3
    
    Assemblies.CVE{s}                 = nan(Ltr,MaxSize);
    Assemblies.Score{s}               = nan(MaxSize,1);
    Assemblies.Score_shuffled{s}      = nan(MaxSize,1);
    Assemblies.Score_shuffled_5pc{s}  = nan(MaxSize,1);
    Assemblies.Score_shuffled_95pc{s} = nan(MaxSize,1);
    Assemblies.ScoreRanks{s}          = nan(MaxSize,10);
    
    for i = 2:MaxSize
        % FIRST ROUND
        fprintf('>Decoding Synthetic %s Assemblies: %d units pooled (%d maximum)...\n',Areas{s},i,MaxSize)
        
        if s<3
            try
                unitIds = RankIn_{s}(1:i);
            catch
                fprintf('>Ran out of units for %s ass. size %d, rank 1\n',Areas{s},i);
                break
            end
        else
            try
                %Assumes equal weighting of assemblies between areas
                %@TODO: Implement sliding ratios
                a  = RankIn_{1}(1:ceil(i/2));
                a_ = ScoreIn_{1}(1:ceil(i/2));
                b  = RankIn_{2}(1:ceil(i/2));
                b_ = ScoreIn_{2}(1:ceil(i/2));
                if rem(i,2)~=0
                    if b_(end)>a_(end)
                        a(end)=[];
                    else
                        b(end)=[];
                    end
                end
                unitIds = [a, b+noUnitsReal(1)];
            catch
                fprintf('>Ran out of units for %s ass. size %d rank 1\n',Areas{s},i);
                break
            end
        end
        if DestroyNoiseCorrs
            [Assemblies.CVE{s}(:,i),...
                Assemblies.Score{s}(i),...
                Assemblies.Score_shuffled{s}(i),...
                Assemblies.Score_shuffled_5pc{s}(i),...
                Assemblies.Score_shuffled_95pc{s}(i)] = RunDecodeGroupsDestroyNoiseCorrs(unitIds,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
        else
            [Assemblies.CVE{s}(:,i),...
                Assemblies.Score{s}(i),...
                Assemblies.Score_shuffled{s}(i),...
                Assemblies.Score_shuffled_5pc{s}(i),...
                Assemblies.Score_shuffled_95pc{s}(i)] = RunDecodeGroups(unitIds,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
        end
        
        % Get rank-ordered next performers
        Assemblies.ScoreRanks{s}(i,1)=Assemblies.Score{s}(i);
        for rank_= 2:10
            if s<3
                idx = (i*(rank_-1)+1):i*rank_;
                if sum(idx>length(RankIn_{s}))>0
                    break
                end
                unitIds = RankIn_{s}(idx);
            else
                offset_ = (rank_-1)*i+1;
                if i==2
                    idx = offset_;
                else
                    idx = offset_:(offset_+ceil(i/2)-1);
                end
                if sum(idx>length(RankIn_{1}))>0 || sum(idx>length(RankIn_{2}))>0
                    fprintf('>Ran out of units for %s ass. size %d rank %d\n',Areas{s},i,rank_);
                    break
                end
                a  = RankIn_{1}(idx);
                a_ = ScoreIn_{1}(idx);
                b  = RankIn_{2}(idx);
                b_ = ScoreIn_{2}(idx);
                % If there's an overhanging unit, trim of the one with the lower score
                if rem(i,2)~=0
                    if b_(end)>a_(end)
                        a(end)=[];
                    else
                        b(end)=[];
                    end
                end
                unitIds = [a, b+noUnitsReal(1)];
            end
            if verbose
                fprintf('>...Ranking next best %d/%d (units: [%s])\n',rank_,10,num2str(unitIds))
            end
            if DestroyNoiseCorrs
                [~,Assemblies.ScoreRanks{s}(i,rank_)] = RunDecodeGroupsDestroyNoiseCorrs(unitIds,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
            else
                [~,Assemblies.ScoreRanks{s}(i,rank_)] = RunDecodeGroups(unitIds,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
            end

        end
    end
end
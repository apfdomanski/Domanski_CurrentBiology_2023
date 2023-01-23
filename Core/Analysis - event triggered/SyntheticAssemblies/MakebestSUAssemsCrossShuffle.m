function Assemblies = MakebestSUAssemsCrossShuffle(RankIn_,ScoreIn_, Ltrials_,Rtrials_,MaxSize,CorrsToDestroy,noTrials,nReps)
%@TODO: add shuffling of trial labels between cells to non- trial-resampled version
Areas = {'PFC','HP','Joint'};
InfoCriteria = 'max';
nReps = 50;
verbose = false;
Ltr = size(Ltrials_{1}{1},1);
noUnitsReal = [size(Ltrials_{1}{1},2),size(Ltrials_{2}{1},2),size(Ltrials_{3}{1},2)];

if  nargin<6
    %     CorrsToDestroy = [1,1,1]; % i.e. shuffle all possible combinaitons
    CorrsToDestroy = [1,1,1]; % i.e. shuffle only half of inter-area combinations
end
if nargin<7
    noTrials = Inf; % Use all trials if noTrials is Inf
end
if nargin<8
    nReps = 10;
end
if isinf(noTrials)
    nReps =1;
end
if sum(CorrsToDestroy)>0
    nReps = 10;
end

for s=find(CorrsToDestroy)
    
    Assemblies.CVE{s}                 = nan(Ltr,MaxSize);
    Assemblies.Score{s}               = nan(MaxSize,1);
    Assemblies.Score_shuffled{s}      = nan(MaxSize,1);
    Assemblies.Score_shuffled_5pc{s}  = nan(MaxSize,1);
    Assemblies.Score_shuffled_95pc{s} = nan(MaxSize,1);
    Assemblies.ScoreRanks{s}          = nan(MaxSize,10);
    
    
    for i = 2:MaxSize
        fprintf('>Constructing synthetic inter-area assemblies, shuffling %s correlations: %d units pooled (%d maximum)...\n',Areas{s},i,MaxSize)
        
        
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
                fprintf('>Ran out of units for %s ass. size %d, rank 1\n',Areas{s},i);
                break
            end
        
        
            
            [Assemblies.CVE{s}(:,i),...
                Assemblies.Score{s}(i),...
                Assemblies.Score_shuffled{s}(i),...
                Assemblies.Score_shuffled_5pc{s}(i),...
                Assemblies.Score_shuffled_95pc{s}(i)] = RunDecodeGroupsDestroyNoiseCorrsCrossArea(unitIds,nReps,Ltrials_,Rtrials_,noTrials,'max',s);
        
        % Get rank-ordered next performers
        Assemblies.ScoreRanks{s}(i,1)=Assemblies.Score{s}(i);
        
        for rank_= 2:10
               
                offset_ = (rank_-1)*i+1;
                if i==2
                    idx = offset_;
                else
                    idx = offset_:(offset_+ceil(i/2)-1);
                end
                if sum(idx>length(RankIn_{1}))>0 || sum(idx>length(RankIn_{2}))>0
                    fprintf('>Ran out of units for %s ass. size %d, rank %d\n',Areas{s},i,rank_);
                    break
                end
                a  = RankIn_{1}(idx);
                a_ = ScoreIn_{1}(idx);
                b  = RankIn_{2}(idx);
                b_ = ScoreIn_{2}(idx);
                % If there's an overhanging unit, trim of the
                % one with the lower score
                if rem(i,2)~=0
                    if b_(end)>a_(end)
                        a(end)=[];
                    else
                        b(end)=[];
                    end
                end
                unitIds = [a, b+noUnitsReal(1)];
                
                if verbose
                    fprintf('>...Ranking next best %d/%d (units: [%s])\n',rank_,10,num2str(unitIds))
                end
                
               [~,Assemblies.ScoreRanks{s}(i,rank_)] = RunDecodeGroupsDestroyNoiseCorrsCrossArea(unitIds,nReps,Ltrials_,Rtrials_,noTrials,'max',s);

        end
    end
end
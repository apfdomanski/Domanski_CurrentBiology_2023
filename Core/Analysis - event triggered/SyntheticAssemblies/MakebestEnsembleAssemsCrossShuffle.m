function Assemblies = MakebestEnsembleAssemsCrossShuffle(UnitsReal, Ltrials_,Rtrials_,MaxSize,CorrsToDestroy,noTrials,nReps)
%@TODO: add shuffling of trial labels between cells to non- trial-resampled version
Areas = {'PFC','HP','Joint'};
InfoCriteria = 'max';
nReps = 50;
noTrials = 10;
verbose = false;
Ltr = size(Ltrials_{1}{1},1);
noUnitsReal = [size(Ltrials_{1}{1},2),size(Ltrials_{2}{1},2),size(Ltrials_{3}{1},2)];

if  nargin<5
    %     CorrsToDestroy = [1,1,1]; % i.e. shuffle all possible combinaitons
    CorrsToDestroy = [1,1,1]; % i.e. shuffle only half of inter-area combinations
end
if nargin<6
    noTrials = Inf; % Use all trials if noTrials is Inf
end
if nargin<8
    nReps = 10;
end
if isinf(noTrials)
    nReps = 1;
end

for s=find(CorrsToDestroy)
    Score          = nan(1,MaxSize);
    bestCombo      = cell(1,MaxSize);
    CVE            = nan(Ltr,MaxSize);
    ScoreRanks     = nan(MaxSize,10);
    bestComboRanks = cell(MaxSize,10);
    CVERanks       = cell(1,MaxSize);
    noUnits = min([length(UnitsReal{s}),MaxSize]);
    
    Flag  = 0;
    if length(UnitsReal{s})==0
        Flag = 1;
    end
   
    if noUnits > 4 & ~Flag
        CVE=zeros(Ltr,noUnits);
        % (1) First round: all to all pairwise combinations
        Poolsize = 2;
        unitIDs  = cell(length(UnitsReal{s}),length(UnitsReal{s}));
        CVE_ = cell(length(UnitsReal{s}),length(UnitsReal{s}));
        Score_ = nan(length(UnitsReal{s}),length(UnitsReal{s}));
        fprintf('Finding optimal Synthetic %s Assembly : first round (all to all) ...\n',Areas{s})
        for ii = 1:length(UnitsReal{s})
            parfor jj = 1:length(UnitsReal{s})
%                 [ii jj]
                    if verbose
                        fprintf('...First round (pair [%d,%d] of %d units total)\n',ii,jj,length(UnitsReal{s}))
                    end
                    
                    Flag2 = 0
                    % Ignore pairs of the same unit and only run once for each pair
                   
                    if ii<jj
                        Flag2 = 1;
                    end
                    % only test inter-area pairs for joint area case
                    if s==3 & Flag2
                        Flag2  = UnitsReal{s}(ii)<=noUnitsReal(1) & UnitsReal{s}(jj)>noUnitsReal(1)
                    end
                    if Flag2
                        unitIDs{ii,jj} = UnitsReal{s}([ii,jj]);
                        [CVE_{ii,jj},Score_(ii,jj)] = RunDecodeGroupsDestroyNoiseCorrsCrossArea(unitIDs{ii,jj},nReps,Ltrials_,Rtrials_,noTrials,'max',s);
                    end
            end
        end
       
        % (2) stats for best performing pair
        idx                 = find(Score_==nanmax(nanmax(Score_)),1);
        Score(Poolsize)     = Score_(idx);
        bestCombo{Poolsize} = unitIDs{idx};
        CVE(:,Poolsize)     = CVE_{idx};
        
        % Get rank-ordered next performers
        ScoreRanks(Poolsize,1)=Score(Poolsize);
        Score__ = Score_; Score__(isnan(Score__))=0;
        [ b, ix ] = sort( Score__(:), 'descend' );
        [ rr, cc ] = ind2sub( size(Score__), ix(1:10));
        
        for rank_ = 1:10
            try
                CVERanks{Poolsize}(:,rank_)    = CVE_{rr(rank_),cc(rank_)};
                bestComboRanks{Poolsize,rank_} = unitIDs{rr(rank_),cc(rank_)};            
                ScoreRanks(Poolsize,rank_)     = Score__(rr(rank_),cc(rank_));

            catch
                CVERanks{Poolsize}(:,rank_)    = nan(Ltr,1);
                bestComboRanks{Poolsize,rank_} = unitIDs{rr(rank_),cc(rank_)};
                ScoreRanks(Poolsize,rank_)     = nan;
            end
        end
        clear b ix rr cc rank_
        
        % (3) Successive rounds: greedy aggregation of contributing units
        fprintf('Ranking sub-optimal synthetic %s assemblies...\n',Areas{s})
        for Poolsize = 3:noUnits
            
            unitIDs = setdiff(UnitsReal{s},bestCombo{Poolsize-1});
            if s==3
                % For joint areas, alternate between PFC and HP on odd and even iterations
                % odd runs: test HP additions, even runs: test PFC additions
                if rem(Poolsize,2)
                    unitIDs(unitIDs>noUnitsReal(1))=[];
                else 
                    unitIDs(unitIDs<=noUnitsReal(1))=[];
                end
            end
            if isempty (unitIDs)
                fprintf('Run out of units on Poolsize=%d.\n ',Poolsize)
                break
            end
            
            % vectors to catch 2(+)nd level deep results
            CVE_2    = zeros(Ltr,length(unitIDs));
            Score_2  = zeros(length(unitIDs),1);
            unitID_2 = zeros(length(unitIDs),1);
            
            parfor ii=1:length(unitIDs)
                if verbose
                    fprintf('...Poolsize %d (of %d), %d Units left to test\n',Poolsize,noUnits,length(unitIDs)-ii)
                end
                unitID_(ii)=unitIDs(ii);
                    [CVE_2(:,ii),Score_2(ii)] = RunDecodeGroupsDestroyNoiseCorrsCrossArea([bestCombo{Poolsize-1},unitID_(ii)],nReps,Ltrials_,Rtrials_,noTrials,'max',s);

               
            end
            idx = find(Score_2 == max(Score_2),1);
            Score(Poolsize)     = Score_2(idx);
            bestCombo{Poolsize} = [bestCombo{Poolsize-1},unitIDs(idx)];
            CVE(:,Poolsize)     = CVE_2(:,idx);
            
            [ b, ix ] = sort( Score_2, 'descend' );
            for rank_ = 1:10
                try
                    ScoreRanks(Poolsize,rank_)     = Score_2(ix(rank_));
                    bestComboRanks{Poolsize,rank_} = [bestCombo{Poolsize-1},unitIDs(ix(rank_))];
                    CVERanks{Poolsize}(:,rank_)    = CVE_2(:,ix(rank_));
                catch
                    ScoreRanks(Poolsize,rank_)     = NaN;
                    bestComboRanks{Poolsize,rank_} = [];
                    CVERanks{Poolsize}(:,rank_)    = nan(Ltr,1);
                end
            end
            clear b ix rank_
        end
    end
    
    Assemblies.Score{s}          = Score;
    Assemblies.CVE{s}            = CVE;
    Assemblies.bestCombo{s}      = bestCombo;
    Assemblies.ScoreRanks{s}     = ScoreRanks;
    Assemblies.CVERanks{s}       = CVERanks;
    Assemblies.bestComboRanks{s} = bestComboRanks;
    
end
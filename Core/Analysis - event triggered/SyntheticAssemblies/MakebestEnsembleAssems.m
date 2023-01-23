function Assemblies = MakebestEnsembleAssems(UnitsReal,Ltrials_,Rtrials_,MaxSize,JointOnly,DestroyNoiseCorrs,noTrials,nReps)
%%
% Compared decoding performance from aggregated cell assemblies.
% Search is two step: first round is an all-to-all pairwise exhaustive
% search. Best pair is chosen, then remaining units are added one-by-one to
% find best triplet, and so on up to maximum assembly size, or until
% available units run out.
%
% Outputs: Assemblies (Structure) with the following fields
% ./Score:          {no Areas}[1,MaxSize]                        Summary classifier performance for best assembly of each size (either mean or max cross-validation error)
% ./CVE:            {no Areas}[Length of trial,MaxSize]          Time-resolved classifier performance for best assembly of each size
% ./BestCombo:      {no Areas}{MaxSize}[1,Assembly size]         IDs of units combined into performance-ranked assemblies
% ./ScoreRanks:     {no Areas}[MaxSize,MaxRank]                  Summary classifier performance for ranked assembly of each size (either mean or max cross-validation error)
% ./CVERanks:       {no Areas}{MaxSize}[Length of trial,MaxRank] Time-resolved classifier performance for all assembly of each size
% ./BestComboRanks: {no Areas}{MaxSize,MaxRank}[1,Assembly size] IDs of units combined into performance-ranked assemblies

%%
if nargin<5 |isempty(JointOnly)
    JointOnly = false;
end
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
    nReps =1;
end

Areas = {'PFC','HP','Joint'};
InfoCriteria = 'max';
verbose = false;
Ltr = size(Ltrials_{1}{1},1);
noUnitsReal = [size(Ltrials_{1}{1},2),size(Ltrials_{2}{1},2),size(Ltrials_{3}{1},2)];

for s=1:3
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
    if ~JointOnly
        if s==3
            Flag =  ~(length(UnitsReal{1}) >0 & length(UnitsReal{2}) >0);
        end
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
                    if verbose
                        fprintf('...First round (pair [%d,%d] of %d units total)\n',ii,jj,length(UnitsReal{s}))
                    end
                    
                    Flag2 = 0;
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
                        if DestroyNoiseCorrs
                            [CVE_{ii,jj},Score_(ii,jj)] = RunDecodeGroupsDestroyNoiseCorrs(unitIDs{ii,jj},nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
                        else
                            [CVE_{ii,jj},Score_(ii,jj)] = RunDecodeGroups(unitIDs{ii,jj},nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
                        end
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
                if DestroyNoiseCorrs
                    [CVE_2(:,ii),Score_2(ii)] = RunDecodeGroupsDestroyNoiseCorrs([bestCombo{Poolsize-1},unitIDs(ii)],nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
                else
                    [CVE_2(:,ii),Score_2(ii)] = RunDecodeGroups([bestCombo{Poolsize-1},unitIDs(ii)],nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');
                end
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
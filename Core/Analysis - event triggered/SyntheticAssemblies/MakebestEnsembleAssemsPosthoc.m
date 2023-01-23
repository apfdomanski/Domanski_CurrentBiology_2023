function Assemblies = MakebestEnsembleAssemsPosthoc(UnitsReal,Ltrials_,Rtrials_,MaxSize,JointOnly,DestroyNoiseCorrs,noTrials,nReps,unitSelection)
    %%
    % ** Takes IDs of groupings previously identified as optimal **
    %
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
    Assemblies= struct;
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

        if noUnits >= 2 && ~Flag

            CVE_   = cell(MaxSize,10); CVE_(:) = {nan(Ltr,1)};
            Score_ = nan(MaxSize,10);

            for iSize = 2:MaxSize
                fprintf('Ass Size = %d\n',iSize)

                parfor iRank = 1:10
                    unitIDs = unitSelection{s}{iSize,iRank};
                    
                    if ~isempty(unitIDs)
                        if DestroyNoiseCorrs
                            [CVE_{iSize,iRank},Score_(iSize,iRank)] = RunDecodeGroupsDestroyNoiseCorrs(unitIDs,nReps,Ltrials_{s},Rtrials_{s},noTrials,InfoCriteria);
                        else
                            [CVE_{iSize,iRank},Score_(iSize,iRank)] = RunDecodeGroups(unitIDs,nReps,Ltrials_{s},Rtrials_{s},noTrials,InfoCriteria);
                        end
                    end
                    
                end
                
            end
            
            Assemblies.Score{s}          = Score_(:,1);
            Assemblies.CVE{s}            = [nan(Ltr,1), cell2mat(CVE_(2:MaxSize,2)')];
            Assemblies.bestCombo{s}      = unitSelection{s}(:,1);
            Assemblies.ScoreRanks{s}     = Score_;
            for iSize = 2:MaxSize
                Assemblies.CVERanks{s}{iSize}= cell2mat(CVE_(iSize,:));
            end
            Assemblies.bestComboRanks{s} = unitSelection{s};

        end
    end
    
    
end


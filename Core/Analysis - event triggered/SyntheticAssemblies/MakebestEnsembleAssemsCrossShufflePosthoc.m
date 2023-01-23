function Assemblies = MakebestEnsembleAssemsCrossShufflePosthoc(UnitsReal, Ltrials_,Rtrials_,MaxSize,CorrsToDestroy,noTrials,nReps,unitSelection)
    Assemblies= struct;

    %@TODO: add shuffling of trial labels between cells to non- trial-resampled version
    InfoCriteria = 'max';
    nReps = 50;
    noTrials = 10;
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

        noUnits = min([length(UnitsReal{3}),MaxSize]);

        Flag  = 0;
        if length(UnitsReal{3})==0
            Flag = 1;
        end

        CVE_   = cell(MaxSize,10); CVE_(:) = {nan(Ltr,1)};
        Score_ = nan(MaxSize,10);
        if noUnits >= 2 && ~Flag

            for iSize = 2:MaxSize
                fprintf('Ass Size = %d\n',iSize)

                parfor iRank = 1:10
                    unitIDs = unitSelection{3}{iSize,iRank};

                    if ~isempty(unitIDs)
                        [CVE_{iSize,iRank},Score_(iSize,iRank)] = RunDecodeGroupsDestroyNoiseCorrsCrossArea(unitIDs,nReps,Ltrials_,Rtrials_,noTrials,InfoCriteria,s);
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
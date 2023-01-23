function Assemblies = MakeAverageAssems(unitIDs,Ltrials_,Rtrials_,MaxSize)

Areas = {'PFC','HP','Joint'};
ResampleTrials = true;
InfoCriteria = 'max';
nReps = 10;   % Draw from trials
nDraws= 10;   % Draw from unit pool
noTrials = 20; % how many trials to use for each decoder
verbose = false;
Ltr = size(Ltrials_{1}{1},1);
noUnitsReal = [size(Ltrials_{1}{1},2),size(Ltrials_{2}{1},2),size(Ltrials_{3}{1},2)];

for s=1:3
    Assemblies.CVE{s}                 = nan(Ltr,MaxSize,nDraws);
    Assemblies.Score{s}               = nan(MaxSize,nDraws);
    Assemblies.Score_shuffled{s}      = nan(MaxSize,nDraws);
    Assemblies.Score_shuffled_5pc{s}  = nan(MaxSize,nDraws);
    Assemblies.Score_shuffled_95pc{s} = nan(MaxSize,nDraws);
    
    if ResampleTrials
        for i = 2:MaxSize
           
            Assemblies.CVE_shuffled{s}{i}        = nan(Ltr,2,nDraws);

            for iDraw = 1:nDraws
                 fprintf('Decoding Synthetic %s Assemblies: %d units pooled (%d maximum) Draw %d of %d...\n',Areas{s},i,MaxSize,iDraw,nDraws) 
                 if s<3
                    try
                        sample_ = randperm(length(unitIDs{s}));
                        unitIds = unitIDs{s}(sample_(1:i));
                    catch
                        fprintf('Ran out of units on size = %d\n',i);
                        break
                    end
                else
                    try
                        %Assumes equal weighting of assemblies between areas
                        %@TODO: Implement sliding ratios
                        sample_ = randperm(length(unitIDs{1}));
                        a  = unitIDs{1}(sample_(1:ceil(i/2)));
                        sample_ = randperm(length(unitIDs{2}));
                        b  = unitIDs{2}(sample_(1:ceil(i/2)));
                        unitIds = [a, b+noUnitsReal(1)];
                        unitIds=unitIds(randperm(length(unitIds)));
                        unitIds=unitIds(1:i);
                        
                    catch
                        fprintf('Ran out of units for %s size = %d',Areas{s},i);
                        break
                    end
                end
                [Assemblies.CVE{s}(:,i,iDraw),...
                    Assemblies.Score{s}(i,iDraw),...
                    Assemblies.Score_shuffled{s}(i,iDraw),...
                    Assemblies.Score_shuffled_5pc{s}(i,iDraw),...
                    Assemblies.Score_shuffled_95pc{s}(i,iDraw),...
                    Assemblies.CVE_shuffled{s}{i}(:,:,iDraw)] = RunDecodeGroups(unitIds,nReps,Ltrials_{s},Rtrials_{s},noTrials,'max');%
                
            end
            Assemblies.CVE_shuffled{s}{i}= nanmean(Assemblies.CVE_shuffled{s}{i},3);
        end
        
        Assemblies.CVE{s}   = squeeze(mean(Assemblies.CVE{s},3));
        Assemblies.Score{s} = nanmean(Assemblies.Score{s},2);
        Assemblies.Score_shuffled{s} = nanmean(Assemblies.Score_shuffled{s},2);
        Assemblies.Score_shuffled_5pc{s} = nanmean(Assemblies.Score_shuffled_5pc{s},2);
        Assemblies.Score_shuffled_95pc{s} = nanmean(Assemblies.Score_shuffled_95pc{s},2);
        
        figure; hold on
        for i =2:MaxSize
            plot(Assemblies.CVE{s}(:,i),'k')
            plot(Assemblies.CVE_shuffled{s}{i},':r')
            plot([1 Ltr],[Assemblies.Score_shuffled_5pc{s}(i) Assemblies.Score_shuffled_5pc{s}(i)],'r')
            plot([1 Ltr],[Assemblies.Score_shuffled_95pc{s}(i) Assemblies.Score_shuffled_95pc{s}(i)],'r')
        end
        
    else
        for i = 1:MaxSize
            fprintf('Decoding Synthetic %s Assemblies from joint area members: %d units pooled (%d maximum)...\n',Areas{s},i,MaxSize)
            
            if s<3
                unitIds = Assemblies.Rank{s}(1:i);
            else
                %Assumes equal weighting of assemblies between areas
                %@TODO: Implement sliding ratios
                a  = Assemblies.Rank{1}(1:ceil(i/2));
                a_ = Assemblies.SortedScore{1}(1:ceil(i/2));
                b  = Assemblies.Rank{2}(1:ceil(i/2));
                b_ = Assemblies.SortedScore{2}(1:ceil(i/2));
                if rem(i,2)~=0
                    if b_(end)>a_(end)
                        a(end)=[];
                    else
                        b(end)=[];
                    end
                end
                unitIds = [a, b+size(FR{1},2)];
            end
            CVE_ = 1-DecodeCVE_mex(FR{s}(:,unitIds),evt0,0.05);
            CVEshuf_ = 1-DecodeCVE_mex(FR{s}(:,unitIds),evt0(randperm(length(evt0))),0.05);
            Assemblies.CVE{s}(:,i)          = CVE_;
            Assemblies.Score{s}(i)          = max(CVE_);
            Assemblies.Score_shuffled{s}(i) = max(CVEshuf_);
        end
    end
end
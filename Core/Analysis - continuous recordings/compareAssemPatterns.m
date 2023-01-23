function P1vsP2 = compareAssemPatterns(P1,P2,plotYN,PermittedTransitions,RankOrder)
plotVerbose = false;
if nargin <3 || isempty(plotYN);
   plotYN = true; 
end
if nargin < 4 || isempty(PermittedTransitions)
    PermittedTransitions = 3;
end
if nargin <5 || isempty(RankOrder);
    RankOrder = true;
end

%% (1a) Change in assembly p(sequence) for global patterns
nGlobalEpoch = length(P1.FSC)+1;
A = P1.ptn.transMatrix{nGlobalEpoch};
B = P2.ptn.transMatrix{nGlobalEpoch};
switch PermittedTransitions
    case 1 % AND:  Ignore any pattern that is insignificant at any point
        A(A<P1.ptn.ciHTransition(nGlobalEpoch)) = NaN; 
        B(B<P2.ptn.ciHTransition(nGlobalEpoch)) = NaN;
    case 2 % OR:   Allow patterns through if they become/cease to be significant
        A(A<P1.ptn.ciHTransition(nGlobalEpoch) & B<P2.ptn.ciHTransition(nGlobalEpoch)) = NaN;
        B(A<P1.ptn.ciHTransition(nGlobalEpoch) & B<P2.ptn.ciHTransition(nGlobalEpoch)) = NaN;
    case 3 % Else: Include all transitions
end
        
temp = B - A; 

P1vsP2.ptn.delta_transMatrix = temp;
temp = temp(1:end);
if RankOrder
    [temp,idx] = sort(temp);
else
    idx = 1:length(temp);
end
labels_ = P1.ptn.labelsGlobal(idx);

if plotYN
    figure('color','w','name',[P1.filename ' to ' P2.filename  ': Change in global assembly activation patterns'])
    title(['Change in assembly sequences between ' P1.filename ' and ' P2.filename ' epochs'])
    hold on
    tempUP = temp;
    tempDOWN = temp;  
    idxUP = idx;
    idxDOWN = idx;
    tempUP(temp<0) = NaN; 
    tempDOWN(temp>0) = NaN; 

    stem(tempUP,'filled',':^b')    
    stem(tempDOWN,'filled',':vr')
    tempx = 1:sum(P1.nAss)^2;
    IDs = find(~isnan(tempUP));% | ~isnan(tempDOWN));
    for i = 1:length(IDs)
        text(tempx(IDs(i)),...
             tempUP(IDs(i)),...
             strcat({' '},labels_{IDs(i)})) %P1.ptn.labelsGlobal{IDs(i)}))
    end
    IDs = find(~isnan(tempDOWN));% | ~isnan(tempDOWN));
    for i = 1:length(IDs)
        text(tempx(IDs(i)),...
             tempDOWN(IDs(i)),...
             strcat({' '},labels_{IDs(i)})) %P1.ptn.labelsGlobal{IDs(i)}))
    end
%     axis([1 length(temp) -0.3 0.3])
%     axis([-Inf Inf -0.3 0.3])
    ylabel('\Delta p(Assembly pattern)');
    if isequal(PermittedTransitions,3)
        xlabel('Assembly activation patterns')
    else
        xlabel('Significant assembly activation patterns')
    end
end
clear i IDs tempx tempUP tempDOWN A B idx idxDOWN idxUP labels_ PermittedTransitions temp
%% (1b) Change in assembly sequence: DURATION and REPEAT TIMING of patterns
for s = 1:length(P1.titles) 
    A = P1.ptn.transMatrix{s};
    B = P2.ptn.transMatrix{s};

    % 1) Ignore any pattern that is insignificant at any point
    % A(A<P1.ptn.ciHTransition(nGlobalEpoch)) = NaN; B(B<P2.ptn.ciHTransition(nGlobalEpoch)) = NaN;

    % 2) Allow patterns through if the become/cease to be significant
    keepers = true(size(A));
    keepers(A<P1.ptn.ciHTransition(s) & B<P2.ptn.ciHTransition(s)) = false;

    P1vsP2.ptn.delta_DurMean{s} = P2.ptn.DurMean{s} - P1.ptn.DurMean{s};
    P1vsP2.ptn.delta_DurMean{s}(~keepers) = NaN;

    P1vsP2.ptn.delta_RepMean{s} = P2.ptn.RepMean{s} - P1.ptn.RepMean{s};
    P1vsP2.ptn.delta_RepMean{s}(~keepers) = NaN;
end
%% (2)  Change in repetitive assembly activation rate
for s = 1:length(P1.titles) 
    for iAss = 1:size(P1.pks.repHist{s},1)
        if plotYN && plotVerbose
        figure;
        hold on

        plot(P1.pks.repHistBins, smooth_hist(P1.pks.repHist{s}(iAss,:)),'b')
        plot(P2.pks.repHistBins, smooth_hist(P2.pks.repHist{s}(iAss,:)),'r')
        legend(P1.filename, P2.filename)
        scatter(P1.pks.repMode{s}(iAss),0,'b')
        scatter(P2.pks.repMode{s}(iAss),0,'r')
        title([P1.titles{s} ' assembly no. ' num2str(iAss)])
        set(gca,'xscale','log')
        xlabel('Activation interval (s)')
        ylabel('Distribution')
        end
        
        P1vsP2.pks.delta_repMode{s}(iAss) = P2.pks.repMode{s}(iAss) - P1.pks.repMode{s}(iAss);
        
    end
end
%% (3)  Change in assembly activation strength
for s = 1:length(P1.titles) 
    for iAss = 1:size(P1.pks.repHist{s},1)
        if plotYN && plotVerbose
            figure;
            hold on
            plot(P1.pks.LogStrengthHist{s}(:,iAss),P1.pks.LogBins,'b')
            plot(P2.pks.LogStrengthHist{s}(:,iAss),P2.pks.LogBins,'r')
            legend(P1.filename, P2.filename)
            scatter(P1.pks.LogStrengthQ3binValue{s}(iAss),P1.pks.LogStrengthQ3bin{s}(iAss),'b')
            scatter(P2.pks.LogStrengthQ3binValue{s}(iAss),P2.pks.LogStrengthQ3bin{s}(iAss),'r')

            scatter(P1.pks.ModeLogStrengthPos{s}(iAss),P1.pks.ModeLogStrength{s}(iAss),'b')
            scatter(P2.pks.ModeLogStrengthPos{s}(iAss),P2.pks.ModeLogStrength{s}(iAss),'r')
                title([P1.titles{s} ' assembly no. ' num2str(iAss)])
                ylabel('Activation strength')
                xlabel('Distribution')
        end
            
        % Change in activation: Mode and upper shoulder (3rd quartile)
        P1vsP2.pks.delta_ModeLogStrength{s}(iAss)  = P2.pks.ModeLogStrengthPos{s}(iAss) - P1.pks.ModeLogStrengthPos{s}(iAss);
        P1vsP2.pks.delta_LogStrengthQ3bin{s}(iAss) = P2.pks.LogStrengthQ3bin{s}(iAss)   - P1.pks.LogStrengthQ3bin{s}(iAss);
        % Change in activation: Upper shoulder (3rd quartile) - Mode 
        P1vsP2.pks.delta_Eccentric{s}(iAss) = (P2.pks.LogStrengthQ3bin{s}(iAss) - P2.pks.ModeLogStrengthPos{s}(iAss)) - (P1.pks.LogStrengthQ3bin{s}(iAss) - P1.pks.ModeLogStrengthPos{s}(iAss));
    end
end

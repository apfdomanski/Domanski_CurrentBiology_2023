function Ass = getCutoutAssemPatterns(Name_,bw,FSC,ciHld,ciHsc,FL,TrialOutcome,plotYN)
% Analyse patterns in assembly activation for repetitive short but
% variable-length delay periods
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
%
% 
% Loads data generated by runFAassem_[epoch] code (i.e. looks at assembly 
% activations that are calculated on a continous timebase, rather than
% cut-outs around behavioural task events).
%
% Wraps data into a structure called 'Ass' ...Cheeky!.
% Ass 
%  |-> Units    Which single units are members of which Assems, how do they overlap between them?
%  |-> Corr     Cross/Auto-correlation between factor activations in time       
%  |-> pks      Times/strengths of assembly activations, histograms of within-Assem repeat interval
%  |-> ptn      Assem activation patterns (1st order Markov process), timing, stats and Bootstrap CI margins
%
% Use to explore the timing and prevalence of assembly activation patterns that occur
% above chance

%% Preamble
if nargin<8
    flags.verbose = false;
else
    flags.verbose = plotYN;
end

Ass.filename  = Name_;
Ass.titles    = {'PFC','HP','Joint'};
Ass.bw        = bw;  % Time bin over which factor is calculated
per_          = 90;  % Percentile to choose for average strenth of assmbly activation 

%% First, analyse each (PFC/HP/Joint) factor type independently
for s = 1:length(FSC)  
%% (0) Basic Stuff                                           
    disp(['(0) Crunching... Basic Stuff (' Ass.titles{s} ')'])  %%% For each factor...
    Ass.nTrials              = numel(FSC{s});                    % Number of latent factors (taken from shuffled-neuron bootstrap)
    Ass.nTimePoints          = cellfun(@(x) size(x,1),FSC{s});
    Ass.nAss(s)              = size(FSC{s}{1},2);
    Ass.FSC{s}               = FSC{s};                          % Factor scores in time [Time x factors]

    Ass.thrsh.ciHld(s)       = ciHld(s,3);                      % Minimum factor loading bound to tag a unit as an assembly member
    Ass.thrsh.ciHsc(s)       = ciHsc(s,2);                      % Minimum factor score bound to tag timepoint as a factor activation
    Ass.LD{s}                = FL{s};                           % Factor loading matrix
    Ass.units.isMem{s}       = Ass.LD{s}>Ass.thrsh.ciHld(s);    % Which units Units tagged as members    
%% Core analysis
if Ass.nAss(s) ~= 0
%% (3) Factor score auto/cross-correlations
    disp(['(3) Crunching... Factor score auto/cross-correlations (' Ass.titles{s} ')'])
  
    for iTrial = 1:Ass.nTrials
        for AssID=1:Ass.nAss(s)
%             [ Ass.Corr.autoCorr{s}{AssID}(iTrial,:) , Ass.Corr.lags]    = xcorr(Ass.FSC{s}{iTrial}(:,AssID),1000,'coeff');
            for AssID_= 1:Ass.nAss(s)
                [Ass.Corr.crossCorr{s}{AssID,AssID_}(iTrial,:),...
                 Ass.Corr.lags] = xcorr(zscore(Ass.FSC{s}{iTrial}(:,AssID)),...
                                        zscore(Ass.FSC{s}{iTrial}(:,AssID_)),1000,'coeff'); %'biased'
            end
        end
    end
    
    % collate autocorrs
    Ass.Corr.autoCorr{s} = Ass.Corr.crossCorr{s}(find(eye(Ass.nAss(s))));
    
    % Trial-sorted averages
    for AssID=1:Ass.nAss(s)
        for  iOutcome=1:2
             Ass.Corr.autoCorr_sortedMean{s}{AssID}(:,iOutcome) = nanmean(Ass.Corr.autoCorr{s}{AssID}(TrialOutcome==iOutcome,:));
             Ass.Corr.autoCorr_sortedSEM{s}{AssID}(:,iOutcome)  = nansem (Ass.Corr.autoCorr{s}{AssID}(TrialOutcome==iOutcome,:));
              
        for AssID_= 1:Ass.nAss(s)
            Ass.Corr.crossCorr_sortedMean{s}{AssID,AssID_}(:,iOutcome) = nanmean(Ass.Corr.crossCorr{s}{AssID,AssID_}(TrialOutcome==iOutcome,:));
            Ass.Corr.crossCorr_sortedSEM{s}{AssID,AssID_}(:,iOutcome)  = nansem (Ass.Corr.crossCorr{s}{AssID,AssID_}(TrialOutcome==iOutcome,:));
        end
        end

    end
    
%%% Plot auto-corrs only
%     figure; c_ = {'r','b'};
%     for AssID=1:Ass.nAss(s)
%         for  iOutcome=1:2
% 
%         subplot(1,Ass.nAss(s),AssID); hold on
%         ciplot(Ass.Corr.autoCorr_sortedMean{s}{AssID}(:,iOutcome)+Ass.Corr.autoCorr_sortedSEM{s}{AssID}(:,iOutcome),...
%                Ass.Corr.autoCorr_sortedMean{s}{AssID}(:,iOutcome)-Ass.Corr.autoCorr_sortedSEM{s}{AssID}(:,iOutcome),...
%                Ass.Corr.lags*Ass.bw,c_{iOutcome})
%         plot(Ass.Corr.lags*Ass.bw,Ass.Corr.autoCorr_sortedMean{s}{AssID}(:,iOutcome),'color',c_{iOutcome})
%         
%         axis([-25 25 -0.1 1])
%         end
%     end
        if flags.verbose
        figure('color','w','NumberTitle','off','Name',[Ass.titles{s},' L/R Assembly sequence xCorr']); c_ = {'r','b'};
         for AssID=1:Ass.nAss(s)^2
            for  iOutcome=1:2
            if ismember(AssID,find(tril(ones(Ass.nAss(s),Ass.nAss(s)),0)))
            subplot(Ass.nAss(s),Ass.nAss(s),AssID); hold on
            ciplot(Ass.Corr.crossCorr_sortedMean{s}{AssID}(:,iOutcome)+Ass.Corr.crossCorr_sortedSEM{s}{AssID}(:,iOutcome),...
                   Ass.Corr.crossCorr_sortedMean{s}{AssID}(:,iOutcome)-Ass.Corr.crossCorr_sortedSEM{s}{AssID}(:,iOutcome),...
                   Ass.Corr.lags*Ass.bw,c_{iOutcome})
            plot(Ass.Corr.lags*Ass.bw,Ass.Corr.crossCorr_sortedMean{s}{AssID}(:,iOutcome),'color',c_{iOutcome})

            axis([-25 25 -0.1 0.2])
            if AssID == Ass.nAss
                %     legend('R trials','L trials'); legend('boxoff')
                text(8,0.15,'L trials','color','r')
                text(8,0.1,'R trials','color','b')
            end
            end
            end
         end
        end
%     
%        colors = [1 0 0;0 0 1];
%         for AssID=5%1:Ass.nAss(s)
%             for AssID_= 4%1:Ass.nAss(s)
%             figure('color','w'); hold on
%                 for iTrial = 1:Ass.nTrials
%                     y1=zscore(Ass.FSC{s}{iTrial}(:,AssID));
%                     y2=zscore(Ass.FSC{s}{iTrial}(:,AssID_));
%                     x = (1:length(y1))*bw;
%                     plot(x,y1 + 5*iTrial,'color',colors(TrialOutcome(iTrial),:))
%                     plot(x,y2 + 5*iTrial,'color',colors(TrialOutcome(iTrial),:))
%                     
%                 end
%             end
%         end
%         axis off
%         plot([30 35],[20 20],'LineWidth',1.5,'Color','k')
%         text(32.5,25,'5s','HorizontalAlignment','center')
%         
%         text(30,6,'L trials','color','r')
%         text(30,1,'R trials','color','b')
    
%% (4) Assembly activation times (factor score peaks: [amp, times]) and average strength
    disp(['(4) Crunching... Assembly activation times and strength (' Ass.titles{s} ')'])
    
    Ass.pks.LogBins=logspace(-1, 1.3,100);
    smooth_flag=0;
    meanshift_flag = 0;
    for AssID=1:Ass.nAss(s)
        for iTrial = 1:Ass.nTrials
            %%% Peak time/strength
            [Ass.pks.peak{s}{iTrial}{AssID},...
                Ass.pks.time{s}{iTrial}{AssID}] = findpeaks(Ass.FSC{s}{iTrial}(:, AssID),'MINPEAKHEIGHT',Ass.thrsh.ciHsc(s));
            %figure; hold on; plot((1:length(Ass.FSC{s}(:, AssID)))*Ass.bw,Ass.FSC{s}(:, AssID),'k'); scatter(Ass.pks.time{s}{AssID}*Ass.bw,Ass.pks.peak{s}{AssID},'xr')

            %%% Average strength (log transformed)
            temp = Ass.FSC{s}{iTrial}(:,AssID);
            if meanshift_flag
               temp = temp-median(temp);
            end
            temp=histc(temp,Ass.pks.LogBins);
            if smooth_flag, temp=smooth_hist(temp); end
            temp=temp./max(temp);
            temp(isnan(temp)) = 0;
            % FSC: mode (histogram peak)
            Ass.pks.LogStrengthHist{s}{AssID}(:,iTrial)=temp;
            Ass.pks.ModeLogStrength{s}(AssID,iTrial) = max(Ass.pks.LogBins(temp==max(temp)));
            Ass.pks.ModeLogStrengthPos{s}(AssID,iTrial)=max(temp);

            % FSC: 75% of cumulative distribution
            [~, per] = min(abs(cumsum(temp) - prctile(cumsum(temp),per_))); % find a percentile specified by per_
            Ass.pks.LogStrengthQ3bin{s}(AssID,iTrial)      = Ass.pks.LogBins(per);
            Ass.pks.LogStrengthQ3binValue{s}(AssID,iTrial) = temp(per);
            % figure;plot(Ass.pks.LogBins,temp,Ass.pks.LogStrengthQ3bin{s}(AssID),Ass.pks.LogStrengthQ3binValue{s}(AssID),'or'
        end  
    end
        
    % Outcome sorted averages
    for AssID=1:Ass.nAss(s)
        for  iOutcome=1:2
            Ass.pks.LogStrengthHist_sortedMean{s}{AssID}(:,iOutcome) = nanmean(Ass.pks.LogStrengthHist{s}{AssID}(:,TrialOutcome==iOutcome),2);
            Ass.pks.LogStrengthHist_sortedSEM{s}{AssID}(:,iOutcome)  = nansem(Ass.pks.LogStrengthHist{s}{AssID}(:,TrialOutcome==iOutcome),2);
            
            Ass.pks.ModeLogStrengthPos_sortedMean{s}(:,iOutcome)  = nanmean(Ass.pks.ModeLogStrengthPos{s}(:,TrialOutcome==iOutcome),2);
            Ass.pks.ModeLogStrengthPos_sortedSEM{s}(:,iOutcome)   = nansem(Ass.pks.ModeLogStrengthPos{s}(:,TrialOutcome==iOutcome),2);
            
            Ass.pks.LogStrengthQ3bin_sortedMean{s}(:,iOutcome)  = nanmean(Ass.pks.LogStrengthQ3bin{s}(:,TrialOutcome==iOutcome),2);
            Ass.pks.LogStrengthQ3bin_sortedSEM{s}(:,iOutcome)   = nansem(Ass.pks.LogStrengthQ3bin{s}(:,TrialOutcome==iOutcome),2);
        end
    end
    
    if flags.verbose
    figure('color','w'); c_ = {'r','b'};
    for AssID=1:Ass.nAss(s)
        subplot(1,Ass.nAss(s),AssID); hold on
        for  iOutcome=1:2
            plot(Ass.pks.LogStrengthHist_sortedMean{s}{AssID}(:,iOutcome),Ass.pks.LogBins,'color',c_{iOutcome})
            plot(0,Ass.pks.LogStrengthQ3bin_sortedMean{s}(AssID,iOutcome),'o','color',c_{iOutcome})
        end
             
    end
    end
%% (5) factor score inter-peak Intervals
    disp(['(5) Crunching... Factor score inter-peak Intervals (' Ass.titles{s} ')'])
    Ass.pks.repHistBins=0:20:1000;
    Ass.pks.repHist{s}=cell(Ass.nAss(s),1);
    Ass.pks.repHist_sortedMean{s}=cell(Ass.nAss(s),1);
    Ass.pks.repHist_sortedSEM{s}=cell(Ass.nAss(s),1);
    for AssID=1:Ass.nAss(s)
        for iTrial = 1:Ass.nTrials
%             [s AssID iTrial]
            if length(Ass.pks.time{s}{iTrial}{AssID})>1
               Ass.pks.repHist{s}{AssID}(iTrial,:) = histc (diff(Ass.pks.time{s}{iTrial}{AssID}),Ass.pks.repHistBins)./(numel(Ass.pks.time{s}{iTrial}{AssID})-1);
               Ass.pks.repMean{s}(AssID,iTrial)    = mean  (diff(Ass.pks.time{s}{iTrial}{AssID}));
               Ass.pks.repMedian{s}(AssID,iTrial)  = median(diff(Ass.pks.time{s}{iTrial}{AssID}));
               Ass.pks.repMode{s}(AssID,iTrial)    = mode  (diff(Ass.pks.time{s}{iTrial}{AssID}));
               Ass.pks.repSEM{s}(AssID,iTrial)     = nansem(diff(Ass.pks.time{s}{iTrial}{AssID}));
               Ass.pks.repSTD{s}(AssID,iTrial)     = nanstd(diff(Ass.pks.time{s}{iTrial}{AssID}));
            else
               Ass.pks.repHist{s}{AssID}(iTrial,:) = NaN*(Ass.pks.repHistBins);
               Ass.pks.repMean{s}(AssID,iTrial)    = NaN;
               Ass.pks.repMedian{s}(AssID,iTrial)  = NaN;
               Ass.pks.repMode{s}(AssID,iTrial)    = NaN;
               Ass.pks.repSEM{s}(AssID,iTrial)     = NaN;
               Ass.pks.repSTD{s}(AssID,iTrial)     = NaN;          
            end
            % plot((0:1000)*Ass.bw*1000,Ass.pks.IEI_hist{s}'); set(gca,'Xscale','log')
        end
        for iOutcome = 1:2
             Ass.pks.repHist_sortedMean{s}{AssID}(:,iOutcome)  = nanmean(Ass.pks.repHist{s}{AssID}(TrialOutcome==iOutcome,:));
             Ass.pks.repHist_sortedSEM {s}{AssID}(:,iOutcome)  = nansem (Ass.pks.repHist{s}{AssID}(TrialOutcome==iOutcome,:));
             Ass.pks.repMean_sortedMean{s}(AssID,iOutcome)         = nanmean(Ass.pks.repMean{s}(AssID,TrialOutcome==iOutcome),2);
             
        end
    end
    if flags.verbose
    figure('color','w'); c_ = {'r','b'};
    for AssID=1:Ass.nAss(s)
        subplot(1,Ass.nAss(s),AssID); hold on
        for  iOutcome=1:2
            plot(Ass.pks.repHistBins*Ass.bw,smooth_hist(Ass.pks.repHist_sortedMean{s}{AssID}(:,iOutcome)),'color',c_{iOutcome})
            plot(Ass.pks.repMean_sortedMean{s}(AssID,iOutcome)*Ass.bw,0,'o','color',c_{iOutcome})
        end
        axis([0 25 0 Inf])
        set(gca,'XScale','log')
    end
    end
%% (6) factor score sequences
    disp(['(6) Crunching... Factor score sequences (' Ass.titles{s} ')'])
    for iTrial = 1:Ass.nTrials
        Ass.ptn.seq{s}{iTrial} = []; % This grows with each turn loop iteration (i.e self-concatenating)
        for AssID=1:Ass.nAss(s)
            Ass.ptn.seq{s}{iTrial} = [Ass.ptn.seq{s}{iTrial}; [Ass.pks.time{s}{iTrial}{AssID},AssID+zeros(size(Ass.pks.time{s}{iTrial}{AssID}))]];
        end
        Ass.ptn.seq{s}{iTrial} = sortrows(Ass.ptn.seq{s}{iTrial}); % rearrange assembly pattern by peak time      end %/ if no assemblies
    end
%% (7) Estimate each Assem->Assem transition probability (~ 1st order Markov process)
    disp(['(7) Crunching... Estimating Assem->Assem transition probability (' Ass.titles{s} ')'])
    Ass.AssLabels{s}=cellfun(@num2str,num2cell(1:Ass.nAss(s)),'UniformOutput',0);
    for iTrial = 1:Ass.nTrials
            seq_temp=Ass.ptn.seq{s}{iTrial}(:,2);

        % Catch the case where a 1>max pattern doesn't exist
        for i = 1:Ass.nAss(s)
            for j = 1:Ass.nAss(s)
                seq_temp = [seq_temp ; i ;j];
            end
        end
        % TODO - fix this ASAP!
        % TODO - higher order Markov processes?
        % [ transitionMatrix columnStates ]=getTransitionMatrix(seq_temp,3);
        % sortrows(transitionMatrix )
        Ass.ptn.transMatrix{s}{iTrial}= zeros(Ass.nAss(s),Ass.nAss(s));
        [Ass.ptn.transMatrix{s}{iTrial},~]=hmmestimate(seq_temp,seq_temp);       
    end
    % Trial average
    for iOutcome=1:2
        temp =  Ass.ptn.transMatrix{s}(TrialOutcome==iOutcome);
        Ass.ptn.transMatrix_sortedMean{s}{iOutcome}=zeros(size(temp{1}));
        for iTrial=1:numel(temp)
            Ass.ptn.transMatrix_sortedMean{s}{iOutcome} = Ass.ptn.transMatrix_sortedMean{s}{iOutcome}+temp{iTrial}; 
        end
        Ass.ptn.transMatrix_sortedMean{s}{iOutcome}=Ass.ptn.transMatrix_sortedMean{s}{iOutcome}./numel(temp);
    end
    
    if flags.verbose & Ass.nAss(s)>1
    figure('color','w'); c_ = {'r','b'}; hold on
    for  iOutcome=1:2
        plot( Ass.ptn.transMatrix_sortedMean{s}{iOutcome}(1:end),'color',c_{iOutcome})
    end
    
    figure('color','w'); c_ = {'r','b'};         titles_ = {'L trials','R trials'};  

    for  iOutcome=1:2
        subplot(2,1,iOutcome)
        title(titles_{iOutcome})
        pcolor(Ass.ptn.transMatrix_sortedMean{s}{iOutcome})
        xlabel('Assembly n')
        ylabel('Assembly n+1')
        
    end
%      colormap parula
%     c = colorbar('Location','manual','Position',[0.75 0.2 0.05 0.2]);
%     c.Label.String = {'Sum decoding power';'of assembly sequence'};

    end
end
end
%% Then, concatenate times from all areas for global transition probability
 % Make Assebbly transiton labels
    Ass.AssLabels{4} = cell(Ass.nTrials,1);
    Ass.labelsGlobal= {};
    for s = 1: numel(Ass.nAss)
        for iAss = 1:Ass.nAss(s)
            Ass.labelsGlobal = [Ass.labelsGlobal;[Ass.titles{s}(1) num2str(iAss)]];
        end
    end
    
    for AssID=1:sum(Ass.nAss)
        for AssID_= 1:sum(Ass.nAss)
            Ass.ptn.labelsGlobal{AssID,AssID_}= [Ass.labelsGlobal{AssID} ' > ' Ass.labelsGlobal{AssID_}];  % record a pattern label
        end
    end
    
    
    disp('(10) Crunching... Assembly Markov model (global)')
    Ass.AssLabels{4} = [];
    for iTrial = 1:Ass.nTrials
        temp = [];    
        for s = 1: numel(Ass.nAss)
            try
                temp{s} =  Ass.ptn.seq{s}{iTrial};
            catch 
                temp{s} = [];
            end
        end
        temp = concatenatePatterns(temp);
        temp = sortrows(temp,1);
        Ass.ptn.seq_global{iTrial} = temp;
        Ass.AssLabels{4}=[ Ass.AssLabels{4} ; unique(temp(:,2))];
            
        seq_temp=temp(:,2);
        % TODO - fix this ASAP!
        % Catch the case where a 1>max pattern doesn't exist
        seq_temp = [seq_temp ; 1 ;max(unique(seq_temp))];
        
        % Catch the case where a 1>max pattern doesn't exist
        for i = 1:sum(Ass.nAss)
            for j = 1:sum(Ass.nAss)
                seq_temp = [seq_temp ; i ;j];
            end
        end
        % TODO - fix this ASAP!
        % [ transitionMatrix columnStates ]=getTransitionMatrix(seq_temp,3);
        [Ass.ptn.transMatrix{4}{iTrial},~]=hmmestimate(seq_temp,seq_temp);    
    end
    Ass.AssLabels{4}=unique( Ass.AssLabels{4});

    
     % Trial average
    for iOutcome=1:2
        temp = Ass.ptn.transMatrix{4}(TrialOutcome==iOutcome);
        Ass.ptn.transMatrix_sortedMean{4}{iOutcome}=zeros(size(temp{1}));
        for iTrial=1:numel(temp)
            Ass.ptn.transMatrix_sortedMean{4}{iOutcome} = Ass.ptn.transMatrix_sortedMean{4}{iOutcome}+temp{iTrial}; 
        end
        Ass.ptn.transMatrix_sortedMean{4}{iOutcome}=Ass.ptn.transMatrix_sortedMean{4}{iOutcome}./numel(temp);
    end
    if flags.verbose
    figure('color','w'); c_ = {'r','b'}; hold on
    for  iOutcome=1:2
        plot( Ass.ptn.transMatrix_sortedMean{4}{iOutcome}(1:end),'color',c_{iOutcome})
    end
    
    figure('color','w'); c_ = {'r','b'};
    for  iOutcome=1:2
        subplot(2,1,iOutcome)
        pcolor(Ass.ptn.transMatrix_sortedMean{4}{iOutcome})
    end
    end
    
end
function globalPattern = concatenatePatterns(pattern)
    globalPattern = [pattern{1}];
    for s = 2:length(pattern)
        try 
            globalPattern = [globalPattern; ...
                             pattern{s}(:,1),pattern{s}(:,2)+max(globalPattern(:,2))]; 
        catch 
            globalPattern = pattern{s};
        end
    end
end
function start = findPattern2(array, pattern)
% Fast lookup: Locate a pattern in an array.
%
%   indices = findPattern2(array, pattern) finds the starting indices of
%   pattern within array.
%
%   Example:
%   a = [0 1 4 9 16 4 9];
%   patt = [4 9];
%   indices = findPattern2(a,patt)
%   indices =
%        3     6

len = length(pattern);
start = find(array==pattern(1));
endVals = start+len-1;
start(endVals>length(array)) = [];
for pattval = 2:len
    locs = pattern(pattval) == array(start+pattval-1);
    start(~locs) = [];
end
end
function [ transitionMatrix, columnStates ] = getTransitionMatrix(markovChain,Norder)
%  Constructs the transition matrix of a markov chain, given the markovChain and the order, i.e., an integer greater or equal to 1.
%
% Inputs:
%    markovChain : the markov chain, in integers.
%    Norder : the order to be analyzed. 
% Ouptuts: 
%    transition matrix: the state-transition matrix (TM), where each value represents
%                       the number of occurrence for a sequence of states, which is the 
%                       previous state (column of TM) followed by the current state (row of TM).
%                       (See references for more info.)
%    columnStates: the sequence for the previous state(column of TM).
%%
% %Example 1:
% %outputs the 1st order transition matrix of the below markov chain.
% markovChain = [1 6 1 6 4 4 4 3 1 2 2 3 4 5 4 5 2 6 2 6 2]; Norder = 1;
% [ transitionMatrix columnStates ] = getTransitionMatrix(markovChain,Norder);
%
%%
% %Example 2:
% %plots the 2nd order transition matrix of a random markov chain with 4 states
% getTransitionMatrix;
%
%%
% ref: http://en.wikipedia.org/wiki/Markov_chain
%      http://stackoverflow.com/questions/11072206/constructing-a-multi-order-markov-chain-transition-matrix-in-matlab
%      
% $ version 1   $ by Pangyu Teng $ 25JUN2012 $ created $ 
% $ version 1.1 $ by Pangyu Teng $ 26Jun2012 $ updated, multi-order support
%

if nargin < 2,
    display(sprintf('2 inputs required (getTransitionMatrix.m)\nusing example data... '));
    %return;        
    markovChain = round(3*rand([1,1000])+1);
    Norder = 2;
end

if Norder < 1,
    display('Norder has to be >1 (getTransitionMatrix.m)');
    return;
end

if numel(markovChain) <= Norder
    display('Need more data for the 1st input. (getTransitionMatrix.m)');
    return;
end

%make markovChain a column
if size(markovChain,1) > 1;
    markovChain = markovChain';
end

%number of states
Nstates = max(markovChain);

if Norder == 1,
    
    %get transition matrix
    transitionMatrix = full(sparse(markovChain(1:end-1),markovChain(2:end),1,Nstates^Norder,Nstates));
    columnStates = cellstr(regexprep(num2str(1:Nstates),'[^\w'']','')');
    
else

    %get Norder-contiguous sequences of the markov chain
    ngrams = [];
    for i = 0:Norder-1
        ngrams = [ngrams, circshift(markovChain,[0 -1*(i)])'];
    end
    ngrams = cellstr(num2str( ngrams));
    ngrams = ngrams(1:end-Norder);

    
    %create  all combinations of Norder-contiguous sequences
    [x{1:Norder}] = ndgrid(1:Nstates);

    %format x to cell
    evalStr = ['xCell = cellstr(num2str(['];
    for i = 1:Norder
        evalStr = [evalStr 'x{' num2str(i) '}(:) '];            
    end
    evalStr = [evalStr ']));'];
    eval(evalStr);

    %map ngrams to numbers
    [gn,trashToo,g]=unique([xCell;ngrams]);
    s1 = g(Nstates^Norder+1:end);

    %states following the ngrams
    s2 = markovChain(Norder+1:end);

    %reordered xCell.
    columnStates = gn(1:Nstates^Norder);    
    
    %get transition matrix
    transitionMatrix = full( sparse(s1,s2,1,Nstates^Norder,Nstates) );

end

if nargout < 1

    %plot the transition matrix
    imagesc(transitionMatrix);
    str= 'Transition Matrix of the Markov Chain:\n';
    str=[str sprintf('%d',markovChain) '...'];
    title(sprintf(str));    
    
    %replace tickLabels with states.
    set(gca,'YTick',1:numel(columnStates));
    set(gca,'YTickLabel',columnStates);
    set(gca,'XTick',1:Nstates);
    set(gca,'XTickLabel',1:Nstates);
    
end
end
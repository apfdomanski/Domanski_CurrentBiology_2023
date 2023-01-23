function Ass = getAssemPatterns(Name_,bw,FSC,ciHld,ciHsc,FL)
% Analyse patterns in assembly activation
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
flags.verbose = false;
Ass.filename  = Name_;
Ass.titles    = {'PFC','HP','Joint'};
Ass.bw        = bw;  % Time bin over which factor is calculated
per_          = 90;  % Percentile to choose for average strenth of assmbly activation 

%% First, analyse each (PFC/HP/Joint) factor type independently
for s = 1:length(FSC)    
    %% (0) Basic Stuff                                           
    disp(['(0) Crunching... Basic Stuff (' Ass.titles{s} ')'])  %%% For each factor...
    [Ass.nTimePoints,...                                        % Number of timepoints
    Ass.nAss(s)]             = size(FSC{s});                    % Number of latent factors (taken from shuffled-neuron bootstrap)
    Ass.FSC{s}               = FSC{s};                          % Factor scores in time [Time x factors]
    Ass.thrsh.ciHld(s)       = ciHld(s,3);                      % Minimum factor loading bound to tag a unit as an assembly member
    Ass.thrsh.ciHsc(s)       = ciHsc(s,2);                      % Minimum factor score bound to tag timepoint as a factor activation
    Ass.LD{s}                = FL{s};                           % Factor loading matrix
    Ass.units.isMem{s}       = Ass.LD{s}>Ass.thrsh.ciHld(s);    % Which units Units tagged as members
    
    if Ass.nAss(s) ~= 0
    %% (1) Promiscuity of units between factors
    disp(['(1) Crunching... Promiscuity of units between factors (' Ass.titles{s} ')'])
    Ass.units.AssOverlap{s}          = sum(Ass.units.isMem{s},2);                           % No. factors each unit is shared between
    Ass.units.FractionSingleAssem(s) = sum(Ass.units.AssOverlap{s}==1)./size(Ass.LD{s},1);  % Fraction of neurons assigned to only one Assem
    Ass.units.AssOverlapHistBins   = 1:10;
    % Ass.units.AssOverlap{s}          = Ass.units.AssOverlap{s}./ Ass.nAss(s);             % Adjust to fraction of detected factors
    % Ass.units.AssOverlapHistBins   = 0:0.1:1;
    Ass.units.AssOverlapHist(:,s)   = histc(Ass.units.AssOverlap{s},Ass.units.AssOverlapHistBins)...
                                       ./sum(Ass.units.AssOverlap{s}>0);                    % distribution of factor overlaps by area, normalised by no. units involved in factors
    %% (2) Which units are members of this factor?
    disp(['(2) Crunching... member units (' Ass.titles{s} ')'])
    for AssID=1:Ass.nAss(s)
        Ass.units.memIDs{s}{AssID}       = find(Ass.units.isMem{s}(:,AssID));
    end 
    %% (3) Factor score auto/cross-correlations
    disp(['(3) Crunching... Factor score auto/cross-correlations (' Ass.titles{s} ')'])
    for AssID=1:Ass.nAss(s)
        [ Ass.Corr.autoCorr{s}(AssID,:) , Ass.Corr.lags]    = xcorr(Ass.FSC{s}(:,AssID),1000,'coeff');
        for AssID_= 1:Ass.nAss(s)
             Ass.Corr.crossCorr{s}{AssID,AssID_}             = xcorr(Ass.FSC{s}(:,AssID), Ass.FSC{s}(:,AssID_),1000,'biased');
        end
    end   
    %% (4) Assembly activation times (factor score peaks: [amp, times]) and average strength
    disp(['(4) Crunching... Assembly activation times and strength (' Ass.titles{s} ')'])
    
    Ass.pks.LogBins=logspace(-1, 1.3,100);
    smooth_flag=0;

    for AssID=1:Ass.nAss(s)
        %%% Peak time/strength
        [Ass.pks.peak{s}{AssID},...
            Ass.pks.time{s}{AssID}] = findpeaks(Ass.FSC{s}(:, AssID),'MINPEAKHEIGHT',Ass.thrsh.ciHsc(s));
        %figure; hold on; plot((1:length(Ass.FSC{s}(:, AssID)))*Ass.bw,Ass.FSC{s}(:, AssID),'k'); scatter(Ass.pks.time{s}{AssID}*Ass.bw,Ass.pks.peak{s}{AssID},'xr')
   
        %%% Average strength (log transformed)
        temp = Ass.FSC{s}(:,AssID);
        temp=histc(temp-median(temp),Ass.pks.LogBins);
        if smooth_flag, temp=smooth_hist(temp); end
        temp=temp./max(temp);
            
        % FSC: mode (histogram peak)
        Ass.pks.LogStrengthHist{s}(:,AssID)=temp;
        Ass.pks.ModeLogStrength{s}(AssID) = max(Ass.pks.LogBins(temp==max(temp)));
        Ass.pks.ModeLogStrengthPos{s}(AssID)=max(temp);
        
        % FSC: 75% of cumulative distribution
        [~, per] = min(abs(cumsum(temp) - prctile(cumsum(temp),per_))); % find a percentile specified by per_
        Ass.pks.LogStrengthQ3bin{s}(AssID)      = Ass.pks.LogBins(per);
        Ass.pks.LogStrengthQ3binValue{s}(AssID) = temp(per);
        % figure;plot(Ass.pks.LogBins,temp,Ass.pks.LogStrengthQ3bin{s}(AssID),Ass.pks.LogStrengthQ3binValue{s}(AssID),'or'
        
    end    
    %% (5) factor score inter-peak Intervals
    disp(['(5) Crunching... Factor score inter-peak Intervals (' Ass.titles{s} ')'])
    for AssID=1:Ass.nAss(s)
        %[s AssID]
       Ass.pks.repHistBins=0:1:1000;
       
        if length(Ass.pks.time{s}{AssID})>1
           Ass.pks.repHist{s}(AssID,:) = histc (diff(Ass.pks.time{s}{AssID}),Ass.pks.repHistBins)./(numel(Ass.pks.time{s}{AssID})-1);
           Ass.pks.repMean{s}(AssID)   = mean  (diff(Ass.pks.time{s}{AssID}));
           Ass.pks.repMedian{s}(AssID) = median(diff(Ass.pks.time{s}{AssID}));
           Ass.pks.repMode{s}(AssID)   = mode  (diff(Ass.pks.time{s}{AssID}));
           Ass.pks.repSEM{s}(AssID)    = nansem(diff(Ass.pks.time{s}{AssID}));
           Ass.pks.repSTD{s}(AssID)    = nanstd(diff(Ass.pks.time{s}{AssID}));
        else
           Ass.pks.repHist{s}(AssID,:)   = NaN*(Ass.pks.repHistBins);
           Ass.pks.repMean{s}(AssID)   = NaN;
           Ass.pks.repMedian{s}(AssID) = NaN;
           Ass.pks.repMode{s}(AssID)   = NaN;
           Ass.pks.repSEM{s}(AssID)    = NaN;
           Ass.pks.repSTD{s}(AssID)    = NaN;          
        end
        
        % plot((0:1000)*Ass.bw*1000,Ass.pks.IEI_hist{s}'); set(gca,'Xscale','log')
    end
    %% (6) factor score sequences
    disp(['(6) Crunching... Factor score sequences (' Ass.titles{s} ')'])
    Ass.ptn.seq{s} = []; % This grows with each turn loop iteration (i.e self-concatenating)
    for AssID=1:Ass.nAss(s)
        Ass.ptn.seq{s} = [Ass.ptn.seq{s}; [Ass.pks.time{s}{AssID},AssID+zeros(size(Ass.pks.time{s}{AssID}))]];
    end
    Ass.ptn.seq{s} = sortrows(Ass.ptn.seq{s}); % rearrange assembly pattern by peak time   
    %% (7) Estimate each Assem->Assem transition probability (~ 1st order Markov process)
    disp(['(7) Crunching... Estimating Assem->Assem transition probability (' Ass.titles{s} ')'])
    Ass.AssLabels{s}=cellfun(@num2str,num2cell(1:Ass.nAss(s)),'UniformOutput',0);
    seq_temp=Ass.ptn.seq{s}(:,2);
    
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
    Ass.ptn.transMatrix{s}= zeros(Ass.nAss(s),Ass.nAss(s));
    [Ass.ptn.transMatrix{s},~]=hmmestimate(seq_temp,seq_temp);    
    %% (8) Bootstrap tranistion matrices on shuffled sequences
    disp(['(8) Crunching... Bootstrap tranistion matrices on shuffled sequences (' Ass.titles{s} ')'])
    alpha = 0.05; % Confidence limit for bootstrap transition probability threshold
    BS_max=500;transMatrixBS=[];
    for BSid=1:BS_max
        if flags.verbose
            disp(['Shuffling... Pattern no. ' num2str(BSid) ' of ' num2str(BS_max)])
        end
        BSperm=seq_temp(randperm(length(seq_temp)));
        [transMatrixBS{BSid},~]=hmmestimate(BSperm,BSperm);
    end
    transMatrixBS=cell2mat(transMatrixBS);
    ws=sort(transMatrixBS(1:end),'ascend');
    Ass.ptn.ciLTransition(s)=ws(max(1,round(length(ws)*alpha)));
    Ass.ptn.ciHTransition(s)=ws(max(1,round(length(ws)*(1-alpha))));    
    %% (9) Pattern searching (duration, delay)
    disp(['(9) Crunching... Pattern searching (' Ass.titles{s} ')'])
    Ass.ptn.DurBins = 1:1:100;  Ass.ptn.RepBins = 1:50:1000; % These bins are in units of Ass.bw, NOT time
    for AssID=1:Ass.nAss(s);
        for AssID_= 1:Ass.nAss(s)

            Ass.ptn.labels{s}{AssID,AssID_}= [num2str(AssID) ' > ' num2str(AssID_)];  % record a pattern label

            % Use exact template search to find pattern
            if flags.verbose
                disp(['Searching for pattern ' num2str(AssID) ' -> ' num2str(AssID_)])
            end
            PtnTimes = findPattern2(Ass.ptn.seq{s}(:,2),[AssID,AssID_]);
            
            if isempty(PtnTimes);   % guard against any non-existent patterns...              
                ptnDur   = NaN;    
                ptnRep   = NaN;    
            else
                ptnDur   = zeros(length(PtnTimes),1);     % Intra-pattern pattern duration
                ptnRep   = zeros(length(PtnTimes)-1,1);   % Inter-pattern pattern repeat

                for Rep_id=1:length(PtnTimes) % For each repeat of this pattern...
                    % Pattern duration
                    ptnDur(Rep_id) = Ass.ptn.seq{s}(PtnTimes(Rep_id)+1,1) - ...
                                     Ass.ptn.seq{s}(PtnTimes(Rep_id),1);
                
                    if Rep_id>1  
                        % Pattern Repeats 
                        ptnRep(Rep_id-1) = Ass.ptn.seq{s}(PtnTimes(Rep_id),1) - Ass.ptn.seq{s}(PtnTimes(Rep_id-1),1);
                    end
                end
%                 ptnRep(ptnRep==0)=[];
            end        
            if isempty(ptnRep) ; ptnRep = NaN; end % guard against any one-off patterns
            
            % How many repeats of this pattern
            Ass.ptn.RepCount{s}(AssID,AssID_) = length(PtnTimes);

            % average duration
            Ass.ptn.DurMedian{s}(AssID,AssID_) = nanmedian(ptnDur);    
            Ass.ptn.DurMean{s}(AssID,AssID_)   = nanmean(ptnDur);      
            Ass.ptn.DurSEM{s}(AssID,AssID_)    = nansem(ptnDur);       
            
            % average repeat interval
            Ass.ptn.RepMedian{s}(AssID,AssID_) = nanmedian(ptnRep);    
            Ass.ptn.RepMean{s}(AssID,AssID_)   = nanmean(ptnRep);     
            Ass.ptn.RepSEM{s}(AssID,AssID_)    = nansem(ptnRep);       
            
            % Histograms
            Ass.ptn.DurHist{s}{AssID,AssID_}  = histc(ptnDur,Ass.ptn.DurBins)./numel(ptnDur);
            Ass.ptn.RepHist{s}{AssID,AssID_}  = histc(ptnRep,Ass.ptn.RepBins)./numel(ptnRep);
        end
    end
    end

end % /if no assemblies
%% Then, concatenate times from all areas for global transition probability
if length(FSC)>1
    
    % Make Assebbly transiton labels
    Ass.labelsGlobal= {};
    for s = 1: numel(Ass.nAss)
        for iAss = 1:Ass.nAss(s)
            Ass.labelsGlobal = [Ass.labelsGlobal;[Ass.titles{s}(1) num2str(iAss)]];
        end
    end
    
    for AssID=1:sum(Ass.nAss);
        for AssID_= 1:sum(Ass.nAss)
            Ass.ptn.labelsGlobal{AssID,AssID_}= [Ass.labelsGlobal{AssID} ' > ' Ass.labelsGlobal{AssID_}];  % record a pattern label
        end
    end
    
    disp('(10) Crunching... Assembly Markov model (global)')
    temp = [];
    temp = concatenatePatterns(Ass.ptn.seq);
    temp = sortrows(temp,1);
    
    Ass.ptn.seq_global = temp;
    Ass.AssLabels{4}=unique(temp(:,2));
    
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
    [Ass.ptn.transMatrix{4},~]=hmmestimate(seq_temp,seq_temp);


    % (8) Bootstrap tranistion matrices on shuffled sequences
    disp(['(8) Crunching... Bootstrap transition matrices on shuffled sequences (Global)'])
    alpha = 0.05; % Confidence limit for bootstrap transition probability threshold
    BS_max=500;transMatrixBS=[];
    for BSid=1:BS_max
        if flags.verbose
        disp(['Shuffling... Pattern no. ' num2str(BSid) ' of ' num2str(BS_max)])
        end
        BSperm=seq_temp(randperm(length(seq_temp)));
        [transMatrixBS{BSid},~]=hmmestimate(BSperm,BSperm);
    end
    transMatrixBS=cell2mat(transMatrixBS);
    ws=sort(transMatrixBS(1:end),'ascend');
    Ass.ptn.ciLTransition(4)=ws(max(1,round(length(ws)*alpha)));
    Ass.ptn.ciHTransition(4)=ws(max(1,round(length(ws)*(1-alpha))));
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
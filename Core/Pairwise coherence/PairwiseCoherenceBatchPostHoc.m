%% Preamble, get file list
clear all; close all
set(0,'DefaultFigureWindowStyle','docked')
flags.plotIndividualChanges = false;
% pat{1} -  % location of the processed assembly memberships
% pat{2} -  % location of raw spike times
if ispc
    %     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - epoch specific FRs\';
    pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - use task mean FRs\';
    pat{2} = 'C:\Analysis\AssemblyAnalysis\raw';          
elseif ismac
    pat{1} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Decoding vs sleep/';
    pat{2} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Task';
else
    %    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    addpath(genpath('/panfs/panasas01/phph/ad15419/MATLAB/SpikeTrainAna'))
    p.home = getenv('HOME');
    cd([p.home '/MATLAB/SpikeTrainAna/'])
    pat{1} = [p.home '/Sleep/DecodingVsSleep/'];
    pat{2} = [p.home '/raw/'];
end
cd(pat{1})
RatList = dir([pat{1} filesep '*_coherence.mat']);

% Uncomment to ignore some specific files...

% Ignore= {'JaroslawLONG1'; ...
%          'MiroslawLONG2'; ...
%          'NorbertLONG2' ; ...
%          'OnufryLONG2'};
% for iList = 1:length(Ignore)
%     RatList = RatList(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
%            {RatList.name},repmat({Ignore(iList)},1,length(RatList)),'UniformOutput',false)));
% end

%% Loop over recordings
for iRat = 1:length(RatList)
    tic

    disp([sprintf('Working on rat %d of %d',iRat, length(RatList)), ' (' RatList(iRat).name ')'])
    load([pat{1} RatList(iRat).name]);
    
    load([pat{1} strtok(RatList(iRat).name,'_') '_Task_DecodingVSsleep.mat'],'P','D','FRtrials','FAtrials');
    load([pat{2} filesep strtok(RatList(iRat).name,'_')],'HPcells','PFCcells')
   
    %% Batch process
    rateThresh = 0; % Minimal joint rate threshold to include cell pair
    
    spikes = [PFCcells(FAtrials.unitIDs{1});...
              HPcells(FAtrials.unitIDs{2})];
          
    % joint firing rate of unit pairs in the task period
    for iUnit = 1:length(spikes)
        temp = spikes{iUnit}.t*1e-4;
        temp(temp<min(FRtrials.EvtTs) | temp>max(FRtrials.EvtTs)) =[];
        avgFR(iUnit)=length(temp)/range(temp);
    end
    for i = 1:length(spikes)
        for j = 1:length(spikes)
            meanRate(i,j) =  sqrt(avgFR(i)*avgFR(j));
        end
    end
    
    len_ = max(unique(cellfun(@length,t_(1:end))));
    tempMembers   = cell(size(coherence)); % index array to follow which in-assembly pairs are included
    tempCrossMembers = cell(size(coherence)); % index array to follow which out-of assembly pairs are included
    SameAssCoherence = cell(length(FAtrials.units{3}),1);     % Same assembly
    DiffAssCoherence = cell(length(FAtrials.units{3}),1);     % Different assembly
    AllAssMembers    = unique(cell2mat(FAtrials.units{3}));   % All member assembly units
    noUnits          = cellfun(@length,FAtrials.unitIDs);
    
    
    for Ass_idx = 1:length(FAtrials.units{3})                 % loop across assemblies

        SameAssCoherence{Ass_idx}=[];
        AssMembers      = FAtrials.units{3}{Ass_idx};         % Member units for this assembly
        OtherAssMembers = setdiff(AllAssMembers,AssMembers);
        
        for i = 1:noUnits(3)
            for j = 1:noUnits(3)


                %%%% Conditional sorting of unit pairs  %%%%

                % (1) Members of same assembly
                if      i~=j && ...                                                    % ...Not the same unit
                        ismember(i,AssMembers) && ismember(j,AssMembers) && ...        % ...Both members of the same assembly
                        i<=noUnits(1) && j > noUnits(1) && ...                         % ...Live in the right areas
                        meanRate(i,j) >= rateThresh                                     % ...Sufficient joint firing rate
                        
                        if isempty(tempMembers{i,j}) %... i.e. only count this pair once in case the two units appear in another assembly
                            SameAssCoherence{Ass_idx}   = cat(3,SameAssCoherence{Ass_idx}, [nan(len_-length(coherence{i,j}),size(params.band,1));coherence{i,j}]);
                        end
                        tempMembers{i,j} = [i,j];

                % (2) Members of different assembly
                elseif  i~=j && ...                                                    % ...Not the same unit
                        ismember(i,AssMembers) && ismember(j,OtherAssMembers) && ...   % ...Members of different assemblies
                        i<=noUnits(1) && j > noUnits(1)  && ...                        % ...Live in right areas
                        meanRate(i,j) >= rateThresh                                     % ...Sufficient joint firing rate
                        
                        if isempty(tempCrossMembers{i,j}) %... i.e. only count this pair once in case the two units appear in another assembly
                            DiffAssCoherence{Ass_idx} = cat(3,DiffAssCoherence{Ass_idx}, [nan(len_-length(coherence{i,j}),size(params.band,1));coherence{i,j}]);
                        end
                        tempCrossMembers{i,j} = [i,j];
                else
                end
            end
        end
    end
        
    idx = cellfun(@isempty,tempMembers) & cellfun(@isempty,tempCrossMembers) ;
    idx(:,1:noUnits(1)) = 0; 
    idx(noUnits(2):noUnits(3),:) = 0;
    idx(meanRate<rateThresh)=0;         % restrict by joint firing rate
    
    tempNonMembers = idx;
    
    tempMembers2  = coherence(idx);    
    NonAssCoherence = [];
    for i = 1:length(tempMembers2)
        NonAssCoherence = cat(3,NonAssCoherence, [nan(len_-length(tempMembers2{i}),size(params.band,1));tempMembers2{i}]);
    end
    clear tempMembers2
    % Strip any empty assemblies (i.e. containing no pairs which make the cut)
    SameAssCoherence(cellfun(@isempty,SameAssCoherence))=[];
    DiffAssCoherence(cellfun(@isempty,DiffAssCoherence))=[];
    
    % Sort by frequency bands
    for iBand = 1:size(params.band,1)
        SameAssCoherence_{iBand,1} = [];
        for iAss=1:length(SameAssCoherence)
            SameAssCoherence_{iBand,1} = cat(2,SameAssCoherence_{iBand},squeeze(SameAssCoherence{iAss}(:,iBand,:)));
            SameAssCoherence_{iBand,1}(:,nanvar(SameAssCoherence_{iBand,1})<0.02)=[];
        end
    end
    clear SameAssCoherence
    
    for iBand = 1:size(params.band,1)
        DiffAssCoherence_{iBand,1} = [];
        for iAss=1:length(DiffAssCoherence)
            DiffAssCoherence_{iBand,1} = cat(2,DiffAssCoherence_{iBand},squeeze(DiffAssCoherence{iAss}(:,iBand,:)));
            DiffAssCoherence_{iBand,1}(:,nanvar(DiffAssCoherence_{iBand,1})<0.02)=[];
        end
    end
    clear DiffAssCoherence
       
    for iBand = 1:size(params.band,1)
        NonAssCoherence_{iBand,1} = [];
            NonAssCoherence_{iBand,1} = squeeze(NonAssCoherence(:,1,:));
            NonAssCoherence(:,1,:)=[];
            NonAssCoherence_{iBand,1}(:,nanvar(NonAssCoherence_{iBand,1})<0.02)=[];
    end
    clear NonAssCoherence

    clear  temp temp2 temp3
    %% Collate histograms + plot
    iBand = 2;

    bins =0:0.01:1;
    AssCoherenceBinned=[];
    DiffAssCoherenceBinned=[];
    NonAssCoherenceBinned=[];
    temp  = SameAssCoherence_{iBand};
    temp2 = DiffAssCoherence_{iBand}; 
    temp3 = NonAssCoherence_{iBand};

    for iAss = 1:size(temp,2)
         AssCoherenceBinned(:,iAss) = histc(temp(:,iAss),bins);
    end
    for iAss = 1:size(temp2,2)
         DiffAssCoherenceBinned(:,iAss) = histc(temp2(:,iAss),bins);
    end
    for iAss = 1:size(temp3,2)    
        NonAssCoherenceBinned(:,iAss) = histc(temp3(:,iAss),bins);
    end

    length(bins)
    nBS = 1000;
    NonAssCoherenceBinnedBS=[];
    AssCoherenceBinnedBS=[];
    DiffAssCoherenceBinnedBS=[];

    % draw a 50% sample of possible pairs for each bootstrap
    k = round(0.5*(min([size(AssCoherenceBinned,2),...
                        size(DiffAssCoherenceBinned,2),...
                        size(NonAssCoherenceBinned,2),...
                        ])));

    for iBS = 1:nBS
       fprintf('Bootsrap draw no. %d of %d\n',iBS,nBS)
       AssCoherenceBinnedBS     (:,iBS) = nanmean(AssCoherenceBinned(:,randperm(size(AssCoherenceBinned,2),k)),2);
       DiffAssCoherenceBinnedBS (:,iBS) = nanmean(DiffAssCoherenceBinned(:,randperm(size(DiffAssCoherenceBinned,2),k)),2);
       NonAssCoherenceBinnedBS  (:,iBS) = nanmean(NonAssCoherenceBinned(:,randperm(size(NonAssCoherenceBinned,2),k)),2);
    end
    AssCoherenceBinnedBS_mean      = nanmean(AssCoherenceBinnedBS');
    NonAssCoherenceBinnedBS_mean   = nanmean(NonAssCoherenceBinnedBS');
    OtherAssCoherenceBinnedBS_mean = nanmean(DiffAssCoherenceBinnedBS');

    AssCoherenceBinnedBS_CI        = nanstd(AssCoherenceBinnedBS');
    NonAssCoherenceBinnedBS_CI     = nanstd(NonAssCoherenceBinnedBS');
    OtherAssCoherenceBinnedBS_CI   = nanstd(DiffAssCoherenceBinnedBS');

    AssCoherenceBinnedBS_mean      = AssCoherenceBinnedBS_mean./length(bins);
    OtherAssCoherenceBinnedBS_mean = OtherAssCoherenceBinnedBS_mean./length(bins);
    NonAssCoherenceBinnedBS_mean   = NonAssCoherenceBinnedBS_mean./length(bins);

    AssCoherenceBinnedBS_CI        = AssCoherenceBinnedBS_CI./length(bins);
    NonAssCoherenceBinnedBS_CI     = NonAssCoherenceBinnedBS_CI./length(bins);
    OtherAssCoherenceBinnedBS_CI   = OtherAssCoherenceBinnedBS_CI./length(bins);

    figure; hold on
    plot(bins,AssCoherenceBinnedBS_mean,'b','LineWidth',2)
    plot(bins,OtherAssCoherenceBinnedBS_mean,'g','LineWidth',2)
    plot(bins,NonAssCoherenceBinnedBS_mean,'k','LineWidth',2)
    legend('Same CA1-mPFC assembly member pairs ',...
           'Different CA1-mPFC assembly member pairs ',...
           'Non-member CA1-mPFC pairs'); legend boxoff
    ciplot(AssCoherenceBinnedBS_mean+AssCoherenceBinnedBS_CI,...
           AssCoherenceBinnedBS_mean-AssCoherenceBinnedBS_CI,...
           bins,'b')
    ciplot(OtherAssCoherenceBinnedBS_mean+OtherAssCoherenceBinnedBS_CI,...
           OtherAssCoherenceBinnedBS_mean-OtherAssCoherenceBinnedBS_CI,...
           bins,'g')   
    ciplot(NonAssCoherenceBinnedBS_mean+NonAssCoherenceBinnedBS_CI,...
           NonAssCoherenceBinnedBS_mean-NonAssCoherenceBinnedBS_CI,...
           bins,'k')
    axis([0 1 0 1])
    view(90, -90)
    xlabel('4-5Hz Spike train Coherence')
    ylabel('Distribution of experiment time')
    %% Cut outs and trial-averaging
    twin_ = 10;
    SameAssCoherenceMean = [];
    DiffAssCoherenceMean = [];
    NonAssCoherenceMean  = [];
     
    % (1) same assembly member pairs
    [cutout_t,cutout]=SelTimesCohereogram(twin_,t_(find(~cellfun(@isempty,tempMembers))),coherence(find(~cellfun(@isempty,tempMembers))),FRtrials.EvtTs);
    bw_ = ceil(mean(diff(cutout_t{1}{1}))); ko = round((twin_+5)/bw_);
    cutoutSel = [];
    for iPair=1:length(cutout_t)
        for iTrial=1:length(cutout_t{iPair})
            try
                cutout_t{iPair}{iTrial} = cutout_t{iPair}{iTrial}([1:ko end-ko+1:end]);
                cutout{iPair}{iTrial}   = cutout  {iPair}{iTrial}([1:ko end-ko+1:end],:)';
            catch
                cutout_t{iPair}{iTrial}=[];
                cutout{iPair}{iTrial}=[];
            end
        end;
        for iBand = 1:size(params.band,1)
            for iTrial=1:length(cutout_t{iPair})
                if ~isempty(cutout_t{iPair}{iTrial})
                    cutoutSel{iPair}{iBand}(iTrial,:) =  cutout{iPair}{iTrial}(iBand,:);
                else
                    cutoutSel{iPair}{iBand}(iTrial,:) = nan(2*ko,1);
                end
            end
            SameAssCoherenceMean{iBand}(iPair,:) =  nanmean(cutoutSel{iPair}{iBand},1);
        end
    end;
    
    % (2) Different assembly member pairs
    [cutout_t,cutout]=SelTimesCohereogram(twin_,t_(find(~cellfun(@isempty,tempCrossMembers))),coherence(find(~cellfun(@isempty,tempCrossMembers))),FRtrials.EvtTs);
    bw_ = ceil(mean(diff(cutout_t{1}{1}))); ko = round((twin_+5)/bw_);
    cutoutSel = [];
    for iPair=1:length(cutout_t)
        for iTrial=1:length(cutout_t{iPair})
            try
                cutout_t{iPair}{iTrial} = cutout_t{iPair}{iTrial}([1:ko end-ko+1:end]);
                cutout{iPair}{iTrial}   = cutout  {iPair}{iTrial}([1:ko end-ko+1:end],:)';
            catch
                cutout_t{iPair}{iTrial}=[];
                cutout{iPair}{iTrial}=[];
            end
        end
        for iBand = 1:size(params.band,1)
            for iTrial=1:length(cutout_t{iPair})
                if ~isempty(cutout_t{iPair}{iTrial})
                    cutoutSel{iPair}{iBand}(iTrial,:) =  cutout{iPair}{iTrial}(iBand,:);
                else
                    cutoutSel{iPair}{iBand}(iTrial,:) = nan(2*ko,1);
                end
            end
            DiffAssCoherenceMean{iBand}(iPair,:) =  nanmean(cutoutSel{iPair}{iBand},1);
        end
    end;
    
    % (3) Non assembly member pairs
    [cutout_t,cutout]=SelTimesCohereogram(twin_,t_(find(tempNonMembers)),coherence(find(tempNonMembers)),FRtrials.EvtTs);
    bw_ = ceil(mean(diff(cutout_t{1}{1}))); ko = round((twin_+5)/bw_);
    cutoutSel = [];
    for iPair=1:length(cutout_t)
        for iTrial=1:length(cutout_t{iPair})
            try
                cutout_t{iPair}{iTrial} = cutout_t{iPair}{iTrial}([1:ko end-ko+1:end]);
                cutout{iPair}{iTrial}   = cutout  {iPair}{iTrial}([1:ko end-ko+1:end],:)';
            catch
                cutout_t{iPair}{iTrial}=[];
                cutout{iPair}{iTrial}=[];
            end
        end
        for iBand = 1:size(params.band,1)
            for iTrial=1:length(cutout_t{iPair})
                if ~isempty(cutout_t{iPair}{iTrial})
                    cutoutSel{iPair}{iBand}(iTrial,:) =  cutout{iPair}{iTrial}(iBand,:);
                else
                    cutoutSel{iPair}{iBand}(iTrial,:) = nan(2*ko,1);
                end
            end
            NonAssCoherenceMean{iBand}(iPair,:) =  nanmean(cutoutSel{iPair}{iBand},1);
        end
    end;

     tb = 1:2*ko;
     
     for iBand = 1:8
         
         figure; hold on
         plot(tb,nanmean(SameAssCoherenceMean{iBand}),'b','LineWidth',1.5);
         plot(tb,nanmean(DiffAssCoherenceMean{iBand}),'g','LineWidth',1.5);
         plot(tb,nanmean(NonAssCoherenceMean{iBand}),'k','LineWidth',1.5);
         
         ciplot(nanmean(SameAssCoherenceMean{iBand}) + nansem(SameAssCoherenceMean{iBand}),...
             nanmean(SameAssCoherenceMean{iBand}) - nansem(SameAssCoherenceMean{iBand}),...
             tb,'b')
         
         ciplot(nanmean(DiffAssCoherenceMean{iBand}) + nansem(DiffAssCoherenceMean{iBand}),...
             nanmean(DiffAssCoherenceMean{iBand}) - nansem(DiffAssCoherenceMean{iBand}),...
             tb,'g')
         
         ciplot(nanmean(NonAssCoherenceMean{iBand})  + nansem(NonAssCoherenceMean{iBand}),...
             nanmean(NonAssCoherenceMean{iBand}) - nansem(NonAssCoherenceMean{iBand}),...
             tb,'k')
         
         axis([min(tb) max(tb) 0.2 0.4])
     end
end




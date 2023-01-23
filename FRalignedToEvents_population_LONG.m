clear
Targets = {'SHORT','MEDIUM','LONG'};
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Dirs = {'Left','Right'};
Target = Targets{3};
bw=0.05;
bw_hist = 0.2;
minFR = 0.1;
DataType = 'SpikeCounts';%'KDE'
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
else
    pat = '/Volumes/HDD2/DNMTP/raw';
end    
cd(pat)
fileList=dir(sprintf('allTimestamps%s*%s*',filesep,Target));
% fileList=dir('allTimestamps\*LONG*.mat');

reject_list={'IreneuszLONG1_Events.mat'}; 
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag
%%
Areas = {'HP','PFC'};

tlims= [-5 10];
length_ = sum(abs([tlims]))/bw;
length_hist = sum(abs([tlims]))/bw_hist;
tb = (1:length_)*bw +tlims(1);
tb_hist = (1:length_hist)*bw_hist +tlims(1);
Events_ = {'CueLight','SamplePress','DelayEnd','NosePoke','ChoicePress','RewardConsume'};
Events__ = {'Cue Light','Sample Press','Delay End','Nose Poke','Choice Press','Reward Collection'};
%%%%%
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    load(fullfile(pat, 'allTimestamps',[fname,'_Events.mat']))
    load(fullfile(pat, [fname,'.mat']))
    
    
    for iArea = 1:length(Areas)
        
        fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})       
        load(fullfile(pat, 'KDE_binsTaskonly',sprintf('%s_%s_iFR50_behavOnly.mat',fname,Areas{iArea})),'iFR','Tmtx');
        
        switch DataType
            case 'KDE'
                % data_{1}.iFR = zscore(iFR);
                data_{1}.iFR =(iFR);
            case 'SpikeCounts'
                nU = eval(sprintf('length(%scells)',Areas{iArea}));
                
                TbaseSpikes = min(Tmtx):bw_hist:max(Tmtx);
                data_{1}.iFR = zeros(length(TbaseSpikes),nU);
                dataZ_{1}.iFR = zeros(length(TbaseSpikes),nU);
                for iUnit =1:nU
                    ST = eval(sprintf('%scells{iUnit}.t*1e-4',Areas{iArea}));
                    data_{1}.iFR(:,iUnit) = eval(sprintf('histc(%scells{iUnit}.t*1e-4,TbaseSpikes);',Areas{iArea}));
                    dataZ_{1}.iFR(:,iUnit) = zscore(data_{1}.iFR(:,iUnit));
                    meanFR(iUnit) = mean(data_{1}.iFR(:,iUnit))./bw_hist; %length(ST)/range(ST);
                    stdFR(iUnit)  = std(data_{1}.iFR(:,iUnit))./bw_hist;
                end
        end

        for iDelay = 1:length(Delays_)
            Delay_ = Delays_{iDelay};
            %%%%%%%%%%%%%%%% Correct trials
            for iEvent = 1:length(Events_)
                times_ = [eval(sprintf('t.%s.%s_LeftCorrect',Delay_,Events_{iEvent})); ...
                          eval(sprintf('t.%s.%s_RightCorrect',Delay_,Events_{iEvent}))];
                      switch DataType
                          case 'KDE'
                              [~,cutout] = SelTimesFRs(tlims,Tmtx,data_,times_*1e-6);
                          case 'SpikeCounts'
                              [~,cutout] = SelTimesFRs(tlims,TbaseSpikes,dataZ_,times_*1e-6);
                      end
                      
                
                cutout_     = cell(size(iFR,2),1);
                cutout_Hist = cell(size(iFR,2),1);
                for iUnit=1:size(iFR,2)
                    eval(sprintf('ST = %scells{iUnit}.t*1e-4;',Areas{iArea}));
                    for iTrial = 1:size(times_,1)
                        bins_ =  ((times_(iTrial)*1e-6)+tlims(1)+bw_hist):bw_hist:(times_(iTrial)*1e-6)+tlims(2);
                        cutout_Hist{iUnit}(iTrial,1:length_hist) = histc(ST,bins_);
                        if ~isempty(cutout{1}{iTrial})
                            cutout_{iUnit}(iTrial,1:length_hist)     =  cutout{1}{iTrial}(1:length_hist,iUnit);
                           %cutout_{iUnit}(iTrial,1:length_)     = zscore(cutout_{iUnit}(iTrial,1:length_) - mean(cutout_{iUnit}(iTrial,1:ceil(0.9*abs(tlims(1))/bw))));
                        else
                            cutout_{iUnit}(iTrial,1:length_hist) = NaN;
                        end
                    end
                end
                FiringRatesMean{iEvent}{iArea}{iDelay}{iFile}     = cell2mat(cellfun(@nanmean,cutout_,'UniformOutput',false))';
                FiringRatesMeanHist{iEvent}{iArea}{iDelay}{iFile} = cell2mat(cellfun(@nanmean,cutout_Hist,'UniformOutput',false))';
                FiringRatesVarHist{iEvent}{iArea}{iDelay}{iFile}  = cell2mat(cellfun(@nanvar,cutout_Hist,'UniformOutput',false))';
            end
            
            times_ = eval(sprintf('t.%s.CueLight_LeftCorrect',Delay_));
            nC_L = length(eval(sprintf('t.%s.CueLight_LeftCorrect',Delay_)));
            nC_R = length(eval(sprintf('t.%s.CueLight_RightCorrect',Delay_)));
            nE_L = length(eval(sprintf('t.%s.CueLight_LeftError',Delay_)));
            nE_R = length(eval(sprintf('t.%s.CueLight_RightError',Delay_)));
            
            for iTrial = 1:nC_L
                times_ =  [eval(sprintf('t.%s.CueLight_LeftCorrect(iTrial)',Delay_)), ...
                          eval(sprintf('t.%s.SamplePress_LeftCorrect(iTrial)',Delay_)), ...
                          eval(sprintf('t.%s.DelayEnd_LeftCorrect(iTrial)',Delay_)), ...
                          eval(sprintf('t.%s.NosePoke_LeftCorrect(iTrial)',Delay_)), ...
                          eval(sprintf('t.%s.ChoicePress_LeftCorrect(iTrial)',Delay_)), ...
                          eval(sprintf('t.%s.RewardConsume_LeftCorrect(iTrial)',Delay_))]*1e-6;
%                 tlims = [-5, times_(end)-times_(1) + 5];
%                   [Tcutout,cutout] = SelTimesFRs(tlims,Tmtx,data_,times_(1));
%                   cutout{1}{1}
            end
            eval(sprintf('nU = length(%scells);',Areas{iArea}));
            clear cutout_Hist cutout_Hist_CL cutout_Hist_CR cutout_Hist_EL cutout_Hist_ER

            for iUnit=1:nU
                
                eval(sprintf('ST = %scells{iUnit}.t*1e-4;',Areas{iArea}));
                
                clear cutout_Hist 
                for iTrial = 1:nC_L 
                    
                      times_ =  [eval(sprintf('t.%s.CueLight_LeftCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.SamplePress_LeftCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.DelayEnd_LeftCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.NosePoke_LeftCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.ChoicePress_LeftCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.RewardConsume_LeftCorrect(iTrial)',Delay_))]*1e-6;
                  
                    
                    bins_ =  (times_(1)-10):1:(times_(end)+10);
                    
                    cutout_Hist{iTrial}(:,iUnit) = zscore(histc(ST,bins_));
                    for i= 1:length(times_)
                        if i<length(times_)
                            tIdx = bins_>=times_(i) & bins_<=times_(i+1);
                        else
                            tIdx = bins_>=times_(i) ;
                        end
                        cutout_Hist_CL{iUnit}(iTrial,i) = nanmean(cutout_Hist{iTrial}(tIdx,iUnit));
                      
                    end
                end
                
                clear cutout_Hist 
                for iTrial = 1:nC_R
                    
                    times_ =  [eval(sprintf('t.%s.CueLight_RightCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.SamplePress_RightCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.DelayEnd_RightCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.NosePoke_RightCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.ChoicePress_RightCorrect(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.RewardConsume_RightCorrect(iTrial)',Delay_))]*1e-6;
                    
                    bins_ =  (times_(1)-10):1:(times_(end)+10);
                    
                    cutout_Hist{iTrial}(:,iUnit) = zscore(histc(ST,bins_));
                    for i= 1:length(times_)
                        if i<length(times_)
                            tIdx = bins_>=times_(i) & bins_<=times_(i+1);
                        else
                            tIdx = bins_>=times_(i) ;
                        end
                        cutout_Hist_CR{iUnit}(iTrial,i) = nanmean(cutout_Hist{iTrial}(tIdx,iUnit));
                      
                    end
                end
                
                clear cutout_Hist 
                for iTrial = 1:nE_L
                    
                    times_ =  [eval(sprintf('t.%s.CueLight_LeftError(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.SamplePress_LeftError(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.DelayEnd_LeftError(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.NosePoke_LeftError(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.ChoicePress_LeftError(iTrial)',Delay_))]*1e-6;
                    
                    bins_ =  (times_(1)-10):1:(times_(end)+10);
                    
                    cutout_Hist{iTrial}(:,iUnit) = zscore(histc(ST,bins_));
                    for i= 1:length(times_)
                        if i<length(times_)
                            tIdx = bins_>=times_(i) & bins_<=times_(i+1);
                        else
                            tIdx = bins_>=times_(i) ;
                        end
                        cutout_Hist_EL{iUnit}(iTrial,i) = nanmean(cutout_Hist{iTrial}(tIdx,iUnit));
                    end
                end
                
                  clear cutout_Hist 
                for iTrial = 1:nE_R
                    
                    times_ =  [eval(sprintf('t.%s.CueLight_RightError(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.SamplePress_RightError(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.DelayEnd_RightError(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.NosePoke_RightError(iTrial)',Delay_)), ...
                        eval(sprintf('t.%s.ChoicePress_RightError(iTrial)',Delay_))]*1e-6;
                    
                    bins_ =  (times_(1)-10):1:(times_(end)+10);
                    
                    cutout_Hist{iTrial}(:,iUnit) = zscore(histc(ST,bins_));
                    for i= 1:length(times_)
                        if i<length(times_)
                            tIdx = bins_>=times_(i) & bins_<=times_(i+1);
                        else
                            tIdx = bins_>=times_(i) ;
                        end
                        cutout_Hist_ER{iUnit}(iTrial,i) = nanmean(cutout_Hist{iTrial}(tIdx,iUnit));
                    end
                end
                
            end
            cutout_Hist_CL = cell2mat(cellfun(@nanmean,cutout_Hist_CL,'UniformOutput',false)');
            cutout_Hist_CR = cell2mat(cellfun(@nanmean,cutout_Hist_CR,'UniformOutput',false)');
            cutout_Hist_EL = cell2mat(cellfun(@nanmean,cutout_Hist_EL,'UniformOutput',false)');
            cutout_Hist_ER = cell2mat(cellfun(@nanmean,cutout_Hist_ER,'UniformOutput',false)');
            
            (cutout_Hist_CL + cutout_Hist_CR)./2;
            (cutout_Hist_EL + cutout_Hist_ER)./2;
            
            %%%%%%%%%%%%%%%% Error trials
            for iEvent = 1:length(Events_)-1
                times_ = [eval(sprintf('t.%s.%s_LeftError',Delay_,Events_{iEvent}))';eval(sprintf('t.%s.%s_RightError',Delay_,Events_{iEvent}))'];
                [~,cutout] = SelTimesFRs(tlims,Tmtx,data_,times_*1e-6);
                
                cutout_ = cell(size(iFR,2),1);
                cutout_Hist = cell(size(iFR,2),1);
                for iUnit=1:size(iFR,2)
                    eval(sprintf('ST = %scells{iUnit}.t*1e-4;',Areas{iArea}));
                    for iTrial = 1:size(times_,1)
                        bins_ =  ((times_(iTrial)*1e-6)+tlims(1)+bw_hist):bw_hist:(times_(iTrial)*1e-6)+tlims(2);
                        cutout_Hist{iUnit}(iTrial,1:length_hist) = histc(ST,bins_);
                        if ~isempty(cutout{1}{iTrial})
                            cutout_{iUnit}(iTrial,1:length_) = cutout{1}{iTrial}(1:length_,iUnit);
                        else
                            cutout_{iUnit}(iTrial,1:length_) = NaN;
                        end
                    end
                end
                %                 if length(times_)>0%size(cutout_{1},1)>1
                
                if size(cutout_{1},1)>1
                    FiringRatesErrorsMean{iEvent}{iArea}{iDelay}{iFile} = cell2mat(cellfun(@nanmean,cutout_,'UniformOutput',false))';
                    FiringRatesErrorsMeanHist{iEvent}{iArea}{iDelay}{iFile} = cell2mat(cellfun(@nanmean,cutout_Hist,'UniformOutput',false))';
                    FiringRatesErrorsVarHist{iEvent}{iArea}{iDelay}{iFile} = cell2mat(cellfun(@nanvar,cutout_Hist,'UniformOutput',false))';
                    %                 elseif length(times_)==1
                    %                     FiringRatesErrorsMean{iEvent}{iArea}{iDelay}{iFile} = NaN(size(cell2mat(cutout_)'));
                    %                 elseif length(times_)==0
                else
                    FiringRatesErrorsMean{iEvent}{iArea}{iDelay}{iFile} = NaN(size(cell2mat(cutout_)'));
                    FiringRatesErrorsMeanHist{iEvent}{iArea}{iDelay}{iFile} = NaN(size(cell2mat(cutout_Hist)'));
                    FiringRatesErrorsVarHist{iEvent}{iArea}{iDelay}{iFile}  = NaN(size(cell2mat(cutout_Hist)'));
                end
            end
            
        end
    end
end

clear cutout cutout_ cutout_Hist data_ Delay_ ERRORtrangeleft_choice ERRORtrangeleft_sample ERRORtrangeright_sample ERRORtrangeright_choice
clear HPcells HPinter HPtonic iArea iDelay iEvent iFile iFR inputData iTrial iUnit PFCcells PFCinter PFCtonic t times_ Tmtx trangeleft_choice trangeleft_sample trangeright_choice trangeright_sample

%% collapse
load('blue_white_red.mat')
color_={'r','b','g'};
for iArea = 1:length(Areas)
    for iDelay = 1:length(Delays_)
        for iEvent = 1:length(Events_)
            FiringRatesMeancollapse{iEvent}{iArea}{iDelay}     = cell2mat(FiringRatesMean{iEvent}{iArea}{iDelay});
            FiringRatesMeanHistcollapse{iEvent}{iArea}{iDelay} = cell2mat(FiringRatesMeanHist{iEvent}{iArea}{iDelay});
            FiringRatesVarHistcollapse{iEvent}{iArea}{iDelay}  = cell2mat(FiringRatesVarHist{iEvent}{iArea}{iDelay});
        end
        for iEvent = 1:length(Events_)-1
            FiringRatesErrorsMeancollapse{iEvent}{iArea}{iDelay}     = cell2mat(FiringRatesErrorsMean{iEvent}{iArea}{iDelay});
            FiringRatesErrorsMeanHistcollapse{iEvent}{iArea}{iDelay} = cell2mat(FiringRatesErrorsMeanHist{iEvent}{iArea}{iDelay});
            FiringRatesErrorsVarHistcollapse{iEvent}{iArea}{iDelay}  = cell2mat(FiringRatesErrorsVarHist{iEvent}{iArea}{iDelay});
        end
    end
end

%% Plot staggered rate activation maps for each event

for iArea = 1:length(Areas)
    for iEvent = 1:length(Events_)
        
        figure('color','w','name',Events_{iEvent})
        for iDelay = 1:length(Delays_)
            [~,i]=max(FiringRatesMeancollapse{iEvent}{iArea}{iDelay});
            [~,i]=sort(i);
            N = max(i);
            subplot(1,length(Delays_),iDelay); hold on
            title(sprintf('%s units: %s delay',Areas{iArea},Delays__{iDelay}))
            imagesc(tb,1:N,FiringRatesMeancollapse{iEvent}{iArea}{iDelay}(:,i)');
            %colormap jet
            colormap(cmap)
            caxis([-3 3])
            
            % Vertical bars
            plot([0 0],[0 N+1],'k')
            axis([min(tb) max(tb) 0 max(i)])
            
            axis off
            
        end
        plot([2, 7],[10,10],'-k','LineWidth',1.5)
        plot([2, 2],[10,60],'-k','LineWidth',1.5)
    end
end
%% Plot staggered spike count histogram activation maps for each event

for iArea = 1:length(Areas)
    for iEvent = 1:length(Events_)
        
        figure('color','w','name',Events_{iEvent})
        for iDelay = 1:length(Delays_)
            [~,i]=max(FiringRatesMeanHistcollapse{iEvent}{iArea}{iDelay});
            [~,i]=sort(i);
            N = max(i);
            subplot(1,length(Delays_),iDelay); hold on
            title(sprintf('%s units: %s delay',Areas{iArea},Delays__{iDelay}))
            imagesc(tb,1:N,FiringRatesMeanHistcollapse{iEvent}{iArea}{iDelay}(:,i)');
            %colormap jet
            colormap(cmap)
            caxis([-3 3])
            
            % Vertical bars
            plot([0 0],[0 N+1],'k')
            axis([min(tb) max(tb) 0 max(i)])
            
            axis off
            
        end
        plot([2, 7],[10,10],'-k','LineWidth',1.5)
        plot([2, 2],[10,60],'-k','LineWidth',1.5)
    end
end
%% Plot average firing rate foreach event
offset = [0.55 0.475 0.4];

for iArea = 1:length(Areas)
    figure('color','w','name',Areas{iArea})
    for iEvent = 1:length(Events_)-1
        subplot(1,length(Events_),iEvent); hold on
        for iDelay = 1:length(Delays_)
            m = nanmean(FiringRatesErrorsMeancollapse{iEvent}{iArea}{iDelay},2);
            s = nansem(FiringRatesErrorsMeancollapse{iEvent}{iArea}{iDelay},2);
            ciplot(m+s,m-s,tb,color_{iDelay},0.2);
        end
    end
    for iEvent = 1:length(Events_)
        subplot(1,length(Events_),iEvent); hold on
        for iDelay = 1:length(Delays_)
            m = nanmean(FiringRatesMeancollapse{iEvent}{iArea}{iDelay},2);
            s = nansem(FiringRatesMeancollapse{iEvent}{iArea}{iDelay},2);
            ciplot(m+s,m-s,tb,color_{iDelay},1);
        end
        title(Events__{iEvent})
        plot([0 0],[-0.1 0.2],'k')
        axis([min(tb) max(tb) -0.5 1])
        axis off
    end
    
    for iEvent = 1:length(Events_)
        subplot(1,length(Events_),iEvent); hold on
        for iDelay = 1:length(Delays_)
            x = FiringRatesMeancollapse{iEvent}{iArea}{iDelay};
            y = FiringRatesErrorsMeancollapse{iEvent}{iArea}{iDelay};
            [p,~] = permtest2vec(x,y,5000,0.05);
            a = nan(size(p));a(p)=1;
            %          a(x==y)=NaN;
            plot(tb,a*offset(iDelay),color_{iDelay},'LineWidth',4)
        end
    end
    plot([2 7],[0.4 0.4],'k','LineWidth',1.5)
    plot([7 7],[0.4 0.6],'k','LineWidth',1.5)
    
end
%% Plot average spike count histogram for each event
offset = [0.55 0.475 0.4];

for iArea = 1:length(Areas)
    figure('color','w','name',Areas{iArea})
    for iEvent = 1:length(Events_)-1
        subplot(1,length(Events_),iEvent); hold on
        for iDelay = 1:length(Delays_)
            x = FiringRatesErrorsMeanHistcollapse{iEvent}{iArea}{iDelay};
            x = zscore(x);
            m = nanmean(x,2);
            s = nansem(x,2);
            ciplot(m+s,m-s,tb_hist,color_{iDelay},0.2);
        end
    end
    for iEvent = 1:length(Events_)
        subplot(1,length(Events_),iEvent); hold on
        for iDelay = 1:length(Delays_)
            x = FiringRatesMeanHistcollapse{iEvent}{iArea}{iDelay};
            x= zscore(x);
            m = nanmean(x,2);
            s = nansem(x,2);
            ciplot(m+s,m-s,tb_hist,color_{iDelay},1);
        end
        title(Events__{iEvent})
        plot([0 0],[-0.1 0.2],'k')
        axis([min(tb_hist) max(tb_hist) -0.5 1])
        axis off
    end
    
    for iEvent = 1:length(Events_)-1
        subplot(1,length(Events_),iEvent); hold on
        for iDelay = 1:length(Delays_)
            x = zscore(FiringRatesMeanHistcollapse{iEvent}{iArea}{iDelay});
            y = zscore(FiringRatesErrorsMeanHistcollapse{iEvent}{iArea}{iDelay});
            [p,~] = permtest2vec(x,y,500,0.05);
            a = nan(size(p));a(p)=1;
            %          a(x==y)=NaN;
            plot(tb_hist,1.5*a*offset(iDelay),color_{iDelay},'LineWidth',4)
        end
    end
    %             plot([2 7],[0.4 0.4],'k','LineWidth',1.5)
    %             plot([7 7],[0.4 0.6],'k','LineWidth',1.5)
    
end
%% Fraction of activated/inactivated cells - correct trials
useHists = false;
SD_thresh = 2;
Dur_thresh = 3;
win_Baseline  = find(tb_hist<0);
win_Analysis  = [0 5];

% Correct trials
for iArea = 1:length(Areas)
    for iEvent = 1:length(Events_)
        for iDelay = 1:length(Delays_)
            if useHists 
                x = FiringRatesMeanHistcollapse{iEvent}{iArea}{iDelay};
            else
                x = FiringRatesMeancollapse{iEvent}{iArea}{iDelay};
            end
%             for i=1:size(x,2)
%                 x(:,i) = sgolayfilt(x(:,i),7,21);
%             end
            x = x - repmat(mean(x(win_Baseline,:)),size(x,1),1);
            x=x./repmat(std(x(win_Baseline,:)),size(x,1),1);
            plot(tb_hist,x)
            % imagesc(x); caxis([0 2])
            sig_UP   = zeros(size(x,2),1);
            sig_DOWN = zeros(size(x,2),1);
            
            for i=1:size(x,2)
                x_ = x(tb_hist<win_Analysis(1) | tb_hist>win_Analysis(2),i);
                [~,~,w] = findpeaks(x_,'MinPeakHeight',SD_thresh);
                if sum ( w>=Dur_thresh ) >0
                    sig_UP(i) = sum ( w>=Dur_thresh );
                end
                
                [~,~,w] = findpeaks(-x_,'MinPeakHeight',SD_thresh);
                if sum ( w>=Dur_thresh ) >0
                    sig_DOWN(i) = sum ( w>=Dur_thresh );
                end
            end
            
            Mod_ = zeros(size(x,2),1);
            
            NonMod = sum([sig_UP, sig_DOWN],2) == 0; 
            
            Mod_(sig_DOWN>sig_UP   & ~NonMod)  = -1;
            Mod_(sig_DOWN<sig_DOWN & ~NonMod)  = 1;
            
            Mod_(min([sig_UP, sig_DOWN],[],2)==0 & sig_UP>0)   = 1;
            Mod_(min([sig_UP, sig_DOWN],[],2)==0 & sig_DOWN>0) = -1;
            
            ModInd{iArea}{iEvent}{iDelay}=[sum(Mod_==-1),sum(Mod_==0),sum(Mod_==1)];
        end
    end
end

% Error trials
for iArea = 1:length(Areas)
    for iEvent = 1:length(Events_)-1
        for iDelay = 1:length(Delays_)
            if useHists
                x = FiringRatesErrorsMeanHistcollapse{iEvent}{iArea}{iDelay};
            else
                x = FiringRatesErrorsMeancollapse{iEvent}{iArea}{iDelay};
            end
%             for i=1:size(x,2)
%                 x(:,i) = sgolayfilt(x(:,i),7,21);
%             end
            x = x - repmat(mean(x(win_Baseline,:)),size(x,1),1);
            x=x./repmat(std(x(win_Baseline,:)),size(x,1),1);
            plot(tb_hist,x)
            % imagesc(x); caxis([0 2])
            sig_UP   = zeros(size(x,2),1);
            sig_DOWN = zeros(size(x,2),1);
            
            for i=1:size(x,2)
                x_ = x(tb_hist<win_Analysis(1) | tb_hist>win_Analysis(2),i);
                [~,~,w] = findpeaks(x_,'MinPeakHeight',SD_thresh);
                if sum ( w>=Dur_thresh ) >0
                    sig_UP(i) = sum ( w>=Dur_thresh );
                end
                
                [~,~,w] = findpeaks(-x_,'MinPeakHeight',SD_thresh);
                if sum ( w>=Dur_thresh ) >0
                    sig_DOWN(i) = sum ( w>=Dur_thresh );
                end
            end
            
            Mod_ = zeros(size(x,2),1);
            
            NonMod = sum([sig_UP, sig_DOWN],2) == 0; 
            
            Mod_(sig_DOWN>sig_UP   & ~NonMod)  = -1;
            Mod_(sig_DOWN<sig_DOWN & ~NonMod)  = 1;
            
            Mod_(min([sig_UP, sig_DOWN],[],2)==0 & sig_UP>0)   = 1;
            Mod_(min([sig_UP, sig_DOWN],[],2)==0 & sig_DOWN>0) = -1;
            
            ModIndErr{iArea}{iEvent}{iDelay}=[sum(Mod_==-1),sum(Mod_==0),sum(Mod_==1)];
        end
    end
end

clear w x_ x sig_UP sig_DOWN Mod_ NonMod iArea iEvent iDelay
%%
x_ = [0.2 ; 0.4; 0.6];
x = x_;
for i=1:length(Events_)-1
    x = [x; x_+i];
end
x=repmat(x,1,3);

% Correct trials
clear y
y{1}=[NaN NaN NaN];y{2}=y{1};
for iArea = 1:length(Areas)
    for iEvent = 1:length(Events_)
        y{iArea} = [y{iArea}; ...
                    ModInd{iArea}{iEvent}{1}./sum(ModInd{iArea}{iEvent}{1}); ...
                    ModInd{iArea}{iEvent}{2}./sum(ModInd{iArea}{iEvent}{2}); ...
                    ModInd{iArea}{iEvent}{3}./sum(ModInd{iArea}{iEvent}{3})];
    end
    y{iArea}(1,:)=[];
end

figure; hold on
for iArea = 1:length(Areas)
   subplot(length(Areas),1,iArea)
   bar(x,y{iArea},'stacked')
   title([Areas{iArea} ' Units (Correct trials)'])
   ylabel('Fraction of responsive units')
end
ticks_ = [.4:1:5.4];
set(gca,'XTick',ticks_,'XTickLabel',Events__,'XTickLabelRotation',45)

% Error trials
clear y
y{1}=[NaN NaN NaN];y{2}=y{1};
for iArea = 1:length(Areas)
    for iEvent = 1:length(Events_)-1
        y{iArea} = [y{iArea}; ...
                    ModIndErr{iArea}{iEvent}{1}./sum(ModIndErr{iArea}{iEvent}{1}); ...
                    ModIndErr{iArea}{iEvent}{2}./sum(ModIndErr{iArea}{iEvent}{2}); ...
                    ModIndErr{iArea}{iEvent}{3}./sum(ModIndErr{iArea}{iEvent}{3})];
    end
    y{iArea}(1,:)=[];
end

figure; hold on
for iArea = 1:length(Areas)
   subplot(length(Areas),1,iArea)
   bar(x(1:end-3,:),y{iArea},'stacked')
   title([Areas{iArea} ' Units (Error trials)'])
   ylabel('Fraction of responsive units')
end
set(gca,'XTick',ticks_(1:end-1),'XTickLabel',Events__(1:end-1),'XTickLabelRotation',45)

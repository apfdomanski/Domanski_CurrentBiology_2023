clear

pat = 'C:\Analysis\AssemblyAnalysis\raw\';

Targets = {'SHORT','MEDIUM','LONG'};
Dirs = {'Left','Right'};
Areas = {'HP','PFC'};
bw=0.05;
tlimsSample = [-5 5];
tlimsChoice = [-5 5];
length_ = sum(abs([tlimsSample, tlimsChoice]))/bw;
tb = (1:length_)*bw;
% reject_list={'MiroslawMEDIUM2_Events.mat','NikodemMEDIUM2_Events.mat'...
%              'NorbertMEDIUM1_Events.mat','NorbertMEDIUM2_Events.mat',...
%              'OnufryMEDIUM1_Events.mat','OnufryMEDIUM2_Events.mat'};
 reject_list={'NorbertMEDIUM1_Events.mat','NorbertMEDIUM2_Events.mat',...
             'OnufryMEDIUM1_Events.mat','OnufryMEDIUM2_Events.mat'};
DelaysNames = {{'No Delay'},{'2s Delay','4s Delay','6s Delay','8s Delay'},{'4s Delay','8s Delay','16s Delay'}};

FanoMinTrialCount = 10;
%% batch process firing rates
for iTarget = 1:length(Targets)
    fileList = dir([pat 'allTimestamps' filesep '*' Targets{iTarget} '*.mat']);
    
    % Remove bad files
    name_flag=zeros(numel(fileList),1);
    for idx=1:numel(fileList)
        fnames_{idx,1}=fileList(idx).name;
        name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
    end
    try
        fileList(find(name_flag))=[];% fnames_(name_flag)=[];
    end
    
    for iFile =1:length(fileList)
        
        
        % Get raw spike data
        fname=strtok(fileList(iFile).name,'_');
        load(sprintf('%s%s.mat',pat,fname));
        % Get event data
        load([pat 'allTimestamps' filesep fileList(iFile).name],'t');
        for iArea = 1:length(Areas)
            % Get convolved spike rates
            load(sprintf('%s\\KDE_binsTaskonly\\%s_%s_iFR50_behavOnly.mat',pat,fname,Areas{iArea}),'iFR','Tmtx');
%             data_{1}.iFR = zscore(iFR);
            data_{1}.iFR = iFR;
            
            % Parse delays
            Delays_ = fieldnames(t);
            nDelays(iTarget) = length(Delays_);
            Delays__{iTarget} = Delays_;
            
            for iDelay = 1:length(Delays_)
                
                fprintf('Analysing run %d/%d %s (%s)...%s\n',iFile,length(fileList),fname,Areas{iArea},Delays_{iDelay})
                
                for iDir=1:2
                    
                    %%%%%%%%%%%%%%%% Correct trials
                    SP = eval(sprintf('t.%s.SamplePress_%sCorrect',Delays_{iDelay},Dirs{iDir}));
                    CP = eval(sprintf('t.%s.ChoicePress_%sCorrect',Delays_{iDelay},Dirs{iDir}));
                    
                    [~,cutout{1}]=SelTimesFRs(tlimsSample,Tmtx,data_,SP*1e-6);
                    [~,cutout{2}]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP*1e-6);
                    
                    cutout{1}{1}(isempty_cell(cutout{1}{1}))=[];
                    cutout{2}{1}(isempty_cell(cutout{2}{1}))=[];
                    
                    cutout_ = cell(size(iFR,2),1);
                    for iUnit=1:size(iFR,2)
                        for iTrial = 1:length(cutout{2}{1})
                            cutout_{iUnit}(iTrial,1:length_) = [cutout{1}{1}{iTrial}(:,iUnit);cutout{2}{1}{iTrial}(:,iUnit)];
                        end
                    end
                    
                    FiringRates{iTarget}{iArea}{iDelay}{iFile}{iDir}     = cutout_;
                    FiringRatesMean{iTarget}{iArea}{iDelay}{iFile}{iDir} = cell2mat(cellfun(@(x) nanmean(x,1),cutout_,'UniformOutput',false))';
                    FiringRatesSEM{iTarget}{iArea}{iDelay}{iFile}{iDir}  = cell2mat(cellfun(@(x) nansem(x,1),cutout_,'UniformOutput',false))';
                    
                    if size(cutout_{1},1)>=FanoMinTrialCount
                        FanoF{iTarget}{iArea}{iDelay}{iFile}{iDir}  = cell2mat(cellfun(@(x) nanvar(x,1),cutout_,'UniformOutput',false))' ./ ...
                                                                      cell2mat(cellfun(@(x) nanmean(x,1),cutout_,'UniformOutput',false))';
                    else
                        FanoF{iTarget}{iArea}{iDelay}{iFile}{iDir}  = nan(length(tb),size(iFR,2));
                    end
                    
                    
                    %%%%%%%%%%%%%%%% Error trials
                    SPe = eval(sprintf('t.%s.SamplePress_%sError',Delays_{iDelay},Dirs{iDir}))';
                    CPe = eval(sprintf('t.%s.ChoicePress_%sError',Delays_{iDelay},Dirs{iDir}))';
                    
                    [~,cutout{1}]=SelTimesFRs(tlimsSample,Tmtx,data_,SPe*1e-6);
                    [~,cutout{2}]=SelTimesFRs(tlimsChoice,Tmtx,data_,CPe*1e-6);
                    
                    cutout{1}{1}(isempty_cell(cutout{1}{1}))=[];
                    cutout{2}{1}(isempty_cell(cutout{2}{1}))=[];

                    cutout_ = cell(size(iFR,2),1);
                    for iUnit=1:size(iFR,2)
                        for iTrial = 1:length(cutout{2}{1})
                            try
                                cutout_{iUnit}(iTrial,1:length_) = [cutout{1}{1}{iTrial}(:,iUnit);cutout{2}{1}{iTrial}(:,iUnit)];
                            catch
                                cutout_{iUnit}(iTrial,1:length_) = NaN;
                            end
                        end
                    end
                    if ~isempty(cutout_{1})
                        FiringRatesError{iTarget}{iArea}{iDelay}{iFile}{iDir}         = cutout_;
                        if size(cutout_{1},1)>0
                            FiringRatesErrorMean{iTarget}{iArea}{iDelay}{iFile}{iDir} = cell2mat(cellfun(@(x) nanmean(x,1),cutout_,'UniformOutput',false))';
                            FiringRatesErrorSEM{iTarget}{iArea}{iDelay}{iFile}{iDir}  = cell2mat(cellfun(@(x) nansem(x,1),cutout_,'UniformOutput',false))';
                            
                            if size(cutout_{1},1)>=FanoMinTrialCount
                                FanoFError{iTarget}{iArea}{iDelay}{iFile}{iDir}       = cell2mat(cellfun(@(x) nanvar(x,1),cutout_,'UniformOutput',false))' ./ ...
                                                                                        cell2mat(cellfun(@(x) nanmean(x,1),cutout_,'UniformOutput',false))';
                            else
                                FanoFError{iTarget}{iArea}{iDelay}{iFile}{iDir}       = nan(length(tb),size(iFR,2));
                            end
                                
                        else
                            FiringRatesErrorMean{iTarget}{iArea}{iDelay}{iFile}{iDir} = cell2mat(cutout_)';
                            FiringRatesErrorSEM{iTarget}{iArea}{iDelay}{iFile}{iDir}  = NaN(size(cell2mat(cutout_)'));
                            FanoFError{iTarget}{iArea}{iDelay}{iFile}{iDir}           = NaN(size(cell2mat(cutout_)'));
                            
                        end
                    end
                    
                end
                %noDirs replaces divide by 2: robust to case where there's
                %L but no R trials and vice versa
                noDirs = sum(~isempty_cell(FiringRatesMean{iTarget}{iArea}{iDelay}{iFile}));
                FiringRatesMeanLR{iTarget}{iArea}{iDelay}{iFile}     = (FiringRatesMean{iTarget}{iArea}{iDelay}{iFile}{1}+ FiringRatesMean{iTarget}{iArea}{iDelay}{iFile}{2})./noDirs;
                FanoFMeanLR{iTarget}{iArea}{iDelay}{iFile}           = (FanoF{iTarget}{iArea}{iDelay}{iFile}{1} + FanoF{iTarget}{iArea}{iDelay}{iFile}{2})./noDirs;
                
                try
                    noDirs = sum(~isempty_cell(FiringRatesErrorMean{iTarget}{iArea}{iDelay}{iFile}));
                    FiringRatesErrorMeanLR{iTarget}{iArea}{iDelay}{iFile} = (FiringRatesErrorMean{iTarget}{iArea}{iDelay}{iFile}{1}+ FiringRatesErrorMean{iTarget}{iArea}{iDelay}{iFile}{2})./noDirs;
                    FanoFMeanErrorLR{iTarget}{iArea}{iDelay}{iFile}       = (FanoFError{iTarget}{iArea}{iDelay}{iFile}{1} + FanoFError{iTarget}{iArea}{iDelay}{iFile}{2})./noDirs;
                catch
                end
                
            end
        end
    end
end


clear iFile fname Delay_ cutout cutout_ ERRORtrangeleft_choice ERRORtrangeleft_sample ERRORtrangeright_choice ERRORtrangeright_sample
clear HPcells HPinter HPtonic iDelay iDir iArea iFile iFR inputData iTrial iUnit PFCcells PFCinter PFCtonic SP t trangeleft_choice trangeleft_sample trangeright_choice trangeright_sample

%% collapse across files

for iTarget = 1:length(Targets)
    for iArea = 1:length(Areas)
        for iDelay = 1:nDelays(iTarget)
            FiringRatesMeanLRcollapse{iTarget}{iArea}{iDelay}       = cell2mat(FiringRatesMeanLR{iTarget}{iArea}{iDelay});
            FanoFMeanLRcollapse{iTarget}{iArea}{iDelay}             = cell2mat(FanoFMeanLR{iTarget}{iArea}{iDelay});
%             FanoFMeanLRcollapse{iTarget}{iArea}{iDelay}(:,nanmax(FanoFMeanLRcollapse{iTarget}{iArea}{iDelay})>10)=NaN;
            
            FiringRatesErrorMeanLRcollapse{iTarget}{iArea}{iDelay}  = cell2mat(FiringRatesErrorMeanLR{iTarget}{iArea}{iDelay});
            FanoFMeanErrorLRcollapse{iTarget}{iArea}{iDelay}        = cell2mat(FanoFMeanErrorLR{iTarget}{iArea}{iDelay});
%             FanoFMeanErrorLRcollapse{iTarget}{iArea}{iDelay}(:,nanmax(FanoFMeanErrorLRcollapse{iTarget}{iArea}{iDelay})>10)=NaN;
        end
    end
end

load('blue_white_red.mat')
color_={'r','b','g'};
color_=[{{'k'}},{num2cell(cool(4),2)},{{'r','b','g'}}];
color_{2}{2} = color_{3}{1};
color_{2}{4} = color_{3}{2};
%% Firing rate maps
for iTarget = 1:length(Targets)
    for iArea = 1:length(Areas)
        figure('color','w')
        for iDelay = 1:nDelays(iTarget)
            if iDelay==1
                [~,i]=max(FiringRatesMeanLRcollapse{iTarget}{iArea}{iDelay});
                [~,i]=sort(i);
                N = max(i);
            end
            subplot(1,nDelays(iTarget),iDelay); hold on
            title(sprintf('%s units: %s delay',Areas{iArea},Delays__{iTarget}{iDelay}))
            imagesc(tb,1:N,FiringRatesMeanLRcollapse{iTarget}{iArea}{iDelay}(:,i)');
            %colormap jet
            colormap(cmap)
            caxis([-1 1])
            
            
            % Vertical bars
            plot(sum(abs(tlimsSample))/2*[1,1],[0 N+1],'color',[0 1 0 0.3],'LineWidth',2)
            plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0 N+1],'color',[1 0 0 0.3],'LineWidth',2)
            plot(sum(abs(tlimsSample))*[1,1],[0 N+1],'k')
            axis([min(tb) max(tb) 0 max(i)])
            
            axis off
            
        end
        plot([17, 19],[10,10],'-k','LineWidth',1.5)
        plot([17, 17],[10,60],'-k','LineWidth',1.5)
    end 
end

%% mean +/- error firing rates
offset = [0.55 0.475 0.4 0.325];
for iTarget = 1:length(Targets)
    
    figure('color','w')
    for iArea = 1:length(Areas)
        subplot(1,2,iArea); hold on

        title(sprintf('%s units: Pop. rate',Areas{iArea}))
        plot([min(tb),max(tb)],[0,0],':k','LineWidth',1,'HandleVisibility','off')
        
        for iDelay = 1:nDelays(iTarget)
            m = nanmean(FiringRatesErrorMeanLRcollapse{iTarget}{iArea}{iDelay},2);
            s = nansem(FiringRatesErrorMeanLRcollapse{iTarget}{iArea}{iDelay},2);
            ciplot(m+s,m-s,tb,color_{iTarget}{iDelay},0.2);
        end
        for iDelay = 1:nDelays(iTarget)
            m = nanmean(FiringRatesMeanLRcollapse{iTarget}{iArea}{iDelay},2);
            s = nansem(FiringRatesMeanLRcollapse{iTarget}{iArea}{iDelay},2);
            ciplot(m+s,m-s,tb,color_{iTarget}{iDelay},1);
        end
        for iDelay = 1:nDelays(iTarget)
            x = FiringRatesMeanLRcollapse{iTarget}{iArea}{iDelay};
            y = FiringRatesErrorMeanLRcollapse{iTarget}{iArea}{iDelay};
            [p,~] = permtest2vec(x,y,5000,0.05);
            a = nan(size(p));a(p)=1;
            %          a(x==y)=NaN;
            plot(tb,a*offset(iDelay),'color',color_{iTarget}{iDelay},'LineWidth',4)
        end
        
        % Vertical bars
        plot(sum(abs(tlimsSample))/2*[1,1],[-0.1 0.1],'color',[0 1 0 0.3],'LineWidth',2)
        plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[-0.1 0.1],'color',[1 0 0 0.3],'LineWidth',2)
        plot(sum(abs(tlimsSample))*[1,1],[-0.1 0.1],'k')
        
        axis off
        
        axis([min(tb) max(tb) -1 1])
        
    end
    plot([17, 19],[-0.6,-0.6],'-k','LineWidth',1.5)
    plot([17, 17],[-0.6 -0.4],'-k','LineWidth',1.5)
    legend(DelaysNames{iTarget});legend boxoff
    

end

%%  mean +/- error Fano factor
offset = [2 2.2 2.4 2.6];
% N.B. rerun the analysis block without zscored rates or this will be junk
for iTarget = 1:length(Targets)

    figure('color','w')
for iArea = 1:length(Areas)
    subplot(1,2,iArea); hold on
    title(sprintf('%s units: Fano factor',Areas{iArea}))
    
    for iDelay = 1:nDelays(iTarget)
        try
        m = nanmean(FanoFMeanErrorLRcollapse{iTarget}{iArea}{iDelay},2);
        s = nansem(FanoFMeanErrorLRcollapse{iTarget}{iArea}{iDelay},2);
        ciplot(m+s,m-s,tb,color_{iTarget}{iDelay},0.2);
        end
    end
    for iDelay = 1:nDelays(iTarget)
        m = nanmean(FanoFMeanLRcollapse{iTarget}{iArea}{iDelay},2);
        s = nansem(FanoFMeanLRcollapse{iTarget}{iArea}{iDelay},2);
        ciplot(m+s,m-s,tb,color_{iTarget}{iDelay},1);
        
    end
    
     for iDelay = 1:nDelays(iTarget)
         try
        x = FiringRatesMeanLRcollapse{iTarget}{iArea}{iDelay};
        y = FiringRatesErrorMeanLRcollapse{iTarget}{iArea}{iDelay};
        [p,~] = permtest2vec(x,y,1000,0.05);
         a = nan(size(p));a(p)=1;
%          a(x==y)=NaN;
         plot(tb,a*offset(iDelay),'color',color_{iTarget}{iDelay},'LineWidth',4)
         end
     end
    
     
    % Vertical bars
    plot(sum(abs(tlimsSample))/2*[1,1],[0.5 1.5],'color',[0 1 0 0.3],'LineWidth',2)
    plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0.5 1.5],'color',[1 0 0 0.3],'LineWidth',2)
    plot(sum(abs(tlimsSample))*[1,1],[0.5 1.5],'k')
    
    %     axis off
    
    axis([min(tb) max(tb) -1 3])
end
    legend(DelaysNames{iTarget});legend boxoff

end

%% Behavioural variable encoding



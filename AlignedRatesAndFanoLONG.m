clear
Targets = {'SHORT','MEDIUM','LONG'};
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Dirs = {'Left','Right'};

bw=0.05;

pat = 'C:\Analysis\AssemblyAnalysis\raw';
cd(pat)
fileList=dir('allTimestamps\*LONG*.mat');
Areas = {'HP','PFC'};

tlimsSample = [-5 5];
tlimsChoice = [-5 5];
length_ = sum(abs([tlimsSample, tlimsChoice]))/bw;
tb = (1:length_)*bw;
%%%%%
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
    load(sprintf('%s\\%s.mat',pat,fname));
    for iArea = 1:length(Areas)
        
        fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
        load(sprintf('%s\\KDE_binsTaskonly\\%s_%s_iFR50_behavOnly.mat',pat,fname,Areas{iArea}),'iFR','Tmtx');
        
        % Spike rates
        data_{1}.iFR = (iFR);
        for iDelay = 1:length(Delays_)
            
            for iDir=1:2
                
                Delay_ = Delays_{iDelay};
                %%%%%%%%%%%%%%%% Correct trials
                SP = eval(sprintf('t.%s.SamplePress_%sCorrect',Delay_,Dirs{iDir}));
                CP = eval(sprintf('t.%s.ChoicePress_%sCorrect',Delay_,Dirs{iDir}));
                
                [~,cutout{1}]=SelTimesFRs(tlimsSample,Tmtx,data_,SP*1e-6);
                [~,cutout{2}]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP*1e-6);
                
                cutout_ = cell(size(iFR,2),1);
                for iUnit=1:size(iFR,2)
                    for iTrial = 1:size(SP,1)
                        cutout_{iUnit}(iTrial,1:length_) = [cutout{1}{1}{iTrial}(:,iUnit);cutout{2}{1}{iTrial}(:,iUnit)];
                    end
                end
                
                FiringRates{iArea}{iDelay}{iFile}{iDir}     = cutout_;
                FiringRatesMean{iArea}{iDelay}{iFile}{iDir} = cell2mat(cellfun(@nanmean,cutout_,'UniformOutput',false))';
                FiringRatesSEM{iArea}{iDelay}{iFile}{iDir}  = cell2mat(cellfun(@nansem,cutout_,'UniformOutput',false))';
                
                FanoF{iArea}{iDelay}{iFile}{iDir}  = cell2mat(cellfun(@nanvar,cutout_,'UniformOutput',false))' ./ cell2mat(cellfun(@nanmean,cutout_,'UniformOutput',false))';

            
                %%%%%%%%%%%%%%%% Error trials
                SPe = eval(sprintf('t.%s.SamplePress_%sError',Delay_,Dirs{iDir}))';
                CPe = eval(sprintf('t.%s.ChoicePress_%sError',Delay_,Dirs{iDir}))';
                
                [~,cutout{1}]=SelTimesFRs(tlimsSample,Tmtx,data_,SPe*1e-6);
                [~,cutout{2}]=SelTimesFRs(tlimsChoice,Tmtx,data_,CPe*1e-6);
                
                cutout_ = cell(size(iFR,2),1);
                for iUnit=1:size(iFR,2)
                    for iTrial = 1:size(SPe,1)
                        try
                            cutout_{iUnit}(iTrial,1:length_) = [cutout{1}{1}{iTrial}(:,iUnit);cutout{2}{1}{iTrial}(:,iUnit)];
                        catch
                            cutout_{iUnit}(iTrial,1:length_) = NaN;
                        end
                    end
                end
                if ~isempty(cutout_{1})
                FiringRatesError{iArea}{iDelay}{iFile}{iDir}         = cutout_;
                if size(cutout_{1},1)>0
                    FiringRatesErrorMean{iArea}{iDelay}{iFile}{iDir} = cell2mat(cellfun(@nanmean,cutout_,'UniformOutput',false))';
                    FiringRatesErrorSEM{iArea}{iDelay}{iFile}{iDir}  = cell2mat(cellfun(@nansem,cutout_,'UniformOutput',false))';
                    FanoFError{iArea}{iDelay}{iFile}{iDir}           = cell2mat(cellfun(@nanvar,cutout_,'UniformOutput',false))' ./ cell2mat(cellfun(@nanmean,cutout_,'UniformOutput',false))';
                else
                    FiringRatesErrorMean{iArea}{iDelay}{iFile}{iDir} = cell2mat(cutout_)';
                    FiringRatesErrorSEM{iArea}{iDelay}{iFile}{iDir}  = NaN(size(cell2mat(cutout_)'));
                    FanoFError{iArea}{iDelay}{iFile}{iDir}           = NaN(size(cell2mat(cutout_)'));

                end
                end

            end
            FiringRatesMeanLR{iArea}{iDelay}{iFile}     = (FiringRatesMean{iArea}{iDelay}{iFile}{1}+ FiringRatesMean{iArea}{iDelay}{iFile}{2})./2;
            FanoFMeanLR{iArea}{iDelay}{iFile}           = (FanoF{iArea}{iDelay}{iFile}{1} + FanoF{iArea}{iDelay}{iFile}{2})./2;

            
            try
            FiringRatesErrorMeanLR{iArea}{iDelay}{iFile} = (FiringRatesErrorMean{iArea}{iDelay}{iFile}{1}+ FiringRatesErrorMean{iArea}{iDelay}{iFile}{2})./2;
            FanoFMeanErrorLR{iArea}{iDelay}{iFile}       = (FanoFError{iArea}{iDelay}{iFile}{1} + FanoFError{iArea}{iDelay}{iFile}{2})./2;
            catch
            end
        end
        
    end
end

clear iFile fname Delay_ cutout cutout_ ERRORtrangeleft_choice ERRORtrangeleft_sample ERRORtrangeright_choice ERRORtrangeright_sample
clear HPcells HPinter HPtonic iDelay iDir iArea iFile iFR inputData iTrial iUnit PFCcells PFCinter PFCtonic SP t trangeleft_choice trangeleft_sample trangeright_choice trangeright_sample

%%
load('blue_white_red.mat')
color_={'r','b','g'};
for iArea = 1:length(Areas)
    for iDelay = 1:length(Delays_)
        FiringRatesMeanLRcollapse{iArea}{iDelay} = cell2mat(FiringRatesMeanLR{iArea}{iDelay});
        FanoFMeanLRcollapse{iArea}{iDelay} = cell2mat(FanoFMeanLR{iArea}{iDelay});
        
        FiringRatesErrorMeanLRcollapse{iArea}{iDelay} = cell2mat(FiringRatesErrorMeanLR{iArea}{iDelay});
        FanoFMeanErrorLRcollapse{iArea}{iDelay} = cell2mat(FanoFMeanErrorLR{iArea}{iDelay});
    end
end
%%

for iArea = 1:length(Areas)
    figure('color','w')
    for iDelay = 1:length(Delays_)
        [~,i]=max(FiringRatesMeanLRcollapse{iArea}{iDelay});
        [~,i]=sort(i);
        N = max(i);
        subplot(1,length(Delays_),iDelay); hold on
        title(sprintf('%s units: %s delay',Areas{iArea},Delays__{iDelay}))
        imagesc(tb,1:N,FiringRatesMeanLRcollapse{iArea}{iDelay}(:,i)');
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
%%
offset = [0.55 0.475 0.4];
figure('color','w')
for iArea = 1:length(Areas)
    subplot(1,2,iArea); hold on
    title(sprintf('%s units: Pop. rate',Areas{iArea}))
    
   
    for iDelay = 1:length(Delays_)
        m = nanmean(FiringRatesErrorMeanLRcollapse{iArea}{iDelay},2);
        s = nansem(FiringRatesErrorMeanLRcollapse{iArea}{iDelay},2);
        ciplot(m+s,m-s,tb,color_{iDelay},0.2);
    end    
    for iDelay = 1:length(Delays_)
        m = nanmean(FiringRatesMeanLRcollapse{iArea}{iDelay},2);
        s = nansem(FiringRatesMeanLRcollapse{iArea}{iDelay},2);
        ciplot(m+s,m-s,tb,color_{iDelay},1);
    end
    for iDelay = 1:length(Delays_)
        x = FiringRatesMeanLRcollapse{iArea}{iDelay};
        y = FiringRatesErrorMeanLRcollapse{iArea}{iDelay};
        [p,~] = permtest2vec(x,y,5000,0.05);
         a = nan(size(p));a(p)=1;
%          a(x==y)=NaN;
         plot(tb,a*offset(iDelay),color_{iDelay},'LineWidth',4)
    end
    
    % Vertical bars
    plot(sum(abs(tlimsSample))/2*[1,1],[-0.1 0.1],'color',[0 1 0 0.3],'LineWidth',2)
    plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[-0.1 0.1],'color',[1 0 0 0.3],'LineWidth',2)
    plot(sum(abs(tlimsSample))*[1,1],[-0.1 0.1],'k')
    axis([min(tb) max(tb) 0 max(i)])
    
    axis off
    
    axis([min(tb) max(tb) -1 1])
    
end
    plot([17, 19],[-0.6,-0.6],'-k','LineWidth',1.5)
    plot([17, 17],[-0.6 -0.4],'-k','LineWidth',1.5)
legend(Delays__);legend boxoff
%%
offset = [2 2.2 2.4];
% N.B. rerun the analysis block without zscored rates or this will be junk
figure('color','w')
for iArea = 1:length(Areas)
    subplot(1,2,iArea); hold on
    title(sprintf('%s units: Fano factor',Areas{iArea}))
    
    for iDelay = 1:length(Delays_)
        m = nanmean(FanoFMeanErrorLRcollapse{iArea}{iDelay},2);
        s = nansem(FanoFMeanErrorLRcollapse{iArea}{iDelay},2);
        ciplot(m+s,m-s,tb,color_{iDelay},0.2);
        
    end
    for iDelay = 1:length(Delays_)
        m = nanmean(FanoFMeanLRcollapse{iArea}{iDelay},2);
        s = nansem(FanoFMeanLRcollapse{iArea}{iDelay},2);
        ciplot(m+s,m-s,tb,color_{iDelay},1);
        
    end
    
     for iDelay = 1:length(Delays_)
        x = FiringRatesMeanLRcollapse{iArea}{iDelay};
        y = FiringRatesErrorMeanLRcollapse{iArea}{iDelay};
        [p,~] = permtest2vec(x,y,1000,0.05);
         a = nan(size(p));a(p)=1;
%          a(x==y)=NaN;
         plot(tb,a*offset(iDelay),color_{iDelay},'LineWidth',4)
     end
    
     
    % Vertical bars
    plot(sum(abs(tlimsSample))/2*[1,1],[0.5 1.5],'color',[0 1 0 0.3],'LineWidth',2)
    plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0.5 1.5],'color',[1 0 0 0.3],'LineWidth',2)
    plot(sum(abs(tlimsSample))*[1,1],[0.5 1.5],'k')
    axis([min(tb) max(tb) 0 max(i)])
    
    %     axis off
    
    axis([min(tb) max(tb) -1 3])
end
legend(Delays__,'Orientation','horizontal' );legend boxoff



%% %%%%%% PREAMBLE %%%%%%
clear 
Delays_ = {'Delay_2','Delay_4','Delay_6','Delay_8',};
Delays__ = {'2s','4s','6s','8s'};
Target = 'MEDIUM';

AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time

tlimsAll = [-5 5];
tlimsDelay_2=[-2 10];
tlimsDelay_4=[-4 10];
tlimsDelay_6=[-6 10];
tlimsDelay_8=[-8 10];
tlimsANOVA = [-2 2];

shift = 0;
plotOnline = false;
runCVEdecoder = true;
bw=0.05;
Nbs = 500;

tbDelay_2=tlimsDelay_2(1):bw:tlimsDelay_2(2);
tbDelay_4=tlimsDelay_4(1):bw:tlimsDelay_4(2);
tbDelay_6=tlimsDelay_6(1):bw:tlimsDelay_6(2);
tbDelay_8=tlimsDelay_8(1):bw:tlimsDelay_8(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
tbANOVA = tlimsANOVA(1):bw:tlimsANOVA(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
pat = 'C:\Analysis\AssemblyAnalysis\raw';
cd(pat)
fileList=dir(sprintf('allTimestamps\\*%s*.mat',Target));
fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileListAss(find(name_flag))=[];% fnames_(name_flag)=[];
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

normaliseFscores = false;
normWin = [-5 -3];
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};

%% %%%%%%%%%%%%%%%%%%%%%%%% ANALYSE UNITS %%%%%%%%%%%%%%%%%%%%%%%%
%% batch import unit data for meta-analysis and plotting
% Import
for iArea = 1:2%length(Areas)
    D_.Delay_2.ANOVA.pThresh{iArea} = [];
    D_.Delay_4.ANOVA.pThresh{iArea} = [];
    D_.Delay_6.ANOVA.pThresh{iArea} = [];
    D_.Delay_8.ANOVA.pThresh{iArea} = [];
    
    for iFile =1:length(fileList)
        try
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity_MEDIUM\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
        load(fnIn ,'D');
        
        for iDelay = 1:length(Delays_)
            
            % L/R: Continuous decoding of Positional information - correct trials
            try
            eval(sprintf('D_temp = D.%s.LR;',Delays_{iDelay}))
            if normaliseFscores
                bp = find(tbAll>normWin(1) & tbAll<normWin(2));
                B  = nanmean(D_temp.Rt2(:,bp)');
                D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
            end
            eval(sprintf('D_.%s.LR.Rt2{iArea}(:,iFile)      = D_temp.Rt2;',Delays_{iDelay}))
            
            if runCVEdecoder
                eval(sprintf('D_.%s.LR.CVE{iArea}(:,iFile)      = D_temp.CVE;',Delays_{iDelay}))
                
                eval(sprintf('D_.%s.LR.CVEsig{iArea}(:,iFile)   = D_temp.cveBSciL>D_temp.CVE;',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
                eval(sprintf('D_.%s.LR.cveBSciH{iArea}(:,iFile) = D_temp.cveBSciH;',Delays_{iDelay}));
                eval(sprintf('D_.%s.LR.cveBSciL{iArea}(:,iFile) = D_temp.cveBSciL;',Delays_{iDelay}));
                eval(sprintf('D_.%s.LR.CVEbs{iArea}(:,iFile)    = D_temp.CVEbs;',Delays_{iDelay}));
            end
            end
            clear D_temp
            
            % L/R: Continuous decoding of Positional information - error trials
            try             
                eval(sprintf('D_temp = D.%s.LR_err;',Delays_{iDelay}))
                if normaliseFscores
                    bp = find(tbAll>normWin(1) & tbAll<-normWin(2));
                    B  = nanmean(D_temp.Rt2(:,bp)');
                    D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
                end
                eval(sprintf('D_.%s.LR_err.Rt2{iArea}(:,iFile)     = D_temp.Rt2;',Delays_{iDelay})),
                if runCVEdecoder
                    eval(sprintf('D_.%s.LR_err.CVE{iArea}(:,iFile)     = D_temp.CVE(1,:);',Delays_{iDelay})),
                    eval(sprintf('D_.%s.LR_err.CVEsig{iArea}(:,iFile)  = D_temp.cveBSciH>D_temp.CVE(1,:);',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
                end
            catch
                eval(sprintf('D_.%s.LR_err.Rt2{iArea}(:,iFile)  = nan(2*length(tbAll),1);',Delays_{iDelay}))
                if runCVEdecoder
                    eval(sprintf('D_.%s.LR_err.CVE{iArea}(:,iFile)  = nan(2*length(tbAll),1);',Delays_{iDelay}))
                    eval(sprintf('D_.%s.LR_err.CVEsig{iArea}(:,iFile)  = nan(2*length(tbAll),1);',Delays_{iDelay}))
                    eval(sprintf('D_.%s.LR_err.CVEbs{iArea}(:,iFile)  = nan(2*length(tbAll),1);',Delays_{iDelay}))
                end
            end
            clear D_temp
            
            % C/E: Continuous decoding of Outcome information
            try
                eval(sprintf('D_tempL = D.%s.CE_L;',Delays_{iDelay}))
                eval(sprintf('D_tempR = D.%s.CE_R;',Delays_{iDelay}))
                if normaliseFscores
                    %eval(sprintf('bp = find(tb%s>normWin(1) & tb%s<(2));',Delays_{iDelay},Delays_{iDelay}));
                    bp = find(tbAll>normWin(1) & tbAll<-normWin(2));

                    B=nanmean(D_tempL.Rt2(:,bp)');
                    D_tempL.Rt2=D_tempL.Rt2./(B'*ones(1,length(D_tempL.Rt2)));  
                    B=nanmean(D_tempR.Rt2(:,bp)');
                    D_tempR.Rt2=D_tempR.Rt2./(B'*ones(1,length(D_tempR.Rt2))); 
                end
                eval(sprintf('D_.%s.CE.Rt2{iArea}(:,iFile)  = mean([D_tempL.Rt2;D_tempR.Rt2]);',Delays_{iDelay}))
            catch
                eval(sprintf('D_.%s.CE.Rt2{iArea}(:,iFile)  = nan(2*length(tbAll),1);',Delays_{iDelay}))
            end
            clear D_tempL D_tempR
            
            % S/C: Continuous decoding of Contextual information - correct trials
            try
                eval(sprintf('D_tempL = D.%s.SC_L;',Delays_{iDelay}))
                eval(sprintf('D_tempR = D.%s.SC_R;',Delays_{iDelay}))
                if normaliseFscores
                    eval(sprintf('bp = find(tb%s>normWin(1) & tb%s<normWin(2));',Delays_{iDelay},Delays_{iDelay}));
                    B=nanmean(D_tempL.Rt2(:,bp)');
                    D_tempL.Rt2=D_tempL.Rt2./(B'*ones(1,length(D_tempL.Rt2)));  
                    B=nanmean(D_tempR.Rt2(:,bp)');
                    D_tempR.Rt2=D_tempR.Rt2./(B'*ones(1,length(D_tempR.Rt2))); 
                end
                eval(sprintf('D_.%s.SC.Rt2{iArea}(:,iFile)  = mean([D_tempL.Rt2;D_tempR.Rt2]);',Delays_{iDelay}))
            catch
                eval(sprintf('D_.%s.SC.Rt2{iArea}(:,iFile)  = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
            end
            clear D_tempL D_tempR
            
            % S/C: Continuous decoding of Contextual information - error trials
            try
                eval(sprintf('D_tempL = D.%s.SCerr_L;',Delays_{iDelay}))
                eval(sprintf('D_tempR = D.%s.SCerr_R;',Delays_{iDelay}))
                if normaliseFscores
                    eval(sprintf('bp = find(tb%s>normWin(1) & tb%s<normWin(2));',Delays_{iDelay},Delays_{iDelay}));
                    B=nanmean(D_tempL.Rt2(:,bp)');
                    D_tempL.Rt2=D_tempL.Rt2./(B'*ones(1,length(D_tempL.Rt2)));  
                    B=nanmean(D_tempR.Rt2(:,bp)');
                    D_tempR.Rt2=D_tempR.Rt2./(B'*ones(1,length(D_tempR.Rt2))); 
                end
                eval(sprintf('D_.%s.SCerr.Rt2{iArea}(:,iFile)  = mean([D_tempL.Rt2;D_tempR.Rt2]);',Delays_{iDelay}))
            catch
                eval(sprintf('D_.%s.SCerr.Rt2{iArea}(:,iFile)  = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
            end
            clear D_tempL D_tempR
            
            % ANOVA: Fractions of tuned cells
            eval(sprintf('D_.%s.ANOVA.prcSig{iArea}(:,iFile)   = D.%s.ANOVA.prcSig;',Delays_{iDelay},Delays_{iDelay}))
            % ANOVA: Fractions of pure and mixed tuned cells
            eval(sprintf('D_.%s.ANOVA.prcTuned{iArea}(:,iFile)  = D.%s.ANOVA.prcTuned;',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.ANOVA.pThresh{iArea}{iFile}     = D.%s.ANOVA.p<0.05;',Delays_{iDelay},Delays_{iDelay}))
            
            % Fano factor: Response selectivity for one condition
            eval(sprintf('D_.%s.Fano.FFR{iArea}{iFile,1}  = D.%s.Fano.FFR;',Delays_{iDelay},Delays_{iDelay}))
            % Fano factor: Trial-to-trial variability in firing rate for each trial type
            eval(sprintf('D_.%s.Fano.FFt{iArea}{iFile}  = D.%s.Fano.FFt;',Delays_{iDelay},Delays_{iDelay}))
            eval(sprintf('D_.%s.Fano.FFtmean{iArea}(:,iFile)  = D.%s.Fano.FFtmean;',Delays_{iDelay},Delays_{iDelay}))
        end
        FFtNames = D.Delay_2.Fano.FFtNames;
        tuning   = D.Delay_2.ANOVA.tuning;
        Factors  = D.Delay_2.ANOVA.Factors;
        clear D
    end
    end
end

% Collapse population means/eerors
for iArea = 1:2%length(Areas)
    %eval(sprintf('',Delays_{iDelay},Delays_{iDelay}))
    for iDelay = 1:length(Delays_)
        eval(sprintf('D_.%s.ANOVA.prcSigMean{iArea}  = nanmean(D_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcSigSEM{iArea}   = nansem(D_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcTunedMean{iArea}  = nanmean(D_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.ANOVA.prcTunedSEM{iArea}  = nansem(D_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFRmean{iArea}  = cellfun(@nanmean,D_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFRsem{iArea}  = cellfun(@nansem,D_.%s.Fano.FFR{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFtmean_popmean{iArea}  = nanmean(D_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
        eval(sprintf('D_.%s.Fano.FFtmean_popSEM{iArea}   = nansem(D_.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}))
    end
    for iFile=1:length(D_.Delay_2.ANOVA.pThresh{iArea})
        Tot_ = D_.Delay_2.ANOVA.pThresh{iArea}{iFile} + D_.Delay_4.ANOVA.pThresh{iArea}{iFile} + D_.Delay_6.ANOVA.pThresh{iArea}{iFile} + D_.Delay_8.ANOVA.pThresh{iArea}{iFile};
        Stability = [sum(Tot_ == 4)./size(Tot_,1)*100;...
                     sum(Tot_ == 3)./size(Tot_,1)*100;...
                     sum(Tot_ == 2)./size(Tot_,1)*100;...
                     sum(Tot_ == 1)./size(Tot_,1)*100];
        Tuningstability_{iArea}(iFile,:) = sum(Stability(1:3,:));
%         Tuningstability_{iArea}(iFile,:) = Stability(1,:);

    end
end

clear D Stability D_temp D_tempL D_tempR B bp 
%% batch import unit data for meta-analysis and plotting - Decoders only
% Import
clear D_
for iArea = 1:2%length(Areas)
 
    for iFile =1:length(fileList)
        
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
        load(fnIn ,'D');
        
        for iDelay = 1:length(Delays_)
            
            % L/R: Continuous decoding of Positional information - correct trials
            eval(sprintf('D_temp = D.%s.LR;',Delays_{iDelay}))
            if normaliseFscores
                bp = find(tbAll>normWin(1) & tbAll<normWin(2));
                B  = nanmean(D_temp.Rt2(:,bp)');
                D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
            end
            eval(sprintf('D_.%s.LR.Rt2{iArea}(:,iFile)      = D_temp.Rt2;',Delays_{iDelay}))
            
            eval(sprintf('D_.%s.LR.CVE{iArea}(:,iFile)      = D_temp.CVE;',Delays_{iDelay}))
            
            eval(sprintf('D_.%s.LR.CVEsig{iArea}(:,iFile)   = nanmean(D_temp.CVE,1)<D_temp.cveBSciL;',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
            eval(sprintf('D_.%s.LR.cveBSciH{iArea}(:,iFile) = D_temp.cveBSciH;',Delays_{iDelay})); 
            eval(sprintf('D_.%s.LR.cveBSciL{iArea}(:,iFile) = D_temp.cveBSciL;',Delays_{iDelay})); 
            eval(sprintf('D_.%s.LR.CVEbs{iArea}(:,iFile)    = D_temp.CVEbs;',Delays_{iDelay})); 

            clear D_temp
            
            % L/R: Continuous decoding of Positional information - error trials
            try             
                eval(sprintf('D_temp = D.%s.LR_err;',Delays_{iDelay}))
                if normaliseFscores
                    bp = find(tbAll>normWin(1) & tbAll<-normWin(2));
                    B  = nanmean(D_temp.Rt2(:,bp)');
                    D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
                end
                eval(sprintf('D_.%s.LR_err.Rt2{iArea}(:,iFile)     = D_temp.Rt2;',Delays_{iDelay})),
                if size(D_temp.CVE,2) ==2*length(tbAll)
                    eval(sprintf('D_.%s.LR_err.CVE{iArea}(:,iFile)     = nanmean(D_temp.CVE,1);',Delays_{iDelay})),
                else
                    eval(sprintf('D_.%s.LR_err.CVE{iArea}(:,iFile)     = nan(2*length(tbAll),1);',Delays_{iDelay}))
                end

                eval(sprintf('D_.%s.LR_err.CVEsig{iArea}(:,iFile)  = D_temp.cveBSciH>nanmean(D_temp.CVE,1);',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
            catch
                eval(sprintf('D_.%s.LR_err.Rt2{iArea}(:,iFile)     = nan(2*length(tbAll),1);',Delays_{iDelay}))
                eval(sprintf('D_.%s.LR_err.CVE{iArea}(:,iFile)     = nan(2*length(tbAll),1);',Delays_{iDelay}))
                eval(sprintf('D_.%s.LR_err.CVEsig{iArea}(:,iFile)  = nan(2*length(tbAll),1);',Delays_{iDelay}))
                eval(sprintf('D_.%s.LR_err.CVEbs{iArea}(:,iFile)   = nan(2*length(tbAll),1);',Delays_{iDelay}))
            end
            clear D_temp
          
            
        end
       
        clear D
    end
    
end


clear D Stability D_temp D_tempL D_tempR B bp 
%% Look for simlarity of tuning across delay length condition

for  iArea = 1:2%length(Areas)
    figure; hold on
   
    points = Tuningstability_{iArea};
    x =repmat(1:length(Factors),size(points,1),1);
    x = x+0.05*randn(size(x));
    scatter(x(:),points(:),5,'Marker','o','MarkerEdgeColor',color_{iArea})

        
    bar(nanmean(Tuningstability_{iArea}),'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
    errorbar(nanmean(Tuningstability_{iArea}),...
             nansem(Tuningstability_{iArea}),'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})   
    
    set(gca,'XTick',1:length(Factors),'XTickLabel',Factors,'XTickLabelRotation',50)
    axis([0 length(Factors)+1 0 100])
    ylabel('% of Units')
    title(sprintf('%s units: Variables encoded for >2 delay lengths',Areas{iArea}))

end
%% Plot Fano FFT for sample/choice - Units

for iArea = 1:2%length(Areas)
    figure('name',['Fano: ' Areas{iArea}]);
    for iDelay=1:length(Delays_)
        subplot(1,length(Delays_),iDelay); hold on
        eval(sprintf('points = D_.%s.Fano.FFtmean{iArea};',Delays_{iDelay}));
%          eval(sprintf('points = cell2mat(D_.%s.Fano.FFt{iArea}'')'';',Delays_{iDelay}));
        bar([1,2],nanmean(points(1:2,:),2),'EdgeColor','k','FaceAlpha',0,'LineWidth',1.5)
        plot([1 2],points(1:2,:),':k')
%         plot([1 2],points(3:4,:),':r')
        errorbar([1,2],nanmean(points(1:2,:),2),nansem(points(1:2,:),2),'k','LineWidth',1.5);
%         errorbar([1,2],nanmean(points(3:4,:),2),nansem(points(3:4,:),2),'r','LineWidth',1.5);
        bar([1,2],nanmean(points(1:2,:),2),'EdgeColor','k','FaceAlpha',0,'LineWidth',1.5)
        set(gca,'XTick',[1,2],'XTickLabel',{'Sample','Choice'},'XTickLabelRotation',50)
        title(Delays_{iDelay})
        axis([0 3 0 1])
        if iDelay==1
           ylabel('Fano Factor') 
        end
    end
end
%% Plot Fano FFR for sample/choice - Units
figure('name','Fano: Tuning');
    for iArea = 1:length(Areas)
        subplot(1,length(Areas),iArea); hold on
        points = [D_.Delay_2.Fano.FFRmean{iArea},...
                  D_.Delay_4.Fano.FFRmean{iArea},...
                  D_.Delay_6.Fano.FFRmean{iArea},...
                  D_.Delay_8.Fano.FFRmean{iArea}];
        bar(1:3,nanmean(points),'EdgeColor','k','FaceAlpha',0,'LineWidth',1.5)
        plot(1:3,points(),':k')
        errorbar(1:3,nanmean(points),nansem(points),'k','LineWidth',1.5);

        set(gca,'XTick',1:3,'XTickLabel',{'Short','Medium','Long'},'XTickLabelRotation',50)
        title(Areas{iArea})
        axis([0 4 0 1])
       if iArea ==1
           ylabel('Response variability (Fano Factor)') 
       end
    end
    
%% Plot proportions of tuned/untuned Units - group by delay
for iArea = 1:2%length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    
    subplot(1,4,1);hold on
    bar(D_.Delay_2.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_.Delay_2.ANOVA.prcTunedMean{iArea},D_.Delay_2.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    ylabel('% of units')
    title('2s delay trials')
    axis([0 6 0 70])
    
    subplot(1,4,2);hold on
    bar(D_.Delay_4.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_.Delay_4.ANOVA.prcTunedMean{iArea},D_.Delay_4.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    title('4s delay trials')
    axis([0 6 0 70])

    subplot(1,4,3);hold on
    bar(D_.Delay_6.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_.Delay_6.ANOVA.prcTunedMean{iArea},D_.Delay_6.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    title('6s delay trials')
    axis([0 6 0 70])
    
    subplot(1,4,4);hold on
    bar(D_.Delay_8.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_.Delay_8.ANOVA.prcTunedMean{iArea},D_.Delay_8.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    title('8s delay trials')
    axis([0 6 0 70])
end
%% Plot proportions of tuned/untuned Units - by delay

for iArea = 1:2%length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]);
    for iType = 1:length(tuning)
        subplot(1,length(tuning),iType);hold on
        m   = [D_.Delay_2.ANOVA.prcTunedMean{iArea}(iType)
            D_.Delay_4.ANOVA.prcTunedMean{iArea}(iType)
            D_.Delay_6.ANOVA.prcTunedMean{iArea}(iType)
            D_.Delay_8.ANOVA.prcTunedMean{iArea}(iType)];
        
        e = [D_.Delay_2.ANOVA.prcTunedSEM{iArea}(iType)
            D_.Delay_4.ANOVA.prcTunedSEM{iArea}(iType)
            D_.Delay_6.ANOVA.prcTunedSEM{iArea}(iType)
            D_.Delay_8.ANOVA.prcTunedSEM{iArea}(iType)];
        
        points = [D_.Delay_2.ANOVA.prcTuned{iArea}(iType,:),...
                  D_.Delay_4.ANOVA.prcTuned{iArea}(iType,:),...
                  D_.Delay_6.ANOVA.prcTuned{iArea}(iType,:),...
                  D_.Delay_8.ANOVA.prcTuned{iArea}(iType,:)];
              
        x = [ones(1,length(D_.Delay_2.ANOVA.prcTuned{iArea}(iType,:))),...
                  2*ones(1,length(D_.Delay_4.ANOVA.prcTuned{iArea}(iType,:))),...
                  3*ones(1,length(D_.Delay_6.ANOVA.prcTuned{iArea}(iType,:))),...
                  4*ones(1,length(D_.Delay_8.ANOVA.prcTuned{iArea}(iType,:)))];
              
        scatter(x+0.05*randn(1,length(x)),points',5,color_{iArea})
        
        b = bar(m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
%         b.LineStyle='none';b.EdgeColor=[0 0 1];b.EdgeColor=[0 0 1];b.FaceAlpha=0.1;
        errorbar(m,e,'Marker','none','color',color_{iArea},'Linewidth',1.2)
        axis([0 4 0 100]);
        set(gca,'XTick',1:length(Delays_),'XTicklabel',Delays__,'XTickLabelRotation',50)
        title(tuning{iType})
        if iType ==1
            ylabel('% of units')
        end

    end
    
end
%% Plot proportions of tuned/untuned Units - by tuning type
for iArea = 1:2%length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    
    for iDelay=1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('b=bar(D_.%s.ANOVA.prcSigMean{iArea},''EdgeColor'',color_{iArea},''EdgeAlpha'',1,''FaceAlpha'',0.2,''FaceColor'',color_{iArea});',Delays_{iDelay}))
        b.LineWidth=1.5;
        eval(sprintf('e=errorbar(D_.%s.ANOVA.prcSigMean{iArea},D_.%s.ANOVA.prcSigSEM{iArea},''color'',color_{iArea},''LineStyle'',''none'');',Delays_{iDelay},Delays_{iDelay}))
        e.LineWidth=1.5;
        set(gca,'XTick',1:length(Factors),'XTicklabel',Factors,'XTickLabelRotation',50)
        title([Delays__{iDelay} ' Delay'])
        axis([0 8 0 50])
        if iDelay==1
            ylabel('% of units')
        end
    
    end
end
%% Plot proportions of tuned/untuned Units - by tuning type AND Delays


for iArea = 1:2%length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    for iType = 1:length(Factors)
        subplot(1,length(Factors),iType);hold on
         m   = [D_.Delay_2.ANOVA.prcSigMean{iArea}(iType)
            D_.Delay_4.ANOVA.prcSigMean{iArea}(iType)
            D_.Delay_6.ANOVA.prcSigMean{iArea}(iType)
            D_.Delay_8.ANOVA.prcSigMean{iArea}(iType)];
          e = [D_.Delay_2.ANOVA.prcSigSEM{iArea}(iType)
            D_.Delay_4.ANOVA.prcSigSEM{iArea}(iType)
            D_.Delay_6.ANOVA.prcSigSEM{iArea}(iType)
            D_.Delay_8.ANOVA.prcSigSEM{iArea}(iType)];
        
        points = [D_.Delay_2.ANOVA.prcSig{iArea}(iType,:),...
                  D_.Delay_4.ANOVA.prcSig{iArea}(iType,:),...
                  D_.Delay_6.ANOVA.prcSig{iArea}(iType,:),...
                  D_.Delay_8.ANOVA.prcSig{iArea}(iType,:)];
              x = [ones(1,length(D_.Delay_2.ANOVA.prcSig{iArea}(iType,:))),...
                  2*ones(1,length(D_.Delay_4.ANOVA.prcSig{iArea}(iType,:))),...
                  3*ones(1,length(D_.Delay_6.ANOVA.prcSig{iArea}(iType,:))),...
                  4*ones(1,length(D_.Delay_8.ANOVA.prcSig{iArea}(iType,:)))];
              
        scatter(x+0.05*randn(1,length(x)),points',5,'Marker','o','MarkerEdgeColor',color_{iArea})
        b = bar(m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
%         b.LineStyle='none';b.EdgeColor=[0 0 1];b.EdgeColor=[0 0 1];b.FaceAlpha=0.1;
        errorbar(m,e,'Marker','none','color',color_{iArea},'Linewidth',1.2)
        axis([0 5 0 40]);
        set(gca,'XTick',1:length(Delays_),'XTicklabel',Delays__,'XTickLabelRotation',50)
        title(Factors{iType})
        if iType ==1
            ylabel('% of units')
        end
    end
end

%% Plot L/R decoding by area
figure('name','L/R decoding');
for iArea =1:2%length(Areas)
    subplot(1,2,iArea);hold on
%     ciplot(mean(D_.Delay_2.LR.Rt2{iArea},2)+nansem(D_.Delay_2.LR.Rt2{iArea},2),...
%            mean(D_.Delay_2.LR.Rt2{iArea},2)-nansem(D_.Delay_2.LR.Rt2{iArea},2),...
%            [0:(length(tbAll)*2-1)]*bw,'r')
%     ciplot(mean(D_.Delay_4.LR.Rt2{iArea},2)+nansem(D_.Delay_4.LR.Rt2{iArea},2),...
%            mean(D_.Delay_4.LR.Rt2{iArea},2)-nansem(D_.Delay_4.LR.Rt2{iArea},2),...
%            [0:(length(tbAll)*2-1)]*bw,'g')    
%     ciplot(mean(D_.Delay_6.LR.Rt2{iArea},2)+nansem(D_.Delay_6.LR.Rt2{iArea},2),...
%            mean(D_.Delay_6.LR.Rt2{iArea},2)-nansem(D_.Delay_6.LR.Rt2{iArea},2),...
%            [0:(length(tbAll)*2-1)]*bw,'b')      
%     ciplot(mean(D_.Delay_8.LR.Rt2{iArea},2)+nansem(D_.Delay_8.LR.Rt2{iArea},2),...
%            mean(D_.Delay_8.LR.Rt2{iArea},2)-nansem(D_.Delay_8.LR.Rt2{iArea},2),...
%            [0:(length(tbAll)*2-1)]*bw,'k')           
       
    plot( [0:(length(tbAll)*2-1)]*bw,mean(D_.Delay_2.LR.Rt2{iArea},2),'r')
    plot( [0:(length(tbAll)*2-1)]*bw,mean(D_.Delay_4.LR.Rt2{iArea},2),'g')
    plot( [0:(length(tbAll)*2-1)]*bw,mean(D_.Delay_6.LR.Rt2{iArea},2),'b')
    plot( [0:(length(tbAll)*2-1)]*bw,mean(D_.Delay_8.LR.Rt2{iArea},2),'k')

    plot([10 10],[0.5 2],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0.1 0.5],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0.1 0.5],'color',[1 0 0 0.3],'LineWidth',1.5)
%     axis([0 20 0 10])
    axis([0 20 0 1])
    legend({'2s delay','4s delay','6s delay','8s delay'})
     if iArea==1
        ylabel('L/R decoding (F-Score)')
    end
end
%% Plot L/R decoding by area - CVE
figure('name','L/R decoding');
alpha_ = 0.6;
alpha2_ = 0.1;
sigcutoff = 0.2;
x_ = [0:(length(tbAll)*2-1)]*bw;
for iArea =1:2%length(Areas)
    subplot(1,2,iArea);hold on
      
    ciplot((mean(1-D_.Delay_2.LR.CVE{iArea},2)+nansem(1-D_.Delay_2.LR.CVE{iArea},2))*100,...
           (mean(1-D_.Delay_2.LR.CVE{iArea},2)-nansem(1-D_.Delay_2.LR.CVE{iArea},2))*100,...
           x_,'r',alpha_)       
    ciplot((mean(1-D_.Delay_4.LR.CVE{iArea},2)+nansem(1-D_.Delay_4.LR.CVE{iArea},2))*100,...
           (mean(1-D_.Delay_4.LR.CVE{iArea},2)-nansem(1-D_.Delay_4.LR.CVE{iArea},2))*100,...
           x_,'g',alpha_)               
    ciplot((mean(1-D_.Delay_6.LR.CVE{iArea},2)+nansem(1-D_.Delay_6.LR.CVE{iArea},2))*100,...
           (mean(1-D_.Delay_6.LR.CVE{iArea},2)-nansem(1-D_.Delay_6.LR.CVE{iArea},2))*100,...
           x_,'b',alpha_)       
    ciplot((mean(1-D_.Delay_8.LR.CVE{iArea},2)+nansem(1-D_.Delay_8.LR.CVE{iArea},2))*100,...
           (mean(1-D_.Delay_8.LR.CVE{iArea},2)-nansem(1-D_.Delay_8.LR.CVE{iArea},2))*100,...
           x_,'k',alpha_)         
       
    ciplot((mean(1-D_.Delay_2.LR.CVEbs{iArea},2)+nansem(1-D_.Delay_2.LR.CVEbs{iArea},2))*100,...
           (mean(1-D_.Delay_2.LR.CVEbs{iArea},2)-nansem(1-D_.Delay_2.LR.CVEbs{iArea},2))*100,...
           x_,'r',alpha2_)
    ciplot((mean(1-D_.Delay_4.LR.CVEbs{iArea},2)+nansem(1-D_.Delay_4.LR.CVEbs{iArea},2))*100,...
           (mean(1-D_.Delay_4.LR.CVEbs{iArea},2)-nansem(1-D_.Delay_4.LR.CVEbs{iArea},2))*100,...
           x_,'g',alpha2_)    
    ciplot((mean(1-D_.Delay_6.LR.CVEbs{iArea},2)+nansem(1-D_.Delay_6.LR.CVEbs{iArea},2))*100,...
           (mean(1-D_.Delay_6.LR.CVEbs{iArea},2)-nansem(1-D_.Delay_6.LR.CVEbs{iArea},2))*100,...
           x_,'b',alpha2_) 
    ciplot((mean(1-D_.Delay_8.LR.CVEbs{iArea},2)+nansem(1-D_.Delay_8.LR.CVEbs{iArea},2))*100,...
           (mean(1-D_.Delay_8.LR.CVEbs{iArea},2)-nansem(1-D_.Delay_8.LR.CVEbs{iArea},2))*100,...
           x_,'k',alpha2_)        
       
%     ciplot(mean(1-D_.Short.LR.cveBSciH{iArea},2)*100,...
%            mean(1-D_.Short.LR.cveBSciL{iArea},2)*100,...
%            x_,'g',alpha2_)       
%     ciplot(mean(1-D_.Medium.LR.cveBSciL{iArea},2)*100,...
%            mean(1-D_.Medium.LR.cveBSciH{iArea},2)*100,...
%            x_,'r',alpha2_)     
%     ciplot(mean(1-D_.Long.LR.cveBSciL{iArea},2)*100,...
%            mean(1-D_.Long.LR.cveBSciH{iArea},2)*100,...
%            x_,'b',alpha2_)     
       
    temp =nansum(D_.Delay_2.LR.CVEsig{iArea},2)./sum(~isnan(sum(D_.Delay_2.LR.CVEsig{iArea},1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-4.5*ones(size(temp,1),1);-3.5*ones(size(temp,1),1)], 100*[temp,temp]')
    rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
    
    temp =nansum(D_.Delay_4.LR.CVEsig{iArea},2)./sum(~isnan(sum(D_.Delay_4.LR.CVEsig{iArea},1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-9.5*ones(size(temp,1),1);-8.5*ones(size(temp,1),1)],100*[temp,temp]')
    rectangle('Position',[0 -10.5 20 3],'LineWidth',1,'EdgeColor','g')

    temp =nansum(D_.Delay_6.LR.CVEsig{iArea},2)./sum(~isnan(sum(D_.Delay_6.LR.CVEsig{iArea},1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-14.5*ones(size(temp,1),1);-13.5*ones(size(temp,1),1)],100*[temp,temp]')
    rectangle('Position',[0 -15.5 20 3],'LineWidth',1,'EdgeColor','b')

    temp =nansum(D_.Delay_8.LR.CVEsig{iArea},2)./sum(~isnan(sum(D_.Delay_8.LR.CVEsig{iArea},1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-14.5*ones(size(temp,1),1);-13.5*ones(size(temp,1),1)],100*[temp,temp]')
    rectangle('Position',[0 -20.5 20 3],'LineWidth',1,'EdgeColor','k')    
   
%     plot(tbShort,mean(D_.Short.LR.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.LR.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.LR.Ft2{iArea},2),'b')
    plot([10 10],[0 100],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 100],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 100],'color',[1 0 0 0.3],'LineWidth',1.5)
    plot([0 20],[50 50],':k')
    axis([0 20 -20 100])
    set(gca,'YTick',[0:25:100])
    
    legend({'2s delay','4s delay','6s delay','8s delay'})
    title(Areas{iArea});
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
     if iArea==1
        ylabel('% Correct decoding')
     end
    xlabel('Time (s)')
    colormap(flipud(gray));caxis([0 100])



end
%% Plot L/R decoding by delay
figure('name','L/R decoding');
for iDelay = 1:length(Delays_)
    subplot(1,length(Delays_),iDelay);hold on
%     title(['L/R trial decoding ' Delays_{iDelay} ' Trials'])
    title([Delays__{iDelay} ' trials'])
    
    for iArea = 1:2%length(Areas)
        eval(sprintf('m=nanmean(D_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e=nansem(D_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,[0:(length(tbAll)*2-1)]*bw,color_{iArea})
    end
    plot([10 10],[0.1 0.5],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0.1 0.5],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0.1 0.5],'color',[1 0 0 0.3],'LineWidth',1.5)
    axis([0 20 0 1])
    if iDelay==1
        ylabel('L/R decoding (F-Score)')
    end
end
legend(Areas{1:2})
%% Plot L/R decoding on correct and errors

for iArea =1:2%length(Areas)
    figure('name',['L/R decoding ',Areas{iArea}]); 
    for iDelay = 1:length(Delays_)
    subplot(1,length(Delays_),iDelay);hold on
    eval(sprintf('m = nanmean(D_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
    eval(sprintf('e = nansem(D_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
    ciplot(m+e,m-e,[0:(length(tbAll)*2-1)]*bw,'k')
    eval(sprintf('m = nanmean(D_.%s.LR_err.Rt2{iArea},2);',Delays_{iDelay}));
    eval(sprintf('e = nansem(D_.%s.LR_err.Rt2{iArea},2);',Delays_{iDelay}));
    ciplot(m+e,m-e,[0:(length(tbAll)*2-1)]*bw,'r')        
    title(Delays_{iDelay})   
    plot([10 10],[0.5 2.5],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 0.5],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 0.5],'color',[1 0 0 0.3],'LineWidth',1.5)
    axis([0 20 0 6])
    if iDelay==1
        ylabel('L/R decoding (F-Score)')
    end
    end
    legend({'Correct','Errors'})
     
end
%% Plot L/R decoding on correct and errors - CVE
alpha_ = 0.6;

for iArea =1:2%length(Areas)
    figure('name',['L/R decoding ',Areas{iArea}]); 
    for iDelay = 1:length(Delays_)
    subplot(1,length(Delays_),iDelay);hold on
    eval(sprintf('m = nanmean(1-D_.%s.LR.CVE{iArea},2)*100;',Delays_{iDelay}));
    eval(sprintf('e = nansem(1-D_.%s.LR.CVE{iArea},2)*100;',Delays_{iDelay}));
    ciplot(m+e,m-e,[0:(length(tbAll)*2-1)]*bw,'k',alpha_)
    eval(sprintf('m = nanmean(1-D_.%s.LR_err.CVE{iArea},2)*100;',Delays_{iDelay}));
    eval(sprintf('e = nansem(1-D_.%s.LR_err.CVE{iArea},2)*100;',Delays_{iDelay}));
    ciplot(m+e,m-e,[0:(length(tbAll)*2-1)]*bw,'r',alpha_)        
    title(Delays__{iDelay})   
    plot([10 10],[0 100],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 100],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 100],'color',[1 0 0 0.3],'LineWidth',1.5)
    plot([0 20],[50 50],':k')

    axis([0 20 0 100])
    if iDelay==1
        ylabel('% Correct decoding')
    end
    xlabel('Time (s)')
    end
    legend({'Correct','Errors'})
     
end
%% Plot S/C decoding

figure('name','S/C decoding');
for iArea =1:2%length(Areas)
    subplot(1,length(Areas),iArea);hold on
    ciplot(mean(D_.Delay_2.SC.Rt2{iArea},2)+nansem(D_.Delay_2.SC.Rt2{iArea},2),...
           mean(D_.Delay_2.SC.Rt2{iArea},2)-nansem(D_.Delay_2.SC.Rt2{iArea},2),...
           tbShort,'r')
    ciplot(mean(D_.Delay_4.SC.Rt2{iArea},2)+nansem(D_.Delay_4.SC.Rt2{iArea},2),...
           mean(D_.Delay_4.SC.Rt2{iArea},2)-nansem(D_.Delay_4.SC.Rt2{iArea},2),...
           tbMedium,'g')    
    ciplot(mean(D_.Delay_6.SC.Rt2{iArea},2)+nansem(D_.Delay_6.SC.Rt2{iArea},2),...
           mean(D_.Delay_6.SC.Rt2{iArea},2)-nansem(D_.Delay_6.SC.Rt2{iArea},2),...
           tbLong,'b')     
    ciplot(mean(D_.Delay_8.SC.Rt2{iArea},2)+nansem(D_.Delay_8.SC.Rt2{iArea},2),...
           mean(D_.Delay_8.SC.Rt2{iArea},2)-nansem(D_.Delay_8.SC.Rt2{iArea},2),...
           tbLong,'k')        
    plot([0 0],[0 0.5],'color',[0 0 0 0.3],'LineWidth',1.5)

           axis([-15 10 0 4])

%     plot(tbShort,mean(D_.Short.SC.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.SC.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.SC.Ft2{iArea},2),'b')
end
legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
%% Plot S/C decoding on correct and errors

for iArea =1:length(Areas)
    figure('name',['S/C decoding ',Areas{iArea}]); 
    for iDelay = 1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('m = nanmean(D_.%s.SC.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e = nansem(D_.%s.SC.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,eval(sprintf('tb%s',Delays_{iDelay})),'k')
        eval(sprintf('m = nanmean(D_.%s.SCerr.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e = nansem(D_.%s.SCerr.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,eval(sprintf('tb%s',Delays_{iDelay})),'r')
        title(Delays_{iDelay})
        plot([0 0],[0 0.5],'color',[0 0 0 0.3],'LineWidth',1.5)

%         axis([-10 10 0 6])

    end
    legend({'Correct','Errors'})
end
%% Plot C/E decoding by area
figure('name','C/E decoding');
for iArea =1:length(Areas)
    subplot(1,length(Areas),iArea);hold on
    title(['C/E trial decoding ' Areas{iArea}])
    ciplot(nanmean(D_.Short.CE.Rt2{iArea},2)+nansem(D_.Short.CE.Rt2{iArea},2),...
           nanmean(D_.Short.CE.Rt2{iArea},2)-nansem(D_.Short.CE.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'r')
    ciplot(nanmean(D_.Medium.CE.Rt2{iArea},2)+nansem(D_.Medium.CE.Rt2{iArea},2),...
           nanmean(D_.Medium.CE.Rt2{iArea},2)-nansem(D_.Medium.CE.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'g')    
    ciplot(nanmean(D_.Long.CE.Rt2{iArea},2)+nansem(D_.Long.CE.Rt2{iArea},2),...
           nanmean(D_.Long.CE.Rt2{iArea},2)-nansem(D_.Long.CE.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'b')           
           legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})

    axis([0 20 0 3])

     plot([10 10],[0 1],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 1],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 1],'color',[1 0 0 0.3],'LineWidth',1.5)
end
%% Plot C/E decoding - by delay
figure('name','C/E decoding');
for iDelay = 1:length(Delays_)
    subplot(1,length(Delays_),iDelay);hold on
    title(['C/E trial decoding ' Delays_{iDelay} ' Trials'])
    
    for iArea = 1:length(Areas)
        eval(sprintf('m=nanmean(D_.%s.CE.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e=nansem(D_.%s.CE.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,...
            [0:(length(tbAll)*2-1)]*bw,color_{iArea})
        
    end
    %     plot(tbShort,mean(D_.Short.CE.Ft2{iArea},2),'r')
    %     plot(tbMedium,mean(D_.Medium.CE.Ft2{iArea},2),'g')
    %     plot(tbLong,mean(D_.Long.CE.Ft2{iArea},2),'b')
    plot([10 10],[0 5],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 1],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 1],'color',[1 0 0 0.3],'LineWidth',1.5)
    axis([0 20 0 3])
end
legend(Areas)

%% %%%%%%%%%%%%%%%%%%% ANALYSE ASSEMBLIES %%%%%%%%%%%%%%%%%%%%%%%%%
%% batch import Assem data for meta-analysis and plotting
for iFile = 1:length(fileListAss)
        
        fname=strtok(fileListAss(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity\\%s_MixedSelectivity_Ass.mat',pat,fname);
        load(fnIn ,'D_Ass');
        
        for iArea = 1:3
            for iDelay = 1:length(Delays_)
                
                % L/R: Continuous decoding of Positional information - correct trials
                flag = 0;
                if eval(sprintf('isfield(D_Ass.%s,''LR'')',Delays_{iDelay})) 
                    if eval(sprintf('~isempty(D_Ass.%s.LR{iArea})',Delays_{iDelay}))
                        flag = 1;
                        eval(sprintf('D_temp = D_Ass.%s.LR{iArea};',Delays_{iDelay}))
                        if normaliseFscores
                            bp = find(tbAll>normWin(1) & tbAll<normWin(2));
                            B  = nanmean(D_temp.Rt2(:,bp)');
                            D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
                        end
                        eval(sprintf(   'D_Ass_.%s.LR.Rt2{iArea}(:,iFile)       = D_temp.Rt2;',Delays_{iDelay}));
                        clear D_temp
                    end
                end
                if ~flag
                	eval(sprintf(   'D_Ass_.%s.LR.Rt2{iArea}(:,iFile)       = nan(2*length(tbAll),1);',Delays_{iDelay}));
                end
                
                % L/R: Continuous decoding of Positional information - error trials
                flag = 0;
                if eval(sprintf('isfield(D_Ass.%s,''LR_err'')',Delays_{iDelay})) 
                    if eval(sprintf('~isempty(D_Ass.%s.LR_err{iArea})',Delays_{iDelay}))
                        flag = 1; 
                        eval(sprintf('D_temp = D_Ass.%s.LR_err{iArea};',Delays_{iDelay}))
                        if normaliseFscores
                            bp = find(tbAll>normWin(1) & tbAll<normWin(2));
                            B  = nanmean(D_temp.Rt2(:,bp)');
                            D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
                        end
                        eval(sprintf(   'D_Ass_.%s.LR_err.Rt2{iArea}(:,iFile)       = D_temp.Rt2;',Delays_{iDelay}));
                        clear D_temp
                    end
                end
                if ~flag
                	eval(sprintf(   'D_Ass_.%s.LR_err.Rt2{iArea}(:,iFile)       = nan(2*length(tbAll),1);',Delays_{iDelay}));
                end
            
                
                % C/E: Continuous decoding of Outcome information
                try
                    eval(sprintf('D_tempL = D_Ass.%s.CE_L{iArea};',Delays_{iDelay}))
                    eval(sprintf('D_tempR = D_Ass.%s.CE_R{iArea};',Delays_{iDelay}))
                    if normaliseFscores
                        %eval(sprintf('bp = find(tb%s>normWin(1) & tb%s<(2));',Delays_{iDelay},Delays_{iDelay}));
                        bp = find(tbAll>normWin(1) & tbAll<-normWin(2));
                        B=nanmean(D_tempL.Rt2(:,bp)');
                        D_tempL.Ft2=D_tempL.Ft2./(B'*ones(1,length(D_tempL.Rt2)));
                        B=nanmean(D_tempR.Rt2(:,bp)');
                        D_tempR.Ft2=D_tempR.Ft2./(B'*ones(1,length(D_tempR.Rt2)));
                    end
                    eval(sprintf(   'D_Ass_.%s.CE.Rt2{iArea}(:,iFile)       = mean([D_tempL.Rt2;D_tempR.Rt2]);',Delays_{iDelay}))
                catch
                    eval(sprintf(   'D_Ass_.%s.CE.Rt2{iArea}(:,iFile)      = nan(2*length(tbAll),1);',Delays_{iDelay}));
                end
                
                % S/C: Continuous decoding of Contextual information - correct trials
                try
                    eval(sprintf('D_tempL = D_Ass.%s.SC_L{iArea};',Delays_{iDelay}))
                    eval(sprintf('D_tempR = D_Ass.%s.SC_R{iArea};',Delays_{iDelay}))
                    if normaliseFscores
                         eval(sprintf('bp = find(tb%s>normWin(1) & tb%s<normWin(2));',Delays_{iDelay},Delays_{iDelay}));
                         B=nanmean(D_tempL.Rt2(:,bp)');
                         D_tempL.Rt2=D_tempL.Rt2./(B'*ones(1,length(D_tempL.Rt2)));
                         B=nanmean(D_tempR.Rt2(:,bp)');
                         D_tempR.Rt2=D_tempR.Rt2./(B'*ones(1,length(D_tempR.Rt2)));
                    end
                    eval(sprintf(   'D_Ass_.%s.SC.Rt2{iArea}(:,iFile)       = mean([D_tempL.Rt2;D_tempR.Rt2]);',Delays_{iDelay}))
                catch
                    eval(sprintf(   'D_Ass_.%s.SC.Rt2{iArea}(:,iFile)       = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                end
                
                % S/C: Continuous decoding of Contextual information - error trials
                try
                    eval(sprintf('D_tempL = D_Ass.%s.SCerr_L{iArea};',Delays_{iDelay}))
                    eval(sprintf('D_tempR = D_Ass.%s.SCerr_R{iArea};',Delays_{iDelay}))
                    if normaliseFscores
                         eval(sprintf('bp = find(tb%s>normWin(1) & tb%s<normWin(2));',Delays_{iDelay},Delays_{iDelay}));
                         B=nanmean(D_tempL.Rt2(:,bp)');
                         D_tempL.Rt2=D_tempL.Rt2./(B'*ones(1,length(D_tempL.Rt2)));
                         B=nanmean(D_tempR.Rt2(:,bp)');
                         D_tempR.Rt2=D_tempR.Rt2./(B'*ones(1,length(D_tempR.Rt2)));
                    end
                    eval(sprintf(   'D_Ass_.%s.SCerr.Rt2{iArea}(:,iFile)       = mean([D_tempL.Rt2;D_tempR.Rt2]);',Delays_{iDelay}))
                catch
                    eval(sprintf(   'D_Ass_.%s.SCerr.Rt2{iArea}(:,iFile)       = nan(length(tb%s),1);',Delays_{iDelay},Delays_{iDelay}))
                end
                
                
                
                        % ANOVA: Fractions of tuned cells
                        eval(sprintf(   'D_Ass_.%s.ANOVA.prcSig{iArea}(:,iFile) = D_Ass.%s.ANOVA.prcSig{iArea};',  Delays_{iDelay},Delays_{iDelay}));
                        
                        % ANOVA: Fractions of pure and mixed tuned cells
                        eval(sprintf(   'D_Ass_.%s.ANOVA.prcTuned{iArea}(:,iFile)   = D_Ass.%s.ANOVA.prcTuned{iArea};',     Delays_{iDelay},Delays_{iDelay}));
                        eval(sprintf(   'D_Ass_.%s.ANOVA.pThresh{iArea}{iFile}      = D_Ass.%s.ANOVA.p{iArea}<0.05;',        Delays_{iDelay},Delays_{iDelay}));
                        
                        % Fano factor: Response selectivity for one condition
                        eval(sprintf(   'D_Ass_.%s.Fano.FFR_Ass{iArea}{iFile}       = abs(D_Ass.%s.Fano.FFR{iArea});',      Delays_{iDelay},Delays_{iDelay}));
                        eval(sprintf(   'D_Ass_.%s.Fano.FFRmean_Ass{iArea}(iFile)   = nanmean(D_Ass.%s.Fano.FFR{iArea});',  Delays_{iDelay},Delays_{iDelay}));
                        
                        % Fano factor: Trial-to-trial variability in firing rate for each trial type
                        eval(sprintf('D_Ass_.%s.Fano.FFtmean_Ass{iArea}{:,iFile}    = abs(D_Ass.%s.Fano.FFtmean{iArea});',  Delays_{iDelay},Delays_{iDelay}));
                    
            end
        end
end

% Collapse population means/errors
for iArea = 1:3
    for iDelay = 1:length(Delays_)

        eval(sprintf(   'D_Ass_.%s.LR.Rt2{iArea}(:,nansum(D_Ass_.%s.LR.Rt2{iArea})==0)=NaN;',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf(   'D_Ass_.%s.CE.Rt2{iArea}(:,nansum(D_Ass_.%s.CE.Rt2{iArea})==0)=NaN;',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf(   'D_Ass_.%s.SC.Rt2{iArea}(:,nansum(D_Ass_.%s.SC.Rt2{iArea})==0)=NaN;',Delays_{iDelay},Delays_{iDelay}));
         
        eval(sprintf(   'D_Ass_.%s.LR.Rt2mean{iArea} = nanmean(D_Ass_.%s.LR.Rt2{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf(   'D_Ass_.%s.CE.Rt2mean{iArea} = nanmean(D_Ass_.%s.CE.Rt2{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf(   'D_Ass_.%s.SC.Rt2mean{iArea} = nanmean(D_Ass_.%s.SC.Rt2{iArea},2);',Delays_{iDelay},Delays_{iDelay}));        
        
        eval(sprintf(   'D_Ass_.%s.LR.Rt2sem{iArea} = nansem(D_Ass_.%s.LR.Rt2{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf(   'D_Ass_.%s.CE.Rt2sem{iArea} = nansem(D_Ass_.%s.CE.Rt2{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
        eval(sprintf(   'D_Ass_.%s.SC.Rt2sem{iArea} = nansem(D_Ass_.%s.SC.Rt2{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
        
         eval(sprintf(   'D_Ass_.%s.ANOVA.prcSigMean{iArea}         = nanmean(D_Ass_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
         eval(sprintf(   'D_Ass_.%s.ANOVA.prcSigSEM{iArea}          = nansem(D_Ass_.%s.ANOVA.prcSig{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
   
    	 eval(sprintf(   'D_Ass_.%s.ANOVA.prcTunedMean{iArea}       = nanmean(D_Ass_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
         eval(sprintf(   'D_Ass_.%s.ANOVA.prcTunedSEM{iArea}        = nansem(D_Ass_.%s.ANOVA.prcTuned{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
         
         eval(sprintf(   'D_Ass_.%s.Fano.FFRmean{iArea}             = cellfun(@nanmean,D_Ass_.%s.Fano.FFR_Ass{iArea});',Delays_{iDelay},Delays_{iDelay}));
         eval(sprintf(   'D_Ass_.%s.Fano.FFRSEM{iArea}              = cellfun(@nansem,D_Ass_.%s.Fano.FFR_Ass{iArea});',Delays_{iDelay},Delays_{iDelay}));
        
         eval(sprintf(   'D_Ass.%s.Fano.FFtmean_AssPopmean{iArea}   = nanmean(D_Ass.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}));
         eval(sprintf(   'D_Ass.%s.Fano.FFtmean_AssPopSEM{iArea}    = nansem(D_Ass.%s.Fano.FFtmean{iArea});',Delays_{iDelay},Delays_{iDelay}));
    end
end
% Stability of tuning across delay lengths
for iArea = 1:3
    for iFile=1:length(D_Ass_.Short.ANOVA.pThresh{iArea})
            Tot_ = D_Ass_.Short.ANOVA.pThresh{iArea}{iFile} + D_Ass_.Medium.ANOVA.pThresh{iArea}{iFile} + D_Ass_.Long.ANOVA.pThresh{iArea}{iFile};
            Stability = [sum(Tot_ == 3)./size(Tot_,1)*100;...
                         sum(Tot_ == 2)./size(Tot_,1)*100;...
                         sum(Tot_ == 1)./size(Tot_,1)*100];
            % Percentage of variables encoded for 2 or more delays 
            Tuningstability_{iArea}(iFile,:) = sum(Stability(1:2,:));
            % Percentage of variables encoded for all delays 
%                     Tuningstability_{iArea}(iFile,:) = Stability(1,:);

    end    
    
end

FFtNames = D_Ass.Long.Fano.FFtNames;
tuning   = D_Ass.Short.ANOVA.tuning;
Factors  = D_Ass.Short.ANOVA.Factors;
%% batch import Assem data for meta-analysis and plotting - decoders only
for iFile = 1:length(fileListAss)
        
        fname=strtok(fileListAss(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity\\%s_MixedSelectivity_Ass.mat',pat,fname);
        load(fnIn ,'D_Ass');
        
        for iArea = 1:3
            for iDelay = 1:length(Delays_)
                
                % L/R: Continuous decoding of Positional information - correct trials
                flag = 0;
                if eval(sprintf('isfield(D_Ass.%s,''LR'')',Delays_{iDelay})) 
                    if eval(sprintf('~isempty(D_Ass.%s.LR{iArea})',Delays_{iDelay}))
                        flag = 1;
                        eval(sprintf('D_temp = D_Ass.%s.LR{iArea};',Delays_{iDelay}))
                        if normaliseFscores
                            bp = find(tbAll>normWin(1) & tbAll<normWin(2));
                            B  = nanmean(D_temp.Rt2(:,bp)');
                            D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
                        end
                        eval(sprintf(   'D_Ass_.%s.LR.Rt2{iArea}(:,iFile)       = D_temp.Rt2;',Delays_{iDelay}));
                        clear D_temp
                    end
                end
                if ~flag
                	eval(sprintf(   'D_Ass_.%s.LR.Rt2{iArea}(:,iFile)       = nan(2*length(tbAll),1);',Delays_{iDelay}));
                end
                
                % L/R: Continuous decoding of Positional information - error trials
                flag = 0;
                if eval(sprintf('isfield(D_Ass.%s,''LR_err'')',Delays_{iDelay})) 
                    if eval(sprintf('~isempty(D_Ass.%s.LR_err{iArea})',Delays_{iDelay}))
                        flag = 1; 
                        eval(sprintf('D_temp = D_Ass.%s.LR_err{iArea};',Delays_{iDelay}))
                        if normaliseFscores
                            bp = find(tbAll>normWin(1) & tbAll<normWin(2));
                            B  = nanmean(D_temp.Rt2(:,bp)');
                            D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
                        end
                        eval(sprintf(   'D_Ass_.%s.LR_err.Rt2{iArea}(:,iFile)       = D_temp.Rt2;',Delays_{iDelay}));
                        clear D_temp
                    end
                end
                if ~flag
                	eval(sprintf(   'D_Ass_.%s.LR_err.Rt2{iArea}(:,iFile)       = nan(2*length(tbAll),1);',Delays_{iDelay}));
                end
                
                % CVE Decoder
                flag = 0;
                if eval(sprintf('isfield(D_Ass.%s,''LR'')',Delays_{iDelay})) 
                    if eval(sprintf('~isempty(D_Ass.%s.LR{iArea})',Delays_{iDelay}))
                        flag = 1;
                        eval(sprintf('D_temp = D_Ass.%s.LR{iArea};',Delays_{iDelay}))
                        if normaliseFscores
                            bp = find(tbAll>normWin(1) & tbAll<normWin(2));
                            B  = nanmean(D_temp.Rt2(:,bp)');
                            D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
                        end
                        eval(sprintf('D_Ass_.%s.LR.Rt2{iArea}(:,iFile)       = D_temp.Rt2;',Delays_{iDelay}));
                        
                         % L/R: Continuous decoding of Positional information - correct trials
                         
                         eval(sprintf('D_Ass_.%s.LR.CVE{iArea}(:,iFile)      = D_temp.CVE;',Delays_{iDelay}))
                         
                        eval(sprintf('D_Ass_.%s.LR.CVEsig{iArea}(:,iFile)   = nanmean(D_temp.CVE,1)<D_temp.cveBSciL;',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
                        eval(sprintf('D_Ass_.%s.LR.cveBSciH{iArea}(:,iFile) = D_temp.cveBSciH;',Delays_{iDelay})); 
                        eval(sprintf('D_Ass_.%s.LR.cveBSciL{iArea}(:,iFile) = D_temp.cveBSciL;',Delays_{iDelay}));
                        eval(sprintf('D_Ass_.%s.LR.CVEbs{iArea}(:,iFile)    = D_temp.CVEbs;',Delays_{iDelay}));

                        
                        
                        
                        clear D_temp
                    end
                end
%                 if ~flag
%                 	eval(sprintf(   'D_Ass_.%s.LR.Rt2{iArea}(:,iFile)       = nan(2*length(tbAll),1);',Delays_{iDelay}));
%                 end        
                
            end
        end
end

% Collapse population means/errors
for iArea = 1:3
    for iDelay = 1:length(Delays_)

        eval(sprintf(   'D_Ass_.%s.LR.Rt2{iArea}(:,nansum(D_Ass_.%s.LR.Rt2{iArea})==0)=NaN;',Delays_{iDelay},Delays_{iDelay}));
       
        eval(sprintf(   'D_Ass_.%s.LR.Rt2mean{iArea} = nanmean(D_Ass_.%s.LR.Rt2{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
     
        
        eval(sprintf(   'D_Ass_.%s.LR.Rt2sem{iArea} = nansem(D_Ass_.%s.LR.Rt2{iArea},2);',Delays_{iDelay},Delays_{iDelay}));
   
        
       
    end
end

%% Look for simlarity of tuning across delay length condition

for  iArea = 1:length(Areas)
    figure; hold on
   
    points = Tuningstability_{iArea};
    x = repmat(1:length(Factors{iArea}),size(points,1),1);
    x = x+0.05*randn(size(x));
    scatter(x(:),points(:),5,'Marker','o','MarkerEdgeColor',color_{iArea})

        
    bar(nanmean(Tuningstability_{iArea}),'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
    errorbar(nanmean(Tuningstability_{iArea}),...
             nansem(Tuningstability_{iArea}),'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})   
    
    set(gca,'XTick',1:length(Factors{iArea}),'XTickLabel',Factors{iArea},'XTickLabelRotation',50)
    axis([0 length(Factors{iArea})+1 0 100])
    ylabel('% of Assemblies')
    title(sprintf('%s assemblies: Variables encoded for >2 delay lengths',Areas{iArea}))

end

%% Plot Fano FFT for sample/choice - Assemblies
% Areas_ = {'HP','PFC','Joint'};
for iArea = 1:length(Areas_)
    figure('name',['Fano: ' Areas_{iArea}]);
    for iDelay=1:length(Delays_)
        subplot(1,length(Delays_),iDelay); hold on
        eval(sprintf('points = D_Ass_.%s.Fano.FFtmean_Ass{iArea};',Delays_{iDelay}));
%          eval(sprintf('points = cell2mat(D_Ass.%s.Fano.FFt{iArea}'')'';',Delays_{iDelay}));
        bar([1,2],nanmean(points(1:2,:),2),'EdgeColor','k','FaceAlpha',0,'LineWidth',1.5)
        plot([1 2],points(1:2,:),':k')
%         plot([1 2],points(3:4,:),':r')
        errorbar([1,2],nanmean(points(1:2,:),2),nansem(points(1:2,:),2),'k','LineWidth',1.5);
%         errorbar([1,2],nanmean(points(3:4,:),2),nansem(points(3:4,:),2),'r','LineWidth',1.5);
        bar([1,2],nanmean(points(1:2,:),2),'EdgeColor','k','FaceAlpha',0,'LineWidth',1.5)
        set(gca,'XTick',[1,2],'XTickLabel',{'Sample','Choice'},'XTickLabelRotation',50)
        title(Delays_{iDelay})
        axis([0 3 0 Inf])
        if iDelay==1
           ylabel('Fano Factor') 
        end
    end
end
%% Plot Fano FFR for sample/choice - Assemblies
figure('name','Fano: Tuning');
% D_Ass.Medium.Fano.FFRmean
    for iArea = 1:length(Areas_)
        subplot(1,length(Areas_),iArea); hold on
        points = [D_Ass.Short.Fano.FFRmean_Ass{iArea}',...
                  D_Ass.Medium.Fano.FFRmean_Ass{iArea}',...
                  D_Ass.Long.Fano.FFRmean_Ass{iArea}'];
        bar(1:3,nanmean(points),'EdgeColor','k','FaceAlpha',0,'LineWidth',1.5)
        plot(1:3,points(),':k')
        errorbar(1:3,nanmean(points),nansem(points),'k','LineWidth',1.5);

        set(gca,'XTick',1:3,'XTickLabel',{'Short','Medium','Long'},'XTickLabelRotation',50)
        title(Areas_{iArea})
        axis([0 4 0 50])
       if iArea ==1
           ylabel('Response variability (Fano Factor)') 
       end
    end
%% Plot proportions of tuned/untuned Assemblies - by delay
color_={'b','r','g'};

for iArea = 1:3%length(Areas)
    figure('name',['Percentage of tuned assemblies: ' Areas{iArea}]);
    for iType = 1:length(tuning{iArea})
        subplot(1,length(tuning{iArea}),iType);hold on
        m   = [D_Ass_.Short.ANOVA.prcTunedMean{iArea}(iType)
            D_Ass_.Medium.ANOVA.prcTunedMean{iArea}(iType)
            D_Ass_.Long.ANOVA.prcTunedMean{iArea}(iType)];
        
        e = [D_Ass_.Short.ANOVA.prcTunedSEM{iArea}(iType)
            D_Ass_.Medium.ANOVA.prcTunedSEM{iArea}(iType)
            D_Ass_.Long.ANOVA.prcTunedSEM{iArea}(iType)];
        
        points = [D_Ass_.Short.ANOVA.prcTuned{iArea}(iType,:),...
                  D_Ass_.Medium.ANOVA.prcTuned{iArea}(iType,:),...
                  D_Ass_.Long.ANOVA.prcTuned{iArea}(iType,:)];
        x = [ones(1,length(D_Ass_.Short.ANOVA.prcTuned{iArea}(iType,:))),...
                  2*ones(1,length(D_Ass_.Medium.ANOVA.prcTuned{iArea}(iType,:))),...
                  3*ones(1,length(D_Ass_.Long.ANOVA.prcTuned{iArea}(iType,:)))];
              
        scatter(x+0.05*randn(1,length(x)),points',5,color_{iArea})
        
        b = bar(m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
%         b.LineStyle='none';b.EdgeColor=[0 0 1];b.EdgeColor=[0 0 1];b.FaceAlpha=0.1;
        errorbar(m,e,'Marker','none','color',color_{iArea},'Linewidth',1.2)
        axis([0 4 0 100]);
        set(gca,'XTick',1:length(Delays_),'XTicklabel',Delays_,'XTickLabelRotation',50)
        title(tuning{iArea}{iType})
        if iType ==1
            ylabel('% of assemblies')
        end

    end
    
end
    
%% Plot proportions of tuned/untuned Assemblies - group by delay
for iArea = 1:length(Areas)
    figure('name',['Percentage of tuned assemblies: ' Areas{iArea}]); 
    
    subplot(1,3,1);hold on
    bar(D_Ass_.Short.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_Ass_.Short.ANOVA.prcTunedMean{iArea},D_Ass_.Short.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    ylabel('% of Assemblies')
    title('Short delay trials')
    axis([0 6 0 100])
    
    subplot(1,3,2);hold on
    bar(D_Ass_.Medium.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_Ass_.Medium.ANOVA.prcTunedMean{iArea},D_Ass_.Medium.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    title('Medium delay trials')
    axis([0 6 0 100])

    subplot(1,3,3);hold on
    bar(D_Ass_.Long.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_Ass_.Long.ANOVA.prcTunedMean{iArea},D_Ass_.Long.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    title('Long delay trials')
    axis([0 6 0 100])
end
%% Plot proportions of tuned/untuned Units - by tuning type
color_={'b','r','g'};
for iArea = 1:3
    figure('name',['Percentage of tuned assemblies: ' Areas{iArea}]); 
    for iDelay=1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('b=bar(D_Ass_.%s.ANOVA.prcSigMean{iArea},''EdgeColor'',color_{iArea},''EdgeAlpha'',1,''FaceAlpha'',0.2,''FaceColor'',color_{iArea});',Delays_{iDelay}))
        b.LineWidth=1.5;
        eval(sprintf('e=errorbar(D_Ass_.%s.ANOVA.prcSigMean{iArea},D_Ass_.%s.ANOVA.prcSigSEM{iArea},''color'',color_{iArea},''LineStyle'',''none'')',Delays_{iDelay},Delays_{iDelay}))
        e.LineWidth=1.5;
        set(gca,'XTick',1:length(Factors{iArea}),'XTicklabel',Factors{iArea},'XTickLabelRotation',50)
        title([Delays_{iDelay} ' Delay'])
        axis([0 8 0 100])
        if iDelay==1
            ylabel('% of Assemblies')
        end
    
    end
end

%% Plot proportions of tuned/untuned Units - by tuning type AND Delays

for iArea = 1%:3%length(Areas)
    figure('name',['Percentage of tuned assemblies: ' Areas{iArea}]); 
    for iType = 1:length(Factors)
        subplot(1,length(Factors),iType);hold on
        m   = [D_Ass_.Short.ANOVA.prcSigMean{iArea}(iType)
            D_Ass_.Medium.ANOVA.prcSigMean{iArea}(iType)
            D_Ass_.Long.ANOVA.prcSigMean{iArea}(iType)];
        e = [D_Ass_.Short.ANOVA.prcSigSEM{iArea}(iType)
            D_Ass_.Medium.ANOVA.prcSigSEM{iArea}(iType)
            D_Ass_.Long.ANOVA.prcSigSEM{iArea}(iType)];
        
        points = [D_Ass_.Short.ANOVA.prcSig{iArea}(iType,:),...
                  D_Ass_.Medium.ANOVA.prcSig{iArea}(iType,:),...
                  D_Ass_.Long.ANOVA.prcSig{iArea}(iType,:)];
              x = [ones(1,length(D_Ass_.Short.ANOVA.prcSig{iArea}(iType,:))),...
                  2*ones(1,length(D_Ass_.Medium.ANOVA.prcSig{iArea}(iType,:))),...
                  3*ones(1,length(D_Ass_.Long.ANOVA.prcSig{iArea}(iType,:)))];
              
        scatter(x+0.05*randn(1,length(x)),points',5,'Marker','o','MarkerEdgeColor',color_{iArea})
        b = bar(m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
%         b.LineStyle='none';b.EdgeColor=[0 0 1];b.EdgeColor=[0 0 1];b.FaceAlpha=0.1;
        errorbar(m,e,'Marker','none','color',color_{iArea},'Linewidth',1.2)
        axis([0.5 3.5 0 100]);
        set(gca,'XTick',1:length(Delays_{iArea}),'XTicklabel',Delays_,'XTickLabelRotation',50)
        title(Factors{iType})
        if iType ==1
            ylabel('% of assemblies')
        end
    end
end

%% Plot L/R decoding by area
figure('name','L/R decoding');
for iArea =1:length(Areas)
    subplot(1,length(Areas),iArea);hold on
    ciplot(nanmean(D_Ass_.Short.LR.Rt2{iArea},2)+nansem(D_Ass_.Short.LR.Rt2{iArea},2),...
           nanmean(D_Ass_.Short.LR.Rt2{iArea},2)-nansem(D_Ass_.Short.LR.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'r')
    ciplot(nanmean(D_Ass_.Medium.LR.Rt2{iArea},2)+nansem(D_Ass_.Medium.LR.Rt2{iArea},2),...
           nanmean(D_Ass_.Medium.LR.Rt2{iArea},2)-nansem(D_Ass_.Medium.LR.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'g')    
    ciplot(nanmean(D_Ass_.Long.LR.Rt2{iArea},2)+nansem(D_Ass_.Long.LR.Rt2{iArea},2),...
           nanmean(D_Ass_.Long.LR.Rt2{iArea},2)-nansem(D_Ass_.Long.LR.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'b')           
       
%     plot(tbShort,mean(D_.Short.LR.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.LR.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.LR.Ft2{iArea},2),'b')
axis([0 20 0 1])
   title(Areas{iArea})
     if iArea==1
        ylabel('L/R decoding (F-Score)')
     end
     plot([10 10],[0 5],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 5],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 5],'color',[1 0 0 0.3],'LineWidth',1.5)
end
 legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
%% Plot L/R decoding by area - CVE
figure('name','L/R decoding');
alpha_ = 0.6;
alpha2_ = 0.1;
sigcutoff = 0.2;
x_ = [0:(length(tbAll)*2-1)]*bw;
for iArea =1:length(Areas)
    subplot(1,3,iArea);hold on
    for iDelay = 1:length(Delays_)  
        eval(sprintf('temp = 1-D_Ass_.%s.LR.CVE{iArea};',Delays_{iDelay}));
        idx{iArea} = sum(temp) == size(temp,1);
%         temp = smooth2a(temp,5,0);
        temp(:,idx{iArea})=[];
        ciplot((mean(temp,2)+nansem(temp,2))*100,...
            (mean(temp,2)-nansem(temp,2))*100,...
            x_,color_{iDelay},alpha_)
    end
    
     for iDelay = 1:length(Delays_)  
        eval(sprintf('temp = 1-D_Ass_.%s.LR.CVEbs{iArea};',Delays_{iDelay}));
%         temp = smooth2a(temp,5,0);
        temp(:,idx{iArea})=[];
        ciplot((mean(temp,2)+nansem(temp,2))*100,...
            (mean(temp,2)-nansem(temp,2))*100,...
            x_,color_{iDelay},alpha2_)
    end
    
       
%     ciplot(mean(1-D_Ass_.Short.LR.cveBSciH{iArea},2)*100,...
%            mean(1-D_Ass_.Short.LR.cveBSciL{iArea},2)*100,...
%            x_,'g',alpha2_)       
%     ciplot(mean(1-D_Ass_.Medium.LR.cveBSciL{iArea},2)*100,...
%            mean(1-D_Ass_.Medium.LR.cveBSciH{iArea},2)*100,...
%            x_,'r',alpha2_)     
%     ciplot(mean(1-D_Ass_.Long.LR.cveBSciL{iArea},2)*100,...
%            mean(1-D_Ass_.Long.LR.cveBSciH{iArea},2)*100,...
%            x_,'b',alpha2_)     

    temp =nansum(D_Ass_.Short.LR.CVEsig{iArea}(:,~idx{iArea}),2)./sum(~isnan(sum(D_Ass_.Short.LR.CVEsig{iArea}(:,~idx{iArea}),1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-4.5*ones(size(temp,1),1);-3.5*ones(size(temp,1),1)], 100*[temp,temp]')
    rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
    
    temp = nansum(D_Ass_.Medium.LR.CVEsig{iArea}(:,~idx{iArea}),2)./sum(~isnan(sum(D_Ass_.Medium.LR.CVEsig{iArea}(:,~idx{iArea}),1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-9.5*ones(size(temp,1),1);-8.5*ones(size(temp,1),1)],100*[temp,temp]')
    rectangle('Position',[0 -10.5 20 3],'LineWidth',1,'EdgeColor','g')

    temp = nansum(D_Ass_.Long.LR.CVEsig{iArea}(:,~idx{iArea}),2)./sum(~isnan(sum(D_Ass_.Long.LR.CVEsig{iArea}(:,~idx{iArea}),1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-14.5*ones(size(temp,1),1);-13.5*ones(size(temp,1),1)],100*[temp,temp]')
    rectangle('Position',[0 -15.5 20 3],'LineWidth',1,'EdgeColor','b')
   
%     plot(tbShort,mean(D_Ass_.Short.LR.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_Ass_.Medium.LR.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_Ass_.Long.LR.Ft2{iArea},2),'b')
    plot([10 10],[0 100],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 100],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 100],'color',[1 0 0 0.3],'LineWidth',1.5)
    plot([0 20],[50 50],':k')
    axis([0 20 -20 100])
    set(gca,'YTick',[0:25:100])
    
    legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'},'Location','northoutside','Orientation','horizontal')
    title(Areas{iArea});
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
     if iArea==1
        ylabel('% Correct decoding')
     end
    xlabel('Time (s)')
    colormap(flipud(gray));caxis([0 100])



end


%% Plot L/R decoding by delay
figure('name','L/R decoding');
for iDelay = 1:length(Delays_)
    subplot(1,length(Delays_),iDelay);hold on
    title(['L/R trial decoding ' Delays_{iDelay} ' Trials'])
    
    for iArea = 1:3
        eval(sprintf('m=nanmean(D_Ass_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e=nansem(D_Ass_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,[0:(length(tbAll)*2-1)]*bw,color_{iArea})
    end
    plot([10 10],[0.5 2.5],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 1],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 1],'color',[1 0 0 0.3],'LineWidth',1.5)
    axis([0 20 0 Inf])
end
legend(Areas)
%% Plot L/R decoding on correct and errors

for iArea =1:3
    figure('name',['L/R decoding ',Areas{iArea}]);
    for iDelay = 1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('m = nanmean(D_Ass_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e = nansem(D_Ass_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,[0:(length(tbAll)*2-1)]*bw,'k')
        eval(sprintf('m = nanmean(D_Ass_.%s.LR_err.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e = nansem(D_Ass_.%s.LR_err.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,[0:(length(tbAll)*2-1)]*bw,'r')
        title(Delays_{iDelay})
        plot([10 10],[0.5 2.5],'color',[0 0 0 0.3],'LineWidth',1.5)
        plot([5 5],[0 0.5],'color',[0 1 0 0.3],'LineWidth',1.5)
        plot([15 15],[0 0.5],'color',[1 0 0 0.3],'LineWidth',1.5)
        axis([0 20 0 20])
        if iDelay==1
            ylabel('L/R decoding (F-Score)')
        end
    end
    legend({'Correct','Errors'})
end
%% Plot S/C decoding

figure('name','S/C decoding');
for iArea =1:3
    subplot(1,length(Areas),iArea);hold on
    title(['Sample/Choice decoding: ' Areas{iArea}])
    ciplot(nanmean(D_Ass_.Short.SC.Rt2{iArea},2)+nansem(D_Ass_.Short.SC.Rt2{iArea},2),...
           nanmean(D_Ass_.Short.SC.Rt2{iArea},2)-nansem(D_Ass_.Short.SC.Rt2{iArea},2),...
           tbShort,'r')
    ciplot(nanmean(D_Ass_.Medium.SC.Rt2{iArea},2)+nansem(D_Ass_.Medium.SC.Rt2{iArea},2),...
           nanmean(D_Ass_.Medium.SC.Rt2{iArea},2)-nansem(D_Ass_.Medium.SC.Rt2{iArea},2),...
           tbMedium,'g')    
    ciplot(nanmean(D_Ass_.Long.SC.Rt2{iArea},2)+nansem(D_Ass_.Long.SC.Rt2{iArea},2),...
           nanmean(D_Ass_.Long.SC.Rt2{iArea},2)-nansem(D_Ass_.Long.SC.Rt2{iArea},2),...
           tbLong,'b')           
    plot([0 0],[0 0.5],'color',[0 0 0 0.3],'LineWidth',1.5)

%            axis([-15 10 0 4])

%     plot(tbShort,mean(D_.Short.SC.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.SC.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.SC.Ft2{iArea},2),'b')
end
legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
%% Plot S/C decoding on correct and errors

for iArea =1:3
    figure('name',['S/C decoding ',Areas{iArea}]); 
    for iDelay = 1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('m = nanmean(D_Ass_.%s.SC.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e = nansem(D_Ass_.%s.SC.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,eval(sprintf('tb%s',Delays_{iDelay})),'k')
        eval(sprintf('m = nanmean(D_Ass_.%s.SCerr.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e = nansem(D_Ass_.%s.SCerr.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,eval(sprintf('tb%s',Delays_{iDelay})),'r')
        title(Delays_{iDelay})
        plot([0 0],[0 0.5],'color',[0 0 0 0.3],'LineWidth',1.5)

%         axis([-10 10 0 6])

    end
    legend({'Correct','Errors'})
end
%% Plot C/E decoding by area
figure('name','C/E decoding');
for iArea =1:length(Areas)
    subplot(1,length(Areas),iArea);hold on
    title(['C/E trial decoding ' Areas{iArea}])
    ciplot(nanmean(D_Ass_.Short.CE.Rt2{iArea},2)+nansem(D_Ass_.Short.CE.Rt2{iArea},2),...
           nanmean(D_Ass_.Short.CE.Rt2{iArea},2)-nansem(D_Ass_.Short.CE.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'r')
    ciplot(nanmean(D_Ass_.Medium.CE.Rt2{iArea},2)+nansem(D_Ass_.Medium.CE.Rt2{iArea},2),...
           nanmean(D_Ass_.Medium.CE.Rt2{iArea},2)-nansem(D_Ass_.Medium.CE.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'g')    
    ciplot(nanmean(D_Ass_.Long.CE.Rt2{iArea},2)+nansem(D_Ass_.Long.CE.Rt2{iArea},2),...
           nanmean(D_Ass_.Long.CE.Rt2{iArea},2)-nansem(D_Ass_.Long.CE.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'b')           
           legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})

    axis([0 20 0 3])

     plot([10 10],[0 1],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 1],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 1],'color',[1 0 0 0.3],'LineWidth',1.5)
end
%% Plot C/E decoding - by delay
figure('name','C/E decoding');
for iDelay = 1:length(Delays_)
    subplot(1,length(Delays_),iDelay);hold on
    title(['C/E trial decoding ' Delays_{iDelay} ' Trials'])
    
    for iArea = 1:3
        eval(sprintf('m=nanmean(D_Ass_.%s.CE.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e=nansem(D_Ass_.%s.CE.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,...
            [0:(length(tbAll)*2-1)]*bw,color_{iArea})
        
    end
    %     plot(tbShort,mean(D_.Short.CE.Ft2{iArea},2),'r')
    %     plot(tbMedium,mean(D_.Medium.CE.Ft2{iArea},2),'g')
    %     plot(tbLong,mean(D_.Long.CE.Ft2{iArea},2),'b')
    plot([10 10],[0 5],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 1],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 1],'color',[1 0 0 0.3],'LineWidth',1.5)
    axis([0 20 0 3])
end
legend(Areas)

%% %%%%%%%%%%%%%%%%%%% ANALYSE ASSEMBLIES %%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the Assembly membership info
% Constrain so that at mothership assembly must have at least one
% significant ANOVA term 
RestrictComparison = true;


switch AssemblyChoice
    case 1
        pat2 = 'C:\Analysis\AssemblyAnalysis\raw\KDE_bins\LONG';
    case 2
        pat2 = 'C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\LONGTaskonly';
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task';
end


for iFile = 1:length(fileListAss)

    fname=strtok(fileListAss(iFile).name,'_');
  
    fprintf('Loading run %d/%d %s ...\n',iFile,length(fileList),fname)
    
    % for assemblies describing  by lever cutouts
    % load(sprintf('%s\\%s_iFR50_Task_FSC.mat',pat2,fname));
    % load(sprintf('%s\\%s_iFR50_AssemRes2.mat',pat2,fname));
    
    
    
switch AssemblyChoice
    case 1
        % for assemblies describing lever cutouts only
        load(sprintf('%s\\%s_iFR50_FSC.mat',pat2,fname),'units','nu');   
        load(sprintf('%s\\%s_iFR50_AssemRes2.mat',pat2,fname),'usel_out');    
    case 2
        % for assemblies describing whole task period (cue  - reward)
        load(sprintf('%s\\%s_iFR50_BehavOnly_FSC.mat',pat2,fname),'units','nu');    
        load(sprintf('%s\\%s_iFR50_BehavOnly_AssemRes2.mat',pat2,fname),'usel_out');    
    case 3
    % for assemblies describing whole task period (continuous)
        load(sprintf('%s\\%s_iFR50_Task_FSC.mat',pat2,fname),'units','nu');    
        %load(sprintf('%s\\%s_iFR50_Task_AssemRes2.mat',pat2,fname),'usel_out');         
        fn = [pat filesep 'KDE_bins' filesep fname '_PFC_iFR50.mat'];
        usel_out = SelCells(fn,0.1,1e6);
end

        
    Ass.usel_out{iFile} = usel_out;
    Ass.nu{iFile} = nu;
    Ass.units{iFile} = units;
    
    clear units usel_out nu
end

% Co-tuning of Assemblies and Single Unit members
for iFile = 1:length(fileListAss)
    fprintf('Analysing run %d/%d %s ...\n',iFile,length(fileList),fname)
    
    for iDelay = 1:length(Delays_)
        for iArea = 1:length(Areas)
            try
                eval(sprintf('p_Ass = D_Ass_.%s.ANOVA.pThresh{iArea}{iFile};',Delays_{iDelay}))
                for iAss = 1:size(p_Ass,1)
                    units_ = Ass.units{iFile}{iArea}{iAss};
                    if iArea < length(Areas)
                        units_ = Ass.usel_out{iFile}{iArea}(units_);
                        eval(sprintf('p_Units = D_.%s.ANOVA.pThresh{iArea}{iFile}(units_,:);',Delays_{iDelay}))
                    else
                        eval(sprintf('L = size(D_.%s.ANOVA.pThresh{1}{iFile},1);',Delays_{iDelay}))
                        %                     usel_out = [Ass.usel_out{iFile}{1},Ass.usel_out{iFile}{2}+ Ass.nu{iFile}(1)];
                        usel_out = [Ass.usel_out{iFile}{1},Ass.usel_out{iFile}{2}+ L];
                        units_   = usel_out(units_);
                        eval(sprintf('pThresh = [D_.%s.ANOVA.pThresh{1}{iFile};D_.%s.ANOVA.pThresh{1}{iFile}];',Delays_{iDelay},Delays_{iDelay}))
                        
                        p_Units = pThresh(units_,:);
                    end
                    
                    if sum(p_Ass(iAss,:))==0 && RestrictComparison
                        Ass.SharedTuning{iFile}{iDelay}{iArea}(iAss,:) = nan(7,1);
                    else
                        Ass.SharedTuning{iFile}{iDelay}{iArea}(iAss,:) = sum(p_Units == p_Ass(iAss,:)) ./ size(p_Units,1);
                    end
                end
            catch
                Ass.SharedTuning{iFile}{iDelay}{iArea} = nan(1,7);
            end
        end
    end
end
for iDelay = 1:length(Delays_)
    for iArea = 1:length(Areas)
        temp = [];
        for iFile = 1:length(fileListAss)
            try
                temp = [temp;Ass.SharedTuning{iFile}{iDelay}{iArea}];
            end
        end
        Ass.SharedTuningCollapsed{iDelay}{iArea} = temp;
    end
end
clear temp fname iArea units_ p_Units p_Ass p_thresh p
 %%

for iArea = 1:length(Areas)
    figure; 
    for iDelay = 1:length(Delays_)


        subplot(1,3,iDelay); hold on
        m = nanmean( Ass.SharedTuningCollapsed{iDelay}{iArea})*100;
        e = nansem( Ass.SharedTuningCollapsed{iDelay}{iArea})*100;
            b = bar(m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
        errorbar(m,e,'Marker','none','color',color_{iArea},'Linestyle','none','LineWidth',1.2)
        title([Areas{iArea} ' (' Delays_{iDelay} ' delay)'])
        axis([0 8 0 100]);
        if iDelay ==1
            ylabel('% Assembly/unit tuning match')
        end
        set(gca,'XTick',1:length(Factors{iArea}),'XTicklabel',Factors{iArea},'XTickLabelRotation',50)

    end
end







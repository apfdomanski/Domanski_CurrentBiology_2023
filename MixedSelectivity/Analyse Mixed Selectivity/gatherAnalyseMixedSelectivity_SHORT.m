%% %%%%%% PREAMBLE %%%%%%
clear 
Delays_ = {'Delay_0'};
Delays__ = {'0s'};
Target = 'SHORT';

tlimsAll = [-5 5];
tlimsDelay_0=[-2 2];
tlimsANOVA = [-2 2];

shift = 0;
plotOnline = false;
runCVEdecoder = true;
bw=0.05;
Nbs = 100;

tbDelay_0=tlimsDelay_0(1):bw:tlimsDelay_0(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));
tbANOVA = tlimsANOVA(1):bw:tlimsANOVA(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];    
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);

fileListAss = fileList;

reject_list={''}; %'ALL_events.mat'
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

AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time
switch AssemblyChoice
    case 1
        pat2 = [pat 'KDE_bins' filesep Target filesep];
    case 2
        pat2 = [pat 'KDE_binsTaskonly' filesep 'LONGTaskonly' filesep];
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end

mkdir([pat 'MixedSelectivity_SHORT'])

%% %%%%%%%%%%%%%%%%%%%%%%%% ANALYSE UNITS %%%%%%%%%%%%%%%%%%%%%%%%
%% batch import unit data for meta-analysis and plotting
% Import
for iArea = 1:2%length(Areas)
    D_.Delay_0.ANOVA.pThresh{iArea} = [];
    
    for iFile =1:length(fileList)
        try
        fname=strtok(fileList(iFile).name,'_');
        fnIn = sprintf('%s\\MixedSelectivity_SHORT\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
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
        FFtNames = D.Delay_0.Fano.FFtNames;
        tuning   = D.Delay_0.ANOVA.tuning;
        Factors  = D.Delay_0.ANOVA.Factors;
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
        title(Delays__{iDelay})
        axis([0 3 0 1])
        if iDelay==1
           ylabel('Fano Factor') 
        end
    end
end
%% Plot Fano FFR for sample/choice - Units
figure('name','Fano: Tuning');
    for iArea = 1:2%length(Areas)
        subplot(1,2,iArea); hold on
        points = [D_.Delay_0.Fano.FFRmean{iArea}];
        bar(1,nanmean(points),'EdgeColor','k','FaceAlpha',0,'LineWidth',1.5)
        plot(1,points(),':k')
        errorbar(1,nanmean(points),nansem(points),'k','LineWidth',1.5);

        set(gca,'XTick',1,'XTickLabel',{'No Delay'},'XTickLabelRotation',50)
        title(Areas{iArea})
        axis([0 2 0 2])
       if iArea ==1
           ylabel('Response variability (Fano Factor)') 
       end
    end
    
%% Plot proportions of tuned/untuned Units - group by delay
for iArea = 1:2%length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    
    hold on
    bar(D_.Delay_0.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_.Delay_0.ANOVA.prcTunedMean{iArea},D_.Delay_0.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    ylabel('% of units')
    title('No delay trials')
    axis([0 6 0 100])
    
   
end
%% Plot proportions of tuned/untuned Units - by delay

for iArea = 1:2%length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]);
    for iType = 1:length(tuning)
        subplot(1,length(tuning),iType);hold on
        m   = [D_.Delay_0.ANOVA.prcTunedMean{iArea}(iType)];
        
        e = [D_.Delay_0.ANOVA.prcTunedSEM{iArea}(iType)];
        
        points = [D_.Delay_0.ANOVA.prcTuned{iArea}(iType,:)];
              
        x = [ones(1,length(D_.Delay_0.ANOVA.prcTuned{iArea}(iType,:)))];
              
        scatter(x+0.05*randn(1,length(x)),points',5,color_{iArea})
        
        b = bar(m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
%         b.LineStyle='none';b.EdgeColor=[0 0 1];b.EdgeColor=[0 0 1];b.FaceAlpha=0.1;
        errorbar(m,e,'Marker','none','color',color_{iArea},'Linewidth',1.2)
        axis([0 2 0 100]);
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
          m = [D_.Delay_0.ANOVA.prcSigMean{iArea}(iType)];
          e = [D_.Delay_0.ANOVA.prcSigSEM{iArea}(iType)];
        
        points = [D_.Delay_0.ANOVA.prcSig{iArea}(iType,:)];
             x = [ones(1,length(D_.Delay_0.ANOVA.prcSig{iArea}(iType,:)))];
              
        scatter(x+0.05*randn(1,length(x)),points',5,'Marker','o','MarkerEdgeColor',color_{iArea})
        b = bar(m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
%         b.LineStyle='none';b.EdgeColor=[0 0 1];b.EdgeColor=[0 0 1];b.FaceAlpha=0.1;
        errorbar(m,e,'Marker','none','color',color_{iArea},'Linewidth',1.2)
        axis([0 2 0 40]);
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
       
    plot( [0:(length(tbAll)*2-1)]*bw,mean(D_.Delay_0.LR.Rt2{iArea},2),'r')

    plot([10 10],[0.5 2],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0.1 0.5],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0.1 0.5],'color',[1 0 0 0.3],'LineWidth',1.5)
%     axis([0 20 0 10])
    axis([0 20 0 1])
    legend({'No delay'})
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
      
    ciplot((mean(1-D_.Delay_0.LR.CVE{iArea},2)+nansem(1-D_.Delay_0.LR.CVE{iArea},2))*100,...
           (mean(1-D_.Delay_0.LR.CVE{iArea},2)-nansem(1-D_.Delay_0.LR.CVE{iArea},2))*100,...
           x_,'r',alpha_)       
    
       
    ciplot((mean(1-D_.Delay_0.LR.CVEbs{iArea},2)+nansem(1-D_.Delay_0.LR.CVEbs{iArea},2))*100,...
           (mean(1-D_.Delay_0.LR.CVEbs{iArea},2)-nansem(1-D_.Delay_0.LR.CVEbs{iArea},2))*100,...
           x_,'r',alpha2_)
   
       
%     ciplot(mean(1-D_.Short.LR.cveBSciH{iArea},2)*100,...
%            mean(1-D_.Short.LR.cveBSciL{iArea},2)*100,...
%            x_,'g',alpha2_)       
%     ciplot(mean(1-D_.Medium.LR.cveBSciL{iArea},2)*100,...
%            mean(1-D_.Medium.LR.cveBSciH{iArea},2)*100,...
%            x_,'r',alpha2_)     
%     ciplot(mean(1-D_.Long.LR.cveBSciL{iArea},2)*100,...
%            mean(1-D_.Long.LR.cveBSciH{iArea},2)*100,...
%            x_,'b',alpha2_)     
       
    temp =nansum(D_.Delay_0.LR.CVEsig{iArea},2)./sum(~isnan(sum(D_.Delay_0.LR.CVEsig{iArea},1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-4.5*ones(size(temp,1),1);-3.5*ones(size(temp,1),1)], 100*[temp,temp]')
    rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
    
    
   
%     plot(tbShort,mean(D_.Short.LR.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.LR.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.LR.Ft2{iArea},2),'b')
    plot([10 10],[0 100],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0 100],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0 100],'color',[1 0 0 0.3],'LineWidth',1.5)
    plot([0 20],[50 50],':k')
    axis([0 20 -20 100])
    set(gca,'YTick',[0:25:100])
    
    legend({'No delay'})
    title(Areas{iArea});

    
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
x_ = [0:(length(tbAll)*2-1)]*bw;

figure('name','S/C decoding');
for iArea =1:2%length(Areas)
    subplot(1,2,iArea);hold on
    ciplot(mean(D_.Delay_0.SC.Rt2{iArea},2)+nansem(D_.Delay_0.SC.Rt2{iArea},2),...
           mean(D_.Delay_0.SC.Rt2{iArea},2)-nansem(D_.Delay_0.SC.Rt2{iArea},2),...
           tbDelay_0,'r')
   
    plot([0 0],[0 0.5],'color',[0 0 0 0.3],'LineWidth',1.5)

           axis([-15 10 0 4])

%     plot(tbShort,mean(D_.Short.SC.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.SC.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.SC.Ft2{iArea},2),'b')
end
legend({'No delay'})
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
    ciplot(nanmean(D_.Delay_0.CE.Rt2{iArea},2)+nansem(D_.Delay_0.CE.Rt2{iArea},2),...
           nanmean(D_.Delay_0.CE.Rt2{iArea},2)-nansem(D_.Delay_0.CE.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'r')
   
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
    
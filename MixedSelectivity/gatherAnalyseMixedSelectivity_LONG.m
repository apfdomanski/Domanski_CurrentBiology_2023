%% %%%%%% PREAMBLE %%%%%%
clear 
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};
Target = 'LONG';

AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time

tlimsAll = [-5 5];
tlimsShort=[-4 10];
tlimsMedium=[-8 10];
tlimsLong=[-16 10];
tlimsANOVA = [-4 4];

shift = 0;
plotOnline = false;
bw=0.05;
Nbs = 500;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbANOVA=tlimsANOVA(1):bw:tlimsANOVA(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
else
%     pat = '/Volumes/HDD2/DNMTP/raw'
    pat = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw/';

end
    cd(pat)
fileList=dir(sprintf('allTimestamps%s*%s*.mat',filesep, Target));
fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
% reject_list={'KrzesimirLONG2_Events.mat','MiroslawLONG1_Events.mat','MiroslawLONG2_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
reject_list={'OnufryLONG1_Events.mat','OnufryLONG2_Events.mat','NorbertLONG1_Events.mat','NorbertLONG2_Events.mat','MiroslawLONG1_Events.mat','MiroslawLONG2_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat',};
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
importCVE = false;
%% %%%%%%%%%%%%%%%%%%%%%%%% ANALYSE UNITS %%%%%%%%%%%%%%%%%%%%%%%%
%% batch import unit data for meta-analysis and plotting
% Import
for iArea = 1:2%length(Areas)
    D_.Short.ANOVA.pThresh{iArea}  = [];
    D_.Medium.ANOVA.pThresh{iArea} = [];
    D_.Long.ANOVA.pThresh{iArea}   = [];
    
    for iFile =1:length(fileList)
            
        fname=strtok(fileList(iFile).name,'_');
%         fnIn = sprintf('%s%sMixedSelectivity_LONG%sLeftVsRight%s%s_%s_MixedSelectivity_Units.mat',pat,filesep,filesep,filesep,fname,Areas{iArea});
        fnIn = sprintf('%s%sMixedSelectivity_LONG%sMixedSelectivity%s%s_%s_MixedSelectivity_Units.mat',pat,filesep,filesep,filesep,fname,Areas{iArea})
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
            if importCVE
                eval(sprintf('D_.%s.LR.CVE{iArea}(:,iFile)      = D_temp.CVE;',Delays_{iDelay}))
                
                eval(sprintf('D_.%s.LR.CVEsig{iArea}(:,iFile)   = D_temp.cveBSciL>D_temp.CVE;',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
                eval(sprintf('D_.%s.LR.cveBSciH{iArea}(:,iFile) = D_temp.cveBSciH;',Delays_{iDelay}));
                eval(sprintf('D_.%s.LR.cveBSciL{iArea}(:,iFile) = D_temp.cveBSciL;',Delays_{iDelay}));
                eval(sprintf('D_.%s.LR.CVEbs{iArea}(:,iFile)    = D_temp.CVEbs;',Delays_{iDelay}));
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
                eval(sprintf('D_.%s.LR_err.Rt2{iArea}(:,iFile)  = D_temp.Rt2;',Delays_{iDelay})),
                if importCVE
                    eval(sprintf('D_.%s.LR_err.CVE{iArea}(:,iFile)  = D_temp.CVE;',Delays_{iDelay})),
                    eval(sprintf('D_.%s.LR_err.CVEsig{iArea}(:,iFile)  = D_temp.cveBSciH>D_temp.CVE;',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
                end
            catch
                eval(sprintf('D_.%s.LR_err.Rt2{iArea}(:,iFile)  = nan(2*length(tbAll),1);',Delays_{iDelay}))
                if importCVE
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
        FFtNames = D.Long.Fano.FFtNames;
        tuning   = D.Short.ANOVA.tuning;
        Factors  = D.Short.ANOVA.Factors;
        clear D
    end
    
end

% Collapse population means/errors
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
    for iFile=1:length(D_.Short.ANOVA.pThresh{iArea})
        Tot_ = D_.Short.ANOVA.pThresh{iArea}{iFile} + D_.Medium.ANOVA.pThresh{iArea}{iFile} + D_.Long.ANOVA.pThresh{iArea}{iFile};
        Stability = [sum(Tot_ == 3)./size(Tot_,1)*100;...
            sum(Tot_ == 2)./size(Tot_,1)*100;...
            sum(Tot_ == 1)./size(Tot_,1)*100];
        Tuningstability_{iArea}(iFile,:) = sum(Stability(1:2,:));
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
%         fnIn = sprintf('%s\\MixedSelectivity\\%s_%s_MixedSelectivity_Units.mat',pat,fname,Areas{iArea});
%         fnIn = sprintf('%s\\MixedSelectivity_%s\\%s_%s_MixedSelectivity_Units.mat',pat,Target,fname,Areas{iArea});

%         fnIn = sprintf('%s\\MixedSelectivity_%s\\LeftvsRight\\%s_%s_MixedSelectivity_Units.mat',pat,Target,fname,Areas{iArea});
        fnIn = sprintf('%s%sMixedSelectivity_%s%sLeftvsRight%s%s_%s_MixedSelectivity_Units.mat',pat,filesep,Target,filesep,filesep,fname,Areas{iArea});

        load(fnIn ,'D');
        
        for iDelay = 1:length(Delays_)
            
            % L/R: Continuous decoding of Positional information - correct trials
            eval(sprintf('D_temp = D.%s.LR;',Delays_{iDelay}))
            if normaliseFscores
                bp = find(tbAll>normWin(1) & tbAll<normWin(2));
                B  = nanmean(D_temp.Rt2(:,bp)');
                D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
                
                B  = nanmean(D_temp.Ft2(:,bp)');
                D_temp.Ft2=D_temp.Ft2./(B'*ones(1,length(D_temp.Ft2)));
                
            end
            

            eval(sprintf('D_.%s.LR.TS{iArea}{iFile}         = D_temp.TS;',Delays_{iDelay}))
            eval(sprintf('D_.%s.LR.avgFR{iArea}(:,iFile)    = D_temp.avgFR;',Delays_{iDelay}))
            eval(sprintf('D_.%s.LR.TSsig{iArea}{iFile}      = D_temp.TSsig;',Delays_{iDelay}))
            eval(sprintf('D_.%s.LR.Rt2{iArea}(:,iFile)      = D_temp.Rt2;',Delays_{iDelay}))
            eval(sprintf('D_.%s.LR.Ft2{iArea}(:,iFile)      = D_temp.Ft2;',Delays_{iDelay}))
            eval(sprintf('D_.%s.LR.Ft2sig{iArea}{iFile}      = [D.%s.LR.Ft2>D.%s.LR.Ft2ciH]'';',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))

            
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
                    
                    B  = nanmean(D_temp.Ft2(:,bp)');
                    D_temp.Ft2=D_temp.Ft2./(B'*ones(1,length(D_temp.Ft2)));
                end
                eval(sprintf('D_.%s.LR_err.TS_err{iArea}{iFile}         = D_temp.TS;',Delays_{iDelay}))
                eval(sprintf('D_.%s.LR_err.avgFR_err{iArea}(:,iFile)    = D_temp.avgFR;',Delays_{iDelay}))
                eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}{iFile}      = D_temp.TSsig;',Delays_{iDelay}))
                eval(sprintf('D_.%s.LR_err.Rt2{iArea}(:,iFile)          = D_temp.Rt2;',Delays_{iDelay})),
                eval(sprintf('D_.%s.LR_err.Ft2{iArea}(:,iFile)          = D_temp.Ft2;',Delays_{iDelay})),
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
D_.Medium.Fano.FFRmean
    for iArea = 1:length(Areas)
        subplot(1,length(Areas),iArea); hold on
        points = [D_.Short.Fano.FFRmean{iArea},...
                  D_.Medium.Fano.FFRmean{iArea},...
                  D_.Long.Fano.FFRmean{iArea}];
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
for iArea = 1:length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    
    subplot(1,3,1);hold on
    bar(D_.Short.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_.Short.ANOVA.prcTunedMean{iArea},D_.Short.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    ylabel('% of units')
    title('Short delay trials')
    axis([0 6 0 70])
    
    subplot(1,3,2);hold on
    bar(D_.Medium.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_.Medium.ANOVA.prcTunedMean{iArea},D_.Medium.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    title('Medium delay trials')
    axis([0 6 0 70])

    subplot(1,3,3);hold on
    bar(D_.Long.ANOVA.prcTunedMean{iArea},'EdgeColor',color_{iArea},'FaceColor',color_{iArea})
    errorbar(D_.Long.ANOVA.prcTunedMean{iArea},D_.Long.ANOVA.prcTunedSEM{iArea},'LineStyle','none','LineWidth',1.5,'Color',color_{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    title('Long delay trials')
    axis([0 6 0 70])
end
%% Plot proportions of tuned/untuned Units - by delay

for iArea = 1:length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]);
    for iType = 1:length(tuning)
        subplot(1,length(tuning),iType);hold on
        m   = [D_.Short.ANOVA.prcTunedMean{iArea}(iType)
            D_.Medium.ANOVA.prcTunedMean{iArea}(iType)
            D_.Long.ANOVA.prcTunedMean{iArea}(iType)];
        
        e = [D_.Short.ANOVA.prcTunedSEM{iArea}(iType)
            D_.Medium.ANOVA.prcTunedSEM{iArea}(iType)
            D_.Long.ANOVA.prcTunedSEM{iArea}(iType)];
        
        points = [D_.Short.ANOVA.prcTuned{iArea}(iType,:),...
                  D_.Medium.ANOVA.prcTuned{iArea}(iType,:),...
                  D_.Long.ANOVA.prcTuned{iArea}(iType,:)];
        x = [ones(1,length(D_.Short.ANOVA.prcTuned{iArea}(iType,:))),...
                  2*ones(1,length(D_.Medium.ANOVA.prcTuned{iArea}(iType,:))),...
                  3*ones(1,length(D_.Long.ANOVA.prcTuned{iArea}(iType,:)))];
              
        scatter(x+0.05*randn(1,length(x)),points',5,color_{iArea})
        
        b = bar(m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
%         b.LineStyle='none';b.EdgeColor=[0 0 1];b.EdgeColor=[0 0 1];b.FaceAlpha=0.1;
        errorbar(m,e,'Marker','none','color',color_{iArea},'Linewidth',1.2)
        axis([0 4 0 100]);
        set(gca,'XTick',1:length(Delays_),'XTicklabel',Delays_,'XTickLabelRotation',50)
        title(tuning{iType})
        if iType ==1
            ylabel('% of units')
        end

    end
    
end
%% Plot proportions of tuned/untuned Units - by tuning type
for iArea = 1:length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    
    for iDelay=1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('b=bar(D_.%s.ANOVA.prcSigMean{iArea},''EdgeColor'',color_{iArea},''EdgeAlpha'',1,''FaceAlpha'',0.2,''FaceColor'',color_{iArea});',Delays_{iDelay}))
        b.LineWidth=1.5;
        eval(sprintf('e=errorbar(D_.%s.ANOVA.prcSigMean{iArea},D_.%s.ANOVA.prcSigSEM{iArea},''color'',color_{iArea},''LineStyle'',''none'')',Delays_{iDelay},Delays_{iDelay}))
        e.LineWidth=1.5;
        set(gca,'XTick',1:length(Factors),'XTicklabel',Factors,'XTickLabelRotation',50)
        title([Delays_{iDelay} ' Delay'])
        axis([0 8 0 100])
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
         m   = [D_.Short.ANOVA.prcSigMean{iArea}(iType)
            D_.Medium.ANOVA.prcSigMean{iArea}(iType)
            D_.Long.ANOVA.prcSigMean{iArea}(iType)];
          e = [D_.Short.ANOVA.prcSigSEM{iArea}(iType)
            D_.Medium.ANOVA.prcSigSEM{iArea}(iType)
            D_.Long.ANOVA.prcSigSEM{iArea}(iType)];
        
        points = [D_.Short.ANOVA.prcSig{iArea}(iType,:),...
                  D_.Medium.ANOVA.prcSig{iArea}(iType,:),...
                  D_.Long.ANOVA.prcSig{iArea}(iType,:)];
              x = [ones(1,length(D_.Short.ANOVA.prcSig{iArea}(iType,:))),...
                  2*ones(1,length(D_.Medium.ANOVA.prcSig{iArea}(iType,:))),...
                  3*ones(1,length(D_.Long.ANOVA.prcSig{iArea}(iType,:)))];
              
        scatter(x+0.05*randn(1,length(x)),points',5,'Marker','o','MarkerEdgeColor',color_{iArea})
        b = bar(m,'EdgeColor',color_{iArea},'FaceColor',color_{iArea},'FaceAlpha',0.2,'LineWidth',1.2);
%         b.LineStyle='none';b.EdgeColor=[0 0 1];b.EdgeColor=[0 0 1];b.FaceAlpha=0.1;
        errorbar(m,e,'Marker','none','color',color_{iArea},'Linewidth',1.2)
        axis([0 4 0 40]);
        set(gca,'XTick',1:length(Delays_),'XTicklabel',Delays_,'XTickLabelRotation',50)
        title(Factors{iType})
        if iType ==1
            ylabel('% of units')
        end
    end
end

%% Get t-scores, plot side by side

for iArea = 1:length(Areas)-1
    figure
    
%     iDelay = 1;
%     TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
%     TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%             [~,idx] = max(TS_,[],2); 
%             [~,idx] = sort(idx);
            
            
%     TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
%     TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%     idxNS = find(nansum(TSsig_,2)==0);
%     [~,idx] = max(TS_,[],2); 
%     [~,idx] = sort(idx);
%     idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];

    for iDelay = 1:length(Delays_)
        TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
        TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%         TS(sum(TSsig_,2)==0,:)=0;
        
        % Sort by peak time
        if iDelay ==1
            [~,idx] = max(TS_,[],2); 
            [~,idx] = sort(idx);
        end
    % Sort by peak time, separate insignificant units
%     if iDelay ==1
%     idxNS = find(nansum(TSsig_,2)==0);
%     [~,idx] = max(TS_,[],2); 
%     [~,idx] = sort(idx);
%     idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];
%     end
    
    TS_ =TS_(idx,:);
    TSsig_ =TSsig_ (idx,:);
%     TS_(~TSsig_ )=0;
    subplot(1,3,iDelay); hold on
    x_  = repmat((1:size(TS_,2)).*bw,size(TS_,1),1);
    y_  = repmat(1:size(TS_,1),size(TS_,2),1);
    imagesc(x_(:),y_(:),TS_);
    set(gca,'YDir','normal')
%     cmap =([1 1 1;flipud(hot)]);
    cmap =([1 1 1;(jet)]);
    cmap =(jet);
    colormap (cmap)
    caxis([0 10])
    plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 0 350])
    if iArea==1 & iDelay ==3
        plot([18 20],[20 20],'k','LineWidth',1.5)
        plot([18 18],[20 70],'k','LineWidth',1.5)
    end
    axis off
    
    end
end
% clear idx TS_ TSsig_
%% Get t-scores, plot side by side for one animal
iFile = 1;
for iArea = 1:length(Areas)-1
    figure
    
%     iDelay = 1;
%     TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
%     TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%             [~,idx] = max(TS_,[],2); 
%             [~,idx] = sort(idx);
            
            
%     TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
%     TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%     idxNS = find(nansum(TSsig_,2)==0);
%     [~,idx] = max(TS_,[],2); 
%     [~,idx] = sort(idx);
%     idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];

    for iDelay = 1:length(Delays_)
        TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}(iFile)',Delays_{iDelay})))';
        TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}(iFile)',Delays_{iDelay})))';
%         TS(sum(TSsig_,2)==0,:)=0;
        
        % Sort by peak time
        if iDelay ==1
            [~,idx] = max(TS_,[],2); 
            [~,idx] = sort(idx);
        end
    % Sort by peak time, separate insignificant units
%     if iDelay ==1
%     idxNS = find(nansum(TSsig_,2)==0);
%     [~,idx] = max(TS_,[],2); 
%     [~,idx] = sort(idx);
%     idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];
%     end
    
    TS_    = TS_(idx,:);
    TSsig_ = TSsig_(idx,:);
%     TS_(~TSsig_ )=0;
    subplot(1,3,iDelay); hold on
    x_  = repmat((1:size(TS_,2)).*bw,size(TS_,1),1);
    y_  = repmat(1:size(TS_,1),size(TS_,2),1);
    imagesc(x_(:),y_(:),TS_);
    set(gca,'YDir','normal')
%     cmap =([1 1 1;flipud(hot)]);
    cmap =([1 1 1;(jet)]);
    cmap =(jet);
    colormap(cmap)
    caxis([0 10])
    plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 0 size(TS_,1)])
    if iArea==1 & iDelay ==3
        plot([18 20],[20 20],'k','LineWidth',1.5)
        plot([18 18],[20 70],'k','LineWidth',1.5)
    end
    axis off
    
    end
end
% clear idx TS_ TSsig_
%% Difference between peak decoding time across delays
clear latency
col_Delay = flipud(gray(3));
figure
for iFile = 1:length(fileList)

for iArea = 1:length(Areas)-1
%     subplot(1,2,iArea); hold on
    for iDelay = 1:length(Delays_)
        TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}(iFile)',Delays_{iDelay})))';
        %         TS_ = cell2mat(eval(sprintf('D_.%s.LR_err.TS_err{iArea}(iFile)',Delays_{iDelay})))';
        
        [~,Peaktime_] = max(TS_,[],2);
        [~,idx] = sort(Peaktime_);
        y = (1:max(idx))./max(idx);
        bins = 1:size(TS_,2);
        for iBin =1:length(bins)
            temp(iBin)=y(min(FindClosestIndex (Peaktime_(idx),bins(iBin))));
        end
        latency{iArea}{iDelay}(iFile,:)=temp;
        
        TS_ = cell2mat(eval(sprintf('D_.%s.LR_err.TS_err{iArea}(iFile)',Delays_{iDelay})))';
        
        [~,Peaktime_] = max(TS_,[],2);
        [~,idx] = sort(Peaktime_);
        y = (1:max(idx))./max(idx);
        bins = 1:size(TS_,2);
        for iBin =1:length(bins)
            temp(iBin)=y(min(FindClosestIndex (Peaktime_(idx),bins(iBin))));
        end
        latencyErr{iArea}{iDelay}(iFile,:)=temp;
    end
end
end

for iArea = 1:length(Areas)-1  
    subplot(1,2,iArea); hold on
    for iDelay = 2:length(Delays_)
      
        m = nanmean(latency{iArea}{iDelay});
        e = nansem(latency{iArea}{iDelay});
%         m = nanmean(latency{iArea}{iDelay}-latency{iArea}{2});
%         e = nansem(latency{iArea}{iDelay}-latency{iArea}{2});
        ciplot(m+e,m-e,bins*bw,col_Delay(iDelay,:),1)
        
        m = nanmean(latencyErr{iArea}{iDelay});
        e = nansem(latencyErr{iArea}{iDelay});
%         m = nanmean(latency{iArea}{iDelay}-latency{iArea}{2});
%         e = nansem(latency{iArea}{iDelay}-latency{iArea}{2});
%         ciplot(m+e,m-e,bins*bw,col_Delay(iDelay,:),1)
        plot(bins*bw,m+e,':','color',col_Delay(iDelay,:),'LineWidth',1.5,'HandleVisibility','off')
        plot(bins*bw,m-e,':','color',col_Delay(iDelay,:),'LineWidth',1.5,'HandleVisibility','off')
        
        legend(Delays_(2:3),'Location','southeast')
    end
    plot(bins*bw,bins./max(bins),':k','LineWidth',1.5,'HandleVisibility','off')
end   
%% Difference between peak firing time across delays
clear latency
col_Delay = flipud(autumn(3));
figure
for iFile = 1:length(fileList)

for iArea = 1:length(Areas)-1
    for iDelay = 1:length(Delays_)
        TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}(iFile)',Delays_{iDelay})))';
        avgFR_  = (cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,iFile)',Delays_{iDelay})))' + ...
            cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,iFile)',Delays_{iDelay})))')./2;
        
        % Sort by peak time
        
        [~,Peaktime_] = max(avgFR_,[],2);
        
        [~,idx] = sort(Peaktime_);
        y = (1:max(idx))./max(idx);
        bins = 1:size(TS_,2);
        for iBin =1:length(bins)
            temp(iBin)=y(min(FindClosestIndex (Peaktime_(idx),bins(iBin))));
        end
        latency{iArea}{iDelay}(iFile,:)=temp;
        
        try
        TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}(iFile)',Delays_{iDelay})))';
              
        avgFR_  = (cell2mat(eval(sprintf('D_.%s.LR_err.avgFR_err{iArea}(1,iFile)',Delays_{iDelay})))' + ...
                   cell2mat(eval(sprintf('D_.%s.LR_err.avgFR_err{iArea}(2,iFile)',Delays_{iDelay})))')./2;               
        
        % Sort by peak time
            [~,Peaktime_] = max(avgFR_,[],2);

            [~,idx] = sort(Peaktime_);
            y = (1:max(idx))./max(idx);
            bins = 1:size(TS_,2);
            for iBin =1:length(bins)
               temp(iBin)=y(min(FindClosestIndex (Peaktime_(idx),bins(iBin))));
            end
            latencyErr{iArea}{iDelay}(iFile,:)=temp;
%         end
        catch
            latencyErr{iArea}{iDelay}(iFile,:)=nan(1,402);
        end
    end
end
end

for iArea = 1:length(Areas)-1  
    subplot(1,2,iArea); hold on
    for iDelay = 2:length(Delays_)
      
        m = nanmean(latency{iArea}{iDelay});
        e = nansem(latency{iArea}{iDelay});
%         m = nanmean(latency{iArea}{iDelay}-latency{iArea}{2});
%         e = nansem(latency{iArea}{iDelay}-latency{iArea}{2});
        ciplot(m+e,m-e,bins*bw,col_Delay(iDelay,:),1)
        
        m = nanmean(latencyErr{iArea}{iDelay});
        e = nansem(latencyErr{iArea}{iDelay});
%         m = nanmean(latency{iArea}{iDelay}-latency{iArea}{2});
%         e = nansem(latency{iArea}{iDelay}-latency{iArea}{2});
%         ciplot(m+e,m-e,bins*bw,col_Delay(iDelay,:),1)
        plot(bins*bw,m+e,':','color',col_Delay(iDelay,:),'LineWidth',1.5,'HandleVisibility','off')
        plot(bins*bw,m-e,':','color',col_Delay(iDelay,:),'LineWidth',1.5,'HandleVisibility','off')
        
        legend(Delays_(2:3),'Location','southeast')
    end
    plot(bins*bw,bins./max(bins),':k','LineWidth',1.5,'HandleVisibility','off')
end   
%% Difference between sort order of peak decoding time across delays
clear Peaktime_ PeakTime PeakTimeErr PeakTimeOut PeakTimeErrOut
col_Delay = flipud(gray(3));
figure
bins =-402:100:402;
for iFile = 1:length(fileList)
    
    for iArea = 1:length(Areas)-1
        %     subplot(1,2,iArea); hold on
        for iDelay = 1:length(Delays_)
            TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}(iFile)',Delays_{iDelay})))';
            [~,Peaktime_] = max(TS_,[],2);
            PeakTime{iArea}{iFile}(iDelay,:) = Peaktime_;
            try
                TS_Err = cell2mat(eval(sprintf('D_.%s.LR_err.TS_err{iArea}(iFile)',Delays_{iDelay})))';
                [~,Peaktime_] = max(TS_Err,[],2);
                PeakTimeErr{iArea}{iFile}(iDelay,:)=Peaktime_;
            catch
                PeakTimeErr{iArea}{iFile}(iDelay,1:size(TS_,1))=nan(size(TS_,1),1);
            end
        end
        PeakTimeErr{iArea}{iFile} = PeakTimeErr{iArea}{iFile}-PeakTimeErr{iArea}{iFile}(2,:);
        PeakTime{iArea}{iFile} = PeakTime{iArea}{iFile}-PeakTime{iArea}{iFile}(2,:);
    end
      for iArea = 1:length(Areas)-1
    for iDelay = 1:length(Delays_)
        PeakTimeErrOut{iArea}{iDelay}(iFile,:) = histc(PeakTimeErr{iArea}{iFile}(iDelay,:),bins)./sum(PeakTimeErr{iArea}{iFile}(iDelay,:));
        PeakTimeOut{iArea}{iDelay}(iFile,:)    = histc(PeakTime{iArea}{iFile}(iDelay,:),bins)./sum(PeakTime{iArea}{iFile}(iDelay,:));
    end
      end
end

for iArea = 1:length(Areas)-1  
    subplot(1,2,iArea); hold on
    for iDelay = 2:length(Delays_)
      
        m = nanmean(PeakTimeOut{iArea}{iDelay});
        e = nansem(PeakTimeOut{iArea}{iDelay});
%         m = nanmean(latency{iArea}{iDelay}-latency{iArea}{2});
%         e = nansem(latency{iArea}{iDelay}-latency{iArea}{2});
        ciplot(m+e,m-e,bins*bw,col_Delay(iDelay,:),1)
        
         m = nanmean(PeakTimeErrOut{iArea}{iDelay});
         e = nansem(PeakTimeErrOut{iArea}{iDelay});
%         m = nanmean(latency{iArea}{iDelay}-latency{iArea}{2});
%         e = nansem(latency{iArea}{iDelay}-latency{iArea}{2});
%         ciplot(m+e,m-e,bins*bw,col_Delay(iDelay,:),1)
        plot(bins*bw,m+e,':','color',col_Delay(iDelay,:),'LineWidth',1.5,'HandleVisibility','off')
        plot(bins*bw,m-e,':','color',col_Delay(iDelay,:),'LineWidth',1.5,'HandleVisibility','off')
        
        legend(Delays_(2:3),'Location','southeast')
    end
end        

%% Get t-scores, plot side by side - Errors

for iArea = 1:length(Areas)-1
    figure
    
%     iDelay = 1;
%     TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
%     TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%             [~,idx] = max(TS_,[],2); 
%             [~,idx] = sort(idx);
            
            
%     TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
%     TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%     idxNS = find(nansum(TSsig_,2)==0);
%     [~,idx] = max(TS_,[],2); 
%     [~,idx] = sort(idx);
%     idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];

    for iDelay = 1:length(Delays_)
        TS_ = cell2mat(eval(sprintf('D_.%s.LR_err.TS_err{iArea}',Delays_{iDelay})))';
        TSsig_ = eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay}))';
        TSsig_(isempty_cell(TSsig_))=[];
        TSsig_ = cell2mat(TSsig_')';
        
        TS(sum(TSsig_,2)==0,:)=0;
        
        % Sort by peak time
%         if iDelay ==1
            [~,idx] = max(TS_,[],2); 
            [~,idx] = sort(idx);
%         end
    % Sort by peak time, separate insignificant units
%     if iDelay ==1
%     idxNS = find(nansum(TSsig_,2)==0);
%     [~,idx] = max(TS_,[],2); 
%     [~,idx] = sort(idx);
%     idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];
%     end
    
    TS_ =TS_(idx,:);
    TSsig_ =TSsig_ (idx,:);
    TS_(~TSsig_ )=0;
    subplot(1,3,iDelay); hold on
    x_  = repmat((1:size(TS_,2)).*bw,size(TS_,1),1);
    y_  = repmat(1:size(TS_,1),size(TS_,2),1);
    imagesc(x_(:),y_(:),TS_);
    set(gca,'YDir','normal')
%     cmap =([1 1 1;flipud(hot)]);
    cmap =([1 1 1;(jet)]);
    colormap (cmap)
    caxis([0 10])
    plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 0 350])
    if iArea==1 & iDelay ==3
        plot([18 20],[20 20],'k','LineWidth',1.5)
        plot([18 18],[20 70],'k','LineWidth',1.5)
    end
%     axis off
    
    end
end
% clear idx TS_ TSsig_
%% Plot mean firing rates next to t-scores
% save('blue_white_redCropped.mat','cmap')
cmap1 = load('blue_white_redCropped.mat','cmap');cmap1 = cmap1.cmap;
cmap2 = ([0.9 0.9 0.9;(jet)]);
% cmap2 = jet;
iDelay = 2;

for iArea = 1:length(Areas)-1
    figure('name',Areas{iArea});
    
    TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
    TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
    % 1) mean FR over L/R trials
%         avgFR_  = (cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))' + ...
%                    cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))')./2;
    % 2) max FR at any point
    %     avgFR_ = zscore(avgFR_,[],2);
    % 3) preferred dirtection only
    avgFR_  = cat(3,zscore(cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))',[],2) ,...
                    zscore(cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))',[],2));
    avgFR_  = cat(3,cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))',...
                    cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))');
          avgFR_=     max(avgFR_,[],3);
%     idx =  squeeze(max(avgFR_,[],2));    idx = (diff(idx,[],2)>0)+1;
%     avgFR__= [];
%     for i = 1:size(avgFR_,1)
%         avgFR__ = [avgFR__;avgFR_(i,:,idx(i))];
%     end
%     avgFR_ = avgFR__; clear avgFR__ idx i            
    % TS(sum(TSsig_,2)==0,:)=0;
        avgFR_ = zscore(avgFR_,[],2);

    % Sort by peak time
%     if iDelay ==1
%         [~,idx] = max(avgFR_,[],2);
    % TS(sum(TSsig_,2)==0,:)=0;

        [~,idx] = max(TS_,[],2);
        [~,idx] = sort(idx);
%     end

    avgFR_ = avgFR_(idx,:);
    TS_    = TS_(idx,:);
    TSsig_ = TSsig_ (idx,:);
    TS_(~TSsig_ )=0;
    
    x_  = repmat((1:size(avgFR_,2)).*bw,size(avgFR_,1),1);
    y_  = repmat(1:size(avgFR_,1),size(avgFR_,2),1);
    
    subplot(1,2,1); hold on
    imagesc(x_(:),y_(:),avgFR_);
    set(gca,'YDir','normal')
    colormap (gca,cmap1)
    caxis(gca,[-2 2])
    plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 0 350])
    axis off

    cb = colorbar; cb.Location='manual'; cb.Position= [0.48,0.12,0.032,0.2]
    subplot(1,2,2); hold on
    imagesc(x_(:),y_(:),TS_);
    set(gca,'YDir','normal')
    colormap (gca,cmap2)
    caxis(gca,[0 10])
    plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 0 350])
%     if iArea==1 & iDelay ==3
        plot([18 20],[20 20],'k','LineWidth',1.5)
        plot([18 18],[20 70],'k','LineWidth',1.5)
%     end
    axis off
    cb = colorbar; cb.Location='manual'; cb.Position= [0.92,0.12,0.032,0.2];

    
end
%% Plot mean firing rates next to t-scores (preferred and non-preferred
% save('blue_white_redCropped.mat','cmap')
cmap1 = load('blue_white_redCropped.mat','cmap');cmap1 = cmap1.cmap;
cmap2 = ([0.9 0.9 0.9;(jet)]);
% cmap2 = jet;
iDelay = 2;

for iArea = 1:length(Areas)-1
     
    figure('name',Areas{iArea});
    
    TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
    TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
    Ltr = size(TS_,2);
    % 3) preferred direction only
%     avgFR_  = cat(3,zscore(cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))',[],2) ,...
%                     zscore(cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))',[],2));
    avgFR_  = cat(3,cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))',...
                    cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))');
%           avgFR_=     max(avgFR_,[],3);
    idx =  squeeze(max(avgFR_,[],2));    idxPref = (diff(idx,[],2)>0)+1; idxNonPref = (diff(idx,[],2)<0)+1;
    avgFR_Pref= [];
    avgFR_NonPref= [];
    for i = 1:size(avgFR_,1)
        avgFR_Pref    = [avgFR_Pref;   avgFR_(i,:,idxPref(i))];
        avgFR_NonPref = [avgFR_NonPref;avgFR_(i,:,idxNonPref(i))];
    end
    clear avgFR_ idx i         
    avgFR_ = [avgFR_Pref,avgFR_NonPref];
    
    avgFR_ = zscore(avgFR_,[],2);
    avgFR_Pref    = avgFR_(:,1:Ltr);
    avgFR_NonPref = avgFR_(:,(Ltr+1):2*Ltr);
    
%     avgFR_Pref = zscore(avgFR_Pref,[],2);
%     avgFR_NonPref = zscore(avgFR_NonPref,[],2);

    % Sort by peak time
% %         TS_(sum(TSsig_,2)==0,:)=0;
%         TS_(TSsig_==0)=0;

    [~,idx] = max(TS_,[],2);
%     [~,idx] = max(avgFR_Pref,[],2);
    [~,idx] = sort(idx);

    avgFR_Pref = avgFR_Pref(idx,:);
    avgFR_NonPref = avgFR_NonPref(idx,:);
    TS_    = TS_(idx,:);
    TSsig_ = TSsig_ (idx,:);
    TS_(~TSsig_ )=0;
    
    x_  = repmat((1:size(avgFR_Pref,2)).*bw,size(avgFR_Pref,1),1);
    y_  = repmat(1:size(avgFR_Pref,1),size(avgFR_Pref,2),1);
    
    subplot(1,3,1); hold on
        imagesc(x_(:),y_(:),avgFR_Pref);
        set(gca,'YDir','normal')
        colormap (gca,cmap1)
        caxis(gca,[-2 2])
        plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
        plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
        plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
        axis([0 20 0 350])
        axis off
    
    subplot(1,3,2); hold on
        imagesc(x_(:),y_(:),avgFR_NonPref);
        set(gca,'YDir','normal')
        colormap (gca,cmap1)
        caxis(gca,[-2 2])
        plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
        plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
        plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
        axis([0 20 0 350])
        axis off
    
    
    
    cb = colorbar; cb.Location='manual'; cb.Position = [0.63,0.12,0.032,0.2];
    subplot(1,3,3); hold on
    imagesc(x_(:),y_(:),TS_);
    set(gca,'YDir','normal')
    colormap (gca,cmap2)
    caxis(gca,[0 10])
    plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 0 350])
%     if iArea==1 & iDelay ==3
        plot([18 20],[20 20],'k','LineWidth',1.5)
        plot([18 18],[20 70],'k','LineWidth',1.5)
%     end
    axis off
    cb = colorbar; cb.Location='manual'; cb.Position= [0.91,0.12,0.032,0.2];

    tempOoutput.FRpreferred{iArea} = avgFR_Pref;
    tempOoutput.FRnonpreferred{iArea} = avgFR_NonPref;
    tempOoutput.TS{iArea} = TS_;
    tempOoutput.TSsig{iArea} = TSsig_;
end
%% Plot mean firing rates (preferred and non-preferred next to t-scores
iDelay = 2;
for iFile =1:length(fileList)
    TS_ = eval(sprintf('D_.%s.LR.TS{iArea}{iFile}',Delays_{iDelay}))';
    TSsig_ = (eval(sprintf('D_.%s.LR.TSsig{iArea}{iFile}',Delays_{iDelay})))';
    Ltr = size(TS_,2);
    % 3) preferred direction only
    
    avgFR_  = cat(3,(eval(sprintf('D_.%s.LR.avgFR{iArea}{1,iFile}',Delays_{iDelay})))',...
                    (eval(sprintf('D_.%s.LR.avgFR{iArea}{2,iFile}',Delays_{iDelay})))');
%           avgFR_=     max(avgFR_,[],3);
    idx =  squeeze(max(avgFR_,[],2));    idxPref = (diff(idx,[],2)>0)+1; idxNonPref = (diff(idx,[],2)<0)+1;
    avgFR_Pref= [];
    avgFR_NonPref= [];
    for i = 1:size(avgFR_,1)
        avgFR_Pref    = [avgFR_Pref;   avgFR_(i,:,idxPref(i))];
        avgFR_NonPref = [avgFR_NonPref;avgFR_(i,:,idxNonPref(i))];
    end
    clear avgFR_ idx i         
    avgFR_ = [avgFR_Pref,avgFR_NonPref];
    
    avgFR_ = zscore(avgFR_,[],2);
    avgFR_Pref    = avgFR_(:,1:Ltr);
    avgFR_NonPref = avgFR_(:,(Ltr+1):2*Ltr);
    
    tempOoutput.FRpreferred{iArea}(iFile,:) = nanmean(avgFR_Pref);
    tempOoutput.FRnonpreferred{iArea}(iFile,:) = nanmean(avgFR_NonPref);
end
    
figure
for iArea = 1:length(Areas)-1
    
    
   
    
    
    
    subplot(1,3,1); hold on
    ciplot(nanmean( tempOoutput.FRpreferred{iArea})+nansem( tempOoutput.FRpreferred{iArea}),...
           nanmean( tempOoutput.FRpreferred{iArea})-nansem( tempOoutput.FRpreferred{iArea}),(1:size( tempOoutput.FRpreferred{iArea},2)).*bw,color_{iArea},0.5);
   
    plot([5 5],0.2+[0 0.5],'color',[0 1 0 0.6],'LineWidth',2)
    plot([10 10],0+[0 0.5],'color',[0 0 0 0.6],'LineWidth',2)
    plot([15 15],-0.15+[0 0.5],'color',[1 0 0 0.6],'LineWidth',2)
    axis([0 20 -0.4 1.6])
    axis off

    subplot(1,3,2); hold on
    ciplot(nanmean(tempOoutput.FRnonpreferred{iArea})+nansem(tempOoutput.FRnonpreferred{iArea}),...
           nanmean(tempOoutput.FRnonpreferred{iArea})-nansem(tempOoutput.FRnonpreferred{iArea}),(1:size(tempOoutput.FRnonpreferred{iArea},2)).*bw,color_{iArea},0.6);
   
    plot([5 5],-0.4+[0 0.5],'color',[0 1 0 0.6],'LineWidth',2)
    plot([10 10],-0.7+[0 0.5],'color',[0 0 0 0.6],'LineWidth',2)
    plot([15 15],-0.7+[0 0.5],'color',[1 0 0 0.6],'LineWidth',2)
    axis([0 20 -1 1])
    axis off
    
    subplot(1,3,3); hold on
%     TSsig_ = tempOoutput.TSsig{iArea};
    TSsig_ = eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay}))';
    TSsig_ = cell2mat(cellfun(@(x) sum(x,2)./size(x,2),TSsig_,'UniformOutput',false)')';
    ciplot(nanmean(TSsig_)+nansem(TSsig_),...
           nanmean(TSsig_)-nansem(TSsig_),(1:size(TSsig_,2)).*bw,color_{iArea},0.5);
    
%     stairs((1:size(TSsig_,2)).*bw,nanmean(TSsig_)+nansem(TSsig_),color_{iArea})
%     stairs((1:size(TSsig_,2)).*bw,nanmean(TSsig_)-nansem(TSsig_),color_{iArea})
%     stairs((1:size(TSsig_,2)).*bw,nanmean(TSsig_),color_{iArea},'LineWidth',1.5)
    plot([5 5],0.08+[0 0.2],'color',[0 1 0 0.6],'LineWidth',2)
    plot([10 10],[0.03 0.13],'color',[0 0 0 0.6],'LineWidth',2)
    plot([15 15],0.08+[0 0.2],'color',[1 0 0 0.6],'LineWidth',2)
    axis([0 20 0 0.6])
%     if iArea==1 & iDelay ==3
%         plot([18 20],[20 20],'k','LineWidth',1.5)
%         plot([18 18],[20 70],'k','LineWidth',1.5)
%     end
    axis off
    
end
%% Plot sort order of FR vs TS
% for iDelay =1:3
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on 
    iDelay = 2;
    TS_    = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
    TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
    
    % 1) mean FR over L/R trials
%         avgFR_  = (cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))' + ...
%                    cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))')./2;
%         avgFR_ = zscore(avgFR_,[],2);

    % 2)
%     avgFR_  = cat(3,zscore(cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))',[],2) ,...
%         zscore(cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))',[],2));
    avgFR_  = cat(3,cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))' ,...
                    cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))');
% 2.1) max FR at any point
%     avgFR_  = zscore(max(avgFR_,[],3),[],2);
%     % 2.2)  preferred direction only
    idx =  squeeze(max(avgFR_,[],2));    
    idxPref = (diff(idx,[],2)>0)+1; 
    idxNonPref = (diff(idx,[],2)<0)+1;
    avgFR_Pref = [];
    avgFR_NonPref = [];
    for i = 1:size(avgFR_,1)
        avgFR_Pref    = [avgFR_Pref;   avgFR_(i,:,idxPref(i))];
        avgFR_NonPref = [avgFR_NonPref;avgFR_(i,:,idxNonPref(i))];
    end
    clear avgFR_ idx i
%     avgFR_ = [avgFR_Pref,avgFR_NonPref];
%     avgFR_ = zscore(avgFR_,[],2);
%     avgFR_Pref    = avgFR_(:,1:Ltr);
%     avgFR_NonPref = avgFR_(:,(Ltr+1):2*Ltr);
%         
    avgFR_Pref = zscore(avgFR_Pref,[],2);
    avgFR_NonPref = zscore(avgFR_NonPref,[],2);

    avgFR_ = avgFR_Pref;
    
   
    % remove any cells with no sig L/R Decoding 
        %     TS_(~TSsig_ )=0;
        %     TS_(sum(TSsig_,2)==0,:)=[];
        %     avgFR_(sum(TSsig_,2)==0,:)=[];
    
    % Threshold by mean FR
        %         idx = max(avgFR_,[],2)<0.5;
        %         TS_(idx,:)=[];
        %         avgFR_(idx,:)=[];

    % Sort by peak time
    [~,idxFR] = max(avgFR_,[],2);
    [~,idxTS] = max(TS_,[],2);
%     idx = idxFR<=1|idxTS<=1;
%     idxTS(idx)=[]; idxFR(idx)=[];


% subplot(1,2,iArea); hold on;

scatter(idxFR*bw,idxTS*bw,10,color_{iArea},'filled')



mdl{iArea} = fitlm(idxFR*bw,idxTS*bw);
pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
%             pFit(iArea)   = coefTest(mdl{iArea});
Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
slope(iArea) = mdl{iArea}.Coefficients{2,1};

h = plot(mdl{iArea});
if pFit(iArea)<0.05
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = color_{iArea};
    h(3).Color = color_{iArea};
    h(4).Color = color_{iArea};

%     h(3).LineStyle='none';
%     h(4).LineStyle='none';
else
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = color_{iArea};
    h(3).Color = color_{iArea};
    h(4).Color = color_{iArea};
%     h(3).LineStyle='none';
%     h(4).LineStyle='none';
end
% a = ezfit(idxFR*bw,idxTS*bw,'poly1'); showfit(a)
legend off
plot([0 20],[0 20],':k')
axis([0 20 0 20])
xlabel('Latency to peak firing rate deflection (s)')
ylabel('Latency to peak cue encoding (s)')
title(Areas{iArea})
if pFit(iArea)<0.01
    text(1,17,{sprintf('R^2=%.1g',Rsq(iArea));'p<0.01'})
else
    text(1,17,{sprintf('R^2=%.1g',Rsq(iArea));sprintf('p=%.2g',pFit(iArea))})
end


%     [~,idxFR] = sort(idxFR);
%     [~,idxTS] = sort(idxTS);
[corr_(iArea),pval_(iArea)] = corr(idxFR,idxTS,'type','spearman');

end
% end
    
%% Plot mean firing rates next to t-scores - Errors
cmap1 = load('blue_white_redCropped.mat','cmap');cmap1 = cmap1.cmap;
cmap2 = ([0.8 0.8 0.8;(jet)]);
% cmap2 = jet;
iDelay = 3;

for iArea = 1:length(Areas)-1
    figure
    
    TS_ = cell2mat(eval(sprintf('D_.%s.LR_err.TS_err{iArea}',Delays_{iDelay})))';
    TSsig_ = cell2mat(cellfun(@double, eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay})), 'UniformOutput', false))';
    % 1) mean FR over L/R trials
        avgFR_  = (cell2mat(eval(sprintf('D_.%s.LR_err.avgFR_err{iArea}(1,:)',Delays_{iDelay})))' + ...
                   cell2mat(eval(sprintf('D_.%s.LR_err.avgFR_err{iArea}(2,:)',Delays_{iDelay})))')./2;
    % 2) max FR at any point
    %     avgFR_ = zscore(avgFR_,[],2);
    % 3) preferred dirtection only
%     avgFR_  = cat(3,zscore(cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))',[],2) ,...
%                     zscore(cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))',[],2));
%     avgFR_  = cat(3,cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(1,:)',Delays_{iDelay})))',...
%                     cell2mat(eval(sprintf('D_.%s.LR.avgFR{iArea}(2,:)',Delays_{iDelay})))');
%           avgFR_=     max(avgFR_,[],3);
%     idx =  squeeze(max(avgFR_,[],2));    idx = (diff(idx,[],2)>0)+1;
%     avgFR__= [];
%     for i = 1:size(avgFR_,1)
%         avgFR__ = [avgFR__;avgFR_(i,:,idx(i))];
%     end
%     avgFR_ = avgFR__; clear avgFR__ idx i            
    % TS(sum(TSsig_,2)==0,:)=0;
        avgFR_ = zscore(avgFR_,[],2);

    % Sort by peak time
%     if iDelay ==1
%         [~,idx] = max(avgFR_,[],2);
        [~,idx] = max(TS_,[],2);
        [~,idx] = sort(idx);
%     end

    avgFR_ = avgFR_(idx,:);
    TS_    = TS_(idx,:);
    TSsig_ = TSsig_ (idx,:);
    TS_(~TSsig_ )=0;
    
    x_  = repmat((1:size(avgFR_,2)).*bw,size(avgFR_,1),1);
    y_  = repmat(1:size(avgFR_,1),size(avgFR_,2),1);
    
    subplot(1,2,1); hold on
    imagesc(x_(:),y_(:),avgFR_);
    set(gca,'YDir','normal')
    colormap (gca,cmap1)
    caxis(gca,[-2 2])
    plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 0 350])
    axis off

    subplot(1,2,2); hold on
    imagesc(x_(:),y_(:),TS_);
    set(gca,'YDir','normal')
    colormap (gca,cmap2)
    caxis(gca,[0 10])
    plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 0 350])
%     if iArea==1 & iDelay ==3
        plot([18 20],[20 20],'k','LineWidth',1.5)
        plot([18 18],[20 70],'k','LineWidth',1.5)
%     end
    axis off
    
end
%% Fraction sig on corr/error trials
iDelay =3;
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    
    TSsig_ = cellfun(@double, eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
    TSsig_(isempty_cell(TSsig_))=[];
    TSsig_ = cell2mat(TSsig_');
    TSsig_ = (nansum(TSsig_,2)./size(TSsig_,2))*100';
    plot((1:size(TSsig_,1)).*bw,TSsig_,'LineWidth',1.5,'Color',color_{iArea})
    
       
    TSsig_ = cellfun(@double, eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
    TSsig_(isempty_cell(TSsig_))=[];
    TSsig_ = cell2mat(TSsig_');
    TSsig_ = (nansum(TSsig_,2)./size(TSsig_,2))*100';
    plot((1:size(TSsig_,1)).*bw,TSsig_,':','LineWidth',1.5,'Color',color_{iArea})
    
end

figure
for iArea = 1:length(Areas)-1
    
    TSsig_CorrM =[];
    TSsig_ErrM =[];
    
    for iDelay = 1:3
        TSsig_ = cellfun(@double, eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
        TSsig_(isempty_cell(TSsig_))=[];
        TSsig_ = cell2mat(TSsig_');
        TSsig_ = (nansum(TSsig_,2)./size(TSsig_,2))*100';
        TSsig_CorrM=[TSsig_CorrM,TSsig_];
        
        TSsig_ = cellfun(@double, eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
        TSsig_(isempty_cell(TSsig_))=[];
        TSsig_ = cell2mat(TSsig_');
        TSsig_ = (nansum(TSsig_,2)./size(TSsig_,2))*100';
        TSsig_ErrM=[TSsig_ErrM,TSsig_];
    end
        subplot(1,2,iArea); hold on

       
%     stairs((1:size(TSsig_CorrM,1)).*bw,nanmean(TSsig_CorrM,2),'LineWidth',1.5,'Color',color_{iArea})
%     stairs((1:size(TSsig_ErrM,1)).*bw,nanmean(TSsig_ErrM,2),':','LineWidth',1.5,'Color',color_{iArea})
    
     stairs((1:size(TSsig_CorrM,1)).*bw,(nanmean(TSsig_ErrM,2)-nanmean(TSsig_CorrM,2))./(nanmean(TSsig_CorrM,2)+nanmean(TSsig_ErrM,2)),'LineWidth',1.5,'Color',color_{iArea})
    
end
%% Fraction sig on corr/err trials - lumped
iDelay = 1;
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    
    TSsig_ = cellfun(@double, eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
    TSsig_(isempty_cell(TSsig_))=[];
    TSsig_ = cell2mat(cellfun(@(x) (nansum(x,2)./size(x,2))*100,TSsig_,'UniformOutput',false)')';
    ciplot(nanmean(TSsig_)+nansem(TSsig_),nanmean(TSsig_)-nansem(TSsig_),...
           (1:size(TSsig_,2)).*bw,color_{iArea},0.9)
       
    TSsig_ = cellfun(@double, eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
    TSsig_(isempty_cell(TSsig_))=[];
    TSsig_ = cell2mat(cellfun(@(x) (nansum(x,2)./size(x,2))*100,TSsig_,'UniformOutput',false)')';
    ciplot(nanmean(TSsig_)+nansem(TSsig_),nanmean(TSsig_)-nansem(TSsig_),...
           (1:size(TSsig_,2)).*bw,color_{iArea},0.3)
end
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    TSsig_CorrM =[];
    TSsig_CorrE =[];
    for iDelay = 1:3
        TSsig_ = cellfun(@double, eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
        TSsig_(isempty_cell(TSsig_))=[];
        TSsig_ = cell2mat(cellfun(@(x) (nansum(x,2)./size(x,2))*100,TSsig_,'UniformOutput',false)')';
        TSsig_CorrM=[TSsig_CorrM;TSsig_];
        
        TSsig_ = cellfun(@double, eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
        TSsig_(isempty_cell(TSsig_))=[];
        TSsig_ = cell2mat(cellfun(@(x) (nansum(x,2)./size(x,2))*100,TSsig_,'UniformOutput',false)')';
        TSsig_ErrM=[TSsig_CorrE;TSsig_];
        
    end
    ciplot(nanmean(TSsig_CorrM)+nansem(TSsig_CorrM),nanmean(TSsig_CorrM)-nansem(TSsig_CorrM),...
        (1:size(TSsig_CorrM,2)).*bw,color_{iArea},0.9)
    ciplot(nanmean(TSsig_ErrM)+nansem(TSsig_ErrM),nanmean(TSsig_ErrM)-nansem(TSsig_ErrM),...
        (1:size(TSsig_ErrM,2)).*bw,color_{iArea},0.3)
end
%%
figure
for iArea = 1:length(Areas)-1
    for iDelay = 1:3
        TSsig_C = cellfun(@double, eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
        TSsig_E = cellfun(@double, eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
        for iFile = 1:length(TSsig_C)
            try
                %TSsig_ratio{iFile} = (TSsig_E{iFile}-TSsig_C{iFile})./(TSsig_E{iFile}+TSsig_C{iFile});
                TSsig_C_ = nansum(TSsig_C{iFile},2)./size(TSsig_C{iFile},2)*100';
                TSsig_E_ = nansum(TSsig_E{iFile},2)./size(TSsig_E{iFile},2)*100';

                TSsig_ratio{iFile} = (TSsig_E_-TSsig_C_)./(TSsig_E_+TSsig_C_);
            catch
                %TSsig_ratio{iFile} = nan(size(TSsig_C{iFile}));
                TSsig_ratio{iFile}  = nan(size(TSsig_C{iFile},1),1);
            end
        end
        ratio_(:,:,iDelay) = cell2mat(TSsig_ratio);
    end
    ratio_M = nanmean(nanmean(ratio_,3),2);
    ratio_E = nansem(nanmean(ratio_,3),2);
%     ratio_M = nanmean(cell2mat(TSsig_ratio),2);
%     ratio_E = nansem(cell2mat(TSsig_ratio),2);
    
    subplot(1,2,iArea); hold on
    ciplot(ratio_M+ratio_E,...
           ratio_M-ratio_E,(1:size(avgFR_,2)).*bw,color_{iArea},0.6);
   
    plot([10 10],[-2 2],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([0 20],[0 0],':k')
    plot([5 5],[-2 2],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[-2 2],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 -2 2])
    axis off
    
    
end
%%
    figure
for iArea = 1:length(Areas)-1
    ratio_ = nan(402,10,3);
    for iDelay = 1:3
        TSsig_C = cellfun(@double, eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
        TSsig_E = cellfun(@double, eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay})), 'UniformOutput', false)';
        for iFile = 1:length(TSsig_C)
            try
%                TSsig_ratio{iFile} = nansum(smooth2a(TSsig_E{iFile},5,0) - smooth2a(TSsig_C{iFile},5,0),2)./size(TSsig_C{iFile},2)*100';
               TSsig_ratio{iFile} = nansum(smooth2a(TSsig_E{iFile}-TSsig_C{iFile},5,0),2)./size(TSsig_C{iFile},2)*100';
            catch
                %TSsig_ratio{iFile} = nan(size(TSsig_C{iFile}));
                TSsig_ratio{iFile}  = nan(size(TSsig_C{iFile},1),1);
            end
        end
        ratio_(:,:,iDelay) = cell2mat(TSsig_ratio);
    end
    ratio_M = nanmean(nanmean(ratio_,3),2);
    ratio_E = nansem(nanmean(ratio_,3),2);
%     ratio_M = nanmean(cell2mat(TSsig_ratio),2);
%     ratio_E = nansem(cell2mat(TSsig_ratio),2);
    
    subplot(1,2,iArea); hold on
    ciplot(ratio_M+ratio_E,...
           ratio_M-ratio_E,(1:size(avgFR_,2)).*bw,color_{iArea},0.6);
   
    plot([10 10],[-2 2],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([0 20],[0 0],':k')
    plot([5 5],[-2 2],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[-2 2],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 -10 10])
    axis off
    
    
end


%% Distributions of peak t-Scores - shaded
plotErr = false;
bins = 0:0.1:10;  
alpha_ = 0.3;
alpha_Err = 0.3;
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    for iDelay = 1:length(Delays_)
       TS_ = eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay}));
       TS_ = cellfun(@(x) max(x,[],1),TS_,'UniformOutput',false);
       for iFile =1:length(TS_)
           maxhist = cumsum(histc(TS_{iFile},bins));
           maxHist(:,iFile) = maxhist./max(maxhist);
       end
       maxHist_Mean = nanmean(maxHist,2);
       maxHist_SEM = nansem(maxHist,2);
       ciplot(maxHist_Mean+maxHist_SEM,maxHist_Mean-maxHist_SEM,...
              bins,color_{iDelay},alpha_)
          
          if plotErr
       stairs(bins,maxHist_Mean,'LineWidth',1.5,'color',color_{iDelay});

       TS_Err = eval(sprintf('D_.%s.LR_err.TS_err{iArea}',Delays_{iDelay}));
       TS_Err(isempty_cell(TS_Err))=[];
       TS_Err = cellfun(@(x) max(x,[],1),TS_Err,'UniformOutput',false);
       for iFile =1:length(TS_Err)
           try
              maxhist = cumsum(histc(TS_Err{iFile},bins));
              maxHist(:,iFile) = maxhist./max(maxhist);
           end
       end
       maxHist_Mean = nanmean(maxHist,2);
       maxHist_SEM  = nansem(maxHist,2);
%        ciplot(maxHist_Mean+maxHist_SEM,maxHist_Mean-maxHist_SEM,...
%               bins,color_{iDelay},alpha_Err)
       stairs(bins,maxHist_Mean,'LineWidth',1.5,'LineStyle',':','color',color_{iDelay});
          end

        set(gca,'Ytick',0:0.5:1,'Xtick',0:5:20)
        xlabel('Peak t-score')
        ylabel('Fraction of units')
        title(Areas{iArea})
    end
end
%% Distributions of peak t-Scores, tim sig. - both areas
plotErr = true;
bins = 0:0.1:10;  
alpha_ = 0.9;
alpha_Err = 0.3;
figure
subplot(1,2,1); hold on
for iArea = 1:length(Areas)-1    
    for iDelay = 2%1:length(Delays_)
        TS_ = eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay}));
        TS_ = cellfun(@(x) max(x,[],1),TS_,'UniformOutput',false);
        for iFile =1:length(TS_)
            maxhist = cumsum(histc(TS_{iFile},bins));
            maxHist(:,iFile) = maxhist./max(maxhist);
        end
        maxHist_Mean = nanmean(maxHist,2);
        maxHist_SEM = nansem(maxHist,2);
        ciplot(maxHist_Mean+maxHist_SEM,maxHist_Mean-maxHist_SEM,...
            bins,color_{iArea},alpha_)
        stairs(bins,maxHist_Mean,'LineWidth',1.5,'color',color_{iArea});
        if plotErr
            
            TS_Err = eval(sprintf('D_.%s.LR_err.TS_err{iArea}',Delays_{iDelay}));
            TS_Err(isempty_cell(TS_Err))=[];
            TS_Err = cellfun(@(x) max(x,[],1),TS_Err,'UniformOutput',false);
            for iFile =1:length(TS_Err)
                try
                    maxhist = cumsum(histc(TS_Err{iFile},bins));
                    maxHist(:,iFile) = maxhist./max(maxhist);
                end
            end
            maxHist_Mean = nanmean(maxHist,2);
            maxHist_SEM  = nansem(maxHist,2);
                   ciplot(maxHist_Mean+maxHist_SEM,maxHist_Mean-maxHist_SEM,...
                          bins,color_{iArea},alpha_Err)
            stairs(bins,maxHist_Mean,'LineWidth',1.5,'LineStyle',':','color',color_{iArea});
        end
    end
end
set(gca,'Ytick',0:0.5:1,'Xtick',0:5:20)
xlabel('Peak t-score')
ylabel('Fraction of units')
title('Peak decoding')

clear timehist timeHist FracSig FracSigErr
Sigthresh = bw;
maxBin = 15;
bins   = Sigthresh:bw:maxBin;
subplot(1,2,2); hold on
for iArea = 1:length(Areas)-1
    for iDelay  = 2%1:length(Delays_)
       TSsig_   = eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay}));
       
%        TSsig_ = cellfun(@(x) sum(x(202:end,:),1).*bw,TSsig_,'UniformOutput',false);   
       TSsig_ = cellfun(@(x) sum(x,1).*bw,TSsig_,'UniformOutput',false);
       FracSig{iArea}(:,iDelay) = cell2mat(cellfun(@(x) sum(sum(x,1).*bw > Sigthresh)./size(x,2),TSsig_,'UniformOutput',false));
       
        for iFile =1:length(TSsig_)

           timehist = cumsum(histc(TSsig_{iFile},bins));
           timeHist(:,iFile) = timehist./max(timehist);
       end
       
       timeHist_Mean = nanmean(timeHist,2);
       timeHist_SEM = nansem(timeHist,2);
       ciplot(timeHist_Mean+timeHist_SEM,timeHist_Mean-timeHist_SEM,...
              bins,color_{iArea},alpha_)
%        plot(bins,timeHist_Mean,'Color',color_{iDelay},'LineWidth',1.5)   
if plotErr
    TSsig_Err   = eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay}));
    emptyYN = isempty_cell(TSsig_Err); 
    TSsig_Err = cellfun(@(x) sum(x,1).*bw,TSsig_Err,'UniformOutput',false);

    for iFile =1:length(emptyYN)
        %        TSsig_ = cellfun(@(x) sum(x(202:end,:),1).*bw,TSsig_,'UniformOutput',false);
        FracSigErr{iArea}(iFile,iDelay) = sum(sum( TSsig_Err{iFile},1).*bw > Sigthresh)./size( TSsig_Err{iFile},2)
    
        timehist = cumsum(histc(TSsig_Err{iFile},bins));
        timeHist(:,iFile) = timehist./max(timehist);
    end
    
    timeHist_Mean = nanmean(timeHist,2);
    timeHist_SEM = nansem(timeHist,2);
    ciplot(timeHist_Mean+timeHist_SEM,timeHist_Mean-timeHist_SEM,...
           bins,color_{iArea},alpha_Err)
    plot(bins,timeHist_Mean,'Color',color_{iArea},'LineWidth',1.5,'LineStyle',':')
end
   
    end
end
axis([0 15 0. 1])
set(gca,'Ytick',0:0.5:1,'Xtick',0:5:maxBin)
xlabel('Time (s)')
% ylabel('Fraction of units')
title('Span of Significant decoding')
%% Distributions of time significant - shaded
plotErr = false;
figure
clear timehist timeHist FracSig FracSigErr
Sigthresh = bw;
maxBin = 15;
bins   = Sigthresh:bw:maxBin;
alpha_ = 0.3;
alpha_Err = 0.2;

for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    for iDelay  = 1:length(Delays_)
       TSsig_   = eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay}));
       
%        TSsig_ = cellfun(@(x) sum(x(202:end,:),1).*bw,TSsig_,'UniformOutput',false);   
       TSsig_ = cellfun(@(x) sum(x,1).*bw,TSsig_,'UniformOutput',false);
       FracSig{iArea}(:,iDelay) = cell2mat(cellfun(@(x) sum(sum(x,1).*bw > Sigthresh)./size(x,2),TSsig_,'UniformOutput',false));
       
        for iFile =1:length(TSsig_)

           timehist = cumsum(histc(TSsig_{iFile},bins));
           timeHist(:,iFile) = timehist./max(timehist);
       end
       
       timeHist_Mean = nanmean(timeHist,2);
       timeHist_SEM = nansem(timeHist,2);
       ciplot(timeHist_Mean+timeHist_SEM,timeHist_Mean-timeHist_SEM,...
              bins,color_{iDelay},alpha_)
%        plot(bins,timeHist_Mean,'Color',color_{iDelay},'LineWidth',1.5)   
if plotErr
    TSsig_Err   = eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay}));
    emptyYN = isempty_cell(TSsig_Err); 
    TSsig_Err = cellfun(@(x) sum(x,1).*bw,TSsig_Err,'UniformOutput',false);

    for iFile =1:length(emptyYN)
        %        TSsig_ = cellfun(@(x) sum(x(202:end,:),1).*bw,TSsig_,'UniformOutput',false);
        FracSigErr{iArea}(iFile,iDelay) = sum(sum( TSsig_Err{iFile},1).*bw > Sigthresh)./size( TSsig_Err{iFile},2)
    
        timehist = cumsum(histc(TSsig_Err{iFile},bins));
        timeHist(:,iFile) = timehist./max(timehist);
    end
    
    timeHist_Mean = nanmean(timeHist,2);
    timeHist_SEM = nansem(timeHist,2);
%     ciplot(timeHist_Mean+timeHist_SEM,timeHist_Mean-timeHist_SEM,...
%            bins,color_{iDelay},alpha_Err)
    plot(bins,timeHist_Mean,'Color',color_{iDelay},'LineWidth',1.5,'LineStyle',':')
end
   title(Areas{iArea})
    end
    axis([0 15 0. 1])
    set(gca,'Ytick',0:0.5:1,'Xtick',0:5:maxBin)
    xlabel('Time (s)')
    ylabel('Fraction of units')
    title(Areas{iArea})
end
%% Fraction significant (error bars)
color_ = {[0 0 0],[0 0 0],[0 0 0]}
alpha_ = 0.9;
figure
for iArea = [2,1]%1:length(Areas)-1
    subplot(1,2,iArea); hold on
    for iDelay  = 1:length(Delays_)

     M = nanmean(FracSig{iArea}(:,iDelay));
     E = nansem(FracSig{iArea}(:,iDelay));
       bar(iDelay,M,...
         'FaceColor','w',...
         'EdgeColor',0.3*[1 1 1],...
         'FaceAlpha',0.3,...
         'LineWidth',1.5)
     scatter(iDelay*ones(size(FracSig{iArea},1),1),FracSig{iArea}(:,iDelay),10, ...
         'MarkerEdgeColor',color_{iDelay}, ...
         'MarkerFaceColor','none', ...
         'LineWidth',1.5,'MarkerEdgeAlpha',alpha_)
     errorbar(iDelay,M,E,'Color',0.3*[1 1 1],'LineWidth',1.5)
   
    end
        plot(FracSig{iArea}','LineStyle',':','Color',0.3*[1 1 1])
        ylabel('Frac. informative units')
        set(gca,'Ytick',0:0.5:1,'Xtick',[1,2,3],'XTickLabel',Delays__,'XTickLabelRotation',90)
        xlabel('Delay')
        title(Areas{iArea})
    axis([0 4 0 1])
end

%%

p = kruskalwallis(FracSig{2})
%% Fraction significant (error bars) - with errors
clear M E FS

for iArea = 1:length(Areas)-1
    M{iArea}=[];
    E{iArea}=[];
    FS{iArea}=[];
    subplot(1,2,iArea); hold on
    for iDelay  = 1:length(Delays_)

     M{iArea} = [M{iArea}, nanmean(FracSig{iArea}(:,iDelay))];
     M{iArea} = [M{iArea}, nanmean(FracSigErr{iArea}(:,iDelay))];
     E{iArea} = [E{iArea}, nansem(FracSig{iArea}(:,iDelay))];
     E{iArea} = [E{iArea}, nansem(FracSigErr{iArea}(:,iDelay))];
     
     FS{iArea} = [FS{iArea},FracSig{iArea}(:,iDelay)];
     FS{iArea} = [FS{iArea},FracSigErr{iArea}(:,iDelay)];
%      scatter(iDelay*ones(size(FracSig{iArea},1),1),FracSig{iArea}(:,iDelay),10, ...
%          'MarkerEdgeColor',color_{iDelay}, ...
%          'MarkerFaceColor','none', ...
%          'LineWidth',1.5,'MarkerEdgeAlpha',0.3)

    end

end
x =[1 2 4 5 7 8];x_ = [1 3 5];
color2_={'b','b','g','g','r','r'};
alpha2_=[0.9 0.2 0.9 0.2 0.9 0.2];
figure
for iArea = 1:length(Areas)-1
        subplot(1,2,iArea); hold on
        
        for i=1:6
            bar(x(i),M{iArea}(i),...
            'FaceColor',color2_{i},...
            'EdgeColor',0.3*[1 1 1],...
            'FaceAlpha',alpha2_(i),...
            'LineWidth',1.5)
        
        scatter(x(i)*ones(size(FS{iArea},1),1),FS{iArea}(:,i),10,...
           'MarkerEdgeColor',color2_{i}, ...
         'MarkerFaceColor','none', ...
         'LineWidth',1.5,'MarkerEdgeAlpha',0.3)
        end
        errorbar(x,M{iArea},E{iArea},'Color',0.3*[1 1 1],'LineWidth',1.5,'LineStyle','none')
        set(gca,'Ytick',0:0.5:1,'Xtick',[1.5,4.5,7.5],'XTickLabel',Delays__,'XTickLabelRotation',90)
        
                for i=1:3
                    plot([x(x_(i)):x(x_(i)+1)],FS{iArea}(:,x_(i):x_(i)+1)','LineStyle',':','Color',0.3*[1 1 1])
                end
        ylabel('Frac. informative units')

        xlabel('Delay')
        title(Areas{iArea})
    axis([0 9 0 1])
end

%% Fractions of significant-encoding neurons
clear TSsig_
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    title(Areas{iArea})
    bar(FracNonSig{iArea})
end
%% Fractions of significant-encoding neurons changing from short trials
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    clear TSsig_

    for iDelay  = 1:length(Delays_)
        x = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})));
       TSsig_(:,iDelay)   = sum(x)';       
    end
    TSsig_=TSsig_+1e-9;
    TSsig_=TSsig_./TSsig_(:,1);
    plot(TSsig_')
end
%% strength of encoding changing from short trials
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    clear TSsig_

    for iDelay  = 1:length(Delays_)
       TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
       TSsig_(:,iDelay)   = max(TS_,[],2);       
    end

    TSsig_=TSsig_-TSsig_(:,1);
    plot((TSsig_)')
end

%% peak strength vs time significant
col_ = {'r','b','g'};
figure;
for iArea = 1:length(Areas)-1
    subplot(1,3,iArea); hold on
    
    for iDelay = 1:length(Delays_)
        TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
        TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';

    max_  = max(TS_,[],2); 
    duration = sum(TSsig_,2).*bw;
    max_(duration==0)=[];
    duration(duration==0)=[];
    
    scatter(duration, max_,col_{iDelay})
    end
    axis([0 20 0 10])
end
clear idx TS_ TSsig_
%%
col_ = {'r','b','g'};
figure;
for iArea = 1:length(Areas)-1
        subplot(1,3,iArea); hold on
    for iDelay = 3:length(Delays_)
        TS_    = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
        TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
        idx    = 1:size(TS_,2);
        idx2   = 1:(abs(tlimsAll(1))/bw);
        idx([idx2,length(idx)-idx2])=[];
        TS_    = TS_(:,idx);
        TSsig_ = TSsig_(:,idx);
        
        clear t_max t_min
        
        duration = sum(TSsig_,2).*bw;
        [~,idx] = sort(duration )
        TSsig_ = sortrows([duration,TSsig_])
%         TSsig_ = TSsig_(idx,:)
        imagesc((TSsig_))
    end 
end
        
        %% Significant span vs. time of sig encoding
col_ = {'r','b','g'};
figure;
for iArea = 1:length(Areas)-1
        subplot(1,3,iArea); hold on

    for iDelay = 3:length(Delays_)
        TS_    = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
        TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
        
        idx    = 1:size(TS_,2);
        idx2   = 1:(abs(tlimsAll(1))/bw);
        idx([idx2,length(idx)-idx2])=[];
        TS_    = TS_(:,idx);
        TSsig_ = TSsig_(:,idx);
        
        clear t_max t_min
        

        duration      = sum(TSsig_,2).*bw;

%         [~,t_max]  = max(TS_,[],2);        
        for iUnit = 1:size(TSsig_,1)
            a = TSsig_(iUnit ,:);
            first_idx = find(a,1,'first'); if isempty(first_idx), first_idx=NaN; end
            last_idx = find(a,1,'last'); if isempty(last_idx), last_idx=NaN; end
            t_min(iUnit)=first_idx*bw;
            t_max(iUnit)=last_idx*bw;
        end
%         t_max(t_max==10.05)=NaN;
        t_mid =t_max-t_min;
%         t_max(duration==0)=[];
%         duration(duration==0)=[];
        scatter(t_max,duration,col_{iDelay})
%         scatter(t_mid,duration,col_{iDelay})
        
    end
    axis([0 20 0 10])
    
end
% cle
% ar idx TS_ TSsig_ 
%% Distributions of peak t-Scores
col_ = {'r','b','g'};
bins = 0:0.1:10;
figure
for iArea = 1:length(Areas)-1
%     figure('name',Areas{iArea})
    subplot(1,2,iArea); hold on
    title(Areas{iArea})
    for iDelay = 1:length(Delays_)
       TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
       max_ = max(TS_,[],2);
       maxHist_ = cumsum(histc(max_,bins));
       maxHist_ =maxHist_ ./max(maxHist_ );
       stairs(bins,maxHist_,'LineWidth',2,'color',col_{iDelay},'LineStyle','-');
       
       TS_ = cell2mat(eval(sprintf('D_.%s.LR_err.TS_err{iArea}',Delays_{iDelay})))';
       max_ = max(TS_,[],2);
       maxHist_ = cumsum(histc(max_,bins));
       maxHist_ =maxHist_ ./max(maxHist_ );
       stairs(bins,maxHist_,'LineWidth',1,'color',col_{iDelay},'LineStyle','-','HandleVisibility','off');
    end
end
legend(Delays_,'location','SouthEast'); legend boxoff
%% Distributions of time significant
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    title(Areas{iArea})
    for iDelay  = 1:length(Delays_)
       TSsig_   = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%        TSsig_=TSsig_(:,201:402);
       bins     = (2*bw:size(TSsig_,2)).*bw;
       max_     = sum(TSsig_,2).*bw ;
       maxHist_ = cumsum(histc(max_,bins));
       maxHist_ = maxHist_ ./max(maxHist_ );
       stairs(bins,maxHist_,'LineWidth',2,'color',col_{iDelay},'LineStyle','-');
       FracNonSig{iArea}(iDelay) = maxHist_(1);

       
       TSsig_   = eval(sprintf('D_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay}));
       TSsig_(isempty_cell(TSsig_))=[];
       TSsig_ = cell2mat(TSsig_)';
       max_     = sum(TSsig_,2).*bw ;
       maxHist_ = cumsum(histc(max_,bins));
       maxHist_ = maxHist_ ./max(maxHist_ );
       stairs(bins,maxHist_,'LineWidth',2,'color',col_{iDelay},'LineStyle',':','HandleVisibility','off');
       FracNonSig_Err{iArea}(iDelay) = maxHist_(1);

    end
    title(Areas{iArea})
       
    axis([0 10 0 1])
end
legend(Delays_,'location','SouthEast'); legend boxoff
%% Fractions of significant-encoding neurons
figure
for iArea = 1:length(Areas)-1
    subplot(1,2,iArea); hold on
    title(Areas{iArea})
    bar(1-FracNonSig{iArea},'FaceColor','b','FaceAlpha',0.5)
    bar(1-FracNonSig_Err{iArea},'FaceColor','w','FaceAlpha',0.5)
end
%% Plot L/R decoding by area - RT2
figure('name','L/R decoding');
for iArea =1:2%length(Areas)
    subplot(1,2,iArea);hold on
    
    ciplot(mean(D_.Short.LR.Rt2{iArea},2)+nansem(D_.Short.LR.Rt2{iArea},2),...
           mean(D_.Short.LR.Rt2{iArea},2)-nansem(D_.Short.LR.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'b',.9)
    ciplot(mean(D_.Medium.LR.Rt2{iArea},2)+nansem(D_.Medium.LR.Rt2{iArea},2),...
           mean(D_.Medium.LR.Rt2{iArea},2)-nansem(D_.Medium.LR.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'g',0.9)    
    ciplot(mean(D_.Long.LR.Rt2{iArea},2)+nansem(D_.Long.LR.Rt2{iArea},2),...
           mean(D_.Long.LR.Rt2{iArea},2)-nansem(D_.Long.LR.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'r',0.9)               
       
%     plot(tbShort,mean(D_.Short.LR.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.LR.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.LR.Ft2{iArea},2),'b')
    plot([10 10],[0.5 2],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0.1 0.5],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0.1 0.5],'color',[1 0 0 0.3],'LineWidth',1.5)
    if normaliseFscores
        axis([0 20 0 10])
    else
        axis([0 20 0 1])
    end

%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'},'Location','northoutside','Orientation','horizontal')
    legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
     if iArea==1
        ylabel('L/R decoding (F-Score)')
    end
end
%% Plot L/R decoding by area - FT2
figure('color','w','name','L/R decoding');
alpha_ = 0.5;
sigcutoff = 0.5;
x_ = [0:(length(tbAll)*2-1)]*bw;
for iArea =1:2%length(Areas)
    subplot(1,2,iArea);hold on
    ciplot(mean(D_.Short.LR.Ft2{iArea},2)+nansem(D_.Short.LR.Ft2{iArea},2),...
           mean(D_.Short.LR.Ft2{iArea},2)-nansem(D_.Short.LR.Ft2{iArea},2),...
           x_,'b',alpha_)
    ciplot(mean(D_.Medium.LR.Ft2{iArea},2)+nansem(D_.Medium.LR.Ft2{iArea},2),...
           mean(D_.Medium.LR.Ft2{iArea},2)-nansem(D_.Medium.LR.Ft2{iArea},2),...
           x_,'g',alpha_)    
    ciplot(mean(D_.Long.LR.Ft2{iArea},2)+nansem(D_.Long.LR.Ft2{iArea},2),...
           mean(D_.Long.LR.Ft2{iArea},2)-nansem(D_.Long.LR.Ft2{iArea},2),...
           x_,'r',alpha_)           
       
       
       
    temp = nansum(cell2mat(D_.Short.LR.Ft2sig{iArea}),2)./sum(~isnan(sum(cell2mat(D_.Short.LR.Ft2sig{iArea}),1)));
%     temp=temp>=sigcutoff;
    imagesc(repmat(x_,1,2),[-0*ones(size(temp,1),1);-0.2*ones(size(temp,1),1)], [temp,temp]')
    rectangle('Position',[0 -0.3 20 0.4],'LineWidth',1,'EdgeColor','b')

    temp =nansum(cell2mat(D_.Medium.LR.Ft2sig{iArea}),2)./sum(~isnan(sum(cell2mat(D_.Medium.LR.Ft2sig{iArea}),1)));
%     temp=temp>=sigcutoff;
    imagesc(repmat(x_,1,2),[-0.5*ones(size(temp,1),1);-0.7*ones(size(temp,1),1)],[temp,temp]')
    rectangle('Position',[0 -0.8 20 0.4],'LineWidth',1,'EdgeColor','g')

    temp =nansum(cell2mat(D_.Long.LR.Ft2sig{iArea}),2)./sum(~isnan(sum(cell2mat(D_.Long.LR.Ft2sig{iArea}),1)));
%     temp=temp>=sigcutoff;
    imagesc(repmat(x_,1,2),[-1*ones(size(temp,1),1);-1.2*ones(size(temp,1),1)],[temp,temp]')
    rectangle('Position',[0 -1.3 20 0.4],'LineWidth',1,'EdgeColor','r')       

       
%     plot(tbShort,mean(D_.Short.LR.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.LR.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.LR.Ft2{iArea},2),'b')
    plot([10 10],[0.2 0.6],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[0.2 0.6],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[0.2 0.6],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    plot([10 10],[-1.8 -1.4],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[-1.8 -1.4],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[-1.8 -1.4],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    axis([0 20 -2 15])
%     axis([0 20 0 1])
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'},'Location','northoutside','Orientation','horizontal')
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
     if iArea==1
        ylabel('L/R decoding (F-Score)')
        plot([2 7], [5 5],'k','LineWidth',1.5)
        plot([2 2], [5 8],'k','LineWidth',1.5)
    end
    colormap(flipud(gray(10)));
    caxis([0 1])
    axis off
end
%% Plot L/R decoding by area - CVE
figure('name','L/R decoding');
alpha_ = 0.6;
alpha2_ = 0.1;
sigcutoff = 0.2;
x_ = [0:(length(tbAll)*2-1)]*bw;
for iArea =1:2%length(Areas)
    subplot(1,2,iArea);hold on
      
    ciplot((mean(1-D_.Short.LR.CVE{iArea},2)+nansem(1-D_.Short.LR.CVE{iArea},2))*100,...
           (mean(1-D_.Short.LR.CVE{iArea},2)-nansem(1-D_.Short.LR.CVE{iArea},2))*100,...
           x_,'b',alpha_)       
    ciplot((mean(1-D_.Medium.LR.CVE{iArea},2)+nansem(1-D_.Medium.LR.CVE{iArea},2))*100,...
           (mean(1-D_.Medium.LR.CVE{iArea},2)-nansem(1-D_.Medium.LR.CVE{iArea},2))*100,...
           x_,'g',alpha_)               
     ciplot((mean(1-D_.Long.LR.CVE{iArea},2)+nansem(1-D_.Long.LR.CVE{iArea},2))*100,...
           (mean(1-D_.Long.LR.CVE{iArea},2)-nansem(1-D_.Long.LR.CVE{iArea},2))*100,...
           x_,'r',alpha_)   
       
    ciplot((mean(1-D_.Short.LR.CVEbs{iArea},2)+nansem(1-D_.Short.LR.CVEbs{iArea},2))*100,...
           (mean(1-D_.Short.LR.CVEbs{iArea},2)-nansem(1-D_.Short.LR.CVEbs{iArea},2))*100,...
           x_,'b',alpha2_)
    ciplot((mean(1-D_.Medium.LR.CVEbs{iArea},2)+nansem(1-D_.Medium.LR.CVEbs{iArea},2))*100,...
           (mean(1-D_.Medium.LR.CVEbs{iArea},2)-nansem(1-D_.Medium.LR.CVEbs{iArea},2))*100,...
           x_,'g',alpha2_)    
    ciplot((mean(1-D_.Long.LR.CVEbs{iArea},2)+nansem(1-D_.Long.LR.CVEbs{iArea},2))*100,...
           (mean(1-D_.Long.LR.CVEbs{iArea},2)-nansem(1-D_.Long.LR.CVEbs{iArea},2))*100,...
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
       
    temp =nansum(D_.Short.LR.CVEsig{iArea},2)./sum(~isnan(sum(D_.Short.LR.CVEsig{iArea},1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-4.5*ones(size(temp,1),1);-3.5*ones(size(temp,1),1)], 100*[temp,temp]')
    rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
    
    temp = nansum(D_.Medium.LR.CVEsig{iArea},2)./sum(~isnan(sum(D_.Medium.LR.CVEsig{iArea},1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-9.5*ones(size(temp,1),1);-8.5*ones(size(temp,1),1)],100*[temp,temp]')
    rectangle('Position',[0 -10.5 20 3],'LineWidth',1,'EdgeColor','g')

    temp = nansum(D_.Long.LR.CVEsig{iArea},2)./sum(~isnan(sum(D_.Long.LR.CVEsig{iArea},1)));
    temp=temp>sigcutoff;
    imagesc(repmat(x_,1,2),[-14.5*ones(size(temp,1),1);-13.5*ones(size(temp,1),1)],100*[temp,temp]')
    rectangle('Position',[0 -15.5 20 3],'LineWidth',1,'EdgeColor','b')
   
%     plot(tbShort,mean(D_.Short.LR.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.LR.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.LR.Ft2{iArea},2),'b')
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
%% Plot L/R decoding by area - CVE HP vs PFC
figure('name','L/R decoding');hold on
alpha_ = 0.6;
alpha2_ = 0.1;
sigcutoff = 0.2;
x_ = [0:(length(tbAll)*2-1)]*bw;
col__ = {[0 0 0],0.5*[1 1 1]};
col__ = {[0 0 1],[1 0 0]};
plot([10 10],[0 100],'color',[0 0 0],'LineWidth',1.5)
plot([5 5],[0 100],'color',[0 1 0],'LineWidth',1.5)
plot([15 15],[0 100],'color',[1 0 0],'LineWidth',1.5)
plot([0 20],[50 50],':k')

Z = 1-[D_.Short.LR.CVE{2},D_.Medium.LR.CVE{2},D_.Long.LR.CVE{2}];
Z_ =[D_.Short.LR.CVEsig{2},D_.Medium.LR.CVEsig{2},D_.Long.LR.CVEsig{2}];
Z_High = nanmean([D_.Short.LR.cveBSciH{2},D_.Medium.LR.cveBSciH{2},D_.Long.LR.cveBSciH{2}],2);
Z_Low = nanmean([D_.Short.LR.cveBSciL{2},D_.Medium.LR.cveBSciL{2},D_.Long.LR.cveBSciL{2}],2);
ciplot(Z_High*100,Z_Low*100,x_,col__{2},alpha_)
ciplot((mean(Z,2)+nansem(Z,2))*100,(mean(Z,2)-nansem(Z,2))*100,x_,col__{2},alpha_)

Y = 1-[D_.Short.LR.CVE{1},D_.Medium.LR.CVE{1},D_.Long.LR.CVE{1}];
Y_ =[D_.Short.LR.CVEsig{1},D_.Medium.LR.CVEsig{1},D_.Long.LR.CVEsig{1}];
Y_High = nanmean([D_.Short.LR.cveBSciH{1},D_.Medium.LR.cveBSciH{1},D_.Long.LR.cveBSciH{1}],2);
Y_Low = nanmean([D_.Short.LR.cveBSciL{1},D_.Medium.LR.cveBSciL{1},D_.Long.LR.cveBSciL{1}],2);
ciplot(Y_High*100,Y_Low*100,x_,col__{1},alpha_)
ciplot((mean(Y,2)+nansem(Y,2))*100,(mean(Y,2)-nansem(Y,2))*100,x_,col__{1},alpha_)

[sig,crit] = permtest2vec(smooth2a(Y,20,0),smooth2a(Z,20,0),100,0.05);

YZ = (nanmean(Y,2)-nanmean(Z,2))>0;
a = nan(size(sig));a(sig)=1; a(YZ)=NaN;
b = nan(size(sig));b(sig)=1; b(~YZ)=NaN;

% enforce significant individual decoding from each area
% Ysig = nansum(Y_,2)./size(Y_,2)>0.5 ;
% Zsig = nansum(Z_,2)./size(Z_,2)>0.5 ;
% YZsig = nansum(Y_,2)./size(Y_,2)>0.5 & nansum(Z_,2)./size(Z_,2)>0.5 ;
% a(~Ysig)=NaN;
% b(~Zsig)=NaN;
% a(~YZsig)=NaN;
% b(~YZsig)=NaN;
%         plot(x_,repmat(105,length(x_),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)
plot(x_,a*102,'color',[col__{2},alpha_],'LineWidth',5)
plot(x_,b*105,'color',[col__{1},alpha_],'LineWidth',5)

axis([0 20 40 110])
set(gca,'YTick',[0:25:100])

ylabel('% Correct decoding')
xlabel('Time (s)')
%% Plot L/R decoding by area - CVE HP vs PFC redux ***
figure('name','L/R decoding');hold on
alpha_ = 0.6;
alpha2_ = 0.1;
sigcutoff = 0.2;
x_ = [0:(length(tbAll)*2-1)]*bw;
col__ = {[0 0 0],0.5*[1 1 1]};
col__ = {[0 0 1],[1 0 0]};
plot([10 10],[0 100],'color',[0 0 0],'LineWidth',1.5)
plot([5 5],[0 100],'color',[0 1 0],'LineWidth',1.5)
plot([15 15],[0 100],'color',[1 0 0],'LineWidth',1.5)
plot([0 20],[50 50],':k','LineWidth',1.5)

Z = [(1-D_.Short.LR.CVE{2})+(1-D_.Medium.LR.CVE{2})+(1-D_.Long.LR.CVE{2})]/3;
Z_ =[nansum(D_.Short.LR.CVEsig{2},2)./size(D_.Short.LR.CVEsig{2},2) + ...
     nansum(D_.Medium.LR.CVEsig{2},2)./size(D_.Medium.LR.CVEsig{2},2) + ...
     nansum(D_.Long.LR.CVEsig{2},2)./size(D_.Long.LR.CVEsig{2},2)] ./3;
% Z_High = nanmean([D_.Short.LR.cveBSciH{2},D_.Medium.LR.cveBSciH{2},D_.Long.LR.cveBSciH{2}],2);
% Z_Low = nanmean([D_.Short.LR.cveBSciL{2},D_.Medium.LR.cveBSciL{2},D_.Long.LR.cveBSciL{2}],2);
% ciplot(Z_High*100,Z_Low*100,x_,col__{2},alpha_)
ciplot((mean(Z,2)+nansem(Z,2))*100,(mean(Z,2)-nansem(Z,2))*100,x_,col__{2},alpha_);

Y = [(1-D_.Short.LR.CVE{1})+(1-D_.Medium.LR.CVE{1})+(1-D_.Long.LR.CVE{1})]/3;
Y_ =[nansum(D_.Short.LR.CVEsig{1},2)./size(D_.Short.LR.CVEsig{1},2) + ...
     nansum(D_.Medium.LR.CVEsig{1},2)./size(D_.Medium.LR.CVEsig{1},2) + ...
     nansum(D_.Long.LR.CVEsig{1},2)./size(D_.Long.LR.CVEsig{1},2)] ./3;
% Y_High = nanmean([D_.Short.LR.cveBSciH{1},D_.Medium.LR.cveBSciH{1},D_.Long.LR.cveBSciH{1}],2);
% Y_Low = nanmean([D_.Short.LR.cveBSciL{1},D_.Medium.LR.cveBSciL{1},D_.Long.LR.cveBSciL{1}],2);
% ciplot(Y_High*100,Y_Low*100,x_,col__{1},alpha_)
ciplot((mean(Y,2)+nansem(Y,2))*100,(mean(Y,2)-nansem(Y,2))*100,x_,col__{1},alpha_);
[sig,crit] = permtest2vec(smooth2a(Y,20,0),smooth2a(Z,20,0),50,0.05);

YZ = (nanmean(Y,2)-nanmean(Z,2))>0;
a = nan(size(sig));a(sig)=1; a(YZ)=NaN;
b = nan(size(sig));b(sig)=1; b(~YZ)=NaN;
c = ones(size(sig)); c(~isnan(a)|~isnan(b))=0; c(c==0)=nan;
% enforce significant individual decoding from each area
% Ysig = nansum(Y_,2)./size(Y_,2)>0.5 ;
% Zsig = nansum(Z_,2)./size(Z_,2)>0.5 ;
% YZsig = nansum(Y_,2)./size(Y_,2)>0.5 & nansum(Z_,2)./size(Z_,2)>0.5 ;
% a(~Ysig)=NaN;
% b(~Zsig)=NaN;
% a(~YZsig)=NaN;
% b(~YZsig)=NaN;
%         plot(x_,repmat(105,length(x_),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)

imagesc(repmat(x_,1,4),repmat(101*ones(1,Ltr),1,4),repmat(Y_,1,4)')

imagesc(repmat(x_,1,4),repmat(106*ones(1,Ltr),1,4),repmat(Z_,1,4)')

rectangle('Position',[0,100.5,max(x_),4],'EdgeColor',[col__{1},alpha_],'LineWidth',2 )
rectangle('Position',[0,105.5,max(x_),4],'EdgeColor',[col__{2},alpha_],'LineWidth',2 )
% imagesc(x_,104*ones(1,Ltr),Z_')
colormap(flipud(gray))
caxis([0 1])
plot(x_,a*113,'color',[col__{2},alpha_],'LineWidth',10)
plot(x_,b*113,'color',[col__{1},alpha_],'LineWidth',10)
plot(x_,c*113,'color',[0.9*ones(1,3),alpha_],'LineWidth',10)

axis([-0.1 20.1 40 115])
set(gca,'YTick',[0:25:100])

ylabel('% Correct decoding')
xlabel('Time (s)')
%% Plot L/R decoding by delay
figure('name','L/R decoding');
for iDelay = 1:length(Delays_)
    subplot(1,length(Delays_),iDelay);hold on
%     title(['L/R trial decoding ' Delays_{iDelay} ' Trials'])
    title([Delays_{iDelay} ' trials'])
    
    for iArea = 1:2%length(Areas)
        eval(sprintf('m=nanmean(D_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
        eval(sprintf('e=nansem(D_.%s.LR.Rt2{iArea},2);',Delays_{iDelay}));
        ciplot(m+e,m-e,[0:(length(tbAll)*2-1)]*bw,color_{iArea})
    end
    plot([10 10],[0.1 0.5],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0.1 0.5],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0.1 0.5],'color',[1 0 0 0.3],'LineWidth',1.5)
    axis([0 20 0 10])
    if iDelay==1
        ylabel('L/R decoding (F-Score)')
    end
end
legend(Areas{1:2})
%% Plot L/R decoding on correct and errors

for iArea =1:length(Areas)
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
    title(Delays_{iDelay})   
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
%% Plot L/R decoding on correct and errors - CVE redux
figure('name','L/R decoding');hold on
alpha_ = 0.6;
alpha2_ = 0.1;
sigcutoff = 0.2;
x_ = [0:(length(tbAll)*2-1)]*bw;
col__ = {[0 0 0],0.5*[1 1 1]};
col__ = {[0 0 1],[1 0 0]};
plot([10 10],[0 100],'color',[0 0 0],'LineWidth',1.5)
plot([5 5],[0 100],'color',[0 1 0],'LineWidth',1.5)
plot([15 15],[0 100],'color',[1 0 0],'LineWidth',1.5)
plot([0 20],[50 50],':k','LineWidth',1.5)

Z = [(1-D_.Short.LR.CVE{2})+(1-D_.Medium.LR.CVE{2})+(1-D_.Long.LR.CVE{2})]/3;
Z_ =[nansum(D_.Short.LR.CVEsig{2},2)./size(D_.Short.LR.CVEsig{2},2) + ...
     nansum(D_.Medium.LR.CVEsig{2},2)./size(D_.Medium.LR.CVEsig{2},2) + ...
     nansum(D_.Long.LR.CVEsig{2},2)./size(D_.Long.LR.CVEsig{2},2)] ./3;
% Z_High = nanmean([D_.Short.LR.cveBSciH{2},D_.Medium.LR.cveBSciH{2},D_.Long.LR.cveBSciH{2}],2);
% Z_Low = nanmean([D_.Short.LR.cveBSciL{2},D_.Medium.LR.cveBSciL{2},D_.Long.LR.cveBSciL{2}],2);
% ciplot(Z_High*100,Z_Low*100,x_,col__{2},alpha_)
ciplot((mean(Z,2)+nansem(Z,2))*100,(mean(Z,2)-nansem(Z,2))*100,x_,col__{2},alpha_);

Y = [(1-D_.Short.LR.CVE{1})+(1-D_.Medium.LR.CVE{1})+(1-D_.Long.LR.CVE{1})]/3;
Y_ =[nansum(D_.Short.LR.CVEsig{1},2)./size(D_.Short.LR.CVEsig{1},2) + ...
     nansum(D_.Medium.LR.CVEsig{1},2)./size(D_.Medium.LR.CVEsig{1},2) + ...
     nansum(D_.Long.LR.CVEsig{1},2)./size(D_.Long.LR.CVEsig{1},2)] ./3;
% Y_High = nanmean([D_.Short.LR.cveBSciH{1},D_.Medium.LR.cveBSciH{1},D_.Long.LR.cveBSciH{1}],2);
% Y_Low = nanmean([D_.Short.LR.cveBSciL{1},D_.Medium.LR.cveBSciL{1},D_.Long.LR.cveBSciL{1}],2);
% ciplot(Y_High*100,Y_Low*100,x_,col__{1},alpha_)
ciplot((mean(Y,2)+nansem(Y,2))*100,(mean(Y,2)-nansem(Y,2))*100,x_,col__{1},alpha_);
[sig,crit] = permtest2vec(smooth2a(Y,20,0),smooth2a(Z,20,0),50,0.05);

YZ = (nanmean(Y,2)-nanmean(Z,2))>0;
a = nan(size(sig));a(sig)=1; a(YZ)=NaN;
b = nan(size(sig));b(sig)=1; b(~YZ)=NaN;
c = ones(size(sig)); c(~isnan(a)|~isnan(b))=0; c(c==0)=nan;
% enforce significant individual decoding from each area
% Ysig = nansum(Y_,2)./size(Y_,2)>0.5 ;
% Zsig = nansum(Z_,2)./size(Z_,2)>0.5 ;
% YZsig = nansum(Y_,2)./size(Y_,2)>0.5 & nansum(Z_,2)./size(Z_,2)>0.5 ;
% a(~Ysig)=NaN;
% b(~Zsig)=NaN;
% a(~YZsig)=NaN;
% b(~YZsig)=NaN;
%         plot(x_,repmat(105,length(x_),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)

imagesc(repmat(x_,1,4),repmat(101*ones(1,Ltr),1,4),repmat(Y_,1,4)')

imagesc(repmat(x_,1,4),repmat(106*ones(1,Ltr),1,4),repmat(Z_,1,4)')

rectangle('Position',[0,100.5,max(x_),4],'EdgeColor',[col__{1},alpha_],'LineWidth',2 )
rectangle('Position',[0,105.5,max(x_),4],'EdgeColor',[col__{2},alpha_],'LineWidth',2 )
% imagesc(x_,104*ones(1,Ltr),Z_')
colormap(flipud(gray))
caxis([0 1])
plot(x_,a*113,'color',[col__{2},alpha_],'LineWidth',10)
plot(x_,b*113,'color',[col__{1},alpha_],'LineWidth',10)
plot(x_,c*113,'color',[0.9*ones(1,3),alpha_],'LineWidth',10)

axis([-0.1 20.1 40 115])
set(gca,'YTick',[0:25:100])

ylabel('% Correct decoding')
xlabel('Time (s)')
%% Plot L/R decoding by area - CVE C vs PFC redux ***
figure('name','L/R decoding');hold on
alpha_ = 0.6;
alpha2_ = 0.1;
sigcutoff = 0.2;
x_ = [0:(length(tbAll)*2-1)]*bw;
col_ = {[0 0 1],[1 0 0]};
Ltr = length(x_);

for iArea =1:2
    subplot(1,2,iArea); hold on
    plot([5 5],[0 100],'color',[0 1 0],'LineWidth',1.5)
    plot([15 15],[0 100],'color',[1 0 0],'LineWidth',1.5)
    plot([0 20],[50 50],':k','LineWidth',1.5)

%     Z = [(1-D_.Short.LR.CVE{iArea})+(1-D_.Medium.LR.CVE{iArea})+(1-D_.Long.LR.CVE{iArea})]/3;
    Z = cat(3,1-D_.Short.LR.CVE{iArea} , ...
         1-D_.Medium.LR.CVE{iArea}, ...
         1-D_.Long.LR.CVE{iArea});
     Z = nanmean(Z,3);
    Z = smooth2a(Z,5,0);
%     Z_ =[nansum(D_.Short.LR.CVEsig{iArea},2)./size(D_.Short.LR.CVEsig{iArea},2) + ...
%          nansum(D_.Medium.LR.CVEsig{iArea},2)./size(D_.Medium.LR.CVEsig{iArea},2) + ...
%          nansum(D_.Long.LR.CVEsig{iArea},2)./size(D_.Long.LR.CVEsig{iArea},2)] ./3;

    ciplot((nanmean(Z,2)+nansem(Z,2))*100,(nanmean(Z,2)-nansem(Z,2))*100,x_,col_{iArea},1);
    
    Y = cat(3,1-D_.Short.LR_err.CVE{iArea} , ...
         1-D_.Medium.LR_err.CVE{iArea}, ...
         1-D_.Long.LR_err.CVE{iArea});
     Y = nanmean(Y,3);
     Y = smooth2a(Y,5,0);
%     Y_ =[nansum(D_.Short.LR_err.CVEsig{iArea},2)./size(D_.Short.LR_err.CVEsig{iArea},2) + ...
%          nansum(D_.Medium.LR_err.CVEsig{iArea},2)./size(D_.Medium.LR_err.CVEsig{iArea},2) + ...
%          nansum(D_.Long.LR_err.CVEsig{iArea},2)./size(D_.Long.LR_err.CVEsig{iArea},2)] ./3;

    ciplot((nanmean(Y,2)+nansem(Y,2))*100,(nanmean(Y,2)-nansem(Y,2))*100,x_,col_{iArea},0.5);
    plot([10 10],[0 100],'color','w','LineWidth',5)
    
      [sig,crit] = permtest2vec(smooth2a(Y,20,0),smooth2a(Z,20,0),50,0.05);

    YZ = (nanmean(Y,2)-nanmean(Z,2))>0;
    a = nan(size(sig));a(sig)=1; a(YZ)=NaN;
    b = nan(size(sig));b(sig)=1; b(~YZ)=NaN;
    c = ones(size(sig)); c(~isnan(a)|~isnan(b))=0; c(c==0)=nan;
    % enforce significant individual decoding from each area
    % Ysig = nansum(Y_,2)./size(Y_,2)>0.5 ;
    % Zsig = nansum(Z_,2)./size(Z_,2)>0.5 ;
    % YZsig = nansum(Y_,2)./size(Y_,2)>0.5 & nansum(Z_,2)./size(Z_,2)>0.5 ;
    % a(~Ysig)=NaN;
    % b(~Zsig)=NaN;
    % a(~YZsig)=NaN;
    % b(~YZsig)=NaN;
    %         plot(x_,repmat(105,length(x_),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)

%     imagesc(repmat(x_,1,4),repmat(101*ones(1,Ltr),1,4),repmat(Y_,1,4)')
% 
%     imagesc(repmat(x_,1,4),repmat(106*ones(1,Ltr),1,4),repmat(Z_,1,4)')

%     rectangle('Position',[0,100.5,max(x_),4],'EdgeColor',[col__{1},alpha_],'LineWidth',2 )
%     rectangle('Position',[0,105.5,max(x_),4],'EdgeColor',[col__{2},alpha_],'LineWidth',2 )
    % imagesc(x_,104*ones(1,Ltr),Z_')
    colormap(flipud(gray))
    caxis([0 1])
    plot(x_,a*103,'color','k','LineWidth',10)
    plot(x_,b*103,'color','k','LineWidth',10)
    plot(x_,c*103,'color',[0.9*ones(1,3),alpha_],'LineWidth',10)

    axis([-0.1 20.1 40 115])
    set(gca,'YTick',[0:25:100])

    ylabel('% Correct decoding')
    xlabel('Time (s)')
end
%%
  
  


%% Plot S/C decoding

figure('name','S/C decoding');
for iArea =1:length(Areas)
    subplot(1,length(Areas),iArea);hold on
    ciplot(mean(D_.Short.SC.Rt2{iArea},2)+nansem(D_.Short.SC.Rt2{iArea},2),...
           mean(D_.Short.SC.Rt2{iArea},2)-nansem(D_.Short.SC.Rt2{iArea},2),...
           tbShort,'r')
    ciplot(mean(D_.Medium.SC.Rt2{iArea},2)+nansem(D_.Medium.SC.Rt2{iArea},2),...
           mean(D_.Medium.SC.Rt2{iArea},2)-nansem(D_.Medium.SC.Rt2{iArea},2),...
           tbMedium,'g')    
    ciplot(mean(D_.Long.SC.Rt2{iArea},2)+nansem(D_.Long.SC.Rt2{iArea},2),...
           mean(D_.Long.SC.Rt2{iArea},2)-nansem(D_.Long.SC.Rt2{iArea},2),...
           tbLong,'b')           
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
%% batch import Assem data for meta-analysis and plotting - decoders only (OLD)
for iFile = 1:length(fileListAss)
        
        fname=strtok(fileListAss(iFile).name,'_');
%         fnIn = sprintf('%s\\MixedSelectivity\\%s_MixedSelectivity_Ass.mat',pat,fname);
        fnIn = sprintf('%s%sMixedSelectivity_%s%sLeftvsRight%s%s_MixedSelectivity_Ass.mat',pat,filesep,Target,filesep,filesep,fname);

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
%                 flag = 0;
%                 if eval(sprintf('isfield(D_Ass.%s,''LR'')',Delays_{iDelay})) 
%                     if eval(sprintf('~isempty(D_Ass.%s.LR{iArea})',Delays_{iDelay}))
%                         flag = 1;
%                         eval(sprintf('D_temp = D_Ass.%s.LR{iArea};',Delays_{iDelay}))
%                         if normaliseFscores
%                             bp = find(tbAll>normWin(1) & tbAll<normWin(2));
%                             B  = nanmean(D_temp.Rt2(:,bp)');
%                             D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
%                         end
%                         eval(sprintf('D_Ass_.%s.LR.Rt2{iArea}(:,iFile)       = D_temp.Rt2;',Delays_{iDelay}));
%                         
%                          % L/R: Continuous decoding of Positional information - correct trials
%                          
%                          eval(sprintf('D_Ass_.%s.LR.CVE{iArea}(:,iFile)      = D_temp.CVE;',Delays_{iDelay}))
%                          
%                         eval(sprintf('D_Ass_.%s.LR.CVEsig{iArea}(:,iFile)   = nanmean(D_temp.CVE,1)<D_temp.cveBSciL;',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
%                         eval(sprintf('D_Ass_.%s.LR.cveBSciH{iArea}(:,iFile) = D_temp.cveBSciH;',Delays_{iDelay})); 
%                         eval(sprintf('D_Ass_.%s.LR.cveBSciL{iArea}(:,iFile) = D_temp.cveBSciL;',Delays_{iDelay}));
%                         eval(sprintf('D_Ass_.%s.LR.CVEbs{iArea}(:,iFile)    = D_temp.CVEbs;',Delays_{iDelay}));
% 
%                         
%                         
%                         
%                         clear D_temp
%                     end
%                 end
% %                 if ~flag
% %                 	eval(sprintf(   'D_Ass_.%s.LR.Rt2{iArea}(:,iFile)       = nan(2*length(tbAll),1);',Delays_{iDelay}));
% %                 end        
                
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
% clear D_Ass  Stability D_temp D_tempL D_tempR B bp 
%% batch import assembly data for meta-analysis and plotting - Decoders only (NEW)
% Import
clear D_Ass

for iFile =1:length(fileListAss)
    
    fname=strtok(fileListAss(iFile).name,'_');
    
    fnIn = sprintf('%s%sMixedSelectivity_%s%sLeftvsRight%s%s_MixedSelectivity_Ass.mat',pat,filesep,Target,filesep,filesep,fname);
    
    load(fnIn ,'D_Ass');
    
    for iDelay = 1:length(Delays_)
        for iArea = 1:length(Areas)
            
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
                        
                        B  = nanmean(D_temp.Ft2(:,bp)');
                        D_temp.Ft2=D_temp.Ft2./(B'*ones(1,length(D_temp.Ft2)));
                        
                    end
                    
                    eval(sprintf('D_Ass_.%s.LR.TS{iArea}{iFile}      = D_temp.TS;',Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.LR.TSsig{iArea}{iFile}   = D_temp.TSsig;',Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.LR.Rt2{iArea}(:,iFile)   = D_temp.Rt2;',Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.LR.Ft2{iArea}(:,iFile)   = D_temp.Ft2;',Delays_{iDelay}))
                    eval(sprintf('D_Ass_.%s.LR.Ft2sig{iArea}{iFile}  = [D_Ass.%s.LR{iArea}.Ft2>D_Ass.%s.LR{iArea}.Ft2ciH]'';',Delays_{iDelay},Delays_{iDelay},Delays_{iDelay}))
                    
                    
                    eval(sprintf('D_Ass_.%s.LR.CVE{iArea}(:,iFile)      = D_temp.CVE;',Delays_{iDelay}))
                    
                    eval(sprintf('D_Ass_.%s.LR.CVEsig{iArea}(:,iFile)   = nanmean(D_temp.CVE,1)<D_temp.cveBSciL;',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
                    eval(sprintf('D_Ass_.%s.LR.cveBSciH{iArea}(:,iFile) = D_temp.cveBSciH;',Delays_{iDelay}));
                    eval(sprintf('D_Ass_.%s.LR.cveBSciL{iArea}(:,iFile) = D_temp.cveBSciL;',Delays_{iDelay}));
                    eval(sprintf('D_Ass_.%s.LR.CVEbs{iArea}(:,iFile)    = D_temp.CVEbs;',Delays_{iDelay}));
                    
                    clear D_temp
                    
                end
            end
            if ~flag
                    eval(sprintf(   'D_Ass_.%s.LR.TS{iArea}{iFile}          = nan(2*length(tbAll),1);',Delays_{iDelay}));
                    eval(sprintf(   'D_Ass_.%s.LR.TSsig{iArea}{iFile}       = false(2*length(tbAll),1);',Delays_{iDelay}));
                    eval(sprintf(   'D_Ass_.%s.LR.Rt2{iArea}(:,iFile)         = nan(2*length(tbAll),1);',Delays_{iDelay}));
                    eval(sprintf(   'D_Ass_.%s.LR.Ft2{iArea}(:,iFile)         = nan(2*length(tbAll),1);',Delays_{iDelay}));
                    eval(sprintf(   'D_Ass_.%s.LR.Ft2sig{iArea}{iFile}      = false(2*length(tbAll),1);',Delays_{iDelay}));
            end
            
            
            % L/R: Continuous decoding of Positional information - error trials
            %             try
            %                 eval(sprintf('D_temp = D_Ass.%s.LR_Err{iArea};',Delays_{iDelay}))
            %                 if normaliseFscores
            %                     bp = find(tbAll>normWin(1) & tbAll<-normWin(2));
            %                     B  = nanmean(D_temp.Rt2(:,bp)');
            %                     D_temp.Rt2=D_temp.Rt2./(B'*ones(1,length(D_temp.Rt2)));
            %
            %                     B  = nanmean(D_temp.Ft2(:,bp)');
            %                     D_temp.Ft2=D_temp.Ft2./(B'*ones(1,length(D_temp.Ft2)));
            %                 end
            %                 eval(sprintf('D_Ass_.%s.LR_err.TS_err{iArea}{iFile}        = D_temp.TS;',Delays_{iDelay}))
            %                 eval(sprintf('D_Ass_.%s.LR_err.TSsig_err{iArea}{iFile}     = D_temp.TSsig;',Delays_{iDelay}))
            %                 eval(sprintf('D_Ass_.%s.LR_err.Rt2{iArea}(:,iFile)         = D_temp.Rt2;',Delays_{iDelay})),
            %                 eval(sprintf('D_Ass_.%s.LR_err.Ft2{iArea}(:,iFile)         = D_temp.Ft2;',Delays_{iDelay})),
            %                 if size(D_temp.CVE,2) ==2*length(tbAll)
            %                     eval(sprintf('D_Ass_.%s.LR_err.CVE{iArea}(:,iFile)     = nanmean(D_temp.CVE,1);',Delays_{iDelay})),
            %                 else
            %                     eval(sprintf('D_Ass_.%s.LR_err.CVE{iArea}(:,iFile)     = nan(2*length(tbAll),1);',Delays_{iDelay}))
            %                 end
            %
            %                 eval(sprintf('D_Ass_.%s.LR_err.CVEsig{iArea}(:,iFile)  = D_temp.cveBSciH>nanmean(D_temp.CVE,1);',Delays_{iDelay})); % Significant decoding (n.b. sign is inverted duer to 1-CVE)
            %             catch
            %                 eval(sprintf('D_Ass_.%s.LR_err.Rt2{iArea}(:,iFile)     = nan(2*length(tbAll),1);',Delays_{iDelay}))
            %                 eval(sprintf('D_Ass_.%s.LR_err.CVE{iArea}(:,iFile)     = nan(2*length(tbAll),1);',Delays_{iDelay}))
            %                 eval(sprintf('D_Ass_.%s.LR_err.CVEsig{iArea}(:,iFile)  = nan(2*length(tbAll),1);',Delays_{iDelay}))
            %                 eval(sprintf('D_Ass_.%s.LR_err.CVEbs{iArea}(:,iFile)   = nan(2*length(tbAll),1);',Delays_{iDelay}))
            %             end
            %             clear D_temp
            
            
        end
        
        clear D
    end
    
end


clear D Stability D_temp D_tempL D_tempR B bp
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
axis([0 20 0 10])
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
clear idx
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
%% Plot L/R decoding by area - CVE HP vs PFC vs Joint
figure('name','L/R decoding');hold on
alpha_ = 1;
alpha2_ = 0.1;
sigcutoff = 0.2;
x_ = [0:(length(tbAll)*2-1)]*bw;
col__ = {[0 0 0],0.5*[1 1 1],[0 1 0]};

plot([10 10],[0 100],'color',[0 0 0],'LineWidth',1.5)
plot([5 5],[0 100],'color',[0 1 0],'LineWidth',1.5)
plot([15 15],[0 100],'color',[1 0 0],'LineWidth',1.5)
plot([0 20],[50 50],':k')

Z = 1-[D_Ass_.Short.LR.CVE{2},D_Ass_.Medium.LR.CVE{2},D_Ass_.Long.LR.CVE{2}];
Z_ =[D_Ass_.Short.LR.CVEsig{2},D_Ass_.Medium.LR.CVEsig{2},D_Ass_.Long.LR.CVEsig{2}];
idx = sum(Z) == size(Z,1);
Z(:,idx)=[];
Z_(:,idx)=[];
% Z_High = nanmean([D_.Short.LR.cveBSciH{2},D_.Medium.LR.cveBSciH{2},D_.Long.LR.cveBSciH{2}],2);
% Z_Low = nanmean([D_.Short.LR.cveBSciL{2},D_.Medium.LR.cveBSciL{2},D_.Long.LR.cveBSciL{2}],2);
% ciplot(Z_High*100,Z_Low*100,x_,col__{2},alpha_)
ciplot((mean(Z,2)+nansem(Z,2))*100,(mean(Z,2)-nansem(Z,2))*100,x_,col__{2},alpha_)
% plot(x_,Z*100,'color',col__{2})

Y = 1-[D_Ass_.Short.LR.CVE{1},D_Ass_.Medium.LR.CVE{1},D_Ass_.Long.LR.CVE{1}];
Y_ =[D_Ass_.Short.LR.CVEsig{1},D_Ass_.Medium.LR.CVEsig{1},D_Ass_.Long.LR.CVEsig{1}];
idx = sum(Y) == size(Y,1);
Y(:,idx)=[];
Y_(:,idx)=[];
% Y_High = nanmean([D_.Short.LR.cveBSciH{1},D_.Medium.LR.cveBSciH{1},D_.Long.LR.cveBSciH{1}],2);
% Y_Low = nanmean([D_.Short.LR.cveBSciL{1},D_.Medium.LR.cveBSciL{1},D_.Long.LR.cveBSciL{1}],2);
% ciplot(Y_High*100,Y_Low*100,x_,col__{1},alpha_)
ciplot((mean(Y,2)+nansem(Y,2))*100,(mean(Y,2)-nansem(Y,2))*100,x_,col__{1},alpha_)
% plot(x_,Y*100,'color',col__{1})

X = 1-[D_Ass_.Short.LR.CVE{3},D_Ass_.Medium.LR.CVE{3},D_Ass_.Long.LR.CVE{3}];
X_ =[D_Ass_.Short.LR.CVEsig{3},D_Ass_.Medium.LR.CVEsig{3},D_Ass_.Long.LR.CVEsig{3}];
idx = sum(X) == size(X,1);
X(:,idx)=[];
X_(:,idx)=[];
% Y_High = nanmean([D_.Short.LR.cveBSciH{1},D_.Medium.LR.cveBSciH{1},D_.Long.LR.cveBSciH{1}],2);
% Y_Low = nanmean([D_.Short.LR.cveBSciL{1},D_.Medium.LR.cveBSciL{1},D_.Long.LR.cveBSciL{1}],2);
% ciplot(Y_High*100,Y_Low*100,x_,col__{1},alpha_)
ciplot((mean(X,2)+nansem(X,2))*100,(mean(X,2)-nansem(X,2))*100,x_,col__{3},alpha_)
% plot(x_,X*100,'color',col__{3})


[sig,crit] = permtest2vec(smooth2a(Y,20,0),smooth2a(Z,20,0),100,0.05);

YZ = (nanmean(Y,2)-nanmean(Z,2))>0;
a = nan(size(sig));a(sig)=1; a(YZ)=NaN;
b = nan(size(sig));b(sig)=1; b(~YZ)=NaN;

% enforce significant individual decoding from each area
% Ysig = nansum(Y_,2)./size(Y_,2)>0.5 ;
% Zsig = nansum(Z_,2)./size(Z_,2)>0.5 ;
% YZsig = nansum(Y_,2)./size(Y_,2)>0.5 & nansum(Z_,2)./size(Z_,2)>0.5 ;
% a(~Ysig)=NaN;

% b(~Zsig)=NaN;
% a(~YZsig)=NaN;
% b(~YZsig)=NaN;
%         plot(x_,repmat(105,length(x_),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',3)
plot(x_,a*105,'color',col__{2},'LineWidth',5)
plot(x_,b*105,'color',col__{1},'LineWidth',5)

axis([0 20 0 110])
set(gca,'YTick',[0:25:100])

ylabel('% Correct decoding')
xlabel('Time (s)')


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

%% Get t-scores, plot side by side

for iArea = 1:length(Areas)
    figure
    
%     iDelay = 1;
%     TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
%     TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%             [~,idx] = max(TS_,[],2); 
%             [~,idx] = sort(idx);
            
            
%     TS_ = cell2mat(eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay})))';
%     TSsig_ = cell2mat(eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
%     idxNS = find(nansum(TSsig_,2)==0);
%     [~,idx] = max(TS_,[],2); 
%     [~,idx] = sort(idx);
%     idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];

    for iDelay = 1:length(Delays_)
        TS_ = cell2mat(eval(sprintf('D_Ass_.%s.LR.TS{iArea}',Delays_{iDelay})))';
        TSsig_ = cell2mat(eval(sprintf('D_Ass_.%s.LR.TSsig{iArea}',Delays_{iDelay})))';
        TS(sum(TSsig_,2)==0,:)=0;
        
        % Sort by peak time
%         if iDelay ==1
            [~,idx] = max(TS_,[],2); 
            [~,idx] = sort(idx);
%         end
    % Sort by peak time, separate insignificant units
%     if iDelay ==1
%     idxNS = find(nansum(TSsig_,2)==0);
%     [~,idx] = max(TS_,[],2); 
%     [~,idx] = sort(idx);
%     idx = [idx(ismember(idx,idxNS)); idx(~ismember(idx,idxNS))];
%     end
    
    TS_ =TS_(idx,:);
    TSsig_ =TSsig_ (idx,:);
    TS_(~TSsig_ )=0;
    subplot(1,3,iDelay); hold on
    x_  = repmat((1:size(TS_,2)).*bw,size(TS_,1),1);
    y_  = repmat(1:size(TS_,1),size(TS_,2),1);
    imagesc(x_(:),y_(:),TS_);
    set(gca,'YDir','normal')
%     cmap =([1 1 1;flipud(hot)]);
    cmap =([1 1 1;(jet)]);
    colormap (cmap)
    caxis([0 10])
    plot([10 10],[1 size(TS_,1)],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[1 size(TS_,1)],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[1 size(TS_,1)],'color',[1 0 0 0.6],'LineWidth',1.5)
    axis([0 20 0 inf])
    if iArea==1 & iDelay ==3
        plot([18 20],[20 20],'k','LineWidth',1.5)
        plot([18 18],[20 70],'k','LineWidth',1.5)
    end
    axis off
    
    end
end
% clear idx TS_ TSsig_
%% Distributions of peak t-Scores - shaded
plotErr = false;
bins = 0:0.1:10;  
alpha_ = 0.3;
alpha_Err = 0.3;
figure
for iArea = 1:length(Areas)
    subplot(1,3,iArea); hold on
    for iDelay = 1:length(Delays_)
       TS_ = eval(sprintf('D_Ass_.%s.LR.TS{iArea}',Delays_{iDelay}));
       TS_ = cellfun(@(x) max(x,[],1),TS_,'UniformOutput',false);
       for iFile =1:length(TS_)
           maxhist = cumsum(histc(TS_{iFile},bins));
           maxHist(:,iFile) = maxhist./max(maxhist);
       end
       maxHist_Mean = nanmean(maxHist,2);
       maxHist_SEM = nansem(maxHist,2);
       ciplot(maxHist_Mean+maxHist_SEM,maxHist_Mean-maxHist_SEM,...
              bins,color_{iDelay},alpha_)
          
          if plotErr
       stairs(bins,maxHist_Mean,'LineWidth',1.5,'color',color_{iDelay});

       TS_Err = eval(sprintf('D_Ass_.%s.LR_err.TS_err{iArea}',Delays_{iDelay}));
       TS_Err(isempty_cell(TS_Err))=[];
       TS_Err = cellfun(@(x) max(x,[],1),TS_Err,'UniformOutput',false);
       for iFile =1:length(TS_Err)
           try
              maxhist = cumsum(histc(TS_Err{iFile},bins));
              maxHist(:,iFile) = maxhist./max(maxhist);
           end
       end
       maxHist_Mean = nanmean(maxHist,2);
       maxHist_SEM  = nansem(maxHist,2);
%        ciplot(maxHist_Mean+maxHist_SEM,maxHist_Mean-maxHist_SEM,...
%               bins,color_{iDelay},alpha_Err)
       stairs(bins,maxHist_Mean,'LineWidth',1.5,'LineStyle',':','color',color_{iDelay});
          end

        set(gca,'Ytick',0:0.5:1,'Xtick',0:5:20)
        xlabel('Peak t-score')
        ylabel('Fraction of Aassemblies')
        title(Areas{iArea})
    end
end
%% Distributions of time significant - shaded
plotErr = false;
figure
clear timehist timeHist FracSig FracSigErr
Sigthresh = bw;
maxBin = 15;
bins   = Sigthresh:bw:maxBin;
alpha_ = 0.3;
alpha_Err = 0.2;

for iArea = 1:length(Areas)
    subplot(1,3,iArea); hold on
    for iDelay  = 1:length(Delays_)
       TSsig_   = eval(sprintf('D_Ass_.%s.LR.TSsig{iArea}',Delays_{iDelay}));
       
%        TSsig_ = cellfun(@(x) sum(x(202:end,:),1).*bw,TSsig_,'UniformOutput',false);   
       TSsig_ = cellfun(@(x) sum(x,1).*bw,TSsig_,'UniformOutput',false);
       FracSig{iArea}(:,iDelay) = cell2mat(cellfun(@(x) sum(sum(x,1).*bw > Sigthresh)./size(x,2),TSsig_,'UniformOutput',false));
       
        for iFile =1:length(TSsig_)

           timehist = cumsum(histc(TSsig_{iFile},bins));
           timeHist(:,iFile) = timehist./max(timehist);
       end
       
       timeHist_Mean = nanmean(timeHist,2);
       timeHist_SEM = nansem(timeHist,2);
       ciplot(timeHist_Mean+timeHist_SEM,timeHist_Mean-timeHist_SEM,...
              bins,color_{iDelay},alpha_)
%        plot(bins,timeHist_Mean,'Color',color_{iDelay},'LineWidth',1.5)   
if plotErr
    TSsig_Err   = eval(sprintf('D_Ass_.%s.LR_err.TSsig_err{iArea}',Delays_{iDelay}));
    emptyYN = isempty_cell(TSsig_Err); 
    TSsig_Err = cellfun(@(x) sum(x,1).*bw,TSsig_Err,'UniformOutput',false);

    for iFile =1:length(emptyYN)
        %        TSsig_ = cellfun(@(x) sum(x(202:end,:),1).*bw,TSsig_,'UniformOutput',false);
        FracSigErr{iArea}(iFile,iDelay) = sum(sum( TSsig_Err{iFile},1).*bw > Sigthresh)./size( TSsig_Err{iFile},2)
    
        timehist = cumsum(histc(TSsig_Err{iFile},bins));
        timeHist(:,iFile) = timehist./max(timehist);
    end
    
    timeHist_Mean = nanmean(timeHist,2);
    timeHist_SEM = nansem(timeHist,2);
%     ciplot(timeHist_Mean+timeHist_SEM,timeHist_Mean-timeHist_SEM,...
%            bins,color_{iDelay},alpha_Err)
    plot(bins,timeHist_Mean,'Color',color_{iDelay},'LineWidth',1.5,'LineStyle',':')
end
   title(Areas{iArea})
    end
    axis([0 15 0. 1])
    set(gca,'Ytick',0:0.5:1,'Xtick',0:5:maxBin)
    xlabel('Time (s)')
    ylabel('Fraction of Assemblies')
    title(Areas{iArea})
end

%% Plot L/R decoding by area - RT2
figure('name','L/R decoding');
for iArea =1:length(Areas)
    subplot(1,3,iArea);hold on
    
    ciplot(nanmean(D_Ass_.Short.LR.Rt2{iArea},2)+nansem(D_Ass_.Short.LR.Rt2{iArea},2),...
           nanmean(D_Ass_.Short.LR.Rt2{iArea},2)-nansem(D_Ass_.Short.LR.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'b',.9)
    ciplot(nanmean(D_Ass_.Medium.LR.Rt2{iArea},2)+nansem(D_Ass_.Medium.LR.Rt2{iArea},2),...
           nanmean(D_Ass_.Medium.LR.Rt2{iArea},2)-nansem(D_Ass_.Medium.LR.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'g',0.9)    
    ciplot(nanmean(D_Ass_.Long.LR.Rt2{iArea},2)+nansem(D_Ass_.Long.LR.Rt2{iArea},2),...
           nanmean(D_Ass_.Long.LR.Rt2{iArea},2)-nansem(D_Ass_.Long.LR.Rt2{iArea},2),...
           [0:(length(tbAll)*2-1)]*bw,'r',0.9)               
       
%     plot(tbShort,mean(D_.Short.LR.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.LR.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.LR.Ft2{iArea},2),'b')
    plot([10 10],[0.5 2],'color',[0 0 0 0.3],'LineWidth',1.5)
    plot([5 5],[0.1 0.5],'color',[0 1 0 0.3],'LineWidth',1.5)
    plot([15 15],[0.1 0.5],'color',[1 0 0 0.3],'LineWidth',1.5)
    if normaliseFscores
        axis([0 20 0 20])
    else
        axis([0 20 0 1])
    end

%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'},'Location','northoutside','Orientation','horizontal')
    legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
     if iArea==1
        ylabel('L/R decoding (F-Score)')
    end
end
%% Plot L/R decoding by area - FT2
figure('color','w','name','L/R decoding');
alpha_ = 0.5;
sigcutoff = 0.5;
x_ = [0:(length(tbAll)*2-1)]*bw;
for iArea =1:length(Areas)
    subplot(1,3,iArea);hold on
    ciplot(nanmean(D_Ass_.Short.LR.Ft2{iArea},2)+nansem(D_Ass_.Short.LR.Ft2{iArea},2),...
           nanmean(D_Ass_.Short.LR.Ft2{iArea},2)-nansem(D_Ass_.Short.LR.Ft2{iArea},2),...
           x_,'b',alpha_)
    ciplot(nanmean(D_Ass_.Medium.LR.Ft2{iArea},2)+nansem(D_Ass_.Medium.LR.Ft2{iArea},2),...
           nanmean(D_Ass_.Medium.LR.Ft2{iArea},2)-nansem(D_Ass_.Medium.LR.Ft2{iArea},2),...
           x_,'g',alpha_)    
    ciplot(nanmean(D_Ass_.Long.LR.Ft2{iArea},2)+nansem(D_Ass_.Long.LR.Ft2{iArea},2),...
           nanmean(D_Ass_.Long.LR.Ft2{iArea},2)-nansem(D_Ass_.Long.LR.Ft2{iArea},2),...
           x_,'r',alpha_)           
       
       
       
    temp = nansum(cell2mat(D_Ass_.Short.LR.Ft2sig{iArea}),2)./sum(~isnan(sum(cell2mat(D_Ass_.Short.LR.Ft2sig{iArea}),1)));
%     temp=temp>=sigcutoff;
    imagesc(repmat(x_,1,2),[-0*ones(size(temp,1),1);-0.2*ones(size(temp,1),1)], [temp,temp]')
    rectangle('Position',[0 -0.3 20 0.4],'LineWidth',1,'EdgeColor','b')

    temp =nansum(cell2mat(D_Ass_.Medium.LR.Ft2sig{iArea}),2)./sum(~isnan(sum(cell2mat(D_Ass_.Medium.LR.Ft2sig{iArea}),1)));
%     temp=temp>=sigcutoff;
    imagesc(repmat(x_,1,2),[-0.5*ones(size(temp,1),1);-0.7*ones(size(temp,1),1)],[temp,temp]')
    rectangle('Position',[0 -0.8 20 0.4],'LineWidth',1,'EdgeColor','g')

    temp =nansum(cell2mat(D_Ass_.Long.LR.Ft2sig{iArea}),2)./sum(~isnan(sum(cell2mat(D_Ass_.Long.LR.Ft2sig{iArea}),1)));
%     temp=temp>=sigcutoff;
    imagesc(repmat(x_,1,2),[-1*ones(size(temp,1),1);-1.2*ones(size(temp,1),1)],[temp,temp]')
    rectangle('Position',[0 -1.3 20 0.4],'LineWidth',1,'EdgeColor','r')       

       
%     plot(tbShort,mean(D_.Short.LR.Ft2{iArea},2),'r')
%     plot(tbMedium,mean(D_.Medium.LR.Ft2{iArea},2),'g')
%     plot(tbLong,mean(D_.Long.LR.Ft2{iArea},2),'b')
    plot([10 10],[0.2 0.6],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[0.2 0.6],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[0.2 0.6],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    plot([10 10],[-1.8 -1.4],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[-1.8 -1.4],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[-1.8 -1.4],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    axis([0 20 -2 50])
%     axis([0 20 0 1])
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'},'Location','northoutside','Orientation','horizontal')
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
     if iArea==1
        ylabel('L/R decoding (F-Score)')
        plot([2 7], [5 5],'k','LineWidth',1.5)
        plot([2 2], [5 8],'k','LineWidth',1.5)
    end
    colormap(flipud(gray(10)));
    caxis([0 1])
    axis off
end

%% %%%%%%%%%%%%%%%%%%% Coding: UNITS vs. ASSEMBLIES %%%%%%%%%%%%%%%%%%%%%%%%%
%% Distributions of peak t-Scores - shaded
plotErr = false;
bins = 0:0.1:10;  
alpha_ = 0.3;
alpha_Err = 0.3;
figure
for iArea = 1:length(Areas)
    subplot(1,3,iArea); hold on
     for iDelay = 1:length(Delays_)
     if iArea<3
       TS_ = eval(sprintf('D_.%s.LR.TS{iArea}',Delays_{iDelay}));
       TS_ = cellfun(@(x) max(x,[],1),TS_,'UniformOutput',false);
       for iFile =1:length(TS_)
           maxhist = cumsum(histc(TS_{iFile},bins));
           maxHist(:,iFile) = maxhist./max(maxhist);
       end
       maxHist_Mean = nanmean(maxHist,2);
       maxHist_SEM = nansem(maxHist,2);
       ciplot(maxHist_Mean+maxHist_SEM,maxHist_Mean-maxHist_SEM,...
              bins,color_{iDelay},alpha_)
     end
       TS_Ass = eval(sprintf('D_Ass_.%s.LR.TS{iArea}',Delays_{iDelay}));
       TS_Ass = cellfun(@(x) max(x,[],1),TS_Ass,'UniformOutput',false);
       for iFile =1:length(TS_Ass)
           maxhist = cumsum(histc(TS_Ass{iFile},bins));
           maxHist(:,iFile) = maxhist./max(maxhist);
       end
       maxHist_Mean = nanmean(maxHist,2);
       maxHist_SEM = nansem(maxHist,2);
       ciplot(maxHist_Mean+maxHist_SEM,maxHist_Mean-maxHist_SEM,...
              bins,color_{iDelay},1)          
     end
 
        set(gca,'Ytick',0:0.5:1,'Xtick',0:5:20)
        xlabel('Peak t-score')
        ylabel('Fraction of units')
        title(Areas{iArea})
    
end
%% Distributions of time significant - shaded
plotErr = false;
figure
clear timehist timeHist FracSig FracSigErr
Sigthresh = bw;
maxBin = 15;
bins   = Sigthresh:bw:maxBin;
alpha_ = 0.3;
alpha_Err = 0.2;

for iArea = 1:length(Areas)
    subplot(1,3,iArea); hold on
    for iDelay  = 1:length(Delays_)
    if iArea<3
        TSsig_   = eval(sprintf('D_.%s.LR.TSsig{iArea}',Delays_{iDelay}));
        
        %        TSsig_ = cellfun(@(x) sum(x(202:end,:),1).*bw,TSsig_,'UniformOutput',false);
        TSsig_ = cellfun(@(x) sum(x,1).*bw,TSsig_,'UniformOutput',false);
        FracSig{iArea}(:,iDelay) = cell2mat(cellfun(@(x) sum(sum(x,1).*bw > Sigthresh)./size(x,2),TSsig_,'UniformOutput',false));
        
        for iFile =1:length(TSsig_)
            timehist = cumsum(histc(TSsig_{iFile},bins));
            timeHist(:,iFile) = timehist./max(timehist);
        end
        
        timeHist_Mean = nanmean(timeHist,2);
        timeHist_SEM = nansem(timeHist,2);
        ciplot(timeHist_Mean+timeHist_SEM,timeHist_Mean-timeHist_SEM,...
            bins,color_{iDelay},alpha_)
        %        plot(bins,timeHist_Mean,'Color',color_{iDelay},'LineWidth',1.5)
    end
    TSsig_Ass_   = eval(sprintf('D_Ass_.%s.LR.TSsig{iArea}',Delays_{iDelay}));
    
    %        TSsig_Ass_ = cellfun(@(x) sum(x(202:end,:),1).*bw,TSsig_Ass_,'UniformOutput',false);
    TSsig_Ass_ = cellfun(@(x) sum(x,1).*bw,TSsig_Ass_,'UniformOutput',false);
    FracSig{iArea}(:,iDelay) = cell2mat(cellfun(@(x) sum(sum(x,1).*bw > Sigthresh)./size(x,2),TSsig_Ass_,'UniformOutput',false));
    
    for iFile =1:length(TSsig_Ass_)
        timehist = cumsum(histc(TSsig_Ass_{iFile},bins));
        timeHist(:,iFile) = timehist./max(timehist);
    end
    
    timeHist_Mean = nanmean(timeHist,2);
    timeHist_SEM = nansem(timeHist,2);
    ciplot(timeHist_Mean+timeHist_SEM,timeHist_Mean-timeHist_SEM,...
        bins,color_{iDelay},1)
    %        plot(bins,timeHist_Mean,'Color',color_{iDelay},'LineWidth',1.5)
    end
    title(Areas{iArea})
    
    axis([0 15 0. 1])
    set(gca,'Ytick',0:0.5:1,'Xtick',0:5:maxBin)
    xlabel('Time (s)')
    ylabel('Fraction of units')
    title(Areas{iArea})
end

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

%     nu(3) = sum(nu(1:2))   
%     usel_out{3} = [usel_out(1),usel_out(1) 
    Ass.usel_out{iFile} = usel_out;
    Ass.nu{iFile} = nu;
    Ass.units{iFile} = units;
%     for iArea =1:2
%          Ass.units{iFile}{iArea} = usel_out{iArea}(units{iArea})
%     end
%     clear units usel_out nu
end
%% Pick example assemblies, show units and assembly decoding side-by-side
iArea =3;
for iFile =1:length(fileList)
iDelay = 3;
for iAss = 1:length(Ass.units{iFile}{iArea})
    try
    eval(sprintf('TS_D_Ass_ = D_Ass_.%s.LR.TS{iArea}{iFile}(:,iAss);',Delays_{iDelay}));
    idx = Ass.units{iFile}{iArea}{iAss};
    if iArea<3
        eval(sprintf('TS_D_Units_ = D_.%s.LR.TS{1}{iFile};',Delays_{iDelay}));
    else
        eval(sprintf('TS_D_Units_ = [D_.%s.LR.TS{1}{iFile},D_.%s.LR.TS{2}{iFile}];',Delays_{iDelay},Delays_{iDelay}));
    end
    TS_D_Units_ = TS_D_Units_(:,idx);
    tb = (1:size(TS_D_Ass_,1))*bw;
    figure; hold on
    plot(tb,TS_D_Ass_,'k','LineWidth',1.5);
    plot(tb,TS_D_Units_,'r','LineWidth',1.5);
    end
end
end
%% ANOVA Co-tuning of Assemblies and Single Unit members
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







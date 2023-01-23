% %%%%%% PREAMBLE %%%%%%
clear

Delays_ = {'Short','Medium','Long'};
Target = 'LONG';
Areas = {'PFC','HP','Joint'};

bw=0.05; tlimsAll = [-5 5]; tbAll = tlimsAll(1):bw:tlimsAll(2);

maxRange = 20;    % Upper limit of assembly size to explore
noTrials = 5;   % How many trials to consider simultanrously
nReps = 50;     % how many repeated draws among trials to perform


InfoCriteria = 'max'; % 'mean','max'

if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
% elseif isunix
    %     pat ='/Volumes/Data/DNMTP/raw';

elseif ismac
    pat = '/Volumes/Aleks Data Portable/AssemblyAnalysis/raw/';
end

fileList = dir([pat 'SyntheticAssemblies' filesep '*' Target '*_Signalcorr_redux8stats.mat']);
fileList = dir([pat 'SyntheticAssemblies' filesep '*' Target '*_Signalcorr_redux8stats.mat']);

%% Batch import
for iFile =1:length(fileList)
    
    %fnIn = fullfile(pat, 'SyntheticAssemblies', fileList(iFile).name);
    fnIn = fullfile(pat, 'SyntheticAssemblies',fileList(iFile).name);
    D{iFile} = load(fnIn,'D_');
    
    %     fnIn = [strtok(fnIn,'.'),'stats.mat'];
    %     temp = load(fnIn,'D_');
    %     D{iFile}.D_.AssReal.Ft2 = temp.D_.AssReal.CVE;
end
clear temp

%% Summary stats
for s=1:3
    for iFile =1:length(D)
        try
            D_.AssReal.unitDistMean{s}{iFile,1} = D{iFile}.D_.AssReal.unitDistMean{s};
            D_.AssReal.peakCVE{s}{iFile,1} = nanmax(D{iFile}.D_.AssReal.CVE{s},1);
            D_.AssReal.peakCVE_FSC{s}{iFile,1} = nanmax(D{iFile}.D_.AssReal.CVE_FSC{s},1);
        end
    end
    idx = isempty_cell(D_.AssReal.unitDistMean{s});
    D_.AssReal.unitDistMeanCollapsed{s} = cell2mat(D_.AssReal.unitDistMean{s}(~idx));
    D_.AssReal.peakCVECollapsed{s} = cell2mat(D_.AssReal.peakCVE{s}(~idx)')';
    D_.AssReal.peakCVE_FSCCollapsed{s} = cell2mat(D_.AssReal.peakCVE_FSC{s}(~idx)')';
    
end

%% no Units vs peak multivariate coding
clear mdl pFit Rsq slope
figure('color','w','name','no Units vs peak multivariate coding');
for s=1:3
    
    subplot(1,3,s);    hold on
    noUnits    = D_.AssReal.unitDistMeanCollapsed{s}(:,1);
    codingDist = D_.AssReal.unitDistMeanCollapsed{s}(:,2);
    peakCoding = D_.AssReal.peakCVECollapsed{s};
    scatter(noUnits,peakCoding)
    
    mdl{s} = fitlm(noUnits,peakCoding,'RobustOpts','on');
    %             mdl{s} = fitlm(codingDist,peakCoding);
    pFit(s)   = mdl{s}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s)   =  mdl{s}.Rsquared.Adjusted;
    slope(s) = mdl{s}.Coefficients{2,1};
    
    h = plot(mdl{s});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'r';
    h(3).Color = 'r';
    h(4).Color = 'r';
    %             axis([2 20 0.4 1])
    legend off
end
%% no Units vs peak Factor score coding
clear mdl pFit Rsq slope
figure('color','w','name','no Units vs peak multivariate coding');
for s=1:3
    
    subplot(1,3,s);    hold on
    noUnits    = D_.AssReal.unitDistMeanCollapsed{s}(:,1);
    codingDist = D_.AssReal.unitDistMeanCollapsed{s}(:,2);
    peakCoding = D_.AssReal.peakCVE_FSCCollapsed{s};
    scatter(noUnits,peakCoding)
    
    mdl{s} = fitlm(noUnits,peakCoding,'RobustOpts','on');
    %             mdl{s} = fitlm(codingDist,peakCoding);
    pFit(s)   = mdl{s}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s)   =  mdl{s}.Rsquared.Adjusted;
    slope(s) = mdl{s}.Coefficients{2,1};
    
    h = plot(mdl{s});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'r';
    h(3).Color = 'r';
    h(4).Color = 'r';
    %             axis([2 20 0.4 1])
    legend off
end
%% coding distances vs. peak multivariate/Factor score coding
clear mdl pFit Rsq slope
figure('color','w','name','coding distances vs. peak multivariate/Factor score coding');
for s=1:3
    
    subplot(1,3,s);    hold on
    noUnits    = D_.AssReal.unitDistMeanCollapsed{s}(:,1);
    codingDist = D_.AssReal.unitDistMeanCollapsed{s}(:,2);
    peakCodingUnits = D_.AssReal.peakCVECollapsed{s};
    peakCodingFSC = D_.AssReal.peakCVE_FSCCollapsed{s};
    scatter(codingDist,peakCodingUnits,'b')
    
    mdl{s} = fitlm(codingDist,peakCodingUnits,'RobustOpts','on');
    PFitUnits(s)   = mdl{s}.Coefficients{2,4};
    Rsq(s)         = mdl{s}.Rsquared.Adjusted;
    slope(s)       = mdl{s}.Coefficients{2,1};
    
    h = plot(mdl{s});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'b';
    h(3).Color = 'b';
    h(4).Color = 'b';
    
    scatter(codingDist,peakCodingFSC,'r')
    mdl{s} = fitlm(codingDist,peakCodingFSC,'RobustOpts','on');
    PFitFSC(s)      = mdl{s}.Coefficients{2,4};
    Rsq(s)          = mdl{s}.Rsquared.Adjusted;
    slope(s)        = mdl{s}.Coefficients{2,1};
    
    h = plot(mdl{s});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'r';
    h(3).Color = 'r';
    h(4).Color = 'r';
    
    %             axis([0 600 0 600])
%     set(gca,'yscale','log')
    
    legend off
end
%% no Units vs. peak multivariate/Factor score coding
clear mdl pFit Rsq slope
figure('color','w','name','no Units vs. peak multivariate/Factor score coding');
for s=1:3
    
    subplot(1,3,s);    hold on
    noUnits    = D_.AssReal.unitDistMeanCollapsed{s}(:,1);
    codingDist = D_.AssReal.unitDistMeanCollapsed{s}(:,2);
    peakCodingUnits = D_.AssReal.peakCVECollapsed{s};
    peakCodingFSC = D_.AssReal.peakCVE_FSCCollapsed{s};
    
    scatter(noUnits,peakCodingUnits,'b')
    mdl{s} = fitlm(noUnits,peakCodingUnits,'RobustOpts','on');
    PFitUnits(s)   = mdl{s}.Coefficients{2,4};
    Rsq(s)   =  mdl{s}.Rsquared.Adjusted;
    slope(s) = mdl{s}.Coefficients{2,1};
    
    h = plot(mdl{s});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'b';
    h(3).Color = 'b';
    h(4).Color = 'b';
    
    scatter(noUnits,peakCodingFSC,'r')
    mdl{s} = fitlm(noUnits,peakCodingFSC,'RobustOpts','on');
    PFitFSC(s)   = mdl{s}.Coefficients{2,4};
    Rsq(s)   =  mdl{s}.Rsquared.Adjusted;
    slope(s) = mdl{s}.Coefficients{2,1};
    
    h = plot(mdl{s},'HandleVisibility','off');
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'r';
    h(3).Color = 'r';
    h(4).Color = 'r';
    
    %             axis([2 20 0.4 1])
    %             set(gca,'yscale','log')
    
    legend off
    
end
%% peak multivariate vs peak Factor score coding
clear mdl pFit Rsq slope
figure('color','w','name','peak multivariate vs peak Factor score coding');
for s=1:3
    
    subplot(1,3,s);    hold on
    peakCodingunits = D_.AssReal.peakCVECollapsed{s};
    peakCodingFSC = D_.AssReal.peakCVE_FSCCollapsed{s};
    scatter(peakCodingunits,peakCodingFSC)
    
    mdl{s} = fitlm(peakCodingunits,peakCodingFSC,'RobustOpts','on');
    %             mdl{s} = fitlm(codingDist,peakCoding);
    pFit(s)   = mdl{s}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s)   =  mdl{s}.Rsquared.Adjusted;
    slope(s) = mdl{s}.Coefficients{2,1};
    
    h = plot(mdl{s});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'r';
    h(3).Color = 'r';
    h(4).Color = 'r';
    %             axis([2 20 0.4 1])
    legend off
    set(gca,'yscale','log','xscale','log')
    
end
%% difference inmultivariate and Factor score coding
clear mdl pFit Rsq slope
figure('color','w','name','peak multivariate vs peak Factor score coding');
hold on
% bins =-1:0.1:1;
bins =-1:0.02:1;
col_ = {'b','r','g'};
for s=1:3
    subplot(1,3,s);    hold on

    peakCodingunits = D_.AssReal.peakCVECollapsed{s};
    peakCodingFSC   = D_.AssReal.peakCVE_FSCCollapsed{s};
    codingDist = D_.AssReal.unitDistMeanCollapsed{s}(:,2);
    noUnits         = D_.AssReal.unitDistMeanCollapsed{s}(:,1);
    delta_ = (peakCodingFSC-peakCodingunits)./(peakCodingFSC+peakCodingunits);
%     delta_ = peakCodingFSC-peakCodingunits;
    temp   = histc(delta_,bins);
    temp   = temp./sum(temp);
%     stairs(bins,temp,'LineWidth',1.5,'Color',col_{s})
%     plot(bins,cumsum(temp),'LineWidth',1.5,'Color',col_{s})
    scatter(codingDist,delta_,col_{s})
%     axis([0 20 0 100])
%     axis([0 20 -1 1])

    mdl{s} = fitlm(codingDist,delta_,'RobustOpts','on');
    %             mdl{s} = fitlm(codingDist,peakCoding);
    pFit(s)   = mdl{s}.Coefficients{2,4};
     h = plot(mdl{s});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'r';
    h(3).Color = 'r';
    h(4).Color = 'r';
    %             axis([2 20 0.4 1])
    legend off
%     set(gca,'yscale','log','xscale','log')
    
end

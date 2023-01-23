%% %%%%%% PREAMBLE %%%%%%

clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};

Target = 'LONG';
RunErrors = false;
tlimsAll = [0 5];
tlimsShort=[0 4];
tlimsMedium=[0 8];
tlimsLong=[0 16];

shift = 0;
plotOnline = false;
bw    = 0.05;
Nbs = 500;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);

fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
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
groupSize = 5;
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
mkdir([pat 'MixedSelectivity'])

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
%% batch import unit data for meta-analysis and plotting
UnitSelection = 'all'; %{'pairs','groups','all','groups_TrialDraw'}


clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
    pat = '/Volumes/Aleks Data Portable/AssemblyAnalysis/raw/'
end

fileList = dir([pat 'MixedSelectivity' filesep 'DelayDecodingVsBehav' filesep '*' Target '*' 'PFC' '*' UnitSelection '*.mat']);

clear D_
for iArea = 1:3
    for iFile = 1:length(fileList)
        try
            
            fname=strtok(fileList(iFile).name,'_');
            fnIn = sprintf('%sMixedSelectivity%s%s%s%s_%s_DelayDecoding_UnitSpan_%s.mat',pat,filesep,'DelayDecodingVsBehav',filesep,fname,Areas{iArea},UnitSelection);
            D_{iArea}{iFile} = load(fnIn);
            fprintf('Done.\n')
        end
    end
end

%% Collapse - Short delays
clear D_Collapsed
iArea=1

for iFile = 1:length(fileList)
    for iDelay =1:3
        D_Collapsed.meanDecoding.Short(iFile,:) = D_{iArea}{iFile}.D.Short.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Medium(iFile,:) = D_{iArea}{iFile}.D.Medium.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Long(iFile,:) = D_{iArea}{iFile}.D.Long.DelaySpan.CVE;
        
        D_Collapsed.meanDecoding.ShortCIL(iFile,:)  = cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledL_raw(:)');
        D_Collapsed.meanDecoding.MediumCIL(iFile,:) = cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledL_raw(:)');
        D_Collapsed.meanDecoding.LongCIL(iFile,:)   = cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledL_raw(:)');
        
        D_Collapsed.meanDecoding.ShortCIH(iFile,:)  = cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledH_raw(:)');
        D_Collapsed.meanDecoding.MediumCIH(iFile,:) = cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledH_raw(:)');
        D_Collapsed.meanDecoding.LongCIH(iFile,:)   = cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledH_raw(:)');
        
        D_Collapsed.meanDecoding.tbShort = D_{iArea}{iFile}.tbShort;
        D_Collapsed.meanDecoding.tbMedium = D_{iArea}{iFile}.tbMedium;
        D_Collapsed.meanDecoding.tbLong = D_{iArea}{iFile}.tbLong;
    end
    
    D_Collapsed.SigDecodeTimeDelayMeans(iFile) = nanmean([D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime ./diff(tlimsShort), ...
        D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeTime ./diff(tlimsMedium), ...
        D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeTime   ./diff(tlimsLong)]);
    D_Collapsed.SigDecodeFracDrawsDelayMeans(iFile) = nanmean([sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsShort) , ...
        sum(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)./diff(tlimsMedium), ...
        sum(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong)]);
    
    D_Collapsed.AccuracyDelayMeans(iFile) = mean(D_{1}{iFile}.D.Accuracy);
    D_Collapsed.pCorrDelayMeans(iFile) = mean(D_{1}{iFile}.D.pCorr);
    
    D_Collapsed.SigDecodeTime(iFile)        = D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime./diff(tlimsLong);
    D_Collapsed.SigDecodeFracDraws(iFile)   = sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong);
    D_Collapsed.LastSigTime(iFile)          = max([0,find(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)])*bw;
    
    D_Collapsed.Accuracy(iFile) = D_{1}{iFile}.D.Accuracy(1);
    D_Collapsed.pCorr(iFile)    = D_{1}{iFile}.D.pCorr(1);
    
    D_Collapsed.AboveChance(iFile,:) = D_{1}{iFile}.D.AboveChance;
end
%% Plot performance vs decoding - Short delays
figure
subplot(1,2,1); hold on

x = D_Collapsed.LastSigTime;
y = D_Collapsed.pCorr;
% y = D_Collapsed.Accuracy;
scatter(x(~D_Collapsed.AboveChance(:,1)),y(~D_Collapsed.AboveChance(:,1)),'b')
scatter(x(D_Collapsed.AboveChance(:,1)),y(D_Collapsed.AboveChance(:,1)),'b','filled')
[R,P]= corrcoef([x;y]')

mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
%     %             pFit(iArea)   = coefTest(mdl{iArea});
Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
slope(iArea) = mdl{iArea}.Coefficients{2,1};
intercept(iArea) = mdl{iArea}.Coefficients{1,1};

plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
[pFit(iArea), Rsq(iArea)]

h = plot(mdl{iArea});
h(1).Marker = 'none';
h(1).Marker = 'none';
h(2).LineWidth=1.5;
h(2).LineStyle='-';
h(2).Color = 'k';
h(3).Color = 'k';
h(4).Color = 'k';
axis([0 5 0.4 1])
legend off
xlabel({'Last significant timepoint (s)'})
ylabel('Fraction of choices correct')
title(Areas{iArea})

subplot(1,2,2); hold on
x = D_Collapsed.SigDecodeTime*4;
% x = D_Collapsed.SigDecodeFracDraws;
y = D_Collapsed.pCorr;
% y = D_Collapsed.Accuracy;
scatter(x(~D_Collapsed.AboveChance(:,2)),y(~D_Collapsed.AboveChance(:,2)),'b')
scatter(x(D_Collapsed.AboveChance(:,2)),y(D_Collapsed.AboveChance(:,2)),'b','filled')
[R,P]= corrcoef([x;y]');

mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
%     %             pFit(iArea)   = coefTest(mdl{iArea});
Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
slope(iArea) = mdl{iArea}.Coefficients{2,1};
intercept(iArea) = mdl{iArea}.Coefficients{1,1};

plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
[pFit(iArea), Rsq(iArea)]

h = plot(mdl{iArea});
h(1).Marker = 'none';
h(1).Marker = 'none';
h(2).LineWidth=1.5;
h(2).LineStyle='-';
h(2).Color = 'k';
h(3).Color = 'k';
h(4).Color = 'k';
axis([0 1 0.4 1])
legend off
xlabel({'Significant cue decoding span (s)'})
ylabel('Fraction of choices correct')
title(Areas{iArea})

%% Collapse - Medium delays
clear D_Collapsed
iArea=1

for iFile = 1:length(fileList)
    for iDelay =1:3
        D_Collapsed.meanDecoding.Short(iFile,:) = D_{iArea}{iFile}.D.Short.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Medium(iFile,:) = D_{iArea}{iFile}.D.Medium.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Long(iFile,:) = D_{iArea}{iFile}.D.Long.DelaySpan.CVE;
        
        D_Collapsed.meanDecoding.ShortCIL(iFile,:)  = cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledL_raw(:)');
        D_Collapsed.meanDecoding.MediumCIL(iFile,:) = cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledL_raw(:)');
        D_Collapsed.meanDecoding.LongCIL(iFile,:)   = cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledL_raw(:)');
        
        D_Collapsed.meanDecoding.ShortCIH(iFile,:)  = cell2mat(D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledH_raw(:)');
        D_Collapsed.meanDecoding.MediumCIH(iFile,:) = cell2mat(D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledH_raw(:)');
        D_Collapsed.meanDecoding.LongCIH(iFile,:)   = cell2mat(D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledH_raw(:)');
        
        D_Collapsed.meanDecoding.tbShort = D_{iArea}{iFile}.tbShort;
        D_Collapsed.meanDecoding.tbMedium = D_{iArea}{iFile}.tbMedium;
        D_Collapsed.meanDecoding.tbLong = D_{iArea}{iFile}.tbLong;
    end
    
    D_Collapsed.SigDecodeTimeDelayMeans(iFile) = nanmean([D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime ./diff(tlimsShort), ...
        D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeTime ./diff(tlimsMedium), ...
        D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeTime   ./diff(tlimsLong)]);
    D_Collapsed.SigDecodeFracDrawsDelayMeans(iFile) = nanmean([sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsShort) , ...
        sum(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)./diff(tlimsMedium), ...
        sum(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong)]);
    
    D_Collapsed.AccuracyDelayMeans(iFile) = mean(D_{1}{iFile}.D.Accuracy);
    D_Collapsed.pCorrDelayMeans(iFile) = mean(D_{1}{iFile}.D.pCorr);
    
    D_Collapsed.SigDecodeTime(iFile)        = D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeTime./diff(tlimsLong);
    D_Collapsed.SigDecodeFracDraws(iFile)   = sum(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong);
    D_Collapsed.LastSigTime(iFile)          = max([0,find(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws>0)])*bw;
    
    D_Collapsed.Accuracy(iFile) = D_{1}{iFile}.D.Accuracy(2);
    D_Collapsed.pCorr(iFile)    = D_{1}{iFile}.D.pCorr(2);
    
    D_Collapsed.AboveChance(iFile,:) = D_{1}{iFile}.D.AboveChance;
end
%% Plot performance vs decoding - Medium delays
figure
subplot(1,2,1); hold on

x = D_Collapsed.LastSigTime;
y = D_Collapsed.pCorr;
% y = D_Collapsed.Accuracy;
scatter(x(~D_Collapsed.AboveChance(:,2)),y(~D_Collapsed.AboveChance(:,2)),'b')
scatter(x(D_Collapsed.AboveChance(:,2)),y(D_Collapsed.AboveChance(:,2)),'b','filled')
[R,P]= corrcoef([x;y]')

mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
%     %             pFit(iArea)   = coefTest(mdl{iArea});
Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
slope(iArea) = mdl{iArea}.Coefficients{2,1};
intercept(iArea) = mdl{iArea}.Coefficients{1,1};

plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
[pFit(iArea), Rsq(iArea)]

h = plot(mdl{iArea});
h(1).Marker = 'none';
h(1).Marker = 'none';
h(2).LineWidth=1.5;
h(2).LineStyle='-';
h(2).Color = 'k';
h(3).Color = 'k';
h(4).Color = 'k';
axis([0 16 0.4 1])
legend off
xlabel({'Last significant timepoint (s)'})
ylabel('Fraction of choices correct')
title(Areas{iArea})

subplot(1,2,2); hold on
x = D_Collapsed.SigDecodeTime*8;
% x = D_Collapsed.SigDecodeFracDraws;
y = D_Collapsed.pCorr;
% y = D_Collapsed.Accuracy;
scatter(x(~D_Collapsed.AboveChance(:,2)),y(~D_Collapsed.AboveChance(:,2)),'b')
scatter(x(D_Collapsed.AboveChance(:,2)),y(D_Collapsed.AboveChance(:,2)),'b','filled')
[R,P]= corrcoef([x;y]');

mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
%     %             pFit(iArea)   = coefTest(mdl{iArea});
Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
slope(iArea) = mdl{iArea}.Coefficients{2,1};
intercept(iArea) = mdl{iArea}.Coefficients{1,1};

plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
[pFit(iArea), Rsq(iArea)]

h = plot(mdl{iArea});
h(1).Marker = 'none';
h(1).Marker = 'none';
h(2).LineWidth=1.5;
h(2).LineStyle='-';
h(2).Color = 'k';
h(3).Color = 'k';
h(4).Color = 'k';
axis([0 6 0.4 1])
legend off
xlabel({'Significant cue decoding span (s)'})
ylabel('Fraction of choices correct')
title(Areas{iArea})

%% Collapse - Long delays
clear D_Collapsed
iArea=1

for iFile = 1:length(fileList)
    for iDelay =1:3
        D_Collapsed.meanDecoding.Short(iFile,:) = D_{iArea}{iFile}.D.Short.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Medium(iFile,:) = D_{iArea}{iFile}.D.Medium.DelaySpan.CVE;
        D_Collapsed.meanDecoding.Long(iFile,:) = D_{iArea}{iFile}.D.Long.DelaySpan.CVE;
        
        D_Collapsed.meanDecoding.ShortCIL(iFile,:)  = nanmean(cell2mat(cellfun(@transpose,D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledL_raw(:),'UniformOutput',false)'),2);
        D_Collapsed.meanDecoding.MediumCIL(iFile,:) = nanmean(cell2mat(cellfun(@transpose,D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledL_raw(:),'UniformOutput',false)'),2);
        D_Collapsed.meanDecoding.LongCIL(iFile,:)   = nanmean(cell2mat(cellfun(@transpose,D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledL_raw(:),'UniformOutput',false)'),2);
        
        D_Collapsed.meanDecoding.ShortCIH(iFile,:)  = nanmean(cell2mat(cellfun(@transpose,D_{iArea}{iFile}.D.Short.DelaySpan.CVE_shuffledH_raw(:),'UniformOutput',false)'),2);
        D_Collapsed.meanDecoding.MediumCIH(iFile,:) = nanmean(cell2mat(cellfun(@transpose,D_{iArea}{iFile}.D.Medium.DelaySpan.CVE_shuffledH_raw(:),'UniformOutput',false)'),2);
        D_Collapsed.meanDecoding.LongCIH(iFile,:)   = nanmean(cell2mat(cellfun(@transpose,D_{iArea}{iFile}.D.Long.DelaySpan.CVE_shuffledH_raw(:),'UniformOutput',false)'),2);
        
        D_Collapsed.meanDecoding.tbShort = D_{iArea}{iFile}.tbShort;
        D_Collapsed.meanDecoding.tbMedium = D_{iArea}{iFile}.tbMedium;
        D_Collapsed.meanDecoding.tbLong = D_{iArea}{iFile}.tbLong;
    end
    
    D_Collapsed.SigDecodeTimeDelayMeans(iFile) = nanmean([D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime ./diff(tlimsShort), ...
        D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeTime ./diff(tlimsMedium), ...
        D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeTime   ./diff(tlimsLong)]);
    D_Collapsed.SigDecodeFracDrawsDelayMeans(iFile) = nanmean([sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsShort) , ...
        sum(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)./diff(tlimsMedium), ...
        sum(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong)]);
    
    D_Collapsed.AccuracyDelayMeans(iFile) = mean(D_{1}{iFile}.D.Accuracy);
    D_Collapsed.pCorrDelayMeans(iFile) = mean(D_{1}{iFile}.D.pCorr);
    
    D_Collapsed.SigDecodeTime(iFile)        = D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeTime./diff(tlimsLong);
    D_Collapsed.SigDecodeFracDraws(iFile)   = sum(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong);
    D_Collapsed.LastSigTime(iFile)          = max([0,find(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws>0.5)])*bw;
    
    D_Collapsed.Accuracy(iFile) = D_{1}{iFile}.D.Accuracy(3);
    D_Collapsed.pCorr(iFile)    = D_{1}{iFile}.D.pCorr(3);
    
    D_Collapsed.AboveChance(iFile,:) = D_{1}{iFile}.D.AboveChance;
end
%% Plot performance vs decoding - Long delays
figure
subplot(1,2,1); hold on

x = D_Collapsed.LastSigTime;
y = D_Collapsed.pCorr;
% y = D_Collapsed.Accuracy;
scatter(x(~D_Collapsed.AboveChance(:,3)),y(~D_Collapsed.AboveChance(:,3)),color_{iArea})
scatter(x(D_Collapsed.AboveChance(:,3)),y(D_Collapsed.AboveChance(:,3)),color_{iArea},'filled')
[R,P]= corrcoef([x;y]');

mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
%     %             pFit(iArea)   = coefTest(mdl{iArea});
Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
slope(iArea) = mdl{iArea}.Coefficients{2,1};
intercept(iArea) = mdl{iArea}.Coefficients{1,1};

plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
[pFit(iArea), Rsq(iArea)]

h = plot(mdl{iArea});
h(1).Marker = 'none';
h(1).Marker = 'none';
h(2).LineWidth=1.5;
h(2).LineStyle='-';
h(2).Color = 'k';
h(3).Color = 'k';
h(4).Color = 'k';
axis([0 17 0.4 1])
legend off
xlabel({'Last significant timepoint (s)'})
ylabel('Fraction of choices correct')
title(Areas{iArea})

subplot(1,2,2); hold on
x = D_Collapsed.SigDecodeTime*16;
% x = D_Collapsed.SigDecodeFracDraws;
y = D_Collapsed.pCorr;
% y = D_Collapsed.Accuracy;
scatter(x(~D_Collapsed.AboveChance(:,3)),y(~D_Collapsed.AboveChance(:,3)),color_{iArea})
scatter(x(D_Collapsed.AboveChance(:,3)),y(D_Collapsed.AboveChance(:,3)),color_{iArea},'filled')
[R,P]= corrcoef([x;y]');

mdl{iArea} = fitlm(x,y,'RobustOpts','on' );
pFit(iArea)   = mdl{iArea}.Coefficients{2,4};
%     %             pFit(iArea)   = coefTest(mdl{iArea});
Rsq(iArea)   =  mdl{iArea}.Rsquared.Adjusted;
slope(iArea) = mdl{iArea}.Coefficients{2,1};
intercept(iArea) = mdl{iArea}.Coefficients{1,1};

plot(mdl{iArea}.VariableInfo.Range{1},(mdl{iArea}.VariableInfo.Range{1})* slope(iArea)+ intercept(iArea),'k','LineWidth',1.5)
[pFit(iArea), Rsq(iArea)]

h = plot(mdl{iArea});
h(1).Marker = 'none';
h(1).Marker = 'none';
h(2).LineWidth=1.5;
h(2).LineStyle='-';
h(2).Color = 'k';
h(3).Color = 'k';
h(4).Color = 'k';
axis([0 6 0.4 1])
legend off
xlabel({'Significant cue decoding span (s)'})
ylabel('Fraction of choices correct')
title(Areas{iArea})

%% Plot mean performance
figure
subplot(3,1,1); hold on
ciplot(nanmean(D_Collapsed.meanDecoding.Short)+nansem(D_Collapsed.meanDecoding.Short),...
    nanmean(D_Collapsed.meanDecoding.Short)-nansem(D_Collapsed.meanDecoding.Short),...
    D_Collapsed.meanDecoding.tbShort,'k')

ciplot(nanmean((D_Collapsed.meanDecoding.ShortCIL+D_Collapsed.meanDecoding.ShortCIH)./2) + nansem((D_Collapsed.meanDecoding.ShortCIL+D_Collapsed.meanDecoding.ShortCIH)./2),...
    nanmean((D_Collapsed.meanDecoding.ShortCIL+D_Collapsed.meanDecoding.ShortCIH)./2) - nansem((D_Collapsed.meanDecoding.ShortCIL+D_Collapsed.meanDecoding.ShortCIH)./2),...
    D_Collapsed.meanDecoding.tbShort,'r')
plot([0 16], [0.5 0.5],':k')
axis([0 16 0.4 0.8])
title('4s  delay trials')
subplot(3,1,2); hold on

ciplot(nanmean(D_Collapsed.meanDecoding.Medium)+nansem(D_Collapsed.meanDecoding.Medium),...
    nanmean(D_Collapsed.meanDecoding.Medium)-nansem(D_Collapsed.meanDecoding.Medium),...
    D_Collapsed.meanDecoding.tbMedium,'k')
% ciplot(nanmean((D_Collapsed.meanDecoding.MediumCIL+D_Collapsed.meanDecoding.MediumCIH)./2) + nansem((D_Collapsed.meanDecoding.MediumCIL+D_Collapsed.meanDecoding.MediumCIH)./2),...
%        nanmean((D_Collapsed.meanDecoding.MediumCIL+D_Collapsed.meanDecoding.MediumCIH)./2) - nansem((D_Collapsed.meanDecoding.MediumCIL+D_Collapsed.meanDecoding.MediumCIH)./2),...
%        D_Collapsed.meanDecoding.tbMedium,'r')
plot([0 16], [0.5 0.5],':k')
axis([0 16 0.4 0.8])
title('8s  delay trials')
ylabel('Fraction correct decoding')

subplot(3,1,3); hold on
ciplot(nanmean(D_Collapsed.meanDecoding.Long)+nansem(D_Collapsed.meanDecoding.Long),...
    nanmean(D_Collapsed.meanDecoding.Long)-nansem(D_Collapsed.meanDecoding.Long),...
    D_Collapsed.meanDecoding.tbLong,'k')

%  ciplot(nanmean((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2) + nansem((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2),...
%        nanmean((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2) - nansem((D_Collapsed.meanDecoding.LongCIL+D_Collapsed.meanDecoding.LongCIH)./2),...
%        D_Collapsed.meanDecoding.tbLong,'r')
plot([0 16], [0.5 0.5],':k')
axis([0 16 0.4 0.8])
title('16s  delay trials')
xlabel('Time (s)')
%%
for iFile = 1:length(fileList)
    figure;
    x = D_Collapsed.meanDecoding.tbLong';
    y=D_Collapsed.meanDecoding.Long(iFile ,:); y(y<0.5) = 0.5;
    f = fit(x,y','exp2');
    plot(f,x,y)
end

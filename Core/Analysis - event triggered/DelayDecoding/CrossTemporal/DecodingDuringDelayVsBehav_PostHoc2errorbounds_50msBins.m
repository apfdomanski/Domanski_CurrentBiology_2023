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
%     pat = '/Volumes/HDD2/DNMTP/raw/';
    pat = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw/';

else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);

fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'KrzysztofLONG1_Events.mat','KrzysztofLONG2_Events.mat','MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
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
%% batch process performance
for iFile =1:length(fileList)
    %% load files
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);
    
    %% Get the assembly membership info
    fname=strtok(fileListAss(iFile).name,'_');
    fprintf('Loading run %d/%d %s ...\n',iFile,length(fileList),fname)
    % for assemblies describing whole task period (cue  - reward)
    switch AssemblyChoice
        case 1
            load(sprintf('%s%s%s_iFR50_FSC.mat',pat2,filesep,fname),'units','nu');
            %             load(sprintf('%s%s%s_iFR50_AssemRes2.mat',pat2,filesep,fname),'usel_out');
            usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
        case 2
            load(sprintf('%s%s%s_iFR50_BehavOnly_FSC.mat',pat2,filesep,fname),'units','nu');
            load(sprintf('%s%s%s_iFR50_BehavOnly_AssemRes2.mat',pat2,filesep,fname),'usel_out');
        case 3
            load(sprintf('%s%s%s_iFR50_Task_FSC.mat',pat2,filesep,fname),'units','nu');
            %             load(sprintf('%s%s%s_iFR50_Task_AssemRes2.mat',pat2,filesep,fname),'usel_out');
            usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
            
    end
    noPFC=nu(1);
    Ass.usel_out = usel_out;
    Ass.usel_out{3} = [Ass.usel_out{1},Ass.usel_out{2}+noPFC];
    Ass.nu = nu; Ass.nu(3) = sum(Ass.nu(1:2));
    Ass.units = units;
    for iArea = 1:3
        Ass.LocalMembers{iArea}    = unique(cell2mat(Ass.units{iArea}));
    end
    
    
    Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3}<=Ass.nu(1)));
    Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3}>Ass.nu(1))-Ass.nu(1));
    for iArea = 1:2
        Ass.JointMembers{iArea} = setdiff(unique(cell2mat(Ass.units{iArea})),Ass.LocalMembers{iArea});
        Ass.NonMembers{iArea} = setdiff(1:Ass.nu(iArea),[Ass.LocalMembers{iArea},Ass.JointMembers{iArea}]);
    end
    for iArea = 1:2
        Ass.NonMembers{iArea}   = Ass.usel_out{iArea}(Ass.NonMembers{iArea});
        Ass.LocalMembers{iArea} = Ass.usel_out{iArea}(Ass.LocalMembers{iArea});
        Ass.JointMembers{iArea} = Ass.usel_out{iArea}(Ass.JointMembers{iArea});
    end
    
    clear units usel_out
    %% Batch process behavioural performance
    L = [length(t.Short.ChoicePress_LeftCorrect)./(length(t.Short.ChoicePress_LeftCorrect)+length(t.Short.ChoicePress_LeftError)) , ...
        length(t.Medium.ChoicePress_LeftCorrect)./(length(t.Medium.ChoicePress_LeftCorrect)+length(t.Medium.ChoicePress_LeftError)) , ...
        length(t.Long.ChoicePress_LeftCorrect)./(length(t.Long.ChoicePress_LeftCorrect)+length(t.Long.ChoicePress_LeftError))];
    
    
    R = [length(t.Short.ChoicePress_RightCorrect)./(length(t.Short.ChoicePress_RightCorrect)+length(t.Short.ChoicePress_RightError)) , ...
        length(t.Medium.ChoicePress_RightCorrect)./(length(t.Medium.ChoicePress_RightCorrect)+length(t.Medium.ChoicePress_RightError)) , ...
        length(t.Long.ChoicePress_RightCorrect)./(length(t.Long.ChoicePress_RightCorrect)+length(t.Long.ChoicePress_RightError))];
    
    D.Accuracy(iFile,:) = (L+R)./2;
    D.AccuracyL(iFile,:) = L;
    D.AccuracyR(iFile,:) = R;
    
    C_ =  [length(t.Short.ChoicePress_LeftCorrect)  + length(t.Short.ChoicePress_RightCorrect),...
        length(t.Medium.ChoicePress_LeftCorrect) + length(t.Medium.ChoicePress_RightCorrect),...
        length(t.Long.ChoicePress_LeftCorrect)   + length(t.Long.ChoicePress_RightCorrect)];
    
    E_ =  [length(t.Short.ChoicePress_LeftError)    + length(t.Short.ChoicePress_RightError),...
        length(t.Medium.ChoicePress_LeftError)   + length(t.Medium.ChoicePress_RightError),...
        length(t.Long.ChoicePress_LeftError)     + length(t.Long.ChoicePress_RightError)];
    
    D.pCorr(iFile,:) = C_./(C_+E_);
    D.AboveChance(iFile,:) = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;

    C_ =  [length(t.Short.ChoicePress_LeftCorrect),...
        length(t.Medium.ChoicePress_LeftCorrect),...
        length(t.Long.ChoicePress_LeftCorrect)];
    
    E_ =  [length(t.Short.ChoicePress_LeftError),...
        length(t.Medium.ChoicePress_LeftError),...
        length(t.Long.ChoicePress_LeftError)];
       
    D.pCorrL(iFile,:) = C_./(C_+E_);
    D.AboveChanceL(iFile,:) = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;


    C_ =  [length(t.Short.ChoicePress_RightCorrect),...
        length(t.Medium.ChoicePress_RightCorrect),...
        length(t.Long.ChoicePress_RightCorrect)];
    
    E_ =  [length(t.Short.ChoicePress_RightError),...
        length(t.Medium.ChoicePress_RightError),...
        length(t.Long.ChoicePress_RightError)];
       
    D.pCorrR(iFile,:) = C_./(C_+E_);
    D.AboveChanceR(iFile,:) = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;
    clear L R C_ E_
    
end
%% batch import unit data for meta-analysis and plotting
UnitSelection = 'ErrorBounds_2500msKernel'; %{'pairs','groups','all','groups_TrialDraw','ErrorBounds'}


clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
%     pat = '/Volumes/HDD2/DNMTP/raw/';
    pat = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw/';

%     KDE_binsTaskonly_50mskernel
end

fileList = dir([pat 'MixedSelectivity' filesep 'DelayDecodingVsBehav' filesep '*' Target '*' 'PFC' '*' UnitSelection '.mat']);

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
        D_Collapsed.meanDecoding.D{iArea}{iDelay}(:,:,iFile)     = D_{iArea}{iFile}.D_{1}{iDelay}(:,:,iFile);
        D_Collapsed.meanDecoding.D_CIh{iArea}{iDelay}(:,:,iFile) = D_{iArea}{iFile}.D_CIh{1}{iDelay}(:,:,iFile);
        D_Collapsed.meanDecoding.D_CIl{iArea}{iDelay}(:,:,iFile) = D_{iArea}{iFile}.D_CIl{1}{iDelay}(:,:,iFile);
      
        D_Collapsed.meanDecoding.tbShort = D_{iArea}{iFile}.tbShort;
        D_Collapsed.meanDecoding.tbMedium = D_{iArea}{iFile}.tbMedium;
        D_Collapsed.meanDecoding.tbLong = D_{iArea}{iFile}.tbLong;
    end
    
%     D_Collapsed.SigDecodeTimeDelayMeans(iFile) = nanmean([D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime ./diff(tlimsShort), ...
%         D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeTime ./diff(tlimsMedium), ...
%         D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeTime   ./diff(tlimsLong)]);
%     D_Collapsed.SigDecodeFracDrawsDelayMeans(iFile) = nanmean([sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsShort) , ...
%         sum(D_{iArea}{iFile}.D.Medium.DelaySpan.SigDecodeFracDraws)./diff(tlimsMedium), ...
%         sum(D_{iArea}{iFile}.D.Long.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong)]);
%     
%     D_Collapsed.AccuracyDelayMeans(iFile) = mean(D_{1}{iFile}.D.Accuracy);
%     D_Collapsed.pCorrDelayMeans(iFile) = mean(D_{1}{iFile}.D.pCorr);
%     
%     D_Collapsed.SigDecodeTime(iFile)        = D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeTime./diff(tlimsLong);
%     D_Collapsed.SigDecodeFracDraws(iFile)   = sum(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)./diff(tlimsLong);
%     D_Collapsed.LastSigTime(iFile)          = max([0,find(D_{iArea}{iFile}.D.Short.DelaySpan.SigDecodeFracDraws)])*bw;
%     
%     D_Collapsed.Accuracy(iFile) = D_{1}{iFile}.D.Accuracy(1);
%     D_Collapsed.pCorr(iFile)    = D_{1}{iFile}.D.pCorr(1);
%     
%     D_Collapsed.AboveChance(iFile,:) = D_{1}{iFile}.D.AboveChance;
end
%%
iDelay = 2;
figure; hold on

A = D_Collapsed.meanDecoding.D{iArea}{iDelay} > D_Collapsed.meanDecoding.D_CIh{iArea}{iDelay};
% % imagesc(sum(A,3)./size(A,3)>0.5)
% imagesc(sum(A,3)./size(A,3))
figure;
for iFile = 1:length(fileList)
   subplot(5,3,iFile)
   contour(A(:,:,iFile),'LineWidth',1.5,'LineColor','k')
   colormap(gray)
   caxis([0 1]);

end
%% plot Correct - mean across sessions
iArea = 1;
figure
% Ticks_ = {[-5 0 4 9],[-5 0 4 8 13],[-5 0 4 8 16 21]}';
Ticks_ = {[0:4:4],[0:4:8],[0:4:16]}';
for iDelay =1:length(Delays_)
    subplot(1,3,iDelay); hold on 
    eval(sprintf('tlims_X = tb%s;',Delays_{iDelay}));
    aboveChance_ = D.AboveChance(:,iDelay);
%     aboveChance_ = ones(size(aboveChance_));
    decoding_ = D_Collapsed.meanDecoding.D{iArea}{iDelay}(:,:,aboveChance_);
    for i=1:size(decoding_,3)
        decoding_(:,:,i) = smooth2a(decoding_(:,:,i) ,5,5);
    end
    decoding_ = nanmean(decoding_,3);
    decoding_ = smooth2a(decoding_ ,5,5);
    imagesc(tlims_X,tlims_X,decoding_);
    plot([0 0],[-5 0],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w','LineWidth',1.5)  
    plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w','LineWidth',1.5);
    plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w','LineWidth',1.5);
    plot([-5 0],[0 0],'w') 
    plot([-5 max(tlims_X)+5],[-5 max(tlims_X)+5],'w','LineWidth',1.5) 

%     rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w','LineWidth',1.5)
%     rectangle('Position',[-5 -5  max(tlims_X)+5 max(tlims_X)+5],'EdgeColor','w','LineWidth',1.5)
    
%     scatter(0,0,50,'k')
%     scatter(max(tlims_X)-5,max(tlims_X)-5,50,'k')
    
    set(gca,'XTick',Ticks_{iDelay},'YTick',Ticks_{iDelay})
    axis([-5 21 -5 21]); 
    colormap(jet)
%     set(gca,'Colormap',jet); 
%     caxis([0.4 0.6])
    axis square
    ylabel('Training time (s)')
    xlabel('Testing time (s)')
end

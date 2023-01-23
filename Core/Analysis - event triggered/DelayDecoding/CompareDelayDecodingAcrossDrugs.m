%% Preamble

clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};

Targets={'LONG','CP55940','AM251','URB597','physostigmine','scopolamine'}; %,'URB597'
Targets_={'Control','CP55940','AM251','URB597','physostigmine','scopolamine'}; %,'URB597'

Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};


tlimsShort=[-5 9];
tlimsMedium=[-5 13];
tlimsLong=[-5 21];

shift = 0;
plotOnline = false;
bw    = 0.05;
Nbs = 500;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
% tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
    patRes = [pat , 'MixedSelectivity\DelayDecodingCrossTemporal\1000bsWithErrors\'];
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
    patRes = [pat , 'MixedSelectivity/DelayDecodingVsBehav/1000bsWithErrors/'];
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end
UnitSelection = 'ErrorBoundsWithErr'; %{'pairs','groups','all','groups_TrialDraw','ErrorBounds','ErrorBoundsWithErr'}
%% Bulk import all conditions
for iTarget=1:length(Targets)
    Target = Targets{iTarget};


    clear Av
    warning ('off')
    if ispc
        pat = 'C:\Analysis\AssemblyAnalysis\raw\';
    elseif ismac
        pat = '/Volumes/HDD2/DNMTP/raw/';

    end

    % % fileList = dir([pat2 '*' Target '*' 'PFC' '*' UnitSelection '.mat']);
    fileList = dir([patRes '*' Target '*' 'PFC' '*' UnitSelection '.mat']);

    clear D_
    for iArea = 1%:3
        for iFile = 1:length(fileList)
            %         try
            
            fname=strtok(fileList(iFile).name,'_');
            %             fnIn = sprintf('%sMixedSelectivity%s%s%s%s_%s_DelayDecoding_UnitSpan_%s.mat',pat,filesep,'DelayDecodingVsBehav',filesep,fname,Areas{iArea},UnitSelection);
            fnIn = sprintf('%s%s_%s_DelayDecoding_UnitSpan_%s.mat',patRes,fname,Areas{iArea},UnitSelection);
            D_{iArea}{iFile} = load(fnIn);
            fprintf('Done.\n')
            %         end
        end
    end
    
    %% Collapse for each delays
    clear D_Collapsed
    
    
    for iFile = 1:length(fileList)
        for iDelay =1:3
            %         temp =  D_{iArea}{iFile}.D_{1}{iDelay};
            %         iFile_ = iFile;
            %         if iFile>1
            %             iFile_ = find(sum(squeeze(sum(temp,1)),1)>0);
            %         end
            %         [iFile ,iFile_]
            %         if iFile_ ~= iFile
            %             iFile_ =
            %         D_Collapsed.meanDecoding.D{iArea}{iDelay}(:,:,iFile_)     = D_{iArea}{iFile}.D_{1}{iDelay}(:,:,iFile_);
            %         D_Collapsed.meanDecoding.D_CIh{iArea}{iDelay}(:,:,iFile_) = D_{iArea}{iFile}.D_CIh{1}{iDelay}(:,:,iFile_);
            %         D_Collapsed.meanDecoding.D_CIl{iArea}{iDelay}(:,:,iFile_) = D_{iArea}{iFile}.D_CIl{1}{iDelay}(:,:,iFile_);
            
            D_Collapsed.meanDecoding.D{iArea}{iDelay}(:,:,iFile)     = D_{iArea}{iFile}.D_{iArea}{iDelay};
            D_Collapsed.meanDecoding.D_CIh{iArea}{iDelay}(:,:,iFile) = D_{iArea}{iFile}.D_CIh{iArea}{iDelay};
            D_Collapsed.meanDecoding.D_CIl{iArea}{iDelay}(:,:,iFile) = D_{iArea}{iFile}.D_CIl{iArea}{iDelay};
            
            D_Collapsed.meanDecoding.D_Err{iArea}{iDelay}(:,:,iFile)    = D_{iArea}{iFile}.D_Err{iArea}{iDelay};
            D_Collapsed.meanDecoding.D_ErrCIh{iArea}{iDelay}(:,:,iFile) = D_{iArea}{iFile}.D_ErrCIh{iArea}{iDelay};
            D_Collapsed.meanDecoding.D_ErrCIl{iArea}{iDelay}(:,:,iFile) = D_{iArea}{iFile}.D_ErrCIl{iArea}{iDelay};
            
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

    D__{iTarget}{iArea} = D_Collapsed;
end
clear D_Collapsed D_ 
%% Plot all conditions - mean across all recordings
smooth_ = 10;
iArea = 1;
figure
% Ticks_ = {[-5 0 4 9],[-5 0 4 8 13],[-5 0 4 8 16 21]}';
Ticks_ = {[0:4:4],[0:4:8],[0:4:16]}';

for iTarget = 1:length(Targets)
    Target = Targets{iTarget};
    D_Collapsed = D__{iTarget};
    for iDelay=1:3
        
        subplot(3,length(Targets),iTarget+length(Targets)*(iDelay-1)); hold on
        if iDelay==1
           title(Targets_{iTarget}) 
        end
        eval(sprintf('tlims_X = tb%s;',Delays_{iDelay}));
        eval(sprintf('tlims_X2 = D_Collapsed{iArea}.meanDecoding.tb%s;',Delays_{iDelay}));
%         aboveChance_ = D.AboveChance(:,iDelay);
        aboveChance_ = ones(length(fileList));
        decoding_    = D_Collapsed{iArea}.meanDecoding.D{iArea}{iDelay};
        
        A =  D_Collapsed{iArea}.meanDecoding.D{iArea}{iDelay};
        B =  (D_Collapsed{iArea}.meanDecoding.D_CIl{iArea}{iDelay} + D_Collapsed{iArea}.meanDecoding.D_CIh{iArea}{iDelay})./2;
        
        
        for i=1:size(A,3)
            A(:,:,i) = smooth2a(A(:,:,i),smooth_,smooth_);
            B(:,:,i) = smooth2a(B(:,:,i),smooth_,smooth_);
        end
        A = A>B;
        %     A            = D_Collapsed.meanDecoding.D{iArea}{iDelay} > D_Collapsed.meanDecoding.D_CIl{iArea}{iDelay};
        
        idxprune = sum(squeeze(sum(decoding_,1)))==0;
        decoding_(:,:,idxprune) = [];
        A(:,:,idxprune) = [];
        aboveChance_(idxprune)  = [];
        
        decoding_ = decoding_(:,:,aboveChance_);
        A         = A(:,:,aboveChance_);
        
        A = sum(double(A),3)./size(A,3);
        %     A=A>0;
        %     A = smooth2a(A,smooth_,smooth_);
        A=A>0.6;
        for i=1:size(decoding_,3)
            decoding_(:,:,i) = smooth2a(decoding_(:,:,i) ,smooth_,smooth_);
        end
        decoding_ = nanmean(decoding_,3);
        %     decoding_ = smooth2a(decoding_ ,smooth_,smooth_);
        imagesc(tlims_X,tlims_X,decoding_,'alphadata',0.3+A);
        contour(tlims_X2,tlims_X2,A,'linewidth',1.5,'linecolor','k')
        plot([0 0],[-5 0],'w','LineWidth',1.5);
        plot([-5 0],[0 0],'w','LineWidth',1.5)
        plot([max(tlims_X)-5,max(tlims_X)-5],[max(tlims_X)-5 max(tlims_X)],'w','LineWidth',1.5);
        plot([max(tlims_X)-5,max(tlims_X)],[max(tlims_X)-5 max(tlims_X)-5],'w','LineWidth',1.5);
        plot([-5 0],[0 0],'w')
        %     plot([-5 max(tlims_X)+5],[-5 max(tlims_X)+5],'w','LineWidth',1.5)
        
        rectangle('Position',[0 0 max(tlims_X)-5 max(tlims_X)-5],'EdgeColor','w','LineWidth',1.5)
        rectangle('Position',[-5 -5  max(tlims_X)+5 max(tlims_X)+5],'EdgeColor','w','LineWidth',1.5)
        
        scatter(0,0,50,'k')
        scatter(max(tlims_X)-5,max(tlims_X)-5,50,'k')
        
        set(gca,'XTick',Ticks_{iDelay},'YTick',Ticks_{iDelay})
        axis([-5 21 -5 21]);
        colormap(jet)
        %     set(gca,'Colormap',jet);
        caxis([0.4 0.6])
        axis square
        if iTarget ==1
            ylabel('Training time (s)')
        end
        if iDelay==3
            xlabel('Testing time (s)') 
        end
        
    end
end
%% Plot the diagonal
clear D
smth=10;
for iTarget =1:length(Targets)
    D_Collapsed =  D__{iTarget}{iArea}.meanDecoding;    
    for iDelay = 1:3
        decoding_  = D_Collapsed.D{iArea}{iDelay};
        idxprune = sum(squeeze(sum(decoding_,1)))==0;
        decoding_(:,:,idxprune) = [];
        D{iDelay}(:,iTarget) = diag(smooth2a(nanmean(decoding_,3),smth,smth));
    end
end

figure;
subplot(1,3,1); hold on 
plot(D_Collapsed.tbShort,D{1}')
plot([-5 21],[0.5 0.5],':k','HandleVisibility','off')
plot([0 0],[0.45 0.55],'g','HandleVisibility','off')
plot([4 4],[0.45 0.55],'r','HandleVisibility','off')
axis([-5 21 0.3 0.7])
subplot(1,3,2); hold on 
plot(D_Collapsed.tbMedium,D{2}')
plot([-5 21],[0.5 0.5],':k','HandleVisibility','off')
plot([0 0],[0.45 0.55],'g','HandleVisibility','off')
plot([8 8],[0.45 0.55],'r','HandleVisibility','off')
axis([-5 21 0.3 0.7])
subplot(1,3,3); hold on 
plot(D_Collapsed.tbLong,D{3}')
plot([-5 21],[0.5 0.5],':k','HandleVisibility','off')
plot([0 0],[0.45 0.55],'g','HandleVisibility','off')
plot([16 16],[0.45 0.55],'r','HandleVisibility','off')
axis([-5 21 0.3 0.7])
legend(Targets,'Location','South','Box','off')
    

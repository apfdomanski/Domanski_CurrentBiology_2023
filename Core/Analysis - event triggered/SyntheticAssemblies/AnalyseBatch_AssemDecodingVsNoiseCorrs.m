% %%%%%% PREAMBLE %%%%%%
clear

Delays_ = {'Short','Medium','Long'};
Target = 'LONG';
Areas = {'PFC','HP','Joint'};

bw=0.05; tlimsAll = [-5 5]; tbAll = tlimsAll(1):bw:tlimsAll(2);

maxRange=20;    % Upper limit of assembly size to explore
noTrials = 5;   % How many trials to consider simultanrously
nReps = 50;     % how many repeated draws among trials to perform


InfoCriteria = 'max'; % 'mean','max'

if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
%         pat = '/Volumes/HDD2/DNMTP/raw/';
        pat = '/Volumes/Aleks Data Portable/AssemblyAnalysis/raw/';
elseif isunix
    pat ='/Volumes/Data/DNMTP/raw';
end

fileList = dir([pat 'SyntheticAssemblies' filesep '*' Target '*_AssDecodingVsNoiseCorrs3.mat']);


%% Batch import
for iFile =1:length(fileList)
    
    %fnIn = fullfile(pat, 'SyntheticAssemblies', fileList(iFile).name);
    fnIn = fullfile(pat, 'SyntheticAssemblies',fileList(iFile).name);
    D{iFile} = load(fnIn,'D_');
    
end
clear temp
%% Summary stats
for s=1:3
    D_.Ass_real.score{s}=[];
    D_.Ass_real.scoreShuf{s}=[];
    D_.Ass_real.score_NoNoiseCorrs{s}=[];
    D_.Ass_real.Deltascore{s}=[];
    
    for iFile =1:length(D)
        D_.Ass_real.score{s}=[ D_.Ass_real.score{s};D{iFile}.D_.AssReal.score{s}];
        D_.Ass_real.scoreShuf{s}=[ D_.Ass_real.scoreShuf{s};D{iFile}.D_.AssReal.score_shuffled{s}];
%         try
            D_.Ass_real.score_NoNoiseCorrs{s}=[ D_.Ass_real.score_NoNoiseCorrs{s};D{iFile}.D_.AssReal_noNoiseCorrs.score{s}];
            
            delta_ = D{iFile}.D_.AssReal_noNoiseCorrs.score{s}(:,2) - D{iFile}.D_.AssReal.score{s}(:,2);
            D_.Ass_real.Deltascore{s}=[ D_.Ass_real.Deltascore{s};...
                [D{iFile}.D_.AssReal_noNoiseCorrs.score{s}(:,1) ,delta_]];
%         end
        
        
    end
end
%% Summary stats
s=3;
for s_=1:3
    D_.Ass_real.CrossShuffleNoNoiseCorrs{s_}=[];
    D_.Ass_real.CrossShuffle_Deltascore{s_}=[];
end
for iFile =1:length(D)
    
    for s_ =1:3
        try
            D_.Ass_real.CrossShuffleNoNoiseCorrs{s_}=[ D_.Ass_real.CrossShuffleNoNoiseCorrs{s_};  D{iFile}.D_.AssReal_noNoiseCorrs.score_crossShuffle{s}{s_}];
            
            delta_ =  D{iFile}.D_.AssReal_noNoiseCorrs.score_crossShuffle{s}{s_}(:,2) - D{iFile}.D_.AssReal.score{s}(:,2);
            D_.Ass_real.CrossShuffle_Deltascore{s_}=[ D_.Ass_real.CrossShuffle_Deltascore{s_};...
                [D{iFile}.D_.AssReal_noNoiseCorrs.score{s}(:,1) ,delta_] ];
        end
    end
    
end

%%
figure
for s=1:3
    subplot(1,3,s); hold on
    scatter(D_.Ass_real.score{s}(:,1),D_.Ass_real.score{s}(:,2),10,'b','filled')
    scatter(D_.Ass_real.scoreShuf{s}(:,1),D_.Ass_real.scoreShuf{s}(:,2),10,'r','filled')
    scatter(D_.Ass_real.score_NoNoiseCorrs{s}(:,1),D_.Ass_real.score_NoNoiseCorrs{s}(:,2),10,'g','filled')
    for n=2:20
        x = D_.Ass_real.score{s};
        m = nanmean(x(x(:,1)==n,2));
        e = nansem(x(x(:,1)==n,2));
        errorbar(n-.1,m,e,'LineWidth',1.5,'Color','b')
        
        x = D_.Ass_real.scoreShuf{s};
        m = nanmean(x(x(:,1)==n,2));
        e = nansem(x(x(:,1)==n,2));
        errorbar(n,m,e,'LineWidth',1.5,'Color','r')
        
        x = D_.Ass_real.score_NoNoiseCorrs{s};
        m = nanmean(x(x(:,1)==n,2));
        e = nansem(x(x(:,1)==n,2));
        errorbar(n+.1,m,e,'LineWidth',1.5,'Color','g')
    end
    x = D_.Ass_real.score{s};
    mdl{s_} = fitlm(x(:,1),x(:,2));
    pFit(s_)   = mdl{s_}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s_)   =  mdl{s_}.Rsquared.Adjusted;
    slope(s_) = mdl{s_}.Coefficients{2,1};
    
    h = plot(mdl{s_});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'b';
    h(3).Color = 'b';
    h(4).Color = 'b';

    x = D_.Ass_real.score_NoNoiseCorrs{s};
    mdl{s_} = fitlm(x(:,1),x(:,2));
    pFit(s_)   = mdl{s_}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s_)   =  mdl{s_}.Rsquared.Adjusted;
    slope(s_) = mdl{s_}.Coefficients{2,1};
    
    h = plot(mdl{s_});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'g';
    h(3).Color = 'g';
    h(4).Color = 'g';
    
    x = D_.Ass_real.scoreShuf{s};
    mdl{s_} = fitlm(x(:,1),x(:,2));
    pFit(s_)   = mdl{s_}.Coefficients{2,4};
    %     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s_)   =  mdl{s_}.Rsquared.Adjusted;
    slope(s_) = mdl{s_}.Coefficients{2,1};
    
    h = plot(mdl{s_});
    h(1).Marker = 'none';
    h(1).Marker = 'none';
    h(2).LineWidth=1.5;
    h(2).LineStyle='-';
    h(2).Color = 'r';
    h(3).Color = 'r';
    h(4).Color = 'r';
    
    axis([1.5 15 0 1])
    legend off
end
%%
figure
for s=1:3
    subplot(1,3,s); hold on
    scatter(D_.Ass_real.Deltascore{s}(:,1),D_.Ass_real.Deltascore{s}(:,2),'b')
    for n=2:20
        x = D_.Ass_real.Deltascore{s};
        m = nanmean(x(x(:,1)==n,2));
        e = nansem(x(x(:,1)==n,2));
        errorbar(n,m,e,'LineWidth',1.5,'Color','b')
    end
    axis([1.5 20 -0.2 0.2])
end
%%
s=3;

figure
for s_=1:3
    subplot(1,3,s_); hold on
    scatter(D_.Ass_real.score{s}(:,1),D_.Ass_real.score{s}(:,2),10,'b','filled')
    scatter(D_.Ass_real.scoreShuf{s}(:,1),D_.Ass_real.scoreShuf{s}(:,2),10,'r','filled')
    scatter(D_.Ass_real.CrossShuffleNoNoiseCorrs{s_}(:,1),D_.Ass_real.CrossShuffleNoNoiseCorrs{s_}(:,2),10,'g','filled')
    for n=2:20
        x = D_.Ass_real.score{s};
        m = nanmean(x(x(:,1)==n,2));
        e = nansem(x(x(:,1)==n,2));
        errorbar(n-.2,m,e,'LineWidth',1.5,'Color','b')
        
        x = D_.Ass_real.scoreShuf{s};
        m = nanmean(x(x(:,1)==n,2));
        e = nansem(x(x(:,1)==n,2));
        errorbar(n,m,e,'LineWidth',1.5,'Color','r')
        
        x = D_.Ass_real.CrossShuffleNoNoiseCorrs{s_};
        m = nanmean(x(x(:,1)==n,2));
        e = nansem(x(x(:,1)==n,2));
        errorbar(n+.2,m,e,'LineWidth',1.5,'Color','g')
    end
    axis([2 Inf 0 1])
end
%%
Color_={'k','k','k'};

figure
for s_=1:3
    subplot(1,3,s_); hold on
    scatter(D_.Ass_real.CrossShuffle_Deltascore{s_}(:,1),D_.Ass_real.CrossShuffle_Deltascore{s_}(:,2),'k')
    for n=2:20
        x = D_.Ass_real.CrossShuffle_Deltascore{s_};
        m = nanmean(x(x(:,1)==n,2));
        e = nansem(x(x(:,1)==n,2));
        errorbar(n,m,e,'LineWidth',1.5,'Color','k')
    end
    mdl{s_} = fitlm(x(:,1),x(:,2));
    pFit(s_)   = mdl{s_}.Coefficients{2,4};
%     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s_)   =  mdl{s_}.Rsquared.Adjusted;
    slope(s_) = mdl{s_}.Coefficients{2,1};
    
    h = plot(mdl{s_});
    if pFit(s_)<0.05
        h(1).Marker = 'none';
        h(2).LineWidth=1.5;
        h(2).LineStyle='-';
        h(2).Color = Color_{s_};
        h(3).Color = Color_{s_};
        h(4).Color = Color_{s_};
        
        %     h(3).LineStyle='none';
        %     h(4).LineStyle='none';
    else
        h(1).Marker = 'none';
        h(1).Marker = 'none';
        h(2).LineWidth=1.5;
        h(2).LineStyle='-';
        h(2).Color = Color_{s_};
        h(3).Color = Color_{s_};
        h(4).Color = Color_{s_};
        %     h(3).LineStyle='none';
        %     h(4).LineStyle='none';
    end
    axis([2 20 -0.2 0.2])
    legend off
end


%%
Color_={'b','r','g'};
figure; hold on
for s_=1:3
    x = D_.Ass_real.CrossShuffle_Deltascore{s_}(:,1);
    y = 100*D_.Ass_real.CrossShuffle_Deltascore{s_}(:,2);
    scatter(x,y,10,Color_{s_},'filled')

    mdl{s_} = fitlm(x,y,'RobustOpts','on');
    pFit(s_)   = mdl{s_}.Coefficients{2,4};
%     %             pFit(iArea)   = coefTest(mdl{iArea});
    Rsq(s_)   =  mdl{s_}.Rsquared.Adjusted;
    slope(s_) = mdl{s_}.Coefficients{2,1};
    
    h = plot(mdl{s_});
    if pFit(s_)<0.05
        h(1).Marker = 'none';
        h(2).LineWidth=1.5;
        h(2).LineStyle='-';
        h(2).Color = Color_{s_};
        h(3).Color = Color_{s_};
        h(4).Color = Color_{s_};
        
        %     h(3).LineStyle='none';
        %     h(4).LineStyle='none';
    else
        h(1).Marker = 'none';
        h(1).Marker = 'none';
        h(2).LineWidth=1.5;
        h(2).LineStyle='-';
        h(2).Color = Color_{s_};
        h(3).Color = Color_{s_};
        h(4).Color = Color_{s_};
        %     h(3).LineStyle='none';
        %     h(4).LineStyle='none';
    end
    % a = ezfit(idxFR*bw,idxTS*bw,'poly1'); showfit(a)
    
end
legend off
set(gca,'XTick',[2:7],'YTick',[-5 0 5])
xlabel('Assembly Size (no. Units)')
ylabel({'\delta % correct  with';'shuffled trial labels'})
% axis([1 8 -6 6])
%% Preamble, get file list
clear all; close all
set(0,'DefaultFigureWindowStyle','docked')
flags.plotIndividualChanges = false;
if ispc
%     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\';       % location of the processed sleep assembly activations
%     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - epoch specific FRs\';
    pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - use task mean FRs\';
    pat{2} = 'C:\Analysis\AssemblyAnalysis\raw';           % location of raw spike times
else
    pat{1} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Decoding vs sleep/';
    pat{2} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Task';
end
pat{3} = [pat{2} filesep 'KDE_bins'];                  % location of calculated firing P.rates for Task period
cd(pat{1})
RatList = dir([pat{1} filesep '*_DecodingVSsleep.mat']);

% Ignore= {'JaroslawLONG1'; ...
%          'MiroslawLONG2'; ...
%          'NorbertLONG2' ; ...
%          'OnufryLONG2'};
% for iList = 1:length(Ignore)
%     RatList = RatList(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
%            {RatList.name},repmat({Ignore(iList)},1,length(RatList)),'UniformOutput',false)));
% end

noClusts  = 6;
ClassNames = {'Pre-Sample',...
              'Sample',...
              'Early Delay',...
              'Late Delay',...
              'Sample/Choice',...
              'Post-Choice'};
ClassColors = {[0.7 0.9 0.9],...
               [0.1 0.9 0.1],...
               [0.9 0.7 0.9],...
               [0.9 0.7 0.9],...
               [0.9 0.6 0.1],...
               [0.9 0.1 0.1]};
%% Collate feature vectors
Group = struct;
for iRat =1:length(RatList)
    disp([sprintf('Working on rat %d of %d',iRat, length(RatList)), ' (' RatList(iRat).name ')'])
    load([pat{1} RatList(iRat).name],'P','D','FRtrials','FAtrials');
    
    for s= 1:3
        for iOutcome = 1:2
            if ~isempty(FAtrials.FSC_Mean{s}) %i.e. as long as there are detected assemblies of this type
                Group.Assem.AllFSCmean{s}{iOutcome}{1,iRat} = FAtrials.FSC_Mean{s}{iOutcome};
                % if ~isnan(nanmax(D.Assem.TS{s}))
                Group.Assem.AllTS{s}{iRat} = D.Assem.TS{s}; 
                % NB, for  Significant decoders: 
                Group.Assem.AllTSsig{s}{iRat} = tpdf(D.Assem.TS{s},numel(FRtrials.evt0)-2);
            else
                Group.Assem.AllFSCmean{s}{iOutcome}{1,iRat} = [];
                Group.Assem.AllTS{s}{iRat} = [];
                Group.Assem.AllTSsig{s}{iRat} = [];
            end
                Group.Units.AllFRmean{s}{iOutcome}{1,iRat}  = FRtrials.iFR0_Mean{s}{iOutcome};
                Group.Units.AllTS{s}{iRat}                  = D.units.TS{s};
                Group.Units.AllTSsig{s}{iRat}               = tpdf(D.units.TS{s},numel(FRtrials.evt0)-2);
        end
    end
end
    
clear iRat s iOutcome D FAtrials FRtrials
%% Collapse L/R specific sequences
for s=1:3
    for iOutcome = 1:2
        Group.Units.AllFRmean{s}{iOutcome}  = cell2mat(Group.Units.AllFRmean{s}{iOutcome});
        Group.Assem.AllFSCmean{s}{iOutcome} = cell2mat(Group.Assem.AllFSCmean{s}{iOutcome});
    end
    Group.Units.AllFRmeanCollapse{s} = cell2mat(Group.Units.AllFRmean{s}');
    Group.Assem.AllFSCmeanCollapse{s} = cell2mat(Group.Assem.AllFSCmean{s}');
end

for s=1:3
    Group.Units.AllTSCollapse{s} =  cell2mat(Group.Units.AllTS{s});
    Group.Assem.AllTSCollapse{s} =  cell2mat(Group.Assem.AllTS{s});
    Group.Units.AllTSsigCollapse{s} = cell2mat(Group.Units.AllTSsig{s});
    Group.Assem.AllTSsigCollapse{s} = cell2mat(Group.Assem.AllTSsig{s});

    % Times when Units/Assems are significant decoders
    Group.Units.AllTSsigCollapse{s} = Group.Units.AllTSsigCollapse{s}<0.05;
    Group.Assem.AllTSsigCollapse{s} = Group.Assem.AllTSsigCollapse{s}<0.05;
    
    % Units/Assems that are significant decoders above some time threshold
    cutoff = 0; % minimum significance cutoff in seconds
    Group.Units.SigDecoders{s} = find(sum(Group.Units.AllTSsigCollapse{s})>cutoff/0.05);
    Group.Assem.SigDecoders{s} = find(sum(Group.Assem.AllTSsigCollapse{s})>cutoff/0.05);
    % bins = 0:600 ;
    % x=cumsum(histc(sum(Group.Units.AllTSsigCollapse{s}),bins));
    % figure; plot(bins*0.05,x./max(x))
end
clear s iOutcome cutoff bins x
%% example plot of significant units/assems
figure
subplot(1,2,1)
temp = [Group.Units.AllTSCollapse{1},Group.Units.AllTSCollapse{2}];
temp(~[Group.Units.AllTSsigCollapse{1},Group.Units.AllTSsigCollapse{2}])=0;
imagesc(rot90(temp))
colormap((hot ))
axis off
caxis([0 20])
subplot(1,2,2)
temp = [Group.Assem.AllTSCollapse{1},Group.Assem.AllTSCollapse{2},Group.Assem.AllTSCollapse{3}];
temp(~[Group.Assem.AllTSsigCollapse{1},Group.Assem.AllTSsigCollapse{2},Group.Assem.AllTSsigCollapse{3}])=0;
imagesc(rot90(temp))
colormap([1 1 1;hot])
axis off
caxis([0 20])
clear temp
%% Make feature vectors
for s=1:3
%     Group.AssemFeatures{s} = [Group.Assem.AllFSCmeanCollapse{s};Group.Assem.AllTSCollapse{s}];
%     Group.UnitFeatures{s} = [Group.Units.AllFRmeanCollapse{s};Group.Units.AllTSCollapse{s}];
%     Group.AssemFeatures{s} = [Group.Assem.AllFSCmeanCollapse{s}];
%     Group.UnitFeatures{s} = [Group.Units.AllFRmeanCollapse{s}];    


    % All units/Assems
    Group.AssemFeatures{s} = Group.Assem.AllTSCollapse{s};
    % Group.UnitFeatures{s}  = Group.Units.AllTSCollapse{s};
    
    % restrict Units to significant decoders only... mask assems later
%     Group.AssemFeatures{s} = Group.Assem.AllTSCollapse{s}(:,Group.Assem.SigDecoders{s});
    Group.UnitFeatures{s}  = Group.Units.AllTSCollapse{s}(:,Group.Units.SigDecoders{s});
%     temp_x = (1:P.Ltr*4)*P.bw;

end
clear s
%% Distribution of decoding scores for each Unit/Assem type
% bins = 0:0.5:20;
bins = 0:0.01:3;
col_ = {'b','r','g'};
figure
% subplot(2,1,1); 
hold on
for s=1:2
%     plot(bins, cumsum(histc(max(Group.UnitFeatures{s}),bins))./size(Group.UnitFeatures{s},2))
%     plot(bins, smooth_hist(histc(max(Group.UnitFeatures{s}),bins))./size(Group.UnitFeatures{s},2))
    plot(bins, cumsum(histc(sum(Group.UnitFeatures{s})./size(Group.UnitFeatures{s},1),bins))./size(Group.UnitFeatures{s},2),'Color',col_{s},'LineWidth',2,'LineStyle',':')
end
legend(P.names{1:2})
title('Distribution of unit decoding scores')

% subplot(2,1,2); hold on
for s=1:3
%     plot(bins, cumsum(histc(max(Group.UnitFeatures{s}),bins))./size(Group.UnitFeatures{s},2))
%     plot(bins, smooth_hist(histc(max(Group.UnitFeatures{s}),bins))./size(Group.UnitFeatures{s},2))
    plot(bins, cumsum(histc(sum(Group.AssemFeatures{s})./size(Group.AssemFeatures{s},1),bins))./size(Group.AssemFeatures{s},2),'Color',col_{s},'LineWidth',2)
end
% legend(P.names{1:3})
legend('PFC units','HP units','PFC Assemblies','HP Assemblies','Joint Assemblies','Location','southeast'); legend  boxoff
title('Distribution of assembly  decoding scores')

xlabel('Decoding capacity (mean t-score)')
ylabel('Fraction of population')
% clear s bins
%% Fraction of time spent significant
SigHist_Units=cell(1,3); SigHist_Assems=cell(1,3);
P.BonfCorrect = false;
bins = 0:0.01:1;
for s= 1:3
    SigHist_Units{s} = [];     SigHist_Assems{s} = [];

    for iRat = 1:length(RatList)
        if P.BonfCorrect
            temp = Group.Units.AllTSsig{s}{iRat}<(0.05/size(Group.Units.AllTSsig{s}{iRat},1));  
        else
            temp = Group.Units.AllTSsig{s}{iRat}<0.05;
        end
        
        temp = sum(temp)./size(Group.Units.AllTSsig{s}{iRat},1);
        SigHist_Units{s} = [SigHist_Units{s} ; cumsum(histc(temp,bins))./size(Group.Units.AllTSsig{s}{iRat},2)];
        if P.BonfCorrect
            temp = Group.Assem.AllTSsig{s}{iRat}<(0.05/size(Group.Units.AllTSsig{s}{iRat},1));  
        else
            temp = Group.Assem.AllTSsig{s}{iRat}<0.05;
        end
        temp = sum(temp)./size(Group.Assem.AllTSsig{s}{iRat},1);
        SigHist_Assems{s} = [SigHist_Assems{s} ; cumsum(histc(temp,bins))./size(Group.Assem.AllTSsig{s}{iRat},2)];
    end
    SigHist_Units{s}(isnan(max(SigHist_Units{s},[],2)),:)   = [];
    SigHist_Assems{s}(isnan(max(SigHist_Assems{s},[],2)),:) = [];
end
bins =bins*P.Ltr*2*P.bw;
col_ = {'r' 'b' 'g'}
figure('color','w')
subplot(1,2,1); hold on
for s=1:2
    plot(bins,nanmean(SigHist_Units{s}),'color',col_{s})
    plot(bins,nanmean(SigHist_Units{s})+nansem(SigHist_Units{s}),'LineStyle',':','color',col_{s})
    plot(bins,nanmean(SigHist_Units{s})-nansem(SigHist_Units{s}),'LineStyle',':','color',col_{s})
    axis([0 max(bins) 0 1])  
end
text(20,0.2,'PFC','color',col_{1})
text(20,0.15,'HP','color',col_{2})
ylabel('Fraction of population')
xlabel('Mean time spent as significant decoders (s)')
title('Units')
subplot(1,2,2); hold on
for s=1:3
    plot(bins,nanmean(SigHist_Assems{s}),'color',col_{s})
    plot(bins,nanmean(SigHist_Assems{s})+nansem(SigHist_Assems{s}),'LineStyle',':','color',col_{s})
    plot(bins,nanmean(SigHist_Assems{s})-nansem(SigHist_Assems{s}),'LineStyle',':','color',col_{s})
    axis([0 max(bins) 0 1])
end
text(20,0.2,'PFC','color',col_{1})
text(20,0.15,'HP','color',col_{2})
text(20,0.1,'Joint PFC-HP','color',col_{3})
title('Assemblies')
xlabel('Mean time spent as significant decoders (s)')
%% PCA on feature vectors - Decoding

s=3
cluster_input = Group.UnitFeatures{s};

[PCA_.coef, PCA_.scores, PCA_.variances, PCA_.t2]= princomp(cluster_input);
PCA_.var_explained =100*(PCA_.variances/sum(PCA_.variances)); 

temp_x = (1:P.Ltr*2)*P.bw

figure('name',P.names{s},'Color','w'); 
subplot(1,2,1)
    plot(1:10,PCA_.var_explained(1:10),'k','LineWidth',2);
    box off
    xlabel('PC no.')
    ylabel('% var explained')

subplot(1,2,2); hold on
    plot([10 10],[-50 0],'g')
    plot([20 20],[-50 0],'r')
    % plot(temp_x,staggerplot(PCA_.coef(1:end-1,1:10),0,-.5),'color',[0 0 0],'linewidth',2);
    for PCid=1:10
        clear temp
        temp = zscore(PCA_.scores(:,PCid));
        p(1) = patchline(temp_x,temp-4*PCid-1,'edgecolor',[0 0 0],'linewidth',4,'edgealpha',PCA_.var_explained(PCid)/100);
    end
%     axis([1000*([-1*Tcorr Tcorr]) -5 1])
    plot([0.8*max(temp_x) 0.8*max(temp_x)+5],[0.5 0.5],'color',[0.5 0.5 0.5],'linewidth',4);
    text(0.8*max(temp_x),2,'100ms')
    box off; set(gca,'YTick',[]); axis off
    title('Pricipal components of feature vector')

clear cluster_input p PCid s temp temp_x
%% Cluster on artificial templates: make prototypes
temp_x = (1:P.Ltr*2)*P.bw;
times =  {[];...
          [5];[10];[15];[20];[25];...
          [5 10];[10 15];[15 20];[20 25];[10 20]};
% times =  {[];...
%           [5];[10];[15];[20];[25];...
%           [10 20]};
bw = 0.05;
sigma    = 2;    % kernel standard deviation (s)
shoulder = 2;    % Time ranges to draw kernel over (sigmas)
tb = (1:2*P.Ltr)*bw;
templates = makeDecodingTemplates(tb,times,bw,sigma,shoulder,false);

clear clusters result param

param.c = length(times)
% Fudge cluster templates:
clusters.result.cluster.v = templates.result
clusters.data.X = Group.UnitFeatures{1}';%./max(max(Group.UnitFeatures{2}))
new =clust_normalize(clusters.data,'var');
% Evaluate data in prototype templates:
clusters.eval = clusteval(new,clusters.result,param);
for i = 1:size(clusters.data.X,1)
%     clusters.eval.f_(i,:) =  clusters.eval.f(i,:) == min(clusters.eval.f(i,:),[],2);
    clusters.eval.f_(i,:) =  clusters.eval.d(i,:) == min(clusters.eval.d(i,:),[],2);
end

offset = 2;

figure('name',' histogram templates','Color','w')
subplot(1,2,1); hold on
title('Prototype decoding templates')
for clID=1:param.c
    plot(temp_x,templates.result(clID,:)+offset*(clID),'k','LineWidth',2)
    axis off
end
subplot(1,2,2); hold on
title('Averaged clustered members')
for clID=1:param.c
    temp=find(clusters.eval.f_(:,clID));
    plot(temp_x,mat2gray(new.X(temp,:))      +offset*clID,'color',[0.8 0.8 0.8],'LineWidth',1);
    plot(temp_x,mat2gray(mean(new.X(temp,:)))+offset*clID,'color','k','LineWidth',2);
    axis off
end
%% Clustering - hierarchical
hierarchical_.distance_ = 'euclidean'
hierarchical_.method_   = 'ward'
s         = 3
iRat      = 3;
cluster_input = zscore(Group.UnitFeatures{s}',[],2); % All units from all rats in area s
% cluster_input = Group.UnitFeatures{s}'; % All units from all rats in area s
% for iClu=1:size(cluster_input,1)
%     cluster_input(iClu,:) = zscore(cluster_input(iClu,:));
% end
% cluster_input = zscore(Group.Units.AllTS{s}{iRat}',[],2);  %Single area, rat
% cluster_input = Group.Units.AllTS{s}{iRat}';




% figure; hold on;
% % plot(staggerplot(cluster_input',0,1),'k','LineWidth',1)
% imagesc(cluster_input)

hierarchical_.eucD = pdist(cluster_input,hierarchical_.distance_); 
hierarchical_.clustTreeEuc = linkage(hierarchical_.eucD,hierarchical_.method_);

hierarchical_.diagnostic.sqEucD=squareform(hierarchical_.eucD);

[hierarchical_.diagnostic.coph hierarchical_.diagnostic.sqCoph]=cophenet(hierarchical_.clustTreeEuc,hierarchical_.eucD);
hierarchical_.diagnostic.sqCoph=squareform(hierarchical_.diagnostic.sqCoph);

hierarchical_.clustTreeEuc(:,3)=hierarchical_.clustTreeEuc(:,3)./(max(hierarchical_.clustTreeEuc(:,3))); % link/link_max

% figure('name','Decoding shape clustering','Color','w'); hold on
[h,nodes] = dendrogram(hierarchical_.clustTreeEuc,0,'colorthreshold',0.5);
set(gca,'TickDir','out',...
    'TickLength',[.002 0],...
    'XTickLabel',[]);
ylabel('linkage / max. linkage','interpreter','none');
set(gcf,'name','Histogram shape clustering','Color','w'); 
view (-90,90)
axis off
set(gcf,'color','w')

% hierarchical plot
% T = clusterdata(clustTreeEuc,...
%                 'maxclust',noClusts,...
%                 'distance',distance_,...
%                 'linkage',method_); % agglomerative
hierarchical_.T = cluster(hierarchical_.clustTreeEuc,noClusts); % 



% figure; silhouette(cluster_input,hierarchical_.T,hierarchical_.distance_)

%%% remap data by clusters peak time (list)
cluster_input_=[];
% list = [4 3 2 5 1];
for clID=1:noClusts
    temp = find(hierarchical_.T==(clID));
    cluster_input_=[cluster_input_ ;cluster_input(temp,:)];
end
hierarchical_.ReorderedInput = cluster_input_;

figure('color','w')
for clID=1:noClusts
    clID
    subplot(ceil(noClusts^0.5),ceil(noClusts^0.5),clID); hold on
    temp = find(hierarchical_.T==(clID));
%     plot(temp_x,cluster_input(temp,:),'color',[0.5 0.5 0.5]);
    imagesc(cluster_input(temp,:))
    plot(5*(1+mean(cluster_input(temp,:))),'k','LineWidth',2);
    plot([0 0],[0 1],':k')
% 	axis([0 30 -1 10])
    axis off
end



% show pairwise similarity plot
m = mat2gray(squareform(pdist(cluster_input_,hierarchical_.distance_)));



% eucD = pdist(cluster_input_,hierarchical_.distance_); 
% clustTreeEuc_ = linkage(hierarchical_.eucD,hierarchical_.method_);
% [~, m]=cophenet(clustTreeEuc_,eucD_);
% m=squareform(m);



ml=triu(logical(ones(size(m))),1);
m(~ml)=nan;
load('BlueWhite.mat')

figure('color','w');
h=pcolor(m);
set(h,'edgealpha',0)
colormap(cmap)
axis off;
axis ij;
% box on


s=3
cluster_input = zscore(Group.UnitFeatures{s}',[],2); % All units from all rats in area s

clob=clustergram(cluster_input,...
    'RowPDist',hierarchical_.distance_,...
    'Cluster','column',...
    'Linkage',hierarchical_.method_,...
    'Dendrogram',1000,... %100
    'DisplayRatio',0.2,...
    'Colormap',parula,...
    'Symmetric',true,...
    'OptimalLeafOrder',true)

% group_markers = struct('GroupNumber', {369,363,287,393,  391,},...
%                       'Annotation', ClassNames,...
%                       'Color', ClassColors);
%                   
%                   clob.RowGroupMarker = group_markers;
                  
clear h nodes s cluster_input iRat temp ml m group_markers cmap cluster_input_ clID clob
%% k-means: how many clusters?
clear clusters result param
clusters.data.X=data';
clusters.data=clust_normalize(clusters.data,'var');

% PC=[];CE=[];SC=[];S=[];XB=[];DI=[];ADI=[];
nbs = 50;
ment=[];
ncmax=15;
for cln=2:ncmax
    param.c=cln;
    param.vis=0;
    for i =1:nbs
        [cln i]
        success = false;
        % Cluster: Keep respawning until it works
        while (~success)
            try
                clusters.result = Kmeans(clusters.data,param);
                success = true;
            catch %err
            end
        end 
        new.X=clusters.data.X;
        clusteval(new,clusters.result,param);
        %validation
        clusters.result=modvalidity(clusters.result,clusters.data,param);
        ment=clusters.result.validity;
        
        model.PC(cln,i) = ment.PC;
        model.CE(cln,i) = ment.CE;
        model.SC(cln,i) = ment.SC;
        model.S(cln,i) = ment.S;
        model.XB(cln,i) = ment.XB;
        model.DI(cln,i) = ment.DI;
        model.ADI(cln,i) = ment.ADI;

    end
end

    figure
    clf
    subplot(7,1,1); hold on
    ciplot(nanmean(model.PC,2)-nansem(model.PC,2),nanmean(model.PC,2)+nansem(model.PC,2),[NaN,2:ncmax],'b')
    plot([NaN,2:ncmax],nanmean(model.PC,2));
    title('Partition Coefficient (PC)')
    
    subplot(7,1,2); hold on
    ciplot(nanmean(model.CE,2)-nansem(model.CE,2),nanmean(model.CE,2)+nansem(model.CE,2),[NaN,2:ncmax],'b')
    plot([NaN,2:ncmax],nanmean(model.CE,2))
    title('Classification Entropy (CE)')
    
    subplot(7,1,3); hold on
    ciplot(nanmean(model.SC,2)-nansem(model.SC,2),nanmean(model.SC,2)+nansem(model.SC,2),[NaN,2:ncmax],'b')
    plot([NaN,2:ncmax],nanmean(model.SC,2))
    title('Partition Index (SC)')
    
    subplot(7,1,4); hold on
    ciplot(nanmean(model.S,2)-nansem(model.S,2),nanmean(model.S,2)+nansem(model.S,2),[NaN,2:ncmax],'b')
    plot([NaN,2:ncmax],nanmean(model.S,2))
    title('Separation Index (S)')
    
    subplot(7,1,5); hold on
    ciplot(nanmean(model.XB,2)-nansem(model.XB,2),nanmean(model.XB,2)+nansem(model.XB,2),[NaN,2:ncmax],'b')
    plot([NaN,2:ncmax],nanmean(model.XB,2))
    title('Xie and Beni Index (XB)')
    
    subplot(7,1,6); hold on
    ciplot(nanmean(model.DI,2)-nansem(model.DI,2),nanmean(model.DI,2)+nansem(model.DI,2),[NaN,2:ncmax],'b')
    plot([NaN,2:ncmax],nanmean(model.DI,2))
    title('Dunn Index (DI)')
    
    
    subplot(7,1,7); hold on
    ciplot(nanmean(model.ADI,2)-nansem(model.ADI,2),nanmean(model.ADI,2)+nansem(model.ADI,2),[NaN,2:ncmax],'b')
    plot([NaN,2:ncmax],nanmean(model.ADI,2))
    title('Alternative Dunn Index (ADI)')

% plot(clusters.result.cluster.v(:,1),clusters.result.cluster.v(:,2),'ro')

% result = validity(clusters.result,data,param);
% result.validity
%% k-means clustering (abonyi toolbox)
clear clusters result param
s = 5

% clusters.data.X=zscore(Group.UnitFeatures{s}',[],2)
clusters.data.X=Group.UnitFeatures{3}';
% clusters.data=clust_normalize(clusters.data,'var');

param.c   = 5; 
param.vis = false;

clusters.result = Kmeans(clusters.data,param);

offset = 5;
figure('name','Derived histogram templates','Color','w')
subplot(1,2,1); hold on
title('Derived histogram templates')
for clID=1:param.c
    plot(temp_x,clusters.result.cluster.v(clID,:)+offset*clID,'k','LineWidth',2)
    axis off
end
subplot(1,2,2); hold on
title('Averaged clustered data')

for clID=1:param.c
    temp=find(clusters.result.data.f(:,clID));
    plot(temp_x,clusters.data.X(temp,:)+offset*clID,'color',[0.8 0.8 0.8],'LineWidth',1);
    plot(temp_x,mean(clusters.data.X(temp,:))+offset*clID,'color','k','LineWidth',2);
    axis off
end
%% k-means clustering (MATLAB)
clear clusters kmeans_
kmeans_.distance_ = 'cosine'

noClusts  = 6;
s=3
cluster_input = zscore(Group.UnitFeatures{s}',[],2); % All units from all rats in area s
% cluster_input = zscore(Group.AssemFeatures{s}',[],2); % All units from all rats in area s
% 
% kmeans_.eval = evalclusters(cluster_input,'linkage','gap','klist',[1:10],'Distance',kmeans_.distance_)
% figure; plot(kmeans_.eval)
% fprintf ('Optimal no. clusters = %d\n',kmeans_.eval.OptimalK)
% noClusts=kmeans_.eval.OptimalK

kmeans_.opts = statset('Display','final');
kmeans_.distBins = 0:0.05:1;
[kmeans_.T,kmeans_.Centroids,kmeans_.sumD,kmeans_.D] = kmeans(cluster_input,noClusts,'Distance',kmeans_.distance_,'Replicates',100,'Options',kmeans_.opts);
% figure; plot(staggerplot(kmeans_.Centroids',0,.1))

% Reorder clusters by peak time
for clID = 1:noClusts
    %   a = peakfinder(zscore(kmeans_.Centroids(clID,:)))
    [~,a] = max(kmeans_.Centroids(clID,:));
    peakLoc(clID) =  a(1);
end
T = zeros(size(kmeans_.T));
D = zeros(size(kmeans_.D));
[~,idx] = sort(peakLoc);
figure; plot(staggerplot(kmeans_.Centroids',0,.1))

for clID = 1:noClusts
     T(kmeans_.T==idx(clID)) = clID;
     D(:,idx(clID)) = kmeans_.D(:,clID);
end

kmeans_.sumD = kmeans_.sumD(idx);
kmeans_.Centroids = kmeans_.Centroids(idx,:) ;
kmeans_.T = T;
kmeans_.D = D;


cluster_input_ = [];
for clID = 1:noClusts
    noUnits(clID) = sum(kmeans_.T==clID);
    cluster_input_=[cluster_input_; cluster_input(find(kmeans_.T==clID),:)];
end
figure; 
    subplot(1,2,1)
    plot(staggerplot(flipud(kmeans_.Centroids'),0,.1))
    subplot(1,2,2)
    imagesc(cluster_input_)
clear a peakLoc T D 

% distance histograms
for clID=1:noClusts
    temp=find(kmeans_.T==clID);
    kmeans_.Dhist(:,clID)=cumsum(histc(kmeans_.D(temp,clID),kmeans_.distBins));  
    kmeans_.Dhist(:,clID)=kmeans_.Dhist(:,clID)./max(kmeans_.Dhist(:,clID))
end

temp_x = (1:P.Ltr*2)*P.bw;
offset = 0.1;
figure('name','Cluster analysis','Color','w')
subplot(1,2,1); hold on
    title('Derived centroids')
    for clID=1:noClusts
        plot(temp_x,kmeans_.Centroids(clID,:)+offset*clID,'color',ClassColors{clID},'LineWidth',2)
        axis off
    end
    plot([10 10],[0 1],'g','LineWidth',2)
    plot([20 20],[0 1],'r','LineWidth',2)

subplot(1,2,2); hold on
    title('Clustered data')
    offset=20
    for clID=1:noClusts
        temp=find(kmeans_.T==clID);
        temp_ = Group.UnitFeatures{s}';
        plot(temp_x,temp_(temp,:)+offset*clID,'color',[0.8 0.8 0.8],'LineWidth',1);
        plot(temp_x,mean(temp_(temp,:))+offset*clID,'color','k','LineWidth',2);
        axis off
    end


% figure;hold on
%     title('distribution of distance from centroid')
%     plot(kmeans_.distBins, kmeans_.Dhist)
%     legend(string(1:noClusts))



cluster_input__=[];
temp=[];
for clID=1:noClusts
    temp = find(kmeans_.T==clID);
    cluster_input__=[cluster_input__ ;cluster_input(temp,:)];
end

m = squareform(pdist(cluster_input__,kmeans_.distance_));

% eucD = pdist(cluster_input_,distance_); 
% clustTreeEuc = linkage(eucD,method_);
% [~, m]=cophenet(clustTreeEuc,eucD);
% m=squareform(m);

ml=triu(logical(ones(size(m))),1);
m(~ml)=nan;

for i=1:size(cluster_input_,1)
    cluster_input__(i,:)=mat2gray(cluster_input_(i,:));
end

figure('color','w');
ax1 = subplot(1,2,1); hold on
x=repmat(temp_x,size(cluster_input_,1),1);
y=repmat(1:size(cluster_input_,1),size(cluster_input_,2),1)';
    imagesc(x(1:end),y(1:end),cluster_input_);
    plot([10 10],[0 size(cluster_input_,1)],'g','LineWidth',2)
    plot([20 20],[0 size(cluster_input_,1)],'r','LineWidth',2)
	colormap(ax1,gray)
	%caxis(ax1,[-5 3])
    for i=1:noClusts
        plot(temp_x,mat2gray(kmeans_.Centroids(i,:))*0.9*noUnits(i)+sum(noUnits(1:i-1)),'LineWidth',4,'color',ClassColors{i})
        plot(temp_x,mat2gray(kmeans_.Centroids(i,:))*0.9*noUnits(i)+sum(noUnits(1:i-1)),'LineWidth',2,'color','k')
        text(31,5+sum(noUnits(1:i-1)),ClassNames{i},'color',ClassColors{i},'HorizontalAlignment','Left','FontSize',12)
        % text(605,5+sum(noUnits(1:i-1)),ClassNames{i},'color','k','HorizontalAlignment','Left','FontSize',12)
    end
        plot([15 15 ],[0 size(cluster_input_,1)],'w','LineWidth',3)

    xlim(ax1,[0 30])
    ylim([0 size(cluster_input_,1)])
    % ylim([0 sum(noUnits)])    
    axis off
%     caxis([0 6])
% colorbar 

ax2 = subplot(1,2,2); hold on
    h=pcolor(rot90(m));
    set(h,'edgealpha',0)
    axis off;
    colormap(ax2,gray)
    noUnits_ = [1 cumsum(fliplr(noUnits))];
    noUnits__ = fliplr(cumsum([1 noUnits]));
    ClassColors_ = fliplr([{[0 0 0]},{[0 0 0]},ClassColors]);
    for i = 2:noClusts+1
        plot([noUnits__(i-1) noUnits__(i)],[noUnits_(i-1) noUnits_(i-1)],'color',(ClassColors_{i-1}),'LineWidth',2)
        plot([noUnits__(i-1) noUnits__(i)],[noUnits_(i-1) noUnits_(i)],'color',(ClassColors_{i-1}),'LineWidth',2)
        plot([noUnits__(i) noUnits__(i)],[noUnits_(i-1) noUnits_(i)],'color',(ClassColors_{i-1}),'LineWidth',2)

    end
    
    caxis([0 1])
    axis ij;
    axis tight
    plot([200 300],[350 350],'k','LineWidth',2)
    plot([200 200],[300 350],'k','LineWidth',2)
    text(220,320,{'50 Units';'/ 5s'},'HorizontalAlignment','Left','FontSize',12)


clear s clID iRat offset temp temp_ temp_x h m ml cluster_input cluster_input_ cluster_input__ ClassColors_
%% Collect unit class assignments

% Units
collapseList = kmeans_.T;
areaList = [];
noUnits_ = [];
Units=[]
Units_Sig=[];
a=[];
for iRat = 1:length(RatList)
    
    listout(iRat).name_    = RatList(iRat).name;
    for s = 1:3
        Units_Sig{s} = find(max(Group.Units.AllTSsig{s}{iRat}<0.05));
        Units{s} = 1:size(Group.Units.AllTS{s}{iRat},2);
        
        noUnits(s) = numel(Units{s});
        noUnits_Sig(s) =  numel(Units_Sig{s});
        b{s} = nan(noUnits(s),1);
    end
    % Cut out units from this experiment from collapsed list
    temp_T = collapseList(1:noUnits_Sig(3)); 
    collapseList(1:noUnits_Sig(3)) = [];
    
    % class assignments
    a{1} = temp_T(1:noUnits_Sig(1));
    a{2} = temp_T(noUnits_Sig(1)+1:noUnits_Sig(1)+noUnits_Sig(2));
    a{3} = temp_T;
    
    for s=1:3
        b{s}(find(ismember(Units{s},Units_Sig{s}))) = a{s};
        listout(iRat).UnitClass{s} =  b{s}
    end
%     listout(iRat).UnitClass{1} =  b
%     listout(iRat).UnitClass{2} = temp_T(noUnits_(1)+1:noUnits_(1)+noUnits_(2));
%     listout(iRat).UnitClass{3} = temp_T;

     areaList = [areaList;ones(noUnits_Sig(1),1);2*ones(noUnits_Sig(2),1)];
end
clear iRat temp_T collapseList s noUnits noUnits_Sig Units Units_Sig  a b 
%% Cast Assems into single unit patterns
clear a
for s=1:3
    cluster_input{s} = zscore(Group.AssemFeatures{s}',[],2); % All units from all rats in area s
    [~,kmeansAssems_.T{s} ,~] =  matchClust(cluster_input{s},kmeans_.Centroids);
    
    % Silence membership of insignificant  decoders
    kmeansAssems_.T{s}(setdiff(1:size(Group.Assem.AllTSCollapse{s},2),Group.Assem.SigDecoders{s})) = NaN;
end

% plot results
s = 2;
temp_x = (1:P.Ltr*2)*P.bw;
offset = 0.1;
figure('name','Cluster analysis','Color','w')
subplot(1,2,1); hold on
title('Unit-derived cluster centroids')
for clID=1:noClusts
    plot(temp_x,kmeans_.Centroids(clID,:)+offset*clID,'k','LineWidth',2)
    axis off
end
subplot(1,2,2); hold on
title('Template-matched assemblies')
offset=20;
for clID=1:noClusts
    temp=find(kmeansAssems_.T{s}==clID);
    temp_ = cluster_input{s};%Group.UnitFeatures{s}';
    plot(temp_x,temp_(temp,:)+offset*clID,'color',[0.8 0.8 0.8],'LineWidth',1);
    plot(temp_x,mean(temp_(temp,:))+offset*clID,'color','k','LineWidth',2);
    axis off
end

clear a
%% Collect Assem class assignments
areaListAssem = [];
for s= 1:3
    collapseList = kmeansAssems_.T{s};
    for iRat = 1:length(RatList)    
        if iRat<=length(Group.Assem.AllTS{s})
            noAssems_ = size(Group.Assem.AllTS{s}{iRat},2);
            disp([s iRat noAssems_]);
            listout(iRat).AssemClass{s} = collapseList(1:noAssems_);
            areaListAssem = [areaListAssem;[collapseList(1:noAssems_),s*ones(noAssems_,1)]];
            
            collapseList(1:noAssems_) = [];

        else
            listout(iRat).AssemClass{s} =[];
        end

    end
end
clear iRat temp_T collapseList s noUnits_
%% Class type pie chart
figure('color','w','NumberTitle','off','name',['Unit decoding type'])
for s = 1:2
    
 
    subplot(1,3,s); hold on
    temp = hist( kmeans_.T(areaList==s),1:noClusts);
    tempLabels = ClassNames(temp~=0);
    tempColors = cell2mat(ClassColors(temp~=0)');
    Hpie = pie(temp,ClassNames);
    scatter(0,0,10000,'w','filled')    
    hp = findobj(Hpie,'Type','patch');
    for iClass= 1:numel(hp)
        set(hp(iClass), 'FaceColor', tempColors(iClass,:),...
                        'EdgeColor', [0.6 0.6 0.6],...
                        'FaceAlpha', 1,...
                        'LineWidth', 1,...
                        'LineStyle','none');
    end
    title({[P.names{s},' Units'];''})
    axis square
    axis off
end


figure('color','w','NumberTitle','off','name',['Assem decoding type'])
for s = 1:3
    subplot(1,3,s); hold on
    areaListAssem(areaListAssem(:,2)==s,1)

    temp = hist(areaListAssem(areaListAssem(:,2)==s,1),1:noClusts);
    tempLabels = ClassNames(temp~=0);
    tempColors = cell2mat(ClassColors(temp~=0)');
    Hpie = pie(temp,ClassNames);
    scatter(0,0,10000,'w','filled')    
    hp = findobj(Hpie,'Type','patch');
    for iClass= 1:numel(hp)
        set(hp(iClass), 'FaceColor', tempColors(iClass,:),...
                        'EdgeColor', [0.6 0.6 0.6],...
                        'FaceAlpha', 1,...
                        'LineWidth', 1,...
                        'LineStyle','none');
    end
    title({P.names{s};,'Assemblies';''})
    axis square
    axis off
end
%% Plot final classified units and assems

temp_x = (1:P.Ltr*2)*P.bw;
offset = 1;


figure('name','Final Classification','color','w') ;
for s = 1:2
    subplot(2,3,s); hold on
    plot([10 10],[0 7],'g','LineWidth',2)
    plot([20 20],[0 7],'r','LineWidth',2)
    temp_y =cell(noClusts,1);
    for iClu = 1:noClusts
        temp_y{iClu} = [];
        for iRat = 1:length(RatList)
            try
            temp_y{iClu} = [temp_y{iClu}, Group.Units.AllTS{s}{iRat}(:,listout(iRat).UnitClass{s}==iClu)]
            end
        end
        temp_y{iClu} = mat2gray(temp_y{iClu})
        plot(temp_x, offset*(iClu-1)+temp_y{iClu},'color',0.6*[1 1 1])
        plot(temp_x, offset*(iClu-1)+mean(temp_y{iClu},2),'color',ClassColors{iClu},'LineWidth',3)
    end
    axis off
    plot([15 15],[0 7],'w','LineWidth',3)
    plot([15 15],[0 7],':k') 
end
for s = 1:3
    subplot(2,3,s+3); hold on
    plot([10 10],[0 7],'g','LineWidth',2)
    plot([20 20],[0 7],'r','LineWidth',2)
    temp_y =cell(noClusts,1);
    for iClu = 1:noClusts
        temp_y{iClu} = [];
        for iRat = 1:length(RatList)
            try
                temp_y{iClu} = [temp_y{iClu}, Group.Assem.AllTS{s}{iRat}(:,listout(iRat).AssemClass{s}==iClu)]
            end
        end
        try
            temp_y{iClu} = mat2gray(temp_y{iClu})        
            plot(temp_x, offset*(iClu-1)+temp_y{iClu},'color',0.6*[1 1 1])
            plot(temp_x, offset*(iClu-1)+mean(temp_y{iClu},2),'color',ClassColors{iClu},'LineWidth',3)
        end

    end
    axis off
	plot([15 15],[0 7],'w','LineWidth',3)   
    plot([15 15],[0 7],':k')


end
% %%
% for s = 1:3
%     temp_y =cell(noClusts,1);
%     for iClu = 1:noClusts
%         temp_y{iClu} = [];
%         for iRat = 1:length(RatList)
%             try
%                 tempa = find(listout(iRat).AssemClass{s}==iClu)
%                 for iA=1:length(tempa)
%                     Assem_ = tempa(iA)
%                     
%             temp_y{iClu} = [temp_y{iClu}, Group.Units.AllTS{s}{iRat}(:,listout(iRat).UnitClass{s}==iClu)]
%             end
%         end
%         temp_y{iClu} = mat2gray(temp_y{iClu})
%         plot(temp_x, offset*(iClu-1)+temp_y{iClu},'color',0.6*[1 1 1])
%         plot(temp_x, offset*(iClu-1)+mean(temp_y{iClu},2),'color',ClassColors{iClu},'LineWidth',3)
%     end
% end
%% save classification results

save([pat{1} 'UnitAssemClassified.mat'],'listout','noClusts','kmeans_','kmeansAssems_','hierarchical_','PCA_','templates')

%% Preamble, load
clear all 
close all
P.target = 'Task';
P.rat    = 'JaroslawLONG1';  %    Post sleep import error
% P.rat    = 'JaroslawLONG2';  %*
% P.rat    = 'KrzesimirLONG1'; %*
% P.rat    = 'KrzesimirLONG2'; %*
% P.rat    = 'KrzysztofLONG1'; %*
% P.rat    = 'KrzysztofLONG2'; %*
% P.rat    = 'MiroslawLONG0';  %*
% P.rat    = 'MiroslawLONG2';  %*   No Post sleep recorded
% P.rat    = 'NorbertLONG1';   %*
% P.rat    = 'NorbertLONG2';   %
% P.rat    = 'OnufryLONG1';    %*
% P.rat    = 'OnufryLONG2';    %      Post sleep error

P.expt   = 'LONG';
P.flag.PlotVerbose = true;
if ispc
%     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\';       % location of the processed sleep assembly activations
%     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - epoch specific FRs\';
    pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - use task mean FRs\';
    pat{2} = 'C:\Analysis\AssemblyAnalysis\raw';           % location of raw spike times
else
    pat{1} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Decoding vs sleep/';
    pat{2} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Task';
end
cd(pat{1})

set(0,'DefaultFigureWindowStyle','docked')

% Load
load([pat{1} P.rat '_' P.target '_DecodingVSsleep.mat'])
%% Specify time ranges 
P.tRangesScheme = 5;

switch P.tRangesScheme
    case 5
    P.tRangesNames= {'Pre-Sample'; 'Sample'; 'Delay'; 'Choice'; 'Post-Sample'};
    % P.tRanges = [2 8;...
    %              9 11;...
    %              13 17;...
    %              19 21;...
    %              24 30];

    P.tRanges = [0     8.99;...
                 9     11;...
                 11.01 18.99;...
                 19    21;...
                 21.01 30];

    P.tRangesColors= [0.7 0.9 0.9;...
                      0.1 0.9 0.1;...
                      0.9 0.7 0.9;...
                      0.9 0.1 0.1;...
                      0.9 0.9 0.7];       
    case 6
    P.tRangesNames= {'Pre-Sample'; 'Sample Press'; 'Early Delay'; 'Late Delay'; 'Choice Press'; 'Post-Sample'};
    % P.tRanges = [2 8;...
    %              9 11;...
    %              13 17;...
    %              19 21;...
    %              24 30];

    P.tRanges = [0     8.99;...
                 9     11;...
                 11.01 14.99;...
                 15    18.99;...
                 19    21;...
                 21.01 30];

    P.tRangesColors= [0.7 0.9 0.9;...
                      0.1 0.9 0.1;...
                      0.9 0.75 0.9;...
                      0.9 0.55 0.9;...
                      0.9 0.1 0.1;...
                      0.9 0.9 0.7];            
    case 2
            
P.tRangesNames= {'Sample Press'; 'Choice Press'};

P.tRanges = [9  11;...
             19 21];
P.tRangesColors= [0.7 0.9 0.7;...
                  0.9 0.7 0.7];        
end

P.tRangesAdjusted= P.tRanges/P.bw;
P.tRangesCentre=P.tRanges(:,1) + (P.tRanges(:,2)-P.tRanges(:,1))/2;
%% Plot decoding results - both units and assemblies
figure('name',P.rat,'color','w')
% Plot decoding
maxH = 10;

for s=1:3
    subplot(2,3,s); hold on
    title(P.names{s})

%     plot([5 5],[0 10],'color',[0.89 0.95 0.89],'LineWidth',20)
%     plot([25 25],[0 10],'color',[0.95 0.89 0.89],'LineWidth',20) 
    plot([10 10],[0 10*maxH],'color',[0.89 0.95 0.89],'LineWidth',20)
    plot([20 20],[0 10*maxH],'color',[0.95 0.89 0.89],'LineWidth',20)
    
%     Y1=[D.Assem.Ft2ciL0{s};D.Assem.Ft2ciH0{s}-D.Assem.Ft2ciL0{s}];
%     h=area((1:P.Ltr*2)*P.bw,Y1');
%     set(h(1),'FaceColor','w'); set(h(2),'FaceColor','w'); hold on
    plot((1:P.Ltr*2)*P.bw,D.Assem.Ft2ciH0{s},'k')
    plot((1:P.Ltr*2)*P.bw,D.units.Ft2ciH0{s},'b')

	plot((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s},'k','LineWidth',1.5)
    plot((1:P.Ltr*2)*P.bw,D.units.Ft2x{s},'b','LineWidth',1.5)
%   area((1:P.Ltr*2)*P.bw,D.units.Ft2x{s},'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
% 	area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s},'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)

% area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s}-D.units.Ft2x{s},'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    plot([15 15],[0 10*maxH],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[0 10*maxH],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    set(gca,'XTick',[])
    
    if s==1
        ylabel({'Decoding score';'(F-value)'})
        axis([0 Inf 0 maxH])

    elseif s==2
%         xlabel('Time (s)')
        text(15,-1,'Single unit decoding','Color',[0 0 1],'HorizontalAlignment','center' )
        text(15,-3,'Assembly decoding','Color',[0 0 0],'HorizontalAlignment','center' )
        axis([0 Inf 0 2*maxH])
    elseif s==3      
        axis([0 Inf 0 2*maxH])
    end
end
% Plot assembly advantage
for s=1:3
    subplot(2,3,s+3); hold on
%     title(P.names{s})

    plot([9 9],[-10 10],'color',[0.89 0.95 0.89],'LineWidth',20)
    plot([20 20],[-10 10],'color',[0.95 0.89 0.89],'LineWidth',20)
    
    temp = D.Assem.Ft2x{s}-D.units.Ft2x{s}; temp(temp>0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    temp = D.Assem.Ft2x{s}-D.units.Ft2x{s}; temp(temp<0)= 0;
    area((1:P.Ltr*2)*P.bw,temp,'FaceColor',[0.6 0.6 0.6],'EdgeColor','k','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
    % area((1:P.Ltr*2)*P.bw,D.Assem.Ft2x{s}-D.units.Ft2x{s},'FaceColor',[0.6 0.6 0.9],'EdgeColor','b','LineWidth',1.5,'FaceAlpha',0.6,'EdgeAlpha',1)
    
%     plot([0 P.Ltr*2*P.bw],[0 0],'color','k','LineWidth',2)
    plot([15 15],[-10 10],'color',[1 1 1],'LineWidth',5, 'LineStyle','-')
    plot([15 15],[-10 10],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    axis([0 Inf -4 4])
    if s==1
%         ylabel({'Assembly decodimg advantage';'Assembly - Unit'})
        ylabel({'Assembly decodimg advantage'})
    elseif s==2
        xlabel('Time (s)')
    end
end
%% Sort units/assems profiles by peak time / amp

% sort unit firing rates
for s=1:3
    % Strength/Time of mean trial-average firing rate increases for each outcome
    [FRtrials.iFRpeak{s}(1,:), FRtrials.iFRpeakLoc{s}(1,:)] = max(FRtrials.iFR0_Mean{s}{1}); 
    [FRtrials.iFRpeak{s}(2,:), FRtrials.iFRpeakLoc{s}(2,:)] = max(FRtrials.iFR0_Mean{s}{2});
    
     % Strength of mean trial-average firing rate increases regardless of outcome
     FRtrials.iFRpeakGlobal{s} = max(FRtrials.iFRpeak{s});
     
     % L or R preference?
     FRtrials.iFRpref{s}       = (FRtrials.iFRpeak{s}(2,:)>=FRtrials.iFRpeak{s}(1,:))+1;
     
     % Position of mean trial-average firing rate increases regardless of outcome
     for uID = 1:length(FRtrials.iFRpref{s})
         FRtrials.iFRlocGlobal{s}(uID)    =  FRtrials.iFRpeakLoc{s}(FRtrials.iFRpref{s}(uID),uID) ;
     end
     
     % 1) Rank units by peak firing time
    [~,FRtrials.iFRPeaktimeSorted{s}] = sort(FRtrials.iFRlocGlobal{s});
    [~,FRtrials.iFRPeakRateSorted{s}] = sort(FRtrials.iFRpeakGlobal{s});
    
    % 2) rank units by time of peak decoding
    [FRtrials.DPeak{s}(1,:), FRtrials.DPeakLoc{s}(1,:)] = max(D.units.TS{s});
    
    [~,FRtrials.DPeakTimeSorted{s}] = sort(FRtrials.DPeakLoc{s});
    [~,FRtrials.DPeakRateSorted{s}] = sort(FRtrials.DPeak{s});

end

% Sort assembly activation
for s=1:3
    if ~isempty(FAtrials.FSC_Mean{s})
    [FAtrials.AssemPeak{s}(1,:), FAtrials.AssemPeakLoc{s}(1,:)] = max(FAtrials.FSC_Mean{s}{1});
    [FAtrials.AssemPeak{s}(2,:), FAtrials.AssemPeakLoc{s}(2,:)] = max(FAtrials.FSC_Mean{s}{2});
    
     FAtrials.AssemPeakGlobal{s} = max(FAtrials.AssemPeak{s});
     FAtrials.Assempref{s}       = (FAtrials.AssemPeak{s}(2,:)>=FAtrials.AssemPeak{s}(1,:))+1;
     
     for AssemID = 1:length(FAtrials.Assempref{s})
         FAtrials.AssemlocGlobal{s}(AssemID)    =  FAtrials.AssemPeakLoc{s}(FAtrials.Assempref{s}(AssemID),AssemID) ;
     end
     
    [~,FAtrials.AssemPeakTimeSorted{s}] = sort(FAtrials.AssemlocGlobal{s});
    [~,FAtrials.AssemPeakRateSorted{s}] = sort(FAtrials.AssemPeakGlobal{s});

    % Sort by decoding power 
    [FAtrials.DPeak{s}(1,:), FAtrials.DPeakLoc{s}(1,:)] = max(D.Assem.TS{s});
    
    [~,FAtrials.DPeakTimeSorted{s}] = sort(FAtrials.DPeakLoc{s});
    [~,FAtrials.DPeakRateSorted{s}] = sort(FAtrials.DPeak{s});
    
    

    end
end
%% Classify units/assems by peak time of activation or by decoding shape
ClassifyMethod = 'Shape'; % {'Shape' or 'Time'}; 
temp_x = (1:2*P.Ltr)*P.bw;
templates = makeDecodingTemplates(temp_x,num2cell(round(P.tRangesCentre)),P.bw,2,2,false);
param.c = length(P.tRangesCentre);        
plotOnline = true;        
ClassTemp = 1:size(P.tRanges,1);
for s=1:3
    % NB classification scheme defaults to 0 if no label is found
    FRtrials.unitClass{s} = zeros(size(FRtrials.DPeakLoc{s}));
    FAtrials.AssemClass{s} = zeros(size(FAtrials.DPeakLoc{s}));
    if      strcmp(ClassifyMethod,'Time')
        for classID = 1:size(P.tRanges,1)
            % Maximum Likelihood classification of units by peak decoding time
            FRtrials.unitClass{s}(FRtrials.DPeakLoc{s}>P.tRangesAdjusted(classID,1) & FRtrials.DPeakLoc{s}<=P.tRangesAdjusted(classID,2)) ...
                = classID;
            % Maximum Likelihood classification of factors by peak decoding time
            FAtrials.AssemClass{s}(FAtrials.DPeakLoc{s}>P.tRangesAdjusted(classID,1) & FAtrials.DPeakLoc{s}<=P.tRangesAdjusted(classID,2)) ...
                = classID;
        end
    elseif  strcmp(ClassifyMethod,'Shape')
        %%%%%%%% Maximum Likelihood classification  by shape of decoding profile
        %%%% (1) units
        clusters = struct;
        clusters.result.cluster.v = templates.result;
        clusters.data.X = D.units.TS{s}';
        clusters.norm    = clust_normalize(clusters.data,'var');
        % Evaluate data in prototype templates:
        clusters.eval = clusteval(clusters.norm,clusters.result,param);
        for i = 1:size(clusters.data.X,1)
            clusters.eval.f_(i,:) =  clusters.eval.d(i,:) == min(clusters.eval.d(i,:),[],2);
            FRtrials.unitClass{s}(:,i) = find(clusters.eval.f_(i,:))
        end
        if plotOnline 
            
            offset = 2;
            figure('name',' histogram templates','Color','w')
            subplot(1,2,1); hold on
            title('Prototype decoding templates')
                boxHeight =offset*(param.c+1);
                for iTrange = 1:size(P.tRanges,1)
                    area(P.tRanges(iTrange,:),[boxHeight boxHeight],...
                      'FaceColor',P.tRangesColors(iTrange,:),...
                      'EdgeColor',P.tRangesColors(iTrange,:),...
                      'LineWidth',1.5,...
                      'FaceAlpha',0.5,...
                      'EdgeAlpha',0)
                end
                plot([15 15],[1 boxHeight],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
                for iTrange=1:param.c
                    plot(temp_x,templates.result(iTrange,:)+offset*(iTrange),'k','LineWidth',2)
                end
                axis([0 30 2 boxHeight]) 
                axis off
                
            subplot(1,2,2); hold on
                title('Averaged clustered members')
                for iTrange = 1:size(P.tRanges,1)
                    area(P.tRanges(iTrange,:),[boxHeight boxHeight],...
                      'FaceColor',P.tRangesColors(iTrange,:),...
                      'EdgeColor',P.tRangesColors(iTrange,:),...
                      'LineWidth',1.5,...
                      'FaceAlpha',0.5,...
                      'EdgeAlpha',0)
                end
                for iTrange=1:param.c
                    temp=find(clusters.eval.f_(:,iTrange));
                    plot(temp_x,mat2gray(clusters.norm.X(temp,:))      +offset*iTrange,'color',[0.8 0.8 0.8],'LineWidth',1);
                    plot(temp_x,mat2gray(mean(clusters.norm.X(temp,:),1))+offset*iTrange,'color','k','LineWidth',2);
                end
                axis([0 30 2 boxHeight]) 
                axis off

        end

        %%%% (2) Assemblies
        clusters = struct;
        clusters.result.cluster.v = templates.result;
        clusters.data.X = D.Assem.TS{s}';
        clusters.norm    = clust_normalize(clusters.data,'var');
        % Evaluate data in prototype templates:
        clusters.eval = clusteval(clusters.norm,clusters.result,param);
        for i = 1:size(clusters.data.X,1)
            clusters.eval.f_(i,:) =  clusters.eval.d(i,:) == min(clusters.eval.d(i,:),[],2);
            FAtrials.AssemClass{s}(:,i) = find(clusters.eval.f_(i,:));
        end
        if plotOnline 
            
            offset = 2;
            figure('name',' histogram templates','Color','w')
            subplot(1,2,1); hold on
            title('Prototype decoding templates')
                boxHeight =offset*(param.c+1);
                for iTrange = 1:size(P.tRanges,1)
                    area(P.tRanges(iTrange,:),[boxHeight boxHeight],...
                      'FaceColor',P.tRangesColors(iTrange,:),...
                      'EdgeColor',P.tRangesColors(iTrange,:),...
                      'LineWidth',1.5,...
                      'FaceAlpha',0.5,...
                      'EdgeAlpha',0)
                end
                plot([15 15],[1 boxHeight],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
                for iTrange=1:param.c
                    plot(temp_x,templates.result(iTrange,:)+offset*(iTrange),'k','LineWidth',2)
                end
                axis([0 30 2 boxHeight]) 
                axis off
                
            subplot(1,2,2); hold on
                title('Averaged clustered members')
                for iTrange = 1:size(P.tRanges,1)
                    area(P.tRanges(iTrange,:),[boxHeight boxHeight],...
                      'FaceColor',P.tRangesColors(iTrange,:),...
                      'EdgeColor',P.tRangesColors(iTrange,:),...
                      'LineWidth',1.5,...
                      'FaceAlpha',0.5,...
                      'EdgeAlpha',0)
                end
                for iTrange=1:param.c
                    temp=find(clusters.eval.f_(:,iTrange));
                    plot(temp_x,mat2gray(clusters.norm.X(temp,:))      +offset*iTrange,'color',[0.8 0.8 0.8],'LineWidth',1);
                    plot(temp_x,mat2gray(mean(clusters.norm.X(temp,:),1))+offset*iTrange,'color','k','LineWidth',2);
                end
                axis([0 30 2 boxHeight]) 
                axis off

        end
        %%%%%%%%
    end       
end
clear uID AssemID
%% Assembly vs unit overlap in LR pref and peak time 
for s=1:3
    if ~isempty(FAtrials.FSC_Mean{s})
        for iAss = 1:length(FAtrials.keepers{s})
            
            % Concordance of Assembly and constituent units' L/R preference
            temp = (FRtrials.iFRpref{s}  == FAtrials.Assempref{s}(iAss));
            temp = temp(intersect(FAcont.Task.units{s}{iAss},find(FRtrials.DPeak{s}>0)));    % Restrict by member units (and optionally decoding score)
            if FAtrials.DPeak{s}(iAss)>0
                FAtrials.AssemUnitConcordLR{s}(iAss) = sum(temp)./numel(temp);
            else
                FAtrials.AssemUnitConcordLR{s}(iAss) = NaN;
            end
            
            % Concordance of Assembly and constituent units' activation classification
            temp = (FRtrials.unitClass{s} == FAtrials.AssemClass{s}(iAss));
            temp = temp(intersect(FAcont.Task.units{s}{iAss},find(FRtrials.DPeak{s}>0)));    % Restrict by member units (and optionally decoding score)
            if FAtrials.DPeak{s}(iAss)>0
                FAtrials.AssemUnitConcordType{s}(iAss) = sum(temp)./numel(temp);
            else
                FAtrials.AssemUnitConcordType{s}(iAss) = NaN;
            end
            
            
            % Fraction of each Assembly's member units in each class
            temp = FRtrials.unitClass{s};
            temp = temp(intersect(FAcont.Task.units{s}{iAss},find(FRtrials.DPeak{s}>0)));    % Restrict by member units (and optionally decoding score)
            FRtrials.unitClasshist{s}(iAss,:) = histc(temp,1:size(P.tRanges,1))./numel(temp);
            FRtrials.unitClasshistRaw{s}(iAss,:) = histc(temp,1:size(P.tRanges,1));
        end
        
        %Collapse 
        FRtrials.unitClasshistAssemType{s} = zeros(numel(P.tRangesNames),numel(P.tRangesNames));
        for iClass= 1:numel(P.tRangesNames)
            % Unit classification for all assemblies of this type
            temp = nansum(FRtrials.unitClasshistRaw{s}(FAtrials.AssemClass{s}==iClass,:));
            if ~isempty(temp)
                if numel(temp)>1
                    FRtrials.unitClasshistAssemType{s}(iClass,:) = temp./nansum(temp); 
                else
                    FRtrials.unitClasshistAssemType{s}(iClass,temp)=1;
                end
            else
                FRtrials.unitClasshistAssemType{s}(iClass,:) = zeros(size(P.tRangesNames,1),1);
            end
        end
    end
end
%% Plot sorted by decoding / FR
if P.flag.PlotVerbose
    figure('name','Firing rate profiles sorted by Peak rate','NumberTitle','off'); hold on
    for s=1:3
        
        subplot(1,3,s); hold on
        plot(staggerplot((FRtrials.iFR0_Mean{s}{1}(:,FRtrials.iFRPeakRateSorted{s})),0,10),'g');
        plot(staggerplot((FRtrials.iFR0_Mean{s}{2}(:,FRtrials.iFRPeakRateSorted{s})),0,10),'r');
    end

    figure('name','Assembly profiles sorted by Peak rate','NumberTitle','off'); hold on
    for s=1:3
        if ~isempty(FAtrials.FSC_Mean{s})
        subplot(1,3,s); hold on
        plot(staggerplot((FAtrials.FSC_Mean{s}{1}(:,FAtrials.AssemPeakRateSorted{s})),0,5),'g');
        plot(staggerplot((FAtrials.FSC_Mean{s}{2}(:,FAtrials.AssemPeakRateSorted{s})),0,5),'r');
        end
    end
end
%% plot activation scores sorted by peak time

figure('color','w','name','Assemblies sorded by activation time','NumberTitle','off' ); hold on
for s=1:3
    subplot(2,3,s); 
    hold on
%     plot(staggerplot((D.units.TS{s}(:,FRtrials.DPeakTimeSorted{s})),0,3),'k');
    x_= repmat((1:P.Ltr*2)*P.bw,size(D.units.TS{s},2),1);
    y_ = repmat(1:size(D.units.TS{s},2),P.Ltr*2,1)';
    z_= zscore(D.units.TS{s}(:,FRtrials.DPeakTimeSorted{s}))';
    imagesc(x_(1:end),y_(1:end),z_)
    colormap(flipud(gray))
    caxis([0 3])
    boxHeight =size(D.units.TS{s},2);
    for iTrange = 1:size(P.tRanges,1)
          area(P.tRanges(iTrange,:),[boxHeight boxHeight],...
                  'FaceColor',P.tRangesColors(iTrange,:),...
                  'EdgeColor',P.tRangesColors(iTrange,:),...
                  'LineWidth',1.5,...
                  'FaceAlpha',0.5,...
                  'EdgeAlpha',0)
    end
    plot([15 15],[1 boxHeight],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    axis([0 30 1 boxHeight]) 
    if s==1
        ylabel('Unit no.')
    elseif s==2
        title('Decoding Strength (z-scored)')
    end
end

for s=1:3
    if ~isempty(FAtrials.FSC_Mean{s})
    subplot(2,3,s+3); 
    hold on
%     plot(staggerplot((D.Assem.TS{s}(:,FAtrials.DPeakTimeSorted{s})),0,3),'k');
    x_= repmat((1:P.Ltr*2)*P.bw,size(D.Assem.TS{s},2),1);
    y_ = repmat(1:size(D.Assem.TS{s},2),P.Ltr*2,1)';
        
    z_= zscore(D.Assem.TS{s}(:,FAtrials.DPeakTimeSorted{s}))';
    imagesc(x_(1:end),y_(1:end),z_)
    % colormap((gray))
    caxis([0 3])
    boxHeight =size(D.Assem.TS{s},2);

    for iTrange = 1:size(P.tRanges,1)
          area(P.tRanges(iTrange,:),[2*boxHeight 2*boxHeight],...
                  'FaceColor',P.tRangesColors(iTrange,:),...
                  'EdgeColor',P.tRangesColors(iTrange,:),...
                  'LineWidth',1.5,...
                  'FaceAlpha',0.5,...
                  'EdgeAlpha',0.)
    end
    plot([15 15],[1 boxHeight],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
    axis([0 30 1 boxHeight])
    end
    title(P.names{s})
    if s==1
        ylabel('Assem. no.')
    elseif s==2
        xlabel('Time (s)')
    end
end
clear x_ y_ z_
%% Plot each assembly with member units
for s=1:3
    if ~isempty(FAcont.Task.keepers{s})
    figure('color','w','name',[P.names{s} ' Assemblies and member units'],'NumberTitle','off'); hold on
    for iAss =1:size(D.Assem.TS{s},2)
        subaxis(ceil(size(D.Assem.TS{s},2)^0.5),...
                ceil(size(D.Assem.TS{s},2)^0.5),iAss,'SpacingVert',0.05,'Margin',0.01)
        hold on
        % Plot shaded background
        % boxHeight =size(D.units.TS{s},2);
        boxHeight =length(FAcont.Task.units{s}{iAss})+2;
        for iTrange = 1:size(P.tRanges,1)
              area(P.tRanges(iTrange,:),[boxHeight boxHeight],...
                      'FaceColor',P.tRangesColors(iTrange,:),...
                      'EdgeColor',P.tRangesColors(iTrange,:),...
                      'LineWidth',1.5,...
                      'FaceAlpha',0.2,...
                      'EdgeAlpha',0.1)
        end
        plot([15 15],[0 boxHeight],'color',[0.6 0.6 0.6],'LineWidth',1.5, 'LineStyle',':')
                
        % Plot member units
        for iUnit = 1:length(FAcont.Task.units{s}{iAss})
            thisUnit = FAcont.Task.units{s}{iAss}(iUnit);
            plot((1:P.Ltr*2)*P.bw,...
                 iUnit + mat2gray( D.units.TS{s}(:,thisUnit)),...
                 'color', P.tRangesColors(FRtrials.unitClass{s}(thisUnit),:),'LineWidth',1)
        end

        % Plot Assembly
        plot((1:P.Ltr*2)*P.bw,...
             mat2gray(D.Assem.TS{s}(:,iAss)),'color', P.tRangesColors(FAtrials.AssemClass{s}(iAss),:),'LineWidth',3)
        plot((1:P.Ltr*2)*P.bw,mat2gray(D.Assem.TS{s}(:,iAss)),'k','LineWidth',1)
        axis([1 P.Ltr*2*P.bw 0 boxHeight])
        end
    end
end
%% Plot unit concordance
for s = 1:3
    if ~isempty(FAcont.Task.keepers{s})
    figure('color','w','name',[P.names{s},' Assembly and unit classification agreement'],'NumberTitle','off');
    subplot(1,2,1)
        temp = sort(FAtrials.AssemUnitConcordLR{s});
        b=bar([temp;1-temp]'*100,'stacked');
        b(1).FaceColor = [0.5 0.5 0.5];
        b(1).LineWidth = 1.5;
        b(2).FaceColor = 'w';
        b(2).LineWidth = 1.5;
        box off; set(gca,'XTick',[]); axis tight
        xlabel('Assembly no.')
        ylabel('Units/Assembly shared pref. (%)')
        title('L/R Preference')
    subplot(1,2,2)
        temp = sort(FAtrials.AssemUnitConcordType{s});
        b=bar([temp;1-temp]'*100,'stacked');
        b(1).FaceColor = [0.5 0.5 0.5];
        b(1).LineWidth = 1.5;
        b(2).FaceColor = 'w';
        b(2).LineWidth = 1.5;

        box off; set(gca,'XTick',[]); axis tight
        xlabel('Assembly no.')
        ylabel('Units/Assembly shared pref. (%)')
        title('Activation time Preference')
    end
end
%% Plot Assem/unit type pie charts
for s=1:3
    if ~isempty(FAcont.Task.keepers{s})
        figure('color','w','NumberTitle','off','name',[P.names{s}, 'Assemblies'])
        subplot(1,3,1);
        temp = hist(FAtrials.AssemClass{s},1:numel(P.tRangesNames));
        tempLabels = P.tRangesNames(temp~=0);
        tempColors = P.tRangesColors(temp~=0,:);
        Hpie = pie(temp,P.tRangesNames);
        
        hp = findobj(Hpie,'Type','patch');
        for iClass= 1:numel(hp)
            set(hp(iClass), 'FaceColor', tempColors(iClass,:),...
                            'EdgeColor', [0.6 0.6 0.6],...
                            'FaceAlpha', 0.4,...
                            'LineWidth', 1.5);
        end
        title('Assembly classification distribution')
        subplot(1,2,2)
        % b= barh(FRtrials.unitClasshist{s},'stacked');
        b= barh(FRtrials.unitClasshistAssemType{s},'stacked','ShowBaseLine','off');
        for iClass= 1:numel(P.tRangesNames)
            b(iClass).FaceColor = P.tRangesColors(iClass,:);
            b(iClass).LineWidth = 1.5;
            b(iClass).FaceAlpha = 0.4;
            b(iClass).EdgeColor = [0.6 0.6 0.6];
            text(-0.05,iClass,P.tRangesNames{iClass},'HorizontalAlignment','right','Color',P.tRangesColors(iClass,:))
        end
%         axis off
        text(-0.3,ceil(numel(P.tRangesNames)/2),'Assembly type','Rotation',90,'FontSize',14,'FontWeight','bold','HorizontalAlignment','center')
        title('Member unit classification distribution')
    end
end

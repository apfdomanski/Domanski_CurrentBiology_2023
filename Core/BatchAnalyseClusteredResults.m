% Preamble

clear all; close all

% ClassNames = {'Pre-Sample','Sample/Choice','Early Delay','Late Delay','Post-Choice'};
% ClassColors = {[0.7 0.9 0.9],[0.9 0.1 0.1],[0.9 0.7 0.9],[0.9 0.5 0.9],[0.9 0.9 0.7]};
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
RatList = dir([pat{1}  '*_DecodingVSsleep.mat']);

% Load classifier results
load([pat{1} 'UnitAssemClassified.mat'])
%% Process class assignments: Loop over rats

for iRat= 1:length(RatList)
    %% Get data
    
    disp([sprintf('Working on rat %d of %d',iRat, length(RatList)), ' (' RatList(iRat).name ')'])
    load([RatList(iRat).folder filesep RatList(iRat).name])
    
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
         % 3) sort units by both peak time and decoding power
         [~,FRtrials.DPeakTimeSorted{s}] = sort(FRtrials.DPeakLoc{s});
         [~,FRtrials.DPeakStrengthSorted{s}] = sort(FRtrials.DPeak{s});

    end

    % Sort assembly activation
    for s=1:3
        if ~isempty(FAtrials.FSC_Mean{s})
        % Strength/Time of mean trial-average Assem activation increases for each outcome
        [FAtrials.AssemPeak{s}(1,:), FAtrials.AssemPeakLoc{s}(1,:)] = max(FAtrials.FSC_Mean{s}{1});
        [FAtrials.AssemPeak{s}(2,:), FAtrials.AssemPeakLoc{s}(2,:)] = max(FAtrials.FSC_Mean{s}{2});
        
        % Strength of mean trial-average Assem activation increases regardless of outcome
        FAtrials.AssemPeakGlobal{s} = max(FAtrials.AssemPeak{s});
        FAtrials.Assempref{s}       = (FAtrials.AssemPeak{s}(2,:)>=FAtrials.AssemPeak{s}(1,:))+1;

        % Position of mean trial-average Assem activation increases regardless of outcome
        for AssemID = 1:length(FAtrials.Assempref{s})
            FAtrials.AssemlocGlobal{s}(AssemID)    =  FAtrials.AssemPeakLoc{s}(FAtrials.Assempref{s}(AssemID),AssemID) ;
        end
        
        % 1) rank Assems by time of peak decoding
        [~,FAtrials.AssemPeakTimeSorted{s}] = sort(FAtrials.AssemlocGlobal{s});
        [~,FAtrials.AssemPeakRateSorted{s}] = sort(FAtrials.AssemPeakGlobal{s});

        % 2) Sort by decoding power 
        [FAtrials.DPeak{s}(1,:), FAtrials.DPeakLoc{s}(1,:)] = max(D.Assem.TS{s});
        % 3) sort units by both peak time and decoding power
        [~,FAtrials.DPeakTimeSorted{s}] = sort(FAtrials.DPeakLoc{s});
        [~,FAtrials.DPeakStrengthSorted{s}] = sort(FAtrials.DPeak{s});
        end
    end

    %% Assembly vs unit overlap in LR pref and activation pattern
    
    for s=1:3
        
        % get classifier results
        idx = find(~cellfun(@isempty,(strfind({listout.name_}',P.rat))));
        FRtrials.unitClass  = listout(idx).UnitClass;
        FAtrials.AssemClass = listout(idx).AssemClass;
        
        if ~isempty(FAtrials.FSC_Mean{s}) % if there are assemblies detected for this rat/area
            for iAss = 1:length(FAtrials.keepers{s})
                
                % Concordance of Assembly and constituent units' L/R preference
                temp = (FRtrials.iFRpref{s}  == FAtrials.Assempref{s}(iAss));
                temp = temp(intersect(FAcont.Task.units{s}{iAss},find(FRtrials.DPeak{s}>0)));    % Restrict by member units (and optionally decoding score)
                if FAtrials.DPeak{s}(iAss)>0
                    FAtrials.AssemUnitConcordLR{s}(iAss) = sum(temp)./numel(temp);               % Fraction of units sharing factor's L/R preference
                else
                    FAtrials.AssemUnitConcordLR{s}(iAss)     = NaN;
                end
                
                % 1) Concordance of Assembly and constituent units' activation pattern classification
                temp = (FRtrials.unitClass{s} == FAtrials.AssemClass{s}(iAss));
                temp = temp(intersect(FAcont.Task.units{s}{iAss},find(FRtrials.DPeak{s}>0)));    % Restrict by member units (and optionally decoding score)
                
                if FAtrials.DPeak{s}(iAss)>0
                    FAtrials.AssemUnitConcordType{s}(iAss) = sum(temp)./numel(temp); % Fraction of units sharing factor's decoding type
                else
                    FAtrials.AssemUnitConcordType{s}(iAss) = NaN;
                end
                
                % Number of each Assembly's member units in each class including un-informative cells
                memberUnits_ = intersect(FAcont.Task.units{s}{iAss},find(FRtrials.DPeak{s}>0));
                temp = FRtrials.unitClass{s}(memberUnits_);
                temp(isnan(temp)) = 0; % Make unclassified clusters = 0
                FRtrials.unitClasshistRawNumbers{s}(iAss,1:noClusts) = histc(temp,1:noClusts);
                FRtrials.unitClasshistRawNumbers{s}(iAss,noClusts+1) = numel(temp(temp==0));
                
                % Fraction of each Assembly's member units in each class
                temp = FRtrials.unitClass{s};
                temp = temp(intersect(FAcont.Task.units{s}{iAss},find(FRtrials.DPeak{s}>0))); % Restrict by member units (and optionally decoding score)
                FRtrials.unitClasshist{s}(iAss,:)    = histc(temp,1:noClusts)./numel(temp);
                FRtrials.unitClasshistRaw{s}(iAss,:) = histc(temp,1:noClusts);
             end                   
                
                %Collapse
                FRtrials.unitClasshistAssemType{s} = zeros(noClusts,noClusts);
                for iClass= 1:noClusts
                    % Unit classification for all assemblies of this type
                    temp_ = FRtrials.unitClasshistRaw{s}(FAtrials.AssemClass{s}==iClass,:);
                    if ~isempty(temp_)
                        temp  = nansum(temp_,1);
                    else
                        temp  = zeros(1,noClusts);
                    end
                    FRtrials.unitClasshistAssemType{s}(iClass,:) = temp./nansum(temp);
                end
            
        end
    end

    %% Output to data archive
    
    listout(iRat).AssemUnitConcordLR      = FAtrials.AssemUnitConcordLR;
    listout(iRat).AssemUnitConcordType    = FAtrials.AssemUnitConcordType;
    listout(iRat).unitClasshistAssemType  = FRtrials.unitClasshistAssemType;
    listout(iRat).unitClasshistRawNumbers = FRtrials.unitClasshistRawNumbers;
    listout(iRat).AssemClass              = FAtrials.AssemClass;
    
end
clear s iRat AssemID temp temp_ nAss 
%% Collate population data 
Collated = struct;

% L/R preference agreement
tempout = cell(1,3);
temp = {listout.AssemUnitConcordLR}';
for iRat = 1:length(RatList)
    tempout = [tempout;temp{iRat}];
end
tempout(1,:)=[];
for s = 1:3
    Collated.AssemUnitConcordLR{s} = cell2mat(cellfun(@transpose,tempout(:,s),'UniformOutput',false));
end

% Mean classification type agreement - average fraction across animals
Collated.AssemUnitConcordType=cell(1,3);
for s = 1:3
    Collated.AssemUnitConcordType{s} = zeros(noClusts,noClusts);
    i = 0;
    for iRat = 1:length(RatList)
        if ~isempty(listout(iRat).unitClasshistAssemType{s})
            temp = listout(iRat).unitClasshistAssemType{s};
            temp(isnan(temp)) = 0;

            Collated.AssemUnitConcordType{s} = Collated.AssemUnitConcordType{s} + temp;
            i=i+1;
        end
    end
    Collated.AssemUnitConcordType{s}=Collated.AssemUnitConcordType{s}./i;
end

% Mean classification type agreement - raw numbers for stats
tempout = cell(1,3);
temp = {listout.unitClasshistRawNumbers}';
for iRat = 1:length(RatList)
    tempout = [tempout;temp{iRat}];
end
tempout(1,:)=[];
for s = 1:3
    Collated.unitClasshistRawNumbers{s} = cell2mat(tempout(:,s));
    for iClass = 1:noClusts
        
    end
end


% Classification type
Collated.Unittype = cell(3,1);
Collated.Assemtype= cell(3,1);
tempUnits = cell(1,3);
tempAssems = cell(1,3);
tempU = {listout.UnitClass}';
tempA = {listout.AssemClass}';
for iRat = 1:length(RatList)
    for s = 1:3
    tempUnits{s}  = [tempUnits{s}; tempU{iRat}{s}];
    tempAssems{s} = [tempAssems{s};tempA{iRat}{s}];
    end
end

for s = 1:3
    Collated.Unittype{s} = tempUnits{s} ;
    Collated.Unittype{s}(:,2) = s*ones(size(Collated.Unittype{s}));
    
    Collated.Assemtype{s} = tempAssems{s};
    Collated.Assemtype{s}(:,2) = s*ones(size(Collated.Assemtype{s}));
end

% Classification type agreement - stats 1
for s = 1:3
    Assemtype_ = Collated.Assemtype{s}(:,1); Assemtype_(isnan(Assemtype_))=0;
    for AssemClass_ = 1:noClusts
        idx = find(Assemtype_ == AssemClass_);
        temphist_ = Collated.unitClasshistRawNumbers{s}(idx,1:noClusts);
            Collated.unitClasshistStats{s}(AssemClass_,1) = sum(temphist_(:,AssemClass_));
            Collated.unitClasshistStats{s}(AssemClass_,2) = sum(sum(temphist_(:,[1:(AssemClass_-1),(AssemClass_+1):noClusts])));
            Collated.unitClasshistStats{s}(AssemClass_,3) = sum(sum(temphist_));
            % Chi-squared test on each Assembly type: No. Units same vs. different class
            [Collated.unitClasshistStats{s}(AssemClass_,4),...  % reject H0?
             Collated.unitClasshistStats{s}(AssemClass_,5),...  % p value
             Collated.unitClasshistStats{s}(AssemClass_,6),...  % Chi-square value
             Collated.unitClasshistStats{s}(AssemClass_,7)]...  % Degrees of freedom
                 = proportionTest([Collated.unitClasshistStats{s}(AssemClass_,1),Collated.unitClasshistStats{s}(AssemClass_,2)],...
                                  [Collated.unitClasshistStats{s}(AssemClass_,3),Collated.unitClasshistStats{s}(AssemClass_,3)],true);
            
    end
    temp = sum(Collated.unitClasshistStats{s}(:,1:3));
    [a(1) a(2) a(3) a(4)] = proportionTest(temp(1:2),[temp(3) temp(3)],false);
    Collated.unitClasshistStatsSum{s}=a;
end

% Classification type agreement - stats 2
for s = 1:3
    Assemtype_ = Collated.Assemtype{s}(:,1); Assemtype_(isnan(Assemtype_))=0;
    for AssemClass_ = 1:noClusts
        idx = find(Assemtype_ == AssemClass_);
        temphist_ = Collated.unitClasshistRawNumbers{s}(idx,1:noClusts);
            Collated.unitClasshistFrac{s}{AssemClass_}(:,1) = (temphist_(:,AssemClass_)./sum(temphist_,2));
            Collated.unitClasshistFrac{s}{AssemClass_}(:,2) = sum(temphist_(:,[1:(AssemClass_-1),(AssemClass_+1):noClusts]),2)./sum(temphist_,2);
            if ~isempty(Collated.unitClasshistFrac{s}{AssemClass_}(:,1))
                Collated.unitClasshistFracStats{s}(AssemClass_) = ranksum(Collated.unitClasshistFrac{s}{AssemClass_}(:,1),Collated.unitClasshistFrac{s}{AssemClass_}(:,2));
            else
                Collated.unitClasshistFracStats{s}(AssemClass_) = NaN;
            end
    end
end

clear a i s iRat tempUnits  tempAssems temp tempU tempA tempout
%% plot cluster centroids

temp_x = (1:P.Ltr*2)*P.bw;
offset=4.5;
figure('color','w'); hold on
for iClass = 1:noClusts
    plot(temp_x,zscore(kmeans_.Centroids(iClass,:))+(iClass-1)*offset,'color',ClassColors{iClass},'LineWidth',3)
    text(32,(iClass-1)*offset,ClassNames{iClass},'color',ClassColors{iClass})
end
plot([10 10],[-5 -2],'g','LineWidth',1.5)
plot([20 20],[-5 -2],'r','LineWidth',1.5)
text(9,-3.5,'Sample','Rotation',90,'HorizontalAlignment','center','color','g')
text(19,-3.5,'Choice','Rotation',90,'HorizontalAlignment','center','color','r')

% plot([15 15],[-5 40],'w','LineWidth',5)
plot([15 15],[-5 -2],':k')
plot([33 35],[-3.5 -3.5],'k','LineWidth',3)
text(34,-2.5,'2s','HorizontalAlignment','center','color','k')

axis([0 40 -5 Inf])

axis off
%% plot types

% (1) Areas as a function of Units/Assems types 
figure('color','w','NumberTitle','off','name','Unit type per area')
for s= 1:2
    subplot(1,3,s); hold on
    temp = histc(Collated.Unittype{s}(:,1),1:noClusts);
    tempLabels = ClassNames(temp~=0);
    tempColors = ClassColors(temp~=0)';
    Hpie = pie(temp,repmat({''},1,numel(tempLabels)));
    hp = findobj(Hpie,'Type','patch');
    for iClass= 1:numel(hp )
        set(hp (iClass), 'FaceColor', tempColors{iClass},...
                        'EdgeColor', [0.6 0.6 0.6],...
                        'FaceAlpha', 0.6,...
                        'LineWidth', 1,...
                        'LineStyle','none');
    end
    
    scatter(0,0,5000,'w','filled') 
    text(0,0,{sprintf('n=%d',numel(Collated.Unittype{s}(:,1)));'Units'},'HorizontalAlignment','center')
    title(P.names{s})
    
    axis square
    axis off
    if s==2
        legend(tempLabels,'Location','SouthOutside','Orientation','horizontal')
        legend('boxoff')
    else
        legend(tempLabels,'Location','SouthOutside','Orientation','horizontal')
        legend('hide')

    end
end

figure('color','w','NumberTitle','off','name','Assembly type per area')
for s= 1:3
    subplot(1,3,s); hold on
    temp = histc(Collated.Assemtype{s}(:,1),1:noClusts);
    tempLabels = ClassNames(temp~=0);
    tempColors = ClassColors(temp~=0)';
    temp(temp==0)=[];
    Hpie = pie(temp,repmat({''},1,numel(tempLabels)));
    hp = findobj(Hpie,'Type','patch');
    for iClass= 1:numel(hp )
        set(hp (iClass), 'FaceColor', tempColors{iClass},...
                        'EdgeColor', [0.6 0.6 0.6],...
                        'FaceAlpha', 0.6,...
                        'LineWidth', 1,...
                        'LineStyle','none');
    end
    
    scatter(0,0,5000,'w','filled') 
    text(0,0,{sprintf('n=%d',numel(Collated.Assemtype{s}(:,1)));'Assemblies'},'HorizontalAlignment','center')
    title(P.names{s})
    
    axis square
    axis off
    if s==2
        legend(tempLabels,'Location','SouthOutside','Orientation','horizontal')
        legend('boxoff')
                legend('hide')

    else
        legend(tempLabels,'Location','SouthOutside','Orientation','horizontal')
        legend('hide')

    end
end

% (2) Units types as function of area
leg_ = P.names(1:2);
col_ = {[0.9 0.6 0.6],[0.6 0.9 0.6]};
figure('color','w','NumberTitle','off','name','Area by Unit type')
D = cell2mat(Collated.Unittype(1:2));
for iClass = 1:noClusts
    subplot(1,noClusts,iClass); hold on
    temp = histc(D(D(:,1)==iClass,2),1:2);
   
    tempLabels = leg_(temp~=0);
    tempColors = col_(temp~=0)';
    Hpie = pie(temp,{'',''});
    hp = findobj(Hpie,'Type','patch');
    for a= 1:numel(hp )
        set(hp (a), 'FaceColor', tempColors{a},...
                        'EdgeColor', [0.6 0.6 0.6],...
                        'FaceAlpha', 0.6,...
                        'LineWidth', 1,...
                        'LineStyle','none');
    end
    scatter(0,0,2000,'w','filled') 
    text(0,0,{sprintf('n=%d',numel(D(D(:,1)==iClass,2)));'Units'},'HorizontalAlignment','center')
    title(ClassNames{iClass},'color',ClassColors{iClass})
    if iClass==noClusts
        legend(leg_,'Location','SouthOutside','Orientation','horizontal')
        legend('boxoff')
    else
        legend(leg_,'Location','SouthOutside','Orientation','horizontal')
        legend('hide')
    end
    axis square
    axis off
end

% (3) Assem types as function of area
leg_ = P.names;
col_ = {[0.9 0.6 0.6],[0.6 0.9 0.6],[0.6 0.6 0.9]};
figure('color','w','NumberTitle','off','name','Area by Assembly type')
D = cell2mat(Collated.Assemtype);
for iClass = 1:noClusts
    subplot(1,noClusts,iClass); hold on
    temp = histc(D(D(:,1)==iClass,2),1:3);
   
    tempLabels = leg_(temp~=0);
    tempColors = col_(temp~=0)';
    Hpie = pie(temp,{'','',''});
    hp = findobj(Hpie,'Type','patch');
    for a= 1:numel(hp )
        set(hp (a), 'FaceColor', tempColors{a},...
                        'EdgeColor', [0.6 0.6 0.6],...
                        'FaceAlpha', 0.6,...
                        'LineWidth', 1,...
                        'LineStyle','none');
    end
    scatter(0,0,2000,'w','filled') 
    text(0,0,{sprintf('n=%d',numel(D(D(:,1)==iClass,2)));'Assems'},'HorizontalAlignment','center')
    title(ClassNames{iClass},'color',ClassColors{iClass})
    if iClass==noClusts
        legend(leg_,'Location','SouthOutside','Orientation','horizontal')
        legend('boxoff')
    else
        legend(leg_,'Location','SouthOutside','Orientation','horizontal')
        legend('hide')
    end
    axis square
    axis off
end; clear iClass
%% plot unit~Assembly LR concordance

col_ = {'k' , 'r'};
leg_ = {'Shared spatial preference','Different spatial preference'};

figure('color','w','NumberTitle','off','name','Assembly~Unit L/R outcome agreement')
for s= 1:3
    subplot(1,3,s); hold on
    temp = mean(Collated.AssemUnitConcordLR{s});
    
% 	Hpie = pie([temp;1-temp],{'Same','Different'});
    Hpie = pie([temp;1-temp],{'',''});
    hp = findobj(Hpie,'Type','patch');
    for iClass= 1:numel(hp)
        set(hp(iClass), 'FaceColor', col_{iClass},...
                        'EdgeColor', [0.6 0.6 0.6],...
                        'FaceAlpha', 1,...
                        'LineWidth', 1,...
                        'LineStyle','none');
    end
   	scatter(0,0,5000,'w','filled') 
    text(0,0,{sprintf('n=%d',numel(Collated.AssemUnitConcordLR{s}));'Assemblies'},'HorizontalAlignment','center')
    title(P.names{s})
    axis square
    axis off
    if s==2
        legend(leg_,'Location','SouthOutside','Orientation','horizontal')
        legend('boxoff')
    else
        legend(leg_,'Location','SouthOutside','Orientation','horizontal')
        legend('hide')

    end
end
%% Plot unit~Assembly type concordance

figure('name','Unit/Assembly type agreement','color','w')
for s=1:3
    subplot(1,3,s); hold on
    title([P.names{s}])
    temp = Collated.AssemUnitConcordType{s};
    temp(isnan(temp)) =0;
    imagesc(temp)
    for iClass = 1:noClusts
       rectangle('Position',[iClass-0.5 iClass-0.5 1 1],'LineWidth',2) 
    end
    rectangle('Position',[0.5 0.5 noClusts noClusts],'LineWidth',2) 

    caxis([0 0.5])
%     axis([1 noClusts 1 noClusts])
    colormap jet
    axis tight; axis square
%     axis off
    if s==1 
        set(gca,'YTick',[])
        for i=1:noClusts
            text(0,i,ClassNames{i},'color',ClassColors{i},'horizontalalignment','Right','FontWeight','bold')
        end
        text(-1,noClusts+2,{'Assembly';'Classification'},'FontWeight','bold','horizontalalignment','Center','FontSize',14)
    elseif s ==2 
            text(ceil(noClusts/2),-3,'Constituent units',...
            'FontWeight','bold',...
            'horizontalalignment','Center',...
            'FontSize',14)
    end
    for i=1:noClusts
            text(i,0,ClassNames{i},'color',ClassColors{i},...
                 'horizontalalignment','Right',...
                 'FontWeight','bold',...
                 'Rotation',45)
     end
     set(gca,'YTick',[],'XTick',[])
end


    



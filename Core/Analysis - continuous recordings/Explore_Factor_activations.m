% Wrapper function for getAssemPatterns
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
%
%
% Provides plots of these data if plot_flag is set. Use to explore
% the timing and prevalence of assembly activation patterns that occur
% above chance
%% Preamble 
clear all
close all
target= 'Task';
plot_flag = true;
if ispc
    home = 'C:\'; %[getenv('HOMEDRIVE') getenv('HOMEPATH')];
    pat = [home 'Analysis\AssemblyAnalysis\Sleep'];
    separator='\';
elseif ismac
    path(path,'/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/MichalDataAna')
    home = getenv('HOME');
    cd ([home '/Documents/Bristol_work/assembly_analysis/MichalDataAna'])
    pat = [home '/Documents/Bristol_work/assembly_analysis/Michal_data/raw/50ms_bins/'];
    separator='/';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/50ms_bins/'];
    separator='/';
end

% run on 1) iFRsc, 2) KDE iFR

%bw=0.01; iFRtyp='iFRsc';
%bw=0.05; iFRtyp='iFRsc';

bw=0.05; iFRtyp='iFR';
A=dir([pat separator '*PFC_iFR' num2str(bw*1e3) strcat('*', target) '*.mat']);


for event_id=1:length(A)
    k=findstr(A(event_id).name,'_'); r=findstr(A(event_id).name,'.');
    Name{event_id}=[A(event_id).name(1:k(1)) A(event_id).name(k(2)+1:r-1)];
end;
clear k
%% Load data
f=1%5% for f=1:length(A) % AKA Krysz long 2
fnOut=[pat separator target separator Name{f} '_AssemRes2']; load(fnOut,'FL')
fnOut=[pat separator target separator Name{f} '__FSCtemp']; load(fnOut,'FSC','ciHld','ciHsc');
%% Process Assemblies
Ass = [];
Ass = getAssemPatterns(Name{f},bw,FSC,ciHld,ciHsc,FL);
%% Optional plotting
if plot_flag
%% Plot soloists vs. choristers 
figure('color','w'); 
cmap = {'r','g','b'};
for s = 1:3
    subplot(1,2,1); hold on 
        bar(s,Ass.units.FractionSingleAssem(s),'EdgeColor', cmap{s},'FaceColor', cmap{s}); 
        set(gca,'XTickLabel',Ass.titles,'XTick',[1 2 3],'Ylim',[0 1]);
        title('Cells assigned to single Assembly')
        ylabel('Fraction')
        set(gca,'ylim',[0 1])
    subplot(1,2,2); hold on
        plot(Ass.units.AssOverlap_hist_bins,Ass.units.AssOverlap_hist(:,s),'color',cmap{s},'LineWidth',2); 
        title('Number of factors units shared among'); 
        ylabel('Fraction'); 
        xlabel('No. Factors');
        legend(Ass.titles); legend('boxoff')
        set(gca,'ylim',[0 1])
end
%% Plot pattern transition probability
try
    load('whiteHot');
catch
    cmap=(hot(64));
end
for s = 1:length(Ass.FSC)      % across all areas/assem types
    figure('Name',[Ass.filename ' : ' Ass.titles{s} ' Assemblies'])
    subplot(2,1,1) % Plot assembly sequence
    plot(Ass.ptn.seq{s}(:,1)*Ass.bw,Ass.ptn.seq{s}(:,2),'k')
    hold on; axis tight
    xlabel('Time (s)')
    ylabel('Assembly no.')
    set(gca, 'ylim',[1 Ass.nAss(s)],...
        'ytick',[1:Ass.nAss(s)],...
        'yticklabel',Ass.AssLabels{s});
    
    subplot(2,2,3); hold on
    temp=Ass.ptn.transMatrix{s}; temp_hist=histc(temp(1:end),0:0.05:1);
    plot([Ass.ptn.ciHTransition(s) Ass.ptn.ciHTransition(s)],[0 max(temp_hist./sum(temp_hist))],'r','LineWidth',1.5)
    plot(0:0.05:1,temp_hist./sum(temp_hist),'k','LineWidth',1.5)
    legend('95th percentile of shuffled patterns') , legend('boxoff')
    xlabel('Transition probability')
    ylabel('Distribution')
    subplot(2,2,4)
    temp(temp<Ass.ptn.ciHTransition(s))=0;
    pcolor(flipud(rot90(temp)))%; set(gca,'YDir','normal')
    set(gca, 'ylim', [1 Ass.nAss(s)],...
        'xtick',[1:Ass.nAss(s)],'xticklabel',Ass.AssLabels{s},...
        'ytick',[1:Ass.nAss(s)],'yticklabel',Ass.AssLabels{s},...
        'clim',[0 1])
    
    colormap(cmap)
    colorbar
    axis square
    title('Assembly sequence transition matrix')
    xlabel('Event (t)')
    ylabel('Event (t+1)')
    % e.g. p(Assem. 3 -> Ass. 2) ...  Ass.pattern.trans_matrix{s}(3,2)
end
%% Plot within-pattern peak histograms
figure('color','w','name',[target ' assembly activation rate']); 
cmap = {'r','g','b'};
hold on
for s=1:3
    plot((Ass.pks.repHistBins*Ass.bw).^-1,mean(Ass.pks.repHist{s}),'color',cmap{s},'LineWidth',2)
end
    set(gca,'XScale','log')       
xlabel('Average Assembly activation rate (Hz)')
ylabel('Distribution')
legend(Ass.titles); legend('boxoff')
%% Plot assembly time autocorrelation
figure('Name',[Ass.filename ' :  Assembly activation time autocorrelation'])
for s = 1:length(Ass.FSC)      % across all areas/assem types
    % Plot autocorrelogram matrix
    subplot(1,3,s), hold on
    title( Ass.titles{s})
    plot(Ass.Corr.lags*Ass.bw,Ass.Corr.autoCorr{s})
    plot([0 0],[0 1],':k')
    % axis([min(Ass.Corr.lags) max(Ass.Corr.lags) 0 1])
    axis([-50 50 0 1])
    xlabel('Time (s)')
    if s==1, ylabel('Autocorrelation.'), end
end
%% plot assembly time cross-correlograms
% for s = 3%:length(Ass.FSC)
%     
% 
%     
%     figure('Color','w','Name',[Ass.filename ' : '  Ass.titles{s} ' :  Assembly activation time cross-correlation'])
%     for plot_id=1:Ass.nAss(s)^2
%         subplot(Ass.nAss(s),Ass.nAss(s),plot_id)
%         hold on
%         plot([0 0],[-.1 1],':k')
%         plot([min(Ass.Corr.lags) max(Ass.Corr.lags)],[0 0],':k')
%         %             text(.1,.8, Ass.pattern.labels{s}{plot_id},'FontSize',6)
%         title(Ass.pattern.labels{s}{plot_id},'FontSize',6)
%         if Ass.pattern.trans_matrix{s}(plot_id)>Ass.pattern.ciH_transition(s)
%             plot(Ass.Corr.lags,Ass.Corr.SC{s}{plot_id},'LineWidth',1.2,'color','b','LineWidth',1.2)
%         else
%             plot(Ass.Corr.lags,Ass.Corr.SC{s}{plot_id},'LineWidth',1.2,'color',[.6 .6 .6])
%         end
%         box off ;axis off
%         set(gca,'Xtick',[]);set(gca,'Ytick',[])
%         axis([min(Ass.Corr.lags) max(Ass.Corr.lags) -.1 1])
%         % xlabel('Time (s)')
%     end
% end
%% plot pattern duration
for s = 1:length(Ass.FSC)      % across all areas/assem types
    figure('Name',[Ass.filename ' : ' Ass.titles{s} ' pattern durations'],...
        'Color','w');
    for plot_id=1:Ass.nAss(s)^2
        subaxis(Ass.nAss(s),Ass.nAss(s),plot_id,'Spacing',0.01,'Margin',0.05)
        hold on
        plot(Ass.ptn.DurBins*Ass.bw,Ass.ptn.DurHist{s}{plot_id},'LineWidth',1.5)
        %             plot([0 0],[0 1],':k')
        set(gca,'Xtick',[]);set(gca,'Ytick',[]); box on;
        axis([min(Ass.ptn.DurBins)*Ass.bw,max(Ass.ptn.DurBins)*Ass.bw, 0 .2])
        %             if plot_id>Ass.nAss(s), title(num2str(plot_id)), end
    end
end
%% plot pattern repeat interval
for s = 1:length(Ass.FSC)      % across all areas/assem types
    figure('Name',[Ass.filename ' : ' Ass.titles{s} ' pattern repeat interval'],...
        'Color','w');
    for plot_id=1:Ass.nAss(s)^2
        %             subplot(Ass.nAss(s),Ass.nAss(s),plot_id)
        subaxis(Ass.nAss(s),Ass.nAss(s),plot_id,'Spacing',0.01,'Margin',0.02)
        hold on
        plot(Ass.ptn.RepBins*Ass.bw,Ass.ptn.RepHist{s}{plot_id},'LineWidth',1.5)
        %             plot([0 0],[0 1],':k')
        set(gca,'Xtick',[]);set(gca,'Ytick',[]); box on;
        axis([min(Ass.ptn.RepBins)*Ass.bw max(Ass.ptn.RepBins)*Ass.bw, 0 .5])
        %             if plot_id>Ass.nAss(s), title(num2str(plot_id)), end
    end
end


% write_matrix_to_pajek(Ass.pattern.trans_matrix{s}, [pat separator Ass.Name '_' Ass.titles{s} '_AssemTransitions'])
% end
%% plot sequence probability vs. duration
for s = 1:length(Ass.FSC)      % across all areas/assem types
    figure('Name',[Ass.filename ' : ' Ass.titles{s} ' pattern probability vs mean duration'],...
           'Color','w');
       subplot(1,2,1)
       scatter(Ass.ptn.DurMedian{s}(1:end)*Ass.bw,Ass.ptn.transMatrix{s}(1:end))
%        set(gca,'XScale','log')
       ylabel('Pattern probability')
       xlabel('Duration (s)')
       axis([0.1 2 0 0.25])
%               set(gca,'XScale','log')
%        set(gca,'YScale','log')
       subplot(1,2,2)
       scatter((Ass.ptn.RepMedian{s}(1:end)*Ass.bw),Ass.ptn.transMatrix{s}(1:end))
       set(gca,'XScale','log')
       set(gca,'YScale','log')
       ylabel('Pattern probability')
       xlabel('Repeat rate (Hz)')
%        axis([1e-3 1 0 0.5])

%     figure
%        stem3(Ass.pattern.PtnDur_median{s}(1:end)*Ass.bw,...
%              Ass.pattern.PtnRep_median{s}(1:end)*Ass.bw,...
%              Ass.pattern.trans_matrix{s}(1:end))
%              xlabel('Repeat rate (s)')
%              ylabel('pattern duration (s)')
%              zlabel('Pattern probability')
% %              set(gca,'YScale','log')

%% Plot as transition probability
figure;
x    = Ass.ptn.transMatrix{s}(1:end);
y    = Ass.ptn.DurMedian{s}(1:end)*Ass.bw*1e3;

Msize = ((Ass.ptn.DurSEM{s}(1:end)*Ass.bw).^-1);
% Msize(x<4)=NaN;
subplot(1,2,1); hold on
scatter(x,y,Msize,Msize,'fill')
%     for label_id=1:Ass.nAss(s)
%         text(Ass.pattern.RepCount{s}(label_id),...
%              Ass.pattern.PtnDur_median{s}(label_id)*Ass.bw,...
%              Ass.pattern.labels{s}{label_id},'FontSize',6)
%     end
xlabel('Pattern tranisiton probability')
ylabel('Average pattern duration (ms)')
title('Size/Colour: Pattern length variability')%;legend('boxoff'); %legend('location','SouthOutside')

subplot(1,2,2)
x    = Ass.ptn.transMatrix{s}(1:end);
y    = 1./(Ass.ptn.RepMedian{s}(1:end)*Ass.bw);
Msize = ((Ass.ptn.RepSEM{s}(1:end)*Ass.bw).^-1)*1e2;
% Msize(x<4)=NaN;

scatter(x,y,Msize,Msize,'fill')
xlabel('Pattern tranisiton probability')
ylabel('Average repeat rate (Patterns per second)')
% set(gca,'XScale','log')       
set(gca,'YScale','log')       
title('Size/Colour: Pattern repeat rhythmicity')
% legend('Size: Pattern duration peridocity');legend('boxoff');legend('location','SouthEast')
%% Plot as pattern count
% s=3
figure;
x    = Ass.ptn.RepCount{s}(1:end);
y    = Ass.ptn.DurMedian{s}(1:end)*Ass.bw*1e3;
Msize = ((Ass.ptn.DurSEM{s}(1:end)*Ass.bw).^-1);
Msize(x<4)=NaN;
subplot(1,2,1); hold on
scatter(x,...
        y,...
        Msize,...
        Msize,'fill')
%     for label_id=1:Ass.nAss(s)
%         text(Ass.pattern.RepCount{s}(label_id),...
%              Ass.pattern.PtnDur_median{s}(label_id)*Ass.bw,...
%              Ass.pattern.labels{s}{label_id},'FontSize',6)
%     end
xlabel('Pattern repeat count')
ylabel('Average pattern duration (ms)')
title('Size/Colour: Pattern length variability')%;legend('boxoff'); %legend('location','SouthOutside')

subplot(1,2,2)
x    = Ass.ptn.RepCount{s}(1:end);
y    = 1./(Ass.ptn.RepMedian{s}(1:end)*Ass.bw);
Msize = ((Ass.ptn.RepSEM{s}(1:end)*Ass.bw).^-1)*1e3;
Msize(x<4)=NaN;

scatter(x,y,Msize,Msize,'fill')
xlabel('Pattern repeat count')
set(gca,'XScale','log')       
ylabel('Average repeat rate (Patterns per second)')
set(gca,'YScale','log')       
title('Size/Colour: Pattern repeat rhythmicity')
% legend('Size: Pattern duration peridocity');legend('boxoff');legend('location','SouthEast')
%% Plot as pattern count ... highlight self-repeating patterns
% s=1;
figure;
x    = Ass.ptn.RepCount{s}(1:end);
y    = Ass.ptn.DurMedian{s}(1:end)*Ass.bw*1e3;
Msize = ((Ass.ptn.DurSEM{s}(1:end)*Ass.bw).^-1);
Msize(x<4)=NaN;
subplot(1,2,1); hold on
scatter(x,...
        y,...
        100,...
        Msize,'fill')
    
set(gca,'XScale','log')      
set(gca,'YScale','log')       
% % Plot only self-self
% diag(size(Ass.pattern.PtnDur_sem{s}))
% x    = diag(Ass.pattern.RepCount{s});
% y    = diag(Ass.pattern.PtnDur_median{s})*Ass.bw*1e3;
% Msize = ((diag(Ass.pattern.PtnDur_sem{s})*Ass.bw).^-1);
% Msize(x<4)=NaN;
% scatter(x,...
%         y,...
%         100,...
%         Msize,'fill')
    
%     for label_id=1:Ass.nAss(s)
%         text(Ass.pattern.RepCount{s}(label_id),...
%              Ass.pattern.PtnDur_median{s}(label_id)*Ass.bw,...
%              Ass.pattern.labels{s}{label_id},'FontSize',6)
%     end
xlabel('Pattern repeat count')
ylabel('Average pattern duration (ms)')
title('Colour: Pattern duration stability')%;legend('boxoff'); %legend('location','SouthOutside')

subplot(1,2,2)
x    = Ass.ptn.RepCount{s}(1:end);
y    = 1./(Ass.ptn.RepMedian{s}(1:end)*Ass.bw);
Msize = ((Ass.ptn.RepSEM{s}(1:end)*Ass.bw).^-1)*1e3;
Msize(x<4)=NaN;

scatter(x,y,100,Msize,'fill')
xlabel('Pattern repeat count')
set(gca,'XScale','log')      
set(gca,'YScale','log')       
ylabel('Average repeat rate (Patterns per second)')
title('Colour: Pattern repeat rhythmicity')
% legend('Size: Pattern duration peridocity');legend('boxoff');legend('location','SouthEast')
%% Plot rhythmicity of each factor
figure; hold on
scatter((Ass.pks.repMean{s}*Ass.bw),Ass.pks.repSTD{s}*Ass.bw,diag(Ass.ptn.RepCount{s})*10)
% plot([1 60],[1 10],':k')
title('Assembly activation rythmicity')
xlabel('Mean activation interval (s)')
ylabel('repeat interval variability (s)')
% axis([0 60 0 10])
set(gca,'XScale','log')
end
end
%% Save a graph layout
r=Ass.ptn.transMatrix{4};
for r_id=1:sum(Ass.nAss)
    if r_id<=Ass.nAss(1)
        r_labels{r_id}=['PFC' num2str(r_id)];
        r_partition(r_id)=1;
    elseif r_id<=Ass.nAss(1)+ Ass.nAss(2)
        r_labels{r_id}=['HP' num2str(r_id-Ass.nAss(1))];
        r_partition(r_id)=2;
    else
       r_labels{r_id}=['Joint' num2str(r_id-(Ass.nAss(1)+Ass.nAss(2)))];
       r_partition(r_id)=3;
    end
    % r_shapes{r_id}='ellipse';
    
end
% r = (r - min(r(1:end))) ./ ( max(r(1:end)) - min(r(1:end)) ); % scale between 0~1

EdgeL=adj2gephilab2('pre_graph',mat2gray(r),r_partition,r_labels);
% EdgeL=adj2gephilab('pre_graph',r,r_partition);

% adj2pajek2(r,'pre_graph',...
%            'nodeNames',r_labels,...
%            'shapes',r_shapes,...
%            'partition',double(r_partition))
clear r r_labels r_partition r_shapes r_id EdgeL
%% Save structures
fnOut=[pat separator target separator Name{f} '_Ass']; save(fnOut,'Ass')
clear temp cmap columnStates f Msize plot_id s target temp_hist transitionMatrix x y separator Name A fnOut pat
clear Rep_id PtnRep PtnDur PtnTimes ws AssID AssID_  k
clear alpha BS_max BSid BSperm bw ciHld ciHsc event_id FL fnOut FSC home iFRtyp  r s   seq_temp  trans_matrixBS plot_flag
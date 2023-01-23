% Takes factor model computed for the task epoch, evaluates assembly activation profile of unseen data by casting another recording period into it
%
% Can take assemblies detected from task dataset either in trials cut-outs or as continuous FSCs

%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
%% Preample, load 
clear all
flags.makeGraph = false;        % save a Gephi ball graph
flags.useEventFSCs = false;     % Use cutout rather than continuous times for factor detection
flags.plotOnline   = false;     % Sanity check: Line up FRs and unit timestamps
p.target = 'Task';
p.rat    = 'JaroslawLONG2';
% p.rat    = 'KrzesimirLONG2';
p.expt   = 'LONG';
p.sep='\';
disp(['*** p.target = ' p.target ' / p.rat = ' p.rat ' / p.expt = ' p.expt ])
                                                    % Where to find the...
p.pat  = 'C:\Analysis\AssemblyAnalysis\Sleep\';     % processed sleep assembly activations
p.pat2 = 'C:\Analysis\AssemblyAnalysis\raw';        % raw spike times
p.pat3 = [p.pat2 p.sep 'KDE_bins'];                 % calculated firing rates for Task period
% p.pat4 = [p.pat3 p.sep p.expt p.sep];             % event cutout factor analysis
% p.pat = '/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/Sleep/';
% cd(p.pat)



A(1) = dir([p.pat '*' p.rat '*PFC' '*PreSleep.mat']);
A(2) = dir([p.pat '*' p.rat '*HP'  '*PreSleep.mat']);
C(1) = dir([p.pat '*' p.rat '*PFC'  '*PostSleep.mat']);
C(2) = dir([p.pat '*' p.rat '*HP'   '*PostSleep.mat']);
p.names={'PFC','HP','Joint HP~PFC'};

disp(['*** loading Continous dataset...'])
for s=1:2
    FR.Pre{s}  = load([p.pat p.sep A(s).name]);
    FR.Post{s} = load([p.pat p.sep C(s).name]);
end

if flags.useEventFSCs 
    
    disp('*** Overriding Task iFRs with event cutout dataset')
    B(1) = dir([p.pat3 p.sep p.rat '*PFC*']);
    B(2) = dir([p.pat3 p.sep p.rat '*HP*']);
    for s=1:2
        FR.Task{s}  =  load([p.pat3 p.sep B(s).name]);
    end
    
else
    
    B(1) = dir([p.pat '*' p.rat '*PFC'  '*Task.mat']);
    B(2) = dir([p.pat '*' p.rat '*HP'   '*Task.mat']);
    for s=1:2
        FR.Task{s}  =  load([p.pat p.sep B(s).name]);
    end
    
end
disp('*** ...Done')
clear A B C s
%% Load FA model from the Task epoch
disp('*** Loading task factor model...')
if flags.useEventFSCs % use cutout
    disp('*** Using trial-cutout factor model')
     fn1=dir([p.pat3 p.sep p.expt p.sep '*' p.rat '*iFR50_AssemRes2.mat']);
     fn2=dir([p.pat3 p.sep p.expt p.sep '*' p.rat '*iFR50__FSCtemp.mat']);
	 A = load([p.pat3 p.sep p.expt p.sep fn1(1).name],'FL','psix','nassem');
     B = load([p.pat3 p.sep p.expt p.sep fn2(1).name],'FSC','ciHld','ciHsc');
else % use continuous
    disp('*** Using continuous factor model')
    fnOut=dir([p.pat p.target p.sep '*' p.rat '*']);
    A = load([p.pat p.target p.sep fnOut(1).name],'FL','psix','nassem','TmtxS');
    B = load([p.pat p.target p.sep fnOut(2).name],'FSCsel','ciHld','ciHsc','units','unit_IDs','keepers');
end
A.TmtxS{3}=A.TmtxS{1};
% These extracted FA model paramters are based on activity in the task period (either continuous or event cutout)
% NB alpha values are in ciV=[0.1 0.05 0.01 5e-3 1e-3]; Daniel uses thresholds >ciHld(s,3)) and >ciHsc(s,2)
 for s=1:3
    Ass.Task.nassem_(s) = A.nassem{s}(3);
    Ass.Task.FL{s}     = A.FL{s}{Ass.Task.nassem(s)};
    Ass.Task.psix{s}   = A.psix{s}{Ass.Task.nassem(s)};
 end
 Ass.Task.units    = B.units;
 Ass.Task.unit_IDs = B.unit_IDs;
 Ass.Task.keepers  = B.keepers;
 Ass.Task.TmtxS    = A.TmtxS;
 Ass.Task.FSC      = B.FSCsel;
 Ass.Task.ciHld = B.ciHld; %(:,3); % >1%
 Ass.Task.ciHsc = B.ciHsc; %(:,2); % >5%

 % 18/10/16 - Corrected to exclude single-unit factors
 Ass.Task.nassem = cell2mat(cellfun(@numel,Ass.Task.units,'UniformOutput',false));
 disp('*** ...Done')
clear fnOut fn1 fn2 A B s
%% Work out which units made into FA 

A=cell2mat((cellfun(@size,Ass.Task.FL,'UniformOutput',false)'));
disp(['*** Factor model No. Units = ', num2str(A(:,1)')])
disp(['*** Factor model No. Factors = ', num2str(A(:,2)')])
clear A

p.bw=0.05;       % Binwidth
p.minFR=0.1;     % minimal acceptable firing srate
p.critCvBW=1e6;  % critical max bound on kernel width drift over time (variance)
p.twin=10;
p.iFRtyp='iFR';  

if flags.useEventFSCs % use cutout: only cells included in Task FR matrix
    A=dir([p.pat3 p.sep  p.rat '*PFC_iFR' num2str(p.bw*1e3) '*.mat']);
    for event_id=1:length(A)
        k=findstr(A(event_id).name,'_'); r=findstr(A(event_id).name,'.');
        fileNames{event_id}=[A(event_id).name(1:k(1)) A(event_id).name(k(2)+1:r-1)];
    end;
    fn=[p.pat3 p.sep A(1).name];
    [~,~,~,~,Ass.Task.unit_IDs]=SelTrialsCells(fn,p.twin,p.minFR,p.critCvBW,p.iFRtyp);

else % use continuous: choice based on whole dataset
    A=dir([p.pat  p.sep p.rat '*PFC_iFR' num2str(p.bw*1e3) strcat('*', p.target) '*.mat']);
    for event_id=1:length(A)
        k=findstr(A(event_id).name,'_'); r=findstr(A(event_id).name,'.');
        fileNames{event_id}=[A(event_id).name(1:k(1)) A(event_id).name(k(2)+1:r-1)];
    end;
    fn=[p.pat A(1).name];
    [~,~,Ass.Task.unit_IDs]=SelCellsSleep(fn,p.twin,p.minFR,p.critCvBW,p.iFRtyp);
    
end

% reconstruct parameters for joint assems
Ass.Task.unit_IDs{3}  = [Ass.Task.unit_IDs{1},Ass.Task.unit_IDs{2}+max(Ass.Task.unit_IDs{1})];
FR.Pre{3}.iFR         = [FR.Pre{1}.iFR,FR.Pre{2}.iFR];
FR.Task{3}.iFR        = [FR.Task{1}.iFR,FR.Task{2}.iFR];
FR.Post{3}.iFR        = [FR.Post{1}.iFR,FR.Post{2}.iFR];
disp(['*** Confirming No. Units = ' num2str(cellfun(@length,Ass.Task.unit_IDs))])
[~, k] = cellfun(@size,Ass.Task.FSC);
disp(['*** Confirming No. Factors = ' num2str(k)])
clear A k r fn event_id  s 
%% load original spike times to recalculate mean firing rates
load([p.pat2 p.sep p.rat '.mat'],'HPcells','PFCcells')
for s=1:2
    avgFR_Pre =[]; avgFR_Task =[]; avgFR_Post =[]; 
    
    for u= Ass.Task.unit_IDs{s}
        disp(strcat('Unit no ', num2str(u)))
        ST=[];
        switch s
            case  1
                ST=PFCcells{u}.t*1e-4;
            case 2
                ST=HPcells{u}.t*1e-4;
        end
        ST_Pre=ST; ST_Task=ST; ST_Post=ST;
        % ST_Pre(ST_Pre<FR.Pre{s}.Tmtx(1) | ST_Pre>FR.Pre{s}.Tmtx(end))= [];
        % ST_Task(ST_Task<FR.Task{s}.Tmtx(1)   | ST_Task>FR.Task{s}.Tmtx(end))= [];
        % ST_Post(ST_Post<FR.Post{s}.Tmtx(1) | ST_Post>FR.Post{s}.Tmtx(end))= [];
        
        % hardcoded time ranges:
        if strcmp(p.rat,'JaroslawLONG1')
            tRanges.Pre  = ([529434636,  3651647307])*1e-6;	
            tRanges.Task = ([3901833652, 8598314309])*1e-6;
            tRanges.Post = ([9044732395, 1234563230])*1e-6;
        elseif strcmp(p.rat,'JaroslawLONG2')
            tRanges.Pre  = ([622141955,  3758558229])*1e-6;	
            tRanges.Task = ([3901833652, 8598314309])*1e-6;
            tRanges.Post = ([9028418876, 12332882982])*1e-6;    
        elseif  strcmp(p.rat,'KrzesimirLONG2')
            tRanges.Pre  = ([882560060,  4186836267])*1e-6;	
            tRanges.Task = ([4326852193, 9038811544])*1e-6;
            tRanges.Post = ([9391117187, 12693345168])*1e-6;
        end
        
        ST_Pre(ST_Pre<tRanges.Pre(1)    | ST_Pre>tRanges.Pre(2))  = [];
        ST_Task(ST_Task<tRanges.Task(1) | ST_Task>tRanges.Task(2))= [];
        ST_Post(ST_Post<tRanges.Post(1) | ST_Post>tRanges.Post(2))= [];
        
        if  ~isempty(ST_Pre)
            avgFR_Pre(u)=length(ST_Pre)/range(ST_Pre);
        else
            avgFR_Pre(u) = 0;
        end
        if  ~isempty(ST_Task)
            avgFR_Task(u)=length(ST_Task)/range(ST_Task);
        else
            avgFR_Task(u) = 0;
        end
        if  ~isempty(ST_Post) 
            avgFR_Post(u)=length(ST_Post)/range(ST_Post);
        else
            avgFR_Post(u)= 0;
        end
                   
    end
    avgFR_Pre(avgFR_Pre==0)  =NaN;
    avgFR_Task(avgFR_Task==0)=NaN;
    avgFR_Post(avgFR_Post==0)=NaN;
    
    FR.Pre{s}.avgFR_Pre   = avgFR_Pre;
    FR.Task{s}.avgFR_Task = avgFR_Task;
    FR.Post{s}.avgFR_Post = avgFR_Post;
end

FR.PFCcells=PFCcells;FR.HPcells=HPcells;
if flags.plotOnline    
    figure; hold on
    tempTimes=PFCcells{1}.t;
    scatter(tempTimes,ones(length(tempTimes),1),'.k')
    plot(tRanges.Pre*1e4,[0.5 0.5])
    plot(tRanges.Task*1e4,[0.5 0.5])
    plot(tRanges.Post*1e4,[0.5 0.5])
    plot(FR.Pre{1}.Tmtx(1:end-1)*1e4,zscore(FR.Pre{1}.iFR(:,1)))
    plot(FR.Task{1}.Tmtx(1:end-1)*1e4,zscore(FR.Task{1}.iFR(:,1)))
    plot(FR.Post{1}.Tmtx(1:end-1)*1e4,zscore(FR.Post{1}.iFR(:,1)))
    axis([min(tempTimes) max(tempTimes) 0 2])
    clear tempTimes
    end
if flags.useEventFSCs 
    
    disp('*** Overriding mean Task FRs with mean from event cutouts')
    %%% catch the mean firing p.rates while we're at it...
    A=dir([p.pat3 p.sep  p.rat '*PFC_iFR' num2str(p.bw*1e3) '*.mat']);
        temp = load([p.pat3 p.sep A.name],'avgFR');
        avgFRtemp{1}=temp.avgFR;
    A=dir([p.pat3 p.sep  p.rat '*HP_iFR' num2str(p.bw*1e3) '*.mat']);
        temp = load([p.pat3 p.sep A.name],'avgFR');
        avgFRtemp{2}=temp.avgFR ;
    avgFRtemp{3} = [avgFRtemp{1}, avgFRtemp{2}]; 
    clear A 
    %%%
    
    for s=1:3
         FR.Task{s}.avgFR_Task = avgFRtemp{s};
    end
    clear avgFRtemp
end

FR.Task{3}.avgFR_Task = [FR.Task{1}.avgFR_Task, FR.Task{2}.avgFR_Task];

clear s u avgFR_Pre avgFR_Task avgFR_Post ST_Pre ST_Task ST_Post ST  HPcells PFCcells temp A
%% demonstrate that firing rate to factor reconstruction works
s=1; iAssem=1; iUnit=1; xRange=[1 10000];


% recreate firing p.rate from FA parameters
% % FR_(model) = FRmean_model + FL_model * FSC_model + Psi_model
% F.nUnits  = length(unit_IDs{s});
% F.nPoints = size(FSC{s},1);
% F.FL   = FL{s}(iUnit,:)'; %F.FL  = repmat(F.FL,1,F.nPoints);
% F.FSC  = FSC{s}(:,1)';    %F.FSC = repmat(F.FSC,size(F.FL,1),1);
% F.mean = mean(FR.Task{s}.iFR(:,unit_IDs{s}(iUnit)))'; %F.mean = repmat(F.mean', 1 , size(F.FSC,1),1);   
% FR_ = F.mean +  F.FL .*  F.FSC  +  repmat(psix{s}(:,1)      ,1,size(FSC{s},1),1);
  
% recreate Factor scores from firing p.rate and FA parameters   
% FSC_eval = (FR_eval - mean_eval - Psi_model) / FL_model 

    F.FR_   = FR.Task{s}.iFR(:,Ass.Task.unit_IDs{s});
    %F.mean_ = repmat(mean(F.FR_,1),length(FR.Task{s}.iFR(:,unit_IDs{s})),1);
    F.mean_ = repmat(mean(F.FR_,1),length(FR.Task{s}.iFR(:,Ass.Task.unit_IDs{s})),1);
    F.psix_ = repmat(Ass.Task.psix{s}',length(F.FR_),1);
    F.FL_   = Ass.Task.FL{s}';
    F.FSC_  = (F.FR_ - F.mean_ -F.psix_)/F.FL_;
    FSC_  = F.FSC_; 
    clear F



figure('color' ,'w','position',[20 500 1400 400]);
% subplot(2,1,1); hold on
%     title('Construct  firing p.rates from factor model')
%     plot(zscore(FR.Task{s}.iFR(:,iUnit)),'b')
%     plot(zscore(FR_(iUnit,:)),'r')
%     axis([xRange -Inf Inf])
% subplot(2,1,2); 
    subplot('position',[0.1 0.1 0.55 0.8])
    hold on
    title('Factor model constructed from firing rates ')
    plot(zscore(Ass.Task.FSC{s}(:,iAssem)),'b')
    plot(zscore(FSC_(:,iAssem)),'r')
    axis([xRange -Inf Inf])
    legend({'Original firing rate','Reconstructed firing rate'}); legend boxoff
    axis off
    subplot('position',[0.7 0.1 0.3 0.8])
    plot(xcorr(zscore(Ass.Task.FSC{s}(:,iAssem)),zscore(FSC_(:,iAssem)),1000),'k')
    axis off; box off
    title('x-corr: original vs. reconstructed')
    clear iUnit iAssem s FSC_ xRange
%% Reconstruct Pre and Post sleep FSCs using Task FA parameters

% Reconstruct data using the following FA model:
% FSC = (FR-mean-Psi)/FL
epoch={'Pre', 'Post'};
for iA = 1:2
    for s = 1:3
        clear F
        eval(cell2mat([ 'F.FR_ = FR.' epoch(iA) '{s}.iFR(:,Ass.Task.unit_IDs{s},1);']))
        F.mean_ = repmat(FR.Task{s}.avgFR_Task(Ass.Task.unit_IDs{s}),size(F.FR_,1),1); 
        % F.mean_ = repmat(mean(F.FR_),size(F.FR_,1),1);
        F.psix_ = repmat(Ass.Task.psix{s}',size(F.FR_,1),1);
        F.FL_   = Ass.Task.FL{s}';
        F.FSC_  = (F.FR_ - F.mean_ -F.psix_)/F.FL_;    
        eval(cell2mat(['Ass.'  epoch(iA)  '.FSCproj{s} = F.FSC_;'])) 
    end
end
clear F epoch iA s
%% Process Reconstructed assemblies 
f=1;
Ass.Pre.patterns   = getAssemPatterns(fileNames{f},p.bw,Ass.Pre.FSCproj,Ass.Task.ciHld,Ass.Task.ciHsc,Ass.Task.FL);
Ass.Task.patterns  = getAssemPatterns(fileNames{f},p.bw,Ass.Task.FSC,Ass.Task.ciHld,Ass.Task.ciHsc,Ass.Task.FL);
Ass.Post.patterns  = getAssemPatterns(fileNames{f},p.bw,Ass.Post.FSCproj,Ass.Task.ciHld,Ass.Task.ciHsc,Ass.Task.FL);
clear f 
%% Example Plot - all assemblies in different colours
if flags.plotOnline
    figure('color','w','name',[p.target 'Assembly activation scores']); 
    hold on
    s=3;
    cmap = jet(Ass.Task.nassem(s));
    ymax = 20;
    a=[]; b=[];
    for i=1:Ass.Task.nassem(s)
        a = [a ; i*ones(length(Ass.Task.patterns.pks.time{s}{i}),1)];
        b = [b ; Ass.Task.patterns.pks.time{s}{i}];
    end
    c = [a,b];
    c = sortrows(c,2);
    a = c(:,1); b = c(:,2);
    % for iPeak = 1:length(a)-1 % plot shaded areas
    %     t1 = b(iPeak)/3600*bw;
    %     t2 = b(iPeak+1)/3600*bw;
    %         h=area([t1 t2],[ymax ymax],'Linewidth',1.5,...
    %                        'FaceColor',cmap(a(iPeak),:),...
    %                        'EdgeColor',cmap(a(iPeak),:),...
    %                        'LineStyle','none');
    % end
    for iPeak = 1:length(a) % plot shaded stripe
        t1 = b(iPeak)/3600*p.bw-1/3600*p.bw;
        t2 = b(iPeak)/3600*p.bw+1/3600*p.bw;
            h=area([t1 t2],[ymax ymax],'Linewidth',1.5,...
                           'FaceColor',cmap(a(iPeak),:),...
                           'EdgeColor',cmap(a(iPeak),:),...
                           'LineStyle','none');
    end               

    y =  Ass.Task.FSC{s};        % Assembly peaks
    x = ((1:length(y))*p.bw)/3600; % get binned times into hours

    for iAss= 1%:Ass.Task.patterns.nAss(s)
        plot(x,y(:,iAss),'color',cmap(iAss,:));
        x_= Ass.Task.patterns.pks.time{s}{iAss}/3600*p.bw;
        y_= Ass.Task.patterns.pks.peak{s}{iAss};

        scatter(x_,y_,15,cmap(iAss,:),'filled')
    end
    plot([min(x) max(x)],[Ass.Task.ciHld(s,3) Ass.Task.ciHsc(s,2)], ':k','Linewidth',1.5)
    axis([min(x) max(x) 0 Inf])
    xlabel('Time (hours)')
    ylabel('Assembly ')
    clear a b c cmap i iAss iPeak s t1 t2 x x_ y y_ ymax h
end
%% example plot - one assembly
figure('color','w','name',[p.target ' assembly activation scores']); 
Assem = [1 2 5];
ymax  = 20;
cmap  = {'r','g','b'}; % Pre-sleep,Task,Post-sleep
Face_alpha=0.05;
for s=1:3
%     sp=subplot(3,1,s); hold on
     sp=subaxis(3,1,s,'Spacing',0.03,'Margin',0.1); hold on

    t1 = (length(Ass.Pre.FSCproj{s}(:,Assem(s))))*p.bw/3600;
    t2 = (length(Ass.Pre.FSCproj{s}(:,Assem(s)))+length(Ass.Task.FSC{s}(:,Assem(s))))*p.bw/3600;
    t3 = (length(Ass.Pre.FSCproj{s}(:,Assem(s)))+length(Ass.Task.FSC{s}(:,Assem(s))) + length(Ass.Post.FSCproj{s}(:,Assem(s))))*p.bw/3600;
    
    
    %plot shaded areas
    h=area([0 t1],[ymax ymax],'Linewidth',1.5,...
                   'FaceColor',cmap{1},...
                   'EdgeColor',cmap{1},...
                   'LineStyle','none');
    drawnow; pause(0.0);  % required to make transparent
    h.Face.ColorType = 'truecoloralpha';
    h.Face.ColorData(4) = 255* Face_alpha;
    h=area([t1 t2],[ymax ymax],'Linewidth',1.5,...
                   'FaceColor',cmap{2},...
                   'EdgeColor',cmap{2},...
                   'LineStyle','none');
    drawnow; pause(0.0);  % required to make transparent
    h.Face.ColorType = 'truecoloralpha';
    h.Face.ColorData(4) = 255* Face_alpha;
    h=area([t2 t3],[ymax ymax],'Linewidth',1.5,...
                   'FaceColor',cmap{3},...
                   'EdgeColor',cmap{3},...
                   'LineStyle','none');
    drawnow; pause(0.0);  % required to make transparent
    h.Face.ColorType = 'truecoloralpha';
    h.Face.ColorData(4) = 255* Face_alpha;
    
    
    % plot traces
    y= [Ass.Pre.FSCproj{s}(:,Assem(s))-median(Ass.Pre.FSCproj{s}(:,Assem(s)));...
        Ass.Task.FSC{s}(:,Assem(s))-median(Ass.Task.FSC{s}(:,Assem(s)));...
        Ass.Post.FSCproj{s}(:,Assem(s))-median(Ass.Post.FSCproj{s}(:,Assem(s)))];
    x = ((1:length(y))*p.bw)/3600; % get binned times into hours
    plot(x,y,'k');
    
    axis([0 Inf -1 ymax])
    
    plot([t1 t1], get(gca,'Ylim')-2,':k')
    plot([t2 t2], get(gca,'Ylim')-2,':k')
    
    plot([0 t3],[Ass.Task.ciHsc(s,3) Ass.Task.ciHsc(s,3)],'--r','LineWidth',1.5)
    
      
    switch s
        case 1
           
            title('Assembly activation strength')
            ylabel('mPFC')
            set(gca,'XTick',[])
        case 2
            text(t1/3.5, ymax-2.5,'Pre-task rest','Color','r')
            text(t1+(t2-t1)/2.5, ymax-2.5,'DNMTP task','Color','g')
            text(t2+(t3-t2)/3, ymax-2.5,'Post-task rest','Color','b')
            ylabel('Assembly activation strength')
            text(t2+(t3-t2)*.6, 2+Ass.Task.ciHsc(s,3),'Threshold','Color','r')
            ylabel('HP')
            set(gca,'XTick',[])
        case 3
            xlabel('Time (hours)')
            ylabel('Joint mPFC~HP')
    end
    box off
end
%

subaxis(3,1,1,'Spacing',0.03,'Margin',0.1)

    s=1
    axesPosition = get(gca,'Position');          %# Get the current axes position
    hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
                'Color','none',...           %#   ... with no background color
                'YAxisLocation','right',...  %#   ... located on the right
                'XTick',[],...               %#   ... with no x tick marks
                'YTick',[0 0.5 1],...               %#   ... with no x tick marks
                'Ycolor','b',...
                'TickLength',[0 0],...
                'Box','off');                %#   ... and no surrounding box
%     ylabel(hNewAxes,{'Av. member unit' ;'Firing Rate (Hz)'},'Color','b','FontSize',10);  
    box off
    hold on
    units_ = Ass.Task.patterns.units.memIDs{s}{Assem(s)};
    Mean_Rate_PFC= [nanmean(FR.Pre{s}.avgFR_Pre(units_)) nanmean(FR.Pre{s}.avgFR_Pre(units_)) ...
                    nanmean(FR.Task{s}.avgFR_Task(units_))   nanmean(FR.Task{s}.avgFR_Task(units_)) ...
                    nanmean(FR.Post{s}.avgFR_Post(units_)) nanmean(FR.Post{s}.avgFR_Post(units_))];
         
        plot(hNewAxes,[0 t1 t1 t2 t2 t3],Mean_Rate_PFC)% add Firing p.rates of HP units
         set(hNewAxes,'Ylim',[-2 1.5])


% subplot(3,1,2)
subaxis(3,1,2,'Spacing',0.03,'Margin',0.1)

        s=2
    axesPosition = get(gca,'Position');          %# Get the current axes position
    hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
                'Color','none',...           %#   ... with no background color
                'YAxisLocation','right',...  %#   ... located on the right
                'XTick',[],...               %#   ... with no x tick marks
                'YTick',[0 0.5 1],...               %#   ... with no x tick marks
                'Ycolor','b',...
                'TickLength',[0 0],...
                'Box','off');                %#   ... and no surrounding box
    ylabel(hNewAxes,'Av. member unit Firing Rate (Hz)','Color','b','FontSize',10);  
    box off
    hold on
        units_=Ass.Pre.patterns.units.memIDs{s}{Assem(s)};
        Mean_p.rate_HP = [nanmean(FR.Pre{s}.avgFR_Pre(units_)) nanmean(FR.Pre{s}.avgFR_Pre(units_)) ...
                        nanmean(FR.Task{s}.avgFR_Task(units_))   nanmean(FR.Task{s}.avgFR_Task(units_)) ...
                        nanmean(FR.Post{s}.avgFR_Post(units_)) nanmean(FR.Post{s}.avgFR_Post(units_))];
        plot(hNewAxes,[0 t1 t1 t2 t2 t3],Mean_p.rate_HP)% add Firing p.rates of HP units
         set(hNewAxes,'Ylim',[-2 1.5])
         
         
subaxis(3,1,3,'Spacing',0.03,'Margin',0.1)

        s=3
    axesPosition = get(gca,'Position');          %# Get the current axes position
    hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
                'Color','none',...           %#   ... with no background color
                'YAxisLocation','right',...  %#   ... located on the right
                'XTick',[],...               %#   ... with no x tick marks
                'YTick',[0 0.5 1],...               %#   ... with no x tick marks
                'Ycolor','b',...
                'TickLength',[0 0],...
                'Box','off');                %#   ... and no surrounding box
    
    box off
    hold on
        units_=Ass.Pre.patterns.units.memIDs{s}{Assem(s)};
        Mean_p.rate_Joint = (Mean_p.rate_HP + Mean_Rate_PFC)/2;
        plot(hNewAxes,[0 t1 t1 t2 t2 t3],Mean_p.rate_Joint)% add Firing p.rates of HP units
         set(hNewAxes,'Ylim',[-2 1.5])   
clear a b c cmap i iAss iPeak s t1 t2 x x_ y y_ ymax h Assem axesPosition Face_alpha Mean_Rate_PFC t3 units_ hNewAxes sp Mean_p
%% Plot factor score distributions
figure('color','w','name',[p.target ' assembly activation scores']); 
% bins=-5:0.5:20;
bins=logspace(-1, 1.3,100);
smooth_flag=0;
cmap = {'r','g','b'}; % Pre-sleep,Task,Post-sleep
% cmap = {[0.3 0.3 0.3]; [0.2 1 0.2]; [0.6 0.6 0.2]};  
Assem=[1 1 1];
% Assem=[5 5 5];
for s=1:3
    subplot(3,1,s);hold on
    plot(Ass.Pre.patterns.pks.LogStrengthHist{s}(:,Assem(s))',Ass.Pre.patterns.pks.LogBins,'Linewidth',2,'color',cmap{1})
    plot(Ass.Task.patterns.pks.LogStrengthHist{s}(:,Assem(s))',Ass.Task.patterns.pks.LogBins,'Linewidth',2,'color',cmap{2})
    plot(Ass.Post.patterns.pks.LogStrengthHist{s}(:,Assem(s))',Ass.Post.patterns.pks.LogBins,'Linewidth',2,'color',cmap{3})
%     scatter(Ass.Pre.patterns.pks.ModeLogStrengthPos{s}(Assem(s)),Ass.Pre.patterns.pks.ModeLogStrength{s}(Assem(s)),...
%         'MarkerEdgecolor',cmap{1},'MarkerFacecolor',cmap{1})
%     scatter(Ass.Task.patterns.pks.ModeLogStrengthPos{s}(Assem(s)),Ass.Task.patterns.pks.ModeLogStrength{s}(Assem(s)),...
%         'MarkerEdgecolor',cmap{2},'MarkerFacecolor',cmap{2})
%     scatter(Ass.Post.patterns.pks.ModeLogStrengthPos{s}(Assem(s)),Ass.Post.patterns.pks.ModeLogStrength{s}(Assem(s)),...
%         'MarkerEdgecolor',cmap{3},'MarkerFacecolor',cmap{3})
     
    scatter(Ass.Pre.patterns.pks.LogStrengthQ3binValue{s}(Assem(s)),Ass.Pre.patterns.pks.LogStrengthQ3bin{s}(Assem(s)),...
        'MarkerEdgecolor',cmap{1},'MarkerFacecolor',cmap{1})
    scatter(Ass.Task.patterns.pks.LogStrengthQ3binValue{s}(Assem(s)),Ass.Task.patterns.pks.LogStrengthQ3bin{s}(Assem(s)),...
        'MarkerEdgecolor',cmap{2},'MarkerFacecolor',cmap{2})
    scatter(Ass.Post.patterns.pks.LogStrengthQ3binValue{s}(Assem(s)),Ass.Post.patterns.pks.LogStrengthQ3bin{s}(Assem(s)),...
        'MarkerEdgecolor',cmap{3},'MarkerFacecolor',cmap{3})    
    
%     figure; plot(Ass.Pre.patterns.pks.LogBins,Ass.Pre.patterns.pks.CumLogStrengthHist{s}(:,:))
    % 
    title(p.names{s})
%     set(gca,'YScale','log') 
%     set(gca,'XScale','log') 
    axis([0 1 0 20])
    switch s 
        case 1
        case 2
        ylabel('Log (Factor activation strength)') 
        case 3
        xlabel('Normalised Distribution')
        legend('Pre-Task Sleep','During DNMTP Task','Post-Task Sleep'); legend('boxoff')
    end
    
end
clear cmap bins  s smooth_flag
%% Plot factor score distributions (areas)
figure('color','w','name',[p.target ' assembly activation scores']); 
% bins=-5:0.5:20;
bins=logspace(-1, 1.3,100);
smooth_flag=0;
cmap = {'r','g','b'}; % Pre-sleep,Task,Post-sleep
% cmap = {[0.3 0.3 0.3]; [0.2 1 0.2]; [0.6 0.6 0.2]};  
Assem=[1 1 1];
% Assem=[5 5 5];
Face_alpha=0.5;
ymax = 20;
for s=1:3
    subaxis(3,1,s,'Spacing',0.03,'Margin',0.1); hold on
    box off
    
    %%%% Pre-sleep
    h=area(Ass.Pre.patterns.pks.LogBins,...
           Ass.Pre.patterns.pks.LogStrengthHist{s}(:,Assem(s))',...
           'Linewidth',1.5,...
           'FaceColor',cmap{1},...
           'LineStyle','-');
    drawnow; pause(0);  % required to make transparent
    h.Face.ColorType = 'truecoloralpha';
    h.Face.ColorData(4) = 255* Face_alpha; 
	
    
    %%%% Task
    i=area(Ass.Task.patterns.pks.LogBins,...
           Ass.Task.patterns.pks.LogStrengthHist{s}(:,Assem(s))',...
           'Linewidth',1.5,...                   
           'FaceColor',cmap{2},...
           'EdgeColor',cmap{2},...
           'LineStyle','-');
    drawnow; pause(0);  % required to make transparent
    i.Face.ColorType = 'truecoloralpha';
    i.Face.ColorData(4) = 255* Face_alpha; 
    
    %%%% Post-sleep
    j=area(Ass.Post.patterns.pks.LogBins,...
           Ass.Post.patterns.pks.LogStrengthHist{s}(:,Assem(s))',...
           'Linewidth',1.5,...                    
           'FaceColor',cmap{3},...
           'EdgeColor',cmap{3},...
           'LineStyle','-');
    drawnow; pause(0);  % required to make transparent
    j.Face.ColorType = 'truecoloralpha';
    j.Face.ColorData(4) = 255* Face_alpha; 
    
    % blank below thresholds?
%         area([Ass.Task.patterns.thrsh.ciHsc(s) Ass.Task.patterns.thrsh.ciHsc(s)],[0 ymax],...
%             'Linewidth',1,...
%             'FaceColor','w',...
%             'EdgeColor','w',...
%             'LineStyle','-');
    
    view([90 -90])
    text(17, 0.05,p.names{s},'FontWeight','bold')
    axis([0 ymax 0 1 ])%Remember x and y axes fipped)
    set(gca,'YTick',[0 0.5 1])
    switch s 
        case 1
            title('Log (Assembly activation strength)')
            set(gca,'YTick',[])
            % legend('Pre-Task Rest','During DNMTP Task','Post-Task Rest','Location','northeast','FontSize',6); legend('boxoff')
            text(14,0.5,'Pre-Task Rest','color','r')
            text(12,0.5,'During DNMTP Task','color','g')
            text(10,0.5,'Post-Task Rest','color','b')
        case 2     
            set(gca,'YTick',[])
        case 3
            ylabel('Normalised Distribution')
    end    
end
clear Assem bins Face_alpha s smooth_flag h i j 
%% Plot factor score peak distributions
figure('color','w','name',[p.target ' assembly activation strength']); 

cmap = {'r','g','b'};
for s=1:3
    bins=Ass.Task.patterns.thrsh.ciHsc(s):0.5:20;
    subplot(1,3,s);hold on
%     temp1=histc(cell2mat(Ass.Pre.patterns.peaks.peak{s}'),bins);
%     temp2=histc(cell2mat(Ass.Task.patterns.peaks.peak{s}'),bins);
%     temp3=histc(cell2mat(Ass.Post.patterns.peaks.peak{s}'),bins);
    temp1=histc((Ass.Pre.patterns.pks.peak{s}{1}),bins);
    temp2=histc((Ass.Task.patterns.pks.peak{s}{1}),bins);
    temp3=histc((Ass.Post.patterns.pks.peak{s}{1}),bins);
    plot(temp1./length(temp1),bins,'Linewidth',2,'color',cmap{1})
    plot(temp2./length(temp2),bins,'Linewidth',2,'color',cmap{2})
    plot(temp3./length(temp3),bins,'Linewidth',2,'color',cmap{3})
    plot([0 1.5],[Ass.Task.patterns.thrsh.ciHsc(s) Ass.Task.patterns.thrsh.ciHsc(s)],':r','LineWidth',2)
    title(p.names{s})
%     set(gca,'YScale','log')       
%     axis([0 10 0 20])
end
clear cmap bins s temp1 temp2 temp3 ymax 
%% Plot assembly repeat timing histograms
figure('color','w','name',[p.target ' assembly activation p.rate']); 
cmap = {'r','g','b'};
for s=1:3
    subplot(1,3,s);hold on
    plot((Ass.Pre.patterns.pks.repHistBins*Ass.Pre.patterns.bw).^-1,   (mean(Ass.Pre.patterns.pks.repHist{s})),'color',cmap{1},'LineWidth',2)
    plot((Ass.Task.patterns.pks.repHistBins*Ass.Task.patterns.bw).^-1, (mean(Ass.Task.patterns.pks.repHist{s})),'color',cmap{2},'LineWidth',2)
    plot((Ass.Post.patterns.pks.repHistBins*Ass.Post.patterns.bw).^-1, (mean(Ass.Post.patterns.pks.repHist{s})),'color',cmap{3},'LineWidth',2)
%     temp1=cell2mat(Ass.Pre.patterns.peaks.peak{s}');  temp1=histc(temp1,Ass.Pre.patterns.peaks.IEI_hist_bins*Ass.Pre.patterns.bw);
%     temp2=cell2mat(Ass.Task.patterns.peaks.peak{s}'); temp2=histc(temp2,Ass.Pre.patterns.peaks.IEI_hist_bins*Ass.Pre.patterns.bw);
%     temp3=cell2mat(Ass.Post.patterns.peaks.peak{s}'); temp3=histc(temp2,Ass.Pre.patterns.peaks.IEI_hist_bins*Ass.Pre.patterns.bw);
%     plot((Ass.Pre.patterns.peaks.IEI_hist_bins*Ass.Pre.patterns.bw).^-1,   temp1,'color',cmap{1},'LineWidth',2)
%     plot((Ass.Task.patterns.peaks.IEI_hist_bins*Ass.Task.patterns.bw).^-1, temp2,'color',cmap{2},'LineWidth',2)
%     plot((Ass.Post.patterns.peaks.IEI_hist_bins*Ass.Post.patterns.bw).^-1, temp3,'color',cmap{3},'LineWidth',2)
    set(gca,'XScale','log')       
%     axis([-Inf Inf 0 .05])
    title(p.names{s})
    xlabel('Average Assembly activation rate (Hz)')
    ylabel('Distribution')
    legend('Pre-sleep','Task','Post-sleep'); legend('boxoff')
end
clear cmap s
%% Plot pattern transition probability
try
    load('whiteHot');
catch
    cmap=(hot(64));
end
for s = 1:length(Ass.Task.patterns.FSC)      % across all areas/assem types
    figure('Name',[Ass.Task.patterns.filename ' : ' Ass.Task.patterns.titles{s} ' Assemblies'])
    subplot(2,1,1) % Plot assembly sequence
    stairs(Ass.Task.patterns.ptn.seq{s}(:,1)*Ass.Task.patterns.bw,Ass.Task.patterns.ptn.seq{s}(:,2),'k')
    hold on; axis tight
    xlabel('Time (s)')
    ylabel('Assembly no.')
    set(gca, 'ylim',[1 Ass.Task.patterns.nAss(s)],...
        'ytick',[1:Ass.Task.patterns.nAss(s)],...
        'yticklabel',Ass.Task.patterns.AssLabels{s});
    
    subplot(2,2,3); hold on
    temp=Ass.Task.patterns.ptn.transMatrix{s}; temp_hist=histc(temp(1:end),0:0.05:1);
    plot([Ass.Task.patterns.ptn.ciHTransition(s) Ass.Task.patterns.ptn.ciHTransition(s)],[0 max(temp_hist./sum(temp_hist))],'r','LineWidth',1.5)
    plot(0:0.05:1,temp_hist./sum(temp_hist),'k','LineWidth',1.5)
    legend('95th percentile of shuffled p.patterns') , legend('boxoff')
    xlabel('Transition probability')
    ylabel('Distribution')
    subplot(2,2,4)
    temp(temp<Ass.Task.patterns.ptn.ciHTransition(s))=0;
    pcolor(flipud(rot90(temp)))%; set(gca,'YDir','normal')
    set(gca, 'ylim', [1 Ass.Task.patterns.nAss(s)],...
        'xtick',[1:Ass.Task.patterns.nAss(s)],'xticklabel',Ass.Task.patterns.AssLabels{s},...
        'ytick',[1:Ass.Task.patterns.nAss(s)],'yticklabel',Ass.Task.patterns.AssLabels{s},...
        'clim',[0 1])
    
    colormap(cmap)
    colorbar
    axis square
    title('Assembly sequence transition matrix')
    xlabel('Event (t)')
    ylabel('Event (t+1)')
    % e.g. p(Assem. 3 -> Ass.Task.patterns. 2) ...  Ass.Task.patterns.p.pattern.trans_matrix{s}(3,2)
end
clear cmap s temp temp_hist 
%% Scatter peak vs p.rate
figure('color','w','name',[p.target ' assembly activation rate vs strength']); 
cmap = {'r','g','b'};
for s=1:3
    subplot(1,3,s);hold on
    try
    scatter((Ass.Pre.patterns.pks.repMean{s}).^-1,cellfun(@mean,Ass.Pre.patterns.pks.peak{s}),...
            'MarkerEdgecolor',cmap{1},'MarkerFacecolor',cmap{1})
    scatter((Ass.Task.patterns.pks.repMean{s}).^-1,cellfun(@mean,Ass.Task.patterns.pks.peak{s}),...
            'MarkerEdgecolor',cmap{2},'MarkerFacecolor',cmap{2})
    scatter((Ass.Post.patterns.pks.repMean{s}).^-1,cellfun(@mean,Ass.Post.patterns.pks.peak{s}),...
            'MarkerEdgecolor',cmap{3},'MarkerFacecolor',cmap{3})   
    catch
    end
end
% % Plot change in 
clear cmap s 
%% Save Graph layouts
if flags.makeGraph
% Make labels and category numbers
for r_id=1:sum(Ass.Task.patterns.nAss)
    if r_id<=Ass.Task.patterns.nAss(1)
        %r_labels{r_id}=['PFC' num2str(r_id)];
        r_labels{r_id}=[num2str(r_id)];
        r_partition(r_id)=1;
    elseif r_id<=Ass.Task.patterns.nAss(1)+ Ass.Task.patterns.nAss(2)
        %r_labels{r_id}=['HP' num2str(r_id-Ass.Task.patterns.nAss(1))];
        r_labels{r_id}=[num2str(r_id-Ass.Task.patterns.nAss(1))];
        r_partition(r_id)=2;
    else
       %r_labels{r_id}=['Joint' num2str(r_id-(Ass.Task.patterns.nAss(1)+Ass.Task.patterns.nAss(2)))];
       r_labels{r_id}=[num2str(r_id-(Ass.Task.patterns.nAss(1)+Ass.Task.patterns.nAss(2)))];
       r_partition(r_id)=3;
    end
    % r_shapes{r_id}='ellipse';
end
% r = (r - min(r(1:end))) ./ ( max(r(1:end)) - min(r(1:end)) ); % scale between 0~1 (now using mat2gray for this)

r_Pre = Ass.Pre.patterns.ptn.transMatrix{4}; r_Pre(r_Pre<Ass.Pre.patterns.ptn.ciHTransition(4)) = 0; % plot significant transitions only
[~]=adj2gephilab2([p.rat '_Pre_SequenceGraph'],mat2gray(r_Pre),r_partition,r_labels);

r_Task = Ass.Task.patterns.ptn.transMatrix{4}; r_Task(r_Task<Ass.Task.patterns.ptn.ciHTransition(4)) = 0; % plot significant transitions only
[~]=adj2gephilab2([p.rat '_Task_SequenceGraph'],mat2gray(r_Task),r_partition,r_labels);

r_Post = Ass.Post.patterns.ptn.transMatrix{4}; r_Post(r_Post<Ass.Post.patterns.ptn.ciHTransition(4)) = 0; % plot significant transitions only
[~]=adj2gephilab2([p.rat '_Post_SequenceGraph'],mat2gray(r_Post),r_partition,r_labels);
change1 = r_Task - r_Pre; change1(~isfinite(change1))=0; change1=mat2gray(change1);
[~]=adj2gephilab2([p.rat '_Pre_vs_Task_SequenceGraph'], change1,r_partition,r_labels);
[~]=adj2gephilab2([p.rat '_Post_vs_Task_SequenceGraph'],r_Post - r_Task,r_partition,r_labels);
end
clear r_labels change1 r_id r_partition r_Post r_Pre r_Task 
%% Task performance vs. reactivation p.rate

% Load task decoding
temp = load([p.pat p.rat '_Task_decoding.mat'],'D');
Decoding = temp.D; clear temp
cmap = {'r','g','b'};
clear x y
%  (1) Scatter peak decoding score vs average p.rate
figure('color','w','name',[p.target ' assembly activation p.rate vs decoding strength']);
for s=1:3
    subplot(3,1,s);hold on
    title(strcat(p.names(s), ' Assemblies'))
    try
        x = max(Decoding.Assem.TS{s});
        y.Pre  = (Ass.Pre.patterns.pks.repMedian{s}).^-1;         fPre = ezfit(x,y.Pre,'a*x+b');
        y.Task = (Ass.Task.patterns.pks.repMedian{s}).^-1;        fTask = ezfit(x,y.Task,'a*x+b');
        y.Post = (Ass.Post.patterns.pks.repMedian{s}).^-1;        fPost = ezfit(x,y.Post,'a*x+b');
        
        plot(x,fPre.m(1)*x+fPre.m(2),'LineWidth',1.5, 'color',cmap{1})
        plot(x,fAss.Task.m(1)*x+fAss.Task.m(2),'LineWidth',1.5, 'color',cmap{2})
        plot(x,fPost.m(1)*x+fPost.m(2),'LineWidth',1.5, 'color',cmap{3})

        
        scatter(x,y.Pre,'MarkerEdgecolor',cmap{1},'MarkerFacecolor',cmap{1})
        scatter(x,y.Task,'MarkerEdgecolor',cmap{2},'MarkerFacecolor',cmap{2})
        scatter(x,y.Post,'MarkerEdgecolor',cmap{3},'MarkerFacecolor',cmap{3})

        text(6,5e-3,strcat('R^2=', sprintf('%3.1g',fPre.r^2)),'color',cmap{1})
        text(6,9e-3,strcat('R^2=', sprintf('%3.1g',fAss.Task.r^2)),'color',cmap{2})
        text(6,13e-3,strcat('R^2=', sprintf('%3.1g',fPost.r^2)),'color',cmap{3})
    catch
    end
%     axis([0 10 0 0.03])
%     set(gca,'yscale','log')
    switch s
        case 1
            legend('Pre-sleep','Task','Post-sleep'); legend('boxoff')
        case 2
            ylabel('Average activation p.rate (Hz)')            
        case 3
            xlabel('Task outcome prediction power (F-score)')
    end
end
clear cmap s x fPre fPost fTask lastfit y
%% Task performance vs. reactivation strength
figure('color','w','name',[p.target ' assembly activation strength vs decoding strength']);
cmap = {'r','g','b'};

for s=1:3
    subplot(3,1,s);hold on
    title(strcat(p.names(s), ' Assemblies'))
    try
        x = max(Decoding.Assem.TS{s});
        y.Pre  = Ass.Pre.patterns.pks.LogStrengthQ3bin{s};         fPre = ezfit(x,y.Pre,'a*x+b');
        y.Task = Ass.Task.patterns.pks.LogStrengthQ3bin{s};        fTask = ezfit(x,y.Task,'a*x+b');
        y.Post = Ass.Post.patterns.pks.LogStrengthQ3bin{s};        fPost = ezfit(x,y.Post,'a*x+b');
        
        plot(x,fPre.m(1)*x+fPre.m(2),'LineWidth',1.5, 'color',cmap{1})
        plot(x,fAss.Task.m(1)*x+fAss.Task.m(2),'LineWidth',1.5, 'color',cmap{2})
        plot(x,fPost.m(1)*x+fPost.m(2),'LineWidth',1.5, 'color',cmap{3})
        
        scatter(x,y.Pre,'MarkerEdgecolor',cmap{1},'MarkerFacecolor',cmap{1})
        scatter(x,y.Task,'MarkerEdgecolor',cmap{2},'MarkerFacecolor',cmap{2})
        scatter(x,y.Post,'MarkerEdgecolor',cmap{3},'MarkerFacecolor',cmap{3})

        text(6,1,strcat('R^2=', sprintf('%3.1g',fPre.r^2)),'color',cmap{1})
        text(6,2,strcat('R^2=', sprintf('%3.1g',fAss.Task.r^2)),'color',cmap{2})
        text(6,3,strcat('R^2=', sprintf('%3.1g',fPost.r^2)),'color',cmap{3})
    catch
    end
    axis([0 10 0 8])
%     set(gca,'yscale','log')
    switch s
        case 1
            legend('Pre-sleep','Task','Post-sleep'); legend('boxoff')
        case 2
            ylabel('Average activation Strength (A.U.)')            
        case 3
            xlabel('Task outcome prediction power (F-score)')
    end
end


%      fnam=([p.pat p.sep p.rat '_reactivation']); % your file name
%      snam='AD1';
%      s=hgexport('readstyle',snam);
%      hgexport(gcf,fnam,s);
clear s x fPost fPre fTask lastfit y 
%% Pre-task vs Post-task activation strength, vs decoding power 
figure('color','w','name',[p.target ' Pre vs Post strength vs decoding strength']);
axlims = [0 8];
clims =[0 5];
for s=1:3
    subplot(1,3,s);hold on

%     title(strcat(p.names(s)))
    text(axlims(2)*0.05,axlims(2)*0.95,p.names(s))
    plot([axlims],[axlims],':k')
        x = Ass.Pre.patterns.pks.LogStrengthQ3bin{s};                 
        y = Ass.Post.patterns.pks.LogStrengthQ3bin{s};         
        z = max(Decoding.Assem.TS{s});
        try
%             f = ezfit(x,y,'a*x+b');
%             plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color',[0.6 0.6 0.6])
%             text(5,0.7,strcat('R^2=', sprintf('%3.1g',f.r^2)),'color',[0.6 0.6 0.6])

            b0=[1,0];
            fit0=@(b,x)(b(1)*x+b(2));
            [f,R,~,~,MSE] = nlinfit(x,y,fit0,b0);% MSE=nanmean(R.^2);
            R2= 1 -nansum(R.^2)/nansum((y-nanmean(y)).^2);
            % or: %C=corrcoef(y,fit0(f,x)); R2=C(1,2)^2 ... but fails with NaN;
            plot(x,fit0(f,x),'LineWidth',1.5, 'color',[0.6 0.6 0.6])
            % text(axlims(2)*0.6,axlims(2)*0.1,strcat('MSE=', sprintf('%3.1g',MSE)),'color',[0.6 0.6 0.6])
            text(axlims(2)*0.6,axlims(2)*0.1,strcat('R^2=', sprintf('%3.1g',R2)),'color',[0.6 0.6 0.6])
        catch
        end
        scatter(x,y,50,z,'filled')
        
        xlabel('Pre-task strength')% activation strength (A.U.)')  
        ylabel('Post-task strength')% activation strength (A.U.)')                    
        switch s
        case 2
            title('Offline assembly activation vs. task outcome decoding score')
        end
        axis([axlims axlims]); axis square; box on
 end
 h=colorbar('EastOutside'); ylabel(h,{'Task-predictive power';'(F-score)'},'fontsize',12,'fontname','Arial');
    set(h,'Position',[0.91,0.4,0.01,0.2300])
    caxis(clims)
%     colormap jet

clear axlims b0 clims f MSE R R2 s x y z fit0 h 
%% Pre-task vs Post-task activation p.rate, vs decoding power 
figure('color','w','name',[p.target ' Pre vs Post p.rate vs decoding strength']);
axlims = [0 0.03    ];
clims =[0 5];
for s=1:3
    subplot(1,3,s);hold on
%     title(strcat(p.names(s)))
    text(axlims(2)*0.05,axlims(2)*0.95,p.names(s))
    plot([axlims],[axlims],':k')
    x = (Ass.Pre.patterns.pks.repMedian{s}).^-1;
    y = (Ass.Post.patterns.pks.repMedian{s}).^-1;
    z = max(Decoding.Assem.TS{s});
    try
%         f = ezfit(x,y,'a*x+b');
%         plot(x,f.m(1)*x+f.m(2),'LineWidth',1.5, 'color',[0.6 0.6 0.6])
%         text(axlims(2)*0.6,axlims(2)*0.1,strcat('R^2=', sprintf('%3.1g',f.r^2)),'color',[0.6 0.6 0.6])
        
        b0=[1,0];
        fit0=@(b,x)(b(1)*x+b(2));
        [f,R,~,~,MSE] = nlinfit(x,y,fit0,b0);% MSE=nanmean(R.^2);
        R2= 1 -nansum(R.^2)/nansum((y-nanmean(y)).^2);
        % or: %C=corrcoef(y,fit0(f,x)); R2=C(1,2)^2 ... but fails with NaN;
        plot(x,fit0(f,x),'LineWidth',1.5, 'color',[0.6 0.6 0.6])
        % text(axlims(2)*0.6,axlims(2)*0.1,strcat('MSE=', sprintf('%3.1g',MSE)),'color',[0.6 0.6 0.6])
        text(axlims(2)*0.6,axlims(2)*0.1,strcat('R^2=', sprintf('%3.1g',R2)),'color',[0.6 0.6 0.6])

    catch
    end
    scatter(x,y,50,z,'filled')
    ylabel('Post-task Rate (Hz)')
    xlabel('Pre-task Rate (Hz)')
    switch s
        case 2
        title('Offline assembly activation vs. task outcome decoding score')
    end
    axis([axlims axlims]); axis square; box on
end

    h = colorbar('EastOutside'); ylabel(h,{'Task-predictive power';'(F-score)'},'fontsize',12,'fontname','Arial');
    set(h,'Position',[0.91,0.4,0.01,0.2300])
    caxis(clims)
%     colormap jet
clear axlims b0 clims f MSE R R2 s x y z fit0 h cmap
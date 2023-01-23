clear 
cd('C:\Analysis\AssemblyAnalysis\LFPanalysis\Aleks chosen LFP segments\CellsOnTetrodes')
load('cellsonTetrodes.mat')
pat = 'C:\Analysis\AssemblyAnalysis\raw';
cd(pat);
useWholeTaskPeriod = true;
restrictToLeverTime = true;
Targets= {'SHORT','MEDIUM','LONG'};
iTarget = 3;

f.bandNames    = {'4-6Hz','theta', 'Spindle', 'low gamma', 'high gamma'};
f.flim = [4 6;...
    6 10; ...
    10 18; ...
    30 45;...
    55 80];
f.fcolors= [0.2 0.2 0.2 ;...
    0.2 0.3 0.4 ; ...
    0.1 0.1 0.8; ...
    0.8 0.4 0.23;...
    0.6 0.3 0.7];


FS = 508;%tb = (1:nSamps)/FS;
% D.lim_lfp    = round([-.1,0.1]*FS);
% D.ThetaBins   = linspace(0, 2*pi, 36);
% D.tb_lfp = (D.lim_lfp(1):D.lim_lfp(2))/FS;

addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\jolato\phasemeth')
% addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\jolato\nlx')

warning('off')
Areas = {'PFC','HP','joint'};

UseResampledResults = true;
if UseResampledResults 
    eval(sprintf('MixedSelectivityfileList_%s = dir([pat ''\\MixedSelectivity\\TrialNosMatched\\*%s*PFC*.mat'']);',Targets{iTarget},Targets{iTarget}))
else
    eval(sprintf('MixedSelectivityfileList_%s = dir([pat ''\\MixedSelectivity\\*%s*PFC*.mat'']);',Targets{iTarget},Targets{iTarget}) )
end
for iFile =1:length(MixedSelectivityfileList_LONG)
    MixedSelectivityfileList_LONG(iFile).name = strtok(MixedSelectivityfileList_LONG(iFile).name,'_');
end

Filelist = intersect({MixedSelectivityfileList_LONG.name},tetrodes.names);
%% How many spikes within task period to bootstrap amongst?
    noDrawnSpikes = 50;

for iDay = 1:length(Filelist)
    disp(sprintf('File %d of %d (%s)',iDay,length(Filelist),Filelist{iDay}))
    % Load cells vs tetrode assignments
    fn = fullfile(pat,sprintf('%s.mat',Filelist{iDay}));
    spikes = load(fn);
    for iArea = 1:2
        eval(sprintf('nU = length(spikes.%scells);',Areas{iArea}))
%          t_ = sort(mean([spikes.trangeleft_choice;spikes.trangeleft_sample;spikes.trangeright_sample;spikes.trangeright_choice],2))*1e-6;
         t_ = sort(mean([spikes.trangeleft_choice;spikes.trangeleft_sample;spikes.trangeright_sample;spikes.trangeright_choice;...
             spikes.ERRORtrangeleft_choice;spikes.ERRORtrangeleft_sample;spikes.ERRORtrangeright_choice;spikes.ERRORtrangeright_sample],2))*1e-6;
         for iUnit =1:nU
             eval(sprintf('st = spikes.%scells{iUnit}.t*1e-4;' ,Areas{iArea}))
             if restrictToLeverTime
                 idx = zeros(size(st));
                 for i=1:length(t_)
                     idx = idx |  st>(t_(i)-4)  & st<(t_(i)+4);
                 end
                 st = st(idx);
             else
                 st(st<=(min(t_)-60) | st>=(max(t_)+60)) = [];
             end
             noSpikes{iArea}{iDay}(iUnit,1)= length(st);
         end
    end
end
figure; hold on ;bins = 0:50:10000;plot(bins,cumsum(histc(cell2mat([noSpikes{1:2}]'),bins)));xlabel('No. Spikes in Task Period');ylabel('Unit count')
plot([noDrawnSpikes noDrawnSpikes],[0 500],':r')
%% Batch process: Compute phase-locking and firing phase preferences
bins = -180:45:180; % for histograms

for iDay = 1:length(Filelist)
    disp(sprintf('File %d of %d (%s)',iDay,length(Filelist),Filelist{iDay}))
    %% Load data
    
    % Load cells vs tetrode assignments
    fn = fullfile(pat,sprintf('%s.mat',Filelist{iDay}));
    spikes = load(fn);
    % Load LFP phase
    fn = fullfile(pat,'LFPphase',sprintf('%s_phase.mat',Filelist{iDay}));
    load(fn)
    %f.phaseData{3} = [f.phaseData{1}(B.usel_out{1});...
    %                  f.phaseData{2}(B.usel_out{2})];
    f.phaseData{3} = [f.phaseData{1};f.phaseData{2}];
    
    
    
    
    % Load mixed selectivity assignments and F-scores
    for iArea = 1:2
        % fidx = find(~cellfun(@isempty,strfind({MixedSelectivityfileList_LONG.name},tetrodes.names{iDay})))
        if UseResampledResults
            fn = [pat '\MixedSelectivity\TrialNosMatched\' Filelist{iDay} '_' Areas{iArea} '_MixedSelectivity_Units.mat'];
        else
            fn = [pat '\MixedSelectivity\' Filelist{iDay} '_' Areas{iArea} '_MixedSelectivity_Units.mat'];
        end
        MixedSelectivity{iArea} = load(fn);
        selectivity.F{iArea}{iDay} = MixedSelectivity{iArea}.D.Long.ANOVA2.F_cor;
        selectivity.p{iArea}{iDay} = MixedSelectivity{iArea}.D.Long.ANOVA2.p_cor;
    end
    %% batch extract phase pref parameters
    for s = 1:2
        for iBand = 1:numel(f.bandNames)
            for iUnit =1:length(f.phaseData{s})
                eval(sprintf('st = spikes.%scells{iUnit}.t*1e-4;' ,Areas{s}))
%                 t_ = mean([spikes.trangeleft_choice;spikes.trangeleft_sample;spikes.trangeright_sample;spikes.trangeright_choice;...
%                     spikes.ERRORtrangeleft_choice;spikes.ERRORtrangeleft_sample;spikes.ERRORtrangeright_choice;spikes.ERRORtrangeright_sample],2)*1e-6;
                t_ = sort(mean([spikes.trangeleft_choice;spikes.trangeleft_sample;spikes.trangeright_sample;spikes.trangeright_choice],2))*1e-6;

                % get spike phase for this band, compute phase locking
                Theta_   = f.phaseData{s}{iUnit}(iBand,:);
                
                % Choose spike time limits
                if restrictToLeverTime
                    idx = zeros(size(st));
                    for i=1:length(t_)
                        idx = idx |  st>=(t_(i)-4) & st<=(t_(i)+4);
                    end
                    Theta_ = Theta_(idx);
                    clear idx i
                else  % restrict spikes to task period only
                    Theta_(st<=(min(t_)-60) | st>=(max(t_)+60) )=[];
                end
                
                clear st t_
                
                fprintf('%s unit %d of %d...\n',Areas{s},iUnit,length(f.phaseData{s}))
                if length(Theta_)>noDrawnSpikes
                    parfor iDraw = 1:100
                        
                        Theta_Draw = Theta_(randsample(length(Theta_),noDrawnSpikes));
                        
                        [Theta__(iDraw,1), ...
                            kappa__(iDraw,1), ...
                            Rbar__(iDraw,1)]   = circStats(Theta_Draw);
                        Phase__(iDraw,1)       = mode(round(rad2deg(Theta_Draw)));
                        [Rp__(iDraw,1),   ...
                            Rbar2__(iDraw,1)]  = rayleigh(Theta_Draw');
                    end
                    Theta(iUnit,1)  = circ_mean(Theta__);
                    Phase(iUnit,1)  = rad2deg(circ_mean(deg2rad(Phase__+180)))-180;
                    kappa(iUnit,1)  = nanmean(kappa__);
                    Rp_(iUnit,1)    = nanmean(Rp__);
                    Rbar_(iUnit,1)  = nanmean(Rbar__);
                    Rbar2_(iUnit,1) = nanmean(Rbar2__);
                else
                    Theta(iUnit,1) = NaN;
                    Phase(iUnit,1) = NaN;
                    kappa(iUnit,1) = NaN;
                    Rp_(iUnit,1)   = NaN;
                    Rbar_(iUnit,1) = NaN;
                    Rbar2_(iUnit,1) = NaN;
                end
                
                
                
                
                %                 [Theta(iUnit) kappa(iUnit), Rbar_(iUnit)]   = circStats(Theta_Draw);
                %
                %
                %                 % Get peak phase preference
                %                 Phase(iUnit)   = mode(round(rad2deg(Theta_)));
                %
                %                 % Rayleigh p value and mean resultant length vector
                %                 [Rp_(iUnit),   Rbar(iUnit)] = rayleigh(Theta_');
                %
                %                 % Rayleigh test only good for n<10, enforce here
                %                 if numel(Theta_)<10, Rp_(iUnit) = 1;end
                
            end
            
            phasestats.Theta{s}{iBand}{iDay} = Theta;
            phasestats.Phase{s}{iBand}{iDay} = Phase;
            phasestats.kappa{s}{iBand}{iDay} = kappa;
            phasestats.Rbar_{s}{iBand}{iDay} = Rbar_;
            phasestats.Rbar2_{s}{iBand}{iDay} = Rbar2_;
            phasestats.Rp_{s}{iBand}{iDay} = Rp_;
        end
        clear Theta Phase kappa Rbar_ Rp_
    end
    
end
%% Get Assembly membership
fname = 'C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\LONGTaskonly\Assemmembertypes.mat';
load(fname,'D')
%% Scatter phase locking on polar axes - mixed selectivity

factor_ = 1;

for s=1:2
    figure('name',[Areas{s},' phase locking']);
    for iBand = 1:numel(f.bandNames)
        ax = subplot(1,numel(f.bandNames),iBand,polaraxes); hold on
        
        % Plot non-sig encoding
        for iDay = 1:length(Filelist)
            Phase = phasestats.Phase{s}{iBand}{iDay}';
            Theta = phasestats.Theta{s}{iBand}{iDay}';
            kappa = phasestats.kappa{s}{iBand}{iDay}';
            Rbar_ = phasestats.Rbar_{s}{iBand}{iDay}';
            Rbar2_ = phasestats.Rbar2_{s}{iBand}{iDay}';
            Rp_ = phasestats.Rp_{s}{iBand}{iDay}';
            p_    = selectivity.p{s}{iDay};
            F_    = selectivity.F{s}{iDay};
            
            if ~sum(isnan(F_)>0)
                idx = p_(:,factor_)>0.05;%;
                %                 idx = idx & Rp_'<0.05;
                pS = polarscatter(ax,...
                    (Theta(idx)),...
                    (Rbar2_(idx)),...
                    30,...
                    [1 1 1]);
                pS.MarkerEdgeColor = 'k';
                pS.MarkerEdgeAlpha = 0.5;
                pS.MarkerFaceAlpha = 0.5;
            end
        end
        % Plot Sig encoding
        for iDay = 1:length(Filelist)
            Phase = phasestats.Phase{s}{iBand}{iDay}';
            Theta = phasestats.Theta{s}{iBand}{iDay}';
            kappa = phasestats.kappa{s}{iBand}{iDay}';
            Rbar_ = phasestats.Rbar_{s}{iBand}{iDay}';
            Rbar2_ = phasestats.Rbar2_{s}{iBand}{iDay}';
            Rp_ = phasestats.Rp_{s}{iBand}{iDay}';
            p_    = selectivity.p{s}{iDay};
            F_    = selectivity.F{s}{iDay};
            
            if ~sum(isnan(F_)>0)
                idx = p_(:,factor_)<0.05;
                %                 idx = idx & Rp_'<0.05;
                pS = polarscatter(ax,...
                    (Theta(idx)),...
                    (Rbar2_(idx)),...
                    30,...
                    F_(idx,factor_),'filled');
                pS.MarkerEdgeAlpha = 0.5;
                pS.MarkerFaceAlpha = 0.9;
            end
        end
        
        %         ax.RLim = [0 0.2];
        title(f.bandNames{iBand})
        colormap(jet)
        caxis([ 0 10])
    end
    
end
%% group averaged phase locking by sig Y/N
for factor_ = 1:3
factors_ = {'Context';'Location';'ContextXLocation'};
for s=1:2
    figure('name',[Areas{s},' phase locking (' factors_{factor_} ')']);
    for iBand = 1:numel(f.bandNames)
        ax = subplot(1,numel(f.bandNames),iBand); hold on
        
        kappa=[];Rbar_=[]; Rp_=[]; p_=[]; F_=[];Rbar2_=[];
         for iDay = 1:length(Filelist)
            kappa = [kappa; phasestats.kappa{s}{iBand}{iDay}];
            Rbar_ = [Rbar_; phasestats.Rbar_{s}{iBand}{iDay}];
            Rbar2_ = [Rbar2_; phasestats.Rbar2_{s}{iBand}{iDay}];
            Rp_    = [Rp_; phasestats.Rp_{s}{iBand}{iDay}];
            p_    = [p_; selectivity.p{s}{iDay}];
            F_    = [F_; selectivity.F{s}{iDay}];
         end
          idx = p_(:,factor_)<0.05;
          %;
%          bar(1,nanmean(kappa(idx))); 
%          bar(2,nanmean(kappa(~idx)))
%          errorbar([1,2],[nanmean(kappa(idx)),nanmean(kappa(~idx))],[nansem(kappa(idx)),nansem(kappa(~idx))])
           
                     
         
           x = Rbar2_(idx); x(isnan(x))=[];
           scatter(0.2*randn(length(x),1)+1,x,'ob','filled','Markerfacealpha',0.5)
           
           x = Rbar2_(~idx); x(isnan(x))=[];
           scatter(0.2*randn(length(x),1)+2,x,'or','filled','Markerfacealpha',0.5) 
         errorbar(1,nanmean(Rbar2_(idx)),nansem(Rbar2_(idx)),'Color','k','LineWidth',1.5)
         errorbar(2,nanmean(Rbar2_(~idx)),nansem(Rbar2_(~idx)),'Color','k','LineWidth',1.5)
         bar(1,nanmean(Rbar2_(idx)),'FaceColor','none','EdgeColor','k','LineWidth',1.5); 
         bar(2,nanmean(Rbar2_(~idx)),'FaceColor','none','EdgeColor','k','LineWidth',1.5)
         title(f.bandNames{iBand})
         set(gca,'Xtick',[1 2],'XtickLabel',{'Encoders','Non-Encoders'},'XTickLabelRotation',-45)
         if iBand==1
             ylabel('Mean resultant vector length')
         end
axis([0 3 0 0.3])
    end
end
end
%% group averaged phase locking by sig Y/N - cum histograms 
exclusive_ = true;
bins = 0:0.005:0.3;
for factor_ = 3;% 1:3
factors_ = {'Context';'Location';'ContextXLocation'};
for s=1:2
    figure('name',[Areas{s},' phase locking (' factors_{factor_} ')']);
    for iBand = 1:numel(f.bandNames)
        ax = subplot(1,numel(f.bandNames),iBand); hold on
        
        kappa=[];Rbar_=[]; Rp_=[]; p_=[]; F_=[];Rbar2_=[];
         for iDay = 1:length(Filelist)
            kappa = [kappa; phasestats.kappa{s}{iBand}{iDay}];
            Rbar_ = [Rbar_; phasestats.Rbar_{s}{iBand}{iDay}];
            Rbar2_ = [Rbar2_; phasestats.Rbar2_{s}{iBand}{iDay}];
            Rp_    = [Rp_; phasestats.Rp_{s}{iBand}{iDay}];
            p_    = [p_; selectivity.p{s}{iDay}];
            F_    = [F_; selectivity.F{s}{iDay}];
         end
         
         if exclusive_
             idx =  sum(p_<0.05 == [0 0 1],2)==3; % NMS for Context * Spatial location only
         else
             idx =  p_(:,3)<0.05==1; %NMS for Context * Spatial location
         end
                        
          
%          bar(1,nanmean(kappa(idx))); 
%          bar(2,nanmean(kappa(~idx)))
%          errorbar([1,2],[nanmean(kappa(idx)),nanmean(kappa(~idx))],[nansem(kappa(idx)),nansem(kappa(~idx))])
           
                     
         
           x = Rbar2_(idx); x(isnan(x))=[];
           x = histc(x,bins);x = x./sum(x);
           x= cumsum(x);
           x = plot(bins,cumsum(x),'b','LineWidth',1.5);
           x = Rbar2_(~idx); x(isnan(x))=[];
           x = histc(x,bins);x = x./sum(x);
           x= cumsum(x);
           x = plot(bins,cumsum(x),':k','LineWidth',1.5);
          
         title(f.bandNames{iBand})
         xlabel('R_b_a_r')
         if iBand==1
            ylabel('Fraction of units')
         end
    axis([0 0.2 0 1])
    end
    legend({'M.S.','No M.S.'}); legend boxoff

end
end
%% group averaged phase locking by sig Y/N
exclusive_ = true;
factors_ = {'none';'Context';'Location';'Linear ContextXLocation';'Nonlinear ContextXLocation'};

for s=1:2
    figure('name',[Areas{s},' phase locking']);
    
    
    for iBand = 1:numel(f.bandNames)
        for factor_ = 1:length(factors_)
        ax = subplot(1,numel(f.bandNames),iBand); hold on
        
        kappa=[];Rbar_=[]; Rp_=[]; p_=[]; F_=[];Rbar2_=[];
         for iDay = 1:length(Filelist)
            kappa = [kappa; phasestats.kappa{s}{iBand}{iDay}];
            Rbar_ = [Rbar_; phasestats.Rbar_{s}{iBand}{iDay}];
            Rbar2_ = [Rbar2_; phasestats.Rbar2_{s}{iBand}{iDay}];
            Rp_    = [Rp_; phasestats.Rp_{s}{iBand}{iDay}];
            p_    = [p_; selectivity.p{s}{iDay}];
            F_    = [F_; selectivity.F{s}{iDay}];
         end
            switch factor_
                    case 1
                            idx =  sum(p_<0.05,2) == 0; % Untuned
                    case 2
                        if exclusive_
                            idx =  sum(p_<0.05 == [1 0 0],2)==3; % Simple context only
                        else
                            idx =  p_(:,1)<0.05; % Simple, context
                        end
                    case 3
                        if exclusive_
                            idx =  sum(p_<0.05 == [0 1 0],2)==3; % Simple spatial location only
                        else
                            idx =  p_(:,2)<0.05; % Simple spatial
                        end
                    case 4
                        if exclusive_
                            idx =  sum(p_<0.05 == [1 1 0],2)==3; % NMS for Context * Spatial location only
                        else
                            idx =  sum(p_(:,1:2)<0.05,2)==2; %NMS for Context * Spatial location
                        end
                    case 5
                        if exclusive_
                            idx =  sum(p_<0.05 == [0 0 1],2)==3; % NMS for Context * Spatial location only
                        else
                            idx =  p_(:,3)<0.05==1; %NMS for Context * Spatial location
                        end
            end
         x = kappa(idx); x(isnan(x))=[];
         scatter(0*randn(length(x),1)+factor_,x,'ob','filled','Markerfacealpha',0.5)
           
         errorbar(factor_,nanmean(kappa(idx)),nansem(kappa(idx)),'Color','k','LineWidth',1.5)
         bar(factor_,nanmean(kappa(idx)),'FaceColor','none','EdgeColor','k','LineWidth',1.5); 
           
                     
         
%            x = Rbar2_(idx); x(isnan(x))=[];
%            scatter(0*randn(length(x),1)+factor_,x,'ob','filled','Markerfacealpha',0.5)
%            
%          errorbar(factor_,nanmean(Rbar2_(idx)),nansem(Rbar2_(idx)),'Color','k','LineWidth',1.5)
%          bar(factor_,nanmean(Rbar2_(idx)),'FaceColor','none','EdgeColor','k','LineWidth',1.5); 
         title(f.bandNames{iBand})
         set(gca,'Xtick',1:length(factors_),'XtickLabel',factors_,'XTickLabelRotation',-45)

axis([0 6 0 1])
    end
end
end
%% Scatter phase locking on Cartesian axes
factor_ = 3;
exclusive_ = true;

for s=1:2
        figure;
    for iBand = 1:numel(f.bandNames)
        ax = subplot(1,numel(f.bandNames),iBand); hold on
        
        % Plot Sig encoding
        for iDay = 1:length(Filelist)
            
            Phase = phasestats.Phase{s}{iBand}{iDay}';
            Theta = phasestats.Theta{s}{iBand}{iDay}';
            kappa = phasestats.kappa{s}{iBand}{iDay}';
            Rbar_ = phasestats.Rbar_{s}{iBand}{iDay}';
            Rbar2_ = phasestats.Rbar2_{s}{iBand}{iDay}';
            Rp_   = phasestats.Rp_{s}{iBand}{iDay}';
            p_    = selectivity.p{s}{iDay};
            F_    = selectivity.F{s}{iDay};
            
            if ~sum(isnan(F_)>0)
                switch factor_
                    case 1
                        if exclusive_
                            idx =  sum(p_<0.05 == [1 0 0],2)==3; % Simple context only
                        else
                            idx =  p_(:,1)<0.05; % Simple, context
                        end
                    case 2
                        if exclusive_
                            idx =  sum(p_<0.05 == [0 1 0],2)==3; % Simple spatial location only
                        else
                            idx =  p_(:,2)<0.05; % Simple spatial
                        end
                    case 3
                        if exclusive_
                            idx =  sum(p_<0.05 == [0 0 1],2)==3; % NMS for Context * Spatial location only
                        else
                            idx =  p_(:,3)<0.05==1; %NMS for Context * Spatial location
                        end
                end
%                 idx=1:size(p_,1);
%                 idx = idx & Rp_<0.05;
                pS = scatter(Rbar2_(idx),F_(idx,factor_),'b','filled');
                pS.MarkerEdgeAlpha = 0.5;
                pS.MarkerFaceAlpha = 0.9;
            end
        end
        
        
        % Plot non-sig encoding
        for iDay = 1:length(Filelist)
            Phase  = phasestats.Phase{s}{iBand}{iDay}';
            Theta  = phasestats.Theta{s}{iBand}{iDay}';
            kappa  = phasestats.kappa{s}{iBand}{iDay}';
            Rbar_  = phasestats.Rbar_{s}{iBand}{iDay}';
            Rbar2_ = phasestats.Rbar2_{s}{iBand}{iDay}';
            Rp_    = phasestats.Rp_{s}{iBand}{iDay}';
            p_     = selectivity.p{s}{iDay};
            F_     = selectivity.F{s}{iDay};
            
            if ~sum(isnan(F_)>0)
                switch factor_
                    case 1
                        if exclusive_
                            idx =  sum(p_<0.05 == [1 0 0],2)==0; % Simple no context only
                        else
                            idx =  p_(:,1)>0.05; % Simple no context
                        end
                    case 2
                        if exclusive_
                            idx =  sum(p_<0.05 == [0 1 0],2)==0; % Simple no spatial location only
                        else
                            idx =  p_(:,2)>0.05; % Simple no spatial
                        end
                    case 3
                        if exclusive_
                            idx =  sum(p_<0.05 == [0 0 1],2)==0; % no NMS for Context * Spatial location only
                        else
                            idx =  p_(:,3)>0.05==1; %NMS for Context * Spatial location
                        end
                end
                pS = scatter(Rbar_(idx),F_(idx,factor_),'w');
                pS.MarkerEdgeColor = 'k';
                pS.MarkerEdgeAlpha = 0.5;
                pS.MarkerFaceAlpha = 0.5;
            end
        end
%           ax.RLim=[0 10]
%     set(ax,'Xscale','log','Yscale','log')
    end
end
%% *** Contextual vs. Spatial information
figure;
for s=1:2
    F_  = cell2mat(selectivity.F{s}');
    p_  = cell2mat(selectivity.p{s}');
    idx = sum(p_<0.05,2)==0;
    subplot(1,2,s); hold on
    scatter3(F_(~idx,1),F_(~idx,2),F_(~idx,3),'b','filled','MarkerEdgeColor','w');
    scatter3(F_(idx,1),F_(idx,2),F_(idx,3),'r','filled','MarkerEdgeColor','r');
%     axis([0 25 0 12])
    view(3)
    xlabel('Contextual information')
    ylabel('Spatial information')
    zlabel('Spatial*Contextual information')
end
%%
figure;
for s=1:2
    p_  = cell2mat(selectivity.p{s}');
    subplot(1,2,s); hold on
    scatter3(p_(:,1),p_(:,2),p_(:,3),'b','filled','MarkerEdgeColor','w');
    
%     axis([0 25 0 12])
end

%% Scatter phase locking on polar axes - Assembly membership

factor_ = 1;
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

for s=1:2
    figure('name',[Areas{s},' phase locking']);
    for iBand = 1:numel(f.bandNames)
        ax = subplot(1,numel(f.bandNames),iBand,polaraxes); hold on
        
        
        Theta=cell(1,3);kappa=cell(1,3);Rbar2_=cell(1,3);idx=cell(1,3);
        
        for iDay = 1:length(Filelist)
            fIDX = find(strcmp(Filelist{iDay},D.fname)); % match the two file structures
            idx{1} = D.nonmemberUnits{fIDX}{s};
            idx{2} = D.localmemberUnits{fIDX}{s};
            idx{3} = D.jointmemberUnits{fIDX}{s};
            for iType =1:3
                Theta{iType} = [Theta{iType}; phasestats.Theta{s}{iBand}{iDay}(idx{iType})];
                kappa{iType} = [kappa{iType}; phasestats.kappa{s}{iBand}{iDay}(idx{iType})];
                Rbar2_{iType} = [Rbar2_{iType}; phasestats.Rbar2_{s}{iBand}{iDay}(idx{iType})];
            end
        end
        
        for iType =1:3
%             pS = polarscatter(ax,Theta{iType},kappa{iType},30,col_{iType},'filled');
            pS = polarscatter(ax,Theta{iType},Rbar2_{iType},30,col_{iType},'filled');
            pS.MarkerEdgeColor = col_{iType};
            pS.MarkerEdgeAlpha = 0.5;
            pS.MarkerFaceAlpha = 0.5;
        end
        
        %         ax.RLim = [0 0.2];
        title(f.bandNames{iBand})
        colormap(jet)
        caxis([ 0 10])
    end
    
end

%% histogram phase locking on polar axes - Assembly membership

factor_ = 1;
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
bins = deg2rad(0:22.5:360);

for s=1:2
    figure('name',[Areas{s},' phase locking']);
    for iBand = 1:numel(f.bandNames)
        ax = subplot(1,numel(f.bandNames),iBand,polaraxes); hold on
        
        
        Theta=cell(1,3);kappa=cell(1,3);Rbar2_=cell(1,3);idx=cell(1,3);
        
        for iDay = 1:length(Filelist)
            fIDX = find(strcmp(Filelist{iDay},D.fname)); % match the two file structures
            idx{1} = D.nonmemberUnits{fIDX}{s};
            idx{2} = D.localmemberUnits{fIDX}{s};
            idx{3} = D.jointmemberUnits{fIDX}{s};
            for iType =1:3
                Theta{iType} = [Theta{iType}; phasestats.Theta{s}{iBand}{iDay}(idx{iType})];
                kappa{iType} = [kappa{iType}; phasestats.kappa{s}{iBand}{iDay}(idx{iType})];
                Rbar2_{iType} = [Rbar2_{iType}; phasestats.Rbar2_{s}{iBand}{iDay}(idx{iType})];
            end
        end
        
        for iType =1:3
                    polarhistogram(ax, Theta{iType},'BinEdges',bins,'FaceColor',col_{iType},'FaceAlpha',.3,'Normalization','probability');
        end
        
        %         ax.RLim = [0 0.2];
        title(f.bandNames{iBand})
        colormap(jet)
        caxis([ 0 10])
    end
    
end
%% histogram phase locking on polar axes - Assembly membership

factor_ = 1;
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
bins = deg2rad(0:22.5:360);

for s=1:2
    figure('name',[Areas{s},' phase locking']);
    for iBand = 1:numel(f.bandNames)
        ax = subplot(1,numel(f.bandNames),iBand,polaraxes); hold on
        
        
        Theta=cell(1,3);kappa=cell(1,3);Rbar2_=cell(1,3);idx=cell(1,3);
        
        for iDay = 1:length(Filelist)
            fIDX = find(strcmp(Filelist{iDay},D.fname)); % match the two file structures
            idx{1} = D.nonmemberUnits{fIDX}{s};
            idx{2} = D.localmemberUnits{fIDX}{s};
            idx{3} = D.jointmemberUnits{fIDX}{s};
            for iType =1:3
                Theta{iType} = [Theta{iType}; phasestats.Theta{s}{iBand}{iDay}(idx{iType})];
                kappa{iType} = [kappa{iType}; phasestats.kappa{s}{iBand}{iDay}(idx{iType})];
                Rbar2_{iType} = [Rbar2_{iType}; phasestats.Rbar2_{s}{iBand}{iDay}(idx{iType})];
            end
        end
        
        for iType =1:3
                    polarhistogram(ax, Theta{iType},'BinEdges',bins,'FaceColor',col_{iType},'FaceAlpha',.3,'Normalization','probability');
        end
        
        %         ax.RLim = [0 0.2];
        title(f.bandNames{iBand})
        colormap(jet)
        caxis([ 0 10])
    end
    
end

%% *** Plot phase pref vs membership as bars
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

for s=1:2
    figure('name',[Areas{s},' phase locking']);
    for iBand = 1:numel(f.bandNames)
        ax = subplot(1,numel(f.bandNames),iBand); hold on
        
        kappa=cell(1,3);Rbar2_=cell(1,3);idx=cell(1,3);
         for iDay = 1:length(Filelist)
            fIDX = find(strcmp(Filelist{iDay},D.fname)); % match the two file structures
            idx{1} = D.nonmemberUnits{fIDX}{s};  
            idx{2} = D.localmemberUnits{fIDX}{s};  
            idx{3} = D.jointmemberUnits{fIDX}{s};  
            for iType =1:3
                kappa{iType} = [kappa{iType}; phasestats.kappa{s}{iBand}{iDay}(idx{iType})];
                Rbar2_{iType} = [Rbar2_{iType}; phasestats.Rbar2_{s}{iBand}{iDay}(idx{iType})];
            end
         end
%          bar(1,nanmean(kappa(idx))); 
%          bar(2,nanmean(kappa(~idx)))
%          errorbar([1,2],[nanmean(kappa(idx)),nanmean(kappa(~idx))],[nansem(kappa(idx)),nansem(kappa(~idx))])
           
                     
          for iType =1:3
            x = Rbar2_{iType}; x(isnan(x))=[];
%             x = kappa{iType}; x(isnan(x))=[];
            scatter(0.2*randn(length(x),1)+iType,x,'filled','MarkerFaceColor',col_{iType},'MarkerEdgeColor',col_{iType},'Markerfacealpha',0.5)
            errorbar(iType,nanmean(x),nansem(x),'Color','k','LineWidth',1.5) 
            bar(iType,nanmean(x),'FaceColor','none','EdgeColor','k','LineWidth',1.5); 
          end
         
       
         title(f.bandNames{iBand})
         set(gca,'Xtick',1:3,'XtickLabel',{'Non-members','Local Members','Joint members'},'XTickLabelRotation',-45)
         if iBand==1
             ylabel('Mean resultant vector length')
         end
        axis([0 4 0 0.3])
    end
end
%% *** 3D Scatter of membership vs info content
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
figure('name','Assembly membership vs. information content');

for s=1:2
        subplot(1,2,s); hold on
        
        F_  = cell2mat(selectivity.F{s}');
        p_  = cell2mat(selectivity.p{s}');
        idx = sum(p_<0.05,2)==0;
    
        idx=cell(1,3);
        p_=cell(1,3); F_=cell(1,3);
        for iDay = 1:length(Filelist)
            fIDX = find(strcmp(Filelist{iDay},D.fname)); % match the two file structures
            idx{1} = D.nonmemberUnits{fIDX}{s};
            idx{2} = D.localmemberUnits{fIDX}{s};
            idx{3} = D.jointmemberUnits{fIDX}{s};
            try
                for iType =1:3
                    F_{iType} = [F_{iType}; selectivity.F{s}{iDay}(idx{iType},:)];
                    p_{iType} = [p_{iType}; selectivity.p{s}{iDay}(idx{iType},:)];
                end
            end
        end
        
    for iType =1:3
        idx_ = sum(p_{iType}<0.05,2)==0;
%         scatter3(F_{iType}(~idx_,1),F_{iType}(~idx_,2),F_{iType}(~idx_,3),'filled','MarkerEdgeColor',col_{iType},'MarkerFaceColor',col_{iType});
%         scatter3(F_{iType}(idx_,1),F_{iType}(idx_,2),F_{iType}(idx_,3),'MarkerEdgeColor',col_{iType});
        
        stem3(F_{iType}(~idx_,1),F_{iType}(~idx_,2),F_{iType}(~idx_,3),'filled','color',col_{iType},'MarkerEdgeColor',col_{iType},'MarkerFaceColor',col_{iType});
        stem3(F_{iType}(idx_,1),F_{iType}(idx_,2),F_{iType}(idx_,3),'color',col_{iType},'MarkerEdgeColor',col_{iType},'MarkerFaceColor','w');

    end
        axis([0 30 0 40 0 Inf])
    view(3); grid on; title([Areas{s} ' units'])
    xlabel('Contextual')
    ylabel('Spatial')
    zlabel('Spatial*Contextual')
end
%% *** Bars of membership vs info content
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
Tunings =     {'Contextual information','Spatial information','Spatial*Contextual information'};

resultsTemp = cell(2,1);

for s=1:2
            figure('name',[Areas{s},' Tuning']);

            
            
        
        F_  = cell2mat(selectivity.F{s}');
        p_  = cell2mat(selectivity.p{s}');
        idx = sum(p_<0.05,2)==0;
        

        idx=cell(1,3);
        p_=cell(1,3); F_=cell(1,3);
        for iDay = 1:length(Filelist)
            fIDX = find(strcmp(Filelist{iDay},D.fname)); % match the two file structures
            idx{1} = D.nonmemberUnits{fIDX}{s};
            idx{2} = D.localmemberUnits{fIDX}{s};
            idx{3} = D.jointmemberUnits{fIDX}{s};
            try
                for iType =1:3
                    F_{iType} = [F_{iType}; selectivity.F{s}{iDay}(idx{iType},:)];
                    p_{iType} = [p_{iType}; selectivity.p{s}{iDay}(idx{iType},:)];
                end
            end
        end
        
        for iTuning =1:3
            subplot(1,3,iTuning); hold on
            for iType =1:3
                idx_ = sum(p_{iType}(:,iTuning)<0.05,2)==0;
                
                resultsTemp{s}{iType}= F_{iType};

                x =F_{iType}(~idx_,iTuning);
                scatter(0.2*randn(length(x),1)+iType,x,'filled','MarkerFaceColor',col_{iType},'MarkerEdgeColor',col_{iType},'Markerfacealpha',0.5)
                x =F_{iType}(idx_,iTuning);
                scatter(0.2*randn(length(x),1)+iType,x,'filled','MarkerFaceColor','w','MarkerEdgeColor',col_{iType},'Markerfacealpha',0.5)
                x = F_{iType}(~idx_,iTuning);
%                 x = F_{iType}(:,iTuning);
                errorbar(iType,nanmean(x),nansem(x),'Color','k','LineWidth',1.5)
                bar(iType,nanmean(x),'FaceColor','none','EdgeColor','k','LineWidth',1.5);
                if iTuning == 1
                    
                    ylabel('F-score')
                end
            end
                axis([0 4 0 25])
                set(gca,'Xtick',1:3,'XtickLabel',{'Non-members','Local Members','Joint members'},'XTickLabelRotation',-45)
                title(Tunings{iTuning})

        end

end

ylim([0 100])

%%
%                 resultsTemp{s}{iType,iTuning}= F_; %{iArea}{Ass type}[no. units x tuning variable]
s       = 1;  %  {PFC CA1}
Tuning  = 3;  %  {'Contextual information','Spatial information','Spatial*Contextual information'};
               [p,tbl,stats] = anova1([resultsTemp{s}{1}(:,Tuning);...
                                       resultsTemp{s}{2}(:,Tuning);...
                                       resultsTemp{s}{3}(:,Tuning)],...
                                    [  ones(length(resultsTemp{s}{1}(:,Tuning)),1);...
                                     2*ones(length(resultsTemp{s}{2}(:,Tuning)),1);...
                                     3*ones(length(resultsTemp{s}{3}(:,Tuning)),1)]);
                
[c,~,~,gnames] = multcompare(stats);



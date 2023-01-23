clear 
cd('C:\Analysis\AssemblyAnalysis\LFPanalysis\Aleks chosen LFP segments\CellsOnTetrodes')
load('cellsonTetrodes.mat')
pat = 'C:\Analysis\AssemblyAnalysis\raw';
cd(pat);
useWholeTaskPeriod = true;
Targets= {'SHORT','MEDIUM','LONG'};
iTarget = 3;

UseResampledResults = true;

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
if UseResampledResults
    eval(sprintf('MixedSelectivityfileList_%s = dir([pat ''\\MixedSelectivity\\TrialNosMatched\\*%s*PFC*.mat'']);',Targets{iTarget},Targets{iTarget}))
else
    eval(sprintf('MixedSelectivityfileList_%s = dir([pat ''MixedSelectivity\\*%s*PFC*.mat'']);',Targets{iTarget},Targets{iTarget}) )
end
%% Batch process: Compute phase-locking and firing phase preferences
bins = -180:45:180; % for histograms
TimeSpan = 1;
threshChoice = 4; % 2=p<0.05 , 3=p<0.01

plotAsswideLockingOnlineYN = true;
kappa_comparison = cell(numel(f.bandNames),1);
phase_comparison = cell(numel(f.bandNames),1);
Rayleighp_comparison = cell(numel(f.bandNames),1);
Rbar_comparison = cell(numel(f.bandNames),1);

kappa_comparisonOutsiders = cell(numel(f.bandNames),1);
phase_comparisonOutsiders = cell(numel(f.bandNames),1);
Rayleighp_comparisonOutsiders = cell(numel(f.bandNames),1);
Rbar_comparisonOutsiders = cell(numel(f.bandNames),1);
for iDay = setdiff(1:length(tetrodes.names),[7,8])
    %% Load data
    
    % Load cells vs tetrode assignments
    fn = fullfile(pat,sprintf('%s.mat',tetrodes.names{iDay}));
    spikes = load(fn);
    
    % Load LFP phase 
    fn = fullfile(pat,'LFPphase',sprintf('%s_phase.mat',tetrodes.names{iDay}));
    load(fn)
     %f.phaseData{3} = [f.phaseData{1}(B.usel_out{1});...
    %                  f.phaseData{2}(B.usel_out{2})];
    f.phaseData{3} = [f.phaseData{1};f.phaseData{2}];    
    
    % Load Assemblies
    if useWholeTaskPeriod
        fn = fullfile(pat,'KDE_binsTaskOnly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_FSC.mat',tetrodes.names{iDay}));
        A  = load(fn);
        fn = fullfile(pat,'KDE_binsTaskOnly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_AssemRes2.mat',tetrodes.names{iDay}));
        B  = load(fn);
        % recalculate included units
        fn = fullfile(pat,'KDE_binsTaskOnly',sprintf('%s_PFC_iFR50_behavOnly.mat',tetrodes.names{iDay}));
        usel_out=SelCells(fn,0.1,1e6);
    else
        fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',tetrodes.names{iDay}));        
        A  = load(fn);
        fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',tetrodes.names{iDay}));        
        B  = load(fn);
        % recalculate included units
        fn = fullfile(pat,'KDE_bins',sprintf('%s_PFC_iFR50.mat',tetrodes.names{iDay}));
        [~,~,~,~,B.usel_out]=SelTrialsCellsWholeTrial(fn,10,0.1,1e8,'iFR');
    end
    spikes.jointcells=[spikes.PFCcells;spikes.HPcells]';
    A.nu(3) =sum(A.nu(1:2));
    B.usel_out{3}= [B.usel_out{1},B.usel_out{2}+max(B.usel_out{1})];
    
    % Load mixed selectivity assignments and F-scores
    for iArea = 1:2
        % fidx = find(~cellfun(@isempty,strfind({MixedSelectivityfileList_LONG.name},tetrodes.names{iDay})))
        if UseResampledResults
            fn = [pat '\MixedSelectivity\TrialNosMatched\' tetrodes.names{iDay} '_' Areas{iArea} '_MixedSelectivity_Units.mat'];
        else
            fn = [pat '\MixedSelectivity\' tetrodes.names{iDay} '_' Areas{iArea} '_MixedSelectivity_Units.mat'];
        end
        MixedSelectivity{iArea} = load(fn);
        selectivity.F{iArea} = MixedSelectivity{iArea}.D.Long.ANOVA2.F_cor;
        selectivity.p{iArea} = MixedSelectivity{iArea}.D.Long.ANOVA2.p_cor;
    end
    
    %% 
    figure; hold on
    for s = 2
        for iBand = 1:numel(f.bandNames)
            for iUnit =1:length(f.phaseData{s})
                
                
                
                % get spike phase for this band, compute phase locking
                Theta_   = f.phaseData{s}{iUnit}(iBand,:);
                [~, kappa_(iUnit), ~]   = circStats(Theta_);
                
                
                % Get peak phase preference
                Phase_(iUnit)   = mode(round(rad2deg(Theta_)));
                
                % Rayleigh p value and mean resultant length vector
                [Rp_(iUnit),   Rbar_(iUnit)] = rayleigh(Theta_');
                
                % Rayleigh test only good for n<10, enforce here
                if numel(Theta_)<10, Rp_(iUnit) = 1;end
                
                
            end
            
        
%              q_ = [sum(sum(p_thresh,2)==0),...              % No tuning
%                         sum(sum(p_thresh == [1 0 0],2)==3),...     % CS for Context
%                         sum(sum(p_thresh == [0 1 0],2)==3),...     % CS for Location
%                         sum(sum(p_thresh == [1 0 0],2)==3) + sum(sum(p_thresh == [0 1 0],2)==3),... % Any CS
%                         sum(sum(p_thresh == [1 1 0],2)==3),...     % LMS for Context x Location
%                         sum(sum(p_thresh(:,3) == 1,2)),...         % NMS for Context x Location
%                         sum(sum(p_thresh,2)>0)];                   % Any tuning
                    
subplot(1,numel(f.bandNames),iBand); hold on
p_ = selectivity.p{s};
x_ = Phase_;
y_ = Rbar_;
%%% Context
% idx =  sum(p_(:,[1 3])<0.05,2)==0; % Simple, no context
%   scatter(x_(idx),y_(idx), 20,selectivity.F{s}(idx,1))
% idx =  sum(p_<0.05 == [1 0 0],2)==3; %Simple, context
%   scatter(x_(idx),y_(idx), 20,selectivity.F{s}(idx,1),'filled')
% idx =  p_(:,3)<0.05==1; %NMS for Context * Spatial location
% %   scatter(x_(idx),y_(idx), 20,selectivity.F{s}(idx,3),'*','filled')


%%% Location
idx =  sum(p_(:,2:3)<0.05,2)==0; % Simple, no Spatial location
  scatter(x_(idx),y_(idx), 20,selectivity.F{s}(idx,2))
idx =  sum(p_<0.05 == [0 1 0],2)==3; %Simple, Spatial location
  scatter(x_(idx),y_(idx), 20,selectivity.F{s}(idx,2),'filled')
idx =  p_(:,3)<0.05==1; %NMS for Context * Spatial location
%   scatter(x_(idx),y_(idx),50,selectivity.F{s}(idx,2),'*')


%%% Mixed
% idx = sum(p_(:,1:2)<0.05,2)==0; % no MS
% idx = p_(:,3)<0.05==0; % no MS
%     scatter(x_(idx),selectivity.F{s}(idx,3),'b')
% idx =  sum(p_<0.05 == [1 1 0],2)==3; %Linear-mixed
% scatter(x_(idx),mean(selectivity.F{s}(idx,1:2),2),'b','filled')
% idx =  p_(:,3)<0.05 == 1; %NMS for Context x Location
% scatter(x_(idx),selectivity.F{s}(idx,3),'*b')

colormap jet
title(f.bandNames{iBand})

        axis([-180 180 -Inf Inf])
switch iBand
    case 1
        ylabel('Phase locking')
    case 3
        xlabel('Mean phase preference (deg)')
end
        end
    end
    
     %% bandpassed phase-locking wrt assembly activation
%         
%    
%     
%     
%     for s =3%:3
%         
%         nAss = length(A.units{s});
%         ActiveThresh = A.ciHsc(s,threshChoice);
%         
%         nonmemberUnits = setdiff(B.usel_out{s},unique(B.usel_out{s}(cell2mat(A.units{s}))));
% 
% 
%         for iAss=1:nAss
%             disp(sprintf('Analysing day %d,area %d,Assembly %d\n',iDay,s,iAss));
%             memberUnits = B.usel_out{s}(A.units{s}{iAss});
%            
%             histActive = cell(length(memberUnits),1);
%             histInactive = cell(length(memberUnits),1);
%             % Stats for member units
%             for iUnit =1:length(memberUnits)
%                 thisUnit = memberUnits(iUnit);
%                 eval(sprintf('times = spikes.%scells{thisUnit}.t/1e4;',Areas{s}));
%                 
%                 isActive = A.FSCsel{s}(:,iAss) >= ActiveThresh;
%                 isActiveTimes = A.Tmtx( isActive );
%                 isInactive = A.FSCsel{s}(:,iAss) < ActiveThresh;
%                 isInactiveTimes = A.Tmtx( isInactive );
%                 
%                 %times(times<min(A.Tmtx) | times>max(A.Tmtx))=[];
%                 times_idx = closest( A.Tmtx , times ); % indices of timebase closest to spike times
%                 times_    = A.Tmtx( times_idx );        % values  of timebase closest to spike times
%                 
%                 % Only keep spikes within the time axis of the trial cutouts
%                 retain = abs( times_ - times ) < TimeSpan;
%                 
%                 
%                 ActiveSpikeIDX   = find( ismember( times_ , A.Tmtx( isActive ) ) & retain);
%                 InactiveSpikeIDX = find( ismember( times_ , A.Tmtx( isInactive ) ) & ~retain);
%                 
%                 % Optional plotting - raw spike times in and of assembly activity
%                 if plotRawOnlineYN
%                     % Plot the raw spike times
%                     figure; hold on
%                     plot( A.Tmtx , A.FSCsel{s}(:,iAss))
%                     % scatter(times,-4*ones(length(times),1))
%                     % scatter(times_,-4*ones(length(times_),1))
%                     plot( [ min(A.Tmtx) max(A.Tmtx) ] ,[ ActiveThresh ActiveThresh ],':r')
% 
%                     T = times( ActiveSpikeIDX );
%                     scatter( T, ones(length(T) , 1 ),'b');
%                     T = times( InactiveSpikeIDX );
%                     scatter( T, ones(length(T) , 1 ),'r');
%                     
%                 end
%                 % Optional plotting - set up a plot for the histograms
%                 if plotOnlineYN | plotAsswideLockingOnlineYN 
%                     figure('position',[101 236 1439 526]);
%                 end
%                 
%                 % Crunch the numbers
%                 %histActive{iBand} = zeros(length(bins),length(memberUnits));
%                 %histInactive{iBand} = zeros(length(bins),length(memberUnits));
%            
%                 for iBand = 1:numel(f.bandNames)
%                    
%                     % get spike phase for this band, compute phase locking
%                     Theta_Active   = f.phaseData{s}{thisUnit}(iBand,ActiveSpikeIDX);
%                     Theta_Inactive = f.phaseData{s}{thisUnit}(iBand,InactiveSpikeIDX);
%                     [~, kappa_Active, ~]   = circStats(Theta_Active);
%                     [~, kappa_Inactive, ~] = circStats(Theta_Inactive);
%                     f.kappa_Active{iDay}{s}{iAss}{thisUnit}(iBand)   = kappa_Active;
%                     f.kappa_Inactive{iDay}{s}{iAss}{thisUnit}(iBand) = kappa_Inactive;
%                     f.kappa_Shift{iDay}{s}{iAss}{thisUnit}(iBand)    = kappa_Active-kappa_Inactive;
%                     
%                     % Get peak phase preference
%                     Phase_Active   = mode(round(rad2deg(Theta_Active)));
%                     Phase_Inactive = mode(round(rad2deg(Theta_Inactive)));
%                     f.peakPhaseActive{iDay}{s}{iAss}{thisUnit}(iBand)   = Phase_Active;
%                     f.peakPhaseInactive{iDay}{s}{iAss}{thisUnit}(iBand) = Phase_Inactive;
%                     f.peakPhaseShift{iDay}{s}{iAss}{thisUnit}(iBand)    = Phase_Active-Phase_Inactive;
%                     
%                     % Rayleigh p value and mean resultant length vector
%                     [p_Active,   Rbar_Active] = rayleigh(Theta_Active');
%                     [p_Inactive, Rbar_Inactive] = rayleigh(Theta_Inactive');
%                     
%                     % Rayleigh test only good for n<10, enforce here
%                     if numel(Theta_Active)<10, p_Active = 1;end
%                     if numel(Theta_Inactive)<10, p_Inactive = 1;end
%                     
%                     f.p_Active{iDay}{s}{iAss}{thisUnit}(iBand) = p_Active;
%                     f.p_Inactive{iDay}{s}{iAss}{thisUnit}(iBand) = p_Inactive;
%                     f.Rbar_Active{iDay}{s}{iAss}{thisUnit}(iBand) = Rbar_Active;
%                     f.Rbar_Inactive{iDay}{s}{iAss}{thisUnit}(iBand) = Rbar_Inactive;
%                     
%                     % cheap and dirty log of all kappas and phase shifts
%                     kappa_comparison{iBand} = [kappa_comparison{iBand};[kappa_Inactive, kappa_Active]];
%                     phase_comparison{iBand} = [phase_comparison{iBand};[Phase_Active-Phase_Inactive]];
%                     Rayleighp_comparison{iBand} = [Rayleighp_comparison{iBand};[p_Inactive, p_Active]];
%                     Rbar_comparison{iBand} = [Rbar_comparison{iBand};[Rbar_Inactive, Rbar_Active]];
% 
%                     h_Active = histc(rad2deg(Theta_Active),bins)';
%                     % h_Active(h_Active==0)=NaN;
%                     % h_Active=h_Active./sum(h_Active);
%                     h_Active=smooth_hist(h_Active);
%                     
%                     h_Active=mat2gray(h_Active);
%                     h_Active(end) = h_Active(1);
%                     
%                     h_Inactive = histc(rad2deg(Theta_Inactive),bins)';
%                     % h_Inactive(h_Inactive==0)=NaN;
%                     % h_Inactive = h_Inactive./sum(h_Inactive);
%                     h_Inactive=smooth_hist(h_Inactive);
%                     h_Inactive=mat2gray(h_Inactive);
%                     h_Inactive(end) = h_Inactive(1);
%                     
%                     histActive{iBand}(:,iUnit)  = h_Active;
%                     histInactive{iBand}(:,iUnit) = h_Inactive;
% 
%                     if plotOnlineYN
%  
%                         subplot(2,numel(f.bandNames),iBand); hold on
%                         
%                         b= bar(bins,(h_Inactive),'FaceColor','r',...
%                             'EdgeColor','r',...
%                             'BaseValue',0,'FaceAlpha',0.8,'BarWidth',1,'Basevalue',0);
%                         b= bar(bins,(h_Active),'FaceColor','b',...
%                             'EdgeColor','b',...
%                             'EdgeAlpha',0,'BaseValue',0,'FaceAlpha',0.8,'BarWidth',1,'Basevalue',0);
% %                         b= bar(bins,(h_Inactive),'FaceColor','r',...
% %                             'EdgeColor','r',...
% %                             'BaseValue',0,'BarWidth',1,'Basevalue',0);
% %                         b= bar(bins,(h_Active),'FaceColor','b',...
% %                             'EdgeColor','b',...
% %                             'BaseValue',0,'BarWidth',1,'Basevalue',0);
% %                         plot(bins,sind(bins)+2,'k','LineWidth',2)
%                         title(f.bandNames{iBand})
%                         
%                         subplot(2,length(f.bandNames),iBand+length(f.bandNames),polaraxes); hold on
%                         pInActive = polarplot(deg2rad(bins),h_Inactive,'r');
%                         pActive = polarplot(deg2rad(bins),h_Active,'b');
%                         pInActive.LineWidth=1.5;
%                         pActive.LineWidth=1.5;
%                         title(sprintf('\\kappa: Active=%1.2f, Inactive=%1.2f',kappa_Active,kappa_Inactive))
%                     end
%                    
%                     
%                 end
%                 
%              
%             end
%             
%             if plotAsswideLockingOnlineYN
%                  col_ =rand(length(memberUnits),3);
%                  for iBand = 1:numel(f.bandNames)
%                         subplot(2,numel(f.bandNames),iBand);hold on
%                         for iUnit =1:length(memberUnits)
%                             b= bar(bins,(histActive{iBand}(:,iUnit)),'FaceColor',col_(iUnit,:),...
%                             'EdgeColor',col_(iUnit,:),...
%                             'FaceAlpha',0.8,'BarWidth',1,'Basevalue',0);
% %                             b = stairs(bins,iUnit*1+mat2gray(histActive{iBand}(:,iUnit)))
%                         end
%                          subplot(2,numel(f.bandNames),iBand+numel(f.bandNames));hold on
%                         for iUnit =1:length(memberUnits)
%                             b= bar(bins,(histInactive{iBand}(:,iUnit)),'FaceColor',col_(iUnit,:),...
%                             'EdgeColor',col_(iUnit,:),...
%                             'FaceAlpha',0.8,'BarWidth',1,'Basevalue',0);
% %                             b = stairs(bins,iUnit*1+mat2gray(histActive{iBand}(:,iUnit)))
%                         end
%                        
%                  end
%              end
%                     
%              
%             % Stats for non-member units
%             for iUnit =1:length(nonmemberUnits)
%                 thisUnit = nonmemberUnits(iUnit);
%                 eval(sprintf('times = spikes.%scells{thisUnit}.t/1e4;',Areas{s}));
%                 
%                 isActive = A.FSCsel{s}(:,iAss)>=ActiveThresh;
%                 isActiveTimes = A.Tmtx(isActive);
%                 isInactive = A.FSCsel{s}(:,iAss)<ActiveThresh;
%                 isInactiveTimes = A.Tmtx(isInactive);
%                 
%                 %times(times<min(A.Tmtx) | times>max(A.Tmtx))=[];
%                 times_idx = closest(A.Tmtx,times); % indices of timebase closest to spike times
%                 times_ = A.Tmtx(times_idx);        % values  of timebase closest to spike times
%                 
%                 % Only keep spikes within the time axis of the trial cutouts
%                 retain = abs(times_-times)<TimeSpan;
%                
%                 ActiveSpikeIDX   = find(ismember(times_,A.Tmtx(isActive)) & retain);
%                 InactiveSpikeIDX = find(ismember(times_,A.Tmtx(isInactive)) & ~retain);
%                 
%               
%                 
%                 % Crunch the numbers
%                 for iBand = 1:numel(f.bandNames)
%                     % get spike phase for this band, compute phase locking
%                     Theta_Active   = f.phaseData{s}{thisUnit}(iBand,ActiveSpikeIDX);
%                     Theta_Inactive = f.phaseData{s}{thisUnit}(iBand,InactiveSpikeIDX);
%                     [~, kappa_Active, ~]   = circStats(Theta_Active);
%                     [~, kappa_Inactive, ~] = circStats(Theta_Inactive);
%                     f.kappa_ActiveOutsiders{iDay}{s}{iAss}{thisUnit}(iBand)   = kappa_Active;
%                     f.kappa_InactiveOutsiders{iDay}{s}{iAss}{thisUnit}(iBand) = kappa_Inactive;
%                     f.kappa_ShiftOutsiders{iDay}{s}{iAss}{thisUnit}(iBand)    = kappa_Active-kappa_Inactive;
%                     
%                     % Get peak phase preference
%                     Phase_Active   = mode(round(rad2deg(Theta_Active)));
%                     Phase_Inactive = mode(round(rad2deg(Theta_Inactive)));
%                     f.peakPhaseActiveOutsiders{iDay}{s}{iAss}{thisUnit}(iBand)   = Phase_Active;
%                     f.peakPhaseInactiveOutsiders{iDay}{s}{iAss}{thisUnit}(iBand) = Phase_Inactive;
%                     f.peakPhaseShiftOutsiders{iDay}{s}{iAss}{thisUnit}(iBand)    = Phase_Active-Phase_Inactive;
%                     
%                     % Rayleigh p value and mean resultant length vector
%                     [p_Active,   Rbar_Active] = rayleigh(Theta_Active');
%                     [p_Inactive, Rbar_Inactive] = rayleigh(Theta_Inactive');
%                     
%                     % Rayleigh test only good for n<10, enforce here
%                     if numel(Theta_Active)<10, p_Active = 1;end
%                     if numel(Theta_Inactive)<10, p_Inactive = 1;end
%                     
%                     f.p_ActiveOutsiders{iDay}{s}{iAss}{thisUnit}(iBand) = p_Active;
%                     f.p_InactiveOutsiders{iDay}{s}{iAss}{thisUnit}(iBand) = p_Inactive;
%                     f.Rbar_ActiveOutsiders{iDay}{s}{iAss}{thisUnit}(iBand) = Rbar_Active;
%                     f.Rbar_InactiveOutsiders{iDay}{s}{iAss}{thisUnit}(iBand) = Rbar_Inactive;
%                     
%                     % cheap and dirty log of all kappas and phase shifts
%                     kappa_comparisonOutsiders{iBand} = [kappa_comparisonOutsiders{iBand};[kappa_Inactive, kappa_Active]];
%                     phase_comparisonOutsiders{iBand} = [phase_comparisonOutsiders{iBand};[Phase_Active-Phase_Inactive]];
%                     Rayleighp_comparisonOutsiders{iBand} = [Rayleighp_comparisonOutsiders{iBand};[p_Inactive, p_Active]];
%                     Rbar_comparisonOutsiders{iBand}      = [Rbar_comparisonOutsiders{iBand};[Rbar_Inactive, Rbar_Active]];
% 
%                     if plotOnlineYN
%                     
% %                     % Optional plotting
% %                     figure('position',[101 236 1439 526]);
% %                         h_Active = histc(rad2deg(Theta_Active),bins)';
% %                         % h_Active(h_Active==0)=NaN;
% %                         % h_Active=h_Active./sum(h_Active);
% %                         h_Active=smooth_hist(h_Active);
% %                         
% %                         h_Active=mat2gray(h_Active);
% %                         h_Active(end) = h_Active(1);
% %                         
% %                         h_Inactive = histc(rad2deg(Theta_Inactive),bins)';
% %                         % h_Inactive(h_Inactive==0)=NaN;
% %                         % h_Inactive = h_Inactive./sum(h_Inactive);
% %                         h_Inactive=smooth_hist(h_Inactive);
% %                         h_Inactive=mat2gray(h_Inactive);
% %                         h_Inactive(end) = h_Inactive(1);
% %                         
% %                         subplot(2,numel(f.bandNames),iBand);hold on
% %                         
% %                         b= bar(bins,(h_Inactive),'FaceColor','r',...
% %                             'EdgeColor','r',...
% %                             'EdgeAlpha',0,'BaseValue',0,'FaceAlpha',0.8,'BarWidth',1,'Basevalue',0);
% %                         b= bar(bins,(h_Active),'FaceColor','b',...
% %                             'EdgeColor','b',...
% %                             'EdgeAlpha',0,'BaseValue',0,'FaceAlpha',0.8,'BarWidth',1,'Basevalue',0);
% %                         plot(bins,sind(bins)+2,'k','LineWidth',2)
% %                         title(f.bandNames{iBand})
% %                         
% %                         subplot(2,length(f.bandNames),iBand+length(f.bandNames),polaraxes); hold on
% %                         pInActive = polarplot(deg2rad(bins),h_Inactive,'r');
% %                         pActive = polarplot(deg2rad(bins),h_Active,'b');
% %                         pInActive.LineWidth=1.5;
% %                         pActive.LineWidth=1.5;
% %                         title(sprintf('\\kappa: Active=%1.2f, Inactive=%1.2f',kappa_Active,kappa_Inactive))
%                     end
%                     
%                 end
%                 
%             end
%            
%         end
%         
%     end
    
end

clear ActiveSpikeIDX ActiveThresh B h_Active h_Inactive iAss iBand iDay InactiveSpikeIDX 
clear isActive isActiveTimes isInactive isInactiveTimes iUnit kappa_Active kappa_Inactive 
clear memberUnits nAss Phase_Active Phase_Inactive plotOnlineYN Rbar retain s 
clear Theta_Active Theta_Inactive thisUnit times times_ times_idx p_Active p_Inactive plotRawOnlineYN Rbar_Active R
% Example plots to show various combinations of unit timestamps, firing rates, behavioural events and factor scores.
% NB Feb 2016 - This is now depreciated and is obsolete, but retained for example plot purposes.
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk


%% Preamble: load files
% start in 'raw'
% cd 'C:\Analysis\AssemblyAnalysis\raw'
% (1) load raw
load('JaroslawLONG2.mat')
% (2) load iFR50
% load('50ms_bins/JaroslawLONG2_PFC_iFR50.mat')
% (3) get aligned trial markers
[TmtxS,iFR0,EvtLs,EvtTs,usel]=SelTrialsCells('KDE_bins/JaroslawLONG2_PFC_iFR50.mat',10,0.1,1e6,'iFR');

load('KDE_bins/LONG/JaroslawLONG2_iFR50_FSC.mat')
%% Example plot: unit timestamps against FR with behavioural events
% Factor participating units
factor_active_units = usel{1}(units{1}{1});
factor_active_iFRs  = units{1}{1};

unit_id=factor_active_units(1);
iFR_id=factor_active_iFRs(1);
temp_firing=PFCcells{(unit_id)}.t/1e4;
figure; hold on
% 
samples_=[trangeleft_sample(:,1)/1e6;trangeright_sample(:,1)/1e6];
scatter(samples_,ones(numel(samples_),1),'or')

choices_=[trangeleft_choice(:,1)/1e6;trangeleft_choice(:,1)/1e6];
scatter(choices_,ones(numel(choices_),1),'ob')

scatter(temp_firing,zeros(numel(temp_firing),1))
for plot_id=1:numel(iFR0{1})
    plot(TmtxS{1}{plot_id},...
         iFR0{1}{plot_id}(:,iFR_id),'b')
end
title('Example plot: firing rate vs. timestamps with  behavioral events overlaid')

legend('sample times','choice times','spikes','firing rate')
axis([0 Inf -1 Inf])
%% Example plot: raster of units that went into analysis pipeline (i.e. chosen by minimum firing rate etc. )
region = 1; % {PFC, HP, joint}
threshold = 1;
plot_win = [4000 8000];
% plot_win = [min(TmtxS{1}{1}), max(TmtxS{1}{1})];
colorscale=flipud(lines(size(FSCsel{region},2)));
lag=0.05;    
figure; 
subplot(2,1,1)
bar(EvtTs,EvtLs)
xlim(plot_win)
ylabel('event type')
subplot(2,1,2); hold on
% plot ALL units going in
for unit_id=1:numel(usel{region})
    this_unit=usel{region}(unit_id);
    temp_firing=PFCcells{this_unit}.t/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
    scatter(temp_firing,0.1*unit_id*ones(length(temp_firing),1),10,'k','filled')
end

% plot ONLY ASSEMBLY units    
for assembly_id=1:size(FSCsel{region},2);
    % factor scores
    this_FSC=FSCsel{region}(:,assembly_id); % plot(-5:0.1:10,histc(this_FSC,-5:0.1:10)) % histogram to choose activation threshold
    % find factor activation times
    [pks,locs]=findpeaks(this_FSC,'MINPEAKHEIGHT',1);
    assembly_times=Tmtx(locs);
    
    %%%% 
    % plot the factor activations
%         scatter(assembly_times,pks,5,colorscale(assembly_id,:),'filled')
        plot(Tmtx,0.1*this_FSC,'color',colorscale(assembly_id,:))
    %%%%
    
    

    % find factor participating units
    factor_active_units = usel{region}(units{region}{assembly_id}); % usel{brain region} and units{brain region}{assembly number}
    factor_active_iFRs  = units{region}{assembly_id};

    for unit_id=1:numel(factor_active_units)
        this_unit=factor_active_units(unit_id);
        temp_firing=PFCcells{this_unit}.t/1e4;
%         temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
        % highlight spike times of assembly units
        if ismember(this_unit,factor_active_units)
            for assem_peak_no=1:numel(assembly_times)
                temp_firing_ids=(temp_firing-lag)<=assembly_times(assem_peak_no) &...
                                (temp_firing+lag)>=assembly_times(assem_peak_no);
                scatter(temp_firing(temp_firing_ids),0.1*find(this_unit==usel{region})*ones(length(temp_firing(temp_firing_ids)),1),40,colorscale(assembly_id,:),'filled')

            end
            
            
        end
    end
    
    
end
% axis([plot_win(1) plot_win(2) 0 0.12*numel(usel{1})])
axis tight;
xlim(plot_win)
xlabel('time')
title('raster highlighting assembly activations and timestamps')
%
%% Example plot: spike times, firing rates and factor scores
unit_id=factor_active_units(1);
figure; hold on
% (1) plot spike times
this_FSC=PFCcells{unit_id}.t/1e4;
y=ones(length(this_FSC),1);
scatter(this_FSC,y,'b')

% (2) plot Firing rates
for trial_id=1:size(iFR0{1},2)
    plot(TmtxS{1}{trial_id},...
         iFR0{1}{trial_id}(:,iFR_id))
end
% plot trial markers
% plot(EvtTs,EvtLs)

% plot factor score
plot(Tmtx,FSCsel{1}(:,1),'r')
%% e.g. plot all cells' firing rates in all trial epochs

%%
% regions are in the order: {PFC,Hippocampus}
figure; hold on
u_id=usel{1}(2);
    for s_id=1:numel(TmtxS{1})
        plot(TmtxS{1}{s_id},iFR0{1}{s_id}(:,u_id))
        
    end
    this_FSC=PFCcells{u_id}.t*1e-4;
    scatter(this_FSC,ones(size(this_FSC)),'r')
    axis([0 Inf 0 50])
    scatter(EvtTs,EvtLs,'filled','b')
    
    
%% Choose number of factors using bootstrap CIs
% NB load the factor score datasets...

s=1; % region
u=1; % factor no.
[ntp,nu]=size(FSC{s});

n=find(abs(FL{s}{nassem{s}(3)}(:,u))>ciHld(s,3)); %Units involved in this factor
subplot(1,2,1), hold off cla
t=(1:ntp)*0.05+0.05/2;
plot(t,FSC{s}(:,u),'b',...  
     t,ciHsc(s,2),'r-','LineWidth',2);
   
set(gca,'FontSize',22), box off, xlim([680 700]);
title(['#units: ' num2str(n')]);

subplot(1,2,2), hold off cla
x=-3:0.01:3; h=histc(FSC{s}(:,u),x)./length(FSC{s}(:,u));
plot(x,h,'b',[ciHsc(s,2) ciHsc(s,2)],[0 0.05],'g','LineWidth',2);
% xlim([-0.5 1]);
set(gca,'FontSize',22), box off
title(['#units: ' num2str(n')]);

%% Choose units using Log-Likelihood
% NB load the factor score datasets...
[ntp,nu]=size(FSC{s});
alpha=0.05
s=1; % region
        z=diff(LL{s});
        Zbs=diff(LLbs{s}')';
        x=2:length(z)+1;
%         hold on, plot(Zbs','r');
        r=round(length(Zbs)*alpha); zz=sort(Zbs,'descend');
        %%%%%% number of significant assemblies
        n=x(find(x(diff(LL{s})> zz(r)),1,'last'))
        %%%%%%
        figure
        subplot(1,2,1), hold off cla, 
        plot(x,z,'b','LineWidth',2);
        hold on, plot(x,zz(r)*ones(1,length(z)),'r--','LineWidth',2); % plot Confidence limit
        set(gca,'FontSize',22), box off
%         title(['(' num2str(f) ',' num2str(s) '): ' num2str(nassem{s})]);
        xlabel('Factor (# assemblies)'), ylabel('log-likelihood-ratio')
        legend('original','trial-permutation'), legend('boxoff')
        subplot(1,2,2), hold off cla
        vs=cell2mat(FLbs{s}{2});
        x=-0.5:0.005:0.5; h=histc(vs(1:end),x)./length(vs(1:end));
        plot(x,h,'b',[ciHld(s,2) ciHld(s,2)],[0 0.05],'g','LineWidth',2);
        hold on, plot(FL{s}{2}(1:end),0.01,'r.'); xlim([min(x) max(x)]);
        set(gca,'FontSize',22), box off
        
%%
s=1;assem_no=2
factor_no=1;
iFR0{s}
(units{s}{factor_no})


plot(FSCsel{s}(:,assem_no))

%% re-wrap output data

% (1) load raw
times=load('JaroslawLONG2.mat');
% (2) load iFR50
% load('50ms_bins/JaroslawLONG2_PFC_iFR50.mat')
% (3) get aligned trial markers
[times.timeAxis,times.FiringRates,times.EventLabels,times.EventTimes,usel]=SelTrialsCells('KDE_bins/JaroslawLONG2_PFC_iFR50.mat',10,0.1,1e6,'iFR');

assemblies=load('KDE_bins/JaroslawLONG2_iFR50_FSC.mat');


% IDs of units going into FA {PFC HP}
data.input_units    = usel;	 
% IDs of units spat out by FA as Assembly members {PFC HP joint}
data.Assembly_units = assemblies.units; 

%re-register assembly units back to original numbering scheme
for region=1:3 % {PFC HP joint}
    for assem_no=1:length(assemblies.units{region})
        if region ~=3
            data.Assembly_units_adjusted{region}{assem_no}=data.input_units{region}(assemblies.units{region}{assem_no});
        else
            all_units=[ones(size(usel{1}')),usel{1}';2*ones(size(usel{2}')),usel{2}'];
            all_units=all_units(assemblies.units{region}{assem_no},:);
            data.Assembly_units_adjusted{region}{assem_no}{1}=all_units(all_units(:,1)==1,2); % PFC units in assembly
            data.Assembly_units_adjusted{region}{assem_no}{2}=all_units(all_units(:,1)==2,2); % HP units in assembly
        end     
    end
end; clear region assem_id all_units

for u_id=1:length(data.input_units{1})
    data.SpikeTimes{1}{u_id}=times.PFCcells{u_id}.t;
end
for u_id=1:length(data.input_units{2})
    data.SpikeTimes{2}{u_id}=times.PFCcells{u_id}.t;
end
data.SpikeTimes{3}=horzcat(data.SpikeTimes{1},data.SpikeTimes{2});

data.FactorScores=assemblies.FSCsel;
data.assemblies = times.timeAxis;

% find the peak assembly activation times
for region=1:3 % {PFC HP joint}
    for assembly_id=1:size(data.FactorScores{region},2);
        % factor scores
        this_FSC=data.FactorScores{region}(:,assembly_id); % plot(-5:0.1:10,histc(this_FSC,-5:0.1:10)) % histogram to choose activation threshold
        % find factor activation times
        [pks,locs]=findpeaks(this_FSC,'MINPEAKHEIGHT',1);
        
        assemblies.assembly_times{region}{assembly_id}=assemblies.Tmtx(locs); %*1e4 TO CONVERT BACK TO NLX TIME
    end
end
clear assem_no assembly_id locs pks region this_FSC u_id usel
%% Make raster of units that went into analysis pipeline (i.e. chosen by minimum firing rate etc. )
region = 2; % {PFC, HP, joint}
threshold = 1;
plot_win = [4000 8000];
% plot_win = [0 Inf];
colorscale=flipud(lines(size(data.FactorScores{region},2)));
lag=0.05;    
figure; 
subplot(2,1,1)
bar(times.EventTimes,times.EventLabels)
xlim(plot_win)

subplot(2,1,2); hold on
% plot ALL units as a spike raster
for unit_id=1:numel(data.input_units{region})
    this_unit=data.input_units{region}(unit_id);
    temp_firing=times.PFCcells{this_unit}.t/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
    scatter(temp_firing,0.1*unit_id*ones(length(temp_firing),1),10,'k','filled')
end

% highlight ONLY ASSEMBLY units    
for assembly_id=1:size(data.FactorScores{region},2);

    %%%% 
    % plot the factor activations
    % scatter(assembly_times,pks,5,colorscale(assembly_id,:),'filled')
	% plot(Tmtx,0.1*this_FSC,'color',colorscale(assembly_id,:))
    %%%%

    % assembly member units... has format data.input_units{brain region} and assemblies.units{brain region}{assembly number}
    factor_active_units = data.input_units{region}(assemblies.units{region}{assembly_id}); % ... has format data.input_units{brain region} and assemblies.units{brain region}{assembly number}
    for unit_id=1:numel(factor_active_units)
        this_unit=factor_active_units(unit_id);
        temp_firing=times.PFCcells{this_unit}.t/1e4;
        
        % Optionally exclude spikes outside plotting window
        % temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
        
        % highlight spike times of assembly units
        if ismember(this_unit,factor_active_units)
            for assem_peak_no=1:numel(assemblies.assembly_times{region}{assembly_id})
                temp_firing_ids=(temp_firing-lag)<=assemblies.assembly_times{region}{assembly_id}(assem_peak_no) &...
                                (temp_firing+lag)>=assemblies.assembly_times{region}{assembly_id}(assem_peak_no);
                scatter(temp_firing(temp_firing_ids),0.1*find(this_unit==data.input_units{region})*ones(length(temp_firing(temp_firing_ids)),1),40,colorscale(assembly_id,:),'filled')
            end
        end
    end
end

axis tight; xlim(plot_win)




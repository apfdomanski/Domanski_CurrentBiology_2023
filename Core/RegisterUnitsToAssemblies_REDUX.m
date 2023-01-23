% Line up assemblies to their constituent units
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% 2016: OBSOLETE!
%% re-wrap output data
% JaroslawLONG2
% KrzesimirLONG2
% KrzysztofLONG2
% MiroslawLONG2
% NorbertLONG1
Animal= 'KrzysztofLONG2';
% (1) load raw
cd('C:\Analysis\AssemblyAnalysis\raw')
times=load([Animal,'.mat']);

% (2) get aligned trial markers, firing rates and included unit IDs
[times.timeAxis,...
 times.FiringRates,...
 times.EventLabels,...
 times.EventTimes,...
 usel]                 = SelTrialsCells(['KDE_bins/',Animal,'_PFC_iFR50.mat'],10,0.1,1e6,'iFR','all');

% (3) load
assemblies=load(['KDE_bins/LONG/',Animal,'_iFR50_FSC.mat']);
thrsh_=load(['KDE_bins/LONG/',Animal,'_iFR50__FSCtemp.mat'],'ciHsc');
threshold = thrsh_.ciHsc(:,2);% Factor score theshold @0.05

% plot example FL loading results
% load(['KDE_bins/LONG/',Animal,'_iFR50_AssemRes2.mat']);
% figure
% subplot(1,3,1); imagesc(FL{1}{5}); axis off
% subplot(1,3,2); imagesc(FL{2}{3}); axis off
% subplot(1,3,3); imagesc(FL{3}{4}); axis off
% colormap hot; colorbar

%% IDs of units going into FA {PFC HP}
data.input_units    = usel;	 
% IDs of units spat out by FA as Assembly members {PFC HP joint}
data.Assembly_units = assemblies.units; 

%re-register assembly units back to original numbering scheme 
% (handle inter-area assemblies as secial case)
for region=1:3 % {PFC HP joint}
    if length(assemblies.units{region})>0
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
    end
end; clear region assem_id all_units
% collate the spike times
for u_id=1:length(data.input_units{1})
    try
        data.SpikeTimes{1}{u_id}=times.PFCcells{u_id}.t;
    catch,end
end
for u_id=1:length(data.input_units{2})
    try
        data.SpikeTimes{2}{u_id}=times.PFCcells{u_id}.t;
	catch,end
end
data.SpikeTimes{3}=horzcat(data.SpikeTimes{1},data.SpikeTimes{2});

data.FactorScores=assemblies.FSCsel;
data.assemblies = times.timeAxis;

%% find the peak assembly activation times
for region=1:3 % {PFC HP joint}
    if ~isempty(data.FactorScores{region})
        for assembly_id=1:size(data.FactorScores{region},2);
            % factor scores
            this_FSC=data.FactorScores{region}(:,assembly_id); % plot(-5:0.1:10,histc(this_FSC,-5:0.1:10)) % histogram to choose activation threshold
            % find factor activation times
            [pks,locs]=findpeaks(this_FSC,'MINPEAKHEIGHT',threshold(region));

            assemblies.assembly_times{region}{assembly_id}=assemblies.Tmtx(locs); %*1e4 TO CONVERT BACK TO NLX TIME
        end
    end
end
clear assem_no assembly_id locs pks region this_FSC u_id usel


IncludedUnits   = data.input_units;
AssemblyUnits   = data.Assembly_units_adjusted;
AssemblyUnitsSUBSET   = data.Assembly_units;

SpikeTimes      = data.SpikeTimes;
TimeAxis        = assemblies.Tmtx;
FactorScores    = data.FactorScores;
FactorTimes     = assemblies.assembly_times;
%% Ignore me
% ID's of neurons chosen (by minimum firing rate)to go into Factor Analysis:
% IncludedUnits{PFC, HP}


% ID's of units showing membership in each assembly:
% AssemblyUnits{PFC, HP, JOINT}{ASSEMBLY NUMBER}(UNIT NUMBERS)
% N.B. data.Assembly_units_adjusted{3} (i.e. Joint HP-PFC assembly units)
% have 2 cells within: {PFC, HP} unit IDs

% Spike timestamps (time in NLX format) of neurons chosen (by minimum firing rate)to go into Factor Analysis:
% SpikeTimes{PFC, HP, JOINT}{UNIT NUMBERS}

% Experiment time axis (e.g. for assembly co-activations)
% TimeAxis
% N.B. Time is in seconds: multiply by 1e4 to correct time to NLX format!

% time-resolved assembly activations (Factor scores):
% FactorScores{PFC, HP, JOINT}(time,Assembly no.)

% Times of assembly activation:
% FactorTimes{PFC, HP, JOINT}{Assembly no.}
% N.B. Time is in seconds: multiply by 1e4 to correct time to NLX format!

% %% Make raster of units higlighting those tagged as 
% region = 2; % {PFC, HP, joint}
% threshold = 1;
% plot_win = [4000 8500];
% % plot_win = [0 Inf];
% 
% colorscale=flipud(lines(size(data.FactorScores{region},2)));
% lag=0.05;    
% figure; 
% subplot(2,1,1)
% bar(times.EventTimes,times.EventLabels)
% xlim(plot_win)
% xlabel('time (s)')
% ylabel('event type')
% 
% subplot(2,1,2); hold on
% % plot ALL units as a spike raster
% for unit_id=1:numel(data.input_units{region})
%         this_unit=data.input_units{region}(unit_id);
%     switch region
%         case 1
%             temp_firing=times.PFCcells{this_unit}.t/1e4;
%         case 2
%             temp_firing=times.HPcells{this_unit}.t/1e4;
%         otherwise
%             % TO DO: plotting for joint HP-PFC dataset
%     end
%         temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
%         scatter(temp_firing,unit_id*ones(length(temp_firing),1),5,'k','filled')
% end
% 
% % highlight ONLY ASSEMBLY units    
% for assembly_id=1:size(data.FactorScores{region},2);
% 
%     %%%% 
%     % plot the factor activations
%     % scatter(assembly_times,pks,5,colorscale(assembly_id,:),'filled')
% 	% plot(Tmtx,0.1*this_FSC,'color',colorscale(assembly_id,:))
%     %%%%
% 
%     % assembly member units... has format data.input_units{brain region} and assemblies.units{brain region}{assembly number}
%     factor_active_units = data.input_units{region}(assemblies.units{region}{assembly_id}); % ... has format data.input_units{brain region} and assemblies.units{brain region}{assembly number}
%     for unit_id=1:numel(factor_active_units)
%         this_unit=factor_active_units(unit_id);
%         switch region
%             case 1
%                 temp_firing=times.PFCcells{this_unit}.t/1e4;
%             case 2
%                 temp_firing=times.HPcells{this_unit}.t/1e4;
%             otherwise
%                 % TO DO: plotting for joint HP-PFC dataset
%         end        
%         % Optionally exclude spikes outside plotting window
%         % temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
%         
%         % highlight spike times of assembly units
%         if ismember(this_unit,factor_active_units)
%             for assem_peak_no=1:numel(assemblies.assembly_times{region}{assembly_id})
%                 temp_firing_ids=(temp_firing-lag)<=assemblies.assembly_times{region}{assembly_id}(assem_peak_no) &...
%                                 (temp_firing+lag)>=assemblies.assembly_times{region}{assembly_id}(assem_peak_no);
%                 scatter(temp_firing(temp_firing_ids),find(this_unit==data.input_units{region})*ones(length(temp_firing(temp_firing_ids)),1),40,colorscale(assembly_id,:))%,'filled')
%             end
%         end
%     end
% end
% xlabel('time (s)')
% ylabel('Unit number')
% 
% axis tight;
% xlim(plot_win)
% clear assembly_id assem_peak_no region colorscale factor_active_units lag plot_win temp_firing temp_firing_ids this_unit threshold unit_id
% %% Make raster of units higlighting those tagged as assembly members
% threshold = 8;
% 
% % adjust the peak assembly activation times... according to 'threshold' parameter
% for region=1:3 % {PFC HP joint}
%     if ~isempty(FactorScores{region})
%         for assembly_id=1:size(FactorScores{region},2);
%             % factor scores
%             this_FSC=FactorScores{region}(:,assembly_id); % plot(-5:0.1:10,histc(this_FSC,-5:0.1:10)) % histogram to choose activation threshold
%             % find factor activation times
%             [pks,locs]=findpeaks(this_FSC,'MINPEAKHEIGHT',threshold);
% 
%            FactorTimes{region}{assembly_id}=TimeAxis(locs); %*1e4 TO CONVERT BACK TO NLX TIME
%         end
%     end
% end
% clear assem_no assembly_id locs pks  this_FSC u_id usel
% 
% region = 1; % {PFC, HP, joint}
% plot_win = [4000 8000];
% 
% colorscale=flipud(lines(size(FactorScores{region},2)));
% lag=0.05;    % window either side of peak activation to include spikes as part of assembly acivation
% figure; hold on
% 
% 
% % Firstly, plot ALL units as a spike raster
% for unit_id=1:numel(IncludedUnits{region})
%     temp_firing=SpikeTimes{region}{unit_id}/1e4;
%     temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
% 	scatter(temp_firing,unit_id*ones(length(temp_firing),1),5,'k','filled')
% end
% 
% % Secondly, highlight ONLY ASSEMBLY units    
% for assembly_id=1:size(FactorScores{region},2);
% 
%     %%%% 
%     % plot the factor activations:
% 	% plot(TimeAxis,FactorScores{region}(:,assembly_id),'color',colorscale(assembly_id,:))
%     %%%%
% 
%     factor_active_units = AssemblyUnitsSUBSET{region}{assembly_id};
%     for unit_id=1:numel(factor_active_units)
%         this_unit=factor_active_units(unit_id);
%         switch region
%             case 1
%                 temp_firing=SpikeTimes{region}{this_unit}/1e4;
%             case 2
%                 temp_firing=SpikeTimes{region}{this_unit}/1e4;
%             otherwise
%                 % TO DO: plotting for joint HP-PFC dataset
%         end        
%         % Optionally exclude spikes outside plotting window
%         % temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
%         
%         % highlight spike times of assembly units
%         if ismember(this_unit,factor_active_units)
%             for assem_peak_no=1:numel(FactorTimes{region}{assembly_id})
%                 temp_firing_ids=(temp_firing-lag)<=FactorTimes{region}{assembly_id}(assem_peak_no) &...
%                                 (temp_firing+lag)>=FactorTimes{region}{assembly_id}(assem_peak_no);
%                             try
%                                 % This gets a try/catch to guard against
%                                 % cases where a cell doesn't fire even
%                                 % though the assembly is active
%                                 scatter(temp_firing(temp_firing_ids),...
%                                         find(this_unit==IncludedUnits{region})*ones(length(temp_firing(temp_firing_ids)),1),40,colorscale(assembly_id,:),'filled')
%                             end
%                 % Vertical lines mark factor times             
%                 plot([FactorTimes{region}{assembly_id}(assem_peak_no) FactorTimes{region}{assembly_id}(assem_peak_no)],[0 numel(IncludedUnits{region})],'color',colorscale(assembly_id,:))
%             end
%         end
%     end
% end
% xlabel('time (s)')
% ylabel('Unit number')
% 
% axis tight;
% xlim(plot_win)
% clear assembly_id assem_peak_no region colorscale factor_active_units lag plot_win temp_firing temp_firing_ids this_unit threshold unit_id
%% (Plot 1) Make raster of units higlighting those tagged as assemblies... spanning all types
lag=0.05;           % Highlight spikes within this value of factor activation peak
assem_dot_size=30;
cmap_stretch = 1;
plot_win = [4910 4960];
% plot_win = [4362 4412];
HP_offset=30;

figure;hold on 
Choice_colour=[0.7 0.7 0.9];
Sample_colour=[0.7 0.9 0.7];
markerpos=[90 110];

% Add correct sample/choice times:
t_list = [times.EventTimes(times.EventLabels==1); times.EventTimes(times.EventLabels==2)];
plot([t_list,t_list],markerpos,'color',Choice_colour,'LineWidth',1.5)
for id=1:numel(t_list)
    text(t_list(id)-1,markerpos(1),'Choice','color',Choice_colour,'Rotation',90);
end
t_list = [times.EventTimes(times.EventLabels==3); times.EventTimes(times.EventLabels==4)];
plot([t_list,t_list],markerpos,'color',Sample_colour,'LineWidth',1.5)
for id=1:numel(t_list)
    text(t_list(id)-1,markerpos(1),'Sample','color',Sample_colour,'Rotation',90);
end

% axis tight;
xlim(plot_win)


% Add error sample/choice times:
% t_list = [times.EventTimes(times.EventLabels==10); times.EventTimes(times.EventLabels==20)];
% plot([t_list,t_list],[-20,100],'color',Choice_colour,'LineWidth',1.5)
% for id=1:numel(t_list)
%     text(t_list(id)-1,80,'Choice','color',Choice_colour,'Rotation',90);
% end
% t_list = [times.EventTimes(times.EventLabels==30); times.EventTimes(times.EventLabels==40)];
% plot([t_list,t_list],[-20,100],'color',Sample_colour,'LineWidth',1.5)
% for id=1:numel(t_list)
%     text(t_list(id)-1,80,'Sample','color',Sample_colour,'Rotation',90);
% end
 
 
 
%%%%% plot ALL units as a spike raster...
%mPFC first...
for unit_id=1:numel(data.input_units{1})
    this_unit=data.input_units{1}(unit_id);
    temp_firing=times.PFCcells{this_unit}.t/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
    scatter(temp_firing,unit_id*ones(length(temp_firing),1),5,[0 0 0],'o','filled')
end
unit_offset=unit_id+HP_offset;
%..then HP
for unit_id=1:numel(data.input_units{2})
    this_unit=data.input_units{2}(unit_id);
    temp_firing=times.HPcells{this_unit}.t/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
    scatter(temp_firing,unit_offset+unit_id*ones(length(temp_firing),1),5,[0 0 0],'o','filled')
end
%
%%%%%




%%%%% plot Factor scores
% mPFC
region = 1;
colorscale=(autumn(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);
    plot(assemblies.Tmtx,...
         data.FactorScores{region}(:,assembly_id)-10,'color',colorscale(assembly_id,:),'LineWidth',2)
end
% HP
region = 2;
colorscale=(winter(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);
    plot(assemblies.Tmtx,...
         unit_offset+data.FactorScores{region}(:,assembly_id)-18,'color',colorscale(assembly_id,:),'LineWidth',2)
end
%Joint
region = 3;
colorscale=flipud(bone(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);
    plot(assemblies.Tmtx,...
         data.FactorScores{region}(:,assembly_id)-30,'color',colorscale(assembly_id,:),'LineWidth',2)
end
%%%%%

% highlight ONLY ASSEMBLY units
% PFC 
region = 1;
colorscale=(autumn(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);    
     factor_active_units = data.input_units{region}(assemblies.units{region}{assembly_id}); % ... has format data.input_units{brain region} and assemblies.units{brain region}{assembly number}
     for unit_id=1:numel(factor_active_units)
        this_unit=factor_active_units(unit_id);
        u_temp(unit_id)=find(this_unit==IncludedUnits{region});
     end
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
                try
                    scatter(temp_firing(temp_firing_ids),find(this_unit==data.input_units{region})*ones(length(temp_firing(temp_firing_ids)),1),assem_dot_size,colorscale(assembly_id,:),'filled','MarkerEdgeColor','k')
                    
%                     % Vertical lines mark factor times          
%                     plot([assemblies.assembly_times{region}{assembly_id}(assem_peak_no) assemblies.assembly_times{region}{assembly_id}(assem_peak_no)],...
%                     [min(u_temp) max(u_temp)],'-',...
%                     'LineWidth',1.5,...
%                     'color',colorscale(assembly_id,:))%[0 numel(IncludedUnits{region})],'-','color',colorscale(assembly_id,:));
                catch
                    warning('no firing for this unit!')
                end
            end
        end
    end
end
% HP
region = 2;
colorscale=(winter(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);
     factor_active_units = data.input_units{region}(assemblies.units{region}{assembly_id}); % ... has format data.input_units{brain region} and assemblies.units{brain region}{assembly number}
     for unit_id=1:numel(factor_active_units)
            this_unit=factor_active_units(unit_id);
            u_temp(unit_id)=find(this_unit==IncludedUnits{region});
     end
     for unit_id=1:numel(factor_active_units)
        this_unit=factor_active_units(unit_id);
        temp_firing=times.PFCcells{this_unit}.t/1e4;
        % temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));         % Optionally exclude spikes outside plotting window

        % highlight spike times of assembly units
        if ismember(this_unit,factor_active_units)
            for assem_peak_no=1:numel(assemblies.assembly_times{region}{assembly_id})
                temp_firing_ids=(temp_firing-lag)<=assemblies.assembly_times{region}{assembly_id}(assem_peak_no) &...
                    (temp_firing+lag)>=assemblies.assembly_times{region}{assembly_id}(assem_peak_no);
                scatter(temp_firing(temp_firing_ids),unit_offset+find(this_unit==data.input_units{region})*ones(length(temp_firing(temp_firing_ids)),1),assem_dot_size,colorscale(assembly_id,:),'filled','MarkerEdgeColor','k')

%             % Vertical lines mark factor times          
%             plot([assemblies.assembly_times{region}{assembly_id}(assem_peak_no) assemblies.assembly_times{region}{assembly_id}(assem_peak_no)],...
%             [unit_offset+min(u_temp) unit_offset+max(u_temp)],'-',...
%             'LineWidth',1.5,...
%             'color',colorscale(assembly_id,:))%[0 numel(IncludedUnits{region})],'-','color',colorscale(assembly_id,:));
            end
        end
        
    end
end
% Joint
colorscale=flipud(bone(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{3},2);
    factor_active_units_PFC = data.Assembly_units_adjusted{3}{assembly_id}{1};%data.input_units{1}(data.Assembly_units_adjusted{3}{assembly_id}{1});
    factor_active_units_HP  = data.Assembly_units_adjusted{3}{assembly_id}{2};%data.input_units{2}(data.Assembly_units_adjusted{3}{assembly_id}{2});
    
    %      for unit_id=1:numel(factor_active_units)
    %             this_unit=factor_active_units(unit_id);
    %             u_temp(unit_id)=find(this_unit==IncludedUnits{region});
    %      end
    if numel(factor_active_units_PFC)>0
        for unit_id=1:numel(factor_active_units_PFC)
            this_unit=factor_active_units_PFC(unit_id);
            temp_firing=times.PFCcells{this_unit}.t/1e4;
            
            % highlight spike times of assembly units
            if ismember(this_unit,factor_active_units_PFC)
                for assem_peak_no=1:numel(assemblies.assembly_times{3}{assembly_id})
                    temp_firing_ids=(temp_firing-lag)<=assemblies.assembly_times{3}{assembly_id}(assem_peak_no) &...
                                    (temp_firing+lag)>=assemblies.assembly_times{3}{assembly_id}(assem_peak_no);
                    try
                        scatter(temp_firing(temp_firing_ids),find(this_unit==data.input_units{1})*ones(length(temp_firing(temp_firing_ids)),1),assem_dot_size,colorscale(assembly_id,:),'filled','MarkerEdgeColor','k')
                    catch
                        warning('no firing for this unit!')
                    end
                    %             % Vertical lines mark factor times
                    %             plot([assemblies.assembly_times{region}{assembly_id}(assem_peak_no) assemblies.assembly_times{region}{assembly_id}(assem_peak_no)],...
                    %             [unit_offset+min(u_temp) unit_offset+max(u_temp)],'-',...
                    %             'LineWidth',1.5,...
                    %             'color',colorscale(assembly_id,:))%[0 numel(IncludedUnits{region})],'-','color',colorscale(assembly_id,:));
                end
            end
        end
    end
    if numel(factor_active_units_HP)>0
        for unit_id=1:numel(factor_active_units_HP)
            this_unit=factor_active_units_HP(unit_id);
            temp_firing=times.HPcells{this_unit}.t/1e4;
            
            % highlight spike times of assembly units
            if ismember(this_unit,factor_active_units_HP)
                for assem_peak_no=1:numel(assemblies.assembly_times{3}{assembly_id})
                    temp_firing_ids=(temp_firing-lag)<=assemblies.assembly_times{3}{assembly_id}(assem_peak_no) &...
                                    (temp_firing+lag)>=assemblies.assembly_times{3}{assembly_id}(assem_peak_no);
                    try
                        scatter(temp_firing(temp_firing_ids),unit_offset+find(this_unit==data.input_units{2})*ones(length(temp_firing(temp_firing_ids)),1),assem_dot_size,colorscale(assembly_id,:),'filled','MarkerEdgeColor','k')
                    catch
                        warning('no firing for this unit!')
                    end
                    %             % Vertical lines mark factor times
                    %             plot([assemblies.assembly_times{region}{assembly_id}(assem_peak_no) assemblies.assembly_times{region}{assembly_id}(assem_peak_no)],...
                    %             [unit_offset+min(u_temp) unit_offset+max(u_temp)],'-',...
                    %             'LineWidth',1.5,...
                    %             'color',colorscale(assembly_id,:))%[0 numel(IncludedUnits{region})],'-','color',colorscale(assembly_id,:));
                end
            end
        end
    end
end


  
xlabel('time (s)')
ylabel('Unit number')
clear assembly_id assem_peak_no region colorscale factor_active_units lag plot_win temp_firing temp_firing_ids this_unit threshold unit_id

% axis off
%% (Plot 2) Make raster of units higlighting those tagged as assemblies: mean activation

lag=0.05;           % Highlight spikes within this value of factor activation peak
assem_dot_size=30;
cmap_stretch = 1;
% plot_win = [4910 4960];
% plot_win = [4362 4412];
 plot_win = [4210 8500];
HP_offset=30;

figure;hold on 
Choice_colour=[0.7 0.7 0.9];
Sample_colour=[0.7 0.9 0.7];
markerpos=[90 110];
%% Add correct sample/choice times:
t_list = [times.EventTimes(times.EventLabels==1); times.EventTimes(times.EventLabels==2)];
plot([t_list,t_list],markerpos,'color',Choice_colour,'LineWidth',1.5)
for id=1:numel(t_list)
    text(t_list(id)-1,markerpos(1),'Choice','color',Choice_colour,'Rotation',90);
end
t_list = [times.EventTimes(times.EventLabels==3); times.EventTimes(times.EventLabels==4)];
plot([t_list,t_list],markerpos,'color',Sample_colour,'LineWidth',1.5)
for id=1:numel(t_list)
    text(t_list(id)-1,markerpos(1),'Sample','color',Sample_colour,'Rotation',90);
end

% axis tight;
xlim(plot_win)
%% plot ALL units as a spike raster...
% mPFC first...
for unit_id=1:numel(data.input_units{1})
    this_unit=data.input_units{1}(unit_id);
    temp_firing=times.PFCcells{this_unit}.t/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
    scatter(temp_firing,unit_id*ones(length(temp_firing),1),5,[0 0 0],'o','filled')
end
unit_offset=unit_id+HP_offset;
%..then HP
for unit_id=1:numel(data.input_units{2})
    this_unit=data.input_units{2}(unit_id);
    temp_firing=times.HPcells{this_unit}.t/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
    scatter(temp_firing,unit_offset+unit_id*ones(length(temp_firing),1),5,[0 0 0],'o','filled')
end
%
%%%%%
%% plot Factor scores
% mPFC
region = 1;
colorscale=(autumn(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);
    plot(assemblies.Tmtx,...
         data.FactorScores{region}(:,assembly_id)-10,'color',colorscale(assembly_id,:),'LineWidth',2)
end
% HP
region = 2;
colorscale=(winter(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);
    plot(assemblies.Tmtx,...
         unit_offset+data.FactorScores{region}(:,assembly_id)-18,'color',colorscale(assembly_id,:),'LineWidth',2)
end
%Joint
region = 3;
colorscale=flipud(bone(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);
    plot(assemblies.Tmtx,...
         data.FactorScores{region}(:,assembly_id)-30,'color',colorscale(assembly_id,:),'LineWidth',2)
end
%% highlight ONLY ASSEMBLY units
% PFC 
region = 1;
colorscale=(autumn(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);    
     factor_active_units = data.input_units{region}(assemblies.units{region}{assembly_id}); % ... has format data.input_units{brain region} and assemblies.units{brain region}{assembly number}
     for unit_id=1:numel(factor_active_units)
        this_unit=factor_active_units(unit_id);
        u_temp(unit_id)=find(this_unit==IncludedUnits{region});
     end
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
                try
                    scatter(temp_firing(temp_firing_ids),find(this_unit==data.input_units{region})*ones(length(temp_firing(temp_firing_ids)),1),assem_dot_size,colorscale(assembly_id,:),'filled','MarkerEdgeColor','k')
                    
%                     % Vertical lines mark factor times          
%                     plot([assemblies.assembly_times{region}{assembly_id}(assem_peak_no) assemblies.assembly_times{region}{assembly_id}(assem_peak_no)],...
%                     [min(u_temp) max(u_temp)],'-',...
%                     'LineWidth',1.5,...
%                     'color',colorscale(assembly_id,:))%[0 numel(IncludedUnits{region})],'-','color',colorscale(assembly_id,:));
                catch
                    warning('no firing for this unit!')
                end
            end
        end
    end
end
% HP
region = 2;
colorscale=(winter(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{region},2);
     factor_active_units = data.input_units{region}(assemblies.units{region}{assembly_id}); % ... has format data.input_units{brain region} and assemblies.units{brain region}{assembly number}
     for unit_id=1:numel(factor_active_units)
            this_unit=factor_active_units(unit_id);
            u_temp(unit_id)=find(this_unit==IncludedUnits{region});
     end
     for unit_id=1:numel(factor_active_units)
        this_unit=factor_active_units(unit_id);
        temp_firing=times.PFCcells{this_unit}.t/1e4;
        % temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));         % Optionally exclude spikes outside plotting window

        % highlight spike times of assembly units
        if ismember(this_unit,factor_active_units)
            for assem_peak_no=1:numel(assemblies.assembly_times{region}{assembly_id})
                temp_firing_ids=(temp_firing-lag)<=assemblies.assembly_times{region}{assembly_id}(assem_peak_no) &...
                    (temp_firing+lag)>=assemblies.assembly_times{region}{assembly_id}(assem_peak_no);
                scatter(temp_firing(temp_firing_ids),unit_offset+find(this_unit==data.input_units{region})*ones(length(temp_firing(temp_firing_ids)),1),assem_dot_size,colorscale(assembly_id,:),'filled','MarkerEdgeColor','k')

%             % Vertical lines mark factor times          
%             plot([assemblies.assembly_times{region}{assembly_id}(assem_peak_no) assemblies.assembly_times{region}{assembly_id}(assem_peak_no)],...
%             [unit_offset+min(u_temp) unit_offset+max(u_temp)],'-',...
%             'LineWidth',1.5,...
%             'color',colorscale(assembly_id,:))%[0 numel(IncludedUnits{region})],'-','color',colorscale(assembly_id,:));
            end
        end
        
    end
end
% Joint
colorscale=flipud(bone(cmap_stretch*size(data.FactorScores{region},2)));
for assembly_id=1:size(data.FactorScores{3},2);
    factor_active_units_PFC = data.Assembly_units_adjusted{3}{assembly_id}{1};%data.input_units{1}(data.Assembly_units_adjusted{3}{assembly_id}{1});
    factor_active_units_HP  = data.Assembly_units_adjusted{3}{assembly_id}{2};%data.input_units{2}(data.Assembly_units_adjusted{3}{assembly_id}{2});
    
    %      for unit_id=1:numel(factor_active_units)
    %             this_unit=factor_active_units(unit_id);
    %             u_temp(unit_id)=find(this_unit==IncludedUnits{region});
    %      end
    if numel(factor_active_units_PFC)>0
        for unit_id=1:numel(factor_active_units_PFC)
            this_unit=factor_active_units_PFC(unit_id);
            temp_firing=times.PFCcells{this_unit}.t/1e4;
            
            % highlight spike times of assembly units
            if ismember(this_unit,factor_active_units_PFC)
                for assem_peak_no=1:numel(assemblies.assembly_times{3}{assembly_id})
                    temp_firing_ids=(temp_firing-lag)<=assemblies.assembly_times{3}{assembly_id}(assem_peak_no) &...
                                    (temp_firing+lag)>=assemblies.assembly_times{3}{assembly_id}(assem_peak_no);
                    try
                        scatter(temp_firing(temp_firing_ids),find(this_unit==data.input_units{1})*ones(length(temp_firing(temp_firing_ids)),1),assem_dot_size,colorscale(assembly_id,:),'filled','MarkerEdgeColor','k')
                    catch
                        warning('no firing for this unit!')
                    end
                    %             % Vertical lines mark factor times
                    %             plot([assemblies.assembly_times{region}{assembly_id}(assem_peak_no) assemblies.assembly_times{region}{assembly_id}(assem_peak_no)],...
                    %             [unit_offset+min(u_temp) unit_offset+max(u_temp)],'-',...
                    %             'LineWidth',1.5,...
                    %             'color',colorscale(assembly_id,:))%[0 numel(IncludedUnits{region})],'-','color',colorscale(assembly_id,:));
                end
            end
        end
    end
    if numel(factor_active_units_HP)>0
        for unit_id=1:numel(factor_active_units_HP)
            this_unit=factor_active_units_HP(unit_id);
            temp_firing=times.HPcells{this_unit}.t/1e4;
            
            % highlight spike times of assembly units
            if ismember(this_unit,factor_active_units_HP)
                for assem_peak_no=1:numel(assemblies.assembly_times{3}{assembly_id})
                    temp_firing_ids=(temp_firing-lag)<=assemblies.assembly_times{3}{assembly_id}(assem_peak_no) &...
                                    (temp_firing+lag)>=assemblies.assembly_times{3}{assembly_id}(assem_peak_no);
                    try
                        scatter(temp_firing(temp_firing_ids),unit_offset+find(this_unit==data.input_units{2})*ones(length(temp_firing(temp_firing_ids)),1),assem_dot_size,colorscale(assembly_id,:),'filled','MarkerEdgeColor','k')
                    catch
                        warning('no firing for this unit!')
                    end
                    %             % Vertical lines mark factor times
                    %             plot([assemblies.assembly_times{region}{assembly_id}(assem_peak_no) assemblies.assembly_times{region}{assembly_id}(assem_peak_no)],...
                    %             [unit_offset+min(u_temp) unit_offset+max(u_temp)],'-',...
                    %             'LineWidth',1.5,...
                    %             'color',colorscale(assembly_id,:))%[0 numel(IncludedUnits{region})],'-','color',colorscale(assembly_id,:));
                end
            end
        end
    end
end


  
xlabel('time (s)')
ylabel('Unit number')
clear assembly_id assem_peak_no region colorscale factor_active_units lag plot_win temp_firing temp_firing_ids this_unit threshold unit_id

% axis off
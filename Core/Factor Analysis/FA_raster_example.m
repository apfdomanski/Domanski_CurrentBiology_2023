% Make raster of units higlighting those tagged as assembly members
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk

threshold = 5;

% adjust the peak assembly activation times... according to 'threshold' parameter
for region=1:3 % {PFC HP joint}
    if ~isempty(FactorScores{region})
        for assembly_id=1:size(FactorScores{region},2);
            % factor scores
            this_FSC=FactorScores{region}(:,assembly_id); % plot(-5:0.1:10,histc(this_FSC,-5:0.1:10)) % histogram to choose activation threshold
            % find factor activation times
            [pks,locs]=findpeaks(this_FSC,'MINPEAKHEIGHT',threshold);

           FactorTimes{region}{assembly_id}=TimeAxis(locs); %*1e4 TO CONVERT BACK TO NLX TIME
        end
    end
end
clear assem_no assembly_id locs pks  this_FSC u_id usel

region = 1; % {PFC, HP, joint}
plot_win = [4000 8000];
colorscale=flipud(lines(size(FactorScores{region},2)));
lag=0.05;    % window either side of peak activation to include spikes as part of assembly acivation

figure; hold on
% Firstly, plot ALL units as a spike raster
for unit_id=1:numel(IncludedUnits{region})
    temp_firing=SpikeTimes{region}{unit_id}/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
	scatter(temp_firing,unit_id*ones(length(temp_firing),1),5,[0.6 0.6 0.6],'filled')
end

% Secondly, highlight ONLY ASSEMBLY units    
for assembly_id=1:size(FactorScores{region},2);
    %%%% 
    % plot the factor activations:
	% plot(TimeAxis,FactorScores{region}(:,assembly_id),'color',colorscale(assembly_id,:))
    %%%%

    factor_active_units = AssemblyUnitsSUBSET{region}{assembly_id};
        for unit_id=1:numel(factor_active_units)
            this_unit=factor_active_units(unit_id);
            u_temp(unit_id)=find(this_unit==IncludedUnits{region});
        end
	for unit_id=1:numel(factor_active_units)
        this_unit=factor_active_units(unit_id);
        switch region
            case 1
                temp_firing=SpikeTimes{region}{this_unit}/1e4;
            case 2
                temp_firing=SpikeTimes{region}{this_unit}/1e4;
            otherwise
                % TO DO: plotting for joint HP-PFC dataset
        end        
        % Optionally exclude spikes outside plotting window
        % temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
        
        % highlight spike times of assembly units
        if ismember(this_unit,factor_active_units)
            for assem_peak_no=1:numel(FactorTimes{region}{assembly_id})
                temp_firing_ids=(temp_firing-lag)<=FactorTimes{region}{assembly_id}(assem_peak_no) &...
                                (temp_firing+lag)>=FactorTimes{region}{assembly_id}(assem_peak_no);
                            try
                                % This gets a try/catch to guard against cases where a cell doesn't fire eventhough the assembly is active
                                scatter(temp_firing(temp_firing_ids),...
                                        find(this_unit==IncludedUnits{region})*ones(length(temp_firing(temp_firing_ids)),1),40,colorscale(assembly_id,:),'filled')
                            end
                % Vertical lines mark factor times      
                
                plot([FactorTimes{region}{assembly_id}(assem_peak_no) FactorTimes{region}{assembly_id}(assem_peak_no)],...
                [min(u_temp) max(u_temp)],'-',...
                'LineWidth',1.5,...
                'color',colorscale(assembly_id,:))%[0 numel(IncludedUnits{region})],'-','color',colorscale(assembly_id,:));
            end
        end
    end
end
xlabel('time (s)')
ylabel('Unit number')

axis tight;
xlim(plot_win)
% clear ans assembly_id assem_peak_no region colorscale factor_active_units lag plot_win temp_firing temp_firing_ids this_unit threshold unit_id


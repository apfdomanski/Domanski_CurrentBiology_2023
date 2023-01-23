function AssemblyRasterPlotAll(Factor,Rates,Units,Events,plot_win)
% Make raster of units higlighting those tagged as assembly members
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk

threshold = 3;
colorscale{1} = zeros(size(Factor.FSC{1},2),3);
colorscale{1}(:,3) = 1-0.5*rand(size(Factor.FSC{1},2),1);

colorscale{2} = zeros(size(Factor.FSC{2},2),3);
colorscale{2}(:,1) = 1-0.5*rand(size(Factor.FSC{2},2),1);

colorscale{3} = zeros(size(Factor.FSC{3},2),3);
colorscale{3}(:,2) = 1-0.5*rand(size(Factor.FSC{3},2),1);
lag        = 1;    % window either side of peak activation to include spikes as part of assembly activation

% Extract assembly peak times
for region = 1:3 % {PFC HP joint}
    if ~isempty(Factor.FSC{region})
        for assembly_id=1:size(Factor.FSC{region},2)
            % factor scores
            this_FSC=Factor.FSC{region}(:,assembly_id); % plot(-5:0.1:10,histc(this_FSC,-5:0.1:10)) % histogram to choose activation threshold
            % find factor activation times. NB ciHsc  p-values are [0.1 0.05 0.01 5e-3 1e-3];
            [~,locs]=findpeaks(this_FSC,'MINPEAKHEIGHT',Factor.ciHsc(region,threshold));
            FactorTimes{region}{assembly_id}=Factor.TmtxS(locs); %*1e4 TO CONVERT BACK TO NLX TIME
            FactorTimes{region}{assembly_id}(FactorTimes{region}{assembly_id}<plot_win(1) | FactorTimes{region}{assembly_id}>plot_win(2))=[];
        end
    end
end
clear assem_no assembly_id locs pks  this_FSC u_id usel



%% (1) plot ALL units as a spike raster
figure('color','w'); hold on
region = 3;
for unit_id=1:numel(Units{region})
    temp_firing=Units{region}{unit_id}.t/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
    switch region
        case 3 % Different colours for different areas
            if unit_id<=length(Factor.unitIDs{1})
                %                 scatter(temp_firing,unit_id*ones(length(temp_firing),1),1,[0.9 0.6 0.6],'.')
                plot(([temp_firing,temp_firing])',(unit_id+[-0.5 0.5].*ones(length(temp_firing),2))','color',[0.6 0.6 0.6])%[0.9 0.6 0.6])
            else
                %                 scatter(temp_firing,unit_id*ones(length(temp_firing),1),1,[0.6 0.6 0.9],'.')
                
                plot(([temp_firing,temp_firing])',(unit_id+[-0.5 0.5].*ones(length(temp_firing),2))','color',[0.6 0.6 0.6])%[0.6 0.6 0.9])
            end
            plot([plot_win(1) plot_win(1)],[0 length(Factor.unitIDs{1})+0.5],'LineWidth',2,'Color',[0 0 1])
            plot([plot_win(1) plot_win(1)],[length(Factor.unitIDs{1})+0.5 length(Factor.unitIDs{3})+0.5],'LineWidth',2,'Color',[1 0 0])
        otherwise
            %             scatter(temp_firing,unit_id*ones(length(temp_firing),1),1,[0.6 0.6 0.6],'.')
            plot(([temp_firing,temp_firing])',(unit_id+[-0.5 0.5].*ones(length(temp_firing),2))','color',[0.6 0.6 0.6])
            
    end
end
%% (2) Plot the assembly traces
for region = 1:3

    for assembly_id=1:size(Factor.FSC{region},2)
        %%%%
        % plot the factor activations:
        idx = Factor.TmtxS>=plot_win(1) & Factor.TmtxS<=plot_win(2);
        plot(Factor.TmtxS(idx),Factor.FSC{region}(idx,assembly_id)-assembly_id/2-10*region,'color',colorscale{region}(assembly_id,:),'LineWidth',1)
        %%%%
    end
    
end
%% (3) plot assembly peak times
% for assembly_id=1:size(Factor.FSC{region},2)
%     if ~isempty(FactorTimes{region}{assembly_id})
% %     plot([FactorTimes{region}{assembly_id}; FactorTimes{region}{assembly_id}],[min(u_temp) max(u_temp)],...
% %         '-','LineWidth',1.5,...
% %         'color',colorscale(assembly_id,:))
%      plot([FactorTimes{region}{assembly_id}; FactorTimes{region}{assembly_id}],[-1 -4],...
%         '-','LineWidth',1.5,...
%         'color',colorscale(assembly_id,:))
%     end
% end
%% (4) highlight ONLY ASSEMBLY units
for region = 1:3
    for assembly_id=1:size(Factor.FSC{region},2)

        factor_active_units = Factor.units{region}{assembly_id};

        %Realign unit numbers
    %     for unit_id=1:numel(factor_active_units)
    %         this_unit=factor_active_units(unit_id);
    %         %u_temp(unit_id)=Factor.unitIDs{region}(this_unit);
    %         u_temp(unit_id)=this_unit;
    %     end
        for unit_id=1:numel(factor_active_units)
            this_unit=factor_active_units(unit_id);
            if region ==2
                this_unit=this_unit+length(Units{1});
            end
            switch region
                case 3
                    temp_firing=Units{region}{this_unit}.t/1e4;
                otherwise
                    temp_firing=Units{region}{this_unit}.t/1e4;
                    % TO DO: plotting for joint HP-PFC dataset
            end

            % highlight spike times of assembly units
            for assem_peak_no=1:numel(FactorTimes{region}{assembly_id})
                temp_firing_ids=(temp_firing-lag)<=FactorTimes{region}{assembly_id}(assem_peak_no) &...
                    (temp_firing+lag)>=FactorTimes{region}{assembly_id}(assem_peak_no) &...
                    (temp_firing-lag)>=plot_win(1) &...
                    (temp_firing-lag)<=plot_win(2);

                try  % This gets a try/catch to guard against cases where a cell doesn't fire eventhough the assembly is active

    %                 scatter(temp_firing(temp_firing_ids),...
    %                     find(this_unit==Factor.unitIDs{region})*ones(length(temp_firing(temp_firing_ids)),1),50,colorscale(assembly_id,:),'s','filled')
                      scatter(temp_firing(temp_firing_ids),...
                        this_unit*ones(length(temp_firing(temp_firing_ids)),1),50,colorscale{region}(assembly_id,:),'s','filled')
                end
            end
        end
    end
end
%% (5) overlay lever press times
% Collapse lever press events
TSample = Events.EvtTs(Events.EvtLs>2);
TChoice = Events.EvtTs(Events.EvtLs<=2);
Outcome = Events.evt0;
Y = length(Units{3})+1;
for iLever = 1:length(TSample)
    if [TSample(iLever) TChoice(iLever)] > plot_win(1) & [TSample(iLever) TChoice(iLever)] < plot_win(2)
        plot([TSample(iLever) TSample(iLever)],[Y Y+2],'-g','LineWidth',2)
        plot([TChoice(iLever) TChoice(iLever)],[Y Y+2],'-r','LineWidth',2)
        text(TSample(iLever),Y+3,'Sample','Rotation',90,'Color','g')
        text(TChoice(iLever),Y+3,'Choice','Rotation',90,'Color','r')
        switch Outcome(iLever)
            case 1
                text((TSample(iLever) + TChoice(iLever))/2,Y+3,'L','Color','k')
            case 2
                text((TSample(iLever) + TChoice(iLever))/2,Y+3,'R','Color','k')
        end
    end
end
%%
xlabel('time (s)')
ylabel('Unit number')
axis tight;
xlim(plot_win)

axis off

text(plot_win(1)-2,length(Factor.unitIDs{1})/2,'mPFC Units','Rotation',90,'Color','b','HorizontalAlignment','center')
text(plot_win(1)-2,length(Factor.unitIDs{1})+(length(Factor.unitIDs{3})-length(Factor.unitIDs{1}))/2,'dCA1 Units','Rotation',90,'Color','r','HorizontalAlignment','center')
plot([plot_win(1)+0 plot_win(1)+20],[-8 -8],'LineWidth',2,'Color','k')
text(plot_win(1)+10,-6,'10s','HorizontalAlignment','center')

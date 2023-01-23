function AssemblyRasterPlotDiscreteJointOnly(Factor,Rates,Units,Events,plot_win)
% Make raster of units higlighting those tagged as assembly members
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk

threshold = 3; % Which level of scrutiny to apply to thresholding
colorscale{1} = zeros(size(Factor.FSC{1},2),3);
colorscale{1}(:,3) = 1-0.5*rand(size(Factor.FSC{1},2),1);

colorscale{2} = zeros(size(Factor.FSC{2},2),3);
colorscale{2}(:,1) = 1-0.5*rand(size(Factor.FSC{2},2),1);

colorscale{3} = zeros(size(Factor.FSC{3},2),3);
colorscale{3}(:,2) = 1-0.5*rand(size(Factor.FSC{3},2),1);

colorscale{3} = jet(size(Factor.FSC{3},2));

lag        = 0.25;    % window either side of peak activation to include spikes as part of assembly activation

% Extract assembly peak times
for iArea = 1:3 % {PFC HP joint}
    if ~isempty(Factor.FSC{iArea})
        for iAss=1:size(Factor.FSC{iArea},2)
            % factor scores
            this_FSC=Factor.FSC{iArea}(:,iAss); % plot(-5:0.1:10,histc(this_FSC,-5:0.1:10)) % histogram to choose activation threshold
            % find factor activation times. NB ciHsc  p-values are [0.1 0.05 0.01 5e-3 1e-3];
            [~,locs]=findpeaks(this_FSC,'MINPEAKHEIGHT',Factor.ciHsc(iArea,threshold));
            FactorTimes{iArea}{iAss}=Factor.TmtxS(locs); %*1e4 TO CONVERT BACK TO NLX TIME
            FactorTimes{iArea}{iAss}(FactorTimes{iArea}{iAss}<plot_win(1) | FactorTimes{iArea}{iAss}>plot_win(2))=[];
        end
    end
end
clear assem_no iAss locs pks  this_FSC u_id usel



%% (1) plot ALL units as a spike raster
figure('color','w'); hold on
iArea = 3;
for unit_id=1:numel(Units{iArea})
    temp_firing=Units{iArea}{unit_id}.t/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
    switch iArea
        case 3 % Different colours for different areas
            if unit_id<=length(Factor.unitIDs{1}) 
                plot(([temp_firing,temp_firing])',...
                     bsxfun(@times,unit_id+[-0.5; 0.5],ones(length(temp_firing),2)'),...
                     'color',[0.6 0.6 0.6])%[0.9 0.6 0.6])
            else
                plot(([temp_firing,temp_firing])',...
                    bsxfun(@times,unit_id+[-0.5; 0.5],ones(length(temp_firing),2)'),...
                     'color',[0.6 0.6 0.6])%[0.6 0.6 0.9])
            end
        otherwise
            %             scatter(temp_firing,unit_id*ones(length(temp_firing),1),1,[0.6 0.6 0.6],'.')
            plot(([temp_firing,temp_firing])',(unit_id+[-0.5 0.5].*ones(length(temp_firing),2))','color',[0.6 0.6 0.6])
            
    end
end
%% Area markers
switch iArea
        case 3 % Different colours for different areas
            
        plot([plot_win(1) plot_win(1)],[0 length(Factor.unitIDs{1})+0.5],'LineWidth',2,'Color',[0 0 1])
            plot([plot_win(1) plot_win(1)],[length(Factor.unitIDs{1})+0.5 length(Factor.unitIDs{3})+0.5],'LineWidth',2,'Color',[1 0 0])
    otherwise        
end
%% (2) Plot the assembly traces
for iArea = 3%1:3

    for iAss=1:size(Factor.FSC{iArea},2)
        %%%%
        % plot the factor activations:
        idx = Factor.TmtxS>=plot_win(1) & Factor.TmtxS<=plot_win(2);
        plot(Factor.TmtxS(idx),Factor.FSC{iArea}(idx,iAss)-3*iArea,'color',colorscale{iArea}(iAss,:),'LineWidth',1.5)
        %%%%
    end
    
end
%% (3) plot assembly peak times
% for iAss=1:size(Factor.FSC{iArea},2)
%     if ~isempty(FactorTimes{iArea}{iAss})
% %     plot([FactorTimes{iArea}{iAss}; FactorTimes{iArea}{iAss}],[min(u_temp) max(u_temp)],...
% %         '-','LineWidth',1.5,...
% %         'color',colorscale(iAss,:))
%      plot([FactorTimes{iArea}{iAss}; FactorTimes{iArea}{iAss}],[-1 -4],...
%         '-','LineWidth',1.5,...
%         'color',colorscale(iAss,:))
%     end
% end
%% (4) highlight ONLY ASSEMBLY units
for iArea = 1:3
    for iAss=1:size(Factor.FSC{iArea},2)

        factor_active_units = Factor.units{iArea}{iAss};

        %Realign unit numbers
    %     for unit_id=1:numel(factor_active_units)
    %         this_unit=factor_active_units(unit_id);
    %         %u_temp(unit_id)=Factor.unitIDs{iArea}(this_unit);
    %         u_temp(unit_id)=this_unit;
    %     end
        for unit_id=1:numel(factor_active_units)
            this_unit=factor_active_units(unit_id);

            switch iArea
                case 3
                    temp_firing=Units{iArea}{this_unit}.t/1e4;
                otherwise
                    temp_firing=Units{iArea}{this_unit}.t/1e4;
                    % TO DO: plotting for joint HP-PFC dataset
            end
            if iArea ==2
                this_unit=this_unit+length(Units{1});
            end
            % highlight spike times of assembly units
            for assem_peak_no=1:numel(FactorTimes{iArea}{iAss})
                temp_firing_ids=(temp_firing-lag)<=FactorTimes{iArea}{iAss}(assem_peak_no) &...
                    (temp_firing+lag)>=FactorTimes{iArea}{iAss}(assem_peak_no) &...
                    (temp_firing-lag)>=plot_win(1) &...
                    (temp_firing-lag)<=plot_win(2);

                try  % This gets a try/catch to guard against cases where a cell doesn't fire eventhough the assembly is active

    %                 scatter(temp_firing(temp_firing_ids),...
    %                     find(this_unit==Factor.unitIDs{iArea})*ones(length(temp_firing(temp_firing_ids)),1),50,colorscale(iAss,:),'s','filled')
                      scatter(temp_firing(temp_firing_ids),...
                        this_unit*ones(length(temp_firing(temp_firing_ids)),1),50,colorscale{iArea}(iAss,:),'s','filled')
                end
            end
        end
    end
end
%% (5) overlay lever press times
% Collapse lever press events
TSample = [Events.Left(:,1);Events.Right(:,1)]/1e6;
TChoice = [Events.Left(:,2);Events.Right(:,2)]/1e6;
Outcome = [ones(size(Events.Left,1),1);2*ones(size(Events.Right,1),1)];
% Y = 3;
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

xlim(plot_win)



text(plot_win(1)-2,length(Factor.unitIDs{1})/2,'mPFC Units','Rotation',90,'Color','b','HorizontalAlignment','center')
text(plot_win(1)-2,length(Factor.unitIDs{1})+(length(Factor.unitIDs{3})-length(Factor.unitIDs{1}))/2,'dCA1 Units','Rotation',90,'Color','r','HorizontalAlignment','center')
plot([plot_win(1)+0 plot_win(1)+20],[-15 -15],'LineWidth',2,'Color','k')
text(plot_win(1)+10,-14,'10s','HorizontalAlignment','center')
axis tight;
axis off

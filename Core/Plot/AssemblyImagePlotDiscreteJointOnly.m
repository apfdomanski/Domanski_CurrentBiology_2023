function AssemblyImagePlotDiscreteJointOnly(Factor,Rates,Units,Events,plot_win)
% Make raster of units higlighting those tagged as assembly members
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
%%
threshold = 3; % Which level of scrutiny to apply to thresholding
colorscale{1} = zeros(size(Factor.FSC{1},2),3);
colorscale{1}(:,3) = 1-0.5*rand(size(Factor.FSC{1},2),1);

colorscale{2} = zeros(size(Factor.FSC{2},2),3);
colorscale{2}(:,1) = 1-0.5*rand(size(Factor.FSC{2},2),1);
 
% colorscale{3} = zeros(size(Factor.FSC{3},2),3);
% colorscale{3}(:,2) = 1-0.5*rand(size(Factor.FSC{3},2),1);

colorscale{3} = (jet(size(Factor.FSC{3},2)));
% colorscale{3} = (colorcube);

% colorscale{3} = 0.7+0.3*rand(size(Factor.FSC{3},2),3);
% colorscale{3}(:,2) = 0.8+0.3*rand(size(Factor.FSC{3},2),1);
% colorscale{3}(colorscale{3}>0.9) = 0.9;
% 
% colorscale{3} = rand(size(Factor.FSC{3},2),3);
colorscale{3} = cbrewer('qual', 'Set1', size(Factor.FSC{3},2), 'cubic');%Set1

lag        = 1;    % window either side of peak activation to include spikes as part of assembly activation

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

figure('color','w'); hold on

%% (1) plot ALL units as a spike raster image

iArea = 3;
tb = plot_win(1):1:plot_win(2);
raster_ = zeros(length(tb),numel(Units{iArea}));
for unit_id=1:numel(Units{iArea})
    temp_firing=Units{iArea}{unit_id}.t/1e4;
    temp_firing=temp_firing(temp_firing>plot_win(1) & temp_firing<plot_win(2));
    raster_(:,unit_id) = histc(temp_firing,tb);
end


imagesc(tb,1:numel(Units{iArea}),raster_');
colormap(flipud(gray));caxis([0,20])
%% (2) Plot the assembly traces
for iArea = 3%1:3
    
%     for iAss=1:size(Factor.FSC{iArea},2)
%         %%%%
%         % plot the factor activations:
%         idx = Factor.TmtxS>=plot_win(1) & Factor.TmtxS<=plot_win(2);
%         plot(Factor.TmtxS(idx),Factor.FSC{iArea}(idx,iAss)-2*iArea - iAss*2,'color',colorscale{iArea}(iAss,:),'LineWidth',1.5)
%         %%%%
%     end
    
    
            idx = Factor.TmtxS>=plot_win(1) & Factor.TmtxS<=plot_win(2);

    for iAss=1:size(Factor.FSC{iArea},2)
        
        % (2) plot the factor activations:
        temp = Factor.FSC{iArea}(idx,iAss);
%         temp = mat2gray(temp);
        temp = temp-2*iArea - iAss*2;
%         temp = temp*nUnit_(iAss)-1;
        
        %     plot(Factor.TmtxS(idx),temp + offset,'color',colorscale{iArea}(iAss,:),'LineWidth',2)
        
        temptb = Factor.TmtxS(idx);
        tJump = [1;find(diff(temptb)>0.1);length(temptb)-1];
        
        for iBlock = 1:length(tJump)-1
            tbIdx = [tJump(iBlock)+1:tJump(iBlock+1)];
            
            plot(temptb(tbIdx),temp(tbIdx),'color',colorscale{iArea}(iAss,:),'LineWidth',1.5)
        end
        
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
nAss_ = cellfun(@(x) size(x,2), Factor.FSC);
for iArea = 3%1:3
    
    i = nAss_(iArea);
    for iAss=1:size(Factor.FSC{iArea},2)

%         factor_active_units = Factor.units{iArea}{iAss};
        factor_active_units = find(Factor.FL{iArea}{i}(:,iAss) >= Factor.ciHld(iArea,threshold));
 
    
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
                        this_unit*ones(length(temp_firing(temp_firing_ids)),1),40,colorscale{iArea}(iAss,:),'s','filled')
                end
            end
        end
    end
end
%% (5) Area markers
plot([plot_win(1) plot_win(1)]-5,[0 length(Factor.unitIDs{1})+0.5],'LineWidth',2,'Color',[0 0 1])
plot([plot_win(1) plot_win(1)]-5,[length(Factor.unitIDs{1})+0.5 length(Factor.unitIDs{3})+0.5],'LineWidth',2,'Color',[1 0 0])
text(plot_win(1)-10,length(Factor.unitIDs{1})/2,'mPFC Units','Rotation',90,'Color','b','HorizontalAlignment','center','VerticalAlignment','bottom')
text(plot_win(1)-10,length(Factor.unitIDs{1})+(length(Factor.unitIDs{3})-length(Factor.unitIDs{1}))/2,'dCA1 Units','Rotation',90,'Color','r','HorizontalAlignment','center','VerticalAlignment','bottom')
%% (6) overlay lever press times
hold on
% Collapse lever press events
TSample = [Events.Left(:,1);Events.Right(:,1)]/1e6;
TChoice = [Events.Left(:,2);Events.Right(:,2)]/1e6;
Outcome = [ones(size(Events.Left,1),1);2*ones(size(Events.Right,1),1)];
tickLength = 5;
Y_top = length(Units{3})+1;
Y_bottom = -25;
for iLever = 1:length(TSample)
    if [TSample(iLever) TChoice(iLever)] > plot_win(1) & [TSample(iLever) TChoice(iLever)] < plot_win(2)
        plot([TSample(iLever) TSample(iLever)],[Y_top Y_top+tickLength],'-g','LineWidth',2)
        plot([TChoice(iLever) TChoice(iLever)],[Y_top Y_top+tickLength],'-r','LineWidth',2)
        
        plot([TSample(iLever) TSample(iLever)],[Y_bottom Y_bottom-tickLength],'-g','LineWidth',2)
        plot([TChoice(iLever) TChoice(iLever)],[Y_bottom Y_bottom-tickLength],'-r','LineWidth',2)
        %plot([TSample(iLever) TChoice(iLever)],(Y_bottom-tickLength/2)*[1,1],':k','LineWidth',2)
        
        text(TSample(iLever),Y_bottom-tickLength,'Sample','Rotation',90,'Color','g','HorizontalAlignment','right')
        text(TChoice(iLever),Y_bottom-tickLength,'Choice','Rotation',90,'Color','r','HorizontalAlignment','right')
        switch Outcome(iLever)
            case 1
                text((TSample(iLever) + TChoice(iLever))/2,Y_bottom-tickLength,'L','Color','k','HorizontalAlignment','center')
            case 2
                text((TSample(iLever) + TChoice(iLever))/2,Y_bottom-tickLength,'R','Color','k','HorizontalAlignment','center')
        end
    end
end
%% (7) Set layout
xlabel('time (s)')
ylabel('Unit number')
xlim(plot_win)
plot([plot_win(1)+0 plot_win(1)+20],[-22 -22],'LineWidth',2,'Color','k')
text(plot_win(1)+10,-20,'10s','HorizontalAlignment','center')
axis tight;
ylim([-45 105])
xlim([5555 plot_win(2)])
 axis off
%% (A) highlight ASSEMBLY units and factors underneath
figure; hold on
iArea = 3;
jump = 20;
nAss_ = cellfun(@(x) size(x,2), Factor.FSC);
   i = nAss_(iArea);
% Work out unit spacing
for iAss=1:size(Factor.FSC{iArea},2)
    nUnit(iAss) = length(Factor.units{iArea}{iAss});
end
nUnit_=nUnit;
nUnit= nUnit+jump;
nUnit = cumsum(nUnit)-nUnit(1);
for iAss=1:size(Factor.FSC{iArea},2)
    % (1) Plot spike times during activations
    %factor_active_units = Factor.units{iArea}{iAss};
	factor_active_units = find(Factor.FL{iArea}{i}(:,iAss) >= Factor.ciHld(iArea,threshold));
    factor_active_units_{iAss} = factor_active_units
    for unit_id=1:numel(factor_active_units)
        % Get delimited spike times
        temp_firing=Units{iArea}{unit_id}.t/1e4;
        FiringTimes=[];
        for assem_peak_no=1:numel(FactorTimes{iArea}{iAss})
            temp_firing_ids=(temp_firing-lag)<=FactorTimes{iArea}{iAss}(assem_peak_no) &...
                (temp_firing+lag)>=FactorTimes{iArea}{iAss}(assem_peak_no) &...
                (temp_firing-lag)>=plot_win(1) &...
                (temp_firing-lag)<=plot_win(2);
            FiringTimes = [FiringTimes;temp_firing(temp_firing_ids)];
        end
            
            % highlight spike times of assembly units

%                 scatter(temp_firing(temp_firing_ids),...
%                     this_unit*ones(length(temp_firing(temp_firing_ids)),1)+nUnit(iAss),50,colorscale{iArea}(iAss,:),'s','filled')
                
                if factor_active_units(unit_id)<=length(Units{1})
                    scatter(FiringTimes,...
                        unit_id*ones(length(FiringTimes),1)+nUnit(iAss),10,[0 0 1],'s','filled','Markerfacealpha',1)
%                     plot([FiringTimes,FiringTimes]',(repmat(unit_id+[-0.5 0.5],length(FiringTimes),1)+nUnit(iAss))','b')
                else
                    scatter(FiringTimes,...
                        unit_id*ones(length(FiringTimes),1)+nUnit(iAss),10,[1 0 0],'s','filled','Markerfacealpha',1)
                    
%                     plot([FiringTimes,FiringTimes]',(repmat(unit_id+[-0.5 0.5],length(FiringTimes),1)+nUnit(iAss))','r');
                end        
    end
    
    % (2) plot the factor activations:
    idx = Factor.TmtxS>=plot_win(1) & Factor.TmtxS<=plot_win(2);
    temp = Factor.FSC{iArea}(idx,iAss);
    temp = mat2gray(temp);
    temp = temp*nUnit_(iAss);
    offset = nUnit(iAss)-nUnit_(iAss);
    %     plot(Factor.TmtxS(idx),temp + offset,'color',colorscale{iArea}(iAss,:),'LineWidth',2)
    
    temptb = Factor.TmtxS(idx);
    tJump = [1;find(diff(temptb)>0.1);length(temptb)-1];
    
    for iBlock = 1:length(tJump)-1
        tbIdx = [tJump(iBlock)+1:tJump(iBlock+1)];

        plot(temptb(tbIdx),temp(tbIdx) + offset,'color',colorscale{iArea}(iAss,:),'LineWidth',1.5)
    end
    
end
%  (5) overlay lever press times
% Collapse lever press events
% Collapse lever press events
TSample = [Events.Left(:,1);Events.Right(:,1)]/1e6;
TChoice = [Events.Left(:,2);Events.Right(:,2)]/1e6;
Outcome = [ones(size(Events.Left,1),1);2*ones(size(Events.Right,1),1)];
tickLength = 5;
Y_top = max(nUnit)+nUnit_(end)+5;
Y_bottom = -20;
for iLever = 1:length(TSample)
    if [TSample(iLever) TChoice(iLever)] > plot_win(1) & [TSample(iLever) TChoice(iLever)] < plot_win(2)
        plot([TSample(iLever) TSample(iLever)],[Y_top Y_top+tickLength],'-g','LineWidth',2)
        plot([TChoice(iLever) TChoice(iLever)],[Y_top Y_top+tickLength],'-r','LineWidth',2)
        
        plot([TSample(iLever) TSample(iLever)],[Y_bottom Y_bottom-tickLength],'-g','LineWidth',2)
        plot([TChoice(iLever) TChoice(iLever)],[Y_bottom Y_bottom-tickLength],'-r','LineWidth',2)
        %plot([TSample(iLever) TChoice(iLever)],(Y_bottom-tickLength/2)*[1,1],':k','LineWidth',2)
        
        text(TSample(iLever),Y_bottom-tickLength,'Sample','Rotation',90,'Color','g','HorizontalAlignment','right')
        text(TChoice(iLever),Y_bottom-tickLength,'Choice','Rotation',90,'Color','r','HorizontalAlignment','right')
        switch Outcome(iLever)
            case 1
                text((TSample(iLever) + TChoice(iLever))/2,Y_bottom-tickLength,'L','Color','k','HorizontalAlignment','center')
            case 2
                text((TSample(iLever) + TChoice(iLever))/2,Y_bottom-tickLength,'R','Color','k','HorizontalAlignment','center')
        end
    end
end
axis([min(plot_win) max(plot_win) -40 Inf]);
% axis([5555 max(plot_win) -40 Inf]);
axis off
%% (B) Plot Factor loading stem plot 
figure; hold on
nAss_ = cellfun(@(x) size(x,2), Factor.FSC);
step_ = 1.2;
for iArea =3% 1:3
    uList = 1:length(Units{iArea});
    i = nAss_(iArea);
    FL_ = Factor.FL{iArea}{i};
    thresh_H = Factor.ciHld(iArea,threshold); 
    thresh_L = Factor.ciLld(iArea,threshold); 
%     subplot(1,3,iArea); hold on
    if iArea==2
        Yoffset=length(Units{1});
    else
        Yoffset=0;
    end
    for iAss = 1:nAss_(iArea)
        sig_ = FL_(:,iAss)>thresh_H;
        FLsig = FL_(sig_,iAss);
        IDs = uList(sig_);
        bump_ = iAss*step_;
        plot([bump_ bump_],Yoffset+[0,length(Units{iArea})],'k')
        scatter(FL_(:,iAss)+bump_,uList+Yoffset,50,[0 0 0],'filled','Markerfacealpha',0.3)
        plot(bump_+[zeros(length(Units{iArea}),1),FL_(:,iAss)]',...
                  repmat((uList+Yoffset),2,1),'color',[0 0 0 0.6])
              
        scatter(FLsig+bump_,IDs+Yoffset,50,colorscale{iArea}(iAss,:),'filled')
        plot(bump_+[zeros(size(FLsig)) FLsig]',...
                  Yoffset+[IDs ;IDs],'color',colorscale{iArea}(iAss,:),'LineWidth',1.5)
%         plot(bump_+thresh_H*[1 1],[1 length(Units{3})],':k','LineWidth',1.5)
%         plot(bump_+thresh_L*[1 1],[1 length(Units{3})],':k','LineWidth',1.5)
%         rectangle('Position',[bump_+thresh_L 0  abs(thresh_L)+abs(thresh_H) length(Units{3})],'LineStyle','none','Facecolor',0.6*[1 1 1],'FaceAlpha',0.2)
        area([bump_+thresh_L bump_+thresh_H],[length(Units{3}) length(Units{3})],'LineStyle','none','Facecolor',0.6*[1 1 1],'FaceAlpha',0.3,'ShowBaseLine','off')
    end
    
plot([0.5 0.5],[0 length(Factor.unitIDs{1})+0.5],'LineWidth',2,'Color',[0 0 1])
plot([0.5 0.5],[length(Factor.unitIDs{1})+0.5 length(Factor.unitIDs{3})+0.5],'LineWidth',2,'Color',[1 0 0])
text(0.3,length(Factor.unitIDs{1})/2,'mPFC Units','Rotation',90,'Color','b','HorizontalAlignment','center')
text(0.3,length(Factor.unitIDs{1})+(length(Factor.unitIDs{3})-length(Factor.unitIDs{1}))/2,'dCA1 Units','Rotation',90,'Color','r','HorizontalAlignment','center')

%     axis([-1.5 nAss_(iArea)+step_ 0 length(Units{3})])
axis([-1.5 10 0 length(Units{3})])
ylim([-45 105])

axis off
end
%% (C) Plot average activation
% 
% TSample = [Events.Left(:,1);Events.Right(:,1)]/1e6;
% TChoice = [Events.Left(:,2);Events.Right(:,2)]/1e6;
tLims = [-5 5];
bw = 0.05;%diff(Factor.TmtxS([1 2]));
LTrials = [];
for iTrial = 1:size(Events.Left,1)
    tSample = FindClosestIndex(Factor.TmtxS,Events.Left(iTrial,1)/1e6)+(tLims/bw);
    tChoice = FindClosestIndex(Factor.TmtxS,Events.Left(iTrial,2)/1e6)+(tLims/bw);
    idx = [tSample(:,1):tSample(:,2),tChoice(:,1):tChoice(:,2)];
    LTrials = cat(3,LTrials,Factor.FSC{3}(idx,:));
end
RTrials = [];
for iTrial = 1:size(Events.Right,1)
    tSample = FindClosestIndex(Factor.TmtxS,Events.Right(iTrial,1)/1e6)+(tLims/bw);
    tChoice = FindClosestIndex(Factor.TmtxS,Events.Right(iTrial,2)/1e6)+(tLims/bw);
    idx = [tSample(:,1):tSample(:,2),tChoice(:,1):tChoice(:,2)];
    RTrials = cat(3,RTrials,Factor.FSC{3}(idx,:));
end
tb=(1:size(RTrials,1))*bw;

figure; hold on
for iAss = 1:size(Factor.FSC{3},2)
    offset = iAss*6;
    lTemp = squeeze(LTrials(:,iAss,:)); rTemp = squeeze(RTrials(:,iAss,:));
    Temp = ([mean(lTemp,2);mean(rTemp,2)]);
    m_ = mean(Temp); s_ = abs(std(Temp));
    lTemp = (lTemp-m_)/s_; 
    rTemp = (rTemp-m_)/s_;
    tempL = nanmean(lTemp,2); tempL_ = nansem(lTemp,2);
    tempR = nanmean(rTemp,2); tempR_ = nansem(rTemp,2);
    plot(tb,offset+tempL,'color',colorscale{3}(iAss,:),'LineStyle',':','LineWidth',1.5)
    plot(tb,offset+tempR,'color',colorscale{3}(iAss,:),'LineStyle','-','LineWidth',1.5)
    ciplot(offset+tempL+tempL_,...
           offset+tempL-tempL_,...
           tb,colorscale{3}(iAss,:),0.3)
       
    ciplot(offset+tempR+tempR_,...
           offset+tempR-tempR_,...
           tb,colorscale{3}(iAss,:),1)       
end
axis off
%% (D) Plot zoomed in spikes of one assembly

figure; hold on
iArea = 3;
jump = 0;
% Work out unit spacing
for iAss=1:size(Factor.FSC{iArea},2)
    nUnit(iAss) = length(Factor.units{iArea}{iAss});
end
nUnit_=nUnit;
nUnit= nUnit+jump;
nAss_ = cellfun(@(x) size(x,2), Factor.FSC);
   i = nAss_(iArea);
 nUnit = cumsum(nUnit)-nUnit(1);
for iAss=1%:size(Factor.FSC{iArea},2)
    % (1) Plot spike times during activations
    factor_active_units = find(Factor.FL{iArea}{i}(:,iAss) >= Factor.ciHld(iArea,threshold));
    %factor_active_units = Factor.units{iArea}{iAss};
    for unit_id=1:numel(factor_active_units)
        % Get all spikes - plot in grey
        temp_firing = Units{iArea}{unit_id}.t/1e4;
        temp_firing_ids = (temp_firing-lag)>=plot_win(1) &...
                          (temp_firing-lag)<=plot_win(2);
        FiringTimes = temp_firing(temp_firing_ids);
%         scatter(FiringTimes,unit_id*ones(length(FiringTimes),1)+nUnit(iAss),10,0.6*[1 1 1],'s','filled','Markerfacealpha',1)
        plot([FiringTimes,FiringTimes]',(repmat(unit_id+[-0.5 0.5],length(FiringTimes),1)+nUnit(iAss))','color',[0 0 0 0.2])
        % Get delimited spike times - show colours
        FiringTimes=[];
        for assem_peak_no=1:numel(FactorTimes{iArea}{iAss})
            temp_firing_ids=(temp_firing-lag)<=FactorTimes{iArea}{iAss}(assem_peak_no) &...
                (temp_firing+lag)>=FactorTimes{iArea}{iAss}(assem_peak_no) &...
                (temp_firing-lag)>=plot_win(1) &...
                (temp_firing-lag)<=plot_win(2);
            FiringTimes = [FiringTimes;temp_firing(temp_firing_ids)];
        end
            
            % highlight spike times of assembly units

                if factor_active_units(unit_id)<=length(Units{1})
%                     scatter(FiringTimes,unit_id*ones(length(FiringTimes),1)+nUnit(iAss),10,[0 0 1],'s','filled','Markerfacealpha',1)
                    plot([FiringTimes,FiringTimes]',(repmat(unit_id+[-0.5 0.5],length(FiringTimes),1)+nUnit(iAss))','color',[0 0 1])
                else
%                     scatter(FiringTimes,unit_id*ones(length(FiringTimes),1)+nUnit(iAss),10,[1 0 0],'s','filled','Markerfacealpha',1)
                    plot([FiringTimes,FiringTimes]',(repmat(unit_id+[-0.5 0.5],length(FiringTimes),1)+nUnit(iAss))','color',[1 0 0])
                end        
    end
    
    % (2) plot the factor activations:
    idx = Factor.TmtxS>=plot_win(1) & Factor.TmtxS<=plot_win(2);
    temp = Factor.FSC{iArea}(idx,iAss);
    temp = mat2gray(temp);
    temp = temp*nUnit_(iAss);
    offset = nUnit(iAss)-nUnit_(iAss);
    %     plot(Factor.TmtxS(idx),temp + offset,'color',colorscale{iArea}(iAss,:),'LineWidth',2)
    
    temptb = Factor.TmtxS(idx);
    tJump = [1;find(diff(temptb)>0.1);length(temptb)-1];
    
    for iBlock = 1:length(tJump)-1
        tbIdx = [tJump(iBlock)+1:tJump(iBlock+1)];

        plot(temptb(tbIdx),temp(tbIdx) + offset,'color',colorscale{iArea}(iAss,:),'LineWidth',2)
    end
    
end
%  (5) overlay lever press times
% Collapse lever press events
% Collapse lever press events
TSample = [Events.Left(:,1);Events.Right(:,1)]/1e6;
TChoice = [Events.Left(:,2);Events.Right(:,2)]/1e6;
Outcome = [ones(size(Events.Left,1),1);2*ones(size(Events.Right,1),1)];
tickLength = 3;
Y_top = max(nUnit)+nUnit_(end)+5;
Y_bottom = -15;
for iLever = 1:length(TSample)
    if [TSample(iLever) TChoice(iLever)] > plot_win(1) & [TSample(iLever) TChoice(iLever)] < plot_win(2)
        
        plot([TSample(iLever) TSample(iLever)],[Y_bottom Y_bottom-tickLength],'-g','LineWidth',2)
        plot([TChoice(iLever) TChoice(iLever)],[Y_bottom Y_bottom-tickLength],'-r','LineWidth',2)
        %plot([TSample(iLever) TChoice(iLever)],(Y_bottom-tickLength/2)*[1,1],':k','LineWidth',2)
        
        text(TSample(iLever),Y_bottom-tickLength,'Sample','Rotation',90,'Color','g','HorizontalAlignment','right')
        text(TChoice(iLever),Y_bottom-tickLength,'Choice','Rotation',90,'Color','r','HorizontalAlignment','right')
        switch Outcome(iLever)
            case 1
                text((TSample(iLever) + TChoice(iLever))/2,Y_bottom-tickLength,'L','Color','k','HorizontalAlignment','center')
            case 2
                text((TSample(iLever) + TChoice(iLever))/2,Y_bottom-tickLength,'R','Color','k','HorizontalAlignment','center')
        end
    end
end
axis([min(plot_win) max(plot_win) -25 Inf]);
axis off;%; grid on

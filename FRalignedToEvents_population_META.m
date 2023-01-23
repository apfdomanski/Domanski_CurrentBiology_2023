clear

Areas = {'HP','PFC'};
Delays_ = {'Naive','Debutant','Expert'};

Events_ = {'CueLight','SamplePress','DelayEnd','NosePoke','ChoicePress','RewardConsume'};
Events__ = {'Cue Light','Sample Press','Delay End','Nose Poke','Choice Press','Reward Collection'};


pat = 'C:\Analysis\AssemblyAnalysis\raw\UnitTuning';
x = load([pat filesep 'UnitTuning_Short.mat'],'ModInd');
ModInd{1} = x.ModInd;
x = load([pat filesep 'UnitTuning_Medium.mat'],'ModInd');
ModInd{2} =x.ModInd;
x = load([pat filesep 'UnitTuning_Long.mat'],'ModInd');
ModInd{3} = x.ModInd;
clear x
%% Plot fraction of neurons tuned to: [cue, sample press]
x = [0.2 ; 0.4; 0.6];
x = repmat(x,1,3);
for iEvent=1:2
    figure('name',Events__{iEvent}); hold on
    for iArea=1:2
        iDelay=2;
        % ModInd{iStage}{iArea}{iEvent}{iDelay};
        
        y = [NaN, NaN, NaN];
        % Early training
        y_ = ModInd{1}{iArea}{iEvent}{1};
        y_ = y_./sum(y_,2);
        y = [y; y_];
        
        % Middle training
        y_= [ModInd{2}{iArea}{iEvent}{1};...
            ModInd{2}{iArea}{iEvent}{2};...
            ModInd{2}{iArea}{iEvent}{3};...
            ModInd{2}{iArea}{iEvent}{4}];
        y_ = mean(y_./sum(y_,2));
        y = [y; y_];
        
        % Expert performance
        y_= [ModInd{3}{iArea}{iEvent}{1};...
            ModInd{3}{iArea}{iEvent}{2};...
            ModInd{3}{iArea}{iEvent}{3}];
        y_ = mean(y_./sum(y_,2));
        y = [y; y_];
        
        
        
        y(1,:)=[];
        
        subplot(length(Areas),1,iArea); hold on
        % bar(x,y,'stacked')
        bar(1:3,-y(:,1),'FaceColor','r')
        bar(1:3,y(:,3),'FaceColor','b')
        title([Areas{iArea} ' Units (Correct trials)'])
        ylabel('Fraction of responsive units')
        
        ticks_ = [1:3];
        set(gca,'XTick',ticks_,'XTickLabel',Delays_,'XTickLabelRotation',45)
        axis([0 4 -1 1])
    end
end
%% Plot fraction of neurons tuned to: [end of delay tone]
x = [0.2 ; 0.4 ; 0.6 ; 0.8 ; 1.2 ; 1.4 ; 1.6];
iEvent=3;
figure('name',Events__{iEvent}); hold on
for iArea=1:2
    % ModInd{iStage}{iArea}{iEvent}{iDelay};
    
    y = [NaN, NaN, NaN];
    
    % Middle training
    y_= [ModInd{2}{iArea}{iEvent}{1};...
        ModInd{2}{iArea}{iEvent}{2};...
        ModInd{2}{iArea}{iEvent}{3};...
        ModInd{2}{iArea}{iEvent}{4}];
    y = [y; y_];
    
    % Expert performance
    y_= [ModInd{3}{iArea}{iEvent}{1};...
        ModInd{3}{iArea}{iEvent}{2};...
        ModInd{3}{iArea}{iEvent}{3}];
    y = [y; y_];
    
    
    
    y(1,:)=[];
    y = y./sum(y,2);
    
    
    subplot(length(Areas),1,iArea); hold on
    % bar(x,y,'stacked')
    bar(x,-y(:,1),'FaceColor','r')
    bar(x,y(:,3),'FaceColor','b')
    title([Areas{iArea} ' Units (Correct trials)'])
    ylabel('Fraction of responsive units')
    if iArea ==1
        set(gca,'XTick',[])
    else
        set(gca,'XTick',x,'XTickLabelRotation',45,...
            'XTickLabel',...
            {'Debutant (2s)',...
            'Debutant (4s)',...
            'Debutant (6s)',...
            'Debutant (8s)',...
            'Expert (4s)',...
            'Expert(8s)',...
            'Expert (16s)'})
    end
    axis([0 2 -1 1])
end

%% Plot fraction of neurons tuned to: [Nose Poke, Choice press, Reward consumation]
x = [-0.2; 0.2 ; 0.4 ; 0.6 ; 0.8 ; 1.2 ; 1.4 ; 1.6];
for iEvent=4:6
    figure('name',Events__{iEvent}); hold on
    
    for iArea=1:2
        
        
        % ModInd{iStage}{iArea}{iEvent}{iDelay};
        
        y = [NaN, NaN, NaN];
        % Early training
        y_ = ModInd{1}{iArea}{iEvent}{1};
        y = [y; y_];
        
        % Middle training
        y_= [ModInd{2}{iArea}{iEvent}{1};...
            ModInd{2}{iArea}{iEvent}{2};...
            ModInd{2}{iArea}{iEvent}{3};...
            ModInd{2}{iArea}{iEvent}{4}];
        y = [y; y_];
        
        % Expert performance
        y_= [ModInd{3}{iArea}{iEvent}{1};...
            ModInd{3}{iArea}{iEvent}{2};...
            ModInd{3}{iArea}{iEvent}{3}];
        y = [y; y_];
        y(1,:)=[];
        y = y./sum(y,2);
        
        subplot(length(Areas),1,iArea); hold on
        % bar(x,y,'stacked')
        bar(x,-y(:,1),'FaceColor','r')
        bar(x,y(:,3),'FaceColor','b')
        title([Areas{iArea} ' Units (Correct trials)'])
        ylabel('Fraction of responsive units')
        if iArea ==1
            set(gca,'XTick',[])
        else
            set(gca,'XTick',x,'XTickLabelRotation',45,...
                'XTickLabel',...
                {'Beginner (0s)',...
                'Debutant (2s)',...
                'Debutant (4s)',...
                'Debutant (6s)',...
                'Debutant (8s)',...
                'Expert (4s)',...
                'Expert(8s)',...
                'Expert (16s)'})
        end
        axis([-0.4 2 -1 1])
    end
end
%% Plot fraction of neurons tuned to all events averaged across delays - Expert stage only
x = 1:length(Events__);
x = repmat(x,1,3);
figure;
for iArea=1:2
    subplot(2,1,iArea); hold on
   
        
        
        % ModInd{iStage}{iArea}{iEvent}{iDelay};
        
        y = [NaN, NaN, NaN];
        for iEvent=1:length(Events__)
        % Expert performance
        y_= [ModInd{3}{iArea}{iEvent}{1};...
            ModInd{3}{iArea}{iEvent}{2};...
            ModInd{3}{iArea}{iEvent}{3}];
        y_ = mean(y_./sum(y_,2));
        
        y = [y; y_];
        end
        
        y(1,:)=[];
        
%         b = bar(y,'stacked')
%         b(1).FaceColor = 'b';
%         b(2).FaceColor = 'w';
%         b(3).FaceColor = 'r';
        bar(1:length(Events__),-y(:,1),'FaceColor','r')
        bar(1:length(Events__),y(:,3),'FaceColor','b')
        title([Areas{iArea} ' Units'])
        ylabel('Fraction of units')
        
        if iArea==1
                set(gca,'XTick',[])
        else
            ticks_ = 1:length(Events__);
            set(gca,'XTick',ticks_,'XTickLabel',Events__,'XTickLabelRotation',45)
        end
    end




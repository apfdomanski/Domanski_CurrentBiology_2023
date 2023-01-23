clear
Delays_ = {'Short','Medium','Long'};

tlimsPoke = [-5 5];

bw=0.05;
warning ('off')
pat = 'C:\Analysis\AssemblyAnalysis\raw';
cd(pat)
%%%%%
fileList=dir('allTimestamps\*.mat');
iFile =6;
Areas = {'HP','PFC'};
iArea =2
%%%%%
fname=strtok(fileList(iFile).name,'_');
load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname),'t');
% load(sprintf('%s\\%s.mat',pat,fname),'PFCcells');

% load(sprintf('%s\\KDE_bins\\%s_%s_iFR50.mat',pat,fname,Areas{iArea}),'iFR','Tmtx');
load(sprintf('%s\\KDE_binsTaskonly\\%s_%s_iFR50_behavOnly.mat',pat,fname,Areas{iArea}),'iFR','Tmtx');
length_ = sum(abs(tlimsPoke))/bw;
tb_  = (1:length_)*bw+tlimsPoke(1);
nUnits = size(iFR,2);
%% loop across units - Correct trials 
data_{1}.iFR = iFR;
for iDelay =1:length(Delays_)
    
    eval(sprintf('t_L = t.%s.NosePoke_LeftCorrect*1e-6;',Delays_{iDelay}))
    eval(sprintf('t_R = t.%s.NosePoke_RightCorrect*1e-6;',Delays_{iDelay}))
    [~,cutoutL{iDelay}]=SelTimesFRs(tlimsPoke,Tmtx,data_,t_L);
    [~,cutoutR{iDelay}]=SelTimesFRs(tlimsPoke,Tmtx,data_,t_R);
end

for iUnit=1:nUnits
    m=[];e=[];
    for iDelay =1:length(Delays_)
        LP = nan(length_,1);  RP = nan(length_,1);
        
        for iTrial = 1:length(cutoutL{iDelay}{1})
            try
                LP= [LP,cutoutL{iDelay}{1}{iTrial}(:,iUnit)];
            end
        end
        for iTrial = 1:length(cutoutR{iDelay}{1})
            try
                RP= [RP,cutoutR{iDelay}{1}{iTrial}(:,iUnit)];
            end
        end
        if size(LP,2)>1
            LP(:,1)=[];
        end
        if size(RP,2)>1
            RP(:,1)=[];
        end   
            
        
        
        m = [m,nanmean(LP,2),nanmean(RP,2')];
        e = [e,nansem(LP,2),nansem(RP,2)];
    end
    m = (m-mean(m))-std(m);e = (e-mean(m))-std(m);
    Lm=[];Rm=[];Le=[];Re=[];
    for iDelay =1:length(Delays_)
        Lm(:,iDelay) = m(1:length_);
        Rm(:,iDelay) = m(length_+1:2*length_);
        Le(:,iDelay) = e(1:length_);
        Re(:,iDelay) = e(length_+1:2*length_);
        m(1:2*length_)=[];e(1:2*length_)=[];
    end
    
    Lpokes.m{iUnit} = Lm;
    Rpokes.m{iUnit} = Rm;
    Lpokes.e{iUnit} = Le;
    Rpokes.e{iUnit} = Re;
end
clear Lm Rm Le Re e m LP RP cutoutL cutoutR 
%% loop across units - Error trials 
data_{1}.iFR = iFR;
for iDelay =1:length(Delays_)
    
    eval(sprintf('t_L = t.%s.NosePoke_LeftError*1e-6;',Delays_{iDelay}))
    eval(sprintf('t_R = t.%s.NosePoke_RightError*1e-6;',Delays_{iDelay}))
    [~,cutoutL{iDelay}]=SelTimesFRs(tlimsPoke,Tmtx,data_,t_L);
    [~,cutoutR{iDelay}]=SelTimesFRs(tlimsPoke,Tmtx,data_,t_R);
end

for iUnit=1:nUnits
    m=[];e=[];
    for iDelay =1:length(Delays_)
        LP = nan(length_,1);  RP = nan(length_,1);
        
        for iTrial = 1:length(cutoutL{iDelay}{1})
            try
                LP= [LP,cutoutL{iDelay}{1}{iTrial}(:,iUnit)];
            end
        end
        for iTrial = 1:length(cutoutR{iDelay}{1})
            try
                RP= [RP,cutoutR{iDelay}{1}{iTrial}(:,iUnit)];
            end
        end
        if size(LP,2)>1
            LP(:,1)=[];
        end
        if size(RP,2)>1
            RP(:,1)=[];
        end   
            
        
        
        m = [m,nanmean(LP,2),nanmean(RP,2')];
        e = [e,nansem(LP,2),nansem(RP,2)];
        
    end
    m = (m-mean(m))-std(m);e = (e-mean(m))-std(m);
    Lm=[];Rm=[];Le=[];Re=[];
    for iDelay =1:length(Delays_)
        Lm(:,iDelay) = m(1:length_);
        Rm(:,iDelay) = m(length_+1:2*length_);
        Le(:,iDelay) = e(1:length_);
        Re(:,iDelay) = e(length_+1:2*length_);
        m(1:2*length_)=[];e(1:2*length_)=[];
    end
    
    LpokesError.m{iUnit} = Lm;
    RpokesError.m{iUnit} = Rm;
    LpokesError.e{iUnit} = Le;
    RpokesError.e{iUnit} = Re;
end
clear Lm Rm Le Re e m LP RP cutoutL cutoutR 

%% Plot correct trials: all delays overlaid, separate by L/R condition
bins =min(tb_):0.5:max(tb_);
Clatencyhist=[];
Platencyhist=[];

SF = 4; offset = -10;
for iDelay =1:length(Delays_)
    eval(sprintf('Clatencyhist.L(iDelay,:)  = histc((t.%s.ChoicePress_LeftCorrect-t.%s.NosePoke_LeftCorrect)*1e-6,bins);',Delays_{iDelay},Delays_{iDelay}))
    eval(sprintf('Clatencyhist.R(iDelay,:)  = histc((t.%s.ChoicePress_RightCorrect-t.%s.NosePoke_RightCorrect)*1e-6,bins);',Delays_{iDelay},Delays_{iDelay}))

    eval(sprintf('Platencyhist.L(iDelay,:)  = histc((t.%s.DelayEnd_LeftCorrect-t.%s.NosePoke_LeftCorrect)*1e-6,bins);',Delays_{iDelay},Delays_{iDelay}))
    eval(sprintf('Platencyhist.R(iDelay,:)  = histc((t.%s.DelayEnd_RightCorrect-t.%s.NosePoke_RightCorrect)*1e-6,bins);',Delays_{iDelay},Delays_{iDelay}))

end



DelCol = {[0.1 0.1 0.7],[0.2 0.7 0.2],[0.7 0.2 0.2]};
AxLims = [min(tb_) max(tb_) -10 10];
for iUnit = 1%:nUnits
    figure('color','w','Position',[200 500 1200 400],'Units','Pixels');
    for iDelay =1:length(Delays_)
        
        subplot(1,2,1);hold on
        ciplot(Lpokes.m{iUnit}(:,iDelay)-Lpokes.e{iUnit}(:,iDelay),...
               Lpokes.m{iUnit}(:,iDelay)+Lpokes.e{iUnit}(:,iDelay),tb_,DelCol{iDelay},0.8);
        axis(AxLims)
        title('Left trials','Color',[0.5 0 0.5])
        
        subplot(1,2,2);hold on
        ciplot(Rpokes.m{iUnit}(:,iDelay)-Rpokes.e{iUnit}(:,iDelay),...
               Rpokes.m{iUnit}(:,iDelay)+Rpokes.e{iUnit}(:,iDelay),tb_,DelCol{iDelay},0.8);
        
        axis(AxLims)
        title('Right trials','Color',[0 0 1])
        
    end
    
    for iDelay =1:length(Delays_)
        
        subplot(1,2,1);hold on
        y = Clatencyhist.L(iDelay,:); y = y./sum(y);
        stairs(bins,SF*y+offset,'color',DelCol{iDelay},'LineWidth',1.5)       
        y = Platencyhist.L(iDelay,:); y = y./sum(y);
        stairs(bins,SF*y+offset,'color',DelCol{iDelay},'LineWidth',1.5,'LineStyle',':')
        
        subplot(1,2,2);hold on
        y = Clatencyhist.R(iDelay,:); y = y./sum(y);
        stairs(bins,SF*y+offset,'color',DelCol{iDelay},'LineWidth',1.5)
        y = Platencyhist.R(iDelay,:); y = y./sum(y);
        stairs(bins,SF*y+offset,'color',DelCol{iDelay},'LineWidth',1.5,'LineStyle',':')
        
    end
    subplot(1,2,1); plot([0 0],[-3 8],'k','LineWidth',2)
    text(0,9,'Nosepoke','HorizontalAlignment','center')
    text(max(tb_),offset,'Choice presses','HorizontalAlignment','right','VerticalAlignment','bottom')
    text(max(tb_),offset+1,'(Delay ends)','HorizontalAlignment','right','VerticalAlignment','bottom')


    subplot(1,2,2); plot([0 0],[-3 8],'k','LineWidth',2)
    
    legend('4s Delay','8s Delay','16s Delay'); legend boxoff
    
end

%% Plot correct trials: all delays overlaid, separate by L/R condition


DelCol = {[0.1 0.1 0.7],[0.2 0.7 0.2],[0.7 0.2 0.2]};
AxLims = [min(tb_) max(tb_) -10 10];
for iUnit = 1:10%nUnits
    figure('color','w','Position',[200 500 1200 400],'Units','Pixels');
        
        subplot(2,3,1);hold on
        ciplot(LpokesError.m{iUnit}(:,1)-LpokesError.e{iUnit}(:,1),...
               LpokesError.m{iUnit}(:,1)+LpokesError.e{iUnit}(:,1),tb_,[1 0 0],0.8);
           
           
        ciplot(Lpokes.m{iUnit}(:,1)-Lpokes.e{iUnit}(:,1),...
               Lpokes.m{iUnit}(:,1)+Lpokes.e{iUnit}(:,1),tb_,[0 0 0],0.8);
           
        axis(AxLims)
        title('Left trials, 4s','Color',[0.5 0 0.5])
        
         subplot(2,3,2);hold on
        ciplot(LpokesError.m{iUnit}(:,2)-LpokesError.e{iUnit}(:,2),...
               LpokesError.m{iUnit}(:,2)+LpokesError.e{iUnit}(:,2),tb_,[1 0 0],0.8);
           
           
        ciplot(Lpokes.m{iUnit}(:,2)-Lpokes.e{iUnit}(:,2),...
               Lpokes.m{iUnit}(:,2)+Lpokes.e{iUnit}(:,2),tb_,[0 0 0],0.8);
           
        axis(AxLims)
        title('Left trials, 8s','Color',[0.5 0 0.5])
        
        subplot(2,3,3);hold on
        ciplot(LpokesError.m{iUnit}(:,3)-LpokesError.e{iUnit}(:,3),...
               LpokesError.m{iUnit}(:,3)+LpokesError.e{iUnit}(:,3),tb_,[1 0 0],0.8);
          
        ciplot(Lpokes.m{iUnit}(:,3)-Lpokes.e{iUnit}(:,3),...
               Lpokes.m{iUnit}(:,3)+Lpokes.e{iUnit}(:,3),tb_,[0 0 0],0.8);
           
        axis(AxLims)
        title('Left trials, 16s','Color',[0.5 0 0.5])
       
        
 subplot(2,3,4);hold on
        ciplot(RpokesError.m{iUnit}(:,1)-RpokesError.e{iUnit}(:,1),...
               RpokesError.m{iUnit}(:,1)+RpokesError.e{iUnit}(:,1),tb_,[1 0 0],0.8);
           
           
        ciplot(Rpokes.m{iUnit}(:,1)-Rpokes.e{iUnit}(:,1),...
               Rpokes.m{iUnit}(:,1)+Rpokes.e{iUnit}(:,1),tb_,[0 0 0],0.8);
           
        axis(AxLims)
        title('Right trials, 4s','Color',[0 0 1])
        
         subplot(2,3,5);hold on
        ciplot(RpokesError.m{iUnit}(:,2)-RpokesError.e{iUnit}(:,2),...
               RpokesError.m{iUnit}(:,2)+RpokesError.e{iUnit}(:,2),tb_,[1 0 0],0.8);
           
           
        ciplot(Rpokes.m{iUnit}(:,2)-Rpokes.e{iUnit}(:,2),...
               Rpokes.m{iUnit}(:,2)+Rpokes.e{iUnit}(:,2),tb_,[0 0 0],0.8);
           
        axis(AxLims)
        title('Right trials, 8s','Color',[0 0 1])
        
        subplot(2,3,6);hold on
        ciplot(RpokesError.m{iUnit}(:,3)-RpokesError.e{iUnit}(:,3),...
               RpokesError.m{iUnit}(:,3)+RpokesError.e{iUnit}(:,3),tb_,[1 0 0],0.8);
          
        ciplot(Lpokes.m{iUnit}(:,3)-Rpokes.e{iUnit}(:,3),...
               Rpokes.m{iUnit}(:,3)+Rpokes.e{iUnit}(:,3),tb_,[0 0 0],0.8);
           
        axis(AxLims)
        title('Right trials, 16s','Color',[0 0 1])    
end




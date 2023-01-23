%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Target = 'LONG';

shift = 0;
plotOnline = false;
bw=0.05;
Nbs = 500;

warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
else
    pat = '/Volumes/HDD2/DNMTP/raw/';
end
cd(pat)
fileList=dir(sprintf('MixedSelectivity_LONG%s*%s*_MixedSelectivity_MembershipsortedUnits.mat',filesep,Target));
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};
MemberClasses ={'LocalMembers','JointMembers','NonMembers'};
MemberClasses_ ={'Local Assembly Members','Inter-area Assembly Members','Non-members'};
col_ ={[0.9 0.6 0],[0.6 0 0.6],[0.6 0.6 0.6]};
normWin = [0 4];
normaliseFscores = true;
%% import
clear D_ D;
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    fnIn = sprintf('%sMixedSelectivity_LONG%s%s_MixedSelectivity_MembershipsortedUnits.mat',pat,filesep,fname);
    load(fnIn ,'D');
    
    for iDelay = 1:length(Delays_)
        for iClass= 1:length(MemberClasses)
            
            for iArea =1:2
                try
                    eval(sprintf('D_.%s.LR.%s{iArea}.TS{iFile}= D.%s.LR.%s{iArea}.TS;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    eval(sprintf('D_.%s.LR.%s{iArea}.Ft2{iFile}= D.%s.LR.%s{iArea}.Ft2;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    eval(sprintf('D_.%s.LR.%s{iArea}.Rt2{iFile}= D.%s.LR.%s{iArea}.Rt2;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    eval(sprintf('D_.%s.LR.%s{iArea}.Ft2ciH{iFile}= D.%s.LR.%s{iArea}.Ft2ciH;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    % eval(sprintf('D_.%s.LR.%s{iArea}.Ft2_drawnTrials{iFile}= D.%s.LR.%s{iArea}.Ft2_drawnTrials;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    
                    eval(sprintf('D_.%s.LR.%s{iArea}.CVE{iFile}= D.%s.LR.%s{iArea}.CVE;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    % eval(sprintf('D_.%s.LR.%s{iArea}.CVEbs{iFile}= D.%s.LR.%s{iArea}.CVEbs;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    disp(fname)
                end
                try
                    
                    eval(sprintf('D_.%s.LR_err.%s{iArea}.TS{iFile}= D.%s.LR_err.%s{iArea}.TS;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    eval(sprintf('D_.%s.LR_err.%s{iArea}.Ft2{iFile}= D.%s.LR_err.%s{iArea}.Ft2;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    eval(sprintf('D_.%s.LR_err.%s{iArea}.Rt2{iFile}= D.%s.LR_err.%s{iArea}.Rt2;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    eval(sprintf('D_.%s.LR_err.%s{iArea}.CVE{iFile}= D.%s.LR_err.%s{iArea}.CVE;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    % eval(sprintf('D_.%s.LR_err.%s{iArea}.CVEbs{iFile}= D.%s.LR_err.%s{iArea}.CVEbs;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    eval(sprintf('D_.%s.LR_err.%s{iArea}.Ft2ciH{iFile}= D.%s.LR_err.%s{iArea}.Ft2ciH;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                    
                    disp(fname)
                end
            end
        end
    end
end
% clear D
%% Plot Fscores

clear y y_ F Fsig
% x = (1:size(y_,2))*bw;
x = (1:402)*bw;
Del_idx = [2:3];
Ypos ={[-0*ones(length(x),1);-0.2*ones(length(x),1)];...
       [-0.5*ones(length(x),1);-0.7*ones(length(x),1)];...
       [-1*ones(length(x),1);-1.2*ones(length(x),1)]};
YposBox = {[0 -0.3 20 0.4];...
           [0 -0.8 20 0.4];...
           [0 -1.3 20 0.4]};
% correct       
for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        clear y y_ yCI_ yCI ysigTF
        
        for iDelay =Del_idx%1:length(Delays_)
            try
                eval(sprintf('y_{iDelay,1}   = D_.%s.LR.%s{iArea}.Ft2;',Delays_{iDelay},MemberClasses{iClass}))
                eval(sprintf('yCI_{iDelay,1} = D_.%s.LR.%s{iArea}.Ft2ciH;',Delays_{iDelay},MemberClasses{iClass}))
            end
        end
        % collapse across delay lengths
        for iFile =1:max(cellfun(@length,y_))
            
            n = 0;
            y{iFile}      = zeros(size(x));
            yCI{iFile}    = zeros(size(x));
            ysigTF{iFile} = cell(size(Del_idx));%nan(size(x));
            
            % Accumulate by delay
            for iDelay = Del_idx%1:length(Delays_)
                ysigTF{iFile}{iDelay} = nan(size(x));
                try
                    temp = y_{iDelay}{iFile};
%                      if normaliseFscores
%                      bp = find( x>normWin(1) & x<normWin(2) );
% % 
%                         temp  = temp - mean(temp(bp));
%                         temp  = temp ./mean(temp(bp));
%                      end
                    y{iFile}      =  y{iFile} + temp;
                    yCI{iFile}    =  yCI{iFile} + yCI_{iDelay}{iFile};
                    ysigTF{iFile}{iDelay} = double(y_{iDelay}{iFile} > yCI_{iDelay}{iFile});
                    n=n+1;
                catch
                    disp('!')
                end
            end
            
            % average across delays
            if  nansum(y{iFile}) == 0
                y{iFile} = nan(size(x));
%                 ysigTF{iFile} = nan(size(x));
            else
                y{iFile} =  y{iFile}./n;
            end
        end
        y = cell2mat(y');
        y(isnan(nansum(y,2)),:)=[];
        
        yCI = cell2mat(yCI');
        yCI(isnan(nansum(y,2)),:)=[];
        Fsig{iArea}{iClass} = cell2mat(cellfun(@(c) cell2mat(c'), ysigTF,'UniformOutput',false)');
        Fsig{iArea}{iClass} = nansum(Fsig{iArea}{iClass})./size(Fsig{iArea}{iClass},1);
        
        if normaliseFscores
            bp = find( x>normWin(1) & x<normWin(2) );
            y = bsxfun(@minus,y,mean(y(:,bp),2));
            y = bsxfun(@rdivide,y,mean(y(:,bp),2));
%             B  = nanmean(y(:,bp)');
%             y = y./(B'*ones(1,length(x)));
        end
            F{iArea}{iClass} = y;
            clear y_ y B bp yCI ysigTF yCI_
    end
end

% error
%for iArea = 1:2
%     for iClass = 1:length(MemberClasses)
%         clear y y_ yCI_ yCI ysigTF
%         
%         for iDelay =Del_idx%1:length(Delays_)
%             try
%                 eval(sprintf('y_{iDelay,1}   = D_.%s.LR_err.%s{iArea}.Ft2;',Delays_{iDelay},MemberClasses{iClass}))
%                 eval(sprintf('yCI_{iDelay,1} = D_.%s.LR_err.%s{iArea}.Ft2ciH;',Delays_{iDelay},MemberClasses{iClass}))
%             end
%         end
%         % collapse across delay lengths
%         for iFile =1:length(fileList)%max(cellfun(@length,y_))
%             
%             n = 0;
%             y{iFile}      = zeros(size(x));
%             yCI{iFile}    = zeros(size(x));
%             ysigTF{iFile} = cell(size(Del_idx));%nan(size(x));
%             
%             % Accumulate by delay
%             for iDelay = Del_idx%1:length(Delays_)
%                 ysigTF{iFile}{iDelay} = nan(size(x));
%                 try
%                     y{iFile}      =  y{iFile} + y_{iDelay}{iFile};
%                     yCI{iFile}    =  yCI{iFile} + yCI_{iDelay}{iFile};
%                     ysigTF{iFile}{iDelay} = double(y_{iDelay}{iFile} > yCI_{iDelay}{iFile});
%                     n=n+1;
%                 end
%             end
%             
%             % average across delays
%             if  nansum(y{iFile}) == 0
%                 y{iFile} = nan(size(x));
% %                 ysigTF{iFile} = nan(size(x));
%             else
%                 y{iFile} =  y{iFile}./n;
%             end
%         end
%         y = cell2mat(y');
%         y(isnan(nansum(y,2)),:)=[];
%         
%         yCI = cell2mat(yCI');
%         yCI(isnan(nansum(y,2)),:)=[];
%         Fsig_err{iArea}{iClass} = cell2mat(cellfun(@(c) cell2mat(c'), ysigTF,'UniformOutput',false)');
%         Fsig_err{iArea}{iClass} = nansum(Fsig{iArea}{iClass})./size(Fsig{iArea}{iClass},1);
%         
%         if normaliseFscores
%             bp = find( x>normWin(1) & x<normWin(2) );
%             B  = nanmean(y(:,bp)');
%             y = y - B'*ones(1,length(x));
%             y = y./(B'*ones(1,length(x)));
%         end
%         if ~isempty(y)
%             F_err{iArea}{iClass} = y;
%         else
%             F_err{iArea}{iClass} = nan(2,length(x));
%         end
%         clear y_ y B bp yCI ysigTF yCI_
%     end
% end




figure
for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        subplot(1,2,iArea); hold on
        title ([Areas{iArea} ' Units'])
        
        ciplot(nanmean( F{iArea}{iClass},1) + nansem( F{iArea}{iClass},1),...
               nanmean( F{iArea}{iClass},1) - nansem( F{iArea}{iClass},1),...
               x, col_{iClass},1);
        
%         ciplot(nanmean( F_err{iArea}{iClass},1) + nansem( F_err{iArea}{iClass},1),...
%                nanmean( F_err{iArea}{iClass},1) - nansem( F_err{iArea}{iClass},1),...
%                x, col_{iClass},0.5);
        
               
%     temp = Fsig{iArea}{iClass};
%     imagesc(repmat(x,1,2),Ypos{iClass}, [temp;temp]);
%     rectangle('Position',YposBox{iClass},'LineWidth',1,'EdgeColor',col_{iClass})


%     plot([10 10],[0.2 0.6],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[0.2 0.6],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[0.2 0.6],'color',[1 0 0 0.6],'LineWidth',1.5)
%     
%     plot([10 10],[-1.8 -1.4],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[-1.8 -1.4],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[-1.8 -1.4],'color',[1 0 0 0.6],'LineWidth',1.5)
%     
    if normaliseFscores
        axis([0 20 -Inf Inf])
    else
        axis([0 20 0 10])
    end
%     axis([0 20 -3 1])
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'},'Location','northoutside','Orientation','horizontal')
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
%      if iArea==1
        ylabel('L/R decoding (F-Score)')
        plot([7 12], [5 5],'k','LineWidth',1.5)
        plot([7 7], [5 8],'k','LineWidth',1.5)
%     end
    colormap(flipud(gray(10)));
    caxis([0 0.5])
    axis off
        
        
        %         plot(x,y,'color',col_{iClass})
    end
    xlabel('Time (s)')
end
% legend(MemberClasses_,'Orientation','horizontal','Location','north'); legend boxoff
%% Plot Fscores all collapsed
clear y y_
% x = (1:size(y_,2))*bw;
x = (1:402)*bw;

% correct
for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        clear y y_
        for iDelay = 1:length(Delays_)
            try
                eval(sprintf('y_{iDelay,1} = cell2mat(D_.%s.LR.%s{iArea}.Ft2'');',Delays_{iDelay},MemberClasses{iClass}))
                eval(sprintf('yCI_{iDelay,1} =  cell2mat(D_.%s.LR.%s{iArea}.Ft2ciH'');',Delays_{iDelay},MemberClasses{iClass}))
            end
        end
        % collapse across delay lengths
        
        y   = cell2mat(y_);
        yCI = cell2mat(yCI_);
                
        yCI(isnan(nansum(y,2)),:)=[];
        y(isnan(nansum(y,2)),:)=[];
        ysigTF = double(y>yCI);
        if normaliseFscores
            bp = find(x>normWin(1) & x<normWin(2));
            B  = nanmean(y(:,bp)');
            y = y./(B'*ones(1,length(x)));
        end
        FCollapse{iArea}{iClass} = y;
        FciHCollapse{iArea}{iClass} = yCI;
        Fsig{iArea}{iClass} = ysigTF;
        Fsig{iArea}{iClass} = nansum(Fsig{iArea}{iClass})./size(Fsig{iArea}{iClass},1);

        clear y_
        
        
    end
end
% error
% for iArea = 1:2
%     for iClass = 1:length(MemberClasses)
%         clear y y_
%         for iDelay = 1:length(Delays_)
%             try
%                 eval(sprintf('y_{iDelay,1} = cell2mat(D_.%s.LR_err.%s{iArea}.Ft2'');',Delays_{iDelay},MemberClasses{iClass}))
%                 %                 y_ = cell2mat(y_');
%             catch
%                 %                 y_ = zeros(12,402);
%             end
%         end
%         % collapse across delay lengths
%         
%         y = cell2mat(y_);
%         y(isnan(nansum(y,2)),:)=[];
%         
%         if normaliseFscores
%             bp = find(x>normWin(1) & x<normWin(2));
%             B  = nanmean(y(:,bp)');
%             %                 y = y - B'*ones(1,length(x));
%             y = y./(B'*ones(1,length(x)));
%         end
%         FCollapse_err{iArea}{iClass} = y;
% 
%         %             y = nanmean(y_,1);
%         clear y_
%         
%         
%     end
% end

figure;
        
for iArea = 1:2
    subplot(1,2,iArea); hold on
    
    for iClass = 1:length(MemberClasses)    
        
        ciplot(nanmean(FCollapse{iArea}{iClass})+nansem(FCollapse{iArea}{iClass}),nanmean(FCollapse{iArea}{iClass})-nansem(FCollapse{iArea}{iClass}),...
            x, col_{iClass},1);
        
%         ciplot(nanmean(FciHCollapse{iArea}{iClass})+nansem(FciHCollapse{iArea}{iClass}),nanmean(FciHCollapse{iArea}{iClass})-nansem(FciHCollapse{iArea}{iClass}),...
%             x, col_{iClass},0.3);
%         ciplot(nanmean(FCollapse_err{iArea}{iClass})+nansem(FCollapse_err{iArea}{iClass}),nanmean(FCollapse_err{iArea}{iClass})-nansem(FCollapse_err{iArea}{iClass}),...
%             x, col_{iClass},0.5);        
    xlim([0 20])
%     plot([10 10],[0.2 0.6],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[0.2 0.6],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[0.2 0.6],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    
%     temp = Fsig{iArea}{iClass};
%     imagesc(repmat(x,1,2),Ypos{iClass}, [temp;temp]);
%     rectangle('Position',YposBox{iClass},'LineWidth',1,'EdgeColor',col_{iClass})
% 
% 
%     plot([10 10],[0.2 0.6],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[0.2 0.6],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[0.2 0.6],'color',[1 0 0 0.6],'LineWidth',1.5)
%     
%     plot([10 10],[-1.8 -1.4],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[-1.8 -1.4],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[-1.8 -1.4],'color',[1 0 0 0.6],'LineWidth',1.5)
    
%     axis([0 20 -2 40])
    axis([0 20 -0.5 20])
%         axis([0 20 -3 1])

%     axis([0 20 0 1])
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'},'Location','northoutside','Orientation','horizontal')
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
%      if iArea==1
        ylabel('L/R decoding (F-Score)')
        plot([7 12], [5 5],'k','LineWidth',1.5)
        plot([7 7], [5 8],'k','LineWidth',1.5)        
         colormap(flipud(gray(10)));
    caxis([0 1])
    axis off
    end
end
%% Plot Fscores sorted by delay length
clear y y_
% x = (1:size(y_,2))*bw;
x = (1:402)*bw;

% correct
for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        clear y y_ yCI_ yCI
        for iDelay = 1:length(Delays_)
            try
                eval(sprintf('y_{iDelay,1}   =  cell2mat(D_.%s.LR.%s{iArea}.Ft2'');',   Delays_{iDelay},MemberClasses{iClass}))
                eval(sprintf('yCI_{iDelay,1} =  cell2mat(D_.%s.LR.%s{iArea}.Ft2ciH'');',Delays_{iDelay},MemberClasses{iClass}))
            end
            yCI_{iDelay,1}(isnan(nansum(y_{iDelay,1},2)),:)=[];
            y_{iDelay,1}(isnan(nansum(y_{iDelay,1},2)),:)=[];
            ysigTF_{iDelay,1} = double(y_{iDelay,1}>yCI_{iDelay,1});
            
        
        
       
        if normaliseFscores
            bp = find(x>normWin(1) & x<normWin(2));
            B  = nanmean(y_{iDelay,1}(:,bp)');
            y_{iDelay,1} = y_{iDelay,1}./(B'*ones(1,length(x)));
        end
        FCollapse_{iArea}{iClass}{iDelay,1} = y_{iDelay,1};
        FciHCollapse_{iArea}{iClass}{iDelay,1} = yCI_{iDelay,1};
        Fsig_{iArea}{iClass}{iDelay,1} = nansum(ysigTF_{iDelay,1})./size(ysigTF_{iDelay,1},1);

%         clear y_
        end
        
    end
end
% error
% for iArea = 1:2
%     for iClass = 1:length(MemberClasses)
%         clear y y_
%         for iDelay = 1:length(Delays_)
%             try
%                 eval(sprintf('y_{iDelay,1} = cell2mat(D_.%s.LR_err.%s{iArea}.Ft2'');',Delays_{iDelay},MemberClasses{iClass}))
%                 %                 y_ = cell2mat(y_');
%             catch
%                 %                 y_ = zeros(12,402);
%             end
%         end
%         % collapse across delay lengths
%         
%         y = cell2mat(y_);
%         y(isnan(nansum(y,2)),:)=[];
%         
%         if normaliseFscores
%             bp = find(x>normWin(1) & x<normWin(2));
%             B  = nanmean(y(:,bp)');
%             %                 y = y - B'*ones(1,length(x));
%             y = y./(B'*ones(1,length(x)));
%         end
%         FCollapse_err{iArea}{iClass} = y;
% 
%         %             y = nanmean(y_,1);
%         clear y_
%         
%         
%     end
% end

figure;
        
for iArea = 1:2
    figure('Name',Areas{iArea})
    
    for iClass = 1:length(MemberClasses)    
        subplot(1,3,iClass ); hold on
         for iDelay = 1:length(Delays_)
             ciplot(nanmean(FCollapse_{iArea}{iClass}{iDelay})+nansem(FCollapse_{iArea}{iClass}{iDelay}),nanmean(FCollapse_{iArea}{iClass}{iDelay})-nansem(FCollapse_{iArea}{iClass}{iDelay}),...
                 x, col_{iClass},1-(iDelay-1)*0.3);
         end
%         ciplot(nanmean(FciHCollapse{iArea}{iClass})+nansem(FciHCollapse{iArea}{iClass}),nanmean(FciHCollapse{iArea}{iClass})-nansem(FciHCollapse{iArea}{iClass}),...
%             x, col_{iClass},0.3);
%         ciplot(nanmean(FCollapse_err{iArea}{iClass})+nansem(FCollapse_err{iArea}{iClass}),nanmean(FCollapse_err{iArea}{iClass})-nansem(FCollapse_err{iArea}{iClass}),...
%             x, col_{iClass},0.5);        
    xlim([0 20])
%     plot([10 10],[0.2 0.6],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[0.2 0.6],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[0.2 0.6],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    
%     temp = Fsig{iArea}{iClass};
%     imagesc(repmat(x,1,2),Ypos{iClass}, [temp;temp]);
%     rectangle('Position',YposBox{iClass},'LineWidth',1,'EdgeColor',col_{iClass})
% 
% 
%     plot([10 10],[0.2 0.6],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[0.2 0.6],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[0.2 0.6],'color',[1 0 0 0.6],'LineWidth',1.5)
%     
%     plot([10 10],[-1.8 -1.4],'color',[0 0 0 0.6],'LineWidth',1.5)
%     plot([5 5],[-1.8 -1.4],'color',[0 1 0 0.6],'LineWidth',1.5)
%     plot([15 15],[-1.8 -1.4],'color',[1 0 0 0.6],'LineWidth',1.5)
    
%     axis([0 20 -2 40])
    axis([0 20 -0.5 20])
%         axis([0 20 -3 1])

%     axis([0 20 0 1])
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'},'Location','northoutside','Orientation','horizontal')
%     legend({'Short delay (4s)','Medium delay (8s)','Long delay (16s)'})
%      if iArea==1
        ylabel('L/R decoding (F-Score)')
        plot([7 12], [5 5],'k','LineWidth',1.5)
        plot([7 7], [5 8],'k','LineWidth',1.5)        
         colormap(flipud(gray(10)));
    caxis([0 1])
    axis off
    end
end

%% Plot correct - CVE
Del_idx = [1:3];
class_idx = [2];
figure;
clear y y_
% x = (1:size(y_,2))*bw;
x = (1:402)*bw;

% correct
for iArea = 1:2
    for iClass = class_idx%1:length(MemberClasses)
        clear y y_
        for iDelay = Del_idx%1:length(Delays_)
            try
                eval(sprintf('y_{iDelay,1} = D_.%s.LR.%s{iArea}.CVE;',Delays_{iDelay},MemberClasses{iClass}))
                %                 y_ = cell2mat(y_');
            catch
                %                 y_ = zeros(12,402);
            end
        end
        % collapse across delay lengths
        for iFile =1:max(cellfun(@length,y_))
            n = 0;
            y{iFile} = zeros(size(x));
            for iDelay = Del_idx%1:length(Delays_)
                try
                    y{iFile} =  y{iFile} + 1-y_{iDelay}{iFile};
                    n=n+1;
                end
            end
            if  nansum(y{iFile}) == 0
                y{iFile} = nan(size(x));
            else
                y{iFile} =  y{iFile}./n;
            end
        end
        y = cell2mat(y');
        y = smooth2a(y,0,5);

        y(isnan(nansum(y,2)),:)=[];
        CVE{iArea}{iClass} = y;
        clear y_
    end
end
% error
for iArea = 1:2
    for iClass = class_idx%1:length(MemberClasses)
        clear y y_
        for iDelay = 1:length(Delays_)
            try
                eval(sprintf('y_{iDelay,1} = D_.%s.LR_err.%s{iArea}.CVE;',Delays_{iDelay},MemberClasses{iClass}))
                %                 y_ = cell2mat(y_');
            catch
                %                 y_ = zeros(12,402);
            end
        end
        % collapse across delay lengths
        for iFile =1:length(fileList)%max(cellfun(@length,y_))
            n = 0;
            y{iFile} = zeros(size(x));
            for iDelay = 1:length(Delays_)
                try
                    y{iFile} =  y{iFile} + 1-y_{iDelay}{iFile};
                    n=n+1;
                end
            end
            % average across delays
            if  nansum(y{iFile}) == 0
                y{iFile} = nan(size(x));
            else
                y{iFile} =  y{iFile}./n;
            end
        end
        y = cell2mat(y');
        y = smooth2a(y,0,5);
        y(isnan(nansum(y,2)),:)=[];
        CVE_err{iArea}{iClass} = y;
        clear y_
    end
end

for iArea = 1:2
    subplot(1,2,iArea); hold on
    for iClass = class_idx%1:length(MemberClasses)        
        ciplot(nanmean(CVE{iArea}{iClass})+nansem(CVE{iArea}{iClass}),nanmean(CVE{iArea}{iClass})-nansem(CVE{iArea}{iClass}),...
            x, col_{iClass},1);
%         ciplot(nanmean(CVE_err{iArea}{iClass})+nansem(CVE_err{iArea}{iClass}),nanmean(CVE_err{iArea}{iClass})-nansem(CVE_err{iArea}{iClass}),...
%             x, col_{iClass},0.5);        
    end
        axis([0 20 0 1])

end
%% Plot correct - CVE all collapsed
Del_idx = [1:3];
class_idx = [1:3];
figure;
clear y y_
% x = (1:size(y_,2))*bw;
x = (1:402)*bw;
% correct
for iArea = 1:2
    for iClass = class_idx%1:length(MemberClasses)
        clear y y_
        for iDelay = Del_idx%1:length(Delays_)
            try
                eval(sprintf('y_{iDelay,1} = 1-cell2mat(D_.%s.LR.%s{iArea}.CVE'');',Delays_{iDelay},MemberClasses{iClass}))
                %                 y_ = cell2mat(y_');
            catch
                %                 y_ = zeros(12,402);
            end
        end
        
        % collapse across delay lengths
        y = cell2mat(y_);
        y = smooth2a(y,0,5);
        y(isnan(nansum(y,2)),:)=[];
        CVEcollapse{iArea}{iClass} = y;
        clear y_
    end
end
% error
for iArea = 1:2
    for iClass = class_idx%1:length(MemberClasses)
        clear y y_
        for iDelay = Del_idx%1:length(Delays_)
            try
                eval(sprintf('y_{iDelay,1} = 1-cell2mat(D_.%s.LR_err.%s{iArea}.CVE'');',Delays_{iDelay},MemberClasses{iClass}))
                %                 y_ = cell2mat(y_');
            catch
                %                 y_ = zeros(12,402);
            end
        end
        
        % collapse across delay lengths
        try
            y = cell2mat(y_);
        y = smooth2a(y,0,5);
        y(isnan(nansum(y,2)),:)=[];
        catch
           y = nan(12,402);
        end
        CVEcollapse_err{iArea}{iClass} = y;
        clear y_
    end
end

for iArea = 1:2
    subplot(1,2,iArea); hold on
    
    for iClass = class_idx%1:length(MemberClasses)
        ciplot(nanmean(CVEcollapse{iArea}{iClass})+nansem(CVEcollapse{iArea}{iClass}),nanmean(CVEcollapse{iArea}{iClass})-nansem(CVEcollapse{iArea}{iClass}),...
               x, col_{iClass},1);
%         ciplot(nanmean(CVEcollapse_err{iArea}{iClass})+nansem(CVEcollapse_err{iArea}{iClass}),nanmean(CVEcollapse_err{iArea}{iClass})-nansem(CVEcollapse_err{iArea}{iClass}),...
%                x, col_{iClass},0.5);           
        xlim([0 20])
         plot([10 10],[0.2 0.6],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot([5 5],[0.2 0.6],'color',[0 1 0 0.6],'LineWidth',1.5)
    plot([15 15],[0.2 0.6],'color',[1 0 0 0.6],'LineWidth',1.5)
    
    
    axis([0 20 0 1])
        
    end
end

%% Plot correct - mean TS correct
trange = 1:402;%100:300;
clear y_ meanTS meanTS_
figure
for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        y_ = nan(3,length(fileList));
        for iDelay = 1:length(Delays_)
            for iFile = 1:length(fileList)
                try
                    eval(sprintf('y_(iDelay,iFile) = nanmean(nanmean(D_.%s.LR.%s{iArea}.TS{iFile}(trange,:)));',Delays_{iDelay},MemberClasses{iClass}))
                end
            end
        end
        meanTS{iArea,iClass} = nanmean(y_);
    end
    
%     meanTS_{iArea}(1,:) = nanmax([meanTS{iArea,1};meanTS{iArea,2}]);
%     meanTS_{iArea}(2,:) = meanTS{iArea,3};    
    meanTS_{iArea}(1,:) = nanmean([meanTS{iArea,1};meanTS{iArea,3}]);
    meanTS_{iArea}(2,:) = meanTS{iArea,2};
    subplot(1,2,iArea)
    x=meanTS_{iArea};
    rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')
    ttest(x(1,:),x(2,:))
    addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\nanmean.m')
    plot(x)
end
x = nanmean(cat(3,meanTS_{1},meanTS_{2}),3);
rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')
ttest(x(1,:),x(2,:))
addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\nanmean.m')
figure
plot(x)

%% Plot correct - mean TS correct
trange = 1:402;%100:300;
clear y_ meanTS meanTS_
figure
for iArea = 1:2
    for iFile = 1:length(fileList)
        for iDelay = 1:3
            x=[];
            try
                eval(sprintf('x = [x , D_.%s.LR.%s{iArea}.TS{iFile}(trange,:)];',Delays_{iDelay},MemberClasses{1}))
            end
            try
                eval(sprintf('x = [x , D_.%s.LR.%s{iArea}.TS{iFile}(trange,:)];',Delays_{iDelay},MemberClasses{3}))
            end
            y_nonmembers(iDelay,iFile) = nanmean(nanmean(x,1));
            
            try
                eval(sprintf('x = D_.%s.LR.%s{iArea}.TS{iFile}(trange,:);',Delays_{iDelay},MemberClasses{2}))
            catch
                x=NaN;
            end
            
            y_members(iDelay,iFile) = nanmean(nanmean(x,1));
        end
    end
    y_{iArea}(1,:) = nanmax(y_members);
    y_{iArea}(2,:) = nanmax(y_nonmembers);
    
     subplot(1,2,iArea)
    x= y_{iArea};
    rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')
    ttest(x(1,:),x(2,:))
    addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\nanmean.m')
    plot(x)
end

    for iClass = 1:length(MemberClasses)
        y_ = nan(3,length(fileList));
        for iDelay = 1:length(Delays_)
            for iFile = 1:length(fileList)
                try
                    eval(sprintf('y_(iDelay,iFile) = nanmean(nanmean(D_.%s.LR.%s{iArea}.TS{iFile}(trange,:)));',Delays_{iDelay},MemberClasses{iClass}))
                end
            end
        end
        meanTS{iArea,iClass} = nanmean(y_);
    end
    
%     meanTS_{iArea}(1,:) = nanmax([meanTS{iArea,1};meanTS{iArea,2}]);
%     meanTS_{iArea}(2,:) = meanTS{iArea,3};    
    meanTS_{iArea}(1,:) = nanmean([meanTS{iArea,1};meanTS{iArea,3}]);
    meanTS_{iArea}(2,:) = meanTS{iArea,2};
    subplot(1,2,iArea)
    x=meanTS_{iArea};
    rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')
    ttest(x(1,:),x(2,:))
    addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\nanmean.m')
    plot(x)
end
x = nanmean(cat(3,meanTS_{1},meanTS_{2}),3);
rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\')
ttest(x(1,:),x(2,:))
addpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite\nanmean.m')
figure
plot(x)



%% Continuous factor model
clear all; close all
set(0,'DefaultFigureWindowStyle','docked')
flags.plotIndividualChanges = true;
if ispc
%     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\';       % location of the processed sleep assembly activations
%     pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - epoch specific FRs\';
    pat{1}  = 'C:\Analysis\AssemblyAnalysis\Sleep\Decoding vs flanking periods - use task mean FRs\';
    pat{2} = 'C:\Analysis\AssemblyAnalysis\raw';           % location of raw spike times
else
    pat{1} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Decoding vs sleep/';
    pat{2} = '/Users/aleksanderdomanski/Dropbox/Assembly analysis/Task';
end
pat{3} = [pat{2} filesep 'KDE_bins'];                  % location of calculated firing P.rates for Task period
cd(pat{1})
Rats = dir([pat{1} filesep '*_DecodingVSsleep.mat']);

Ignore= {'JaroslawLONG1'; ...
         'MiroslawLONG2'; ...
         'NorbertLONG2' ; ...
         'OnufryLONG2'};
     
% Ignore= {};

for iList = 1:length(Ignore)
    Rats = Rats(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
           {Rats.name},repmat({Ignore(iList)},1,length(Rats)),'UniformOutput',false)));
end
clear Ignore iList 
%  Process assembly activation averages
Group = struct;
Group.Assem.FSC = {};
Group.Assem.FSC_Mean = cell(3,1);
Group.Assem.FSC_SEM = cell(3,1);
for iRat= 1:length(Rats)
    load([pat{1} Rats(iRat).name],'P','D','FAtrials');
    disp(['Working on ' Rats(iRat).name '...'])
    for s=1:length(P.names) 
        if ~isempty(FAtrials.FSC_{s}(1))
    temp = cell2mat(FAtrials.FSC_{s}'); temp_m=[]; temp_s=[];
    for i = 1:size(temp,2)
        temp_m (:,i) = nanmean(reshape(temp(:,i),[600,length(temp(:,i))/600])');
        temp_s (:,i) = nansem(reshape(temp(:,i),[600,length(temp(:,i))/600])');
    end
    Group.Assem.FSC_Mean{s} = [Group.Assem.FSC_Mean{s}; temp_m'];
    Group.Assem.FSC_SEM{s}  = [Group.Assem.FSC_SEM{s}; temp_m'];
    
        end
    end
end
clear iRat s i temp_m temp_s temp

% Plot each assembly

figure('color','w'); hold on
tb = (1:P.Ltr*2)*P.bw;
clr2={'b','g','r'};

for s= 1:3
    temp_m  = Group.Assem.FSC_Mean{s};
    temp_s  = Group.Assem.FSC_SEM{s};
    clr_    = rand(size(temp_m,1),3);
    for i=1:size(temp_m,1)
        subplot(1,3,s); hold on
        title(P.names{s})
        
        %         ciplot(5*i+temp_m(i,:)-temp_s(i,:),...
        %             5*i+temp_m(i,:)+temp_s(i,:),...
        %             tb,clr_(i,:));
        plot(tb,5*i+temp_m(i,:),'color',clr_(i,:),'LineWidth',1.5)
        %     plot(tb,staggerplot(temp_m',0, 5))
        %     axis([min(tb) max(tb), -10 10])
        
        
    end
    plot([15 15],[-1.5 max(get(gca,'YLim'))],':k')
    plot([10 10],[-1.5 max(get(gca,'YLim'))],'color','g','LineWidth',2)
    plot([20 20],[-1.5 max(get(gca,'YLim'))],'r','LineWidth',2)
    set(gca,'children',flipud(get(gca,'children')))
    axis([0 30 -10 140])
    
    if s==2
        plot([17 17],[60 65],'k','LineWidth',1.5)
        text(18,62,'5 A.U.')
    end
    axis off
end

% Plot mean assembly activation

figure('color','w'); hold on
tb = (1:P.Ltr*2)*P.bw;
clr2={'b','g','r'};

for s= 1:3
   temp_m  = Group.Assem.FSC_Mean{s};
   temp_s  = Group.Assem.FSC_SEM{s};


   ciplot(mean(temp_m)+nansem(temp_m),...
          mean(temp_m)-nansem(temp_m),...
          tb,clr2{s},0.6); 
%    plot(tb,mean(temp_m),'color',clr2{s},'LineWidth',1.5)
    axis([min(tb) max(tb), -2 5])
   
end
legend(P.names); legend('boxoff')
% xlabel('Time (s)')
set(gca, 'XTick',[5:5:25],'XTickLabel',{'Pre-sample','Sample','Delay','Choice','Post-choice'});

ylabel({'Assembly activation';'(Mean +/- SEM Factor score)'})

plot([15 15],[-1.5 2],':k')
plot([10 10],[-1.5 2],'g','LineWidth',3)
plot([20 20],[-1.5 2],'r','LineWidth',3)
clear clr2 clr_ temp_m temp_s tb i s

%% Discrete factor models
if ispc
%   
    pat{1} = 'C:\Analysis\AssemblyAnalysis\raw\KDE_bins\LONG\';           % location of raw spike times
else
    pat{1} =  '/Volumes/HDD2/DNMTP/raw/KDE_binsTaskonly/LONGTaskonly/';

end
cd(pat{1})
Rats = dir([pat{1} filesep '*FSC.mat']);
% Rats = dir([pat{1} filesep '*_iFR50__FSCtemp.mat']);

%  Process assembly activation averages
Group = struct;
Group.Assem.FSC = {};
Group.Assem.FSC_Mean = cell(3,1);
Group.Assem.FSC_SEM = cell(3,1);
for iRat= 1:length(Rats)
    load([pat{1} Rats(iRat).name],'FSCsel');
%     load([pat{1} Rats(iRat).name],'FSC');
%     FSCsel = FSC;
    disp(['Working on ' Rats(iRat).name '...'])
    
    for s=1:3
        if ~isempty(FSCsel{s})
            temp = FSCsel{s}; temp_m=[]; temp_s=[];
            for i = 1:size(temp,2)
                
                temp_m (:,i) = nanmean(reshape(temp(:,i),[600,length(temp(:,i))/600])');
                temp_s (:,i) = nansem(reshape(temp(:,i),[600,length(temp(:,i))/600])');
            end
            Group.Assem.FSC_Mean{s} = [Group.Assem.FSC_Mean{s}; temp_m'];
            Group.Assem.FSC_SEM{s}  = [Group.Assem.FSC_SEM{s}; temp_m'];
            
        end
    end
end
clear iRat s i temp_m temp_s temp

% Plot each assembly
figure('color','w'); hold on
tb = (1:600)*0.05-0.05/2;
clr2={'b','g','r'};
for s= 1:3
    
    temp_m  = Group.Assem.FSC_Mean{s};
    temp_s  = Group.Assem.FSC_SEM{s};
    clr_    = rand(size(temp_m,1),3);
    for i=1:size(temp_m,1)
        subplot(1,3,s); hold on
        title(P.names{s})
%         ciplot(5*i+temp_m(i,:)-temp_s(i,:),...
%             5*i+temp_m(i,:)+temp_s(i,:),...
%             tb,clr_(i,:));
        plot(tb,5*i+temp_m(i,:),'color',clr_(i,:),'LineWidth',1.5)
        %     plot(tb,staggerplot(temp_m',0, 5))
        %     axis([min(tb) max(tb), -10 10])
    end
plot([15 15],[-1.5 max(get(gca,'YLim'))],':k')
plot([10 10],[-1.5 max(get(gca,'YLim'))],'color','g','LineWidth',2)
plot([20 20],[-1.5 max(get(gca,'YLim'))],'r','LineWidth',2)
set(gca,'children',flipud(get(gca,'children')))
axis([0 30 -10 140])
    if s==2
        plot([17 17],[60 65],'k','LineWidth',1.5)
        text(18,62,'5 A.U.')
    end
axis off
end

% Plot mean assembly activation
figure('color','w'); hold on
tb = (1:600)*0.05;
clr2={'b','g','r'};

for s= 1:3
   temp_m  = Group.Assem.FSC_Mean{s};
   temp_s  = Group.Assem.FSC_SEM{s};

   ciplot(mean(temp_m)+nansem(temp_m),...
          mean(temp_m)-nansem(temp_m),...
          tb,clr2{s},0.6); 
%    plot(tb,mean(temp_m),'color',clr2{s},'LineWidth',1.5)
    axis([min(tb) max(tb), -2 5])
   
end
legend(P.names); legend('boxoff')
% xlabel('Time (s)')
% ylabel('Mean +/- SEM assembly activation strength')
set(gca, 'XTick',[5:5:25],'XTickLabel',{'Pre-sample','Sample','Delay','Choice','Post-choice'});
ylabel({'Assembly activation';'(Mean +/- SEM Factor score)'})

plot([15 15],[-1.5 2],':k')
plot([10 10],[-1.5 2],'g','LineWidth',3)
plot([20 20],[-1.5 2],'r','LineWidth',3)

clear clr2 clr_ temp_m temp_s tb i ss
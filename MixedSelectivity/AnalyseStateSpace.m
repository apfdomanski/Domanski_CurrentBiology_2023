clear 
Delays_ = {'Short','Medium','Long'};
Delays__ = [4,8,16];
Target = 'LONG';

offset = [-5 5];
plotOnline = false;
bw=0.05;
clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
else
    pat = '/Volumes/HDD2/DNMTP/raw';
end
cd(pat)
fileList=dir(sprintf('allTimestamps%s*%s*.mat',filesep,Target));
fileListAss = fileList;

reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat' }; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileListAss(find(name_flag))=[];% fnames_(name_flag)=[];
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag



Areas = {'HP','PFC'};
%% Batch process units
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s%sallTimestamps%s%s_Events.mat',pat,filesep,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname));
    for iArea = 1:length(Areas)
        fprintf('Analysing run %d/%d %s (%s)...',iFile,length(fileList),fname,Areas{iArea})
        load(sprintf('%s%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep,filesep,fname,Areas{iArea}));
        
        iFR_=iFR;
        iFR_=zscore(iFR_);
        
        
%% PCA on whole trials concatenated the lever presses into one vector
% 
% Ltrials = [];                   Rtrials = [];
% nL = 0;                         nR = 0;
% LtrialsTimes = [];              RtrialsTimes = [];
% LtrialsLengths = [];            RtrialsLengths = [];
% 
% LeftTrials =[t.Long.CueLight_LeftCorrect,t.Long.ChoicePress_LeftCorrect];
% LeftTrialsEvents = [t.Long.CueLight_LeftCorrect,    ...
%     t.Long.SamplePress_LeftCorrect, ...
%     t.Long.DelayEnd_LeftCorrect,    ...
%     t.Long.ChoicePress_LeftCorrect];
% 
% 
% RightTrials =[t.Long.CueLight_RightCorrect,t.Long.ChoicePress_RightCorrect];
% RightTrialsEvents = [t.Long.CueLight_RightCorrect,    ...
%                      t.Long.SamplePress_RightCorrect, ...
%                      t.Long.DelayEnd_RightCorrect,    ...
%                      t.Long.ChoicePress_RightCorrect];
%                  
% LjPCA=[];
% for iTrial =1:size(LeftTrials,1)
%     try
%         tlimsLong_ = LeftTrials(iTrial,:)/1e6 + offset;
%         tlimsLong_ = closest(Tmtx,tlimsLong_);
%         FR = iFR_(tlimsLong_(1):tlimsLong_(2),:);
%         Ltrials = [Ltrials;FR];
%         LjPCA(iTrial).A=FR;
%         LjPCA(iTrial).times = (0:size(FR,1)-1)*mean(diff(Tmtx))'+offset(1);
%         LjPCA(iTrial).times=LjPCA(iTrial).times*1000;
% %         LjPCA(iTrial).A=zscore(FR(1:650,:));
% %         LjPCA(iTrial).times =  ((0:399)*mean(diff(Tmtx))+offset(1))';
% %         LjPCA(iTrial).times =  1000*((0:649)*mean(diff(Tmtx)))';
%         LtrialsLengths(iTrial) = diff(tlimsLong_)+1;
%         tlimsLong_= LeftTrialsEvents(iTrial,:)/1e6 + [offset(1), 0 0 0];
%         tlimsLong_= closest(Tmtx,tlimsLong_ );
%         LtrialsTimes(iTrial,:) = tlimsLong_-tlimsLong_(1);
%         nL=nL+1;
%         
%         
%     end
% end
% 
% RjPCA=[];
% for iTrial =1:size(RightTrials,1)
%     try
%         tlimsLong_ = RightTrials(iTrial,:)/1e6 + offset;
%         tlimsLong_ = closest(Tmtx,tlimsLong_);
%         FR = iFR_(tlimsLong_(1):tlimsLong_(2),:);
%         Rtrials = [Rtrials;FR];
%         RjPCA(iTrial).A=FR;
%         RjPCA(iTrial).times = (0:size(FR,1)-1)*mean(diff(Tmtx))'+offset(1);
%         RjPCA(iTrial).times=RjPCA(iTrial).times*1000;
% %         RjPCA(iTrial).A=FR(1:200,:);
% %         RjPCA(iTrial).times =  ((0:79)*mean(diff(Tmtx))+offset(1))';
%         
%         RtrialsLengths(iTrial) = diff(tlimsLong_)+1;
%         
%         tlimsLong_= RightTrialsEvents(iTrial,:)/1e6 + [offset(1), 0 0 0];
%         tlimsLong_= closest(Tmtx,tlimsLong_ );
%         RtrialsTimes(iTrial,:) = tlimsLong_-tlimsLong_(1);
%         nR=nR+1;
%         
%         
%     end
% end
%% PCA on delay period concatenated the lever presses into one vector

iDelay = 3;
Ltrials = [];                   Rtrials = [];
nL = 0;                         nR = 0;
LtrialsTimes = [];              RtrialsTimes = [];
LtrialsLengths = [];            RtrialsLengths = [];

% LeftTrials =[t.Long.SamplePress_LeftCorrect,t.Long.ChoicePress_LeftCorrect];
% RightTrials =[t.Long.SamplePress_RightCorrect,t.Long.ChoicePress_RightCorrect];
LeftTrials =[t.Long.SamplePress_LeftCorrect,t.Long.SamplePress_LeftCorrect+Delays__(iDelay)*1e6];
RightTrials =[t.Long.SamplePress_RightCorrect,t.Long.SamplePress_RightCorrect+Delays__(iDelay)*1e6];

LjPCA=[];
for iTrial =1:size(LeftTrials,1)
    try
        tlimsLong_ = LeftTrials(iTrial,:)/1e6 + offset;
        tlimsLong_ = closest(Tmtx,tlimsLong_);
        FR = iFR_(tlimsLong_(1):tlimsLong_(2),:);
        Ltrials = [Ltrials;FR];
        LjPCA(iTrial).A=FR;
        LjPCA(iTrial).times = (0:size(FR,1)-1)*bw'+offset(1);
        LjPCA(iTrial).times=LjPCA(iTrial).times*1000;
%         LjPCA(iTrial).A=zscore(FR(1:650,:));
%         LjPCA(iTrial).times =  ((0:399)*mean(diff(Tmtx))+offset(1))';
%         LjPCA(iTrial).times =  1000*((0:649)*mean(diff(Tmtx)))';
        LtrialsLengths(iTrial) = diff(tlimsLong_)+1;
        tlimsLong_= LeftTrialsEvents(iTrial,:)/1e6 + [offset(1), 0 0 0];
        tlimsLong_= closest(Tmtx,tlimsLong_ );
        LtrialsTimes(iTrial,:) = tlimsLong_-tlimsLong_(1);
        nL=nL+1;
        
        
    end
end

RjPCA=[];
for iTrial =1:size(RightTrials,1)
    try
        tlimsLong_ = RightTrials(iTrial,:)/1e6 + offset;
        tlimsLong_ = closest(Tmtx,tlimsLong_);
        FR = iFR_(tlimsLong_(1):tlimsLong_(2),:);
        Rtrials = [Rtrials;FR];
        RjPCA(iTrial).A=FR;
        RjPCA(iTrial).times = (0:size(FR,1)-1)*bw'+offset(1);
        RjPCA(iTrial).times=RjPCA(iTrial).times*1000;
%         RjPCA(iTrial).A=FR(1:200,:);
%         RjPCA(iTrial).times =  ((0:79)*mean(diff(Tmtx))+offset(1))';
        
        RtrialsLengths(iTrial) = diff(tlimsLong_)+1;
        
        tlimsLong_= RightTrialsEvents(iTrial,:)/1e6 + [offset(1), 0 0 0];
        tlimsLong_= closest(Tmtx,tlimsLong_ );
        RtrialsTimes(iTrial,:) = tlimsLong_-tlimsLong_(1);
        nR=nR+1;
        
        
    end
end

%% Run PCA
% % input = [Ltrials;Rtrials];
% %     input=zscore(input);    
%         
%         
%         [COEFF,SCORE,~,~,var_explained]=pca(input);
%         percent_var_explained = 100*var_explained/sum(var_explained);
%         PC_L=cell(length(LtrialsLengths),1);
%         PC_R=cell(length(RtrialsLengths),1);
% for iTrial = 1:length(LtrialsLengths)
%     PC_L{iTrial} = SCORE(1:LtrialsLengths(1),1:3);
%     SCORE(1:LtrialsLengths(1),:)=[];
%     LtrialsLengths(1)=[];
% end
% 
% for iTrial = 1:length(RtrialsLengths)
%     PC_R{iTrial} = SCORE(1:RtrialsLengths(1),1:3);
%     SCORE(1:RtrialsLengths(1),:)=[];
%     RtrialsLengths(1)=[];
% end
% 
% figure;hold on
% for  iTrial = 1:length(PC_L)
%      plot(PC_L{iTrial}(:,1),PC_L{iTrial}(:,2),'color',[0 0 1 0.2])
%      scatter(PC_L{iTrial}(LtrialsTimes(iTrial,[4]),1),PC_L{iTrial}(LtrialsTimes(iTrial,[4]),2),'ob')
% end
% 
% % 
% for  iTrial = 1:length(PC_R)
%     plot(PC_R{iTrial}(:,1),PC_R{iTrial}(:,2),'color',[1 0 0 0.2])
%     scatter(PC_R{iTrial}(RtrialsTimes(iTrial,[4]),1),PC_R{iTrial}(RtrialsTimes(iTrial,[4]),2),'or')
% end        
%% Plot PCA results
% figure;hold on
% set(gcf,'pos',[1197 610 484 369],'Color','w');
% for  iTrial = 1:length(PC_L)
%     [X Y Z] =tubeplot(PC_L{iTrial}(:,1),...
%         PC_L{iTrial}(:,2),...
%         PC_L{iTrial}(:,3),...
%         0.1);
%     
%     fig_=surf(X,Y,Z);
%     set(fig_,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
% end
% for  iTrial = 1:length(PC_R)
%     [X Y Z] =tubeplot(PC_R{iTrial}(:,1),...
%         PC_R{iTrial}(:,2),...
%         PC_R{iTrial}(:,3),...
%         0.1);
%     
%     fig_=surf(X,Y,Z);
%     set(fig_,'EdgeColor','none','FaceColor',[1 0 0],'FaceLighting','phong')% axis([0 3 -2 2 -1 2]);%axis([0 2 -1 1 -1 1])
% end
%     %   set(WT_fig,'FaceColor','interp')
% grid on
% view(3)
% axis vis3d
% camlight
% xlabel('PC 1');ylabel('PC 2');zlabel('PC 3')
%% Run jPCA

jPCA_params.softenNorm = 5;  % how each neuron's rate is normized, see below
jPCA_params.suppressBWrosettes = true;  % these are useful sanity plots, but lets ignore them for now
jPCA_params.suppressHistograms = true;  % these are useful sanity plots, but lets ignore them for now
jPCA_params.numPCs = 10;  % default anyway, but best to be specific

% times = 0:50:32450;  % 50 ms before 'neural movement onset' until 150 ms after
times = -5000:50:5000;  % 50 ms before 'neural movement onset' until 150 ms after
% times = (offset(1):bw:Delays__(iDelay)+offset(2))*1000;  % 50 ms before 'neural movement onset' until 150 ms after


[D.ProjectionL{iArea}, D.SummaryL{iArea}] = jPCA(LjPCA, times, jPCA_params);
phaseSpace(D.ProjectionL{iArea}, D.SummaryL{iArea});  % makes the plot
[D.ProjectionR{iArea}, D.SummaryR{iArea}] = jPCA(RjPCA, times, jPCA_params);
phaseSpace(D.ProjectionR{iArea}, D.SummaryR{iArea});  % makes the plot
% plotParams.planes2plot = [1 2 3];
plotParams.crossCondMean=1;

% phaseSpace(D.ProjectionL{iArea}, D.SummaryL{iArea}, plotParams);  % makes all three plots
% phaseMovie(D.ProjectionL{iArea}, D.SummaryL{iArea});
figure; hold on 
plot(D.SummaryL{1}.crossCondMean(:,1),D.SummaryL{1}.crossCondMean(:,2),'r')
plot(D.SummaryR{1}.crossCondMean(:,1),D.SummaryR{1}.crossCondMean(:,2),'b')

D.LjPCA{iArea} = LjPCA;
D.RjPCA{iArea} = RjPCA;
        %% Save results
        fnOut = sprintf('%s\\jPCAresults\\%s_%s_jPCA.mat',pat,fname,Areas{iArea});
        save(fnOut,'D');
        fprintf('Done.\n')
        
        
        iFR_=iFR;
        % iFR_=zscore(iFR_);
        
    end
end
% clearvars -except Delays_ tbShort tbMedium tbLong tbANOVA  pat fileList fileListAss Areas tlimsANOVA shift

clear 
Targets = {'SHORT','MEDIUM','LONG'};
Delays_ = {'Short','Medium','Long'};
tlimsShort=[-10 5];
tlimsMedium=[-10 5];
tlimsLong=[-10 5];
tlimsANOVA = [-4 4];
shift = 0;
plotOnline = false;
bw=0.05;
tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbANOVA=tlimsANOVA(1):bw:tlimsANOVA(2);
clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
else ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
end
cd(pat)
fileList=dir(sprintf('allTimestamps%s*LONG*.mat',filesep));



%%%%%
Areas = {'PFC','HP'};
iFile = 6%:length(fileList)
iArea = 1%:length(Areas)
%% %%%%% Load an example
fname=strtok(fileList(iFile).name,'_');

% load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
% load(sprintf('%s\\%s.mat',pat,fname));
% 
% fprintf('Analysing run %d/%d %s (%s)...',iFile,length(fileList),fname,Areas{iArea})
% % load(sprintf('%s\\KDE_bins\\%s_%s_iFR50.mat',pat,fname,Areas{iArea}),'iFR','Tmtx');
% load(sprintf('%s\\KDE_binsTaskonly\\%s_%s_iFR50_behavOnly.mat',pat,fname,Areas{iArea}),'iFR','Tmtx');

load(sprintf('%s%sallTimestamps%s%s_Events.mat',pat,filesep,filesep,fname));
load(sprintf('%s%s%s.mat',pat,filesep,fname));

fprintf('Analysing run %d/%d %s (%s)...',iFile,length(fileList),fname,Areas{iArea})
% load(sprintf('%s\\KDE_bins\\%s_%s_iFR50.mat',pat,fname,Areas{iArea}),'iFR','Tmtx');
load(sprintf('%s%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep,filesep,fname,Areas{iArea}),'iFR','Tmtx');

iFR_=iFR;

iFR_=iFR;
% iFR_=zscore(iFR_);

% Sanity check firing rates
% iUnit = 2;
% figure; hold on
% eval(sprintf('ST = %scells{iUnit}.t*1e-4;',Areas{iArea}))
% scatter(ST,zeros(size(ST)),'.b')
% plot(Tmtx(1:end-1),mat2gray(iFR_(:,iUnit)),'b')
% stairs(Tmtx,mat2gray(histc(ST,Tmtx)),'r')

%% ANOVA on L/R S/C C/E
q_names = {'Untuned','Pure only','Mixed only','Mixed and Pure','Any tuning'};
for iDelay =1:length(Delays_)
    eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
    Trials = [t_.SamplePress_LeftCorrect;...
        t_.SamplePress_LeftError';...
        t_.SamplePress_RightCorrect;...
        t_.SamplePress_RightError';...
        t_.ChoicePress_LeftCorrect;...
        t_.ChoicePress_LeftError';...
        t_.ChoicePress_RightCorrect;...
        t_.ChoicePress_RightError'];
    SC = [repmat({'Sample'},length(t_.SamplePress_LeftCorrect),1);...
        repmat({'Sample'},length(t_.SamplePress_LeftError),1);...
        repmat({'Sample'},length(t_.SamplePress_RightCorrect),1);...
        repmat({'Sample'},length(t_.SamplePress_RightError),1);...
        repmat({'Choice'},length(t_.ChoicePress_LeftCorrect),1);...
        repmat({'Choice'},length(t_.ChoicePress_LeftError),1);...
        repmat({'Choice'},length(t_.ChoicePress_RightCorrect),1);...
        repmat({'Choice'},length(t_.ChoicePress_RightError),1)];
    
    CE = [repmat({'Correct'},length(t_.SamplePress_LeftCorrect),1);...
        repmat({'Error'},length(t_.SamplePress_LeftError),1);...
        repmat({'Correct'},length(t_.SamplePress_RightCorrect),1);...
        repmat({'Error'},length(t_.SamplePress_RightError),1);...
        repmat({'Correct'},length(t_.ChoicePress_LeftCorrect),1);...
        repmat({'Error'},length(t_.ChoicePress_LeftError),1);...
        repmat({'Correct'},length(t_.ChoicePress_RightCorrect),1);...
        repmat({'Error'},length(t_.ChoicePress_RightError),1)];
    
    LR = [repmat({'Left'},length(t_.SamplePress_LeftCorrect),1);...
        repmat({'Left'},length(t_.SamplePress_LeftError),1);...
        repmat({'Right'},length(t_.SamplePress_RightCorrect),1);...
        repmat({'Right'},length(t_.SamplePress_RightError),1);...
        repmat({'Left'},length(t_.ChoicePress_LeftCorrect),1);...
        repmat({'Left'},length(t_.ChoicePress_LeftError),1);...
        repmat({'Right'},length(t_.ChoicePress_RightCorrect),1);...
        repmat({'Right'},length(t_.ChoicePress_RightError),1)];
    FR =[];
    for iTrial =1:length(Trials)
        try
            tlims_ =Trials(iTrial)/1e6+tlimsANOVA + shift;
            tlims_ = closest(Tmtx,tlims_);
            FR(iTrial,1:size(iFR_,2))= mean(iFR_(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
        catch
            FR(iTrial,1:size(iFR_,2))= nan(1,size(iFR_,2));
        end
    end
    SC(isnan(sum(FR,2)))=[];
    CE(isnan(sum(FR,2)))=[];
    LR(isnan(sum(FR,2)))=[];
    FR(isnan(sum(FR,2)),:)=[];
    p_ = zeros(size(iFR_,2),7); F_=p_;
    
    
    for iUnit=1:size(iFR_,2)
        [p_(iUnit,:),tbl_] = anovan(FR(:,iUnit),{SC CE LR},'model','full',...
            'varnames',{'Context','Outcome','Position'},'display','off');
        F_(iUnit,:)=[tbl_{2:8,6}];
        % [results,~,~,gnames] = multcompare(stats,'Dimension',[1 2 3])
    end
    p_thresh = p_<0.05;
    q_ = [sum(sum(p_thresh,2)==0),...                                      % No tuning
        sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)==0),...   % Pure tuning only
        sum(sum(p_thresh(:,1:3),2)==0 & sum(p_thresh(:,4:7),2)>0),...    % Mixed tuning only
        sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)>0),...    % Mixed and pure tuning
        sum(sum(p_thresh,2)>0)];                                          % Any tuning
    %             q_ = [sum((sum(p_'<0.05)<1)), ...        % No tuning
    %                 sum((sum(p_(:,1:3)'<0.05)>1)), ... % Pure tuning only
    %                 sum((sum(p_(:,4:7)'<0.05)>1))];    % Mixed tuning
    eval(sprintf('D.%s.ANOVA.p=p_;',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.F=F_;',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.Factors=transpose({tbl_{2:8,1}});',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.prcSig=sum(p_<0.05)./size(iFR_,2)*100;',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.prcTuned=q_./size(iFR_,2)*100;',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.tuning=q_names;',Delays_{iDelay}));
    
    
end

% clear p_ q_ q_names FR SC LR CE p_thresh
%% Make raster + SDF ( Medium and long combined)

D.Long.ANOVA.p<0.05;
% iUnit  = 27 %43%39% 33; %20
iUnit  = 2% 27 %43%39% 33; %20

offset=20;
SF=5;

%%% SAMPLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tlimsLong = [-5 10];



      tlims_ = [trangeleft_sample(:,1) ; ...
          trangeright_sample(:,1)]*1e-6 ;
tlims_ = bsxfun(@plus,tlims_,+5 + tlimsLong);
figure('color','w'); hold on

plot([5 5],[-offset length(tlims_)],'color',[0 1 0 0.6],'LineWidth',2)
plot([25 25],[-offset length(tlims_)],'color',[1 0 0 0.6],'LineWidth',2)
plot([15 15],[-offset length(tlims_)],'k')

for iTrial =1:length(tlims_)
    if strcmp(Areas{iArea},'PFC')
        times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) ; 
    elseif strcmp(Areas{iArea},'HP')
        times_ = restrictFRrange(HPcells{iUnit}.t*1e-4,tlims_(iTrial,:))'  - tlims_(iTrial,1) ; 
    end
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
end


%%%  CHOICE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tlimsLong = [-10 5];
tlims_ = [trangeleft_choice(:,1) ; ...
          trangeright_choice(:,1)]*1e-6 +5 + tlimsLong;

for iTrial =1:length(tlims_)
    if strcmp(Areas{iArea},'PFC')
        times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) +15; 
    elseif strcmp(Areas{iArea},'HP')
        times_ = restrictFRrange(HPcells{iUnit}.t*1e-4,tlims_(iTrial,:))'  - tlims_(iTrial,1) +15; 
    end
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
end

axis off
plot([31 31],[1 size(trangeleft_choice,1)],'color',[0.5 0 0.5],'LineWidth',1.5,'LineWidth',4)
plot([31 31],[1 size(trangeright_choice,1)]+size(trangeleft_choice,1),'color',[0 0 1],'LineWidth',4)

text(31.5,size(trangeleft_choice,1)/2,'L','color',[0.5 0 0.5])
text(31.5,size(trangeright_choice,1)+size(trangeleft_choice,1)-size(trangeright_choice,1)/2,'R','color',[0 0 1])


data_{1}.iFR = iFR_;
[~,cutout]=SelTimesFRs([-5 10],Tmtx,data_,trangeleft_sample(:,1)*1e-6 +5 );
LS = nan(300,1);
for iTrial = 1:size(trangeleft_sample,1)
    try
    LS= [LS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs([-10 5],Tmtx,data_,trangeleft_choice(:,1)*1e-6 +5 );
LC = nan(300,1);
for iTrial = 1:size(trangeleft_choice,1)
    try
    LC= [LC,cutout{1}{iTrial}(:,iUnit)];
    end
end


[~,cutout]=SelTimesFRs([-5 10],Tmtx,data_,trangeright_sample(:,1)*1e-6 +5 );
RS = nan(300,1);
for iTrial = 1:size(trangeright_sample,1)
    try
    RS= [RS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs([-10 5],Tmtx,data_,trangeright_choice(:,1)*1e-6 +5 );
RC = nan(300,1);
for iTrial = 1:size(trangeright_choice,1)
    try
    RC= [RC,cutout{1}{iTrial}(:,iUnit)];
    end
end


Lm = [nanmean(LS'),nanmean(LC')]*SF;Le =[nansem(LS'),nansem(LC')]*SF;
Rm = [nanmean(RS'),nanmean(RC')]*SF;Re =[nansem(RS'),nansem(RC')]*SF;
tb_  = (1:600)*0.05;
% figure; hold on
% ciplot(Lm+Le-mean([Lm,Rm])-offset,Lm-Le-mean([Lm,Rm])-offset,tb_,[0.5 0 0.5],0.8);
% ciplot(Rm+Re-mean([Lm,Rm])-offset,Rm-Re-mean([Lm,Rm])-offset,tb_,[0 0 1],0.8);
ciplot(Lm+Le-offset,Lm-Le-offset,tb_,[0.5 0 0.5],0.8);
ciplot(Rm+Re-offset,Rm-Re-offset,tb_,[0 0 1],0.8);
D.Long.ANOVA.p<0.05;
%% Make raster + SDF loop across units - no errors
offset=30;
offsetpeak = 20;
SF=2;
NBS = 500;
tlimsSample = [-5 5];
tlimsChoice = [-5 5];
% file 6 / iUnit 11
close all     
figure('color','w','Position',  [99 687 1329 243]); 

 for iUnit=3%1:size(iFR_,2)
     clf
     %% Make raster + SDF  - ordered by Choice press latency
     for iDelay = 1:length(Delays_)
         
         subplot(1,3,iDelay);hold on
         Delay_ = Delays_{iDelay};
         % Cue light 
         CL_L = eval(sprintf('t.%s.CueLight_LeftCorrect',Delay_));
         CL_R = eval(sprintf('t.%s.CueLight_RightCorrect',Delay_));
         % Sample Press
         SP_L = eval(sprintf('t.%s.SamplePress_LeftCorrect',Delay_));
         SP_R = eval(sprintf('t.%s.SamplePress_RightCorrect',Delay_));
         % Nosepoke
         NP_L = eval(sprintf('t.%s.NosePoke_LeftCorrect',Delay_));
         NP_R = eval(sprintf('t.%s.NosePoke_RightCorrect',Delay_));
         % Choice Press
         CP_L = eval(sprintf('t.%s.ChoicePress_LeftCorrect',Delay_));
         CP_R = eval(sprintf('t.%s.ChoicePress_RightCorrect',Delay_));
         
         % Sort all trials by nosepoke latency
         d_Lpoke = CP_L-NP_L;
         d_Rpoke = CP_R-NP_R;
         [~,idxL]=sort(d_Lpoke);[~,idxR]=sort(d_Rpoke);
         NP_offset = [d_Lpoke(idxL);d_Rpoke(idxR)]*1e-6;
         NP_offset(NP_offset>abs(tlimsChoice(1)))=NaN;
         
         d_Lcue = SP_L-CL_L;
         d_Rcue = SP_R-CL_R;
         CL_offset = [d_Lcue(idxL);d_Rcue(idxR)]*1e-6;
         CL_offset(CL_offset>abs(tlimsSample(1)))=NaN;

         
         % Vertical bars
         plot(sum(abs(tlimsSample))/2*[1,1],[0 length(NP_offset)+1],'color',[0 1 0 0.6],'LineWidth',2)
         plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0 length(NP_offset)+1],'color',[1 0 0 0.6],'LineWidth',2)
         plot(sum(abs(tlimsSample))*[1,1],[0 length(NP_offset)+1],'k')
         
         plot(sum(abs(tlimsSample))/2*[1,1],[-offset -offsetpeak],'color',[0 1 0 0.6],'LineWidth',2)
         plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[-offset -offsetpeak],'color',[1 0 0 0.6],'LineWidth',2)
         plot(sum(abs(tlimsSample))*[1,1],[-offset -offsetpeak],'k')
         
         %%%  rasters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         eval(sprintf('ST = %scells{iUnit}.t*1e-4;',Areas{iArea}))
         
         %%% SAMPLE %%%

         tlims_ = [SP_L(idxL) ; ...
                   SP_R(idxR)]*1e-6 + tlimsSample;
         
         for iTrial =1:length(tlims_)
             times_ = restrictFRrange(ST,tlims_(iTrial,:))' - tlims_(iTrial,1) ;
             plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
             plot(sum(abs(tlimsSample))/2-[CL_offset(iTrial)';CL_offset(iTrial)'],[iTrial iTrial]+[-0.5; 0.5],'color',[0.9 0.4 0.2],'LineWidth',1.5)

         end
         
         %%% CHOICE %%%

         
         tlims_ = [CP_L(idxL) ; ...
                   CP_R(idxR)]*1e-6 + tlimsChoice;
         
         for iTrial =1:length(tlims_)
             times_ = restrictFRrange(ST,tlims_(iTrial,:))' - tlims_(iTrial,1) + sum(abs(tlimsSample));
             plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
             plot(-[NP_offset(iTrial)';NP_offset(iTrial)']+15,[iTrial iTrial]+[-0.5; 0.5],'color',[0.9 0.4 0.2],'LineWidth',1.5)
         end
         
         % Trial markers
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_L,1)],'color',[0.5 0 0.5],'LineWidth',1.5,'LineWidth',4)
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_R,1)]+size(CP_L,1),'color',[0 0 1],'LineWidth',4)
         % Trial markers text
         text(1.5+sum(abs([tlimsSample tlimsChoice])),size(CP_L,1)/2,'L','color',[0.5 0 0.5])
         text(1.5+sum(abs([tlimsSample tlimsChoice])),size(CP_L,1)+size(CP_R,1)/2,'R','color',[0 0 1])
         
         % Spike rates
         data_{1}.iFR = iFR_;
         
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_L*1e-6 );
         LS = nan(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_L,1)
             try
                 LS= [LS,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_L*1e-6 );
         LC = nan(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(CP_L,1)
             try
                 LC= [LC,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         
         
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_R*1e-6 );
         RS  = nan(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_R,1)
             try
                 RS= [RS,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_R*1e-6);
         RC= nan(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(SP_R,1)
             try
                 RC= [RC,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         Lm = [nanmean(LS'),nanmean(LC')]*SF;Le =[nansem(LS'),nansem(LC')]*SF;
         Rm = [nanmean(RS'),nanmean(RC')]*SF;Re =[nansem(RS'),nansem(RC')]*SF;
         Lm = (Lm-mean(Lm))-std(Lm);Le = (Le-mean(Lm))-std(Lm);
         Rm = (Rm-mean(Rm))-std(Rm);Re = (Re-mean(Rm))-std(Rm);
         tb_  = (1:sum(abs([tlimsSample tlimsChoice]))/0.05)*0.05;
         fr{1} = [LS(:,2:end);LC(:,2:end)];
         fr{2} = [RS(:,2:end);RC(:,2:end)];
         
         [p,~] = permtest2vec(fr{1},fr{2},NBS,0.05);

        % Add statistical significance ticks
       
         fr_{1}=[];fr_{2}=[];
         for iTrial=1:size(fr{1},2)
             fr_{1} = [fr_{1};fr{1}(:,iTrial)];
         end
         for iTrial=1:size(fr{2},2)
             fr_{2} = [fr_{2};fr{2}(:,iTrial)];
         end                 
         [~,~,~,~,~,~,TS,~,~]=DecodeStats([fr_{1};fr_{2}],[ones(size(fr{1},2),1);1+ones(size(fr{2},2),1)],0.05);
         p = tpdf(TS,(size(fr{1},2)+size(fr{1},2))-2);
         p = p<0.05;   
  
         a = nan(size(p));a(p)=1;
         a(Lm==Rm)=NaN;
      

         ciplot(Lm+Le-offset,Lm-Le-offset,tb_,[0.5 0 0.5],0.8);
         ciplot(Rm+Re-offset,Rm-Re-offset,tb_,[0 0 1],0.8);
         
         plot(tb_,a - offsetpeak,'-k','LineWidth',2.5)

         axis tight;
         % axis off
     end
     
     %%% Axis adjustments
     for iDelay = 1:length(Delays_)
         subplot(1,length(Delays_),iDelay);
         ylims(iDelay,:) = get(gca,'YLim');
     end
     ylims = [floor(min(ylims(:,1))),ceil(max(ylims(:,2)))];
     for iDelay = 1:length(Delays_)
         subplot(1,length(Delays_),iDelay);
         set(gca,'YLim',ylims);
         axis off
     end
%      subplot(1,length(Delays_),iDelay);
     plot([19 20],[-15 -15],'-k','LineWidth',1.5)
     plot([20 20],[-15 -10],'-k','LineWidth',1.5)
%      printFigs(1, pwd, '-djpeg', sprintf('%s_%s_unit%d',fname,Areas{iArea},iUnit))
     
     
%      clear SP_L SP_R  NP_L  NP_R CP_L CP_R d_L d_R LC RC RS RC fr fr_ TS
%      clear idxL idxR NP_offset tlims_ times_ iDelay iTrial ylims
%      clear Lm Le Rm Re tb_ a
 end
 
%% Make raster + SDF loop across units - with errors
errGap = 1;
offset=20;
offsetpeak = 10;
SF=2;
NBS = 500;
tlimsSample = [-5 5];
tlimsChoice = [-5 5];   
SigBars = [-41 -47 -53];
SigBarsWidth = 3;
SigBarsWidthAlpha = 1;
% file 6 / iUnit 11
close all     
figure('color','w','Position',  [99 300 1329 300]); 

 for iUnit=3%1:size(iFR_,2)
     spikes_ = eval(sprintf('%scells{iUnit}.t*1e-4;',Areas{iArea}));

     clf
     %% Make raster + SDF  - ordered by Choice press latency
     for iDelay = 1:length(Delays_)
         %%
%          subplot(1,3,iDelay);hold on
         subaxis(1,3,iDelay,'SpacingHoriz',0.01,'Margin',0.01);hold on
         Delay_ = Delays_{iDelay};
         Cue_L= eval(sprintf('t.%s.CueLight_LeftCorrect',Delay_));
         Cue_R= eval(sprintf('t.%s.CueLight_RightCorrect',Delay_));
         SP_L = eval(sprintf('t.%s.SamplePress_LeftCorrect',Delay_));
         SP_R = eval(sprintf('t.%s.SamplePress_RightCorrect',Delay_));
         NP_L = eval(sprintf('t.%s.NosePoke_LeftCorrect',Delay_));
         NP_R = eval(sprintf('t.%s.NosePoke_RightCorrect',Delay_));
         CP_L = eval(sprintf('t.%s.ChoicePress_LeftCorrect',Delay_));
         CP_R = eval(sprintf('t.%s.ChoicePress_RightCorrect',Delay_));
         Rew_L = eval(sprintf('t.%s.RewardConsume_LeftCorrect',Delay_));
         Rew_R = eval(sprintf('t.%s.RewardConsume_RightCorrect',Delay_));
         
         Cue_Le= eval(sprintf('t.%s.CueLight_LeftError',Delay_))';
         Cue_Re= eval(sprintf('t.%s.CueLight_RightError',Delay_))';         
         SP_Le = eval(sprintf('t.%s.SamplePress_LeftError',Delay_))';
         SP_Re = eval(sprintf('t.%s.SamplePress_RightError',Delay_))';
         NP_Le = eval(sprintf('t.%s.NosePoke_LeftError',Delay_))';
         NP_Re = eval(sprintf('t.%s.NosePoke_RightError',Delay_))';
         CP_Le = eval(sprintf('t.%s.ChoicePress_LeftError',Delay_))';
         CP_Re = eval(sprintf('t.%s.ChoicePress_RightError',Delay_))';
         
         % Sample to press latency (correct)
         d_Lsample = SP_L-Cue_L;
         d_Rsample = SP_R-Cue_R;
         
         % nosepoke to choice press latency (correct)
         d_Lpoke = CP_L-NP_L;
         d_Rpoke = CP_R-NP_R;
         
         % Sample to press latency (correct)
         d_Lpoke_e = CP_Le-NP_Le;
         d_Rpoke_e = CP_Re-NP_Re;
         
         % Sample to press latency (correct)
         d_Lsample_e = SP_Le-Cue_Le;
         d_Rsample_e = SP_Re-Cue_Re;
         
         % nosepoke to choice press latency (correct)
         d_LRew = Rew_L-CP_L;
         d_RRew = Rew_R-CP_R;
          
%          % Sort trials on nosepoke to choice latency
         [~,idxL]=sort(d_Lpoke);[~,idxR]=sort(d_Rpoke);
         [~,idxLe]=sort(d_Lpoke_e);[~,idxRe]=sort(d_Rpoke_e);
         
         % Sort trials on nosepoke to choice latency
%          [~,idxL]=sort(d_Lsample);[~,idxR]=sort(d_Rsample);
%          [~,idxLe]=sort(d_Lsample_e);[~,idxRe]=sort(d_Rsample_e);
         
         
         Sampl_offset = [d_Lsample(idxL);d_Rsample(idxR);d_Lsample_e(idxLe);d_Rsample_e(idxRe)]*1e-6;
         Sampl_offset(Sampl_offset>abs(tlimsSample(1)))=NaN;

         NP_offset = [d_Lpoke(idxL);d_Rpoke(idxR);d_Lpoke_e(idxLe);d_Rpoke_e(idxRe)]*1e-6;
         NP_offset(NP_offset>abs(tlimsChoice(1)))=NaN;
         
         Rew_offset = [d_LRew(idxL);d_RRew(idxR)]*1e-6;
         Rew_offset(Rew_offset>abs(tlimsChoice(2)))=NaN;
         
         
         nTrials = [0,length(SP_L), length(SP_R), length(SP_Le), length(SP_Re)];

         tlims_ = [SP_L(idxL) ; SP_R(idxR) ; SP_Le(idxLe) ; SP_Re(idxRe)]*1e-6 + tlimsSample;
         
         
         % Raster: sample trials
         for iTrial =1:sum(nTrials)
             times_ = restrictFRrange(spikes_,tlims_(iTrial,:))' - tlims_(iTrial,1) ;
             if iTrial<=sum(nTrials(1:3))
                 plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
             else
                 plot([times_;times_],errGap + iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0.6 0.6 0.6 0.6],'LineWidth',1.5)
             end
             plot(abs(tlimsSample(1))-[Sampl_offset(iTrial)';Sampl_offset(iTrial)'],[iTrial iTrial]+[-0.5; 0.5],'color',[0 1 0 0.6],'LineWidth',1.5)
            
         end
         
         tlims_ = [CP_L(idxL) ; CP_R(idxR) ; CP_Le(idxLe) ; CP_Re(idxRe)]*1e-6 + tlimsChoice;
         
         for iTrial =1:sum(nTrials)
             times_ = restrictFRrange(spikes_,tlims_(iTrial,:))' - tlims_(iTrial,1) + sum(abs(tlimsSample));
             if iTrial<=sum(nTrials(1:3))
             	plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
                plot(-[NP_offset(iTrial)';NP_offset(iTrial)']+15,[iTrial iTrial]+[-0.5; 0.5],'color',[0.9 0.4 0.2],'LineWidth',1.5)
                plot([Rew_offset(iTrial)';Rew_offset(iTrial)']+15,[iTrial iTrial]+[-0.5; 0.5],'color',[0.9 0.4 0.2],'LineWidth',1.5)
             else
                plot([times_;times_],errGap + iTrial*ones(2,length(times_))+[-0.5; 0.5],   'color',[0.6 0.6 0.6 0.6],'LineWidth',1.5)            
                plot(-[NP_offset(iTrial)';NP_offset(iTrial)']+15,errGap + [iTrial iTrial]+[-0.5; 0.5],'color',[0.9 0.4 0.2],'LineWidth',1.5)
             end
         end
         
         % bars for correct rasters
         plot(sum(abs(tlimsSample))/2*[1,1],[0.5 sum(nTrials(1:3))+0.5],'color',[0 1 0 0.6],'LineWidth',2)
         plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0.5 sum(nTrials(1:3))+0.5],'color',[1 0 0 0.6],'LineWidth',2)
         plot(sum(abs(tlimsSample))*[1,1],[0.5 sum(nTrials(1:3))+0.5],'color',[0 0 0 0.8],'LineWidth',1)
         
         % bars for error rasters
         plot(sum(abs(tlimsSample))/2*[1,1],...
             errGap+[sum(nTrials(1:3))-0.5+1 sum(nTrials)+0.5],'color',[0 1 0 0.6],'LineWidth',2)
         plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],...
             errGap+[sum(nTrials(1:3))-0.5+1 sum(nTrials)+0.5],'color',[1 0 0 0.6],'LineWidth',2)
         plot(sum(abs(tlimsSample))*[1,1],...
             errGap+[sum(nTrials(1:3))-0.5+1 sum(nTrials)+0.5],'color',[0 0 0 0.8],'LineWidth',1)
         
         % bars for correct spike rates
         plot(sum(abs(tlimsSample))/2*[1,1],[-offset-10 -offsetpeak],'color',[0 1 0 0.6],'LineWidth',2)
         plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[-offset-10 -offsetpeak],'color',[1 0 0 0.6],'LineWidth',2)
         plot(sum(abs(tlimsSample))*[1,1],[-offset-10 -offsetpeak],'color',[0 0 0 0.8],'LineWidth',1)
         
         
         
         
         
         
           % Trial markers
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_L,1)]+[-0.5 0.5],'color',[0.5 0 0.5],'LineWidth',1.5,'LineWidth',4)
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_R,1)]+size(CP_L,1)+[-0.5 0.5],'color',[0 0 1],'LineWidth',4)
         
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_Le,1)]+sum(nTrials(1:3))+errGap+[-0.5 0.5],'color',[0.5 0 0.5 0.6],'LineWidth',1.5,'LineWidth',4)
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_Re,1)]+sum(nTrials(1:4))+errGap+[-0.5 0.5],'color',[0 0 1 0.6],'LineWidth',4)
         
         % Trial markers text
         if iDelay==length(Delays_)
         text(1.5+sum(abs([tlimsSample tlimsChoice])),size(CP_L,1)/2,'L correct','color',[0.5 0 0.5])
         text(1.5+sum(abs([tlimsSample tlimsChoice])),size(CP_L,1)+size(CP_R,1)/2,'R correct','color',[0 0 1])   
                  
         text(1.5+sum(abs([tlimsSample tlimsChoice])),sum(nTrials(1:3))+errGap+nTrials(4)/2,'L errors','color',[0.5 0 0.5])
         text(1.5+sum(abs([tlimsSample tlimsChoice])),sum(nTrials(1:4))+errGap+nTrials(5)/2,'R errors','color',[0 0 1])
         end
         
         data_{1}.iFR = iFR_;
         % Spike rates (correct trials)
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_L*1e-6 );
         LS = zeros(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_L,1)
             try
                 LS= [LS,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_L*1e-6 );
         LC = zeros(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(CP_L,1)
             try
                 LC= [LC,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         
         
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_R*1e-6 );
         RS  = zeros(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_R,1)
             try
                 RS= [RS,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_R*1e-6);
         RC= zeros(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(SP_R,1)
             try
                 RC= [RC,cutout{1}{iTrial}(:,iUnit)];
             end
         end
        
         % Spike rates (error trials)
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_Le*1e-6 );    
         LSe = zeros(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_Le,1)
             if ~isempty(cutout{1}{iTrial}(:,iUnit))
                 LSe= [LSe,cutout{1}{iTrial}(:,iUnit)];
             else
                 LSe= [LSe,zeros(sum(abs([tlimsSample]))/0.05,1)];
             end
         end
         
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_Le*1e-6 );
         LCe = zeros(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(CP_Le,1)
             if ~isempty(cutout{1}{iTrial}(:,iUnit))
                 LCe= [LCe,cutout{1}{iTrial}(:,iUnit)];
             else
                 LCe= [LCe,zeros(sum(abs([tlimsSample]))/0.05,1)];
             end
         end
         
         
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_Re*1e-6 );
         RSe  = zeros(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_Re,1)
             if ~isempty(cutout{1}{iTrial}(:,iUnit))
                 RSe= [RSe,cutout{1}{iTrial}(:,iUnit)];
             else
                 RSe= [RSe,zeros(sum(abs([tlimsSample]))/0.05,1)];
             end
         end
         
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_Re*1e-6);
         RCe= zeros(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(CP_Re,1)
             if ~isempty(cutout{1}{iTrial}(:,iUnit))
                 RCe= [RCe,cutout{1}{iTrial}(:,iUnit)];
             else
                 RCe= [RCe,zeros(sum(abs([tlimsSample]))/0.05,1)];
             end
         end         
         
         LS(isnan(LS))=0; RS(isnan(RS))=0;
         LC(isnan(LC))=0; RC(isnan(RC))=0;
         LSe(isnan(LSe))=0; RSe(isnan(RSe))=0;
         LCe(isnan(LCe))=0; RCe(isnan(RCe))=0;
         
         Lm = [nanmean(LS,2);nanmean(LC,2);nanmean(LSe,2);nanmean(LCe,2)]*SF;
         Le = [nansem(LS,2);nansem(LC,2);nansem(LSe,2);nansem(LCe,2)]*SF;
         Rm = [nanmean(RS,2);nanmean(RC,2);nanmean(RSe,2);nanmean(RCe,2)]*SF;
         Re = [nansem(RS,2);nansem(RC,2);nansem(RSe,2);nansem(RCe,2)]*SF;
         Lm = (Lm-mean(Lm))-std(Lm);Le = (Le-mean(Lm))-std(Lm);
         Rm = (Rm-mean(Rm))-std(Rm);Re = (Re-mean(Rm))-std(Rm);
         tb_  = (1:sum(abs([tlimsSample tlimsChoice]))/0.05)*0.05;
         
         % Spike rates (errors)
         if size(LCe,2)>1
            plot(tb_,Lm((length(tb_)+1):length(tb_)*2)-offset,'color',[0.5 0 0.5],'LineStyle',':','LineWidth',1.2);
         end
         if size(RCe,2)>1
            plot(tb_,Rm((length(tb_)+1):length(tb_)*2)-offset,'color',[0 0 1],'LineStyle',':','LineWidth',1.2);
         end
         % Spike rates (correct)
         ciplot(Lm(1:length(tb_))+Le(1:length(tb_))-offset,Lm(1:length(tb_))-Le(1:length(tb_))-offset,tb_,[0.5 0 0.5],0.8);
         ciplot(Rm(1:length(tb_))+Re(1:length(tb_))-offset,Rm(1:length(tb_))-Re(1:length(tb_))-offset,tb_,[0 0 1],0.8);

         

         
         
         % Add statistical significance ticks - L vs R
         fr{1} = [LS(:,2:end);LC(:,2:end)];
         fr{2} = [RS(:,2:end);RC(:,2:end)];
                        
         fr_{1}=[];fr_{2}=[];
         for iTrial=1:size(fr{1},2)
             fr_{1} = [fr_{1};fr{1}(:,iTrial)];
         end
         for iTrial=1:size(fr{2},2)
             fr_{2} = [fr_{2};fr{2}(:,iTrial)];
         end                 
         [~,~,~,~,~,~,TS,~,~]=DecodeStats([fr_{1};fr_{2}],[ones(size(fr{1},2),1);1+ones(size(fr{2},2),1)],0.05);
         a=[];p=[];
         p = tpdf(TS,(size(fr{1},2)+size(fr{2},2))-2);
         p = p<0.05;   
  
         a = nan(size(p));a(p)=1;
         a(Lm==Rm)=NaN;                  
%          plot(tb_,a -45,'-k','LineWidth',4)
         plot(tb_,repmat(SigBars(1),length(tb_),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',SigBarsWidth)

         sig1 =  (Lm(1:length(tb_))>Rm(1:length(tb_)) & p); 
         sig2 =  (Lm(1:length(tb_))<Rm(1:length(tb_)) & p); 
         a = nan(size(sig1));a(sig1)=1;
         plot(tb_,a*SigBars(1),'color',[0.5 0 0.5],'LineWidth',SigBarsWidth)
         a = nan(size(sig2));a(sig2)=1;
        plot(tb_,a*SigBars(1),'color',[0 0 1],'LineWidth',SigBarsWidth)
    
        
          % Add statistical significance ticks - C vs E
         fr{1} = [[LS(:,2:end);LC(:,2:end)],[RS(:,2:end);RC(:,2:end)]];
         fr{2} = [[LSe(:,2:end);LCe(:,2:end)],[RSe(:,2:end);RCe(:,2:end)]];
                        
         fr_{1}=[];fr_{2}=[];
         for iTrial=1:size(fr{1},2)
             fr_{1} = [fr_{1};fr{1}(:,iTrial)];
         end
         for iTrial=1:size(fr{2},2)
             fr_{2} = [fr_{2};fr{2}(:,iTrial)];
         end                 
         [~,~,~,~,~,~,TS,~,~]=DecodeStats([fr_{1};fr_{2}],[ones(size(fr{1},2),1);1+ones(size(fr{2},2),1)],0.05);
         a=[];p=[];
         p = tpdf(TS,(size(fr{1},2)+size(fr{2},2))-2);
         p = p<0.05;   
  
         a = nan(size(p));a(p)=1;
%          a(Lm==Rm)=NaN;                  
         plot(tb_,repmat(SigBars(2),length(tb_),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',SigBarsWidth)         
         plot(tb_,a*SigBars(2),'color',[0 0 0 SigBarsWidthAlpha],'LineWidth',SigBarsWidth)

         
          % Add statistical significance ticks - Sample vs Choice
         fr{1} = [LS(:,2:end),RS(:,2:end)];
         fr{2} = [LC(:,2:end),RC(:,2:end)];
%          fr{1} = [LS(:,2:end),RS(:,2:end),LSe(:,2:end),RSe(:,2:end)];
%          fr{2} = [LC(:,2:end),RC(:,2:end),LCe(:,2:end),RCe(:,2:end)];                        
         fr_{1}=[];fr_{2}=[];
         for iTrial=1:size(fr{1},2)
             fr_{1} = [fr_{1};fr{1}(:,iTrial)];
         end
         for iTrial=1:size(fr{2},2)
             fr_{2} = [fr_{2};fr{2}(:,iTrial)];
         end                 
         [~,~,~,~,~,~,TS,~,~]=DecodeStats([fr_{1};fr_{2}],[ones(size(fr{1},2),1);1+ones(size(fr{2},2),1)],0.05);
         TS = [TS;TS];
         a=[];p=[];
         p = tpdf(TS,(size(fr{1},2)+size(fr{2},2))-2);
         p = p<0.05;   
  
         a = nan(size(p));a(p)=1;
         a(Lm==Rm)=NaN;                  
         plot(tb_,repmat(SigBars(3),length(tb_),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',SigBarsWidth)
%          plot(tb_,a*SigBars(3),'color',[0 0 0 SigBarsWidthAlpha],'LineWidth',SigBarsWidth)

         fr =cellfun(@mean,cellfun(@transpose,fr,'UniformOutput',0),'UniformOutput',0);
         
         sig1 =  (repmat(fr{1}>fr{2},1,2)' & p); 
         sig2 =  (repmat(fr{1}<fr{2},1,2)' & p); 
         a = nan(size(sig1));a(sig1)=1;
         plot(tb_,a*SigBars(3),'color',[0 1 0 SigBarsWidthAlpha],'LineWidth',SigBarsWidth)
         a = nan(size(sig2));a(sig2)=1;
         plot(tb_,a*SigBars(3),'color',[1 0 0 SigBarsWidthAlpha],'LineWidth',SigBarsWidth)
         plot([round(median(tb_)) round(median(tb_))],[1 -1] + SigBars([1,3]),'color',[0 0 0 0.8],'LineWidth',1)
     end
     
     % Match Y scale across all delay lengths
     for iDelay = 1:length(Delays_)
%          subplot(1,length(Delays_),iDelay);
        subaxis(1,3,iDelay,'SpacingHoriz',0.01,'Margin',0.01);hold on
            ylims(iDelay,:) = get(gca,'YLim');
     end
     ylims = [floor(min(ylims(:,1))),ceil(max(ylims(:,2)))];
     for iDelay = 1:length(Delays_)
%          subplot(1,length(Delays_),iDelay);
            subaxis(1,3,iDelay,'SpacingHoriz',0.01,'Margin',0.01);hold on
            set(gca,'YLim',ylims);
         axis off
     end
%      subplot(1,length(Delays_),iDelay);
     plot([19 20],[-15 -15],'-k','LineWidth',1.5)
     plot([20 20],[-15 -10],'-k','LineWidth',1.5)
%      printFigs(1, pwd, '-djpeg', sprintf('%s_%s_unit%d_withErorrs',fname,Areas{iArea},iUnit))
 end
 
 %% Make raster + SDF loop across units - with errors but no ANOVA bars
errGap = 1;
offset=20;
offsetpeak = 10;
SF=1;
NBS = 500;
tlimsSample = [-5 5];
tlimsChoice = [-5 5];   
SigBars = [-41 -47 -53];
SigBarsWidth = 3;
SigBarsWidthAlpha = 1;
behavTickalpha = 1;
col_Cue = [0 1 1 behavTickalpha];
col_NP  = [1 0.5 0 behavTickalpha];
col_Rew = [1 0 1 behavTickalpha];
% file 6 / iUnit 11
close all     
figure('color','w','Position',  [99 300 1329 300]); 

 for iUnit=1:size(iFR_,2)
     spikes_ = eval(sprintf('%scells{iUnit}.t*1e-4;',Areas{iArea}));

     clf
     %% Make raster + SDF  - ordered by Choice press latency
     for iDelay = 1:length(Delays_)
         %%
         subaxis(1,3,iDelay,'SpacingHoriz',0.01,'Margin',0.01);hold on
         Delay_ = Delays_{iDelay};
         Cue_L= eval(sprintf('t.%s.CueLight_LeftCorrect',Delay_));
         Cue_R= eval(sprintf('t.%s.CueLight_RightCorrect',Delay_));
         SP_L = eval(sprintf('t.%s.SamplePress_LeftCorrect',Delay_));
         SP_R = eval(sprintf('t.%s.SamplePress_RightCorrect',Delay_));
         NP_L = eval(sprintf('t.%s.NosePoke_LeftCorrect',Delay_));
         NP_R = eval(sprintf('t.%s.NosePoke_RightCorrect',Delay_));
         CP_L = eval(sprintf('t.%s.ChoicePress_LeftCorrect',Delay_));
         CP_R = eval(sprintf('t.%s.ChoicePress_RightCorrect',Delay_));
         Rew_L = eval(sprintf('t.%s.RewardConsume_LeftCorrect',Delay_));
         Rew_R = eval(sprintf('t.%s.RewardConsume_RightCorrect',Delay_));
         
         Cue_Le= eval(sprintf('t.%s.CueLight_LeftError',Delay_))';
         Cue_Re= eval(sprintf('t.%s.CueLight_RightError',Delay_))';         
         SP_Le = eval(sprintf('t.%s.SamplePress_LeftError',Delay_))';
         SP_Re = eval(sprintf('t.%s.SamplePress_RightError',Delay_))';
         NP_Le = eval(sprintf('t.%s.NosePoke_LeftError',Delay_))';
         NP_Re = eval(sprintf('t.%s.NosePoke_RightError',Delay_))';
         CP_Le = eval(sprintf('t.%s.ChoicePress_LeftError',Delay_))';
         CP_Re = eval(sprintf('t.%s.ChoicePress_RightError',Delay_))';
         
         % Sample to press latency (correct)
         d_Lsample = SP_L-Cue_L;
         d_Rsample = SP_R-Cue_R;
         
         % nosepoke to choice press latency (correct)
         d_Lpoke = CP_L-NP_L;
         d_Rpoke = CP_R-NP_R;
         
         % Sample to press latency (correct)
         d_Lpoke_e = CP_Le-NP_Le;
         d_Rpoke_e = CP_Re-NP_Re;
         
         % Sample to press latency (correct)
         d_Lsample_e = SP_Le-Cue_Le;
         d_Rsample_e = SP_Re-Cue_Re;
         
         % nosepoke to choice press latency (correct)
         d_LRew = Rew_L-CP_L;
         d_RRew = Rew_R-CP_R;
          
%          % Sort trials on nosepoke to choice latency
         [~,idxL]=sort(d_Lpoke);[~,idxR]=sort(d_Rpoke);
         [~,idxLe]=sort(d_Lpoke_e);[~,idxRe]=sort(d_Rpoke_e);
         
        % Sort trials on cue to sample latency
        % [~,idxL]=sort(d_Lsample);[~,idxR]=sort(d_Rsample);
        % [~,idxLe]=sort(d_Lsample_e);[~,idxRe]=sort(d_Rsample_e);
         
         
         Sampl_offset = [d_Lsample(idxL);d_Rsample(idxR);d_Lsample_e(idxLe);d_Rsample_e(idxRe)]*1e-6;
         Sampl_offset(Sampl_offset>abs(tlimsSample(1)))=NaN;

         NP_offset = [d_Lpoke(idxL);d_Rpoke(idxR);d_Lpoke_e(idxLe);d_Rpoke_e(idxRe)]*1e-6;
         NP_offset(NP_offset>abs(tlimsChoice(1)))=NaN;
         
         Rew_offset = [d_LRew(idxL);d_RRew(idxR)]*1e-6;
         Rew_offset(Rew_offset>abs(tlimsChoice(2)))=NaN;
         
         
         nTrials = [0,length(SP_L), length(SP_R), length(SP_Le), length(SP_Re)];

         % Raster: sample trials

         tlims_ = [SP_L(idxL) ;   ...
                   SP_R(idxR) ;   ...
                   SP_Le(idxLe) ; ...
                   SP_Re(idxRe)]*1e-6 + tlimsSample;
         
         
         for iTrial =1:sum(nTrials)
             times_ = restrictFRrange(spikes_,tlims_(iTrial,:))' - tlims_(iTrial,1) ;
             if iTrial<=sum(nTrials(1:3))
                 plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
             else
                 plot([times_;times_],errGap + iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0.6 0.6 0.6 0.6],'LineWidth',1.5)
             end
             plot(abs(tlimsSample(1))-[Sampl_offset(iTrial)';Sampl_offset(iTrial)'],[iTrial iTrial]+[-0.5; 0.5],'color',col_Cue,'LineWidth',1.5)
            
         end
         
         % Raster: choice trials
         
         tlims_ = [CP_L(idxL) ;   ...
                   CP_R(idxR) ;   ...
                   CP_Le(idxLe) ; ...
                   CP_Re(idxRe)]*1e-6 + tlimsChoice;
         
         for iTrial =1:sum(nTrials)
             times_ = restrictFRrange(spikes_,tlims_(iTrial,:))' - tlims_(iTrial,1) + sum(abs(tlimsSample));
             if iTrial<=sum(nTrials(1:3))
             	plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
                plot(-[NP_offset(iTrial)';NP_offset(iTrial)']+15,[iTrial iTrial]+[-0.5; 0.5],'color',col_NP,'LineWidth',1.5)
                plot([Rew_offset(iTrial)';Rew_offset(iTrial)']+15,[iTrial iTrial]+[-0.5; 0.5],'color',col_Rew,'LineWidth',1.5)
             else
                plot([times_;times_],errGap + iTrial*ones(2,length(times_))+[-0.5; 0.5],   'color',[0.6 0.6 0.6 0.6],'LineWidth',1.5)            
                plot(-[NP_offset(iTrial)';NP_offset(iTrial)']+15,errGap + [iTrial iTrial]+[-0.5; 0.5],'color',col_NP,'LineWidth',1.5)
             end
         end
         
         % bars for correct rasters
         plot(sum(abs(tlimsSample))/2*[1,1],[0.5 sum(nTrials(1:3))+0.5],'color',[0 1 0 0.6],'LineWidth',2)
         plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[0.5 sum(nTrials(1:3))+0.5],'color',[1 0 0 0.6],'LineWidth',2)
         plot(sum(abs(tlimsSample))*[1,1],[0.5 sum(nTrials(1:3))+0.5],'color',[0 0 0 0.8],'LineWidth',1)
         
         % bars for error rasters
         plot(sum(abs(tlimsSample))/2*[1,1],...
             errGap+[sum(nTrials(1:3))-0.5+1 sum(nTrials)+0.5],'color',[0 1 0 0.6],'LineWidth',2)
         plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],...
             errGap+[sum(nTrials(1:3))-0.5+1 sum(nTrials)+0.5],'color',[1 0 0 0.6],'LineWidth',2)
         plot(sum(abs(tlimsSample))*[1,1],...
             errGap+[sum(nTrials(1:3))-0.5+1 sum(nTrials)+0.5],'color',[0 0 0 0.8],'LineWidth',1)
         
         % bars for correct spike rates
         plot(sum(abs(tlimsSample))/2*[1,1],[-offset-10 -offsetpeak],'color',[0 1 0 0.6],'LineWidth',2)
         plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[-offset-10 -offsetpeak],'color',[1 0 0 0.6],'LineWidth',2)
         plot(sum(abs(tlimsSample))*[1,1],[-offset-10 -offsetpeak],'color',[0 0 0 0.8],'LineWidth',1)
         
         % Trial markers
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_L,1)]+[-0.5 0.5],'color',[0.5 0 0.5],'LineWidth',1.5,'LineWidth',4)
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_R,1)]+size(CP_L,1)+[-0.5 0.5],'color',[0 0 1],'LineWidth',4)
         
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_Le,1)]+sum(nTrials(1:3))+errGap+[-0.5 0.5],'color',[0.5 0 0.5 0.6],'LineWidth',1.5,'LineWidth',4)
         plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(CP_Re,1)]+sum(nTrials(1:4))+errGap+[-0.5 0.5],'color',[0 0 1 0.6],'LineWidth',4)
         
         % Trial markers text
         if iDelay==length(Delays_)
             text(1.5+sum(abs([tlimsSample tlimsChoice])),size(CP_L,1)/2,'L correct','color',[0.5 0 0.5])
             text(1.5+sum(abs([tlimsSample tlimsChoice])),size(CP_L,1)+size(CP_R,1)/2,'R correct','color',[0 0 1])
             
             text(1.5+sum(abs([tlimsSample tlimsChoice])),sum(nTrials(1:3))+errGap+nTrials(4)/2,'L errors','color',[0.5 0 0.5])
             text(1.5+sum(abs([tlimsSample tlimsChoice])),sum(nTrials(1:4))+errGap+nTrials(5)/2,'R errors','color',[0 0 1])
         end
         
         % Spike rates
         data_{1}.iFR = iFR_;
         
         % Spike rates:  (correct trials)
         % Left sample
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_L*1e-6 );
         LS = zeros(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_L,1)
             try
                 LS= [LS,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         % Left choice
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_L*1e-6 );
         LC = zeros(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(CP_L,1)
             try
                 LC= [LC,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         % Right sample
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_R*1e-6 );
         RS  = zeros(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_R,1)
             try
                 RS= [RS,cutout{1}{iTrial}(:,iUnit)];
             end
         end
         % Right choice
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_R*1e-6);
         RC= zeros(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(SP_R,1)
             try
                 RC= [RC,cutout{1}{iTrial}(:,iUnit)];
             end
         end
        
         % Spike rates: (error trials)
         % Left sample
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_Le*1e-6 );    
         LSe = zeros(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_Le,1)
             if ~isempty(cutout{1}{iTrial}(:,iUnit))
                 LSe= [LSe,cutout{1}{iTrial}(:,iUnit)];
             else
                 LSe= [LSe,zeros(sum(abs([tlimsSample]))/0.05,1)];
             end
         end
         % Left choice
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_Le*1e-6 );
         LCe = zeros(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(CP_Le,1)
             if ~isempty(cutout{1}{iTrial}(:,iUnit))
                 LCe= [LCe,cutout{1}{iTrial}(:,iUnit)];
             else
                 LCe= [LCe,zeros(sum(abs([tlimsSample]))/0.05,1)];
             end
         end
         % Right sample
         [~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,SP_Re*1e-6 );
         RSe  = zeros(sum(abs([tlimsSample]))/0.05,1);
         for iTrial = 1:size(SP_Re,1)
             if ~isempty(cutout{1}{iTrial}(:,iUnit))
                 RSe= [RSe,cutout{1}{iTrial}(:,iUnit)];
             else
                 RSe= [RSe,zeros(sum(abs([tlimsSample]))/0.05,1)];
             end
         end
         % Right choice
         [~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,CP_Re*1e-6);
         RCe= zeros(sum(abs([tlimsChoice]))/0.05,1);
         for iTrial = 1:size(CP_Re,1)
             if ~isempty(cutout{1}{iTrial}(:,iUnit))
                 RCe= [RCe,cutout{1}{iTrial}(:,iUnit)];
             else
                 RCe= [RCe,zeros(sum(abs([tlimsSample]))/0.05,1)];
             end
         end         
         
         LS(isnan(LS))=0; RS(isnan(RS))=0;
         LC(isnan(LC))=0; RC(isnan(RC))=0;
         LSe(isnan(LSe))=0; RSe(isnan(RSe))=0;
         LCe(isnan(LCe))=0; RCe(isnan(RCe))=0;
         
         Lm = [nanmean(LS,2);nanmean(LC,2);nanmean(LSe,2);nanmean(LCe,2)]*SF;
         Le = [nansem(LS,2);nansem(LC,2);nansem(LSe,2);nansem(LCe,2)]*SF;
         Rm = [nanmean(RS,2);nanmean(RC,2);nanmean(RSe,2);nanmean(RCe,2)]*SF;
         Re = [nansem(RS,2);nansem(RC,2);nansem(RSe,2);nansem(RCe,2)]*SF;
         
%          Lm = (Lm-mean([Lm;Rm]));%-std([Lm;Rm]);
%          Le = (Le-mean([Lm;Rm]));%-std([Lm;Rm]);
%          Rm = (Rm-mean([Lm;Rm]));%-std([Lm;Rm]);
%          Re = (Re-mean([Lm;Rm]));%-std([Lm;Rm]);
         tb_  = (1:sum(abs([tlimsSample tlimsChoice]))/0.05)*0.05;
         
         % Spike rates (errors)
         if size(LCe,2)>1
            plot(tb_,Lm((length(tb_)+1):length(tb_)*2)-offset,'color',[0.5 0 0.5],'LineStyle',':','LineWidth',1.2);
         end
         if size(RCe,2)>1
            plot(tb_,Rm((length(tb_)+1):length(tb_)*2)-offset,'color',[0 0 1],'LineStyle',':','LineWidth',1.2);
         end
         % Spike rates (correct)
         hold on
         ciplot(Lm(1:length(tb_))+Le(1:length(tb_))-offset, ...
                Lm(1:length(tb_))-Le(1:length(tb_))-offset, ...
                tb_ , [0.5 0 0.5],0.8);
         ciplot(Rm(1:length(tb_))+Re(1:length(tb_))-offset, ...
                Rm(1:length(tb_))-Re(1:length(tb_))-offset, ...
                tb_ , [0 0 1]   ,0.8);

      
         
         
         % Add statistical significance ticks - L vs R
         fr{1} = [LS(:,2:end);LC(:,2:end)];
         fr{2} = [RS(:,2:end);RC(:,2:end)];
                        
         fr_{1}=[];fr_{2}=[];
         for iTrial=1:size(fr{1},2)
             fr_{1} = [fr_{1};fr{1}(:,iTrial)];
         end
         for iTrial=1:size(fr{2},2)
             fr_{2} = [fr_{2};fr{2}(:,iTrial)];
         end                 
         [~,~,~,~,~,~,TS,~,~]=DecodeStats([fr_{1};fr_{2}],[ones(size(fr{1},2),1);1+ones(size(fr{2},2),1)],0.05);
         a=[];p=[];
         p = tpdf(TS,(size(fr{1},2)+size(fr{2},2))-2);
         p = p<0.05;   
  
         a = nan(size(p));a(p)=1;
         a(Lm(1:length(tb_))==Rm(1:length(tb_)))=NaN;                  
         plot(tb_,a - offsetpeak+1,'-k','LineWidth',2.5)
         
% %          plot(tb_,a -45,'-k','LineWidth',4)
%          plot(tb_,repmat(SigBars(1),length(tb_),1),'color',[0.6 0.6 0.6 0.6],'LineWidth',SigBarsWidth)
% 
%          sig1 =  (Lm(1:length(tb_))>Rm(1:length(tb_)) & p); 
%          sig2 =  (Lm(1:length(tb_))<Rm(1:length(tb_)) & p); 
%          a = nan(size(sig1));a(sig1)=1;
%          plot(tb_,a*SigBars(1),'color',[0.5 0 0.5],'LineWidth',SigBarsWidth)
%          a = nan(size(sig2));a(sig2)=1;
%          plot(tb_,a*SigBars(1),'color',[0 0 1],'LineWidth',SigBarsWidth)
    
     end
     
     % Match Y scale across all delay lengths
     for iDelay = 1:length(Delays_)
%          subplot(1,length(Delays_),iDelay);
        subaxis(1,3,iDelay,'SpacingHoriz',0.01,'Margin',0.01);hold on
            ylims(iDelay,:) = get(gca,'YLim');
     end
%      ylims = [floor(min(ylims(:,1))),ceil(max(ylims(:,2)))];
     ylims = [-22,ceil(max(ylims(:,2)))];
     for iDelay = 1:length(Delays_)
%          subplot(1,length(Delays_),iDelay);
            subaxis(1,3,iDelay,'SpacingHoriz',0.01,'Margin',0.01);hold on
            set(gca,'YLim',ylims);
         axis off
     end
%      subplot(1,length(Delays_),iDelay);
     plot([19 20],[-13 -13],'-k','LineWidth',1.5)
     plot([20 20],[-13 -8],'-k','LineWidth',1.5)
     printFigs(1, pwd, '-djpeg', sprintf('%s_%s_unit%d_withErorrs',fname,Areas{iArea},iUnit))
 end
 
 
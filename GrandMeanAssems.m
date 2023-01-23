%% %%%%%% PREAMBLE %%%%%%
clear
close all
Target = 'LONG';
bw=0.05;
clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
%     pat = '/Volumes/HDD2/DNMTP/raw/';
    pat = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw/';

else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);

fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat','Ireneuszscopolamine_Events.mat','Krzysztofscopolamine_Events.mat',...
            'Ireneuszphysostigmine_Events.mat','IreneuszAM251_Events.mat','KrzesimirAM251_Events.mat'}; %'ALL_events.mat'
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
% pat2 = sprintf('C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\%sTaskonly',Target);
pat2 = fullfile(pat,'KDE_binsTaskonly',[Target 'Taskonly']);
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

% Thresholds=[0.1 0.05 0.01 5e-3 1e-3];
ThreshDurS = 1 %
threshold_ = 3;  % FL
threshold__ = 3; %FSCs
AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time
switch AssemblyChoice
    case 1
        pat2 = [pat 'KDE_bins' filesep Target filesep];
    case 2
        %pat2 = [pat 'KDE_binsTaskonly' filesep 'LONGTaskonly' filesep];
        pat2 = fullfile(pat,'KDE_binsTaskonly',[Target 'Taskonly']);
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end
%% Load files
AssOutput=[];
for iFile =1:length(fileList)
    %% Get the single unit files
    fname=strtok(fileList(iFile).name,'_');
    
    for iArea = 1:2
        fprintf('Getting units: %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
        load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
        
        iFR_{iArea} = iFR;
        % iFR_ = zscore(iFR_);
    end
    fprintf('Getting units: %d/%d %s (Joint)...\n',iFile,length(fileList),fname)
    iFR_{3} = [iFR_{1},iFR_{2}];
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    load(sprintf('%s%s%s.mat',pat,filesep,fname),'HPcells','PFCcells');
    %% Get the assembly membership info
    fname=strtok(fileListAss(iFile).name,'_');
    fprintf('Loading run %d/%d %s ...\n',iFile,length(fileList),fname)
    % for assemblies describing whole task period (cue - reward)
    load(sprintf('%s%s%s_iFR50_BehavOnly_FSC.mat',pat2,filesep,fname));
    load(sprintf('%s%s%s_iFR50_BehavOnly_AssemRes2.mat',pat2,filesep,fname));
    usel_out{3} = [usel_out{1},usel_out{2}+length(PFCcells);];
    Ass.usel_out = usel_out;
    Ass.nu = nu; Ass.nu(3) = sum(Ass.nu(1:2));
    Ass.units = units;
    for iArea = 1:3
        Ass.LocalMembers{iArea}    = unique(cell2mat(Ass.units{iArea}));
    end
    
    Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3}<=Ass.nu(1)));
    Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3}>Ass.nu(1))-Ass.nu(1));
    for iArea = 1:2
        Ass.JointMembers{iArea} = setdiff(unique(cell2mat(Ass.units{iArea})),Ass.LocalMembers{iArea});
        Ass.NonMembers{iArea}   = setdiff(1:Ass.nu(iArea),[Ass.LocalMembers{iArea},Ass.JointMembers{iArea}]);
    end
    
     Ass.LocalMembers
    
    clear units usel_out nu
    %% Prune any factors with insignificant activation
%     nAss = cellfun(@(x) size(x,2), FSCsel);
%     for iArea=1:3
%         if nAss(iArea)>0
%             i = nAss(iArea);
%             thresh_H  = ciHsc(iArea,threshold__);
%             FSCthresh = FSCsel{iArea}>thresh_H;
% %             idx =sum(FSCthresh,1)<1;
%             idx =sum(FSCthresh,1)<(ThreshDurS*bw^-1);
%            
%             FL{iArea}{i}(:,idx)=[];
%             FSCsel{iArea}(:,idx)=[];
%             Ass.units{iArea}(idx)=[];
%         end
%     end
    
    %% Prune any factors with less than 2 units (or single area factors for dCA1-mPFC)
    nAss = cellfun(@(x) size(x,2), FSCsel);
    for iArea=1:3
        if nAss(iArea)>0
            i = nAss(iArea);
            FL_ = FL{iArea}{i};
            thresh_H = ciHld(iArea,threshold_);
            FLthresh= FL_>thresh_H;
            idx =sum(FLthresh,1)<2;
            
            if iArea==3 % remove assems with single area loading
                idx =idx | (sum(FLthresh(1:length(Ass.usel_out{1}),:))<1 | ...
                    sum(FLthresh((length(Ass.usel_out{1})+1):end,:))<1);
            end
            FL{iArea}{i}(:,idx)=[];
            FSCsel{iArea}(:,idx)=[];
            Ass.units{iArea}(idx)=[];
        end
    end 

    %% Prepare datastructure for plotting 
    clear Factor
    Factor.FSC = FSCsel;
    Factor.ciHsc = ciHsc;
    Factor.TmtxS = Tmtx;
    Factor.units = Ass.units;
    Factor.unitIDs = Ass.usel_out;
    Factor.ciLld = ciLld;
    Factor.ciHld = ciHld;

    Factor.FL = FL;
    %% Remove double counted tb entries
    clear FSC_ TmtxS_
        
    [TmtxS_,ia]=unique(Factor.TmtxS);
    for iArea=1:3
        if size(Factor.FSC{iArea},2)>0
            FSC_{iArea} = Factor.FSC{iArea}(ia,:); 
            
        else
            FSC_{iArea} = [];
        end
    end
    
        Factor.FSC = FSC_;
        Factor.TmtxS = TmtxS_;
    
%% Get event times 
LC = [t.Short.SamplePress_LeftCorrect, t.Short.ChoicePress_LeftCorrect;...
      t.Medium.SamplePress_LeftCorrect,t.Medium.ChoicePress_LeftCorrect;...
      t.Long.SamplePress_LeftCorrect,  t.Long.ChoicePress_LeftCorrect];

RC = [t.Short.SamplePress_RightCorrect,t.Short.ChoicePress_RightCorrect;...
      t.Medium.SamplePress_RightCorrect,t.Medium.ChoicePress_RightCorrect;...
      t.Long.SamplePress_RightCorrect,t.Long.ChoicePress_RightCorrect];

LE = [t.Short.SamplePress_LeftError', t.Short.ChoicePress_LeftError';...
      t.Medium.SamplePress_LeftError',t.Medium.ChoicePress_LeftError';...
      t.Long.SamplePress_LeftError',  t.Long.ChoicePress_LeftError'];

RE = [t.Short.SamplePress_RightError',t.Short.ChoicePress_RightError';...
      t.Medium.SamplePress_RightError',t.Medium.ChoicePress_RightError';...
      t.Long.SamplePress_RightError',t.Long.ChoicePress_RightError'];

% LC = [t.Medium.SamplePress_LeftCorrect,t.Medium.ChoicePress_LeftCorrect;...
%       t.Long.SamplePress_LeftCorrect,  t.Long.ChoicePress_LeftCorrect];
% 
% RC = [t.Medium.SamplePress_RightCorrect,t.Medium.ChoicePress_RightCorrect;...
%       t.Long.SamplePress_RightCorrect,t.Long.ChoicePress_RightCorrect];
% 
% LE = [t.Medium.SamplePress_RightError',t.Medium.ChoicePress_RightError';...
%       t.Long.SamplePress_LeftError',  t.Long.ChoicePress_LeftError'];
% 
% RE = [t.Medium.SamplePress_RightError',t.Medium.ChoicePress_RightError';...
%       t.Long.SamplePress_RightError',t.Long.ChoicePress_RightError'];
  
  LC(sum(isnan(LC),2)>0,:)=[];
  RC(sum(isnan(RC),2)>0,:)=[];
  LE(sum(isnan(LE),2)>0,:)=[];
  RE(sum(isnan(RE),2)>0,:)=[];
  
  Events.Left = LC;
  Events.Right = RC;
  Events.LeftErr = LE;
  Events.RightErr = RE;
%%

tLims = [-5 5];
bw = 0.05;%diff(Factor.TmtxS([1 2]));
for iArea = 1:3
    LTrials = [];
    RTrials = [];
    LTrialsE = [];
    RTrialsE = [];
    if size(Factor.FSC{iArea},1)
    
    for iTrial = 1:size(Events.Left,1)
        tSample = FindClosestIndex(Factor.TmtxS,Events.Left(iTrial,1)/1e6)+(tLims/bw);
        tChoice = FindClosestIndex(Factor.TmtxS,Events.Left(iTrial,2)/1e6)+(tLims/bw);
        idx = [tSample(:,1):tSample(:,2),tChoice(:,1):tChoice(:,2)];
        LTrials = cat(3,LTrials,Factor.FSC{iArea}(idx,:));
    end
    for iTrial = 1:size(Events.Right,1)
        tSample = FindClosestIndex(Factor.TmtxS,Events.Right(iTrial,1)/1e6)+(tLims/bw);
        tChoice = FindClosestIndex(Factor.TmtxS,Events.Right(iTrial,2)/1e6)+(tLims/bw);
        idx = [tSample(:,1):tSample(:,2),tChoice(:,1):tChoice(:,2)];
        RTrials = cat(3,RTrials,Factor.FSC{iArea}(idx,:));
    end
    for iTrial = 1:size(Events.LeftErr,1)
        tSample = FindClosestIndex(Factor.TmtxS,Events.LeftErr(iTrial,1)/1e6)+(tLims/bw);
        tChoice = FindClosestIndex(Factor.TmtxS,Events.LeftErr(iTrial,2)/1e6)+(tLims/bw);
        idx = [tSample(:,1):tSample(:,2),tChoice(:,1):tChoice(:,2)];
        LTrialsE = cat(3,LTrialsE,Factor.FSC{iArea}(idx,:));
    end
    for iTrial = 1:size(Events.RightErr,1)
        tSample = FindClosestIndex(Factor.TmtxS,Events.RightErr(iTrial,1)/1e6)+(tLims/bw);
        tChoice = FindClosestIndex(Factor.TmtxS,Events.RightErr(iTrial,2)/1e6)+(tLims/bw);
        idx = [tSample(:,1):tSample(:,2),tChoice(:,1):tChoice(:,2)];
        RTrialsE = cat(3,RTrialsE,Factor.FSC{iArea}(idx,:));
    end
    
    tb=(1:size(RTrials,1))*bw;
    
    LTrials = nanmean(LTrials,3);
    RTrials = nanmean(RTrials,3);
    LTrialsE = nanmean(LTrialsE,3);
    RTrialsE = nanmean(RTrialsE,3);
    
    ciHsc(iArea,threshold_)
    
    for iAss =1:size(LTrials,2)
        if  sum(LTrials(:,iAss))>sum(RTrials(:,iAss));
%                if sum(LTrials(1:max(find((tb<15))),iAss))>sum(RTrials(max(find((tb<15))),iAss));
            AssOutput.Pref{iArea}{iFile}(:,iAss) = LTrials(:,iAss);
            AssOutput.PrefE{iArea}{iFile}(:,iAss) = LTrialsE(:,iAss);
            AssOutput.NonPref{iArea}{iFile}(:,iAss)  = RTrials(:,iAss);
            AssOutput.NonPrefE{iArea}{iFile}(:,iAss) = RTrialsE(:,iAss);
            
            AssOutput.PrefSig{iArea}{iFile}(:,iAss) = LTrials(:,iAss)>ciHsc(iArea,threshold__);
            AssOutput.PrefESig{iArea}{iFile}(:,iAss) = LTrialsE(:,iAss)>ciHsc(iArea,threshold__);
            AssOutput.NonPrefSig{iArea}{iFile}(:,iAss)  = RTrials(:,iAss)>ciHsc(iArea,threshold__);
            AssOutput.NonPrefESig{iArea}{iFile}(:,iAss) = RTrialsE(:,iAss)>ciHsc(iArea,threshold__);
            
        else
            AssOutput.Pref{iArea}{iFile}(:,iAss) = RTrials(:,iAss);
            AssOutput.PrefE{iArea}{iFile}(:,iAss) = RTrialsE(:,iAss);
            AssOutput.NonPref{iArea}{iFile}(:,iAss) = LTrials(:,iAss);
            AssOutput.NonPrefE{iArea}{iFile}(:,iAss) = LTrialsE(:,iAss);
            
            
            AssOutput.PrefSig{iArea}{iFile}(:,iAss) = RTrials(:,iAss)>ciHsc(iArea,threshold__);
            AssOutput.PrefESig{iArea}{iFile}(:,iAss) = RTrialsE(:,iAss)>ciHsc(iArea,threshold__);
            AssOutput.NonPrefSig{iArea}{iFile}(:,iAss)  = LTrials(:,iAss)>ciHsc(iArea,threshold__);
            AssOutput.NonPrefESig{iArea}{iFile}(:,iAss) = LTrialsE(:,iAss)>ciHsc(iArea,threshold__);
        end
        
        
    end
    end
end
end
%% Plot Assembly activation - average over all directions
figure('color','w');hold on
for iArea =1:3

%     P = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     P = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),P,'UniformOutput',false));
%     nP = AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea}));
%     nP = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),nP,'UniformOutput',false));
%     C = P';
    P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
    nP = cell2mat(AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea})));
%    P=zscore(P);nP=zscore(nP);
    P = (P+nP)./2;
    C = P';
%    C=zscore(P)';
   
   PSig = cell2mat(AssOutput.PrefSig{iArea}(~isempty_cell(AssOutput.PrefSig{iArea})));
   PSig= sum(PSig,2)./size(PSig,2)*100; 
   
   
   ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb',color_{iArea},0.8)
   
%    imagesc(repmat(tb,1,2),iArea+[-0.125*ones(size(PSig,1),1);0.25*ones(size(PSig,1),1)], 1*[PSig,PSig]')
%     rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
%     plot(tb,LTrials,'b');plot(tb,LTrialsE,'r')
end

% for iArea =1:3
%    P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
%    C=zscore(P)';
%    plot(tb,nanmean(C),color_{iArea},'LineWidth',3.5)
% end
colormap(flipud(gray))
caxis([0 100])
axis off
%% Plot Assembly activation - Preferred
figure('color','w');hold on
for iArea =1:3

%     P = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     P = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),P,'UniformOutput',false));
%     nP = AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea}));
%     nP = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),nP,'UniformOutput',false));
%     C = P';
   P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
   nP = cell2mat(AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea})));
%    P=zscore(P);nP=zscore(nP);
%     P = (P+nP)./2;
    C = P';
%    C=zscore(P)';
   
   PSig = cell2mat(AssOutput.PrefSig{iArea}(~isempty_cell(AssOutput.PrefSig{iArea})));
   PSig= sum(PSig,2)./size(PSig,2)*100; 
   
   
   ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb',color_{iArea},0.8)
   
%    imagesc(repmat(tb,1,2),iArea+[-0.125*ones(size(PSig,1),1);0.25*ones(size(PSig,1),1)], 1*[PSig,PSig]')
%     rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
%     plot(tb,LTrials,'b');plot(tb,LTrialsE,'r')
end

% for iArea =1:3
%    P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
%    C=zscore(P)';
%    plot(tb,nanmean(C),color_{iArea},'LineWidth',3.5)
% end
colormap(flipud(gray))
caxis([0 100])
axis off
%% Plot Assembly activation - Preferred 2
figure('color','w');hold on
for iArea =1:3

%     P = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     P = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),P,'UniformOutput',false));
%     nP = AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea}));
%     nP = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),nP,'UniformOutput',false));
%     C = P';
   P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
   nP = cell2mat(AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea})));
%    P=zscore(P);nP=zscore(nP);
%    P = (P+nP)./2;
    
%    C=zscore(P)';
   C=(P)';
   
   PSig = cell2mat(AssOutput.PrefSig{iArea}(~isempty_cell(AssOutput.PrefSig{iArea})));
   PSig = sum(PSig,2)./size(PSig,2)*100; 
   col__ = color_{iArea};
   if iArea<3
    col__ = col_{2};
   else
    col__ = col_{3};
   end
   plot(tb,iArea*2.5-1+ones(size(tb)),'--','LineWidth',1.5,'color',0.6*[1 1 1])   
   if iArea ==1
       
       SampleAss = mean(C(:,tb>0 & tb<=10),2)>mean(C(:,tb>10 & tb<=15),2);
       ChoiceAss = mean(C(:,tb>0 & tb<=10),2)<mean(C(:,tb>10 & tb<=15),2);
       sum([SampleAss,ChoiceAss])./length(SampleAss);
       Sample = C(SampleAss,:);
       Choice = C(ChoiceAss,:);
       
       ciplot(iArea*2.5+nanmean(Sample)+nansem(Sample),iArea*2.5+nanmean(Sample)-nansem(Sample),tb',col__,0.8)
       plot(tb,iArea*2.5+nanmean(Sample)+nansem(Sample),'-','LineWidth',2.5,'color',0*[1 1 1])
       plot(tb,iArea*2.5+nanmean(Sample)-nansem(Sample),'-','LineWidth',2.5,'color',0*[1 1 1])
       
       ciplot(iArea*2.5+nanmean(Choice)+nansem(Choice),iArea*2.5+nanmean(Choice)-nansem(Choice),tb',col__,0.8)
       plot(tb,iArea*2.5+nanmean(Choice)+nansem(Choice),'-','LineWidth',3,'color',0*[1 1 1])
       plot(tb,iArea*2.5+nanmean(Choice)-nansem(Choice),'-','LineWidth',3,'color',0*[1 1 1])
   else
       ciplot(iArea*2.5 + nanmean(C)+nansem(C),iArea*2.5 + nanmean(C)-nansem(C),tb',col__,1)
       plot(tb,iArea*2.5+nanmean(C)+nansem(C),'-','LineWidth',2.5,'color',0*[1 1 1])   
       plot(tb,iArea*2.5+nanmean(C)-nansem(C),'-','LineWidth',2.5,'color',0*[1 1 1]) 
%    ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb',col__,1)
   end
   
   
   
%    imagesc(repmat(tb,1,2),iArea*3+[-0.125*ones(size(PSig,1),1);0.25*ones(size(PSig,1),1)], 1*[PSig,PSig]')
%     rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
%     plot(tb,LTrials,'b');plot(tb,LTrialsE,'r')
end

% for iArea =1:3
%    P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
%    C=zscore(P)';
%    plot(tb,nanmean(C),color_{iArea},'LineWidth',3.5)
% end
% colormap(flipud(gray))
% colormap(jet)
caxis([0 50])
axis off
%% Plot Assembly activation - Preferred vs non-preferred 
figure('color','w');hold on
for iArea =1:3

%     P = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     P = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),P,'UniformOutput',false));
%     nP = AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea}));
%     nP = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),nP,'UniformOutput',false));
%     C = P';
   P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})))';
   nP = cell2mat(AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea})))';
%    P=zscore(P);nP=zscore(nP);
   
   PSig = cell2mat(AssOutput.PrefSig{iArea}(~isempty_cell(AssOutput.PrefSig{iArea})));
   PSig = sum(PSig,2)./size(PSig,2)*100; 
   col__ = color_{iArea};

   plot(tb,iArea*2.5-1+ones(size(tb)),'--','LineWidth',1.5,'color',0.6*[1 1 1])   
   ciplot(iArea*2.5 + nanmean(P)+nansem(P),iArea*2.5 + nanmean(P)-nansem(P),tb',col__,1)
   plot(tb,iArea*2.5+nanmean(P)+nansem(P),'-','LineWidth',2.5,'color',0*[1 1 1])   
   plot(tb,iArea*2.5+nanmean(P)-nansem(P),'-','LineWidth',2.5,'color',0*[1 1 1]) 

   ciplot(iArea*2.5 + nanmean(nP)+nansem(nP),iArea*2.5 + nanmean(nP)-nansem(nP),tb',col__,0.6)
   plot(tb,iArea*2.5+nanmean(nP)+nansem(nP),':','LineWidth',2.5,'color',0*[1 1 1])   
   plot(tb,iArea*2.5+nanmean(nP)-nansem(nP),':','LineWidth',2.5,'color',0*[1 1 1]) 
   
   %    ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb',col__,1)
%    imagesc(repmat(tb,1,2),iArea*3+[-0.125*ones(size(PSig,1),1);0.25*ones(size(PSig,1),1)], 1*[PSig,PSig]')
%     rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
%     plot(tb,LTrials,'b');plot(tb,LTrialsE,'r')
end

% for iArea =1:3
%    P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
%    C=zscore(P)';
%    plot(tb,nanmean(C),color_{iArea},'LineWidth',3.5)
% end
% colormap(flipud(gray))
% colormap(jet)
caxis([0 50])
axis off
%% Plot Assembly activation - C vs E
figure
for iArea =1:3
    subplot(1,3,iArea);hold on
    C = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
    E = cell2mat(AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
    C = zscore(C,1);
    E = zscore(E,1);
%     C = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     E = AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     C = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),C,'UniformOutput',false));
%     E = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),E,'UniformOutput',false));
%     
   temp = [C;E];
%     temp = zscore(temp);
   C = temp(1:length(tb),:)';
   E = temp(length(tb)+1:end,:)';
%    ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb','k',0.8)
%    ciplot(nanmean(E)+nansem(E),nanmean(E)-nansem(E),tb','r',0.8)
% ciplot(nanmean((E-C))+nansem((E-C)),nanmean((E-C))-nansem((E-C)),tb','g',1)

   ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb',color_{iArea},0.8)
   ciplot(nanmean(E)+nansem(E),nanmean(E)-nansem(E),tb',color_{iArea},0.3)
   
 sig1 = permtest2vec(smooth2a(C,0,20)',smooth2a(E,0,20)',50,0.05);
                    a = nan(size(sig1));a(sig1)=1;
%     plot(f,maxY*a-0.1,'-k','LineWidth',2.5)
    a(isnan(a))=0;
    [~,Xpk,Wpk] = findpeaks(a);
    plot([min(tb),max(tb)],1.2*[1 1],'color',0.6*[1 1 1],'LineWidth',4) 
    for iPeak=1:length(Xpk)
       plot(tb([Xpk(iPeak) Xpk(iPeak)+Wpk(iPeak)]),ones(1,2)*1.2,'-k','LineWidth',4) 
    end
%     plot([5 5],[-1 .6],'g','LineWidth',1.5)
%     plot([15 15],[-1 0.6],'r','LineWidth',1.5)
axis([0 20 -1.5 2.2])
axis off
end
%% Plot Assembly activation - C vs E 2
figure; hold on
for iArea =1:3
    
    C = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
    E = cell2mat(AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
%     C = zscore(C,1);
%     E = zscore(E,1);
%     C = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     E = AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     C = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),C,'UniformOutput',false));
%     E = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),E,'UniformOutput',false));
%     
   temp = [C;E];
%     temp = zscore(temp);
   C = temp(1:length(tb),:)';
   E = temp(length(tb)+1:end,:)';
%    ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb','k',0.8)
%    ciplot(nanmean(E)+nansem(E),nanmean(E)-nansem(E),tb','r',0.8)
% ciplot(nanmean((E-C))+nansem((E-C)),nanmean((E-C))-nansem((E-C)),tb','g',1)
   col__ = color_{iArea};
   if iArea<3
    col__ = col_{2};
   else
    col__ = col_{3};
   end
   plot(tb,iArea*2.5-1+ones(size(tb)),'--','LineWidth',2.5,'color',0.6*[1 1 1])   
   
   
   
   
   ciplot(iArea*2.5+nanmean(C)+nansem(C),iArea*2.5+nanmean(C)-nansem(C),tb',col__,0.6)
   plot(tb,iArea*2.5+nanmean(C)+nansem(C),'-','LineWidth',2.5,'color',col__)   
   plot(tb,iArea*2.5+nanmean(C)-nansem(C),'-','LineWidth',2.5,'color',col__)   
   
   
   ciplot(iArea*2.5+nanmean(E)+nansem(E),iArea*2.5+nanmean(E)-nansem(E),tb','r',1)
   plot(tb,iArea*2.5+nanmean(E)+nansem(E),'-','LineWidth',2.5,'color','k')   
   plot(tb,iArea*2.5+nanmean(E)-nansem(E),'-','LineWidth',2.5,'color','k')   
   
   
   sig1 = permtest2vec(smooth2a(C,0,20)',smooth2a(E,0,20)',50,0.05);
                    a = nan(size(sig1));a(sig1)=1;
%     plot(f,maxY*a-0.1,'-k','LineWidth',2.5)
    a(isnan(a))=0;
    [~,Xpk,Wpk] = findpeaks(a);
%     plot([min(tb),max(tb)],1.2*[1 1],'color',0.6*[1 1 1],'LineWidth',4) 
%     for iPeak=1:length(Xpk)
%        plot(tb([Xpk(iPeak) Xpk(iPeak)+Wpk(iPeak)]),ones(1,2)*1.2,'-k','LineWidth',4) 
%     end
%     plot([5 5],[-1 .6],'g','LineWidth',1.5)
%     plot([15 15],[-1 0.6],'r','LineWidth',1.5)
% axis([0 20 -1.5 2.2])
axis off
end
%% Plot Assembly activation - C - E 
figure
for iArea =1:3
    subplot(3,1,iArea);hold on
    C = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
    E = cell2mat(AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
%     C = zscore(C,1);
%     E = zscore(E,1);
%     C = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     E = AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     C = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),C,'UniformOutput',false));
%     E = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),E,'UniformOutput',false));
%     
   temp = [C;E];
%     temp = zscore(temp);
   C = temp(1:length(tb),:)';
   E = temp(length(tb)+1:end,:)';
%    ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb','k',0.8)
%    ciplot(nanmean(E)+nansem(E),nanmean(E)-nansem(E),tb','r',0.8)
% ciplot(nanmean((E-C))+nansem((E-C)),nanmean((E-C))-nansem((E-C)),tb','g',1)
   col__ = color_{iArea};
   if iArea<3
    col__ = col_{2};
   else
    col__ = col_{3};
   end
   plot(tb,iArea*2.5-1+ones(size(tb)),'--','LineWidth',2.5,'color',0.6*[1 1 1])   
   
   
   
   
   ciplot(iArea*2.5+nanmean(E-C)+nansem(C-E),iArea*2.5+nanmean(E-C)-nansem(C-E),tb',col__,0.6)
   plot(tb,iArea*2.5+nanmean(E-C)+nansem(E-C),'-','LineWidth',2.5,'color',col__)   
   plot(tb,iArea*2.5+nanmean(E-C)-nansem(E-C),'-','LineWidth',2.5,'color',col__)   
   
   
   sig1 = permtest2vec(smooth2a(C,0,20)',smooth2a(E,0,20)',50,0.05);
                    a = nan(size(sig1));a(sig1)=1;
%     plot(f,maxY*a-0.1,'-k','LineWidth',2.5)
    a(isnan(a))=0;
    [~,Xpk,Wpk] = findpeaks(a);
%     plot([min(tb),max(tb)],1.2*[1 1],'color',0.6*[1 1 1],'LineWidth',4) 
%     for iPeak=1:length(Xpk)
%        plot(tb([Xpk(iPeak) Xpk(iPeak)+Wpk(iPeak)]),ones(1,2)*1.2,'-k','LineWidth',4) 
%     end
%     plot([5 5],[-1 .6],'g','LineWidth',1.5)
%     plot([15 15],[-1 0.6],'r','LineWidth',1.5)
% axis([0 20 -1.5 2.2])
axis off
end

%% Plot Assembly activation - Preferred C vs E, sample 
figure
for iArea =1:3
    subplot(1,3,iArea);hold on
    C = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})))';
    E = cell2mat(AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea})))';
    SampleAss = mean(C(:,tb>0 & tb<=10),2)>mean(C(:,tb>10 & tb<=15),2);
    ChoiceAss = mean(C(:,tb>0 & tb<=10),2)<mean(C(:,tb>10 & tb<=15),2);
    sum([SampleAss,ChoiceAss])./length(SampleAss);
    C = C(SampleAss,:);E = E(SampleAss,:);
%     C = C(ChoiceAss,:);E = E(ChoiceAss,:);
%     C = zscore(C,1);
%     E = zscore(E,1);
%     C = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     E = AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     C = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),C,'UniformOutput',false));
%     E = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),E,'UniformOutput',false));
%     
%    temp = [C;E];
% %     temp = zscore(temp);
%    C = temp(1:length(tb),:);
%    E = temp(length(tb)+1:end,:);
%    ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb','k',0.8)
%    ciplot(nanmean(E)+nansem(E),nanmean(E)-nansem(E),tb','r',0.8)
% ciplot(nanmean((E-C))+nansem((E-C)),nanmean((E-C))-nansem((E-C)),tb','g',1)

   ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb',color_{iArea},0.8)
   ciplot(nanmean(E)+nansem(E),nanmean(E)-nansem(E),tb',color_{iArea},0.3)
   
 sig1 = permtest2vec(smooth2a(C,0,20)',smooth2a(E,0,20)',100,0.05);
                    a = nan(size(sig1));a(sig1)=1;
%     plot(f,maxY*a-0.1,'-k','LineWidth',2.5)
    a(isnan(a))=0;
    [~,Xpk,Wpk] = findpeaks(a);
    plot([min(tb),max(tb)],2.5*[1 1],'color',0.6*[1 1 1],'LineWidth',4) 
    for iPeak=1:length(Xpk)
       plot(tb([Xpk(iPeak) Xpk(iPeak)+Wpk(iPeak)]),ones(1,2)*2.5,'-k','LineWidth',4) 
    end
    plot([5 5],[-1 .6],'g','LineWidth',1.5)
    plot([15 15],[-1 0.6],'r','LineWidth',1.5)
    plot([10 10],[-1 0.6],'w','LineWidth',3)
    plot([10 10],[-1 0.6],':k','LineWidth',1.5)
axis([0 20 -1.5 2.7])
axis off
end
%% Plot Assembly activation - Preferred C vs E, choice
figure
for iArea =1:3
    subplot(1,3,iArea);hold on
    C = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})))';
    E = cell2mat(AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea})))';
    SampleAss = mean(C(:,tb>0 & tb<=10),2)>mean(C(:,tb>10 & tb<=15),2);
    ChoiceAss = mean(C(:,tb>0 & tb<=10),2)<mean(C(:,tb>10 & tb<=15),2);
    sum([SampleAss,ChoiceAss])./length(SampleAss);
%     C = C(SampleAss,:);E = E(SampleAss,:);
    C = C(ChoiceAss,:);E = E(ChoiceAss,:);
%     C = zscore(C,1);
%     E = zscore(E,1);
%     C = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     E = AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     C = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),C,'UniformOutput',false));
%     E = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),E,'UniformOutput',false));
%     
%    temp = [C;E];
% %     temp = zscore(temp);
%    C = temp(1:length(tb),:);
%    E = temp(length(tb)+1:end,:);
%    ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb','k',0.8)
%    ciplot(nanmean(E)+nansem(E),nanmean(E)-nansem(E),tb','r',0.8)
% ciplot(nanmean((E-C))+nansem((E-C)),nanmean((E-C))-nansem((E-C)),tb','g',1)

   ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb',color_{iArea},0.8)
   ciplot(nanmean(E)+nansem(E),nanmean(E)-nansem(E),tb',color_{iArea},0.3)
   
 sig1 = permtest2vec(smooth2a(C,0,20)',smooth2a(E,0,20)',100,0.05);
                    a = nan(size(sig1));a(sig1)=1;
%     plot(f,maxY*a-0.1,'-k','LineWidth',2.5)
    a(isnan(a))=0;
    [~,Xpk,Wpk] = findpeaks(a);
    plot([min(tb),max(tb)],2.5*[1 1],'color',0.6*[1 1 1],'LineWidth',4) 
    for iPeak=1:length(Xpk)
       plot(tb([Xpk(iPeak) Xpk(iPeak)+Wpk(iPeak)]),ones(1,2)*2.5,'-k','LineWidth',4) 
    end
    plot([5 5],[-1 .6],'g','LineWidth',1.5)
    plot([15 15],[-1 0.6],'r','LineWidth',1.5)
    axis([0 20 -1.5 2.7])
    axis off
end
%% Plot Assembly activation - , sort by sample or choice
figure('color','w');hold on
for iArea =1:3

%     P = AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
%     P = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),P,'UniformOutput',false));
%     nP = AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea}));
%     nP = cell2mat(cellfun(@(x) nanmean(zscore(x,1),2),nP,'UniformOutput',false));
%     C = P';
   P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})))';
   nP = cell2mat(AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea})))';
%    P=zscore(P);nP=zscore(nP);
   
   PSig = cell2mat(AssOutput.PrefSig{iArea}(~isempty_cell(AssOutput.PrefSig{iArea})));
   PSig = sum(PSig,2)./size(PSig,2)*100; 
   col__ = color_{iArea};

   plot(tb,iArea*2.5-1+ones(size(tb)),'--','LineWidth',1.5,'color',0.6*[1 1 1])   
   ciplot(iArea*2.5 + nanmean(P)+nansem(P),iArea*2.5 + nanmean(P)-nansem(P),tb',col__,1)
   plot(tb,iArea*2.5+nanmean(P)+nansem(P),'-','LineWidth',2.5,'color',0*[1 1 1])   
   plot(tb,iArea*2.5+nanmean(P)-nansem(P),'-','LineWidth',2.5,'color',0*[1 1 1]) 

   ciplot(iArea*2.5 + nanmean(nP)+nansem(nP),iArea*2.5 + nanmean(nP)-nansem(nP),tb',col__,0.6)
   plot(tb,iArea*2.5+nanmean(nP)+nansem(nP),':','LineWidth',2.5,'color',0*[1 1 1])   
   plot(tb,iArea*2.5+nanmean(nP)-nansem(nP),':','LineWidth',2.5,'color',0*[1 1 1]) 
   
   %    ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb',col__,1)
%    imagesc(repmat(tb,1,2),iArea*3+[-0.125*ones(size(PSig,1),1);0.25*ones(size(PSig,1),1)], 1*[PSig,PSig]')
%     rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
%     plot(tb,LTrials,'b');plot(tb,LTrialsE,'r')
end

% for iArea =1:3
%    P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
%    C=zscore(P)';
%    plot(tb,nanmean(C),color_{iArea},'LineWidth',3.5)
% end
% colormap(flipud(gray))
% colormap(jet)
caxis([0 50])
axis off
%% Plot PFC, sample and choice active assemblies

figure;hold on
iArea =1;
A_ = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})))';

SampleAss = mean(A_(:,tb>0 & tb<=10),2)>mean(A_(:,tb>10 & tb<=15),2);
ChoiceAss = mean(A_(:,tb>0 & tb<=10),2)<mean(A_(:,tb>10 & tb<=15),2);
sum([SampleAss,ChoiceAss])./length(SampleAss);
S = A_(SampleAss,:);
C = A_(ChoiceAss,:);

ciplot(nanmean(S)+nansem(S),nanmean(S)-nansem(S),tb',color_{iArea},0.8)
plot(tb,nanmean(S)+nansem(S),'-','LineWidth',2.5,'color',0*[1 1 1])   
plot(tb,nanmean(S)-nansem(S),'-','LineWidth',2.5,'color',0*[1 1 1]) 

ciplot(nanmean(C)+nansem(C),nanmean(C)-nansem(C),tb',color_{iArea},0.8)
plot(tb,nanmean(C)+nansem(C),'-','LineWidth',3,'color',0*[1 1 1])   
plot(tb,nanmean(C)-nansem(C),'-','LineWidth',3,'color',0*[1 1 1]) 
   

plot([5 5],[-1 .6],'g','LineWidth',1.5)
plot([15 15],[-1 0.6],'r','LineWidth',1.5)
plot([10 10],[-1 0.6],'w','LineWidth',3)
plot([10 10],[-1 0.6],':k','LineWidth',1.5)
axis([0 20 -1.5 2.7])
axis off
%% Plot Assembly activation - Preferred 2 ****
figure('color','w')
offset = 3;
for iArea =1:3

   P = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
   
   C=(P)';
   
   
   col__ = color_{iArea};
   if iArea<3
    col__ = col_{2};
   else
    col__ = col_{3};
   end
   
   subplot(1,2,1); hold on 
   plot(tb,iArea*offset*ones(size(tb)),'--','LineWidth',1.5,'color',0.6*[1 1 1])
   SampleAss = mean(C(:,tb>0 & tb<=10),2)>mean(C(:,tb>10 & tb<=15),2);
   Sample = C(SampleAss,:);
   ciplot(iArea*offset+nanmean(Sample)+nansem(Sample),iArea*offset+nanmean(Sample)-nansem(Sample),tb',col__,0.8)
   plot(tb,iArea*offset+nanmean(Sample)+nansem(Sample),'-','LineWidth',2.5,'color',0*[1 1 1])
   plot(tb,iArea*offset+nanmean(Sample)-nansem(Sample),'-','LineWidth',2.5,'color',0*[1 1 1])
   
   subplot(1,2,2); hold on 
   plot(tb,iArea*offset*ones(size(tb)),'--','LineWidth',1.5,'color',0.6*[1 1 1])
   ChoiceAss = mean(C(:,tb>0 & tb<=10),2)<mean(C(:,tb>10 & tb<=15),2);
   Choice = C(ChoiceAss,:);
   ciplot(iArea*offset+nanmean(Choice)+nansem(Choice),iArea*offset+nanmean(Choice)-nansem(Choice),tb',col__,0.8)
   plot(tb,iArea*offset+nanmean(Choice)+nansem(Choice),'-','LineWidth',3,'color',0*[1 1 1])
   plot(tb,iArea*offset+nanmean(Choice)-nansem(Choice),'-','LineWidth',3,'color',0*[1 1 1])
   
 
   
   
   
%    imagesc(repmat(tb,1,2),iArea*3+[-0.125*ones(size(PSig,1),1);0.25*ones(size(PSig,1),1)], 1*[PSig,PSig]')
%     rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
%     plot(tb,LTrials,'b');plot(tb,LTrialsE,'r')
end

   subplot(1,2,1); axis off
   subplot(1,2,2); axis off
   
%    imagesc(repmat(tb,1,2),iArea*3+[-0.125*ones(size(PSig,1),1);0.25*ones(size(PSig,1),1)], 1*[PSig,PSig]')
%     rectangle('Position',[0 -5.5 20 3],'LineWidth',1,'EdgeColor','r')
%     plot(tb,LTrials,'b');plot(tb,LTrialsE,'r')

%% Plot Assembly activation - C Vs E ****
figure('color','w')
offset = 4;
smooth_ = 50;
subplot(1,2,1); hold on
plot([5 5],[2.5 14],'linewidth',2,'color',[0 1 0 0.3])
plot([15 15],[2.5 14],'linewidth',2,'color',[1 0 0 0.3])
axis off
subplot(1,2,2); hold on
plot([5 5],[2.5 14],'linewidth',2,'color',[0 1 0 0.3])
plot([15 15],[2.5 14],'linewidth',2,'color',[1 0 0 0.3])
axis off
for iArea =1:3
    col__ = color_{iArea};
    if iArea<3
        col__ = col_{2};
    else
        col__ = col_{3};
    end
    
    C = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})))';
    E = cell2mat(AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea})))';
    
    temp = [C,E];
%     temp = zscore(temp')';
   C = temp(:,1:length(tb));
   E = temp(:,length(tb)+1:end);
   
   
    subplot(1,2,1); hold on
    SampleAss = mean(C(:,tb>0 & tb<=10),2)>mean(C(:,tb>10 & tb<=15),2);
    Sample = C(SampleAss,:);
    SampleErr = E(SampleAss,:);
    
    plot(tb,iArea*offset*ones(size(tb)),'--','LineWidth',1.5,'color',0.6*[1 1 1])
    
    ciplot(iArea*offset+nanmean(Sample)+nansem(Sample),iArea*offset+nanmean(Sample)-nansem(Sample),tb',col__,0.8)
    plot(tb,iArea*offset+nanmean(Sample)+nansem(Sample),'-','LineWidth',2,'color',0*[1 1 1])
    plot(tb,iArea*offset+nanmean(Sample)-nansem(Sample),'-','LineWidth',2,'color',0*[1 1 1])
    
    ciplot(iArea*offset+nanmean(SampleErr)+nansem(SampleErr),iArea*offset+nanmean(SampleErr)-nansem(SampleErr),tb','r',0.8)
    plot(tb,iArea*offset+nanmean(SampleErr)+nansem(SampleErr),'-','LineWidth',2,'color',[0 0 0 1])
    plot(tb,iArea*offset+nanmean(SampleErr)-nansem(SampleErr),'-','LineWidth',2,'color',[0 0 0 1])
    
    sig1 = permtest2vec(smooth2a(Sample,0,smooth_)',smooth2a(SampleErr,0,smooth_)',100,0.05);
    a = nan(size(sig1));a(sig1)=1;
    a(isnan(a))=0;
    [~,Xpk,Wpk] = findpeaks(a);
    for iPeak=1:length(Xpk)
        plot(tb([Xpk(iPeak) Xpk(iPeak)+Wpk(iPeak)]),iArea*offset+ones(1,2),'-k','LineWidth',4)
    end
    
    subplot(1,2,2); hold on
    ChoiceAss = mean(C(:,tb>0 & tb<=10),2)<mean(C(:,tb>10 & tb<=15),2);
    Choice    = C(ChoiceAss,:);
    ChoiceErr = E(ChoiceAss,:);
    
    plot(tb,iArea*offset*ones(size(tb)),'--','LineWidth',1.5,'color',0.6*[1 1 1])
    
    ciplot(iArea*offset+nanmean(Choice)+nansem(Choice),iArea*offset+nanmean(Choice)-nansem(Choice),tb',col__,0.8)
    plot(tb,iArea*offset+nanmean(Choice)+nansem(Choice),'-','LineWidth',2,'color',0*[1 1 1])
    plot(tb,iArea*offset+nanmean(Choice)-nansem(Choice),'-','LineWidth',2,'color',0*[1 1 1])
    
    ciplot(iArea*offset+nanmean(ChoiceErr)+nansem(ChoiceErr),iArea*offset+nanmean(ChoiceErr)-nansem(ChoiceErr),tb','r',0.8)
    plot(tb,iArea*offset+nanmean(ChoiceErr)+nansem(ChoiceErr),'-','LineWidth',2,'color',[0 0 0 1])
    plot(tb,iArea*offset+nanmean(ChoiceErr)-nansem(ChoiceErr),'-','LineWidth',2,'color',[0 0 0 1])
    
    sig1 = permtest2vec(smooth2a(Choice,0,smooth_)',smooth2a(ChoiceErr,0,smooth_)',100,0.05);
    a = nan(size(sig1));a(sig1)=1;
    a(isnan(a))=0;
    [~,Xpk,Wpk] = findpeaks(a);
    for iPeak=1:length(Xpk)
        plot(tb([Xpk(iPeak) Xpk(iPeak)+Wpk(iPeak)]),iArea*offset+ones(1,2),'-k','LineWidth',4)
    end
    
end

subplot(1,2,1); hold on
plot([10 10],[2.5 14],'linewidth',4,'color','w')
axis off
subplot(1,2,2); hold on
plot([10 10],[2.5 14],'linewidth',4,'color','w')

%% Plot Assembly activation - C vs E 2
figure; hold on
for iArea =1:3
    
    C = cell2mat(AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
    E = cell2mat(AssOutput.PrefE{iArea}(~isempty_cell(AssOutput.Pref{iArea})));
  
   temp = [C;E];
   C = temp(1:length(tb),:)';
   E = temp(length(tb)+1:end,:)';

   if iArea==1
       SampleAss = mean(C(:,tb>0 & tb<=10),2)>mean(C(:,tb>10 & tb<=15),2);
       ChoiceAss = mean(C(:,tb>0 & tb<=10),2)<mean(C(:,tb>10 & tb<=15),2);
       sum([SampleAss,ChoiceAss])./length(SampleAss);
          C = C(SampleAss,:);E = E(SampleAss,:);
%        C = C(ChoiceAss,:);E = E(ChoiceAss,:);
   end
   col__ = color_{iArea};
   if iArea<3
    col__ = col_{2};
   else
    col__ = col_{3};
   end
   plot(tb,iArea*2.5-1+ones(size(tb)),'--','LineWidth',2.5,'color',0.6*[1 1 1])   
   
   
   
   
   ciplot(iArea*2.5+nanmean(C)+nansem(C),iArea*2.5+nanmean(C)-nansem(C),tb',col__,0.6)
   plot(tb,iArea*2.5+nanmean(C)+nansem(C),'-','LineWidth',2.5,'color',col__)   
   plot(tb,iArea*2.5+nanmean(C)-nansem(C),'-','LineWidth',2.5,'color',col__)   
   
   
   ciplot(iArea*2.5+nanmean(E)+nansem(E),iArea*2.5+nanmean(E)-nansem(E),tb','r',1)
   plot(tb,iArea*2.5+nanmean(E)+nansem(E),'-','LineWidth',2.5,'color','k')   
   plot(tb,iArea*2.5+nanmean(E)-nansem(E),'-','LineWidth',2.5,'color','k')   
   
   
   sig1 = permtest2vec(smooth2a(C,0,20)',smooth2a(E,0,20)',50,0.05);
                    a = nan(size(sig1));a(sig1)=1;
%     plot(f,maxY*a-0.1,'-k','LineWidth',2.5)
    a(isnan(a))=0;
    [~,Xpk,Wpk] = findpeaks(a);
%     plot([min(tb),max(tb)],iArea*2.5+1.2*[1 1],'color',0.6*[1 1 1],'LineWidth',4) 
    for iPeak=1:length(Xpk)
       plot(tb([Xpk(iPeak) Xpk(iPeak)+Wpk(iPeak)]),iArea*2.5+ones(1,2)*1.2,'-k','LineWidth',4) 
    end
%     plot([5 5],[-1 .6],'g','LineWidth',1.5)
%     plot([15 15],[-1 0.6],'r','LineWidth',1.5)
% axis([0 20 -1.5 2.2])
axis off
end
%% Fractions of Sample and Choice mPFC Assem
iArea =1;

A_= AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
FracSample=[];
for i=1:length(A_)
    C = A_{i}';
    SampleAss = mean(C(:,tb>0 & tb<=10),2)>mean(C(:,tb>10 & tb<=15),2);
    ChoiceAss = mean(C(:,tb>0 & tb<=10),2)<mean(C(:,tb>10 & tb<=15),2);
    FracSample(i) = (sum(SampleAss)./sum(ChoiceAss))./length(SampleAss);
end
FracSample(isinf(FracSample))=1;
figure; hold on
b = bar(1,nanmean(FracSample));
b.FaceColor = col_{2};
b.EdgeColor = 'k';
b.LineWidth = 1.5;
errorbar(1,nanmean(FracSample),nansem(FracSample),'k','LineWidth',1.5)
set(gca,'xtick',[])
axis([0.5 1.5 0 1])
view([90,-90])

%% Fractions of Sample and Choice assems - discrete ****
figure; hold on
for iArea =1:3
    
    A_= AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
    A_nP= AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea}));
    FracSample=[];
    for i=1:length(A_)
        C = A_{i}';
%         C = (A_{i}+A_nP{i})';
        if size(C,1)>1
            SampleAss = -1*(nanmean(C(:,tb>0 & tb<=10),2)>nanmean(C(:,tb>10 & tb<=15),2));
            ChoiceAss = +1*(nanmean(C(:,tb>0 & tb<=10),2)<nanmean(C(:,tb>10 & tb<=15),2));
            FracSample(i) = nanmean(SampleAss+ChoiceAss);
        else
            FracSample(i) = NaN;
        end
    end
%     FracSample(isinf(FracSample))=1;

    if iArea<3
        col__ = col_{2};
    else
        col__ = col_{3};
    end
    FracSample_{iArea}=FracSample;
    barh(iArea,nanmean(FracSample),'BaseValue',0.,'EdgeColor','k','FaceColor',col__,'LineWidth',1.5,'showbaseline','off') 
    errorbar_x(nanmean(FracSample),iArea,nansem(FracSample),'k')
end
plot([0 0],[0.5 3.5],':k','LineWidth',2.5) 
set(gca,'ytick',[1:3],'YTickLabel',{'Within-mPFC','Within-dCA1','dCA1-mPFC'})
set(gca,'xtick',[-1:1],'xTickLabel',[0:0.5:1])
% set(gca,'xtick',[-1:1],'xTickLabel',{'More ''Sample'' Assemblies','50:50','More ''Choice'' Assemblies'})
axis([-1 1 0.5 3.5])
%% Counts of Sample and Choice assems - discrete ****
figure; hold on
for iArea =1:3
    
    A_= AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
    A_nP= AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea}));
    FracSample=[];
    for i=1:length(A_)
        C = A_{i}';
%         C = (A_{i}+A_nP{i})';
        
            SampleAss = -1*(nanmean(C(:,tb>0 & tb<=10),2)>nanmean(C(:,tb>10 & tb<=15),2));
            ChoiceAss = +1*(nanmean(C(:,tb>0 & tb<=10),2)<nanmean(C(:,tb>10 & tb<=15),2));
            FracSample = [FracSample ; -1*SampleAss];
    end
    sum(FracSample)
%     FracSample(isinf(FracSample))=1;

    if iArea<3
        col__ = col_{2};
    else
        col__ = col_{3};
    end
    FracSample_{iArea}=FracSample;
    barh(iArea,nanmean(FracSample),'BaseValue',0.,'EdgeColor','k','FaceColor',col__,'LineWidth',1.5,'showbaseline','off') 
    errorbar_x(nanmean(FracSample),iArea,nansem(FracSample),'k')
end
plot([0 0],[0.5 3.5],':k','LineWidth',2.5) 
set(gca,'ytick',[1:3],'YTickLabel',{'Within-mPFC','Within-dCA1','dCA1-mPFC'})
set(gca,'xtick',[-1:1],'xTickLabel',[0:0.5:1])
% set(gca,'xtick',[-1:1],'xTickLabel',{'More ''Sample'' Assemblies','50:50','More ''Choice'' Assemblies'})
axis([-1 1 0.5 3.5])

%%
x = []
for iArea =1:3
    x=[x;[FracSample_{iArea}',iArea*ones(length(FracSample_{iArea}'),1)]]
end
anova1(x(:,1),x(:,2))
kruskalwallis(x(:,1),x(:,2))
%%
% rmpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\functions\nansuite');
rmpath('/Users/domansa/Dropbox/MATLAB/PC MATLAB Path/functions/nansuite')

clear h p ci stats
for iArea =1:3
    [h(iArea),p(iArea),ci{iArea},stats{iArea}] =ttest(FracSample_{iArea},0.5)
end
%% Fractions of Sample and Choice assems - discrete **** Bars 
fnames=[];
for iFile =1:length(fileList)
    %% Get the single unit files
    fname=strtok(fileList(iFile).name,'_');
    i = strfind (fname,'LONG');
    fname(i:i+3)=[];
    fnames{iFile} = fname;
end
figure
for iArea =1:3;
    subplot(1,3,iArea); hold on
    A_= AssOutput.Pref{iArea};
    A_nP= AssOutput.NonPref{iArea};
    FracSample=[];
    for i=1:length(A_)
        C = A_{i}';
%         C = (A_{i}+A_nP{i})';
        if size(C,1)>1
        SampleAss = (mean(C(:,tb>0 & tb<=10),2)>mean(C(:,tb>10 & tb<=15),2));
        ChoiceAss = (mean(C(:,tb>0 & tb<=10),2)<mean(C(:,tb>10 & tb<=15),2));
        FracSample(i,:) =[sum(SampleAss),sum(ChoiceAss)];
        else
            FracSample(i,:) = [NaN,NaN];
        end
    end
%     FracSample(isinf(FracSample))=1;
    b = bar(FracSample,'stacked');
    b(1).FaceColor='g';
    b(2).FaceColor='r';
    b(1).EdgeColor='none';
    b(2).EdgeColor='none';
    axis([0.5 length(A_)+0.5 0 10])
    set(gca,'XTick',1:length(fileList),'XTickLabel',fnames,'XTickLabelRotation',45)
    if iArea==2
        xlabel('Recording no.')
    elseif iArea==3
        legend('Sample-active','Choice-active'); legend boxoff
    end
    ylabel('No. Assemblies')
end
 

%% Fractions of Sample and Choice assems - continuous
figure; hold on
for iArea =1:3;
    
    A_= AssOutput.Pref{iArea}(~isempty_cell(AssOutput.Pref{iArea}));
    A_nP= AssOutput.NonPref{iArea}(~isempty_cell(AssOutput.NonPref{iArea}));
    FracSample=[];
    for i=1:length(A_)
        C = A_{i}';
%         C = (A_{i}+A_nP{i})';
        SampleAss  = mean(C(:,tb>0 & tb<=10),2);
        ChoiceAss  = mean(C(:,tb>10 & tb<=15),2);
        FracSample = [FracSample;((ChoiceAss-SampleAss)./(SampleAss+ChoiceAss))];
    end
%     FracSample(isinf(FracSample))=1;

    if iArea<3
        col__ = col_{2};
    else
        col__ = col_{3};
    end
    barh(iArea,nanmean(FracSample),'BaseValue',0.,'EdgeColor','k','FaceColor',col__,'LineWidth',1.5,'showbaseline','off') 
    errorbar_x(nanmean(FracSample),iArea,nansem(FracSample),'k')
end
plot([0 0],[0.5 3.5],':k','LineWidth',2.5) 
set(gca,'ytick',[1:3],'YTickLabel',{'Within-mPFC','Within-dCA1','dCA1-mPFC'})
% set(gca,'xtick',[-1:1],'xTickLabel',[0:0.5:1])
set(gca,'xtick',[-1:1],'xTickLabel',{'More ''Sample'' Assemblies','50:50','More ''Choice'' Assemblies'})
% axis([-1 1 0.5 3.5])

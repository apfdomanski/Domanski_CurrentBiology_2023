%% %%%%%% PREAMBLE %%%%%%
clear
Target = 'LONG';
bw=0.05;
clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);

fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
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
pat2 = 'C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\LONGTaskonly';
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};
threshold_ = 3;
AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time
switch AssemblyChoice
    case 1
        pat2 = [pat 'KDE_bins' filesep Target filesep];
    case 2
        pat2 = [pat 'KDE_binsTaskonly' filesep 'LONGTaskonly' filesep];
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end
%% Load files
iFile =1;
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
    %% Prune any factors with less than 2 units (or single area factors for dCA1-mPFC
    nAss = cellfun(@(x) size(x,2), FSCsel);
    for iArea=1:3
        
        i = nAss(iArea);
        FL_ = FL{iArea}{i};
        thresh_H = ciHld(iArea,threshold_);
        FLthresh= FL_>thresh_H;
        idx =sum(FLthresh,1)<2;
        if iArea==3
            idx =idx | (sum(FLthresh(1:length(Ass.usel_out{1}),:))<1 | ...
            sum(FLthresh((length(Ass.usel_out{1})+1):end,:))<1);
        end
        FL{iArea}{i}(:,idx)=[];
        FSCsel{iArea}(:,idx)=[];
        Ass.units{iArea}(idx)=[];
    end
    %% Prepare datastructure for plotting 
Factor.FSC = FSCsel;
Factor.ciHsc = ciHsc;
Factor.TmtxS = Tmtx;
Factor.units = Ass.units;
Factor.unitIDs = Ass.usel_out;
Factor.ciLld = ciLld;
Factor.ciHld = ciHld;

Factor.FL = FL;
%% Remove double counted tb entries

[TmtxS_,ia]=unique(Factor.TmtxS);
for iArea=1:3
    FSC_{iArea} = Factor.FSC{iArea}(ia,:); 
end
Factor.FSC = FSC_;
Factor.TmtxS = TmtxS_;

%% Plot Grand spike raster figure
spikes = {PFCcells(Ass.usel_out{1}),HPcells(Ass.usel_out{2})};
spikes{3} = [spikes{1};spikes{2}];
t_.Left = [t.Short.SamplePress_LeftCorrect,t.Short.ChoicePress_LeftCorrect;...
           t.Medium.SamplePress_LeftCorrect,t.Medium.ChoicePress_LeftCorrect;...    
           t.Long.SamplePress_LeftCorrect,t.Long.ChoicePress_LeftCorrect];

t_.Right = [t.Short.SamplePress_RightCorrect,t.Short.ChoicePress_RightCorrect;...
            t.Medium.SamplePress_RightCorrect,t.Medium.ChoicePress_RightCorrect;...    
            t.Long.SamplePress_RightCorrect,t.Long.ChoicePress_RightCorrect];
% plot_win = [0 Inf];
firstSample = min(min([t_.Left;t_.Right]))/1e6;
% plot_win = [firstSample-10 firstSample+1200];
plot_win = [4652 4772]; % Animal 2
AssemblyRasterPlotDiscrete(Factor,[],spikes,t_,plot_win)          
%% Plot Grand spike raster figure  Joint only
spikes = {PFCcells(Ass.usel_out{1}),HPcells(Ass.usel_out{2})};
spikes{3} = [spikes{1};spikes{2}];
t_.Left = [t.Short.SamplePress_LeftCorrect,t.Short.ChoicePress_LeftCorrect;...
           t.Medium.SamplePress_LeftCorrect,t.Medium.ChoicePress_LeftCorrect;...    
           t.Long.SamplePress_LeftCorrect,t.Long.ChoicePress_LeftCorrect];

t_.Right = [t.Short.SamplePress_RightCorrect,t.Short.ChoicePress_RightCorrect;...
            t.Medium.SamplePress_RightCorrect,t.Medium.ChoicePress_RightCorrect;...    
            t.Long.SamplePress_RightCorrect,t.Long.ChoicePress_RightCorrect];
% plot_win = [0 Inf];
firstSample = min(min([t_.Left;t_.Right]))/1e6;
% plot_win = [firstSample-10 firstSample+1200];
plot_win = [4652 4772]; % Animal 2
AssemblyRasterPlotDiscreteJointOnly(Factor,[],spikes,t_,plot_win)          
%% Plot Grand spike Image figure  Joint only
spikes = {PFCcells(Ass.usel_out{1}),HPcells(Ass.usel_out{2})};
spikes{3} = [spikes{1};spikes{2}];
t_.Left = [t.Short.SamplePress_LeftCorrect,t.Short.ChoicePress_LeftCorrect;...
           t.Medium.SamplePress_LeftCorrect,t.Medium.ChoicePress_LeftCorrect;...    
           t.Long.SamplePress_LeftCorrect,t.Long.ChoicePress_LeftCorrect];

t_.Right = [t.Short.SamplePress_RightCorrect,t.Short.ChoicePress_RightCorrect;...
            t.Medium.SamplePress_RightCorrect,t.Medium.ChoicePress_RightCorrect;...    
            t.Long.SamplePress_RightCorrect,t.Long.ChoicePress_RightCorrect];
% plot_win = [0 Inf];
firstSample = min(min([t_.Left;t_.Right]))/1e6;
plot_win = [firstSample-10 firstSample+1200];
% plot_win = [4000 5000]; % Animal 2
% plot_win = [4650 4768]; % Animal 2
% plot_win = [5395 5700]; % Animal 2
AssemblyImagePlotDiscreteJointOnly(Factor,iFR_,spikes,t_,plot_win)          

%% Plot oscillating FSCs with event markers

 AssemblyPlotDiscrete(Factor,[],spikes,t_,plot_win)          
 
%% Plot Factor loading stem plot - Horizontal 

nAss_ = cellfun(@(x) size(x,2), Factor.FSC);
nAss
for iArea =1:3
    i = nAss(iArea);
    i_ = nAss_(iArea);
    FL_ = FL{iArea}{i};
    thresh_H = ciHld(iArea,threshold_); 
    thresh_L = ciLld(iArea,threshold_); 
    subplot(1,3,iArea); hold on
    stem(FL_(:,1));
    plot([1 length(FL_)],thresh_H*[1 1],':k','LineWidth',1.5)
    plot([1 length(FL_)],thresh_L*[1 1],':k','LineWidth',1.5)
end
%% Plot Factor loading stem plot - Horizontal Vertical
threshold_ = 3;
nAss_ = cellfun(@(x) size(x,2), Factor.FSC);
for iArea =1:3
    uList = 1:length(spikes{iArea});
    i = nAss(iArea);
    FL_ = FL{iArea}{i};
    thresh_H = ciHld(iArea,threshold_); 
    thresh_L = ciLld(iArea,threshold_); 
    subplot(1,3,iArea); hold on
    if iArea==2
        Yoffset=length(spikes{1});
    else
        Yoffset=0;
    end
    for iAss = 1:nAss_(iArea)
        sig_ = FL_(:,iAss)>thresh_H;
        FLsig = FL_(sig_,iAss);
        IDs = uList(sig_);
        bump_ = iAss;
        plot([bump_ bump_],Yoffset+[0,length(spikes{iArea})],'k')
        scatter(FL_(:,iAss)+bump_,uList+Yoffset,5,'k')
        plot(bump_+[zeros(length(spikes{iArea}),1),FL_(:,iAss)]',...
                  repmat((uList+Yoffset),2,1),'k')
              
        scatter(FLsig+bump_,IDs+Yoffset,10,'g','filled')
        plot(bump_+[zeros(size(FLsig)) FLsig]',...
                  Yoffset+[IDs ;IDs],'g')
%         plot(iAss+thresh_H*[1 1],[1 length(spikes{3})],':k','LineWidth',1.5)
%         plot(iAss+thresh_L*[1 1],[1 length(spikes{3})],':k','LineWidth',1.5)
    end
    axis([0 max(nAss_)+1 0 length(spikes{3})])
end

%% Plot Assembly activation

function DecodingDuringDelayCrossTemporalBatch2(iFile)
if isstr(iFile)
    iFile = str2num(iFile);
end
%% %%%%%% PREAMBLE %%%%%%

Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};

Target = 'LONG';
RunErrors = false;
tlimsShort=[0 4];
tlimsMedium=[0 8];
tlimsLong=[0 16];
tlimsShort=[-5 9];
tlimsMedium=[-5 13];
tlimsLong=[-5 21];
plotOnline = false;
bw    = 0.05;
Nbs = 500;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);

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

normaliseFscores = false;
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};

AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time
switch AssemblyChoice
    case 1
        pat2 = [pat 'KDE_bins' filesep Target filesep];
    case 2
        pat2 = [pat 'KDE_binsTaskOnly' filesep 'LONGTaskonly' filesep];
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end
mkdir([pat 'MixedSelectivity'])

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};

TrialsSubsample = Inf;
if isinf(TrialsSubsample)
    nReps = 1;
else
    nReps =500;
end
%% load files
fname=strtok(fileList(iFile).name,'_');

load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
load(sprintf('%s%s%s.mat',pat,filesep,fname));
nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu);

%% Batch process units
for iArea = 1%1:length(Areas)
    %% Get the files
    fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
    if iArea < 3
        % load(sprintf('%sKDE_bins%s%s_%s_iFR50.mat',pat,filesep,fname,Areas{iArea}));
        load(sprintf('%sKDE_binsTaskOnly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
        iFR_ = iFR;
    else
        U_{1} = load(sprintf('%sKDE_binsTaskOnly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{1}));
        U_{2} = load(sprintf('%sKDE_binsTaskOnly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{2}));
        iFR_ = [U_{1}.iFR,U_{2}.iFR];
        Tmtx = U_{1}.Tmtx;
        clear U_
    end
    
    % iFR_ = zscore(iFR_);
    %% Delay period decoding for units
    shuffleCorrect = true;
    for iDelay =1:length(Delays_)
        % run correct trials
        eval(sprintf('LeftTrials = [t.%s.SamplePress_LeftCorrect];',Delays_{iDelay}));
        eval(sprintf('RightTrials = [t.%s.SamplePress_RightCorrect];',Delays_{iDelay}));
        eval(sprintf('tlims_X = tlims%s;',Delays_{iDelay}));
        length_ = sum(abs(tlims_X))/bw;
        eval(sprintf('tb_ = tb%s;',Delays_{iDelay}))
        
        Ltrials = [];nL = 0; Ltrials_={};
        for iTrial =1:size(LeftTrials,1)
            try
                tlims_  = LeftTrials(iTrial,1)/1e6+tlims_X(1);
                tlims_  = closest(Tmtx,tlims_);
                tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                Ltrials = [Ltrials;iFR_(tlims_,:)];
                Ltrials_{iTrial} = iFR_(tlims_,:);
                nL=nL+1;
            end
        end
        Rtrials = [];nR = 0;  Rtrials_={};
        for iTrial =1:size(RightTrials,1)
            try
                tlims_  = RightTrials(iTrial,1)/1e6+tlims_X(1);
                tlims_  = closest(Tmtx,tlims_);
                tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                Rtrials = [Rtrials;iFR_(tlims_,:)];
                Rtrials_{iTrial} = iFR_(tlims_,:);
                nR=nR+1;
            end
        end
        
        flag_ = true;
        if ~isinf(TrialsSubsample)
            if sum([nL,nR]<=TrialsSubsample)==2
                flag_ = false;
            end
        end
        clear PE PE_CIl PE_CIh
        if flag_
            
            FR   = [Ltrials;Rtrials];
            evt0 = [ones(nL,1);2*ones(nR,1)];
            reg  = 0.05;
            [PE,PE_CIl,PE_CIh] = DecodePredErrCrossTemporal_CV_GPU(FR,evt0,reg,true);%(:,Ass.JointMembers{iArea})
            D_{iArea}{iDelay}    = PE;
            D_CIh{iArea}{iDelay} = PE_CIh;
            D_CIl{iArea}{iDelay} = PE_CIl;
        else
            D_{iArea}{iDelay}    = nan(length_+1,length_+1);
            D_CIh{iArea}{iDelay} = nan(length_+1,length_+1);
            D_CIl{iArea}{iDelay} = nan(length_+1,length_+1);
        end
        
        % Run error trials
                eval(sprintf('LeftTrialsE = [t.%s.SamplePress_LeftError''];',Delays_{iDelay}));
                eval(sprintf('RightTrialsE = [t.%s.SamplePress_RightError''];',Delays_{iDelay}));
                LtrialsE = [];nLe = 0; LtrialsE_={};
                for iTrial =1:size(LeftTrialsE,1)
                    try
                        tlims_  = LeftTrialsE(iTrial,1)/1e6+tlims_X(1);
                        tlims_ = closest(Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                        LtrialsE = [LtrialsE;iFR_(tlims_,:)];
                        LtrialsE_{iTrial} = iFR_(tlims_,:);
                        
                        nLe=nLe+1;
                    end
                end
                RtrialsE = [];nRe = 0; RtrialsE_={};
                for iTrial =1:size(RightTrialsE,1)
                    try
                        tlims_  = RightTrialsE(iTrial,1)/1e6+tlims_X(1) ;
                        tlims_ = closest(Tmtx,tlims_);
                        tlims_  = [tlims_(1):(tlims_(1)+length(tb_)-1)];
                        RtrialsE = [RtrialsE;iFR_(tlims_,:)];
                        RtrialsE_{iTrial} = iFR_(tlims_,:);
                        nRe=nRe+1;
                    end
                end
                flag_ = true;
                if ~isinf(TrialsSubsample)
                    if sum([nLe,nRe]<=TrialsSubsample)==2
                        flag_ = false;
                    end
                end
                if flag_
                    FR = [LtrialsE;RtrialsE];
                    evt0 = [ones(nLe,1);2*ones(nRe,1)];
                    reg = 0.05;
                    [PE,PE_CIl,PE_CIh]      = DecodePredErrCrossTemporal_CV_GPU(FR,evt0,reg,true);%(:,Ass.JointMembers{iArea})
                    D_Err{iArea}{iDelay}    = PE;
                    D_ErrCIh{iArea}{iDelay} = PE_CIh;
                    D_ErrCIl{iArea}{iDelay} = PE_CIl;
                else
                    D_Err{iArea}{iDelay}    = nan(length_+1,length_+1);
                    D_ErrCIh{iArea}{iDelay} = nan(length_+1,length_+1);
                    D_ErrCIl{iArea}{iDelay} = nan(length_+1,length_+1);
                end
                    
                   
        
    end
    %% Save results
    if exist('D_')
        fnOut = sprintf('%sMixedSelectivity%s%s_%s_DelayDecoding_UnitSpan_ErrorBoundsWithErr.mat',pat,filesep,fname,Areas{iArea});
        save(fnOut,'D_','D_CIh','D_CIl','D_Err','D_ErrCIh','D_ErrCIl','tbShort','tbMedium','tbLong','-v7.3');
        fprintf('Done.\n')
        clear PE PE_CIl PE_CIh
    end
    
end

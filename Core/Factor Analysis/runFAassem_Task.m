clear all
close all
target= 'Task';
if ispc
    home = 'C:\'; % [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    pat = [home 'Analysis\AssemblyAnalysis\Sleep'];
elseif ismac
    path(path,'/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/MichalDataAna')
    home = getenv('HOME');
    cd ([home '/Documents/Bristol_work/assembly_analysis/MichalDataAna'])
    pat = [home '/Documents/Bristol_work/assembly_analysis/Michal_data/raw/50ms_bins/'];
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/Sleep'];
end

if ~exist([pat filesep target])
    mkdir([pat filesep target])
end
% run on 1) iFRsc, 2) KDE iFR

%bw=0.01; iFRtyp='iFRsc';
%bw=0.05; iFRtyp='iFRsc';

bw=0.05; iFRtyp='iFR';
A=dir([pat  filesep '*PFC_iFR' num2str(bw*1e3) strcat('*', target) '*.mat']);

for event_id=1:length(A)
    k=findstr(A(event_id).name,'_'); r=findstr(A(event_id).name,'.');
    Name{event_id}=[A(event_id).name(1:k(1)) A(event_id).name(k(2)+1:r-1)];
end;
minFR=0.1;     % minimal acceptable firing rate
critCvBW=1e6;  % critical max bound on kernel width drift over time (variance)
alpha=0.01;    % significance threshold
kmax=30;       % number of factors to check up to
Nbs=500;       % number of bootstrap draws
twin=10;       % Time window (s)
twin = twin/bw;
%% run Factor analysis
for f=1:length(A)
    % output vector constructors:
    % Factor analysis variables
    FL     = cell(1,3);     FSC    = cell(1,3);     psix    = cell(1,3);
    % Metrics for number of factors choices
    nassem = cell(1,3); 
    LL     = cell(1,3); AIC = cell(1,3); BIC=cell(1,3); Pr = cell(1,3); Chi2 = cell(1,3); 
    % Bootstrapped variables
    perm   = cell(1,3); 
    FLbs   = cell(1,3);     FSCbs  = cell(1,3);     psixBS  = cell(1,3); 
    LLbs   = cell(1,3);     prBS   = cell(1,3);     Chi2bs  = cell(1,3); nuBS=cell(1,3);  
    
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fn=[pat filesep A(f).name];
    [TmtxS,SCM,unit_IDs]=SelCellsSleep(fn,twin,minFR,critCvBW,'iFR');
%   [TmtxS,FiringRates,EvtLs,EvtTs,unit_IDs] = SelCellsSleep(fn,twin,minFR,critCvBW,iFRtype)
    
    for area_id=1:2  % ...loop across all brain areas
        for event_id=1:length(TmtxS{area_id}) % ...loop across all events
            TmtxS{area_id}{event_id} = TmtxS{area_id}{event_id}(1:end);
            SCM{area_id}{event_id}=SCM{area_id}{event_id}(1:end,:)';
        end;
    end;
    % (1) run FA-based assembly detection algo separately for each area
    for s=1:2
        [nassem{s},...
         FL{s},~,...
         LL{s},AIC{s},BIC{s},Pr{s},Chi2{s}, ...
         prBS{s},LLbs{s},Chi2bs{s},psix{s},FLbs{s},~,psixBS{s},perm{s}]  =   FAassem_(SCM{s},kmax,Nbs,alpha);
    end;
    
    % (2) combine PFC and HP, run Factor analysis again
    s=3; SCMcomb=cell(1,length(SCM{1}));
    for event_id=1:length(SCM{1}), SCMcomb{event_id}=[SCM{1}{event_id};SCM{2}{event_id}]; end;
        [nassem{s},FL{s},~,LL{s},AIC{s},BIC{s},Pr{s},Chi2{s}, ...
                prBS{s},LLbs{s},Chi2bs{s},psix{s},FLbs{s},~,psixBS{s},perm{s}]      =   FAassem_(SCMcomb,kmax,Nbs,alpha);
  
 fnOut=[pat filesep target filesep Name{f} '_AssemRes2'];
    save(fnOut,'TmtxS','nassem','FL','LL','AIC','BIC', ...
        'Pr','prBS','LLbs','FLbs','Chi2bs','Chi2','psix','psixBS','perm');
end
%% compute factor scores and BS confidence limits
for f=1:length(A)
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fn=[pat filesep A(f).name];
    [TmtxS,SCM,unit_IDs]=SelCellsSleep(fn,twin,minFR,critCvBW,'iFR');
    fnOut=[pat filesep target filesep Name{f} '_AssemRes2'];
    [FSC,FSCbs,ciLld,ciHld,ciLsc,ciHsc]=BSciSleep(TmtxS,SCM,fnOut);
    k=findstr(fnOut,'_');
    save([fnOut(1:k(end)) '_FSCtemp']','FSC','FSCbs','ciLld','ciHld','ciLsc','ciHsc');
end;

%% extract assemblies (choose BS as reasonable but relatively conservative)
% extract assembly activations as f(t) through factor scores
for f=1:length(A)
    fnOut=[pat filesep target filesep Name{f} '_AssemRes2']; load(fnOut);
    fprintf([Name{f}, ' AssemRes2: %G, %G, %G, factors.\n'],nassem{1}(3),nassem{2}(3),nassem{3}(3))
    fnOut = [pat filesep target filesep Name{f} '__FSCtemp']; load(fnOut);    
    units=cell(1,length(nassem)); FSCsel=units; Tmtx=FSCsel;
    
    % Process units
    for s=1:length(nassem)
        numf=nassem{s}(3);
        for n=1:numf, 
            units{s}{n}=find(abs(FL{s}{numf}(:,n))>ciHld(s,3))'; 
        end;
        % Weed out factors with only one unit
        k=find(cellfun(@length,units{s})>1);
        keepers{s} = k; % Indices of factors to keep
        units{s}=units{s}(k);
        FSCsel{s}=FSC{s}(:,k);
    end;

    nu=[]; 
    for s=1:3 
        nu(s)=size(FL{s}{1},1); 
    end;

    TmtxS{3} = TmtxS{1};
    Tmtx=cell2mat(TmtxS{1})';

    minFR=0.1;     % minimal acceptable firing rate
    critCvBW=1e6;  % critical max bound on kernel width drift over time (variance)
    twin=10;       % Time window (s)
    twin = twin/bw;
	fn=[pat filesep A(f).name];
    [~,~,unit_IDs]=SelCellsSleep(fn,twin,minFR,critCvBW,'iFR');
    
    k=findstr(fnOut,'_');
    save([fnOut(1:k(end-1)) 'FSC'],'units','FSCsel','Tmtx','nassem', ...
        'ciLld','ciHld','ciLsc','ciHsc','nu','unit_IDs','keepers');
end

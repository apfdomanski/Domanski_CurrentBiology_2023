clear all
close all
target= 'PreSleep';
if ispc
    home = 'C:\'; % [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    pat = [home 'Analysis\AssemblyAnalysis\Sleep'];
    separator='\';
elseif ismac
    path(path,'/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/MichalDataAna')
    home = getenv('HOME');
    cd ([home '/Documents/Bristol_work/assembly_analysis/MichalDataAna'])
    pat = [home '/Documents/Bristol_work/assembly_analysis/Michal_data/raw/50ms_bins/'];
    separator='/';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB/MichalDataAna')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/Sleep'];
    separator='/';
end

mkdir([pat separator target])
% run on 1) iFRsc, 2) KDE iFR

%bw=0.01; iFRtyp='iFRsc';
%bw=0.05; iFRtyp='iFRsc';

bw=0.05; iFRtyp='iFR';
A=dir([pat  separator '*PFC_iFR' num2str(bw*1e3) strcat('*', target) '*.mat']);

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

%% compute factor scores and BS confidence limits
for f=1:length(A)
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fn=[pat separator A(f).name];
    [TmtxS,SCM,unit_IDs]=SelCellsSleep(fn,twin,minFR,critCvBW,'iFR');
    fnOut=[pat separator target separator Name{f} '_AssemRes2'];
    [FSC,FSCbs,ciLld,ciHld,ciLsc,ciHsc]=BSciSleep(TmtxS,SCM,fnOut);
    k=findstr(fnOut,'_');
    save([fnOut(1:k(end)) '_FSCtemp']','FSC','FSCbs','ciLld','ciHld','ciLsc','ciHsc');
end;

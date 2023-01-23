% Code to run spike train Xcorr function in parallel on Blue
% Crystal.

% APFD UoB March 2016

% Dependencies: restrictTranges.m , catpad.m, CompCCMtx.m, CCMtx (MEX)

clear all
close all

p.rat = 'JaroslawLONG2';

if ispc
    p.home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    p.pat = 'C:\Analysis\AssemblyAnalysis\Sleep\';
    cd ([p.home '\Documents\MATLAB\Default dependencies\AssemblyCode\'])
elseif ismac
%     path(path,'/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/MichalDataAna')
%     home = getenv('HOME');
%     cd ([home '/Documents/Bristol_work/assembly_analysis/MichalDataAna'])
%     pat = [home '/Documents/Bristol_work/assembly_analysis/Michal_data/raw/KDE_bins/'];
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    p.home = getenv('HOME');
    cd ([p.home '/MATLAB/SpikeTrainAna/'])
    p.pat = [p.home '/Sleep/'];    
end

mkdir([p.pat '/' 'SpikeCorrs'])
load ([p.pat p.rat '_analysis.mat'],'p','FR','tRanges')

%% Set up for cross-correlograms
% concatenate all unit spike times
FR.Jointcells = [FR.PFCcells;FR.HPcells];
FR.iCells ={1:length(FR.PFCcells),(1:length(FR.HPcells))+length(FR.PFCcells)};
FR.unitLocs = [ones(length(FR.PFCcells),1) ; 2*ones(length(FR.HPcells),1)];
epochs ={'Pre','Task','Post'};
STMtx{1}  = restrictTranges(FR.Jointcells,tRanges.Pre)';
STMtx{2} = restrictTranges(FR.Jointcells,tRanges.Task)';
STMtx{3} = restrictTranges(FR.Jointcells,tRanges.Post)';
% [no. Units, no. Timestamps]

% Can restrict to equal spike no's with this:
% mxsize=min([size(STMtxIN,1) size(STMtxOUT,1)]);
% if size(STMtxIN,1)>=mxsize;    STMtxIN  = STMtxIN(1:mxsize,:);  end


%% Run cross-correlation code
clear UUcorr

p.xcorr.Tcorr = 0.05;
p.xcorr.bw    = 0.001;
p.xcorr.optf  = 2;  
%  0: raw CCTH, 
%  1: phase histo, 
%  2: cc/(n1*n2), 
%  3: true corr 

Tinv{1}  = [min(STMtx{1}(STMtx{1}>0)),max(max(STMtx{1}))] ;
Tinv{2} = [min(STMtx{2}(STMtx{2}>0)),max(max(STMtx{2}))] ;
Tinv{3} = [min(STMtx{3}(STMtx{3}>0)),max(max(STMtx{3}))] ;

% Optional: check for low FR neurons to exclude
% Ncrit=100; unit_mask = ones(size(FR.Jointcells));
% for j=1:size(STMtxPre,1)
%     k=find(STMtxPre(j,:)>=Tinv(1) & STMtxPre(j,:)<=Tinv(2));
%     if length(k)<Ncrit
%         unit_mask(j,:)=0;
%         disp(['excluded unit# ' num2str(j)])
%     end;
% end;

parpool local

parfor iEpoch = 1:length(epochs)

    [UUcorr{iEpoch}.Tb,...
     UUcorr{iEpoch}.CC,...
     UUcorr{iEpoch}.CCMax,...
     UUcorr{iEpoch}.CCDel,...
     UUcorr{iEpoch}.CCSig,...
     UUcorr{iEpoch}.CCSum,...
     UUcorr{iEpoch}.ccmxSurr,...
     UUcorr{iEpoch}.CCSurr]=CompCCMtx(STMtx{iEpoch}, Tinv{iEpoch}, p.xcorr.Tcorr, p.xcorr.bw, p.xcorr.optf);

end

   fnOut=[p.pat 'SpikeCorrs\' p.rat '_Corrs.mat'];
   save(fnOut,'UUcorr');



% matlabpool close
%%
plot(UUcorr{2}.Tb,UUcorr{2}.CC(1:2000,:))
ylim([0 0.0005])

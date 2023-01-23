% Decode Left/Right preference of assemblies.
%
% Reports decoding metrics of both each assemblies and assebly population
% with both F-score and Cross-validation error metrics. Bootsrapped
% F-decoding error also calculated.
%
% (c) 2014 Daniel Durstewitz, Bernstein Center for Computational Neuroscience, CIMH/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk


clear all
close all
target= 'LONG';
if ispc
    home = 'C:/'; % [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    pat = [home 'Analysis\AssemblyAnalysis\raw\KDE_bins\' target];
    cd(pat)
elseif ismac
    path(path,'/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/MichalDataAna')
    home = getenv('HOME');
    cd ([home '/Documents/Bristol_work/assembly_analysis/MichalDataAna'])
    pat = [home '/Documents/Bristol_work/assembly_analysis/Michal_data/raw/KDE_bins/'];
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/KDE_bins/'];    
end
dirOut = [pat filesep target];
if ~isdir(dirOut)
    mkdir([pat filesep target])
end

bw=0.05; 
iFRtyp='iFR';
A=dir([pat '/' strcat('*', target,'*') '*FSC.mat']);
% remove error trial datasets
Ignore= 'inc';
A = A(cellfun(@isempty,cellfun(@(x,a) strfind(x,a),...
    {A.name},repmat({Ignore},1,length(A)),'UniformOutput',false)));



disp(['...Working path: ' pat '.'])
disp(['...Output path: ' dirOut '.'])
fprintf('...Working on %G recordings.\n',length(A))
fprintf(['...Ignoring files with "' Ignore '" in the title.\n'])
%%
reg=0.05;
twin=10; 
bw=0.05;
Nbs=1000; 

prFt2=cell(1,3); 
Rt2=cell(1,3); 
Ft2=cell(1,3);
Ft2ciL=cell(1,3); 
Ft2ciH=cell(1,3); 
TS=cell(1,3); 
dfnum=cell(1,3); 
dfden=cell(1,3);
Ft2bs=cell(1,3); 
Ft2BSciH=cell(1,3); 
Ft2BSciL=cell(1,3);
CVE=cell(1,3); 
CVEbs=cell(1,3); 
cveBSciH=cell(1,3); 
cveBSciL=cell(1,3);
for f=1:length(A)
    fn=[pat '/' A(f).name];
    fprintf(['...Working on file %G of %G: ',A(f).name '\n'], f , length(A))
    [TmtxS,FSCs,EvtLs,EvtTs]=SelTrialsFSCs(fn,twin);
    ko=round((twin+5)/bw);
    FSCsel=cell(1,3); n_assem=zeros(1,3);
    % extract equal length segments for each trial
    for s=1:3 % areas
        for i=1:length(TmtxS{s}) % trials
            TmtxS{s}{i}=TmtxS{s}{i}([1:ko end-ko+1:end]);
            FSCs{s}{i}=FSCs{s}{i}([1:ko end-ko+1:end],:)';
        end;
        FSCsel{s}=cell2mat(FSCs{s})'; n_assem(s)=size(FSCsel{s},2);
    end;
    Tmtx=cell2mat(TmtxS{1})';

    ntr=length(EvtTs)/2; Ltr=round(size(FSCsel{1},1)/ntr);
    evt0=EvtLs(((1:ntr)-1)*2+1)-2;
    for s=1:3       
        if n_assem(s)~=0
        %% (1) decoding in terms of F-score
        [~,~,Ft2x,Rt2x,Ft2ciL0,Ft2ciH0,TS{s}{f},dfnum{s}(f),dfd]=DecodeStats(FSCsel{s},evt0,reg);
        prFt2{s}(f,:)=1-fcdf(Ft2x,dfnum{s}(f),dfd);
        Rt2{s}(f,:)=Rt2x;
        Ft2{s}(f,:)=Ft2x;
        Ft2ciL{s}(f,:)=Ft2ciL0;
        Ft2ciH{s}(f,:)=Ft2ciH0;
        dfden{s}(f,:)=dfd;       
        %% (2) decoding in terms of CVE
        if n_assem(s)==min(n_assem(n_assem>0)), nrep=1; else nrep=10; end;
        rs=cell(1,nrep);
        pe0=zeros(nrep,Ltr);
        for i=1:nrep
            rs{i}=randsample(n_assem(s),min(n_assem(n_assem>0)));
            pe0(i,:)=DecodeCVE(FSCsel{s}(:,rs{i}),evt0,reg);
        end;
        CVE{s}(f,:)=mean(pe0,1);        
        %% (3) bootstraps (scramble assignments of trajectories to conditions)
        Ft2BSciH{s}=zeros(1,Ltr); Ft2BSciL{s}=Ft2BSciH{s}; Ft2bs{s}=[];
        cveBSciH{s}=zeros(1,Ltr); cveBSciL{s}=zeros(1,Ltr);
        CVEbs{s}=[];
        for b=1:Nbs
            disp(b);
            k=randperm(length(evt0)); evt1=evt0(k);
            while true
                try
                    [~,~,Ft2bs{s}(b,:)]=DecodeStats(FSCsel{s},evt1,reg);
                    break;
                catch
                    fprintf('no convergence, retry...')
                    k=randperm(length(evt0)); evt1=evt0(k);
                end
            end
            
            %cve0=zeros(nrep,Ltr);
            %for i=1:nrep, cve0(i,:)=DecodeCVE(iFRsel{s}(:,rs{i}),evt1,reg); end;
            %CVEbs{s}(b,:)=mean(cve0,1);
        end;       
        for t=1:Ltr
            Fts=sort(Ft2bs{s}(:,t),'ascend');
            Ft2BSciH{s}(t)=Fts(round(0.95*Nbs));
            Ft2BSciL{s}(t)=Fts(round(0.05*Nbs));
            %cves=sort(CVEbs{s}(:,t),'ascend');
            %cveBSciH{s}(t)=cves(round(0.95*Nbs));
            %cveBSciL{s}(t)=cves(round(0.05*Nbs));
        end
        %%
        end
    end
     k=findstr(A(f).name,'_');
        save([pat '\' target '\' 'ASSEMPopDecodeBS_' A(f).name(1:k) iFRtyp num2str(bw*1e3)], ...
            'Ft2bs','Ft2BSciH','Ft2BSciL','CVEbs','cveBSciH','cveBSciL');
end
%%
Ft2Z=cell(1,3); for s=1:3, Ft2Z{s}=zscore(Ft2{s}')'; end;
Tv=(1:size(Rt2{1},2))*bw-bw/2;
save([pat '\' target '\' 'ASSEMPopDecodeStats_' iFRtyp num2str(bw*1e3)],'prFt2','Rt2','Ft2', ...
    'Ft2Z','Ft2ciL','Ft2ciH','TS','Tv','dfnum','dfden','CVE');


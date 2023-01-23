% Wrapper code: single unit and population decoding between two behavioural outcomes. 
%
% (c) 2014 Daniel Durstewitz, Bernstein Center for Computational Neuroscience, CIMH/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Iterates across all animals/files
% Uses both F-test based LDA decoder (DecodeStats), and leave-one out
% cross-validation error decoder (decodeCVE). Treats error trials as
% special case of the latter.

clear all
close all
% path(path,'/home/dd/bin/MLana')

% run on 1) iFRsc, 2) KDE iFR

% pat='data/';
pat='C:\Analysis\AssemblyAnalysis\raw\KDE_bins\';
% pat=pwd;
%bw=0.01; iFRtyp='iFRsc';
%bw=0.05; iFRtyp='iFRsc';
bw=0.05; iFRtyp='iFR';
A=dir([pat '/*PFC_iFR' num2str(bw*1e3) '*.mat']);
%A=dir(['*.mat']);

reg=0.05;
twin=10; minFR=0.1; critCvBW=1e8;
prFt2=cell(1,2); Rt2=cell(1,2); Ft2=cell(1,2);
Ft2ciL=cell(1,2); Ft2ciH=cell(1,2); TS=cell(1,2); dfnum=cell(1,2); dfden=cell(1,2);
Nbs=1000; Ft2bs=cell(1,2); Ft2BSciH=cell(1,2); Ft2BSciL=cell(1,2);
CVE=cell(1,2); CVEbs=cell(1,2); cveBSciH=cell(1,2); cveBSciL=cell(1,2);
for f=1:length(A)
    % Choose trials/neurons: cut out equally-sized trial periods, kick out bad-behaved neurons
    fn=[pat '\' A(f).name];
    [TmtxS,iFR0,EvtLs,EvtTs]=SelTrialsCells(fn,twin,minFR,critCvBW,iFRtyp);
    ko=round((twin+5)/bw);
    iFRsel=cell(1,2); nu=zeros(1,2);
    for s=1:2
        for i=1:length(TmtxS{s})
            TmtxS{s}{i}=TmtxS{s}{i}([1:ko end-ko+1:end]);
            iFR0{s}{i}=iFR0{s}{i}([1:ko end-ko+1:end],:)';
        end;
        iFRsel{s}=cell2mat(iFR0{s})'; nu(s)=size(iFRsel{s},2);
    end;
    Tmtx=cell2mat(TmtxS{1})';
    
    ntr=length(EvtTs)/2; Ltr=round(size(iFRsel{1},1)/ntr);
    evt0=EvtLs(((1:ntr)-1)*2+1)-2;
    %Tv=(1:Ltr)*bw-bw/2;
    for s=1:2
        %%% (1) decoding in terms of F-score
            % avgFR     {trial type}(time , unit no.) Average firing rate for each event type
            % seFR      {trial type}(time , unit no.) Firing rate squared error for each event type
            % Ft2x      Hotelling's T^2 statistic: discriminabiltiy of 2 classes
            % Rt2x      Adjusted T^2
            % Ft2ciL0   Low 95% confidence limits on T^2
            % Ft2ciH0   High 95% confidence limits on T^2
            % TS        T-score
            % dfnum     F-distribution numerator DofF
            % dfd       F-distribution denominator DofF
            
        [~,~,Ft2x,Rt2x,Ft2ciL0,Ft2ciH0,TS{s}{f},dfnum{s}(f),dfd]=DecodeStats(iFRsel{s},evt0,reg);
        
        prFt2 {s}(f,:)=1-fcdf(Ft2x,dfnum{s}(f),dfd); % p value for Ft2x (T^2 values)
        Rt2   {s}(f,:)=Rt2x;
        Ft2   {s}(f,:)=Ft2x;
        Ft2ciL{s}(f,:)=Ft2ciL0;
        Ft2ciH{s}(f,:)=Ft2ciH0;
        dfden {s}(f,:)=dfd;
        
        %%% (2) decoding in terms of CVE
        if nu(s)==min(nu), nrep=1; else nrep=10; end;
        rs=cell(1,nrep);
        pe0=zeros(nrep,Ltr);
        for i=1:nrep
            rs{i}=randsample(nu(s),min(nu)); % random draw of units
            
            pe0(i,:)=DecodeCVE(iFRsel{s}(:,rs{i}),evt0,reg); % run leave-one out CVE on this draw
            
        end;
        CVE{s}(f,:)=mean(pe0,1); % Collapse results across all draws

        
        %%% (3) bootstraps (scramble assignments of trajectories to conditions)
            % shuffle trial outcomes and run DecodeStats and decodeCVE for
            % each draw (only results from DecodeStats are used here)
        Ft2BSciH{s}=zeros(1,Ltr); Ft2BSciL{s}=Ft2BSciH{s}; Ft2bs{s}=[];
        cveBSciH{s}=zeros(1,Ltr); cveBSciL{s}=zeros(1,Ltr);
        CVEbs{s}=[];
        for b=1:Nbs
            disp(b);
            k=randperm(length(evt0)); evt1=evt0(k);
            [~,~,Ft2bs{s}(b,:)]=DecodeStats(iFRsel{s},evt1,reg);
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
        end;
        k=findstr(A(f).name,'_');
        save([pat 'PopDecodeBS_' A(f).name(1:k) iFRtyp num2str(bw*1e3)], ...
            'Ft2bs',...
            'Ft2BSciH',...
            'Ft2BSciL',...
            'CVEbs',...
            'cveBSciH',...
            'cveBSciL');
    end;
    
end;
Ft2Z=cell(1,2); for s=1:2, Ft2Z{s}=zscore(Ft2{s}')'; end;
Tv=(1:size(Rt2{1},2))*bw-bw/2;
save([pat 'PopDecodeStats_' iFRtyp num2str(bw*1e3)],...
    'prFt2',...  % p value for Ft2x (T^2 values)
    'Rt2',...    % Rt2x      Adjusted T^2
    'Ft2', ...   % Ft2x      Hotelling's T^2 statistic: discriminabiltiy of 2 classes
    'Ft2Z',...   % Z-scored T^2 values
    'Ft2ciL',... % Low CI on T^2 value
    'Ft2ciH',... % High CI on T^2 value
    'TS',...     % T-score
    'Tv',...     % Timebase
    'dfnum',...  % F-distribution numerator DofF
    'dfden',...  % F-distribution denominator DofF
    'CVE');      % CVE decoding results



%% ----- run decoding on error trials ---------------------------------
% As for correct trials, except error trials are signalled by a label of outcome x10
% in evtLs
reg=0.05;
twin=10; minFR=0.1; critCvBW=1e8;
PE=cell(1,2);
for f=1:length(A)
    
    % cut out equally-sized trial periods, kick out bad-behaved neurons
    fn=[pat A(f).name];
    [TmtxS,iFR0,EvtLs,EvtTs]=SelTrialsCells(fn,twin,minFR,critCvBW,iFRtyp,'all');
    ko=round((twin+5)/bw);
    iFRsel=cell(1,2); nu=zeros(1,2);
    for s=1:2
        for i=1:length(TmtxS{s})
            TmtxS{s}{i}=TmtxS{s}{i}([1:ko end-ko+1:end]);
            iFR0{s}{i}=iFR0{s}{i}([1:ko end-ko+1:end],:)';
        end;
        iFRsel{s}=cell2mat(iFR0{s})'; nu(s)=size(iFRsel{s},2);
    end;
    ntr=length(EvtTs)/2; Ltr=round(size(iFRsel{1},1)/ntr);
    
    for s=1:2
        if nu(s)==min(nu), nrep=1; else nrep=10; end;
        rs=cell(1,nrep);
        pe0=zeros(nrep,Ltr);
        for i=1:nrep
            rs{i}=randsample(nu(s),min(nu));
            pe0(i,:)=DecodePredErr(iFRsel{s}(:,rs{i}),EvtLs(2:2:end),reg);
        end;
        PE{s}(f,:)=mean(pe0,1);
    end;
    
end;
Tv=(1:size(PE{1},2))*bw-bw/2;
save([pat 'ErrTrDecodeStats_' iFRtyp num2str(bw*1e3)],'Tv','PE');

% (c) 2014, Daniel Durstewitz, Dept. Theoretical Neuroscience, CIMH
% Mannheim/ Heidelberg University

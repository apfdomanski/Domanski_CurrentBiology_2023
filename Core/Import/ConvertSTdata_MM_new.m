% Prototype function for calculating smoothed estimates of firing rate by
% optimised kernel density smoothing. Also contains a bloc of code to calculate FR for error trials.
%
% (c) 2014 Daniel Durstewitz, Bernstein Center for Computational Neuroscience, CIMH/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% For each unit's spike train in parallel...
% (1) firing rate kernel is optimised using UniKDE:  This routine implements univariate nonparametric density 
% estimation with bandwidth selection, starting with initial guess at mean
% firing rate of that unit
% (2) Then calls a MEX'ed C routine (STtoGfAvg) to convert spike train matrix into Gaussian-filtered activity matrix
% 
% NB this routine cuts out +/- 5s windows around each behavioural event,
% and creates a timebase with gaps in it. If you want the whole experiment
% as continuous data use something like ConvertSTdata_MM_pre_sleep.m 
%
% NB edit and run for each of HP and PFC brain areas. 
%
% Variables (*=exported):
% avgFR*                  Average firing rate
% pat    (str)            data path
% A      {no. cell types} list of recording type ending in HP/PFC cells
% B      {no. neurons}    Cell of spike each neuron's spike times
% i                       recording type pointer
% j                       Condition/trial pointer
% ST     (no. spikes)     temporary spike times vector
% bw*                     KDE bandwidth
% EvtT*                   Event times (samples or choices)
% EvtL*                   Event labels (direction/sample or choice conditions)
% Tmtx*                   Time axis with deltaT = KDE bandwidth
% H0                      optimization initial guess of (mean firing rate)^2
% hucvS  {no. neurons}    optimised kernel width for each trial individually
% hucv*  (no. neurons)    optimized kernel width 
% MISE*  (no. neurons)    mean integrated squared error of optimization
% CvBW*  (no. neurons)    variance in KDE bandwidth 12 trials (used to be choice-specific)
% iFR0   {no. neurons}    Gaussian-convolved spike rates
% iFR*   (time,neurons}   Gaussian-convolved spike rates
% TL     (epochs)         time vector cut up into 12 equal portions (repeated trials)
%

clear all
% path(path,'/home/dd/bin/MLana')
% path(path,'/home/dd/bin/MLana/MUAana')

% generate iFR and event matrices
load JaroslawLONG1.mat; ext='JaroslawLONG1';
load load('D:\Aleks\OnufryLONG2.m.mat; ext='JaroslawLONG1';
%load data/LONG2.mat; ext='L2';g

bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
%A=who('*PFCcells'); br='PFC';
A=who('*HPcells'); br='HP';


parpool local
for i=1:length(A)
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr)    % for each direction/sample condition
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5;

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    for u=1:length(B)
        [i u]
        % collect ST
        ST=B{u}.t*1e-4;
        avgFR(u)=length(ST)/range(ST);
        
        % convert STMtx into iFR via KDE
        h0=mean(diff(ST))^2;
        [hucv(u),MISE(u)]=UniKDE(ST,h0,3);
        if ~isnan(hucv(u)), iFR0{u}=STtoGfAvg(ST',Tmtx,sqrt(hucv(u)))';
        else iFR0{u}=zeros(length(Tmtx)-1,1); end;
        
        % compute var. in bw-estimates by dividing spike times across  12 repeated trials
        %tL=sort(EvtT(EvtL==3 | EvtL==4),'ascend'); %(used to be across sample trials)
        tL=Tmtx(1):range(Tmtx)/11:Tmtx(end); % cut up time vector into 12 equal portions
        for j=1:length(tL)-1
            kk=find(ST>=tL(j) & ST<tL(j+1)); % spike indices for this trial
            if length(kk)>1, hucvS{u}(j)=UniKDE(ST(kk),h0,3);
            else hucvS{u}(j)=0; end;
        end;
        CvBW(u)=std(hucvS{u})/mean(hucvS{u});
        
        % OR, enforce common bin-width of 0.5ms (?)
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)'; 
        CvBW(u)=0;
    end;
    iFR=cell2mat(iFR0);
    
    save(['data\' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],...
                    'iFR',...  % KDE optimised firing rates (size; time vector,no.units)
                    'Tmtx',... % time vector
                    'EvtT',... % event times (samples and choices)
                    'EvtL',... % event labels (left/right, sample/choice) 
                    'avgFR',...% each neuron's average firing rate
                    'CvBW',... % variance in optimal kernel width over repeated trials
                    'hucv',... % each neuron's optimised kernel width
                    'MISE',... % mean-squared-error for each neuron's optimal kernel width
                    'hucvS');  % optimised kernel for each trial individually
end;
delete(gcp)
%% add t-info on error trials
pat='data\';
ext='';
 A=dir('*.mat')

brL={'PFC','HP'}; 
Vstr={'left_choice','right_choice','left_sample','right_sample'};
for i=1:length(A)
    
    load([A(i).name]);
    EvtTinc=[]; EvtLinc=[];
    for j=1:length(Vstr)
        f=['ERRORtrange' Vstr{j}];
        t=eval(f); EvtTinc=[EvtTinc;(t(:,1)*1e-6+5)];
        EvtLinc(end+1:end+length(t))=j*10;
    end;
    
    k=findstr('.mat',A(i).name);
    for b=1:2
        br=brL{b};
        fn=[pat A(i).name(1:k-1) '_' br '_iFR50' ext]
        save(fn,'EvtTinc','EvtLinc','-append');
    end;
    
end;



% (c) 2013 Daniel Durstewitz, Bernstein Center for Computational
% Neuroscience, CIMH/ Heidelberg University

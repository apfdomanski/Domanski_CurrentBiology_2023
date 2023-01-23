% Function for calculating smoothed estimates of firing rate by optimised kernel density smoothing. Wrapped up to run on as a single Blue Crystal job for each of HP and PFC brain areas. Uses Save-parfor to save interim data on each unit as it finishes, rather than waiting until all are finished.
%
% (c) 2014 Daniel Durstewitz, Bernstein Center for Computational Neuroscience, CIMH/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% **** run KDE_collapse_interim_results.m once finished to collate all units ****
%
% For each unit's spike train in parallel...
% (1) firing rate kernel is optimised using UniKDE:  This routine implements univariate nonparametric density 
% estimation with bandwidth selection, starting with initial guess at mean
% firing rate of that unit
% (2) Then calls a MEX'ed C routine (STtoGfAvg) to convert spike train matrix into Gaussian-filtered activity matrix
% 
% NB this routine creates a timebase for the whole experiment
% as continuous data, ignoring behavioural timestamps. 
%
% **** NB doens't contain a block of code to calculate FR for error trials ****
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

%% setup
% Converts single unit spike times to smoothed firing rates using optimised
% kernel bins
clear all

load JaroslawLONG1.mat; ext='JaroslawLONG1';

bw=50e-3;
% Vstr={'_PreSleep','_Task','_PostSleep'};
Vstr={'_PreSleep'};
trange_PreSleep  = [529434636, 3651647307];
% trange_Task       = [3793118224, 8725684155];
% trange_PostSleep = [9044732395, 1234563230];
%% Instantiate parallel pool
parpool local
poolobj = gcp;
addAttachedFiles(poolobj,{'UniKDE.m','gfilter.mexa64'})
%% process HP units
A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    
    k=findstr('_',A{i}); % collect event times & types
    
    EvtT=[]; EvtT2=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6)]; EvtT2=[EvtT2;(t(:,2)*1e-6)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT):bw:max(EvtT2); % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
%     parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
%         parfor_progress
        % collect ST
        %ST=Data(B{u})*1e-4;%  
        ST=B{u}.t*1e-4;
        avgFR(u)=length(ST)/range(ST);
        
        % convert STMtx into iFR via KDE
        h0=mean(diff(ST))^2;
        [hucv(u),MISE(u)]=UniKDE(ST,h0,3);
        if ~isnan(hucv(u)), iFR0{u}=STtoGfAvg(ST',Tmtx,sqrt(hucv(u)))';
        else iFR0{u}=zeros(length(Tmtx)-1,1); end;
        
        % compute var. in bw-estimates across trials
        %tL=sort(EvtT(EvtL==3 | EvtL==4),'ascend');
        tL=Tmtx(1):range(Tmtx)/11:Tmtx(end);    
        for j=1:length(tL)-1
            kk=find(ST>=tL(j) & ST<tL(j+1));
            if length(kk)>1, hucvS{u}(j)=UniKDE(ST(kk),h0,3);
            else hucvS{u}(j)=0; end;
        end;
        CvBW(u)=std(hucvS{u})/mean(hucvS{u});
        % OR, enforce common bin-width of 0.5ms (?)
%             iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)'; 
%             CvBW(u)=0;
        iFR_u   = iFR0{u};
        avgFR_u = avgFR(u);
        CvBW_u  = CvBW(u);
        hucv_u  = hucv(u);
        MISE_u  = MISE(u) ;
        hucvS_u = hucvS{u};
        fname= [ 'data/'  A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3) '_PreSleep_HPu_' num2str(u) 'of' num2str(length(B))];
        save_parfor(fname, iFR_u);   % KDE optimised firing rates (size; time vector,no.units)
        save_parfor(fname, Tmtx);    % time vector
        save_parfor(fname, EvtT);    % event times (samples and choices)
        save_parfor(fname, EvtL);    % event labels (left/right, sample/choice) 
        save_parfor(fname, avgFR_u); % each neuron's average firing rate
        save_parfor(fname, CvBW_u);  % variance in optimal kernel width over repeated trials
        save_parfor(fname, hucv_u);  % each neuron's optimised kernel width
        save_parfor(fname, MISE_u);  % mean-squared-error for each neuron's optimal kernel width
        save_parfor(fname, hucvS_u); % optimised kernel for each trial individually
    end;
    %     parfor_progress(0);
    iFR=cell2mat(iFR0);

    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3) '_PreSleep'],...
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
%% process PFC units
A=who('*PFCcells'); br='PFC';
for i=1:length(A) % over all datasets from this brain region
    
    k=findstr('_',A{i}); % collect event times & types
    
    EvtT=[]; EvtT2=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6)]; EvtT2=[EvtT2;(t(:,2)*1e-6)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT):bw:max(EvtT2); % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
%     parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
%         parfor_progress
        % collect ST
        %ST=Data(B{u})*1e-4;%  
        ST=B{u}.t*1e-4;
        avgFR(u)=length(ST)/range(ST);
        
        % convert STMtx into iFR via KDE
        h0=mean(diff(ST))^2;
        [hucv(u),MISE(u)]=UniKDE(ST,h0,3);
        if ~isnan(hucv(u)), iFR0{u}=STtoGfAvg(ST',Tmtx,sqrt(hucv(u)))';
        else iFR0{u}=zeros(length(Tmtx)-1,1); end;
        
        % compute var. in bw-estimates across trials
        %tL=sort(EvtT(EvtL==3 | EvtL==4),'ascend');
        tL=Tmtx(1):range(Tmtx)/11:Tmtx(end);    
        for j=1:length(tL)-1
            kk=find(ST>=tL(j) & ST<tL(j+1));
            if length(kk)>1, hucvS{u}(j)=UniKDE(ST(kk),h0,3);
            else hucvS{u}(j)=0; end;
        end;
        CvBW(u)=std(hucvS{u})/mean(hucvS{u});
        % OR, enforce common bin-width of 0.5ms (?)
%             iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)'; 
%             CvBW(u)=0;
        iFR_u   = iFR0{u};
        avgFR_u = avgFR(u);
        CvBW_u  = CvBW(u);
        hucv_u  = hucv(u);
        MISE_u  = MISE(u) ;
        hucvS_u = hucvS{u};
        fname= [ 'data/'  A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3) '_PreSleep_PFCu_' num2str(u) 'of' num2str(length(B))];
        save_parfor(fname, iFR_u);   % KDE optimised firing rates (size; time vector,no.units)
        save_parfor(fname, Tmtx);    % time vector
        save_parfor(fname, EvtT);    % event times (samples and choices)
        save_parfor(fname, EvtL);    % event labels (left/right, sample/choice) 
        save_parfor(fname, avgFR_u); % each neuron's average firing rate
        save_parfor(fname, CvBW_u);  % variance in optimal kernel width over repeated trials
        save_parfor(fname, hucv_u);  % each neuron's optimised kernel width
        save_parfor(fname, MISE_u);  % mean-squared-error for each neuron's optimal kernel width
        save_parfor(fname, hucvS_u); % optimised kernel for each trial individually
    end;
    %     parfor_progress(0);
    iFR=cell2mat(iFR0);

    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3) '_PreSleep'],...
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
%% clean up
delete(poolobj);


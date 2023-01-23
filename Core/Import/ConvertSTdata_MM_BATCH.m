% Function for calculating smoothed estimates of firing rate by optimised kernel density smoothing. 
%
% (c) 2014 Daniel Durstewitz, Bernstein Center for Computational Neuroscience, CIMH/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% A variant of ConvertSTdata_MM.m ....run sequentially for each of many
% experimental conditions.
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
% Variables (*=exported):
% avgFR* Average firing rate
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
%%
clear all
% generate iFR and event matrices
load JaroslawCP55940.mat; ext='JaroslawCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
A=who('*PFCcells'); br='PFC';
% A=who('*HPcells'); br='HP';
parpool local
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;
%%
clear all
% generate iFR and event matrices
load JaroslawCP55940.mat; ext='JaroslawCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
% A=who('*PFCcells'); br='PFC';
A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;
%%
clear all
% generate iFR and event matrices
load KrzesimirCP55940.mat; ext='KrzesimirCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
 A=who('*PFCcells'); br='PFC';
%A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;
%%
clear all
% generate iFR and event matrices
load KrzesimirCP55940.mat; ext='KrzesimirCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
% A=who('*PFCcells'); br='PFC';
A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;

%%
clear all
% generate iFR and event matrices
load KrzysztofCP55940.mat; ext='KrzysztofCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
A=who('*PFCcells'); br='PFC';
% A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;
%%
clear all
% generate iFR and event matrices
load KrzysztofCP55940.mat; ext='KrzysztofCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
% A=who('*PFCcells'); br='PFC';
A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;
%%
clear all
% generate iFR and event matrices
load MiroslawCP55940.mat; ext='MiroslawCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
A=who('*PFCcells'); br='PFC';
% A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;
%%
clear all
% generate iFR and event matrices
load MiroslawCP55940.mat; ext='MiroslawCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
% A=who('*PFCcells'); br='PFC';
A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;
%%
clear all
% generate iFR and event matrices
load NorbertCP55940.mat; ext='NorbertCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
A=who('*PFCcells'); br='PFC';
% A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;
%%
clear all
% generate iFR and event matrices
load NorbertCP55940.mat; ext='NorbertCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
% A=who('*PFCcells'); br='PFC';
A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;

%%
clear all
% generate iFR and event matrices
load OnufryCP55940.mat; ext='OnufryCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
A=who('*PFCcells'); br='PFC';
% A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;
%%
clear all
% generate iFR and event matrices
load OnufryCP55940.mat; ext='OnufryCP55940';
bw=50e-3;
Vstr={'left_choice','right_choice','left_sample','right_sample'};
% A=who('*PFCcells'); br='PFC';
A=who('*HPcells'); br='HP';
for i=1:length(A) % over all datasets from this brain region
    % collect event times & types
    k=findstr('_',A{i});
    
    EvtT=[]; EvtL=[];
    for j=1:length(Vstr) % over all choice and conditions
        f=[A{i}(1:k) 'trange' Vstr{j}];
        t=eval(f); EvtT=[EvtT;(t(:,1)*1e-6+5)];
        EvtL(end+1:end+length(t))=j;
    end;
    Tmtx=min(EvtT)-5:bw:max(EvtT)+5; % cut out (t) +/-5s

    B=eval(A{i}); ST=[]; avgFR=[];
    hucv=[]; MISE=[]; CvBW=[]; iFR0=cell(1,length(B));
    hucvS=cell(1,length(B));
    parfor_progress(length(B))
    parfor u=1:length(B)
        [i u]
        parfor_progress
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
        
        iFR0{u}=STtoGfAvg(ST',Tmtx,mean(diff(Tmtx))/100)';
        CvBW(u)=0;
    end;
    parfor_progress(0);
    iFR=cell2mat(iFR0);
    
    save(['data/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)],'iFR','Tmtx','EvtT','EvtL','avgFR','CvBW','hucv','MISE','hucvS');
end;


delete(gcp)

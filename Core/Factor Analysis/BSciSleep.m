function [FSC,FSCbs,ciLld,ciHld,ciLsc,ciHsc]=BSciSleep(TmtxS,SCM,fnA,twin)
% Calculates bootstrapped confidence intervals on loading and activation of factor scores
%
% Deals with long time-series typical of sleep recordings rather than short single-trial behavioural data.
% 
% Core part of Factor analysis Assembly detection workflow
% NB Run FAassem_.m on unit firing rate data via runFAassem_.m first.
%
% Aleksander PF Domanski PhD UoB 2015 modified from 
% 2014 Daniel Durstewitz, Bernstein Center for Computational Neuroscience, CIMH/ Heidelberg University
% aleks.domanski@bristol.ac.uk
%

      
% cut out equally-sized trial periods, kick out bad-behaved neurons
bw=mean(diff(TmtxS{1}{1}));
for br=1:2
    for i=1:length(TmtxS{br})
%         TmtxS{br}{i}=TmtxS{br}{i}([1:ko end-ko+1:end]);
        SCM{br}{i}=SCM{br}{i}';
    end;
end;
SCM{3}=cell(1,length(SCM{1}));
for i=1:length(SCM{1}), SCM{3}{i}=[SCM{1}{i};SCM{2}{i}]; end;


% Combine HP and PFC spike rasters
SCM{3}=cell(1,length(SCM{1}));
for i=1:length(SCM{1}), SCM{3}{i}=[SCM{1}{i};SCM{2}{i}]; end;

ciV=[0.1 0.05 0.01 5e-3 1e-3];
disp(['loading... ' fnA])
load(fnA);
% parpool local
for s=1:length(SCM)
    try
        na=nassem{s}(3);

    V=cell2mat(FLbs{s}{2}); vs=sort(V(1:end),'ascend');
    for i=1:length(ciV)
        ciLld(s,i)=vs(round(length(vs)*ciV(i)));
        ciHld(s,i)=vs(round(length(vs)*(1-ciV(i))));
    end;

    SCM0   = cell2mat(SCM{s});                          % [no. units x no. time points]
    SCM0   = SCM0'-ones(size(SCM0,2),1)*mean(SCM0');    % [no. time points x no. units]
    psix0  = psix{s}{na};                               % [no. units x 1]
    FL0    = FL{s}{na};                                 % [no. units x no. Assems]
    Wgt    = FL0'*(FL0*FL0'+diag(psix0))^-1;            % [no. Assems x no. units]
    FSC{s} = SCM0*Wgt';                                 % {Area no.}[no. time points x no. Assems]
    [n,...                                              % n = no. time points
     p]    = size(SCM0);                                % p = no. units
    Nbs    = length(perm{s});                           % Nbs = no. bootstrap draws
    X      = cell(Nbs,1);                               
    parfor b=1:Nbs                                      % loop over all bootstrap draws
        disp(['[Bootstraping Factor activation thresholds...  Assem area no. ' num2str(s) ' of ' num2str(length(SCM)) ' , Draw no. ' num2str(b) ' of ' num2str(Nbs) ']'])
        SCM0 = zeros(n,p);                              % [no. time points x no. units]
        for i = 1:p                                     % loop over all units
            SCMsh=cell2mat(SCM{s}(perm{s}{b}(i,:)));    % [no. units x no. time points] // perm has size {no. areas}{no. BS draws}(no. units x no events)
            SCM0(:,i)=SCMsh(i,:)';                      % [no. time points x no. units]
        end;
        SCM0=SCM0-ones(size(SCM0,1),1)*mean(SCM0);
        psix0=psixBS{s}{2}{b}; 
        FL0=FLbs{s}{2}{b};
        Wgt=FL0'*(FL0*FL0'+diag(psix0))^-1;
        Sc=SCM0*Wgt';
        X{b}=Sc(1:end);
    end;
    FSCbs{s}=cell2mat(X);
    
    ws=sort(FSCbs{s}(1:end),'ascend');
    for i=1:length(ciV)
        ciLsc(s,i)=ws(max(1,round(length(ws)*ciV(i))));
        ciHsc(s,i)=ws(max(1,round(length(ws)*(1-ciV(i)))));
    end;
    catch;
    disp(['BSciSleep>>error on area ' num2str(s)])
    end
end;
% % parpool close
% delete(gcp)
% (c) 2014 Daniel Durstewitz, Bernstein Center for Computational
% Neuroscience, CIMH/ Heidelberg University

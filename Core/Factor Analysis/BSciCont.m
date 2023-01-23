function [Ass]=BSciCont(P,D,Ass)
% Calculates bootstrapped confidence intervals on loading and activation of factor scores
% Deals with long time-series typical of sleep recordings rather than short single-trial behavioural data.
% 
% Core part of Factor analysis Assembly detection workflow
% NB Run FAassem_.m on unit firing rate data via runFAassem_.m first.
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk


% Ass = getAssemPatterns(fName,P.bw,FSC,ciHld,ciHsc,FL);


SCM = D.iFR;
      
% cut out equally-sized  periods
%     for i=1:length(D.iFRtb_)
%         SCM{i}=SCM{i}';
%     end;


ciV=[0.1 0.05 0.01 5e-3 1e-3];
% parpool local
% try
    na=Ass.nassem(3);
    V=cell2mat(Ass.FLbs{2}); vs=sort(V(1:end),'ascend');
for i=1:length(ciV)
    Ass.ciLld{1}(i)=vs(round(length(vs)*ciV(i)));
    Ass.ciHld{1}(i)=vs(round(length(vs)*(1-ciV(i))));
end;

SCM0   = cell2mat(SCM);                          % [no. units x no. time points]
SCM0   = SCM0'-ones(size(SCM0,2),1)*mean(SCM0');    % [no. time points x no. units]
psix0  = Ass.psix{na};                               % [no. units x 1]
FL0    = Ass.FL{1}{na};                                 % [no. units x no. Assems]
Wgt    = FL0'*(FL0*FL0'+diag(psix0))^-1;            % [no. Assems x no. units]
Ass.FSC{1} = SCM0*Wgt';                                 % {Area no.}[no. time points x no. Assems]
[n,...                                              % n = no. time points
 p]    = size(SCM0);                                % p = no. units
Nbs    = length(Ass.perm);                           % Nbs = no. bootstrap draws
X      = cell(Nbs,1);                               
parfor b=1:Nbs                                      % loop over all bootstrap draws
    disp(['[Bootstraping Factor activation thresholds...   Draw no. ' num2str(b) ' of ' num2str(Nbs) ']'])
    SCM0 = zeros(n,p);                              % [no. time points x no. units]
    for i = 1:p                                     % loop over all units
        SCMsh=cell2mat(SCM(Ass.perm{b}(i,:)));    % [no. units x no. time points] // perm has size {no. areas}{no. BS draws}(no. units x no events)
        SCM0(:,i)=SCMsh(i,:)';                      % [no. time points x no. units]
    end;
    SCM0=SCM0-ones(size(SCM0,1),1)*mean(SCM0);
    psix0=Ass.psixBS{2}{b}; 
    FL0=Ass.FLbs{2}{b};
    Wgt=FL0'*(FL0*FL0'+diag(psix0))^-1;
    Sc=SCM0*Wgt';
    X{b}=Sc(1:end);
end;
Ass.FSCbs=cell2mat(X);

ws=sort(Ass.FSCbs(1:end),'ascend');
for i=1:length(ciV)
    Ass.ciLsc(1,i)=ws(max(1,round(length(ws)*ciV(i))));
    Ass.ciHsc(1,i)=ws(max(1,round(length(ws)*(1-ciV(i)))));
end;


% catch;
% disp(['BSciCont>>error !'])
% end


% % parpool close
% delete(gcp)
% (c) 2014 Daniel Durstewitz, Bernstein Center for Computational
% Neuroscience, CIMH/ Heidelberg University

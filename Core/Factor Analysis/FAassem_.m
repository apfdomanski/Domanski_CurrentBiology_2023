function [nassem,FL,FSC,LL,AIC,BIC,Pr,Chi2,prBS,LLbs,Chi2bs,psix,FLbs,FSCbs,psixBS,perm]=FAassem_(SCM,kmax,Nbs,alpha)
% Factor analysis assembly detection algorithm
%
% Original 2014 Daniel Durstewitz, Bernstein Center for Computational Neuroscience, CIMH/ Heidelberg University
% Updates by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Uses MATLAB's factoran function. Can specify maximum mdoel complexity (i.e. expected no. factor) to test up to as kmax.
% Reports various measures of confidence on parameters and model
% complexity. Time series is also cut up and bootstrapped into factoran to
% provide non-parametric CIs.
%
% Input variables:
% SCM       cell array of equally sized p x N spike count matrices 
% kmax      number of potential factors to evlauate up to
% Nbs       number of bootstrap draws to perform
% alpha     significance threshold for Chi^2 stat
%
% Output variables:
% nassem       4 metrics of number of detected assemblies, order: AIC, BIC, BS, pChi2
% ...for each value of potential no. factors (counter = k):
% FL{k}        Factor loading matrix
% FSC{k}       Factor scores (predictions of common factors)    
% psix{k}      ML estimates of psi(specific variances)
% LL(k)        Log-likelihood ratio
% AIC(k)       Akaike Information criterion
% BIC(k)       Bayesian Information criterion
% Pr(k)        Right tail p-value for null hypothesis 
% Chi2(k)      Ch^2 stat for null hopethesis
% ... as above, bootstrapped (bootstrap trial counter = b) versions for each number of k's:
% prBS(k)
% LLbs(b,k)
% Chi2bs(b,k)
% FLb{b}{k}
% FSCbs(k)
% psixBS(k)
% perm         Bootstrap-drawn spike times (shuffled neuron-trial assignments)


if iscell(SCM)
    SCM0=cell2mat(SCM);
else
    SCM0=SCM';
end
[p,N]=size(SCM0); % p=no. units, N=no. data points

s=floor((2*p+1-sqrt(8*p+1))/2); % optional default no. params to check up to

if nargin<2 || isempty(kmax), kmax=s;
else kmax=min(kmax,s); end;
if nargin<3 || isempty(Nbs), Nbs=1000; end;
if nargin<4 || isempty(alpha), alpha=0.01; end;
nrep=0;
Pr=zeros(1,kmax)+nan; AIC=Pr; BIC=Pr; FL=cell(1,kmax); Chi2=Pr; LL=Pr; psix=FL; FSC=FL;
for k=1:kmax
    disp(['nfac = ' num2str(k)]);
    if nrep>0, 
        [FL{k},psix{k},~,stats,FSC{k}]=factoran(SCM0',k,'maxit',10000,'delta',1e-4,'start',nrep);
    else
        [FL{k},psix{k},~,stats,FSC{k}]=factoran(SCM0',k,'maxit',10000,'delta',1e-4); 
    end;
    if isfield(stats,'p') 
        Pr(k)   = stats.p;
        Chi2(k) = stats.chisq; 
    end;
    Cpred   = FL{k}*FL{k}' + diag(psix{k});
    LL(k)   = -N/2*(p*log(2*pi) + log(det(Cpred)) + trace(Cpred^-1*corr(SCM0'))); % log-likelihood ratio
    Chi2(k) = (N-(2*p+11)/6-2*k/3)*(log(det(Cpred)) - log(det(corr(SCM0'))));     % Chi^2
    np      = k*p+p-k*(k-1)/2;                                                    % number of parameters
    AIC(k)  = -2*LL(k)+2*np;
    BIC(k)  = -2*LL(k)+log(N)*np;
end;


%order: AIC, BIC, BS, pChi2
nassem=zeros(1,4);
try
    [~,nassem(1)]=min(AIC);
catch
    nassem(1)=NaN;
    disp('caution; no AIC reported! ...No Assemblies detected?!')
end
try
    [~,nassem(2)]=min(BIC);
catch
    nassem(2)=NaN;
    disp('caution; no BIC reported! ...No Assemblies detected?!')

end
r=find(Pr>alpha,1);
if ~isempty(r), nassem(4)=r-1; else nassem(4)=nan; end;

% Bootstrap distribution (shuffle neuron-trial assignments)
kmax=2;
Chi2bs=zeros(Nbs,kmax)+nan; LLbs=Chi2bs; perm=cell(Nbs,1);
for b=1:Nbs
    disp(['BS# ' num2str(b)]);
    SCM0bs=zeros(p,N);
    for iUnit=1:p
        perm{b}(iUnit,:)=randperm(length(SCM));
        SCMsh=cell2mat(SCM(perm{b}(iUnit,:)));
        SCM0bs(iUnit,:)=SCMsh(iUnit,:);
    end;
    
    parfor k=1:kmax
        if nrep>0, [FLbs0,psixBS0,~,~,FSCbs0]=factoran(SCM0bs',k,'maxit',10000,'delta',1e-4,'start',nrep);
        else [FLbs0,psixBS0,~,~,FSCbs0]=factoran(SCM0bs',k,'maxit',10000,'delta',1e-4); end;
        Cpred=FLbs0*FLbs0'+diag(psixBS0);
        LLbs(b,k)=-N/2*(p*log(2*pi)+log(det(Cpred))+trace(Cpred^-1*corr(SCM0bs')));
        FLbs{k}{b}=FLbs0; FSCbs{k}{b}=FSCbs0; psixBS{k}{b}=psixBS0;
        Chi2bs(b,k)=(N-(2*p+11)/6-2*k/3)*(log(det(Cpred))-log(det(corr(SCM0bs'))));
    end;
end;
z=diff(LL); Zbs=diff(LLbs')';
prBS=zeros(1,length(z))+nan;
for i=1:length(z)
    k=find(Zbs(:,min(i,kmax-1))>=z(i));
    prBS(i)=length(k)/length(Zbs(:,min(i,kmax-1)));
end;
k=find(prBS>alpha,1);
if ~isempty(k), nassem(3)=k; else nassem(3)=1; end;


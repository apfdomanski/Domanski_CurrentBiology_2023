function FracCorr=DecodePredErrCrossTemporal(FR,evt0,reg)
% CVE decoding on erroroneous behavioural trials
%
% (c) 2014, Daniel Durstewitz, Dept. Theoretical Neuroscience, CIMH
% Mannheim/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Input Variables:
% FR:   Firing rate matrix [units,time x no. trials]
% evt0: Event category labels (either 1 or 2) [no. trials]
% reg:  Regularization parameter
%
% Output Variables:
% PE:  Decoding result 
%
if nargin<3, reg=0; end;

Ltr=round(size(FR,1)/length(evt0));
PE=zeros(Ltr,Ltr); N=PE;

parfor t=1:Ltr % Loop across each time bin
    X=FR(t:Ltr:end,:);
    Xtrain1=X(evt0==1,:); k1= ~isnan(Xtrain1(:,1)); Xtrain1=Xtrain1(k1,:); % Firing during correct responses to class 1
    Xtrain2=X(evt0==2,:); k2= ~isnan(Xtrain2(:,1)); Xtrain2=Xtrain2(k2,:); % Firing during correct responses to class 2
    n1=length(k1); n2=length(k2);
    Cpool=((n1-1)*cov(Xtrain1)+(n2-1)*cov(Xtrain2))./(n1+n2-2);            % Pooled covariance matrix of 2 response sets
    Creg=Cpool+reg*eye(size(Cpool));                                       % Regularized (pooled) covariance matrix
    m1=mean(Xtrain1); m2=mean(Xtrain2);                                    % Mean at this time point
    for t_=1:Ltr % Loop across each time bin
       
        X_ = FR(t_:Ltr:end,:);
        for s=1:2 % Loop across each outcome type
            
            Xtest=X_(evt0==s,:); k= ~isnan(Xtest(:,1)); Xtest=Xtest(k,:);
            n=length(k);
            D1=diag((Xtest-ones(n,1)*m1)*Creg^-1*(Xtest-ones(n,1)*m1)');       % Mah. distance from response 1 mean at this time point
            D2=diag((Xtest-ones(n,1)*m2)*Creg^-1*(Xtest-ones(n,1)*m2)');       % Mah. distance from response 2 mean at this time point
            cpred=(D1>D2)+1;                                                   % Predicted class at this time point
            PE(t,t_)=PE(t,t_)+sum(cpred~=s);                                   % Classifcation error at this time point
            N(t,t_)=N(t,t_)+n;

            
        end
      
    end
end
PE=PE./N;                                                                  % Normalise by length
FracCorr=1-PE;

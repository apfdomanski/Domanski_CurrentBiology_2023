function CVE=DecodeCVE(FR,evt0,reg)
% Leave-one out cross-validation error (CVE) Decoder 
%
% (c) 2014, Daniel Durstewitz, Dept. Theoretical Neuroscience, CIMH
% Mannheim/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Reports binary classification strength between two outcomes across time
%
% Input Variables:
% FR:   Firing rate matrix [time x no. trials]
% evt0: Event category labels (either 1 or 2) [no. trials]
% reg:  Regularization parameter
%
% Output Variables:
% CVE:  Decoding strength 
%
if nargin<3, reg=0; end;
% b=waitbar(0,'Decoding...');

ntr=length(evt0);
Ltr=round(size(FR,1)/ntr);
CVE=zeros(1,Ltr);
for t=1:Ltr % Loop across each time bin
    %     waitbar(t/Ltr,b,sprintf('Decoding timestep: %d / %d\n',t,Ltr))
    X=FR(t:Ltr:end,:);
    if sum(sum(X))==0
        CVE(t) = ntr;
    else
        for i=1:ntr % Loop across each trial
            evt0i=evt0([1:i-1 i+1:end]);                            % Event outcomes (leave this one out)
            Xi=X([1:i-1 i+1:end],:);                                % Firing rates (leave this one out)
            X1=Xi(evt0i==1,:); k1= ~isnan(X1(:,1)); X1=X1(k1,:);    % Sort responses into outcome 1 or 2
            X2=Xi(evt0i==2,:); k2= ~isnan(X2(:,1)); X2=X2(k2,:);
            n1=length(k1); n2=length(k2);                           % No of each type of response
            Cpool=((n1-1)*cov(X1)+(n2-1)*cov(X2))./(n1+n2-2);       % Pooled covariance matrix
            Creg=Cpool+reg*eye(size(Cpool));                        % Regularized (pooled) covariance matrix
            m1=mean(X1); m2=mean(X2);                               % Mean at each time point
            d1=(X(i,:)-m1)*Creg^-1*(X(i,:)-m1)';                    % Mah. distance from class 1 mean at this time point
            d2=(X(i,:)-m2)*Creg^-1*(X(i,:)-m2)';                    % Mah. distance from class 2 mean at this time point
            cpred=(d1>d2)+1;                                        % Predicted class at this time point
            CVE(t)=CVE(t)+(cpred~=evt0(i));                         % Cross-validation error at this time point
        end;
    end
end;
CVE=CVE./ntr;                                                   % normalise byt no. trials



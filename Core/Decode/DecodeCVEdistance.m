function [CVE, pClass1,CVEunits,pClassunits] =DecodeCVEdistance(FR,evt0,reg)
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

ntr=length(evt0);
Ltr=round(size(FR,1)/ntr);
CVE=zeros(1,Ltr);
printf('Decoding group...')
for t=1:Ltr % Loop across each time bin
    X=FR(t:Ltr:end,:); %[no. Trials,no. Units]
    for i=1:ntr % Loop across each trial

        evt0i=evt0([1:i-1 i+1:end]);                            % Event outcomes (leave this one out)
        Xi=X([1:i-1 i+1:end],:);                                % Firing rates (leave this one out)
        X1=Xi(evt0i==1,:); k1= ~isnan(X1(:,1)); X1=X1(k1,:);    % Firing rates outcome 1 [N(k=1)-this trial, No. units]
        X2=Xi(evt0i==2,:); k2= ~isnan(X2(:,1)); X2=X2(k2,:);    % Firing rates outcome 2 [N(k=2)-this trial, No. units]
        n1=length(k1); n2=length(k2);                           % No of each type of response
        Cpool=((n1-1)*cov(X1)+(n2-1)*cov(X2))./(n1+n2-2);       % Pooled covariance matrix
        Creg=Cpool+reg*eye(size(Cpool));                        % Regularized (pooled) covariance matrix
        m1=mean(X1); m2=mean(X2);                               % Mean at this time point
        d1=(X(i,:)-m1)*Creg^-1*(X(i,:)-m1)';                    % Mah. distance from class 1 mean at this time point
        d2=(X(i,:)-m2)*Creg^-1*(X(i,:)-m2)';                    % Mah. distance from class 2 mean at this time point
        cpred=(d1>d2)+1;                                        % Predicted class at this time point
        CVE(t)=CVE(t)+(cpred~=evt0(i));                         % Cross-validation error at this time point
        
        pClass1(t,i)= 1-(d1/(d1+d2));                           % p(class 1) at time=t for trial i
        %pClass2(t,i)= 1-(d2/(d1+d2));                           % p(class 2) at time=t for trial i
        
    end;
end;
CVE=CVE./ntr;                                                   % normalise by no. trials



% pClass1_mean = nanmean(pClass1(:,evt0==1),2);
% pClass2_mean = nanmean(pClass2(:,evt0==2),2);
% pClass1_SEM  = nansem(pClass1(:,evt0==1),2);
% pClass2_SEM  = nansem(pClass2(:,evt0==2),2);


% figure; 
% subplot(2,1,1);hold on
% ciplot(pClass1_mean+pClass1_SEM,...
%        pClass1_mean-pClass1_SEM,...
%        1:t,'b')
% ciplot(pClass2_mean+pClass2_SEM,...
%        pClass2_mean-pClass2_SEM,...
%        1:t,'r')   
% plot(1:t,pClass1_mean,'b','LineWidth',1.5)
% plot(1:t,pClass2_mean,'r','LineWidth',1.5)
% legend('Left trials','Right trials')
% plot([0 t],[0.5 0.5],':k')
% axis([0 t 0 1])
% ylabel('p(Correct trial classification)')
% 
% subplot(2,1,2);hold on
% plot(1-CVE)
% plot([0 t],[0.5 0.5],':k')
% axis([0 t 0 1])
% ylabel('Decoder accuracy')




CVEunits = zeros(size(FR,2),Ltr);
pClassunits = cell(size(FR,2),1);
for iUnit =1:size(FR,2)
    fprintf('decoding unit %d of %d\n',iUnit,size(FR,2))
    for t=1:Ltr % Loop across each time bin
        X=FR(t:Ltr:end,iUnit); %[no. Trials,no. Units]
        for i=1:ntr % Loop across each trial
            
            evt0i=evt0([1:i-1 i+1:end]);                            % Event outcomes (leave this one out)
            Xi=X([1:i-1 i+1:end],:);                                % Firing rates (leave this one out)
            X1=Xi(evt0i==1,:); k1= ~isnan(X1(:,1)); X1=X1(k1,:);    % Firing rates outcome 1 [N(k=1)-this trial, No. units]
            X2=Xi(evt0i==2,:); k2= ~isnan(X2(:,1)); X2=X2(k2,:);    % Firing rates outcome 2 [N(k=2)-this trial, No. units]
            n1=length(k1); n2=length(k2);                           % No of each type of response
            Cpool=((n1-1)*cov(X1)+(n2-1)*cov(X2))./(n1+n2-2);       % Pooled covariance matrix
            Creg=Cpool+reg*eye(size(Cpool));                        % Regularized (pooled) covariance matrix
            m1=mean(X1); m2=mean(X2);                               % Mean at this time point
            d1=(X(i,:)-m1)*Creg^-1*(X(i,:)-m1)';                    % Mah. distance from class 1 mean at this time point
            d2=(X(i,:)-m2)*Creg^-1*(X(i,:)-m2)';                    % Mah. distance from class 2 mean at this time point
            cpred=(d1>d2)+1;                                        % Predicted class at this time point
            CVEunits(iUnit,t)=CVEunits(iUnit,t)+(cpred~=evt0(i));   % Cross-validation error at this time point
            
            pClassunits{iUnit}(t,i)= 1-(d1/(d1+d2));                % p(class 1) at time=t for trial i
            %pClass2(t,i)= 1-(d2/(d1+d2));                           % p(class 1) at time=t for trial i
            
        end;
    end;
    CVEunits(iUnit,:)=CVEunits(iUnit,:)./ntr;
end



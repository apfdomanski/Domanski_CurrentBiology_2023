function [TS,avgFR,sqFR]=DecodeStatsTSonly(FR,evt0,reg)
% Fisher F-test based classifier/decoder
%
% (c) 2014, Daniel Durstewitz, Dept. Theoretical Neuroscience, CIMH
% Mannheim/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Input Variables:
% FR:   Firing rate matrix [time x no. Units]
% evt0: Event category labels (either 1 or 2) [no. trials]
% reg:  Regularization parameter
%
% Output Variables:
% avgFR     {trial type}(time , unit no.) Average firing rate for each event type
% seFR      {trial type}(time , unit no.) Firing rate squared error for each event type
% Ft2x      Hotelling's T^2 statistic: discriminabiltiy of 2 classes
% Rt2x      Adjusted T^2
% Ft2ciL0   Low 95% confidence limits on T^2
% Ft2ciH0   High 95% confidence limits on T^2
% TS        T-score
% dfnum     F-distribution numerator DofF
% dfd       F-distribution denominator DofF

if nargin<3, reg=0; end;

ntr=length(evt0);
Ltr=round(size(FR,1)/ntr);
nu=size(FR,2);
% s is the trial outcome, either 1 or 2
avgFR{1}=zeros(Ltr,nu); avgFR{2}=avgFR{1};
sqFR{1}=zeros(Ltr,nu); sqFR{2}=sqFR{1};
N{1}=zeros(Ltr,nu); N{2}=N{1};

for i=1:ntr
    s=evt0(i);
    FR0=FR((i-1)*Ltr+(1:Ltr),:);
    k=find(~isnan(FR0));
    FR0(isnan(FR0))=0;
    avgFR{s}=FR0+avgFR{s};
    sqFR{s}=FR0.^2+sqFR{s};
    N{s}(k)=N{s}(k)+1;
end;

for s=1:2
    avgFR{s}=avgFR{s}./N{s};
    seFR{s}=(sqFR{s}./N{s}-avgFR{s}.^2)./sqrt(N{s});
end;

% if nargout>2
%     dfnum=size(FR,2);
%     Ft2=[]; Ft2ciL=[]; Ft2ciH=[]; dfden=[]; Rt2=[]; Rt2b=[];
%     for i=1:Ltr
%         i
%         X=FR(i:Ltr:end,:);
%         X1=X(evt0==1,:);
%         k1= ~isnan(X1(:,1)); X1=X1(k1,:);
%         X2=X(evt0==2,:);
%         k2= ~isnan(X2(:,1)); X2=X2(k2,:);
%         n1=length(k1); n2=length(k2);
%         Cpool=((n1-1)*cov(X1)+(n2-1)*cov(X2))./(n1+n2-2);
%         Creg=Cpool+reg*eye(size(Cpool));
%         U=sqrtm(Cpool+1e-8*eye(size(Cpool)));
%         dfeff=real(trace(U*Creg^-1*U'));
%         dfden(i)=n1+n2-dfeff-1;
%         dm=avgFR{1}(i,:)-avgFR{2}(i,:);
%         HT2=(dm*Creg^-1*dm')*n1*n2/(n1+n2);
%         Ft2(i)=HT2*dfden(i)/((n1+n2-2)*dfeff);
%         if nargout>4
%             Ft2ciL(i)=finv(0.05,dfnum,dfden(i));
%             Ft2ciH(i)=finv(0.95,dfnum,dfden(i));
%         end;
%         Rt2(i)=Ft2(i)/(Ft2(i)+dfden(i));
%     end;
% end;

    md=abs(avgFR{1}-avgFR{2});
    ntr1=length(find(evt0==1)); ntr2=length(find(evt0==2));
    v1=seFR{1}*sqrt(ntr1)*ntr1;
    v2=seFR{2}*sqrt(ntr2)*ntr2;
    spool=sqrt((v1+v2)./(ntr1+ntr2-2));
    TS=md./(spool*sqrt(1/ntr1+1/ntr2));


% (c) 2014, Daniel Durstewitz, Dept. Theoretical Neuroscience, CIMH
% Mannheim/ Heidelberg University

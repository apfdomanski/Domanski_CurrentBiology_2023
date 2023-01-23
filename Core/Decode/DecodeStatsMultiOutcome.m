function [avgFR,seFR,Ft2,Rt2,Ft2ciL,Ft2ciH]=DecodeStatsMultiOutcome(FR,reg)
% Fisher F-test based classifier/decoder for multi-way comparisons (>2 Outcomes)
% Aleksander PF Domanski PhD UoB 2017
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

CompOption = 3;
% 1: average of one-vs-one comparisons
% 2: average of one-vs-all comparisons
% 3: class mean vs total mean
if nargin<2, reg=0; end;
% evt0 = 1:length(FR);
nu   = length(FR{1});
Ltr  = length(FR{1}{1});

% s is the trial outcome, either 1 or 2
for s=1:length(FR)
    avgFR{s} = zeros(Ltr,nu);
    sqFR{s}  = zeros(Ltr,nu);
    N{s}     = zeros(Ltr,nu); 
end
for s = 1:length(FR)
    for iUnit = 1:nu
        for iPos = 1:Ltr
            FR_ = FR{s}{iUnit}{iPos};
            avgFR_ = nanmean(FR_);
            N_   = numel(FR_(~isnan(FR_)));
            
            sqFR_ = nansum(FR_.^2);
            seFR_ = (sqFR_./N_ - avgFR_.^0.5)./N_^0.5;
            
            seFR{s}(iPos,iUnit)  = seFR_;
            avgFR{s}(iPos,iUnit) = avgFR_;
        end
    end
end
   
if nargout>2
    dfnum=nu;
    Ft2=[]; Ft2ciL=[]; Ft2ciH=[]; dfden=[]; Rt2=[]; Rt2b=[];
    for i=1:Ltr
        X = cell(length(FR),1);
        for s = 1:length(FR)
            for iUnit=1:nu
                X{s}(:,iUnit) =  FR{s}{iUnit}{i};
            end
            k= ~isnan(X{s}(:,1)); X{s} =  X{s}(k,:);
            n(s) = sum(k);
        end
        Cpool = zeros(nu,nu);
        for s= 1:length(FR)
            Cpool = Cpool + (n(s)-1)*cov(X{s});
        end
        Cpool=Cpool./(sum(n)-length(FR));
        Creg=Cpool+reg*eye(size(Cpool));
        U=sqrtm(Cpool+1e-8*eye(size(Cpool)));
        dfeff=real(trace(U*Creg^-1*U'));
        dfden(i)=sum(n)-dfeff - length(FR); % Check this... -no classes or -1
        
       
        switch CompOption
            case 1  % calculate difference in means as average of one vs. one
                
                dm = zeros(1,nu);
                for i_ = 1:length(FR)
                    for j_ = 1:length(FR)
                        if i_ ~= j_
                            dm = dm+avgFR{j_}(i,:) - avgFR{i_}(i,:);
                        end
                    end
                end
                dm = dm./(length(FR)^2 - length(FR));
                
            case 2 % calculate difference in means as average of one vs. all
                
                dm = zeros(1,nu);
                for i_ = 1:length(FR)
                    m_ = [];
                    for j_ = setdiff(1:length(FR),i_)
                        m_(j_,:) = avgFR{j_}(i,:);
                    end
                    dm = dm + avgFR{i_}(i,:) - mean(m_(setdiff(1:length(FR),i_),:));
                end
                dm = dm/length(FR);
                
            case 3 % class mean vs. total mean
                dm = zeros(1,nu);
                m_ = zeros(length(FR),nu);
                for i_ = 1:length(FR)
                    m_(i_,:) = avgFR{i_}(i,:);
                end
                
            for i_ = 1:length(FR)
                    dm(i_,:) = avgFR{i_}(i,:) - mean(m_);
            end
            dm = mean(dm);
            
        end
        HT2=(dm*Creg^-1*dm')*prod(n)/(sum(n));
        Ft2(i)=HT2*dfden(i)/((sum(n)-length(FR))*dfeff);
        Rt2(i)=Ft2(i)/(Ft2(i)+dfden(i));
        if nargout>4
            Ft2ciL(i)=finv(0.05,dfnum,dfden(i));
            Ft2ciH(i)=finv(0.95,dfnum,dfden(i));
        end;
    end;
end;


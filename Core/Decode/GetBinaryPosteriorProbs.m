function [Accuracy, pCorrect,Accuracyunits,pCorrectunits] = GetBinaryPosteriorProbs(FR,evt0)
% Leave-one out cross-validation error (CVE) Decoder - linear decision boundaries

% Aleksander Domanski PhD 2016 UoB

% Input Variables:
% FR:   Firing rate matrix [trial time x no. trials,noUnits]
% evt0: Event category labels (either 1 or 2) [no. trials]

% Output Variables:
% Accuracy: Fraction of decoders (leave-one-out) that are correct for each time-step
% pCorrect: Posterior probability that each decoder is correct on a trial-by-trial basis

b=waitbar(0,'Decoding group...');
% setappdata(b,'canceling',0)
    
if nargout<3
    DecodeUnivariate=false;
else
    DecodeUnivariate=true;
end
ntr=length(evt0);
Ltr=round(size(FR,1)/ntr);
Accuracy=zeros(1,Ltr);
for t=1:Ltr 
	
    waitbar(t/Ltr,b,sprintf('Decoding timestep: %d / %d\n',t,Ltr))       
    X=FR(t:Ltr:end,:); %[no. Trials,no. Units]                  % Data for classification 
    for i=1:ntr % Loop across each trial
%         if getappdata(b,'canceling')
%             break
%         end
        evt0i=evt0([1:i-1 i+1:end]);                            % Event outcomes (leave this one out)
        Xi=X([1:i-1 i+1:end],:);                                % Firing rates (leave this one out)
        try
        [class,~,Posterior] = classify(X(i,:),Xi,evt0i');       % Linear decision boundary classifier
        
        pCorrect(t,i)= Posterior(evt0(i));                      % Posterior probability (pCorrect) at this time for this trial
        Accuracy(t) = Accuracy(t) + (class==evt0(i));           % Accumulate decoder accuracy
        catch
            pCorrect(t,i) = NaN;
            Accuracy(t)   = NaN;
        end
    end;
end;
Accuracy = Accuracy./ntr;

if DecodeUnivariate
    
    waitbar(0,b,'Decoding units...')
    Accuracyunits = zeros(size(FR,2),Ltr);
    pCorrectunits = cell(size(FR,2),1);
    for iUnit =1:size(FR,2)
        for t=1:Ltr % Loop across each time bin    
            waitbar(t/Ltr,b,sprintf('Unit %d / %d, decoding timestep: %d / %d\n',iUnit,size(FR,2),t,Ltr)  )       
            X=FR(t:Ltr:end,iUnit); %[no. Trials,no. Units]
            for i=1:ntr % Loop across each trial
    %             if getappdata(b,'canceling')
    %                 break
    %             end
                evt0i=evt0([1:i-1 i+1:end]);                            % Event outcomes (leave this one out)
                Xi=X([1:i-1 i+1:end],:);                                % Firing rates (leave this one out)

                [class,~,Posterior] = classify(X(i,:),Xi,evt0i');       % Linear decision boundary classifier

                pCorrectunits{iUnit}(t,i)= Posterior(evt0(i));                      % Posterior probability (pCorrect) at this time for this trial
                Accuracyunits(iUnit,t) = Accuracyunits(iUnit,t) + (class==evt0(i));           % Accumulate decoder accuracy


            end;
        end;
        Accuracyunits(iUnit,:)=Accuracyunits(iUnit,:)./ntr;
    end
end
delete(b)


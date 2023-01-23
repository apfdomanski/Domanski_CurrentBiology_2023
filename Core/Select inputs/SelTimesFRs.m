function [Tcutout,cutout]=SelTimesFRs(twin,Tmtx,data,tRange)
% Select specified time range around trial periods from continuous data in
% workspace
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk

if length(twin)==1;twin =[-twin twin]; end

for s=1:length(data) % across brain regions: PFC, HP, joint
    Tmtx=Tmtx(1:end-1); 
    cutout{s}=cell(1,length(tRange)); Tcutout{s}=cutout{s};
    
    for iEvent=1:length(tRange)
        k=find(Tmtx>=tRange(iEvent)+twin(1) & Tmtx<=tRange(iEvent)+twin(2));
        cutout{s}{iEvent} = (data{s}.iFR(k,:));
        Tcutout{s}{iEvent}= Tmtx(k);
    end;
end;

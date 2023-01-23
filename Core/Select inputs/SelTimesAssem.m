function [Tcutout,cutout]=SelTimesAssem(twin,Tmtx,data,tRange)
% Select specified time range around trial periods from continuous data in
% workspace
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk

if size(data{1},1)~=length(Tmtx)
        Tmtx=Tmtx(1:end-1); 
end
for s=1:length(data) % across brain regions: PFC, HP, joint
    
    cutout{s}=cell(1,length(tRange)); Tcutout{s}=cutout{s};
    for iEvent=1:length(tRange)-1
        k=find(Tmtx>=tRange(iEvent)-twin & Tmtx<=tRange(iEvent)+twin);
        cutout{s}{iEvent} = data{s}(k,:);
        Tcutout{s}{iEvent}= Tmtx(k);
    end
end

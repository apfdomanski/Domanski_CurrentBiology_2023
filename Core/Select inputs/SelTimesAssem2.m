function [cutout,Tcutout]=SelTimesAssem2(twin,Tmtx,data,tRange,bw)
% Select specified time range around trial periods from continuous data in
% workspace
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk
nTrial = size(tRange,1);
LTr = length([-twin:bw:twin]);
if size(data{1},1)~=length(Tmtx)
    Tmtx=Tmtx(1:end-1);
end
for s=1:length(data) % across brain regions: PFC, HP, joint
    nAss = size(data{s},2);
    cutout{s}=cell(nAss,1); Tcutout=zeros(nTrial,2*LTr);
    for iEvent=1:nTrial
        idx = [FindClosestIndex(Tmtx,tRange(iEvent,1)-twin);FindClosestIndex(Tmtx,tRange(iEvent,2)-twin)];
        k = [idx(1):(idx(1)+LTr-1),idx(2):(idx(2)+LTr-1)];
        for iAss = 1:nAss
            cutout{s}{iAss}  = [cutout{s}{iAss},data{s}(k,iAss)];
            Tcutout(iEvent,1:2*LTr) = Tmtx(k);
        end
    end
end

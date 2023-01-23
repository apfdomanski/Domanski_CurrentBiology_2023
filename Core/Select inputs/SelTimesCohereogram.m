function [Tcutout,cutout]=SelTimesCohereogram(twin,Tmtx,data,tRange)
% Select specified time range around trial periods from continuous data in
% workspace
%
% Aleksander PF Domanski PhD UoB 2016
% aleks.domanski@bristol.ac.uk


% Tmtx=Tmtx(1:end-1); 
cutout=cell(1,length(data)); Tcutout=cutout;
%len_ = [];
for iPair = 1:length(data)
    for iEvent=1:2:length(tRange)-1
        k=find(Tmtx{iPair}>=tRange(iEvent)-twin & Tmtx{iPair}<=tRange(iEvent+1)+twin);
        cutout{iPair}{ceil(iEvent/2)} = data{iPair}(k,:);
        Tcutout{iPair}{ceil(iEvent/2)}= Tmtx{iPair}(k);%-tRange(iEvent);
        %len_=[len_ , cellfun(@length, Tcutout{iPair})];
    end
end

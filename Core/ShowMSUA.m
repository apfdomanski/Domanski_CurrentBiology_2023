function [evtB,Clr]=ShowMSUA(Tmtx,iFR,CL,iFRred,ret,TrCol,sc)
% Plot a heatmap timeseries showing the firing rates of all single units 
%
% (c) 2014 Daniel Durstewitz, Bernstein Center for Computational Neuroscience, CIMH/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
if nargin<4, iFRred=[]; end;
if nargin<5, ret=0; end;
if nargin<6 || isempty(TrCol), TrCol='k'; end;
if nargin<7, sc=0; end;

Font_Size = 12;


% scale all FR to [0,1]
if sc
    V=ones(size(iFR,1),1);
    iFR=(iFR-V*min(iFR))./(V*range(iFR));
end;

cmap=colormap;
for i=1:max(CL), Clr{i}=cmap(35+(i-1)*20,:); end;
cl0=[0 CL 0]; evtB=find(diff(cl0))-1;
k=~(CL(evtB(evtB>0))==0 | CL(evtB+1)==0);
evtB=sort([evtB evtB(k)]);

% plot spike density mtx
if ~isempty(Tmtx)
    figure(ret+2), hold off cla;
    imagesc(Tmtx,1:size(iFR,2),iFR');
    colormap('jet'); colorbar;
    xlabel('Time (s)','FontSize',Font_Size); ylabel('unit #','FontSize',Font_Size);
    set(gca,'FontSize',Font_Size);
    hold on
    dt=mean(diff(Tmtx));
    for i=1:2:length(evtB)
        k0=evtB(i)+1; k1=evtB(i+1); 
        c=CL(k0);
        plot([Tmtx(k0)-dt/2 Tmtx(k0)-dt/2],[0 size(iFR,2)+1],'Color',Clr{c},'LineWidth',3);
        plot([Tmtx(k1)+dt/2 Tmtx(k1)+dt/2],[0 size(iFR,2)+1],'Color',Clr{c},'LineWidth',3);
        plot([Tmtx(k0)-dt/2 Tmtx(k1)+dt/2],[0.5 0.5],'Color',Clr{c},'LineWidth',3);
        plot([Tmtx(k0)-dt/2 Tmtx(k1)+dt/2],[size(iFR,2)+0.5 size(iFR,2)+0.5],'Color',Clr{c},'LineWidth',3);
    end;
    ylim([0 size(iFR,2)+1]);
end;

% plot 3D MSUA space representation    
if ~isempty(iFRred)
    s=size(iFR,1)-size(iFRred,1);
    figure(1); if ret==0, hold off cla; end;
    plot3(iFRred(:,1),iFRred(:,2),iFRred(:,3),'Color',TrCol,'LineWidth',3), hold on
    for i=1:2:length(evtB)
        k0=evtB(i)+1; k1=evtB(i+1); 
        c=CL(k0);
        r=max(1,k0-s):k1-s;
        if r(end)>0
            plot3(iFRred(r,1),iFRred(r,2),iFRred(r,3),'Color',Clr{c},'LineWidth',3)
            hold on
        end;
    end;
end;

% (c) 2013 Daniel Durstewitz, Bernstein Center for Computational
% Neuroscience, CIMH/ Heidelberg University

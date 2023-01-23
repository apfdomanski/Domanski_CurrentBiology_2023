function [ActT,Za]=ShowAssemByTrialG(Za,...
                                     Tmtx,...
                                     STPat,...
                                     EvtT,EvtL,...
                                     Lag,crit,show)
% Plot assembly activation during a trial
%
% (c) 2014, Daniel Durstewitz, Dept. Theoretical Neuroscience, CIMH
% Mannheim/ Heidelberg University
                                 
if nargin<6 || isempty(Lag), Lag=zeros(size(Za,2)); end;
if nargin<7 || isempty(crit), crit=1; end;
if nargin<8 || isempty(show), show=1; end;

L=cellfun(@length,STPat);
k= L>1; STPat=STPat(k);
% bw=mean(diff(Tmtx{1}(1:100)));
bw=mean(diff(Tmtx(1:100)));

if show
    clr=colormap('gray');
    clr2=colormap('hot');
    clr(30,:)=clr2(28,:); %red
    clr2=colormap('summer');
    clr(32,:)=clr2(45,:);   %green
    clr(31,:)=clr2(end,:);  %yellow
    clr2=colormap('cool'); 
    clr(33,:)=clr2(10,:);    %light blue
    clr(34,:)=clr2(1,:);    %purple
    clr2=colormap('jet');
    clr(35,:)=clr2(1,:);    %dark blue
    ccL={'w','y','b','g','r','c','m','k'};
end;
ActT=cell(1,length(Tmtx));
for i=1:length(Tmtx)
    hold off cla
    Z0=28*Za{i};
    [N,p]=size(Z0);
    for j=1:length(STPat) %j=unit_id?
        r=mod(j-1,6)+1;
        
        L0=Lag(STPat{j},STPat{j}); %this_lag
        s0=abs(min(L0(1:end))); sx=max(L0(1:end));
        X=zeros(N+sx+s0,length(STPat{j}));
        u1=STPat{j}(1);
        X(s0+1:s0+N,1)=Za{i}(:,u1);
        for k2=2:length(STPat{j})
            u2=STPat{j}(k2);
            l=Lag(u1,u2);
            X(s0+1+l:s0+N+l,k2)=Za{i}(:,u2);
        end;
        xf=sum(X(s0+1:s0+N,:)')./length(STPat{j}); w=find(xf>=crit);
        for k2=1:length(STPat{j})
            u2=STPat{j}(k2);
            w2=w+Lag(u2,u1); w2=w2(w2>0 & w2<N+1);
            Z0(w2,u2)=(29+r)*sign(Z0(w2,u2));
        end;
        ActT{i}{j}=[Tmtx{i}(w)-s0*bw;Tmtx{i}(w)+sx*bw]';
    end;
    
    if show
        %image(Tmtx{i},1:size(Z0,2),Z0'); colormap(clr);
        
        Zx=Za{i}; k=find(Zx==0); Zx(k)=nan;
        [nt,nu]=size(Zx); Zx=Zx.*repmat(1:nu,nt,1);
        
        hold off cla; plot(Tmtx{i},Zx,'.','Color',clr(28,:),'MarkerSize',15); hold on
        k=find(Z0~=0 & Z0~=28);
        for kk=1:length(k)
            plot(Tmtx{i}(mod(k(kk)-1,nt)+1),Zx(k(kk)),'.','Color',clr(Z0(k(kk)),:),'MarkerSize',15);
            plot(Tmtx{i}(mod(k(kk)-1,nt)+1),Zx(k(kk)),'o','Color',clr(Z0(k(kk)),:),'MarkerSize',15,'LineWidth',2);
        end;
        
        axis([min(Tmtx{i}) max(Tmtx{i}) 0 nu+1]), box off
        set(gca,'FontSize',22); ylabel('unit #'); xlabel('time (s)');
        hold on
        for t=1:length(EvtT{i})
            cc=ccL{EvtL{i}(t)+1};
            plot([EvtT{i}(t) EvtT{i}(t)],[0 p+1],cc,'LineWidth',2);
        end;
        %q=input([num2str(i) '...']);
    end;
    
    Za{i}=Z0;
end;


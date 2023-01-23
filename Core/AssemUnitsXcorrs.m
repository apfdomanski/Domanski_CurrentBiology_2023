%% Examine spike-spike cross-correlations in and out of assemblies
rat = 'JaroslawLONG2'
%% Set up for cross-correlograms
% concatenate all unit spike times
FR.Jointcells = [FR.PFCcells;FR.HPcells];
FR.iCells ={1:length(FR.PFCcells),(1:length(FR.HPcells))+length(FR.PFCcells)};
FR.unitLocs = [ones(length(FR.PFCcells),1) ; 2*ones(length(FR.HPcells),1)];
STMtxPre  = restrictTranges(FR.Jointcells,tRanges.Pre)';
STMtxTask = restrictTranges(FR.Jointcells,tRanges.Task)';
STMtxPost = restrictTranges(FR.Jointcells,tRanges.Post)';
% [no. Units, no. Timestamps]

% Can restrict to equal spike no's with this:
% mxsize=min([size(STMtxIN,1) size(STMtxOUT,1)]);
% if size(STMtxIN,1)>=mxsize;    STMtxIN  = STMtxIN(1:mxsize,:);  end
%% Run cross-correlation code
clear UUcorr

UUcorr.p.Tcorr = 0.05;
UUcorr.p.bw    = 0.001;
UUcorr.p.optf  = 2;  
%  0: raw CCTH, 
%  1: phase histo, 
%  2: cc/(n1*n2), 
%  3: true corr 

TinvPre  = [min(STMtxPre(STMtxPre>0)),max(max(STMtxPre))] ;
TinvTask = [min(STMtxTask(STMtxTask>0)),max(max(STMtxTask))] ;
TinvPost = [min(STMtxPost(STMtxPost>0)),max(max(STMtxPost))] ;

% Optional: check for low FR neurons to exclude
% Ncrit=100; unit_mask = ones(size(FR.Jointcells));
% for j=1:size(STMtxPre,1)
%     k=find(STMtxPre(j,:)>=Tinv(1) & STMtxPre(j,:)<=Tinv(2));
%     if length(k)<Ncrit
%         unit_mask(j,:)=0;
%         disph(['excluded unit# ' num2str(j)])
%     end;
% end;
UUcorr.unitmask=ones(length(FR.unitLocs),1);

% [UUcorr.Pre.Tb,UUcorr.Pre.CC]   = CompCCMtxBasic (STMtxPre, TinvPre, UUcorr.p.Tcorr,UUcorr.p.bw,UUcorr.p.optf); 
% [UUcorr.Task.Tb,UUcorr.Task.CC] = CompCCMtxBasic (STMtxTask,TinvTask,UUcorr.p.Tcorr,UUcorr.p.bw,UUcorr.p.optf); 
% [UUcorr.Post.Tb,UUcorr.Post.CC] = CompCCMtxBasic (STMtxPost,TinvPost,UUcorr.p.Tcorr,UUcorr.p.bw,UUcorr.p.optf); 
% [UUcorr.Pre.ccmx,UUcorr.Pre.ccidx]=max(UUcorr.Pre.CC');
% [UUcorr.Task.ccmx,UUcorr.Task.ccidx]=max(UUcorr.Task.CC');
% [UUcorr.Post.ccmx,UUcorr.Post.ccidx]=max(UUcorr.Post.CC');


[UUcorr.Task.Tb,...
 UUcorr.Task.CC,...
 UUcorr.Task.CCMax,...
 UUcorr.Task.CCDel,...
 UUcorr.Task.CCSig,...
 UUcorr.Task.CCSum,...
 UUcorr.Task.ccmxSurr,...
 UUcorr.Task.CCSurr]=CompCCMtx(STMtxTask, TinvTask, UUcorr.p.Tcorr, UUcorr.p.bw, UUcorr.p.optf);

% [UUcorr.Pre.Tb,...
%  UUcorr.Pre.CCbase,...
%  UUcorr.Pre.CCMaxbase,...
%  UUcorr.Pre.CCDelbase,...
%  UUcorr.Pre.CCSigbase,...
%  UUcorr.Pre.CCSumbase,...
%  UUcorr.Pre.ccmxSurrbase,...
%  UUcorr.Pre.CCSurrbase]=CompCCMtx(STMtxOUT,Tinv,Tcorr,bw,optf);
%% Plot xcorr histogram matrix
% temp_x=1000*LRcorr.Tb-bw/2;
% temp_x(:,[1:10 end-10:end])=[];

IDs=reshape(1:size(STMtxTask,1)^2,size(STMtxTask,1),size(STMtxTask,1))'; % IDs(diag(IDs)) = 0 ;
IDdiag=IDs(diag(IDs));
IDs=triu(IDs,0); %2nd arg 0 to include autocorrs
IDs=sort(IDs(IDs>0));

figure('name','Task','Color','w'); hold on
for iPlot=1:length(IDs)%size(STMtxTask,1)^2
    if ismember(IDs(iPlot),IDdiag)
%         disp([i IDs(i)])
        subaxis(size(STMtxTask,1),size(STMtxTask,1),IDs(iPlot),'spacing',0.001); hold on %IDs(i)
        

            
            temp=UUcorr.Task.CC(iPlot,1:end-1);
            temp=smooth_hist1(temp);
            temp=zscore(temp);
            %tempout(:,iPlot)=[iPlot;temp'];
            % tempout(:,sum(tempout(2:end,:),1)==0)=[];
            
    if isfield(UUcorr.Task,'CCSig')
            temp_Surr=UUcorr.ccmxSurr(:,iPlot);
            temp_Surr=smooth_hist1(temp_Surr);
            temp_Surr=zscore(temp_Surr);
            
        if UUcorr.Pre.CCSig(iPlot)
            ph(1) = patchline(LRcorr.Pre.Tb(1:end-1)-UUcorr.p.bw/2,temp_Surr,'edgecolor','r','linewidth',1,'edgealpha',1);       
            ph(1) = patchline(LRcorr.Pre.Tb-UUcorr.p.bw/2,temp ,'edgecolor','k','linewidth',1.5,'edgealpha',1);    
            plot([0 0],[-2 5],':k')

        else
%             ph(1) = patchline(UUcorr.Pre.Tb-bw/2,temp ,'edgecolor',[0.3 0.3 0.3],'linewidth',1,'edgealpha',.8);    
        end
    else
       ph(1) = patchline(UUcorr.Task.Tb-UUcorr.p.bw/2,temp ,'edgecolor','k','linewidth',1.2,'edgealpha',.8);
    end
	axis off; axis([-.01 .01 -2 3])
    end
end
clear iPlot IDs IDdiag temp 
%% work out which xcorrs are are which
% (1) Indicies of interactions in the x-corr matrix
% IDs=reshape(1:size(STMtxIN,1)^2,size(STMtxIN,1),size(STMtxIN,1))'; % IDs(diag(IDs)) = 0 ;
% IDs=triu(IDs,0);% IDs=triu(IDs,1); % to remove autocorrs

% (2) Indices related to xcorr matrix output
unitLocMasked=FR.unitLocs; unitLocMasked(~UUcorr.unitmask)=[];
IDs = triuIndex(length(unitLocMasked),0); 

% PFC<>PFC correlations
areaMask= zeros(size(IDs)); areaMask((unitLocMasked==1),(unitLocMasked==1))=1;
UUcorr.mask.PFCcorrs = IDs; UUcorr.mask.PFCcorrs(~areaMask)=0;
% UUcorr.mask.PFCcorrs = sort(UUcorr.mask.PFCcorrs(1:end),'descend'); UUcorr.mask.PFCcorrs=UUcorr.mask.PFCcorrs(UUcorr.mask.PFCcorrs>0);

% HP<>HP correlations
areaMask= zeros(size(IDs)); areaMask((unitLocMasked==2),(unitLocMasked==2))=1;
UUcorr.mask.HP_HPcorrs=IDs; UUcorr.mask.HP_HPcorrs(~areaMask)=0;
% UUcorr.mask.HP_HPcorrs = sort(UUcorr.mask.HP_HPcorrs(1:end),'descend'); UUcorr.mask.HP_HPcorrs = UUcorr.mask.HP_HPcorrs(UUcorr.mask.HP_HPcorrs>0);


% PFC<>HP correlations
areaMask= zeros(size(IDs)); areaMask((unitLocMasked==1),(unitLocMasked==2))=1;
UUcorr.mask.PFC_HPcorrs=IDs; UUcorr.mask.PFC_HPcorrs(~areaMask)=0;
% UUcorr.mask.PFC_HPcorrs = sort(UUcorr.mask.PFC_HPcorrs(1:end),'descend'); UUcorr.mask.PFC_HPcorrs = UUcorr.mask.PFC_HPcorrs(UUcorr.mask.PFC_HPcorrs>0);
%% Plot peak vs delay
TF_base     = triu(LRcorr.CCSigbase,1);
TF_stim     = TF_base;%triu(LRcorr.CCSigstim,1);
delay_stim  = triu(LRcorr.CCDelstim,1);              delay_base  = triu(LRcorr.CCDelbase,1);
corr_stim   = triu(LRcorr.CCMaxstim,1);              corr_base   = triu(LRcorr.CCMaxbase,1);
delay_stim  = delay_stim(TF_stim==1 & LRcons>1);     delay_base  = delay_base(TF_base==1 & LRcons>1);

corr_stim   = corr_stim(TF_stim==1 & LRcons>1)./(LRcorr.CCSumstim(TF_base==1 & LRcons>1));
corr_base   = corr_base(TF_base==1 & LRcons>1)./(LRcorr.CCSumbase(TF_base==1 & LRcons>1));
    
figure; hold on 
subplot(1,2,1)
    scatter(delay_base(1:end),corr_base(1:end),'k','filled')
    set(gca,'YScale','Log')
    axis([-Inf Inf 0 1e5])
    title('Baseline')
    xlabel('Peak Time-lag (s)')
    ylabel('Peak Correlation')
subplot(1,2,2)
    scatter(delay_stim(1:end),corr_stim(1:end),'g','filled')
    set(gca,'YScale','Log')
    axis([-Inf Inf 0 1e5])
    delay_delta=delay_stim-delay_base;
    corr_delta=corr_stim-corr_base;
    title('Illumination')   
    xlabel('Peak Time-lag (s)')
    ylabel('Peak Correlation')


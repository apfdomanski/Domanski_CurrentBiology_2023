% quality control for factor analysis detected assemblies
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk

clear all
close all
target= 'PreSleep';
if ispc
    home = 'C:\'; % [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    pat = [home 'Analysis\AssemblyAnalysis\Sleep'];
    separator='\';
elseif ismac
    path(path,'/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/MichalDataAna')
    home = getenv('HOME');
    cd ([home '/Documents/Bristol_work/assembly_analysis/MichalDataAna'])
    pat = [home '/Documents/Bristol_work/assembly_analysis/Michal_data/raw/50ms_bins/'];
	separator='/';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/50ms_bins/'];
    separator='/';
end

mkdir([pat separator target])
% run on 1) iFRsc, 2) KDE iFR

%bw=0.01; iFRtyp='iFRsc';
%bw=0.05; iFRtyp='iFRsc';

bw=0.05; iFRtyp='iFR';
A=dir([pat separator '*PFC_iFR' num2str(bw*1e3) strcat('*', target) '*.mat']);


for event_id=1:length(A)
    k=findstr(A(event_id).name,'_'); r=findstr(A(event_id).name,'.');
    Name{event_id}=[A(event_id).name(1:k(1)) A(event_id).name(k(2)+1:r-1)];
end;
minFR=0.1;     % minimal acceptable firing rate
critCvBW=1e6;  % critical max bound on kernel width drift over time (variance)
alpha=0.01;    % significance threshold
kmax=30;       % number of factors to check up to
Nbs=500;       % number of bootstrap draws
twin=10;       % Time window (s)
twin = twin/bw;

%% visualize results
for f=1%:length(A)
   
fnOut=[pat separator target separator Name{f} '_AssemRes2']; load(fnOut,'FL','LL','LLbs','nassem','FLbs')
fnOut=[pat separator target separator Name{f} '__FSCtemp']; load(fnOut,'FSC','ciHld','ciHsc','FSCbs');
    titles={'PFC','HP','Joint'};
    for s = 1:length(FSC)      % across all areas/assem types
        figure('NumberTitle','off','Name',[Name{f} ' : ' titles{s} ' Assemblies'],'Position',[20 570 1100 420])
        
        %%% (1) plot estimated number of Assemblies
        z   = diff(LL{s});     % gradient between successive Log-likelihood measures of factor numbers
        Zbs = diff(LLbs{s}')'; % gradient between successive Log-likelihood measures of Bootrapped factor numbers
        x_noAssems   = 2:length(z)+1;
        r=round(length(Zbs)*alpha);
        zz=sort(Zbs,'descend');
        subplot(1,3,1); hold on
            plot(x_noAssems,z,'b','LineWidth',2);
            %plot(Zbs','r');
            plot(x_noAssems,zz(r)*ones(1,length(z)),'r--','LineWidth',2);
            title(['Assembly count prediction'])
            xlabel({'No. Factors ' ;...
                ['(Estimated counts: AIC: ' num2str(nassem{s}(1)) ' / BIC: ' num2str(nassem{s}(2)) ' / BS: ' num2str(nassem{s}(3)) ' / pChi2: ' num2str(nassem{s}(4)) ')' ]})
            ylabel('Log-Likelihood Ratio')
            legend('Original','Trial-permutation'),legend('boxoff')
            set(gca,'FontSize',8)
        
        %%% (2) plot factor loading distributions
        x_FLdist=-0.5:0.1:1;
        vs = cell2mat(FLbs{s}{2});                         % Bootstrapped factor loadings
        h  = histc(vs(1:end),x_FLdist)./length(vs(1:end)); % Distribution of BS factor loadings
        temp = FL{s}{2}(1:end);                            % Original factor loadings
        subplot(1,3,2); hold on
            plot(x_FLdist,h,'b','LineWidth',2)                          % Distribution of BS factor loadings
            plot([ciHld(s,2) ciHld(s,2)],[0 0.05],'g','LineWidth',2);   % Confidence limit for significant factor loading
            plot(temp(temp>ciHld(s,2)),0.01,'go');
            plot(temp,0.01,'r.');
            plot(temp(temp>ciHld(s,2)),0.01,'g.');
            xlim([min(x_FLdist) 1.2*max([max(x_FLdist), ciHld(s,2)])]);
            legend('Factor loading distribution',...
                '95th percentile of B/Strapped range',...
                'FL>95% of B/Strapped range') , legend('boxoff')
            title('Strength of assembly membership')
            xlabel('Factor loading strength')
            ylabel('Distribution')
            set(gca,'FontSize',8)
        
        %%% (3) plot factor score distributions
        ws=sort(FSCbs{s}(1:end),'ascend');
        ciV=[0.1 0.05 0.01 5e-3 1e-3];
        x_FSCdist=-3:0.1:10;%ceil(max(FSC{s})); 
        h=histc(ws,x_FSCdist)./length(ws);
        temp = FSC{s}(find(rand(size(FSC{s}))>.999)); %only plot 0.1% to keep the plotting overhead down
        for i=1:length(ciV)
            ciLsc(s,i)=ws(round(length(ws)*ciV(i)));
            ciHsc(s,i)=ws(round(length(ws)*(1-ciV(i))));
        end;
        subplot(1,3,3); hold on
            plot(x_FSCdist,h,'b',[ciHsc(s,2) ciHsc(s,2)],[0 0.05],'g','LineWidth',2);
            plot(FSC{s}(FSC{s}>ciHsc(s,2)),0.01,'go')
            plot(temp,0.01,'r.'); 
            plot(FSC{s}(FSC{s}>ciHsc(s,2)),0.01,'g.')
            xlim([min(x_FSCdist) max(x_FSCdist)]);
            legend('Factor activation distribution',...
                '95th percentile of B/Strapped range',...
                'Scores>95% of B/Strapped range') , legend('boxoff')
            title('Strength of assembly activation')
            xlabel('Factor score')
            ylabel('Distribution')
            set(gca,'FontSize',8)
    end
end
%% visualize factor scores/ activations in time & distribution
for f=1:length(A)
    fnOut=[pat '/' target '/' Name{f} '_AssemRes2']; load(fnOut);
    fnOut=[pat '/' target '/' Name{f} '__FSCtemp']; load(fnOut);
    titles={'PFC','HP','Joint'};
    %%% (1) plot factor scores versus time to show confidence thresholds on top
    for s=1:length(FSC)
        [nTimePoints,nAssems]=size(FSC{s});
        for Assem_id=1%:nAssems
            figure('NumberTitle','off','Name',[Name{f} ' : ' titles{s} ' assembly no. ' num2str(Assem_id) ' of ' num2str(nAssems)],'Position',[20 570 1100 420])
            plot_factor = 1; % downsampling factor
            n = find(abs(FL{s}{nassem{s}(3)}(:,Assem_id))>ciHld(s,3)); % Units involved
            t = (1:ceil(nTimePoints/plot_factor))*plot_factor*0.05+0.05/2;
%             t = decimate(cell2mat(TmtxS{s}),plot_factor);
            subplot('Position' , [0.1 0.1 0.7 0.8]); hold on
                plot(t,ciHsc(s,2),'g','LineWidth',2);
                plot(t,decimate(FSC{s}(:,Assem_id),plot_factor),'b','LineWidth',2);
                set(gca,'FontSize',12)
                xlabel('Time (s)')
                ylabel('Assembly activation strength') 
                xlim([min(t) max(t)]);
%                 title(['No. units involved: ' num2str(n')]);
                title([titles{s} ' assembly no. ' num2str(Assem_id) ' of ' num2str(nAssems)])
                xlims=get(gca,'ylim');
                
            x=-3:0.01:ceil(max(FSC{s}(:,Assem_id))); 
            h=histc(FSC{s}(:,Assem_id),x)./length(FSC{s}(:,Assem_id));
            subplot('Position' , [0.82 0.1 0.15 0.8]); hold on
                plot([0 0.05],[ciHsc(s,2) ciHsc(s,2)],'g','LineWidth',2);
                plot(h,x,'b','LineWidth',2);
                xlim([0 1.2*max(h)])
                ylim(xlims);
                set(gca,'YTick',[])
                set(gca,'FontSize',12)
                xlabel('Assembly activation strength')
                title('Distribution') 
                legend('95% B/strap CI') , legend('boxoff')
%             title(['#units: ' num2str(n')]);
        end;
    end;
    %%% (1) plot factor scores versus time for all with 
    figure('NumberTitle','off','Name',[Name{f} ' : All assemblies'],'Position',[20 570 1100 420]); 
    hold on
    colours = hsv(length(FSC)); % one per colour per  area
    for s=1:length(FSC)
        [nTimePoints,nAssems]=size(FSC{s});
        plot_factor = 1; % downsampling factor
        for Assem_id=1%:nAssems
            t = (1:ceil(nTimePoints/plot_factor))*plot_factor*0.05+0.05/2;
            plot(t,decimate(FSC{s}(:,Assem_id),plot_factor),'color',colours(s,:),'LineWidth',1.2);
            legend(titles) , legend('boxoff')
             xlabel('Time (s)')
             ylabel('Assembly activation strength') 
        end
    end      
    
    %%% (1) plot factor scores versus time for all overlaid 
    figure('NumberTitle','off','Name',[Name{f} ' : All assemblies'],'Position',[20 570 1100 420]); 
    hold on
    colours = hsv(length(FSC)); % one per colour per area
    for s=1:length(FSC)
        [nTimePoints,nAssems]=size(FSC{s});
        plot_factor = 1; % downsampling factor
        for Assem_id=1:nAssems
            t = (1:ceil(nTimePoints/plot_factor))*plot_factor*0.05+0.05/2;
%                         t = decimate(cell2mat(TmtxS{s}),plot_factor);

            plot(t,decimate(FSC{s}(:,Assem_id),plot_factor),'color',colours(s,:),'LineWidth',2);
%             legend(titles) , legend('boxoff')
             xlabel('Time (s)')
             ylabel('Assembly activation strength') 
        end
    end      
end;


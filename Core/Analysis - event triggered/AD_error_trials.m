% Loads and plots the sequence of correct and error trials from behavioural data and calculates the (Markov-esque) 1st orde transition matrix (sequence probabily matrix)
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk


%% Preamble: load files
clear
% start in 'raw'
% cd 'C:\Analysis\AssemblyAnalysis\raw'
% (1) load raw
% load('JaroslawLONG2.mat')

% (2) load iFR50
load('50ms_bins/JaroslawLONG2_PFC_iFR50.mat')
%include errors?:
[EvtTs_,s]=sort([EvtT;EvtTinc]);
EvtLs_ = [EvtL';EvtLinc'/10];
%otherwise...
% [EvtTs_,s]=sort(EvtT);
% EvtLs_ = [EvtL];
% 
EvtLs_ = EvtLs_(s);
Error_times=EvtTinc;


% (3) get aligned trial markers
% [TmtxS,iFR0,EvtLs,EvtTs,usel]=SelTrialsCells('KDE_bins/JaroslawLONG2_PFC_iFR50.mat',...
%                                              10,...         % time window (s)
%                                              0.1,...        % minimum acceptable firing rate
%                                              1e6,...        % maxmium permitted variance in BW over time 
%                                              'iFR',...      % data type: Firing rate or spike counts
%                                              'all');        % include errors? ("corr" N or "all" Y)

% Remove error trial markup (were multiplied by 10)
% EvtTs_=EvtTs;
% EvtLs_=EvtLs;
% EvtLs_(rem(EvtLs,10)==0)=EvtLs(rem(EvtLs,10)==0)/10;
% Error_times=EvtTs(rem(EvtLs,10)==0);


% load('KDE_bins/JaroslawLONG2_iFR50_FSC.mat')
%% Plot sample/choice sequence along with error trials, add transition matrix

% test:
% EvtLs__=repmat([1:4],1,5)
% options={'one','two','three','four'};

options={'Left Choice','Right Choice','Left Sample','Right Sample'};
trans_matrix=hmmestimate(EvtLs_,EvtLs_);

figure
subplot(1,2,1); hold on
plot([Error_times,Error_times],[0,5],'r')
plot(EvtTs_,EvtLs_,'-ob')
axis tight
legend({'Errors'})
xlabel('Time (s)')
ylabel('Event sequence')
set(gca,'ytick',[1:4],'yticklabel',options)

subplot(1,2,2)
imagesc(flipud(rot90(trans_matrix)))%; set(gca,'YDir','normal')
set(gca,'xtick',[1:4],'xticklabel',options)
set(gca,'ytick',[1:4],'yticklabel',options)
colorbar
axis square
title('Choice probability transition matrix')
xlabel('Event (t)')
ylabel('Event (t+1)')
% e.g. p(Left Sample -> Right choice)
trans_matrix(3,2)

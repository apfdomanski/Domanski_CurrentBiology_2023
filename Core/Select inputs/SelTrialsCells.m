function [TmtxS,iFRs,EvtLs,EvtTs,usel_out]=SelTrialsCells(fn,twin,minFR,critCvBW,iFRtype,trials)
% Select cells from PFC and HC data sets & cut out event periods (trial case)
%
% (c) 2014 Daniel Durstewitz, Bernstein Center for Computational Neuroscience, CIMH/ Heidelberg University
% Edited by Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Corrects for mismatched unit numbers if using data generated on Blue Crystal.
% Cuts up time series into equal length epochs for use in downstream
% bootstrapping.
%
% outputs
% Tmtxs {no. trials}                                         time axes for selected trials
% iFRs  {no. brain areas}{no. trials}(timerange,no. neurons) spike FR matrix
% EvtLs (no. events)                                         event labels       
% EvtTs (no. events)                                         event times
if nargin<2 || isempty(twin),     twin=10;       end; % time window (s)
if nargin<3 || isempty(minFR),    minFR=0.5;     end; % minimum acceptable firing rate
if nargin<4 || isempty(critCvBW), critCvBW=1;    end; % maxmium permitted variance in BW over time
if nargin<5 || isempty(iFRtype),  iFRtype='iFR'; end; % data type: Firing rate or spike counts
if nargin<6 || isempty(trials),   trials='corr'; end; % include errors? ("corr" N or "all" Y)


fn0=cell(1,2);
fn0{1}=fn; k=strfind(fn,'PFC');
fn0{2}=[fn(1:k-1) 'HP' fn(k+3:end)];
    
T=cell(1,2); iFRs=T; TmtxS=T;
for f=1:length(fn0) % across brain regions: PFC, HP
    load(fn0{f}); iFR0=eval(iFRtype);
    if strncmp(trials,'all',3)
        k=find(EvtTinc>min(EvtT) & EvtTinc<max(EvtT));
        EvtT=[EvtT;EvtTinc(k)]; 
        EvtL=[EvtL,EvtLinc(k)];
    end;
    try
        usel=find(avgFR>=minFR & sum(iFR0)>0 & CvBW<=critCvBW); % indices of cells to keep
    catch
        usel=find(avgFR>=minFR & sum(iFR0)>0 ); % indices of cells to keep... equal bin case
    end    
    usel_out{f}=usel;
    Tmtx=Tmtx(1:end-1); T{f}=Tmtx; iFR0=iFR0(:,usel);
    [EvtTs,s]=sort(EvtT); % event times
    EvtLs=EvtL(s);        % event types
    % {'left_choice','right_choice','left_sample','right_sample'};
    iFRs{f}=cell(1,round(length(EvtTs)/2)); TmtxS{f}=iFRs{f};
    for i=1:2:length(EvtTs)-1
        k=find(Tmtx>=EvtTs(i)-twin & Tmtx<=EvtTs(i+1)+twin);
        iFRs{f}{ceil(i/2)}=iFR0(k,:);
        TmtxS{f}{ceil(i/2)}=Tmtx(k);
    end;
    
    %%%% AD 9/2/16 edit - add this condition to handle jumbled unit IDs
    if exist('files')
        % This means that units' firing rates were added in piecemeal fashion as they
        % were calculated on Blue Crystal, thus not necessarily in the
        % right order.
        Name ={}; k=[]; r=[];
        for file_ID=1:length(files)
            k=findstr(files(file_ID).name,'_'); r=findstr(files(file_ID).name,'.');
            Name{file_ID}=[files(file_ID).name(k(end)+1:r-1)];
            Name{file_ID}=str2num(strtok(Name{file_ID}, 'of'));
        end;
        Name=cell2mat(Name);
        temp=[];
        for u=1:length(usel_out{f})
            temp(u)=Name(usel_out{f}(u));
        end
        usel_out{f}=sort(temp);
    else
        % otherwise keep as is since the units are already in the right order
    end
    
end;
if ~isempty(find(T{1}~=T{2},1)), warning('T1 and T2 do not match!'); end;

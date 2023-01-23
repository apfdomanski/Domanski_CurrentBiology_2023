function [TmtxS,iFRs,unit_IDs] = SelCellsSleep(fn,twin,minFR,critCvBW,iFRtype)
% select cells from PFC and HC data sets (Continuous data case)
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Corrects for mismatched unit numbers if using data generated on Blue Crystal.
% Cuts up time series into equal length epochs for use in downstream
% bootstrapping
%
%
%
%

if nargin<2 || isempty(minFR),    minFR=0.05;    end; % minimum acceptable firing rate
if nargin<3 || isempty(critCvBW), critCvBW=1;    end; % maxmium permitted variance in BW over time
if nargin<4 || isempty(iFRtype),  iFRtype='iFR'; end; % data type: Firing rate or spike counts
%outputs
% Tmtxs {no. trials}                                         time axes for selected trials
% iFRs  {no. brain areas}{no. trials}(timerange,no. neurons) spike FR matrix
% EvtLs (no. events)                                         event labels       
% EvtTs (no. events)                                         event times

fn0=cell(1,2);
fn0{1}=fn; k=strfind(fn,'PFC');
fn0{2}=[fn(1:k-1) 'HP' fn(k+3:end)];
iFRs=cell(1,2); TmtxS=cell(1,2);

for f_id=1:length(fn0) % across brain regions: PFC, HP
    disp(['Loading... ' fn0{f_id}])
    load(fn0{f_id}); tempRate=eval(iFRtype);
    % get useful cells 
    try
        unit_IDs{f_id}=find(avgFR>=minFR & sum(tempRate)>0 & CvBW<=critCvBW); % indices of cells to keep
    catch
        unit_IDs{f_id}=find(avgFR>=minFR & sum(tempRate)>0 ); % indices of cells to keep... equal bin case
    end    
   
    % get time and units and subdivide into epochs (for bootstrap analysis)
    %%% need multiple segments to mimic trials, or could collapse all into one mega-trial
    % tTmtxS{f}{1}  = Tmtx(1:end-1);  
    % iFRs{f}{1}= tempRate(:,tempUnits);
    noEpochs = floor((length(Tmtx)-1)/twin);
    for epoch_id=1:noEpochs
        TmtxS{f_id}{epoch_id} = Tmtx(     (twin*(epoch_id-1)+1) : twin*epoch_id);  
        iFRs{f_id}{epoch_id}  = tempRate(((twin*(epoch_id-1)+1) : twin*epoch_id),unit_IDs{f_id});
    end
    
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
        for u = 1:length(unit_IDs{f_id})
                temp(u) = Name(unit_IDs{f_id}(u));
        end
        unit_IDs{f_id}=sort(temp);
    else
        % otherwise keep as is since the units are already in the right order
    end
    
end
disp('Loading... Done.')


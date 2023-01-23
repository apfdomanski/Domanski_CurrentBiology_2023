function usel_out=SelCells(fn,minFR,critCvBW)
% Select cells from PFC and HC data sets & cut out event periods (trial case)

% outputs
% Tmtxs {no. trials}                                         time axes for selected trials
% iFRs  {no. brain areas}{no. trials}(timerange,no. neurons) spike FR matrix
% EvtLs (no. events)                                         event labels       
% EvtTs (no. events)                                         event times
if nargin<3 || isempty(minFR),    minFR=0.5;     end; % minimum acceptable firing rate
if nargin<4 || isempty(critCvBW), critCvBW=1;    end; % maxmium permitted variance in BW over time


fn0=cell(1,2);
fn0{1}=fn; k=strfind(fn,'PFC');
fn0{2}=[fn(1:k-1) 'HP' fn(k+3:end)];
    
for f=1:length(fn0) % across brain regions: PFC, HP
    load(fn0{f}); iFR0=eval('iFR');

    try
        usel=find(avgFR>=minFR & sum(iFR0)>0 & CvBW<=critCvBW); % indices of cells to keep
    catch
        usel=find(avgFR>=minFR & sum(iFR0)>0 ); % indices of cells to keep... equal bin case
    end    
    usel_out{f}=usel;
  
end

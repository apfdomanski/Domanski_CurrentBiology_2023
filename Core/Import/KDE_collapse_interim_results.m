% Collect interim single cell data from KDE convert code into HP and PFC archives
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
% 
% You need to run this after big batch jobs on Blue Crystal to gather all
% individually saved single unit files into one population .MAT file.
%% Collate HP cells
clear
files=dir(['*HP*' '*.mat']);

for file_id = 1:size(files,1)
    disp(['Reading ' files(file_id).name])
    temp=load(files(file_id).name);
    avgFR(1,file_id) = temp.avgFR_u;
    hucv(1,file_id)	 = temp.hucv_u;
    MISE(1,file_id)  = temp.MISE_u;
    CvBW(1,file_id)  = temp.CvBW_u;
    hucvS{file_id}   = temp.hucvS_u;
    iFR(:,file_id) = temp.iFR_u;    
    delete(files(file_id).name)
end
EvtT= temp.EvtT;
EvtL= temp.EvtL;
Tmtx= temp.Tmtx;

fname=strsplit(files(1).name,'_HPu'); fname=[fname{1} '.mat'];
save(fname,...
    'files',...% File list just in case of mix ups
    'iFR',...  % KDE optimised firing rates (size; time vector,no.units)
    'Tmtx',... % time vector
    'EvtT',... % event times (samples and choices)
    'EvtL',... % event labels (left/right, sample/choice)
    'avgFR',...% each neuron's average firing rate
    'CvBW',... % variance in optimal kernel width over repeated trials
    'hucv',... % each neuron's optimised kernel width
    'MISE',... % mean-squared-error for each neuron's optimal kernel width
    'hucvS');  % optimised kernel for each trial individually
%% Collate PFC cells
clear
files=dir(['*PFC*' '*.mat']);

for file_id = 1:size(files,1)
    disp(['Reading ' files(file_id).name])
    temp=load(files(file_id).name);
    avgFR(1,file_id) = temp.avgFR_u;
    hucv(1,file_id)	 = temp.hucv_u;
    MISE(1,file_id)  = temp.MISE_u;
    CvBW(1,file_id)  = temp.CvBW_u;
    hucvS{file_id}   = temp.hucvS_u;
    iFR(:,file_id)   = temp.iFR_u;      
    delete(files(file_id).name)
  
end
EvtT= temp.EvtT;
EvtL= temp.EvtL;
Tmtx= temp.Tmtx;

fname=strsplit(files(1).name,'_PFCu'); fname=[fname{1} '.mat'];
save(fname,...
    'files',...% File list just in case of mix ups
    'iFR',...  % KDE optimised firing rates (size; time vector,no.units)
    'Tmtx',... % time vector
    'EvtT',... % event times (samples and choices)
    'EvtL',... % event labels (left/right, sample/choice)
    'avgFR',...% each neuron's average firing rate
    'CvBW',... % variance in optimal kernel width over repeated trials
    'hucv',... % each neuron's optimised kernel width
    'MISE',... % mean-squared-error for each neuron's optimal kernel width
    'hucvS');  % optimised kernel for each trial individually

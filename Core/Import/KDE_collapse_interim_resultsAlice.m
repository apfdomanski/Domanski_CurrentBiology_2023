% Collect interim single cell data from KDE convert code into HP and PFC archives
% No event times pass through this version
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
% 
% You need to run this after big batch jobs on Blue Crystal to gather all
% individually saved single unit files into one population .MAT file.
%% Collate data
clear
targets = {'hipp', 'amygd'};
for iTarget=1:length(targets);
    avgFR = []; hucv  = []; MISE  = [];
    CvBW  = []; hucvS = []; iFR   = [];   
    
    files=dir(['*' targets{iTarget} '*.mat']);
    fname_ = strsplit(files(1).name,'u_');
    
    for file_id = 1:size(files,1)
       
        fname = [fname_{1} 'u_' num2str(file_id) 'of' num2str(size(files,1)) '.mat'];
        disp(['Reading ' fname])
        temp=load(fname);
        avgFR(1,file_id) = temp.avgFR_u;
        hucv(1,file_id)	 = temp.hucv_u;
        MISE(1,file_id)  = temp.MISE_u;
        CvBW(1,file_id)  = temp.CvBW_u;
        hucvS{file_id}   = temp.hucvS_u;
        iFR(:,file_id) = temp.iFR_u;    
        delete(fname)
    end
    
    Tmtx= temp.Tmtx;

     fname=[fname_{1} '_all.mat'];
    save(fname,...
        'files',...% File list just in case of mix ups
        'iFR',...  % KDE optimised firing rates (size; time vector,no.units)
        'Tmtx',... % time vector
        'avgFR',...% each neuron's average firing rate
        'CvBW',... % variance in optimal kernel width over repeated trials
        'hucv',... % each neuron's optimised kernel width
        'MISE',... % mean-squared-error for each neuron's optimal kernel width
        'hucvS');  % optimised kernel for each trial individually
end


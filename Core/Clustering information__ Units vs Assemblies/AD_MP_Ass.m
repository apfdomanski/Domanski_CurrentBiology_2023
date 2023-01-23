function [rEnergy, f, t] = AD_MP_Ass(data,params)
% Temporary workspace for MP decomposition
folderName = ['/Users/aleksanderdomanski/Documents/MATLAB/MPtest' '/data/'];
tag = 'test/';

Max_iterations = 1000; % number of iterations
% Import the data


L_orig = size(data,1);           % Unpadded signal length
L  = 2^nextpow2(size(data,1));   % Padded signal length
Ntrials = size(data,2);          % No. trials

% Pad data up to nearest power of 2 length
inputSignal=zeros(L,Ntrials,1);
inputSignal(1:size(data,1),:)=data;

importData(inputSignal,folderName,tag,[1 L],params.Fs);

% perform Gabor decomposition
runGabor(folderName,tag,L , Max_iterations);
params.gaborInfo = getGaborData(folderName,tag,1);


f = 0:(params.Fs/params.downres_factor)/L:(params.Fs/params.downres_factor)/2;
t  = (0:L-1)/(params.Fs/params.downres_factor);

% construct spectrogram and optionally reconsruct input signal
wrap=1;
atomList=[]; % all atoms, or optionally choose which to include in reconstruction
rEnergy = zeros(length(f),L,Ntrials);
for iTrial=1:Ntrials
    
    if isempty(atomList)
        disp(['Reconstructing trial ' num2str(iTrial) ', all atoms']);
    else
        disp(['Reconstructing trial ' num2str(iTrial) ', atoms ' num2str(atomList(1)) ':' nusm2str(atomList(end))]);
    end
    % reconstruct signal
    %     rSignal(:,iTrial)   = reconstructSignalFromAtomsMPP(params.gaborInfo{iTrial}.gaborData,L,wrap,atomList);
    % reconstruct energy
    rEnergy(:,:,iTrial) = reconstructEnergyFromAtomsMPP(params.gaborInfo{iTrial}.gaborData,L,wrap,atomList);
end
clear iTrial wrap tag folderName atomList
rEnergy = rEnergy(:,1:L_orig,:);


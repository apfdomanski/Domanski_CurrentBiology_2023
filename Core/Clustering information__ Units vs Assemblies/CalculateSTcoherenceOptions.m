function [coh_,f] = CalculateSTcoherenceOptions(testspikes,STmultiplier,WL,tlims)
% Calulates a signal-noise-ration based measurement of single unit
% frequency locking, across a spectrum. Uses windowed Welch's method for
% calculation and has an optional adjustment for Poisson firing rate
% assumption based on refractory period effects.
% (c) Aleks PFC Domanski (UoB) 2019
% Inspired by Matzer & Ben-Gad  PLOS Computational Biology 2015 https://doi.org/10.1371/journal.pcbi.1004252
%
% Inputs:
% ST:           Spike times (struct array with times in seconds stored in ST{iUnit}.t)
% STmultiplier: Optional correction factor to bring spike times to seconds (typically 1e-4)
% WL:           Hamming window length (Optional)
% tlims:        Optional time lims to trim the spikes between (Optional)
%
% Outputs:
% ModIndex:           Spectral modulation index for each queried neuron [freq x no. Units]
% f:                  Returned frequency vector in Hz
% powWelch:           Uncorrected spike train power spectrum (Welch's method represented as SNR vs 100~500Hz mean power)
% powWelchCorrected:  Power spectrum corrected for window power and baseline firing rate
% r0:                 Baseline firing rate for each neurons
%
% NB choose window length to cover units' estimated freq following jitter, i.e. +/- [X]Hz around central freq assuming 1ms bin width (500Hz Nyquist):
% 100 points ~ +/-5Hz
% 500 points ~ +/-2.5Hz
% 1000 points ~ +/-0.5Hz
%
%%%% Example calling script: 
% Fs = 1e3;         
% t_max = 60*Fs;
% t  = 0:Fs^-1:t_max/Fs; % timebase (s)
% Rate0 = sin(2*pi*15 *t) + sin(2*pi*66 *t); % Example underlying rate comprised of uncoupled 15Hz, 66Hz oscillations 
% % Generate Poisson(ish) spikes to follow Rate0:
% a = (Rate0>0); y = zeros(size(a));
% y(a==0) = rand(nnz(a==0),1)<0; % floor of p(spike)
% y(a==1) = rand(nnz(a==1),1)<0.2; % ceiling of p(spike)
% testspikes{1}.t = find(y==1)/Fs';
% [ModIndex,f,powWelch] = CalculateFrequencyModulation(testspikes,STmutliplier);
% figure;
% subplot(2,1,1);hold on
% plot(f,mean(powWelch,2),'r')
% set(gca,'Xscale','log')
% subplot(2,1,2);hold on
% plot(f,mean(modIndex,2),'b')
% set(gca,'Xscale','log')
%%%%

if nargin<2 | isempty(STmultiplier)
    STmultiplier = 1;
end
if nargin<3 | isempty(WL)
    WL = 2;
end
if nargin<4 | isempty(tlims)
    tlims =[-Inf Inf];
end


bw = 2e-2; % Spike binning vector size
Fs = 1/bw; % Sampling frequency
WL = WL*Fs; % window length

win_ = hamming(WL);
nWins = 1;
NFFT = 2^nextpow2(WL);
f = Fs*(0:(NFFT/2))/NFFT; % frequency vector
f = 0:0.25:25;
nOverlap = 0;
winpow_ = win_ .* conj(win_);
winCorrectionFactor = max(( mean(win_) * WL) / winpow_ );

NormIdx = f>=100 & f<=500; % 'Noise' range of spectrum for SNR calculation
noUnits = length(testspikes);


NormIdx = f>=100 & f<=500; % 'Noise' range of spectrum for SNR calculation
noUnits = length(testspikes);
Spctm            = nan(length(f),noUnits);
% powWelchCorrected   = nan(length(f),noUnits);
% modIndex            = nan(length(f),noUnits);


% iUnit = 3; jUnit =5;
coh_ = cell(noUnits);
for iUnit = 40%1:noUnits
    for jUnit = 1:noUnits    
        if jUnit>iUnit
            
            fprintf('Spike train coherence: Unit [#%d,%d/%d]...\n',iUnit,jUnit,noUnits)
            ST_i = testspikes{iUnit}.t*STmultiplier;
            ST_j = testspikes{jUnit}.t*STmultiplier;
            
            ST_i(ST_i<tlims(1) | ST_i>tlims(2)) = [];
            ST_j(ST_j<tlims(1) | ST_j>tlims(2)) = [];
            
     if ~isempty(ST_i) & ~isempty(ST_j)
    tb = min([ST_i;ST_j]):bw:max([ST_i;ST_j]);
    STraster_i = histc(ST_i,tb);
    STraster_j = histc(ST_j,tb);
    
    r0_i = (sum(STraster_i)./length(STraster_i))*Fs; % Estimated mean firing rate, Poisson assumption from binned data
    r0_j = (sum(STraster_j)./length(STraster_j))*Fs; % Estimated mean firing rate, Poisson assumption from binned data
%     r0(iUnit) = r0_i;
    STraster_i = STraster_i(:)./r0_i;
    STraster_j = STraster_j(:)./r0_j;
%     Si = pwelch(STraster_i,win_,nOverlap,f,Fs,)/r0_i;
%     Sj = pwelch(STraster_j,win_,nOverlap,f,Fs)/r0_j;
%     tb=(1:10000)/1000;
%     STraster_i = sin(2*pi*8*tb')+sin(2*pi*12*tb');
%     plot(tb,STraster_i)
%     f = 1000*(0:(NFFT/2))/NFFT; % frequency vector
% 
%     Si = localComputeSpectra(STraster_i,win_,200,NFFT);
%     plot(f,10*log10(abs(Si)))
    
%     Si = localComputeSpectra(STraster_i,win_,200,NFFT)./r0_i;
%     Sj = localComputeSpectra(STraster_j,win_,200,NFFT)./r0_j;
% 
%     figure;hold on
%     plot(f,10*log10(abs(Si))); 
%     plot(f,10*log10(abs(Sj)));
%     
%     Sij = conj(Si) .* Sj;
%     Sii = conj(Si) .* Si;
%     Sjj = conj(Sj) .* Sj;
%     
%     figure;hold on
%     plot(f,10*log10(abs(Sii))); 
%     plot(f,10*log10(abs(Sjj)))
%     
%     Cij=(abs(Sij).^2)./( Sii .* Sjj );
%     Cij = abs(Sij)./(sqrt(Sii.*Sjj));
% 
% 
%     Coh = abs(Cij);
%     plot(f,Cij)
    if length(STraster_i)>WL & length(STraster_j)>WL 
%         try
            [Cxy,F] = mscohere(STraster_i,STraster_j,win_,0,f,Fs);
            % figure;plot(F,Cxy)
            coh_{iUnit,jUnit} = Cxy;
%         catch
%             warning(sprintf('Spike train frequency modulation :Insufficient spikes in window for Unit [#%d, %d/%d]!',iUnit,jUnit,noUnits))
%         end
    else
        warning(sprintf('Spike train frequency modulation :Insufficient spikes in window for Unit [#%d, %d/%d]!',iUnit,jUnit,noUnits))
    end
%     powWelchCorrected(:,iUnit) = powWelch(:,iUnit)*Fs/2 * nWins-(nWins-1) * winCorrectionFactor;


     end
    else
        warning(sprintf('Spike train frequency modulation :No spikes in window for Unit [#%d, %d/%d]!',iUnit,jUnit,noUnits))
    end
    end
end
end

function [S] = localComputeSpectra(sig,win_,nOverlap,NFFT)
    N = length(sig);
    L = length(win_);
    % Pad up to full length of multiples of window length
    sig = [sig;zeros(L-rem(N,L),1)];
    N = length(sig);
    nSeg = (N/L);
  
    
    

    sigBuff = buffer(sig,L,nOverlap);
    sigBuff = bsxfun(@times,sigBuff,win_);
    S = (sum(fft(sigBuff,NFFT),2)./nSeg);   
    S = S(1:NFFT/2+1);
%     plot(10*log10(abs(S)))
end
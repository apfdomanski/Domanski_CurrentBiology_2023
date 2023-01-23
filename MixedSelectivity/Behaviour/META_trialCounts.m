clear

pat = 'C:\Analysis\AssemblyAnalysis\raw\';

fileList = dir([pat 'allTimestamps' filesep '*SHORT*.mat']);
N_SHORT = [];N_SHORTerr = [];
for iFile =1:length(fileList)
    load([pat 'allTimestamps' filesep fileList(iFile).name],'t');
    Delays_ = fieldnames(t);
    N = zeros(length(Delays_),2);
    Ne = zeros(length(Delays_),2);
    for iDelay = 1:length(Delays_)
        eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
        N(iDelay,:) = [length(t_.CueLight_LeftCorrect),length(t_.CueLight_RightCorrect)];
        Ne(iDelay,:) = [length(t_.CueLight_LeftError),length(t_.CueLight_RightError)];
    end
    N_SHORT= [N_SHORT;N];
    N_SHORTerr=[N_SHORTerr;Ne];
end 



fileList = dir([pat 'allTimestamps' filesep '*MEDIUM*.mat']);
N_MEDIUM = [];N_MEDIUMerr = [];
for iFile =1:length(fileList)
    load([pat 'allTimestamps' filesep fileList(iFile).name],'t');
    Delays_ = fieldnames(t);
    N = zeros(length(Delays_),2);
    Ne = zeros(length(Delays_),2);
    for iDelay = 1:length(Delays_)
        eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
        N(iDelay,:) = [length(t_.CueLight_LeftCorrect),length(t_.CueLight_RightCorrect)];
        Ne(iDelay,:) = [length(t_.CueLight_LeftError),length(t_.CueLight_RightError)];
    end
    N_MEDIUM= [N_MEDIUM;N];
    N_MEDIUMerr=[N_MEDIUMerr;Ne];
end 

fileList = dir([pat 'allTimestamps' filesep '*LONG*.mat']);
N_LONG = [];N_LONGerr = [];
for iFile =1:length(fileList)
    load([pat 'allTimestamps' filesep fileList(iFile).name],'t');
    Delays_ = fieldnames(t);
    N = zeros(length(Delays_),2);
    Ne = zeros(length(Delays_),2);
    for iDelay = 1:length(Delays_)
        eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
        N(iDelay,:) = [length(t_.CueLight_LeftCorrect),length(t_.CueLight_RightCorrect)];
        Ne(iDelay,:) = [length(t_.CueLight_LeftError),length(t_.CueLight_RightError)];
    end
    N_LONG= [N_LONG;N];
    N_LONGerr=[N_LONGerr;Ne];
end 

bins = 1:50;
figure; hold on

plot(bins,cumsum(histc(N_SHORT(:),bins))./numel(N_SHORT),'b')
plot(bins,cumsum(histc(N_MEDIUM(:),bins))./numel(N_MEDIUM),'r')
plot(bins,cumsum(histc(N_LONG(:),bins))./numel(N_LONG),'k')
plot(bins,cumsum(histc(N_SHORTerr(:),bins))./numel(N_SHORTerr),':b')
plot(bins,cumsum(histc(N_MEDIUMerr(:),bins))./numel(N_MEDIUMerr),':r')
plot(bins,cumsum(histc(N_LONGerr(:),bins))./numel(N_LONGerr),':k')
legend({'Early (Correct)','Mid (Correct)','Late (Correct)','Early (Error)','Mid (Error)','Late (Error)'},'Location','southeast');legend boxoff
xlabel('No. Trials')
ylabel('Fraction of experiments')
clear 
Targets = {'SHORT','MEDIUM','LONG'};
Delays_ = {'Short','Medium','Long'};

bw=0.05;


warning ('off')
pat = 'C:\Analysis\AssemblyAnalysis\raw';
cd(pat)
fileList=dir('allTimestamps\*.mat');
Areas = {'HP','PFC'};



%%%%%
iFile =1%:length(fileList)
iArea =1% 1:length(Areas)
%%%%%
%% Load an example
fname=strtok(fileList(iFile).name,'_');

load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
load(sprintf('%s\\%s.mat',pat,fname));

fprintf('Analysing run %d/%d %s (%s)...',iFile,length(fileList),fname,Areas{iArea})
% load(sprintf('%s\\KDE_bins\\%s_%s_iFR50.mat',pat,fname,Areas{iArea}));
load(sprintf('%s\\KDE_binsTaskonly\\%s_%s_iFR50_behavOnly.mat',pat,fname,Areas{iArea}));

iFR_=zscore(iFR);
% iFR_=zscore(iFR_);

%% Long ordered by Choice press latency

tlimsSample = [-5 5];
tlimsChoice = [-5 5];

length_ = (sum(abs(tlimsSample))+sum(abs(tlimsChoice)))/bw;
tb_  = (1:length_)*bw;

data_{1}.iFR = iFR_;
d_LR = [t.Medium.ChoicePress_LeftCorrect;t.Medium.ChoicePress_RightCorrect] - ...
       [t.Medium.NosePoke_LeftCorrect;t.Medium.NosePoke_RightCorrect];
[~,idxLR]=sort(d_LR);
NP_offset = d_LR(idxLR)*1e-6;
NP_offset(NP_offset>abs(tlimsChoice(1)))=NaN;

times_S = [t.Medium.SamplePress_LeftCorrect;t.Medium.SamplePress_RightCorrect]*1e-6;
times_C = [t.Medium.ChoicePress_LeftCorrect;t.Medium.ChoicePress_RightCorrect]*1e-6;
[~,cutoutS]=SelTimesFRs(tlimsSample,Tmtx,data_,times_S(idxLR));
[~,cutoutC]=SelTimesFRs(tlimsChoice,Tmtx,data_,times_C(idxLR));

%%
figure; hold on

% units_ = [1 9 11 31 34];
units_ =1:length(PFCcells);
units_(avgFR<1)=[];


cmap_ = jet(length(units_));
offset=15;
for iTrial = 1:length(cutoutS{1})
    try
        FR_ = [cutoutS{1}{iTrial};cutoutC{1}{iTrial}];
        FR_=zscore(FR_);
%             m =nanmean(nanmean(FR_));e =nanstd(nanstd(FR_));
%             FR_ = (FR_-m)./e;   
for iUnit=1:length(units_)
    plot(tb_,FR_(1:length_,units_(iUnit))+offset*(iTrial-1),'color',cmap_(iUnit,:),'LineWidth',1.5)
end
        plot((sum(abs(tlimsChoice))/2+sum(abs(tlimsSample)))-NP_offset(iTrial)*[1,1],[-offset/3 +offset/3]+offset*(iTrial-1),'color',[0.9 0.4 0.2],'LineWidth',2)
    end
end

plot(sum(abs(tlimsSample))/2*[1,1],[-offset +offset*length(cutoutS{1})],'color',[0 1 0 0.6],'LineWidth',2)
plot((sum(abs(tlimsChoice))/2+sum(abs(tlimsSample)))*[1,1],[-offset +offset*length(cutoutS{1})],'color',[1 0 0 0.6],'LineWidth',2)
plot(sum(abs(tlimsSample))*[1,1],[-offset +offset*length(cutoutS{1})],'k')



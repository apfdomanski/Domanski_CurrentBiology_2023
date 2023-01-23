%% %%%%%% PREAMBLE %%%%%%
clear 
TaskOnly = true;

Epochs = {'PreSleep','Task','PostSleep'};

Epoch_ = 'PreSleep';
Target = 'LONG';

warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
    addpath(genpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\AssemblyCode\Eleonora\Mic\programs'))
elseif ismac
    pat = '/Volumes/Data/DNMTP/';    
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    addpath(genpath('/panfs/panasas01/phph/ad15419/MATLAB/Mic'))
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];    
end

fileList = dir([pat '*' Target '*.mat']);


% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={};%{'MiroslawLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileListAss(find(name_flag))=[];% fnames_(name_flag)=[];
    fileList(find(name_flag))=[];   % fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};

Targets = {'SHORT','MEDIUM','LONG'};
Delays_ = {'Short','Medium','Long'};

tlimsSample = [-5 5];
tlimsChoice = [-5 5];
bw = 0.05;
    cd ([pat 'MICmetaAnalysis'])

BinSizes=[0.005 0.01 0.0150    0.0200    0.0300    0.0500    0.0800    0.1200    0.2000    0.3500    0.5000   0.7000    1.0000];
MaxLags=[10   10   10   10   10    10    10    10    10     10     10     10     10];


%%
PreSleep    = load('PreSleepAssembliesMeta');
Task        = load('TaskAssembliesMeta');
PostSleep   = load('PostSleepAssembliesMeta');
comparison = @(A,B) (B-A)./(B+A);


for s=1:3
    n = min ([size(Task.Ass_{s}.Duration,2),size(PostSleep.Ass_{s}.Duration,2)]);
    
    subplot(2,3,1); hold on % raw bins
        
        A = PostSleep.Ass_{s}.TimeBinHist(:,1:n);
        B = Task.Ass_{s}.TimeBinHist(:,1:n);
        
        X = comparison(A,B);

        ciplot(nanmean(X,2)+nansem(X,2),nanmean(X,2)-nansem(X,2),BinSizes,color_{s},'0.5')
        % stairs(BinSizes,Ass_{s}.TimeBinHist,'Color',color_{s},'LineWidth',1.5)
        set(gca,'Xscale','log','Xlim',[min(BinSizes) max(BinSizes)],'Ylim',[-Inf Inf])
        xlabel('Binwidth (s)')
        ylabel('Task vs Post-task Sleep')
        
        
     subplot(2,3,4); hold on % bin-adjusted total pattern duration
        A = PostSleep.Ass_{s}.DurationHist(:,1:n);
        B = Task.Ass_{s}.DurationHist(:,1:n);
        X = comparison(A,B);
        
        a_ = nanmean(X,2)+nansem(X,2);
        b_ = nanmean(X,2)-nansem(X,2);
        
        x_ = logspace(-1.3,1.3,10);
        idx = isnan(a_) | isnan(b_);
        x_(idx) = [];
        a_(idx)=[];
        b_(idx)=[];
        
        ciplot(a_,b_,x_,color_{s},'0.5')
        %stairs(logspace(-1.3,1.3,10),Ass_{s}.DurationHist,'Color',color_{s},'LineWidth',1.5)
        set(gca,'Xscale','log','Ylim',[-Inf Inf])
        xlabel('Pattern duration (s)')
        ylabel('Fraction of Assemblies')
        
    subplot(2,3,2); hold on
        A = PostSleep.Ass_{s}.membercountHist(:,1:n);
        B = Task.Ass_{s}.membercountHist(:,1:n);
        X = comparison(A,B);
        
        a_ = nanmean(X,2)+nansem(X,2);
        b_ = nanmean(X,2)-nansem(X,2);
        
        x_ = 0:20;
        idx = isnan(a_) | isnan(b_);
        x_(idx) = [];
        a_(idx)=[];
        b_(idx)=[];
        
        ciplot(a_,b_,x_,color_{s},'0.5')
        set(gca,'Xlim',[0 20],'Ylim',[-Inf Inf])
        xlabel('No. interacting units')
        
    subplot(2,3,5); hold on
        A = PostSleep.Ass_{s}.membercountHistAjusted(:,1:n);
        B = Task.Ass_{s}.membercountHistAjusted(:,1:n);
        X = comparison(A,B);
        
        a_ = nanmean(X,2)+nansem(X,2);
        b_ = nanmean(X,2)-nansem(X,2);
        
        x_ = 0:0.05:0.5;
        idx = isnan(a_) | isnan(b_);
        x_(idx) = [];
        a_(idx)=[];
        b_(idx)=[];
        
        ciplot(a_,b_,x_,color_{s},'0.5')
        set(gca,'Xlim',[0 0.5],'Ylim',[-Inf Inf])
        xlabel('Fraction of recorded units interacting ')      
        
        
    subplot(2,3,3); hold on
        A = PostSleep.Ass_{s}.membersharingHist(:,1:n);
        B = Task.Ass_{s}.membersharingHist(:,1:n);
        X = comparison(A,B);
        
        a_ = nanmean(X,2)+nansem(X,2);
        b_ = nanmean(X,2)-nansem(X,2);
        
        x_ = 0:20;
        idx = isnan(a_) | isnan(b_);
        x_(idx) = [];
        a_(idx)=[];
        b_(idx)=[];
        
        ciplot(a_,b_,x_,color_{s},'0.5')
        xlabel({'Promiscuity of units between assemblies'})
        set(gca,'Xlim',[0 20],'Ylim',[-Inf Inf])        
        
    subplot(2,3,6); hold on

        A = PostSleep.Ass_{s}.membersharingHistAjusted(:,1:n);
        B = Task.Ass_{s}.membersharingHistAjusted(:,1:n);
        X = comparison(A,B);
        
        a_ = nanmean(X,2)+nansem(X,2);
        b_ = nanmean(X,2)-nansem(X,2);
        
        x_ = 0:0.05:0.5;
        idx = isnan(a_) | isnan(b_);
        x_(idx) = [];
        a_(idx)=[];
        b_(idx)=[];
        
        ciplot(a_,b_,x_,color_{s},'0.5')
        xlabel({'Promiscuity of units between assemblies';'(As Fraction of total Assemblies)'})
        set(gca,'Xlim',[0 0.5],'Ylim',[-Inf Inf])        
    
end
    
%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
DelayTimes = [4,8,16];
Target = 'LONG';

AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time

tlimsAll = [-5 5];
tlimsShort=[-4 10];
tlimsMedium=[-8 10];
tlimsLong=[-16 10];
tlimsANOVA = [-4 4];

shift = 0;
plotOnline = false;
bw=0.05;
Nbs = 500;

tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbANOVA=tlimsANOVA(1):bw:tlimsANOVA(2);
tbAll = tlimsAll(1):bw:tlimsAll(2);%0:bw:2*sum(abs(tlimsAll));

clear Av
warning ('off')
pat = 'C:\Analysis\AssemblyAnalysis\raw';
cd(pat)
fileList=dir(sprintf('allTimestamps\\*%s*.mat',Target));
fileListAss = fileList;

% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={'MiroslawLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileListAss(find(name_flag))=[];% fnames_(name_flag)=[];
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

normaliseFscores = false;
normWin = [-5 -3];
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};


%% %%%%%%%%%%%%%%%%%%%%%%%% ANALYSE UNITS %%%%%%%%%%%%%%%%%%%%%%%%
%% Batch process units
for iFile =1:length(fileList)
    fname=strtok(fileList(iFile).name,'_');
    load(sprintf('%s\\allTimestamps\\%s_Events.mat',pat,fname));
    load(sprintf('%s\\%s.mat',pat,fname));
    for iArea = 1:2%length(Areas)
        %% Get the files
        fprintf('Analysing run %d/%d %s (%s)...',iFile,length(fileList),fname,Areas{iArea})
        %         load(sprintf('%s\\KDE_bins\\%s_%s_iFR50.mat',pat,fname,Areas{iArea}));
        load(sprintf('%s\\KDE_binsTaskonly\\%s_%s_iFR50_behavOnly.mat',pat,fname,Areas{iArea}));
        
                iFR_ = iFR;
%         iFR_ = zscore(iFR);
        
        %% Decode Left/Right on correct trials
        for iDelay =1:length(Delays_)
            
            eval(sprintf('LeftTrials  = t.%s.SamplePress_LeftError'';',Delays_{iDelay}));
            eval(sprintf('RightTrials = t.%s.SamplePress_RightError'';',Delays_{iDelay}));
            
%             LeftTrials(1:(floor(0.75*length(LeftTrials))))=[];
%             RightTrials(1:(floor(0.75*length(RightTrials))))=[];
            
            
            Ltrials = cell(size(iFR,2),1);nL = 0;
            Rtrials = cell(size(iFR,2),1);nR = 0;
            
            for iTrial =1:size(LeftTrials,1)
                for iUnit =1:size(iFR,2)
                    try
                        tlims_  = LeftTrials(iTrial,1)/1e6;
                        tlims_  = closest(Tmtx,tlims_);
                        tlims_  = tlims_(1):(tlims_(1)+DelayTimes(iDelay)/bw);
                        Ltrials{iUnit} = [Ltrials{iUnit};iFR_(tlims_,iUnit)'];
                        nL=nL+1;
                    end
                end
            end
            
            
            for iTrial =1:size(RightTrials,1)
                for iUnit =1:size(iFR,2)
                    try
                        tlims_  = RightTrials(iTrial,1)/1e6;
                        tlims_  = closest(Tmtx,tlims_);
                        tlims_  = tlims_(1):(tlims_(1)+DelayTimes(iDelay)/bw);
                        Rtrials{iUnit} = [Rtrials{iUnit};iFR_(tlims_,iUnit)'];
                        nR=nR+1;
                    end
                end
            end
            
            D_.Ml  = cell2mat(cellfun(@nanmean,Ltrials,'UniformOutput',false))';
            D_.Mr  = cell2mat(cellfun(@nanmean,Rtrials,'UniformOutput',false))';
            D_.El  = cell2mat(cellfun(@nansem,Ltrials,'UniformOutput',false))';
            D_.Er  = cell2mat(cellfun(@nansem,Rtrials,'UniformOutput',false))';
            D_.FFl = (cell2mat(cellfun(@nanvar,Ltrials,'UniformOutput',false))./cell2mat(cellfun(@nanmean,Ltrials,'UniformOutput',false)))';
            D_.FFr = (cell2mat(cellfun(@nanvar,Rtrials,'UniformOutput',false))./cell2mat(cellfun(@nanmean,Rtrials,'UniformOutput',false)))';
            eval(sprintf('D{iArea}.%s = D_;',Delays_{iDelay}));
            
        end
        
        
        
        clear FF Ltrials nL Rtrials nR LeftTrials RightTrials iTrial D_
    end
    
    %% Plot Fano factors
    figure
    for iArea =1:length(Areas)-1
        subplot(length(Areas)-1,1,iArea);hold on
        
        
        for iDelay = 1:length(Delays_)
            eval(sprintf('temp = D{iArea}.%s.FFl;',Delays_{iDelay}))
%             temp = zscore(temp);
            x = (1:size(temp,1))*bw;
            for iUnit =1:size(temp,2)
                plot(x,temp(:,iUnit)+2*iUnit)
            end
        end
    end
        
        
    
 %% Plot Fano factors
    figure
    for iArea =1:length(Areas)-1
        subplot(1,length(Areas)-1,iArea);hold on
        
        title([Areas{iArea} ' units'])
        for iDelay = 1:length(Delays_)
            eval(sprintf('tempM = D{iArea}.%s.Mr;',Delays_{iDelay}))
            eval(sprintf('tempE = D{iArea}.%s.Er;',Delays_{iDelay}))
%             temp = zscore(temp);
            x = (1:size(tempM,1))*bw;
            for iUnit =1:size(tempM,2)
                ciplot(tempM(:,iUnit)+tempE(:,iUnit)+2*iUnit,tempM(:,iUnit)-tempE(:,iUnit)+2*iUnit,x,'b');
            end
        end
    end    
end




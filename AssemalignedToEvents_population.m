clear
Targets = {'SHORT','MEDIUM','LONG'};
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Dirs = {'Left','Right'};

bw=0.05;
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
    cd(pat) 
    fileList=dir('allTimestamps\*.mat');
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';    
    cd(pat) 
    fileList=dir('allTimestamps/*LONG*.mat');
end
    
% reject_list={'OnufryLONG1_Events.mat',...
%              'OnufryLONG2_Events.mat',...
%              'NorbertLONG1_Events.mat',...
%              'NorbertLONG2_Events.mat',...
%              'MiroslawLONG1_Events.mat',...
%              'MiroslawLONG2_Events.mat',...
%              'IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
reject_list={'IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
         
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end


    
% Areas = {'HP','PFC', 'HP-PFC'};
Areas = {'HP','PFC'};
p.col_ = {'r','b','g'};

tlims= [-5 10];
length_ = sum(abs([tlims]))/bw;
tb = (1:length_)*bw +tlims(1);
Events_ = {'CueLight','SamplePress','DelayEnd','NosePoke','ChoicePress','RewardConsume'};
Events__ = {'Cue Light','Sample Press','Delay End','Nose Poke','Choice Press','Reward Collect'};
%% %%%
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
    A = load(sprintf('%s%s.mat',pat,fname));
    spikeTimes{1}=A.PFCcells;
    spikeTimes{2}=A.HPcells;
    spikeTimes{3}=[spikeTimes{1};spikeTimes{2}];
    
    for iArea = 1:length(Areas)
        
        fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
        load(sprintf('%s\\KDE_binsTaskonly\\%s_%s_iFR50_behavOnly.mat',pat,fname,Areas{iArea}),'iFR','Tmtx');
        
        % Spike rates
        data_{1}.iFR = zscore(iFR);
        for iDelay = 1:length(Delays_)
            Delay_ = Delays_{iDelay};
            %%%%%%%%%%%%%%%% Correct trials
            for iEvent = 1:length(Events_)
                times_ = [eval(sprintf('t.%s.%s_LeftCorrect',Delay_,Events_{iEvent}));eval(sprintf('t.%s.%s_RightCorrect',Delay_,Events_{iEvent}))];
                [~,cutout]=SelTimesFRs(tlims,Tmtx,data_,times_*1e-6);
                
                cutout_ = cell(size(iFR,2),1);
                for iUnit=1:size(iFR,2)
                    for iTrial = 1:size(times_,1)
                        if ~isempty(cutout{1}{iTrial})
                            cutout_{iUnit}(iTrial,1:length_) = cutout{1}{iTrial}(1:length_,iUnit);
                        else
                            cutout_{iUnit}(iTrial,1:length_) = NaN;
                        end
                    end
                end
                FiringRatesMean{iEvent}{iArea}{iDelay}{iFile} = cell2mat(cellfun(@nanmean,cutout_,'UniformOutput',false))';
            end
            
            %%%%%%%%%%%%%%%% Error trials
            for iEvent = 1:length(Events_)-1
                times_ = [eval(sprintf('t.%s.%s_LeftError',Delay_,Events_{iEvent}))';eval(sprintf('t.%s.%s_RightError',Delay_,Events_{iEvent}))'];
                [~,cutout]=SelTimesFRs(tlims,Tmtx,data_,times_*1e-6);
                
                cutout_ = cell(size(iFR,2),1);
                for iUnit=1:size(iFR,2)
                    for iTrial = 1:size(times_,1)
                        if ~isempty(cutout{1}{iTrial})
                            cutout_{iUnit}(iTrial,1:length_) = cutout{1}{iTrial}(1:length_,iUnit);
                        else
                            cutout_{iUnit}(iTrial,1:length_) = NaN;
                        end
                    end
                end
%                 if length(times_)>0%size(cutout_{1},1)>1

                if size(cutout_{1},1)>1
                    FiringRatesErrorsMean{iEvent}{iArea}{iDelay}{iFile} = cell2mat(cellfun(@nanmean,cutout_,'UniformOutput',false))';
%                 elseif length(times_)==1
%                     FiringRatesErrorsMean{iEvent}{iArea}{iDelay}{iFile} = NaN(size(cell2mat(cutout_)'));
%                 elseif length(times_)==0
                else
                    FiringRatesErrorsMean{iEvent}{iArea}{iDelay}{iFile} = NaN(size(cell2mat(cutout_)'));
                end
            end
            
            
        end
    end
end

clear cutout cutout_ data_ Delay_ ERRORtrangeleft_choice ERRORtrangeleft_sample ERRORtrangeright_sample ERRORtrangeright_choice
clear HPcells HPinter HPtonic iArea iDelay iEvent iFile iFR inputData iTrial iUnit PFCcells PFCinter PFCtonic t times_ Tmtx trangeleft_choice trangeleft_sample trangeright_choice trangeright_sample

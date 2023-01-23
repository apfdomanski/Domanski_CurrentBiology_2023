%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Target = 'LONG';

shift = 0;
plotOnline = false;
bw=0.05;
Nbs = 500;

warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
else
    pat = '/Volumes/HDD2/DNMTP/raw/';
end
cd(pat)
fileList=dir(sprintf('MixedSelectivity_LONG%s*%s*_MixedSelectivity_MembershipsortedUnitsEventAnalysis_Fano.mat',filesep,Target));
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};
MemberClasses ={'LocalMembers','JointMembers','NonMembers'};
MemberClasses_ ={'Local Assembly Members','Inter-area Assembly Members','Non-members'};
col_ ={[0.9 0.6 0],[0.6 0 0.6],[0.6 0.6 0.6]};
normaliseFscores = false;
Outcome = {'Correct','Error'};
Events_ = {'CueLight','SamplePress','DelayEnd','NosePoke','ChoicePress','RewardConsume'};
Events__ = {'Cue Light','Sample Press','Delay End','Nose Poke','Choice Press','Reward Collection'};
no_Events_ = [6,5];
normWin = [-2 0];

%% Gather Files - delay means
clear D_ D;
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    fnIn = sprintf('%sMixedSelectivity_LONG%s%s_MixedSelectivity_MembershipsortedUnitsEventAnalysis_Fano.mat',pat,filesep,fname);
    load(fnIn ,'D','tbAll');
    
    % D{iOutcome}{iEvent}{iDelay}{iArea}{iClass} = D_;
    
    for iOutcome = 1
        for iArea =1:2
            for iEvent = 1:no_Events_(iOutcome)
                
                for iClass= 1:length(MemberClasses)
                    x{iClass} = nan(length(tbAll),length(Delays_));
                    xFano{iClass} = nan(length(tbAll),length(Delays_));
                    xMean{iClass} = nan(length(tbAll),length(Delays_));
                    for iDelay = 1:length(Delays_)
                        try
                            x{iClass}(:,iDelay) = D{iOutcome}{iEvent}{iDelay}{iArea}{iClass}.Ft2;
                        end
                        try
                            xFano{iClass}(:,iDelay) = nanmean(D{iOutcome}{iEvent}{iDelay}{iArea}{iClass}.Fano,2);
                        end
                        try
                            temp = D{iOutcome}{iEvent}{iDelay}{iArea}{iClass}.Mean;
%                             idx = find(tbAll>=normWin(1) & tbAll<=normWin(2));                      
%                             tempMean = repmat(nanmean(temp(idx,:),1),length(tbAll),1);
%                             temp=temp./tempMean;
% %                             temp=zscore(temp);
                            xMean{iClass}(:,iDelay) = nanmean(temp,2);
                        end
                    end
                    
                end
                D_{iArea}{iFile}{iEvent} = cell2mat(cellfun(@(X) nanmean(X,2) ,x,'UniformOutput',false));
                Fano_{iArea}{iFile}{iEvent} = cell2mat(cellfun(@(X) nanmean(X,2) ,xFano,'UniformOutput',false));
                Mean_{iArea}{iFile}{iEvent} = cell2mat(cellfun(@(X) nanmean(X,2) ,xMean,'UniformOutput',false));
                
            end
            
            
        end
        
    end
end
%% Collapse and plot Decoding
for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        for iEvent = 1:no_Events_(iOutcome)
            D_Collapse{iArea}{iEvent}{iClass}= [];
            for iFile = 1:length(fileList)
                D_Collapse{iArea}{iEvent}{iClass} = [D_Collapse{iArea}{iEvent}{iClass} , D_{iArea}{iFile}{iEvent}(:,iClass)];
            end
            D_CollapseMean{iArea}{iEvent} = cell2mat(cellfun(@(x)nanmean(x,2),D_Collapse{iArea}{iEvent},'UniformOutput',false));
            D_CollapseSEM{iArea}{iEvent}  = cell2mat(cellfun(@(x)nansem(x,2),D_Collapse{iArea}{iEvent},'UniformOutput',false));
        end
    end
end

for iArea = 1:2
    figure('name',Areas{iArea})
    for iEvent = 1:no_Events_(iOutcome)
        subplot(1,no_Events_(iOutcome),iEvent); hold on
        for iClass = 1:length(MemberClasses)
           m = D_CollapseMean{iArea}{iEvent}(:,iClass);
           e = D_CollapseSEM{iArea}{iEvent}(:,iClass);
           ciplot(m+e,m-e,tbAll,col_{iClass},1) 
           plot([0 0],[0 20],':k')
           axis([min(tbAll),max(tbAll),0 10])
           title(Events__{iEvent})
        end

    end
end
%% Collapse and plot Fano

for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        for iEvent = 1:no_Events_(iOutcome)
            Fano_Collapse{iArea}{iEvent}{iClass}= [];
            for iFile = 1:length(fileList)
                Fano_Collapse{iArea}{iEvent}{iClass} = [Fano_Collapse{iArea}{iEvent}{iClass} , Fano_{iArea}{iFile}{iEvent}(:,iClass)];
            end
            Fano_CollapseMean{iArea}{iEvent} = cell2mat(cellfun(@(x)nanmean(x,2),Fano_Collapse{iArea}{iEvent},'UniformOutput',false));
            Fano_CollapseSEM{iArea}{iEvent}  = cell2mat(cellfun(@(x)nansem(x,2),Fano_Collapse{iArea}{iEvent},'UniformOutput',false));
        end
    end
end

for iArea = 1:2
    figure('name',Areas{iArea})
    for iEvent = 1:no_Events_(iOutcome)
        subplot(1,no_Events_(iOutcome),iEvent); hold on
        for iClass = 1:length(MemberClasses)
           m = Fano_CollapseMean{iArea}{iEvent}(:,iClass);
           e = Fano_CollapseSEM{iArea}{iEvent}(:,iClass);
           ciplot(m+e,m-e,tbAll,col_{iClass},1) 
           plot([0 0],[0 20],':k')
           axis([min(tbAll),max(tbAll),0 10])
           title(Events__{iEvent})
        end

    end
end
%% Collapse and plot Mean

for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        for iEvent = 1:no_Events_(iOutcome)
            Mean_Collapse{iArea}{iEvent}{iClass}= [];
            for iFile = 1:length(fileList)
                Mean_Collapse{iArea}{iEvent}{iClass} = [Mean_Collapse{iArea}{iEvent}{iClass} , Mean_{iArea}{iFile}{iEvent}(:,iClass)];
            end
            Mean_CollapseMean{iArea}{iEvent} = cell2mat(cellfun(@(x)nanmean(x,2),Mean_Collapse{iArea}{iEvent},'UniformOutput',false));
            Mean_CollapseSEM{iArea}{iEvent}  = cell2mat(cellfun(@(x)nansem(x,2),Mean_Collapse{iArea}{iEvent},'UniformOutput',false));
        end
    end
end

for iArea = 1:2
    figure('name',Areas{iArea})
    for iEvent = 1:no_Events_(iOutcome)
        subplot(1,no_Events_(iOutcome),iEvent); hold on
        for iClass = 1:length(MemberClasses)
           m = Mean_CollapseMean{iArea}{iEvent}(:,iClass);
           e = Mean_CollapseSEM{iArea}{iEvent}(:,iClass);
           ciplot(m+e,m-e,tbAll,col_{iClass},1) 
           plot([0 0],[0 20],':k')
           axis([min(tbAll),max(tbAll),-0.5 2])
           title(Events__{iEvent})
        end

    end
end

%% Gather Files - collapsed
clear D_ Fano_ D Mean_ x xMean
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    fnIn = sprintf('%sMixedSelectivity_LONG%s%s_MixedSelectivity_MembershipsortedUnitsEventAnalysis_Fano.mat',pat,filesep,fname);
    load(fnIn ,'D','tbAll');
    
    % D{iOutcome}{iEvent}{iDelay}{iArea}{iClass} = D_;
    
    for iOutcome = 1
        for iArea =1:2
            for iEvent = 1:no_Events_(iOutcome)
                
                for iClass= 1:length(MemberClasses)
                    D_{iArea}{iFile}{iEvent}{iClass}    = nan(length(tbAll),1); %%%%
                    Fano_{iArea}{iFile}{iEvent}{iClass} = nan(length(tbAll),1); %%%%
                    Mean_{iArea}{iFile}{iEvent}{iClass} = nan(length(tbAll),1); %%%%
                    
                    x{iClass} = nan(length(tbAll),length(Delays_));
                    xFano{iClass} = nan(length(tbAll),length(Delays_));
                    xMean{iClass} = nan(length(tbAll),length(Delays_));
                    for iDelay = 1:length(Delays_)
                        try
                            x{iClass}(:,iDelay) = D{iOutcome}{iEvent}{iDelay}{iArea}{iClass}.Ft2;                
                        end
                        try
                            xFano{iClass}(:,iDelay) = nanmean(D{iOutcome}{iEvent}{iDelay}{iArea}{iClass}.Fano,2);
                        end
                        try
                            temp = D{iOutcome}{iEvent}{iDelay}{iArea}{iClass}.Mean;
%                             idx = find(tbAll>=normWin(1) & tbAll<=normWin(2));                      
%                             tempMean = repmat(nanmean(temp(idx,:),1),length(tbAll),1);
%                             temp=temp./tempMean;
% %                             temp=zscore(temp);
                            xMean{iClass}(:,iDelay) = nanmean(temp,2);
                            clear temp
                            
                        end
                    end
                    D_{iArea}{iEvent}{iClass}{iFile} = x{iClass};
                    Fano_{iArea}{iEvent}{iClass}{iFile} = xFano{iClass};
                    Mean_{iArea}{iEvent}{iClass}{iFile} = xMean{iClass};

                end
                %D_{iArea}{iFile}{iEvent}    = cell2mat(cellfun(@(X) nanmean(X,2) ,x,'UniformOutput',false));
                %Fano_{iArea}{iFile}{iEvent} = cell2mat(cellfun(@(X) nanmean(X,2) ,xFano,'UniformOutput',false));
                
%                 Fano_{iArea}{iFile}{iEvent} = cell2mat(cellfun(@(X) nanmean(X,2) ,xFano,'UniformOutput',false));
            end
            
        end
        
    end
end
%% Collapse and plot Decoding
for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        for iEvent = 1:no_Events_(iOutcome)
            D_Collapse{iArea}{iEvent}{iClass} = cell2mat(D_{iArea}{iEvent}{iClass});
            D_CollapseMean{iArea}{iEvent} = cell2mat(cellfun(@(x)nanmean(x,2),D_Collapse{iArea}{iEvent},'UniformOutput',false));
            D_CollapseSEM{iArea}{iEvent}  = cell2mat(cellfun(@(x)nansem(x,2),D_Collapse{iArea}{iEvent},'UniformOutput',false));
        end
    end
end

for iArea = 1:2
    figure('name',Areas{iArea})
    for iEvent = 1:no_Events_(iOutcome)
        subplot(1,no_Events_(iOutcome),iEvent); hold on
        for iClass = 1:length(MemberClasses)
           m = D_CollapseMean{iArea}{iEvent}(:,iClass);
           e = D_CollapseSEM{iArea}{iEvent}(:,iClass);
           ciplot(m+e,m-e,tbAll,col_{iClass},1) 
           plot([0 0],[0 20],':k')
           axis([min(tbAll),max(tbAll),0 Inf])
           title(Events__{iEvent})
        end

    end
end
%% Collapse and plot Fano

for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        for iEvent = 1:no_Events_(iOutcome)
            Fano_Collapse{iArea}{iEvent}{iClass} = cell2mat(Fano_{iArea}{iEvent}{iClass});
            Fano_CollapseMean{iArea}{iEvent} = cell2mat(cellfun(@(x)nanmean(x,2),Fano_Collapse{iArea}{iEvent},'UniformOutput',false));
            Fano_CollapseSEM{iArea}{iEvent}  = cell2mat(cellfun(@(x)nansem(x,2),Fano_Collapse{iArea}{iEvent},'UniformOutput',false));
        end
    end
end

for iArea = 1:2
    figure('name',Areas{iArea})
    for iEvent = 1:no_Events_(iOutcome)
        subplot(1,no_Events_(iOutcome),iEvent); hold on
        for iClass = 1:length(MemberClasses)
           m = Fano_CollapseMean{iArea}{iEvent}(:,iClass);
           e = Fano_CollapseSEM{iArea}{iEvent}(:,iClass);
           ciplot(m+e,m-e,tbAll,col_{iClass},1) 
           plot([0 0],[0 20],':k')
           axis([min(tbAll),max(tbAll),0 10])
           title(Events__{iEvent})
        end

    end
end
%% Collapse and plot Mean

for iArea = 1:2
    for iClass = 1:length(MemberClasses)
        for iEvent = 1:no_Events_(iOutcome)    
            x = cell2mat(Mean_{iArea}{iEvent}{iClass});
            x(isinf(x))=NaN;
            Mean_Collapse{iArea}{iEvent}{iClass} =x;
            Mean_CollapseMean{iArea}{iEvent} = cell2mat(cellfun(@(x)nanmean(x,2),Mean_Collapse{iArea}{iEvent},'UniformOutput',false));
            Mean_CollapseSEM{iArea}{iEvent}  = cell2mat(cellfun(@(x)nansem(x,2),Mean_Collapse{iArea}{iEvent},'UniformOutput',false));
        end
    end
end

for iArea = 1:2
    figure('name',Areas{iArea})
    for iEvent = 1:no_Events_(iOutcome)
        subplot(1,no_Events_(iOutcome),iEvent); hold on
        for iClass = 1:length(MemberClasses)
           m = Mean_CollapseMean{iArea}{iEvent}(:,iClass);
           e = Mean_CollapseSEM{iArea}{iEvent}(:,iClass);
           ciplot(m+e,m-e,tbAll,col_{iClass},1) 
           plot([0 0],[-1 1],':k')
           axis([min(tbAll),max(tbAll),-1 1])
           title(Events__{iEvent})
        end

    end
end




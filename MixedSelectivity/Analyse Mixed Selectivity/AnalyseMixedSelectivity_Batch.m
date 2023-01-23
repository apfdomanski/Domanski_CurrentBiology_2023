function AnalyseMixedSelectivity_Batch()
%% %%%%%% PREAMBLE %%%%%%
clear
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

mkdir([pat 'MixedSelectivity'])
TargetList = {'SHORT','MEDIUM','LONG'};

AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time
for iTarget=1:length(TargetList)
    switch AssemblyChoice
        case 1
            pat2{iTarget} = [pat 'KDE_bins' filesep TargetList{iTarget} filesep];
        case 2
            pat2{iTarget} = [pat 'KDE_binsTaskonly' filesep TargetList{iTarget} 'Taskonly' filesep];
        case 3
            pat2{iTarget} = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
    end
end
DelaysList_ = {{'Delay_0'};{'Delay_2','Delay_4','Delay_6','Delay_8'};{'Short','Medium','Long'}};
DelaysList__ = {{'0s'};{'2s','4s','6s','8s'};{'4s','8s','16s'}};
reject_listList={{''};...
    {'NorbertMEDIUM1_Events.mat','NorbertMEDIUM2_Events.mat','OnufryMEDIUM1_Events.mat','OnufryMEDIUM2_Events.mat'};...
    {'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}};

Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};
bw = 0.05;
tlimsANOVA = [-2 2];
tbANOVA = tlimsANOVA(1):bw:tlimsANOVA(2);
shift = 0;

% Re-sample a subset of trials to keep consistent between experiments?
resampleTrials = true;
nDraws = 1000;
nDrawnTrials = 5;

ProcessUnits = false;
ProcessAssemblies = true;
PlotOnline = false;
%% Batch analyse
for iTarget =1:length(TargetList)
    
    Target = TargetList{iTarget};
    
    Delays_     = DelaysList_{iTarget};
    Delays__    = DelaysList__{iTarget};
    reject_list = reject_listList{iTarget};
    
    fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);
    fileListAss = fileList;
    
    % Strip troublemakers
    name_LowTrialFlag=zeros(numel(fileList),1);
    for idx=1:numel(fileList)
        fnames_{idx,1}=fileList(idx).name;
        name_LowTrialFlag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical LowTrialFlag for testing
    end
    try
        fileListAss(find(name_LowTrialFlag))=[];% fnames_(name_LowTrialFlag)=[];
        fileList(find(name_LowTrialFlag))=[];% fnames_(name_LowTrialFlag)=[];
    end
    clear  reject_list idx fnames_ name_LowTrialFlag
    
    %% Batch process units
    if ProcessUnits
        for iFile =1:length(fileList)
            fname=strtok(fileList(iFile).name,'_');
            load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
            load(sprintf('%s%s%s.mat',pat,filesep,fname));
            
            for iArea = 1:2
                %% Get the files
                fprintf('Analysing run %d/%d %s (%s)...\n',iFile,length(fileList),fname,Areas{iArea})
                %         load(sprintf('%sKDE_bins%s%s_%s_iFR50.mat',pat,filesep,fname,Areas{iArea}));
                load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep, fname,Areas{iArea}));
                
                iFR_ = iFR;
                % iFR_ = zscore(iFR_);
                %% ANOVA on L/R S/C C/E
                q_names = {'Untuned','Pure only','Mixed only','Mixed and Pure','Any tuning'};
                q2_names = {'Untuned','CS','LMS','NMS',...
                    'LMS CxO','LMS CxP','LMS OxP','LMS CxOxP',...
                    'NMS CxO','NMS CxP','NMS OxP','NMS CxOxP'};
                
                for iDelay =1:length(Delays_)
                    eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
                    Trials = [t_.SamplePress_LeftCorrect;...
                        t_.SamplePress_LeftError';...
                        t_.SamplePress_RightCorrect;...
                        t_.SamplePress_RightError';...
                        t_.ChoicePress_LeftCorrect;...
                        t_.ChoicePress_LeftError';...
                        t_.ChoicePress_RightCorrect;...
                        t_.ChoicePress_RightError'];
                    SC = [repmat({'Sample'},length(t_.SamplePress_LeftCorrect),1);...
                        repmat({'Sample'},length(t_.SamplePress_LeftError),1);...
                        repmat({'Sample'},length(t_.SamplePress_RightCorrect),1);...
                        repmat({'Sample'},length(t_.SamplePress_RightError),1);...
                        repmat({'Choice'},length(t_.ChoicePress_LeftCorrect),1);...
                        repmat({'Choice'},length(t_.ChoicePress_LeftError),1);...
                        repmat({'Choice'},length(t_.ChoicePress_RightCorrect),1);...
                        repmat({'Choice'},length(t_.ChoicePress_RightError),1)];
                    
                    CE = [repmat({'Correct'},length(t_.SamplePress_LeftCorrect),1);...
                        repmat({'Error'},length(t_.SamplePress_LeftError),1);...
                        repmat({'Correct'},length(t_.SamplePress_RightCorrect),1);...
                        repmat({'Error'},length(t_.SamplePress_RightError),1);...
                        repmat({'Correct'},length(t_.ChoicePress_LeftCorrect),1);...
                        repmat({'Error'},length(t_.ChoicePress_LeftError),1);...
                        repmat({'Correct'},length(t_.ChoicePress_RightCorrect),1);...
                        repmat({'Error'},length(t_.ChoicePress_RightError),1)];
                    
                    LR = [repmat({'Left'},length(t_.SamplePress_LeftCorrect),1);...
                        repmat({'Left'},length(t_.SamplePress_LeftError),1);...
                        repmat({'Right'},length(t_.SamplePress_RightCorrect),1);...
                        repmat({'Right'},length(t_.SamplePress_RightError),1);...
                        repmat({'Left'},length(t_.ChoicePress_LeftCorrect),1);...
                        repmat({'Left'},length(t_.ChoicePress_LeftError),1);...
                        repmat({'Right'},length(t_.ChoicePress_RightCorrect),1);...
                        repmat({'Right'},length(t_.ChoicePress_RightError),1)];
                    FR =[];
                    for iTrial =1:length(Trials)
                        try
                            tlims_ =Trials(iTrial)/1e6+tlimsANOVA + shift;
                            tlims_ = closest(Tmtx,tlims_);
                            FR(iTrial,1:size(iFR_,2))= mean(iFR_(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
                        catch
                            FR(iTrial,1:size(iFR_,2))= nan(1,size(iFR_,2));
                        end
                    end
                    SC(isnan(sum(FR,2)))=[];
                    CE(isnan(sum(FR,2)))=[];
                    LR(isnan(sum(FR,2)))=[];
                    FR(isnan(sum(FR,2)),:)=[];
                    p_ = zeros(size(iFR_,2),7); F_=p_;
                    
                    
                    for iUnit=1:size(iFR_,2)
                        [p_(iUnit,:),tbl_] = anovan(FR(:,iUnit),{SC CE LR},'model','full',...
                            'varnames',{'Context','Outcome','Position'},'display','off');
                        F_(iUnit,:)=[tbl_{2:8,6}];
                        % [results,~,~,gnames] = multcompare(stats,'Dimension',[1 2 3])
                    end
                    p_thresh = p_<0.05;
                    
                    %Fractions of units with different tunings:
                    q_ = [sum(sum(p_thresh,2)==0),...                                      % No tuning
                        sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)==0),...   % Pure tuning only
                        sum(sum(p_thresh(:,1:3),2)==0 & sum(p_thresh(:,4:7),2)>0),...    % Mixed tuning only
                        sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)>0),...    % Mixed and pure tuning
                        sum(sum(p_thresh,2)>0)];                                          % Any tuning
                    
                    q2_= [sum(sum(p_thresh,2)==0),...                                       % No tuning
                        sum(sum(p_thresh(:,1:3),2)==1 & sum(p_thresh(:,4:7),2)==0),...    % CS:  Classical selectivity for one variable
                        sum(sum(p_thresh(:,1:3),2)>1  & sum(p_thresh(:,4:7),2)==0),...    % LMS: Linear mixed selectivity for 2+ variables
                        sum(sum(p_thresh(:,4:7),2)>1),...                                      % NMS: Non-linear mixed selectivity
                        sum(sum(p_thresh == [1 0 0 0 0 0 0],2)==7),...     CS for Context
                        sum(sum(p_thresh == [0 1 0 0 0 0 0],2)==7),...     CS for Outcome
                        sum(sum(p_thresh == [0 0 1 0 0 0 0],2)==7),...     CS for Position
                        sum(sum(p_thresh == [1 1 0 0 0 0 0],2)==7),...     LMS for Context X Outcome
                        sum(sum(p_thresh == [1 0 1 0 0 0 0],2)==7),...     LMS for Context X Position
                        sum(sum(p_thresh == [0 1 1 0 0 0 0],2)==7),...     LMS for Outcome X Position
                        sum(sum(p_thresh == [1 1 1 0 0 0 0],2)==7),...     LMS for Context X Outcome X Position
                        sum(sum(p_thresh(:,[4 7])==[1 0],2)==2),...        NMS for Context X Outcome
                        sum(sum(p_thresh(:,[5 7])==[1 0],2)==2),...        NMS for Context X Position
                        sum(sum(p_thresh(:,[6 7])==[1 0],2)==2),...        NMS for Outcome X Position
                        sum(sum(p_thresh(:,7) == 1 ,2)==2),...             NMS for Outcome X Position X Position
                        ];
                    %             q_ = [sum((sum(p_'<0.05)<1)), ...        % No tuning
                    %                 sum((sum(p_(:,1:3)'<0.05)>1)), ... % Pure tuning only
                    %                 sum((sum(p_(:,4:7)'<0.05)>1))];    % Mixed tuning
                    eval(sprintf('D.%s.ANOVA3.p=p_;',Delays_{iDelay}));
                    eval(sprintf('D.%s.ANOVA3.F=F_;',Delays_{iDelay}));
                    eval(sprintf('D.%s.ANOVA3.Factors=transpose({tbl_{2:8,1}});',Delays_{iDelay}));
                    eval(sprintf('D.%s.ANOVA3.prcSig=sum(p_<0.05)./size(iFR_,2)*100;',Delays_{iDelay}));
                    eval(sprintf('D.%s.ANOVA3.prcTuned=q_./size(iFR_,2)*100;',Delays_{iDelay}));
                    eval(sprintf('D.%s.ANOVA3.prcTuned2=q2_./size(iFR_,2)*100;',Delays_{iDelay}));
                    eval(sprintf('D.%s.ANOVA3.tuning=q_names;',Delays_{iDelay}));
                    eval(sprintf('D.%s.ANOVA3.tuning2=q2_names;',Delays_{iDelay}));
                    
                end
                
                clear p_ q_ q2_ q_names q2_names FR SC LR CE p_thresh
                %% ANOVA on L/R S/C correct
                q_names = {'Untuned','CS (Context)','CS (Location)','Any CS','LMS','NMS','Any Tuning'};
                for iDelay =1:length(Delays_)
                    
                    eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
                    
                    if length(t_.SamplePress_LeftCorrect)>=nDrawnTrials && length(t_.SamplePress_RightCorrect)>=nDrawnTrials
                        Trials = [t_.SamplePress_LeftCorrect;...
                            t_.ChoicePress_LeftCorrect;...
                            t_.SamplePress_RightCorrect;...
                            t_.ChoicePress_RightCorrect];
                        
                        SC = [repmat({'Sample'},length(t_.SamplePress_LeftCorrect),1);...
                            repmat({'Choice'},length(t_.ChoicePress_LeftCorrect),1);...
                            repmat({'Sample'},length(t_.SamplePress_RightCorrect),1);...
                            repmat({'Choice'},length(t_.ChoicePress_RightCorrect),1)];
                        
                        
                        LR = [repmat({'Left'},length(t_.SamplePress_LeftCorrect),1);...
                            repmat({'Left'},length(t_.ChoicePress_LeftCorrect),1);...
                            repmat({'Right'},length(t_.SamplePress_RightCorrect),1);...
                            repmat({'Right'},length(t_.ChoicePress_RightCorrect),1)];
                        
                        nL = length(t_.SamplePress_LeftCorrect);
                        nR = length(t_.SamplePress_RightCorrect);
                        
                        FR = ExtractFiringRates(Trials,tlimsANOVA,Tmtx,iFR_,tbANOVA);
                        
                        Stats_ = CalculateANOVA(FR,LR,SC,nDraws,nDrawnTrials,PlotOnline,resampleTrials);
                        %% Collapse results and write out
                        if resampleTrials
                            p_       = Stats_.p_DrawnMetaP;
                            p_thresh = Stats_.p_DrawnThresh;
                            F_       = Stats_.F_DrawnMetaP;
                        else
                            p_       = Stats_.p_Total;
                            p_thresh = Stats_.p_TotalThresh;
                            F_       = Stats_.F_Total;
                        end
                        %Fractions of units with different tunings:
                        q_ = [sum(sum(p_thresh,2)==0),...              % No tuning
                            sum(sum(p_thresh == [1 0 0],2)==3),...     % CS for Context
                            sum(sum(p_thresh == [0 1 0],2)==3),...     % CS for Location
                            sum(sum(p_thresh == [1 0 0],2)==3) + sum(sum(p_thresh == [0 1 0],2)==3),... % Any CS
                            sum(sum(p_thresh == [1 1 0],2)==3),...     % LMS for Context x Location
                            sum(sum(p_thresh(:,3) == 1,2)),...         % NMS for Context x Location
                            sum(sum(p_thresh,2)>0)];                   % Any tuning
                        
                        eval(sprintf('D.%s.ANOVA2.p_cor=p_;',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.F_cor=F_;',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.prcSig_cor=sum(p_<0.05)./size(iFR_,2)*100;',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.prcTuned_cor=q_./size(iFR_,2)*100;',Delays_{iDelay}));
                        
                    else
                        eval(sprintf('D.%s.ANOVA2.p_cor=nan(1,3);',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.F_cor=nan(1,3);',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.prcSig_cor=nan(1,3);',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.prcTuned_cor=nan(1,length(q_names));',Delays_{iDelay}));
                    end
                    eval(sprintf('D.%s.ANOVA2.tuning_cor=q_names;',Delays_{iDelay}));
                    eval(sprintf('D.%s.ANOVA2.Factors_cor={''Context'',''Position'',''Context*Position''};',Delays_{iDelay}));
                    
                end
                clear p_ q_ q_names  FR SC LR CE p_thresh
                %% ANOVA on L/R S/C error
                q_names = {'Untuned','CS (Context)','CS (Location)','Any CS','LMS','NMS','Any Tuning'};
                for iDelay =1:length(Delays_)
                    
                    eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
                    
                    if length(t_.SamplePress_LeftError)>=nDrawnTrials && length(t_.SamplePress_RightError)>=nDrawnTrials
                        Trials = [t_.SamplePress_LeftError';...
                            t_.SamplePress_RightError';...
                            t_.ChoicePress_LeftError';...
                            t_.ChoicePress_RightError'];
                        SC = [repmat({'Sample'},length(t_.SamplePress_LeftError),1);...
                            repmat({'Sample'},length(t_.SamplePress_RightError),1);...
                            repmat({'Choice'},length(t_.ChoicePress_LeftError),1);...
                            repmat({'Choice'},length(t_.ChoicePress_RightError),1)];
                        
                        LR = [repmat({'Left'},length(t_.SamplePress_LeftError),1);...
                            repmat({'Right'},length(t_.SamplePress_RightError),1);...
                            repmat({'Left'},length(t_.ChoicePress_LeftError),1);...
                            repmat({'Right'},length(t_.ChoicePress_RightError),1)];
                        
                        FR = ExtractFiringRates(Trials,tlimsANOVA,Tmtx,iFR_,tbANOVA);
                        
                        Stats_ = CalculateANOVA(FR,LR,SC,nDraws,nDrawnTrials,PlotOnline,resampleTrials);
                        %% Collapse results and write out
                        if resampleTrials
                            p_       = Stats_.p_DrawnMetaP;
                            p_thresh = Stats_.p_DrawnThresh;
                            F_       = Stats_.F_DrawnMetaP;
                        else
                            p_       = Stats_.p_Total;
                            p_thresh = Stats_.p_TotalThresh;
                            F_       = Stats_.F_Total;
                        end
                        %Fractions of units with different tunings:
                        q_ = [sum(sum(p_thresh,2)==0),...              % No tuning
                            sum(sum(p_thresh == [1 0 0],2)==3),...     % CS for Context
                            sum(sum(p_thresh == [0 1 0],2)==3),...     % CS for Location
                            sum(sum(p_thresh == [1 0 0],2)==3) + sum(sum(p_thresh == [0 1 0],2)==3),... % Any CS
                            sum(sum(p_thresh == [1 1 0],2)==3),...     % LMS for Context x Location
                            sum(sum(p_thresh(:,3) == 1,2)),...         % NMS for Context x Location
                            sum(sum(p_thresh,2)>0)];                   % Any tuning
                        
                        eval(sprintf('D.%s.ANOVA2.p_err=p_;',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.F_err=F_;',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.prcSig_err=sum(p_<0.05)./size(iFR_,2)*100;',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.prcTuned_err=q_./size(iFR_,2)*100;',Delays_{iDelay}));
                        
                    else
                        eval(sprintf('D.%s.ANOVA2.p_err=nan(1,3);',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.F_err=nan(1,3);',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.prcSig_err=nan(1,3);',Delays_{iDelay}));
                        eval(sprintf('D.%s.ANOVA2.prcTuned_err=nan(1,length(q_names));',Delays_{iDelay}));
                    end
                    eval(sprintf('D.%s.ANOVA2.tuning_err=q_names;',Delays_{iDelay}));
                    eval(sprintf('D.%s.ANOVA2.Factors_err={''Context'',''Position'',''Context*Position''};',Delays_{iDelay}));
                    
                end
                clear p_ q_ q_names  FR SC LR CE p_thresh
                %% Save results
                if resampleTrials
                    fnOut = sprintf('%sMixedSelectivity%sTrialNosMatched%s%s_%s_MixedSelectivity_Units.mat',pat,filesep,filesep,fname,Areas{iArea});
                else
                    fnOut = sprintf('%sMixedSelectivity%s%s_%s_MixedSelectivity_Units.mat',pat,filesep,fname,Areas{iArea});
                end
                save(fnOut,'D','tlimsANOVA','tbANOVA','-v7.3');
                
                fprintf('Done.\n')
            end
        end
    end
    %% Batch process assemblies
    if ProcessAssemblies
        for iFile =1:length(fileListAss)
            %% Get the files
            fname=strtok(fileListAss(iFile).name,'_');
            load(sprintf('%sallTimestamps%s%s_Events.mat',pat,filesep,fname));
            fprintf('Analysing run %d/%d %s ...\n',iFile,length(fileList),fname)
            switch AssemblyChoice
                case 1
                    Ass = load(sprintf('%s%s_iFR50_FSC.mat',pat2{iTarget},fname));
                case 2
                    Ass = load(sprintf('%s%s_iFR50_BehavOnly_FSC.mat',pat2{iTarget},fname));
                case 3
                    Ass = load(sprintf('%s%s_iFR50_Task_FSC.mat',pat2{iTarget},fname));
            end
            for iArea = 1:length(Ass.FSCsel)
                if ~isempty(Ass.FSCsel{iArea})
                    %% ANOVA on L/R S/C C/E
                    Factors = {'Context';'Outcome';'Position';'Context*Outcome';'Context*Position';'Outcome*Position';'Context*Outcome*Position'};
                    q_names = {'Untuned','Pure only','Mixed only','Mixed and Pure','Any tuning'};
                    q2_names = {'Untuned','CS','LMS','NMS',...
                        'LMS CxO','LMS CxP','LMS OxP','LMS CxOxP',...
                        'NMS CxO','NMS CxP','NMS OxP','NMS CxOxP'};
                    for iDelay =1:length(Delays_)
                        eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
                        Trials = [t_.SamplePress_LeftCorrect;...
                            t_.SamplePress_LeftError';...
                            t_.SamplePress_RightCorrect;...
                            t_.SamplePress_RightError';...
                            t_.ChoicePress_LeftCorrect;...
                            t_.ChoicePress_LeftError';...
                            t_.ChoicePress_RightCorrect;...
                            t_.ChoicePress_RightError'];
                        
                        SC = [repmat({'Sample'},length(t_.SamplePress_LeftCorrect),1);...
                            repmat({'Sample'},length(t_.SamplePress_LeftError),1);...
                            repmat({'Sample'},length(t_.SamplePress_RightCorrect),1);...
                            repmat({'Sample'},length(t_.SamplePress_RightError),1);...
                            repmat({'Choice'},length(t_.ChoicePress_LeftCorrect),1);...
                            repmat({'Choice'},length(t_.ChoicePress_LeftError),1);...
                            repmat({'Choice'},length(t_.ChoicePress_RightCorrect),1);...
                            repmat({'Choice'},length(t_.ChoicePress_RightError),1)];
                        
                        CE = [repmat({'Correct'},length(t_.SamplePress_LeftCorrect),1);...
                            repmat({'Error'},length(t_.SamplePress_LeftError),1);...
                            repmat({'Correct'},length(t_.SamplePress_RightCorrect),1);...
                            repmat({'Error'},length(t_.SamplePress_RightError),1);...
                            repmat({'Correct'},length(t_.ChoicePress_LeftCorrect),1);...
                            repmat({'Error'},length(t_.ChoicePress_LeftError),1);...
                            repmat({'Correct'},length(t_.ChoicePress_RightCorrect),1);...
                            repmat({'Error'},length(t_.ChoicePress_RightError),1)];
                        
                        LR = [repmat({'Left'},length(t_.SamplePress_LeftCorrect),1);...
                            repmat({'Left'},length(t_.SamplePress_LeftError),1);...
                            repmat({'Right'},length(t_.SamplePress_RightCorrect),1);...
                            repmat({'Right'},length(t_.SamplePress_RightError),1);...
                            repmat({'Left'},length(t_.ChoicePress_LeftCorrect),1);...
                            repmat({'Left'},length(t_.ChoicePress_LeftError),1);...
                            repmat({'Right'},length(t_.ChoicePress_RightCorrect),1);...
                            repmat({'Right'},length(t_.ChoicePress_RightError),1)];
                        FSC = [];
                        nAss = size(Ass.FSCsel{iArea},2);
                        for iTrial =1:length(Trials)
                            try
                                tlims_ = Trials(iTrial)/1e6+tlimsANOVA + shift;
                                tlims_ = closest(Ass.Tmtx,tlims_);
                                FSC(iTrial,1:nAss) = mean(Ass.FSCsel{iArea}(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
                            catch
                                FSC(iTrial,1:nAss)= nan(1,nAss);
                            end
                        end
                        SC(isnan(sum(FSC,2)))=[];
                        CE(isnan(sum(FSC,2)))=[];
                        LR(isnan(sum(FSC,2)))=[];
                        FSC(isnan(sum(FSC,2)),:)=[];
                        p_ = zeros(nAss,7); F_=p_;
                        
                        for iAss=1:nAss
                            [p_(iAss,:),tbl_] = anovan(FSC(:,iAss),{SC CE LR},'model','full',...
                                'varnames',{'Context','Outcome','Position'},'display','off');
                            F_(iAss,:)=[tbl_{2:8,6}];
                            % [results,~,~,gnames] = multcompare(stats,'Dimension',[1 2 3])
                        end
                        p_thresh = p_<0.05;
                        
                        %Fractions of units with different tunings:
                        q_ = [sum(sum(p_thresh,2)==0),...                                    % No tuning
                            sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)==0),...   % Pure tuning only
                            sum(sum(p_thresh(:,1:3),2)==0 & sum(p_thresh(:,4:7),2)>0),...    % Mixed tuning only
                            sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)>0),...    % Mixed and pure tuning
                            sum(sum(p_thresh,2)>0)];                                         % Any tuning
                        
                        q2_= [sum(sum(p_thresh,2)==0),...                                    % No tuning
                            sum(sum(p_thresh(:,1:3),2)==1 & sum(p_thresh(:,4:7),2)==0),...   % CS:  Classical selectivity for one variable
                            sum(sum(p_thresh(:,1:3),2)>1  & sum(p_thresh(:,4:7),2)==0),...   % LMS: Linear mixed selectivity for 2+ variables
                            sum(sum(p_thresh(:,4:7),2)>1),...                                % NMS: Non-linear mixed selectivity
                            sum(sum(p_thresh == [1 0 0 0 0 0 0],2)==7),...     CS for Context
                            sum(sum(p_thresh == [0 1 0 0 0 0 0],2)==7),...     CS for Outcome
                            sum(sum(p_thresh == [0 0 1 0 0 0 0],2)==7),...     CS for Position
                            sum(sum(p_thresh == [1 1 0 0 0 0 0],2)==7),...     LMS for Context X Outcome
                            sum(sum(p_thresh == [1 0 1 0 0 0 0],2)==7),...     LMS for Context X Position
                            sum(sum(p_thresh == [0 1 1 0 0 0 0],2)==7),...     LMS for Outcome X Position
                            sum(sum(p_thresh == [1 1 1 0 0 0 0],2)==7),...     LMS for Context X Outcome X Position
                            sum(sum(p_thresh(:,[4 7])==[1 0],2)==2),...        NMS for Context X Outcome
                            sum(sum(p_thresh(:,[5 7])==[1 0],2)==2),...        NMS for Context X Position
                            sum(sum(p_thresh(:,[6 7])==[1 0],2)==2),...        NMS for Outcome X Position
                            sum(sum(p_thresh(:,7) == 1 ,2)==2),...             NMS for Outcome X Position X Position
                            ];
                        
                        %             q_ = [sum((sum(p_'<0.05)<1)), ...        % No tuning
                        %                 sum((sum(p_(:,1:3)'<0.05)>1)), ... % Pure tuning only
                        %                 sum((sum(p_(:,4:7)'<0.05)>1))];    % Mixed tuning
                        eval(sprintf('D_Ass.%s.ANOVA3.p=p_;',Delays_{iDelay}));
                        eval(sprintf('D_Ass.%s.ANOVA3.F=F_;',Delays_{iDelay}));
                        %eval(sprintf('D_Ass.%s.ANOVA3.Factors=transpose({tbl_{2:8,1}});',Delays_{iDelay}));
                        eval(sprintf('D_Ass.%s.ANOVA3.Factors=Factors;',Delays_{iDelay}));
                        eval(sprintf('D_Ass.%s.ANOVA3.prcSig=sum(p_<0.05)./size(FSC,2)*100;',Delays_{iDelay}));
                        eval(sprintf('D_Ass.%s.ANOVA3.prcTuned=q_./size(FSC,2)*100;',Delays_{iDelay}));
                        eval(sprintf('D_Ass.%s.ANOVA3.prcTuned2=q2_./size(FSC,2)*100;',Delays_{iDelay}));
                        eval(sprintf('D_Ass.%s.ANOVA3.tuning=q_names;',Delays_{iDelay}));
                        eval(sprintf('D_Ass.%s.ANOVA3.tuning2=q2_names;',Delays_{iDelay}));
                    end %/delays
                    %% ANOVA on L/R S/C correct
                    q_names = {'Untuned','CS (Context)','CS (Location)','Any CS','LMS','NMS','Any Tuning'};
                    for iDelay =1:length(Delays_)
                        
                        eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
                        
                        if length(t_.SamplePress_LeftCorrect)>=nDrawnTrials && length(t_.SamplePress_RightCorrect)>=nDrawnTrials
                            Trials = [t_.SamplePress_LeftCorrect;...
                                t_.ChoicePress_LeftCorrect;...
                                t_.SamplePress_RightCorrect;...
                                t_.ChoicePress_RightCorrect];
                            
                            SC = [repmat({'Sample'},length(t_.SamplePress_LeftCorrect),1);...
                                repmat({'Choice'},length(t_.ChoicePress_LeftCorrect),1);...
                                repmat({'Sample'},length(t_.SamplePress_RightCorrect),1);...
                                repmat({'Choice'},length(t_.ChoicePress_RightCorrect),1)];
                            
                            
                            LR = [repmat({'Left'},length(t_.SamplePress_LeftCorrect),1);...
                                repmat({'Left'},length(t_.ChoicePress_LeftCorrect),1);...
                                repmat({'Right'},length(t_.SamplePress_RightCorrect),1);...
                                repmat({'Right'},length(t_.ChoicePress_RightCorrect),1)];
                            
                            nL = length(t_.SamplePress_LeftCorrect);
                            nR = length(t_.SamplePress_RightCorrect);
                            
                            FSC = ExtractFiringRates(Trials,tlimsANOVA,Ass.Tmtx,Ass.FSCsel{iArea},tbANOVA);
                            
                            Stats_ = CalculateANOVA(FSC,LR,SC,nDraws,nDrawnTrials,PlotOnline,resampleTrials);
                            %% Collapse results and write out
                            if resampleTrials
                                p_       = Stats_.p_DrawnMetaP;
                                p_thresh = Stats_.p_DrawnThresh;
                                F_       = Stats_.F_DrawnMetaP;
                            else
                                p_       = Stats_.p_Total;
                                p_thresh = Stats_.p_TotalThresh;
                                F_       = Stats_.F_Total;
                            end
                            %Fractions of units with different tunings:
                            q_ = [sum(sum(p_thresh,2)==0),...              % No tuning
                                sum(sum(p_thresh == [1 0 0],2)==3),...     % CS for Context
                                sum(sum(p_thresh == [0 1 0],2)==3),...     % CS for Location
                                sum(sum(p_thresh == [1 0 0],2)==3) + sum(sum(p_thresh == [0 1 0],2)==3),... % Any CS
                                sum(sum(p_thresh == [1 1 0],2)==3),...     % LMS for Context x Location
                                sum(sum(p_thresh(:,3) == 1,2)),...         % NMS for Context x Location
                                sum(sum(p_thresh,2)>0)];                   % Any tuning
                            
                            eval(sprintf('D_Ass.%s.ANOVA2.p_cor=p_;',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.F_cor=F_;',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.prcSig_cor=sum(p_<0.05)./size(FSC,2)*100;',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.prcTuned_cor=q_./size(FSC,2)*100;',Delays_{iDelay}));
                            
                        else
                            eval(sprintf('D_Ass.%s.ANOVA2.p_cor=nan(1,3);',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.F_cor=nan(1,3);',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.prcSig_cor=nan(1,3);',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.prcTuned_cor=nan(1,length(q_names));',Delays_{iDelay}));
                        end
                        eval(sprintf('D_Ass.%s.ANOVA2.tuning_cor=q_names;',Delays_{iDelay}));
                        eval(sprintf('D_Ass.%s.ANOVA2.Factors_cor={''Context'',''Position'',''Context*Position''};',Delays_{iDelay}));
                        
                    end
                    clear p_ q_ q_names FSC SC LR CE p_thresh
                    %% ANOVA on L/R S/C error
                    q_names = {'Untuned','CS (Context)','CS (Location)','Any CS','LMS','NMS','Any Tuning'};
                    for iDelay =1:length(Delays_)
                        
                        eval(sprintf('t_ = t.%s;',Delays_{iDelay}))
                        
                        if length(t_.SamplePress_LeftError)>=nDrawnTrials && length(t_.SamplePress_RightError)>=nDrawnTrials
                            Trials = [t_.SamplePress_LeftError';...
                                t_.SamplePress_RightError';...
                                t_.ChoicePress_LeftError';...
                                t_.ChoicePress_RightError'];
                            SC = [repmat({'Sample'},length(t_.SamplePress_LeftError),1);...
                                repmat({'Sample'},length(t_.SamplePress_RightError),1);...
                                repmat({'Choice'},length(t_.ChoicePress_LeftError),1);...
                                repmat({'Choice'},length(t_.ChoicePress_RightError),1)];
                            
                            LR = [repmat({'Left'},length(t_.SamplePress_LeftError),1);...
                                repmat({'Right'},length(t_.SamplePress_RightError),1);...
                                repmat({'Left'},length(t_.ChoicePress_LeftError),1);...
                                repmat({'Right'},length(t_.ChoicePress_RightError),1)];
                            
                            FSC = ExtractFiringRates(Trials,tlimsANOVA,Ass.Tmtx,Ass.FSCsel{iArea},tbANOVA);
                            idx = find(isnan(FSC(:,1)));
                            if ~isempty(idx)
                                if idx<=length(LR)/2
                                    idx = [idx,idx+length(LR)/2];
                                else
                                    idx = [idx-length(LR)/2,idx];
                                end
                                SC(idx)=[];
                                LR(idx)=[];
                                FSC(idx,:)=[];
                            end
                           
                            Stats_ = CalculateANOVA(FSC,LR,SC,nDraws,nDrawnTrials,PlotOnline,resampleTrials);
                            %% Collapse results and write out
                            if resampleTrials
                                p_       = Stats_.p_DrawnMetaP;
                                p_thresh = Stats_.p_DrawnThresh;
                                F_       = Stats_.F_DrawnMetaP;
                            else
                                p_       = Stats_.p_Total;
                                p_thresh = Stats_.p_TotalThresh;
                                F_       = Stats_.F_Total;
                            end
                            %Fractions of units with different tunings:
                            q_ = [sum(sum(p_thresh,2)==0),...              % No tuning
                                sum(sum(p_thresh == [1 0 0],2)==3),...     % CS for Context
                                sum(sum(p_thresh == [0 1 0],2)==3),...     % CS for Location
                                sum(sum(p_thresh == [1 0 0],2)==3) + sum(sum(p_thresh == [0 1 0],2)==3),... % Any CS
                                sum(sum(p_thresh == [1 1 0],2)==3),...     % LMS for Context x Location
                                sum(sum(p_thresh(:,3) == 1,2)),...         % NMS for Context x Location
                                sum(sum(p_thresh,2)>0)];                   % Any tuning
                            
                            eval(sprintf('D_Ass.%s.ANOVA2.p_err=p_;',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.F_err=F_;',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.prcSig_err=sum(p_<0.05)./size(FSC,2)*100;',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.prcTuned_err=q_./size(FSC,2)*100;',Delays_{iDelay}));
                        else
                            eval(sprintf('D_Ass.%s.ANOVA2.p_err=nan(1,3);',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.F_err=nan(1,3);',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.prcSig_err=nan(1,3);',Delays_{iDelay}));
                            eval(sprintf('D_Ass.%s.ANOVA2.prcTuned_err=nan(1,length(q_names));',Delays_{iDelay}));
                        end
                        eval(sprintf('D_Ass.%s.ANOVA2.tuning_err=q_names;',Delays_{iDelay}));
                        eval(sprintf('D_Ass.%s.ANOVA2.Factors_err={''Context'',''Position'',''Context*Position''};',Delays_{iDelay}));
                    end %/Delays
                    clear p_ q_ q_names  FSC SC LR CE p_thresh
                    %% Save results
                    if resampleTrials
                        fnOut = sprintf('%sMixedSelectivity%sTrialNosMatched%s%s_%s_MixedSelectivity_Ass.mat',pat,filesep,filesep,fname,Areas{iArea});
                    else
                        fnOut = sprintf('%sMixedSelectivity%s%s_%s_MixedSelectivity_Ass.mat',pat,filesep,fname,Areas{iArea});
                    end
                    save(fnOut,'D_Ass','tlimsANOVA','tbANOVA','-v7.3');
                    
                    fprintf('Done.\n')
                end
            end %/areas
        end %/files
    end
end%/Targets

end

function FR = ExtractFiringRates(Trials,tlimsANOVA,Tmtx,iFR_,tbANOVA)
    FR =[];
    for iTrial =1:length(Trials)
        try
            tlims_ = Trials(iTrial)/1e6+tlimsANOVA;
            tlims_ = closest(Tmtx,tlims_);
            FR(iTrial,1:size(iFR_,2))= mean(iFR_(tlims_(1):tlims_(1)+length(tbANOVA)-1,:));
        catch
            FR(iTrial,1:size(iFR_,2))= nan(1,size(iFR_,2));
        end
    end
end

function Stats_ = CalculateANOVA(FR,LR,SC,nDraws,nDrawnTrials,PlotOnline,resampleTrials)
    %% (0) break out trials matrices to two cells, one for each direction:
    
    % How to report the meta-p-values?
    StatChoice = 3;
    % 1 - Thresholded by distribution of null p-values
    % 2 - Meta-p-value Using Fisher's Method
    % 3 - Meta-p-value Using Stouffer's Method
    % 4 - Approximate p-value approach from Parthasarathy 2017: No. of overlapping points divided by no.draws
    % 5 - Report the lowest returned p-value from all the draws
            
    nL = sum(strcmp(LR,'Left'))./2;
    nR = sum(strcmp(LR,'Right'))./2;

    FR_{1} = FR(1:2*nL,:);
    FR_{2} = FR((2*nL+1):2*(nL+nR),:);

    LR_{1} = LR(1:2*nL,:);
    LR_{2} = LR((2*nL+1):2*(nL+nR),:);

    SC_{1} = SC(1:2*nL,:);
    SC_{2} = SC((2*nL+1):2*(nL+nR),:); 

    SC(isnan(sum(FR,2)))=[];
    LR(isnan(sum(FR,2)))=[];
    FR(isnan(sum(FR,2)),:)=[];
    
    for i=1:2
        idx = isnan(sum(FR_{i},2));
        LR_{i}(idx)=[];
        SC_{i}(idx)=[];
        FR_{i}(idx,:)=[];
    end; clear i idx
    nL = length(LR_{1})/2;
    nR = length(LR_{2})/2;

    LowTrialFlag = 0;
    if min([nL nR])<=nDrawnTrials
        LowTrialFlag = 1;
    end

                    
    %% (1) Loop across units

    for iUnit=1:size(FR,2)
        disp(sprintf('Bootstrapping item no %d/%d...',iUnit,size(FR,2)))
        
            
        % (0) H1_total: Run on all trials:
        [Stats_.p_Total(iUnit,:),tbl_] = anovan(FR(:,iUnit),{SC LR},'model','full','varnames',{'Context','Position'},'display','off');
         Stats_.F_Total(iUnit,:)=[tbl_{2:4,6}];
         
        
        Stats_.p_Drawn{iUnit}     = nan(nDraws,3);
        Stats_.p_DrawnPerm{iUnit} = nan(nDraws,3);
        Stats_.F_Drawn{iUnit}     = nan(nDraws,3);
        Stats_.F_DrawnPerm{iUnit} = nan(nDraws,3);
        
        parfor iDraw = 1:nDraws
            
            % (1) All trials but with permuted labels %%%%%%%%%%%%%%%%%%%%%
            idxL = 1:nL;            idxL = [idxL,idxL+nL]; % N.B. Double up to include both sample
            idxR = 1:nR;            idxR = [idxR,idxR+nR]; % and choice epochs due to stacking
            
            FR__ = [FR_{1}(idxL,iUnit); FR_{2}(idxR,iUnit)];
            LR__ = [LR_{1}(idxL);       LR_{2}(idxR)];
            SC__ = [SC_{1}(idxL);       SC_{2}(idxR)];
            
            [p_perm_Total(iDraw,:),tbl_] = anovan(FR__(randperm(length(FR__))),{SC__(randperm(length(SC__))), LR__(randperm(length(LR__)))},'model','full','varnames',{'Context','Position'},'display','off');
            F_Totaldraw(iDraw,:)=[tbl_{2:4,6}];
            if resampleTrials
            if ~LowTrialFlag
                % (2) Random subsample of trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                idxL = randsample(nL,nDrawnTrials);       idxL = [idxL;idxL+nL];
                idxR = randsample(nR,nDrawnTrials);       idxR = [idxR;idxR+nR];
                
                FR__ = [FR_{1}(idxL,iUnit); FR_{2}(idxR,iUnit)];
                LR__ = [LR_{1}(idxL);       LR_{2}(idxR)];
                SC__ = [SC_{1}(idxL);       SC_{2}(idxR)];
                
                % H1_Subsample
                [p_(iDraw,:),tbl_] = anovan(FR__,{SC__ LR__},'model','full','varnames',{'Context','Position'},'display','off');
                F_(iDraw,:)=[tbl_{2:4,6}];
                
                % H0_subsample
                [p_perm_draw(iDraw,:),tbl_] = anovan(FR__(randperm(length(FR__))),{SC__(randperm(length(SC__))) LR__(randperm(length(LR__)))},'model','full','varnames',{'Context','Position'},'display','off');
                F_draw(iDraw,:)=[tbl_{2:4,6}];
            end
            end
            
        end
        
        % 5th percentile, i.e. lower 0.05 bound of distribution of randomly-drawn p-values... should average out to 0.05
        Stats_.p_TotalPerm(iUnit,:) = prctile(p_perm_Total,5);
        %Stats_.F_TotalPerm{iUnit}   = F_Totaldraw; Don't really need this!
        if ~LowTrialFlag && resampleTrials
            Stats_.p_Drawn{iUnit}     = p_;
            Stats_.p_DrawnPerm{iUnit} = p_perm_draw;
            Stats_.F_Drawn{iUnit}     = F_;
            Stats_.F_DrawnPerm{iUnit} = F_draw;
        end
        clear p_ p_perm_draw F_ F_draw F_Totaldraw p_perm_Total tbl_
    end
    
    
    %% Threshold the drawn distributions
    % (1) ANOVA on all trials...
    Stats_.p_TotalThresh = Stats_.p_Total<Stats_.p_TotalPerm;
    Stats_.F_TotalThresh = Stats_.F_Total; Stats_.F_TotalThresh(~Stats_.p_TotalThresh) = NaN;
    % (2) ANOVA on random draws
    if resampleTrials
        if ~LowTrialFlag
            
            bins = 0:0.01:1;
            % calculate summary stats on random draws
            for iUnit=1:size(FR,2)
                switch StatChoice
                    case 1 % Thresholded by distribution of null p-values
                        Stats_.p_DrawnMetaP(iUnit,:)  = mean(Stats_.p_Drawn{iUnit});
                        Stats_.p_DrawnThresh(iUnit,:) = mean(Stats_.p_Drawn{iUnit})<prctile(Stats_.p_DrawnPerm{iUnit},5);
                        
                    case 2 % Meta-p-value Using Fisher's Method
                        for itest=1:3
                            pvals = Stats_.p_Drawn{iUnit}(:,itest);
                            group_pval(itest) = Fisher(pvals);
                        end; clear itest
                        Stats_.p_DrawnMetaP(iUnit,:)  = group_pval;
                        Stats_.p_DrawnThresh(iUnit,:) = group_pval<0.05;
                        
                    case 3 % Meta-p-value Using Stouffer's Method
                        for itest=1:3 
                            pvals = Stats_.p_Drawn{iUnit}(:,itest);
                            group_pval(itest) = stouffer(pvals);
                        end; clear itest
                        Stats_.p_DrawnMetaP(iUnit,:)  = group_pval;
                        Stats_.p_DrawnThresh(iUnit,:) = group_pval<0.05;
                        
                    case 4 % Approximate p-value approach from Parthasarathy 2017 No. of overlapping points divided by no.draws.
                       % N.B. Need to investigate effect of no. bins
                        for itest =1:3
                            A = Stats_.p_Drawn{iUnit}(:,itest);
                            B = Stats_.p_DrawnPerm{iUnit}(:,itest);
                            % A(A>prctile(A,95))=[];B(B<prctile(B,5))=[];
                            x = histc(A,bins) - histc(B,bins);
                            x(x<=0 | isnan(x)) = [];
                            Stats_.p_DrawnMetaP(iUnit,itest)  = (1+nDraws-sum(x))./(nDraws+1);
                            Stats_.p_DrawnThresh(iUnit,itest) = Stats_.p_DrawnMetaP(iUnit,itest)<0.05;
                        end; clear itest A B x
                        
                    case 5 % Report the lowest returned p-value from all the draws
                        Stats_.p_DrawnMetaP(iUnit,:) = min(Stats_.p_Drawn{iUnit});
                        Stats_.p_DrawnThresh(iUnit,:)= min(Stats_.p_Drawn{iUnit})<= min(Stats_.p_DrawnPerm{iUnit}); % prctile(Stats_.p_DrawnPerm{iUnit},5); 
                end
                Stats_.F_DrawnMetaP(iUnit,:) = mean(Stats_.F_Drawn{iUnit});

            end
            if PlotOnline
                for iUnit=1:size(FR,2)
                    figure('name',sprintf('Unit %d',iUnit));
                    for itest = 1:3
                        pThresh = prctile(Stats_.p_Drawn{iUnit}(:,itest),[5 95]);
                        p_permThresh = prctile(Stats_.p_DrawnPerm{iUnit}(:,itest),[5 95]);

                        pMeta(itest)  = Stats_.p_DrawnMetaP(iUnit,itest); % This is Bonferoni correction
                        pThresh(itest) = Stats_.p_DrawnThresh(iUnit,itest); % This is Bonferoni correction
                        
                        subplot(1,3,itest); hold on
                        
                        area(bins,histc(Stats_.p_DrawnPerm{iUnit}(:,itest),bins)./nDraws,'EdgeColor','r','FaceColor','r','FaceAlpha',0.5)
                        area(bins,histc(Stats_.p_Drawn{iUnit}(:,itest),bins)./nDraws,'EdgeColor','k','FaceColor','k','FaceAlpha',0.5)
                        
                        scatter(Stats_.p_Total(iUnit,itest),0.1,'*k')
                        scatter(Stats_.p_TotalPerm(iUnit,itest),0.1,'*r')
                        title(sprintf('Sig?=%d(p~%0.4f)',pThresh(itest),pMeta(itest)))
                        axis([0 1 0 0.3]); %set(gca,'XScale','log')
                        clear pThresh p_permThresh
                    end
                    drawnow;%pause(0.01)
                end
            end

        else
            Stats_.p_DrawnThresh = nan(size(FR_,2),3);
            Stats_.p_DrawnMetaP  = nan(size(FR_,2),3);
            Stats_.F_DrawnMetaP  = nan(size(FR_,2),3);
        end
    end
    
end

function pcomb = stouffer(p)
    % Stouffer et al's (1949) unweighted method for combination of independent p-values via z's
    if length(p)==0
        error('pfast was passed an empty array of p-values')
        pcomb=1;
    else
        pcomb = (1-erf(sum(sqrt(2) * erfinv(1-2*p))/sqrt(2*length(p))))/2;
    end
end

function pcomb = Fisher(p)
    % Fisher's method for combining p-values
    if length(p)==0
        error('pfast was passed an empty array of p-values')
        pcomb=1;
    else
        pcomb = 1 - chi2cdf(sum(-2.*log(p)),2*length(p));
    end
end
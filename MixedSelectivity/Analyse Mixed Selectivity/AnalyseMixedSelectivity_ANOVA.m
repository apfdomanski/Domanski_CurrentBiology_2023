clear 
Targets = {'SHORT','MEDIUM','LONG'};
Delays_ = {'Short','Medium','Long'};
tlimsShort  = [-4 10];
tlimsMedium = [-8 10];
tlimsLong   = [-16 10];
tlimsANOVA  = [-5 5];
shift = 0;
plotOnline = false;
bw = 0.05;
tbShort  = tlimsShort(1) : bw : tlimsShort(2);
tbMedium = tlimsMedium(1): bw : tlimsMedium(2);
tbLong   = tlimsLong(1)  : bw : tlimsLong(2);
tbANOVA  = tlimsANOVA(1) : bw : tlimsANOVA(2);
clear Av
warning ('off')
pat = 'C:\Analysis\AssemblyAnalysis\raw';
cd(pat)
fileList=[];
for i=1:length(Targets)
    fileList=[fileList;eval(sprintf('dir(''allTimestamps\\*%s*.mat'');',Targets{i}))];
end
Areas = {'HP','PFC'};
%% Batch process
for iFile =1:length(fileList)
    
    fname = strtok( fileList(iFile).name , '_' );
    
    load( sprintf( '%s\\allTimestamps\\%s_Events.mat' , pat , fname ) );
    load( sprintf('%s\\%s.mat' , pat , fname ) );
    for iArea = 1:length( Areas )
        fprintf( 'Analysing run %d/%d %s (%s)...' , iFile , length(fileList) , fname , Areas{iArea} )
        load( sprintf( '%s\\KDE_bins\\%s_%s_iFR50.mat' , pat , fname , Areas{iArea} ) );
        
        iFR_ = iFR;
        % iFR_=zscore(iFR_);
        
        %% sanity check: plot original and imported timestamps
        % figure; hold on
        % t_ = [trangeleft_sample(:,1)]+5000000;
        % plot(t_ ,ones(length(t_)),'.b')
        % t_ = [t.Medium.SamplePress_LeftCorrect;t.Long.SamplePress_LeftCorrect];
        % plot(t_,ones(length(t_),1),'or')     
        %% ANOVA on L/R S/C C/E 
        
        q_names = {'Untuned','Pure only','Mixed only','Mixed and Pure','Any tuning'}; 
        
        for iDelay = 1 : length( Delays_ )
            
            eval( sprintf( 't_ = t.%s;' , Delays_{iDelay} ) )
            
            % Trial times
            
            Trials   = [t_.SamplePress_LeftCorrect;...
                        t_.SamplePress_LeftError;...
                        t_.SamplePress_RightCorrect;...
                        t_.SamplePress_RightError;...
                        t_.ChoicePress_LeftCorrect;...
                        t_.ChoicePress_LeftError;...
                        t_.ChoicePress_RightCorrect;...
                        t_.ChoicePress_RightError];
            
            % Mark up each trial type with at string for the ANOVA
            
            SC       = [repmat( {'Sample'} , length(t_.SamplePress_LeftCorrect)  , 1);...
                        repmat( {'Sample'} , length(t_.SamplePress_LeftError)    , 1);...
                        repmat( {'Sample'} , length(t_.SamplePress_RightCorrect) , 1);...
                        repmat( {'Sample'} , length(t_.SamplePress_RightError)   , 1);...
                        repmat( {'Choice'} , length(t_.ChoicePress_LeftCorrect)  , 1);...
                        repmat( {'Choice'} , length(t_.ChoicePress_LeftError)    , 1);...
                        repmat( {'Choice'} , length(t_.ChoicePress_RightCorrect) , 1);...
                        repmat( {'Choice'} , length(t_.ChoicePress_RightError)   , 1)];
            
            CE       = [repmat( {'Correct'} , length(t_.SamplePress_LeftCorrect)  , 1);...
                        repmat( {'Error'}   , length(t_.SamplePress_LeftError)    , 1);...
                        repmat( {'Correct'} , length(t_.SamplePress_RightCorrect) , 1);...
                        repmat( {'Error'}   , length(t_.SamplePress_RightError)   , 1);...
                        repmat( {'Correct'} , length(t_.ChoicePress_LeftCorrect)  , 1);...
                        repmat( {'Error'}   , length(t_.ChoicePress_LeftError)    , 1);...
                        repmat( {'Correct'} , length(t_.ChoicePress_RightCorrect) , 1);...
                        repmat( {'Error'}   , length(t_.ChoicePress_RightError)   , 1)];
            
            LR       = [repmat( {'Left'}    , length(t_.SamplePress_LeftCorrect)  , 1);...
                        repmat( {'Left'}    , length(t_.SamplePress_LeftError)    , 1);...
                        repmat( {'Right'}   , length(t_.SamplePress_RightCorrect) , 1);...
                        repmat( {'Right'}   , length(t_.SamplePress_RightError)   , 1);...
                        repmat( {'Left'}    , length(t_.ChoicePress_LeftCorrect)  , 1);...
                        repmat( {'Left'}    , length(t_.ChoicePress_LeftError)    , 1);...
                        repmat( {'Right'}   , length(t_.ChoicePress_RightCorrect) , 1);...
                        repmat( {'Right'}   , length(t_.ChoicePress_RightError)   , 1)];
            
            FR =[];
            for iTrial =1 : length(Trials)
                
                try
                    
                    tlims_ = Trials(iTrial) / 1e6 + tlimsANOVA + shift;
                    tlims_ = closest(Tmtx,tlims_);
                    
                    FR( iTrial  , 1:size(iFR_,2) ) = mean( iFR_(tlims_(1) : tlims_(1) + length(tbANOVA)-1 , : ) ...
                                                     );
                catch
                    
                    FR( iTrial  , 1:size(iFR_,2) ) =  nan(1,size(iFR_,2));
                    
                end
                
            end
           
            SC( isnan(sum(FR,2) ))     = [];
            CE( isnan(sum(FR,2) ))     = [];
            LR( isnan(sum(FR,2) ))     = [];
            FR( isnan(sum(FR,2) ),:)   = [];
            
            p_ = zeros( size(iFR_,2) , 7);
            F_ = p_;
            
            % Run multi-way ANOVA for each single unit
            for iUnit = 1:size(iFR_,2)
                
                [ p_(iUnit,:) , tbl_ ] = anovan( FR(:,iUnit) , {SC CE LR},...
                                                'model','full', ...
                                                'varnames',{'Context','Outcome','Position'}, ...
                                                'display','off');
                                            
                F_(iUnit,:) = [tbl_{2:8,6}];
                % [results,~,~,gnames] = multcompare(stats,'Dimension',[1 2 3])
                
            end
            
            % Gather significnant simple/interaction terms
            p_thresh = p_ < 0.05;
            
            q_ = [sum( sum(p_thresh,2) == 0)                                    ,...  % No tuning
                  sum( sum(p_thresh(:,1:3),2)>0  & sum( p_thresh(:,4:7),2)==0 ) ,...  % Pure tuning only
                  sum( sum(p_thresh(:,1:3),2)==0 & sum( p_thresh(:,4:7),2)>0 )  ,...  % Mixed tuning only
                  sum( sum(p_thresh(:,1:3),2)>0  & sum( p_thresh(:,4:7),2)>0 )  ,...  % Mixed and pure tuning
                  sum( sum(p_thresh,2) > 0) ...                                       % Any tuning
                  ];                                         
            % q_ = [sum( (sum(p_'<0.05)<1))        , ...    % No tuning
            %       sum( (sum(p_(:,1:3)'<0.05)>1)) , ...    % Pure tuning only
            %       sum( (sum(p_(:,4:7)'<0.05)>1)) ...      % Mixed tuning
            %       ]
            
            % Write results out to data structure 
            
            eval( sprintf('D.%s.ANOVA.p = p_;' ,                                  Delays_{iDelay}));
            eval( sprintf('D.%s.ANOVA.F = F_;' ,                                  Delays_{iDelay}));
            eval( sprintf('D.%s.ANOVA.Factors = transpose({tbl_{2:8,1}});' ,      Delays_{iDelay}));
            eval( sprintf('D.%s.ANOVA.prcSig = sum(p_<0.05)./size(iFR_,2)*100;' , Delays_{iDelay}));
            eval( sprintf('D.%s.ANOVA.prcTuned = q_./size(iFR_,2)*100;' ,         Delays_{iDelay}));
            eval( sprintf('D.%s.ANOVA.tuning = q_names;' ,                        Delays_{iDelay}));
            
            
        end
        
        clear p_ q_ q_names FR SC LR CE p_thresh
        %% Save results
        
        fnOut = sprintf( '%s\\MixedSelectivity\\%s_%s_MixedSelectivity.mat' , pat , fname , Areas{iArea});
        save( fnOut , 'D' , 'tbShort' , 'tbMedium' , 'tbLong');
        fprintf('Done.\n')
        
        iFR_ = iFR;
        % iFR_=zscore(iFR_);
        
    end
end
clearvars -except Delays_ tbShort tbMedium tbLong tbANOVA  pat fileList Areas 
%% batch import for plotting
D_.Short.LR.Ft2         = cell(length(Areas),1);
D_.Medium.LR.Ft2        = cell(length(Areas),1);
D_.Long.LR.Ft2          = cell(length(Areas),1);
D_.Short.CE.Ft2         = cell(length(Areas),1);
D_.Medium.CE.Ft2        = cell(length(Areas),1);
D_.Long.CE.Ft2          = cell(length(Areas),1);
D_.Short.SC.Ft2         = cell(length(Areas),1);
D_.Medium.SC.Ft2        = cell(length(Areas),1);
D_.Long.SC.Ft2          = cell(length(Areas),1);
D_.Short.ANOVA.prcSig   = cell(length(Areas),1);
D_.Short.ANOVA.tuning   = cell(length(Areas),1);

for iArea = 1 : length(Areas)
    
    for iFile = 1 : length(fileList)
        
        fname = strtok(fileList(iFile).name , '_' );
        fnIn  = sprintf('%s\\MixedSelectivity\\%s_%s_MixedSelectivity.mat' , pat , fname , Areas{iArea} );
        load(fnIn ,'D');
        
        D_.Short.LR.Ft2{iArea}(:,iFile)  = D.Short.LR.Ft2;
        D_.Medium.LR.Ft2{iArea}(:,iFile) = D.Medium.LR.Ft2;
        D_.Long.LR.Ft2{iArea}(:,iFile)   = D.Long.LR.Ft2;
        
        D_.Short.CE.Ft2{iArea}(:,iFile)  = D.Short.CE.Ft2;
        D_.Medium.CE.Ft2{iArea}(:,iFile) = D.Medium.CE.Ft2;
        D_.Long.CE.Ft2{iArea}(:,iFile)   = D.Long.CE.Ft2;
        
        D_.Short.SC.Ft2{iArea}(:,iFile)  = D.Short.SC.Ft2;
        D_.Medium.SC.Ft2{iArea}(:,iFile) = D.Medium.SC.Ft2;
        D_.Long.SC.Ft2{iArea}(:,iFile)   = D.Long.SC.Ft2;
        
        D_.Short.ANOVA.prcSig{iArea}(:,iFile)  = D.Short.ANOVA.prcSig;
        D_.Medium.ANOVA.prcSig{iArea}(:,iFile) = D.Medium.ANOVA.prcSig;
        D_.Long.ANOVA.prcSig{iArea}(:,iFile)   = D.Long.ANOVA.prcSig;
        
        D_.Short.ANOVA.prcTuned{iArea}(:,iFile)  = D.Short.ANOVA.prcTuned;
        D_.Medium.ANOVA.prcTuned{iArea}(:,iFile) = D.Medium.ANOVA.prcTuned;
        D_.Long.ANOVA.prcTuned{iArea}(:,iFile)   = D.Long.ANOVA.prcTuned;

    end
    
end

for iArea = 1:length(Areas)
    D_.Short.ANOVA.prcSigMean{iArea}    = nanmean(D_.Short.ANOVA.prcSig{iArea},2);
    D_.Medium.ANOVA.prcSigMean{iArea}   = nanmean(D_.Medium.ANOVA.prcSig{iArea},2);
    D_.Long.ANOVA.prcSigMean{iArea}     = nanmean(D_.Long.ANOVA.prcSig{iArea},2);
    
    D_.Short.ANOVA.prcSigSEM{iArea}     = nansem(D_.Short.ANOVA.prcSig{iArea},2);
    D_.Medium.ANOVA.prcSigSEM{iArea}    = nansem(D_.Medium.ANOVA.prcSig{iArea},2);
    D_.Long.ANOVA.prcSigSEM{iArea}      = nansem(D_.Long.ANOVA.prcSig{iArea},2);
    
    D_.Short.ANOVA.prcTunedMean{iArea}  = nanmean(D_.Short.ANOVA.prcTuned{iArea},2);
    D_.Medium.ANOVA.prcTunedMean{iArea} = nanmean(D_.Medium.ANOVA.prcTuned{iArea},2);
    D_.Long.ANOVA.prcTunedMean{iArea}   = nanmean(D_.Long.ANOVA.prcTuned{iArea},2);
    
    D_.Short.ANOVA.prcTunedSEM{iArea}   = nansem(D_.Short.ANOVA.prcTuned{iArea},2);
    D_.Medium.ANOVA.prcTunedSEM{iArea}  = nansem(D_.Medium.ANOVA.prcTuned{iArea},2);
    D_.Long.ANOVA.prcTunedSEM{iArea}    = nansem(D_.Long.ANOVA.prcTuned{iArea},2);
end

tuning  = D.Short.ANOVA.tuning;
Factors = D.Short.ANOVA.Factors;
clear D

%% Plot proportions of tuned/untuned cells - group by delay
for iArea = 1:length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    
    subplot(1,3,1);hold on
    bar(D_.Short.ANOVA.prcTunedMean{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    ylabel('% of units')
    title('Short delay trials')
    
    subplot(1,3,2);hold on
    bar(D_.Medium.ANOVA.prcTunedMean{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    title('Medium delay trials')
    
    subplot(1,3,3);hold on
    bar(D_.Long.ANOVA.prcTunedMean{iArea})
    set(gca,'XTick',1:length(tuning),'XTicklabel',tuning,'XTickLabelRotation',50)
    title('Long delay trials')
end
%% Plot proportions of tuned/untuned cells - by delay

for iArea = 1:length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]);
    for iType = 1:length(tuning)
        subplot(1,length(tuning),iType);hold on
        m   = [D_.Short.ANOVA.prcTunedMean{iArea}(iType)
            D_.Medium.ANOVA.prcTunedMean{iArea}(iType)
            D_.Long.ANOVA.prcTunedMean{iArea}(iType)];
        
        e = [D_.Short.ANOVA.prcTunedSEM{iArea}(iType)
            D_.Medium.ANOVA.prcTunedSEM{iArea}(iType)
            D_.Long.ANOVA.prcTunedSEM{iArea}(iType)];
        
        
        b = bar(m);
        b.LineStyle='none';b.EdgeColor=[0 0 1];
        errorbar(m,e,'.b','Linewidth',1.5)
        axis([0 4 0 100]);
        set(gca,'XTick',1:length(Delays_),'XTicklabel',Delays_,'XTickLabelRotation',50)
        title(tuning{iType})
        if iType ==1
            ylabel('% of units')
        end

    end
    
end
%% Plot proportions of tuned/untuned cells - by tuning type
for iArea = 1:length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    for iDelay=1:length(Delays_)
        subplot(1,length(Delays_),iDelay);hold on
        eval(sprintf('b=bar(D_.%s.ANOVA.prcSigMean{iArea});',Delays_{iDelay}))
        b.LineStyle='none';b.EdgeColor=[0 0 1];
        eval(sprintf('e=errorbar(D_.%s.ANOVA.prcSigMean{iArea},D_.%s.ANOVA.prcSigSEM{iArea},''color'',''b'',''LineStyle'',''none'')',Delays_{iDelay},Delays_{iDelay}))
        e.LineWidth=1.5;
        set(gca,'XTick',1:length(Factors),'XTicklabel',Factors,'XTickLabelRotation',50)
        title([Delays_{iDelay} ' Delay'])

        if iDelay==1
            ylabel('% of units')
        end
    
    end
end
%% Plot proportions of tuned/untuned cells - by tuning type AND Delays
for iArea = 2%1:length(Areas)
    figure('name',['Percentage of tuned cells: ' Areas{iArea}]); 
    for iType = 1:length(Factors)
        subplot(1,length(Factors),iType);hold on
         m   = [D_.Short.ANOVA.prcSigMean{iArea}(iType)
            D_.Medium.ANOVA.prcSigMean{iArea}(iType)
            D_.Long.ANOVA.prcSigMean{iArea}(iType)];
          e = [D_.Short.ANOVA.prcSigSEM{iArea}(iType)
            D_.Medium.ANOVA.prcSigSEM{iArea}(iType)
            D_.Long.ANOVA.prcSigSEM{iArea}(iType)];
        
        b = bar(m);
        b.LineStyle='none';b.EdgeColor=[0 0 1];
        errorbar(m,e,'.b','Linewidth',1.5)
        axis([0 4 0 50]);
        set(gca,'XTick',1:length(Delays_),'XTicklabel',Delays_,'XTickLabelRotation',50)
        title(Factors{iType})
        if iType ==1
            ylabel('% of units')
        end
    end
end



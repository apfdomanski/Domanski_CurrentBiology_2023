clear 
Targets = {'SHORT','MEDIUM','LONG'};
Delays_ = {'Short','Medium','Long'};
tlimsShort=[-10 5];
tlimsMedium=[-10 5];
tlimsLong=[-10 5];
tlimsANOVA = [-5 5];
shift = 0;
plotOnline = false;
bw=0.05;
tbShort=tlimsShort(1):bw:tlimsShort(2);
tbMedium=tlimsMedium(1):bw:tlimsMedium(2);
tbLong=tlimsLong(1):bw:tlimsLong(2);
tbANOVA=tlimsANOVA(1):bw:tlimsANOVA(2);
clear Av
warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
else
    pat = '/Volumes/HDD2/DNMTP/raw';
end
cd(pat)
fileList=dir(['allTimestamps',filesep,'*LONG*.mat']);
Areas = {'HP','PFC'};



%%%%%
iFile =6%:length(fileList)
iArea =2% 1:length(Areas)
%%%%%
%% Load an example
fname=strtok(fileList(iFile).name,'_');

load(sprintf('%s%sallTimestamps%s%s_Events.mat',pat,filesep,filesep,fname));
load(sprintf('%s%s%s.mat',pat,filesep,fname));

fprintf('Analysing run %d/%d %s (%s)...',iFile,length(fileList),fname,Areas{iArea})
% load(sprintf('%s%sKDE_bins%s%s_%s_iFR50.mat',pat,filesep,filesep,fname,Areas{iArea}));
load(sprintf('%s%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',pat,filesep,filesep,fname,Areas{iArea}));

iFR_=iFR;
% iFR_=zscore(iFR_);
%% ANOVA on L/R S/C C/E
q_names = {'Untuned','Pure only','Mixed only','Mixed and Pure','Any tuning'};
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
    q_ = [sum(sum(p_thresh,2)==0),...                                      % No tuning
        sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)==0),...   % Pure tuning only
        sum(sum(p_thresh(:,1:3),2)==0 & sum(p_thresh(:,4:7),2)>0),...    % Mixed tuning only
        sum(sum(p_thresh(:,1:3),2)>0  & sum(p_thresh(:,4:7),2)>0),...    % Mixed and pure tuning
        sum(sum(p_thresh,2)>0)];                                          % Any tuning
    %             q_ = [sum((sum(p_'<0.05)<1)), ...        % No tuning
    %                 sum((sum(p_(:,1:3)'<0.05)>1)), ... % Pure tuning only
    %                 sum((sum(p_(:,4:7)'<0.05)>1))];    % Mixed tuning
    eval(sprintf('D.%s.ANOVA.p=p_;',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.F=F_;',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.Factors=transpose({tbl_{2:8,1}});',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.prcSig=sum(p_<0.05)./size(iFR_,2)*100;',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.prcTuned=q_./size(iFR_,2)*100;',Delays_{iDelay}));
    eval(sprintf('D.%s.ANOVA.tuning=q_names;',Delays_{iDelay}));
    
    
end

% clear p_ q_ q_names FR SC LR CE p_thresh
%% Make raster + SDF ( Medium and long combined)

D.Long.ANOVA.p<0.05;
% iUnit  = 27 %43%39% 33; %20
iUnit  = 7% 27 %43%39% 33; %20
tlimsLong = [-5 10];

offset=20;
SF=5;



tlims_ = [trangeleft_sample(:,1) ; trangeright_sample(:,1)]*1e-6 +5 + tlimsLong;

figure('color','w'); hold on

plot([5 5],[-offset length(tlims_)],'color',[0 1 0 0.6],'LineWidth',2)
plot([25 25],[-offset length(tlims_)],'color',[1 0 0 0.6],'LineWidth',2)
plot([15 15],[-offset length(tlims_)],'k')

for iTrial =1:length(tlims_)
    times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) ; 
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
end

tlimsLong = [-10 5];

tlims_ = [trangeleft_choice(:,1) ; trangeright_choice(:,1)]*1e-6 +5 + tlimsLong;

for iTrial =1:length(tlims_)
    times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) +15; 
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
end

axis off
plot([31 31],[1 size(trangeleft_choice,1)],'color',[0.5 0 0.5],'LineWidth',1.5,'LineWidth',4)
plot([31 31],[1 size(trangeright_choice,1)]+size(trangeleft_choice,1),'color',[0 0 1],'LineWidth',4)

text(31.5,size(trangeleft_choice,1)/2,'L','color',[0.5 0 0.5])
text(31.5,size(trangeright_choice,1)+size(trangeleft_choice,1)-size(trangeleft_choice,1)/2,'R','color',[0 0 1])


data_{1}.iFR = iFR_;
[~,cutout]=SelTimesFRs([-5 10],Tmtx,data_,trangeleft_sample(:,1)*1e-6 +5 );
LS = nan(300,1);
for iTrial = 1:size(trangeleft_sample,1)
    try
    LS= [LS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs([-5 10],Tmtx,data_,trangeleft_choice(:,1)*1e-6 +5 );
LC = nan(300,1);
for iTrial = 1:size(trangeleft_choice,1)
    try
    LC= [LC,cutout{1}{iTrial}(:,iUnit)];
    end
end


[~,cutout]=SelTimesFRs([-5 10],Tmtx,data_,trangeright_sample(:,1)*1e-6 +5 );
RS = nan(300,1);
for iTrial = 1:size(trangeright_sample,1)
    try
    RS= [RS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs([-5 10],Tmtx,data_,trangeright_choice(:,1)*1e-6 +5 );
RC = nan(300,1);
for iTrial = 1:size(trangeright_choice,1)
    try
    RC= [RC,cutout{1}{iTrial}(:,iUnit)];
    end
end


Lm = [nanmean(LS'),nanmean(LC')]*SF;Le =[nansem(LS'),nansem(LC')]*SF;
Rm = [nanmean(RS'),nanmean(RC')]*SF;Re =[nansem(RS'),nansem(RC')]*SF;
tb_  = (1:600)*0.05;
% figure; hold on
% ciplot(Lm+Le-mean([Lm,Rm])-offset,Lm-Le-mean([Lm,Rm])-offset,tb_,[0.5 0 0.5],0.8);
% ciplot(Rm+Re-mean([Lm,Rm])-offset,Rm-Re-mean([Lm,Rm])-offset,tb_,[0 0 1],0.8);
ciplot(Lm+Le-offset,Lm-Le-offset,tb_,[0.5 0 0.5],0.8);
ciplot(Rm+Re-offset,Rm-Re-offset,tb_,[0 0 1],0.8);
D.Long.ANOVA.p<0.05;
%% 
% file 6 / iUnit 11
iUnit  = 9 %43%39% 33; %20
% iUnit  = 1% 27 %43%39% 33; %20
offset=20;
SF=2;

tlimsSample = [-5 5];
tlimsChoice = [-5 5];

%% Make raster + SDF  - Short ordered by Choice press latency

D.Short.ANOVA.p<0.05;
d_L = t.Short.ChoicePress_LeftCorrect-t.Short.DelayEnd_LeftCorrect;
d_R = t.Short.ChoicePress_RightCorrect-t.Short.DelayEnd_RightCorrect;
% d_L = t.Short.ChoicePress_LeftCorrect-t.Short.NosePoke_LeftCorrect;
% d_R = t.Short.ChoicePress_RightCorrect-t.Short.NosePoke_RightCorrect;
[~,idxL]=sort(d_L);[~,idxR]=sort(d_R);
NP_offset = [d_L(idxL);d_R(idxR)]*1e-6;
NP_offset(NP_offset>abs(tlimsChoice(1)))=NaN;

tlims_ = [t.Short.SamplePress_LeftCorrect(idxL) ; t.Short.SamplePress_RightCorrect(idxR) ]*1e-6 + tlimsSample;

figure('color','w'); hold on

plot(sum(abs(tlimsSample))/2*[1,1],[-offset length(tlims_)+1],'color',[0 1 0 0.6],'LineWidth',2)
plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[-offset length(tlims_)+1],'color',[1 0 0 0.6],'LineWidth',2)
plot(sum(abs(tlimsSample))*[1,1],[-offset length(tlims_)+1],'k')

for iTrial =1:length(tlims_)
    times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) ; 
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
end

tlims_ = [t.Short.ChoicePress_LeftCorrect(idxL) ; t.Short.ChoicePress_RightCorrect(idxR) ]*1e-6+ tlimsChoice;

for iTrial =1:length(tlims_)
    times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) + sum(abs(tlimsSample)); 
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot(-[NP_offset(iTrial)';NP_offset(iTrial)']+15,[iTrial iTrial]+[-0.5; 0.5],'color',[0.9 0.4 0.2],'LineWidth',1.5)

end

axis off
plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(t.Short.ChoicePress_LeftCorrect,1)],'color',[0.5 0 0.5],'LineWidth',1.5,'LineWidth',4)
plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(t.Short.ChoicePress_RightCorrect,1)]+size(t.Short.ChoicePress_LeftCorrect,1),'color',[0 0 1],'LineWidth',4)

text(1.5+sum(abs([tlimsSample tlimsChoice])),size(t.Short.ChoicePress_LeftCorrect,1)/2,'L','color',[0.5 0 0.5])
text(1.5+sum(abs([tlimsSample tlimsChoice])),size(t.Short.ChoicePress_LeftCorrect,1)+size(t.Short.ChoicePress_RightCorrect,1)/2,'R','color',[0 0 1])

data_{1}.iFR = iFR_;
[~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,t.Short.SamplePress_LeftCorrect*1e-6 );
LS = nan(sum(abs([tlimsSample]))/0.05,1);
for iTrial = 1:size(t.Short.SamplePress_LeftCorrect,1)
    try
    LS= [LS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,t.Short.ChoicePress_LeftCorrect*1e-6 );
LC = nan(sum(abs([tlimsChoice]))/0.05,1);
for iTrial = 1:size(t.Short.ChoicePress_LeftCorrect,1)
    try
    LC= [LC,cutout{1}{iTrial}(:,iUnit)];
    end
end


[~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,t.Short.SamplePress_RightCorrect*1e-6 );
RS  = nan(sum(abs([tlimsSample]))/0.05,1);
for iTrial = 1:size(t.Short.SamplePress_RightCorrect,1)
    try
    RS= [RS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,t.Short.ChoicePress_RightCorrect*1e-6);
RC = nan(sum(abs([tlimsChoice]))/0.05,1);
for iTrial = 1:size(t.Short.ChoicePress_RightCorrect,1)
    try
    RC= [RC,cutout{1}{iTrial}(:,iUnit)];
    end
end


Lm = [nanmean(LS'),nanmean(LC')]*SF;Le =[nansem(LS'),nansem(LC')]*SF;
Rm = [nanmean(RS'),nanmean(RC')]*SF;Re =[nansem(RS'),nansem(RC')]*SF;
tb_  = (1:sum(abs([tlimsSample tlimsChoice]))/0.05)*0.05;
% figure; hold on
% ciplot(Lm+Le-mean([Lm,Rm])-offset,Lm-Le-mean([Lm,Rm])-offset,tb_,[0.5 0 0.5],0.8);
% ciplot(Rm+Re-mean([Lm,Rm])-offset,Rm-Re-mean([Lm,Rm])-offset,tb_,[0 0 1],0.8);
ciplot(Lm+Le-offset,Lm-Le-offset,tb_,[0.5 0 0.5],0.8);
ciplot(Rm+Re-offset,Rm-Re-offset,tb_,[0 0 1],0.8);
%% Make raster + SDF  - Medium ordered by Choice press latency

D.Medium.ANOVA.p<0.05;
% d_L = t.Medium.ChoicePress_LeftCorrect-t.Medium.DelayEnd_LeftCorrect;
% d_R = t.Medium.ChoicePress_RightCorrect-t.Medium.DelayEnd_RightCorrect;
d_L = t.Medium.ChoicePress_LeftCorrect-t.Medium.NosePoke_LeftCorrect;
d_R = t.Medium.ChoicePress_RightCorrect-t.Medium.NosePoke_RightCorrect;
[~,idxL]=sort(d_L);[~,idxR]=sort(d_R);
NP_offset = [d_L(idxL);d_R(idxR)]*1e-6;
NP_offset(NP_offset>abs(tlimsChoice(1)))=NaN;
tlims_ = [t.Medium.SamplePress_LeftCorrect(idxL) ; t.Medium.SamplePress_RightCorrect(idxR)]*1e-6 + tlimsSample;
figure('color','w'); hold on

plot(sum(abs(tlimsSample))/2*[1,1],[-offset length(tlims_)+1],'color',[0 1 0 0.6],'LineWidth',2)
plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[-offset length(tlims_)+1],'color',[1 0 0 0.6],'LineWidth',2)
plot(sum(abs(tlimsSample))*[1,1],[-offset length(tlims_)+1],'k')

for iTrial =1:length(tlims_)
    times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) ; 
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
end


tlims_ = [t.Medium.ChoicePress_LeftCorrect(idxL) ; t.Medium.ChoicePress_RightCorrect(idxR) ]*1e-6+ tlimsChoice;

for iTrial =1:length(tlims_)
    times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) + sum(abs(tlimsSample)); 
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot(-[NP_offset(iTrial)';NP_offset(iTrial)']+15,[iTrial iTrial]+[-0.5; 0.5],'color',[0.9 0.4 0.2],'LineWidth',1.5)

end

axis off
plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(t.Medium.ChoicePress_LeftCorrect,1)],'color',[0.5 0 0.5],'LineWidth',1.5,'LineWidth',4)
plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(t.Medium.ChoicePress_RightCorrect,1)]+size(t.Medium.ChoicePress_LeftCorrect,1),'color',[0 0 1],'LineWidth',4)

text(1.5+sum(abs([tlimsSample tlimsChoice])),size(t.Medium.ChoicePress_LeftCorrect,1)/2,'L','color',[0.5 0 0.5])
text(1.5+sum(abs([tlimsSample tlimsChoice])),size(t.Medium.ChoicePress_LeftCorrect,1)+size(t.Medium.ChoicePress_RightCorrect,1)/2,'R','color',[0 0 1])


data_{1}.iFR = iFR_;
[~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,t.Medium.SamplePress_LeftCorrect*1e-6 );
LS = nan(sum(abs([tlimsSample]))/0.05,1);
for iTrial = 1:size(t.Medium.SamplePress_LeftCorrect,1)
    try
    LS= [LS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,t.Medium.ChoicePress_LeftCorrect*1e-6 );
LC = nan(sum(abs([tlimsChoice]))/0.05,1);
for iTrial = 1:size(t.Medium.ChoicePress_LeftCorrect,1)
    try
    LC= [LC,cutout{1}{iTrial}(:,iUnit)];
    end
end


[~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,t.Medium.SamplePress_RightCorrect*1e-6 );
RS  = nan(sum(abs([tlimsSample]))/0.05,1);
for iTrial = 1:size(t.Medium.SamplePress_RightCorrect,1)
    try
    RS= [RS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,t.Medium.ChoicePress_RightCorrect*1e-6);
RC= nan(sum(abs([tlimsChoice]))/0.05,1);
for iTrial = 1:size(t.Medium.ChoicePress_RightCorrect,1)
    try
    RC= [RC,cutout{1}{iTrial}(:,iUnit)];
    end
end


Lm = [nanmean(LS'),nanmean(LC')]*SF;Le =[nansem(LS'),nansem(LC')]*SF;
Rm = [nanmean(RS'),nanmean(RC')]*SF;Re =[nansem(RS'),nansem(RC')]*SF;
tb_  = (1:sum(abs([tlimsSample tlimsChoice]))/0.05)*0.05;
% figure; hold on
% ciplot(Lm+Le-mean([Lm,Rm])-offset,Lm-Le-mean([Lm,Rm])-offset,tb_,[0.5 0 0.5],0.8);
% ciplot(Rm+Re-mean([Lm,Rm])-offset,Rm-Re-mean([Lm,Rm])-offset,tb_,[0 0 1],0.8);
ciplot(Lm+Le-offset,Lm-Le-offset,tb_,[0.5 0 0.5],0.8);
ciplot(Rm+Re-offset,Rm-Re-offset,tb_,[0 0 1],0.8);
%% Make raster + SDF  - Long ordered by Choice press latency

D.Long.ANOVA.p<0.05;
% d_L = t.Long.ChoicePress_LeftCorrect-t.Long.DelayEnd_LeftCorrect;
% d_R = t.Long.ChoicePress_RightCorrect-t.Long.DelayEnd_RightCorrect;
d_L = t.Long.ChoicePress_LeftCorrect-t.Long.NosePoke_LeftCorrect;
d_R = t.Long.ChoicePress_RightCorrect-t.Long.NosePoke_RightCorrect;
[~,idxL]=sort(d_L);[~,idxR]=sort(d_R);
NP_offset = [d_L(idxL);d_R(idxR)]*1e-6;
NP_offset(NP_offset>abs(tlimsChoice(1)))=NaN;


tlims_ = [t.Long.SamplePress_LeftCorrect(idxL) ; t.Long.SamplePress_RightCorrect(idxR)]*1e-6 + tlimsSample;

figure('color','w'); hold on

plot(sum(abs(tlimsSample))/2*[1,1],[-offset length(tlims_)+1],'color',[0 1 0 0.6],'LineWidth',2)
plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],[-offset length(tlims_)+1],'color',[1 0 0 0.6],'LineWidth',2)
plot(sum(abs(tlimsSample))*[1,1],[-offset length(tlims_)+1],'k')

for iTrial =1:length(tlims_)
    times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) ; 
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
end

tlims_ = [t.Long.ChoicePress_LeftCorrect(idxL) ; t.Long.ChoicePress_RightCorrect(idxR) ]*1e-6+ tlimsChoice;

for iTrial =1:length(tlims_)
    times_ = restrictFRrange(PFCcells{iUnit}.t*1e-4,tlims_(iTrial,:))' - tlims_(iTrial,1) + sum(abs(tlimsSample)); 
    plot([times_;times_],iTrial*ones(2,length(times_))+[-0.5; 0.5],'color',[0 0 0 0.6],'LineWidth',1.5)
    plot(-[NP_offset(iTrial)';NP_offset(iTrial)']+15,[iTrial iTrial]+[-0.5; 0.5],'color',[0.9 0.4 0.2],'LineWidth',1.5)
end

axis off
plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(t.Long.ChoicePress_LeftCorrect,1)],'color',[0.5 0 0.5],'LineWidth',1.5,'LineWidth',4)
plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 size(t.Long.ChoicePress_RightCorrect,1)]+size(t.Long.ChoicePress_LeftCorrect,1),'color',[0 0 1],'LineWidth',4)

text(1.5+sum(abs([tlimsSample tlimsChoice])),size(t.Long.ChoicePress_LeftCorrect,1)/2,'L','color',[0.5 0 0.5])
text(1.5+sum(abs([tlimsSample tlimsChoice])),size(t.Long.ChoicePress_LeftCorrect,1)+size(t.Long.ChoicePress_RightCorrect,1)/2,'R','color',[0 0 1])


data_{1}.iFR = iFR_;
[~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,t.Long.SamplePress_LeftCorrect*1e-6 );
LS = nan(sum(abs([tlimsSample]))/0.05,1);
for iTrial = 1:size(t.Long.SamplePress_LeftCorrect,1)
    try
    LS= [LS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,t.Long.ChoicePress_LeftCorrect*1e-6 );
LC = nan(sum(abs([tlimsChoice]))/0.05,1);
for iTrial = 1:size(t.Long.ChoicePress_LeftCorrect,1)
    try
    LC= [LC,cutout{1}{iTrial}(:,iUnit)];
    end
end


[~,cutout]=SelTimesFRs(tlimsSample,Tmtx,data_,t.Long.SamplePress_RightCorrect*1e-6 );
RS  = nan(sum(abs([tlimsSample]))/0.05,1);
for iTrial = 1:size(t.Long.SamplePress_RightCorrect,1)
    try
    RS= [RS,cutout{1}{iTrial}(:,iUnit)];
    end
end

[~,cutout]=SelTimesFRs(tlimsChoice,Tmtx,data_,t.Long.ChoicePress_RightCorrect*1e-6);
RC= nan(sum(abs([tlimsChoice]))/0.05,1);
for iTrial = 1:size(t.Long.ChoicePress_RightCorrect,1)
    try
    RC= [RC,cutout{1}{iTrial}(:,iUnit)];
    end
end


Lm = [nanmean(LS'),nanmean(LC')]*SF;Le =[nansem(LS'),nansem(LC')]*SF;
Rm = [nanmean(RS'),nanmean(RC')]*SF;Re =[nansem(RS'),nansem(RC')]*SF;
tb_  = (1:sum(abs([tlimsSample tlimsChoice]))/0.05)*0.05;
% figure; hold on
% ciplot(Lm+Le-mean([Lm,Rm])-offset,Lm-Le-mean([Lm,Rm])-offset,tb_,[0.5 0 0.5],0.8);
% ciplot(Rm+Re-mean([Lm,Rm])-offset,Rm-Re-mean([Lm,Rm])-offset,tb_,[0 0 1],0.8);
ciplot(Lm+Le-offset,Lm-Le-offset,tb_,[0.5 0 0.5],0.8);
ciplot(Rm+Re-offset,Rm-Re-offset,tb_,[0 0 1],0.8);
%% Make SDF by delay
delayTime = [4 8 16];
figure('color','w','name',sprintf('Unit %d',iUnit))
for iDelay =1:length(Delays_)
    subplot(2,1,1); hold on

    eval(sprintf('tlims_  = t.%s.SamplePress_LeftCorrect*1e-6;',Delays_{iDelay}));
    tlims_ = [closest(Tmtx,tlims_-5)];
    
    length_ = (1/0.05)*(delayTime(iDelay)+5);
    tb_ = (1:length_)*0.05-5;
    
    cutout = nan(length_,1);
    for iTrial=1:length(tlims_)
        cutout = [cutout, iFR_(tlims_(iTrial):tlims_(iTrial)+length_-1,iUnit)];
    end
    
    lm = nanmean(cutout,2);le = nansem(cutout,2);
    ciplot(lm+le,lm-le,tb_,'r')
%     plot(tb_,lm,'k','LineWidth',1.5)
    
    subplot(2,1,2); hold on

    eval(sprintf('tlims_  = t.%s.SamplePress_RightCorrect*1e-6;',Delays_{iDelay}));
    tlims_ = [closest(Tmtx,tlims_-5)];
    
    length_ = (1/0.05)*(delayTime(iDelay)+5);
    tb_ = (1:length_)*0.05-5;
    
    cutout = nan(length_,1);
    for iTrial=1:length(tlims_)
        cutout = [cutout, iFR_(tlims_(iTrial):tlims_(iTrial)+length_-1,iUnit)];
    end
    
    lm = nanmean(cutout,2);le = nansem(cutout,2);
    ciplot(lm+le,lm-le,tb_,'r')
%     plot(tb_,lm,'k','LineWidth',1.5)
end





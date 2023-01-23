
clear; iFile = 1;
p.target= 'LONG';

p.Delays_ = {'Short','Medium','Long'};
p.bw=0.05;

clear Av
warning ('off')
if ispc
    p.pat = 'C:\Analysis\AssemblyAnalysis\raw\';
else ismac
    p.pat = '/Volumes/HDD2/DNMTP/raw/';
end

cd(p.pat)
fileList=dir(['allTimestamps' filesep '*.mat']);
p.Areas = {'PFC','HP','HP-PFC'};
p.col_ = {'b','r','g'};

%%%%%

%%%%%%% Load an example
fname=strtok(fileList(iFile).name,'_');

load(sprintf('%sallTimestamps%s%s_Events.mat',p.pat,filesep,fname),'t');
A=load(sprintf('%s%s.mat',p.pat,fname));
spikeTimes{1}=A.PFCcells;
spikeTimes{2}=A.HPcells;
spikeTimes{3}=[spikeTimes{1};spikeTimes{2}];
%%
for iArea =1:2
    % load(sprintf('%s\\KDE_bins\\%s_%s_iFR50.mat',p.pat,fname,p.Areas{iArea}));
    A=load(sprintf('%sKDE_binsTaskonly%s%s_%s_iFR50_behavOnly.mat',p.pat,filesep,fname,p.Areas{iArea}));
    iFR{iArea}=A.iFR; 
%     iFR{iArea}=zscore(A.iFR);
end
iFR{3}=[iFR{1},iFR{2}];
Tmtx = A.Tmtx';
clear iArea
% figure; hold on
% plot(Tmtx(1:end-1),zscore(iFR{1}))
%%
fn=[p.pat 'KDE_binsTaskonly' filesep p.target  'Taskonly' filesep fname '_iFR50_behavOnly_'];       
load([fn 'AssemRes2'],'usel_out')
Ass = load([fn 'FSC'],'FSCsel','Tmtx','ciHsc','nu','units');
Ass.nu(3) = sum(Ass.nu(1:2));
Ass.usel_out = usel_out;
Ass.usel_out{3} = [Ass.usel_out{1},Ass.usel_out{2}+length(spikeTimes{1})];
UnitOffset = [0 Ass.nu(1) 0];
clear usel_out fn
%%
% plot(Ass.Tmtx,Ass.FSCsel{1})
%% Strip overlapping segments
[Ass.Tmtx,i]= unique(Ass.Tmtx);
for iArea=1:length(Ass.FSCsel)
    Ass.FSCsel{iArea} = Ass.FSCsel{iArea}(i,:);
end
% j = setdiff(1:length(Ass.Tmtx),i);
clear i iArea
%% mark up L and R trials, Choices/errors

t.LR = [ones(size(t.Short.CueLight_LeftCorrect));...
      2*ones(size(t.Short.CueLight_RightCorrect));...
        ones(size(t.Short.CueLight_LeftError'));...
      2*ones(size(t.Short.CueLight_RightError'));...
        ones(size(t.Medium.CueLight_LeftCorrect));...
      2*ones(size(t.Medium.CueLight_RightCorrect));...
        ones(size(t.Medium.CueLight_LeftError'));...
      2*ones(size(t.Medium.CueLight_RightError'));...
        ones(size(t.Long.CueLight_LeftCorrect));...
      2*ones(size(t.Long.CueLight_RightCorrect));...
        ones(size(t.Long.CueLight_LeftError'));...
      2*ones(size(t.Long.CueLight_RightError'))];
    
t.CE = [ones(size(t.Short.CueLight_LeftCorrect));...
        ones(size(t.Short.CueLight_RightCorrect));...
        zeros(size(t.Short.CueLight_LeftError'));...
        zeros(size(t.Short.CueLight_RightError'));...
        ones(size(t.Medium.CueLight_LeftCorrect));...
        ones(size(t.Medium.CueLight_RightCorrect));...
        zeros(size(t.Medium.CueLight_LeftError'));...
        zeros(size(t.Medium.CueLight_RightError'));...
        ones(size(t.Long.CueLight_LeftCorrect));...
        ones(size(t.Long.CueLight_RightCorrect));...
        zeros(size(t.Long.CueLight_LeftError'));...
        zeros(size(t.Long.CueLight_RightError'))];  
%% Plot each trial aligned - Long ordered by Choice press latency
iArea = 3;

offset=0;
SF=20;
tlimsSample = [-5 5];
tlimsChoice = [-5 5];
% subplot(1,3,3);
figure('color','w'); hold on

% Choice press latency
% d_L = t.Long.ChoicePress_LeftCorrect-t.Long.DelayEnd_LeftCorrect;
% d_R = t.Long.ChoicePress_RightCorrect-t.Long.DelayEnd_RightCorrect;
d_L = t.Long.ChoicePress_LeftCorrect-t.Long.NosePoke_LeftCorrect;
d_R = t.Long.ChoicePress_RightCorrect-t.Long.NosePoke_RightCorrect;
[~,idxL]=sort(d_L);[~,idxR]=sort(d_R);
NP_offset = [d_L(idxL);d_R(idxR)]*1e-6;
NP_offset(NP_offset>abs(tlimsChoice(1)))=NaN;
% Sample press latency
d_L = t.Long.SamplePress_LeftCorrect-t.Long.CueLight_LeftCorrect;
d_R = t.Long.SamplePress_RightCorrect-t.Long.CueLight_RightCorrect;
SP_offset = [d_L(idxL);d_R(idxR)]*1e-6;
% SP_offset(SP_offset>abs(tlimsSample(1)))=NaN;

tlims_ = [t.Long.SamplePress_LeftCorrect(idxL) ; t.Long.SamplePress_RightCorrect(idxR)]*1e-6 ;
% tlims_ = [t.Medium.SamplePress_LeftCorrect(idxL) ; t.Medium.SamplePress_RightCorrect(idxR)]*1e-6 ;
tlims_ =repmat(tlims_,1,2) + repmat(tlimsSample,length(tlims_),1);

plot(sum(abs(tlimsSample))/2*[1,1],SF*[-offset length(tlims_)+1],'color',[0 1 0 0.6],'LineWidth',2)
plot((sum(abs(tlimsChoice))/2+ sum(abs(tlimsSample)))*[1,1],SF*[-offset length(tlims_)+1],'color',[1 0 0 0.6],'LineWidth',2)
plot(sum(abs(tlimsSample))*[1,1],SF*[-offset length(tlims_)+1],'k')
idx=[];

for iTrial =1:length(tlims_)
    try
        idx = closest(Ass.Tmtx,tlims_(iTrial,:));
    catch
        for i=1:size(tlims_,2)
            idx(i) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,i));
        end
    end
    idx = idx(1):idx(2);      
%     idx = idx(1):(idx(1)+199);
    times_ = Ass.Tmtx(idx) - tlims_(iTrial,1); 
     
    for iAss = 1:size(Ass.FSCsel{iArea},2)
        plot(times_,iTrial*SF+Ass.FSCsel{iArea}(idx,iAss),'color',[0 0 0 0.6],'LineWidth',1.5)
    end
end

tlims_ = [t.Long.ChoicePress_LeftCorrect(idxL) ; t.Long.ChoicePress_RightCorrect(idxR)]*1e-6 ;
tlims_ =repmat(tlims_,1,2) + repmat(tlimsChoice,length(tlims_),1);

for iTrial =1:length(tlims_)
    try
        idx = closest(Ass.Tmtx,tlims_(iTrial,:));
    catch
        for i=1:size(tlims_,2)
            idx(i) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,i));
        end
    end
    idx = idx(1):idx(2);      
%     idx = idx(1):(idx(1)+199);
    times_ = Ass.Tmtx(idx)-tlims_(iTrial,1)+ sum(abs(tlimsSample));
%     times_ = times_-times_(1)
    for iAss = 1:size(Ass.FSCsel{iArea},2)
        plot(times_,iTrial*SF+Ass.FSCsel{iArea}(idx,iAss),'color',[0 0 0 0.6],'LineWidth',1.5)
    end
    % Add cue times
    plot(5-[SP_offset(iTrial)';SP_offset(iTrial)'],...
         SF*[iTrial iTrial]+[-0.5, 0.5]*SF,'color',[0.1 0.9 0.1],'LineWidth',1.5)
     
    % Add Nose poke times
    plot(-[NP_offset(iTrial)';NP_offset(iTrial)']+15,...
         SF*[iTrial iTrial]+[-0.5, 0.5]*SF,'color',[0.9 0.4 0.2],'LineWidth',1.5)
    
end


plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 SF*size(t.Long.ChoicePress_LeftCorrect,1)],'color',[0.5 0 0.5],'LineWidth',1.5,'LineWidth',4)
plot(1+sum(abs([tlimsSample tlimsChoice]))*[1,1],[1 SF*size(t.Long.ChoicePress_RightCorrect,1)]+SF*size(t.Long.ChoicePress_LeftCorrect,1),'color',[0 0 1],'LineWidth',4)

text(1.5+sum(abs([tlimsSample tlimsChoice])),SF*size(t.Long.ChoicePress_LeftCorrect,1)/2,'L','color',[0.5 0 0.5])
text(1.5+sum(abs([tlimsSample tlimsChoice])),SF*size(t.Long.ChoicePress_LeftCorrect,1)+SF*size(t.Long.ChoicePress_RightCorrect,1)/2,'R','color',[0 0 1])
plot([22 24],SF*[1 1],'LineWidth',1.5,'color','k')
plot([22 22],SF*[1 2],'LineWidth',1.5,'color','k')
axis([0 sum(abs([tlimsChoice, tlimsSample]))+5 0 SF*length(tlims_)+50])
axis off

clear NP_offset d_L d_R iArea iAss idxL idxR iTrial offset SF times_ tlims_ tlimsChoice tlimsSample

%% Plot all activations on continuous timesale with units overlaid
SF = 20;
bounds_ = 10;
iArea = 3;
 tlims_ = [...
        t.Short.CueLight_LeftCorrect,  t.Short.SamplePress_LeftCorrect,    t.Short.ChoicePress_LeftCorrect;...
        t.Short.CueLight_RightCorrect, t.Short.SamplePress_RightCorrect,   t.Short.ChoicePress_RightCorrect;...
        t.Short.CueLight_LeftError',   t.Short.SamplePress_LeftError',     t.Short.ChoicePress_LeftError';...
        t.Short.CueLight_RightError',  t.Short.SamplePress_RightError',    t.Short.ChoicePress_RightError';...
        t.Medium.CueLight_LeftCorrect, t.Medium.SamplePress_LeftCorrect,   t.Medium.ChoicePress_LeftCorrect;...
        t.Medium.CueLight_RightCorrect,t.Medium.SamplePress_RightCorrect,  t.Medium.ChoicePress_RightCorrect;...
        t.Medium.CueLight_LeftError',  t.Medium.SamplePress_LeftError',    t.Medium.ChoicePress_LeftError';...
        t.Medium.CueLight_RightError', t.Medium.SamplePress_RightError',   t.Medium.ChoicePress_RightError';...
        t.Long.CueLight_LeftCorrect,   t.Long.SamplePress_LeftCorrect,     t.Long.ChoicePress_LeftCorrect;...
        t.Long.CueLight_RightCorrect,  t.Long.SamplePress_RightCorrect,    t.Long.ChoicePress_RightCorrect;...
        t.Long.CueLight_LeftError',    t.Long.SamplePress_LeftError',      t.Long.ChoicePress_LeftError';...
        t.Long.CueLight_RightError',   t.Long.SamplePress_RightError',     t.Long.ChoicePress_RightError']*1e-6;
 
%         tlims_=sortrows(tlims_,1);  
    
    

figure; hold on
% shaded blue regions to mark trial times
%         for iTrial =1:size(tlims_,1) 
%             try
% %             r= rectangle('position',[tlims_(iTrial,1)-10,5 tlims_(iTrial,3)-tlims_(iTrial,1)+20 size(Ass.FSCsel{iArea},2)*SF  ])
% %             r.FaceColor= 'b';
% %             r.LineStyle='none';
%             patch([tlims_(iTrial,1)-bounds_,tlims_(iTrial,3)+bounds_,tlims_(iTrial,3)+bounds_,tlims_(iTrial,1)-bounds_],...
%                    [5 5 size(Ass.FSCsel{iArea},2)*SF+5  size(Ass.FSCsel{iArea},2)*SF+5],[0 0 1],'FaceAlpha',0.6,'LineStyle','none')
%             catch
%             end
%         end
        
for iAss=1:size(Ass.FSCsel{iArea},2)
    units_=Ass.usel_out{iArea}(Ass.units{iArea}{iAss});
    iFR_ = iFR{iArea}(:,units_);
% plot(Tmtx(1:end-1),zscore(iFR_)+iAss*SF,'color',0.6*[1 1 1],'LineWidth',1);
try
    plot(Tmtx(1:end-1),zscore(iFR_(:,units_<=Ass.nu(1)))+iAss*SF,'color',[0 0 1],'LineWidth',1);
end
try
    plot(Tmtx(1:end-1),zscore(iFR_(:,units_>Ass.nu(1)))+iAss*SF,'color',[1 0 0],'LineWidth',1);
end 
for iTrial =1:size(tlims_,1) 
    try
    idx=[];
    try
        idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
    catch
        idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
        idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
    end
    idx = idx(1):idx(2);      
    plot(Ass.Tmtx(idx),Ass.FSCsel{iArea}(idx,iAss)+iAss*SF,'color','k','LineWidth',2)
    plot([tlims_(iTrial,2);tlims_(iTrial,2)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[0 1 0 0.6],'LineWidth',2)
    plot([tlims_(iTrial,3);tlims_(iTrial,3)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[1 0 0 0.6],'LineWidth',2)
    catch
        disp(sprintf('No timestamps found for trial %G',iTrial))
    end
end
end

% max_ = SF*size(Ass.FSCsel{iArea},2);

set(gca,'xtick',[min(tlims_)-bounds_:60:max(tlims_)+bounds_])
grid on

% plot(Ass.Tmtx(:),Ass.FSCsel{iArea}(:,iAss),'color','k','LineWidth',1.2)
% plot(Ass.Tmtx(j),Ass.FSCsel{iArea}(j,iAss)-10,'color',0.6*[1 0 0],'LineWidth',1.2)
 clear bounds_ max_ iTrial iArea tlims_ SF iAss idx units_
%% Plot all activations on continuous timesale with units overlaid - one example assembly showing non-member units 
SF = 1.2;
SF_Ass = 5;
bounds_ = 10;
Offset = 4;
plotNonMembers = false;
tExtreme = [5240 5490];
% tExtreme = [5000 8000];
shading_ = [0.3 0.3];
TrialSkip = 2;


for iArea = 3%1:3;
 tlims_ = [...
        t.Short.CueLight_LeftCorrect,  t.Short.SamplePress_LeftCorrect,    t.Short.ChoicePress_LeftCorrect;...
        t.Short.CueLight_RightCorrect, t.Short.SamplePress_RightCorrect,   t.Short.ChoicePress_RightCorrect;...
        t.Short.CueLight_LeftError',   t.Short.SamplePress_LeftError',     t.Short.ChoicePress_LeftError';...
        t.Short.CueLight_RightError',  t.Short.SamplePress_RightError',    t.Short.ChoicePress_RightError';...
        t.Medium.CueLight_LeftCorrect, t.Medium.SamplePress_LeftCorrect,   t.Medium.ChoicePress_LeftCorrect;...
        t.Medium.CueLight_RightCorrect,t.Medium.SamplePress_RightCorrect,  t.Medium.ChoicePress_RightCorrect;...
        t.Medium.CueLight_LeftError',  t.Medium.SamplePress_LeftError',    t.Medium.ChoicePress_LeftError';...
        t.Medium.CueLight_RightError', t.Medium.SamplePress_RightError',   t.Medium.ChoicePress_RightError';...
        t.Long.CueLight_LeftCorrect,   t.Long.SamplePress_LeftCorrect,     t.Long.ChoicePress_LeftCorrect;...
        t.Long.CueLight_RightCorrect,  t.Long.SamplePress_RightCorrect,    t.Long.ChoicePress_RightCorrect;...
        t.Long.CueLight_LeftError',    t.Long.SamplePress_LeftError',      t.Long.ChoicePress_LeftError';...
        t.Long.CueLight_RightError',   t.Long.SamplePress_RightError',     t.Long.ChoicePress_RightError']*1e-6;
           
        [tlims_ ,I]=sortrows(tlims_,1);  
        LR_ = t.LR(I');

LR_(min(tlims_,[],2)<tExtreme(1) | max(tlims_,[],2)>tExtreme(2)) = [];
tlims_(tlims_<tExtreme(1) | tlims_>tExtreme(2)) = NaN;
tlims_(sum(isnan(tlims_),2)>0,:)=[];



idx_ = closest(Tmtx,tExtreme);
Tmtx_ = Tmtx(idx_(1):idx_(2));  


AssCol = {[0.081 0.09 0.631],[0.761 0.505 0.046],[0.2,0.67,0.4]};
% [0.16,0.36,0.19]
% [0.2,0.67,0.4]
for iAss=1%:size(Ass.FSCsel{iArea},2)
    figure('color','w'); hold on

    %%%%%%%%%%%% plot the assemblies
    for iTrial =1:TrialSkip:size(tlims_,1)
        try
            idx=[];
            try
                idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
            catch
                idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
                idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
            end
            idx = idx(1):idx(2);
            temp_ = Ass.FSCsel{iArea}(:,iAss);
            temp_ = mat2gray(temp_);
            temp_ = temp_(idx);
            temp_([1,end]) =0;%min(temp_([1,end]));
            
%             patch(Ass.Tmtx(idx),Offset*temp_,p.col_{iArea},'EdgeColor','none','FaceAlpha',0.8)
            patch(Ass.Tmtx(idx),SF_Ass*temp_,AssCol{iArea},'EdgeColor',AssCol{iArea},'FaceAlpha',0.8)
            
            %     %%%%%%%%%%%% plot the lever presses
            %     plot([tlims_(iTrial,2);tlims_(iTrial,2)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[0 1 0 0.6],'LineWidth',2)
            %     plot([tlims_(iTrial,3);tlims_(iTrial,3)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[1 0 0 0.6],'LineWidth',2)
        catch
            disp(sprintf('No timestamps found for trial %G',iTrial))
        end
    end
    
    
    %%%%%%%%%%%% Plot the units

        
        iFR_ = iFR{3};%(:,units_Real);
        iFR_ = zscore(iFR_);
        
        iFR_ = bsxfun(@minus,iFR_,iFR_(1,:));
        iFR_ = iFR_(idx_(1):idx_(2),:);
        
        units_ = Ass.units{iArea}{iAss};
        %
        switch iArea
            case 1
                orderedUnits = [setdiff(1:Ass.nu(1),units_),...
                                units_ ];
            case 2
                orderedUnits = [units_, ...
                                setdiff((Ass.nu(1)+1):Ass.nu(3),units_)];
            case 3
                orderedUnits = [setdiff(1:Ass.nu(1),units_),...
                                units_ , ...
                                setdiff((Ass.nu(1)+1):Ass.nu(3),units_)];
        end
        %units_Real = Ass.usel_out{3}(Ass.units{iArea}{iAss}+UnitOffset(iArea));
        units_Real = Ass.usel_out{3}(orderedUnits);
        for iUnit =1:length(orderedUnits)
            temp_ = mat2gray(iFR_(:,units_Real(iUnit)));
            temp_([1,end]) = 0;%min(temp_([1,end]));
             ShadingLevel = shading_(ismember(orderedUnits(iUnit),units_)+1);
             if ismember(orderedUnits(iUnit),units_)
            if units_Real(iUnit) <= Ass.usel_out{1}(Ass.nu(1))
                patch(Tmtx_,Offset+temp_+iUnit*SF,[0.081 0.09 0.631],'EdgeColor',[0.081 0.09 0.631],'FaceAlpha',ShadingLevel)
            else
                patch(Tmtx_,Offset+temp_+iUnit*SF,[0.761 0.505 0.046],'EdgeColor',[0.761 0.505 0.046],'FaceAlpha',ShadingLevel)
            end
             else
                patch(Tmtx_,Offset+temp_+iUnit*SF,[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6],'FaceAlpha',ShadingLevel)
             end
        end

%%%%%%%%%%%% plot the lever presses
for iTrial =1:TrialSkip:size(tlims_,1) 
    try
    idx=[];
    try
        idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
    catch
        idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
        idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
    end
    idx = idx(1):idx(2);      
    
 
    plot([tlims_(iTrial,2);tlims_(iTrial,2)],[-1 -0.2]*SF,'color',[0 1 0 0.6],'LineWidth',2)
    plot([tlims_(iTrial,3);tlims_(iTrial,3)],[-1 -0.2]*SF,'color',[1 0 0 0.6],'LineWidth',2)
    switch LR_(iTrial)
        case 1
            text(mean(tlims_(iTrial,2:3)),-2,'L','HorizontalAlignment','center')
        case 2
            text(mean(tlims_(iTrial,2:3)),-2,'R','HorizontalAlignment','center')
    end
    
    end
end
% max_ = SF*size(Ass.FSCsel{iArea},2);

set(gca,'xtick',[min(tlims_)-bounds_:60:max(tlims_)+bounds_])
grid on
axis([tExtreme -2 Inf])
axis off
end
end
%% Plot all activations on continuous timesale with units overlaid - one example assembly 
SF = 1.2;
SF_Ass = 2;
bounds_ = 10;
Offset = 1;
plotNonMembers = false;

AssCol = {[0.081 0.09 0.631],[0.761 0.505 0.046],[0.2,0.67,0.4]};
% [0.16,0.36,0.19]
% [0.2,0.67,0.4]
tlims_ = [...
    t.Short.CueLight_LeftCorrect,  t.Short.SamplePress_LeftCorrect,    t.Short.ChoicePress_LeftCorrect;...
    t.Short.CueLight_RightCorrect, t.Short.SamplePress_RightCorrect,   t.Short.ChoicePress_RightCorrect;...
    t.Short.CueLight_LeftError',   t.Short.SamplePress_LeftError',     t.Short.ChoicePress_LeftError';...
    t.Short.CueLight_RightError',  t.Short.SamplePress_RightError',    t.Short.ChoicePress_RightError';...
    t.Medium.CueLight_LeftCorrect, t.Medium.SamplePress_LeftCorrect,   t.Medium.ChoicePress_LeftCorrect;...
    t.Medium.CueLight_RightCorrect,t.Medium.SamplePress_RightCorrect,  t.Medium.ChoicePress_RightCorrect;...
    t.Medium.CueLight_LeftError',  t.Medium.SamplePress_LeftError',    t.Medium.ChoicePress_LeftError';...
    t.Medium.CueLight_RightError', t.Medium.SamplePress_RightError',   t.Medium.ChoicePress_RightError';...
    t.Long.CueLight_LeftCorrect,   t.Long.SamplePress_LeftCorrect,     t.Long.ChoicePress_LeftCorrect;...
    t.Long.CueLight_RightCorrect,  t.Long.SamplePress_RightCorrect,    t.Long.ChoicePress_RightCorrect;...
    t.Long.CueLight_LeftError',    t.Long.SamplePress_LeftError',      t.Long.ChoicePress_LeftError';...
    t.Long.CueLight_RightError',   t.Long.SamplePress_RightError',     t.Long.ChoicePress_RightError']*1e-6;
        
[tlims_ ,I]=sortrows(tlims_,1);
LR_ = t.LR(I');
CE_ = t.CE(I');

tExtreme = [min(min(tlims_)) max(max(tlims_))];
% tExtreme = [5240 5490];

LR_(min(tlims_,[],2)<tExtreme(1) | max(tlims_,[],2)>tExtreme(2)) = [];
CE_(min(tlims_,[],2)<tExtreme(1) | max(tlims_,[],2)>tExtreme(2)) = [];


tlims_(tlims_<tExtreme(1) | tlims_>tExtreme(2)) = NaN;
tlims_(sum(isnan(tlims_),2)>0,:)=[];


idx_ = closest(Tmtx,tExtreme);
Tmtx_ = Tmtx(idx_(1):idx_(2));


TrialSkip = 1;

for iArea = 3%1:3;

for iAss=1:size(Ass.FSCsel{iArea},2)
    figure('color','w'); hold on

    %%%%%%%%%%%% plot the assemblies
    for iTrial =1:TrialSkip:size(tlims_,1)
        try
            idx=[];
            try
                idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
            catch
                idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
                idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
            end
            idx = idx(1):idx(2);
            temp_ = Ass.FSCsel{iArea}(:,iAss);
            temp_ = mat2gray(temp_);
            temp_ = temp_(idx);
            temp_([1,end]) =0;%min(temp_([1,end]));
            
%             patch(Ass.Tmtx(idx),Offset*temp_,p.col_{iArea},'EdgeColor','none','FaceAlpha',0.8)
            patch(Ass.Tmtx(idx),SF_Ass*temp_,AssCol{iArea},'EdgeColor',AssCol{iArea},'FaceAlpha',1)
%             plot(Ass.Tmtx(idx),SF_Ass*temp_,'Color',AssCol{iArea},'Linewidth',1.5)
            %     %%%%%%%%%%%% plot the lever presses
            %     plot([tlims_(iTrial,2);tlims_(iTrial,2)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[0 1 0 0.6],'LineWidth',2)
            %     plot([tlims_(iTrial,3);tlims_(iTrial,3)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[1 0 0 0.6],'LineWidth',2)
        catch
            disp(sprintf('No timestamps found for trial %G',iTrial))
        end
    end
    
    
    %%%%%%%%%%%% Plot the units

        units_Real = Ass.usel_out{3}(Ass.units{iArea}{iAss}+UnitOffset(iArea));
        iFR_ = iFR{3}(:,units_Real);
        iFR_ = zscore(iFR_);
        
        iFR_ = bsxfun(@minus,iFR_,iFR_(1,:));
        iFR_ = iFR_(idx_(1):idx_(2),:);
        
        units_ = Ass.units{iArea}{iAss};

        for iUnit =1:length(units_)
            temp_ = mat2gray(iFR_(:,iUnit));
            temp_([1,end]) = 0;%min(temp_([1,end]));
            if units_Real(iUnit) <= Ass.usel_out{1}(Ass.nu(1))
                patch(Tmtx_,Offset+temp_+iUnit*SF,[0.081 0.09 0.631],'EdgeColor',[0.081 0.09 0.631],'FaceAlpha',0.3)
            else
                patch(Tmtx_,Offset+temp_+iUnit*SF,[0.761 0.505 0.046],'EdgeColor',[0.761 0.505 0.046],'FaceAlpha',0.3)
            end
        end

%%%%%%%%%%%% plot the lever presses
for iTrial =1:TrialSkip:size(tlims_,1) 
    try
    idx=[];
    try
        idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
    catch
        idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
        idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
    end
    idx = idx(1):idx(2);      
    
 
    plot([tlims_(iTrial,2);tlims_(iTrial,2)],[-1 -0.2]*SF,'color',[0 1 0 0.6],'LineWidth',2)
    plot([tlims_(iTrial,3);tlims_(iTrial,3)],[-1 -0.2]*SF,'color',[1 0 0 0.6],'LineWidth',2)
    
     plot([tlims_(iTrial,2);mean(tlims_(iTrial,2:3))],[-1 -1.5]*SF,'k')
     plot([mean(tlims_(iTrial,2:3)),tlims_(iTrial,3)],[-1.5 -1]*SF,'k')
   
    switch CE_(iTrial)
        case 0
            text(mean(tlims_(iTrial,2:3)),-2,'Correct','HorizontalAlignment','center')
        case 1
            text(mean(tlims_(iTrial,2:3)),-2,'*Error*','HorizontalAlignment','center')
    end
    switch LR_(iTrial)
        case 1
            text(mean(tlims_(iTrial,2:3)),-2.8,'(L)','HorizontalAlignment','center')
        case 2
            text(mean(tlims_(iTrial,2:3)),-2.8,'(R)','HorizontalAlignment','center')
    end
    end
end
% max_ = SF*size(Ass.FSCsel{iArea},2);

set(gca,'xtick',[min(tlims_)-bounds_:60:max(tlims_)+bounds_])
grid on
axis([tExtreme -3.5 22])
% axis off
end
end
%% Plot all activations on continuous timescale with units overlaid - one example assembly with non-member units included
SF = 0.1;
bounds_ = 10;
Offset = 4;
plotNonMembers = false;
iArea = 3;
 tlims_ = [...
        t.Short.CueLight_LeftCorrect,  t.Short.SamplePress_LeftCorrect,    t.Short.ChoicePress_LeftCorrect;...
        t.Short.CueLight_RightCorrect, t.Short.SamplePress_RightCorrect,   t.Short.ChoicePress_RightCorrect;...
        t.Short.CueLight_LeftError',   t.Short.SamplePress_LeftError',     t.Short.ChoicePress_LeftError';...
        t.Short.CueLight_RightError',  t.Short.SamplePress_RightError',    t.Short.ChoicePress_RightError';...
        t.Medium.CueLight_LeftCorrect, t.Medium.SamplePress_LeftCorrect,   t.Medium.ChoicePress_LeftCorrect;...
        t.Medium.CueLight_RightCorrect,t.Medium.SamplePress_RightCorrect,  t.Medium.ChoicePress_RightCorrect;...
        t.Medium.CueLight_LeftError',  t.Medium.SamplePress_LeftError',    t.Medium.ChoicePress_LeftError';...
        t.Medium.CueLight_RightError', t.Medium.SamplePress_RightError',   t.Medium.ChoicePress_RightError';...
        t.Long.CueLight_LeftCorrect,   t.Long.SamplePress_LeftCorrect,     t.Long.ChoicePress_LeftCorrect;...
        t.Long.CueLight_RightCorrect,  t.Long.SamplePress_RightCorrect,    t.Long.ChoicePress_RightCorrect;...
        t.Long.CueLight_LeftError',    t.Long.SamplePress_LeftError',      t.Long.ChoicePress_LeftError';...
        t.Long.CueLight_RightError',   t.Long.SamplePress_RightError',     t.Long.ChoicePress_RightError']*1e-6;
 
%         tlims_=sortrows(tlims_,1);  
    
tExtreme = [4500 4600];
tlims_(tlims_<tExtreme(1) | tlims_>tExtreme(2)) = NaN;
tlims_(sum(isnan(tlims_),2)>0,:)=[];

idx_ = closest(Tmtx,tExtreme);
Tmtx_ = Tmtx(idx_(1):idx_(2));

figure; hold on
% shaded blue regions to mark trial times
%         for iTrial =1:size(tlims_,1) 
%             try
% %             r= rectangle('position',[tlims_(iTrial,1)-10,5 tlims_(iTrial,3)-tlims_(iTrial,1)+20 size(Ass.FSCsel{iArea},2)*SF  ])
% %             r.FaceColor= 'b';
% %             r.LineStyle='none';
%             patch([tlims_(iTrial,1)-bounds_,tlims_(iTrial,3)+bounds_,tlims_(iTrial,3)+bounds_,tlims_(iTrial,1)-bounds_],...
%                    [5 5 size(Ass.FSCsel{iArea},2)*SF+5  size(Ass.FSCsel{iArea},2)*SF+5],[0 0 1],'FaceAlpha',0.6,'LineStyle','none')
%             catch
%             end
%         end

for iAss=1%:size(Ass.FSCsel{iArea},2)
    
    %%%%%%%%%%%% plot the assemblies
    for iTrial =1:size(tlims_,1)
        try
            idx=[];
            try
                idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
            catch
                idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
                idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
            end
            idx = idx(1):idx(2);
            temp_ = Ass.FSCsel{iArea}(:,iAss);
            temp_ = mat2gray(temp_);
            temp_ = temp_(idx);
            temp_([1,end]) =0;%min(temp_([1,end]));
            
            patch(Ass.Tmtx(idx),Offset*temp_,p.col_{iArea},'EdgeColor','none','FaceAlpha',0.8)
            
            %     %%%%%%%%%%%% plot the lever presses
            %     plot([tlims_(iTrial,2);tlims_(iTrial,2)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[0 1 0 0.6],'LineWidth',2)
            %     plot([tlims_(iTrial,3);tlims_(iTrial,3)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[1 0 0 0.6],'LineWidth',2)
        catch
            disp(sprintf('No timestamps found for trial %G',iTrial))
        end
    end
    
    
    %%%%%%%%%%%% Plot the units
    if plotNonMembers
        units_=Ass.usel_out{iArea}(Ass.units{iArea}{iAss});
        iFR_ = iFR{iArea}(:,:);
        iFR_ = zscore(iFR_);
        
        iFR_ = bsxfun(@minus,iFR_,iFR_(1,:));
        iFR_ = iFR_(idx_(1):idx_(2),:);
        iFR_ = mat2gray(iFR_);
        for iUnit =1:size(iFR_,2)
            temp_ = iFR_(:,iUnit);
            temp_([1,end]) =0;%min(temp_([1,end]));
            if ismember(iUnit,units_)
                if iUnit<=Ass.nu(1)
                     plot(Tmtx_,Offset-2+temp_+iUnit*SF,'color',[0.081 0.09 0.631],'LineWidth',1.5)
                else
                     plot(Tmtx_,Offset-2+temp_+iUnit*SF,'color',[0.761 0.505 0.046],'LineWidth',1.5)
                end
            else                
            	plot(Tmtx_,Offset-2+temp_+iUnit*SF,'color',0.8*[1 1 1 1])
            end
        end
    else
        units_ = Ass.usel_out{iArea}(Ass.units{iArea}{iAss});
        iFR_ = iFR{iArea}(:,units_);
        iFR_ = zscore(iFR_);
        
        iFR_ = bsxfun(@minus,iFR_,iFR_(1,:));
        iFR_ = iFR_(idx_(1):idx_(2),:);

        for iUnit =1:length(units_)
            temp_ = mat2gray(iFR_(:,iUnit));
            temp_([1,end]) =0;%min(temp_([1,end]));
            if units_(iUnit)<=Ass.nu(1)
                patch(Tmtx_,Offset-2+temp_+iUnit*SF,[0.081 0.09 0.631],'EdgeColor','none','FaceAlpha',1)%+iAss*SF
            else
                patch(Tmtx_,Offset-2+temp_+iUnit*SF,[0.761 0.505 0.046],'EdgeColor','none','FaceAlpha',1)%+iAss*SF
            end
        end
    end
    
end
%%%%%%%%%%%% plot the lever presses
for iTrial =1:size(tlims_,1) 
    try
    idx=[];
    try
        idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
    catch
        idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
        idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
    end
    idx = idx(1):idx(2);      
    
 
    plot([tlims_(iTrial,2);tlims_(iTrial,2)],[-1 0]*SF,'color',[0 1 0 0.6],'LineWidth',2)
    plot([tlims_(iTrial,3);tlims_(iTrial,3)],[-1 0]*SF,'color',[1 0 0 0.6],'LineWidth',2)
    end
end
% max_ = SF*size(Ass.FSCsel{iArea},2);

set(gca,'xtick',[min(tlims_)-bounds_:60:max(tlims_)+bounds_])
grid on

% plot(Ass.Tmtx(:),Ass.FSCsel{iArea}(:,iAss),'color','k','LineWidth',1.2)
% plot(Ass.Tmtx(j),Ass.FSCsel{iArea}(j,iAss)-10,'color',0.6*[1 0 0],'LineWidth',1.2)
%  clear bounds_ max_ iTrial iArea tlims_ SF iAss idx units_
%% Plot all activations on continuous timescale with units overlaid - all assemblies
SF = 1;
bounds_ = 10;
Offset = 4;
plotNonMembers = true;
 tlims_ = [...
        t.Short.CueLight_LeftCorrect,  t.Short.SamplePress_LeftCorrect,    t.Short.ChoicePress_LeftCorrect;...
        t.Short.CueLight_RightCorrect, t.Short.SamplePress_RightCorrect,   t.Short.ChoicePress_RightCorrect;...
        t.Short.CueLight_LeftError',   t.Short.SamplePress_LeftError',     t.Short.ChoicePress_LeftError';...
        t.Short.CueLight_RightError',  t.Short.SamplePress_RightError',    t.Short.ChoicePress_RightError';...
        t.Medium.CueLight_LeftCorrect, t.Medium.SamplePress_LeftCorrect,   t.Medium.ChoicePress_LeftCorrect;...
        t.Medium.CueLight_RightCorrect,t.Medium.SamplePress_RightCorrect,  t.Medium.ChoicePress_RightCorrect;...
        t.Medium.CueLight_LeftError',  t.Medium.SamplePress_LeftError',    t.Medium.ChoicePress_LeftError';...
        t.Medium.CueLight_RightError', t.Medium.SamplePress_RightError',   t.Medium.ChoicePress_RightError';...
        t.Long.CueLight_LeftCorrect,   t.Long.SamplePress_LeftCorrect,     t.Long.ChoicePress_LeftCorrect;...
        t.Long.CueLight_RightCorrect,  t.Long.SamplePress_RightCorrect,    t.Long.ChoicePress_RightCorrect;...
        t.Long.CueLight_LeftError',    t.Long.SamplePress_LeftError',      t.Long.ChoicePress_LeftError';...
        t.Long.CueLight_RightError',   t.Long.SamplePress_RightError',     t.Long.ChoicePress_RightError']*1e-6;
 
        tlims_=sortrows(tlims_,1);
tExtreme = [4580 5200];
% tExtreme = [0 Inf];
tlims_(tlims_<tExtreme(1) | tlims_>tExtreme(2)) = NaN;
tlims_(sum(isnan(tlims_),2)>0,:)=[];

idx_ = closest(Tmtx,tExtreme);
Tmtx_ = Tmtx(idx_(1):idx_(2));

GlobalMembers = unique([cell2mat(Ass.units{1}),cell2mat(Ass.units{1})+Ass.nu(1),cell2mat(Ass.units{3})]);
GlobalNonMembers = setdiff(1:Ass.nu(3),GlobalMembers);

colOrder = [3 1 2];

figure('color','w'); hold on

for iArea=1:length(Ass.FSCsel)
    for iAss=1:3%size(Ass.FSCsel{iArea},2)
        col_ = min(rand(1,3),0.1*[1 1 1]); col_(colOrder(iArea))=0.9;
        Ass_ = mat2gray(Ass.FSCsel{iArea}(:,iAss));

        %%%%%%%%%%%% plot the assemblies
        for iTrial =1:2:size(tlims_,1)
            try
                idx=[];
                try
                    idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
                catch
                    idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
                    idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
                end
                idx = idx(1):idx(2);
                temp_ = Ass_(idx) ;

                temp_([1 end]) =0;%min(temp_([1,end]));
                patch(Ass.Tmtx(idx),temp_+iArea-1,col_,'EdgeColor',col_,'FaceAlpha',0.1);
%                 plot(Ass.Tmtx(idx),temp_+iArea-1,'color',col_,'LineWidth',1.2);

            end
        end

        
        
        %%%%%%%%%%%% Plot the member units
        units_ = Ass.usel_out{iArea}(Ass.units{iArea}{iAss});
        
        iFR_ = iFR{iArea}(:,units_);
        iFR_ = zscore(iFR_);
        iFR_ = bsxfun(@minus,iFR_,iFR_(1,:));
        iFR_ = mat2gray(iFR_);
        iFR_ = iFR_(idx_(1):idx_(2),:);
        
        units_ = Ass.units{iArea}{iAss};
        
        for iUnit =1:size(iFR_,2)
            temp_ = (iFR_(:,iUnit));
            temp_([1,end]) =0;
            patch(Tmtx_,Offset + temp_ + (UnitOffset(iArea)+units_(iUnit))*SF,col_,'EdgeColor','none','FaceAlpha',0.6)
        end
        
    end
end
axis off
 %%%%%%%%%% Plot the global nonmember units
units_ = Ass.usel_out{3}(GlobalNonMembers);
iFR_ = iFR{3}(:,units_);
iFR_ = zscore(iFR_);
iFR_ = bsxfun(@minus,iFR_,iFR_(1,:));
iFR_ = iFR_(idx_(1):idx_(2),:);

for iUnit =1:size(iFR_,2)
    temp_ = mat2gray(iFR_(:,iUnit));
    temp_([1,end]) =0;
    patch(Tmtx_,Offset+temp_ +units_(iUnit)*SF,[0.6 0.6 0.6],'EdgeColor','none','FaceAlpha',0.6)
end


%%%%%%%%%% plot the lever presses
for iTrial =1:size(tlims_,1) 
    try
    idx=[];
    try
        idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
    catch
        idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
        idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
    end
    idx = idx(1):idx(2);      
    
 
    plot([tlims_(iTrial,2);tlims_(iTrial,2)],[-1 -0.2]*SF,'color',[0 1 0 0.6],'LineWidth',2)
    plot([tlims_(iTrial,3);tlims_(iTrial,3)],[-1 -0.2]*SF,'color',[1 0 0 0.6],'LineWidth',2)
    end
end
% max_ = SF*size(Ass.FSCsel{iArea},2);

set(gca,'xtick',[min(tlims_)-bounds_:60:max(tlims_)+bounds_])
grid on

% plot(Ass.Tmtx(:),Ass.FSCsel{iArea}(:,iAss),'color','k','LineWidth',1.2)
% plot(Ass.Tmtx(j),Ass.FSCsel{iArea}(j,iAss)-10,'color',0.6*[1 0 0],'LineWidth',1.2)
%  clear bounds_ max_ iTrial iArea tlims_ SF iAss idx units_

%% FINAL - Plot all activations on continuous timesale with units overlaid - one example assembly showing non-member units 
SF = 1.2;
SF_Ass = 5;
bounds_ = 10;
Offset = 4;
plotNonMembers = false;
tExtreme = [5240 5490];
% tExtreme = [5000 8000];
shading_ = [0.3 0.3];
TrialSkip = 2;


for iArea = 3%1:3;
 tlims_ = [...
        t.Short.CueLight_LeftCorrect,  t.Short.SamplePress_LeftCorrect,    t.Short.ChoicePress_LeftCorrect;...
        t.Short.CueLight_RightCorrect, t.Short.SamplePress_RightCorrect,   t.Short.ChoicePress_RightCorrect;...
        t.Short.CueLight_LeftError',   t.Short.SamplePress_LeftError',     t.Short.ChoicePress_LeftError';...
        t.Short.CueLight_RightError',  t.Short.SamplePress_RightError',    t.Short.ChoicePress_RightError';...
        t.Medium.CueLight_LeftCorrect, t.Medium.SamplePress_LeftCorrect,   t.Medium.ChoicePress_LeftCorrect;...
        t.Medium.CueLight_RightCorrect,t.Medium.SamplePress_RightCorrect,  t.Medium.ChoicePress_RightCorrect;...
        t.Medium.CueLight_LeftError',  t.Medium.SamplePress_LeftError',    t.Medium.ChoicePress_LeftError';...
        t.Medium.CueLight_RightError', t.Medium.SamplePress_RightError',   t.Medium.ChoicePress_RightError';...
        t.Long.CueLight_LeftCorrect,   t.Long.SamplePress_LeftCorrect,     t.Long.ChoicePress_LeftCorrect;...
        t.Long.CueLight_RightCorrect,  t.Long.SamplePress_RightCorrect,    t.Long.ChoicePress_RightCorrect;...
        t.Long.CueLight_LeftError',    t.Long.SamplePress_LeftError',      t.Long.ChoicePress_LeftError';...
        t.Long.CueLight_RightError',   t.Long.SamplePress_RightError',     t.Long.ChoicePress_RightError']*1e-6;
           
        [tlims_ ,I]=sortrows(tlims_,1);  
        LR_ = t.LR(I');

LR_(min(tlims_,[],2)<tExtreme(1) | max(tlims_,[],2)>tExtreme(2)) = [];
tlims_(tlims_<tExtreme(1) | tlims_>tExtreme(2)) = NaN;
tlims_(sum(isnan(tlims_),2)>0,:)=[];



idx_ = closest(Tmtx,tExtreme);
Tmtx_ = Tmtx(idx_(1):idx_(2));  


AssCol = {[0.081 0.09 0.631],[0.761 0.505 0.046],[0.2,0.67,0.4]};
% [0.16,0.36,0.19]
% [0.2,0.67,0.4]
for iAss=1%:size(Ass.FSCsel{iArea},2)
    figure('color','w'); hold on

    %%%%%%%%%%%% plot the assemblies
    for iTrial =1:TrialSkip:size(tlims_,1)
        try
            idx=[];
            try
                idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
            catch
                idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
                idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
            end
            idx = idx(1):idx(2);
            temp_ = Ass.FSCsel{iArea}(:,iAss);
            temp_ = mat2gray(temp_);
            temp_ = temp_(idx);
            temp_([1,end]) =0;%min(temp_([1,end]));
            
%             patch(Ass.Tmtx(idx),Offset*temp_,p.col_{iArea},'EdgeColor','none','FaceAlpha',0.8)
            patch(Ass.Tmtx(idx),SF_Ass*temp_,AssCol{iArea},'EdgeColor',AssCol{iArea},'FaceAlpha',0.8)
            
            %     %%%%%%%%%%%% plot the lever presses
            %     plot([tlims_(iTrial,2);tlims_(iTrial,2)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[0 1 0 0.6],'LineWidth',2)
            %     plot([tlims_(iTrial,3);tlims_(iTrial,3)],SF*[iAss iAss]+[-0.5,0.5]*SF,'color',[1 0 0 0.6],'LineWidth',2)
        catch
            disp(sprintf('No timestamps found for trial %G',iTrial))
        end
    end
    
    
    %%%%%%%%%%%% Plot the units

        
        iFR_ = iFR{3};%(:,units_Real);
        iFR_ = zscore(iFR_);
        
        iFR_ = bsxfun(@minus,iFR_,iFR_(1,:));
        iFR_ = iFR_(idx_(1):idx_(2),:);
        
        units_ = Ass.units{iArea}{iAss};
        %
        switch iArea
            case 1
                orderedUnits = [setdiff(1:Ass.nu(1),units_),...
                                units_ ];
            case 2
                orderedUnits = [units_, ...
                                setdiff((Ass.nu(1)+1):Ass.nu(3),units_)];
            case 3
                orderedUnits = [setdiff(1:Ass.nu(1),units_),...
                                units_ , ...
                                setdiff((Ass.nu(1)+1):Ass.nu(3),units_)];
        end
        %units_Real = Ass.usel_out{3}(Ass.units{iArea}{iAss}+UnitOffset(iArea));
        units_Real = Ass.usel_out{3}(orderedUnits);
        for iUnit =1:length(orderedUnits)
            temp_ = mat2gray(iFR_(:,units_Real(iUnit)));
            temp_([1,end]) = 0;%min(temp_([1,end]));
             ShadingLevel = shading_(ismember(orderedUnits(iUnit),units_)+1);
             if ismember(orderedUnits(iUnit),units_)
            if units_Real(iUnit) <= Ass.usel_out{1}(Ass.nu(1))
                patch(Tmtx_,Offset+temp_+iUnit*SF,[0.081 0.09 0.631],'EdgeColor',[0.081 0.09 0.631],'FaceAlpha',ShadingLevel)
            else
                patch(Tmtx_,Offset+temp_+iUnit*SF,[0.761 0.505 0.046],'EdgeColor',[0.761 0.505 0.046],'FaceAlpha',ShadingLevel)
            end
             else
                patch(Tmtx_,Offset+temp_+iUnit*SF,[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6],'FaceAlpha',ShadingLevel)
             end
        end

%%%%%%%%%%%% plot the lever presses
for iTrial =1:TrialSkip:size(tlims_,1) 
    try
    idx=[];
    try
        idx = closest(Ass.Tmtx,tlims_(iTrial,[1 3])+[-10 10]);
    catch
        idx(1) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,1)-10);
        idx(2) = FindClosestIndex(Ass.Tmtx,tlims_(iTrial,3)+10);
    end
    idx = idx(1):idx(2);      
    
 
    plot([tlims_(iTrial,2);tlims_(iTrial,2)],[-1 -0.2]*SF,'color',[0 1 0 0.6],'LineWidth',2)
    plot([tlims_(iTrial,3);tlims_(iTrial,3)],[-1 -0.2]*SF,'color',[1 0 0 0.6],'LineWidth',2)
    switch LR_(iTrial)
        case 1
            text(mean(tlims_(iTrial,2:3)),-2,'L','HorizontalAlignment','center')
        case 2
            text(mean(tlims_(iTrial,2:3)),-2,'R','HorizontalAlignment','center')
    end
    
    end
end
% max_ = SF*size(Ass.FSCsel{iArea},2);

set(gca,'xtick',[min(tlims_)-bounds_:60:max(tlims_)+bounds_])
grid on
axis([tExtreme -2 Inf])
axis off
end
end

%% Get cutouts and averages/variances for Assems and FRs 

SF = 10;
bounds_ = 10;
tlimsSample = [-5 5];
tlimsChoice = [-5 5];
KOs = sum(abs(tlimsSample))/p.bw;
KOc = sum(abs(tlimsChoice))/p.bw;

tb_ = (1:(KOs + KOc ))*p.bw;

tlims_ = cell(1,2);

for iDelay = 1:length(p.Delays_)
    
    tlims_{1} = eval(sprintf('[t.%s.SamplePress_LeftCorrect, t.%s.ChoicePress_LeftCorrect]*1e-6;',p.Delays_{iDelay},p.Delays_{iDelay}));
    tlims_{2} = eval(sprintf('[t.%s.SamplePress_RightCorrect, t.%s.ChoicePress_RightCorrect]*1e-6;',p.Delays_{iDelay},p.Delays_{iDelay}));

    tlims_{3} = eval(sprintf('[t.%s.SamplePress_LeftError, t.%s.ChoicePress_LeftError]*1e-6;',p.Delays_{iDelay},p.Delays_{iDelay}));
    tlims_{4} = eval(sprintf('[t.%s.SamplePress_RightError, t.%s.ChoicePress_RightError]*1e-6;',p.Delays_{iDelay},p.Delays_{iDelay}));
    
    % cutout indices for Assemblies
    idx = cell(1,length(tlims_));
    
        for iType =1:length(tlims_)
            if ~isempty(tlims_{iType})
                try
                    idx{iType} = [closest(Ass.Tmtx,tlims_{iType}(:,1)+tlimsSample(1)),closest(Ass.Tmtx,tlims_{iType}(:,2)+tlimsChoice(1))];
                catch
                    for iTrial = 1:length(tlims_{iType})
                        idx{iType}(iTrial,:) = [FindClosestIndex(Ass.Tmtx,tlims_{iType}(iTrial,1)+tlimsSample(1)),FindClosestIndex(Ass.Tmtx,tlims_{iType}(iTrial,2)+tlimsChoice(1))];
                    end
                end
            end
        end
        
    
    % cutout indices for Units 
    idx2 = cell(1,length(tlims_));
    for iType =1:length(tlims_)
        if ~isempty(tlims_{iType})
            try
                idx2{iType} = [closest(Tmtx,tlims_{iType}(:,1)+tlimsSample(1)),closest(Tmtx,tlims_{iType}(:,2)+tlimsChoice(1))];
            catch
                for iTrial = 1:length(tlims_{iType})
                    idx2{iType}(iTrial,:) = [FindClosestIndex(Tmtx,tlims_{iType}(iTrial,1)+tlimsSample(1)),FindClosestIndex(Tmtx,tlims_{iType}(iTrial,2)+tlimsChoice(1))];
                end
            end
        end
    end
    
            
    
    for iType =1:length(tlims_)
        if ~isempty(tlims_{iType})
            for iTrial =1:size(tlims_{iType},1)
                % Assemblies
                idx_ = [idx{iType}(iTrial,1):(idx{iType}(iTrial,1)+KOs-1) , ...
                    idx{iType}(iTrial,2):(idx{iType}(iTrial,2)+KOs-1)];
                for iArea = 1:3
                    for iAss = 1:size(Ass.FSCsel{iArea},2)
                        Ass.Cutouts{iArea}{iDelay}{iAss}{iType}(:,iTrial)=Ass.FSCsel{iArea}(idx_,iAss);
                    end
                end
                % Units
                idx2_ = [idx2{iType}(iTrial,1):(idx2{iType}(iTrial,1)+KOs-1) , ...
                    idx2{iType}(iTrial,2):(idx2{iType}(iTrial,2)+KOs-1)];
                for iArea = 1:3
                    for iUnit = 1:size(iFR{iArea},2)
                        Units.Cutouts{iArea}{iDelay}{iUnit}{iType}(:,iTrial)=iFR{iArea}(idx_,iUnit);
                    end
                end
                
            end
            
            for iArea = 1:3
                for iAss = 1:size(Ass.FSCsel{iArea},2)
                    Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,iType) = nanmean(Ass.Cutouts{iArea}{iDelay}{iAss}{iType},2);
                    Ass.CutoutsSEM{iArea}{iDelay}{iAss}(:,iType)  = nansem(Ass.Cutouts{iArea}{iDelay}{iAss}{iType},2);
                end
                for iUnit = 1:size(iFR{iArea},2)
                    Units.CutoutsMean{iArea}{iDelay}{iUnit}(:,iType) = nanmean(Units.Cutouts{iArea}{iDelay}{iUnit}{iType},2);
                    Units.CutoutsSEM{iArea}{iDelay}{iUnit}(:,iType)  = nansem(Units.Cutouts{iArea}{iDelay}{iUnit}{iType},2);
                end
            end
        end
    end
end
%% Plot all Assemblies staggered for L/R conditions on average
SF = 10  ;  
for iArea = 1:3
figure('name',p.Areas{iArea}); hold on

for iDelay = 1:length(p.Delays_)
    
    for iAss = 1:size(Ass.FSCsel{iArea},2)
        ciplot(iAss*SF+[Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,1)+Ass.CutoutsSEM{iArea}{iDelay}{iAss}(:,1)],...
            iAss*SF+[Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,1)-Ass.CutoutsSEM{iArea}{iDelay}{iAss}(:,1)],...
            tb_+(max(tb_)+2)*(iDelay-1),'r');
        ciplot(iAss*SF+[Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,2)+Ass.CutoutsSEM{iArea}{iDelay}{iAss}(:,2)],...
            iAss*SF+[Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,2)-Ass.CutoutsSEM{iArea}{iDelay}{iAss}(:,2)],...
            tb_+(max(tb_)+2)*(iDelay-1),'g');
    end
end
end
%% Choose a condition/Assembly and plot
iDelay = 3;
iArea = 2; 
iAss = 3;

col_ = fliplr({[0.5 0 0.5],[0 0 1]});
figure('color','w','name',p.Areas{iArea}); hold on
plot(tb_,Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,3),'color',col_{1},'LineStyle',':','LineWidth',1.2);    
plot(tb_,Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,4),'color',col_{2},'LineStyle',':','LineWidth',1.2);        

ciplot([Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,1)+Ass.CutoutsSEM{iArea}{iDelay}{iAss}(:,1)],...
            [Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,1)-Ass.CutoutsSEM{iArea}{iDelay}{iAss}(:,1)],...
            tb_,col_{1},0.8);
ciplot([Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,2)+Ass.CutoutsSEM{iArea}{iDelay}{iAss}(:,2)],...
            [Ass.CutoutsMean{iArea}{iDelay}{iAss}(:,2)-Ass.CutoutsSEM{iArea}{iDelay}{iAss}(:,2)],...
            tb_,col_{2},0.8);
        
        L  = Ass.Cutouts{iArea}{iDelay}{iAss}{1};
        R  = Ass.Cutouts{iArea}{iDelay}{iAss}{2};
        [p_,~] = permtest2vec(L,R,1000,0.05);
        
a = nan(size(p_));a(p_)=1;
% a(x==y)=NaN;
plot(tb_,15*a,'k','LineWidth',4)

        plot(abs(tlimsSample(1))*[1 1],[-10 -5],'g','LineWidth',2)
        plot((sum(abs(tlimsSample))+abs(tlimsSample(1)))*[1 1],[-10 -5],'r','LineWidth',2)
                plot(sum(abs(tlimsSample))*[1 1],[-5 5],':k','LineWidth',1.5)
axis([0 20 -10 15])
      axis off
%% Assembly average with individual member units overlaid
for iArea = 1:3
    % Sort Assemblies by size
    
    I = cellfun(@length,Ass.units{iArea});
    [~,I] = sort(I,'descend');
    units_ = Ass.units{iArea}(I);
    
    
    figure('name',p.Areas{iArea}); hold on
    
    for iDelay = 3%1:length(p.Delays_)
        
        for iUnit = 1:size(iFR{iArea},2)
            ciplot(iUnit*SF+[Units.CutoutsMean{iArea}{iDelay}{iUnit}(:,1)+Units.CutoutsSEM{iArea}{iDelay}{iUnit}(:,1)],...
                iUnit*SF+[Units.CutoutsMean{iArea}{iDelay}{iUnit}(:,1)-Units.CutoutsSEM{iArea}{iDelay}{iUnit}(:,1)],...
                tb_+(max(tb_)+2)*(iDelay-1),'r');
            ciplot(iUnit*SF+[Units.CutoutsMean{iArea}{iDelay}{iUnit}(:,2)+Units.CutoutsSEM{iArea}{iDelay}{iUnit}(:,2)],...
                iUnit*SF+[Units.CutoutsMean{iArea}{iDelay}{iUnit}(:,2)-Units.CutoutsSEM{iArea}{iDelay}{iUnit}(:,2)],...
                tb_+(max(tb_)+2)*(iDelay-1),'g');
        end
    end
    
end


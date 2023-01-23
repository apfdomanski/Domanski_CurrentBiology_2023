%% %%%%%% PREAMBLE %%%%%%
clear
warning ('off')
Target = 'LONG';
bw=0.05;
minFR = 0.5;
clear Av

if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
elseif ismac
    pat = '/Volumes/HDD2/DNMTP/raw/';
    home = getenv('HOME');
else
    
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    addpath(genpath('/panfs/panasas01/phph/ad15419/MATLAB'))
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir(fullfile(pat,'KDE_binsTaskonly',sprintf('*%s*PFC*.mat',Target)));

reject_list={'MiroslawLONG1_PFC_iFR50_behavOnly.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

Areas = {'PFC','HP','Joint'};
Areas_ = {'dCA1','mPFC','dCA1-mPFC'};
color_={'b','r','g'};

AssemblyChoice = 2;
% 1 - lever press cutouts
% 2 - trial cutouts (cue to reward inclusive)
% 3 - continuous time
switch AssemblyChoice
    case 1
        pat2 = [pat 'KDE_bins' filesep Target filesep];
    case 2
        pat2 = [pat 'KDE_binsTaskonly' filesep 'LONGTaskonly' filesep];
    case 3
        pat2 = 'C:\Analysis\AssemblyAnalysis\Sleep\Task\';
end

MemberClasses ={'NonMembers','LocalMembers','JointMembers'};
MemberClasses_ ={'Non-members','Local Assembly Members','Inter-area Assembly Members'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
%% Batch process units
for iFile = 1:length(fileList)
   
    %% Get the assembly files
    fname=strtok(fileList(iFile).name,'_');
    fprintf('Analysing run %d/%d %s...\n',iFile,length(fileList),fname)

    load(sprintf('%s%s%s.mat',pat,filesep,fname),'PFCcells','HPcells');
    nu = [length(PFCcells),length(HPcells)];nu(3)= sum(nu); clear PFCcells HPcells
    %% Get the assembly membership info
           
            fprintf('Loading Assemblies %d/%d %s ...\n',iFile,length(fileList),fname)
            % for assemblies describing whole task period (cue  - reward)
            switch AssemblyChoice
                case 1
                    load(sprintf('%s%s%s_iFR50_FSC.mat',pat2,filesep,fname),'units','nu');
                    %             load(sprintf('%s%s%s_iFR50_AssemRes2.mat',pat2,filesep,fname),'usel_out');
                    usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
                case 2
                    load(sprintf('%s%s_iFR50_BehavOnly_FSC.mat',pat2,fname),'units','nu','FSCsel','Tmtx');
                    load(sprintf('%s%s_iFR50_BehavOnly_AssemRes2.mat',pat2,fname),'usel_out');
                case 3
                    load(sprintf('%s%s%s_iFR50_Task_FSC.mat',pat2,filesep,fname),'units','nu');
                    %             load(sprintf('%s%s%s_iFR50_Task_AssemRes2.mat',pat2,filesep,fname),'usel_out');
                    usel_out=SelCells(sprintf('%sKDE_binsTaskonly%s%s_PFC_iFR50_behavOnly.mat',pat,filesep, fname));
            end
            Ass.usel_out    = usel_out;
            Ass.FSC         = FSCsel;
            Ass.Tmtx        = Tmtx;
            Ass.usel_out{3} = [Ass.usel_out{1},Ass.usel_out{2}+nu(1)];
            Ass.units       = units;
            Ass.nu = nu; Ass.nu(3) = sum(Ass.nu(1:2));
            % Check that inter-area assemblies actually span the two areas
            if ~isempty(Ass.FSC{3})
                nAss = length(Ass.units{3});
                idx = false(nAss,1);
                for iAss = 1:nAss
                    U_ = Ass.usel_out{3}(Ass.units{3}{iAss});
                    
                    if  min(U_)>max(Ass.usel_out{1})
                        idx(iAss) = true;
                    end
                    Ass.units{3}(idx)=[];
                    Ass.FSC{3}(:,idx)=[];
                    
                end
                
            end
            for iArea_ = 1:3
                Ass.LocalMembers{iArea_}    = unique(cell2mat(Ass.units{iArea_}));
            end
            
            
            Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3}<=Ass.nu(1)));
            Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3}>Ass.nu(1))-Ass.nu(1));
            for iArea_ = 1:2
                Ass.JointMembers{iArea_} = setdiff(unique(cell2mat(Ass.units{iArea_})),Ass.LocalMembers{iArea_});
                Ass.NonMembers{iArea_}   = setdiff(1:Ass.nu(iArea_),[Ass.LocalMembers{iArea_},Ass.JointMembers{iArea_}]);
            end
            for iArea_ = 1:2
                Ass.NonMembers{iArea_}   = Ass.usel_out{iArea_}(Ass.NonMembers{iArea_});
                Ass.LocalMembers{iArea_} = Ass.usel_out{iArea_}(Ass.LocalMembers{iArea_});
                Ass.JointMembers{iArea_} = Ass.usel_out{iArea_}(Ass.JointMembers{iArea_});
            end
            
            for iArea_ = 1:2
                %Membership_{iArea_} = -ones( Ass.nu(iArea_),1);
                Membership_{iArea_} = -ones(nu(iArea_),1);
                Membership_{iArea_}(Ass.NonMembers{iArea_})=0;
                Membership_{iArea_}(Ass.LocalMembers{iArea_})=1;
                Membership_{iArea_}(Ass.JointMembers{iArea_})=2;
                Membership_{iArea_}(setdiff(1:nu(iArea_),Ass.usel_out{iArea_}))=[];
            end

            for iArea_ = 1:2
                % Pad this out to ensure that there's an entry even if there's no detected assembly
                LocalMembersMatrix_{iArea_} = nan(Ass.nu(iArea_),length(Ass.units{iArea_})+1);
                if ~isempty(Ass.units{iArea_})
                    for iAss=1:length(Ass.units{iArea_})
                        LocalMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_),Ass.units{iArea_}{iAss});
                    end
                end
                
                JointMembersMatrix_{iArea_} = nan(Ass.nu(iArea_),length(Ass.units{3})+1);
                for iAss=1:length(Ass.units{3})
                    if ~isempty(Ass.units{3})
                        units_ = Ass.units{3}{iAss};
                        if iArea_==1
                            JointMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_), units_(units_<=Ass.nu(1)));
                        else
                            JointMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_), units_(units_>Ass.nu(1))-Ass.nu(1));
                        end
                    end
                end
            end
            iArea_ = 3;
            JointMembersMatrix_{iArea_} = nan(Ass.nu(iArea_),length(Ass.units{3})+1);
            for iAss=1:length(Ass.units{3})
                if ~isempty(Ass.units{3})
                    units_ = Ass.units{3}{iAss};
                    JointMembersMatrix_{iArea_}(:,iAss+1)=ismember(1:Ass.nu(iArea_), units_);
                end
            end
            clear units usel_out nu FSCsel Tmtx  iFR noPFC idx nAss U_ iAss
            nu = cellfun(@length,Membership_);nu(3) = sum(nu);
            
    
    nonMemberUnits = cell(3,1);
    for iArea =1:2
        nonMemberUnits{iArea} = find(Membership_{iArea}==0);
    end
    nonMemberUnits{2} = nonMemberUnits{2}+nu(1);
    nonMemberUnits{3} = find([Membership_{1};Membership_{2}]==0);
     
    % Remove local assemblies that are subsets of joint assemblies
    for iArea=1:2
        toremove = zeros(1,size(LocalMembersMatrix_{iArea},2));
        for iAss = 2:size(LocalMembersMatrix_{iArea},2)
            members_=find(LocalMembersMatrix_{iArea}(:,iAss)>0);
            toremove_= zeros(1,size(JointMembersMatrix_{iArea},2));
            for jAss=2:size(JointMembersMatrix_{iArea},2)
                toremove_(jAss) = sum( sum(ismember(members_,find(JointMembersMatrix_{iArea}(:,jAss)))) == length(members_));
            end
            toremove(iAss) = sum(toremove_)>0;
        end
       LocalMembersMatrix_{iArea}(:,toremove==1) =[];
       
       Membership_{iArea} = zeros(size(LocalMembersMatrix_{iArea},1),1);
       Membership2_{iArea} = zeros(size(LocalMembersMatrix_{iArea},1),1);
       
       % 0=nonmember, 1=localmember, 2 = jointmember | local&jointmember
       Membership_{iArea} = max([ nansum(LocalMembersMatrix_{iArea},2),...
                                  2*(nansum(JointMembersMatrix_{iArea},2)>0)...
                                   ],[],2);
       
       % 0=nonmember, 1=localmember, 2 = jointmember, 3=local&jointmember
       Membership2_{iArea} = max([1*(nansum(LocalMembersMatrix_{iArea},2)==1 & nansum(JointMembersMatrix_{iArea},2)==0),... 
                             2*(nansum(LocalMembersMatrix_{iArea},2)<1 & nansum(JointMembersMatrix_{iArea},2)>0),... 
                             3*(nansum(LocalMembersMatrix_{iArea},2)>0 & nansum(JointMembersMatrix_{iArea},2)>0)],[],2); 
    end
    
    
    clear iAss jAss toremove toremove_ 
    %% Get the firing rate info
    for iArea=1:2
        fn = fullfile(pat,'KDE_binsTaskonly',sprintf('%s_%s_iFR50_behavOnly.mat',fname,Areas{iArea}));
        load(fn,'avgFR','hucv','MISE','CvBW');
        D.KDEsigma{iArea}{iFile} = hucv(Ass.usel_out{iArea}');
        D.avgFR{iArea}{iFile}    = avgFR(Ass.usel_out{iArea}');
        D.MISE{iArea}{iFile}     = MISE(Ass.usel_out{iArea}');
        D.CvBW{iArea}{iFile}     = CvBW(Ass.usel_out{iArea}');
    end
    %% Assembly statistics
    D.nAss{iFile} = [size(LocalMembersMatrix_{1},2),...
                     size(LocalMembersMatrix_{2},2),...
                     size(JointMembersMatrix_{1},2)]-1;
                 
    D.AssSizeLocal{iFile} = cellfun(@(x) nansum(x(:,2:end)),LocalMembersMatrix_,'UniformOutput',false);
    D.AssSizeJoint{iFile} = cellfun(@(x) nansum(x(:,2:end)),JointMembersMatrix_,'UniformOutput',false);
    for iArea = 1:2
        
        D.FracCellsEachClass{iArea}{iFile} = [sum(Membership_{iArea}==0),...
                                     sum(Membership_{iArea}==1),...
                                     sum(Membership_{iArea}==2),...
                                     sum(Membership2_{iArea}==2),...
                                     sum(Membership2_{iArea}==3)]./length(Membership_{iArea});
%        D.MembershipClass{iArea}{iFile} = Membership2_{iArea};
       D.MembershipClass{iArea}{iFile} = Membership_{iArea};
       D.NodeDegree{iArea}{iFile} = nansum([LocalMembersMatrix_{iArea},JointMembersMatrix_{iArea}],2);
       
    end
    
    % Balance between dCA1 and mPFC in joint assems
    for iAss=2:size(JointMembersMatrix_{3},2)
        members_ = find(JointMembersMatrix_{3}(:,iAss));
        D.JointAssemFracCA1{iFile}(iAss-1) = sum(members_>Ass.nu(1))./length(members_);
    end
    
    % equivalence of dCA1 and mPFC unit counts
    D.UnitCountFracCA1(iFile) = (Ass.nu(1)-Ass.nu(2))./sum(Ass.nu(1:2));
%     D.UnitCountFracCA1(iFile) = Ass.nu(2)/Ass.nu(1);
end
clearvars -except D Areas Areas_ MemberClasses MemberClasses_ color_ col_ fileList 
%% Meta-analysis

%% 0) FR and KDE dCA1 vs mPFC
clear avFR KDEsigma MISE CvBW FR_ KDE_ MISE_ CvBW_
avFR{1} = cell2mat(D.avgFR{1}); avFR{2} = cell2mat(D.avgFR{2});
KDEsigma{1} = cell2mat(D.KDEsigma{1}); KDEsigma{2} = cell2mat(D.KDEsigma{2});
MISE{1} = cell2mat(D.MISE{1}); MISE{2} = cell2mat(D.MISE{2});
CvBW{1} = cell2mat(D.CvBW{1}); CvBW{2} = cell2mat(D.CvBW{2});
for iArea = 1:2
    FR_  = avFR{iArea};
    KDE_ = KDEsigma{iArea};
    MISE_ = MISE{iArea};
    CvBW_ = CvBW{iArea};

    idx = FR_<0.1;
    
    FR_(idx)=[];
    KDE_(idx)=[];
    MISE_(idx)=[];
    CvBW_(idx)=[];
    
    avFR{iArea}=FR_;
    KDEsigma{iArea}=KDE_.^0.5;
    MISE{iArea} = MISE_;
    CvBW{iArea} = CvBW_;
end
c_ = {[0 0 1 0.2],[1 0 0 0.2]};

figure; 
subplot(1,2,1);hold on
bins = 0:0.25:10;
for iArea=1:2
    for iFile=1:length(fileList)
%         D.KDEsigma{iArea}{iFile};
        y1 = histc(D.avgFR{iArea}{iFile},bins); y1=y1./sum(y1);
        plot(bins,cumsum(y1),'color',c_{iArea},'LineWidth',1,'HandleVisibility','off')
    end
end
y1 = histc(avFR{1},bins); y1=y1./sum(y1);
y2 = histc(avFR{2},bins); y2=y2./sum(y2);
stairs(bins,cumsum(y1),'b','LineWidth',1.5)
stairs(bins,cumsum(y2),'r','LineWidth',1.5)
[~,p(1),stat(1)] = kstest2(avFR{1},avFR{2})
% p_(1) = ranksum(avFR{1},avFR{2})
xlabel('Mean Firing rate (Hz)')
ylabel('Fraction of units')
axis([min(bins) max(bins) 0 1])
legend('mPFC','dCA1','Location','SouthEast'); legend boxoff

subplot(1,2,2);hold on
bins = 0:0.25:10;
for iArea=1:2
    for iFile=1:length(fileList)
%         D.KDEsigma{iArea}{iFile};
        y1 = histc(D.KDEsigma{iArea}{iFile},bins); y1=y1./sum(y1);
        plot(bins,cumsum(y1),'color',c_{iArea},'LineWidth',1,'HandleVisibility','off')
    end
end
y1 = histc(KDEsigma{1},bins); y1=y1./sum(y1);
y2 = histc(KDEsigma{2},bins); y2=y2./sum(y2);
stairs(bins,cumsum(y1),'b','LineWidth',1.5)
stairs(bins,cumsum(y2),'r','LineWidth',1.5)
[p(2),~,stat(2)] =kstest2(KDEsigma{1},KDEsigma{2})
% p_(2) = ranksum(KDEsigma{1},KDEsigma{2})
xlabel('Rate smoothing kernel sigma (ms)')
axis([min(bins) max(bins) 0 1])
%% 0.1) Fit resduals and KDE estimation stability 
figure; 
subplot(1,2,1);hold on
for iArea=1:2
    scatter3(avFR{iArea},KDEsigma{iArea},-MISE{iArea},20,color_{iArea},'Filled','MarkerFaceAlpha',0.3)
end
box on; grid on
view(50,45)
subplot(1,2,2);hold on
for iArea=1:2
     scatter3(avFR{iArea},KDEsigma{iArea},CvBW{iArea},20,color_{iArea},'Filled','MarkerFaceAlpha',0.3)
end
axis([0 10 0 10 0 3])
box on; grid on
view(50,45)
%% 0.2) FR vs MISE and CvBW
figure;
subplot(1,3,1); hold on
for iArea=1:2
    scatter(avFR{iArea},MISE{iArea},20,color_{iArea},'Filled','MarkerFaceAlpha',0.3)
end
axis([0 10 -Inf Inf])
xlabel('Mean Firing rate (Hz)')
ylabel('KDE estimation residual')

[~,p(1) stat(2)] = kstest2(MISE{1},MISE{2})
% p_(1) = ranksum(MISE{1},MISE{2})

subplot(1,3,2); hold on
for iArea=1:2
	scatter(avFR{iArea},CvBW{iArea},20,color_{iArea},'Filled','MarkerFaceAlpha',0.3)
end
axis([0 10 -Inf Inf])
xlabel('Mean Firing rate (Hz)')
ylabel('KDE estimation variance over time')
subplot(1,3,3); hold on
for iArea=1:2
	scatter(MISE{iArea},CvBW{iArea},20,color_{iArea},'Filled','MarkerFaceAlpha',0.3)
end

[~,p(1), stat(2)] = kstest2(CvBW{1},CvBW{2})
% p_(1) = ranksum(CvBW{1},CvBW{2})
%% 1) no. Ass in each recording
y = cell2mat(D.nAss');
ym = nanmean(y(:,[2,1,3]));
ys = nansem(y(:,[2,1,3]));
figure; hold on
b(1) = bar(1,ym(1),'EdgeColor','r','FaceColor','w','LineWidth',1.5);
    b(2) = bar(2,ym(2),'EdgeColor','b','FaceColor','w','LineWidth',1.5);
    b(3) = bar(3,ym(3),'EdgeColor','g','FaceColor','w','LineWidth',1.5);
    
    
    errorbar(1,ym(1),ys(1),'LineStyle','none','color','r','LineWidth',1.5)
    errorbar(2,ym(2),ys(2),'LineStyle','none','color','b','LineWidth',1.5)
    errorbar(3,ym(3),ys(3),'LineStyle','none','color','g','LineWidth',1.5)
    scatter(1+0.05*randn(size(y,1),1),y(:,1),40,'r','filled','MarkerFaceAlpha',0.3)
    scatter(2+0.05*randn(size(y,1),1),y(:,2),40,'b','filled','MarkerFaceAlpha',0.3)
    scatter(3+0.05*randn(size(y,1),1),y(:,3),40,'g','filled','MarkerFaceAlpha',0.3)
    
% bar(1:3,ym,'EdgeColor','k','LineWidth',1.5,'FaceColor','w')
% errorbar(1:3,ym,ys,'LineStyle','none','color','k','LineWidth',1.5)
set(gca,'XTick',1:3,'YTick',0:5:15,'XTickLabel',Areas_,'XTickLabelRotation',45)
ylabel('No. Assemblies per session')
axis([0 4 0 15])
% [p,tbl,stats] = kruskalwallis(y);
% multcompare(stats)

%% 2) Fraction of units in each type of assem
% figure;
% for iArea=1:2
%     subplot(1,2,iArea); hold on
%     y = cell2mat(D.FracCellsEachClass{iArea}');
%     ym = nanmean(y);
%     ys = nansem(y);
%     
%     bar(1:3,ym(1:3),'EdgeColor','k','LineWidth',1.5,'FaceColor','w')
%     errorbar(1:3,ym(1:3),ys(1:3),'LineStyle','none','color','k','LineWidth',1.5)
%     set(gca,'XTick',1:3,'YTick',1:6,'XTickLabel',MemberClasses_,'XTickLabelRotation',45)
%     ylabel('Fraction of units in each class')
%     axis([0 4 0 1])
% end

figure; hold on
y1 = cell2mat(D.FracCellsEachClass{1}');
y2 = cell2mat(D.FracCellsEachClass{2}');
    ymPFC = nanmean(y1);
    ysPFC = nansem(y1);
    ymCA1 = nanmean(y2);
    ysCA1 = nansem(y2);
    ym = [ymCA1(1),ymPFC(1),ymCA1(2),ymPFC(2),ymCA1(3),ymPFC(3)];
    ys = [ysCA1(1),ysPFC(1),ysCA1(2),ysPFC(2),ysCA1(3),ysPFC(3)];

    ym = [ymCA1(1),ymPFC(1);ymCA1(2),ymPFC(2);ymCA1(3),ymPFC(3)];
    b = bar(ym,'FaceColor','w','LineWidth',1.5);
    b(1).EdgeColor=color_{2};   
    b(2).EdgeColor=color_{1};   
    set(gca,'XTick',1:3,'YTick',0:0.2:1,'XTickLabel',MemberClasses_,'XTickLabelRotation',45)
    errorbar((1:3)-0.15,ymCA1(1:3),ysCA1(1:3),'LineStyle','none','color',color_{2},'LineWidth',1.5)
    errorbar((1:3)+0.15,ymPFC(1:3),ysPFC(1:3),'LineStyle','none','color',color_{1},'LineWidth',1.5)
    legend('dCA1','mPFC'); legend boxoff
    axis([0.5 3.5 0 1])
    ylabel('Frac. neurons in Assemblies')
[p,tbl,stats] = anova1([y1(:,1),y2(:,1),y1(:,2),y2(:,2),y1(:,3),y2(:,3)]);
multcompare(stats)
%% 3) Size of assemblies

AssSize=cell(1,3);
for iFile =1:length(fileList)
    for iArea=1:2
        AssSize{iArea}{iFile,1} =D.AssSizeLocal{iFile}{iArea};
    end
    AssSize{3}{iFile,1} =D.AssSizeJoint{iFile}{3};
end
for iArea=1:3
    a = cellfun(@transpose,AssSize{iArea}(~isempty_cell(AssSize{iArea})),'UniformOutput',false)
    AssSize_{iArea}= cell2mat(a);
    rnge_{iArea}(1:2) = [min(AssSize_{iArea}),max(AssSize_{iArea})];
    mode_(iArea) = mode(AssSize_{iArea});
end


y = [cellfun(@mean,AssSize{2}),cellfun(@mean,AssSize{1}),cellfun(@mean,AssSize{3})];
ym = nanmean(y);
ys = nansem(y);

figure; hold on

    b(1) = bar(1,ym(1),'EdgeColor','r','FaceColor','w','LineWidth',1.5);
    b(2) = bar(2,ym(2),'EdgeColor','b','FaceColor','w','LineWidth',1.5);
    b(3) = bar(3,ym(3),'EdgeColor','g','FaceColor','w','LineWidth',1.5);
    
    
    errorbar(1,ym(1),ys(1),'LineStyle','none','color','r','LineWidth',1.5)
    errorbar(2,ym(2),ys(2),'LineStyle','none','color','b','LineWidth',1.5)
    errorbar(3,ym(3),ys(3),'LineStyle','none','color','g','LineWidth',1.5)
    scatter(ones(size(AssSize_{1}))+0.05*randn(size(AssSize_{1})),AssSize_{1},50,'r','filled','MarkerFaceAlpha',0.3)
    scatter(1+ones(size(AssSize_{2}))+0.05*randn(size(AssSize_{2})),AssSize_{2},50,'b','filled','MarkerFaceAlpha',0.3)
    scatter(2+ones(size(AssSize_{3}))+0.05*randn(size(AssSize_{3})),AssSize_{3},50,'g','filled','MarkerFaceAlpha',0.3)
%     legend('dCA1','mPFC'); legend boxoff
    axis([0 4 0 15])
    set(gca,'XTick',1:3,'YTick',0:5:15,'XTickLabel',Areas_,'XTickLabelRotation',45)
    
[p,tbl,stats] = anova1(y);
multcompare(stats)
%% 4) Node Degree
figure; hold on
x = 0:5;
y = cellfun(@(x) cell2mat(x'), D.NodeDegree,'UniformOutput',false');
for iArea=1:2
%     y_{iArea} = y{iArea}(y{iArea}>0);
	y_{iArea} = y{iArea};
    yHist = (histc(y_{iArea},x));
%     yHist = yHist./max(yHist);
    yHist = yHist./sum(yHist);
    stairs(x,yHist,'color',color_{iArea},'LineWidth',1.5)
    
end
set(gca,'XTick',0.5:1:4.5,'XTicklabel',[0:4],'YTick',[0:0.3:0.6])
for iArea=1:2
    y_{iArea}(y_{iArea}<1)=[];
    FracOverlap(iArea) = sum(y_{iArea}==1)./length(y_{iArea});
end
% cellfun(@nanmean,y)
% cellfun(@nansem,y)
% [h,p,stat] = kstest2(y{1},y{2})
%% Fraction of units in assemblies of any type
y = cellfun(@(x) cell2mat(x'), D.NodeDegree,'UniformOutput',false');

n1 = sum(y{1}==0); N1 = length(y{1});
n2 = sum(y{2}==0); N2 = length(y{2});
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2)
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
[h,p,stats] = chi2gof([1 2 3 4],'freq',observed,'expected',expected,'ctrs',[1 2 3 4],'nparams',2)
%% Node degree vs.firing rate
clear mdl
y = cellfun(@(x) cell2mat(x'), D.NodeDegree,'UniformOutput',false');
figure
for iArea = 1:2
    clear ym ys yVals xVals y_
    subplot(1,2,iArea); hold on

    xVals = y{iArea};
    yVals = avFR{iArea};
    idx = isnan(xVals);
    xVals(idx)=[];yVals(idx)=[];
    
    mdl{iArea} = fitlm(xVals,yVals);
    
    [Rcorr{iArea},pcorr{iArea}] = corrcoef(xVals,yVals);
    Rcorr{iArea}=Rcorr{iArea}.^2;
    
    h = plot(mdl{iArea});
    h(1).Marker = 'o';
    h(1).MarkerEdgeColor = 0.6*[1 1 1];
    h(2).LineWidth=1.5;
    h(3).LineWidth=1.5;
    h(4).LineWidth=1.5;
    for i=1:4
        ym(i) = nanmean(yVals(xVals==i-1));
        ys(i) = nansem(yVals(xVals==i-1));
    end
    errorbar(0:3,ym,ys,'color','k','LineWidth',2,'LineStyle','-')
%     set(gca,'YScale','log')
    axis([-0.5 3.5 0 10.5])
    set(gca,'XTick',0:3,'YTick',0:5:10)
    
    
    p(iArea) = ranksum(yVals(xVals==0),yVals(xVals>0));
    
    legend off
    xlabel(''); ylabel('')
    title('')
    
end
%% Node degree vs. KDE kernels
clear mdl R p 
figure
for iArea = 1:2
    clear ym ys
    subplot(1,2,iArea); hold on
    
    xVals = y{iArea};
    yVals = KDEsigma{iArea}.^0.5;
    idx = isnan(xVals);
    xVals(idx)=[];yVals(idx)=[];
    mdl{iArea} = fitlm(xVals,yVals);
    [Rcorr{iArea},pcorr{iArea}] = corrcoef(xVals,yVals);
    Rcorr{iArea}=Rcorr{iArea}.^2;
    h = plot(mdl{iArea});
    h(1).Marker = 'o';
    h(1).MarkerEdgeColor = 0.6*[1 1 1];
    h(2).LineWidth=1.5;
    h(3).LineWidth=1.5;
    h(4).LineWidth=1.5;
    for i=1:4
        ym(i) = nanmean(yVals(xVals==i-1));
        ys(i) = nansem(yVals(xVals==i-1));
    end
    errorbar(0:3,ym,ys,'color','k','LineWidth',2,'LineStyle',':')
%     set(gca,'YScale','log')
    axis([-0.5 3.5 0 10])
    set(gca,'XTick',0:3,'YTick',0:5:20)
        p(iArea) = ranksum(yVals(xVals==0),yVals(xVals>0));

    
    legend off
    xlabel(''); ylabel('')
    title('')
    
end
%% Unit balance in joint assems
% balance of recorded units between dCA1 and mPFC
bins = -1:0.2:1;
% bins = 0:0.1:2;
x = D.UnitCountFracCA1;
y_ = histc(x,bins);
figure
subplot(1,2,1); hold on
stairs(bins,y_,'k','LineWidth',1.5)
% errorbar_x(nanmean(x),0,nansem(x),nansem(x),'b')
% scatter(nanmean(x),0,'b','filled')
errorbar_x(median(x),0,iqr(x),iqr(x),'b')
scatter(median(x),0,'b','filled')
plot([0 0],[0 6],':k')
xlabel({'Balance of [dCA1,mPFC] unit counts';'in inter-area recordings'})
ylabel('No. Sessions')
axis([min(bins) max(bins) -0.1 6])
[~,ks_p] = kstest(x);
[~,lillie_p] = lillietest(x);
[~,jb_p] = jbtest(x);
[~,ad_p] = adtest(x);
[~,sw_p] = swtest(x);


bins = 0:0.05:1;
x = cell2mat(D.JointAssemFracCA1(~isempty_cell(D.JointAssemFracCA1)));
y_ = histc(x,bins);
subplot(1,2,2); hold on
stairs(bins,y_,'k','LineWidth',1.5)
% errorbar_x(nanmean(x),0,nansem(x),nansem(x),'b')
% scatter(nanmean(x),0,'b','filled')
errorbar_x(median(x),0,iqr(x),iqr(x),'b')
scatter(median(x),0,'b','filled')
plot([0.5 0.5],[0 18],':k')
xlabel({'dCA1/mPFC Balance of units','membership in inter-area assemblies'})
ylabel('No. Assemblies')
axis([min(bins) max(bins) -0.2 15])
[~,ks_p] = kstest(x);
[~,lillie_p] = lillietest(x);
[~,jb_p] = jbtest(x);
[~,ad_p] = adtest(x);
[~,sw_p] = swtest(x);


%% Firing rate vs different member types  
clear p stats tbl yHist
% D.FracCellsEachClass{iArea}{iFile} = [sum(Membership_{iArea}==0),...
%                                      sum(Membership_{iArea}==1),...
%                                      sum(Membership_{iArea}==2),...
%                                      sum(Membership2_{iArea}==2),...
%                                      sum(Membership2_{iArea}==3)]./length(Membership_{iArea});
% 0=nonmember, 1=localmember, 2 = jointmember, 3=local&jointmember
bins=0:0.05:10.5;

y = cellfun(@(x) cell2mat(x'), D.MembershipClass,'UniformOutput',false');
figure
for iArea = 1:2
    clear ym ys yVals xVals y_
    subplot(1,2,iArea); hold on

    xVals = y{iArea};
    yVals = avFR{iArea};
%     yVals = KDEsigma{iArea}.^0.5;

    idx = isnan(xVals);
    xVals(idx)=[];yVals(idx)=[];
    
%     scatter
    for i=1:3
        ym(i) = nanmean(yVals(xVals==i-1));
        ys(i) = nansem(yVals(xVals==i-1));
        yHist{iArea}(i,:) = histc(yVals(xVals==i-1),bins);
        yHist{iArea}(i,:)=yHist{iArea}(i,:)./sum(yHist{iArea}(i,:));
%         ym(i) = nanmedian(yVals(xVals==i-1));
%         ys(i) = iqr(yVals(xVals==i-1));
        bar(i-1,ym(i),'FaceColor','none','EdgeColor',col_{i},'LineWidth',2)
        scatter(i-1+0.05*randn(size(yVals(xVals==i-1))),yVals(xVals==i-1),50,col_{i},'filled','MarkerFaceAlpha',0.3)
        errorbar(i-1,ym(i),ys(i),'color',col_{i},'LineWidth',2,'LineStyle','-')
        
    end
axis([-0.5 2.5 0 10.5])
    set(gca,'XTick',0:2,'XTicklabels',MemberClasses_,'YTick',1:11,'XTicklabelrotation',45)
    set(gca,'XTick',[],'YTick',1:11)
%     [p{iArea},tbl{iArea},stats{iArea}] =kruskalwallis(yVals,xVals+1)
%     multcompare(stats{iArea})
    legend off
    xlabel(''); ylabel('')
    title('')
    
end
figure
for iArea = 1:2
    subplot(1,2,iArea); hold on
    for i=1:3
        stairs(bins,cumsum(yHist{iArea}(i,:)),'color',col_{i},'LineWidth',2)
    end
    axis([0 10 0 1])
end
    

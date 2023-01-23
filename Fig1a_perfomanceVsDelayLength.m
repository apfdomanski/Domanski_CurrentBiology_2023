clear
Target = 'LONG';
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
else
    pat = '/Volumes/HDD2/DNMTP/raw';
end    
cd(pat)
fileList=dir(sprintf('allTimestamps%s*%s*Events.mat',filesep,Target));
% reject_list={'IreneuszLONG1_Events.mat',...
%              'NorbertMEDIUM1_Events.mat',...
%              'NorbertMEDIUM2_Events.mat',...
%              'OnufryMEDIUM1_Events.mat',...
%              'OnufryMEDIUM2_Events.mat'}; % These only have 5s delays
reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
% reject_list={}; %'ALL_events.mat'

name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag
% fileList([1,7])=[]
%%
clear L R C_ E_ Accuracy pout
for iFile =1:length(fileList)
    clear L R C_ E_
    fname=strtok(fileList(iFile).name,'_');
    
    load(sprintf('%s%sallTimestamps%s%s_Events.mat',pat,filesep,filesep,fname));
    
    L = [length(t.Short.ChoicePress_LeftCorrect)./(length(t.Short.ChoicePress_LeftCorrect)+length(t.Short.ChoicePress_LeftError)) , ...
        length(t.Medium.ChoicePress_LeftCorrect)./(length(t.Medium.ChoicePress_LeftCorrect)+length(t.Medium.ChoicePress_LeftError)) , ...
        length(t.Long.ChoicePress_LeftCorrect)./(length(t.Long.ChoicePress_LeftCorrect)+length(t.Long.ChoicePress_LeftError))];
    
    
    R = [length(t.Short.ChoicePress_RightCorrect)./(length(t.Short.ChoicePress_RightCorrect)+length(t.Short.ChoicePress_RightError)) , ...
        length(t.Medium.ChoicePress_RightCorrect)./(length(t.Medium.ChoicePress_RightCorrect)+length(t.Medium.ChoicePress_RightError)) , ...
        length(t.Long.ChoicePress_RightCorrect)./(length(t.Long.ChoicePress_RightCorrect)+length(t.Long.ChoicePress_RightError))];
    
    Accuracy(iFile,:) = (L+R)./2;
    
    C_ =  [length(t.Short.ChoicePress_LeftCorrect)  + length(t.Short.ChoicePress_RightCorrect),...
           length(t.Medium.ChoicePress_LeftCorrect) + length(t.Medium.ChoicePress_RightCorrect),...         
           length(t.Long.ChoicePress_LeftCorrect)   + length(t.Long.ChoicePress_RightCorrect)];
  
    E_ =  [length(t.Short.ChoicePress_LeftError)    + length(t.Short.ChoicePress_RightError),...
           length(t.Medium.ChoicePress_LeftError)   + length(t.Medium.ChoicePress_RightError),...         
           length(t.Long.ChoicePress_LeftError)     + length(t.Long.ChoicePress_RightError)];
    
    pout(iFile,:) = myBinomTest(C_,C_+E_,0.5*ones(1,3),'two, equal counts')<0.05;
     


end
%%
figure;  hold on
plot([0.1 0.5],[50 50],':r','LineWidth',1.5)
plot([3.5 3.9],[50 50],':r','LineWidth',1.5)
plot(1:3,100*Accuracy(pout(:,3),:),'k')
plot(1:3,100*Accuracy(~pout(:,3),:),':k')
errorbar(1:3,100*nanmean(Accuracy),100*nansem(Accuracy),'k','LineWidth',2.5)
axis([0 4 25 100])

ylabel('% Choices correct')
set(gca,'xTick',1:3,'xTickLabel',{'4s','8s','16s'},'XTickLabelRotation',0,'TickDir','out')
get(gca)
B = bar(1:3,100*nanmean(Accuracy));
    B.EdgeColor = 'k';
    B.FaceColor = 'none';
    B.LineWidth = 1.5;
% [p,tbl,stats]=anova1(Accuracy)
% [c,~,~,gnames] = multcompare(stats,'bonferroni');

% tmcomptest(Accuracy,
% [h,p, chi2stat,df] = proportionTest(X , N, correct)

% pout=myBinomTest(1,204,'two','two') 

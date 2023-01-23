% Analyse interplay between SWRs and assembly activation times

% Run Project_Sleep_into_Task first (... including DecodeAssemFromSleepProjections)
temp = load([p.pat 'Ripples' p.sep p.rat 'Ripp.mat']);

rip.PreTimes = reshape([temp.Jlong2rip5sd_R1.tsabsmax],size(temp.Jlong2rip5sd_R1)); 
rip.TaskTimes = rip.PreTimes;
rip.PostTimes = rip.PreTimes;
% Cut out specific ripple times
rip.PreTimes(rip.PreTimes<tRanges.Pre(1)    | rip.PreTimes>tRanges.Pre(2))  = [];
rip.TaskTimes(rip.TaskTimes<tRanges.Task(1) | rip.TaskTimes>tRanges.Task(2))= [];
rip.PostTimes(rip.PostTimes<tRanges.Post(1) | rip.PostTimes>tRanges.Post(2))= [];
clear temp
%% Check that everything lines up

figure; hold on
tempTimes=FR.PFCcells{1}.t*1e-4;
scatter(tempTimes/3600,0.5+ones(length(tempTimes),1),'.k')
plot(tRanges.Pre/3600,[0.5 0.5], 'r')
plot(tRanges.Task/3600,[0.5 0.5], 'g')
plot(tRanges.Post/3600,[0.5 0.5], 'b')
plot(FR.Pre{1}.Tmtx(1:end-1)/3600,2+(FR.Pre{1}.iFR(:,1)), 'r')
plot(FR.Task{1}.Tmtx(1:end-1)/3600,2+(FR.Task{1}.iFR(:,1)), 'g')
plot(FR.Post{1}.Tmtx(1:end-1)/3600,2+(FR.Post{1}.iFR(:,1)), 'b')

scatter(rip.PreTimes/3600,ones(length(rip.PreTimes),1), 'g')
scatter(rip.TaskTimes/3600,ones(length(rip.TaskTimes),1), 'g')
scatter(rip.PostTimes/3600,ones(length(rip.PostTimes),1), 'g')

% axis([min(tempTimes) max(tempTimes) 0 Inf])
clear tempTimes
xlabel('Experiment time (h)')
%% Cutouts around Ripple times - FSCs and FRs
epochs={'Pre','Task','Post'};
cutout=round((p.twin+5)/p.bw); 
rip.cutout.Ltr=cutout;

for epochID=1:length(epochs) %
    if  ~strcmp(epochs(epochID),'Task')
        eval([ '[Tcutout,Dcutout]=SelTimesAssem(p.twin+5,FR.' epochs{epochID} '{1}.Tmtx, Ass.' epochs{epochID} '.FSCproj,rip.' epochs{epochID} 'Times);'])
        % [Tcutout,Dcutout]=SelTimesAssem(p.twin+5,FR.Pre{1}.Tmtx, Ass.Pre.FSCproj,rip.PreTimes);
    else
        eval([ '[Tcutout,Dcutout]=SelTimesAssem(p.twin+5,cell2mat(Ass.' epochs{epochID} '.TmtxS{1}), Ass.' epochs{epochID} '.FSC,rip.' epochs{epochID} 'Times);'])
    end
    
    eval(['[Tcutout_U,Dcutout_U] = SelTimesFRs(p.twin+5,FR.' epochs{epochID} '{1}.Tmtx,FR.'  epochs{epochID} ',rip.'  epochs{epochID} 'Times);'])
    
    % Excise an equal duration section for each ripple
    eval(['rip.cutout.TmtxS_' epochs{epochID} '=cell(1,3);'])
    eval(['rip.cutout.FSC_' epochs{epochID} '=cell(1,3);'])
    eval(['rip.cutout.FR_' epochs{epochID} '=cell(1,3);'])
    for s=1:3
        for iTrial= 1:length(Tcutout{1})-1
            try
                eval(['rip.cutout.TmtxS_' epochs{epochID} '{s}{iTrial}       = Tcutout{1}{iTrial}([1:cutout end-cutout+1:end]);'])
                eval(['rip.cutout.FSC_' epochs{epochID}   '{s}{iTrial}       = Dcutout{s}{iTrial}  ([1:cutout end-cutout+1:end],:);'])
            catch 
                warning(['Overhanging data on '     epochs{epochID} ' recording (Area ' num2str(s) ', Ripple # ' num2str(iTrial) ')'] )
            end
        end;    
        for iTrial= 1:length(Tcutout_U{1})-1
            try
                eval(['rip.cutout.TmtxS_FR_' epochs{epochID} '{s}{iTrial}       = Tcutout_U{1}{iTrial}([1:cutout end-cutout+1:end]);'])
                eval(['rip.cutout.FR_' epochs{epochID}   '{s}{iTrial}       = Dcutout_U{s}{iTrial}  ([1:cutout end-cutout+1:end],:);'])
            catch 
                warning(['Overhanging data on '     epochs{epochID} ' recording (Area ' num2str(s) ', Ripple # ' num2str(iTrial) ')'] )
            end
            
        end
    end;
    % rip.Tmtx=cell2mat( rip.TmtxS{1})';
end
clear Tcutout Dcutout Tcutout_U Dcutout_U cutout s iTrial

figure; hold on

for iCutout=1:length(rip.cutout.TmtxS_Pre{1})
    plot(rip.cutout.TmtxS_Pre{1}{iCutout}/3600,rip.cutout.FSC_Pre{1}{iCutout}(:,:),'r')
end
for iCutout=1:length(rip.cutout.TmtxS_Task{1})
    plot(rip.cutout.TmtxS_Task{1}{iCutout}/3600,rip.cutout.FSC_Task{1}{iCutout}(:,:),'g')
end
for iCutout=1:length(rip.cutout.TmtxS_Post{1})
    plot(rip.cutout.TmtxS_Post{1}{iCutout}/3600,rip.cutout.FSC_Post{1}{iCutout}(:,:),'b')
end
plot(rip.PreTimes/3600,-2*ones(length(rip.PreTimes),1),'or',...
     rip.TaskTimes/3600,-2*ones(length(rip.TaskTimes),1),'og',...
     rip.PostTimes/3600,-2*ones(length(rip.PostTimes),1),'ob')
clear iCutout Tcutout Dcutout epochID 
%% Ripple triggered average Assem activity
clear temp
for s=1:3
    for iAss=1:size(rip.cutout.FSC_Pre{s}{1},2)
        for iRip=1:size(rip.cutout.FSC_Pre{s},2)
            temp{iAss}(iRip,:)=zscore(rip.cutout.FSC_Pre{s}{iRip}(:,iAss));
        end
    end
    rip.cutout.Av.PreMean{s}=cell2mat(cellfun(@mean,temp,'UniformOutput' ,false)');
    rip.cutout.Av.PreSEM{s}=cell2mat(cellfun(@nansem,temp,'UniformOutput' ,false)');
    clear temp
    
    for iAss=1:size(rip.cutout.FSC_Task{s}{1},2)
        for iRip=1:size(rip.cutout.FSC_Task{s},2)
            temp{iAss}(iRip,:)=zscore(rip.cutout.FSC_Task{s}{iRip}(:,iAss));
        end
    end
    rip.cutout.Av.TaskMean{s}=cell2mat(cellfun(@mean,temp,'UniformOutput' ,false)');
	rip.cutout.Av.TaskSEM{s} =cell2mat(cellfun(@nansem,temp,'UniformOutput' ,false)');
    clear temp
    
    for iAss=1:size(rip.cutout.FSC_Post{s}{1},2)
        for iRip=1:size(rip.cutout.FSC_Post{s},2)
            temp{iAss}(iRip,:)=zscore(rip.cutout.FSC_Post{s}{iRip}(:,iAss));
        end
    end
    rip.cutout.Av.PostMean{s}=cell2mat(cellfun(@mean,temp,'UniformOutput' ,false)');
	rip.cutout.Av.PostSEM{s}=cell2mat(cellfun(@nansem,temp,'UniformOutput' ,false)');
    clear temp
end
clear iRip iAss s i Trial

tb=((1:size(rip.cutout.Av.PreMean{1},2))*p.bw);tb=tb-max(tb)/2;
for s = 1:3
    figure('name',['Ripple-triggered average Assembly activity: ' p.names{s}]) ;
    for iAss=1:size(rip.cutout.Av.PreMean{s},1)
        subaxis(ceil(size(rip.cutout.Av.PreMean{s},1)^0.5),ceil(size(rip.cutout.Av.PreMean{s},1)^0.5),iAss,'Spacing',0.01); hold on
%         subaxis(size(rip.cutout.Av.PreMean{s},1),1,iAss,'SpacingVert',0); hold on
        ciplot(rip.cutout.Av.PreMean{s}(iAss,:)-rip.cutout.Av.PreSEM{s}(iAss,:),...
               rip.cutout.Av.PreMean{s}(iAss,:)+rip.cutout.Av.PreSEM{s}(iAss,:),...
               tb,'r');
      if iAss == size(rip.cutout.Av.PostMean{s},1)
            plot([5 5], [0.5 0.9], '-k',  [5 10], [0.5 0.5], '-k', 'LineWidth', 1)
            text(4,0.7, '0.5 SD', 'HorizontalAlignment','right')
            text(7.5,0.35, '5s', 'HorizontalAlignment','center')
       end
       axis off
       axis([-15 15 -1 1])
       
    end
    for iAss=1:size(rip.cutout.Av.TaskMean{s},1)
        subaxis(ceil(size(rip.cutout.Av.TaskMean{s},1)^0.5),ceil(size(rip.cutout.Av.TaskMean{s},1)^0.5),iAss,'Spacing',0.01); hold on
%         subaxis(size(rip.cutout.Av.TaskMean{s},1),1,iAss,'SpacingVert',0); hold on
        ciplot(rip.cutout.Av.TaskMean{s}(iAss,:)-rip.cutout.Av.TaskSEM{s}(iAss,:),...
               rip.cutout.Av.TaskMean{s}(iAss,:)+rip.cutout.Av.TaskSEM{s}(iAss,:),...
               tb,'g');
       if iAss == size(rip.cutout.Av.TaskMean{s},1)
            plot([5 5], [0.5 0.9], '-k',  [5 10], [0.5 0.5], '-k', 'LineWidth', 1)
            text(4,0.7, '0.5 SD', 'HorizontalAlignment','right')
            text(7.5,0.35, '5s', 'HorizontalAlignment','center')
       end           
       axis off
       axis([-15 15 -1 1])
    end
    for iAss=1:size(rip.cutout.Av.PostMean{s},1)
        subaxis(ceil(size(rip.cutout.Av.PostMean{s},1)^0.5),ceil(size(rip.cutout.Av.PostMean{s},1)^0.5),iAss,'Spacing',0.01); hold on
%         subaxis(size(rip.cutout.Av.PostMean{s},1),1,iAss,'SpacingVert',0); hold on
        ciplot(rip.cutout.Av.PostMean{s}(iAss,:)-rip.cutout.Av.PostSEM{s}(iAss,:),...
               rip.cutout.Av.PostMean{s}(iAss,:)+rip.cutout.Av.PostSEM{s}(iAss,:),...
               tb,'b');
       if iAss == size(rip.cutout.Av.PostMean{s},1)
            plot([5 5], [0.5 0.9], '-k',  [5 10], [0.5 0.5], '-k', 'LineWidth', 1)
            text(4,0.7, '0.5 SD', 'HorizontalAlignment','right')
            text(7.5,0.35, '5s', 'HorizontalAlignment','center')
       end
       axis off
       axis([-15 15 -1 1])
       if iAss ==1
           legend('Pre-task rest','Task','Post-task rest'); legend('boxoff')
       end
       
       
       plot([0 0],[-0.5 0.5],':k')
    end
end
%% Ripple triggered average Firing rates
clear temp
for s=1:3
    for ineuron=1:size(rip.cutout.FR_Pre{s}{1},2)
        for iRip=1:size(rip.cutout.FR_Pre{s},2)
            temp{ineuron}(iRip,:)=zscore(rip.cutout.FR_Pre{s}{iRip}(:,ineuron));
        end
    end
    rip.cutout.Av.PreMeanFR{s}=cell2mat(cellfun(@mean,temp,'UniformOutput' ,false)');
    rip.cutout.Av.PreSEMFR{s}=cell2mat(cellfun(@nansem,temp,'UniformOutput' ,false)');
    clear temp
    
    for ineuron=1:size(rip.cutout.FR_Task{s}{1},2)
        for iRip=1:size(rip.cutout.FR_Task{s},2)
            temp{ineuron}(iRip,:)=zscore(rip.cutout.FR_Task{s}{iRip}(:,ineuron));
        end
    end
    rip.cutout.Av.TaskMeanFR{s}=cell2mat(cellfun(@mean,temp,'UniformOutput' ,false)');
    rip.cutout.Av.TaskSEMFR{s}=cell2mat(cellfun(@nansem,temp,'UniformOutput' ,false)');
    clear temp
    
    for ineuron=1:size(rip.cutout.FR_Post{s}{1},2)
        for iRip=1:size(rip.cutout.FR_Post{s},2)
            temp{ineuron}(iRip,:)=zscore(rip.cutout.FR_Post{s}{iRip}(:,ineuron));
        end
    end
    rip.cutout.Av.PostMeanFR{s}=cell2mat(cellfun(@mean,temp,'UniformOutput' ,false)');
    rip.cutout.Av.PostSEMFR{s}=cell2mat(cellfun(@nansem,temp,'UniformOutput' ,false)');
    clear temp
    
end
clear iRip ineuron s 

tb=((1:size(rip.cutout.Av.PreMeanFR{1},2))*p.bw);tb=tb-max(tb)/2;
for s = 1:2
    figure('name',['Ripple-triggered average Unit activity: ' p.names{s}]) ;
    for ineuron=1:size(rip.cutout.Av.PreMeanFR{s},1)
        subaxis(ceil(size(rip.cutout.Av.PreMeanFR{s},1)^0.5),ceil(size(rip.cutout.Av.PreMeanFR{s},1)^0.5),ineuron,'Spacing',0.01); hold on
%         subaxis(size(rip.cutout.Av.PreMeanFR{s},1),1,ineuron,'SpacingVert',0); hold on
        ciplot(rip.cutout.Av.PreMeanFR{s}(ineuron,:)-rip.cutout.Av.PreSEMFR{s}(ineuron,:),...
               rip.cutout.Av.PreMeanFR{s}(ineuron,:)+rip.cutout.Av.PreSEMFR{s}(ineuron,:),...
               tb,'r');
      if ineuron == size(rip.cutout.Av.PostMeanFR{s},1)
            plot([5 5], [0.5 0.9], '-k',  [5 10], [0.5 0.5], '-k', 'LineWidth', 1)
            text(4,0.7, '0.5 SD', 'HorizontalAlignment','right')
            text(7.5,0.35, '5s', 'HorizontalAlignment','center')
       end
       axis off
       axis([-15 15 -1 1])
       
    end
    for ineuron=1:size(rip.cutout.Av.TaskMeanFR{s},1)
        subaxis(ceil(size(rip.cutout.Av.TaskMeanFR{s},1)^0.5),ceil(size(rip.cutout.Av.TaskMeanFR{s},1)^0.5),ineuron,'Spacing',0.01); hold on
%         subaxis(size(rip.cutout.Av.TaskMeanFR{s},1),1,ineuron,'SpacingVert',0); hold on
        ciplot(rip.cutout.Av.TaskMeanFR{s}(ineuron,:)-rip.cutout.Av.TaskSEMFR{s}(ineuron,:),...
               rip.cutout.Av.TaskMeanFR{s}(ineuron,:)+rip.cutout.Av.TaskSEMFR{s}(ineuron,:),...
               tb,'g');
       if ineuron == size(rip.cutout.Av.TaskMeanFR{s},1)
            plot([5 5], [0.5 0.9], '-k',  [5 10], [0.5 0.5], '-k', 'LineWidth', 1)
            text(4,0.7, '0.5 SD', 'HorizontalAlignment','right')
            text(7.5,0.35, '5s', 'HorizontalAlignment','center')
       end           
       axis off
       axis([-15 15 -1 1])
    end
    for ineuron=1:size(rip.cutout.Av.PostMeanFR{s},1)
        subaxis(ceil(size(rip.cutout.Av.PostMeanFR{s},1)^0.5),ceil(size(rip.cutout.Av.PostMeanFR{s},1)^0.5),ineuron,'Spacing',0.01); hold on
%         subaxis(size(rip.cutout.Av.PostMeanFR{s},1),1,ineuron,'SpacingVert',0); hold on
        ciplot(rip.cutout.Av.PostMeanFR{s}(ineuron,:)-rip.cutout.Av.PostSEMFR{s}(ineuron,:),...
               rip.cutout.Av.PostMeanFR{s}(ineuron,:)+rip.cutout.Av.PostSEMFR{s}(ineuron,:),...
               tb,'b');
       if ineuron == size(rip.cutout.Av.PostMeanFR{s},1)
            plot([5 5], [0.5 0.9], '-k',  [5 10], [0.5 0.5], '-k', 'LineWidth', 1)
            text(4,0.7, '0.5 SD', 'HorizontalAlignment','right')
            text(7.5,0.35, '5s', 'HorizontalAlignment','center')
       end
       axis off
       axis([-15 15 -1 1])
       if ineuron ==1
           legend('Pre-task rest','Task','Post-task rest'); legend('boxoff')
       end
       
       
       plot([0 0],[-0.5 0.5],':k')
    end
end
clear epochID ineuron s 
%% Rank-ordered FRs
clims=[-1 1];

load('blue_white_red.mat')
figure('name','Ripple modulation of units')
for s=1:2
    temp=rip.cutout.Av.PreMeanFR{s};
    temp=sortrows(temp,find(tb==0));
    subplot(2,3,3*s-2)
    imagesc(tb,1:size(temp,1),flipud(temp))
    ylabel(strcat(p.names(s), ' Unit no.'))
    colormap(cmap)
    caxis(clims)
    if s==1, title(strcat(epochs(1), '-task rest'),'Fontsize',12); end
end
for s=1:2
    temp=rip.cutout.Av.TaskMeanFR{s};
    temp=sortrows(temp,find(tb==0));
    subplot(2,3,3*s-1)
    imagesc(tb,1:size(temp,1),flipud(temp))
    colormap(cmap)
    caxis(clims)
    if s==1, title(epochs(2),'Fontsize',12); end
    if s==2, xlabel('Peri-ripple time (s)'); end
    % ylabel('Unit no.')
end
for s=1:2
    temp=rip.cutout.Av.PostMeanFR{s};
    temp=sortrows(temp,find(tb==0));
    subplot(2,3,3*s)
    imagesc(tb,1:size(temp,1),flipud(temp))
    if s==1, title(strcat(epochs(3), '-task rest'),'Fontsize',12); end
    colormap(cmap)
    caxis(clims)
end    
%% Rank-ordered Assems
load('blue_white_red.mat')
figure('name','Ripple modulation of assemblies')
for s=1:3
    temp=rip.cutout.Av.PreMean{s};
    temp=sortrows(temp,find(tb==0));
    subplot(3,3,s+2*s-2)
    imagesc(tb,1:size(temp,1),flipud(temp))
    colormap(cmap)
    caxis(clims)
    if s==1, title(strcat(epochs(1), '-task rest'),'Fontsize',12); end
    ylabel(strcat(p.names{s}, ' Assem no.'))
end
for s=1:3
    temp=rip.cutout.Av.TaskMean{s};
    temp=sortrows(temp,find(tb==0));
    subplot(3,3,s+2*s-1)
    imagesc(tb,1:size(temp,1),flipud(temp))
    colormap(cmap)
    caxis(clims)
    if s==1, title(epochs(2),'Fontsize',12); end
    if s==3,  xlabel('Peri-ripple time (s)'); end
end
% figure
for s=1:3
    temp=rip.cutout.Av.PostMean{s};
    temp=sortrows(temp,find(tb==0));
    subplot(3,3,s+2*s)
    imagesc(tb,1:size(temp,1),flipud(temp))
    colormap(cmap)
    caxis(clims)
    if s==1, title(strcat(epochs(3), '-task rest'),'Fontsize',12); end
end
%% Zero-lag ripple modulation before/during/after
figure('name','Assemblies')
for s=1:3
    subplot(3,1,s)
    temp=[rip.cutout.Av.PreMean{s}(:,find(tb==0)),...
          rip.cutout.Av.TaskMean{s}(:,find(tb==0)),...   
          rip.cutout.Av.PostMean{s}(:,find(tb==0))];
%       temp=temp./repmat(temp(:,1),1,3);
    plot(repmat([1,2,3],size(temp,1),1)',temp','-ob')
   set(gca,'XTick',[1 2 3]) 
   set(gca,'XTickLabel',{'Pre-Task', 'Task ', 'Post-Task'}) 
   ylabel('Ripple Firing-rate modulation (S.D.)')
   axis ([0.5 3.5 -0.5 Inf])
   title(strcat(p.names(s)),'Fontsize',12)

end

figure('name','Units')
for s=1:2
    subplot(2,1,s)
    temp=[rip.cutout.Av.PreMeanFR{s}(:,find(tb==0)),...
          rip.cutout.Av.TaskMeanFR{s}(:,find(tb==0)),...   
          rip.cutout.Av.PostMeanFR{s}(:,find(tb==0))];
%       temp=temp./repmat(temp(:,1),1,3);
    plot(repmat([1,2,3],size(temp,1),1)',temp','-ob')
   set(gca,'XTick',[1 2 3]) 
   set(gca,'XTickLabel',{'Pre-Task', 'Task ', 'Post-Task'}) 
   ylabel('Ripple Firing-rate modulation (S.D.)')
   axis ([0.5 3.5 -0.5 Inf])
   title(strcat(p.names(s)),'Fontsize',12)

end
%% Ripple modulation vs. decoding power - Units
figure('color','w','name','Unit ripple modulation vs decoding strength');
cmap = {'r','g','b'};
clear xPre xTask xPost y fPre fTask fPost
lims=[-0.6 1 0 6];
for s=1:2
    xPre  = rip.cutout.Av.PreMeanFR{s}(Ass.Task.unit_IDs{s},(tb==0));
    xTask = rip.cutout.Av.TaskMeanFR{s}(Ass.Task.unit_IDs{s},(tb==0));
    xPost = rip.cutout.Av.PostMeanFR{s}(Ass.Task.unit_IDs{s},(tb==0));        
    y    = max(Decoding.units.TS{s},[],1);
    
    fPre = ezfit(xPre,y,'a*x+b');
    fTask = ezfit(xTask,y,'a*x+b');
    fPost = ezfit(xPost,y,'a*x+b');
    
    
    subplot(2,3,3*s-2); hold on
        plot(xPre,fPre.m(1)*xPre+fPre.m(2),'LineWidth',1.5, 'color',cmap{1})
        scatter(xPre,y,cmap{1})
        text(0.5,5.5,strcat('R^2=', sprintf('%3.1g',fPre.r^2)),'color',cmap{1})
        ylabel(strcat(p.names(s), ' Unit decoding strength (t-score)'))
        axis(lims)
        if s==1, title(strcat(epochs(1), '-task rest'),'Fontsize',12); end

    subplot(2,3,3*s-1); hold on
        plot(xTask,fTask.m(1)*xTask+fTask.m(2),'LineWidth',1.5, 'color',cmap{2})
        scatter(xTask,y,cmap{2})
        text(0.5,5.5,strcat('R^2=', sprintf('%3.1g',fTask.r^2)),'color',cmap{2})
        if s==2, xlabel('Ripple FR modulation (S.D.)'); end
        axis(lims)
        if s==1, title(epochs(2),'Fontsize',12); end

    subplot(2,3,3*s); hold on
        plot(xPost,fPost.m(1)*xPost+fPost.m(2),'LineWidth',1.5, 'color',cmap{3})
        scatter(xPost,y,cmap{3})
        text(0.5,5.5,strcat('R^2=', sprintf('%3.1g',fPost.r^2)),'color',cmap{3})
        axis(lims)
        if s==1, title(strcat(epochs(3), '-task rest'),'Fontsize',12); end

    
end
%% Ripple modulation vs. decoding power - Assemblies

figure('color','w','name','Assembly ripple modulation vs decoding strength');
cmap = {'r','g','b'};
clear xPre xTask xPost y fPre fTask fPost
lims=[-1 1 0 6];
for s=1:3
    xPre  = rip.cutout.Av.PreMean{s}(:,(tb==0));
    xTask = rip.cutout.Av.TaskMean{s}(:,(tb==0));
    xPost = rip.cutout.Av.PostMean{s}(:,(tb==0));        
    y    = max(Decoding.Assem.TS{s},[],1);
    
    fPre = ezfit(xPre,y,'a*x+b');
    fTask = ezfit(xTask,y,'a*x+b');
    fPost = ezfit(xPost,y,'a*x+b');
    
    subplot(3,3,s+2*s-2); hold on
        plot(xPre,fPre.m(1)*xPre+fPre.m(2),'LineWidth',1.5, 'color',cmap{1})
        scatter(xPre,y,cmap{1})
        text(0.5,5,strcat('R^2=', sprintf('%3.1g',fPre.r^2)),'color',cmap{1})
        axis(lims)
        ylabel(strcat(p.names(s), ' decoding (t-score)'))
        if s==1, title(strcat(epochs(1), '-task rest'),'Fontsize',12); end

    subplot(3,3,s+2*s-1); hold on
        plot(xTask,fTask.m(1)*xTask+fTask.m(2),'LineWidth',1.5, 'color',cmap{2})
        scatter(xTask,y,cmap{2})
        text(0.5,5,strcat('R^2=', sprintf('%3.1g',fTask.r^2)),'color',cmap{2})  
        axis(lims)
        if s==1, title(epochs(2),'Fontsize',12); end
        if s==3, xlabel('Ripple activation modulation (S.D.)'); end
        
    subplot(3,3,s+2*s); hold on
        plot(xPost,fPost.m(1)*xPost+fPost.m(2),'LineWidth',1.5, 'color',cmap{3})
        scatter(xPost,y,cmap{3})
        text(0.5,5,strcat('R^2=', sprintf('%3.1g',fPost.r^2)),'color',cmap{3})
        axis(lims)
        if s==1, title(strcat(epochs(3), '-task rest'),'Fontsize',12); end
        
end
%% Change in ripple modulation vs decoding score - Units

figure('color','w','name','Unit ripple modulation vs decoding strength');
cmap = {'r','g','b'};
clear xPreTask xTaskPost y fPreTask fTaskPost
lims=[-0.5 1 0 6];
for s=1:2
    xPreTask  = rip.cutout.Av.PreMeanFR{s}(Ass.Task.unit_IDs{s},(tb==0)) -...
                rip.cutout.Av.TaskMeanFR{s}(Ass.Task.unit_IDs{s},(tb==0));
    xTaskPost = rip.cutout.Av.PostMeanFR{s}(Ass.Task.unit_IDs{s},(tb==0)) -...
                rip.cutout.Av.TaskMeanFR{s}(Ass.Task.unit_IDs{s},(tb==0));
    y    = max(Decoding.units.TS{s},[],1);
    
    fPreTask  = ezfit(xPreTask,y,'a*x+b');
    fTaskPost = ezfit(xTaskPost,y,'a*x+b');
    
    
    subplot(2,2,2*s-1); hold on
        plot(xPreTask,fPreTask.m(1)*xPreTask+fPreTask.m(2),'LineWidth',1.5, 'color','k')
        scatter(xPreTask,y,'k')
        text(0.5,5.5,strcat('R^2=', sprintf('%3.1g',fPreTask.r^2)),'color','k')
        ylabel(strcat(p.names(s), ' Unit decoding strength (t-score)'))
        axis(lims)
%         if s==2, xlabel('Change in ripple FR modulation (S.D.)'); end
        if s==1, title('Pre-Task to Task change','Fontsize',12); end

    subplot(2,2,2*s); hold on
        plot(xTaskPost,fTaskPost.m(1)*xTaskPost+fTaskPost.m(2),'LineWidth',1.5, 'color','k')
        scatter(xTaskPost,y,'k')
        text(0.5,5.5,strcat('R^2=', sprintf('%3.1g',fTaskPost.r^2)),'color','k')
        if s==2, xlabel('Change in ripple FR modulation (S.D.)', 'HorizontalAlignment','right'); end
        if s==1, title('Task to Post-Task change','Fontsize',12); end
        axis(lims)
end
%% Change in ripple modulation vs. decoding power - Assemblies
figure('color','w','name','Assembly ripple modulation vs decoding strength');
cmap = {'r','g','b'};
clear xPre xTask xPost y fPre fTask fPost
lims=[-1 1 0 6];
for s=1:3
    xPreTask  = rip.cutout.Av.PreMean{s}(:,(tb==0))-rip.cutout.Av.TaskMean{s}(:,(tb==0));
    xTaskPost = rip.cutout.Av.PostMean{s}(:,(tb==0))-rip.cutout.Av.TaskMean{s}(:,(tb==0));
    y    = max(Decoding.Assem.TS{s},[],1);
    
    fPreTask  = ezfit(xPreTask,y,'a*x+b');
    fTaskPost = ezfit(xTaskPost,y,'a*x+b');
    
    subplot(3,2,2*s-1); hold on
        plot(xPreTask,fPreTask.m(1)*xPreTask+fPreTask.m(2),'LineWidth',1.5, 'color','k')
        scatter(xPreTask,y,'k')
        text(0.5,5.5,strcat('R^2=', sprintf('%3.1g',fPreTask.r^2)),'color','k')
        axis(lims)
        ylabel(strcat(p.names(s), ' decoding (t-score)'))
        if s==1, title('Pre-Task to Task change','Fontsize',12); end

    subplot(3,2,2*s); hold on
        plot(xTaskPost,fTaskPost.m(1)*xTaskPost+fTaskPost.m(2),'LineWidth',1.5, 'color','k')
        scatter(xTaskPost,y,'k')
        text(0.5,5.5,strcat('R^2=', sprintf('%3.1g',fTaskPost.r^2)),'color','k')
        axis(lims)
        if s==1, title('Task to Post-Task change','Fontsize',12); end
        if s==3, xlabel('Change in Assembly ripple modulation (S.D.)', 'HorizontalAlignment','right'); end
        
   
        
end
clear lims s tb temp xPre xPost xPreTask xTask xTaskPost lastfit y iRip iAss fPre fPost fPreTask fTask fTaskPost cmap
%% Load and process behavioural timestamps
behav = load([p.pat2 p.sep p.rat],'*trange*');

tRanges.trange_sample = [behav.trangeleft_sample;behav.trangeright_sample];
tRanges.trange_choice = [behav.trangeleft_choice;behav.trangeright_choice];
tRanges.trange_sampleError =[behav.ERRORtrangeleft_sample;behav.ERRORtrangeright_sample];
tRanges.trange_choiceError =[behav.ERRORtrangeleft_choice;behav.ERRORtrangeright_choice];

tRanges.trange_sample = tRanges.trange_sample(:,1)*1e-6+5;
tRanges.trange_choice = tRanges.trange_choice(:,1)*1e-6+5;
tRanges.trange_sampleError = tRanges.trange_sampleError(:,1)*1e-6+5;
tRanges.trange_choiceError = tRanges.trange_choiceError(:,1)*1e-6+5;
%% make a peri-event histogram of all ripple events
rip.perievent.bins =(-20:2:20)';
rip.perievent.samplehist=[];
rip.perievent.choicehist=[];
rip.perievent.Errorsamplehist=[];
rip.perievent.Errorchoicehist=[];
for iEvent=1:length(tRanges.trange_sample)
    rip.perievent.samplehist(:,iEvent) =  histc(rip.TaskTimes-tRanges.trange_sample(iEvent),rip.perievent.bins); 
end
for iEvent=1:length(tRanges.trange_choice)
    rip.perievent.choicehist(:,iEvent) =  histc(rip.TaskTimes-tRanges.trange_choice(iEvent),rip.perievent.bins); 
end
for iEvent=1:length(tRanges.trange_sampleError)
    rip.perievent.Errorsamplehist(:,iEvent) =  histc(rip.TaskTimes-tRanges.trange_sampleError(iEvent),rip.perievent.bins); 
end
for iEvent=1:length(tRanges.trange_choiceError)
    rip.perievent.Errorchoicehist(:,iEvent) =  histc(rip.TaskTimes-tRanges.trange_choiceError(iEvent),rip.perievent.bins); 
end

% Count how many ripples take place specifically before and after the evets
rip.perievent.preSampleCount  = sum(rip.perievent.samplehist(rip.perievent.bins<0,:));
rip.perievent.postSampleCount = sum(rip.perievent.samplehist(rip.perievent.bins>0,:));
rip.perievent.preChoiceCount  = sum(rip.perievent.choicehist(rip.perievent.bins<0,:));
rip.perievent.postChoiceCount = sum(rip.perievent.choicehist(rip.perievent.bins>0,:));

rip.perievent.preSampleCountError  = sum(rip.perievent.Errorsamplehist(rip.perievent.bins<0,:));
rip.perievent.postSampleCountError = sum(rip.perievent.Errorsamplehist(rip.perievent.bins>0,:));
rip.perievent.preChoiceCountError  = sum(rip.perievent.Errorchoicehist(rip.perievent.bins<0,:));
rip.perievent.postChoiceCountError = sum(rip.perievent.Errorchoicehist(rip.perievent.bins>0,:));

% Collapse to aveage histogram
rip.perievent.samplehist      = sum(rip.perievent.samplehist,2);
rip.perievent.choicehist      = sum(rip.perievent.choicehist,2);
rip.perievent.Errorsamplehist = sum(rip.perievent.Errorsamplehist,2);
rip.perievent.Errorchoicehist = sum(rip.perievent.Errorchoicehist,2);

rip.perievent.samplehist      = rip.perievent.samplehist./sum(rip.perievent.samplehist);
rip.perievent.choicehist      = rip.perievent.choicehist./sum(rip.perievent.choicehist);
rip.perievent.Errorsamplehist = rip.perievent.Errorsamplehist./sum(rip.perievent.Errorsamplehist);
rip.perievent.Errorchoicehist = rip.perievent.Errorchoicehist./sum(rip.perievent.Errorchoicehist);

% No. ripples in delay period (using sample-choice times as proxy)
for iEvent = 1:length(tRanges.trange_sample)
    temp=rip.TaskTimes;
    temp(temp<tRanges.trange_sample(iEvent) | temp>tRanges.trange_choice(iEvent))=[];
    if ~isempty(temp)
        rip.perievent.delaycountCorrect(iEvent) =numel(temp)/(tRanges.trange_choice(iEvent)- tRanges.trange_sample(iEvent));
    else
        rip.perievent.delaycountCorrect(iEvent) = 0;
    end
end
for iEvent = 1:length(tRanges.trange_sampleError)
    temp=rip.TaskTimes;
    temp(temp<tRanges.trange_sampleError(iEvent) | temp>tRanges.trange_choiceError(iEvent))=[];
    if ~isempty(temp)
        rip.perievent.delaycountError(iEvent) =numel(temp)/(tRanges.trange_choiceError(iEvent)- tRanges.trange_sampleError(iEvent));
    else
        rip.perievent.delaycountError(iEvent) = 0;
    end
end
%% Plot ripple histogram
figure
subplot(2,1,1); hold on
    stairs(rip.perievent.bins,(rip.perievent.samplehist),'b','LineWidth',1.5)
    stairs(rip.perievent.bins,(rip.perievent.Errorsamplehist),'r','LineWidth',1.5)
    axis([min(rip.perievent.bins) max(rip.perievent.bins) 0 0.2])
    plot([0 0],[0 0.2],':k','LineWidth',1.5)
    title('Ripples around sample presses')
	legend('Correct Trials','Error Trials','Location','NorthEast'); legend('boxoff')
    ylabel('Norm. ripple count')
subplot(2,1,2); hold on
    stairs(rip.perievent.bins,(rip.perievent.choicehist),'b','LineWidth',1.5)
    stairs(rip.perievent.bins,(rip.perievent.Errorchoicehist),'r','LineWidth',1.5)
    axis([min(rip.perievent.bins) max(rip.perievent.bins) 0 0.2])
    plot([0 0],[0 0.2],':k','LineWidth',1.5)
    xlabel('Time around lever press (s)')
    ylabel('Norm. ripple count')
    title('Ripples around choice presses')
%% test ripple rate in delay period 
    figure; hold on
    y=[mean(rip.perievent.delaycountCorrect);mean(rip.perievent.delaycountError)];
    e=[nansem(rip.perievent.delaycountCorrect);nansem(rip.perievent.delaycountError)];
    bar(y,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0 0 0],'LineWidth',1.5)
        errorbar(y,e,'LineStyle','none','LineWidth',1.5,'Color',[0 0 0 ])

    set(gca,'XTick',[1 2],'XTicklabel',{'Correct trials', 'Error trials'})
    ylabel('Norm. rate'); title('Ripple rate in delay period')
    pValue=ranksum(rip.perievent.delaycountCorrect,rip.perievent.delaycountError)
    if pValue<0.05;
       text(1.5,1.6*max(y),strcat('* p=', sprintf('%.2g',pValue)), 'HorizontalAlignment','center')
       plot([1 2],[1.5*max(y) 1.5*max(y)],'-k' ,'LineWidth',1.5)
    end
    
    axis([0.5 2.5 0 2*max(y)])    
%% test ripple rate in pre-sample period 
    figure; hold on
    y=[mean(rip.perievent.preSampleCount);mean(rip.perievent.preSampleCountError)];
    e=[nansem(rip.perievent.preSampleCount);nansem(rip.perievent.preSampleCountError)];
    bar(y,'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0 0 0],'LineWidth',1.5)
        errorbar(y,e,'LineStyle','none','LineWidth',1.5,'Color',[0 0 0 ])

    set(gca,'XTick',[1 2],'XTicklabel',{'Correct trials', 'Error trials'})
    set(gca,'YTick',[0 1 2 3])
    ylabel('Ripple count /20s'); title('Ripple count in run-up to sample press')
    pValue=ranksum(rip.perievent.preSampleCount,rip.perievent.preSampleCountError)
           plot([1 2],[1.5*max(y) 1.5*max(y)],'-k' ,'LineWidth',1.5)

    if pValue<0.05;
       text(1.5,1.6*max(y),strcat('* p=', sprintf('%.2g',pValue)), 'HorizontalAlignment','center')
    else
       text(1.5,1.6*max(y),'N.S.', 'HorizontalAlignment','center')
    end
    
    axis([0.5 2.5 0 2*max(y)])    
    
    
    clear temp y pValue iEvent epochs e clims 
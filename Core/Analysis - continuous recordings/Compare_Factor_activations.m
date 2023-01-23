% Explores behaviour of Assemblies in Pre-sleep/Task/Post-sleep epochs
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% NB Run Explore_Factor_activations.m first for each of Pre/Task/Post datasets.
% 
% NB currently this code assumes that FA Assembly detection is run for each
% experimental epoch independently thus unit and assembly assignments may
% (will?) change between periods
% 
% Contains code to plot...
% (1) Relationship of single units to assemblies during Pre/Task/Post periods.
% (2) Promiscuity of units between factors as a function of area/task phase
% (3) Differences in Assembly sequence transition probabilities during each task phase

%%%%%%%%%
%% Preamble: load data
clear 
pat = 'C:\Analysis\AssemblyAnalysis\Sleep\';
cd(pat);
load([pat ,'PreSleep\' ,'KrzysztofLONG2_iFR50_PreSleep_Ass.mat'],'Ass')
Kryzs.Pre = Ass; clear Ass

load([pat ,'Task\' ,'KrzysztofLONG2_iFR50_Task_Ass.mat'],'Ass')
Kryzs.Task = Ass; clear Ass

load([pat ,'PostSleep\' ,'KrzysztofLONG2_iFR50_PostSleep_Ass.mat'],'Ass')
Kryzs.Post = Ass; clear Ass

%% Plot unit membership raster
tempCM=[1 1 1; 0.3 0.3 0.8];
figure('name','Unit membership graph (sorted to loading to factor 1)','color','w');
 subplot(3,3,1)
 temp=[Kryzs.Pre.LD{1}(:,1),Kryzs.Pre.units.isMem{1}]; temp=sortrows(temp,1); temp(:,1)=[];
 imagesc(flipud(temp));set(gca,'YDir','normal'); colormap(tempCM)
 ylabel('PFC Unit no.'); title('Pre-sleep')
 subplot(3,3,2)
 temp=[Kryzs.Task.LD{1}(:,1),Kryzs.Task.units.isMem{1}]; temp=sortrows(temp,1); temp(:,1)=[];
 imagesc(flipud(temp));set(gca,'YDir','normal'); colormap(tempCM)
 title('Task')
 subplot(3,3,3)
 temp=[Kryzs.Post.LD{1}(:,1),Kryzs.Post.units.isMem{1}]; temp=sortrows(temp,1); temp(:,1)=[];
 imagesc(flipud(temp));set(gca,'YDir','normal'); colormap(tempCM)
 title('Post-sleep')

 subplot(3,3,4)
 temp=[Kryzs.Pre.LD{2}(:,1),Kryzs.Pre.units.isMem{2}]; temp=sortrows(temp,1); temp(:,1)=[];
 imagesc(flipud(temp));set(gca,'YDir','normal'); colormap(tempCM) 
 ylabel('HP Unit no.');
 subplot(3,3,5)
 temp=[Kryzs.Task.LD{2}(:,1),Kryzs.Task.units.isMem{2}]; temp=sortrows(temp,1); temp(:,1)=[];
 imagesc(flipud(temp));set(gca,'YDir','normal'); colormap(tempCM) 
 subplot(3,3,6)
 temp=[Kryzs.Post.LD{2}(:,1),Kryzs.Post.units.isMem{2}]; temp=sortrows(temp,1); temp(:,1)=[];
 imagesc(flipud(temp));set(gca,'YDir','normal'); colormap(tempCM)
 
 subplot(3,3,7)
 temp=[Kryzs.Pre.LD{3}(:,1),Kryzs.Pre.units.isMem{3}]; temp=sortrows(temp,1); temp(:,1)=[];
 imagesc(flipud(temp));set(gca,'YDir','normal'); colormap(tempCM)  
 ylabel('PFC/HP Unit no.');
 subplot(3,3,8)
 temp=[Kryzs.Task.LD{3}(:,1),Kryzs.Task.units.isMem{3}]; temp=sortrows(temp,1); temp(:,1)=[];
 imagesc(flipud(temp));set(gca,'YDir','normal'); colormap(tempCM) 
 xlabel('Assembly no.'); 
 subplot(3,3,9)
 temp=[Kryzs.Post.LD{3}(:,1),Kryzs.Post.units.isMem{3}]; temp=sortrows(temp,1); temp(:,1)=[];
 imagesc(flipud(temp));set(gca,'YDir','normal'); colormap(tempCM) 
%% Plot soloists vs. choristers 
cmap = {'r','g','b'};
labels= {'Pre','Task','Post'}
figure('color','w','name','Fraction of cells assigned to a single assembly'); 
for s = 1:3
    subplot(1,3,s); hold on 
    title(Kryzs.Pre.titles{s})
        bar(1,Kryzs.Pre.units.FractionSingleAssem(s),'EdgeColor', cmap{s},'FaceColor', cmap{s},'LineWidth',2); 
        bar(2,Kryzs.Task.units.FractionSingleAssem(s),'EdgeColor', cmap{s},'FaceColor', cmap{s},'LineWidth',2); 
        bar(3,Kryzs.Post.units.FractionSingleAssem(s),'EdgeColor', cmap{s},'FaceColor', cmap{s},'LineWidth',2);
        set(gca,'XTickLabel',labels,'XTick',[1 2 3],'Ylim',[0 1]);
        if s==1,ylabel('Fraction of units assigned to a single assembly'), end
        set(gca,'ylim',[0 1])
end
%
figure('color','w','name','Number of factors each unit shared among'); 
for s = 1:3
    subplot(1,3,s); hold on
         title(Kryzs.Pre.titles{s})
        plot(Kryzs.Pre.units.AssOverlap_hist_bins,  Kryzs.Pre.units.AssOverlap_hist(:,s), 'color',cmap{s},'LineWidth',2,'LineStyle',':'); 
        plot(Kryzs.Task.units.AssOverlap_hist_bins, Kryzs.Task.units.AssOverlap_hist(:,s),'color',cmap{s},'LineWidth',2,'LineStyle','-'); 
        plot(Kryzs.Post.units.AssOverlap_hist_bins, Kryzs.Post.units.AssOverlap_hist(:,s),'color',cmap{s},'LineWidth',2,'LineStyle','-.'); 
        
        if s==1, ylabel('Fraction of units'); end
        if s==2, xlabel('No. assemblies units shared by '); end
        set(gca,'ylim',[0 1])
end
        legend(labels); legend('boxoff')
%% Deal with tranistion matrices
% subplot(3,3,1); 
r=Kryzs.Pre.pattern.trans_matrix{1}
r=Ass.pattern.trans_matrix{4};
% normr_data = (r - min(min(r))) ./ ( max(max(r)) - min(min(r)) );
schemaball(r,[],[0 0 0;1 1 0],[1 1 0])
adj2pajek2(r*1e4,'pre_graph',[])
% graph_to_dot(r,'filename','pre_graph.dot','directed',1)
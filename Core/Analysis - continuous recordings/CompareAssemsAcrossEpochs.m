% Loads and plots factor loading (assembly membership of units into assemblies during each experimental phase
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% NB assumes factors have been calculated on each area independently, thus
% units may not line up based on exlusion criteria not being met in each
% epoch...

%% Preamble: Load data
clear
pat= 'C:\Analysis\AssemblyAnalysis\Sleep'
cd(pat)
% vars={'FL','nassem'}
Kry1.Pre  = load([pat     '\PreSleep\KrzysztofLONG1_iFR50_PreSleep_AssemRes2.mat'], 'FL','nassem')
Kry1.Task = load([pat       '\Task\KrzysztofLONG1_iFR50_Task_AssemRes2.mat'], 'FL','nassem')
Kry1.Post = load([pat '\PostSleep\KrzysztofLONG1_iFR50_PostSleep_AssemRes2.mat'], 'FL','nassem')


%% Plot loading of cells into an example 
BrainStructure=1;
figure; 
nAss=Kry1.Pre.nassem{BrainStructure}(3);
for AssemID=1:nAss;

subaxis(ceil(nAss^0.5),ceil(nAss^0.5),AssemID,'SpacingVert',0.02); hold on
plot(Kry1.Pre.FL{BrainStructure}{nAss}(:,AssemID),'r')
plot(Kry1.Task.FL{BrainStructure}{nAss}(:,AssemID),'g')
plot(Kry1.Post.FL{BrainStructure}{nAss}(:,AssemID),'b')
set(gca,'Xtick',[])
if AssemID==1
    title(['Assembly no. ' num2str(AssemID)])
end
if AssemID==nAss
    title(['Assembly no. ' num2str(AssemID)])
end
end
xlabel('Unit ID')
ylabel('Assembly loading ')
legend('Before','During','After'); legend('Boxoff')
% Fig 2 - proportion of units in assemblies

% Delays_ = {'Short','Medium','Long'};
Areas = {'HP','PFC'};

pat = 'C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\LONGTaskonly\';
cd(pat)
%%%%%
fileList=dir('*_FSC.mat');
%%
for iFile =1:length(fileList)
    %%%%%
    % fname=strtok(fileList(iFile).name,'_');
    fname=fileList(iFile).name;
    load(fullfile(pat,fname),'units','nu');
    
    for iArea = 1:2
        JointAssUnits = unique(cell2mat(units{3}));
        AssUnits = unique(cell2mat(units{iArea}));
        nonAssUnits = setdiff(1:nu(iArea),AssUnits);
        
        % (1) Fraction of neurons that are assembly members
        if ~isempty(AssUnits)
            fracAssUnits(iArea,iFile) = length(AssUnits)./nu(iArea);
        else
            fracAssUnits(iArea,iFile) = NaN;
        end
        % (2) Fraction of neurons bridging local assemblies
        
        A = cell2mat(units{iArea});
        [unqA,~,id] = unique(A);
        B = unqA(histc(id,1:max(id))>1); if isempty(B), B = 0; end
        fracLocalNodes(iArea,iFile) = length(B)./length(A);
        
        % (3) Fraction of inter-area assemblies that bridge both local and long range assemblies
        
        switch iArea
            case 1
                JointAssUnits(JointAssUnits>nu(1))=[];
            case 2
                JointAssUnits = JointAssUnits-nu(1);
                JointAssUnits(JointAssUnits<0)=[];
        end
        bridgeUnits = intersect(AssUnits,JointAssUnits);
        fracBridgeNeurons(iArea,iFile) = length(bridgeUnits)./length(JointAssUnits);
        
        fracJointNeurons(iArea,iFile) = length(JointAssUnits)./nu(iArea);

                
                
        clear A unqA id B 
    end
    clear nu units
end
%% Plot fraction of neurons in assemblies - local

% Ignore datasets with no detected assemblies
% fracAssUnits(fracAssUnits==0)=NaN;

m_ = nanmean(fracAssUnits,2);
e_ = nansem(fracAssUnits,2);
figure; hold on
B = bar([1,2],m_*100);
B.FaceColor = [1 1 1];
B.EdgeColor = [0 0 0];
B.LineWidth = 1.5;

EB = errorbar([1,2],m_*100,e_*100,'.b');
EB.Color = [0 0 0];
EB.LineWidth = 1.5;
ylim([0 50])

set(gca,'XTick',[1,2],'XTickLabel',Areas)
% title('Fraction of units in detected assemblies')
ylabel('Fraction of units in assembly type (%)');
%% Plot fraction of neurons in assemblies - joint

% Ignore datasets with no detected assemblies
fracJointNeurons(fracJointNeurons==0)=NaN;

m_ = nanmean(fracJointNeurons,2);
e_ = nansem(fracJointNeurons,2);
figure; hold on
B = bar([1,2],m_*100);
B.FaceColor = [1 1 1];
B.EdgeColor = [0 0 0];
B.LineWidth = 1.5;

EB = errorbar([1,2],m_*100,e_*100,'.b');
EB.Color = [0 0 0];
EB.LineWidth = 1.5;
ylim([0 50])

set(gca,'XTick',[1,2],'XTickLabel',Areas)
% title('Fraction of units in detected joint assemblies')
ylabel('Fraction of units in assembly type (%)');
%% Plot fraction of neurons acting as local nodes

m_ = nanmean(fracBridgeNeurons,2);
e_ = nansem(fracBridgeNeurons,2);
figure; hold on
B = bar([1,2],m_*100);
B.FaceColor = [1 1 1];
B.EdgeColor = [0 0 0];
B.LineWidth = 1.5;

EB = errorbar([1,2],m_*100,e_*100,'.b');
EB.Color = [0 0 0];
EB.LineWidth = 1.5;
ylim([0 50])
set(gca,'XTick',[1,2],'XTickLabel',Areas)
% title('Fraction of units bridging assemblies')
ylabel('Fraction of units %');
%% Plot fraction of neurons acting as inter-area bridges

% Ignore datasets with no detected assemblies
fracLocalNodes(isinf(fracLocalNodes))=NaN; 

m_ = nanmean(fracLocalNodes,2);
e_ = nansem(fracLocalNodes,2);
figure; hold on
B = bar([1,2],m_*100);
B.FaceColor = [1 1 1];
B.EdgeColor = [0 0 0];
B.LineWidth = 1.5;

EB = errorbar([1,2],m_*100,e_*100,'.b');
EB.Color = [0 0 0];
EB.LineWidth = 1.5;
ylim([0 50])
set(gca,'XTick',[1,2],'XTickLabel',Areas)
% title('Fraction of units bridging assemblies')
ylabel('Fraction of units %');
%%


m_ =  [nanmean(fracAssUnits,2);nanmean(fracJointNeurons,2);nanmean(fracBridgeNeurons,2);nanmean(fracLocalNodes,2);];
e_ =  [nansem(fracAssUnits,2);nansem(fracJointNeurons,2);nansem(fracBridgeNeurons,2);nansem(fracLocalNodes,2);];
figure; hold on
B = bar([1:8],m_*100);
B.FaceColor = [1 1 1];
B.EdgeColor = [0 0 0];
B.LineWidth = 1.5;

EB = errorbar([1,2],m_*100,e_*100,'.b');
EB.Color = [0 0 0];
EB.LineWidth = 1.5;
ylim([0 50])

set(gca,'XTick',[1,2],'XTickLabel',Areas)
% title('Fraction of units in detected joint assemblies')
ylabel('%')
clear
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
else
    pat = '/Volumes/Akasa/Bristol/AssemblyAnalysis/raw/';
end

Target = 'LONG';
cd(pat);
fileList = dir([pat 'allTimestamps' filesep '*' Target '*.mat']);

reject_list={'MiroslawLONG1_Events.mat','IreneuszLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idxJoint=1:numel(fileList)
    fnames_{idxJoint,1}=fileList(idxJoint).name;
    name_flag(idxJoint,1)= logical(sum(ismember(reject_list,fileList(idxJoint).name))); %AD inverted logical flag for testing
end
if sum(name_flag)>0
    fileList(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag
useWholeTaskPeriod = true;
resampleSpikes = true;
noDrawnSpikes  = 400;
nBS = 500;
Delays_ = {'Short','Medium','Long'};

%% Batch process
TimeSpan = 4;
threshChoice = 4; % 2=p<0.05 , 3=p<0.01
warning('off')
Areas = {'PFC','HP','joint'};
binsize = 1;
viLag = (-1000:binsize:1000)/1000;

len = length(viLag);
nfft = 2^nextpow2(2*len-1);
Fs = (1/binsize)*1000;
f = Fs*(0:(nfft/2))/nfft;
win_ = hamming(0.1*Fs);

plotOnlineYN = false;
plotRawOnlineYN = false;
plotAsswideLockingOnlineYN = false;

for iDay = 1:length(fileList)
    %% Get spike times
    fname = strtok(fileList(iDay).name,'_');
    fn = fullfile(pat,sprintf('%s.mat',fname ));
    spikes = load(fn);
    noPFC = length(spikes.PFCcells);
    %% Get assembly information
    
    if useWholeTaskPeriod
        fn = fullfile(pat,'KDE_binsTaskOnly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_FSC.mat',fname));
        A  = load(fn);
        fn = fullfile(pat,'KDE_binsTaskOnly','LONGTaskonly',sprintf('%s_iFR50_behavOnly_AssemRes2.mat',fname));
        B  = load(fn);
        % recalculate included units
%         fn = fullfile(pat,'KDE_binsTaskOnly',sprintf('%s_PFC_iFR50_behavOnly.mat',fname));
%         B.usel_out=SelCells(fn,0.1,1e6);
    else
        fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',fname));
        A  = load(fn);
        fn = fullfile(pat,'KDE_bins','LONG',sprintf('%s_iFR50_FSC.mat',fname));
        B  = load(fn);
        % recalculate included units
        fn = fullfile(pat,'KDE_bins',sprintf('%s_PFC_iFR50.mat',fname));
        [~,~,~,~,B.usel_out]=SelTrialsCellsWholeTrial(fn,10,0.1,1e8,'iFR');
    end
    
    
    spikes.jointcells=[spikes.PFCcells;spikes.HPcells]';
    A.nu(3) =sum(A.nu(1:2));
%     B.usel_out{3}= [B.usel_out{1},B.usel_out{2}+max(B.usel_out{1})];
    B.usel_out{3}= [B.usel_out{1},B.usel_out{2}+noPFC];

   %% Sort units by assembly membership types
   
    if ~isempty(A.FSCsel{3})
        
        nAssJoint = length(A.units{3});

        % 1) Check that inter-area assems actually span the two areas
        
        idxJoint = false(nAssJoint,1);
        for iJointAss = 1:nAssJoint
            U_ = B.usel_out{3}(A.units{3}{iJointAss});
            
            if  min(U_)>max(B.usel_out{1})
                idxJoint(iJointAss) = true;
            end
            A.units{3}(idxJoint)=[];
            A.FSCsel{3}(:,idxJoint)=[];
        end
        
        
        % 2) Remove local assems that are subsets of larger inter-area assems
        
        for iArea = 1:2
            if ~isempty(A.FSCsel{iArea})
                
                nAssLocal = length(A.units{iArea});
                idxLocal  = false(nAssLocal,1);
                
                for iJointAss = 1:nAssJoint
                    Jointunits = A.units{3}{iJointAss};
                    
                    for iLocalAss = 1:nAssLocal
                        Localunits = A.units{iArea}{iLocalAss};
                        
                        if iArea==2
                            Localunits = Localunits + A.nu(1);
                        end
                        
                        if length(Localunits)==length(intersect(Localunits,Jointunits))
                            idxLocal(iLocalAss)  = true;
                        end
                        
                    end
                end
                
                A.units{iArea}(idxLocal)=[];
                A.FSCsel{iArea}(:,idxLocal)=[];
                
            end
        end
        
    end   
    
    
    
    % 3) Prune out inter-area members from local members
    for iArea = 1:3
        Ass.LocalMembers{iArea}    = unique(cell2mat(A.units{iArea}));
    end
    Ass.LocalMembers{1} = setdiff(Ass.LocalMembers{1}, Ass.LocalMembers{3}(Ass.LocalMembers{3} <= A.nu(1)));
    Ass.LocalMembers{2} = setdiff(Ass.LocalMembers{2}, Ass.LocalMembers{3}(Ass.LocalMembers{3} > A.nu(1)) - A.nu(1));
    
    % 4) Get inter-area members and nonmembers
    for iArea = 1:2
        Ass.JointMembers{iArea} = setdiff(unique(cell2mat(A.units{iArea})),Ass.LocalMembers{iArea});
        Ass.NonMembers{iArea}   = setdiff(1:A.nu(iArea),[Ass.LocalMembers{iArea},Ass.JointMembers{iArea}]);
    end
    for iArea = 1:2
        Ass.NonMembers{iArea}   = B.usel_out{iArea}(Ass.NonMembers{iArea});
        Ass.LocalMembers{iArea} = B.usel_out{iArea}(Ass.LocalMembers{iArea});
        Ass.JointMembers{iArea} = B.usel_out{iArea}(Ass.JointMembers{iArea});
    end
 
    % 4) Create [Units x Assembly] assignment matrices
     for iArea = 1:2
        % Pad this out to ensure that there's an entry even if there's no detected assembly
        LocalMembersMatrix_{iArea} = nan(A.nu(iArea),length(A.units{iArea})+1);
        if ~isempty(A.units{iArea})
            for iAss=1:length(A.units{iArea})
                LocalMembersMatrix_{iArea}(:,iAss+1)=ismember(1:A.nu(iArea),A.units{iArea}{iAss});
            end
        end
        
        JointMembersMatrix_{iArea} = nan(A.nu(iArea),length(A.units{3})+1);
        for iAss=1:length(A.units{3})
            if ~isempty(A.units{3})
                units_ = A.units{3}{iAss};
                if iArea==1
                    JointMembersMatrix_{iArea}(:,iAss+1)=ismember(1:A.nu(iArea), units_(units_<=A.nu(1)));
                else
                    JointMembersMatrix_{iArea}(:,iAss+1)=ismember(1:A.nu(iArea), units_(units_>A.nu(1))-A.nu(1));
                end
            end
            
        end
     end
    
    iArea = 3;
    JointMembersMatrix_{iArea} = nan(A.nu(iArea),length(A.units{3})+1);
    for iAss=1:length(A.units{3})
        if ~isempty(A.units{3})
            units_ = A.units{3}{iAss};
            JointMembersMatrix_{iArea}(:,iAss+1)=ismember(1:A.nu(iArea), units_);
        end
    end
    
    
    for iArea = 1:2
        D.fname{iDay} = fname;
        D.nonmemberUnits{iDay}{iArea}   = Ass.NonMembers{iArea};
        D.localmemberUnits{iDay}{iArea} = Ass.LocalMembers{iArea};
        D.jointmemberUnits{iDay}{iArea} = Ass.JointMembers{iArea};
    end

    D.LocalMembersMatrix{iDay} = LocalMembersMatrix_;
    D.JointMembersMatrix{iDay} = JointMembersMatrix_;
    D.usel_out{iDay} = B.usel_out;
    D.units{iDay} = A.units;
    D.FSCsel{iDay} = A.FSCsel;
end
%% Save
% fname = 'C:\Analysis\AssemblyAnalysis\raw\KDE_binsTaskonly\LONGTaskonly\Assemmembertypes.mat';
fname = fullfile(pat,'KDE_binsTaskonly','LONGTaskonly','Assemmembertypes.mat');


save(fname,'D')

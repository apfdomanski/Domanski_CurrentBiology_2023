%% %%%%%% PREAMBLE %%%%%%
clear

Target = 'LONG';

warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw\';
    addpath(genpath('C:\Users\ad15419\Dropbox\MATLAB\PC MATLAB Path\AssemblyCode\Eleonora\Mic\programs'))
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    addpath(genpath('/panfs/panasas01/phph/ad15419/MATLAB/Mic'))
    
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];
end

fileList = dir([pat '*' Target '*.mat']);


% reject_list={'MiroslawLONG1_Events.mat', 'NorbertLONG2_Events.mat', 'OnufryLONG1_Events.mat' , 'OnufryLONG2_Events.mat' }; %'ALL_events.mat'
reject_list={};%{'MiroslawLONG1_Events.mat'}; %'ALL_events.mat'
name_flag=zeros(numel(fileList),1);
for idx=1:numel(fileList)
    fnames_{idx,1}=fileList(idx).name;
    name_flag(idx,1)= logical(sum(ismember(reject_list,fileList(idx).name))); %AD inverted logical flag for testing
end
try
    fileListAss(find(name_flag))=[];% fnames_(name_flag)=[];
    fileList(find(name_flag))=[];   % fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag

Areas = {'PFC','HP','Joint'};

mkdir([pat 'MICresults'])


BinSizes=[0.005 0.01 0.0150    0.0200    0.0300    0.0500    0.0800    0.1200    0.2000    0.3500    0.5000   0.7000    1.0000];
MaxLags=[10   10   10   10   10    10    10    10    10     10     10     10     10];
alph=0.05; %Sig level
Dc=100;    % Bin length for variance #abba calculation
No_th=10;  % Minimal no. activations
O_th = Inf;   % Max assembly order
bytelimit = 200e9; %Memory limit

display='ordunit'; % 'clustered' 'ordunit' 'raw'
criteria = 'distance'; % 'biggest' or 'distance';
act_count = 'combined'; % 'full' 'partial' 'combined'
lagChoice = 'duration'; % 'duration' 'beginning'
%% Batch process units
for iFile =1:length(fileList)
    % Get the spike times
    fname=fileList(iFile).name;
    load(sprintf('%s%s%s',pat,filesep,fname),'PFCcells','HPcells');
    
    % Wrangle data
    for i =1:length(PFCcells)
        ST{1}{i} =  PFCcells{i}.t*1e-4;
        N_{1}(i) = length(ST{1}{i}); % No spikes
        M_{1}(i) = max(ST{1}{i});    % last spike time
    end
    for i =1:length(HPcells)
        ST{2}{i} =  HPcells{i}.t*1e-4;
        N_{2}(i)  = length(ST{2}{i}); % No spikes
        M_{2}(i)  = max(ST{2}{i});    % last spike time
    end
    
    ST{3} = [ST{1},ST{2}];
    N_{3} = [N_{1},N_{2}];
    M_{3} = [M_{1},M_{2}];
    
    epoch = [0, max(cell2mat(M_))];
    
    for s =1:3
        try
            mkdir([pat 'temp_' strtok(fname,'.') '_' Areas{s}])
            
            cd([pat 'temp_' strtok(fname,'.') '_' Areas{s}])
            noNeurons = length(ST{s});
            longest_  = max(N_{s});
            
            spM{s}=nan(noNeurons,longest_);
            
            for iUnit = 1:noNeurons
                temp = ST{s}{iUnit};
                temp(temp<epoch(1))=[];
                temp(temp>epoch(2))=[];
                spM{s}(iUnit,1:length(temp))  = temp;
            end
            spM{s}(:,isnan(nansum(spM{s},1)))=[];
            
            % Find assemblies
            [assemblies]=Main_assemblies_detection(spM{s}, MaxLags, BinSizes,alph,Dc,No_th,O_th,bytelimit);
            
            fname_ = [pat 'MICresults' filesep strtok(fname,'.') '_' Areas{s} '_MICresults.mat'];
            save(fname_,'assemblies','epoch','noNeurons','-v7.3')
            
            % Visualise
            [As_across_bins,As_across_bins_index]=assemblies_across_bins(assemblies, BinSizes);
            [Amatrix,Binvector,Unit_order, As_order]=assembly_assignment_matrix(As_across_bins, noNeurons, BinSizes, display);
            
            % Pruning
            if strcmp(criteria,'biggest')
                [As_across_bins_pr, As_across_bins_index_pr]=pruning_across_bins(As_across_bins,As_across_bins_index,noNeurons,criteria);
            elseif strcmp(criteria,'distance')
                th=0.2;
                style = 'pvalue'; % 'pvalue' or 'occ'
                [As_across_bins_pr, As_across_bins_index_pr]=pruning_across_bins(As_across_bins,As_across_bins_index,noNeurons,criteria,th,style);
            end
            [Amatrix,Binvector,Unit_order, As_order]=assembly_assignment_matrix(As_across_bins_pr, noNeurons, BinSizes, display);
            
            % assembly activity
            assembly_activity = assembly_activity_function(As_across_bins_pr, assemblies, spM{s}, BinSizes, lagChoice, act_count);
            
            % Save
            fname_ = [pat 'MICresults' filesep strtok(fname,'.') '_' Areas{s} '_Assemblies.mat'];
            save(fname_,'assembly_activity','assemblies','As_across_bins','As_across_bins_index','Amatrix','Binvector','Unit_order','As_order','As_across_bins_pr','As_across_bins_index_pr','-v7.3')
            
            % tidy up
            clear assemblies As_across_bins As_across_bins_index Amatrix Binvector Unit_order As_order As_across_bins_pr As_across_bins_index_pr
        catch ME
            disp(['failed on ' strtok(fname,'.') '_' Areas{s} ': ' ME.message])
            
        end
    end
    
end




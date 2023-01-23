% Add information on erroneous behavioural outcomes for import scripts if missing
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
% 
clear all
close all
target= 'physostigmine';

if ispc
    home = 'C:\'; % [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    pat = [home 'Analysis\AssemblyAnalysis\raw\'];
elseif ismac
    path(path,'/Users/aleksanderdomanski/Documents/Bristol_work/assembly_analysis/MichalDataAna')
    home = getenv('HOME');
    cd ([home '/Documents/Bristol_work/assembly_analysis/MichalDataAna'])
    pat = [home '/Documents/Bristol_work/assembly_analysis/Michal_data/raw/'];
else
    path(path,'/panfs/panasas01/phph/ad15419/MATLAB')
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw/'];    
end

mkdir([pat '/' target])
bw=0.05; iFRtyp='iFR';
reject_list={'ALL_events.mat'}; %'ALL_events.mat'
fnames   = dir([pat '*.MAT']);
name_flag=zeros(numel(fnames),1);
for idx=1:numel(fnames)
fnames_{idx,1}=fnames(idx).name;
name_flag(idx,1)= logical(sum(ismember(reject_list,fnames(idx).name))); %AD inverted logical flag for testing
end
try
    fnames(find(name_flag))=[];% fnames_(name_flag)=[];
end
clear  reject_list idx fnames_ name_flag
%% add t-info on error trials
for file_id=1:length(fnames)
    load([pat fnames(file_id).name])
    ext = strtok(fnames(file_id).name, '.');
    bw=50e-3;
    Vstr={'left_choice','right_choice','left_sample','right_sample'};
   
    
    brL={'PFC','HP'}; 
    EvtTinc=[]; EvtLinc=[];
    A=who('*HPcells'); br='HP';
    k=1; i=1;
    k=findstr('_',A{i});
    for j=1:length(Vstr) % loop across trial types
        f=['ERRORtrange' Vstr{j}];
        t=eval(f); EvtTinc=[EvtTinc;(t(:,1)*1e-6+5)];
        EvtLinc(end+1:end+length(t))=j*10;
    end;
    for b=1:2            % loop across brain regions
        br=brL{b};
        try
        fn=[pat 'KDE_bins/' A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3)];
        save(fn,'EvtTinc','EvtLinc','-append');
        catch
            disp([A{i}(1:k) ext '_' br '_iFR' num2str(bw*1e3) ' missing!'])
        end
    end;
end;
disp done!

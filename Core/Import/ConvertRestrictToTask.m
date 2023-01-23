clear 
Target = 'CP55940';
if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
    cd(pat)
else
    home = getenv('HOME');
    cd ([home '/MATLAB/MichalDataAna'])
    pat = [home '/raw'];
end

bw=0.05;        
offset = [-60 60];
Areas = {'HP','PFC'};

warning ('off')
f = [pat filesep 'allTimestamps' filesep '*' Target '*' '.mat']
fileList=dir(f);fileList

%% Instantiate parallel pool
delete(gcp('nocreate'))
parpool local
poolobj = gcp;
addAttachedFiles(poolobj,{'UniKDE.m','gfilter.mexa64'})
%%

for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    try
    load(sprintf('%s%sallTimestamps%s%s_Events.mat',pat,filesep,filesep,fname),'t');
    load(sprintf('%s%s%s.mat',pat,filesep,fname),'PFCcells','HPcells');
    
    t_ = [ collapseStructure(t.Short);collapseStructure(t.Medium);collapseStructure(t.Long)]*1e-6;
    tlims_ = [min(t_)+offset(1), max(t_)+offset(2)];
    %%
    for iArea =1:length(Areas)
        
        fprintf('Analysing run %d/%d %s (%s)...',iFile,length(fileList),fname,Areas{iArea})
        
        % Timebase confined to behavioural timestamps +/- a bit of buffer:
        Tmtx=tlims_(1):bw:tlims_(2);
        
        eval(sprintf('B = %scells;',Areas{iArea}));

        ST=[]; avgFR=[];
        hucv=[]; 
        MISE=[]; 
        CvBW=[]; 
        iFR0=cell(1,length(B));
        hucvS=cell(1,length(B));
        
        parfor u=1:length(B)
            [u]
            % collect ST
            ST = restrictFRrange( B{u}.t*1e-4 , tlims_ );
            avgFR(u)=length(ST)/range(ST);
            % convert STMtx into iFR via KDE
            h0=mean(diff(ST))^2;
            %[hucv(u),MISE(u)]=UniKDE(ST,h0,3);
            
            [hucv(u),MISE(u),CvBW(u)]=UniKDE_AD(ST,h0,3,[],100,12,0);
            if ~isnan(hucv(u))
                iFR0{u}=STtoGfAvg(ST',Tmtx,sqrt(hucv(u)))';
%                 iFR0{u} = KDens(Tmtx(1:end-1)',ST,(hucv(u)))';
            else
                iFR0{u}=zeros(length(Tmtx)-1,1);
            end
            
            % Sanity check firing rates
%             figure; hold on
%             scatter(ST,zeros(size(ST)),'.b')
%             plot(Tmtx(1:end-1),mat2gray(iFR0{u}),'b')
%             stairs(Tmtx,mat2gray(histc(ST,Tmtx)),'r')
        end
   
        iFR=cell2mat(iFR0);

        save([pat filesep 'KDE_binsTaskonly' filesep fname '_' Areas{iArea} '_iFR' num2str(bw*1e3) '_behavOnly'],...
                'iFR',...  % KDE optimised firing rates (size; time vector,no.units)
                'Tmtx',... % time vector
                'avgFR',...% each neuron's average firing rate
                'CvBW',... % variance in optimal kernel width over repeated trials
                'hucv',... % each neuron's optimised kernel width
                'MISE'); % mean-squared-error for each neuron's optimal kernel width
    
    end
    catch 
        disp(['no unit data for ' fname '!'])
    end
end


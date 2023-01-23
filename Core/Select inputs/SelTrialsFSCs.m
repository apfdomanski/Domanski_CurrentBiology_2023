function [TmtxS,FSCs,EvtLs,EvtTs]=SelTrialsFSCs(fn,twin,trials)
% Select Assemblies from PFC/HC/Joint PFC~HP data sets & cut out event periods (trial case)
%
% Aleksander PF Domanski PhD UoB 2015
% aleks.domanski@bristol.ac.uk
%
% Cuts up time series into equal length epochs for use in downstream bootstrapping

if nargin<2 || isempty(twin),     twin=10;       end; % time window (s)
if nargin<3 || isempty(trials),   trials='corr'; end; % include errors? ("corr" N or "all" Y)

%outputs
% Tmtxs {no. trials}                                         time axes for selected trials
% iFRs  {no. brain areas}{no. trials}(timerange,no. neurons) spike FR matrix
% EvtLs (no. events)                                         event labels       
% EvtTs (no. events)                                         event times
load(fn)
hFig = figure('name','Trial import progress');
T=cell(1,2); FSCs=T; TmtxS=T;
for f=1:length(FSCsel) % across assem types (regions)
    
    % NB including error trials not yet implemented...
    %     if strncmp(trials,'all',3)
    %         k=find(EvtTinc>min(EvtT) & EvtTinc<max(EvtT));
    %         EvtT=[EvtT;EvtTinc(k)]; 
    %         EvtL=[EvtL,EvtLinc(k)];
    %     end;
    
    usel=1:size(FSCsel{f},2);
    
%     Tmtx=Tmtx(1:end-1);

    T{f}=Tmtx; FSCtemp=FSCsel{f}(:,usel);
    
    [EvtTs,s]=sort(EvtTs); % event times
    EvtLs=EvtLs(s);        % event types
    % {'left_choice','right_choice','left_sample','right_sample'};
    FSCs{f}=cell(1,round(length(EvtTs)/2)); TmtxS{f}=FSCs{f};
    for i=1:2:length(EvtTs)-1
        k=find(Tmtx>=EvtTs(i)-twin & Tmtx<=EvtTs(i+1)+twin); % trial cutout indices
%         if max(k)>max(Tmtx) break; end
        FSCs{f} {ceil(i/2)} = FSCtemp(k,:);
        TmtxS{f}{ceil(i/2)} = Tmtx(k);
        
        clf(hFig); hold on
        scatter(EvtTs,EvtLs); 
        rectangle('Position',[min( Tmtx(k)),0,max( Tmtx(k))-min( Tmtx(k)),40]);
        title(['Importing trial: ', num2str(i)])
        drawnow
    end;
    
end;
vstr= {'L Choice','R Choice','L Sample','R Sample'};
figure; hold on
for i=1:length(TmtxS{f})
plot(TmtxS{f}{i},i+0*TmtxS{f}{i})
end
for i=1:length(EvtTs)
    text(EvtTs(i),round(i/2),vstr{EvtLs(i)},'Rotation',90)
end
xlabel('Time')
ylabel('Trial')
if ~isempty(find(T{1}~=T{2},1)), warning('T1 and T2 do not match!'); end;
close(hFig)
% (c) 2013 Daniel Durstewitz, Bernstein Center for Computational
% Neuroscience, CIMH/ Heidelberg University

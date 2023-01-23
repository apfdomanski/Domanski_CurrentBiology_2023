%% %%%%%% PREAMBLE %%%%%%
clear
Delays_ = {'Short','Medium','Long'};
Delays__ = {'4s','8s','16s'};
Target = 'LONG';

shift = 0;
plotOnline = false;
bw=0.05;
Nbs = 500;

warning ('off')
if ispc
    pat = 'C:\Analysis\AssemblyAnalysis\raw';
else
    pat = '/Volumes/HDD2/DNMTP/raw/';
end
cd(pat)
fileList=dir(sprintf('MixedSelectivity_LONG%s*%s*_MixedSelectivity_MembershipsortedUnits.mat',filesep,Target));
normWin = [-5 -3];
Areas = {'PFC','HP','Joint'};
color_={'b','r','g'};
MemberClasses ={'LocalMembers','JointMembers','NonMembers'};
col_ ={[0.6 0.6 0.6],[0.9 0.6 0],[0.6 0 0.6]};
normWin = [0 2];
normaliseFscores = true;

%% import correct
clear D_ D;
for iFile =1:length(fileList)
    
    fname=strtok(fileList(iFile).name,'_');
    fnIn = sprintf('%sMixedSelectivity_LONG%s%s_MixedSelectivity_MembershipsortedUnits.mat',pat,filesep,fname);
    load(fnIn ,'D');
    
    for iDelay =1:length(Delays_)
        for iClass= 1:length(MemberClasses)
            
            for iArea =1:2  
                eval(sprintf('D_.%s.LR.%s{iArea}.Ft2{iFile}= D.%s.LR.%s{iArea}.Ft2;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                eval(sprintf('D_.%s.LR.%s{iArea}.Rt2{iFile}= D.%s.LR.%s{iArea}.Rt2;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                eval(sprintf('D_.%s.LR.%s{iArea}.CVE{iFile}= D.%s.LR.%s{iArea}.CVE;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                eval(sprintf('D_.%s.LR.%s{iArea}.CVEbs{iFile}= D.%s.LR.%s{iArea}.CVEbs;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))

            end
        end
    end
end
%% import errors
clear D;
for iFile =1:length(fileList)
        fname=strtok(fileList(iFile).name,'_');
    try
    fnIn = sprintf('%sMixedSelectivity_LONG%s%s_MixedSelectivity_MembershipsortedUnits_Err.mat',pat,filesep,fname);
    load(fnIn ,'D');
    end
    for iDelay =1:length(Delays_)
        for iClass= 1:length(MemberClasses)
            
            for iArea =1:2
       
                    try
                eval(sprintf('D_.%s.LR_err.%s{iArea}.Ft2{iFile}= D.%s.LR_err.%s{iArea}.Ft2;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                eval(sprintf('D_.%s.LR_err.%s{iArea}.Rt2{iFile}= D.%s.LR_err.%s{iArea}.Rt2;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                eval(sprintf('D_.%s.LR_err.%s{iArea}.CVE{iFile}= D.%s.LR_err.%s{iArea}.CVE;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                eval(sprintf('D_.%s.LR_err.%s{iArea}.CVEbs{iFile}= D.%s.LR_err.%s{iArea}.CVEbs;',Delays_{iDelay},MemberClasses{iClass},Delays_{iDelay},MemberClasses{iClass}))
                disp(fname)
                    end
            end
        end
    
    end
    
end

%% Plot correct
figure;
clear y y_


for iArea = 1:2
    subplot(2,1,iArea); hold on
    
    for iClass = 1:length(MemberClasses)
        for iDelay = 1:length(Delays_)
            try
                eval(sprintf('y_ = D_.%s.LR.%s{iArea}.Ft2;',Delays_{iDelay},MemberClasses{iClass}))
                y_ = cell2mat(y_');
            catch 
                y_ = zeros(12,402);
            end
            x = (1:size(y_,2))*bw;
            if normaliseFscores
                bp = find(x>normWin(1) & x<normWin(2));
                B  = nanmean(y_(:,bp)');
%                 y_ = y_ - B'*ones(1,length(x));
%                 y_ = y_./(B'*ones(1,length(x)));
            end
            y{iDelay} = y_;
            clear y_
        end
        
        
        y_ = zeros(size(y{1}));
        for iDelay = 1:length(Delays_)
            y_=y_+y{iDelay};
        end
        y_ = y_./length(Delays_);
        
        ciplot(nanmean(y_)+nansem(y_),nanmean(y_)-nansem(y_),...
            x, col_{iClass},0.5);
        
    end
end
%% Plot errors
figure;
clear y y_


for iArea = 1:2
    subplot(2,1,iArea); hold on
    
    for iClass = 1:length(MemberClasses)
        for iDelay = 1:length(Delays_)
            try
                eval(sprintf('y_ = D_.%s.LR_err.%s{iArea}.Ft2;',Delays_{iDelay},MemberClasses{iClass}))
                y_ = cell2mat(y_');
            catch 
                y_ = zeros(1,402);
            end
            x = (1:size(y_,2))*bw;
            if normaliseFscores
                bp = find(x>normWin(1) & x<normWin(2));
                B  = nanmean(y_(:,bp)');
%                 y_ = y_ - B'*ones(1,length(x));
%                 y_ = y_./(B'*ones(1,length(x)));
            end
            y{iDelay} = y_;
            clear y_
        end
        
        
        y_ = zeros(size(y{1}));
        for iDelay = 1:length(Delays_)
            y_=y_+y{iDelay};
        end
        y_ = y_./length(Delays_);
        
        ciplot(nanmean(y_)+nansem(y_),nanmean(y_)-nansem(y_),...
            x, col_{iClass},0.5);
        
    end
end
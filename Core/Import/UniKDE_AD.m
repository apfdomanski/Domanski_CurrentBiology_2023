function [h,MISE,CvBW,F]=UniKDE_AD(X,h0,Ef,Z,ub,noBlocks,usePSO)
% Fit spike train rate kernel: Univariate nonparametric density estimation with bandwidth selection
%
% Input arguements
% X  - spike times
% h0 - inital guess for kernel density width
% Ef - choice of error function (1:3)
% Z  - Optional time vector for production of smoothed firing rate
% ub - upper bound to constrain fitting
% noBlocks - no. chunks to process spikes over
% usePSO - optional flag to use particle swarm instead of fmincon
% 
% Output arguements
% h     - Optimised kernel sigma
% MISE  - Mean integrated Squared error on firnal iteration
% F     - Optional spike density estimate projected onto Z
% CvBW  - variance in kernel estimate over blocks
% 
% This routine implements univariate nonparametric density estimation
% with bandwidth selection based on minimizing the mean integrated
% squared error (MISE), E[int((f^(x,H)-f(x))^2)], estimated through
% unbiased cross-validation (Ef=1) or bootstrap (Ef=2); based on Taylor
% (1989), Biometrika 76, pp 705-12.
%
% (c) 2011 Daniel Durstewitz, CIMH & BCCN Heidelberg-Mannheim, Germany
% AD edited Dec 2015 to adjust 'optimoptions' for MATLAB 2015a


if nargin<3, Ef=1; end;
if Ef==1, [h,MISE]=fminsearch(@UCVErr1,h0,[],X); end;   % unbiased cross-validation
if Ef==2, [h,MISE]=fminsearch(@BSErr,h0,[],X); end;   % bootstrap
if Ef==3
    if nargin<5, ub=100; end;
    
    t1 = tic;
    if usePSO
        options = optimoptions('particleswarm','Display','iter','SwarmSize',10,'PlotFcns',@pswplotbestf,'UseParallel',false);
    else
        %options=optimoptions('fmincon','Display','iter','GradObj','on','PlotFcns',@optimplotfval);
%         options=optimoptions('fmincon','Display','iter','GradObj','on');
        options=optimoptions('fmincon','GradObj','on');
    end
    
    
    % Run on all spike times first
    disp('Fitting all spikes')

    if usePSO
        [h,MISE,exitFlag,output] = particleswarm(@(h0)UCVErr2(h0,X),1,0.1*h0,10*h0,options);
    else
        [h,MISE,exitflag,output] = fmincon(@(h)UCVErr2(h,X),h0,[],[],[],[],0.1*h0,10*h0,[],options);
    end
    
    
    for iBlock = 1:noBlocks
        t2 = tic;
        times = X;
        times_= prctile(times',[iBlock-1,iBlock]*100/noBlocks);
        times(times<=floor(times_(1)) | times>ceil(times_(2)))=[];
        h0(iBlock) = mean(diff(times))^2;
        try
        if usePSO
            [h_(iBlock),MISE_(iBlock),exitFlag,output] = particleswarm(@(h)UCVErr2(h,times),1,0.1*h0(iBlock),10*h0(iBlock),options);
        else
            [h_(iBlock),MISE_(iBlock),exitflag,output]=fmincon(@(h)UCVErr2(h,times),h0(iBlock),[],[],[],[],0.1*h0(iBlock),10*h0(iBlock),[],options);
        end
        catch
            h_(iBlock) = NaN;
            MISE_(iBlock) = NaN;
            output.iterations = 0;
        end
        
        elaps(iBlock) = toc(t2);
        iterations(iBlock) = output.iterations;
       
%         disp(sprintf('Fitting result: Block %d of %d // %d iterations // MISE = %f // optimum kernel = %d(took %.2fs)\n',iBlock,noBlocks,output.iterations,MISE(iBlock),h(iBlock),toc(t2)))
        %             Z = min(times):0.05:max(times);
        %             F=KDens(Z,times,h(iBlock));
        %             H = histc(times,Z);
        %             close all
        %             figure('Position',[100,400,1500,400] ,'name',sprintf('Fitting result: Block %d of %d // %d iterations // MISE = %f // optimum kernel = %d(took %.2fs)',iBlock,noBlocks,output.iterations,MISE(iBlock),h(iBlock),toc(t2))); hold on
        %             area(Z,mat2gray(F));stairs(Z,mat2gray(H),'r');
        %             pause(0.01);drawnow
        
    end
    
%     figure('Position',[100,400,1500,400] ,'name',sprintf('Fitting result: Block %d of %d // %d iterations // MISE = %f // optimum kernel = %d(took %.2fs)',iBlock,noBlocks,output.iterations,MISE(iBlock),h(iBlock),toc(t2))); hold on
    CvBW=std(h_)/mean(h_);
    
%     h = mean(h);
%     MISE = mean(MISE);
    
       


% evaluate density at grid points specified in Z
if ~isempty(Z) | nargout>3
    %F=KDens(Z,X,mean(h));
    idx =1:200;
    clear F
    noBlocks = ceil(length(Z)/length(idx));
    for iBlock = 1:noBlocks
        offset = (iBlock-1)*length(idx);
        idx_ = idx+offset;
        try
            X_  =X(X>=Z(idx_(1)) & X<=Z(idx_(end)));
            F{iBlock}=KDens(Z(idx_),X_,h);
        catch
            X_  =X(X>=Z(idx_(1)));
            F{iBlock}=KDens(Z(idx_(1):end),X_,h);
        end
    end
    F = cell2mat(F);
end
end
%%
%     H = mat2gray(histc(X,Z));
%     idx =1:200;
%     noBlocks = ceil(length(Z)/length(idx));
% range_ = 0.1:0.1:5;
% for i=1:length(range_)
% i     
%     clear F
%     
%     for iBlock = 1:noBlocks
%         
%         offset = (iBlock-1)*length(idx);
%         idx_ = idx+offset;
%         try
%             X_  =X(X>=Z(idx_(1)) & X<=Z(idx_(end)));
%             F{iBlock}=KDens(Z(idx_),X_,mean(h)*range_(i));
%         catch
%             X_  =X(X>=Z(idx_(1)));
%             F{iBlock}=KDens(Z(idx_(1):end),X_,mean(h)*range_(i));
%         end
%     end
%     F = cell2mat(F);
%     A(i,:) = [range_(i),mean(mat2gray(F)-H')];
%     
% end
% 
% figure; plot(A(:,1)*mean(h),A(:,2))
%%
function [err,derr,dderr]=UCVErr2(h,X)
N=length(X);
k0=1/(2*N^2*sqrt(2*pi));
k1=sqrt(2);
k2=4*N/(N-1);
k3=4*N^2/(N-1);
B=0; C=0;
dB=0; dC=0;
ddB=0; ddC=0;
for i=1:N-1
    for j=i+1:N
        z=(X(i)-X(j))^2;
        y1=exp(-z/(4*h));
        y2=exp(-z/(2*h));
        B=B+y1;
        C=C+y2;
        dB=dB+z/(4*h^2)*y1;
        dC=dC+z/(2*h^2)*y2;
        ddB=ddB+(z^2/(16*h^4)-z/(2*h^3))*y1;
        ddC=ddC+(z^2/(4*h^4)-z/h^3)*y2;
    end;
end;
B=2*B+N; C=2*C+N; A=1/sqrt(h);
dB=2*dB; dC=2*dC; dA=-1/2*h^(-3/2);
ddB=2*ddB; ddC=2*ddC; ddA=3/4*h^(-5/2);
err=k0*A*(k1*B-k2*C+k3);
derr=k0*dA*(k1*B-k2*C+k3)+k0*A*(k1*dB-k2*dC);
dderr=k0*ddA*(k1*B-k2*C+k3)+k0*dA*(k1*dB-k2*dC)+ ...
    k0*dA*(k1*dB-k2*dC)+k0*A*(k1*ddB-k2*ddC);


function err=UCVErr1(h,X)
N=length(X);
a1=1/(N^2*sqrt(4*h*pi));
a2=2/(N*(N-1)*sqrt(2*h*pi));
a3=2/((N-1)*sqrt(2*h*pi));
err1=0; err2=0;
for i=1:N-1
    for j=i+1:N
        b=(X(i)-X(j))^2;
        err1=err1+exp(-b/(4*h));
        err2=err2+exp(-b/(2*h));
    end;
end;
err=a1*(N+2*err1)-a2*(N+2*err2)+a3;


function err=BSErr(h,X)
N=length(X);
a=1/(2*N^2*sqrt(2*h*pi));
b=N*sqrt(2)/2;
err1=0; err2=0; err3=0;
for i=1:N-1
    for j=i+1:N   % exclude i==j
        z=(X(i)-X(j))^2;
        err1=err1+exp(-z/(8*h));
        err2=err2+exp(-z/(6*h));
        err3=err3+exp(-z/(4*h));
    end;
end;
%err=2*a*(err1-4/sqrt(3)*err2+sqrt(2)*err3+b);   % something wrong here!!!###
err=2*a*(err1-4/sqrt(3)*err2-sqrt(2)*err3+b);   % must be wrong as well ...
% There must be a bug in the original eqns. given in Taylor (1989)


function F=KDens(Z,X,h)
N=length(Z);
F=zeros(1,N);
a=1/sqrt(2*pi*h);
for i=1:N
    for j=1:length(X)
        dt=Z(i)-X(j);
        F(i)=F(i)+exp(-dt^2/(2*h));
    end;
end;
F=a*F./length(X);


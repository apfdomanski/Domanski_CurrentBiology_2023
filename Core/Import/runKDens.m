function F = runKDens(Z,X,h)

    idx =1:1000;
    clear F
    noBlocks = ceil(length(Z)/length(idx));
    for iBlock = 1:noBlocks
        offset = (iBlock-1)*length(idx);
        idx_ = idx+offset;
        try
            X_  =X(X>=Z(idx_(1)) & X<=Z(idx_(end)));
            F{iBlock}=KDens(Z(idx_),X_,mean(h));
        catch
            X_  =X(X>=Z(idx_(1)));
            F{iBlock}=KDens(Z(idx_(1):end),X_,mean(h));
        end
    end
    F = cell2mat(F);

function F=KDens(Z,X,h)
% Kernel density convolution
% Z = time bins
% X = spike times
% h = kernel S.D.
% F = convolved output

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

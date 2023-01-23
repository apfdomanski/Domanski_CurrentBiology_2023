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
% F=a*F./length(X);
F=a*F;%./length(X);
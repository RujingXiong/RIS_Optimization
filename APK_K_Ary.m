function [w,optvalue]=APK_K_Ary(h,bit)
% method from Yaowen Zhang in "Configuring intelligent reflecting surface with performance
% guarantees: Optimal beamforming"
% z is a complex column vector. The software is to solve max |W^T Z|, where W is
% on unit circle uniformly. 
%

%   Author(s): Jialong Lu, Rujing Xiong
%   Date: 02-10-2022

%% check input argument

if ~iscolumn(h)
    error('H must be a column vector.')
end
% initialization data
N=length(h);
R=h*h';
h=conj(h);
Phi=zeros(N,1);
W=zeros(N-1,3);
W_max=zeros(N-1,1);
g=zeros(1,3);
% record the phase
for i=1:N
    Phi(i,1)=angle(h(i,1));
end
%calculate w
for j=1:3
    Phi_centre=angle(exp(1i*(Phi(N,1)+(j-2)*pi/(2^bit))));
    for k=1:N-1
        for l=0:(2^bit)-1
            if Phi_centre-pi/(2^bit)<=angle(exp(1i*(Phi(k)+l*pi/(2^(bit-1)))))&&angle(exp(1i*(Phi(k)+l*pi/(2^(bit-1)))))<=Phi_centre+pi/(2^bit)
                W(k,j)=exp(1i*l*pi/(2^(bit-1)));
            end
        end
    end
end
% choose the optimal w
for j=1:3
    d=h(N,1);
    for k=1:N-1
        d=d+h(k,1)*W(k,j);
    end
    g(j)=abs(d);   
end
[~,p]=max(g);
for j=1:N-1
    W_max(j,1)=W(j,p);
end
w=[W_max;1];
optvalue=abs(w'*h);


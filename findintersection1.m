function[intersection] = findintersection1(A,B)
[~,column] = size(A);
C=(diag(B)^-1)*A;
intersection = null([real(C),imag(C)]);
intersection = 1i*intersection(1:column,:) + intersection(column+1:end,:);
clear C

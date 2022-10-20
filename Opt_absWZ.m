function [W, valueOpt] = Opt_absWZ(Z,numBits)
% Syntex: [W, valueOpt] = Opt_absWZ(Z,numBits)
% Z is a complex vector. The software is to solve max |W^H Z|, where W is
% on unit circle uniformly. 


%   Author(s): Tiebin Mi, Rujing Xiong
%   Date: 09-29-2022

%% check input argument
if ~exist('Z','var')
    error('The Z is empty.')
end

if ~exist('numBits','var')
    warning('The numBits is empty. Set numBits=4 by default.')
    numBits = 4;
end

if ~isvector(Z)
    error('Z must be a vector.')
end

%%
numElements = length(Z);

angleWeightRange = 2*pi/2^numBits;

%%
thetaZ = mod(angle(Z),2*pi);

thetaZNormalized = mod(thetaZ,angleWeightRange);

[~, indexThetaZ] = sort(thetaZNormalized,'ascend');

%%
thetaZShift = round((thetaZ-thetaZNormalized)/angleWeightRange);
thetaZShifting = repmat(thetaZShift.',2^numBits*numElements,1);

%%
angleWeightNormalized = zeros(numElements*2^numBits,numElements);

for iPart = 1:1:2^numBits
    angleWeightNormalized((iPart-1)*numElements+indexThetaZ, indexThetaZ) = (iPart-1)*ones(numElements,numElements) + triu(ones(numElements,numElements)).';
end

angleWeightNormalized = angleWeightNormalized - thetaZShifting;
angleWeightNormalized = mod(angleWeightNormalized,2^numBits);

%%
weightChecked = exp(1i*angleWeightRange.*angleWeightNormalized);

[valueOpt, indexOpt] = max(abs(weightChecked*Z));

%%
W = weightChecked(indexOpt,:)';

end
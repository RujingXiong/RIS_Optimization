function[Vector]=finds2(a,B)
Omega = 2*pi/2^B;
[n,d] = size(a);
numGroup = nchoosek(n, 2*d-1);
indexGroup = nchoosek(1:n, 2*d-1);
centerSectors = permn(0:2^(B-1)-1, 2*d-1);
boundarySectors = centerSectors+0.5;
[number,~] = size(boundarySectors);
Vector = zeros(n,1);
% parpool(4)
for i = 1:(numGroup*number)
    disp(i)
    disp(numGroup*number)
    intersection = findintersection1(a(indexGroup(floor((i-0.5)/number)+1,:),:),exp(1i*Omega*boundarySectors(mod((i-0.5), number)+0.5,:).'));
    [~, size1] = size(intersection);
    for j =1:size1
        Vector = [Vector,findvector1(intersection(:,j),a,indexGroup(floor((i-0.5)/number)+1,:),n,B)];
    end
end
% delete(gcp('nocreate'))
Vector(:,all(Vector == 0)) = [];
clear indexGroup
clear centerSectors
clear boundarySectors
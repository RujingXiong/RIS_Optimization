function[Vector]=findvector1(intersection,A,indexGroup,n,B)
[~,d]=size(A);
group = permn(2:3,2*d-1).';
vector = zeros(n,2^(2*d-1));
for i=1:2*d-1
    group1 = group(i,:);
    group1(group1==2) = exp(1i*(angle(A(indexGroup(i),:)*intersection)+pi/2^B));
    group1(group1==3) = exp(1i*(angle(A(indexGroup(i),:)*intersection)-pi/2^B));
    group(i,:) = group1;
end
for i=1:1:n
    if ismember(i,indexGroup)
        [~,l] = find(indexGroup==i);
        vector(i,:) = group(l,:);
    else
        for j=0:1:2^B-1
            if distance(angle(A(i,:)*intersection),j*pi/2^(B-1)) <= pi/2^B
                vector(i,:) = exp(1i*j*pi/2^(B-1));
            end
        end
    end
end
Vector = vector;
clear group
clear group1

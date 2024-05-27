function [permutationGroup] = permutation(orginGroup)
a = (orginGroup~=0);
if sum(a(:)) > 2
    [maxval1] = max(orginGroup);
    orginGroup(maxval1) = 0;
    permutationGroup = [orginGroup;permutation(orginGroup)];
    [maxval2] = max(orginGroup);
    orginGroup(maxval2) = 0;
    permutationGroup = [permutationGroup;orginGroup;permutation(orginGroup)];
    orginGroup(maxval1) = maxval1;
    permutationGroup = [permutationGroup;orginGroup;permutation(orginGroup)];
else    
    permutationGroup = [];
end





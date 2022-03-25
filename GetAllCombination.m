function [AllCombination] = GetAllCombination (speciesToConsider)

speciesToConsider(1) = [];
AllCombination = [];
k = 1;
AllCombination(k,1) = 1;
k = k + 1;
for i = 1:length(speciesToConsider)
   Comnination = combntns(speciesToConsider,i);
   count = size(Comnination);
   for j = 1:count(1)
       AllCombination(k,1:(count(2)+1)) = [1  Comnination(j,:)];
       k = k + 1;
   end
end
%AllCombination(AllCombination == 0) =NaN; 
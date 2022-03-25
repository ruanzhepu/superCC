function [modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversibleModel(model)

%model = models(i).COBRAmodel;
NotEXRxn = {};
AllRxn = model.rxns;
for j=1:length(AllRxn)
    temout1 = regexp(AllRxn(j),'EX_');
    temout2 = regexp(AllRxn(j),'bio1');    
    if ((isempty(temout1{1,1})) & (isempty(temout2{1,1})))
        NotEXRxn(length(NotEXRxn)+1,1) = AllRxn(j);
    end
end
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model,'sRxns',NotEXRxn,'orderReactions', true);
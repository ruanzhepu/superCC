function [EssensialReaction, NotEssensialReaction, NotEssensialReactionButFuntion] = getEssensialReaction (model1, model2, Rxnlist)
%% find the essensial reactions that in Rxnlist
%% model1: base models
%% model2: model with model.rxns that include Rxnlist, used in combineTwoModel that provide information of folulat of
%% reactions in Rxnlist
tef = 1e-6; %biomass cutoff
NotEssensialReaction = {};
NotEssensialReactionButFuntion = {};
EssensialReaction = {};
[Addedmodel] = combineTwoModel(model1,model2,Rxnlist);
Solution = optimizeCbModel(Addedmodel);
if (Solution.f > tef)
    for i = 1:length(Rxnlist)
        tempmodel = removeRxns(Addedmodel,Rxnlist(i));
        tempS = optimizeCbModel(tempmodel);
        if (tempS.f < tef)
            EssensialReaction(length(EssensialReaction) + 1) = Rxnlist(i);             
        elseif (abs(Solution.f - tempS.f) < tef)
            NotEssensialReaction(length(NotEssensialReaction) +1 ) = Rxnlist(i);
        else
            NotEssensialReactionButFuntion(length(NotEssensialReactionButFuntion)+1) = Rxnlist(i);
        end
    end
end
end

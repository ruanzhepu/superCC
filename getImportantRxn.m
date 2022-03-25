initCobraToolbox(false);
changeCobraSolver ('gurobi', 'all', 1);

models(1).COBRAmodel = readCbModel('D2Addedmodel20210516.mat');
models(2).COBRAmodel = readCbModel('D2model20210519.mat');
models(3).COBRAmodel = readCbModel('X1model20210619.mat');
tempmodel = GetMultiSpeciesModel([1:2], {'D2model20210519.mat','X1model20210619.mat'}, 'D:\models\supersoldier\');
models(4).COBRAmodel = tempmodel(2).COBRAmodel; 
%load('X1D1Joint.mat');
%clear X1D1Joint
mediumtxt = 'ComMedium_XXXBJ.txt';
medium = readtable(mediumtxt); %Read the culture conditions.  It is recommended to store the exchange response conditions for all models in one file and use this file for the rest.
exnames = strcat('EX_',medium.compounds,'[u]'); %for joint models
exnamesB = strcat('EX_',medium.compounds,'_e0'); %for the first model, which is basic model and control model

difRxn = setdiff(models(3).COBRAmodel.rxns, models(1).COBRAmodel.rxns);
%difRxn(strmatch('EX_',difRxn))=[]; %remove the EX reactions

for i = 1:length(models)
    tempEX = strmatch('EX_',models(i).COBRAmodel.rxns);
    temout = regexp(models(i).COBRAmodel.rxns(tempEX(1)),'[u]');
    if (isempty(temout{1,1}))
        models(i).COBRAmodel = changeRxnBounds(models(i).COBRAmodel, exnamesB, medium.minFlux, 'l');
        models(i).COBRAmodel = changeRxnBounds(models(i).COBRAmodel, exnamesB, medium.maxFlux, 'u');
    else
        models(i).COBRAmodel = changeRxnBounds(models(i).COBRAmodel, exnames, medium.minFlux, 'l');
        models(i).COBRAmodel = changeRxnBounds(models(i).COBRAmodel, exnames, medium.maxFlux, 'u');
    end
    models(i).Solution = optimizeCbModel(models(i).COBRAmodel);
end

%%Identify essentail reactions: perform a gene knocked-out analysis%%
models(i).COBRAmodel.genes = models(i).COBRAmodel.rxns;
[grRatio] = singleGeneDeletion(models(i).COBRAmodel, 'FBA');
EssentialRxns = models(i).COBRAmodel.rxns(find(grRatio<1));
%%no EssentialRxns was found. use fva of .
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversibleModel(models(i).COBRAmodel);
[MinF,MaxF] = fluxVariability(modelIrrev,100);

MinFwithFlux = intersect(find(MinF>0),strmatch('X1',modelIrrev.rxns));
MinFwithFluxRxn = modelIrrev.rxns(MinFwithFlux);

for i = 1:length(MinFwithFluxRxn)
   MinFwithFluxRxn(i)  = regexprep(MinFwithFluxRxn(i), 'X1model20210617I', '');
   MinFwithFluxRxn(i)  = regexprep(MinFwithFluxRxn(i), 'X1model20210617', '');
   MinFwithFluxRxn(i)  = regexprep(MinFwithFluxRxn(i), '_f', '');
   MinFwithFluxRxn(i)  = regexprep(MinFwithFluxRxn(i), '_b', '');
   MinFwithFluxRxn(i)  = regexprep(MinFwithFluxRxn(i), '\[u\]tr', '_e0');
end
difRxnwithMinF = intersect(MinFwithFluxRxn,difRxn);
printRxnFormula(models(3).COBRAmodel,difRxnwithMinF);

%%The above reactions found were added to D2added to test whether biomass could be increased.
[Addedmodel] = combineTwoModel(models(2).COBRAmodel,models(3).COBRAmodel,difRxnwithMinF);
Solution = optimizeCbModel(Addedmodel);

%%%We did it the other way around. We added the reactions to X-1.

MinFwithFlux2 = intersect(find(MinF>0),strmatch('D2',modelIrrev.rxns));
MinFwithFluxRxn2 = modelIrrev.rxns(MinFwithFlux2);

for i = 1:length(MinFwithFluxRxn2)
   MinFwithFluxRxn2(i)  = regexprep(MinFwithFluxRxn2(i), 'D2model20210519I', '');
   MinFwithFluxRxn2(i)  = regexprep(MinFwithFluxRxn2(i), 'D2model20210519', '');
   MinFwithFluxRxn2(i)  = regexprep(MinFwithFluxRxn2(i), '_f', '');
   MinFwithFluxRxn2(i)  = regexprep(MinFwithFluxRxn2(i), '_b', '');
   MinFwithFluxRxn2(i)  = regexprep(MinFwithFluxRxn2(i), '\[u\]tr', '_e0');
end
difRxn2 = setdiff(models(2).COBRAmodel.rxns, models(3).COBRAmodel.rxns);
difRxnwithMinF2 = intersect(MinFwithFluxRxn2,difRxn2);
printRxnFormula(models(2).COBRAmodel,difRxnwithMinF2);

%%The above reactions found were added to D2added to test whether biomass could be increased.
[Addedmodel2] = combineTwoModel(models(3).COBRAmodel,models(2).COBRAmodel,difRxnwithMinF2);
Solution2 = optimizeCbModel(Addedmodel);











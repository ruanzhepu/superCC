function [EXwithC] = validateModel (model, modelName)

tempmodel = model;
modelname = modelName;
EX = tempmodel.rxns(strmatch('EX_',tempmodel.rxns));
EXwithC = {};
%i = 43;
for i = 1:length(EX); %Find all the exchange reactions starting with [c]
	temmet = tempmodel.mets(find(tempmodel.S(:,findRxnIDs(tempmodel,EX(i)))));
	tempvalue = strfind(temmet,'[c]');
	if (~(isempty(tempvalue{1})))
		EXwithC(length(EXwithC)+1) = temmet;
	end
end
EXwithC = strcat(modelname,EXwithC);


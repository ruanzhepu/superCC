function [models] = GetMultiSpeciesModel(scenarioID, speciesToConsider, filenames, directory,PercentOfSpeciesBio)
%speciesToConsider = [1:7];
%directory = 'D:\models\supersoldier\';
%PercentOfSpeciesBio is the ratio of artificial biomass of each strain, and the number of its numerical value must be the number of models, and must be greater than 0.
[AllCombination] = GetAllCombination (speciesToConsider);
SizeAllCombination = size(AllCombination);
for i = 1:SizeAllCombination(1)
	templine = AllCombination(i,:);
	templine(templine == 0) = [];
	if (length(templine) == 1)
		tempfile = strcat(directory, char(filenames(templine(i))));
		[filepath,name,ext] = fileparts(tempfile);
		models(i).COBRAmodel = readCbModel(tempfile);
		models(i).Name = name;
		rxnBiomassId = findRxnIDs(models(i).COBRAmodel, 'bio1'); % ids
		models(i).COBRAmodel.c(:) = 0;
		models(i).COBRAmodel.c(rxnBiomassId) = 1;%Set objective function
	else
		tempNameList = {};
		EXrxnAdd = {};
		for j = 1:length(templine)
			tempfile = strcat(directory, char(filenames(templine(j))));
			[filepath,name,ext] = fileparts(tempfile);
			tempNameList{length(tempNameList)+1,1} = name;
			tempmodel = readCbModel(tempfile);
			tempmodel = changeObjective(tempmodel,'bio1');
			tempmodel.lb(findRxnIDs(tempmodel,'bio1')) = 0;
			tempmodels{j,1}=tempmodel;
			EXwithC = validateModel(tempmodel, name);
			EXrxnAdd = [EXrxnAdd, EXwithC];
		end		
		models(i).COBRAmodel = createMultipleSpeciesModel(tempmodels, 'nameTagsModels',tempNameList);
		models(i).Name = tempNameList;
		if (~isempty(EXrxnAdd))
			templb = 0 * ones(length(EXrxnAdd), 1);
			tempub = 1000 * ones(length(EXrxnAdd), 1);
			models(i).COBRAmodel = addExchangeRxn(models(i).COBRAmodel, EXrxnAdd,templb,tempub);
		end
		rxnBiomass = strcat(tempNameList, 'bio1'); % biomass reaction names
		rxnBiomassId = findRxnIDs(models(i).COBRAmodel, rxnBiomass); % ids
		models(i).COBRAmodel.c(:) = 0;
		if (scenarioID == 1)
			models(i).COBRAmodel.c(rxnBiomassId) = 1;%All strains were set as objective functions
		elseif (scenarioID == 2)
			models(i).COBRAmodel.c(rxnBiomassId(1)) = 1; %Set the first strain as the objective function
		elseif (scenarioID == 3)                        %Add a biomass reaction, the total biomass is the sum of the biomass of all strains, and the biomass ratio of all strains is 1.
			metlist = strcat(models(i).Name,'cpd11416[c]')';
			metlist(end +1) = {'cpd11416Joint[c]'};
			stoiList = [-1*ones(length(models(i).Name), 1)' length(models(i).Name)];
			stoiList = stoiList/length(models(i).Name);
			models(i).COBRAmodel = addReaction(models(i).COBRAmodel,'bio1Joint','metaboliteList', metlist,'stoichCoeffList', stoiList);
			models(i).COBRAmodel = addExchangeRxn(models(i).COBRAmodel,'cpd11416Joint[c]');
			models(i).COBRAmodel.c(:) = 0;
			models(i).COBRAmodel.c(findRxnIDs(models(i).COBRAmodel,'bio1Joint')) = 1;
		elseif (scenarioID == 4) %Add a biomass reaction, and the total biomass is the sum of the biomass of all strains, which is different from 3, where the proportion of the biomass of all strains in the total biomass is artificially specified.
%                                 The proportion value will be extracted from PercentOfSpeciesBio according to the composition of each flora model, after which it will be summed and normalized to 1, i.e. the coefficient on the right side of the biomass reaction is 1, and the coefficient on the left is "-1* coefficients /sum(coefficients)".
			metlist = strcat(models(i).Name,'cpd11416[c]')';			
			tempSpeciesList = AllCombination(i,:);
			tempPercentOfSpeciesBio = PercentOfSpeciesBio(tempSpeciesList(find(tempSpeciesList)));
			if (length(metlist) == length(tempPercentOfSpeciesBio))				
				metlist(end +1) = {'cpd11416Joint[c]'};
				stoiList = -1*tempPercentOfSpeciesBio/sum(tempPercentOfSpeciesBio);
				stoiList(end + 1) = 1;
				models(i).COBRAmodel = addReaction(models(i).COBRAmodel,'bio1Joint','metaboliteList', metlist,'stoichCoeffList', stoiList);
				models(i).COBRAmodel = addExchangeRxn(models(i).COBRAmodel,'cpd11416Joint[c]');
				models(i).COBRAmodel.c(:) = 0;
				models(i).COBRAmodel.c(findRxnIDs(models(i).COBRAmodel,'bio1Joint')) = 1;
			else 
				error ('Numbers in PercentOfSpeciesBio are not equal to those in speciesToConsider. Aborting')
			end
		else
			error('Unknown model scenario. Aborting')
		end			
	end
end




function [BioAndDegradation, models, medium] = GetEnhancedCombination (scenarioID, directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compound, PercentOfSpeciesBio)
%scenarioID = 1; The objective function of the bacterial community model is the sum of the biomass of all strains, and the biomass proportion of all strains is unlimited (the biomass of all strains is allowed to be 0). 
%scenarioID = 2; The objective function of the bacterial community model was the biomass of the first strain, which was used to find other strains that could promote the biomass of the first strain.
%scenarioID = 3; The objective function of the bacterial community model was the sum of the biomass of all strains, and the total biomass was the sum of the biomass of all strains, and the biomass ratio of all strains was 1.  
%scenarioID = 4; The objective function of the bacterial community model is the sum of the biomass of all strains, and the total biomass is the sum of the biomass of all strains, and the biomass ratio of all strains is artificially specified as PercentOfSpeciesBio.  
%PercentOfSpeciesBio is the ratio of artificial biomass of each strain, and the number of its numerical value must be the number of models, and must be greater than 0.
%PercentOfSpeciesBio must be specified when scenarioID is 4 and not if scenarioID is 1-3.  
%AllCombinationModel: the multispeciesModel cover all the combination of all species
%medium: Medium contains all the EX_rxn[u] of all model, which is suitable for all multispeciesModel.
%directory = 'D:\models\supersoldier\';
%CNsourcesName = {'cpd00013','cpd00027'};
%CNsourcesName = {'cpd50000'};
%CNsourcesNumber = [-100];
%CNsourcesNumber = [-50,-50];
%compound = {'cpd02568','cpd00013','cpd50000'};
%compound = {'cpd02568','cpd50000'};
%mmMedium = 'mmMedium.txt'
%foldchange = 1.2

%initCobraToolbox(false);
%changeCobraSolver ('gurobi', 'all', 1);

%get all model file
warning off

if exist('filenames') == 0
	filenames = { 
	'D2model20210519.mat' ...
	'X_1model20210602.mat' ...
	'E44model20210430.mat' ...
	'Lm5model20210602.mat' ...
	'NRRLmodel20210526.mat' ...
	'Ag3model20210617.mat' ...
	'Ba1model20210428.mat' ...
    };
	speciesToConsider = [1:7];
end	


switch scenarioID
	case 1
		[models] = GetMultiSpeciesModel(scenarioID, speciesToConsider, filenames, directory, PercentOfSpeciesBio); %The objective function of the bacterial community model is the sum of the biomass of all strains, and the biomass proportion of all strains is unlimited (the biomass of all strains is allowed to be 0). 
	case 2
		[models] = GetMultiSpeciesModel(scenarioID, speciesToConsider, filenames, directory, PercentOfSpeciesBio); %The objective function of the bacterial community model was the biomass of the first strain, which was used to find other strains that could promote the biomass of the first strain.
	case 3
		[models] = GetMultiSpeciesModel(scenarioID, speciesToConsider, filenames, directory, PercentOfSpeciesBio); %The objective function of the bacterial community model was the sum of the biomass of all strains, and the total biomass was the sum of the biomass of all strains, and the biomass ratio of all strains was 1. 
	case 4
		[models] = GetMultiSpeciesModel(scenarioID, speciesToConsider, filenames, directory, PercentOfSpeciesBio); %The objective function of the bacterial community model is the sum of the biomass of all strains, and the total biomass is the sum of the biomass of all strains, and the biomass ratio of all strains is artificially specified as PercentOfSpeciesBio.  
	otherwise
        error('Unknown model scenario. Aborting')
end
	
[medium] = GetMedium(directory,filenames, mmMedium, CNsourcesName, CNsourcesNumber);%Based on EX reactions of all models, the inorganic elements provided by MM medium, carbon and nitrogen sources provided, etc., are required to obtain culture conditions suitable for all bacterial combinations.

exnames = strcat('EX_',medium.compounds,'[u]'); %for joint models
exnamesB = strcat('EX_',medium.compounds,'_e0'); %for the first model, which is basic model and control model
compoundJoint = strcat('EX_',compound,'[u]'); 
compoundB = strcat('EX_',compound,'_e0'); 

k = 1; %Control the number of lines of output
BioAndDegradation= {}; %Save the result as the final output. The second column is the total biomass, followed by the biomass of all strains.
colnames = [{'Combination','TotalBiomass'},filenames,compound,'biofoldchange'];
countfilenames = length(filenames);
countcompound = length(compound);
for (i = 1:length(colnames))
    BioAndDegradation(k,i)= colnames(i); %The first line writes the column name.
end
clear i
k = k +1; %The second line begins writing the result.
for i = 1:length(models)
	if (i == 1)
		BioAndDegradation{k,1} =  models(i).Name;
		models(i).COBRAmodel = changeRxnBounds(models(i).COBRAmodel, exnamesB, medium.minFlux, 'l');
		models(i).COBRAmodel = changeRxnBounds(models(i).COBRAmodel, exnamesB, medium.maxFlux, 'u');
		solution = optimizeCbModel(models(i).COBRAmodel);
		%[minF, maxF] = fluxVariability(models(i).COBRAmodel, 100);%too slow, use optimizeCbModel.
		if(solution.stat == 1)
			%%%print resluts%%%			
			BioAndDegradation{k,2} = solution.f;
			%%
			models(i).COBRAmodel = changeRxnBounds(models(i).COBRAmodel, 'bio1', solution.f, 'l');%Set biomass to the maximum
			models(i).COBRAmodel.c(:) = 0;
			temIDs = findRxnIDs(models(i).COBRAmodel,compoundB);
			temIDsList = temIDs;
			temIDs(temIDs == 0) = [];
			temIDCount = find(temIDsList);
			models(i).COBRAmodel.c(temIDs) = 1;%Set the objective function, take the degradation amount of the target pollutant as the objective function, and find the minimum value.
			solutionMin = optimizeCbModel(models(i).COBRAmodel,'min');
			models(i).MinCompoundSolution = solutionMin;
			if (solutionMin.stat==1)
				for h = 1:length(temIDs)
					tempn=temIDCount(h);
					BioAndDegradation{k,2+countfilenames+tempn} = solutionMin.x(temIDs(h));
				end
			else
				BioAndDegradation{k,2+countfilenames} = 'NoFBAsolutionMin';
			end			
			BioAndDegradation{k,end} = 1;
		else
			BioAndDegradation{k,2} = 'NoFBAsolution';
		end
		models(i).COBRAmodel.c(:) = 0;
		models(i).COBRAmodel = changeObjective(models(i).COBRAmodel,'bio1'); %set objet to bio1
		k = k + 1;
	else
		ComName = char(models(i).Name(1)); 
		for j = 2:length(models(i).Name)
			ComName = strcat(ComName, char(models(i).Name(j))); %names of combined species
		end
		BioAndDegradation{k,1} = ComName;
		models(i).COBRAmodel = changeRxnBounds(models(i).COBRAmodel, exnames, medium.minFlux, 'l');
		models(i).COBRAmodel = changeRxnBounds(models(i).COBRAmodel, exnames, medium.maxFlux, 'u');
		solution = optimizeCbModel(models(i).COBRAmodel);
		if(solution.stat == 1)
			%%%print resluts%%%	
			BioAndDegradation{k,2} = solution.f;	
			bioIDs = find(models(i).COBRAmodel.c);
			for m = 1:length(bioIDs)
				tempnum=strmatch(strcat(models(i).Name(m),'.mat'),filenames);
				BioAndDegradation{k,2+tempnum} = solution.x(bioIDs(m));
			end
			%%
			models(i).COBRAmodel.lb(bioIDs) = solution.x(bioIDs);%Set biomass to the maximum
			models(i).COBRAmodel.c(:) = 0;
			temIDs = findRxnIDs(models(i).COBRAmodel,compoundJoint);
			temIDsList = temIDs;
			temIDs(temIDs == 0) = [];
			temIDCount = find(temIDsList);
			models(i).COBRAmodel.c(temIDs) = 1;%Set the objective function, take the degradation amount of the target pollutant as the objective function, and find the minimum value.  
			solutionMin = optimizeCbModel(models(i).COBRAmodel,'min');
			models(i).MinCompoundSolution = solutionMin;
			if (solutionMin.stat==1)
				for n = 1:length(temIDs)
					tempnum2=temIDCount(n);
					BioAndDegradation{k,2+countfilenames+tempnum2} = solutionMin.x(temIDs(n));
				end
			else
				BioAndDegradation{k,2+countfilenames} = 'NoFBAsolutionMin';
			end	
			BioAndDegradation{k,end} = BioAndDegradation{k,2}/BioAndDegradation{2,2};
		else
			BioAndDegradation{k,2} = 'NoFBAsolution';
		end
		models(i).COBRAmodel.c(:) = 0;
		models(i).COBRAmodel.c(bioIDs) = 1;
		k = k+1;
	end
end













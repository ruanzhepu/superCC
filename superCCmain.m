function [BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell)
%find best combination and reations 
%PercentOfSpeciesBio is the ratio of artificial biomass of each strain, and the number of its numerical value must be the number of models, and must be greater than 0.
%scenarioID = 1; The objective function of the bacterial community model is the sum of the biomass of all strains, and the biomass proportion of all strains is unlimited (the biomass of all strains is allowed to be 0). 
%scenarioID = 2; The objective function of the bacterial community model was the biomass of the first strain, which was used to find other strains that could promote the biomass of the first strain.
%scenarioID = 3; The objective function of the bacterial community model was the sum of the biomass of all strains, and the total biomass was the sum of the biomass of all strains, and the biomass ratio of all strains was 1.  
%scenarioID = 4; The objective function of the bacterial community model is the sum of the biomass of all strains, and the total biomass is the sum of the biomass of all strains, and the biomass ratio of all strains is artificially specified as PercentOfSpeciesBio.  
%PercentOfSpeciesBio must be specified when scenarioID is 4 and not if scenarioID is 1-3.  
%for syntheticCell = 1,  design the singel synthetic Cell based on the best comnination.
%for syntheticCell = 0, NOT desigen the singel synthetic Cell based on the best comnination. by default,syntheticCell = 0.  

%filenames = {'D2model20210519.mat','X1model20210619.mat','Ba1modelAdded20210630.mat','BRAmodel20210602.mat' ...
%    'Lm5model20210602.mat', 'E44model20210430.mat', 'WWmodel20210430.mat'};
%speciesToConsider = [1:7];
%CNsourcesName = {'cpd50000','cpd00027','cpd00013','cpd00209'};
%CNsourcesNumber = [-50,-50,-50,-50];
%mmMedium = 'mmMedium.txt';
%directory = 'D:\models\supersoldier\XXXBJ_models\R\';
%compoundtest = {'cpd50000','cpd00027','cpd00013','cpd00209'};
%scenarioID = 2;
%syntheticCell = 0;
%PercentOfSpeciesBio=[1 1];
%PercentOfSpeciesBio=[1 1 1 1 1 1 1];

if ((exist('PercentOfSpeciesBio') == 0) & (scenarioID == 4))
	error ('PercentOfSpeciesBio must be defined when scenarioID is 4. Aborting')
end

if exist('PercentOfSpeciesBio') == 0
	PercentOfSpeciesBio = [];
end

if exist('syntheticCell') == 0
	PercentOfSpeciesBio = 0;
end

[BioAndDegradationAll, modelsAll, mediumAll] = GetEnhancedCombination (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio);


%find the best combination with least number of models
%for syntheticCell = 1, do as the follows to desigen the singel synthetic cell based on the best comnination.
%for syntheticCell = 0, NOT do as the follows to desigen the singel synthetic cell based on the best comnination.

if (syntheticCell == 1)
	biomass=cell2mat(BioAndDegradationAll(2:end,2));
	biomass = roundn(biomass,-4);
	maxdex=min(find(biomass == max(biomass)));

	maxmodel.model = modelsAll(maxdex).COBRAmodel;
	maxmodel.Name = modelsAll(maxdex).Name;
	exnames = strcat('EX_',mediumAll.compounds,'[u]'); %for joint models
	exnamesB = strcat('EX_',mediumAll.compounds,'_e0'); %for single model

	maxmodel.modelIrrev = convertToIrreversibleModel(maxmodel.model);
	[maxmodel.MinF,maxmodel.MaxF] = fluxVariability(maxmodel.modelIrrev,100);
	maxmodel.EssentialRxn = maxmodel.modelIrrev.rxns(find(maxmodel.MinF>0 | maxmodel.MaxF<0));
	for j = 1:length(maxmodel.EssentialRxn)
		maxmodel.EssentialRxn(j)  = regexprep(maxmodel.EssentialRxn(j), '_f', '');
		maxmodel.EssentialRxn(j)  = regexprep(maxmodel.EssentialRxn(j), '_b', '');
		maxmodel.EssentialRxn(j)  = regexprep(maxmodel.EssentialRxn(j), '\[u\]tr', '_e0');
		for n = 1:length(maxmodel.Name)
           maxmodel.EssentialRxn(j)  = regexprep(maxmodel.EssentialRxn(j), char(maxmodel.Name(n)), '');
		end
		maxmodel.EssentialRxn(j)  = regexprep(maxmodel.EssentialRxn(j), 'IEX', 'EX');
		maxmodel.EssentialRxn(j)  = regexprep(maxmodel.EssentialRxn(j), '\[u\]', '_e0');
	end
	maxmodel.uniqueEssentialRxn = unique(maxmodel.EssentialRxn);

	%read sigle models
	for i = 1:length(maxmodel.Name)
		tempfile = strcat(directory, char(maxmodel.Name(i)),'.mat');
		[filepath,name,ext] = fileparts(tempfile);
		Singlemodels(i).COBRAmodel = readCbModel(tempfile);
		Singlemodels(i).Name = name;
		rxnBiomassId = findRxnIDs(Singlemodels(i).COBRAmodel, 'bio1'); % ids
		Singlemodels(i).COBRAmodel.c(:) = 0;
		Singlemodels(i).COBRAmodel.c(rxnBiomassId) = 1;%Set objective function
		Singlemodels(i).COBRAmodel = changeRxnBounds(Singlemodels(i).COBRAmodel, exnamesB, mediumAll.minFlux, 'l');
		Singlemodels(i).COBRAmodel = changeRxnBounds(Singlemodels(i).COBRAmodel, exnamesB, mediumAll.maxFlux, 'u');    
	end

	%find different rxns
	EssentialRxn = maxmodel.uniqueEssentialRxn;
	for i = 2:length(Singlemodels)
        Singlemodels(i).difRxn = setdiff(Singlemodels(i).COBRAmodel.rxns, Singlemodels(1).COBRAmodel.rxns); %Differential reactions, 7D-2 in the first place, X-1 in the second.
        Singlemodels(i).difRxnEssensial = intersect(Singlemodels(i).difRxn, EssentialRxn); %The difference in X-1 has to be reflected, and in general X-1 goes to the second place.
        EssentialRxn = setdiff(EssentialRxn, Singlemodels(i).difRxnEssensial);
	end

	%add EssentialRxn test if biomass produced ,then test which rxn is essential by deletion one by one
	MinMixmodel = Singlemodels(1).COBRAmodel; %first model as the initial MinMixmodel
	Solution = optimizeCbModel(MinMixmodel);
	Singlemodels(1).Bio = Solution.f;
	for i = 2:length(Singlemodels)
		[Singlemodels(i).EssensialReaction, Singlemodels(i).NotEssensialReaction, Singlemodels(i).NotEssensialReactionButFuntion] = getEssensialReaction (MinMixmodel, Singlemodels(i).COBRAmodel, Singlemodels(i).difRxnEssensial);
		if (~isempty(Singlemodels(i).EssensialReaction))
			MinMixmodel = combineTwoModel(MinMixmodel, Singlemodels(i).COBRAmodel,Singlemodels(i).EssensialReaction);
			Solution = optimizeCbModel(MinMixmodel);
			Singlemodels(i).EssensialReactionBio = Solution.f;
		end
		if (~isempty(Singlemodels(i).NotEssensialReactionButFuntion))
			MinMixmodel = combineTwoModel(MinMixmodel, Singlemodels(i).COBRAmodel,Singlemodels(i).NotEssensialReactionButFuntion);
			Solution = optimizeCbModel(MinMixmodel);
			Singlemodels(i).FuntionBio = Solution.f;
		end
	end
	superModel = MinMixmodel;
else
	Singlemodels = {};
	superModel = {};
	disp("Analyze best combination only without single cell design!");	
end
end


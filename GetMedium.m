function [medium] = GetMedium(directory,filenames, mmMedium, CNsourcesName, CNsourcesNumber)
%Based on EX reactions of all models, the inorganic elements provided by MM medium, carbon and nitrogen sources provided, etc., are required to obtain culture conditions suitable for all bacterial combinations.  
%CNsourcesName = {'cpd00013','cpd00027'};
%CNsourcesNumber = [-50,-50];
%mmMedium = 'mmMedium.txt'
EXrxn = {};
for i = 1:length(filenames)
	tempfile = strcat(directory, char(filenames(i)));
	tempmodel = readCbModel(tempfile);
	IDs = find(strncmp(tempmodel.rxns, 'EX_', 3));
	temEXrxn = tempmodel.rxns(IDs);
	EXrxn = [EXrxn, temEXrxn'];
end
EXrxn = unique(EXrxn);

EXcpd = {};
for j = 1:length(EXrxn)
	tempcpd = strsplit(char(EXrxn(j)),'_');
	EXcpd(length(EXcpd)+1) = tempcpd(2);
end
EXcpd = unique(EXcpd)';

mmmedium = readtable(strcat(directory, char(mmMedium))); 
cpddiff = setdiff(EXcpd,mmmedium.compounds);
cpddiff = setdiff(cpddiff,CNsourcesName);
mmcount = size(mmmedium);
mmcount = mmcount(1);

for i = 1:length(cpddiff)
	mmmedium.compounds(mmcount+i) = cpddiff(i);
	mmmedium.concentration(mmcount+i) = 0.01;
	mmmedium.minFlux(mmcount+i) = 0;
	mmmedium.maxFlux(mmcount+i) = 1000;
end
mmcount = size(mmmedium);
mmcount = mmcount(1);
for i = 1:length(CNsourcesName)
	mmmedium.compounds(mmcount+i) = CNsourcesName(i);
	mmmedium.concentration(mmcount+i) = 0.01;
	mmmedium.minFlux(mmcount+i) = CNsourcesNumber(i);
	mmmedium.maxFlux(mmcount+i) = 1000;
end
medium =  mmmedium;
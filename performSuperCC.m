%new analysis
%cpd00013 (NH4) and cpd00209(NO3) should be defined, if not in the medium, then use 0
initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');%The 'glpk' model is prone to error.
filenames = {'D2model20211129.mat','X1model20211129.mat','Ba1modelAdded20211129.mat','BRAmodel20211129.mat' ...
    'Lm5model20211129.mat', 'E44model20211129.mat', 'WWmodel20211129.mat'};
speciesToConsider = [1:7];
PercentOfSpeciesBio = [];
%PercentOfSpeciesBio = [0.1 0.1 0.1 0.1 0.1 0.1 0.1];
CNsourcesName = {'cpd50000','cpd00027','cpd00013','cpd00209'};
%CNsourcesNumber = [-50,-50,-50,-50];
CNsourcesNumber = [-100,0,0,0];
mmMedium = 'mmMedium.txt';
directory = 'D:\models\supersoldier\supersoldier\CO2revisedModel\';
compoundtest = {'cpd50000','cpd00027','cpd00013','cpd00209'};
scenarioID = 1;
%scenarioID = 3;
syntheticCell = 0;
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);
save('R_BO_100_S1.mat')
clear BioAndDegradationAll modelsAll mediumAll Singlemodels superModel

CNsourcesNumber = [-100,-100,0,0];
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);
save('R_BO_G_100_100_S1.mat')
clear BioAndDegradationAll modelsAll mediumAll Singlemodels superModel

CNsourcesNumber = [-100,0,-100,0];
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);
save('R_BO_NH4_100_100_S1.mat')
clear BioAndDegradationAll modelsAll mediumAll Singlemodels superModel

CNsourcesNumber = [-100,0,0,-100];
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);
save('R_BO_NO3_100_100_S1.mat')
clear BioAndDegradationAll modelsAll mediumAll Singlemodels superModel

CNsourcesNumber = [-100,-100,-100,0];
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);
save('R_BO_G_NH4_100_100_100_S1.mat')
clear BioAndDegradationAll modelsAll mediumAll Singlemodels superModel

CNsourcesNumber = [-100,-100,0,-100];
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);
save('R_BO_G_NO3_100_100_100_S1.mat')
clear BioAndDegradationAll modelsAll mediumAll Singlemodels superModel

CNsourcesNumber = [-100,-100,-100,-100];
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);
save('R_BO_G_NH4_NO3_100_100_100_100_S1.mat')
clear BioAndDegradationAll modelsAll mediumAll Singlemodels superModel


%more complicated medium
CNsourcesName = {'cpd50000','cpd00027','cpd00013','cpd00209','cpd00033', 'cpd00035', 'cpd00039', 'cpd00041', 'cpd00051', 'cpd00054', 'cpd00060', 'cpd00065', 'cpd00066', 'cpd00069', 'cpd00084', 'cpd00107', 'cpd00119', 'cpd00129', 'cpd00156', 'cpd00161', 'cpd00322', 'cpd00381'};
compoundtest = CNsourcesName;
CNsourcesNumber = [-100,0,0,0, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10];
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);


scenarioID = 2;
syntheticCell = 1;
CNsourcesNumber = [-100,0,0,0, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10];
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);

CNsourcesNumber = [-100, 0, 0, 0, 0, -10, 0, 0, -10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];




CNsourcesName = {'cpd50000','cpd00027','cpd00013','cpd00209', 'cpd00092'};
scenarioID = 1;
syntheticCell = 0;
CNsourcesNumber = [-100,0,0,0, -10];
compoundtest = CNsourcesName;
[BioAndDegradationAll, modelsAll, mediumAll, Singlemodels, superModel] = superCCmain (scenarioID,directory,filenames, speciesToConsider,mmMedium, CNsourcesName, CNsourcesNumber,compoundtest,PercentOfSpeciesBio,syntheticCell);





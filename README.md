# Super Community Combinations (SuperCC)

**Overview:**

A newly developed metabolic modeling pipeline, Super Community Combinations (SuperCC), which could simulate the performances of different microbiomes.



**Main function:**

superCCmain: find the best combination with least number of models, and desigen the singel synthetic cell based on the best combination.

Input parameters: 
scenarioID, PercentOfSpeciesBio, syntheticCell, filenames, speciesToConsider, CNsourcesName, CNsourcesNumber, mmMedium, directory, compoundtest.



**Other functions：**

_[For MultiSpecies]_

validateModel: Validate models to standard format.

convertToIrreversibleModel: Convert models to irreversible models.

GetAllCombination: Based on models provided for each species, the multispeciesModel cover all the combination of all species.

GetMedium: Based on EX reactions of all models, the inorganic elements provided by MM medium, carbon and nitrogen sources provided, etc., are required to obtain culture conditions suitable for all bacterial combinations. 

GetEnhancedCombination: Calculate the biomass of each combination (each multispeciesModel).

GetMultiSpeciesModel: Based on GetAllCombination results, merge species into models and all reaction IDs from different databases were converted into for consistency.

performSuperCC: Perform SuperCC under different simple/complicated media.



_[For syntheticCell]_

combineTwoModel: Combination of base models and another selected model. 

getEssensialReaction: Find the essensial reactions that in Rxnlist based on combineTwoModel results.

getImportantRxn: Find the important reactions that could increase the biomass of the species.

See superCC/doc for example documentation.



**System Requirements:**

The SuperCC framwork is tested on Windows operating systems. 

The matlab code should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the SuperCC framwork code, users should have MATLAB version 2014b or higher, and COBRA Toolbox should be installed according to the instructions from "https://opencobra.github.io/cobratoolbox/stable/installation.html".



**Package and solver dependencies:**

Packages: git, curl, COBRA Toolbox 

Solver: IBM CPLEX 12.10, GUROBI 9.1.1, MOSEK 9.2



**Models and media dependencies:**

Microbial models and the growth media should be prepared before using this code. The detailed instructions on how to get microbial models and modified them can be found in the reference: "Thiele, I., Palsson, B. A protocol for generating a high-quality genome-scale metabolic reconstruction. Nat Protoc 5, 93–121 (2010). https://doi.org/10.1038/nprot.2009.203". 



**DHP-Com:**

A novel and rational paradigm for consortium design is proposed using SuperCC, designated as DHP-Com (Degrader-Helper-Potentiator). Among the components, Potentiators refer to core species that maintain stable habitat abundance. Notably, synthetic microbiomes constructed under this framework not only exhibit enhanced ecological fitness (encompassing survival and growth) but also—most importantly—deliver substantially improved remediation performance across diverse co-contamination scenarios.

The modified single-strain metabolic models of the six representative strains involved in DHP-Com, along with their corresponding culture media, are detailed in "DHP-Com (six representative species)".



**Reference:**

Ruan, Z., Chen, K., Cao, W. et al. Engineering natural microbiomes toward enhanced bioremediation by microbiome modeling. Nat Commun 15, 4694 (2024). https://doi.org/10.1038/s41467-024-49098-z


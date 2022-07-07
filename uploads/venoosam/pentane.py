# Data sources
database(
    thermoLibraries = ['BurkeH2O2','Klippenstein_Glarborg2016','primaryThermoLibrary','thermo_DFT_CCSDTF12_BAC','DFT_QCI_thermo','CBS_QB3_1dHR','CurranPentane','JetSurF2.0'],
    reactionLibraries = ['BurkeH2O2inN2','Klippenstein_Glarborg2016','CurranPentane'],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)


# Constraints on generated species
generatedSpeciesConstraints(
    allowed=['input species','seed mechanisms','reaction libraries','thermoLibraries'],
    maximumCarbonAtoms=5,
    maximumRadicalElectrons=2,
    maximumSingletCarbenes=1,
    maximumCarbeneRadicals=0,
    allowSingletO2 = True,
)


# List of species
species(
    label='pentane',
    reactive=True,
    structure=SMILES("CCCCC"),
)


species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]"),
)

species(
    label='N2',
    reactive=False,
    structure=InChI("InChI=1/N2/c1-2"),
)

# Reaction

simpleReactor(
	temperature=[(750,'K'),(1350,'K')],
	nSims=12,
	pressure=[(1,'atm'),(10,'atm')],
	initialMoleFractions={
	"pentane":[1.2,4.99],
	"O2":19.97,
	"N2":75.05
},
	terminationConversion={
		'pentane':0.9,
},
	terminationTime=(1e6,'s'),
    sensitivity=['pentane',"O2"],
    sensitivityTemperature = (1150,'K'),
    sensitivityPressure = (10.0,'atm'),
    sensitivityMoleFractions = {"pentane":4.99,"O2":19.97},
    sensitivityThreshold=0.01,
    balanceSpecies = "N2",
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0001,
    toleranceMoveToCore=0.05,
    maximumEdgeSpecies=100000,
    toleranceInterruptSimulation=0.05,
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,2100,'K',8),
    pressures=(0.001,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
    maximumAtoms=20,
    )
    
#uncertainty(
#	localAnalysis=True,
#	globalAnalysis=True,
#	uncorrelated=True,
#	correlated=True,
#	localNumber=10,
#	globalNumber=5,
#	pceRunTime=1800,
#	pceErrorTol=None,
#	pceMaxEvals=None,
#	)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=True,
   )
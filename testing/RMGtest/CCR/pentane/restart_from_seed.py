restartFromSeed(path='seed')

# Data sources
database(
    thermoLibraries = ['BurkeH2O2','DFT_QCI_thermo','primaryThermoLibrary'],
    reactionLibraries = ['BurkeH2O2inN2',],
    seedMechanisms = [],
    kineticsDepositories = ['training'], 
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
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
    label='Ar',
    reactive=False,
    structure=SMILES("[Ar]"),
)

# Reaction systems
simpleReactor(
    temperature=(1600,'K'),
    pressure=(10.13,'bar'),
    initialMoleFractions={
        "pentane": 1,
        "O2": 8,
        "Ar": 1,
    },
    terminationConversion={
        'O2': 0.5,
    },
    terminationTime=(1e6,'s'),
    sensitivity=['pentane'],
    sensitivityThreshold=0.01,
)

simulator(
    atol=1e-16,
    rtol=1e-8,
)

model(
    toleranceKeepInEdge=0.0001,
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=100000
)

pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,3000,'K',8),
    pressures=(0.001,100,'bar',5),
    interpolation=('Chebyshev', 6, 4),
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=True,
)

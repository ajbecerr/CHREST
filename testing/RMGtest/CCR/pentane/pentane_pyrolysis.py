database(
thermoLibraries = ["primaryThermoLibrary"],
reactionLibraries = [],
seedMechanisms = [],
kineticsDepositories = ["training"],
kineticsFamilies = "default",
kineticsEstimator = "rate rules",
)

species(
label="n-pentane",
reactive=True,
structure=SMILES("CCCCC"),
)

species(
label="Ar",
reactive=False,
structure=SMILES("[Ar]"),
)
generatedSpeciesConstraints(
    #maximumRadicalElectrons = 1,
    maximumCarbonAtoms = 7,
    )
simpleReactor(
temperature=(1000,"K"),
pressure=(400,"Pa"),
initialMoleFractions={
"n-pentane": 0.02,
"Ar": 0.98,
},
terminationConversion={
"n-pentane": 0.99,
},
terminationTime=(1e6,"s"),
)

simpleReactor(
temperature=(1500,"K"),
pressure=(400,"Pa"),
initialMoleFractions={
"n-pentane": 0.02,
"Ar": 0.98,
},
terminationConversion={
"n-pentane": 0.99,
},
terminationTime=(1e6,"s"),
)

simpleReactor(
temperature=(2000,"K"),
pressure=(400,"Pa"),
initialMoleFractions={
"n-pentane": 0.02,
"Ar": 0.98,
},
terminationConversion={
"n-pentane": 0.99,
},
terminationTime=(1e6,"s"),
)

simulator(
atol=1e-16,
rtol=1e-8,
)

model(
toleranceMoveToCore=0.02,
toleranceInterruptSimulation=0.02,
)

pressureDependence(
method="modified strong collision",
maximumGrainSize=(0.5,"kcal/mol"),
minimumNumberOfGrains=250,
temperatures=(300,3000,"K",8),
pressures=(0.001,100,"bar",5),
interpolation=("Chebyshev", 6, 4),
)

options(
    units='si',
    generateOutputHTML=True,
    generatePlots=True,
)


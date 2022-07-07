species(
    label = 'C=[C]OC[C]=CO[O](21047)',
    structure = SMILES('C=[C]OC[C]=CO[O]'),
    E0 = (469.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1670,1700,300,440,296.758,296.758,296.758,296.759],'cm^-1')),
        HinderedRotor(inertia=(0.197423,'amu*angstrom^2'), symmetry=1, barrier=(12.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197423,'amu*angstrom^2'), symmetry=1, barrier=(12.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197423,'amu*angstrom^2'), symmetry=1, barrier=(12.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197423,'amu*angstrom^2'), symmetry=1, barrier=(12.3376,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.399396,0.0774331,-9.26817e-05,5.63556e-08,-1.33811e-11,56580.4,32.9743], Tmin=(100,'K'), Tmax=(1036.25,'K')), NASAPolynomial(coeffs=[15.7335,0.0182428,-7.00292e-06,1.23508e-09,-8.31418e-14,53402.3,-41.5502], Tmin=(1036.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'ketene(1375)',
    structure = SMILES('C=C=O'),
    E0 = (-60.9858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48302,0.0083299,5.93417e-06,-1.33939e-08,5.73133e-12,-7313.53,6.51898], Tmin=(100,'K'), Tmax=(954.863,'K')), NASAPolynomial(coeffs=[5.88246,0.00580173,-1.91266e-06,3.35917e-10,-2.37398e-14,-8114.73,-6.74202], Tmin=(954.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.9858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C=CO[O](16806)',
    structure = SMILES('C=C=CO[O]'),
    E0 = (250.435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.895395,'amu*angstrom^2'), symmetry=1, barrier=(20.5869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3584.1,'J/mol'), sigma=(5.7752,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.83 K, Pc=42.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2049,0.0363638,-3.70353e-05,1.95556e-08,-4.06022e-12,30187.6,17.323], Tmin=(100,'K'), Tmax=(1179.28,'K')), NASAPolynomial(coeffs=[9.983,0.00998085,-3.47674e-06,5.84074e-10,-3.83218e-14,28353.2,-21.4843], Tmin=(1179.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = 'O2(S)(666)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,10302.3,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'C#CCO[C]=C(15760)',
    structure = SMILES('C#CCO[C]=C'),
    E0 = (329.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2175,525,1685,370,750,770,3400,2100,247.375,247.554,247.578],'cm^-1')),
        HinderedRotor(inertia=(0.530786,'amu*angstrom^2'), symmetry=1, barrier=(23.1871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900264,'amu*angstrom^2'), symmetry=1, barrier=(39.1818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.531721,'amu*angstrom^2'), symmetry=1, barrier=(23.1792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.832267,0.0556815,-5.42555e-05,2.63749e-08,-4.82793e-12,39785.3,23.9073], Tmin=(100,'K'), Tmax=(1543.81,'K')), NASAPolynomial(coeffs=[15.2665,0.00987207,-1.57431e-06,9.6601e-11,-1.05828e-15,36330.8,-48.7517], Tmin=(1543.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=CJO)"""),
)

species(
    label = 'H2CC(T)(1341)',
    structure = SMILES('[C]=C'),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9483,-0.00181651,2.07124e-05,-2.3669e-08,8.45456e-12,72206.8,5.3155], Tmin=(100,'K'), Tmax=(953.668,'K')), NASAPolynomial(coeffs=[4.20274,0.00473419,-1.57302e-06,2.859e-10,-2.08638e-14,71811.9,2.2838], Tmin=(953.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O]C[C]=CO[O](21018)',
    structure = SMILES('[O]C[C]=CO[O]'),
    E0 = (384.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,180,3343.48],'cm^-1')),
        HinderedRotor(inertia=(0.364939,'amu*angstrom^2'), symmetry=1, barrier=(8.39067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364275,'amu*angstrom^2'), symmetry=1, barrier=(8.3754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00409,0.0524605,-9.36391e-05,9.11787e-08,-3.32926e-11,46277.3,24.0788], Tmin=(100,'K'), Tmax=(871.249,'K')), NASAPolynomial(coeffs=[3.50389,0.0253808,-1.22496e-05,2.29726e-09,-1.54913e-13,46782.4,21.4483], Tmin=(871.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'O(4)',
    structure = SMILES('[O]'),
    E0 = (243.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29226.7,5.11107], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29226.7,5.11107], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=[C]OC[C]=C[O](21789)',
    structure = SMILES('C=[C]OC[C]=C[O]'),
    E0 = (334.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1670,1700,300,440,284.954,284.954,284.954,284.954,284.954],'cm^-1')),
        HinderedRotor(inertia=(0.354779,'amu*angstrom^2'), symmetry=1, barrier=(20.4426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354779,'amu*angstrom^2'), symmetry=1, barrier=(20.4426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354779,'amu*angstrom^2'), symmetry=1, barrier=(20.4426,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.798051,0.0592176,-3.71581e-05,-6.86068e-09,1.05288e-11,40333.5,28.1505], Tmin=(100,'K'), Tmax=(934.851,'K')), NASAPolynomial(coeffs=[19.1455,0.00895831,-1.83533e-06,2.68657e-10,-2.0521e-14,35668.9,-65.7305], Tmin=(934.851,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[C]=CO[O](21346)',
    structure = SMILES('[C]=CO[O]'),
    E0 = (656.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.353101,'amu*angstrom^2'), symmetry=1, barrier=(8.11848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89582,0.0288275,-5.64855e-05,5.40338e-08,-1.90645e-11,78935.3,14.7462], Tmin=(100,'K'), Tmax=(877.636,'K')), NASAPolynomial(coeffs=[5.05128,0.00860812,-4.16047e-06,7.90553e-10,-5.33044e-14,78957.3,6.90957], Tmin=(877.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]O[C]=C(2363)',
    structure = SMILES('[CH2]O[C]=C'),
    E0 = (283.793,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,607.529,608.083],'cm^-1')),
        HinderedRotor(inertia=(0.0541714,'amu*angstrom^2'), symmetry=1, barrier=(14.1917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0541016,'amu*angstrom^2'), symmetry=1, barrier=(14.1952,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.71885,0.0199755,1.11287e-05,-3.28501e-08,1.5309e-11,34186,17.0662], Tmin=(100,'K'), Tmax=(932.145,'K')), NASAPolynomial(coeffs=[10.6432,0.00688551,-1.46255e-06,2.25502e-10,-1.75111e-14,31800.1,-25.481], Tmin=(932.145,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.793,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COCJ) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OCC1=CO[O](21790)',
    structure = SMILES('C=C1OCC1=CO[O]'),
    E0 = (110.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.813278,0.0579205,-2.69312e-05,-1.56965e-08,1.27376e-11,13437.8,23.4018], Tmin=(100,'K'), Tmax=(952.958,'K')), NASAPolynomial(coeffs=[18.2427,0.0141207,-4.20134e-06,7.3189e-10,-5.36583e-14,8782.75,-66.8398], Tmin=(952.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.682,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]OCC1=COO1(21791)',
    structure = SMILES('C=[C]OCC1=COO1'),
    E0 = (259.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827228,0.0543694,-1.17035e-05,-3.35318e-08,1.90942e-11,31374.1,29.4294], Tmin=(100,'K'), Tmax=(965.073,'K')), NASAPolynomial(coeffs=[19.8474,0.0124603,-3.95714e-06,7.63337e-10,-6.01799e-14,25983.4,-70.5648], Tmin=(965.073,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OCC1=CO1(21792)',
    structure = SMILES('C=[C]OCC1=CO1'),
    E0 = (271.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.917558,0.0549065,-2.41326e-05,-2.06741e-08,1.54394e-11,32788.3,25.7875], Tmin=(100,'K'), Tmax=(936.513,'K')), NASAPolynomial(coeffs=[19.429,0.00830402,-1.48565e-06,2.18009e-10,-1.84129e-14,27897.5,-69.9059], Tmin=(936.513,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(C=CJO)"""),
)

species(
    label = 'C=COC=C=CO[O](19338)',
    structure = SMILES('C=COC=C=CO[O]'),
    E0 = (178.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.143844,0.0780796,-8.42529e-05,4.39951e-08,-8.70261e-12,21580.8,31.4349], Tmin=(100,'K'), Tmax=(1354.79,'K')), NASAPolynomial(coeffs=[20.8427,0.00999642,-2.09564e-06,2.32345e-10,-1.16848e-14,16456,-74.1126], Tmin=(1354.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]OC=C=COO(19789)',
    structure = SMILES('C=[C]OC=C=COO'),
    E0 = (265.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.262234,0.0780174,-8.75381e-05,4.89496e-08,-1.06626e-11,32111.2,32.3466], Tmin=(100,'K'), Tmax=(1127.78,'K')), NASAPolynomial(coeffs=[17.374,0.0173252,-6.81432e-06,1.23107e-09,-8.45258e-14,28251.6,-52.2655], Tmin=(1127.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.839,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=CJO)"""),
)

species(
    label = 'H(3)',
    structure = SMILES('[H]'),
    E0 = (211.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25472.7,-0.459566], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25472.7,-0.459566], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.792,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=[C]OC=C=CO[O](19790)',
    structure = SMILES('C=[C]OC=C=CO[O]'),
    E0 = (417.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,540,610,2055,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.756338,'amu*angstrom^2'), symmetry=1, barrier=(17.3897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.75701,'amu*angstrom^2'), symmetry=1, barrier=(17.4051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.756312,'amu*angstrom^2'), symmetry=1, barrier=(17.3891,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.577668,0.0738171,-9.01471e-05,5.5656e-08,-1.3389e-11,50379.5,31.974], Tmin=(100,'K'), Tmax=(1024.01,'K')), NASAPolynomial(coeffs=[15.2172,0.0166325,-6.3821e-06,1.12252e-09,-7.53738e-14,47381.3,-39.0006], Tmin=(1024.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OCC#CO[O](21793)',
    structure = SMILES('C=[C]OCC#CO[O]'),
    E0 = (456.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,2100,2250,500,550,1685,370,291.811,291.813,291.813,291.814],'cm^-1')),
        HinderedRotor(inertia=(0.25733,'amu*angstrom^2'), symmetry=1, barrier=(15.55,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257333,'amu*angstrom^2'), symmetry=1, barrier=(15.55,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.440239,'amu*angstrom^2'), symmetry=1, barrier=(26.6025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.22894,'amu*angstrom^2'), symmetry=1, barrier=(74.2601,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.535257,0.0731705,-8.75135e-05,5.32644e-08,-1.25992e-11,55075.5,31.2338], Tmin=(100,'K'), Tmax=(1043.76,'K')), NASAPolynomial(coeffs=[15.4216,0.0161226,-5.53099e-06,9.0195e-10,-5.76674e-14,51967.8,-41.2223], Tmin=(1043.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCt) + group(O2s-OsH) + group(Cs-CtOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'C#COC[C]=CO[O](21200)',
    structure = SMILES('C#COC[C]=CO[O]'),
    E0 = (438.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.87178,'amu*angstrom^2'), symmetry=1, barrier=(20.0439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.869843,'amu*angstrom^2'), symmetry=1, barrier=(19.9994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.870616,'amu*angstrom^2'), symmetry=1, barrier=(20.0172,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.869477,'amu*angstrom^2'), symmetry=1, barrier=(19.991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.140502,0.081911,-0.000104074,6.42071e-08,-1.5212e-11,52932.4,28.2303], Tmin=(100,'K'), Tmax=(1045.82,'K')), NASAPolynomial(coeffs=[18.2747,0.0125524,-4.59441e-06,7.92813e-10,-5.29961e-14,49139.3,-60.069], Tmin=(1045.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CO[O](16807)',
    structure = SMILES('C=[C][CH]O[O]'),
    E0 = (437.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,274.987],'cm^-1')),
        HinderedRotor(inertia=(0.170957,'amu*angstrom^2'), symmetry=1, barrier=(9.18865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171244,'amu*angstrom^2'), symmetry=1, barrier=(9.18809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2197,0.041819,-5.77934e-05,4.42117e-08,-1.35872e-11,52657.9,18.5046], Tmin=(100,'K'), Tmax=(795.764,'K')), NASAPolynomial(coeffs=[7.632,0.0146118,-6.50507e-06,1.24122e-09,-8.65948e-14,51796.6,-6.36981], Tmin=(795.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36761e-06,-4.93131e-09,1.45956e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.7,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00338e-07,1.59031e-10,-1.14892e-14,-1048.44,6.08305], Tmin=(1087.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2][C]=O(1376)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,623.763,623.847],'cm^-1')),
        HinderedRotor(inertia=(0.000434039,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44947,0.00889162,4.91274e-06,-1.14995e-08,4.62261e-12,19378.1,9.99239], Tmin=(100,'K'), Tmax=(1030.77,'K')), NASAPolynomial(coeffs=[6.10651,0.00625992,-2.43247e-06,4.78594e-10,-3.54632e-14,18422.4,-4.88573], Tmin=(1030.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=[C]CO[C]=C(15762)',
    structure = SMILES('[CH]=[C]CO[C]=C'),
    E0 = (648.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,352.873,352.906,352.974],'cm^-1')),
        HinderedRotor(inertia=(0.175753,'amu*angstrom^2'), symmetry=1, barrier=(15.5287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175629,'amu*angstrom^2'), symmetry=1, barrier=(15.5263,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175652,'amu*angstrom^2'), symmetry=1, barrier=(15.5289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3208.06,'J/mol'), sigma=(5.53191,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=501.09 K, Pc=43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28575,0.0546868,-5.51644e-05,2.8389e-08,-5.73388e-12,78126.7,25.5344], Tmin=(100,'K'), Tmax=(1210.97,'K')), NASAPolynomial(coeffs=[13.4853,0.01439,-5.2497e-06,9.09901e-10,-6.09175e-14,75172,-35.6567], Tmin=(1210.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OC=[C][CH]O[O](19792)',
    structure = SMILES('C=[C]O[CH][C]=CO[O]'),
    E0 = (580.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1670,1700,300,440,324.022,324.06,324.067,324.085],'cm^-1')),
        HinderedRotor(inertia=(0.00160538,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166305,'amu*angstrom^2'), symmetry=1, barrier=(12.395,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395431,'amu*angstrom^2'), symmetry=1, barrier=(29.4704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166377,'amu*angstrom^2'), symmetry=1, barrier=(12.395,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.705759,0.073884,-9.31857e-05,6.08554e-08,-1.56935e-11,69909.5,32.8151], Tmin=(100,'K'), Tmax=(950.461,'K')), NASAPolynomial(coeffs=[13.4227,0.020365,-8.72331e-06,1.6125e-09,-1.10865e-13,67492.1,-27.8908], Tmin=(950.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC[C]=[C]O[O](21794)',
    structure = SMILES('C=[C]OC[C]=[C]O[O]'),
    E0 = (709.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,273.345,273.492,273.724,3207.87],'cm^-1')),
        HinderedRotor(inertia=(0.244816,'amu*angstrom^2'), symmetry=1, barrier=(12.9765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.24469,'amu*angstrom^2'), symmetry=1, barrier=(12.98,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244487,'amu*angstrom^2'), symmetry=1, barrier=(12.9793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244979,'amu*angstrom^2'), symmetry=1, barrier=(12.9778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.60859,0.0792058,-0.000119692,9.55207e-08,-2.99114e-11,85401.4,35.352], Tmin=(100,'K'), Tmax=(854.942,'K')), NASAPolynomial(coeffs=[11.4108,0.0226856,-1.00356e-05,1.83123e-09,-1.22588e-13,83772.9,-13.7916], Tmin=(854.942,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(709.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C#COC[C][CH]O[O](20836)',
    structure = SMILES('C#COC[C][CH]O[O]'),
    E0 = (700.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.213438,0.0960346,-0.000141243,1.02438e-07,-2.88156e-11,84401.7,27.7941], Tmin=(100,'K'), Tmax=(878.626,'K')), NASAPolynomial(coeffs=[16.8922,0.018155,-8.27783e-06,1.54272e-09,-1.05376e-13,81396,-52.5161], Tmin=(878.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(CCsJOOH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=[C]OC=C[CH]O[O](19794)',
    structure = SMILES('C=[C]O[CH]C=CO[O]'),
    E0 = (342.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.590974,0.0710983,-7.50214e-05,4.04641e-08,-8.59866e-12,41312.4,32.7477], Tmin=(100,'K'), Tmax=(1149.05,'K')), NASAPolynomial(coeffs=[15.4169,0.0194872,-7.64674e-06,1.37398e-09,-9.3766e-14,37905.2,-40.8386], Tmin=(1149.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJ(O)C) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OCC=[C]O[O](21795)',
    structure = SMILES('C=[C]OCC=[C]O[O]'),
    E0 = (471.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.757914,0.0731459,-8.93517e-05,5.8333e-08,-1.52097e-11,56793.1,34.3482], Tmin=(100,'K'), Tmax=(936.216,'K')), NASAPolynomial(coeffs=[12.3488,0.023627,-1.00182e-05,1.84469e-09,-1.26502e-13,54622.7,-20.8082], Tmin=(936.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC[C]=CO[O](21135)',
    structure = SMILES('[CH]=COC[C]=CO[O]'),
    E0 = (476.697,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.748832,'amu*angstrom^2'), symmetry=1, barrier=(17.2171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.748274,'amu*angstrom^2'), symmetry=1, barrier=(17.2043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.749471,'amu*angstrom^2'), symmetry=1, barrier=(17.2318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.748499,'amu*angstrom^2'), symmetry=1, barrier=(17.2095,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.244014,0.0834878,-9.76658e-05,5.49579e-08,-1.17316e-11,57494.7,32.5542], Tmin=(100,'K'), Tmax=(1222.67,'K')), NASAPolynomial(coeffs=[21.2054,0.00943736,-2.06124e-06,2.34824e-10,-1.18865e-14,52539.4,-74.0533], Tmin=(1222.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.697,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=COC=[C][CH]O[O](19799)',
    structure = SMILES('C=CO[CH][C]=CO[O]'),
    E0 = (340.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,324.233,324.546,324.686,324.724],'cm^-1')),
        HinderedRotor(inertia=(0.249236,'amu*angstrom^2'), symmetry=1, barrier=(18.6304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249321,'amu*angstrom^2'), symmetry=1, barrier=(18.6299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248925,'amu*angstrom^2'), symmetry=1, barrier=(18.626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248806,'amu*angstrom^2'), symmetry=1, barrier=(18.623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.267037,0.0749271,-7.65563e-05,3.59832e-08,-5.6645e-12,41098.2,31.2537], Tmin=(100,'K'), Tmax=(993.565,'K')), NASAPolynomial(coeffs=[18.3112,0.0149472,-5.12311e-06,8.81418e-10,-6.0175e-14,36887.5,-58.8281], Tmin=(993.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(340.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJ(O)C) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]OC[C]=[C]OO(21796)',
    structure = SMILES('C=[C]OC[C]=[C]OO'),
    E0 = (557.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.492071,0.0809229,-0.00010777,7.58656e-08,-2.12809e-11,67124.7,35.0204], Tmin=(100,'K'), Tmax=(872.172,'K')), NASAPolynomial(coeffs=[12.7186,0.0248473,-1.13259e-05,2.14443e-09,-1.48839e-13,64992.1,-22.2932], Tmin=(872.172,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OCC=CO[O](10675)',
    structure = SMILES('[CH]=[C]OCC=CO[O]'),
    E0 = (478.599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,229.122,229.123,229.123],'cm^-1')),
        HinderedRotor(inertia=(0.364336,'amu*angstrom^2'), symmetry=1, barrier=(13.5726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364336,'amu*angstrom^2'), symmetry=1, barrier=(13.5726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364336,'amu*angstrom^2'), symmetry=1, barrier=(13.5726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364336,'amu*angstrom^2'), symmetry=1, barrier=(13.5726,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.224687,0.0779635,-9.02814e-05,5.2014e-08,-1.15854e-11,57702.5,33.5284], Tmin=(100,'K'), Tmax=(1109.09,'K')), NASAPolynomial(coeffs=[17.668,0.0150534,-5.19805e-06,8.70938e-10,-5.73037e-14,53833.3,-52.4314], Tmin=(1109.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=COC[C]=[C]O[O](21797)',
    structure = SMILES('C=COC[C]=[C]O[O]'),
    E0 = (469.345,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1670,1700,300,440,296.758,296.758,296.758,296.759],'cm^-1')),
        HinderedRotor(inertia=(0.197423,'amu*angstrom^2'), symmetry=1, barrier=(12.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197423,'amu*angstrom^2'), symmetry=1, barrier=(12.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197423,'amu*angstrom^2'), symmetry=1, barrier=(12.3376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197423,'amu*angstrom^2'), symmetry=1, barrier=(12.3376,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.399396,0.0774331,-9.26817e-05,5.63556e-08,-1.33811e-11,56580.4,32.9743], Tmin=(100,'K'), Tmax=(1036.25,'K')), NASAPolynomial(coeffs=[15.7335,0.0182428,-7.00292e-06,1.23508e-09,-8.31418e-14,53402.3,-41.5502], Tmin=(1036.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC=[C][CH]OO(19797)',
    structure = SMILES('C=[C]O[CH][C]=COO'),
    E0 = (428.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.457341,0.0773134,-8.79758e-05,5.09122e-08,-1.16436e-11,51638.3,32.9461], Tmin=(100,'K'), Tmax=(1066.51,'K')), NASAPolynomial(coeffs=[15.3635,0.021409,-9.3511e-06,1.76601e-09,-1.23668e-13,48458.7,-39.9276], Tmin=(1066.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(428.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OC[C]=COO(21798)',
    structure = SMILES('[CH]=[C]OC[C]=COO'),
    E0 = (564.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0811323,0.0842997,-0.000103669,6.30211e-08,-1.48618e-11,68028.8,33.762], Tmin=(100,'K'), Tmax=(1045.28,'K')), NASAPolynomial(coeffs=[17.7207,0.0167974,-6.80091e-06,1.23917e-09,-8.52414e-14,64341.2,-52.1196], Tmin=(1045.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,-745,4.3663], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,-745,4.3663], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.19426,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'Ne',
    structure = SMILES('[Ne]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (20.1797,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1235.53,'J/mol'), sigma=(3.758e-10,'m'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with fixed Lennard Jones Parameters. This is the fallback method! Try improving transport databases!"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,3.35532], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ne""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'N2',
    structure = SMILES('N#N'),
    E0 = (-8.69489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61263,-0.00100893,2.49899e-06,-1.43376e-09,2.58636e-13,-1051.1,2.6527], Tmin=(100,'K'), Tmax=(1817.03,'K')), NASAPolynomial(coeffs=[2.9759,0.00164142,-7.19726e-07,1.25378e-10,-7.91532e-15,-1025.84,5.5376], Tmin=(1817.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.69489,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: BurkeH2O2"""),
)

transitionState(
    label = 'TS1',
    E0 = (469.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (469.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (984.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (577.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (995.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (477.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (477.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (600.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (562.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (494.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (642.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (683.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (665.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (469.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (469.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (598.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (640.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (792.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (920.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (912.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (608.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (703.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (582.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (611.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (694.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (787.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (778.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (534.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (715.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['ketene(1375)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['O2(S)(666)', 'C#CCO[C]=C(15760)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H2CC(T)(1341)', '[O]C[C]=CO[O](21018)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', 'C=[C]OC[C]=C[O](21789)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[C]=CO[O](21346)', '[CH2]O[C]=C(2363)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['C=C1OCC1=CO[O](21790)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['C=[C]OCC1=COO1(21791)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['O(4)', 'C=[C]OCC1=CO1(21792)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO;Y_rad_intra;OOJ] for rate rule [R2OO_D;Cd_rad_in/NonDeC;OOJ]
Euclidian distance = 3.1622776601683795
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['C=COC=C=CO[O](19338)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.08e+10,'s^-1'), n=-0.305, Ea=(93.094,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3;Y_rad_De;XH_Rrad_De] for rate rule [R3radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['C=[C]OC=C=COO(19789)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', 'C=[C]OC=C=CO[O](19790)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', 'C=[C]OCC#CO[O](21793)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', 'C#COC[C]=CO[O](21200)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['ketene(1375)', '[CH2][C]=CO[O](16807)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(93.0203,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond
Ea raised from 87.8 to 93.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['O2(2)', 'C#CCO[C]=C(15760)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.3016e+08,'cm^3/(mol*s)'), n=1.628, Ea=(148.22,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-Cs;OJ] for rate rule [Ct-H_Ct-Cs;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -1.9 to 148.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]=O(1376)', '[CH2][C]=CO[O](16807)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', '[CH]=[C]CO[C]=C(15762)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.12762e+06,'m^3/(mol*s)'), n=-0.0535546, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0261245836817, var=0.43760691139, Tref=1000.0, N=5, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R
    Total Standard Deviation in ln(k): 1.39180927994
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', 'C=[C]OC=[C][CH]O[O](19792)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', 'C=[C]OC[C]=[C]O[O](21794)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', 'C#COC[C][CH]O[O](20836)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['C=[C]OC=C[CH]O[O](19794)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['C=[C]OCC=[C]O[O](21795)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=COC[C]=CO[O](21135)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['C=COC=[C][CH]O[O](19799)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=[C]OC[C]=[C]OO(21796)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.58963e+06,'s^-1'), n=1.78157, Ea=(137.44,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]OCC=CO[O](10675)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=COC[C]=[C]O[O](21797)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_single;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out_singleNd;Cd_H_out_doubleC]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]OC[C]=CO[O](21047)'],
    products = ['C=[C]OC=[C][CH]OO(19797)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(23000,'s^-1'), n=2.11, Ea=(64.7265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]OC[C]=COO(21798)'],
    products = ['C=[C]OC[C]=CO[O](21047)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R8Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4727',
    isomers = [
        'C=[C]OC[C]=CO[O](21047)',
    ],
    reactants = [
        ('ketene(1375)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4727',
    Tmin = (300,'K'),
    Tmax = (3000,'K'),
    Tcount = 8,
    Tlist = ([302.617,324.619,374.997,470.374,649.057,1000.02,1706.11,2761.25],'K'),
    Pmin = (0.001,'bar'),
    Pmax = (100,'bar'),
    Pcount = 5,
    Plist = ([0.00132544,0.0107284,0.316228,9.32101,75.4469],'bar'),
    maximumGrainSize = (0.5,'kcal/mol'),
    minimumGrainCount = 250,
    method = 'modified strong collision',
    interpolationModel = ('Chebyshev', 6, 4),
    activeKRotor = True,
    activeJRotor = True,
    rmgmode = True,
)


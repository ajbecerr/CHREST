species(
    label = '[O]C=C([O])C[C]=O(11162)',
    structure = SMILES('[O]C=C([O])C[C]=O'),
    E0 = (-77.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,338.167,338.779,339.1],'cm^-1')),
        HinderedRotor(inertia=(0.217869,'amu*angstrom^2'), symmetry=1, barrier=(17.9781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218111,'amu*angstrom^2'), symmetry=1, barrier=(17.9671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24419,0.0477059,-2.08267e-05,-1.8493e-08,1.30808e-11,-9181.61,25.2029], Tmin=(100,'K'), Tmax=(971.866,'K')), NASAPolynomial(coeffs=[19.0515,0.00477417,-1.42205e-06,3.38529e-10,-3.1611e-14,-14076.6,-67.575], Tmin=(971.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCCJ=O)"""),
)

species(
    label = 'OCHCO(3676)',
    structure = SMILES('O=[C]C=O'),
    E0 = (-75.5464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,180,525.376,1512.41,1512.65,1513.37],'cm^-1')),
        HinderedRotor(inertia=(0.00619061,'amu*angstrom^2'), symmetry=1, barrier=(0.260399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25386,0.0154074,-1.14326e-05,4.29104e-09,-6.58698e-13,-9058.53,11.1539], Tmin=(100,'K'), Tmax=(1504.29,'K')), NASAPolynomial(coeffs=[6.43153,0.00695778,-3.00696e-06,5.56991e-10,-3.81281e-14,-10014.6,-5.47402], Tmin=(1504.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=C([O])C([O])=C[O](11158)',
    structure = SMILES('C=C([O])C([O])=C[O]'),
    E0 = (-125.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0729095,0.0724912,-8.7986e-05,4.83801e-08,-9.6461e-12,-14943,25.3649], Tmin=(100,'K'), Tmax=(1456.33,'K')), NASAPolynomial(coeffs=[21.1206,-0.00117263,3.8032e-06,-9.24506e-10,6.85246e-14,-19477.2,-79.2224], Tmin=(1456.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[C]=O(2355)',
    structure = SMILES('[C]=O'),
    E0 = (438.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3323.79],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09152,0.00193295,-1.59153e-05,2.47563e-08,-1.11287e-11,52770.5,4.46624], Tmin=(100,'K'), Tmax=(866.029,'K')), NASAPolynomial(coeffs=[1.05092,0.00526657,-3.13864e-06,6.40624e-10,-4.47939e-14,53698.8,21.0169], Tmin=(866.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'C=C([O])[CH][O](2850)',
    structure = SMILES('[CH2]C([O])=C[O]'),
    E0 = (20.9566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,217.215,217.577],'cm^-1')),
        HinderedRotor(inertia=(0.665078,'amu*angstrom^2'), symmetry=1, barrier=(22.287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00575,0.032418,-6.0508e-07,-3.82059e-08,2.20969e-11,2603.17,17.8437], Tmin=(100,'K'), Tmax=(887.044,'K')), NASAPolynomial(coeffs=[16.9262,-0.00266727,4.27963e-06,-9.58521e-10,6.70145e-14,-1310.55,-59.4906], Tmin=(887.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.9566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=C(O)CJ)"""),
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
    label = '[O]C=[C]C[C]=O(9880)',
    structure = SMILES('[O]C=[C]C[C]=O'),
    E0 = (237.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,229.413,229.899],'cm^-1')),
        HinderedRotor(inertia=(0.549517,'amu*angstrom^2'), symmetry=1, barrier=(20.398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.544492,'amu*angstrom^2'), symmetry=1, barrier=(20.3927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8503,0.0374587,-1.49937e-05,-1.21738e-08,8.00458e-12,28630.6,21.945], Tmin=(100,'K'), Tmax=(1027.37,'K')), NASAPolynomial(coeffs=[14.6713,0.00877593,-4.11938e-06,8.88192e-10,-6.9517e-14,24875.6,-45.709], Tmin=(1027.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C([O])C[C]=O(2677)',
    structure = SMILES('[CH]=C([O])C[C]=O'),
    E0 = (237.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1855,455,950,3120,650,792.5,1650,457.716],'cm^-1')),
        HinderedRotor(inertia=(0.0811444,'amu*angstrom^2'), symmetry=1, barrier=(11.9485,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0819876,'amu*angstrom^2'), symmetry=1, barrier=(12.0127,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3870.22,'J/mol'), sigma=(6.17892,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.52 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73346,0.0430516,-3.78896e-05,1.47612e-08,-1.97815e-12,28611.6,22.588], Tmin=(100,'K'), Tmax=(1203.61,'K')), NASAPolynomial(coeffs=[13.8763,0.00944907,-4.42733e-06,8.87773e-10,-6.46617e-14,25199.4,-40.2764], Tmin=(1203.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[O]C=C1CC(=O)O1(11816)',
    structure = SMILES('[O]C=C1CC(=O)O1'),
    E0 = (-312.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70382,0.0279883,5.14183e-05,-1.03518e-07,4.63734e-11,-37435.9,21.4778], Tmin=(100,'K'), Tmax=(923.172,'K')), NASAPolynomial(coeffs=[21.7966,-0.000855149,3.69222e-06,-7.43295e-10,4.29234e-14,-43626.5,-87.288], Tmin=(923.172,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-312.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + ring(4-Methylene-2-oxetanone) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]CC1=COO1(11817)',
    structure = SMILES('O=[C]CC1=COO1'),
    E0 = (162.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90375,0.0322893,1.17627e-05,-4.08215e-08,1.75547e-11,19670.2,23.1389], Tmin=(100,'K'), Tmax=(1029.84,'K')), NASAPolynomial(coeffs=[15.4189,0.0122108,-6.20672e-06,1.37546e-09,-1.08606e-13,15167.6,-50.8067], Tmin=(1029.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1=COC(=O)C1(11818)',
    structure = SMILES('[O]C1=COC(=O)C1'),
    E0 = (-356.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45301,0.00747765,0.000100194,-1.44442e-07,5.75343e-11,-42854.4,17.8733], Tmin=(100,'K'), Tmax=(951.494,'K')), NASAPolynomial(coeffs=[19.4853,0.00335783,3.04212e-07,8.38653e-11,-2.35668e-14,-49150.4,-79.5032], Tmin=(951.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-356.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C=C(O)C=C=O(11819)',
    structure = SMILES('[O]C=C(O)C=C=O'),
    E0 = (-260.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0154503,0.0695631,-8.05317e-05,4.21143e-08,-8.0157e-12,-31189.6,24.7627], Tmin=(100,'K'), Tmax=(1516.41,'K')), NASAPolynomial(coeffs=[21.335,-0.00080472,2.97186e-06,-7.06486e-10,5.11201e-14,-36049.5,-81.8038], Tmin=(1516.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-260.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]CC(O)=C=O(11820)',
    structure = SMILES('O=[C]CC(O)=C=O'),
    E0 = (-183.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.975557,0.061689,-6.90927e-05,3.27154e-08,-4.21888e-12,-21936.9,24.1894], Tmin=(100,'K'), Tmax=(906.846,'K')), NASAPolynomial(coeffs=[16.7336,0.00641745,-1.21558e-06,1.26038e-10,-6.74914e-15,-25380.2,-53.5206], Tmin=(906.846,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-183.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-OdCsH) + missing(Cdd-CdO2d) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(C=C=O)=CO(11821)',
    structure = SMILES('[O]C(C=C=O)=CO'),
    E0 = (-264.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0746075,0.0711196,-8.71351e-05,4.84059e-08,-9.78638e-12,-31635.8,24.5201], Tmin=(100,'K'), Tmax=(1423.53,'K')), NASAPolynomial(coeffs=[20.7844,-0.000934215,3.39493e-06,-8.31092e-10,6.18645e-14,-36127.5,-77.7731], Tmin=(1423.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-264.319,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C(=C=O)CC=O(11822)',
    structure = SMILES('O=[C]C(=O)CC=O'),
    E0 = (-220.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68385,0.0552103,-7.32913e-05,5.58911e-08,-1.76511e-11,-26487,23.9151], Tmin=(100,'K'), Tmax=(766.232,'K')), NASAPolynomial(coeffs=[7.71155,0.0237425,-1.16868e-05,2.28974e-09,-1.61888e-13,-27410.7,-3.55989], Tmin=(766.232,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-220.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = '[O][CH]C([O])[CH][C]=O(11823)',
    structure = SMILES('[O][CH]C([O])[CH][C]=O'),
    E0 = (436.825,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,1855,455,950,217.56,217.593,1770.49,1770.49],'cm^-1')),
        HinderedRotor(inertia=(0.00356154,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148565,'amu*angstrom^2'), symmetry=1, barrier=(4.99199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148632,'amu*angstrom^2'), symmetry=1, barrier=(4.99201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08904,0.0757754,-0.000143319,1.36777e-07,-4.86444e-11,52631.5,30.6782], Tmin=(100,'K'), Tmax=(875.492,'K')), NASAPolynomial(coeffs=[5.69292,0.0274747,-1.38485e-05,2.61558e-09,-1.75971e-13,52870.4,15.0471], Tmin=(875.492,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.825,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCJCO) + radical(CCsJOH) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1O[C]1C[C]=O(11824)',
    structure = SMILES('[O]C1O[C]1C[C]=O'),
    E0 = (155.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54656,0.0623681,-0.000105033,9.41407e-08,-3.11796e-11,18786.4,25.1004], Tmin=(100,'K'), Tmax=(950.197,'K')), NASAPolynomial(coeffs=[4.99118,0.0248198,-9.37523e-06,1.49951e-09,-8.92093e-14,19172.3,14.1331], Tmin=(950.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(155.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1(C[C]=O)[CH]O1(11825)',
    structure = SMILES('[O]C1(C[C]=O)[CH]O1'),
    E0 = (153.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538246,0.065494,-7.94957e-05,4.67175e-08,-1.00568e-11,18595.2,25.5637], Tmin=(100,'K'), Tmax=(1357.34,'K')), NASAPolynomial(coeffs=[15.945,0.00616031,1.46927e-06,-6.10434e-10,5.29224e-14,15696.1,-48.7449], Tmin=(1357.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CCsJO) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1CC(=O)C1[O](11826)',
    structure = SMILES('[O][C]1CC(=O)C1[O]'),
    E0 = (192.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06288,0.0320229,5.26158e-06,-3.19207e-08,1.48459e-11,23224.4,24.9021], Tmin=(100,'K'), Tmax=(984.823,'K')), NASAPolynomial(coeffs=[12.7262,0.0137913,-5.16764e-06,9.97101e-10,-7.44584e-14,19907.9,-32.5535], Tmin=(984.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[O][CH]C1([O])CC1=O(11827)',
    structure = SMILES('[O][CH]C1([O])CC1=O'),
    E0 = (286.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15241,0.070336,-0.000119576,1.05788e-07,-3.55655e-11,34513.1,24.3311], Tmin=(100,'K'), Tmax=(888.436,'K')), NASAPolynomial(coeffs=[7.75515,0.0231882,-1.05617e-05,1.91612e-09,-1.26198e-13,34027.4,-2.87306], Tmin=(888.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CC(C)(C=O)OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'CO(2039)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(815.652,'J/mol'), sigma=(3.65,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.95,'angstroms^3'), rotrelaxcollnum=1.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35309e-07,1.51269e-10,-9.88873e-15,-14292.7,6.51158], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
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
    label = '[O]C=C([O])C=C=O(11380)',
    structure = SMILES('[O]C=C([O])C=C=O'),
    E0 = (-122.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10808,'amu*angstrom^2'), symmetry=1, barrier=(25.4769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.729912,0.0612206,-7.32087e-05,4.05217e-08,-8.30047e-12,-14649.2,23.5818], Tmin=(100,'K'), Tmax=(1359.99,'K')), NASAPolynomial(coeffs=[18.0567,0.00205243,1.10232e-06,-3.4264e-10,2.70638e-14,-18603.2,-62.5473], Tmin=(1359.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(=C=O)C[C]=O(11828)',
    structure = SMILES('O=[C]CC(=O)[C]=O'),
    E0 = (-60.8595,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1850,1860,440,470,900,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.478851,'amu*angstrom^2'), symmetry=1, barrier=(11.0097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.478701,'amu*angstrom^2'), symmetry=1, barrier=(11.0063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.478388,'amu*angstrom^2'), symmetry=1, barrier=(10.9991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17211,0.0701346,-0.000125684,1.12856e-07,-3.84676e-11,-7225.45,24.3442], Tmin=(100,'K'), Tmax=(867.39,'K')), NASAPolynomial(coeffs=[8.47541,0.019855,-1.0027e-05,1.89959e-09,-1.28288e-13,-7867.94,-6.25129], Tmin=(867.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.8595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + group(Cds-O2d(Cds-O2d)H) + radical(C=OCCJ=O) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]=C[O](9592)',
    structure = SMILES('[O][CH][C]=O'),
    E0 = (204.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,343.122],'cm^-1')),
        HinderedRotor(inertia=(0.528453,'amu*angstrom^2'), symmetry=1, barrier=(44.0956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.99411,0.0153541,4.13359e-06,-1.97935e-08,9.27968e-12,24622.8,13.8337], Tmin=(100,'K'), Tmax=(987.313,'K')), NASAPolynomial(coeffs=[10.2587,0.00223056,-7.0477e-07,2.03549e-10,-2.00526e-14,22393.5,-25.1462], Tmin=(987.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=OCOJ) + radical(OCJC=O) + radical(CsCJ=O)"""),
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
    label = '[O][C]=CC([O])=C[O](11388)',
    structure = SMILES('[O]C=C([O])[CH][C]=O'),
    E0 = (122.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,380.143,380.433,380.493],'cm^-1')),
        HinderedRotor(inertia=(0.00116362,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20124,'amu*angstrom^2'), symmetry=1, barrier=(20.6549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43853,0.0424997,-1.2256e-05,-2.78613e-08,1.68493e-11,14855,27.0016], Tmin=(100,'K'), Tmax=(957.826,'K')), NASAPolynomial(coeffs=[19.5475,0.000592525,5.68507e-07,-3.47189e-11,-6.39784e-15,9839.26,-67.6573], Tmin=(957.826,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]=C([O])C[C]=O(11829)',
    structure = SMILES('[O][C]=C([O])C[C]=O'),
    E0 = (162.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1855,455,950,1685,370,377.747,378.659,378.953],'cm^-1')),
        HinderedRotor(inertia=(0.0011805,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14693,'amu*angstrom^2'), symmetry=1, barrier=(14.8814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32315,0.0508629,-5.18798e-05,2.45824e-08,-4.4553e-12,19645.2,28.057], Tmin=(100,'K'), Tmax=(1359.44,'K')), NASAPolynomial(coeffs=[16.3931,0.00652122,-2.9534e-06,5.88978e-10,-4.29284e-14,15547.9,-49.2744], Tmin=(1359.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCCJ=O) + radical(C=CJO)"""),
)

species(
    label = '[O]C=CC([O])=C[O](11298)',
    structure = SMILES('[O]C=CC([O])=C[O]'),
    E0 = (-116.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.729844,0.0535171,-1.11373e-05,-5.09621e-08,3.15912e-11,-13915.9,23.0184], Tmin=(100,'K'), Tmax=(895.886,'K')), NASAPolynomial(coeffs=[26.111,-0.00884055,7.9371e-06,-1.65676e-09,1.12664e-13,-20508.9,-108.055], Tmin=(895.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-116.827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C(=[C]O)C[C]=O(11830)',
    structure = SMILES('[O]C(=[C]O)C[C]=O'),
    E0 = (21.0196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,1855,455,950,345.862,345.868],'cm^-1')),
        HinderedRotor(inertia=(0.16334,'amu*angstrom^2'), symmetry=1, barrier=(13.8649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163342,'amu*angstrom^2'), symmetry=1, barrier=(13.8649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163327,'amu*angstrom^2'), symmetry=1, barrier=(13.8649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02888,0.0566005,-5.17188e-05,1.48911e-08,1.24373e-12,2642.89,27.6938], Tmin=(100,'K'), Tmax=(979.793,'K')), NASAPolynomial(coeffs=[18.2492,0.00506838,-1.56193e-06,3.1571e-10,-2.60955e-14,-1632.53,-59.6305], Tmin=(979.793,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.0196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=CJO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=C(O)[CH][C]=O(11831)',
    structure = SMILES('[O]C=C(O)[CH][C]=O'),
    E0 = (-15.1647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,354.262,354.993],'cm^-1')),
        HinderedRotor(inertia=(0.220581,'amu*angstrom^2'), symmetry=1, barrier=(19.7456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221115,'amu*angstrom^2'), symmetry=1, barrier=(19.7607,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.221627,'amu*angstrom^2'), symmetry=1, barrier=(19.7419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07843,0.0463968,-4.43722e-06,-4.54194e-08,2.51232e-11,-1702.33,26.7932], Tmin=(100,'K'), Tmax=(945.573,'K')), NASAPolynomial(coeffs=[23.2398,-0.00268921,2.58106e-06,-4.16316e-10,1.82003e-14,-7890,-89.4401], Tmin=(945.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-15.1647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]=C(O)C[C]=O(11832)',
    structure = SMILES('O=[C]C[C](O)[C]=O'),
    E0 = (16.7646,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,1850,1860,440,470,900,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.203228,'amu*angstrom^2'), symmetry=1, barrier=(4.6726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202966,'amu*angstrom^2'), symmetry=1, barrier=(4.66658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20294,'amu*angstrom^2'), symmetry=1, barrier=(4.66599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202911,'amu*angstrom^2'), symmetry=1, barrier=(4.66532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14341,0.0680473,-0.000111172,9.30476e-08,-3.0092e-11,2114.21,29.7969], Tmin=(100,'K'), Tmax=(864.759,'K')), NASAPolynomial(coeffs=[10.0074,0.0179038,-8.33492e-06,1.54168e-09,-1.03317e-13,923.025,-9.70204], Tmin=(864.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.7646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(CCCJ=O) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]=C([O])CC=O(11833)',
    structure = SMILES('[O][C]=C([O])CC=O'),
    E0 = (2.52164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2782.5,750,1395,475,1775,1000,1685,370,572.052,572.176,572.364],'cm^-1')),
        HinderedRotor(inertia=(0.0605313,'amu*angstrom^2'), symmetry=1, barrier=(14.0637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0605645,'amu*angstrom^2'), symmetry=1, barrier=(14.0681,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47806,0.0469749,-3.38107e-05,6.66205e-09,1.22777e-12,401.276,26.6794], Tmin=(100,'K'), Tmax=(1110.47,'K')), NASAPolynomial(coeffs=[14.966,0.0123392,-5.86745e-06,1.19817e-09,-8.87812e-14,-3454.37,-43.6781], Tmin=(1110.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.52164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C([CH][C]=O)=CO(11834)',
    structure = SMILES('[O]C([CH]O)=C[C]=O'),
    E0 = (-69.3373,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,322.886,323.079],'cm^-1')),
        HinderedRotor(inertia=(0.219257,'amu*angstrom^2'), symmetry=1, barrier=(16.1945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.219108,'amu*angstrom^2'), symmetry=1, barrier=(16.1947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218982,'amu*angstrom^2'), symmetry=1, barrier=(16.196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09142,0.0638388,-7.78407e-05,4.60078e-08,-1.0592e-11,-8234.39,26.1226], Tmin=(100,'K'), Tmax=(1065.18,'K')), NASAPolynomial(coeffs=[14.6085,0.013079,-6.35973e-06,1.26967e-09,-9.18479e-14,-11114,-39.9427], Tmin=(1065.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-69.3373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(C=CCJ=O)"""),
)

species(
    label = 'O=[C]CC(=O)C=O(11154)',
    structure = SMILES('O=[C]CC(=O)C=O'),
    E0 = (-220.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26473,0.0669541,-0.000110029,9.81424e-08,-3.42775e-11,-26466.6,23.1924], Tmin=(100,'K'), Tmax=(827.885,'K')), NASAPolynomial(coeffs=[7.504,0.0249903,-1.25841e-05,2.43063e-09,-1.68035e-13,-27094.6,-3.2839], Tmin=(827.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-220.82,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + group(Cds-O2d(Cds-O2d)H) + radical(C=OCCJ=O)"""),
)

species(
    label = '[O]C([C]=O)C[C]=O(11149)',
    structure = SMILES('[O]C([C]=O)C[C]=O'),
    E0 = (83.8699,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1850,1860,440,470,900,1000,237.783,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00298127,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249441,'amu*angstrom^2'), symmetry=1, barrier=(10.0051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.249547,'amu*angstrom^2'), symmetry=1, barrier=(10.0042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5941,0.0556513,-7.49427e-05,5.44764e-08,-1.584e-11,10171.4,29.5406], Tmin=(100,'K'), Tmax=(841.41,'K')), NASAPolynomial(coeffs=[9.4815,0.0181555,-8.09868e-06,1.51499e-09,-1.04211e-13,8844.05,-7.14967], Tmin=(841.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.8699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH][O](1548)',
    structure = SMILES('[CH][O]'),
    E0 = (431.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([806.798,806.798,2696.84],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86882,-0.000695042,1.53508e-05,-2.00576e-08,7.65929e-12,51865.1,6.76618], Tmin=(100,'K'), Tmax=(973.285,'K')), NASAPolynomial(coeffs=[6.03251,-0.000301449,4.32946e-07,-3.67045e-11,-1.25373e-15,51004.1,-5.87326], Tmin=(973.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(H3COJ) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = 'O=[C]C[C]=O(2682)',
    structure = SMILES('O=[C]C[C]=O'),
    E0 = (47.1985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1850,1860,440,470,900,1000],'cm^-1')),
        HinderedRotor(inertia=(0.556552,'amu*angstrom^2'), symmetry=1, barrier=(12.7962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.55598,'amu*angstrom^2'), symmetry=1, barrier=(12.7831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96235,0.0530441,-0.000104388,9.76067e-08,-3.34753e-11,5742.3,16.165], Tmin=(100,'K'), Tmax=(905.297,'K')), NASAPolynomial(coeffs=[6.1662,0.0143356,-6.89083e-06,1.24225e-09,-7.98349e-14,5806.21,0.858947], Tmin=(905.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.1985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCCJ=O) + radical(C=OCCJ=O)"""),
)

species(
    label = '[O]C1C(=O)CC1=O(11165)',
    structure = SMILES('[O]C1C(=O)CC1=O'),
    E0 = (-129.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26986,0.0294484,3.67018e-06,-2.41288e-08,1.05904e-11,-15453.7,22.1265], Tmin=(100,'K'), Tmax=(1028.96,'K')), NASAPolynomial(coeffs=[10.4233,0.0173242,-7.18614e-06,1.39007e-09,-1.0096e-13,-18167.6,-22.4763], Tmin=(1028.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-129.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(C=OCOJ)"""),
)

species(
    label = '[O]CC(=O)C=C=O(11835)',
    structure = SMILES('[O]CC(=O)C=C=O'),
    E0 = (-109.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28204,0.04581,-3.13634e-05,-3.55203e-08,5.42643e-11,-13083.3,21.8162], Tmin=(100,'K'), Tmax=(468.372,'K')), NASAPolynomial(coeffs=[6.2237,0.0261576,-1.32941e-05,2.62577e-09,-1.85695e-13,-13606.2,4.1492], Tmin=(468.372,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(C=OCOJ)"""),
)

species(
    label = '[O][CH][C]1CC(=O)O1(11836)',
    structure = SMILES('[O]C=C1C[C]([O])O1'),
    E0 = (108.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56939,0.0352037,2.52036e-05,-7.61981e-08,3.74659e-11,13112.5,21.0056], Tmin=(100,'K'), Tmax=(903.128,'K')), NASAPolynomial(coeffs=[20.8176,-0.00118555,4.48746e-06,-9.99525e-10,6.67338e-14,7643.05,-80.927], Tmin=(903.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(CCOJ) + radical(C=COJ) + radical(Cs_P)"""),
)

species(
    label = 'O=[C]C[C]1[CH]OO1(11837)',
    structure = SMILES('O=[C]C[C]1[CH]OO1'),
    E0 = (375.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60947,0.064332,-0.000118408,1.19727e-07,-4.51905e-11,45215,21.674], Tmin=(100,'K'), Tmax=(858.072,'K')), NASAPolynomial(coeffs=[2.01768,0.0344981,-1.74291e-05,3.33777e-09,-2.28181e-13,46173.2,25.7588], Tmin=(858.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxetane) + radical(C2CsJOO) + radical(CCsJOO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]CC(=O)[CH][C]=O(11838)',
    structure = SMILES('[O]CC([O])=C[C]=O'),
    E0 = (39.0727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.21494,'amu*angstrom^2'), symmetry=1, barrier=(4.94189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214919,'amu*angstrom^2'), symmetry=1, barrier=(4.94141,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1399,0.0735239,-0.000135665,1.30253e-07,-4.74458e-11,4792.1,26.2043], Tmin=(100,'K'), Tmax=(843.303,'K')), NASAPolynomial(coeffs=[5.7907,0.0286486,-1.5262e-05,2.98729e-09,-2.06773e-13,4818.97,9.36959], Tmin=(843.303,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.0727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=CCJ=O)"""),
)

species(
    label = 'HCO(1372)',
    structure = SMILES('[CH]=O'),
    E0 = (32.5964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1116.11,1837.47,2716.03],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06159,-0.0018195,9.37895e-06,-8.24742e-09,2.33851e-12,3918.59,3.98389], Tmin=(100,'K'), Tmax=(1105.12,'K')), NASAPolynomial(coeffs=[3.05335,0.00410464,-1.74962e-06,3.28532e-10,-2.28985e-14,4002.52,8.32038], Tmin=(1105.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.5964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O][C]C[C]=O(11839)',
    structure = SMILES('[O][C]C[C]=O'),
    E0 = (483.503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,180,180,2670.18],'cm^-1')),
        HinderedRotor(inertia=(0.388352,'amu*angstrom^2'), symmetry=1, barrier=(8.92898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.388437,'amu*angstrom^2'), symmetry=1, barrier=(8.93093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2072,0.0492927,-0.000101434,1.02484e-07,-3.78192e-11,58207.2,18.9402], Tmin=(100,'K'), Tmax=(872.783,'K')), NASAPolynomial(coeffs=[3.38492,0.0201932,-1.06874e-05,2.05291e-09,-1.39197e-13,58904.3,18.5903], Tmin=(872.783,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CH2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1(C=O)CC1=O(11840)',
    structure = SMILES('[O]C1(C=O)CC1=O'),
    E0 = (-49.6726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78025,0.0443951,-3.81929e-05,1.71885e-08,-3.11048e-12,-5890.65,23.0823], Tmin=(100,'K'), Tmax=(1324.42,'K')), NASAPolynomial(coeffs=[11.1217,0.0161818,-6.23876e-06,1.10361e-09,-7.42119e-14,-8365,-24.6092], Tmin=(1324.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.6726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(cyclopropanone) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[O]C(C=O)C=C=O(11841)',
    structure = SMILES('[O]C(C=O)C=C=O'),
    E0 = (-119.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37368,0.0644429,-0.000105589,9.32465e-08,-3.17858e-11,-14321.3,24.2537], Tmin=(100,'K'), Tmax=(865.923,'K')), NASAPolynomial(coeffs=[7.20175,0.0237032,-1.1081e-05,2.05781e-09,-1.3842e-13,-14812.5,-0.0330156], Tmin=(865.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-OdCsH) + missing(Cdd-CdO2d) + radical(C=OCOJ)"""),
)

species(
    label = '[O]C=CC(=O)C=O(11296)',
    structure = SMILES('[O]C=CC(=O)C=O'),
    E0 = (-249.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37035,0.0505955,-4.35754e-05,1.68869e-08,-2.16249e-12,-29954,23.0866], Tmin=(100,'K'), Tmax=(1144.57,'K')), NASAPolynomial(coeffs=[14.3293,0.0135016,-5.70176e-06,1.0822e-09,-7.66855e-14,-33457.3,-43.5274], Tmin=(1144.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-249.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)H) + radical(C=COJ)"""),
)

species(
    label = '[O][C]1[CH]OC(=O)C1(11842)',
    structure = SMILES('[O][C]1CC([O])=CO1'),
    E0 = (25.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94221,0.0296169,2.80331e-05,-7.21042e-08,3.47037e-11,3148.94,19.9941], Tmin=(100,'K'), Tmax=(895.37,'K')), NASAPolynomial(coeffs=[17.5518,0.00218763,3.11111e-06,-7.77351e-10,5.38905e-14,-1342.13,-63.0579], Tmin=(895.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = 'C=C([O])C=O(2859)',
    structure = SMILES('[CH2]C(=O)C=O'),
    E0 = (-74.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.43166,'amu*angstrom^2'), symmetry=1, barrier=(9.92472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231641,'amu*angstrom^2'), symmetry=1, barrier=(31.0079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.689,0.0301409,-2.32871e-05,9.30138e-09,-1.56938e-12,-8891.4,16.525], Tmin=(100,'K'), Tmax=(1335.25,'K')), NASAPolynomial(coeffs=[7.34752,0.0161853,-7.60952e-06,1.4738e-09,-1.03803e-13,-10135.4,-7.29647], Tmin=(1335.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJCC=O)"""),
)

species(
    label = '[O]C([CH][C]=O)C=O(11843)',
    structure = SMILES('[O]C([CH][C]=O)C=O'),
    E0 = (123.811,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1855,455,950,599.381,4000],'cm^-1')),
        HinderedRotor(inertia=(0.560574,'amu*angstrom^2'), symmetry=1, barrier=(12.8887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131115,'amu*angstrom^2'), symmetry=1, barrier=(3.33972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0132243,'amu*angstrom^2'), symmetry=1, barrier=(3.3931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93226,0.046651,-4.84424e-05,2.70509e-08,-6.17321e-12,14964.6,30.0038], Tmin=(100,'K'), Tmax=(1049.88,'K')), NASAPolynomial(coeffs=[9.29519,0.0185981,-8.36162e-06,1.59945e-09,-1.12558e-13,13418.6,-5.87634], Tmin=(1049.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(CCCJ=O)"""),
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
    E0 = (-77.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (168.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (515.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (756.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (756.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (-68.9776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (162.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-70.1492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-13.8618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (-13.8618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (-52.2887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (-52.2887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (459.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (156.937,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (155.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (192.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (288.768,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (-73.4704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (99.9279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (189.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (188.593,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (131.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (365.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (334.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (374.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (34.6636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (183.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (139.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (210.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (46.8302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (69.5733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (-77.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (247.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (533.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (-68.9776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (-13.8618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (164.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (376.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (168.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (571.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (-49.6726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (-54.4006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (-49.4091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (25.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (364.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (282.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['OCHCO(3676)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(2355)', 'C=C([O])[CH][O](2850)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[O]C=[C]C[C]=O(9880)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH]=C([O])C[C]=O(2677)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C=C1CC(=O)O1(11816)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['O=[C]CC1=COO1(11817)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(240.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 237.9 to 240.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C1=COC(=O)C1(11818)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C=C(O)C=C=O(11819)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['O=[C]CC(O)=C=O(11820)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C(C=C=O)=CO(11821)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C(=C=O)CC=O(11822)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][CH]C([O])[CH][C]=O(11823)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C1O[C]1C[C]=O(11824)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(234.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C1(C[C]=O)[CH]O1(11825)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(232.443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O][C]1CC(=O)C1[O](11826)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(269.702,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 267.4 to 269.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O][CH]C1([O])CC1=O(11827)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(366.03,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CO(2039)', 'C=C([O])[CH][O](2850)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(754.069,'m^3/(mol*s)'), n=1.0822, Ea=(24.7919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [COm;C_pri_rad] for rate rule [COm;C_rad/H2/Cd]
Euclidian distance = 1.0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[O]C=C([O])C=C=O(11380)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-OneDeH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[O]C(=C=O)C[C]=O(11828)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][C]=C[O](9592)', 'ketene(1375)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.560775,'m^3/(mol*s)'), n=2.01066, Ea=(45.2043,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ck;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['OCHCO(3676)', '[CH2][C]=O(1376)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O][C]=C[O](9592)', '[CH2][C]=O(1376)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[O][C]=C([O])C[C]=O(11829)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C=CC([O])=C[O](11298)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]C(=[C]O)C[C]=O(11830)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]C=C(O)[CH][C]=O(11831)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(52326.9,'s^-1'), n=2.1859, Ea=(154.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O][C]=C(O)C[C]=O(11832)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O][C]=C([O])CC=O(11833)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out;XH_out] for rate rule [R4H_DSS;Cd_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]C([CH][C]=O)=CO(11834)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(8.08094e+07,'s^-1'), n=1.1965, Ea=(138.911,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;XH_out] for rate rule [R4H_SDS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['O=[C]CC(=O)C=O(11154)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C([C]=O)C[C]=O(11149)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH][O](1548)', 'O=[C]C[C]=O(2682)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C1C(=O)CC1=O(11165)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]CC(=O)C=C=O(11835)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O][CH][C]1CC(=O)O1(11836)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['O=[C]C[C]1[CH]OO1(11837)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(454.203,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra_O] for rate rule [R4_linear;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O]CC(=O)[CH][C]=O(11838)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(11891.5,'s^-1'), n=2.58622, Ea=(129.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['HCO(1372)', '[O][C]C[C]=O(11839)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C1(C=O)CC1=O(11840)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(27.5894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination
Ea raised from 26.7 to 27.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C(C=O)C=C=O(11841)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O]C=CC(=O)C=O(11296)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['[O][C]1[CH]OC(=O)C1(11842)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(102.705,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;carbonyl_intra_H;radadd_intra_CO] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_CO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 100.4 to 102.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction45',
    reactants = ['[C]=O(2355)', 'C=C([O])C=O(2859)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CO_birad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction46',
    reactants = ['[O]C([CH][C]=O)C=O(11843)'],
    products = ['[O]C=C([O])C[C]=O(11162)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2830',
    isomers = [
        '[O]C=C([O])C[C]=O(11162)',
    ],
    reactants = [
        ('OCHCO(3676)', 'ketene(1375)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2830',
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


species(
    label = 'C[CH]CC[CH][O](1427)',
    structure = SMILES('C[CH]CC[CH][O]'),
    E0 = (272.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,180,1634.92,1634.97,1634.99],'cm^-1')),
        HinderedRotor(inertia=(0.227599,'amu*angstrom^2'), symmetry=1, barrier=(5.23294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227523,'amu*angstrom^2'), symmetry=1, barrier=(5.23121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227601,'amu*angstrom^2'), symmetry=1, barrier=(5.233,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227555,'amu*angstrom^2'), symmetry=1, barrier=(5.23193,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35387,0.0677678,-9.77581e-05,9.50642e-08,-3.64608e-11,32864.7,27.1936], Tmin=(100,'K'), Tmax=(841.762,'K')), NASAPolynomial(coeffs=[1.0997,0.0470954,-2.19303e-05,4.12968e-09,-2.82507e-13,33682.7,32.9805], Tmin=(841.762,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CCsJOH)"""),
)

species(
    label = 'vinoxy(1351)',
    structure = SMILES('[CH2]C=O'),
    E0 = (6.53237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1127.76,1127.76,1127.77,1945],'cm^-1')),
        HinderedRotor(inertia=(0.00342653,'amu*angstrom^2'), symmetry=1, barrier=(3.09254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51396,0.0045106,2.69938e-05,-3.71579e-08,1.42465e-11,808.756,8.86322], Tmin=(100,'K'), Tmax=(957.493,'K')), NASAPolynomial(coeffs=[6.49564,0.00770518,-2.52913e-06,4.69083e-10,-3.5128e-14,-479.661,-9.13857], Tmin=(957.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.53237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""vinoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C3H6(27)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36193e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35107e-06,1.16619e-09,-8.2762e-14,-487.139,-4.54469], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH]CCC=O(1424)',
    structure = SMILES('C[CH]CCC=O'),
    E0 = (-56.4374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.67194,0.0433391,-1.95106e-05,3.62425e-09,-2.37546e-13,-6754.83,21.3081], Tmin=(100,'K'), Tmax=(2727.9,'K')), NASAPolynomial(coeffs=[29.3419,0.00991848,-4.26031e-06,6.61388e-10,-3.60429e-14,-23421.1,-138.001], Tmin=(2727.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.4374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C)C[CH][O](1428)',
    structure = SMILES('[CH2]C(C)C[CH][O]'),
    E0 = (277.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,1187.06,1194.76],'cm^-1')),
        HinderedRotor(inertia=(0.181577,'amu*angstrom^2'), symmetry=1, barrier=(4.17481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183737,'amu*angstrom^2'), symmetry=1, barrier=(4.22447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181979,'amu*angstrom^2'), symmetry=1, barrier=(4.18404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183609,'amu*angstrom^2'), symmetry=1, barrier=(4.22153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3655.69,'J/mol'), sigma=(6.48678,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.01 K, Pc=30.39 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19047,0.0665903,-7.97859e-05,6.4446e-08,-2.19921e-11,33469.6,26.5501], Tmin=(100,'K'), Tmax=(826.775,'K')), NASAPolynomial(coeffs=[5.3147,0.0393165,-1.70223e-05,3.12742e-09,-2.12286e-13,33037.8,8.95072], Tmin=(826.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([O])C[CH]C(1417)',
    structure = SMILES('[CH2]C([O])C[CH]C'),
    E0 = (295.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,250.144,1331.88,1331.91],'cm^-1')),
        HinderedRotor(inertia=(0.00269429,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00269472,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188466,'amu*angstrom^2'), symmetry=1, barrier=(8.36508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188445,'amu*angstrom^2'), symmetry=1, barrier=(8.36471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46999,0.058291,-4.76951e-05,2.35457e-08,-5.16521e-12,35670.4,26.4723], Tmin=(100,'K'), Tmax=(1040.28,'K')), NASAPolynomial(coeffs=[7.28477,0.0359322,-1.5455e-05,2.88417e-09,-1.99774e-13,34460.6,-1.81019], Tmin=(1040.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJC) + radical(CJCO)"""),
)

species(
    label = 'CHCH3(T)(21)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438698,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82365,-0.000909765,3.21389e-05,-3.73492e-08,1.33096e-11,41371.4,7.10941], Tmin=(100,'K'), Tmax=(960.802,'K')), NASAPolynomial(coeffs=[4.3048,0.00943081,-3.27566e-06,5.95138e-10,-4.27321e-14,40709.2,1.84242], Tmin=(960.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C[CH][O](1549)',
    structure = SMILES('[CH2]C[CH][O]'),
    E0 = (338.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1524.15,1525.95],'cm^-1')),
        HinderedRotor(inertia=(0.310045,'amu*angstrom^2'), symmetry=1, barrier=(7.12854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00431249,'amu*angstrom^2'), symmetry=1, barrier=(7.11449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81289,0.029873,-3.56716e-05,3.591e-08,-1.53006e-11,40735.5,17.4078], Tmin=(100,'K'), Tmax=(773.473,'K')), NASAPolynomial(coeffs=[2.08215,0.0257045,-1.21752e-05,2.37389e-09,-1.67404e-13,41086.3,22.2823], Tmin=(773.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(RCCJ) + radical(CCsJOH)"""),
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
    label = '[CH2]C[CH]C(36)',
    structure = SMILES('[CH2]C[CH]C'),
    E0 = (255.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3120.74],'cm^-1')),
        HinderedRotor(inertia=(0.20781,'amu*angstrom^2'), symmetry=1, barrier=(4.77797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.207841,'amu*angstrom^2'), symmetry=1, barrier=(4.77868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104176,'amu*angstrom^2'), symmetry=1, barrier=(71.9779,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98997,0.0287412,-9.51473e-06,4.19256e-10,1.9052e-13,30780.1,16.8971], Tmin=(100,'K'), Tmax=(2154.58,'K')), NASAPolynomial(coeffs=[12.4234,0.0182237,-7.06297e-06,1.16765e-09,-7.1179e-14,25091.2,-39.6233], Tmin=(2154.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CH3(17)',
    structure = SMILES('[CH3]'),
    E0 = (136.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([532.913,1391.12,1391.12,2779.21,3448.45,3448.46],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (15.0345,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.948,0.0008276,8.34932e-06,-9.82634e-09,3.80104e-12,16425.4,0.336655], Tmin=(100,'K'), Tmax=(660.467,'K')), NASAPolynomial(coeffs=[3.2217,0.00522646,-1.64125e-06,2.58225e-10,-1.62579e-14,16521.3,3.53938], Tmin=(660.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""CH3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]CC[CH][O](1551)',
    structure = SMILES('[CH]CC[CH][O]'),
    E0 = (550.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,267.058,267.102,267.117,267.119,2495.83,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0270636,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.633872,'amu*angstrom^2'), symmetry=1, barrier=(32.087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00725894,'amu*angstrom^2'), symmetry=1, barrier=(32.0875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76252,0.0556378,-8.30439e-05,7.63546e-08,-2.78809e-11,66234.5,21.7518], Tmin=(100,'K'), Tmax=(834.654,'K')), NASAPolynomial(coeffs=[4.24691,0.0306058,-1.44688e-05,2.73997e-09,-1.88e-13,66277,12.9539], Tmin=(834.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOH) + radical(CCJ2_triplet)"""),
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
    label = 'C[C]CC[CH][O](1607)',
    structure = SMILES('C[C]CC[CH][O]'),
    E0 = (526.306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09038,0.0707393,-9.81414e-05,8.63858e-08,-3.10856e-11,63398.3,25.3106], Tmin=(100,'K'), Tmax=(813.376,'K')), NASAPolynomial(coeffs=[5.20804,0.0389531,-1.82468e-05,3.46397e-09,-2.38961e-13,63110.1,8.64175], Tmin=(813.376,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCsJOH) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH]CC[C][O](1608)',
    structure = SMILES('C[CH]CC[C][O]'),
    E0 = (553.228,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39791,0.0651243,-9.22518e-05,8.68014e-08,-3.28923e-11,66624,26.2182], Tmin=(100,'K'), Tmax=(822.298,'K')), NASAPolynomial(coeffs=[2.80234,0.0418568,-1.98268e-05,3.77682e-09,-2.60695e-13,66948.7,23.0962], Tmin=(822.298,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CH2_triplet)"""),
)

species(
    label = 'CC1CCC1[O](1609)',
    structure = SMILES('CC1CCC1[O]'),
    E0 = (38.9416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01798,0.0241676,6.78132e-05,-1.03324e-07,4.05756e-11,4772.47,19.7233], Tmin=(100,'K'), Tmax=(965.457,'K')), NASAPolynomial(coeffs=[13.6643,0.024601,-8.5014e-06,1.60472e-09,-1.20214e-13,254.659,-47.805], Tmin=(965.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(38.9416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC=CCC[O](1610)',
    structure = SMILES('CC=CCC[O]'),
    E0 = (13.3303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31296,0.0451672,-2.11346e-05,4.09792e-09,-2.77547e-13,1655.89,21.4808], Tmin=(100,'K'), Tmax=(2471.4,'K')), NASAPolynomial(coeffs=[23.6483,0.0160869,-6.7931e-06,1.12179e-09,-6.67736e-14,-10554.5,-104.122], Tmin=(2471.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(13.3303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ)"""),
)

species(
    label = 'CCCC=C[O](1611)',
    structure = SMILES('CCCC=C[O]'),
    E0 = (-107.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29527,0.0448541,1.22627e-05,-5.07265e-08,2.34352e-11,-12834.6,23.2996], Tmin=(100,'K'), Tmax=(962.097,'K')), NASAPolynomial(coeffs=[15.4358,0.021716,-7.24852e-06,1.31039e-09,-9.51746e-14,-17205.5,-52.9492], Tmin=(962.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-107.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'C=CCCC[O](1431)',
    structure = SMILES('C=CCCC[O]'),
    E0 = (25.5756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07901,0.0470706,-2.35534e-05,5.22033e-09,-4.44166e-13,3139.76,22.7635], Tmin=(100,'K'), Tmax=(2690.68,'K')), NASAPolynomial(coeffs=[21.2435,0.0185803,-7.67055e-06,1.28505e-09,-7.85235e-14,-7173.27,-88.6627], Tmin=(2690.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.5756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ)"""),
)

species(
    label = 'CC=CC[CH][O](1612)',
    structure = SMILES('CC=CC[CH][O]'),
    E0 = (193.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,253.595,254.2,1965.9],'cm^-1')),
        HinderedRotor(inertia=(0.141536,'amu*angstrom^2'), symmetry=1, barrier=(7.80578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169182,'amu*angstrom^2'), symmetry=1, barrier=(7.83736,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183541,'amu*angstrom^2'), symmetry=1, barrier=(7.82155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53383,0.0596908,-7.22149e-05,6.2919e-08,-2.36218e-11,23371.6,23.8497], Tmin=(100,'K'), Tmax=(781.527,'K')), NASAPolynomial(coeffs=[3.69953,0.0402207,-1.87509e-05,3.58326e-09,-2.49322e-13,23289.2,15.5737], Tmin=(781.527,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH]CC=C[O](1613)',
    structure = SMILES('C[CH]CC=C[O]'),
    E0 = (86.8142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,301.661,904.776,909.759],'cm^-1')),
        HinderedRotor(inertia=(0.568771,'amu*angstrom^2'), symmetry=1, barrier=(13.0772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0225057,'amu*angstrom^2'), symmetry=1, barrier=(13.0941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569629,'amu*angstrom^2'), symmetry=1, barrier=(13.097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52126,0.0447421,-6.09408e-06,-2.3494e-08,1.22034e-11,10539.3,25.224], Tmin=(100,'K'), Tmax=(982.154,'K')), NASAPolynomial(coeffs=[12.1696,0.0238802,-8.60428e-06,1.54067e-09,-1.07718e-13,7362.14,-31.4825], Tmin=(982.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.8142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CCC=O(1483)',
    structure = SMILES('[CH2][CH]CCC=O'),
    E0 = (148.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1550.9,1552.3],'cm^-1')),
        HinderedRotor(inertia=(0.144291,'amu*angstrom^2'), symmetry=1, barrier=(3.31754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144856,'amu*angstrom^2'), symmetry=1, barrier=(3.33053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144612,'amu*angstrom^2'), symmetry=1, barrier=(3.32491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144802,'amu*angstrom^2'), symmetry=1, barrier=(3.32928,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.13066,0.0398929,-1.8153e-05,3.39906e-09,-2.30727e-13,17906.2,21.3822], Tmin=(100,'K'), Tmax=(3053.92,'K')), NASAPolynomial(coeffs=[39.1437,-0.00225523,8.25472e-08,-4.33358e-11,6.99894e-15,-6431.62,-196.399], Tmin=(3053.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][O](1556)',
    structure = SMILES('[CH2][CH][O]'),
    E0 = (367.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1667.93],'cm^-1')),
        HinderedRotor(inertia=(0.00517725,'amu*angstrom^2'), symmetry=1, barrier=(10.2207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22654,0.0212781,-3.59354e-05,3.91027e-08,-1.59281e-11,44278.5,12.2199], Tmin=(100,'K'), Tmax=(830.699,'K')), NASAPolynomial(coeffs=[2.17156,0.016018,-7.76586e-06,1.51127e-09,-1.05387e-13,44810.5,19.2612], Tmin=(830.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO) + radical(CCsJOH)"""),
)

species(
    label = 'C3H6(T)(28)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238388,'amu*angstrom^2'), symmetry=1, barrier=(5.48101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0090964,'amu*angstrom^2'), symmetry=1, barrier=(22.1004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93779,0.019099,4.26879e-06,-1.44878e-08,5.7496e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.8,'K')), NASAPolynomial(coeffs=[5.93905,0.0171893,-6.69156e-06,1.21547e-09,-8.39803e-14,33151.2,-4.14862], Tmin=(1046.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH][CH]C[CH][O](1614)',
    structure = SMILES('C[CH][CH]C[CH][O]'),
    E0 = (466.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1832.14,1832.57,1833.1,1833.72],'cm^-1')),
        HinderedRotor(inertia=(0.184148,'amu*angstrom^2'), symmetry=1, barrier=(4.23393,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185422,'amu*angstrom^2'), symmetry=1, barrier=(4.26321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184368,'amu*angstrom^2'), symmetry=1, barrier=(4.23898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184685,'amu*angstrom^2'), symmetry=1, barrier=(4.24628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6055,0.0682665,-0.000121744,1.31137e-07,-5.16432e-11,56237.9,28.8539], Tmin=(100,'K'), Tmax=(868.954,'K')), NASAPolynomial(coeffs=[-3.15893,0.0510339,-2.4391e-05,4.5789e-09,-3.09864e-13,58544.5,59.6783], Tmin=(868.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH]C[CH][CH][O](1615)',
    structure = SMILES('C[CH]C[CH][CH][O]'),
    E0 = (472.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,180,1751.95,1752.31,1752.69],'cm^-1')),
        HinderedRotor(inertia=(0.107734,'amu*angstrom^2'), symmetry=1, barrier=(2.47701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107747,'amu*angstrom^2'), symmetry=1, barrier=(2.47732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108273,'amu*angstrom^2'), symmetry=1, barrier=(2.48942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107464,'amu*angstrom^2'), symmetry=1, barrier=(2.47082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53489,0.0627085,-8.96356e-05,8.61415e-08,-3.28018e-11,56901.9,29.0405], Tmin=(100,'K'), Tmax=(840.959,'K')), NASAPolynomial(coeffs=[1.71363,0.0427141,-1.98249e-05,3.72938e-09,-2.55054e-13,57548.8,32.2341], Tmin=(840.959,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(472.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CCJCO) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]CC[CH][O](1487)',
    structure = SMILES('[CH2][CH]CC[CH][O]'),
    E0 = (477.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,1700.57,1701.02,1701.76],'cm^-1')),
        HinderedRotor(inertia=(0.13324,'amu*angstrom^2'), symmetry=1, barrier=(3.06346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13226,'amu*angstrom^2'), symmetry=1, barrier=(3.04093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132646,'amu*angstrom^2'), symmetry=1, barrier=(3.04979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132376,'amu*angstrom^2'), symmetry=1, barrier=(3.0436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37842,0.0686508,-0.000108717,1.07913e-07,-4.11036e-11,57547.9,28.8918], Tmin=(100,'K'), Tmax=(853.943,'K')), NASAPolynomial(coeffs=[0.965186,0.0446752,-2.10876e-05,3.96924e-09,-2.704e-13,58563.2,36.3519], Tmin=(853.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = 'CC[CH]C[CH][O](1589)',
    structure = SMILES('CC[CH]C[CH][O]'),
    E0 = (272.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,180,1766.9,1769.77,1771.53],'cm^-1')),
        HinderedRotor(inertia=(0.216375,'amu*angstrom^2'), symmetry=1, barrier=(4.97489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216127,'amu*angstrom^2'), symmetry=1, barrier=(4.9692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.217102,'amu*angstrom^2'), symmetry=1, barrier=(4.99161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.215685,'amu*angstrom^2'), symmetry=1, barrier=(4.95903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39325,0.068218,-0.000102858,1.03326e-07,-4.0249e-11,32863.5,26.8804], Tmin=(100,'K'), Tmax=(844.753,'K')), NASAPolynomial(coeffs=[0.0134669,0.0490355,-2.31335e-05,4.3723e-09,-2.99314e-13,34014.1,38.7351], Tmin=(844.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH]C[CH]C[O](1616)',
    structure = SMILES('C[CH]C[CH]C[O]'),
    E0 = (292.142,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,180,1612.9,1614.45,1614.57],'cm^-1')),
        HinderedRotor(inertia=(0.204571,'amu*angstrom^2'), symmetry=1, barrier=(4.70348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205189,'amu*angstrom^2'), symmetry=1, barrier=(4.7177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205215,'amu*angstrom^2'), symmetry=1, barrier=(4.71829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00254928,'amu*angstrom^2'), symmetry=1, barrier=(4.71154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.815,0.0538447,-5.80713e-05,5.36507e-08,-2.17343e-11,35209.6,28.4864], Tmin=(100,'K'), Tmax=(796.49,'K')), NASAPolynomial(coeffs=[0.995246,0.0460902,-2.11107e-05,4.00171e-09,-2.77137e-13,35716.7,34.6186], Tmin=(796.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO) + radical(RCCJC)"""),
)

species(
    label = '[CH2]CCC[CH][O](1617)',
    structure = SMILES('[CH2]CCC[CH][O]'),
    E0 = (283.338,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,180,1501.22,1502.37],'cm^-1')),
        HinderedRotor(inertia=(0.159852,'amu*angstrom^2'), symmetry=1, barrier=(3.6753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158992,'amu*angstrom^2'), symmetry=1, barrier=(3.65553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159028,'amu*angstrom^2'), symmetry=1, barrier=(3.65637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158864,'amu*angstrom^2'), symmetry=1, barrier=(3.65259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19726,0.0682176,-8.84199e-05,7.8232e-08,-2.89254e-11,34172.2,26.8081], Tmin=(100,'K'), Tmax=(802.525,'K')), NASAPolynomial(coeffs=[4.00257,0.0429161,-1.99724e-05,3.79698e-09,-2.62747e-13,34086.4,16.1621], Tmin=(802.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.338,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = 'CCC[CH][CH][O](1618)',
    structure = SMILES('CCC[CH][CH][O]'),
    E0 = (277.993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,342.943,342.944,1713.23,1713.24],'cm^-1')),
        HinderedRotor(inertia=(0.00143345,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.080103,'amu*angstrom^2'), symmetry=1, barrier=(6.68525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0801029,'amu*angstrom^2'), symmetry=1, barrier=(6.68521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0800997,'amu*angstrom^2'), symmetry=1, barrier=(6.68515,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25784,0.0474353,1.22184e-05,-1.25321e-07,1.21215e-10,33488,23.9309], Tmin=(100,'K'), Tmax=(431.675,'K')), NASAPolynomial(coeffs=[4.53947,0.0413495,-1.89525e-05,3.61724e-09,-2.52566e-13,33150.7,13.2152], Tmin=(431.675,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH][CH]CC[O](1619)',
    structure = SMILES('C[CH][CH]CC[O]'),
    E0 = (286.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1903.83,1904.96,1904.97,1905.5],'cm^-1')),
        HinderedRotor(inertia=(0.17247,'amu*angstrom^2'), symmetry=1, barrier=(3.96543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172705,'amu*angstrom^2'), symmetry=1, barrier=(3.97083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17266,'amu*angstrom^2'), symmetry=1, barrier=(3.96979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172112,'amu*angstrom^2'), symmetry=1, barrier=(3.95718,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85492,0.0597785,-9.15327e-05,1.00383e-07,-4.12641e-11,34546.8,28.4088], Tmin=(100,'K'), Tmax=(851.798,'K')), NASAPolynomial(coeffs=[-3.7299,0.0541488,-2.55216e-05,4.81376e-09,-3.28789e-13,36653.9,61.2401], Tmin=(851.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C[CH][CH]O(1620)',
    structure = SMILES('C[CH]C[CH][CH]O'),
    E0 = (246.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23988,0.0645945,-7.39894e-05,5.61691e-08,-1.83474e-11,29771,29.2434], Tmin=(100,'K'), Tmax=(795.041,'K')), NASAPolynomial(coeffs=[6.24778,0.0365777,-1.58075e-05,2.91855e-09,-1.99344e-13,29063.9,6.79249], Tmin=(795.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCJCO) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH][CH]C[CH]O(1621)',
    structure = SMILES('C[CH][CH]C[CH]O'),
    E0 = (241.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26177,0.0707789,-0.000108543,1.04723e-07,-3.889e-11,29109,29.2282], Tmin=(100,'K'), Tmax=(867.965,'K')), NASAPolynomial(coeffs=[1.50951,0.0446569,-2.02297e-05,3.73313e-09,-2.51197e-13,30007,33.4885], Tmin=(867.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]CCC[O](1622)',
    structure = SMILES('[CH2][CH]CCC[O]'),
    E0 = (297.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,1513.31,1513.37,1513.38],'cm^-1')),
        HinderedRotor(inertia=(0.203947,'amu*angstrom^2'), symmetry=1, barrier=(4.68915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203961,'amu*angstrom^2'), symmetry=1, barrier=(4.68947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00288566,'amu*angstrom^2'), symmetry=1, barrier=(4.68973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203956,'amu*angstrom^2'), symmetry=1, barrier=(4.68934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.46601,0.0404004,-1.59913e-05,1.98405e-09,-4.38285e-15,35768.9,21.6905], Tmin=(100,'K'), Tmax=(2682.42,'K')), NASAPolynomial(coeffs=[36.4398,0.00376259,-2.51164e-06,3.75637e-10,-1.68064e-14,13570.2,-178.329], Tmin=(2682.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[CH]O(1623)',
    structure = SMILES('[CH2][CH]CC[CH]O'),
    E0 = (252.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04882,0.0709834,-9.48247e-05,8.05107e-08,-2.78854e-11,30418.5,29.2162], Tmin=(100,'K'), Tmax=(839.54,'K')), NASAPolynomial(coeffs=[5.59304,0.0383703,-1.69692e-05,3.13385e-09,-2.12609e-13,30041.8,10.3886], Tmin=(839.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
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
    E0 = (272.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (272.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (437.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (453.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (737.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (742.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (742.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (738.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (765.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (280.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (335.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (335.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (297.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (411.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (306.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (363.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (386.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (327.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (652.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (678.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (684.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (689.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (431.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (429.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (431.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (401.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (402.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (428.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (395.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (371.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (366.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]CC[CH][O](1427)'],
    products = ['vinoxy(1351)', 'C3H6(27)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[CH]CC[CH][O](1427)'],
    products = ['C[CH]CCC=O(1424)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C)C[CH][O](1428)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([O])C[CH]C(1417)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHCH3(T)(21)', '[CH2]C[CH][O](1549)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][O](1548)', '[CH2]C[CH]C(36)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', '[CH]CC[CH][O](1551)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C[C]CC[CH][O](1607)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C[CH]CC[C][O](1608)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]CC[CH][O](1427)'],
    products = ['CC1CCC1[O](1609)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]CC[CH][O](1427)'],
    products = ['CC=CCC[O](1610)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH]CC[CH][O](1427)'],
    products = ['CCCC=C[O](1611)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH]CC[CH][O](1427)'],
    products = ['C=CCCC[O](1431)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'CC=CC[CH][O](1612)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'C[CH]CC=C[O](1613)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2][CH]CCC=O(1483)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH][O](1556)', 'C3H6(27)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00620445,'m^3/(mol*s)'), n=2.46568, Ea=(12.4666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['vinoxy(1351)', 'C3H6(T)(28)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH][O](1556)', 'C3H6(T)(28)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', 'C[CH][CH]C[CH][O](1614)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', 'C[CH]C[CH][CH][O](1615)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][CH]CC[CH][O](1487)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC[CH]C[CH][O](1589)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C[CH]C[CH]C[O](1616)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]CCC[CH][O](1617)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CCC[CH][CH][O](1618)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.3913e+06,'s^-1'), n=1.66106, Ea=(123.033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C[CH][CH]CC[O](1619)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[CH]CC[CH][O](1427)'],
    products = ['C[CH]C[CH][CH]O(1620)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.30814e+09,'s^-1'), n=1.19923, Ea=(155.469,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[CH]CC[CH][O](1427)'],
    products = ['C[CH][CH]C[CH]O(1621)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.04154e+06,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]CCC[O](1622)'],
    products = ['C[CH]CC[CH][O](1427)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[CH]CC[CH][O](1427)'],
    products = ['[CH2][CH]CC[CH]O(1623)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.68e+09,'s^-1'), n=0, Ea=(93.5124,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R6Hall;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #473',
    isomers = [
        'C[CH]CC[CH][O](1427)',
    ],
    reactants = [
        ('vinoxy(1351)', 'C3H6(27)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #473',
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


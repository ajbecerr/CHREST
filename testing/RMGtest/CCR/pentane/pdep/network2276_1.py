species(
    label = 'C[CH]C[CH][C]=O(2434)',
    structure = SMILES('C[CH]C[CH][C]=O'),
    E0 = (271.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,2191.4,2191.84],'cm^-1')),
        HinderedRotor(inertia=(0.00214683,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109702,'amu*angstrom^2'), symmetry=1, barrier=(6.30235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21024,'amu*angstrom^2'), symmetry=1, barrier=(67.5894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113548,'amu*angstrom^2'), symmetry=1, barrier=(6.32349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23108,0.0440156,-2.01448e-05,-2.09479e-08,2.52039e-11,32658.9,23.6185], Tmin=(100,'K'), Tmax=(539.722,'K')), NASAPolynomial(coeffs=[4.99773,0.0326275,-1.38305e-05,2.54716e-09,-1.74683e-13,32227.4,10.7471], Tmin=(539.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'HCCO(2227)',
    structure = SMILES('[CH]=C=O'),
    E0 = (165.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,231.114,1089.61,3388.99],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47174,0.0113766,-1.1308e-05,6.13951e-09,-1.35087e-12,19887.1,6.87336], Tmin=(100,'K'), Tmax=(1096.2,'K')), NASAPolynomial(coeffs=[5.39217,0.00436893,-1.71893e-06,3.07751e-10,-2.08706e-14,19466,-2.56797], Tmin=(1096.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][CH][C]=O(9632)',
    structure = SMILES('[CH2][CH][C]=O'),
    E0 = (334.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.0670133,'amu*angstrom^2'), symmetry=1, barrier=(7.51319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00540518,'amu*angstrom^2'), symmetry=1, barrier=(7.51488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.78792,0.0268854,-2.93724e-05,1.83573e-08,-4.61895e-12,40303.1,15.8795], Tmin=(100,'K'), Tmax=(968.231,'K')), NASAPolynomial(coeffs=[6.80976,0.01027,-3.63122e-06,6.33225e-10,-4.24838e-14,39524.3,-3.39371], Tmin=(968.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJCHO) + radical(CJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=[C][O](6861)',
    structure = SMILES('[CH][C]=O'),
    E0 = (425.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1491.4,1492.35],'cm^-1')),
        HinderedRotor(inertia=(0.0736325,'amu*angstrom^2'), symmetry=1, barrier=(5.08901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62981,0.00805288,-4.35162e-06,1.03741e-09,-8.43183e-14,51134.2,10.4631], Tmin=(100,'K'), Tmax=(1941.7,'K')), NASAPolynomial(coeffs=[6.62042,0.00292938,-1.19496e-06,2.28726e-10,-1.56225e-14,49777.3,-6.4529], Tmin=(1941.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet) + radical(CsCJ=O)"""),
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
    label = '[CH]C[CH][C]=O(9634)',
    structure = SMILES('[CH]C[CH][C]=O'),
    E0 = (548.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1855,455,950,180,1003.55,1004.27,1004.31],'cm^-1')),
        HinderedRotor(inertia=(0.0732029,'amu*angstrom^2'), symmetry=1, barrier=(1.68308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0736242,'amu*angstrom^2'), symmetry=1, barrier=(1.69277,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0735852,'amu*angstrom^2'), symmetry=1, barrier=(1.69187,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17502,0.0387035,-3.76272e-05,2.01297e-08,-4.36699e-12,66048.5,19.7699], Tmin=(100,'K'), Tmax=(1112.39,'K')), NASAPolynomial(coeffs=[8.76008,0.0150243,-5.6968e-06,9.93278e-10,-6.62104e-14,64583.5,-12.7006], Tmin=(1112.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(548.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = 'C[C]C[CH][C]=O(9682)',
    structure = SMILES('C[C]C[CH][C]=O'),
    E0 = (524.821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,423.136,423.176,423.377],'cm^-1')),
        HinderedRotor(inertia=(0.0825782,'amu*angstrom^2'), symmetry=1, barrier=(10.5006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00603184,'amu*angstrom^2'), symmetry=1, barrier=(10.501,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828279,'amu*angstrom^2'), symmetry=1, barrier=(10.5013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.560572,'amu*angstrom^2'), symmetry=1, barrier=(71.293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56358,0.0530671,-4.99918e-05,2.62761e-08,-5.69795e-12,63209.6,23.112], Tmin=(100,'K'), Tmax=(1100.41,'K')), NASAPolynomial(coeffs=[9.7471,0.0233194,-9.44117e-06,1.70864e-09,-1.16412e-13,61408.6,-17.1516], Tmin=(1100.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C[C][C]=O(9683)',
    structure = SMILES('C[CH]C[C][C]=O'),
    E0 = (551.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,180,180,1259.64],'cm^-1')),
        HinderedRotor(inertia=(0.174771,'amu*angstrom^2'), symmetry=1, barrier=(4.01832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174661,'amu*angstrom^2'), symmetry=1, barrier=(4.01579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174744,'amu*angstrom^2'), symmetry=1, barrier=(4.01772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174703,'amu*angstrom^2'), symmetry=1, barrier=(4.01676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57772,0.059909,-8.80875e-05,7.94099e-08,-2.83201e-11,66439.6,26.0342], Tmin=(100,'K'), Tmax=(852.578,'K')), NASAPolynomial(coeffs=[4.50643,0.032167,-1.46452e-05,2.72007e-09,-1.84299e-13,66449.1,15.3563], Tmin=(852.578,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'CC1CC1[C]=O(9668)',
    structure = SMILES('CC1CC1[C]=O'),
    E0 = (55.5993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10499,0.0277926,3.45362e-05,-6.39914e-08,2.62384e-11,6767.77,20.2815], Tmin=(100,'K'), Tmax=(967.749,'K')), NASAPolynomial(coeffs=[12.2118,0.0203299,-7.07977e-06,1.31421e-09,-9.65707e-14,3204.88,-36.4482], Tmin=(967.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.5993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'CCC=C[C]=O(9684)',
    structure = SMILES('CCC=C[C]=O'),
    E0 = (20.6116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83641,0.0525271,-4.2117e-05,1.82991e-08,-3.48325e-12,2552.66,20.4505], Tmin=(100,'K'), Tmax=(1168.55,'K')), NASAPolynomial(coeffs=[7.92868,0.0316729,-1.53478e-05,3.02698e-09,-2.15933e-13,1128.84,-9.8901], Tmin=(1168.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.6116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC=CC[C]=O(2431)',
    structure = SMILES('CC=CC[C]=O'),
    E0 = (30.7963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,387.023],'cm^-1')),
        HinderedRotor(inertia=(0.112079,'amu*angstrom^2'), symmetry=1, barrier=(11.8552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111689,'amu*angstrom^2'), symmetry=1, barrier=(11.8567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111223,'amu*angstrom^2'), symmetry=1, barrier=(11.856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83002,0.0398038,-9.11535e-06,-1.05206e-08,4.83654e-12,3788.49,22.0947], Tmin=(100,'K'), Tmax=(1207.67,'K')), NASAPolynomial(coeffs=[11.1225,0.0255305,-1.18873e-05,2.32622e-09,-1.6555e-13,340.438,-29.4728], Tmin=(1207.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.7963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CCC[C]=O(2291)',
    structure = SMILES('C=CCC[C]=O'),
    E0 = (36.8592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,302.475,302.475],'cm^-1')),
        HinderedRotor(inertia=(0.110264,'amu*angstrom^2'), symmetry=1, barrier=(7.15875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110264,'amu*angstrom^2'), symmetry=1, barrier=(7.15876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110265,'amu*angstrom^2'), symmetry=1, barrier=(7.15876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88306,0.0471705,-3.34342e-05,1.2911e-08,-2.1467e-12,4508.58,22.7686], Tmin=(100,'K'), Tmax=(1342.11,'K')), NASAPolynomial(coeffs=[8.35823,0.0278721,-1.18655e-05,2.1972e-09,-1.50995e-13,2770.5,-10.3756], Tmin=(1342.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.8592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=C[CH][C]=O(9685)',
    structure = SMILES('C[CH]C=C[C]=O'),
    E0 = (161.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,203.658],'cm^-1')),
        HinderedRotor(inertia=(0.567983,'amu*angstrom^2'), symmetry=1, barrier=(16.7708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.5604,'amu*angstrom^2'), symmetry=1, barrier=(16.7664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564368,'amu*angstrom^2'), symmetry=1, barrier=(16.7593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99911,0.0478573,-3.55847e-05,1.33834e-08,-2.13733e-12,19519.3,18.7486], Tmin=(100,'K'), Tmax=(1379.98,'K')), NASAPolynomial(coeffs=[9.08184,0.0273273,-1.32692e-05,2.6028e-09,-1.84291e-13,17564.5,-17.7026], Tmin=(1379.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_S) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2][CH]CC=C=O(2614)',
    structure = SMILES('[CH2][CH]CC=C=O'),
    E0 = (274.229,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2120,512.5,787.5,1362.85,1365.22],'cm^-1')),
        HinderedRotor(inertia=(0.191838,'amu*angstrom^2'), symmetry=1, barrier=(4.41074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191184,'amu*angstrom^2'), symmetry=1, barrier=(4.3957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191816,'amu*angstrom^2'), symmetry=1, barrier=(4.41023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02067,0.0470356,-5.03443e-05,3.98626e-08,-1.42782e-11,33050.1,24.6875], Tmin=(100,'K'), Tmax=(760.726,'K')), NASAPolynomial(coeffs=[4.20178,0.0322327,-1.45815e-05,2.75999e-09,-1.91573e-13,32814.7,15.3956], Tmin=(760.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.229,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CH2CHCO(3668)',
    structure = SMILES('C=C[C]=O'),
    E0 = (83.3963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.27992,'amu*angstrom^2'), symmetry=1, barrier=(29.4278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08334,0.0153205,6.54259e-06,-1.77535e-08,7.39369e-12,10067.5,11.8951], Tmin=(100,'K'), Tmax=(1002.99,'K')), NASAPolynomial(coeffs=[6.9784,0.0111885,-4.32948e-06,8.06749e-10,-5.75264e-14,8712.67,-9.76694], Tmin=(1002.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.3963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH][CH][CH][C]=O(9686)',
    structure = SMILES('C[CH][CH]C=[C][O]'),
    E0 = (467.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,1685,370,308.955,1572.09,3520.03],'cm^-1')),
        HinderedRotor(inertia=(0.253637,'amu*angstrom^2'), symmetry=1, barrier=(93.8127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181237,'amu*angstrom^2'), symmetry=1, barrier=(0.120391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.6314,'amu*angstrom^2'), symmetry=1, barrier=(42.9139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82747,0.0422666,-2.64809e-05,8.50649e-09,-1.12213e-12,56330.6,26.1604], Tmin=(100,'K'), Tmax=(1724.78,'K')), NASAPolynomial(coeffs=[10.9601,0.0210867,-8.06124e-06,1.38686e-09,-9.01712e-14,53180.2,-22.8776], Tmin=(1724.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]C[CH][C]=O(2616)',
    structure = SMILES('[CH2][CH]C[CH][C]=O'),
    E0 = (476.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,1652.74,1652.79],'cm^-1')),
        HinderedRotor(inertia=(0.1661,'amu*angstrom^2'), symmetry=1, barrier=(8.24258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0333532,'amu*angstrom^2'), symmetry=1, barrier=(64.6535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00425184,'amu*angstrom^2'), symmetry=1, barrier=(8.24363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16632,'amu*angstrom^2'), symmetry=1, barrier=(8.24369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8346,0.0511755,-6.13261e-05,4.91033e-08,-1.65214e-11,57360.1,26.7553], Tmin=(100,'K'), Tmax=(839.036,'K')), NASAPolynomial(coeffs=[5.14519,0.0296716,-1.2654e-05,2.30336e-09,-1.55372e-13,57005.9,12.5647], Tmin=(839.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(476.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[CH][CH][C]=O(9665)',
    structure = SMILES('CC[CH]C=[C][O]'),
    E0 = (273.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,180,180,912.471],'cm^-1')),
        HinderedRotor(inertia=(0.821508,'amu*angstrom^2'), symmetry=1, barrier=(18.8881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00550086,'amu*angstrom^2'), symmetry=1, barrier=(3.25338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0318527,'amu*angstrom^2'), symmetry=1, barrier=(19.0465,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69094,0.0415104,-5.97338e-06,-2.0245e-08,1.02188e-11,32952.6,23.9027], Tmin=(100,'K'), Tmax=(1012.71,'K')), NASAPolynomial(coeffs=[11.5756,0.0229622,-8.85563e-06,1.63504e-09,-1.15547e-13,29899.6,-29.0987], Tmin=(1012.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = 'C[CH][CH]C[C]=O(2433)',
    structure = SMILES('C[CH][CH]C[C]=O'),
    E0 = (303.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,2057.55,2058.05],'cm^-1')),
        HinderedRotor(inertia=(0.00216421,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106956,'amu*angstrom^2'), symmetry=1, barrier=(5.89274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106314,'amu*angstrom^2'), symmetry=1, barrier=(5.90618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108095,'amu*angstrom^2'), symmetry=1, barrier=(5.90437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8548,0.051267,-5.85796e-05,4.84354e-08,-1.73425e-11,36566.9,27.6321], Tmin=(100,'K'), Tmax=(806.443,'K')), NASAPolynomial(coeffs=[4.01809,0.0342252,-1.51418e-05,2.8214e-09,-1.93431e-13,36423.3,18.9333], Tmin=(806.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJC) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]CC[CH][C]=O(2618)',
    structure = SMILES('[CH2]CC[CH][C]=O'),
    E0 = (281.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,776.85,2892.87],'cm^-1')),
        HinderedRotor(inertia=(0.00816995,'amu*angstrom^2'), symmetry=1, barrier=(3.5861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.571391,'amu*angstrom^2'), symmetry=1, barrier=(13.1374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153428,'amu*angstrom^2'), symmetry=1, barrier=(3.52762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.44186,'amu*angstrom^2'), symmetry=1, barrier=(79.1352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6661,0.0506142,-4.05288e-05,1.83483e-08,-3.52333e-12,33983.6,24.6233], Tmin=(100,'K'), Tmax=(1208.74,'K')), NASAPolynomial(coeffs=[8.90296,0.0266656,-1.08096e-05,1.95699e-09,-1.33155e-13,32234.1,-11.6623], Tmin=(1208.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][CH][CH]C=O(9687)',
    structure = SMILES('C[CH][CH]C=C[O]'),
    E0 = (227.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72947,0.0393994,3.26722e-06,-3.24953e-08,1.53866e-11,27504.5,23.3701], Tmin=(100,'K'), Tmax=(973.332,'K')), NASAPolynomial(coeffs=[12.192,0.0214294,-7.60793e-06,1.37046e-09,-9.70404e-14,24282.3,-32.9126], Tmin=(973.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(Allyl_S)"""),
)

species(
    label = '[CH2][CH]CC[C]=O(2285)',
    structure = SMILES('[CH2][CH]CC[C]=O'),
    E0 = (308.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180,1901.01],'cm^-1')),
        HinderedRotor(inertia=(0.197881,'amu*angstrom^2'), symmetry=1, barrier=(4.54968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198896,'amu*angstrom^2'), symmetry=1, barrier=(4.57302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198328,'amu*angstrom^2'), symmetry=1, barrier=(4.55996,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00176133,'amu*angstrom^2'), symmetry=1, barrier=(4.53194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3428.95,'J/mol'), sigma=(5.98513,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.59 K, Pc=36.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68036,0.0574379,-7.85385e-05,7.14535e-08,-2.62199e-11,37213.7,27.5467], Tmin=(100,'K'), Tmax=(841.119,'K')), NASAPolynomial(coeffs=[3.32794,0.0360821,-1.63422e-05,3.04616e-09,-2.07501e-13,37414.8,22.7263], Tmin=(841.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]CC=C[O](1849)',
    structure = SMILES('[CH2][CH]CC=C[O]'),
    E0 = (292.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,402.327,406.001,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00103691,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.198951,'amu*angstrom^2'), symmetry=1, barrier=(23.1322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416619,'amu*angstrom^2'), symmetry=1, barrier=(49.087,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55837,0.0454769,-1.65489e-05,-1.12387e-08,7.76625e-12,35222,26.8773], Tmin=(100,'K'), Tmax=(990.008,'K')), NASAPolynomial(coeffs=[11.9568,0.0215968,-7.842e-06,1.39948e-09,-9.72242e-14,32274.4,-27.6729], Tmin=(990.008,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CC=C=O(2432)',
    structure = SMILES('C[CH]CC=C=O'),
    E0 = (68.9827,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,246.668,1625.9],'cm^-1')),
        HinderedRotor(inertia=(0.115588,'amu*angstrom^2'), symmetry=1, barrier=(5.00744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11461,'amu*angstrom^2'), symmetry=1, barrier=(5.00193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114333,'amu*angstrom^2'), symmetry=1, barrier=(4.96075,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14409,0.0442309,-3.17307e-05,1.54891e-08,-3.86259e-12,8360.58,22.4694], Tmin=(100,'K'), Tmax=(855.785,'K')), NASAPolynomial(coeffs=[4.10428,0.0350693,-1.56735e-05,2.98095e-09,-2.088e-13,8025.06,13.3177], Tmin=(855.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.9827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C)[CH][C]=O(2395)',
    structure = SMILES('[CH2]C(C)[CH][C]=O'),
    E0 = (272.647,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,1708.32],'cm^-1')),
        HinderedRotor(inertia=(0.174097,'amu*angstrom^2'), symmetry=1, barrier=(7.16066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00289944,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0028914,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.65508,'amu*angstrom^2'), symmetry=1, barrier=(68.2877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3452.89,'J/mol'), sigma=(6.01829,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=539.33 K, Pc=35.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38508,0.0522742,-4.37688e-05,2.08043e-08,-4.04955e-12,32890.4,25.3467], Tmin=(100,'K'), Tmax=(1232.68,'K')), NASAPolynomial(coeffs=[10.557,0.0225114,-7.55133e-06,1.21675e-09,-7.69674e-14,30629.2,-20.8212], Tmin=(1232.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.647,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'CC1CC=C1[O](9688)',
    structure = SMILES('CC1C[CH]C1=O'),
    E0 = (16.6585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43631,0.0154637,7.56211e-05,-1.0735e-07,4.14556e-11,2076.94,18.869], Tmin=(100,'K'), Tmax=(963.514,'K')), NASAPolynomial(coeffs=[12.6444,0.0208549,-7.13958e-06,1.3686e-09,-1.04449e-13,-2107.55,-41.5062], Tmin=(963.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(16.6585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O)"""),
)

species(
    label = 'C[CH]C=CC=O(3830)',
    structure = SMILES('CC=CC=C[O]'),
    E0 = (-9.43439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38932,0.0415479,1.478e-05,-5.78685e-08,2.77931e-11,-1026.12,20.4409], Tmin=(100,'K'), Tmax=(936.101,'K')), NASAPolynomial(coeffs=[17.682,0.0121592,-2.59283e-06,4.14127e-10,-3.31253e-14,-5839.1,-66.5008], Tmin=(936.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.43439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'C=CCC=C[O](6911)',
    structure = SMILES('C=CCC=C[O]'),
    E0 = (21.4197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61317,0.0368327,2.28023e-05,-6.1666e-08,2.7844e-11,2676.41,22.6599], Tmin=(100,'K'), Tmax=(948.605,'K')), NASAPolynomial(coeffs=[16.117,0.0147804,-4.16448e-06,7.44299e-10,-5.71671e-14,-1834.76,-55.8216], Tmin=(948.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.4197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]CC#C[O](9689)',
    structure = SMILES('C[CH]C[C]=C=O'),
    E0 = (271.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1685,370,2120,512.5,787.5,306.07,1989.85],'cm^-1')),
        HinderedRotor(inertia=(0.973019,'amu*angstrom^2'), symmetry=1, barrier=(64.6092,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0813198,'amu*angstrom^2'), symmetry=1, barrier=(5.40218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0812879,'amu*angstrom^2'), symmetry=1, barrier=(5.40149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32678,0.0419342,-2.07042e-05,-1.62898e-08,2.14865e-11,32680.8,21.5771], Tmin=(100,'K'), Tmax=(529.338,'K')), NASAPolynomial(coeffs=[4.60751,0.0322238,-1.45091e-05,2.76093e-09,-1.9327e-13,32333.9,11.029], Tmin=(529.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(RCCJC) + radical(CCCJ=C=O)"""),
)

species(
    label = 'C[CH]C[C]=C[O](9690)',
    structure = SMILES('C[CH]C[C]=C[O]'),
    E0 = (324.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,207.595,786.49,1689.13],'cm^-1')),
        HinderedRotor(inertia=(0.601732,'amu*angstrom^2'), symmetry=1, barrier=(13.835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602133,'amu*angstrom^2'), symmetry=1, barrier=(13.8442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0315643,'amu*angstrom^2'), symmetry=1, barrier=(13.8306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44674,0.0498142,-3.24375e-05,7.62536e-09,5.33056e-13,39144.4,25.9661], Tmin=(100,'K'), Tmax=(1065.61,'K')), NASAPolynomial(coeffs=[11.2955,0.0228573,-8.58564e-06,1.5206e-09,-1.03372e-13,36477,-24.8415], Tmin=(1065.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = 'CCC[C]=[C][O](9691)',
    structure = SMILES('CCC[C][C]=O'),
    E0 = (357.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,315.532,315.532,315.533],'cm^-1')),
        HinderedRotor(inertia=(0.00169323,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140811,'amu*angstrom^2'), symmetry=1, barrier=(9.94828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14081,'amu*angstrom^2'), symmetry=1, barrier=(9.94828,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140811,'amu*angstrom^2'), symmetry=1, barrier=(9.94828,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53741,0.0576467,-6.05113e-05,3.87996e-08,-1.06929e-11,43057.9,23.4555], Tmin=(100,'K'), Tmax=(860.784,'K')), NASAPolynomial(coeffs=[7.30743,0.0308332,-1.3785e-05,2.60978e-09,-1.81893e-13,42064.6,-3.51637], Tmin=(860.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C[C]=[C]O(9692)',
    structure = SMILES('C[CH]C[C]=[C]O'),
    E0 = (422.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1670,1700,300,440,352.869,361.531],'cm^-1')),
        HinderedRotor(inertia=(0.00133564,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0898391,'amu*angstrom^2'), symmetry=1, barrier=(8.08814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0884066,'amu*angstrom^2'), symmetry=1, barrier=(8.09899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0884277,'amu*angstrom^2'), symmetry=1, barrier=(8.11771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40356,0.0566658,-5.61511e-05,3.17025e-08,-7.35288e-12,50961.5,27.8407], Tmin=(100,'K'), Tmax=(1037.18,'K')), NASAPolynomial(coeffs=[9.78373,0.0243466,-9.41003e-06,1.65864e-09,-1.11147e-13,49223.2,-12.8946], Tmin=(1037.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJC) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C[CH][CH]C=[C]O(9693)',
    structure = SMILES('C[CH][CH]C=[C]O'),
    E0 = (326.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,1685,370,272.403,2105.9],'cm^-1')),
        HinderedRotor(inertia=(0.323363,'amu*angstrom^2'), symmetry=1, barrier=(17.049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00226844,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00227028,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701904,'amu*angstrom^2'), symmetry=1, barrier=(36.941,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50904,0.0483553,-2.78423e-05,1.1728e-09,3.42838e-12,39329.2,25.8794], Tmin=(100,'K'), Tmax=(988.148,'K')), NASAPolynomial(coeffs=[11.412,0.0216858,-7.72604e-06,1.3425e-09,-9.10989e-14,36717,-25.0935], Tmin=(988.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJC) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]CC=[C]O(9694)',
    structure = SMILES('[CH2][CH]CC=[C]O'),
    E0 = (390.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,763.521,3850.7],'cm^-1')),
        HinderedRotor(inertia=(0.00964503,'amu*angstrom^2'), symmetry=1, barrier=(3.98877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.593253,'amu*angstrom^2'), symmetry=1, barrier=(13.64,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.593198,'amu*angstrom^2'), symmetry=1, barrier=(13.6388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.75275,'amu*angstrom^2'), symmetry=1, barrier=(86.2832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30653,0.0547918,-4.88375e-05,2.38044e-08,-4.69571e-12,47048,29.4997], Tmin=(100,'K'), Tmax=(1220.57,'K')), NASAPolynomial(coeffs=[11.5397,0.0212557,-7.62329e-06,1.29326e-09,-8.48713e-14,44550,-21.909], Tmin=(1220.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJC) + radical(RCCJ) + radical(C=CJO)"""),
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
    E0 = (271.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (734.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (765.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (740.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (736.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (763.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (278.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (293.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (293.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (296.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (384.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (488.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (431.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (427.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (679.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (688.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (430.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (461.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (433.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (424.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (439.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (365.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (271.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (432.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (278.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (334.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (296.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (498.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (472.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (528.793,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (498.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (573.283,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (410.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (482.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['HCCO(2227)', 'C3H6(27)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CHCH3(T)(21)', '[CH2][CH][C]=O(9632)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C][O](6861)', 'C3H6(T)(28)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH3(17)', '[CH]C[CH][C]=O(9634)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C[C]C[CH][C]=O(9682)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C[CH]C[C][C]=O(9683)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['CC1CC1[C]=O(9668)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['CCC=C[C]=O(9684)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['CC=CC[C]=O(2431)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['C=CCC[C]=O(2291)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', 'CC=C[CH][C]=O(9685)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH2][CH]CC=C=O(2614)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C][O](6861)', 'C3H6(27)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.06e+08,'cm^3/(mol*s)'), n=1.58, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.4 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['CHCH3(T)(21)', 'CH2CHCO(3668)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'C[CH][CH][CH][C]=O(9686)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2][CH]C[CH][C]=O(2616)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC[CH][CH][C]=O(9665)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[CH][CH]C[C]=O(2433)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC[CH][C]=O(2618)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['C[CH][CH][CH]C=O(9687)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.21412e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]CC[C]=O(2285)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]CC=C[O](1849)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5Hall;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['C[CH]CC=C=O(2432)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C)[CH][C]=O(2395)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['CC1CC=C1[O](9688)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['C[CH]C=CC=O(3830)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C[CH]C[CH][C]=O(2434)'],
    products = ['C=CCC=C[O](6911)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', 'C[CH]CC#C[O](9689)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['HCCO(2227)', 'C3H6(T)(28)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C[CH]C[C]=C[O](9690)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CCC[C]=[C][O](9691)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C[CH]C[C]=[C]O(9692)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[CH][CH]C=[C]O(9693)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]CC=[C]O(9694)'],
    products = ['C[CH]C[CH][C]=O(2434)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.54044,'s^-1'), n=3.06876, Ea=(92.3667,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_2H;O_H_out] + [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6Hall;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2276',
    isomers = [
        'C[CH]C[CH][C]=O(2434)',
    ],
    reactants = [
        ('HCCO(2227)', 'C3H6(27)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2276',
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


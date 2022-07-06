species(
    label = 'C[CH]C(C)[CH][C]=O(2720)',
    structure = SMILES('C[CH]C(C)[CH][C]=O'),
    E0 = (241.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,1855,455,950,180,461.499],'cm^-1')),
        HinderedRotor(inertia=(0.0789565,'amu*angstrom^2'), symmetry=1, barrier=(1.81537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0789769,'amu*angstrom^2'), symmetry=1, barrier=(1.81583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0789442,'amu*angstrom^2'), symmetry=1, barrier=(1.81508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0120114,'amu*angstrom^2'), symmetry=1, barrier=(1.81535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0120106,'amu*angstrom^2'), symmetry=1, barrier=(1.81526,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25858,0.0601113,-4.25222e-05,1.64913e-08,-2.74771e-12,29165.2,28.8169], Tmin=(100,'K'), Tmax=(1346.5,'K')), NASAPolynomial(coeffs=[9.69046,0.035063,-1.46185e-05,2.67593e-09,-1.82655e-13,26894.5,-14.3707], Tmin=(1346.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(CCCJ=O)"""),
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
    label = 'butene2t(396)',
    structure = SMILES('CC=CC'),
    E0 = (-28.2774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.464594,'amu*angstrom^2'), symmetry=1, barrier=(10.6819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.460586,'amu*angstrom^2'), symmetry=1, barrier=(10.5898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52,0.0276811,9.44413e-07,-1.25016e-08,4.67173e-12,-3343.75,13.0195], Tmin=(100,'K'), Tmax=(1157.77,'K')), NASAPolynomial(coeffs=[6.04212,0.0260075,-1.04845e-05,1.90884e-09,-1.30603e-13,-4862.7,-7.52639], Tmin=(1157.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.2774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene2t""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH2(S)(11)',
    structure = SMILES('[CH2]'),
    E0 = (418.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1381.44,2681.05,3148.03],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08586,-0.00146571,5.70225e-06,-3.92867e-09,8.88752e-13,50365.3,-0.343718], Tmin=(100,'K'), Tmax=(1285.36,'K')), NASAPolynomial(coeffs=[2.65113,0.0038016,-1.38109e-06,2.30893e-10,-1.47438e-14,50667.9,6.68035], Tmin=(1285.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(S)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3428.95,'J/mol'), sigma=(5.98513,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.59 K, Pc=36.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23108,0.0440156,-2.01448e-05,-2.09479e-08,2.52039e-11,32658.9,23.6185], Tmin=(100,'K'), Tmax=(539.722,'K')), NASAPolynomial(coeffs=[4.99773,0.0326275,-1.38305e-05,2.54716e-09,-1.74683e-13,32227.4,10.7471], Tmin=(539.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'CC(C)[CH][CH][C]=O(10143)',
    structure = SMILES('CC(C)[CH]C=[C][O]'),
    E0 = (240.403,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,180,609.204,613.915],'cm^-1')),
        HinderedRotor(inertia=(0.00608863,'amu*angstrom^2'), symmetry=1, barrier=(1.6318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00759793,'amu*angstrom^2'), symmetry=1, barrier=(1.99326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.858157,'amu*angstrom^2'), symmetry=1, barrier=(20.7528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.902453,'amu*angstrom^2'), symmetry=1, barrier=(20.7512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00764,0.0543627,-1.18433e-05,-2.17695e-08,1.17227e-11,29031.6,27.6268], Tmin=(100,'K'), Tmax=(1013.58,'K')), NASAPolynomial(coeffs=[13.7861,0.028898,-1.11024e-05,2.04256e-09,-1.4396e-13,25158.8,-40.5202], Tmin=(1013.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = 'C[CH][CH]C(C)[C]=O(10144)',
    structure = SMILES('C[CH][CH]C(C)[C]=O'),
    E0 = (272.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,1855,455,950,439.411,1678.93],'cm^-1')),
        HinderedRotor(inertia=(0.000873088,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0548464,'amu*angstrom^2'), symmetry=1, barrier=(7.51479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0548461,'amu*angstrom^2'), symmetry=1, barrier=(7.51479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0548466,'amu*angstrom^2'), symmetry=1, barrier=(7.51479,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00375684,'amu*angstrom^2'), symmetry=1, barrier=(7.51479,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56626,0.057089,-3.97769e-05,1.63007e-08,-3.07651e-12,32872.3,31.1149], Tmin=(100,'K'), Tmax=(1143.62,'K')), NASAPolynomial(coeffs=[6.44177,0.0400361,-1.74097e-05,3.26186e-09,-2.26146e-13,31757.1,6.93915], Tmin=(1143.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJC) + radical(CC(C)CJ=O)"""),
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
    label = 'C[CH][CH][C]=O(9638)',
    structure = SMILES('C[CH]C=[C][O]'),
    E0 = (297.005,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,484.767,484.8],'cm^-1')),
        HinderedRotor(inertia=(0.000717709,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163781,'amu*angstrom^2'), symmetry=1, barrier=(27.3151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43094,0.0256217,1.17561e-05,-3.33743e-08,1.45333e-11,35785.9,19.0015], Tmin=(100,'K'), Tmax=(974.935,'K')), NASAPolynomial(coeffs=[10.1488,0.0154088,-5.53624e-06,1.01957e-09,-7.35997e-14,33261.5,-23.2652], Tmin=(974.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(C=CJO)"""),
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
    label = 'C[CH][CH]C(1186)',
    structure = SMILES('C[CH][CH]C'),
    E0 = (244.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2670.04],'cm^-1')),
        HinderedRotor(inertia=(0.341695,'amu*angstrom^2'), symmetry=1, barrier=(7.85623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00155291,'amu*angstrom^2'), symmetry=1, barrier=(7.85531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136081,'amu*angstrom^2'), symmetry=1, barrier=(68.8251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62536,0.022997,-1.74547e-06,-3.25078e-09,7.61424e-13,29450.7,14.8448], Tmin=(100,'K'), Tmax=(2137.09,'K')), NASAPolynomial(coeffs=[11.4859,0.0201914,-8.13358e-06,1.34908e-09,-8.16544e-14,23371.9,-35.4092], Tmin=(2137.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC)"""),
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
    label = '[CH]C(C)[CH][C]=O(9667)',
    structure = SMILES('[CH]C(C)[CH][C]=O'),
    E0 = (515.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,383.587,383.702,1646.93,2343.56],'cm^-1')),
        HinderedRotor(inertia=(0.13425,'amu*angstrom^2'), symmetry=1, barrier=(14.0993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0011367,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00112417,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.711103,'amu*angstrom^2'), symmetry=1, barrier=(74.7345,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42311,0.0523556,-4.62266e-05,2.20041e-08,-4.23731e-12,62130.4,24.4335], Tmin=(100,'K'), Tmax=(1247.21,'K')), NASAPolynomial(coeffs=[11.4677,0.0201404,-7.48114e-06,1.29332e-09,-8.58157e-14,59624.9,-26.2448], Tmin=(1247.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = 'C[C]C(C)[CH][C]=O(10145)',
    structure = SMILES('C[C]C(C)[CH][C]=O'),
    E0 = (495.347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82903,0.0665415,-5.80793e-05,2.76138e-08,-5.38313e-12,59693.3,27.7114], Tmin=(100,'K'), Tmax=(1220.09,'K')), NASAPolynomial(coeffs=[12.4565,0.0284208,-1.12121e-05,2.00466e-09,-1.35629e-13,56856,-30.6969], Tmin=(1220.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C(C)[C][C]=O(10146)',
    structure = SMILES('C[CH]C(C)[C][C]=O'),
    E0 = (522.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10512,0.0685392,-7.50787e-05,5.19402e-08,-1.56318e-11,62925.2,29.5305], Tmin=(100,'K'), Tmax=(788.628,'K')), NASAPolynomial(coeffs=[7.05989,0.0383359,-1.76304e-05,3.37605e-09,-2.36536e-13,61986,2.21625], Tmin=(788.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'CC1C(C)C1[C]=O(10147)',
    structure = SMILES('CC1C(C)C1[C]=O'),
    E0 = (22.7776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42597,0.0405962,2.88252e-05,-6.57005e-08,2.78097e-11,2846.54,23.9901], Tmin=(100,'K'), Tmax=(972.688,'K')), NASAPolynomial(coeffs=[14.3907,0.0263183,-9.35642e-06,1.72871e-09,-1.25559e-13,-1522.28,-47.691], Tmin=(972.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.7776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'CCC(C)=C[C]=O(9735)',
    structure = SMILES('CCC(C)=C[C]=O'),
    E0 = (-18.4434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02287,0.0719152,-6.94429e-05,3.99452e-08,-1.01529e-11,-2116.47,24.146], Tmin=(100,'K'), Tmax=(910.307,'K')), NASAPolynomial(coeffs=[7.79752,0.0421468,-2.03912e-05,4.02242e-09,-2.87469e-13,-3349.88,-7.90122], Tmin=(910.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.4434,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC=C(C)C[C]=O(2716)',
    structure = SMILES('CC=C(C)C[C]=O'),
    E0 = (-8.25877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,237.634],'cm^-1')),
        HinderedRotor(inertia=(0.0028318,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2945,'amu*angstrom^2'), symmetry=1, barrier=(11.9195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288366,'amu*angstrom^2'), symmetry=1, barrier=(11.9321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339467,'amu*angstrom^2'), symmetry=1, barrier=(11.9145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995304,0.059308,-3.62871e-05,1.01549e-08,-1.11022e-12,-879.273,25.8764], Tmin=(100,'K'), Tmax=(2114.34,'K')), NASAPolynomial(coeffs=[20.9521,0.0215535,-9.50297e-06,1.70984e-09,-1.11695e-13,-9318.54,-85.3467], Tmin=(2114.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.25877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC(C)C[C]=O(2718)',
    structure = SMILES('C=CC(C)C[C]=O'),
    E0 = (5.10705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,294.781,294.781],'cm^-1')),
        HinderedRotor(inertia=(0.00193998,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158066,'amu*angstrom^2'), symmetry=1, barrier=(9.74687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158069,'amu*angstrom^2'), symmetry=1, barrier=(9.74687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158069,'amu*angstrom^2'), symmetry=1, barrier=(9.74686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03891,0.0615754,-4.38931e-05,1.6588e-08,-2.61263e-12,723.49,27.8013], Tmin=(100,'K'), Tmax=(1454.03,'K')), NASAPolynomial(coeffs=[12.0242,0.0313552,-1.27177e-05,2.29427e-09,-1.55041e-13,-2471.13,-29.3089], Tmin=(1454.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.10705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=C(C)[CH][C]=O(10148)',
    structure = SMILES('C[CH]C(C)=C[C]=O'),
    E0 = (122.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.614999,'amu*angstrom^2'), symmetry=1, barrier=(14.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.615229,'amu*angstrom^2'), symmetry=1, barrier=(14.1453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.614949,'amu*angstrom^2'), symmetry=1, barrier=(14.1389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0656154,'amu*angstrom^2'), symmetry=1, barrier=(14.1434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26792,0.0662243,-5.92347e-05,3.02767e-08,-6.8298e-12,14847,22.154], Tmin=(100,'K'), Tmax=(1010,'K')), NASAPolynomial(coeffs=[8.09383,0.0391911,-1.90864e-05,3.77623e-09,-2.70289e-13,13468.1,-10.845], Tmin=(1010,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_S) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2][CH]C(C)C=C=O(4176)',
    structure = SMILES('[CH2][CH]C(C)C=C=O'),
    E0 = (241.131,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2120,512.5,787.5,344.134,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0106282,'amu*angstrom^2'), symmetry=1, barrier=(26.2918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0919169,'amu*angstrom^2'), symmetry=1, barrier=(7.72529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0919256,'amu*angstrom^2'), symmetry=1, barrier=(7.72529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312855,'amu*angstrom^2'), symmetry=1, barrier=(26.2917,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4928,0.058788,-4.53704e-05,1.9742e-08,-3.79184e-12,29088.5,28.2399], Tmin=(100,'K'), Tmax=(1162.17,'K')), NASAPolynomial(coeffs=[8.00864,0.0363616,-1.64249e-05,3.13769e-09,-2.20012e-13,27574,-4.17453], Tmin=(1162.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'CC=C[C]=O(9636)',
    structure = SMILES('C[CH]C=C=O'),
    E0 = (39.4293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,1069.02],'cm^-1')),
        HinderedRotor(inertia=(0.163668,'amu*angstrom^2'), symmetry=1, barrier=(3.76304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30969,'amu*angstrom^2'), symmetry=1, barrier=(30.1123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47551,0.0278456,-3.33541e-06,-1.12635e-08,5.08824e-12,4801.91,15.1824], Tmin=(100,'K'), Tmax=(1090.48,'K')), NASAPolynomial(coeffs=[7.95588,0.0195027,-8.03539e-06,1.49913e-09,-1.04907e-13,2907.47,-14.9378], Tmin=(1090.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.4293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Allyl_S)"""),
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
    label = 'C[CH][C](C)[CH][C]=O(10149)',
    structure = SMILES('C[CH][C](C)C=[C][O]'),
    E0 = (426.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,349.54,349.653,349.662],'cm^-1')),
        HinderedRotor(inertia=(4.10431e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00137927,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00137888,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.46318,'amu*angstrom^2'), symmetry=1, barrier=(40.1616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21737,0.0533312,-2.70727e-05,1.70212e-09,2.02869e-12,51448.9,29.5223], Tmin=(100,'K'), Tmax=(1134.45,'K')), NASAPolynomial(coeffs=[11.4077,0.0295571,-1.17118e-05,2.12105e-09,-1.45249e-13,48354.6,-24.3734], Tmin=(1134.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(Cs_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([CH]C)[CH][C]=O(4240)',
    structure = SMILES('[CH2]C([CH]C)[CH][C]=O'),
    E0 = (446.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,957.83,958.279],'cm^-1')),
        HinderedRotor(inertia=(0.0744504,'amu*angstrom^2'), symmetry=1, barrier=(1.71176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749299,'amu*angstrom^2'), symmetry=1, barrier=(1.72279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749341,'amu*angstrom^2'), symmetry=1, barrier=(1.72288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0744698,'amu*angstrom^2'), symmetry=1, barrier=(1.71221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0748947,'amu*angstrom^2'), symmetry=1, barrier=(1.72198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2049,0.0625364,-5.72642e-05,3.16551e-08,-7.47846e-12,53832.1,30.8], Tmin=(100,'K'), Tmax=(999.296,'K')), NASAPolynomial(coeffs=[8.4989,0.0333397,-1.3438e-05,2.41688e-09,-1.63698e-13,52374.4,-4.38418], Tmin=(999.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C(C)[CH][C]=O(4178)',
    structure = SMILES('[CH2][CH]C(C)[CH][C]=O'),
    E0 = (446.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,850.486,850.567],'cm^-1')),
        HinderedRotor(inertia=(0.122248,'amu*angstrom^2'), symmetry=1, barrier=(2.81072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122387,'amu*angstrom^2'), symmetry=1, barrier=(2.81391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00547187,'amu*angstrom^2'), symmetry=1, barrier=(2.80925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122372,'amu*angstrom^2'), symmetry=1, barrier=(2.81356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122311,'amu*angstrom^2'), symmetry=1, barrier=(2.81218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37781,0.0597816,-4.89364e-05,2.32743e-08,-4.82459e-12,53844.7,30.1838], Tmin=(100,'K'), Tmax=(1105.4,'K')), NASAPolynomial(coeffs=[8.29099,0.0347652,-1.49891e-05,2.80038e-09,-1.94078e-13,52316.4,-3.86095], Tmin=(1105.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[C](C)[CH][C]=O(9738)',
    structure = SMILES('CC[C](C)C=[C][O]'),
    E0 = (232.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08131,0.0539849,-1.41161e-05,-1.78653e-08,1.0211e-11,28058.2,27.4495], Tmin=(100,'K'), Tmax=(1006.89,'K')), NASAPolynomial(coeffs=[12.6996,0.0299354,-1.12204e-05,2.0215e-09,-1.40376e-13,24597.9,-34.2466], Tmin=(1006.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(C=CJO)"""),
)

species(
    label = 'C[CH][C](C)C[C]=O(2719)',
    structure = SMILES('C[CH][C](C)C[C]=O'),
    E0 = (259.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,1855,455,950,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.102702,'amu*angstrom^2'), symmetry=1, barrier=(2.36132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102692,'amu*angstrom^2'), symmetry=1, barrier=(2.36109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102669,'amu*angstrom^2'), symmetry=1, barrier=(2.36057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102644,'amu*angstrom^2'), symmetry=1, barrier=(2.35999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102694,'amu*angstrom^2'), symmetry=1, barrier=(2.36114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05689,0.0753737,-0.000112452,1.09476e-07,-4.18903e-11,31314.4,30.916], Tmin=(100,'K'), Tmax=(836.988,'K')), NASAPolynomial(coeffs=[1.40105,0.049758,-2.35859e-05,4.47624e-09,-3.07506e-13,32096.4,34.3328], Tmin=(836.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]CC(C)[CH][C]=O(4182)',
    structure = SMILES('[CH2]CC(C)[CH][C]=O'),
    E0 = (252.378,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,1003.3,1003.99],'cm^-1')),
        HinderedRotor(inertia=(0.147649,'amu*angstrom^2'), symmetry=1, barrier=(3.39473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.640907,'amu*angstrom^2'), symmetry=1, barrier=(14.7357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.641053,'amu*angstrom^2'), symmetry=1, barrier=(14.7391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000611004,'amu*angstrom^2'), symmetry=1, barrier=(3.38286,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.51184,'amu*angstrom^2'), symmetry=1, barrier=(80.744,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.897881,0.0644572,-4.97697e-05,2.099e-08,-3.68566e-12,30468.8,29.3453], Tmin=(100,'K'), Tmax=(1328,'K')), NASAPolynomial(coeffs=[11.8998,0.0313182,-1.23378e-05,2.19839e-09,-1.4801e-13,27546.7,-26.8531], Tmin=(1328,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH][C]=O)CC(2468)',
    structure = SMILES('[CH2]C([CH][C]=O)CC'),
    E0 = (252.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,327.303,2711.63],'cm^-1')),
        HinderedRotor(inertia=(1.01456,'amu*angstrom^2'), symmetry=1, barrier=(77.1271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00157363,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168326,'amu*angstrom^2'), symmetry=1, barrier=(12.7961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00245239,'amu*angstrom^2'), symmetry=1, barrier=(12.7961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01457,'amu*angstrom^2'), symmetry=1, barrier=(77.1272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79021,0.0664689,-5.56502e-05,2.64471e-08,-5.20762e-12,30453.4,29.7261], Tmin=(100,'K'), Tmax=(1206.53,'K')), NASAPolynomial(coeffs=[11.5566,0.0307744,-1.12727e-05,1.92592e-09,-1.26606e-13,27855.4,-24.2368], Tmin=(1206.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH]C)C[C]=O(2467)',
    structure = SMILES('[CH2]C([CH]C)C[C]=O'),
    E0 = (279.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,413.61,2621.97],'cm^-1')),
        HinderedRotor(inertia=(0.0627707,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627708,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627705,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627707,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0141566,'amu*angstrom^2'), symmetry=1, barrier=(69.0623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05704,0.0686662,-7.37601e-05,5.27733e-08,-1.65495e-11,33685.5,31.5724], Tmin=(100,'K'), Tmax=(762.453,'K')), NASAPolynomial(coeffs=[6.46712,0.0401417,-1.73637e-05,3.21791e-09,-2.20793e-13,32864.7,6.96619], Tmin=(762.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][C](C)[CH]C=O(10150)',
    structure = SMILES('C[CH][C](C)C=C[O]'),
    E0 = (187.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08161,0.0507928,2.0476e-06,-3.91745e-08,1.87473e-11,22624.6,26.8755], Tmin=(100,'K'), Tmax=(978.556,'K')), NASAPolynomial(coeffs=[14.6381,0.0269164,-9.69625e-06,1.76145e-09,-1.25196e-13,18461.5,-45.9482], Tmin=(978.556,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(187.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C([CH]C)[CH]C=O(4241)',
    structure = SMILES('[CH2]C([CH]C)C=C[O]'),
    E0 = (260.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,658.043,658.729],'cm^-1')),
        HinderedRotor(inertia=(0.106045,'amu*angstrom^2'), symmetry=1, barrier=(2.43818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00793951,'amu*angstrom^2'), symmetry=1, barrier=(2.45055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.762056,'amu*angstrom^2'), symmetry=1, barrier=(17.5212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0630873,'amu*angstrom^2'), symmetry=1, barrier=(19.5191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.831755,0.057764,-1.55562e-05,-2.43502e-08,1.46598e-11,31424.5,31.3003], Tmin=(100,'K'), Tmax=(959.338,'K')), NASAPolynomial(coeffs=[15.5669,0.0247168,-8.27676e-06,1.44056e-09,-1.00385e-13,27290.8,-45.9857], Tmin=(959.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(C)C[C]=O(2721)',
    structure = SMILES('[CH2][CH]C(C)C[C]=O'),
    E0 = (279.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,1904.05,1904.05],'cm^-1')),
        HinderedRotor(inertia=(0.0839111,'amu*angstrom^2'), symmetry=1, barrier=(8.91661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00112577,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0839118,'amu*angstrom^2'), symmetry=1, barrier=(8.91662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0839121,'amu*angstrom^2'), symmetry=1, barrier=(8.91661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0839115,'amu*angstrom^2'), symmetry=1, barrier=(8.91662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23781,0.0656711,-6.39225e-05,4.1548e-08,-1.23157e-11,33698,30.9379], Tmin=(100,'K'), Tmax=(784.457,'K')), NASAPolynomial(coeffs=[5.8029,0.0423938,-1.94138e-05,3.72327e-09,-2.61537e-13,32981.8,10.0222], Tmin=(784.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C(C)[CH]C=O(4181)',
    structure = SMILES('[CH2][CH]C(C)C=C[O]'),
    E0 = (260.404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,461.019,461.028,461.031],'cm^-1')),
        HinderedRotor(inertia=(0.000793129,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000793147,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110754,'amu*angstrom^2'), symmetry=1, barrier=(16.7042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110748,'amu*angstrom^2'), symmetry=1, barrier=(16.7042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.835648,0.0569854,-1.40792e-05,-2.38866e-08,1.35378e-11,31444.5,31.2917], Tmin=(100,'K'), Tmax=(994.216,'K')), NASAPolynomial(coeffs=[15.6861,0.025602,-9.52352e-06,1.75347e-09,-1.24994e-13,27089.8,-47.3167], Tmin=(994.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]C(C)C=C=O(2717)',
    structure = SMILES('C[CH]C(C)C=C=O'),
    E0 = (35.8846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,308.584,308.584],'cm^-1')),
        HinderedRotor(inertia=(0.00177037,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00177033,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127888,'amu*angstrom^2'), symmetry=1, barrier=(8.64178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127889,'amu*angstrom^2'), symmetry=1, barrier=(8.64176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36562,0.0592537,-3.95434e-05,1.37604e-08,-2.04111e-12,4409.05,26.8974], Tmin=(100,'K'), Tmax=(1476.41,'K')), NASAPolynomial(coeffs=[10.128,0.0355141,-1.54246e-05,2.86966e-09,-1.96994e-13,1821.67,-18.79], Tmin=(1476.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.8846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Cs_S)"""),
)

species(
    label = 'CC1C=C([O])C1C(10151)',
    structure = SMILES('CC1[CH]C(=O)C1C'),
    E0 = (-16.1633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75879,0.0282498,6.99705e-05,-1.09135e-07,4.30574e-11,-1844.35,22.5723], Tmin=(100,'K'), Tmax=(967.355,'K')), NASAPolynomial(coeffs=[14.8149,0.0268575,-9.42434e-06,1.78502e-09,-1.33595e-13,-6831.15,-52.7017], Tmin=(967.355,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.1633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O)"""),
)

species(
    label = 'CC=C(C)C=C[O](10152)',
    structure = SMILES('CC=C(C)C=C[O]'),
    E0 = (-47.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547615,0.0615355,-1.57712e-05,-3.07043e-08,1.82956e-11,-5518.45,24.0622], Tmin=(100,'K'), Tmax=(954.268,'K')), NASAPolynomial(coeffs=[18.5934,0.0207736,-6.52603e-06,1.14053e-09,-8.19117e-14,-10550.7,-70.4746], Tmin=(954.268,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-47.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'C=CC(C)C=C[O](6952)',
    structure = SMILES('C=CC(C)C=C[O]'),
    E0 = (-12.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.667065,0.0613493,-2.36052e-05,-1.60036e-08,1.13216e-11,-1330.69,25.642], Tmin=(100,'K'), Tmax=(986.338,'K')), NASAPolynomial(coeffs=[16.217,0.025028,-9.03419e-06,1.63365e-09,-1.15387e-13,-5698.9,-55.7571], Tmin=(986.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]C(C)C#C[O](10153)',
    structure = SMILES('C[CH]C(C)[C]=C=O'),
    E0 = (241.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,1685,370,2120,512.5,787.5,180,2703.86],'cm^-1')),
        HinderedRotor(inertia=(0.510122,'amu*angstrom^2'), symmetry=1, barrier=(11.7287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.51009,'amu*angstrom^2'), symmetry=1, barrier=(11.728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186307,'amu*angstrom^2'), symmetry=1, barrier=(11.7266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254272,'amu*angstrom^2'), symmetry=1, barrier=(36.8237,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40698,0.0618312,-5.36684e-05,2.84962e-08,-6.8459e-12,29156.3,26.6047], Tmin=(100,'K'), Tmax=(946.817,'K')), NASAPolynomial(coeffs=[6.71004,0.0394278,-1.81759e-05,3.50573e-09,-2.47429e-13,28152.1,1.31032], Tmin=(946.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Cs_S) + radical(CC(C)CJ=C=O)"""),
)

species(
    label = 'C[CH]C(C)[C]=C[O](10154)',
    structure = SMILES('C[CH]C(C)[C]=C[O]'),
    E0 = (293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,455.838,455.846,455.849],'cm^-1')),
        HinderedRotor(inertia=(0.0908835,'amu*angstrom^2'), symmetry=1, barrier=(13.4013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000811247,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0908838,'amu*angstrom^2'), symmetry=1, barrier=(13.4013,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0908817,'amu*angstrom^2'), symmetry=1, barrier=(13.4013,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.744991,0.0610812,-2.91563e-05,-6.01617e-09,6.69874e-12,35366.1,30.3049], Tmin=(100,'K'), Tmax=(1025.82,'K')), NASAPolynomial(coeffs=[14.8934,0.0270798,-1.039e-05,1.90319e-09,-1.33488e-13,31349.6,-43.7423], Tmin=(1025.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S) + radical(Cds_S)"""),
)

species(
    label = 'CCC(C)[C]=[C][O](9744)',
    structure = SMILES('CCC(C)[C][C]=O'),
    E0 = (327.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961983,0.069288,-6.24035e-05,3.23947e-08,-7.18678e-12,39534.7,27.4815], Tmin=(100,'K'), Tmax=(1053.12,'K')), NASAPolynomial(coeffs=[9.54499,0.0366875,-1.5969e-05,2.99963e-09,-2.08657e-13,37726.9,-14.3707], Tmin=(1053.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C(C)[C]=[C]O(10155)',
    structure = SMILES('C[CH]C(C)[C]=[C]O'),
    E0 = (391.281,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,1670,1700,300,440,180,355.99],'cm^-1')),
        HinderedRotor(inertia=(0.00014563,'amu*angstrom^2'), symmetry=1, barrier=(1.65349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183459,'amu*angstrom^2'), symmetry=1, barrier=(16.4638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183093,'amu*angstrom^2'), symmetry=1, barrier=(16.4638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.716043,'amu*angstrom^2'), symmetry=1, barrier=(16.4632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183017,'amu*angstrom^2'), symmetry=1, barrier=(16.4633,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.508362,0.0701937,-6.0611e-05,2.77377e-08,-5.11582e-12,47191.5,32.8741], Tmin=(100,'K'), Tmax=(1300.24,'K')), NASAPolynomial(coeffs=[14.7984,0.0262323,-9.89557e-06,1.73456e-09,-1.16124e-13,43475.4,-39.819], Tmin=(1300.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.281,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cs_S) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C[CH][C](C)C=[C]O(10156)',
    structure = SMILES('C[CH][C](C)C=[C]O'),
    E0 = (285.42,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,427.41,433.713],'cm^-1')),
        HinderedRotor(inertia=(0.000933728,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00084773,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0621925,'amu*angstrom^2'), symmetry=1, barrier=(8.42365,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0672785,'amu*angstrom^2'), symmetry=1, barrier=(8.52456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181448,'amu*angstrom^2'), symmetry=1, barrier=(23.6528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.859777,0.0597646,-2.91134e-05,-5.44709e-09,6.76775e-12,34449.4,29.3898], Tmin=(100,'K'), Tmax=(992.149,'K')), NASAPolynomial(coeffs=[13.8685,0.0271556,-9.80455e-06,1.7312e-09,-1.19066e-13,30891.7,-38.1875], Tmin=(992.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.42,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_T) + radical(Cs_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([CH]C)C=[C]O(9747)',
    structure = SMILES('[CH2]C([CH]C)C=[C]O'),
    E0 = (358.522,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0288564,'amu*angstrom^2'), symmetry=1, barrier=(13.2704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288454,'amu*angstrom^2'), symmetry=1, barrier=(13.2713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0082897,'amu*angstrom^2'), symmetry=1, barrier=(84.9638,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288512,'amu*angstrom^2'), symmetry=1, barrier=(13.2715,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.577207,'amu*angstrom^2'), symmetry=1, barrier=(13.2711,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622726,0.066583,-4.61787e-05,8.67914e-09,2.97541e-12,43248.7,33.7689], Tmin=(100,'K'), Tmax=(957.56,'K')), NASAPolynomial(coeffs=[14.7381,0.0250562,-8.44276e-06,1.42392e-09,-9.53821e-14,39746,-37.8919], Tmin=(957.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cs_S) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]C(C)C=[C]O(10157)',
    structure = SMILES('[CH2][CH]C(C)C=[C]O'),
    E0 = (358.686,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,998.812,1004.5],'cm^-1')),
        HinderedRotor(inertia=(0.781027,'amu*angstrom^2'), symmetry=1, barrier=(17.9573,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777723,'amu*angstrom^2'), symmetry=1, barrier=(17.8814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.781821,'amu*angstrom^2'), symmetry=1, barrier=(17.9756,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176733,'amu*angstrom^2'), symmetry=1, barrier=(4.06344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00158428,'amu*angstrom^2'), symmetry=1, barrier=(17.988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.595824,0.0661633,-4.59252e-05,1.06633e-08,1.24193e-12,43270.1,33.8709], Tmin=(100,'K'), Tmax=(1031.18,'K')), NASAPolynomial(coeffs=[15.0359,0.0256426,-9.51918e-06,1.69693e-09,-1.16704e-13,39468.3,-40.2317], Tmin=(1031.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.686,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Cs_S) + radical(RCCJ) + radical(C=CJO)"""),
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
    E0 = (241.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (689.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (435.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (471.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (696.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (725.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (707.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (707.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (734.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (249.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (263.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (264.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (266.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (345.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (455.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (400.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (383.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (325.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (604.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (638.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (658.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (658.712,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (378.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (417.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (400.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (399.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (397.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (394.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (391.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (410.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (334.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (241.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (249.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (305.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (266.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (468.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (432.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (497.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (469.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (541.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (369.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (432.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (451.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['HCCO(2227)', 'butene2t(396)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', 'C[CH]C[CH][C]=O(2434)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CC(C)[CH][CH][C]=O(10143)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.11838e+10,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;CH3]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C[CH][CH]C(C)[C]=O(10144)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.37952e+08,'s^-1'), n=1.37167, Ea=(199.228,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHCH3(T)(21)', 'C[CH][CH][C]=O(9638)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C][O](6861)', 'C[CH][CH]C(1186)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', '[CH]C(C)[CH][C]=O(9667)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C[C]C(C)[CH][C]=O(10145)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C[CH]C(C)[C][C]=O(10146)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['CC1C(C)C1[C]=O(10147)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['CCC(C)=C[C]=O(9735)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['CC=C(C)C[C]=O(2716)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['C=CC(C)C[C]=O(2718)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'CC=C(C)[CH][C]=O(10148)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2][CH]C(C)C=C=O(4176)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C][O](6861)', 'butene2t(396)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.01161,'m^3/(mol*s)'), n=1.87483, Ea=(3.76713,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds;Y_1centerbirad] + [Cds-CsH_Cds-CsH;YJ] for rate rule [Cds-CsH_Cds-CsH;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CHCH3(T)(21)', 'CC=C[C]=O(9636)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(17)', 'CC=C[CH][C]=O(9685)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00632077,'m^3/(mol*s)'), n=2.46811, Ea=(26.7526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH3(17)', 'C[CH][CH][CH][C]=O(9686)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.64e+07,'m^3/(mol*s)'), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', 'C[CH][C](C)[CH][C]=O(10149)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2]C([CH]C)[CH][C]=O(4240)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][CH]C(C)[CH][C]=O(4178)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['CC[C](C)[CH][C]=O(9738)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(956916,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C[CH][C](C)C[C]=O(2719)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]CC(C)[CH][C]=O(4182)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['C[CH][C](C)[CH]C=O(10150)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.10706e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C([CH]C)[CH]C=O(4241)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_2;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH]C(C)C[C]=O(2721)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]C(C)[CH]C=O(4181)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5Hall;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['C[CH]C(C)C=C=O(2717)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['CC1C=C([O])C1C(10151)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['CC=C(C)C=C[O](10152)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C[CH]C(C)[CH][C]=O(2720)'],
    products = ['C=CC(C)C=C[O](6952)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', 'C[CH]C(C)C#C[O](10153)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['HCCO(2227)', 'C[CH][CH]C(1186)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.04713,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C[CH]C(C)[C]=C[O](10154)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['CCC(C)[C]=[C][O](9744)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C[CH]C(C)[C]=[C]O(10155)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C[CH][C](C)C=[C]O(10156)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([CH]C)C=[C]O(9747)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][CH]C(C)C=[C]O(10157)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.54044,'s^-1'), n=3.06876, Ea=(92.3667,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_2H;O_H_out] + [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6Hall;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2361',
    isomers = [
        'C[CH]C(C)[CH][C]=O(2720)',
    ],
    reactants = [
        ('HCCO(2227)', 'butene2t(396)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2361',
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


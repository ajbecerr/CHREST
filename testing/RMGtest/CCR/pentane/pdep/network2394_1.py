species(
    label = '[O][CH]C[CH][C]=O(2758)',
    structure = SMILES('[O][CH]C[CH][C]=O'),
    E0 = (346.562,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,1855,455,950,180,1271.58,1272.67],'cm^-1')),
        HinderedRotor(inertia=(0.189104,'amu*angstrom^2'), symmetry=1, barrier=(4.34788,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187414,'amu*angstrom^2'), symmetry=1, barrier=(4.30901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188327,'amu*angstrom^2'), symmetry=1, barrier=(4.33001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62877,0.0609896,-0.000105677,9.96975e-08,-3.55197e-11,41758.7,23.8395], Tmin=(100,'K'), Tmax=(877.42,'K')), NASAPolynomial(coeffs=[4.49343,0.0273091,-1.2845e-05,2.37828e-09,-1.59104e-13,42049.8,14.9171], Tmin=(877.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.562,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCHO) + radical(CCsJOH) + radical(CCCJ=O)"""),
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
    label = 'O=[C][CH]CC=O(2652)',
    structure = SMILES('O=[C][CH]CC=O'),
    E0 = (49.9599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1855,455,950,2417.03],'cm^-1')),
        HinderedRotor(inertia=(0.098458,'amu*angstrom^2'), symmetry=1, barrier=(6.65102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0984391,'amu*angstrom^2'), symmetry=1, barrier=(6.65088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0984336,'amu*angstrom^2'), symmetry=1, barrier=(6.65069,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05098,0.0477841,-6.88496e-05,6.1783e-08,-2.23783e-11,6074.21,23.9783], Tmin=(100,'K'), Tmax=(820.155,'K')), NASAPolynomial(coeffs=[4.69391,0.0258595,-1.22273e-05,2.32583e-09,-1.60384e-13,5944.55,13.6041], Tmin=(820.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.9599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
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
    label = '[O][CH]C[C][C]=O(10176)',
    structure = SMILES('[O][CH]C[C][C]=O'),
    E0 = (627.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1855,455,950,180,180,180,1262.75],'cm^-1')),
        HinderedRotor(inertia=(0.20261,'amu*angstrom^2'), symmetry=1, barrier=(4.65841,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202786,'amu*angstrom^2'), symmetry=1, barrier=(4.66245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203329,'amu*angstrom^2'), symmetry=1, barrier=(4.67494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41641,0.0703576,-0.000142472,1.41589e-07,-5.14024e-11,75520.7,24.7459], Tmin=(100,'K'), Tmax=(883.036,'K')), NASAPolynomial(coeffs=[3.64405,0.0275186,-1.40731e-05,2.65368e-09,-1.77533e-13,76404,21.5054], Tmin=(883.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCsJOH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]C[CH][C]=O(10177)',
    structure = SMILES('[O][C]C[CH][C]=O'),
    E0 = (627.253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1855,455,950,180,180,914.653,919.741],'cm^-1')),
        HinderedRotor(inertia=(0.17601,'amu*angstrom^2'), symmetry=1, barrier=(4.04682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1809,'amu*angstrom^2'), symmetry=1, barrier=(4.15924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182383,'amu*angstrom^2'), symmetry=1, barrier=(4.19335,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66506,0.0584358,-0.00010046,9.1737e-08,-3.20248e-11,75518.4,22.892], Tmin=(100,'K'), Tmax=(863.403,'K')), NASAPolynomial(coeffs=[6.25246,0.0219716,-1.06831e-05,2.01138e-09,-1.36111e-13,75293.2,4.71764], Tmin=(863.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCHO) + radical(CH2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1CC1[C]=O(10178)',
    structure = SMILES('[O]C1CC1[C]=O'),
    E0 = (146.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20641,0.0287139,1.11353e-05,-3.77954e-08,1.71666e-11,17671.5,20.119], Tmin=(100,'K'), Tmax=(966.403,'K')), NASAPolynomial(coeffs=[12.4567,0.0124414,-4.20235e-06,7.89399e-10,-5.93484e-14,14469,-35.3014], Tmin=(966.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=[C]C[CH]C=O(2651)',
    structure = SMILES('[O]C=CC[C]=O'),
    E0 = (-0.508181,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,189.324,190.206],'cm^-1')),
        HinderedRotor(inertia=(0.840078,'amu*angstrom^2'), symmetry=1, barrier=(21.5362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.85679,'amu*angstrom^2'), symmetry=1, barrier=(21.5267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88818,0.0328175,9.85638e-06,-4.13901e-08,1.88828e-11,27.0453,21.3343], Tmin=(100,'K'), Tmax=(992.314,'K')), NASAPolynomial(coeffs=[15.7297,0.00949103,-3.96267e-06,8.67221e-10,-7.04814e-14,-4318.53,-53.3907], Tmin=(992.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-0.508181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]CC=C[C]=O(10179)',
    structure = SMILES('[O]C[CH]C=C=O'),
    E0 = (81.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84174,0.0523848,-7.25779e-05,6.41991e-08,-2.36911e-11,9828.37,19.3234], Tmin=(100,'K'), Tmax=(770.557,'K')), NASAPolynomial(coeffs=[5.08642,0.0290028,-1.43332e-05,2.79491e-09,-1.96264e-13,9522.45,5.77495], Tmin=(770.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(81.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O][C]=CC=C[O](9609)',
    structure = SMILES('[O]C=C[CH][C]=O'),
    E0 = (161.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1855,455,950,242.424,244.574],'cm^-1')),
        HinderedRotor(inertia=(0.688074,'amu*angstrom^2'), symmetry=1, barrier=(29.9476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.692779,'amu*angstrom^2'), symmetry=1, barrier=(29.9454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99504,0.025578,3.51523e-05,-7.39183e-08,3.2058e-11,19527.4,21.3351], Tmin=(100,'K'), Tmax=(963.766,'K')), NASAPolynomial(coeffs=[18.6873,0.0023414,-3.43071e-07,2.0471e-10,-2.785e-14,14171.6,-69.6729], Tmin=(963.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCJC=O) + radical(CCCJ=O)"""),
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
    label = '[O][CH][CH][CH][C]=O(10180)',
    structure = SMILES('[O][C]=C[CH][CH][O]'),
    E0 = (518.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,1685,370,362.805,362.891,362.917,1140.43],'cm^-1')),
        HinderedRotor(inertia=(0.436047,'amu*angstrom^2'), symmetry=1, barrier=(40.7139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.671228,'amu*angstrom^2'), symmetry=1, barrier=(62.69,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51144,0.0589901,-8.8389e-05,7.27965e-08,-2.40753e-11,62505,23.7227], Tmin=(100,'K'), Tmax=(790.373,'K')), NASAPolynomial(coeffs=[8.466,0.0207726,-1.01245e-05,1.94546e-09,-1.34937e-13,61500,-7.59598], Tmin=(790.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(518.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(C=CJO)"""),
)

species(
    label = '[O][CH][CH]C[C]=O(2757)',
    structure = SMILES('[O][CH][CH]C[C]=O'),
    E0 = (378.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,1855,455,950,180,1965.54,1966.61],'cm^-1')),
        HinderedRotor(inertia=(0.204619,'amu*angstrom^2'), symmetry=1, barrier=(4.70459,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205405,'amu*angstrom^2'), symmetry=1, barrier=(4.72267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204055,'amu*angstrom^2'), symmetry=1, barrier=(4.69162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66112,0.0621186,-0.000114462,1.12651e-07,-4.13139e-11,45649.3,26.4585], Tmin=(100,'K'), Tmax=(870.975,'K')), NASAPolynomial(coeffs=[3.28609,0.0293461,-1.44328e-05,2.72199e-09,-1.83884e-13,46326.3,24.3544], Tmin=(870.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCO) + radical(CCsJOH) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C[CH][CH][C]=O(10181)',
    structure = SMILES('[O][C]=C[CH]C[O]'),
    E0 = (338.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96711,0.0479086,-4.8249e-05,2.76754e-08,-6.76196e-12,40805.1,22.5476], Tmin=(100,'K'), Tmax=(962.957,'K')), NASAPolynomial(coeffs=[7.64513,0.0243233,-1.15106e-05,2.24136e-09,-1.58965e-13,39711.6,-4.6313], Tmin=(962.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.687,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(C=CCJCO) + radical(C=CJO)"""),
)

species(
    label = '[O][CH][CH][CH]C=O(10182)',
    structure = SMILES('[O][CH][CH]C=C[O]'),
    E0 = (279.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30309,0.0571973,-6.11819e-05,3.30141e-08,-7.02987e-12,33683.8,20.6488], Tmin=(100,'K'), Tmax=(1143.62,'K')), NASAPolynomial(coeffs=[13.1786,0.0156617,-6.70392e-06,1.25717e-09,-8.78156e-14,30967.5,-38.2375], Tmin=(1143.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(C=CCJCO) + radical(CCsJOH)"""),
)

species(
    label = 'O=[C][CH][CH][CH]O(10183)',
    structure = SMILES('[O][C]=C[CH][CH]O'),
    E0 = (293.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16839,0.0614964,-7.51252e-05,4.61072e-08,-1.10452e-11,35376.1,24.094], Tmin=(100,'K'), Tmax=(1026.67,'K')), NASAPolynomial(coeffs=[13.359,0.0139994,-5.72838e-06,1.04309e-09,-7.15028e-14,32873,-35.039], Tmin=(1026.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([O])[CH][C]=O(2671)',
    structure = SMILES('[CH2]C([O])[CH][C]=O'),
    E0 = (402.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,851.202,2647.18],'cm^-1')),
        HinderedRotor(inertia=(0.229533,'amu*angstrom^2'), symmetry=1, barrier=(5.27742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0274428,'amu*angstrom^2'), symmetry=1, barrier=(14.1083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.613626,'amu*angstrom^2'), symmetry=1, barrier=(14.1085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63814,0.0543197,-7.05236e-05,4.96507e-08,-1.40152e-11,48461.1,26.2344], Tmin=(100,'K'), Tmax=(865.41,'K')), NASAPolynomial(coeffs=[9.44214,0.0182491,-8.00344e-06,1.48877e-09,-1.02264e-13,47110.3,-10.2874], Tmin=(865.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCJCO) + radical(CJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1=CCC1[O](10184)',
    structure = SMILES('[O]C1C[CH]C1=O'),
    E0 = (120.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89758,0.00235898,9.66653e-05,-1.31212e-07,5.12136e-11,14498.5,20.7005], Tmin=(100,'K'), Tmax=(945.99,'K')), NASAPolynomial(coeffs=[14.6694,0.00841195,-1.45678e-06,3.22674e-10,-3.45183e-14,9773.23,-48.6414], Tmin=(945.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ) + radical(CCJC=O)"""),
)

species(
    label = '[O]C=CC=C[O](6971)',
    structure = SMILES('[O]C=CC=C[O]'),
    E0 = (-40.7388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42358,0.0345724,3.5106e-05,-9.2736e-08,4.46261e-11,-4786.29,19.09], Tmin=(100,'K'), Tmax=(909.421,'K')), NASAPolynomial(coeffs=[24.1225,-0.0068109,6.94716e-06,-1.41385e-09,9.17745e-14,-11332.1,-101.554], Tmin=(909.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-40.7388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ)"""),
)

species(
    label = '[O][C]=[C]CC=O(9881)',
    structure = SMILES('O=[C][C]CC=O'),
    E0 = (298.273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.184169,'amu*angstrom^2'), symmetry=1, barrier=(4.23441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184237,'amu*angstrom^2'), symmetry=1, barrier=(4.23596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.572776,'amu*angstrom^2'), symmetry=1, barrier=(13.1692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78933,0.0562281,-9.75861e-05,9.16398e-08,-3.28252e-11,35946.2,22.326], Tmin=(100,'K'), Tmax=(858.014,'K')), NASAPolynomial(coeffs=[5.13849,0.0238797,-1.17777e-05,2.2359e-09,-1.52217e-13,35987.5,10.2708], Tmin=(858.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(298.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][CH]C[C]=C[O](10185)',
    structure = SMILES('[O][CH]C[C]=C[O]'),
    E0 = (400.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,401.746,402.028,402.948,402.97],'cm^-1')),
        HinderedRotor(inertia=(0.0481801,'amu*angstrom^2'), symmetry=1, barrier=(5.54541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139359,'amu*angstrom^2'), symmetry=1, barrier=(16.0596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5425,0.0571376,-7.55146e-05,5.47417e-08,-1.60098e-11,48214.5,23.7625], Tmin=(100,'K'), Tmax=(833.909,'K')), NASAPolynomial(coeffs=[9.22672,0.0202792,-9.21582e-06,1.73972e-09,-1.20385e-13,46932.9,-11.9138], Tmin=(833.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[O][C]=[C]CC[O](10186)',
    structure = SMILES('[O]CC[C][C]=O'),
    E0 = (446.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1855,455,950,180,180,180,1510.1],'cm^-1')),
        HinderedRotor(inertia=(0.193147,'amu*angstrom^2'), symmetry=1, barrier=(4.44083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193026,'amu*angstrom^2'), symmetry=1, barrier=(4.43805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193058,'amu*angstrom^2'), symmetry=1, barrier=(4.43879,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66183,0.0619139,-0.000112394,1.10954e-07,-4.1037e-11,53829.8,24.3152], Tmin=(100,'K'), Tmax=(863.757,'K')), NASAPolynomial(coeffs=[3.10747,0.0305737,-1.51687e-05,2.88017e-09,-1.95757e-13,54499.4,22.8746], Tmin=(863.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][CH]C[C]=[C]O(10187)',
    structure = SMILES('[O][CH]C[C]=[C]O'),
    E0 = (498.447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1670,1700,300,440,363.082,363.092,363.134],'cm^-1')),
        HinderedRotor(inertia=(0.0863502,'amu*angstrom^2'), symmetry=1, barrier=(8.07906,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0863547,'amu*angstrom^2'), symmetry=1, barrier=(8.07918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0863502,'amu*angstrom^2'), symmetry=1, barrier=(8.07926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21173,0.0674879,-0.000111941,9.59357e-08,-3.14591e-11,60043.9,26.6614], Tmin=(100,'K'), Tmax=(886.464,'K')), NASAPolynomial(coeffs=[8.83377,0.0198527,-8.92956e-06,1.61437e-09,-1.06247e-13,59212.9,-6.25732], Tmin=(886.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[O][C]=[C]C[CH]O(10188)',
    structure = SMILES('O=[C][C]C[CH]O'),
    E0 = (401.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1855,455,950,230.814,230.824,230.825],'cm^-1')),
        HinderedRotor(inertia=(0.00316445,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172622,'amu*angstrom^2'), symmetry=1, barrier=(6.52588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172623,'amu*angstrom^2'), symmetry=1, barrier=(6.52596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400436,'amu*angstrom^2'), symmetry=1, barrier=(15.1387,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06613,0.0729545,-0.000129605,1.15673e-07,-3.88938e-11,48392.1,25.1432], Tmin=(100,'K'), Tmax=(890.534,'K')), NASAPolynomial(coeffs=[8.32345,0.0211218,-9.89984e-06,1.80501e-09,-1.18621e-13,47862.2,-4.74539], Tmin=(890.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCsJOH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][CH][CH]C=[C]O(10189)',
    structure = SMILES('[O]C=C[CH][C]O'),
    E0 = (334.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.31894,'amu*angstrom^2'), symmetry=1, barrier=(30.3251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32079,'amu*angstrom^2'), symmetry=1, barrier=(30.3676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.32121,'amu*angstrom^2'), symmetry=1, barrier=(30.3773,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.27501,0.065554,-7.17726e-05,3.59548e-08,-6.66266e-12,40346.2,23.3617], Tmin=(100,'K'), Tmax=(1510.38,'K')), NASAPolynomial(coeffs=[20.667,0.00203578,7.56998e-07,-2.29302e-10,1.68241e-14,35271.4,-79.8346], Tmin=(1510.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJCO) + radical(CH2_triplet)"""),
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
    E0 = (346.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (346.562,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (848.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (821.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (839.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (839.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (349.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (369.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (369.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (381.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (514.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (431.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (730.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (537.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (504.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (499.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (457.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (559.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (354.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (409.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (525.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (555.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (604.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (588.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (648.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (549.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (418.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['HCCO(2227)', 'vinoxy(1351)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['O=[C][CH]CC=O(2652)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH][O](1556)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH][O](1548)', '[CH2][CH][C]=O(9632)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[O][CH]C[C][C]=O(10176)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[O][C]C[CH][C]=O(10177)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['[O]C1CC1[C]=O(10178)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['O=[C]C[CH]C=O(2651)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['[O]CC=C[C]=O(10179)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[O][C]=CC=C[O](9609)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH][O](1548)', 'CH2CHCO(3668)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C][O](6861)', 'vinoxy(1351)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[O][CH][CH][CH][C]=O(10180)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][CH][CH]C[C]=O(2757)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['[O]C[CH][CH][C]=O(10181)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['[O][CH][CH][CH]C=O(10182)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.21412e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['O=[C][CH][CH][CH]O(10183)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.47378e+09,'s^-1'), n=1.09705, Ea=(110.443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;XH_out] for rate rule [R3HJ;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([O])[CH][C]=O(2671)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['[O]C1=CCC1[O](10184)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O][CH]C[CH][C]=O(2758)'],
    products = ['[O]C=CC=C[O](6971)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[O][C]=[C]CC=O(9881)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HCCO(2227)', '[CH2][CH][O](1556)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O][CH]C[C]=C[O](10185)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O][C]=[C]CC[O](10186)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][CH]C[C]=[C]O(10187)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O][C]=[C]C[CH]O(10188)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1194.19,'s^-1'), n=2.56519, Ea=(148.442,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4Hall;Y_rad_out;O_H_out] + [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][CH][CH]C=[C]O(10189)'],
    products = ['[O][CH]C[CH][C]=O(2758)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2394',
    isomers = [
        '[O][CH]C[CH][C]=O(2758)',
    ],
    reactants = [
        ('HCCO(2227)', 'vinoxy(1351)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2394',
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


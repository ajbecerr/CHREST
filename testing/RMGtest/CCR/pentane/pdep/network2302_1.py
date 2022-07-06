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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79021,0.0664689,-5.56502e-05,2.64471e-08,-5.20762e-12,30453.4,29.7261], Tmin=(100,'K'), Tmax=(1206.53,'K')), NASAPolynomial(coeffs=[11.5566,0.0307744,-1.12727e-05,1.92592e-09,-1.26606e-13,27855.4,-24.2368], Tmin=(1206.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = 'butene1(35)',
    structure = SMILES('C=CCC'),
    E0 = (-16.4325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,247.029],'cm^-1')),
        HinderedRotor(inertia=(0.178818,'amu*angstrom^2'), symmetry=1, barrier=(7.72187,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177126,'amu*angstrom^2'), symmetry=1, barrier=(7.72778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.176,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58774,0.0232777,1.93416e-05,-3.55502e-08,1.36908e-11,-1918.73,14.575], Tmin=(100,'K'), Tmax=(1007.28,'K')), NASAPolynomial(coeffs=[7.20513,0.0236362,-9.03154e-06,1.65394e-09,-1.1602e-13,-3797.32,-12.4424], Tmin=(1007.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.4325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene1""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][CH]C([C]=O)CC(9730)',
    structure = SMILES('[CH2][CH]C([C]=O)CC'),
    E0 = (283.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,1496.8,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00533532,'amu*angstrom^2'), symmetry=1, barrier=(8.48228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0781593,'amu*angstrom^2'), symmetry=1, barrier=(8.48228,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227741,'amu*angstrom^2'), symmetry=1, barrier=(24.7179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227764,'amu*angstrom^2'), symmetry=1, barrier=(24.7179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0781606,'amu*angstrom^2'), symmetry=1, barrier=(8.48227,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0976,0.0612738,-4.35914e-05,1.63638e-08,-2.56755e-12,34192.9,31.8431], Tmin=(100,'K'), Tmax=(1451.41,'K')), NASAPolynomial(coeffs=[11.7661,0.0318723,-1.32061e-05,2.40735e-09,-1.63633e-13,31096,-23.601], Tmin=(1451.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'CCC[CH][CH][C]=O(9731)',
    structure = SMILES('CCC[CH]C=[C][O]'),
    E0 = (249.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.933137,0.0575857,-2.42316e-05,-6.63272e-09,5.79628e-12,30120.1,28.8692], Tmin=(100,'K'), Tmax=(1064.96,'K')), NASAPolynomial(coeffs=[13.2192,0.0301624,-1.19775e-05,2.2049e-09,-1.53777e-13,26441.5,-36.1621], Tmin=(1064.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(C=CJO)"""),
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
    label = '[CH2][CH]CC(37)',
    structure = SMILES('[CH2][CH]CC'),
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
    label = 'CH2(T)(20)',
    structure = SMILES('[CH2]'),
    E0 = (382.412,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1035.69,2892.37,3631.64],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02765,-0.00110501,7.85495e-06,-9.33167e-09,3.86549e-12,45993,0.552247], Tmin=(100,'K'), Tmax=(735.837,'K')), NASAPolynomial(coeffs=[3.71316,0.00171632,-1.62522e-07,-1.4674e-11,2.48393e-15,46009.2,1.76854], Tmin=(735.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.412,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C([C][C]=O)CC(9732)',
    structure = SMILES('[CH2]C([C][C]=O)CC'),
    E0 = (532.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788995,0.0732268,-8.28009e-05,5.51873e-08,-1.52288e-11,64206.3,29.8836], Tmin=(100,'K'), Tmax=(875.177,'K')), NASAPolynomial(coeffs=[9.62657,0.0328358,-1.35754e-05,2.45637e-09,-1.6627e-13,62659.4,-11.5746], Tmin=(875.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C([CH][C]=O)CC(9733)',
    structure = SMILES('[CH]C([CH][C]=O)CC'),
    E0 = (495.347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82903,0.0665415,-5.80793e-05,2.76138e-08,-5.38313e-12,59693.3,28.81], Tmin=(100,'K'), Tmax=(1220.09,'K')), NASAPolynomial(coeffs=[12.4565,0.0284208,-1.12121e-05,2.00466e-09,-1.35629e-13,56856,-29.5983], Tmin=(1220.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC1CC1[C]=O(9734)',
    structure = SMILES('CCC1CC1[C]=O'),
    E0 = (31.8191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37253,0.0435983,1.70645e-05,-5.11379e-08,2.20112e-11,3934.14,25.1553], Tmin=(100,'K'), Tmax=(988.716,'K')), NASAPolynomial(coeffs=[13.5704,0.0279964,-1.04631e-05,1.94458e-09,-1.39742e-13,-127.369,-41.8955], Tmin=(988.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.8191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=C(CC)C[C]=O(2464)',
    structure = SMILES('C=C(CC)C[C]=O'),
    E0 = (5.12134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,403.965,405.398],'cm^-1')),
        HinderedRotor(inertia=(0.00103039,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116029,'amu*angstrom^2'), symmetry=1, barrier=(13.4637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11592,'amu*angstrom^2'), symmetry=1, barrier=(13.4662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11581,'amu*angstrom^2'), symmetry=1, barrier=(13.4677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.958161,0.0564238,-2.27169e-05,-5.86096e-09,4.51892e-12,734.097,27.0856], Tmin=(100,'K'), Tmax=(1176.87,'K')), NASAPolynomial(coeffs=[14.2068,0.0303204,-1.35698e-05,2.62263e-09,-1.86094e-13,-3694.99,-44.5577], Tmin=(1176.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.12134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
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
    label = '[CH2]C(=C[C]=O)CC(9736)',
    structure = SMILES('[CH2]C(=C[C]=O)CC'),
    E0 = (133.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.713916,'amu*angstrom^2'), symmetry=1, barrier=(16.4143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712993,'amu*angstrom^2'), symmetry=1, barrier=(16.3931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712789,'amu*angstrom^2'), symmetry=1, barrier=(16.3884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712997,'amu*angstrom^2'), symmetry=1, barrier=(16.3932,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99177,0.0699981,-6.37894e-05,3.13321e-08,-6.4891e-12,16108,24.5225], Tmin=(100,'K'), Tmax=(1122.41,'K')), NASAPolynomial(coeffs=[10.9256,0.0345967,-1.64789e-05,3.23176e-09,-2.30219e-13,13878,-24.5495], Tmin=(1122.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = 'CCC=C[C]=O(9684)',
    structure = SMILES('CCC=C[C]=O'),
    E0 = (20.6116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.66804,'amu*angstrom^2'), symmetry=1, barrier=(15.3596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.66818,'amu*angstrom^2'), symmetry=1, barrier=(15.3628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.667854,'amu*angstrom^2'), symmetry=1, barrier=(15.3553,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83641,0.0525271,-4.2117e-05,1.82991e-08,-3.48325e-12,2552.66,20.4505], Tmin=(100,'K'), Tmax=(1168.55,'K')), NASAPolynomial(coeffs=[7.92868,0.0316729,-1.53478e-05,3.02698e-09,-2.15933e-13,1128.84,-9.8901], Tmin=(1168.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.6116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'C2H5(14)',
    structure = SMILES('C[CH2]'),
    E0 = (109.359,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,180,1148.63,1148.78,3699.3,3699.53],'cm^-1')),
        HinderedRotor(inertia=(0.00634458,'amu*angstrom^2'), symmetry=1, barrier=(5.93959,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.0611,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2097.75,'J/mol'), sigma=(4.302,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.69276,0.00187276,3.1195e-05,-3.71212e-08,1.32027e-11,13168.3,7.07155], Tmin=(100,'K'), Tmax=(967.563,'K')), NASAPolynomial(coeffs=[4.21409,0.0122658,-4.37068e-06,7.87987e-10,-5.5604e-14,12480.1,1.53827], Tmin=(967.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.359,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""C2H5""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=CC=[C][O](9589)',
    structure = SMILES('[CH2]C=C[C]=O'),
    E0 = (194.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.05258,'amu*angstrom^2'), symmetry=1, barrier=(24.201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05258,'amu*angstrom^2'), symmetry=1, barrier=(24.201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33777,0.0398735,-3.70418e-05,1.81131e-08,-3.74289e-12,23480.8,16.0206], Tmin=(100,'K'), Tmax=(1118.07,'K')), NASAPolynomial(coeffs=[7.9573,0.0197688,-1.00692e-05,2.03004e-09,-1.46689e-13,22224.2,-11.7174], Tmin=(1118.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(194.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2][CH][CH][C]=O(9637)',
    structure = SMILES('[CH2][CH]C=[C][O]'),
    E0 = (502.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,458.276,458.348],'cm^-1')),
        HinderedRotor(inertia=(0.000816836,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000816843,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47043,0.0263284,1.39799e-06,-2.1241e-08,1.01461e-11,60468.5,21.3394], Tmin=(100,'K'), Tmax=(980.457,'K')), NASAPolynomial(coeffs=[9.92325,0.0131468,-4.7862e-06,8.81255e-10,-6.33441e-14,58179.2,-18.6908], Tmin=(980.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]([CH][C]=O)CC(9737)',
    structure = SMILES('[CH2][C](C=[C][O])CC'),
    E0 = (437.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,180,759.973,761.477],'cm^-1')),
        HinderedRotor(inertia=(0.0959315,'amu*angstrom^2'), symmetry=1, barrier=(2.20565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0953191,'amu*angstrom^2'), symmetry=1, barrier=(2.19157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00520167,'amu*angstrom^2'), symmetry=1, barrier=(2.19709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87885,'amu*angstrom^2'), symmetry=1, barrier=(66.1904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12649,0.0553595,-2.558e-05,-6.64603e-09,7.12075e-12,52720.3,29.7616], Tmin=(100,'K'), Tmax=(960.613,'K')), NASAPolynomial(coeffs=[12.2939,0.0268901,-9.28177e-06,1.5839e-09,-1.06634e-13,49742.8,-27.9966], Tmin=(960.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(Isobutyl) + radical(C=CJO)"""),
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
    label = '[CH2]C([CH2])[CH][C]=O(2582)',
    structure = SMILES('[CH2]C([CH2])[CH][C]=O'),
    E0 = (477.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,318.813],'cm^-1')),
        HinderedRotor(inertia=(0.0281889,'amu*angstrom^2'), symmetry=1, barrier=(69.0787,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00379562,'amu*angstrom^2'), symmetry=1, barrier=(9.32132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281264,'amu*angstrom^2'), symmetry=1, barrier=(69.0081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133159,'amu*angstrom^2'), symmetry=1, barrier=(9.3129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41759,0.0538488,-5.61778e-05,3.36165e-08,-7.9844e-12,57553.1,26.3147], Tmin=(100,'K'), Tmax=(1127.64,'K')), NASAPolynomial(coeffs=[9.9417,0.0197888,-5.78525e-06,8.17611e-10,-4.62739e-14,55873.7,-14.7554], Tmin=(1127.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2]CC([CH2])[CH][C]=O(5300)',
    structure = SMILES('[CH2]CC([CH2])[CH][C]=O'),
    E0 = (457.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,338.159,338.16],'cm^-1')),
        HinderedRotor(inertia=(0.00147421,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00147422,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00147419,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0014742,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00147422,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922897,0.0660645,-6.20976e-05,3.35322e-08,-7.47196e-12,55132,31.0378], Tmin=(100,'K'), Tmax=(1076.22,'K')), NASAPolynomial(coeffs=[10.7763,0.0294415,-1.10529e-05,1.91188e-09,-1.26599e-13,53011.1,-17.2228], Tmin=(1076.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C](CC)C[C]=O(2466)',
    structure = SMILES('[CH2][C](CC)C[C]=O'),
    E0 = (270.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3100,440,815,1455,1000,1855,455,950,180,3046.75],'cm^-1')),
        HinderedRotor(inertia=(0.151978,'amu*angstrom^2'), symmetry=1, barrier=(3.49427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152017,'amu*angstrom^2'), symmetry=1, barrier=(3.49516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151921,'amu*angstrom^2'), symmetry=1, barrier=(3.49297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151989,'amu*angstrom^2'), symmetry=1, barrier=(3.49452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151961,'amu*angstrom^2'), symmetry=1, barrier=(3.49389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.84503,0.0786801,-0.000114575,1.04244e-07,-3.72684e-11,32591.1,30.9046], Tmin=(100,'K'), Tmax=(866.786,'K')), NASAPolynomial(coeffs=[3.66509,0.0448084,-1.98638e-05,3.63798e-09,-2.44169e-13,32885.8,22.2224], Tmin=(866.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2][C]([CH]C=O)CC(9739)',
    structure = SMILES('[CH2][C](C=C[O])CC'),
    E0 = (197.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932166,0.0534101,2.01471e-06,-4.64217e-08,2.37989e-11,23898.7,27.3311], Tmin=(100,'K'), Tmax=(930.339,'K')), NASAPolynomial(coeffs=[16.3548,0.0228998,-6.51328e-06,1.05075e-09,-7.2452e-14,19479.7,-54.2873], Tmin=(930.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(Isobutyl)"""),
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
    label = '[CH2]CC([CH2])C[C]=O(2469)',
    structure = SMILES('[CH2]CC([CH2])C[C]=O'),
    E0 = (289.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,229.445,1259.08],'cm^-1')),
        HinderedRotor(inertia=(0.00320221,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114217,'amu*angstrom^2'), symmetry=1, barrier=(4.26686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114217,'amu*angstrom^2'), symmetry=1, barrier=(4.26687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114215,'amu*angstrom^2'), symmetry=1, barrier=(4.26686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114217,'amu*angstrom^2'), symmetry=1, barrier=(4.26685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925285,0.0703288,-7.16022e-05,4.48254e-08,-1.19616e-11,34979,31.2773], Tmin=(100,'K'), Tmax=(892.658,'K')), NASAPolynomial(coeffs=[8.39977,0.0368361,-1.53227e-05,2.79457e-09,-1.90489e-13,33644.5,-3.93412], Tmin=(892.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
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
    label = '[CH2]CC([CH2])[CH]C=O(5535)',
    structure = SMILES('[CH2]CC([CH2])C=C[O]'),
    E0 = (270.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,489.462,489.463,489.47],'cm^-1')),
        HinderedRotor(inertia=(0.000703655,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.080465,'amu*angstrom^2'), symmetry=1, barrier=(13.6792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0804634,'amu*angstrom^2'), symmetry=1, barrier=(13.6792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.703656,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690581,0.0595632,-1.40415e-05,-3.1106e-08,1.85148e-11,32718.4,31.0378], Tmin=(100,'K'), Tmax=(938.115,'K')), NASAPolynomial(coeffs=[17.3253,0.0217135,-6.41292e-06,1.05961e-09,-7.36308e-14,28141.7,-55.9105], Tmin=(938.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C=C=O)CC(2465)',
    structure = SMILES('[CH2]C(C=C=O)CC'),
    E0 = (52.4284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2120,512.5,787.5,262.196,262.199],'cm^-1')),
        HinderedRotor(inertia=(0.222037,'amu*angstrom^2'), symmetry=1, barrier=(10.8319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00245211,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222039,'amu*angstrom^2'), symmetry=1, barrier=(10.8319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222032,'amu*angstrom^2'), symmetry=1, barrier=(10.8319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.910917,0.0706994,-6.56844e-05,3.52187e-08,-8.04121e-12,6414.62,26.3981], Tmin=(100,'K'), Tmax=(1027.93,'K')), NASAPolynomial(coeffs=[9.65419,0.0366769,-1.60377e-05,3.02056e-09,-2.10457e-13,4617.11,-16.0241], Tmin=(1027.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.4284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O)"""),
)

species(
    label = 'CC[CH]C[CH][C]=O(2531)',
    structure = SMILES('CC[CH]C[CH][C]=O'),
    E0 = (247.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,318.557,951.52,1897.58],'cm^-1')),
        HinderedRotor(inertia=(0.0721619,'amu*angstrom^2'), symmetry=1, barrier=(3.22971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721619,'amu*angstrom^2'), symmetry=1, barrier=(3.22971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721619,'amu*angstrom^2'), symmetry=1, barrier=(3.22971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721619,'amu*angstrom^2'), symmetry=1, barrier=(3.22971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721619,'amu*angstrom^2'), symmetry=1, barrier=(3.22971,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97333,0.0532059,-4.21998e-06,-8.51042e-08,8.31634e-11,29805.5,26.7174], Tmin=(100,'K'), Tmax=(469.408,'K')), NASAPolynomial(coeffs=[4.95805,0.042757,-1.87154e-05,3.49004e-09,-2.40422e-13,29360.2,12.8166], Tmin=(469.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC1C=C([O])C1(9740)',
    structure = SMILES('CCC1[CH]C(=O)C1'),
    E0 = (-10.4027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75774,0.0302003,6.02813e-05,-9.5727e-08,3.73775e-11,-1153.36,23.4076], Tmin=(100,'K'), Tmax=(980.327,'K')), NASAPolynomial(coeffs=[13.5886,0.029356,-1.09976e-05,2.09729e-09,-1.54764e-13,-5752.04,-45.0583], Tmin=(980.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.4027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O)"""),
)

species(
    label = 'C=C(C=C[O])CC(9741)',
    structure = SMILES('C=C(C=C[O])CC'),
    E0 = (-33.6449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.557549,0.0576784,2.95888e-06,-5.55944e-08,2.85158e-11,-3906.08,25.1333], Tmin=(100,'K'), Tmax=(940.5,'K')), NASAPolynomial(coeffs=[20.7726,0.0169643,-4.29411e-06,7.1665e-10,-5.42759e-14,-9710.32,-81.7947], Tmin=(940.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.6449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C#C[O])CC(9742)',
    structure = SMILES('[CH2]C([C]=C=O)CC'),
    E0 = (258.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,196.899,196.913],'cm^-1')),
        HinderedRotor(inertia=(0.33143,'amu*angstrom^2'), symmetry=1, barrier=(9.11886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.331451,'amu*angstrom^2'), symmetry=1, barrier=(9.11883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00434867,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0681,'amu*angstrom^2'), symmetry=1, barrier=(29.388,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.772416,0.0756731,-8.9207e-05,6.31827e-08,-1.88213e-11,31168.7,26.7285], Tmin=(100,'K'), Tmax=(807.034,'K')), NASAPolynomial(coeffs=[8.66698,0.0365475,-1.64917e-05,3.11958e-09,-2.16683e-13,29894.4,-9.66634], Tmin=(807.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O) + radical(CC(C)CJ=C=O)"""),
)

species(
    label = '[CH2]C([C]=C[O])CC(9743)',
    structure = SMILES('[CH2]C([C]=C[O])CC'),
    E0 = (303.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,402.109,402.317,402.585],'cm^-1')),
        HinderedRotor(inertia=(0.00104255,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128396,'amu*angstrom^2'), symmetry=1, barrier=(14.7563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128654,'amu*angstrom^2'), symmetry=1, barrier=(14.7562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171023,'amu*angstrom^2'), symmetry=1, barrier=(19.6272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.617124,0.0634729,-2.85612e-05,-1.37927e-08,1.18273e-11,36639.2,29.9883], Tmin=(100,'K'), Tmax=(944.223,'K')), NASAPolynomial(coeffs=[16.3637,0.0234709,-7.43753e-06,1.24614e-09,-8.51449e-14,32475.1,-51.3798], Tmin=(944.223,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Isobutyl) + radical(Cds_S)"""),
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
    label = '[CH2]C([C]=[C]O)CC(9745)',
    structure = SMILES('[CH2]C([C]=[C]O)CC'),
    E0 = (401.821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1670,1700,300,440,373.424,373.428],'cm^-1')),
        HinderedRotor(inertia=(0.00120886,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113003,'amu*angstrom^2'), symmetry=1, barrier=(11.1822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113005,'amu*angstrom^2'), symmetry=1, barrier=(11.1822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169378,'amu*angstrom^2'), symmetry=1, barrier=(16.7613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112999,'amu*angstrom^2'), symmetry=1, barrier=(11.1821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0693357,0.0762695,-7.2998e-05,3.69915e-08,-7.35743e-12,48478.1,33.673], Tmin=(100,'K'), Tmax=(1290.9,'K')), NASAPolynomial(coeffs=[16.778,0.0217356,-6.42347e-06,9.53615e-10,-5.7441e-14,44394.3,-50.3123], Tmin=(1290.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Isobutyl) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](C=[C]O)CC(9746)',
    structure = SMILES('[CH2][C](C=[C]O)CC'),
    E0 = (295.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,377.749,382.315],'cm^-1')),
        HinderedRotor(inertia=(0.0010801,'amu*angstrom^2'), symmetry=1, barrier=(0.119834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114773,'amu*angstrom^2'), symmetry=1, barrier=(13.1898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00666869,'amu*angstrom^2'), symmetry=1, barrier=(13.3073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132245,'amu*angstrom^2'), symmetry=1, barrier=(13.3161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.703974,'amu*angstrom^2'), symmetry=1, barrier=(72.2279,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738112,0.0620367,-2.78535e-05,-1.45018e-08,1.26533e-11,35722.3,29.7469], Tmin=(100,'K'), Tmax=(914.046,'K')), NASAPolynomial(coeffs=[15.5093,0.0232708,-6.69876e-06,1.03892e-09,-6.78626e-14,31941.1,-46.1006], Tmin=(914.046,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_T) + radical(Isobutyl) + radical(C=CJO)"""),
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
    label = '[CH2]CC([CH2])C=[C]O(9748)',
    structure = SMILES('[CH2]CC([CH2])C=[C]O'),
    E0 = (369.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,266.132,266.987],'cm^-1')),
        HinderedRotor(inertia=(0.00236745,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00237577,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00236553,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279213,'amu*angstrom^2'), symmetry=1, barrier=(14.0651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.516393,'amu*angstrom^2'), symmetry=1, barrier=(26.0867,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.128828,0.0754535,-6.88374e-05,3.25471e-08,-5.94914e-12,44569.4,35.7049], Tmin=(100,'K'), Tmax=(1480.17,'K')), NASAPolynomial(coeffs=[18.1221,0.0192926,-4.9929e-06,6.6981e-10,-3.77926e-14,39915.7,-56.9718], Tmin=(1480.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Isobutyl) + radical(RCCJ) + radical(C=CJO)"""),
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
    E0 = (252.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (691.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (446.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (409.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (736.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (711.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (744.692,'kJ/mol'),
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
    E0 = (257.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (275.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (275.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (352.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (441.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (408.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (328.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (611.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (649.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (614.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (658.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (669.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (428.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (393.702,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (423.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (399.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (405.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (378.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (365.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (335.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (344.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (252.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (412.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (260.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (315.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (485.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (443.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (507.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (479.697,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (552.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (379.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (385.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (442.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['HCCO(2227)', 'butene1(35)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2]C(C)[CH][C]=O(2395)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]C([C]=O)CC(9730)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['CCC[CH][CH][C]=O(9731)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH]CC(37)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', 'CC[CH][CH][C]=O(9665)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2]C([C][C]=O)CC(9732)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C([CH][C]=O)CC(9733)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['CCC1CC1[C]=O(9734)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['C=C(CC)C[C]=O(2464)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['CCC(C)=C[C]=O(9735)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH2]C(=C[C]=O)CC(9736)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2(T)(20)', 'CCC=C[C]=O(9684)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C][O](6861)', 'butene1(35)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['C2H5(14)', 'C=CC=[C][O](9589)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.0941959,'m^3/(mol*s)'), n=1.97652, Ea=(24.8542,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C2H5(14)', '[CH2][CH][CH][C]=O(9637)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.30526e+07,'m^3/(mol*s)'), n=-0.283333, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-6.50053689359e-11, var=0.305422193575, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C
    Total Standard Deviation in ln(k): 1.10791715097
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2][C]([CH][C]=O)CC(9737)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(17)', '[CH2]C([CH2])[CH][C]=O(2582)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.26105e+08,'m^3/(mol*s)'), n=-0.283333, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-6.50053689359e-11, var=0.305422193575, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C
    Total Standard Deviation in ln(k): 1.10791715097
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2]C([CH]C)[CH][C]=O(4240)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2]CC([CH2])[CH][C]=O(5300)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][C](CC)C[C]=O(2466)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['CC[C](C)[CH][C]=O(9738)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.17661e+06,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['[CH2][C]([CH]C=O)CC(9739)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.10706e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH]C)[CH]C=O(4241)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(30253,'s^-1'), n=2.05523, Ea=(118.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R4HJ_2;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]CC([CH2])C[C]=O(2469)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(927.918,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]CC(C)[CH][C]=O(4182)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]CC([CH2])[CH]C=O(5535)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['[CH2]C(C=C=O)CC(2465)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['CC[CH]C[CH][C]=O(2531)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['CCC1C=C([O])C1(9740)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([CH][C]=O)CC(2468)'],
    products = ['C=C(C=C[O])CC(9741)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', '[CH2]C(C#C[O])CC(9742)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['HCCO(2227)', '[CH2][CH]CC(37)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([C]=C[O])CC(9743)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['CCC(C)[C]=[C][O](9744)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([C]=[C]O)CC(9745)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][C](C=[C]O)CC(9746)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([CH]C)C=[C]O(9747)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(0.0378492,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R5HJ_3;C_rad_out_H/NonDeC;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]CC([CH2])C=[C]O(9748)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_4;C_rad_out_2H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2302',
    isomers = [
        '[CH2]C([CH][C]=O)CC(2468)',
    ],
    reactants = [
        ('HCCO(2227)', 'butene1(35)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2302',
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


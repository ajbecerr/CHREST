species(
    label = '[CH2]C=CC[CH][C]=O(5329)',
    structure = SMILES('[CH2]C=CC[CH][C]=O'),
    E0 = (319.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.778956,'amu*angstrom^2'), symmetry=1, barrier=(17.9097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403324,'amu*angstrom^2'), symmetry=1, barrier=(50.6156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779054,'amu*angstrom^2'), symmetry=1, barrier=(17.912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968321,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907921,0.0595326,-4.70532e-05,1.94918e-08,-3.25348e-12,38588.9,28.1278], Tmin=(100,'K'), Tmax=(1429.35,'K')), NASAPolynomial(coeffs=[14.097,0.0226237,-8.32024e-06,1.42646e-09,-9.37883e-14,34818.5,-40.2132], Tmin=(1429.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Allyl_P) + radical(CCCJ=O)"""),
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
    label = 'butadiene13(1350)',
    structure = SMILES('C=CC=C'),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.30711,'amu*angstrom^2'), symmetry=1, barrier=(30.0531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80606,0.0102576,6.17291e-05,-9.01684e-08,3.59136e-11,11658.5,12.0619], Tmin=(100,'K'), Tmax=(946.037,'K')), NASAPolynomial(coeffs=[12.4692,0.0100558,-2.41229e-06,4.5713e-10,-3.93205e-14,8010.87,-43.6362], Tmin=(946.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]C=C(3743)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.210055,'amu*angstrom^2'), symmetry=1, barrier=(25.2323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779472,'amu*angstrom^2'), symmetry=1, barrier=(93.4341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56316,0.022343,1.87062e-05,-3.93092e-08,1.63979e-11,33100.5,13.4098], Tmin=(100,'K'), Tmax=(974.267,'K')), NASAPolynomial(coeffs=[9.83,0.0151965,-5.22268e-06,9.67646e-10,-7.07852e-14,30607.7,-26.9852], Tmin=(974.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC([CH2])[C]=O(10002)',
    structure = SMILES('[CH2]C=CC([CH2])[C]=O'),
    E0 = (354.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,274.38],'cm^-1')),
        HinderedRotor(inertia=(0.185804,'amu*angstrom^2'), symmetry=1, barrier=(9.92957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185813,'amu*angstrom^2'), symmetry=1, barrier=(9.92918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00223893,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.585306,'amu*angstrom^2'), symmetry=1, barrier=(31.278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12294,0.0662939,-6.68565e-05,3.81396e-08,-9.12615e-12,42758.8,27.4808], Tmin=(100,'K'), Tmax=(991.425,'K')), NASAPolynomial(coeffs=[9.75709,0.0314585,-1.41515e-05,2.69894e-09,-1.89347e-13,41046.8,-14.0996], Tmin=(991.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(Allyl_P) + radical(CC(C)CJ=O)"""),
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
    label = '[CH]=CC[CH][C]=O(10003)',
    structure = SMILES('[CH]=CC[CH][C]=O'),
    E0 = (451.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,3120,650,792.5,1650,464.284],'cm^-1')),
        HinderedRotor(inertia=(0.00228971,'amu*angstrom^2'), symmetry=1, barrier=(19.1528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.576136,'amu*angstrom^2'), symmetry=1, barrier=(13.2465,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.62747,'amu*angstrom^2'), symmetry=1, barrier=(83.4027,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73837,0.0470818,-4.22717e-05,2.05779e-08,-4.07462e-12,54384.7,23.4548], Tmin=(100,'K'), Tmax=(1209.25,'K')), NASAPolynomial(coeffs=[10.2657,0.0188748,-7.28266e-06,1.28827e-09,-8.66922e-14,52322.4,-19.3048], Tmin=(1209.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Cds_P) + radical(CCCJ=O)"""),
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
    label = '[CH2]C=CC[C][C]=O(10004)',
    structure = SMILES('[CH2]C=CC[C][C]=O'),
    E0 = (600.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950,201.083,893.444,893.971],'cm^-1')),
        HinderedRotor(inertia=(0.558185,'amu*angstrom^2'), symmetry=1, barrier=(12.8338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558438,'amu*angstrom^2'), symmetry=1, barrier=(12.8396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558281,'amu*angstrom^2'), symmetry=1, barrier=(12.836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.91041,'amu*angstrom^2'), symmetry=1, barrier=(89.908,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21223,0.0628526,-6.29322e-05,3.46186e-08,-7.86405e-12,72328.4,27.1781], Tmin=(100,'K'), Tmax=(1049.51,'K')), NASAPolynomial(coeffs=[10.5612,0.0272212,-1.20072e-05,2.27064e-09,-1.58656e-13,70366,-18.3773], Tmin=(1049.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C=CC[CH][C]=O(10005)',
    structure = SMILES('[CH]C=CC[CH][C]=O'),
    E0 = (539.047,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22475,0.0583617,-4.32797e-05,1.75529e-08,-3.00774e-12,64934.4,27.809], Tmin=(100,'K'), Tmax=(1340.45,'K')), NASAPolynomial(coeffs=[10.4507,0.0308308,-1.24719e-05,2.23074e-09,-1.5009e-13,62461,-19.4039], Tmin=(1340.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(539.047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'O=[C]C1CC=CC1(10006)',
    structure = SMILES('O=[C]C1CC=CC1'),
    E0 = (68.3456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03543,0.0271839,4.70065e-05,-7.74264e-08,3.04754e-11,8305.1,22.3191], Tmin=(100,'K'), Tmax=(982.238,'K')), NASAPolynomial(coeffs=[12.4669,0.0243262,-9.13853e-06,1.74932e-09,-1.29373e-13,4344.47,-37.5495], Tmin=(982.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.3456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'CC=CC=C[C]=O(10007)',
    structure = SMILES('CC=C[CH]C=C=O'),
    E0 = (67.3711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02962,0.0580588,-4.29012e-05,1.62743e-08,-2.50267e-12,8215.76,23.0833], Tmin=(100,'K'), Tmax=(1525.46,'K')), NASAPolynomial(coeffs=[13.9687,0.0241309,-9.53981e-06,1.69461e-09,-1.13308e-13,4268.13,-44.8043], Tmin=(1525.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.3711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(C=CCJC=C=O)"""),
)

species(
    label = 'C=C=CCC[C]=O(5327)',
    structure = SMILES('C=C=CCC[C]=O'),
    E0 = (177.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950,191.718,192.02],'cm^-1')),
        HinderedRotor(inertia=(0.385958,'amu*angstrom^2'), symmetry=1, barrier=(10.518,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.391139,'amu*angstrom^2'), symmetry=1, barrier=(10.4994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.39649,'amu*angstrom^2'), symmetry=1, barrier=(10.5014,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3804,0.0609701,-5.7359e-05,3.1645e-08,-7.52253e-12,21432.4,25.447], Tmin=(100,'K'), Tmax=(983.611,'K')), NASAPolynomial(coeffs=[8.23605,0.0330905,-1.48426e-05,2.82834e-09,-1.98312e-13,20083.7,-7.51422], Tmin=(983.611,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C[CH][CH][C]=O(10008)',
    structure = SMILES('[CH2][CH]C[CH]C=[C][O]'),
    E0 = (649.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,320.938,1161.84,1161.86,1161.87],'cm^-1')),
        HinderedRotor(inertia=(0.00535115,'amu*angstrom^2'), symmetry=1, barrier=(48.6604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28324,'amu*angstrom^2'), symmetry=1, barrier=(93.7698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00163705,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0699164,'amu*angstrom^2'), symmetry=1, barrier=(5.10982,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52939,0.054353,-3.98316e-05,1.59854e-08,-2.74662e-12,78161.9,31.246], Tmin=(100,'K'), Tmax=(1315.24,'K')), NASAPolynomial(coeffs=[9.25591,0.0308549,-1.30328e-05,2.40188e-09,-1.64692e-13,76129.4,-8.14725], Tmin=(1315.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJC) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = 'O=[C][CH]CC1[CH]C1(10009)',
    structure = SMILES('O=[C][CH]CC1[CH]C1'),
    E0 = (432.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57865,0.0439913,-6.43065e-06,-2.0472e-08,1.03444e-11,52088.2,27.0165], Tmin=(100,'K'), Tmax=(1009.77,'K')), NASAPolynomial(coeffs=[11.3167,0.0258391,-9.80373e-06,1.78445e-09,-1.24787e-13,49080.4,-25.2142], Tmin=(1009.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(cyclopropane) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1[CH]CC1[C]=O(9973)',
    structure = SMILES('[CH2]C1[CH]CC1[C]=O'),
    E0 = (420.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85077,0.0367537,1.2942e-05,-3.86013e-08,1.62746e-11,50639.5,27.582], Tmin=(100,'K'), Tmax=(996.266,'K')), NASAPolynomial(coeffs=[10.2984,0.0275406,-1.03818e-05,1.89585e-09,-1.33377e-13,47730.3,-19.2943], Tmin=(996.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Isobutyl) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C1[CH]C[CH]C1=O(10010)',
    structure = SMILES('[CH2]C1[CH]CC=C1[O]'),
    E0 = (305.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83143,0.0327085,3.50915e-05,-7.08122e-08,3.03492e-11,36844,25.4419], Tmin=(100,'K'), Tmax=(941.017,'K')), NASAPolynomial(coeffs=[13.5268,0.0207751,-6.10933e-06,1.04159e-09,-7.4947e-14,32970.1,-39.1586], Tmin=(941.017,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(C=C(C)OJ) + radical(cyclopentene-4) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1CC1[C]=O(10011)',
    structure = SMILES('[CH2][CH]C1CC1[C]=O'),
    E0 = (431.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5886,0.0432148,-4.93183e-06,-2.08235e-08,9.88487e-12,52005.7,28.7249], Tmin=(100,'K'), Tmax=(1048.08,'K')), NASAPolynomial(coeffs=[11.6117,0.0260685,-1.06e-05,1.99666e-09,-1.41734e-13,48745.4,-25.6319], Tmin=(1048.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Cs_S) + radical(RCCJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2][CH]C1C[CH]C1=O(10012)',
    structure = SMILES('[CH2][CH]C1C[CH]C1=O'),
    E0 = (398.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94415,0.0266983,5.65855e-05,-9.19716e-08,3.66849e-11,47962.2,27.1629], Tmin=(100,'K'), Tmax=(968.923,'K')), NASAPolynomial(coeffs=[14.2169,0.022172,-7.83596e-06,1.49995e-09,-1.13182e-13,43418.2,-42.8348], Tmin=(968.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O) + radical(CCJCC=O) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC=C[C]=O(10013)',
    structure = SMILES('[CH2]C=CC=C[C]=O'),
    E0 = (213.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,3000,3100,440,815,1455,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.17025,'amu*angstrom^2'), symmetry=1, barrier=(26.9063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.17104,'amu*angstrom^2'), symmetry=1, barrier=(26.9245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16951,'amu*angstrom^2'), symmetry=1, barrier=(26.8894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00471,0.0617231,-5.34235e-05,2.3453e-08,-4.14143e-12,25739.3,24.7144], Tmin=(100,'K'), Tmax=(1346.54,'K')), NASAPolynomial(coeffs=[14.2019,0.0225202,-9.75339e-06,1.83242e-09,-1.27383e-13,22185.2,-42.8815], Tmin=(1346.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CCJ=O)"""),
)

species(
    label = 'C=C=CC[CH][C]=O(10014)',
    structure = SMILES('C=C=CC[CH][C]=O'),
    E0 = (344.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950,371.154,371.257],'cm^-1')),
        HinderedRotor(inertia=(0.0909989,'amu*angstrom^2'), symmetry=1, barrier=(8.90382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0909397,'amu*angstrom^2'), symmetry=1, barrier=(8.90337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.269233,'amu*angstrom^2'), symmetry=1, barrier=(26.3334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29967,0.0576537,-5.12447e-05,2.47978e-08,-4.93277e-12,41588.6,25.4864], Tmin=(100,'K'), Tmax=(1193.82,'K')), NASAPolynomial(coeffs=[11.0631,0.0249405,-1.01415e-05,1.84454e-09,-1.26084e-13,39257.5,-23.3461], Tmin=(1193.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CCJCHO) + radical(CCCJ=O)"""),
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
    label = '[CH]C=C(3739)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,230.296,230.449,230.657],'cm^-1')),
        HinderedRotor(inertia=(1.34138,'amu*angstrom^2'), symmetry=1, barrier=(50.5187,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36193e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35107e-06,1.16619e-09,-8.2762e-14,44095,-3.44607], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
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
    label = '[CH2]C=C[CH][CH][C]=O(10015)',
    structure = SMILES('[CH2]C=C[CH]C=[C][O]'),
    E0 = (478.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,602.609,602.83,602.961],'cm^-1')),
        HinderedRotor(inertia=(0.291655,'amu*angstrom^2'), symmetry=1, barrier=(75.1758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106574,'amu*angstrom^2'), symmetry=1, barrier=(27.4797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106612,'amu*angstrom^2'), symmetry=1, barrier=(27.4803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49026,0.039342,2.13071e-05,-6.21794e-08,2.85287e-11,57638.6,27.9866], Tmin=(100,'K'), Tmax=(942.386,'K')), NASAPolynomial(coeffs=[16.4744,0.0155637,-4.23045e-06,7.26991e-10,-5.48589e-14,53046.1,-52.7961], Tmin=(942.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJC=C) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=[C]C[CH][C]=O(10016)',
    structure = SMILES('[CH2]C=[C]C[CH][C]=O'),
    E0 = (557.704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1855,455,950,285.427,4000],'cm^-1')),
        HinderedRotor(inertia=(0.35198,'amu*angstrom^2'), symmetry=1, barrier=(20.3507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5297,'amu*angstrom^2'), symmetry=1, barrier=(88.4401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0020688,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52925,'amu*angstrom^2'), symmetry=1, barrier=(88.4403,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27781,0.0594976,-5.62055e-05,2.93436e-08,-6.31326e-12,67174.6,27.2673], Tmin=(100,'K'), Tmax=(1108.13,'K')), NASAPolynomial(coeffs=[10.5949,0.0258656,-1.06798e-05,1.95442e-09,-1.34056e-13,65109.8,-18.6387], Tmin=(1108.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Allyl_P) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=CC[CH][C]=O(10017)',
    structure = SMILES('[CH2][C]=CC[CH][C]=O'),
    E0 = (557.704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,1855,455,950,285.475,4000],'cm^-1')),
        HinderedRotor(inertia=(0.352017,'amu*angstrom^2'), symmetry=1, barrier=(20.3507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52945,'amu*angstrom^2'), symmetry=1, barrier=(88.4404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00206908,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.52966,'amu*angstrom^2'), symmetry=1, barrier=(88.4405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27781,0.0594976,-5.62055e-05,2.93436e-08,-6.31326e-12,67174.6,27.2673], Tmin=(100,'K'), Tmax=(1108.13,'K')), NASAPolynomial(coeffs=[10.5949,0.0258656,-1.06798e-05,1.95442e-09,-1.34056e-13,65109.8,-18.6387], Tmin=(1108.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Allyl_P) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=C[CH]C[C]=O(5328)',
    structure = SMILES('[CH2]C=C[CH]C[C]=O'),
    E0 = (352.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950,354.573,1044.56],'cm^-1')),
        HinderedRotor(inertia=(0.0840196,'amu*angstrom^2'), symmetry=1, barrier=(7.49586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0840197,'amu*angstrom^2'), symmetry=1, barrier=(7.49587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0840179,'amu*angstrom^2'), symmetry=1, barrier=(7.49591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0264725,'amu*angstrom^2'), symmetry=1, barrier=(20.497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22306,0.0574204,-4.48442e-05,1.84806e-08,-3.13231e-12,42466.9,29.7248], Tmin=(100,'K'), Tmax=(1378.77,'K')), NASAPolynomial(coeffs=[12.0671,0.0259601,-1.06173e-05,1.93092e-09,-1.31482e-13,39476.7,-26.0742], Tmin=(1378.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = 'C[C]=CC[CH][C]=O(10018)',
    structure = SMILES('C[C]=CC[CH][C]=O'),
    E0 = (406.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37421,0.0607494,-5.99697e-05,3.60017e-08,-9.30553e-12,48947.1,26.649], Tmin=(100,'K'), Tmax=(913.589,'K')), NASAPolynomial(coeffs=[7.72852,0.0329288,-1.42928e-05,2.67107e-09,-1.84954e-13,47786.1,-3.43277], Tmin=(913.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=[C]CC[C]=O(5330)',
    structure = SMILES('[CH2]C=[C]CC[C]=O'),
    E0 = (390.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1855,455,950,357.601,357.601],'cm^-1')),
        HinderedRotor(inertia=(0.0702931,'amu*angstrom^2'), symmetry=1, barrier=(6.37876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702931,'amu*angstrom^2'), symmetry=1, barrier=(6.37876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702932,'amu*angstrom^2'), symmetry=1, barrier=(6.37876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179162,'amu*angstrom^2'), symmetry=1, barrier=(16.258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27283,0.0638107,-6.57226e-05,4.04754e-08,-1.06761e-11,47022.1,27.5359], Tmin=(100,'K'), Tmax=(897.349,'K')), NASAPolynomial(coeffs=[8.06495,0.0335346,-1.51138e-05,2.87698e-09,-2.01323e-13,45803.1,-4.49649], Tmin=(897.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=[C]C[CH][C]=O(10019)',
    structure = SMILES('CC=[C]C[CH][C]=O'),
    E0 = (406.205,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,1855,455,950,345.32,345.328],'cm^-1')),
        HinderedRotor(inertia=(0.00141359,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0957173,'amu*angstrom^2'), symmetry=1, barrier=(8.09951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0957222,'amu*angstrom^2'), symmetry=1, barrier=(8.09984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0957154,'amu*angstrom^2'), symmetry=1, barrier=(8.09961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37421,0.0607494,-5.99697e-05,3.60017e-08,-9.30553e-12,48947.1,26.649], Tmin=(100,'K'), Tmax=(913.589,'K')), NASAPolynomial(coeffs=[7.72852,0.0329288,-1.42928e-05,2.67107e-09,-1.84954e-13,47786.1,-3.43277], Tmin=(913.589,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=C[CH][CH]C=O(10020)',
    structure = SMILES('[CH2]C=C[CH]C=C[O]'),
    E0 = (238.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27925,0.0375948,4.81722e-05,-1.00984e-07,4.47846e-11,28817.6,25.6155], Tmin=(100,'K'), Tmax=(930.859,'K')), NASAPolynomial(coeffs=[20.6022,0.0114585,-1.3953e-06,1.7802e-10,-1.93605e-14,22755.2,-79.4629], Tmin=(930.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]C[CH]C=O(10021)',
    structure = SMILES('[CH2]C=[C]CC=C[O]'),
    E0 = (374.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,1685,370,310.129,310.767,311.506],'cm^-1')),
        HinderedRotor(inertia=(0.313605,'amu*angstrom^2'), symmetry=1, barrier=(21.4174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31402,'amu*angstrom^2'), symmetry=1, barrier=(21.4075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311814,'amu*angstrom^2'), symmetry=1, barrier=(21.4138,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916072,0.0547191,-1.47179e-05,-2.639e-08,1.57386e-11,45193.1,27.7221], Tmin=(100,'K'), Tmax=(964.164,'K')), NASAPolynomial(coeffs=[17.308,0.0178921,-5.9287e-06,1.071e-09,-7.8044e-14,40583,-58.2762], Tmin=(964.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC[C]=O(5331)',
    structure = SMILES('[CH2][C]=CCC[C]=O'),
    E0 = (390.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1855,455,950,357.601,357.601],'cm^-1')),
        HinderedRotor(inertia=(0.0702931,'amu*angstrom^2'), symmetry=1, barrier=(6.37876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702931,'amu*angstrom^2'), symmetry=1, barrier=(6.37876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0702932,'amu*angstrom^2'), symmetry=1, barrier=(6.37876,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179162,'amu*angstrom^2'), symmetry=1, barrier=(16.258,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27283,0.0638107,-6.57226e-05,4.04754e-08,-1.06761e-11,47022.1,27.5359], Tmin=(100,'K'), Tmax=(897.349,'K')), NASAPolynomial(coeffs=[8.06495,0.0335346,-1.51138e-05,2.87698e-09,-2.01323e-13,45803.1,-4.49649], Tmin=(897.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=C[CH][CH][C]=O(10022)',
    structure = SMILES('C[CH]C=C[CH][C]=O'),
    E0 = (311.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48667,0.039376,2.1295e-05,-5.40161e-08,2.21245e-11,37556.2,25.5063], Tmin=(100,'K'), Tmax=(1028.65,'K')), NASAPolynomial(coeffs=[15.2159,0.0235917,-1.05212e-05,2.14113e-09,-1.60646e-13,32742.2,-50.787], Tmin=(1028.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_S) + radical(CCJC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=CC[CH]C=O(10023)',
    structure = SMILES('[CH2][C]=CCC=C[O]'),
    E0 = (374.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,1685,370,310.129,310.767,311.506],'cm^-1')),
        HinderedRotor(inertia=(0.313605,'amu*angstrom^2'), symmetry=1, barrier=(21.4174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31402,'amu*angstrom^2'), symmetry=1, barrier=(21.4075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311814,'amu*angstrom^2'), symmetry=1, barrier=(21.4138,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916072,0.0547191,-1.47179e-05,-2.639e-08,1.57386e-11,45193.1,27.7221], Tmin=(100,'K'), Tmax=(964.164,'K')), NASAPolynomial(coeffs=[17.308,0.0178921,-5.9287e-06,1.071e-09,-7.8044e-14,40583,-58.2762], Tmin=(964.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C2H3(60)',
    structure = SMILES('[CH]=C'),
    E0 = (286.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,728.586,728.586,3521.64],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8307,-0.00147628,3.09002e-05,-3.80476e-08,1.43171e-11,34524.2,5.61949], Tmin=(100,'K'), Tmax=(933.662,'K')), NASAPolynomial(coeffs=[5.36086,0.00527345,-1.3196e-06,2.21564e-10,-1.68768e-14,33658.6,-4.76326], Tmin=(933.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=CC1CC1[C]=O(9970)',
    structure = SMILES('C=CC1CC1[C]=O'),
    E0 = (160.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6935,0.0353824,2.86298e-05,-6.34663e-08,2.69839e-11,19421.1,24.5116], Tmin=(100,'K'), Tmax=(968.637,'K')), NASAPolynomial(coeffs=[14.3093,0.0209246,-7.26853e-06,1.35703e-09,-1.00401e-13,15211.3,-45.0651], Tmin=(968.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=CCC=C[C]=O(10024)',
    structure = SMILES('C=CC[CH]C=C=O'),
    E0 = (125.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46376,0.0557229,-4.04046e-05,1.5268e-08,-2.42406e-12,15226.1,25.6131], Tmin=(100,'K'), Tmax=(1420.92,'K')), NASAPolynomial(coeffs=[10.6538,0.0298518,-1.30932e-05,2.45388e-09,-1.69477e-13,12614.4,-21.9518], Tmin=(1420.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(125.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(CCJC(C)=C=O)"""),
)

species(
    label = 'C=CC=CC[C]=O(5341)',
    structure = SMILES('C=CC=CC[C]=O'),
    E0 = (118.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.781505,'amu*angstrom^2'), symmetry=1, barrier=(17.9683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.782468,'amu*angstrom^2'), symmetry=1, barrier=(17.9905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.784775,'amu*angstrom^2'), symmetry=1, barrier=(18.0435,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19676,0.048709,-7.29039e-06,-2.55696e-08,1.26834e-11,14375.7,25.1849], Tmin=(100,'K'), Tmax=(1042.89,'K')), NASAPolynomial(coeffs=[15.5436,0.0219839,-9.55904e-06,1.9032e-09,-1.40484e-13,9844.12,-52.012], Tmin=(1042.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]C1C[CH][CH]C1(10025)',
    structure = SMILES('O=[C]C1C[CH][CH]C1'),
    E0 = (326.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5005,0.0267355,1.88816e-05,-2.97497e-08,9.79727e-12,39320.3,27.4212], Tmin=(100,'K'), Tmax=(1128.67,'K')), NASAPolynomial(coeffs=[4.46115,0.0358762,-1.46488e-05,2.68556e-09,-1.84658e-13,37853,13.1851], Tmin=(1128.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-OdCsH) + ring(Cyclopentane) + radical(cyclopentane) + radical(cyclopentane) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=C1[CH]C[CH][CH]C1(10026)',
    structure = SMILES('[O]C1=CC[CH][CH]C1'),
    E0 = (290.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26562,0.0292527,2.27953e-05,-4.01958e-08,1.47788e-11,35056.5,24.1237], Tmin=(100,'K'), Tmax=(1041.2,'K')), NASAPolynomial(coeffs=[6.87779,0.032835,-1.30527e-05,2.40578e-09,-1.67926e-13,32941.5,-3.85819], Tmin=(1041.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = 'C=CC[CH][CH][C]=O(9967)',
    structure = SMILES('C=CC[CH]C=[C][O]'),
    E0 = (377.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,626.551,626.552,626.559,626.569],'cm^-1')),
        HinderedRotor(inertia=(0.154469,'amu*angstrom^2'), symmetry=1, barrier=(3.55155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0750306,'amu*angstrom^2'), symmetry=1, barrier=(20.9024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.909118,'amu*angstrom^2'), symmetry=1, barrier=(20.9024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27843,0.0493894,-1.33895e-05,-1.7535e-08,1.00144e-11,45477,28.1022], Tmin=(100,'K'), Tmax=(1016.14,'K')), NASAPolynomial(coeffs=[13.486,0.0239851,-9.32375e-06,1.73372e-09,-1.23164e-13,41826.8,-36.7416], Tmin=(1016.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]CC[CH][C]=O(10027)',
    structure = SMILES('C=[C]CC[CH][C]=O'),
    E0 = (418.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,342.751,342.752,3415.49],'cm^-1')),
        HinderedRotor(inertia=(0.158313,'amu*angstrom^2'), symmetry=1, barrier=(13.198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158315,'amu*angstrom^2'), symmetry=1, barrier=(13.198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00143496,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.949821,'amu*angstrom^2'), symmetry=1, barrier=(79.1829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30612,0.0608379,-5.66796e-05,3.05298e-08,-6.95928e-12,50423.6,27.3273], Tmin=(100,'K'), Tmax=(1034.61,'K')), NASAPolynomial(coeffs=[9.08846,0.0307499,-1.30573e-05,2.42108e-09,-1.67159e-13,48813.3,-10.4827], Tmin=(1034.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CCC[CH][C]=O(10028)',
    structure = SMILES('[CH]=CCC[CH][C]=O'),
    E0 = (427.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,3120,650,792.5,1650,307.777,1345.71],'cm^-1')),
        HinderedRotor(inertia=(0.00166985,'amu*angstrom^2'), symmetry=1, barrier=(15.1971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.661025,'amu*angstrom^2'), symmetry=1, barrier=(15.1983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.662067,'amu*angstrom^2'), symmetry=1, barrier=(15.2222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.64933,'amu*angstrom^2'), symmetry=1, barrier=(83.9052,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13788,0.0613322,-5.43057e-05,2.63723e-08,-5.28219e-12,51545.4,27.8551], Tmin=(100,'K'), Tmax=(1184.09,'K')), NASAPolynomial(coeffs=[11.2527,0.027163,-1.10201e-05,2.00148e-09,-1.36682e-13,49150,-22.6519], Tmin=(1184.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C=CCC[C]=O(5318)',
    structure = SMILES('[CH]C=CCC[C]=O'),
    E0 = (371.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39138,0.06072,-4.63156e-05,2.07975e-08,-4.22267e-12,44774.4,27.4582], Tmin=(100,'K'), Tmax=(1095.05,'K')), NASAPolynomial(coeffs=[7.13312,0.0397467,-1.75867e-05,3.30744e-09,-2.29725e-13,43516.9,-0.763788], Tmin=(1095.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH]C[CH]C=O(10029)',
    structure = SMILES('[CH]C=CCC=C[O]'),
    E0 = (356.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,232.993,232.993,232.993,232.993,232.993,232.994],'cm^-1')),
        HinderedRotor(inertia=(1.26969,'amu*angstrom^2'), symmetry=1, barrier=(48.9114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26969,'amu*angstrom^2'), symmetry=1, barrier=(48.9114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26969,'amu*angstrom^2'), symmetry=1, barrier=(48.9114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948504,0.052543,1.91086e-06,-4.29321e-08,2.10012e-11,42949.3,27.9604], Tmin=(100,'K'), Tmax=(967.507,'K')), NASAPolynomial(coeffs=[16.1467,0.0245388,-8.67182e-06,1.56882e-09,-1.12261e-13,38378.3,-53.2844], Tmin=(967.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=CCC=C=O(5326)',
    structure = SMILES('[CH2]C=CCC=C=O'),
    E0 = (117.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,2120,512.5,787.5,180,531.315],'cm^-1')),
        HinderedRotor(inertia=(0.60249,'amu*angstrom^2'), symmetry=1, barrier=(13.8524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.406524,'amu*angstrom^2'), symmetry=1, barrier=(9.3468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02516,'amu*angstrom^2'), symmetry=1, barrier=(23.5704,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23481,0.0597159,-4.88265e-05,2.14386e-08,-3.91701e-12,14261.1,24.0971], Tmin=(100,'K'), Tmax=(1275.07,'K')), NASAPolynomial(coeffs=[11.1768,0.028527,-1.21357e-05,2.25487e-09,-1.55699e-13,11725.8,-26.2831], Tmin=(1275.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Allyl_P)"""),
)

species(
    label = '[O]C1=CCC=CC1(9599)',
    structure = SMILES('[O]C1=CCC=CC1'),
    E0 = (17.0638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86576,0.027947,5.57732e-05,-9.46644e-08,3.88509e-11,2146.41,21.8715], Tmin=(100,'K'), Tmax=(952.465,'K')), NASAPolynomial(coeffs=[15.2611,0.0193966,-5.88984e-06,1.0815e-09,-8.26406e-14,-2569.19,-53.4603], Tmin=(952.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.0638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cyclohexadiene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C=CCC=C[O](9600)',
    structure = SMILES('C=C=CCC=C[O]'),
    E0 = (161.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990037,0.0522398,-7.46629e-06,-3.39513e-08,1.84062e-11,19604.9,25.7561], Tmin=(100,'K'), Tmax=(962.78,'K')), NASAPolynomial(coeffs=[17.4007,0.0175999,-5.75338e-06,1.0465e-09,-7.7133e-14,14890.4,-60.8666], Tmin=(962.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(161.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = '[CH2][CH]CC[C]=[C][O](10030)',
    structure = SMILES('[CH2][CH]CC[C][C]=O'),
    E0 = (733.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918776,0.0760353,-0.000114671,1.03004e-07,-3.64756e-11,88287.1,32.4301], Tmin=(100,'K'), Tmax=(850.283,'K')), NASAPolynomial(coeffs=[5.379,0.0380123,-1.75319e-05,3.27185e-09,-2.22158e-13,88144.6,15.2578], Tmin=(850.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C[CH]C[C]=[C][O](10031)',
    structure = SMILES('[CH2]C[CH]C[C][C]=O'),
    E0 = (733.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.958964,0.0764755,-0.000119734,1.11216e-07,-4.02421e-11,88285.8,32.114], Tmin=(100,'K'), Tmax=(852.672,'K')), NASAPolynomial(coeffs=[4.28949,0.0399583,-1.87386e-05,3.5153e-09,-2.39036e-13,88477.4,21.0308], Tmin=(852.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1[CH]CC=[C]O1(10032)',
    structure = SMILES('[CH2]C1[CH]CC=[C]O1'),
    E0 = (459.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61415,0.0362184,3.36223e-05,-7.66234e-08,3.46756e-11,55347.6,23.3408], Tmin=(100,'K'), Tmax=(917.107,'K')), NASAPolynomial(coeffs=[15.8372,0.0164634,-3.21803e-06,4.2438e-10,-3.0102e-14,50960.8,-53.7402], Tmin=(917.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(CCJCO) + radical(CJC(C)OC) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]C1CC=[C]O1(10033)',
    structure = SMILES('[CH2][CH]C1CC=[C]O1'),
    E0 = (456.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85181,0.0284013,5.5881e-05,-1.01127e-07,4.39799e-11,54958.9,24.9131], Tmin=(100,'K'), Tmax=(911.71,'K')), NASAPolynomial(coeffs=[16.3144,0.0141473,-1.61203e-06,1.01788e-10,-8.27644e-15,50277.1,-54.7375], Tmin=(911.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CCJCO) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=CCC#C[O](10034)',
    structure = SMILES('[CH2]C=CC[C]=C=O'),
    E0 = (320.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,235.804,236.674],'cm^-1')),
        HinderedRotor(inertia=(0.451036,'amu*angstrom^2'), symmetry=1, barrier=(18.0444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.456385,'amu*angstrom^2'), symmetry=1, barrier=(18.0413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.41384,'amu*angstrom^2'), symmetry=1, barrier=(96.287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24023,0.0605004,-5.51439e-05,2.71453e-08,-5.51439e-12,38589.1,23.7892], Tmin=(100,'K'), Tmax=(1164.81,'K')), NASAPolynomial(coeffs=[11.1167,0.0265843,-1.1468e-05,2.14791e-09,-1.49272e-13,36288.2,-25.366], Tmin=(1164.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Allyl_P) + radical(CCCJ=C=O)"""),
)

species(
    label = '[CH2]C=CC[C]=C[O](10035)',
    structure = SMILES('[CH2]C=CC[C]=C[O]'),
    E0 = (374.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,1685,370,310.129,310.767,311.506],'cm^-1')),
        HinderedRotor(inertia=(0.313605,'amu*angstrom^2'), symmetry=1, barrier=(21.4174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31402,'amu*angstrom^2'), symmetry=1, barrier=(21.4075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311814,'amu*angstrom^2'), symmetry=1, barrier=(21.4138,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.916072,0.0547191,-1.47179e-05,-2.639e-08,1.57386e-11,45193.1,27.7221], Tmin=(100,'K'), Tmax=(964.164,'K')), NASAPolynomial(coeffs=[17.308,0.0178921,-5.9287e-06,1.071e-09,-7.8044e-14,40583,-58.2762], Tmin=(964.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC[C]=[C]O(10036)',
    structure = SMILES('[CH2]C=CC[C]=[C]O'),
    E0 = (473.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,326.095,326.095],'cm^-1')),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703702,0.0635789,-4.54885e-05,6.83891e-09,3.96591e-12,57017.5,30.2025], Tmin=(100,'K'), Tmax=(966.332,'K')), NASAPolynomial(coeffs=[16.491,0.0182112,-6.08291e-06,1.05156e-09,-7.28078e-14,53033.4,-50.2487], Tmin=(966.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=C[CH]C=[C]O(10037)',
    structure = SMILES('[CH2]C=C[CH]C=[C]O'),
    E0 = (336.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,304.987,305.034],'cm^-1')),
        HinderedRotor(inertia=(0.479964,'amu*angstrom^2'), symmetry=1, barrier=(31.6994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.480287,'amu*angstrom^2'), symmetry=1, barrier=(31.6983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00181318,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.480025,'amu*angstrom^2'), symmetry=1, barrier=(31.6983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07595,0.0463379,1.78596e-05,-6.84306e-08,3.3341e-11,40641.6,28.064], Tmin=(100,'K'), Tmax=(923.081,'K')), NASAPolynomial(coeffs=[19.778,0.0117919,-1.55846e-06,1.60819e-10,-1.43186e-14,35208,-71.3954], Tmin=(923.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CCJC=C) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'CC=CC[C]=[C][O](10038)',
    structure = SMILES('CC=CC[C][C]=O'),
    E0 = (449.049,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,318.548,318.548,318.549],'cm^-1')),
        HinderedRotor(inertia=(0.14232,'amu*angstrom^2'), symmetry=1, barrier=(10.2481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142319,'amu*angstrom^2'), symmetry=1, barrier=(10.2481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142319,'amu*angstrom^2'), symmetry=1, barrier=(10.248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142319,'amu*angstrom^2'), symmetry=1, barrier=(10.2481,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23752,0.0649474,-6.96586e-05,4.5146e-08,-1.25269e-11,54103.9,26.8144], Tmin=(100,'K'), Tmax=(856.174,'K')), NASAPolynomial(coeffs=[7.85775,0.0340181,-1.54715e-05,2.95305e-09,-2.06772e-13,52970.3,-4.09639], Tmin=(856.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=[C]CC=[C]O(10039)',
    structure = SMILES('[CH2]C=[C]CC=[C]O'),
    E0 = (473.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,326.095,326.095],'cm^-1')),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703702,0.0635789,-4.54885e-05,6.83891e-09,3.96591e-12,57017.5,30.2025], Tmin=(100,'K'), Tmax=(966.332,'K')), NASAPolynomial(coeffs=[16.491,0.0182112,-6.08291e-06,1.05156e-09,-7.28078e-14,53033.4,-50.2487], Tmin=(966.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=CCC=[C]O(10040)',
    structure = SMILES('[CH2][C]=CCC=[C]O'),
    E0 = (473.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,326.095,326.095],'cm^-1')),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216948,'amu*angstrom^2'), symmetry=1, barrier=(16.3708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.703702,0.0635789,-4.54885e-05,6.83891e-09,3.96591e-12,57017.5,30.2025], Tmin=(100,'K'), Tmax=(966.332,'K')), NASAPolynomial(coeffs=[16.491,0.0182112,-6.08291e-06,1.05156e-09,-7.28078e-14,53033.4,-50.2487], Tmin=(966.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([CH][C]=O)C=C(5308)',
    structure = SMILES('[CH2]C([CH][C]=O)C=C'),
    E0 = (377.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,491.276,3175.32],'cm^-1')),
        HinderedRotor(inertia=(0.0484541,'amu*angstrom^2'), symmetry=1, barrier=(8.30035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.466265,'amu*angstrom^2'), symmetry=1, barrier=(79.8763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0484603,'amu*angstrom^2'), symmetry=1, barrier=(8.29861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.465793,'amu*angstrom^2'), symmetry=1, barrier=(79.8784,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3598.66,'J/mol'), sigma=(6.12618,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.10 K, Pc=35.52 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896554,0.0607427,-5.26077e-05,2.49297e-08,-4.75266e-12,45547.1,29.8552], Tmin=(100,'K'), Tmax=(1269.15,'K')), NASAPolynomial(coeffs=[13.033,0.022491,-7.397e-06,1.18052e-09,-7.43624e-14,42466.6,-31.5884], Tmin=(1269.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC1CC=C1[O](10041)',
    structure = SMILES('C=CC1C[CH]C1=O'),
    E0 = (119.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96641,0.0269733,5.41839e-05,-8.74625e-08,3.44868e-11,14464.7,22.6067], Tmin=(100,'K'), Tmax=(977.337,'K')), NASAPolynomial(coeffs=[13.4639,0.0239953,-8.89615e-06,1.71215e-09,-1.27921e-13,10112.2,-43.3682], Tmin=(977.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone) + radical(CCJC=O)"""),
)

species(
    label = 'C=CC=CC=C[O](10042)',
    structure = SMILES('[CH2]C=CC=CC=O'),
    E0 = (52.4633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19817,0.051949,-2.18837e-05,-5.89876e-09,4.9049e-12,6418.88,24.4996], Tmin=(100,'K'), Tmax=(1111.19,'K')), NASAPolynomial(coeffs=[13.2992,0.0262162,-1.12126e-05,2.13734e-09,-1.51482e-13,2628.95,-40.109], Tmin=(1111.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.4633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ)"""),
)

species(
    label = '[C]1=CC[CH][CH]CO1(10043)',
    structure = SMILES('[C]1=CC[CH][CH]CO1'),
    E0 = (470.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06046,-0.0126817,0.000240267,-3.45209e-07,1.44926e-10,56661,27.815], Tmin=(100,'K'), Tmax=(899.122,'K')), NASAPolynomial(coeffs=[38.1028,-0.0278263,2.32967e-05,-4.72389e-09,3.13902e-13,44310.6,-174.874], Tmin=(899.122,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(RCCJCC) + radical(CCJCO) + radical(C=CJO)"""),
)

species(
    label = 'C=CCC[C]=[C][O](10044)',
    structure = SMILES('C=CCC[C][C]=O'),
    E0 = (461.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,223.44,223.441,223.442,223.444],'cm^-1')),
        HinderedRotor(inertia=(0.00337662,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00337653,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423628,'amu*angstrom^2'), symmetry=1, barrier=(15.0084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423629,'amu*angstrom^2'), symmetry=1, barrier=(15.0084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21031,0.0645332,-6.45213e-05,3.71526e-08,-9.04734e-12,55578.7,27.3476], Tmin=(100,'K'), Tmax=(971.765,'K')), NASAPolynomial(coeffs=[9.1073,0.0320268,-1.4344e-05,2.72852e-09,-1.91084e-13,54043.9,-10.5244], Tmin=(971.765,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C[CH]CC=[C]O(10045)',
    structure = SMILES('[CH]C=CCC=[C]O'),
    E0 = (454.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734479,0.0614227,-2.89307e-05,-9.61004e-09,9.1886e-12,54773.8,30.4467], Tmin=(100,'K'), Tmax=(971.389,'K')), NASAPolynomial(coeffs=[15.3367,0.024846,-8.81912e-06,1.54774e-09,-1.06889e-13,50825.7,-45.2964], Tmin=(971.389,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(AllylJ2_triplet)"""),
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
    E0 = (319.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (511.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (755.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (889.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (812.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (750.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (326.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (344.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (344.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (671.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (545.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (445.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (396.146,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (431.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (408.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (432.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (568.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (496.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (711.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (690.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (769.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (769.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (547.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (527.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (532.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (568.126,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (463.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (524.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (434.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (478.898,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (407.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (890.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (327.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (342.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (347.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (418.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (380.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (555.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (489.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (576.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (572.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (404.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (507.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (319.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (328.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (344.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (796.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (758.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (459.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (456.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (546.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (493.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (578.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (623.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (442.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (510.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (624.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (547.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (537.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (327.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (398.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (470.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (585.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (605.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['HCCO(2227)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C=CC([CH2])[C]=O(10002)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(20)', '[CH]=CC[CH][C]=O(10003)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C=CC[C][C]=O(10004)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]C=CC[CH][C]=O(10005)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['O=[C]C1CC=CC1(10006)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['CC=CC=C[C]=O(10007)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['C=C=CCC[C]=O(5327)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]C[CH][CH][C]=O(10008)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['O=[C][CH]CC1[CH]C1(10009)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[CH2]C1[CH]CC1[C]=O(9973)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.66065e+07,'s^-1'), n=1.25778, Ea=(125.799,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[CH2]C1[CH]C[CH]C1=O(10010)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.72761e+09,'s^-1'), n=0.568593, Ea=(76.2838,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_pri;radadd_intra] for rate rule [R5_linear;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[CH2][CH]C1CC1[C]=O(10011)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.86885e+09,'s^-1'), n=0.611527, Ea=(111.745,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 109.7 to 111.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[CH2][CH]C1C[CH]C1=O(10012)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.68391e+08,'s^-1'), n=0.691643, Ea=(88.1932,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra;radadd_intra] for rate rule [R5;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]C=CC=C[C]=O(10013)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', 'C=C=CC[CH][C]=O(10014)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH2CHCO(3668)', '[CH]C=C(3739)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH][C]=O(9632)', '[CH]C=C(3739)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.30526e+07,'m^3/(mol*s)'), n=-0.283333, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-6.50053689359e-11, var=0.305422193575, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C
    Total Standard Deviation in ln(k): 1.10791715097
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2]C=C[CH][CH][C]=O(10015)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2]C=[C]C[CH][C]=O(10016)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][C]=CC[CH][C]=O(10017)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C=C[CH]C[C]=O(5328)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.01225e+10,'s^-1'), n=0.845153, Ea=(195.403,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['C[C]=CC[CH][C]=O(10018)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C=[C]CC[C]=O(5330)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CC=[C]C[CH][C]=O(10019)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[CH2]C=C[CH][CH]C=O(10020)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.34586e+06,'s^-1'), n=1.99734, Ea=(143.275,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3HJ;CO_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C=[C]C[CH]C=O(10021)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.13341e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_Cd;XH_out] for rate rule [R4HJ_2;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=CCC[C]=O(5331)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['CC=C[CH][CH][C]=O(10022)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(987669,'s^-1'), n=1.76395, Ea=(159.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]=CC[CH]C=O(10023)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_3;Cd_rad_out;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C2H3(60)', '[CH]C[CH][C]=O(9634)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['C=CC1CC1[C]=O(9970)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['C=CCC=C[C]=O(10024)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['C=CC=CC[C]=O(5341)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['O=[C]C1C[CH][CH]C1(10025)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(7.16776e+08,'s^-1'), n=0.66239, Ea=(99.0943,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['O=C1[CH]C[CH][CH]C1(10026)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(60.8864,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=[C][O](6861)', 'butadiene13(1350)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.18583,'m^3/(mol*s)'), n=2.36967, Ea=(33.8986,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds-CdH;YJ] for rate rule [Cds-HH_Cds-CdH;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=CC[CH][CH][C]=O(9967)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.42705e+08,'s^-1'), n=1.50595, Ea=(111.926,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C]CC[CH][C]=O(10027)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=CCC[CH][C]=O(10028)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]C=CCC[C]=O(5318)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C[CH]C[CH]C=O(10029)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[CH2]C=CCC=C=O(5326)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[O]C1=CCC=CC1(9599)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['C=C=CCC=C[O](9600)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2][CH]CC[C]=[C][O](10030)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C[CH]C[C]=[C][O](10031)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[CH2]C1[CH]CC=[C]O1(10032)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.18543e+10,'s^-1'), n=0.209443, Ea=(139.485,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 137.0 to 139.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[CH2][CH]C1CC=[C]O1(10033)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(4.34604e+09,'s^-1'), n=0.547362, Ea=(136.303,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra;radadd_intra_O] + [R6;doublebond_intra;radadd_intra] for rate rule [R6;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 132.7 to 136.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction51',
    reactants = ['H(3)', '[CH2]C=CCC#C[O](10034)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction52',
    reactants = ['HCCO(2227)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(0.113192,'m^3/(mol*s)'), n=2.418, Ea=(54.0165,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CsJ-CdHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2]C=CC[C]=C[O](10035)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C=CC[C]=[C]O(10036)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C=C[CH]C=[C]O(10037)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.78e+06,'s^-1'), n=1.75, Ea=(105.855,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/Cd;XH_out] for rate rule [R4HJ_2;C_rad_out_H/Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['CC=CC[C]=[C][O](10038)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(770185,'s^-1'), n=1.81245, Ea=(61.132,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5H_RSMS;Cd_rad_out;Cs_H_out_2H] + [R5H_SSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R5H_SSMS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C=[C]CC=[C]O(10039)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2][C]=CCC=[C]O(10040)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_4;Cd_rad_out;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['C=CC1CC=C1[O](10041)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['C=CC=CC=C[O](10042)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH2]C=CC[CH][C]=O(5329)'],
    products = ['[C]1=CC[CH][CH]CO1(10043)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(2.76476e+07,'s^-1'), n=0.815689, Ea=(150.23,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 142.1 to 150.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction63',
    reactants = ['C=CCC[C]=[C][O](10044)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.29711e+07,'s^-1'), n=1.52333, Ea=(124.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH]=C[CH]CC=[C]O(10045)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2471',
    isomers = [
        '[CH2]C=CC[CH][C]=O(5329)',
    ],
    reactants = [
        ('HCCO(2227)', 'butadiene13(1350)'),
        ('HCCO(2227)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2471',
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


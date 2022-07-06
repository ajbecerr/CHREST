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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.896554,0.0607427,-5.26077e-05,2.49297e-08,-4.75266e-12,45547.1,29.8552], Tmin=(100,'K'), Tmax=(1269.15,'K')), NASAPolynomial(coeffs=[13.033,0.022491,-7.397e-06,1.18052e-09,-7.43624e-14,42466.6,-31.5884], Tmin=(1269.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2][CH]C([C]=O)C=C(9966)',
    structure = SMILES('[CH2][CH]C([C]=O)C=C'),
    E0 = (410.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,278.169,278.181],'cm^-1')),
        HinderedRotor(inertia=(0.0021786,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00217833,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00217876,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00217862,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46945,0.0557083,-4.3267e-05,1.81267e-08,-3.20571e-12,49409.7,30.3143], Tmin=(100,'K'), Tmax=(1295.61,'K')), NASAPolynomial(coeffs=[10.0415,0.029243,-1.26262e-05,2.35995e-09,-1.63314e-13,47188.5,-13.2607], Tmin=(1295.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCC=O) + radical(RCCJ) + radical(CC(C)CJ=O)"""),
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
    label = '[CH2]C=C[CH][C]=O(9826)',
    structure = SMILES('[CH2]C=C[CH][C]=O'),
    E0 = (344.424,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,579.211],'cm^-1')),
        HinderedRotor(inertia=(0.109865,'amu*angstrom^2'), symmetry=1, barrier=(26.367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109836,'amu*angstrom^2'), symmetry=1, barrier=(26.387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110789,'amu*angstrom^2'), symmetry=1, barrier=(26.3527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96036,0.0299846,2.4079e-05,-5.41179e-08,2.23924e-11,41511.2,22.2791], Tmin=(100,'K'), Tmax=(1008.12,'K')), NASAPolynomial(coeffs=[14.799,0.0148331,-6.62881e-06,1.40459e-09,-1.09437e-13,37104,-48.7837], Tmin=(1008.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(Allyl_P) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([C][C]=O)C=C(9968)',
    structure = SMILES('[CH2]C([C][C]=O)C=C'),
    E0 = (658.405,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,259.233,259.241,259.252],'cm^-1')),
        HinderedRotor(inertia=(1.49162,'amu*angstrom^2'), symmetry=1, barrier=(71.1421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00250822,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00250827,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.214458,'amu*angstrom^2'), symmetry=1, barrier=(10.2276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04946,0.0657575,-7.39816e-05,4.65634e-08,-1.18775e-11,79293.3,29.4546], Tmin=(100,'K'), Tmax=(952.031,'K')), NASAPolynomial(coeffs=[10.7116,0.0251651,-1.00308e-05,1.78533e-09,-1.19952e-13,77453.4,-16.6857], Tmin=(952.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(658.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C([CH][C]=O)C=C(9969)',
    structure = SMILES('[CH]C([CH][C]=O)C=C'),
    E0 = (620.852,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.934732,0.0608243,-5.50719e-05,2.61424e-08,-4.9466e-12,74787.1,28.9413], Tmin=(100,'K'), Tmax=(1279.35,'K')), NASAPolynomial(coeffs=[13.9642,0.0200865,-7.30793e-06,1.2527e-09,-8.28519e-14,71453.3,-37.1281], Tmin=(1279.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.852,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = 'C=CC(=C)C[C]=O(5305)',
    structure = SMILES('C=CC(=C)C[C]=O'),
    E0 = (117.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1855,455,950,231.882,231.882],'cm^-1')),
        HinderedRotor(inertia=(0.485105,'amu*angstrom^2'), symmetry=1, barrier=(18.5096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.485105,'amu*angstrom^2'), symmetry=1, barrier=(18.5096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.485106,'amu*angstrom^2'), symmetry=1, barrier=(18.5096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0234,0.0528587,-1.63228e-05,-1.81531e-08,1.05442e-11,14193.4,24.6015], Tmin=(100,'K'), Tmax=(1042.87,'K')), NASAPolynomial(coeffs=[16.2801,0.0210796,-9.07354e-06,1.79863e-09,-1.32534e-13,9557.21,-56.6148], Tmin=(1042.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC(C)=C[C]=O(9971)',
    structure = SMILES('C=C[C](C)C=C=O'),
    E0 = (64.6108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08946,0.0584489,-4.37233e-05,1.69285e-08,-2.67624e-12,7880.11,24.1875], Tmin=(100,'K'), Tmax=(1476.6,'K')), NASAPolynomial(coeffs=[13.2129,0.025607,-1.03605e-05,1.86545e-09,-1.25901e-13,4299.87,-39.0258], Tmin=(1476.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.6108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCd(CCO)H) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=CCJ(C)C=C=O)"""),
)

species(
    label = '[CH2]C[C]([CH2])[CH][C]=O(9972)',
    structure = SMILES('[CH2]C[C]([CH2])C=[C][O]'),
    E0 = (642.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,432.458,432.47,432.472],'cm^-1')),
        HinderedRotor(inertia=(0.000901346,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000901331,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000901373,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000901359,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16563,0.0560714,-3.59609e-05,5.52499e-09,2.71332e-12,77402.9,31.4076], Tmin=(100,'K'), Tmax=(963.201,'K')), NASAPolynomial(coeffs=[12.0659,0.024632,-8.53392e-06,1.44608e-09,-9.64179e-14,74661.7,-24.1009], Tmin=(963.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(642.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_T) + radical(RCCJ) + radical(Isobutyl) + radical(C=CJO)"""),
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
    label = 'O=[C][CH]C1[CH]CC1(6406)',
    structure = SMILES('O=[C][CH]C1[CH]CC1'),
    E0 = (389.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3025,407.5,1350,352.5,1855,455,950,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8805,0.0348006,2.04557e-05,-4.70744e-08,1.92804e-11,46968.9,25.9474], Tmin=(100,'K'), Tmax=(996.263,'K')), NASAPolynomial(coeffs=[10.8218,0.0270763,-1.03352e-05,1.91658e-09,-1.36547e-13,43789,-24.173], Tmin=(996.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(cyclobutane) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1[CH]CC(=O)[CH]1(9974)',
    structure = SMILES('[CH2]C1[CH]CC([O])=C1'),
    E0 = (305.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83143,0.0327085,3.50915e-05,-7.08122e-08,3.03492e-11,36844,25.4419], Tmin=(100,'K'), Tmax=(941.017,'K')), NASAPolynomial(coeffs=[13.5268,0.0207751,-6.10933e-06,1.04159e-09,-7.4947e-14,32970.1,-39.1586], Tmin=(941.017,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(C=C(C)OJ) + radical(cyclopentene-4) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1C([CH2])C1[C]=O(9975)',
    structure = SMILES('[CH2]C1C([CH2])C1[C]=O'),
    E0 = (432.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49144,0.0435942,5.27953e-06,-4.29088e-08,2.17155e-11,52172,27.3198], Tmin=(100,'K'), Tmax=(916.226,'K')), NASAPolynomial(coeffs=[13.9068,0.0196852,-5.17213e-06,7.81993e-10,-5.22046e-14,48625.4,-38.4296], Tmin=(916.226,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Isobutyl) + radical(Isobutyl) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C1CC1[CH][C]=O(9976)',
    structure = SMILES('[CH2]C1CC1[CH][C]=O'),
    E0 = (402.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51254,0.0417552,1.233e-05,-5.06746e-08,2.43654e-11,48501.7,25.7153], Tmin=(100,'K'), Tmax=(925.185,'K')), NASAPolynomial(coeffs=[14.4236,0.0192295,-5.12952e-06,8.03484e-10,-5.54264e-14,44687.7,-43.2699], Tmin=(925.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1[CH]C(=O)C1[CH2](9977)',
    structure = SMILES('[CH2]C1[CH]C(=O)C1[CH2]'),
    E0 = (399.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51482,0.0396604,2.27954e-05,-6.04507e-08,2.68339e-11,48143.4,25.5807], Tmin=(100,'K'), Tmax=(952.694,'K')), NASAPolynomial(coeffs=[14.7294,0.0213268,-6.83006e-06,1.21083e-09,-8.77144e-14,43939.6,-46.3798], Tmin=(952.694,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O) + radical(Isobutyl) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=C)=C[C]=O(9978)',
    structure = SMILES('[CH2]C=C([CH2])C=C=O'),
    E0 = (226.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2120,512.5,787.5,280.801],'cm^-1')),
        HinderedRotor(inertia=(0.0021379,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.931075,'amu*angstrom^2'), symmetry=1, barrier=(52.0977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.931087,'amu*angstrom^2'), symmetry=1, barrier=(52.0976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29052,0.056049,-4.44762e-05,1.89338e-08,-3.32721e-12,27336,25.1916], Tmin=(100,'K'), Tmax=(1333.84,'K')), NASAPolynomial(coeffs=[11.4265,0.0256529,-1.02937e-05,1.84913e-09,-1.2507e-13,24632,-26.6283], Tmin=(1333.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(C=C(CJ)C=C=O) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=C[C]=O(9979)',
    structure = SMILES('C=C[CH]C=C=O'),
    E0 = (103.397,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,2120,512.5,787.5,292.373,293.234],'cm^-1')),
        HinderedRotor(inertia=(0.48653,'amu*angstrom^2'), symmetry=1, barrier=(29.5014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.48292,'amu*angstrom^2'), symmetry=1, barrier=(29.5,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78547,0.0413552,-1.91969e-05,-5.0012e-09,4.91136e-12,12521.9,18.5543], Tmin=(100,'K'), Tmax=(1025,'K')), NASAPolynomial(coeffs=[11.5734,0.0180947,-7.0149e-06,1.29196e-09,-9.09464e-14,9730.73,-32.7359], Tmin=(1025,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=CCJC=C=O)"""),
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
    label = '[CH2][C]([CH][C]=O)C=C(9980)',
    structure = SMILES('[CH2]C=C([CH2])[CH][C]=O'),
    E0 = (456.868,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,462.808],'cm^-1')),
        HinderedRotor(inertia=(0.000787065,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220497,'amu*angstrom^2'), symmetry=1, barrier=(33.5137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.220498,'amu*angstrom^2'), symmetry=1, barrier=(33.5139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22049,'amu*angstrom^2'), symmetry=1, barrier=(33.5138,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12051,0.0474388,2.31123e-06,-4.07993e-08,1.92158e-11,55066.2,26.3309], Tmin=(100,'K'), Tmax=(1011.16,'K')), NASAPolynomial(coeffs=[17.96,0.0174688,-7.58992e-06,1.56851e-09,-1.20288e-13,49787.4,-64.3601], Tmin=(1011.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.868,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(Allyl_P) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([C]=C)[CH][C]=O(9981)',
    structure = SMILES('[CH2]C([C]=C)[CH][C]=O'),
    E0 = (615.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,481.102,1638.17],'cm^-1')),
        HinderedRotor(inertia=(0.0445487,'amu*angstrom^2'), symmetry=1, barrier=(7.31983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0449476,'amu*angstrom^2'), symmetry=1, barrier=(7.30331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0445294,'amu*angstrom^2'), symmetry=1, barrier=(7.31231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.413001,'amu*angstrom^2'), symmetry=1, barrier=(67.4536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13988,0.0621042,-6.61971e-05,3.99119e-08,-9.74294e-12,74138.5,29.4551], Tmin=(100,'K'), Tmax=(995.703,'K')), NASAPolynomial(coeffs=[10.6197,0.0240206,-8.82404e-06,1.49743e-09,-9.76888e-14,72250.7,-16.2385], Tmin=(995.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(615.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Isobutyl) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CC([CH2])[CH][C]=O(9982)',
    structure = SMILES('[CH]=CC([CH2])[CH][C]=O'),
    E0 = (624.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,3120,650,792.5,1650,350.924],'cm^-1')),
        HinderedRotor(inertia=(0.088817,'amu*angstrom^2'), symmetry=1, barrier=(7.10444,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000163,'amu*angstrom^2'), symmetry=1, barrier=(0.120051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0890014,'amu*angstrom^2'), symmetry=1, barrier=(7.03873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.825322,'amu*angstrom^2'), symmetry=1, barrier=(67.0541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.977476,0.0624918,-6.33101e-05,3.49601e-08,-7.69708e-12,75260.1,29.965], Tmin=(100,'K'), Tmax=(1108.81,'K')), NASAPolynomial(coeffs=[12.5021,0.0209171,-7.06768e-06,1.14456e-09,-7.27737e-14,72704.4,-26.8247], Tmin=(1108.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(624.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Isobutyl) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C](C=C)C[C]=O(5307)',
    structure = SMILES('[CH2]C=C([CH2])C[C]=O'),
    E0 = (294.74,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,453.35],'cm^-1')),
        HinderedRotor(inertia=(0.168961,'amu*angstrom^2'), symmetry=1, barrier=(24.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000807629,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167628,'amu*angstrom^2'), symmetry=1, barrier=(24.4438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164811,'amu*angstrom^2'), symmetry=1, barrier=(24.4416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.983791,0.0549561,-2.35716e-05,-8.13187e-09,6.27952e-12,35567.2,26.4417], Tmin=(100,'K'), Tmax=(1099.76,'K')), NASAPolynomial(coeffs=[15.591,0.0236658,-1.06797e-05,2.10917e-09,-1.53017e-13,31033.7,-51.4227], Tmin=(1099.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.74,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C[C](C)[CH][C]=O(9983)',
    structure = SMILES('C=C[C](C)[CH][C]=O'),
    E0 = (304.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34258,0.0514931,-2.72118e-05,2.73287e-10,3.40644e-12,36738.7,24.7245], Tmin=(100,'K'), Tmax=(1018,'K')), NASAPolynomial(coeffs=[11.1854,0.0262021,-9.66737e-06,1.69886e-09,-1.15331e-13,34041.2,-26.3433], Tmin=(1018,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.617,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Allyl_T) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([C]=C)C[C]=O(5309)',
    structure = SMILES('[CH2]C([C]=C)C[C]=O'),
    E0 = (448.031,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,303.082,303.094],'cm^-1')),
        HinderedRotor(inertia=(0.0018351,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0910627,'amu*angstrom^2'), symmetry=1, barrier=(5.93722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.091079,'amu*angstrom^2'), symmetry=1, barrier=(5.93738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0910517,'amu*angstrom^2'), symmetry=1, barrier=(5.93697,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06792,0.0673097,-7.9288e-05,5.62682e-08,-1.65742e-11,53988.7,29.957], Tmin=(100,'K'), Tmax=(821.713,'K')), NASAPolynomial(coeffs=[8.55184,0.0308768,-1.27776e-05,2.30438e-09,-1.55232e-13,52758.8,-4.67863], Tmin=(821.713,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(448.031,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C(C)[CH][C]=O(9984)',
    structure = SMILES('C=[C]C(C)[CH][C]=O'),
    E0 = (410.478,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,306.149,307.839],'cm^-1')),
        HinderedRotor(inertia=(0.09463,'amu*angstrom^2'), symmetry=1, barrier=(6.3952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00175608,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0017969,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322315,'amu*angstrom^2'), symmetry=1, barrier=(21.7966,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17824,0.059762,-5.13926e-05,2.43467e-08,-4.76423e-12,49472.6,27.5349], Tmin=(100,'K'), Tmax=(1210.06,'K')), NASAPolynomial(coeffs=[11.1148,0.0269156,-1.06761e-05,1.91453e-09,-1.29724e-13,47067.8,-22.2979], Tmin=(1210.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]([CH]C=O)C=C(9985)',
    structure = SMILES('[CH2]C=C([CH2])C=C[O]'),
    E0 = (222.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.846672,0.0498815,1.71286e-05,-7.18134e-08,3.53048e-11,26895.7,24.8286], Tmin=(100,'K'), Tmax=(924.862,'K')), NASAPolynomial(coeffs=[21.4789,0.0108943,-1.13216e-06,9.16872e-11,-1.05856e-14,20930.4,-84.716], Tmin=(924.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C([C]=C)[CH]C=O(9986)',
    structure = SMILES('[CH2]C([C]=C)C=C[O]'),
    E0 = (430.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,323.449,323.449,323.449],'cm^-1')),
        HinderedRotor(inertia=(0.185161,'amu*angstrom^2'), symmetry=1, barrier=(13.7464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347412,'amu*angstrom^2'), symmetry=1, barrier=(25.792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185162,'amu*angstrom^2'), symmetry=1, barrier=(13.7464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.663462,0.0675076,-6.0486e-05,2.52748e-08,-3.05743e-12,51935.5,27.9098], Tmin=(100,'K'), Tmax=(970.674,'K')), NASAPolynomial(coeffs=[14.7366,0.0212902,-7.26319e-06,1.21925e-09,-8.08429e-14,48648.7,-42.4238], Tmin=(970.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C[C]=O(5310)',
    structure = SMILES('[CH]=CC([CH2])C[C]=O'),
    E0 = (457.286,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,3120,650,792.5,1650,298.412],'cm^-1')),
        HinderedRotor(inertia=(0.162686,'amu*angstrom^2'), symmetry=1, barrier=(10.2738,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162169,'amu*angstrom^2'), symmetry=1, barrier=(10.2668,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162769,'amu*angstrom^2'), symmetry=1, barrier=(10.2692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0112125,'amu*angstrom^2'), symmetry=1, barrier=(76.5607,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05006,0.0659547,-7.01265e-05,4.29e-08,-1.07991e-11,55104,29.9509], Tmin=(100,'K'), Tmax=(958.085,'K')), NASAPolynomial(coeffs=[10.0328,0.0284516,-1.14103e-05,2.04303e-09,-1.37881e-13,53382.8,-13.0008], Tmin=(958.085,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.286,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC(C)[CH][C]=O(9987)',
    structure = SMILES('[CH]=CC(C)[CH][C]=O'),
    E0 = (419.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,3120,650,792.5,1650,366.985],'cm^-1')),
        HinderedRotor(inertia=(0.086055,'amu*angstrom^2'), symmetry=1, barrier=(8.21253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0861469,'amu*angstrom^2'), symmetry=1, barrier=(8.21971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00125512,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281637,'amu*angstrom^2'), symmetry=1, barrier=(26.8508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.914305,0.0612669,-5.20521e-05,2.3511e-08,-4.28047e-12,50598.8,28.4145], Tmin=(100,'K'), Tmax=(1317.64,'K')), NASAPolynomial(coeffs=[13.5254,0.022983,-8.46962e-06,1.46013e-09,-9.66833e-14,47275.4,-35.9053], Tmin=(1317.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CC([CH2])[CH]C=O(9988)',
    structure = SMILES('[CH]=CC([CH2])C=C[O]'),
    E0 = (440.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,3120,650,792.5,1650,219.725,219.738],'cm^-1')),
        HinderedRotor(inertia=(0.47454,'amu*angstrom^2'), symmetry=1, barrier=(16.2592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.474553,'amu*angstrom^2'), symmetry=1, barrier=(16.2594,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.788357,'amu*angstrom^2'), symmetry=1, barrier=(27.0135,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.216587,0.0711368,-6.84218e-05,3.3655e-08,-6.40662e-12,53069.7,29.4479], Tmin=(100,'K'), Tmax=(1385.39,'K')), NASAPolynomial(coeffs=[17.5797,0.0166066,-4.61859e-06,6.60648e-10,-3.91238e-14,48680.8,-58.4561], Tmin=(1385.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(C=C)C=C=O(5306)',
    structure = SMILES('[CH2]C(C=C)C=C=O'),
    E0 = (179.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.582259,'amu*angstrom^2'), symmetry=1, barrier=(13.3873,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.581618,'amu*angstrom^2'), symmetry=1, barrier=(13.3725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.581965,'amu*angstrom^2'), symmetry=1, barrier=(13.3805,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898514,0.0708851,-7.87956e-05,4.90713e-08,-1.25483e-11,21671.2,25.0677], Tmin=(100,'K'), Tmax=(942.233,'K')), NASAPolynomial(coeffs=[10.6572,0.0294572,-1.28436e-05,2.40749e-09,-1.67056e-13,19832.3,-21.4314], Tmin=(942.233,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCd(CCO)H) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O)"""),
)

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
    label = 'C=CC1C=C([O])C1(9989)',
    structure = SMILES('C=CC1[CH]C(=O)C1'),
    E0 = (118.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07338,0.0220468,7.16335e-05,-1.07791e-07,4.22444e-11,14333.8,22.7831], Tmin=(100,'K'), Tmax=(966.974,'K')), NASAPolynomial(coeffs=[14.3586,0.022232,-7.77305e-06,1.50271e-09,-1.14843e-13,9573.39,-48.4031], Tmin=(966.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone) + radical(CCJC=O)"""),
)

species(
    label = 'C=CC(=C)C=C[O](9990)',
    structure = SMILES('[CH2]C(C=C)=CC=O'),
    E0 = (84.3418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.687475,0.0630829,-4.71674e-05,1.54348e-08,-1.51128e-12,10271.3,23.8557], Tmin=(100,'K'), Tmax=(1227.37,'K')), NASAPolynomial(coeffs=[16.4738,0.021936,-9.46969e-06,1.79659e-09,-1.26139e-13,5620.32,-58.6993], Tmin=(1227.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.3418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])[C]=[C][O](9991)',
    structure = SMILES('[CH2]CC([CH2])[C][C]=O'),
    E0 = (738.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.749271,0.0749362,-9.70031e-05,7.28256e-08,-2.2215e-11,88892.3,31.8082], Tmin=(100,'K'), Tmax=(840.646,'K')), NASAPolynomial(coeffs=[9.61166,0.0302019,-1.26051e-05,2.265e-09,-1.5155e-13,87492.9,-8.87031], Tmin=(840.646,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([C]=[C][O])[CH]C(9992)',
    structure = SMILES('[CH2]C([C][C]=O)[CH]C'),
    E0 = (727.442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.908236,0.0729429,-9.79424e-05,7.90518e-08,-2.5975e-11,87597.7,32.0061], Tmin=(100,'K'), Tmax=(837.771,'K')), NASAPolynomial(coeffs=[7.69459,0.0334748,-1.46244e-05,2.6827e-09,-1.81351e-13,86708.5,1.94712], Tmin=(837.771,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1[CH]CO[C]=C1(9993)',
    structure = SMILES('[CH2]C1[CH]CO[C]=C1'),
    E0 = (458.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04526,0.0239336,6.75221e-05,-1.12136e-07,4.79351e-11,55243.9,23.9242], Tmin=(100,'K'), Tmax=(904.165,'K')), NASAPolynomial(coeffs=[15.229,0.0156248,-1.66967e-06,6.12916e-11,-3.65804e-15,50815.4,-49.6574], Tmin=(904.165,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(CCJCO) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1C=[C]OC1[CH2](9994)',
    structure = SMILES('[CH2]C1C=[C]OC1[CH2]'),
    E0 = (458.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35036,0.0438349,1.62141e-05,-6.59687e-08,3.39562e-11,55270.5,23.03], Tmin=(100,'K'), Tmax=(876.935,'K')), NASAPolynomial(coeffs=[17.2976,0.0121574,1.59863e-07,-3.66805e-10,3.1495e-14,50894.7,-60.8142], Tmin=(876.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(Isobutyl) + radical(CJC(C)OC) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(C#C[O])C=C(9995)',
    structure = SMILES('[CH2]C([C]=C=O)C=C'),
    E0 = (385.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.626945,'amu*angstrom^2'), symmetry=1, barrier=(14.4147,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.62767,'amu*angstrom^2'), symmetry=1, barrier=(14.4314,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.02055,'amu*angstrom^2'), symmetry=1, barrier=(23.4644,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.926056,0.0735462,-9.20031e-05,5.89555e-08,-1.25295e-11,46418.1,24.8214], Tmin=(100,'K'), Tmax=(650.788,'K')), NASAPolynomial(coeffs=[9.84556,0.0290224,-1.31186e-05,2.46372e-09,-1.69701e-13,45039,-16.0541], Tmin=(650.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCd(CCO)H) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O) + radical(CC(C)CJ=C=O)"""),
)

species(
    label = '[CH2]C([C]=C[O])C=C(9996)',
    structure = SMILES('[CH2]C([C]=C[O])C=C'),
    E0 = (430.774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,323.449,323.449,323.449],'cm^-1')),
        HinderedRotor(inertia=(0.185161,'amu*angstrom^2'), symmetry=1, barrier=(13.7464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.347412,'amu*angstrom^2'), symmetry=1, barrier=(25.792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185162,'amu*angstrom^2'), symmetry=1, barrier=(13.7464,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.663462,0.0675076,-6.0486e-05,2.52748e-08,-3.05743e-12,51935.5,27.9098], Tmin=(100,'K'), Tmax=(970.674,'K')), NASAPolynomial(coeffs=[14.7366,0.0212902,-7.26319e-06,1.21925e-09,-8.08429e-14,48648.7,-42.4238], Tmin=(970.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = 'C=CC(C)[C]=[C][O](9997)',
    structure = SMILES('C=CC(C)[C][C]=O'),
    E0 = (453.322,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,288.98,289.036,289.071],'cm^-1')),
        HinderedRotor(inertia=(0.00201806,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203706,'amu*angstrom^2'), symmetry=1, barrier=(12.0635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339522,'amu*angstrom^2'), symmetry=1, barrier=(20.114,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203684,'amu*angstrom^2'), symmetry=1, barrier=(12.0639,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14321,0.062784,-5.70781e-05,2.8431e-08,-5.86876e-12,54625,27.3343], Tmin=(100,'K'), Tmax=(1145.79,'K')), NASAPolynomial(coeffs=[10.9903,0.0284074,-1.20745e-05,2.24621e-09,-1.55523e-13,52368.4,-21.5126], Tmin=(1145.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.322,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([C]=[C]O)C=C(9998)',
    structure = SMILES('[CH2]C([C]=[C]O)C=C'),
    E0 = (529.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,305.28,307.112],'cm^-1')),
        HinderedRotor(inertia=(0.178706,'amu*angstrom^2'), symmetry=1, barrier=(11.8629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178228,'amu*angstrom^2'), symmetry=1, barrier=(11.8624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177245,'amu*angstrom^2'), symmetry=1, barrier=(11.8605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241504,'amu*angstrom^2'), symmetry=1, barrier=(16.0365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.543127,0.0752182,-8.69295e-05,5.24275e-08,-1.20188e-11,63756,30.0643], Tmin=(100,'K'), Tmax=(867.281,'K')), NASAPolynomial(coeffs=[13.6573,0.0220677,-7.68686e-06,1.26434e-09,-8.10159e-14,61205.5,-32.9265], Tmin=(867.281,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](C=C)C=[C]O(9999)',
    structure = SMILES('[CH2]C=C([CH2])C=[C]O'),
    E0 = (320.811,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.646821,0.058579,-1.29983e-05,-3.95448e-08,2.40055e-11,38719.6,27.265], Tmin=(100,'K'), Tmax=(913.908,'K')), NASAPolynomial(coeffs=[20.6552,0.0112274,-1.29553e-06,7.4594e-11,-5.55629e-15,33382.7,-76.6524], Tmin=(913.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.811,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(Allyl_P) + radical(C=CC=CCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C([C]=C)C=[C]O(10000)',
    structure = SMILES('[CH2]C([C]=C)C=[C]O'),
    E0 = (529.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,305.28,307.112],'cm^-1')),
        HinderedRotor(inertia=(0.178706,'amu*angstrom^2'), symmetry=1, barrier=(11.8629,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178228,'amu*angstrom^2'), symmetry=1, barrier=(11.8624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177245,'amu*angstrom^2'), symmetry=1, barrier=(11.8605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.241504,'amu*angstrom^2'), symmetry=1, barrier=(16.0365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.543127,0.0752182,-8.69295e-05,5.24275e-08,-1.20188e-11,63756,30.0643], Tmin=(100,'K'), Tmax=(867.281,'K')), NASAPolynomial(coeffs=[13.6573,0.0220677,-7.68686e-06,1.26434e-09,-8.10159e-14,61205.5,-32.9265], Tmin=(867.281,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=CC([CH2])C=[C]O(10001)',
    structure = SMILES('[CH]=CC([CH2])C=[C]O'),
    E0 = (538.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,215.75],'cm^-1')),
        HinderedRotor(inertia=(0.406802,'amu*angstrom^2'), symmetry=1, barrier=(13.437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.406808,'amu*angstrom^2'), symmetry=1, barrier=(13.4371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.406812,'amu*angstrom^2'), symmetry=1, barrier=(13.4371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506409,'amu*angstrom^2'), symmetry=1, barrier=(16.7274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.474045,0.074491,-8.00932e-05,4.2312e-08,-7.76378e-12,64873.6,30.2405], Tmin=(100,'K'), Tmax=(892.235,'K')), NASAPolynomial(coeffs=[15.1808,0.0195718,-6.27977e-06,9.93781e-10,-6.29147e-14,61810.8,-41.4909], Tmin=(892.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(C=CJO) + radical(Cds_P)"""),
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
    E0 = (377.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (573.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (622.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (755.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (782.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (870.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (832.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (383.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (400.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (400.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (665.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (517.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (559.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (446.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (433.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (414.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (460.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (446.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (524.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (537.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (494.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (789.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (668.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (827.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (836.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (535.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (519.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (590.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (562.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (530.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (580.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (501.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (464.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (473.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (377.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (537.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (386.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (441.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (801.546,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (752.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (458.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (476.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (611.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (494.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (634.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (605.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (679.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (500.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (680.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (689.831,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['HCCO(2227)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([C]=O)C=C(9966)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=CC[CH][CH][C]=O(9967)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', '[CH2]C=C[CH][C]=O(9826)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2]C([C][C]=O)C=C(9968)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH]C([CH][C]=O)C=C(9969)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['C=CC1CC1[C]=O(9970)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['C=CC(=C)C[C]=O(5305)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['C=CC(C)=C[C]=O(9971)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C[C]([CH2])[CH][C]=O(9972)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C1[CH]CC1[C]=O(9973)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.61918e+08,'s^-1'), n=0.930343, Ea=(139.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['O=[C][CH]C1[CH]CC1(6406)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.18363e+11,'s^-1'), n=0.209288, Ea=(182.141,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 125 used for R4_S_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C1[CH]CC(=O)[CH]1(9974)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.61878e+08,'s^-1'), n=0.712108, Ea=(68.5201,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_linear;doublebond_intra_pri_2H;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C1C([CH2])C1[C]=O(9975)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.65211e+08,'s^-1'), n=0.77, Ea=(55.5948,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C1CC1[CH][C]=O(9976)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C1[CH]C(=O)C1[CH2](9977)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.89869e+09,'s^-1'), n=0.700948, Ea=(83.1192,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra_2H_pri;radadd_intra] for rate rule [R5;doublebond_intra_2H_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2]C(C=C)=C[C]=O(9978)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CH2(T)(20)', 'C=CC=C[C]=O(9979)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C][O](6861)', 'butadiene13(1350)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.10641,'m^3/(mol*s)'), n=2.065, Ea=(15.9938,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds-HH;YJ] for rate rule [Cds-CdH_Cds-HH;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C2H3(60)', 'C=CC=[C][O](9589)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0131003,'m^3/(mol*s)'), n=2.40999, Ea=(12.7705,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C2H3(60)', '[CH2][CH][CH][C]=O(9637)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.30526e+07,'m^3/(mol*s)'), n=-0.283333, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-6.50053689359e-11, var=0.305422193575, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C
    Total Standard Deviation in ln(k): 1.10791715097
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][C]([CH][C]=O)C=C(9980)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C([C]=C)[CH][C]=O(9981)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH]=CC([CH2])[CH][C]=O(9982)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2][C](C=C)C[C]=O(5307)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['C=C[C](C)[CH][C]=O(9983)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([C]=C)C[C]=O(5309)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]C(C)[CH][C]=O(9984)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2][C]([CH]C=O)C=C(9985)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.10706e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([C]=C)[CH]C=O(9986)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.13341e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_Cd;XH_out] for rate rule [R4HJ_2;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=CC([CH2])C[C]=O(5310)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=CC(C)[CH][C]=O(9987)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=CC([CH2])[CH]C=O(9988)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C(C=C)C=C=O(5306)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['C=CC1C=C([O])C1(9989)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['C=CC(=C)C=C[O](9990)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]CC([CH2])[C]=[C][O](9991)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([C]=[C][O])[CH]C(9992)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C1[CH]CO[C]=C1(9993)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(80.8712,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 76.6 to 80.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C1C=[C]OC1[CH2](9994)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(6.82166e+09,'s^-1'), n=0.527281, Ea=(99.1325,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H_pri;radadd_intra_O] + [R6;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(3)', '[CH2]C(C#C[O])C=C(9995)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['HCCO(2227)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(2.08e-06,'m^3/(mol*s)'), n=3.05, Ea=(54.8104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CdsJ=Cdd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([C]=C[O])C=C(9996)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=CC(C)[C]=[C][O](9997)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C([C]=[C]O)C=C(9998)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2][C](C=C)C=[C]O(9999)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(520772,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C([C]=C)C=[C]O(10000)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(4.81182e+08,'s^-1'), n=1.25566, Ea=(151.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cd;XH_out] for rate rule [R5HJ_3;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=CC([CH2])C=[C]O(10001)'],
    products = ['[CH2]C([CH][C]=O)C=C(5308)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_4;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2472',
    isomers = [
        '[CH2]C([CH][C]=O)C=C(5308)',
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
    label = 'PDepNetwork #2472',
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


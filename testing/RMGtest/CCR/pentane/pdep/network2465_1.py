species(
    label = '[CH]=[C]OC([CH2])C=C(5360)',
    structure = SMILES('[CH]=[C]OC([CH2])C=C'),
    E0 = (579.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,270.496,270.898,271.049],'cm^-1')),
        HinderedRotor(inertia=(0.332271,'amu*angstrom^2'), symmetry=1, barrier=(17.2435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332851,'amu*angstrom^2'), symmetry=1, barrier=(17.2473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332778,'amu*angstrom^2'), symmetry=1, barrier=(17.244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332976,'amu*angstrom^2'), symmetry=1, barrier=(17.2483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.153209,0.0784731,-8.11336e-05,4.09389e-08,-7.90232e-12,69808.4,32.4073], Tmin=(100,'K'), Tmax=(1347.35,'K')), NASAPolynomial(coeffs=[20.8927,0.0123646,-3.4967e-06,5.25983e-10,-3.29628e-14,64466.5,-74.1798], Tmin=(1347.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(579.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]O[CH]CC=C(10046)',
    structure = SMILES('[CH]=[C]O[CH]CC=C'),
    E0 = (575.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.110381,0.0792971,-8.37804e-05,4.34564e-08,-8.6559e-12,69372,31.4085], Tmin=(100,'K'), Tmax=(1275.08,'K')), NASAPolynomial(coeffs=[20.3705,0.0132469,-3.96093e-06,6.15793e-10,-3.91685e-14,64295.4,-71.8028], Tmin=(1275.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[C]=[CH](9646)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([554.803,1738.79,3454.47],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97621,0.00212915,-8.08978e-08,-3.83305e-10,9.76908e-14,101881,6.00119], Tmin=(100,'K'), Tmax=(1982.31,'K')), NASAPolynomial(coeffs=[5.05695,0.00131032,-4.91873e-07,1.01502e-10,-7.16167e-15,101185,-0.627232], Tmin=(1982.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])C=C(1523)',
    structure = SMILES('[CH2]C([O])C=C'),
    E0 = (252.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,529.328,529.355],'cm^-1')),
        HinderedRotor(inertia=(0.0692523,'amu*angstrom^2'), symmetry=1, barrier=(13.7797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0692738,'amu*angstrom^2'), symmetry=1, barrier=(13.7796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95358,0.036167,-6.30157e-06,-1.92918e-08,1.03178e-11,30464.9,21.1322], Tmin=(100,'K'), Tmax=(988.725,'K')), NASAPolynomial(coeffs=[12.0826,0.0154802,-5.70142e-06,1.06015e-09,-7.65757e-14,27470.1,-32.6352], Tmin=(988.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = '[CH]=[C]OC=C[CH2](9859)',
    structure = SMILES('[CH]=[C]O[CH]C=C'),
    E0 = (521.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,751.257,751.599],'cm^-1')),
        HinderedRotor(inertia=(0.769937,'amu*angstrom^2'), symmetry=1, barrier=(17.7024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0962829,'amu*angstrom^2'), symmetry=1, barrier=(38.6301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.769836,'amu*angstrom^2'), symmetry=1, barrier=(17.7,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53123,0.0477042,-3.52209e-05,9.58444e-09,2.4435e-13,62856.4,25.1157], Tmin=(100,'K'), Tmax=(1055.07,'K')), NASAPolynomial(coeffs=[12.6255,0.0165414,-6.41032e-06,1.16984e-09,-8.15187e-14,59908.8,-31.8771], Tmin=(1055.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]C(C=C)O[C]=[CH](10047)',
    structure = SMILES('[CH]C(C=C)O[C]=[CH]'),
    E0 = (816.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.536157,0.0642651,-4.19037e-05,-3.15593e-09,8.8876e-12,98372.8,31.1792], Tmin=(100,'K'), Tmax=(958.951,'K')), NASAPolynomial(coeffs=[19.8728,0.0116503,-3.46804e-06,6.18437e-10,-4.64663e-14,93374.8,-68.0215], Tmin=(958.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(816.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCJ2_triplet) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[C]=[C]OC([CH2])C=C(10048)',
    structure = SMILES('[C]=[C]OC([CH2])C=C'),
    E0 = (890.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,340.039,340.041,340.041,340.041],'cm^-1')),
        HinderedRotor(inertia=(0.193946,'amu*angstrom^2'), symmetry=1, barrier=(15.9134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193943,'amu*angstrom^2'), symmetry=1, barrier=(15.9134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19394,'amu*angstrom^2'), symmetry=1, barrier=(15.9134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193943,'amu*angstrom^2'), symmetry=1, barrier=(15.9134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259021,0.0769851,-8.68513e-05,4.82655e-08,-1.03831e-11,107193,30.8977], Tmin=(100,'K'), Tmax=(1145.45,'K')), NASAPolynomial(coeffs=[18.0154,0.014979,-5.65314e-06,1.00752e-09,-6.89383e-14,103125,-57.178], Tmin=(1145.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(890.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(C=CJO) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C1CC(C=C)O1(10049)',
    structure = SMILES('[CH]=C1CC(C=C)O1'),
    E0 = (259.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30301,0.0352865,5.78179e-05,-1.16883e-07,5.26719e-11,31282.2,21.2205], Tmin=(100,'K'), Tmax=(911.922,'K')), NASAPolynomial(coeffs=[22.3383,0.00599984,2.3944e-06,-6.30451e-10,3.92556e-14,24826.9,-92.6817], Tmin=(911.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Cds_P)"""),
)

species(
    label = '[CH]=COC(=C)C=C(10050)',
    structure = SMILES('[CH]=COC(=C)C=C'),
    E0 = (263.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605453,0.0650352,-4.59651e-05,6.20191e-09,4.25506e-12,31832.4,25.6455], Tmin=(100,'K'), Tmax=(976.469,'K')), NASAPolynomial(coeffs=[17.1696,0.0181789,-6.24145e-06,1.10252e-09,-7.73563e-14,27596.5,-58.9984], Tmin=(976.469,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.582,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]O[C]([CH2])C[CH2](10051)',
    structure = SMILES('[CH]=[C]O[C]([CH2])C[CH2]'),
    E0 = (857.102,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3120,650,792.5,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0106316,0.0920762,-0.000131253,9.68937e-08,-2.81201e-11,103226,33.7978], Tmin=(100,'K'), Tmax=(848.452,'K')), NASAPolynomial(coeffs=[14.4862,0.0237218,-1.03906e-05,1.91366e-09,-1.29941e-13,100767,-33.7565], Tmin=(848.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(857.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C2CsJOC(O)) + radical(RCCJ) + radical(CJC(C)OC) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC1[CH]CC1(6425)',
    structure = SMILES('[CH]=[C]OC1[CH]CC1'),
    E0 = (599.209,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,1685,370,3120,650,792.5,1650,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35738,0.0365789,4.56408e-05,-9.46495e-08,4.13059e-11,72183.1,27.9087], Tmin=(100,'K'), Tmax=(943.979,'K')), NASAPolynomial(coeffs=[19.8247,0.0126713,-2.72533e-06,4.95252e-10,-4.35331e-14,66275.2,-72.9458], Tmin=(943.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(599.209,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCJCO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C[CH]C([CH2])O1(10052)',
    structure = SMILES('[CH]=C1C[CH]C([CH2])O1'),
    E0 = (471.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30304,0.0356884,5.37783e-05,-1.08864e-07,4.81921e-11,56807.1,25.3478], Tmin=(100,'K'), Tmax=(928.838,'K')), NASAPolynomial(coeffs=[21.7693,0.00820421,2.1322e-07,-1.15581e-10,-1.47576e-16,50388.7,-85.9637], Tmin=(928.838,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCJCO) + radical(CJC(C)OC) + radical(Cds_P)"""),
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
    label = '[CH]=[C]OC1CC1[CH2](10053)',
    structure = SMILES('[CH]=[C]OC1CC1[CH2]'),
    E0 = (599.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953007,0.0493903,1.05466e-05,-6.24826e-08,3.17037e-11,72261.4,27.6595], Tmin=(100,'K'), Tmax=(920.201,'K')), NASAPolynomial(coeffs=[20.4451,0.0105506,-9.46116e-07,3.8044e-11,-5.38544e-15,66731.1,-75.3148], Tmin=(920.201,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(599.767,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1OC([CH2])C1[CH2](10054)',
    structure = SMILES('[CH]=C1OC([CH2])C1[CH2]'),
    E0 = (539.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.09623,0.0777981,-7.63913e-05,3.66236e-08,-6.36852e-12,65073.2,29.0748], Tmin=(100,'K'), Tmax=(1744,'K')), NASAPolynomial(coeffs=[18.2386,0.0104973,1.23689e-06,-5.97956e-10,4.90069e-14,61820.1,-64.9496], Tmin=(1744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(539.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Isobutyl) + radical(CJC(C)OC) + radical(Cds_P)"""),
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
    label = '[CH]=[C]OC(=C)C=C(10055)',
    structure = SMILES('C#COC([CH2])=C[CH2]'),
    E0 = (411.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2175,525,750,770,3400,2100,261.437,261.624],'cm^-1')),
        HinderedRotor(inertia=(0.00245885,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.635993,'amu*angstrom^2'), symmetry=1, barrier=(30.9299,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.633361,'amu*angstrom^2'), symmetry=1, barrier=(30.9334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.633932,'amu*angstrom^2'), symmetry=1, barrier=(30.9287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684429,0.0678694,-6.75806e-05,3.44492e-08,-6.93968e-12,49644.4,25.1782], Tmin=(100,'K'), Tmax=(1207.93,'K')), NASAPolynomial(coeffs=[15.2599,0.0196032,-7.64362e-06,1.36927e-09,-9.3238e-14,46123.2,-47.8935], Tmin=(1207.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=C(O)CJ) + radical(Allyl_P)"""),
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
    label = 'C#CO[CH][CH2](3100)',
    structure = SMILES('C#CO[CH][CH2]'),
    E0 = (415.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19325,'amu*angstrom^2'), symmetry=1, barrier=(27.4351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19703,'amu*angstrom^2'), symmetry=1, barrier=(27.522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20177,'amu*angstrom^2'), symmetry=1, barrier=(27.631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.563675,0.0617773,-7.43297e-05,4.03021e-08,-7.97839e-12,50152.6,19.2452], Tmin=(100,'K'), Tmax=(1440.47,'K')), NASAPolynomial(coeffs=[19.2396,-0.00109988,2.61778e-06,-6.19395e-10,4.5176e-14,45915.1,-73.7043], Tmin=(1440.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCsJOCs) + radical(CJCO)"""),
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
    label = '[CH]=[C]O[CH][CH2](9651)',
    structure = SMILES('[CH]=[C]O[CH][CH2]'),
    E0 = (706.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,203.909,203.912],'cm^-1')),
        HinderedRotor(inertia=(0.602174,'amu*angstrom^2'), symmetry=1, barrier=(17.7676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602165,'amu*angstrom^2'), symmetry=1, barrier=(17.7676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602172,'amu*angstrom^2'), symmetry=1, barrier=(17.7676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.946535,0.0601409,-7.68027e-05,4.53754e-08,-9.98375e-12,85130.4,23.972], Tmin=(100,'K'), Tmax=(1232.6,'K')), NASAPolynomial(coeffs=[17.1138,0.00205631,7.20679e-07,-2.52539e-10,2.07817e-14,81571.7,-55.6752], Tmin=(1232.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(706.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(CJCO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]O[C]([CH2])C=C(10056)',
    structure = SMILES('[CH]=[C]OC([CH2])=C[CH2]'),
    E0 = (689.253,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.153615,'amu*angstrom^2'), symmetry=1, barrier=(16.2489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00562003,'amu*angstrom^2'), symmetry=1, barrier=(16.2504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0122803,'amu*angstrom^2'), symmetry=1, barrier=(16.2455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.25355,'amu*angstrom^2'), symmetry=1, barrier=(74.8055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.985437,0.0655472,-6.80989e-05,3.79369e-08,-8.51768e-12,83007.2,30.1749], Tmin=(100,'K'), Tmax=(1077.04,'K')), NASAPolynomial(coeffs=[12.2843,0.0235844,-9.65655e-06,1.76195e-09,-1.20789e-13,80573.4,-25.1739], Tmin=(1077.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.253,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC([CH2])[C]=C(10057)',
    structure = SMILES('[CH]=[C]OC([CH2])[C]=C'),
    E0 = (816.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,282.509,282.514,282.521],'cm^-1')),
        HinderedRotor(inertia=(0.288072,'amu*angstrom^2'), symmetry=1, barrier=(16.3137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288025,'amu*angstrom^2'), symmetry=1, barrier=(16.3137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288016,'amu*angstrom^2'), symmetry=1, barrier=(16.3136,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.287998,'amu*angstrom^2'), symmetry=1, barrier=(16.3136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.176784,0.078888,-9.17375e-05,5.24586e-08,-1.15721e-11,98395.9,31.691], Tmin=(100,'K'), Tmax=(1120.07,'K')), NASAPolynomial(coeffs=[18.2978,0.0141743,-5.07296e-06,8.75966e-10,-5.88932e-14,94336.5,-57.7874], Tmin=(1120.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(816.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC([CH2])C=[CH](10058)',
    structure = SMILES('[CH]=[C]OC([CH2])C=[CH]'),
    E0 = (826.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3115,3125,620,680,785,800,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.762988,'amu*angstrom^2'), symmetry=1, barrier=(17.5426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.763664,'amu*angstrom^2'), symmetry=1, barrier=(17.5581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.763938,'amu*angstrom^2'), symmetry=1, barrier=(17.5644,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.763094,'amu*angstrom^2'), symmetry=1, barrier=(17.545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.133765,0.0778463,-8.37856e-05,4.09099e-08,-6.72424e-12,99512.4,31.7743], Tmin=(100,'K'), Tmax=(980.772,'K')), NASAPolynomial(coeffs=[19.5636,0.0121132,-3.91542e-06,6.64167e-10,-4.56554e-14,95051.3,-64.8989], Tmin=(980.772,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(826.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(C=CJO) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]O[C](C)C=C(10059)',
    structure = SMILES('[CH]=[C]OC(C)=C[CH2]'),
    E0 = (530.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19351,0.0629527,-5.76002e-05,2.91025e-08,-6.14826e-12,63884.7,28.7402], Tmin=(100,'K'), Tmax=(1114.18,'K')), NASAPolynomial(coeffs=[10.4031,0.0298897,-1.30882e-05,2.46902e-09,-1.72253e-13,61832.5,-16.6865], Tmin=(1114.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=COC([CH2])=C[CH2](7153)',
    structure = SMILES('[CH]=COC([CH2])=C[CH2]'),
    E0 = (449.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3120,650,792.5,1650,311.116,314.047],'cm^-1')),
        HinderedRotor(inertia=(0.00175406,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300249,'amu*angstrom^2'), symmetry=1, barrier=(21.5283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317514,'amu*angstrom^2'), symmetry=1, barrier=(21.5056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.401386,'amu*angstrom^2'), symmetry=1, barrier=(27.8361,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.591939,0.0660051,-4.92092e-05,9.89423e-09,2.94168e-12,54194.1,28.455], Tmin=(100,'K'), Tmax=(973.314,'K')), NASAPolynomial(coeffs=[16.8117,0.0187956,-6.4258e-06,1.11934e-09,-7.75137e-14,50115.5,-54.0894], Tmin=(973.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(449.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(C)[C]=C(10060)',
    structure = SMILES('[CH]=[C]OC(C)[C]=C'),
    E0 = (606.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,321.932,322.265,322.908],'cm^-1')),
        HinderedRotor(inertia=(0.227159,'amu*angstrom^2'), symmetry=1, barrier=(16.6855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226199,'amu*angstrom^2'), symmetry=1, barrier=(16.6911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227426,'amu*angstrom^2'), symmetry=1, barrier=(16.69,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226842,'amu*angstrom^2'), symmetry=1, barrier=(16.6889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538158,0.0660139,-4.78839e-05,7.33844e-09,3.92281e-12,73068.9,30.6416], Tmin=(100,'K'), Tmax=(986.873,'K')), NASAPolynomial(coeffs=[17.9311,0.0169856,-5.99519e-06,1.08514e-09,-7.73733e-14,68590.5,-58.3364], Tmin=(986.873,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(606.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C](C=C)O[C]=C(5357)',
    structure = SMILES('[CH2]C=C([CH2])O[C]=C'),
    E0 = (442.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00759795,'amu*angstrom^2'), symmetry=1, barrier=(5.27335,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0213457,'amu*angstrom^2'), symmetry=1, barrier=(14.8257,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0213531,'amu*angstrom^2'), symmetry=1, barrier=(14.8251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.87926,'amu*angstrom^2'), symmetry=1, barrier=(89.1919,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949344,0.0633398,-5.6079e-05,2.65248e-08,-5.0955e-12,53292.1,29.8991], Tmin=(100,'K'), Tmax=(1243.42,'K')), NASAPolynomial(coeffs=[12.8404,0.0250871,-9.93299e-06,1.78335e-09,-1.21042e-13,50335,-30.0591], Tmin=(1243.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC([CH2])[C]=C(10061)',
    structure = SMILES('[CH]=COC([CH2])[C]=C'),
    E0 = (577.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.889318,'amu*angstrom^2'), symmetry=1, barrier=(20.4472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.889649,'amu*angstrom^2'), symmetry=1, barrier=(20.4548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.891212,'amu*angstrom^2'), symmetry=1, barrier=(20.4907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.891575,'amu*angstrom^2'), symmetry=1, barrier=(20.4991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.658154,0.0843748,-8.96475e-05,4.51447e-08,-8.51931e-12,69602.3,31.567], Tmin=(100,'K'), Tmax=(1456.93,'K')), NASAPolynomial(coeffs=[23.7057,0.00785439,-9.50616e-07,2.18808e-11,1.94832e-15,63525,-91.6358], Tmin=(1456.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(C)C=[CH](10062)',
    structure = SMILES('[CH]=[C]OC(C)C=[CH]'),
    E0 = (615.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700,1685,370,228.214,233.505],'cm^-1')),
        HinderedRotor(inertia=(0.460648,'amu*angstrom^2'), symmetry=1, barrier=(17.8669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.487594,'amu*angstrom^2'), symmetry=1, barrier=(17.8653,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.455077,'amu*angstrom^2'), symmetry=1, barrier=(17.8808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.485872,'amu*angstrom^2'), symmetry=1, barrier=(17.872,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.507264,0.0648206,-3.93485e-05,-5.07092e-09,9.194e-12,74184.8,30.6819], Tmin=(100,'K'), Tmax=(965.326,'K')), NASAPolynomial(coeffs=[19.2921,0.0147684,-4.74977e-06,8.52932e-10,-6.24646e-14,69263.5,-65.987], Tmin=(965.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(615.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([C]=C)O[C]=C(5358)',
    structure = SMILES('[CH2]C([C]=C)O[C]=C'),
    E0 = (569.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,340.956,340.956,340.956,340.956],'cm^-1')),
        HinderedRotor(inertia=(0.194194,'amu*angstrom^2'), symmetry=1, barrier=(16.0199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194194,'amu*angstrom^2'), symmetry=1, barrier=(16.0199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194194,'amu*angstrom^2'), symmetry=1, barrier=(16.0199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194194,'amu*angstrom^2'), symmetry=1, barrier=(16.0199,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.108991,0.0770016,-8.0628e-05,4.19904e-08,-8.47411e-12,68682.3,31.5328], Tmin=(100,'K'), Tmax=(1220.43,'K')), NASAPolynomial(coeffs=[18.772,0.0158324,-5.44584e-06,9.2136e-10,-6.12131e-14,64127,-62.223], Tmin=(1220.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=COC([CH2])C=[CH](10063)',
    structure = SMILES('[CH]=COC([CH2])C=[CH]'),
    E0 = (586.438,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,3115,3125,620,680,785,800,1600,1700,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.945758,'amu*angstrom^2'), symmetry=1, barrier=(21.7448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.946317,'amu*angstrom^2'), symmetry=1, barrier=(21.7577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.945691,'amu*angstrom^2'), symmetry=1, barrier=(21.7433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.946045,'amu*angstrom^2'), symmetry=1, barrier=(21.7514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.940294,0.086048,-9.07347e-05,4.47206e-08,-8.17347e-12,70729.4,32.5158], Tmin=(100,'K'), Tmax=(1542.23,'K')), NASAPolynomial(coeffs=[24.8707,0.00579301,2.68607e-07,-2.14038e-10,1.77344e-14,64351.1,-98.0577], Tmin=(1542.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.438,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])O[C]=C(5359)',
    structure = SMILES('[CH]=CC([CH2])O[C]=C'),
    E0 = (579.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,270.496,270.898,271.049],'cm^-1')),
        HinderedRotor(inertia=(0.332271,'amu*angstrom^2'), symmetry=1, barrier=(17.2435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332851,'amu*angstrom^2'), symmetry=1, barrier=(17.2473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332778,'amu*angstrom^2'), symmetry=1, barrier=(17.244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.332976,'amu*angstrom^2'), symmetry=1, barrier=(17.2483,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.153209,0.0784731,-8.11336e-05,4.09389e-08,-7.90232e-12,69808.4,32.4073], Tmin=(100,'K'), Tmax=(1347.35,'K')), NASAPolynomial(coeffs=[20.8927,0.0123646,-3.4967e-06,5.25983e-10,-3.29628e-14,64466.5,-74.1798], Tmin=(1347.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(579.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CJC(C)OC) + radical(C=CJO) + radical(Cds_P)"""),
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
    E0 = (579.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (737.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1099.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (959.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1028.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1101.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (587.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (642.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (879.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (703.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (635.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (591.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (616.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (666.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (592.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (632.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (579.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (711.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (699.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (993.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (901.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1028.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1037.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (706.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (732.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (758.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (788.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (678.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (659.981,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (888.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (619.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (770.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['HCCO(2227)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH]=[C]O[CH]CC=C(10046)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.95289e+09,'s^-1'), n=1, Ea=(158.627,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CsJ-HH;C] + [cCs(-HR!H)CJ;CsJ;C] for rate rule [cCs(-HR!H)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=[CH](9646)', '[CH2]C([O])C=C(1523)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(20)', '[CH]=[C]OC=C[CH2](9859)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH]C(C=C)O[C]=[CH](10047)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[C]=[C]OC([CH2])C=C(10048)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH]=C1CC(C=C)O1(10049)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH]=COC(=C)C=C(10050)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]O[C]([CH2])C[CH2](10051)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH]=[C]OC1[CH]CC1(6425)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH]=C1C[CH]C([CH2])O1(10052)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH2]C1[CH]CC=[C]O1(10032)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.16959e+10,'s^-1'), n=0.31, Ea=(12.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH]=[C]OC1CC1[CH2](10053)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH]=C1OC([CH2])C1[CH2](10054)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.65009e+08,'s^-1'), n=1.00067, Ea=(87.7294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH2]C1C=[C]OC1[CH2](9994)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.4226e+09,'s^-1'), n=0.3735, Ea=(13.1942,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6;doublebond_intra_2H_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH]=[C]OC(=C)C=C(10055)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(87.5179,'m^3/(mol*s)'), n=1.66467, Ea=(9.4546,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-OneDe_Cds;HJ] for rate rule [Cds-CdOs_Cds;HJ]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HCCO(2227)', '[CH2][CH]C=C(3743)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.6e+11,'cm^3/(mol*s)'), n=0, Ea=(139.182,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_rad/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 135.7 to 139.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C2H3(60)', 'C#CO[CH][CH2](3100)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00977588,'m^3/(mol*s)'), n=2.40996, Ea=(8.29121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CdsJ-H] for rate rule [Cds-OsH_Cds;CdsJ-H]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH]C=C(3743)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.206e+07,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C2H3(60)', '[CH]=[C]O[CH][CH2](9651)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.30526e+07,'m^3/(mol*s)'), n=-0.283333, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-6.50053689359e-11, var=0.305422193575, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C
    Total Standard Deviation in ln(k): 1.10791715097
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH]=[C]O[C]([CH2])C=C(10056)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH]=[C]OC([CH2])[C]=C(10057)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.28426e+17,'m^3/(mol*s)'), n=-3.05017, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=3.08884453463, var=70.5675449198, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R
    Total Standard Deviation in ln(k): 24.6015908735
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH]=[C]OC([CH2])C=[CH](10058)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH]=[C]O[C](C)C=C(10059)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.08414e+09,'s^-1'), n=1.12201, Ea=(127.544,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;C_rad_out_single;Cs_H_out_Cd] + [R2H_S;C_rad_out_2H;Cs_H_out_OneDe] for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cd]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH]=COC([CH2])=C[CH2](7153)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.38739e+08,'s^-1'), n=1.16185, Ea=(153.379,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Cd_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]OC(C)[C]=C(10060)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH2][C](C=C)O[C]=C(5357)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(868165,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=COC([CH2])[C]=C(10061)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]OC(C)C=[CH](10062)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C]OC([CH2])C=C(5360)'],
    products = ['[CH2]C([C]=C)O[C]=C(5358)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=COC([CH2])C=[CH](10063)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=CC([CH2])O[C]=C(5359)'],
    products = ['[CH]=[C]OC([CH2])C=C(5360)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R6HJ_4;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2465',
    isomers = [
        '[CH]=[C]OC([CH2])C=C(5360)',
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
    label = 'PDepNetwork #2465',
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


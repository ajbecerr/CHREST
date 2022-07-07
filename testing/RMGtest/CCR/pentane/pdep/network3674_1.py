species(
    label = '[CH2]C([CH2])=C[C]=O(9813)',
    structure = SMILES('[CH2]C([CH2])=C[C]=O'),
    E0 = (307.2,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.12597,'amu*angstrom^2'), symmetry=1, barrier=(25.8882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12599,'amu*angstrom^2'), symmetry=1, barrier=(25.8887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.126,'amu*angstrom^2'), symmetry=1, barrier=(25.889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50488,0.0572345,-5.84403e-05,3.09075e-08,-6.68244e-12,37035.6,20.0483], Tmin=(100,'K'), Tmax=(1099.74,'K')), NASAPolynomial(coeffs=[11.0466,0.022529,-1.11032e-05,2.21137e-09,-1.58996e-13,34936.9,-26.8919], Tmin=(1099.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.2,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CCJ=O)"""),
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
    label = 'allene(458)',
    structure = SMILES('C=C=C'),
    E0 = (175.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37447,0.0070462,2.78306e-05,-3.99445e-08,1.55729e-11,21188.6,7.62048], Tmin=(100,'K'), Tmax=(949.705,'K')), NASAPolynomial(coeffs=[6.79956,0.00959979,-3.02068e-06,5.37827e-10,-3.92606e-14,19772.3,-12.7582], Tmin=(949.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""allene""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][C]=C[C]=O(9897)',
    structure = SMILES('[CH2][C]=C[C]=O'),
    E0 = (432.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.00399,'amu*angstrom^2'), symmetry=1, barrier=(23.0837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00442,'amu*angstrom^2'), symmetry=1, barrier=(23.0936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33986,0.0429161,-4.87447e-05,1.23757e-08,1.41321e-11,52082.9,16.5479], Tmin=(100,'K'), Tmax=(512.282,'K')), NASAPolynomial(coeffs=[6.6704,0.019476,-1.04855e-05,2.11536e-09,-1.51263e-13,51503.1,-2.77645], Tmin=(512.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(Cds_S) + radical(C=CCJ=O)"""),
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
    label = '[CH]=C([CH2])[CH2](16419)',
    structure = SMILES('[CH]C([CH2])=C'),
    E0 = (488.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,412.019,412.344,412.414],'cm^-1')),
        HinderedRotor(inertia=(0.422651,'amu*angstrom^2'), symmetry=1, barrier=(50.9185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.421927,'amu*angstrom^2'), symmetry=1, barrier=(50.9131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54708,0.0237852,1.61126e-05,-3.4475e-08,1.40391e-11,58778.2,14.8721], Tmin=(100,'K'), Tmax=(988.925,'K')), NASAPolynomial(coeffs=[8.49332,0.0194138,-7.10753e-06,1.30164e-09,-9.24087e-14,56639.8,-18.6144], Tmin=(988.925,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([CH2])=C[C]=O(16420)',
    structure = SMILES('[CH]C([CH2])=C[C]=O'),
    E0 = (526.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.10045,'amu*angstrom^2'), symmetry=1, barrier=(48.2934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09767,'amu*angstrom^2'), symmetry=1, barrier=(48.2295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10075,'amu*angstrom^2'), symmetry=1, barrier=(48.3004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52118,0.0593111,-6.47653e-05,4.05501e-08,-1.0837e-11,63394.7,20.8288], Tmin=(100,'K'), Tmax=(884.972,'K')), NASAPolynomial(coeffs=[8.03996,0.0298466,-1.48236e-05,2.92786e-09,-2.08896e-13,62241,-9.82393], Tmin=(884.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = 'O=[C]C=C1CC1(16421)',
    structure = SMILES('O=C=C[C]1CC1'),
    E0 = (172.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38052,0.0246252,2.46765e-05,-4.76415e-08,1.92892e-11,20799.7,18.3224], Tmin=(100,'K'), Tmax=(986.41,'K')), NASAPolynomial(coeffs=[10.429,0.0187548,-7.1002e-06,1.34454e-09,-9.81767e-14,17909.7,-26.9976], Tmin=(986.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + ring(Cyclopropane) + radical(CCJ(C)C=C=O)"""),
)

species(
    label = 'C=C1C=C([O])C1(16199)',
    structure = SMILES('C=C1C=C([O])C1'),
    E0 = (142.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00816,0.0242072,5.14884e-05,-9.63447e-08,4.23359e-11,17285.3,17.2213], Tmin=(100,'K'), Tmax=(920.7,'K')), NASAPolynomial(coeffs=[18.3514,0.00397859,1.72241e-06,-4.11476e-10,2.26176e-14,12123.8,-71.9624], Tmin=(920.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.97,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C(C)=C=C=O(16422)',
    structure = SMILES('C=C(C)[C]=C=O'),
    E0 = (174.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.15306,0.044099,-4.27717e-05,2.75623e-08,-8.13148e-12,21001.7,19.6058], Tmin=(100,'K'), Tmax=(785.885,'K')), NASAPolynomial(coeffs=[5.18647,0.0286605,-1.33064e-05,2.56847e-09,-1.81115e-13,20524.9,5.70208], Tmin=(785.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)CJ=C=O)"""),
)

species(
    label = '[CH2][C]1CC1[C]=O(16423)',
    structure = SMILES('[CH2][C]1CC1[C]=O'),
    E0 = (446.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23889,0.0365384,-2.34607e-05,8.38496e-09,-1.28835e-12,53719,23.4101], Tmin=(100,'K'), Tmax=(1454.62,'K')), NASAPolynomial(coeffs=[7.58124,0.0218477,-8.31167e-06,1.44201e-09,-9.50874e-14,52164.8,-4.36571], Tmin=(1454.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Tertalkyl) + radical(Isobutyl) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH2]C1([CH2])C=C1[O](16215)',
    structure = SMILES('[CH2]C1([CH2])[CH]C1=O'),
    E0 = (518.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.47933,0.0651742,-7.00093e-05,3.63809e-08,-7.11558e-12,62512.4,21.8731], Tmin=(100,'K'), Tmax=(1403.99,'K')), NASAPolynomial(coeffs=[18.0549,0.00792516,-1.17883e-06,5.72547e-11,5.66031e-16,58284.4,-66.364], Tmin=(1403.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(518.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsCs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CCJC=O) + radical(CJC(C)2C=O) + radical(CJC(C)2C=O)"""),
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
    label = '[CH2]C([CH2])=C=C=O(16424)',
    structure = SMILES('[CH2]C(=C)[C]=C=O'),
    E0 = (332.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,2120,512.5,787.5,405.755],'cm^-1')),
        HinderedRotor(inertia=(0.0143527,'amu*angstrom^2'), symmetry=1, barrier=(28.5283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.71656,'amu*angstrom^2'), symmetry=1, barrier=(62.459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17443,0.038353,-2.71919e-05,9.92722e-09,-1.50259e-12,40015.5,20.7715], Tmin=(100,'K'), Tmax=(1504.1,'K')), NASAPolynomial(coeffs=[9.3411,0.019294,-8.18494e-06,1.50274e-09,-1.0234e-13,37859.6,-16.7289], Tmin=(1504.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(CJ)C=C=O) + radical(C=C(C)CJ=C=O)"""),
)

species(
    label = 'C=C=C[C]=O(9908)',
    structure = SMILES('C=C=C[C]=O'),
    E0 = (213.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.13978,'amu*angstrom^2'), symmetry=1, barrier=(26.2058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00684,0.0284132,-1.9291e-05,6.01948e-09,-7.90367e-13,25757,13.5223], Tmin=(100,'K'), Tmax=(1551.73,'K')), NASAPolynomial(coeffs=[6.76337,0.0187298,-9.93048e-06,1.99793e-09,-1.42457e-13,24591.2,-6.25135], Tmin=(1551.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]C([CH2])=[C][C]=O(16425)',
    structure = SMILES('[CH2]C([CH2])=[C][C]=O'),
    E0 = (551.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.04859,'amu*angstrom^2'), symmetry=1, barrier=(24.1091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04797,'amu*angstrom^2'), symmetry=1, barrier=(24.0948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04673,'amu*angstrom^2'), symmetry=1, barrier=(24.0664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37854,0.0618654,-7.68575e-05,4.8691e-08,-1.24091e-11,66380.3,21.2148], Tmin=(100,'K'), Tmax=(949.276,'K')), NASAPolynomial(coeffs=[11.3264,0.0199475,-1.06207e-05,2.17348e-09,-1.58306e-13,64491.7,-26.2601], Tmin=(949.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C([CH2])=[C]C=O(16426)',
    structure = SMILES('[CH2]C(=C)[C]=C[O]'),
    E0 = (339.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.66404,'amu*angstrom^2'), symmetry=1, barrier=(38.2596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.66303,'amu*angstrom^2'), symmetry=1, barrier=(38.2363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15684,0.0494782,-1.36073e-05,-3.29225e-08,2.09859e-11,40946.5,19.9492], Tmin=(100,'K'), Tmax=(904.914,'K')), NASAPolynomial(coeffs=[19.1692,0.00529734,8.8255e-07,-3.18907e-10,2.2043e-14,36235.5,-73.1677], Tmin=(904.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(339.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_P) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C(C)=[C][C]=O(16427)',
    structure = SMILES('[CH2]C(C)=[C][C]=O'),
    E0 = (399.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.81042,'amu*angstrom^2'), symmetry=1, barrier=(18.6332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799855,'amu*angstrom^2'), symmetry=1, barrier=(18.3902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.860507,'amu*angstrom^2'), symmetry=1, barrier=(19.7847,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28961,0.065411,-8.9151e-05,6.71869e-08,-2.08349e-11,48160.6,21.2538], Tmin=(100,'K'), Tmax=(780.327,'K')), NASAPolynomial(coeffs=[8.94064,0.0261916,-1.37611e-05,2.77849e-09,-1.9995e-13,46966.6,-13.7602], Tmin=(780.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
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
    label = '[CH2][C]=C(459)',
    structure = SMILES('[CH2][C]=C'),
    E0 = (395.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.110742,'amu*angstrom^2'), symmetry=1, barrier=(23.6192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28575,0.0103679,1.6796e-05,-2.70041e-08,1.05227e-11,47575.3,10.4038], Tmin=(100,'K'), Tmax=(980.899,'K')), NASAPolynomial(coeffs=[6.52804,0.0104286,-3.60844e-06,6.68464e-10,-4.85509e-14,46300.2,-8.43297], Tmin=(980.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C1CC1[C]=O(16162)',
    structure = SMILES('C=C1CC1[C]=O'),
    E0 = (232.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05078,0.0352586,-7.3581e-06,-1.43847e-08,7.60992e-12,27999,18.4145], Tmin=(100,'K'), Tmax=(1022.34,'K')), NASAPolynomial(coeffs=[10.4813,0.0187541,-7.32326e-06,1.36095e-09,-9.64156e-14,25413.9,-26.6562], Tmin=(1022.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=[C][CH][C]1CC1(16428)',
    structure = SMILES('O=[C][CH][C]1CC1'),
    E0 = (415.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07545,0.0366533,-2.2294e-05,7.04581e-09,-9.16778e-13,50057.3,21.7906], Tmin=(100,'K'), Tmax=(1750.39,'K')), NASAPolynomial(coeffs=[10.001,0.0185416,-6.773e-06,1.1343e-09,-7.24532e-14,47282.7,-20.8823], Tmin=(1750.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(415.584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Tertalkyl) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]1[CH]C(=O)C1(16429)',
    structure = SMILES('[CH2][C]1C=C([O])C1'),
    E0 = (373.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38279,0.0200451,5.00902e-05,-8.65424e-08,3.73737e-11,44937.2,19.3924], Tmin=(100,'K'), Tmax=(908.176,'K')), NASAPolynomial(coeffs=[13.867,0.0101016,-6.06221e-07,-5.68593e-11,3.12843e-15,41175.4,-44.1325], Tmin=(908.176,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Allyl_T) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C(=C)C[C]=O(15749)',
    structure = SMILES('[CH]C(=C)C[C]=O'),
    E0 = (398.451,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,409.495,409.594,409.708,409.814],'cm^-1')),
        HinderedRotor(inertia=(0.412128,'amu*angstrom^2'), symmetry=1, barrier=(49.0855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.412124,'amu*angstrom^2'), symmetry=1, barrier=(49.0851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.41201,'amu*angstrom^2'), symmetry=1, barrier=(49.0817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64491,0.0440108,-1.81928e-05,-3.18729e-09,2.77897e-12,48013.7,22.8183], Tmin=(100,'K'), Tmax=(1230.57,'K')), NASAPolynomial(coeffs=[11.7847,0.0247584,-1.14335e-05,2.2026e-09,-1.54967e-13,44480.4,-32.4211], Tmin=(1230.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(C)[CH][C]=O(16430)',
    structure = SMILES('[CH]C(C)=C[C]=O'),
    E0 = (374.887,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.08095,'amu*angstrom^2'), symmetry=1, barrier=(47.8452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07997,'amu*angstrom^2'), symmetry=1, barrier=(47.8226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07981,'amu*angstrom^2'), symmetry=1, barrier=(47.819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31184,0.0645295,-8.42334e-05,7.05664e-08,-2.53515e-11,45180.1,21.2842], Tmin=(100,'K'), Tmax=(733.831,'K')), NASAPolynomial(coeffs=[5.9949,0.0354621,-1.758e-05,3.43834e-09,-2.42454e-13,44588.1,0.790077], Tmin=(733.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]=C([CH2])[CH]C=O(16431)',
    structure = SMILES('[CH]C(=C)C=C[O]'),
    E0 = (359.685,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.04994,'amu*angstrom^2'), symmetry=1, barrier=(47.1322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04968,'amu*angstrom^2'), symmetry=1, barrier=(47.1261,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23547,0.0453703,7.13725e-06,-5.25549e-08,2.66821e-11,43373.9,20.8973], Tmin=(100,'K'), Tmax=(926.731,'K')), NASAPolynomial(coeffs=[18.256,0.0115273,-2.21611e-06,3.08098e-10,-2.39635e-14,38517.8,-69.1015], Tmin=(926.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)C=C=O(15565)',
    structure = SMILES('[CH2]C(=C)C=C=O'),
    E0 = (144.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,2120,512.5,787.5,196.042],'cm^-1')),
        HinderedRotor(inertia=(0.548013,'amu*angstrom^2'), symmetry=1, barrier=(14.9464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29539,'amu*angstrom^2'), symmetry=1, barrier=(35.3181,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79035,0.0478152,-4.18093e-05,1.97259e-08,-3.84494e-12,17449.8,19.3405], Tmin=(100,'K'), Tmax=(1209.73,'K')), NASAPolynomial(coeffs=[9.77577,0.0214113,-9.06985e-06,1.68365e-09,-1.16369e-13,15517.8,-20.705], Tmin=(1209.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(CJ)C=C=O)"""),
)

species(
    label = 'C=[C]CC=[C][O](15568)',
    structure = SMILES('C=[C]C[CH][C]=O'),
    E0 = (442.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,488.96,3187.22],'cm^-1')),
        HinderedRotor(inertia=(0.0736111,'amu*angstrom^2'), symmetry=1, barrier=(12.6387,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.547723,'amu*angstrom^2'), symmetry=1, barrier=(12.5932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.629662,'amu*angstrom^2'), symmetry=1, barrier=(78.6091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3431.29,'J/mol'), sigma=(5.7779,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.96 K, Pc=40.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92383,0.0464097,-4.41284e-05,2.41862e-08,-5.55934e-12,53262.2,22.8634], Tmin=(100,'K'), Tmax=(1032.8,'K')), NASAPolynomial(coeffs=[8.11657,0.0224254,-9.29441e-06,1.70107e-09,-1.16554e-13,51983,-7.21289], Tmin=(1032.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH2])[C][C]=O(9809)',
    structure = SMILES('[CH2]C([CH2])[C][C]=O'),
    E0 = (758.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,180,1164.67],'cm^-1')),
        HinderedRotor(inertia=(0.0061557,'amu*angstrom^2'), symmetry=1, barrier=(5.93723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257638,'amu*angstrom^2'), symmetry=1, barrier=(5.9236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257948,'amu*angstrom^2'), symmetry=1, barrier=(5.93074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.74058,'amu*angstrom^2'), symmetry=1, barrier=(63.0114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33432,0.0616237,-8.70494e-05,6.72732e-08,-2.0086e-11,91309.5,26.7625], Tmin=(100,'K'), Tmax=(958.22,'K')), NASAPolynomial(coeffs=[8.88883,0.0203607,-7.22876e-06,1.14499e-09,-6.90815e-14,90308.3,-7.03094], Tmin=(958.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(758.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]1C=[C]OC1(16432)',
    structure = SMILES('[CH2][C]1C=[C]OC1'),
    E0 = (437.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68383,0.0111066,7.64078e-05,-1.16919e-07,4.98467e-11,52638.2,18.8652], Tmin=(100,'K'), Tmax=(886.709,'K')), NASAPolynomial(coeffs=[14.0045,0.00752205,2.14566e-06,-6.93116e-10,5.05976e-14,48763.9,-44.9155], Tmin=(886.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1([CH2])C=[C]O1(16192)',
    structure = SMILES('[CH2]C1([CH2])C=[C]O1'),
    E0 = (574.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.137687,0.0672495,-7.20302e-05,3.57176e-08,-6.42737e-12,69259.5,25.8046], Tmin=(100,'K'), Tmax=(1640.53,'K')), NASAPolynomial(coeffs=[19.6426,0.00249588,2.28575e-06,-6.22404e-10,4.61543e-14,64993.1,-72.6382], Tmin=(1640.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=C)[C]=[C]O(16433)',
    structure = SMILES('[CH2]C(=C)[C]=[C]O'),
    E0 = (437.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(1.34356,'amu*angstrom^2'), symmetry=1, barrier=(30.891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34107,'amu*angstrom^2'), symmetry=1, barrier=(30.8338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34099,'amu*angstrom^2'), symmetry=1, barrier=(30.8319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.248036,0.0665984,-7.3572e-05,3.87681e-08,-7.52519e-12,52801.2,24.9257], Tmin=(100,'K'), Tmax=(1484.17,'K')), NASAPolynomial(coeffs=[18.3036,0.0053864,9.76894e-07,-4.15542e-10,3.47649e-14,48823.9,-64.6549], Tmin=(1484.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([CH2])C=[C]O(16434)',
    structure = SMILES('[CH]C(=C)C=[C]O'),
    E0 = (457.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.05717,'amu*angstrom^2'), symmetry=1, barrier=(47.2984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05804,'amu*angstrom^2'), symmetry=1, barrier=(47.3183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.05405,'amu*angstrom^2'), symmetry=1, barrier=(47.2266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03736,0.0540461,-2.29091e-05,-2.0398e-08,1.54337e-11,55197.7,23.3275], Tmin=(100,'K'), Tmax=(912.153,'K')), NASAPolynomial(coeffs=[17.4268,0.01187,-2.3851e-06,2.92351e-10,-1.90469e-14,50972.4,-61.0069], Tmin=(912.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH2](15952)',
    structure = SMILES('[CH2][C][CH2]'),
    E0 = (738.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,618.688],'cm^-1')),
        HinderedRotor(inertia=(0.00838111,'amu*angstrom^2'), symmetry=1, barrier=(2.27602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0989004,'amu*angstrom^2'), symmetry=1, barrier=(2.27392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.05819,0.0178835,-7.35728e-06,-7.95768e-10,1.01791e-12,88882.1,14.263], Tmin=(100,'K'), Tmax=(1175.48,'K')), NASAPolynomial(coeffs=[6.73198,0.0102413,-3.80615e-06,7.07005e-10,-4.9648e-14,87682.7,-5.48288], Tmin=(1175.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]([CH2])[C]=C[O](16435)',
    structure = SMILES('[CH2]C([CH2])=[C][CH][O]'),
    E0 = (661.748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,1685,370,362.972,363.651],'cm^-1')),
        HinderedRotor(inertia=(0.00125664,'amu*angstrom^2'), symmetry=1, barrier=(0.119762,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00130752,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00129714,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62256,0.0461737,-3.37856e-05,1.23243e-08,-1.80886e-12,79680.5,24.5224], Tmin=(100,'K'), Tmax=(1595.85,'K')), NASAPolynomial(coeffs=[12.7599,0.0182578,-7.54624e-06,1.36274e-09,-9.16585e-14,76125.8,-34.4145], Tmin=(1595.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH2])C#C[O](9818)',
    structure = SMILES('[CH2]C([CH2])[C]=C=O'),
    E0 = (493.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.362801,'amu*angstrom^2'), symmetry=1, barrier=(8.34151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.362475,'amu*angstrom^2'), symmetry=1, barrier=(8.33402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.50701,'amu*angstrom^2'), symmetry=1, barrier=(103.625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02895,0.0729765,-0.000120849,1.06643e-07,-3.61612e-11,59404.3,23.265], Tmin=(100,'K'), Tmax=(872.549,'K')), NASAPolynomial(coeffs=[7.77751,0.0258499,-1.20021e-05,2.21523e-09,-1.48194e-13,58842.9,-4.84133], Tmin=(872.549,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O) + radical(CJC(C)C=C=O) + radical(CC(C)CJ=C=O)"""),
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
    E0 = (341.124,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (870.4,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (982.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (738.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (315.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (315.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (370.601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (446.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (518.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (389.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (560.625,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (610.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (762.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (502.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (555.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (875.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (312.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (415.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (436.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (615.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (590.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (504.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (568.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (307.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (612.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (781.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (437.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (574.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (583.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (630.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (491.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (959.287,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (684.609,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (670.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['HCCO(2227)', 'allene(458)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(33.9238,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission
Ea raised from 0.0 to 33.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(20)', '[CH2][C]=C[C]=O(9897)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(2355)', '[CH]=C([CH2])[CH2](16419)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['O=[C]C=C1CC1(16421)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['C=C1C=C([O])C1(16199)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['[CH2]C(C)=C=C=O(16422)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['[CH2][C]1CC1[C]=O(16423)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.36786e+12,'s^-1'), n=-0.105173, Ea=(138.904,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 135.1 to 138.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['[CH2]C1([CH2])C=C1[O](16215)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(211.456,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CO(2039)', '[CH]=C([CH2])[CH2](16419)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.51e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2]C([CH2])=C=C=O(16424)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Cdd;HJ] for rate rule [Ca_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['CH2(T)(20)', 'C=C=C[C]=O(9908)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.70446,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH2]C([CH2])=[C][C]=O(16425)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH2])=[C]C=O(16426)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C)=[C][C]=O(16427)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C][O](6861)', '[CH2][C]=C(459)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['C=C1CC1[C]=O(16162)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.18842e+14,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['O=[C][CH][C]1CC1(16428)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(108.383,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 104.4 to 108.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['[CH2][C]1[CH]C(=O)C1(16429)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.5455e+09,'s^-1'), n=0.708429, Ea=(129.065,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;doublebond_intra;radadd_intra] for rate rule [R4_linear;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C][O](6861)', 'allene(458)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.70446,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C(=C)C[C]=O(15749)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C(C)[CH][C]=O(16430)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C([CH2])[CH]C=O(16431)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(868165,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_2;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['[CH2]C(=C)C=C=O(15565)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]CC=[C][O](15568)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH2])[C][C]=O(9809)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['[CH2][C]1C=[C]OC1(16432)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.10072e+11,'s^-1'), n=0.23641, Ea=(129.924,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 123.3 to 129.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C([CH2])=C[C]=O(9813)'],
    products = ['[CH2]C1([CH2])C=[C]O1(16192)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(267.247,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 266.6 to 267.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['HCCO(2227)', '[CH2][C]=C(459)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(=C)[C]=[C]O(16433)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C([CH2])C=[C]O(16434)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['HCCO(2227)', '[CH2][C][CH2](15952)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][C]([CH2])[C]=C[O](16435)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH2])C#C[O](9818)'],
    products = ['[CH2]C([CH2])=C[C]=O(9813)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3674',
    isomers = [
        '[CH2]C([CH2])=C[C]=O(9813)',
    ],
    reactants = [
        ('HCCO(2227)', 'allene(458)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #3674',
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


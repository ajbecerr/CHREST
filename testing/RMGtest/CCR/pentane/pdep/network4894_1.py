species(
    label = '[CH]=[C]OC(=C)[CH]O[O](20826)',
    structure = SMILES('[CH]=[C]OC(=C)[CH]O[O]'),
    E0 = (602.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,1011.24,2529.22],'cm^-1')),
        HinderedRotor(inertia=(0.75715,'amu*angstrom^2'), symmetry=1, barrier=(17.4084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755273,'amu*angstrom^2'), symmetry=1, barrier=(17.3652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.757444,'amu*angstrom^2'), symmetry=1, barrier=(17.4151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162372,'amu*angstrom^2'), symmetry=1, barrier=(69.6184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.957147,0.0682026,-7.7276e-05,4.56742e-08,-1.08186e-11,72576.4,33.3776], Tmin=(100,'K'), Tmax=(1023.59,'K')), NASAPolynomial(coeffs=[12.6843,0.0223757,-1.01207e-05,1.93633e-09,-1.36267e-13,70175.6,-23.4728], Tmin=(1023.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = 'C=C=CO[O](16806)',
    structure = SMILES('C=C=CO[O]'),
    E0 = (250.435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.895395,'amu*angstrom^2'), symmetry=1, barrier=(20.5869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3584.1,'J/mol'), sigma=(5.7752,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.83 K, Pc=42.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2049,0.0363638,-3.70353e-05,1.95556e-08,-4.06022e-12,30187.6,17.323], Tmin=(100,'K'), Tmax=(1179.28,'K')), NASAPolynomial(coeffs=[9.983,0.00998085,-3.47674e-06,5.84074e-10,-3.83218e-14,28353.2,-21.4843], Tmin=(1179.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
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
    label = 'C=C([O])[CH]O[O](21019)',
    structure = SMILES('[CH2]C(=O)[CH]O[O]'),
    E0 = (157.104,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,375,552.5,462.5,1710,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,299.468],'cm^-1')),
        HinderedRotor(inertia=(0.0771649,'amu*angstrom^2'), symmetry=1, barrier=(4.91325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310443,'amu*angstrom^2'), symmetry=1, barrier=(19.7863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05278,'amu*angstrom^2'), symmetry=1, barrier=(67.1902,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45237,0.0587496,-8.09575e-05,5.69462e-08,-1.57319e-11,18984.6,21.1846], Tmin=(100,'K'), Tmax=(889.554,'K')), NASAPolynomial(coeffs=[11.2448,0.0147155,-6.70326e-06,1.29549e-09,-9.13534e-14,17242.5,-24.9116], Tmin=(889.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.104,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(OCJC=O) + radical(C2JC=O)"""),
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
    label = '[CH]=[C]OC(=C)[CH][O](23372)',
    structure = SMILES('[CH]=[C]OC([CH2])=C[O]'),
    E0 = (506.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,180,572.121,572.21],'cm^-1')),
        HinderedRotor(inertia=(0.139671,'amu*angstrom^2'), symmetry=1, barrier=(17.5669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0756225,'amu*angstrom^2'), symmetry=1, barrier=(17.5667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.764047,'amu*angstrom^2'), symmetry=1, barrier=(17.5669,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596282,0.0661455,-7.45336e-05,4.08059e-08,-8.49436e-12,61041.8,30.7505], Tmin=(100,'K'), Tmax=(1262.53,'K')), NASAPolynomial(coeffs=[17.5756,0.0087312,-2.01977e-06,2.44892e-10,-1.30213e-14,57042.9,-53.9808], Tmin=(1262.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(O)CJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36761e-06,-4.93131e-09,1.45956e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.7,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00338e-07,1.59031e-10,-1.14892e-14,-1048.44,6.08305], Tmin=(1087.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH]C(=C)O[C]=[CH](16185)',
    structure = SMILES('[CH]C(=C)O[C]=[CH]'),
    E0 = (785.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,437.732,437.738,437.738,437.743,437.744],'cm^-1')),
        HinderedRotor(inertia=(0.382557,'amu*angstrom^2'), symmetry=1, barrier=(52.0177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.382551,'amu*angstrom^2'), symmetry=1, barrier=(52.0178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.382544,'amu*angstrom^2'), symmetry=1, barrier=(52.0175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3230.54,'J/mol'), sigma=(5.55692,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=504.60 K, Pc=42.72 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80405,0.0501869,-4.64286e-05,2.49013e-08,-5.72833e-12,94557,25.5735], Tmin=(100,'K'), Tmax=(1017.5,'K')), NASAPolynomial(coeffs=[7.77914,0.0266971,-1.17992e-05,2.21162e-09,-1.53347e-13,93341.1,-3.35641], Tmin=(1017.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(AllylJ2_triplet) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]O[O](21387)',
    structure = SMILES('[CH]O[O]'),
    E0 = (471.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.569664,'amu*angstrom^2'), symmetry=1, barrier=(13.0977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (45.0174,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.33696,0.0153854,-2.48163e-05,1.94592e-08,-5.80272e-12,56680.3,11.4536], Tmin=(100,'K'), Tmax=(834.588,'K')), NASAPolynomial(coeffs=[6.14825,0.00191155,-5.99818e-07,1.15174e-10,-8.22083e-15,56211.1,-1.6009], Tmin=(834.588,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=[C]O[C]=C(9939)',
    structure = SMILES('[CH]=[C]O[C]=C'),
    E0 = (696.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650,180,1606.5],'cm^-1')),
        HinderedRotor(inertia=(0.521695,'amu*angstrom^2'), symmetry=1, barrier=(11.9948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.774936,'amu*angstrom^2'), symmetry=1, barrier=(17.8173,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47358,0.0352897,-4.15989e-05,2.86663e-08,-8.17618e-12,83810.8,23.1837], Tmin=(100,'K'), Tmax=(847.016,'K')), NASAPolynomial(coeffs=[6.63847,0.0156221,-6.77105e-06,1.25548e-09,-8.62002e-14,83105.2,3.78175], Tmin=(847.016,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(696.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]OC(=C)[C]O[O](23373)',
    structure = SMILES('[CH]=[C]OC([CH2])=[C]O[O]'),
    E0 = (881.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,242.965,2278.37],'cm^-1')),
        HinderedRotor(inertia=(0.186149,'amu*angstrom^2'), symmetry=1, barrier=(7.80608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186156,'amu*angstrom^2'), symmetry=1, barrier=(7.80755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186078,'amu*angstrom^2'), symmetry=1, barrier=(7.80642,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1861,'amu*angstrom^2'), symmetry=1, barrier=(7.8062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.636769,0.0834317,-0.000147666,1.31056e-07,-4.37774e-11,106100,37.1266], Tmin=(100,'K'), Tmax=(895.753,'K')), NASAPolynomial(coeffs=[9.11009,0.0236817,-1.09168e-05,1.9705e-09,-1.2851e-13,105461,2.08752], Tmin=(895.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(881.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=C(O)CJ) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[C]=[C]OC(=C)[CH]O[O](23374)',
    structure = SMILES('[C]=[C]OC(=C)[CH]O[O]'),
    E0 = (913.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,355.666,355.698,355.728,355.732],'cm^-1')),
        HinderedRotor(inertia=(0.476464,'amu*angstrom^2'), symmetry=1, barrier=(42.7799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165726,'amu*angstrom^2'), symmetry=1, barrier=(14.8825,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00133252,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00416032,'amu*angstrom^2'), symmetry=1, barrier=(42.7806,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.958771,0.0714232,-9.88854e-05,7.29314e-08,-2.15785e-11,109979,33.3502], Tmin=(100,'K'), Tmax=(825.423,'K')), NASAPolynomial(coeffs=[10.9032,0.0232317,-1.13081e-05,2.19688e-09,-1.54414e-13,108337,-12.7177], Tmin=(825.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(913.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CJO) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=[C]OC(=C)C1OO1(23375)',
    structure = SMILES('C#CO[C]([CH2])C1OO1'),
    E0 = (412.763,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.275746,0.0885077,-0.000105953,5.27991e-08,-7.14234e-12,49803.2,27.194], Tmin=(100,'K'), Tmax=(872.731,'K')), NASAPolynomial(coeffs=[22.4882,0.00549712,7.0249e-08,-1.94297e-10,1.81805e-14,45017.8,-84.1827], Tmin=(872.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + ring(dioxirane) + radical(C2CsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH]=C1OC(=C)C1O[O](23376)',
    structure = SMILES('[CH]=C1OC(=C)C1O[O]'),
    E0 = (364.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30528,0.047913,-1.61337e-05,-1.60887e-08,9.8506e-12,43946.1,27.7166], Tmin=(100,'K'), Tmax=(1024.83,'K')), NASAPolynomial(coeffs=[15.3976,0.0173768,-7.25097e-06,1.42902e-09,-1.05641e-13,39772.8,-46.885], Tmin=(1024.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]O[C]1CC1O[O](23377)',
    structure = SMILES('[CH]=[C]O[C]1CC1O[O]'),
    E0 = (678.762,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.254669,0.075716,-8.73488e-05,4.85562e-08,-1.03492e-11,81776.7,31.7958], Tmin=(100,'K'), Tmax=(1161.63,'K')), NASAPolynomial(coeffs=[19.0136,0.0111207,-3.93753e-06,6.85922e-10,-4.67587e-14,77418.6,-61.5157], Tmin=(1161.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(678.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(ROOJ) + radical(C2CsJOC(O)) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C[C]([CH]O[O])O1(23378)',
    structure = SMILES('[CH]=C1C[C]([CH]O[O])O1'),
    E0 = (607.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.459572,0.0745725,-8.45171e-05,4.3036e-08,-6.52162e-12,73246.1,26.4943], Tmin=(100,'K'), Tmax=(862.315,'K')), NASAPolynomial(coeffs=[17.0016,0.0128563,-3.28339e-06,4.28272e-10,-2.38171e-14,69834.9,-54.0981], Tmin=(862.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(ROOJ) + radical(C2CsJOC(O)) + radical(CCsJOOH) + radical(Cds_P)"""),
)

species(
    label = '[O]O[CH][C]1CC=[C]O1(23379)',
    structure = SMILES('[O]O[CH][C]1CC=[C]O1'),
    E0 = (527.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11937,0.0613531,-6.34301e-05,3.26435e-08,-5.66722e-12,63522.5,26.8875], Tmin=(100,'K'), Tmax=(871.475,'K')), NASAPolynomial(coeffs=[12.3698,0.0189149,-6.2206e-06,9.93239e-10,-6.29272e-14,61212.3,-27.846], Tmin=(871.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(ROOJ) + radical(C2CsJOC(O)) + radical(CCsJOOH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]O[C]1[CH]OOC1(23380)',
    structure = SMILES('[CH]=[C]O[C]1[CH]OOC1'),
    E0 = (647.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.156962,0.0755053,-8.42276e-05,4.5854e-08,-9.21652e-12,78060.1,31.7185], Tmin=(100,'K'), Tmax=(1439.77,'K')), NASAPolynomial(coeffs=[18.7055,0.00857242,6.41743e-07,-4.52393e-10,4.10522e-14,74134.4,-60.9277], Tmin=(1439.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.67,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(C2CsJOC(O)) + radical(CCsJOOC) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1OC1([CH2])[CH]O[O](23381)',
    structure = SMILES('[CH]=C1OC1([CH2])[CH]O[O]'),
    E0 = (643.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.237074,0.0884162,-0.000109498,5.87403e-08,-9.93122e-12,77494.5,27.472], Tmin=(100,'K'), Tmax=(878.967,'K')), NASAPolynomial(coeffs=[21.8944,0.00612958,-5.21837e-07,-5.99559e-11,8.32108e-15,72892.1,-80.4942], Tmin=(878.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(ROOJ) + radical(CCsJOOH) + radical(CJC(C)OC) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1([CH]O[O])C=[C]O1(23382)',
    structure = SMILES('[CH2]C1([CH]O[O])C=[C]O1'),
    E0 = (616.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.1635,0.0816593,-9.74554e-05,5.5498e-08,-1.19203e-11,74270.4,30.4997], Tmin=(100,'K'), Tmax=(1239.08,'K')), NASAPolynomial(coeffs=[21.0054,0.00772726,-1.18227e-06,5.58345e-11,1.03341e-15,69453.9,-74.4327], Tmin=(1239.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(616.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(ROOJ) + radical(CCsJOOH) + radical(CJC(C)OC) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OC1([CH2])[CH]OO1(23383)',
    structure = SMILES('[CH]=[C]OC1([CH2])[CH]OO1'),
    E0 = (695.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.50316,0.090055,-0.000100879,5.03798e-08,-9.09743e-12,83815.6,33.1974], Tmin=(100,'K'), Tmax=(1621.46,'K')), NASAPolynomial(coeffs=[27.3415,-0.00238406,4.32364e-06,-9.68878e-10,6.76386e-14,77259.1,-111.276], Tmin=(1621.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(695.013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(CCsJOO) + radical(CJCOOH) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(=C)C=O(23384)',
    structure = SMILES('C#COC([CH2])=C[O]'),
    E0 = (228.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23465,'amu*angstrom^2'), symmetry=1, barrier=(28.3869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23475,'amu*angstrom^2'), symmetry=1, barrier=(28.3893,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23388,'amu*angstrom^2'), symmetry=1, barrier=(28.3694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.239099,0.0690469,-7.57124e-05,3.9158e-08,-7.58081e-12,27681.7,25.9616], Tmin=(100,'K'), Tmax=(1417.01,'K')), NASAPolynomial(coeffs=[20.0008,0.00562349,-4.88127e-07,-3.75954e-11,5.59142e-15,22848.1,-73.5592], Tmin=(1417.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=COJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2][C]=CO[O](16807)',
    structure = SMILES('C=[C][CH]O[O]'),
    E0 = (437.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,274.987],'cm^-1')),
        HinderedRotor(inertia=(0.170957,'amu*angstrom^2'), symmetry=1, barrier=(9.18865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171244,'amu*angstrom^2'), symmetry=1, barrier=(9.18809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2197,0.041819,-5.77934e-05,4.42117e-08,-1.35872e-11,52657.9,18.5046], Tmin=(100,'K'), Tmax=(795.764,'K')), NASAPolynomial(coeffs=[7.632,0.0146118,-6.50507e-06,1.24122e-09,-8.65948e-14,51796.6,-6.36981], Tmin=(795.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S)"""),
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
    label = '[CH][C]([CH]O[O])OC#C(20815)',
    structure = SMILES('[CH]=[C]OC(=[CH])[CH]O[O]'),
    E0 = (849.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3115,3125,620,680,785,800,1600,1700,1685,370,464.458,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0439339,'amu*angstrom^2'), symmetry=1, barrier=(15.7675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.68555,'amu*angstrom^2'), symmetry=1, barrier=(15.7621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.6545,'amu*angstrom^2'), symmetry=1, barrier=(38.0403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.68219,'amu*angstrom^2'), symmetry=1, barrier=(61.6687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852592,0.0720528,-9.49292e-05,6.41835e-08,-1.71651e-11,102298,34.1585], Tmin=(100,'K'), Tmax=(915.711,'K')), NASAPolynomial(coeffs=[12.8345,0.0197132,-9.19254e-06,1.76396e-09,-1.23678e-13,100103,-22.5921], Tmin=(915.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(849.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(=[CH])CO[O](19429)',
    structure = SMILES('[CH][C](CO[O])OC#C'),
    E0 = (679.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,2175,525,750,770,3400,2100,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3951.63,'J/mol'), sigma=(6.62301,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=617.24 K, Pc=30.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.408825,0.102984,-0.000168676,1.35803e-07,-4.20311e-11,81937.9,29.0397], Tmin=(100,'K'), Tmax=(871.728,'K')), NASAPolynomial(coeffs=[15.8439,0.0192798,-8.9394e-06,1.63154e-09,-1.07953e-13,79451.1,-45.1502], Tmin=(871.728,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(C2CsJOCs) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=COC(=[CH])[CH]O[O](20827)',
    structure = SMILES('[CH]=COC(=[CH])[CH]O[O]'),
    E0 = (609.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700,180,749.269],'cm^-1')),
        HinderedRotor(inertia=(0.75372,'amu*angstrom^2'), symmetry=1, barrier=(17.3295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0434747,'amu*angstrom^2'), symmetry=1, barrier=(17.3251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.257498,'amu*angstrom^2'), symmetry=1, barrier=(17.3255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51573,'amu*angstrom^2'), symmetry=1, barrier=(34.8495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.319521,0.0742207,-8.22501e-05,4.43711e-08,-9.23429e-12,73490.3,32.9343], Tmin=(100,'K'), Tmax=(1185.3,'K')), NASAPolynomial(coeffs=[18.4278,0.0131107,-4.91443e-06,8.73592e-10,-5.98262e-14,69197.6,-57.5057], Tmin=(1185.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(609.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH]O[O])O[C]=C(19780)',
    structure = SMILES('[CH]=C([CH]O[O])O[C]=C'),
    E0 = (602.532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,180,1011.24,2529.22],'cm^-1')),
        HinderedRotor(inertia=(0.75715,'amu*angstrom^2'), symmetry=1, barrier=(17.4084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755273,'amu*angstrom^2'), symmetry=1, barrier=(17.3652,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.757444,'amu*angstrom^2'), symmetry=1, barrier=(17.4151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162372,'amu*angstrom^2'), symmetry=1, barrier=(69.6184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.957147,0.0682026,-7.7276e-05,4.56742e-08,-1.08186e-11,72576.4,33.3776], Tmin=(100,'K'), Tmax=(1023.59,'K')), NASAPolynomial(coeffs=[12.6843,0.0223757,-1.01207e-05,1.93633e-09,-1.36267e-13,70175.6,-23.4728], Tmin=(1023.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]([CH]OO)OC#C(20818)',
    structure = SMILES('[CH]=[C]OC(=[CH])[CH]OO'),
    E0 = (697.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3025,407.5,1350,352.5,3115,3125,620,680,785,800,1600,1700,1685,370,694.194,694.674],'cm^-1')),
        HinderedRotor(inertia=(0.893184,'amu*angstrom^2'), symmetry=1, barrier=(20.5361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.92568,'amu*angstrom^2'), symmetry=1, barrier=(67.2671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.39661,'amu*angstrom^2'), symmetry=1, barrier=(32.1109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0597733,'amu*angstrom^2'), symmetry=1, barrier=(20.5258,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105693,'amu*angstrom^2'), symmetry=1, barrier=(2.4301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.642209,0.075034,-8.81575e-05,5.22222e-08,-1.2256e-11,84024.7,34.1531], Tmin=(100,'K'), Tmax=(1037.39,'K')), NASAPolynomial(coeffs=[14.6767,0.0209196,-9.91177e-06,1.93866e-09,-1.38212e-13,81112.8,-34.0706], Tmin=(1037.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(697.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(Cds_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C][CH]O[O](22548)',
    structure = SMILES('[CH2][C][CH]O[O]'),
    E0 = (793.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,1493.49],'cm^-1')),
        HinderedRotor(inertia=(0.167737,'amu*angstrom^2'), symmetry=1, barrier=(3.85659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168521,'amu*angstrom^2'), symmetry=1, barrier=(3.87462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16787,'amu*angstrom^2'), symmetry=1, barrier=(3.85967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07407,0.0474498,-7.98325e-05,7.1518e-08,-2.45548e-11,95458.7,21.5901], Tmin=(100,'K'), Tmax=(855.952,'K')), NASAPolynomial(coeffs=[6.67897,0.0159696,-7.20999e-06,1.35978e-09,-9.27035e-14,95035.3,2.22186], Tmin=(855.952,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(793.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCsJOOH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]1O[C]=CC1O[O](23385)',
    structure = SMILES('[CH2][C]1O[C]=CC1O[O]'),
    E0 = (536.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796786,0.0661372,-6.7254e-05,2.82091e-08,-1.88038e-12,64612.6,28.2128], Tmin=(100,'K'), Tmax=(870.727,'K')), NASAPolynomial(coeffs=[15.9644,0.012651,-3.00698e-06,3.7526e-10,-2.06324e-14,61357.4,-46.3872], Tmin=(870.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(ROOJ) + radical(C2CsJOC(O)) + radical(CJC(C)OC) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]1[CH]OOC=[C]O1(23386)',
    structure = SMILES('C=C1[CH]OO[CH][C]O1'),
    E0 = (545.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.817033,0.00996696,0.00020345,-3.27127e-07,1.42393e-10,65827.3,25.111], Tmin=(100,'K'), Tmax=(905.192,'K')), NASAPolynomial(coeffs=[48.3536,-0.0422362,2.83675e-05,-5.52221e-09,3.61384e-13,50754.1,-235.214], Tmin=(905.192,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cycloheptane) + radical(C=CCJO) + radical(CCsJOOC) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=C1O[C]([CH2])C1O[O](23387)',
    structure = SMILES('[CH][C]1OC(=C)C1O[O]'),
    E0 = (607.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.341097,0.0761923,-8.31717e-05,4.28059e-08,-8.12092e-12,73235,29.689], Tmin=(100,'K'), Tmax=(1510.01,'K')), NASAPolynomial(coeffs=[21.2916,0.00513612,1.07377e-06,-4.19364e-10,3.40558e-14,68269.6,-78.4002], Tmin=(1510.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(ROOJ) + radical(C2CsJOC(O)) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=C1OO[CH][C]([CH2])O1(23388)',
    structure = SMILES('[CH][C]1OO[CH]C(=C)O1'),
    E0 = (486.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.985682,0.0365292,6.33796e-05,-1.20686e-07,5.03885e-11,58628.6,24.4088], Tmin=(100,'K'), Tmax=(980.827,'K')), NASAPolynomial(coeffs=[26.4005,0.00642343,-3.04648e-06,9.08014e-10,-8.87519e-14,50105.7,-115.744], Tmin=(980.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(C=CCJO) + radical(Cs_P) + radical(CCJ2_triplet)"""),
)

species(
    label = 'HC2(2881)',
    structure = SMILES('[C]#C'),
    E0 = (556.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01451,0.0139907,-3.08144e-05,3.10836e-08,-1.10946e-11,66983,5.75946], Tmin=(100,'K'), Tmax=(918.723,'K')), NASAPolynomial(coeffs=[3.14385,0.00498487,-2.32624e-06,4.08757e-10,-2.55787e-14,67315.6,7.08554], Tmin=(918.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.809,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""HC2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][C]([O])[CH]O[O](22833)',
    structure = SMILES('[CH2][C]([O])[CH]O[O]'),
    E0 = (577.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,2362.58],'cm^-1')),
        HinderedRotor(inertia=(0.284489,'amu*angstrom^2'), symmetry=1, barrier=(6.54095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.285315,'amu*angstrom^2'), symmetry=1, barrier=(6.55994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.59577,'amu*angstrom^2'), symmetry=1, barrier=(59.6819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (87.0541,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05836,0.0757261,-0.00014534,1.35191e-07,-4.65652e-11,69594.6,24.5151], Tmin=(100,'K'), Tmax=(893.702,'K')), NASAPolynomial(coeffs=[7.28479,0.0216789,-1.06872e-05,1.96768e-09,-1.29125e-13,69527.1,1.02502], Tmin=(893.702,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CsJOH) + radical(CCsJOOH) + radical(CJCO)"""),
)

species(
    label = '[C]#COC([CH2])[CH]O[O](23389)',
    structure = SMILES('[C]#COC([CH2])[CH]O[O]'),
    E0 = (787.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2175,525,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.636369,0.1089,-0.000183561,1.48992e-07,-4.57975e-11,94852.6,29.5495], Tmin=(100,'K'), Tmax=(912.212,'K')), NASAPolynomial(coeffs=[16.4989,0.0178448,-7.661e-06,1.31084e-09,-8.17785e-14,92388.6,-47.9137], Tmin=(912.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(787.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(CCsJOOH) + radical(CJC(C)OC) + radical(Acetyl)"""),
)

species(
    label = '[C]#CO[C]([CH2])CO[O](23390)',
    structure = SMILES('[C]#CO[C]([CH2])CO[O]'),
    E0 = (779.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,2175,525,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.722345,0.116683,-0.000216833,1.89824e-07,-6.16889e-11,93902.6,30.3476], Tmin=(100,'K'), Tmax=(916.287,'K')), NASAPolynomial(coeffs=[13.7485,0.0219302,-1.00179e-05,1.73545e-09,-1.07629e-13,92576.4,-30.9665], Tmin=(916.287,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(779.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(C2CsJOCs) + radical(CJC(C)OC) + radical(Acetyl)"""),
)

species(
    label = '[C]#CO[C](C)[CH]O[O](23391)',
    structure = SMILES('[C]#CO[C](C)[CH]O[O]'),
    E0 = (757.509,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,2175,525,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.355644,0.106523,-0.000189687,1.64458e-07,-5.35644e-11,91254,28.7356], Tmin=(100,'K'), Tmax=(903.783,'K')), NASAPolynomial(coeffs=[12.6479,0.0239088,-1.09763e-05,1.9512e-09,-1.25109e-13,89927.1,-27.0209], Tmin=(903.783,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(757.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(C2CsJOCs) + radical(CCsJOOH) + radical(Acetyl)"""),
)

species(
    label = '[C]#CO[C]([CH2])[CH]OO(23392)',
    structure = SMILES('[C]=[C]OC(=C)[CH]OO'),
    E0 = (761.533,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.844407,0.0732004,-8.75721e-05,5.46035e-08,-1.3735e-11,91701.7,33.0053], Tmin=(100,'K'), Tmax=(961.627,'K')), NASAPolynomial(coeffs=[12.4394,0.0249693,-1.2338e-05,2.44574e-09,-1.75146e-13,89471.7,-22.4802], Tmin=(961.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(761.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJO) + radical(C=CJO) + radical(CdCdJ2_triplet)"""),
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
    label = '[CH]=[C]O[C]=CO[O](23393)',
    structure = SMILES('C#CO[C][CH]O[O]'),
    E0 = (737.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13704,'amu*angstrom^2'), symmetry=1, barrier=(26.1429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13767,'amu*angstrom^2'), symmetry=1, barrier=(26.1574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1373,'amu*angstrom^2'), symmetry=1, barrier=(26.1489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13762,'amu*angstrom^2'), symmetry=1, barrier=(26.1562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.286213,0.083999,-0.000132567,9.73025e-08,-2.72058e-11,88781,22.5801], Tmin=(100,'K'), Tmax=(889.33,'K')), NASAPolynomial(coeffs=[17.2753,0.00758616,-3.68419e-06,6.88275e-10,-4.65187e-14,85759.2,-57.3897], Tmin=(889.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(737.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(CCsJOOH) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=C1CC(=CO[O])O1(23394)',
    structure = SMILES('[CH]=C1CC(=CO[O])O1'),
    E0 = (386.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35053,0.0451431,-4.73712e-06,-3.27279e-08,1.73711e-11,46596.1,29.0009], Tmin=(100,'K'), Tmax=(965.774,'K')), NASAPolynomial(coeffs=[16.6272,0.0133933,-4.38395e-06,8.24613e-10,-6.27732e-14,42175.2,-51.7792], Tmin=(965.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.53,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC1=COOC1(19203)',
    structure = SMILES('C#CO[C]1[CH]OOC1'),
    E0 = (350.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.326737,0.0785116,-8.83054e-05,4.82626e-08,-9.72091e-12,42342.9,25.393], Tmin=(100,'K'), Tmax=(1439.71,'K')), NASAPolynomial(coeffs=[19.3618,0.0083387,9.25447e-07,-5.20346e-10,4.6164e-14,38277.2,-71.1995], Tmin=(1439.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.644,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + ring(12dioxolane) + radical(C2CsJOCs) + radical(CCsJOOC)"""),
)

species(
    label = '[CH]=C1OOC=C([CH2])O1(23395)',
    structure = SMILES('[CH]=C1OOC=C([CH2])O1'),
    E0 = (365.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3327,0.0481898,-1.54729e-05,-1.11762e-08,6.16827e-12,44046.4,23.1165], Tmin=(100,'K'), Tmax=(1150.72,'K')), NASAPolynomial(coeffs=[14.0422,0.0247154,-1.18625e-05,2.36824e-09,-1.71337e-13,39750.6,-45.9399], Tmin=(1150.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(C=C(O)CJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC1=COC1(20824)',
    structure = SMILES('C#CO[C]1[CH]OC1'),
    E0 = (372.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0462133,0.074024,-8.90885e-05,5.12126e-08,-1.06942e-11,44991.4,22.7431], Tmin=(100,'K'), Tmax=(1414.97,'K')), NASAPolynomial(coeffs=[18.1027,0.00503959,2.78288e-06,-9.03028e-10,7.37919e-14,41625.1,-64.8608], Tmin=(1414.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.761,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Ct-CtOs) + group(Ct-CtH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CCsJOCs)"""),
)

species(
    label = '[CH]=C1OC=C([CH2])O1(23396)',
    structure = SMILES('[CH]=C1O[CH]C(=C)O1'),
    E0 = (240.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10765,0.0269136,3.28132e-05,-6.42288e-08,2.69218e-11,29033.3,21.4427], Tmin=(100,'K'), Tmax=(965.431,'K')), NASAPolynomial(coeffs=[13.7948,0.0148891,-5.05647e-06,9.73113e-10,-7.47903e-14,25080.4,-43.3151], Tmin=(965.431,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=CCJ(O)C) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]O[C]([CH2])C1OO1(23397)',
    structure = SMILES('[CH]=[C]O[C]([CH2])C1OO1'),
    E0 = (709.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.146954,0.0860101,-0.000103788,5.30906e-08,-7.90115e-12,85522.2,33.665], Tmin=(100,'K'), Tmax=(877.607,'K')), NASAPolynomial(coeffs=[21.7319,0.00587096,-2.82777e-07,-1.11848e-10,1.19788e-14,80927.9,-73.3277], Tmin=(877.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(709.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(dioxirane) + radical(C2CsJOC(O)) + radical(CJC(C)OC) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC([CH2])=[C]OO(23398)',
    structure = SMILES('[CH]=[C]OC([CH2])=[C]OO'),
    E0 = (729.226,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,3120,650,792.5,1650,322.851,3323.43],'cm^-1')),
        HinderedRotor(inertia=(0.530517,'amu*angstrom^2'), symmetry=1, barrier=(39.2214,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125034,'amu*angstrom^2'), symmetry=1, barrier=(9.24917,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125037,'amu*angstrom^2'), symmetry=1, barrier=(9.24845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125106,'amu*angstrom^2'), symmetry=1, barrier=(9.24871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61828,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.383243,0.0868959,-0.00014251,1.21205e-07,-3.98345e-11,87828.8,37.278], Tmin=(100,'K'), Tmax=(864.974,'K')), NASAPolynomial(coeffs=[10.7523,0.0252492,-1.18536e-05,2.1983e-09,-1.47558e-13,86547.4,-8.28072], Tmin=(864.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.226,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC(C)=[C]O[O](23399)',
    structure = SMILES('C#CO[C](C)[C]O[O]'),
    E0 = (692.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.547039,0.105602,-0.000173293,1.37946e-07,-4.22259e-11,83479.5,27.2643], Tmin=(100,'K'), Tmax=(857.588,'K')), NASAPolynomial(coeffs=[17.1488,0.0171294,-8.16679e-06,1.51157e-09,-1.00946e-13,80662.6,-54.1171], Tmin=(857.588,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(692.771,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(ROOJ) + radical(C2CsJOCs) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=COC([CH2])=[C]O[O](23400)',
    structure = SMILES('[CH]=COC([CH2])=[C]O[O]'),
    E0 = (641.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,335.326,335.326],'cm^-1')),
        HinderedRotor(inertia=(0.11061,'amu*angstrom^2'), symmetry=1, barrier=(8.826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110611,'amu*angstrom^2'), symmetry=1, barrier=(8.8261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110611,'amu*angstrom^2'), symmetry=1, barrier=(8.826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110608,'amu*angstrom^2'), symmetry=1, barrier=(8.8259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.436501,0.0815426,-0.000120226,9.13781e-08,-2.70847e-11,77278.3,34.7179], Tmin=(100,'K'), Tmax=(875.187,'K')), NASAPolynomial(coeffs=[13.231,0.0195912,-8.09139e-06,1.42405e-09,-9.32345e-14,75171.8,-24.5424], Tmin=(875.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=C(O)CJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OC(=C)[C]O[O](21821)',
    structure = SMILES('[CH2]C(=[C]O[O])O[C]=C'),
    E0 = (634.134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,228.221,228.231,2393.66],'cm^-1')),
        HinderedRotor(inertia=(0.00323691,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222854,'amu*angstrom^2'), symmetry=1, barrier=(8.23631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222826,'amu*angstrom^2'), symmetry=1, barrier=(8.23632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22283,'amu*angstrom^2'), symmetry=1, barrier=(8.23649,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723329,0.0797573,-0.000130483,1.13021e-07,-3.76066e-11,76379.4,36.4129], Tmin=(100,'K'), Tmax=(881.308,'K')), NASAPolynomial(coeffs=[8.79188,0.0266486,-1.20288e-05,2.18781e-09,-1.44923e-13,75597.5,2.13894], Tmin=(881.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(634.134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=C(O)CJ) + radical(C=CJO) + radical(C=CJO)"""),
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
    E0 = (602.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1004.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1025.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (832.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1222.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1093.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1125.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (610.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (610.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (833.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (728.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (634.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (661.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (658.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (656.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (696.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (602.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (633.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (862.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1061.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (825.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (654.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (793.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (730.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (958.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (660.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (733.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (686.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (779.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (713.913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (1134.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (865.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (844.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (799.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (855.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1174.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (610.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (609.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (610.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (685.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (695.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (712.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (871.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (857.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (685.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (816.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['HCCO(2227)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=[CH](9646)', 'C=C([O])[CH]O[O](21019)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH]=[C]OC(=C)[CH][O](23372)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O2(2)', '[CH]C(=C)O[C]=[CH](16185)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]O[O](21387)', '[CH]=[C]O[C]=C(9939)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]=[C]OC(=C)[C]O[O](23373)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[C]=[C]OC(=C)[CH]O[O](23374)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=[C]OC(=C)C1OO1(23375)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_H/OneDe]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=C1OC(=C)C1O[O](23376)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_H/NonDeO]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=[C]O[C]1CC1O[O](23377)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_csHO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=C1C[C]([CH]O[O])O1(23378)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[O]O[CH][C]1CC=[C]O1(23379)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(32.3534,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=[C]O[C]1[CH]OOC1(23380)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.80789e+10,'s^-1'), n=0.280222, Ea=(59.4174,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R5_linear;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=C1OC1([CH2])[CH]O[O](23381)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(56.4092,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH2]C1([CH]O[O])C=[C]O1(23382)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=[C]OC1([CH2])[CH]OO1(23383)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.20137e+08,'s^-1'), n=0.864, Ea=(93.7721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5;doublebond_intra_2H;radadd_intra_O] + [R5;doublebond_intra_2H_secNd;radadd_intra] for rate rule [R5;doublebond_intra_2H_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O(4)', '[CH]=[C]OC(=C)C=O(23384)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(373000,'cm^3/(mol*s)'), n=2.53, Ea=(130.591,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Od_CO-CdH;YJ] for rate rule [Od_CO-CdH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond
Ea raised from 125.5 to 130.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['HCCO(2227)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.41198,'m^3/(mol*s)'), n=2.03228, Ea=(31.3669,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CJ] + [Od_R;YJ] for rate rule [Od_R;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C][O](6861)', '[CH2][C]=CO[O](16807)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH][C]([CH]O[O])OC#C(20815)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C]OC(=[CH])CO[O](19429)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=COC(=[CH])[CH]O[O](20827)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C([CH]O[O])O[C]=C(19780)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH][C]([CH]OO)OC#C(20818)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HCCO(2227)', '[CH2][C][CH]O[O](22548)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH2][C]1O[C]=CC1O[O](23385)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(58.1576,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_csHNd] for rate rule [R5_linear;triplebond_intra_H;radadd_intra_csHO]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH2][C]1[CH]OOC=[C]O1(23386)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(2.07e+10,'s^-1'), n=0.124, Ea=(130.708,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;triplebond_intra_H;radadd_intra] for rate rule [R7_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=C1O[C]([CH2])C1O[O](23387)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(84.4457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;triplebond_intra_H;radadd_intra_cs] for rate rule [R5;triplebond_intra_H;radadd_intra_csHNd]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=C1OO[CH][C]([CH2])O1(23388)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.05e+09,'s^-1'), n=0.155, Ea=(177.192,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;triplebond_intra_H;radadd_intra] for rate rule [R7;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['HC2(2881)', 'C=C([O])[CH]O[O](21019)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5e+12,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1300,'K'), comment="""Estimated using template [R_R;CtJ_Ct] for rate rule [Od_R;CtJ_Ct]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['HC2(2881)', '[CH2][C]([O])[CH]O[O](22833)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(308.407,'m^3/(mol*s)'), n=0.967216, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0049236047651, var=0.0302625193734, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.361117102774
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[C]#COC([CH2])[CH]O[O](23389)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(10470.7,'s^-1'), n=2.29593, Ea=(78.4716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;XH_out] for rate rule [R4H_TSS;Ct_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[C]#CO[C]([CH2])CO[O](23390)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(26847.4,'s^-1'), n=2.09453, Ea=(64.7206,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5HJ_3;Ct_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[C]#CO[C](C)[CH]O[O](23391)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(505536,'s^-1'), n=1.7378, Ea=(41.5716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;Ct_rad_out;Cs_H_out_2H]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[C]#CO[C]([CH2])[CH]OO(23392)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1494.33,'s^-1'), n=2.17569, Ea=(93.5537,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;O_H_out] + [R7Hall;Y_rad_out;XH_out] for rate rule [R7Hall;Ct_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(T)(20)', '[CH]=[C]O[C]=CO[O](23393)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=C1CC(=CO[O])O1(23394)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=[C]OC1=COOC1(19203)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=C1OOC=C([CH2])O1(23395)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Ypri_rad_out] for rate rule [R6_SSDSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['O(4)', '[CH]=[C]OC1=COC1(20824)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;C_pri_rad_intra;OO] for rate rule [R3OO_SD;C_pri_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['O(4)', '[CH]=C1OC=C([CH2])O1(23396)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(9.14131e+09,'s^-1'), n=0.54, Ea=(93.1987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnOO;Y_rad_intra;OOJ] + [R4OO;Y_rad_intra;OO] for rate rule [R4OO;Y_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=[C]O[C]([CH2])C1OO1(23397)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(109.664,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=[C]OC([CH2])=[C]OO(23398)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    products = ['[CH]=[C]OC(C)=[C]O[O](23399)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.91e+11,'s^-1'), n=0.63, Ea=(255.224,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 199 used for R3H_SD;C_rad_out_2H;Cd_H_out_singleNd
Exact match found for rate rule [R3H_SD;C_rad_out_2H;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=COC([CH2])=[C]O[O](23400)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_single;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=[C]OC(=C)[C]O[O](21821)'],
    products = ['[CH]=[C]OC(=C)[CH]O[O](20826)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.74e+08,'s^-1'), n=1.713, Ea=(181.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleNd;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out_singleNd;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4894',
    isomers = [
        '[CH]=[C]OC(=C)[CH]O[O](20826)',
    ],
    reactants = [
        ('HCCO(2227)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4894',
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


species(
    label = '[O]C=[C]OC[C]=O(11163)',
    structure = SMILES('[O]C=[C]OC[C]=O'),
    E0 = (74.9892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,204.976,205.081,205.127,205.223],'cm^-1')),
        HinderedRotor(inertia=(0.813285,'amu*angstrom^2'), symmetry=1, barrier=(24.2916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813077,'amu*angstrom^2'), symmetry=1, barrier=(24.2911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813048,'amu*angstrom^2'), symmetry=1, barrier=(24.2914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995888,0.0523368,-2.47565e-05,-2.24611e-08,1.70709e-11,9139.99,26.7528], Tmin=(100,'K'), Tmax=(927.563,'K')), NASAPolynomial(coeffs=[21.2353,0.000483082,1.80928e-06,-3.79357e-10,2.17044e-14,3861.35,-77.5836], Tmin=(927.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.9892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = 'OCHCO(3676)',
    structure = SMILES('O=[C]C=O'),
    E0 = (-75.5464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,180,525.376,1512.41,1512.65,1513.37],'cm^-1')),
        HinderedRotor(inertia=(0.00619061,'amu*angstrom^2'), symmetry=1, barrier=(0.260399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25386,0.0154074,-1.14326e-05,4.29104e-09,-6.58698e-13,-9058.53,11.1539], Tmin=(100,'K'), Tmax=(1504.29,'K')), NASAPolynomial(coeffs=[6.43153,0.00695778,-3.00696e-06,5.56991e-10,-3.81281e-14,-10014.6,-5.47402], Tmin=(1504.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'ketene(1375)',
    structure = SMILES('C=C=O'),
    E0 = (-60.9858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48302,0.0083299,5.93417e-06,-1.33939e-08,5.73133e-12,-7313.53,6.51898], Tmin=(100,'K'), Tmax=(954.863,'K')), NASAPolynomial(coeffs=[5.88246,0.00580173,-1.91266e-06,3.35917e-10,-2.37398e-14,-8114.73,-6.74202], Tmin=(954.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.9858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[C]=C[O](6859)',
    structure = SMILES('[C]C=O'),
    E0 = (491.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53714,0.00730776,2.32451e-06,-8.52039e-09,3.75032e-12,59149.2,7.56143], Tmin=(100,'K'), Tmax=(1018.09,'K')), NASAPolynomial(coeffs=[6.47558,0.00262438,-8.8465e-07,2.00876e-10,-1.68063e-14,58195.3,-8.41394], Tmin=(1018.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJ3)"""),
)

species(
    label = '[O]C[C]=O(2320)',
    structure = SMILES('[O]C[C]=O'),
    E0 = (65.8932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,298.422],'cm^-1')),
        HinderedRotor(inertia=(0.460283,'amu*angstrom^2'), symmetry=1, barrier=(28.9433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01454,0.0123699,2.06297e-05,-3.93769e-08,1.66626e-11,7969.02,14.0082], Tmin=(100,'K'), Tmax=(965.917,'K')), NASAPolynomial(coeffs=[10.932,0.0028129,-6.04041e-07,1.77093e-10,-1.91371e-14,5355.8,-29.5241], Tmin=(965.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.8932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=OCOJ) + radical(CsCJ=O)"""),
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
    label = '[CH2]O[C]=C[O](11426)',
    structure = SMILES('[CH2]O[C]=C[O]'),
    E0 = (216.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,400.773,400.953,401.443],'cm^-1')),
        HinderedRotor(inertia=(0.189956,'amu*angstrom^2'), symmetry=1, barrier=(21.6878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18896,'amu*angstrom^2'), symmetry=1, barrier=(21.685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09104,0.028776,1.00821e-05,-4.80262e-08,2.47341e-11,26115.5,20.59], Tmin=(100,'K'), Tmax=(909.214,'K')), NASAPolynomial(coeffs=[17.1662,-0.00241134,3.5706e-06,-7.50851e-10,4.89516e-14,21922,-58.6903], Tmin=(909.214,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COCJ) + radical(C=CJO)"""),
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
    label = '[CH]=[C]OC[C]=O(2658)',
    structure = SMILES('[CH]=[C]OC[C]=O'),
    E0 = (389.415,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1685,370,3120,650,792.5,1650,295.126,295.365],'cm^-1')),
        HinderedRotor(inertia=(0.321176,'amu*angstrom^2'), symmetry=1, barrier=(19.8305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.320376,'amu*angstrom^2'), symmetry=1, barrier=(19.8356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.321556,'amu*angstrom^2'), symmetry=1, barrier=(19.8342,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3435.15,'J/mol'), sigma=(5.65948,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.56 K, Pc=43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53967,0.0471702,-4.06961e-05,1.0405e-08,1.667e-12,46930.7,23.9343], Tmin=(100,'K'), Tmax=(966.948,'K')), NASAPolynomial(coeffs=[15.0526,0.00677811,-2.09294e-06,3.75429e-10,-2.80054e-14,43592.5,-44.5523], Tmin=(966.948,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.415,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCJ=O) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'O=[C]COC1=CO1(11807)',
    structure = SMILES('O=[C]COC1=CO1'),
    E0 = (67.7969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.95168,0.0582611,-5.1932e-05,1.41939e-08,1.65913e-12,8271.74,22.0513], Tmin=(100,'K'), Tmax=(972.475,'K')), NASAPolynomial(coeffs=[18.1804,0.00666557,-2.07159e-06,3.89236e-10,-3.0313e-14,4009.64,-65.2716], Tmin=(972.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.7969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C=C1OCC1=O(11808)',
    structure = SMILES('[O]C=C1OCC1=O'),
    E0 = (-198.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62764,0.0268973,5.88805e-05,-1.12923e-07,4.94062e-11,-23817.5,18.3616], Tmin=(100,'K'), Tmax=(938.61,'K')), NASAPolynomial(coeffs=[23.732,-0.0029222,3.64726e-06,-6.14366e-10,2.81291e-14,-30802.9,-101.986], Tmin=(938.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(C=COJ)"""),
)

species(
    label = 'O=[C]COC=C=O(11167)',
    structure = SMILES('O=[C]COC=C=O'),
    E0 = (-124.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.336904,0.0627665,-6.67288e-05,3.25598e-08,-5.8629e-12,-14855.8,25.3974], Tmin=(100,'K'), Tmax=(1572.81,'K')), NASAPolynomial(coeffs=[19.952,0.00247767,6.90932e-07,-2.22907e-10,1.64316e-14,-19739.2,-74.0262], Tmin=(1572.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C=COC=C=O(11302)',
    structure = SMILES('[O]C=COC=C=O'),
    E0 = (-124.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01811,0.0513189,-1.99387e-05,-2.98534e-08,2.06899e-11,-14893.2,23.6764], Tmin=(100,'K'), Tmax=(911.124,'K')), NASAPolynomial(coeffs=[21.6465,-0.000988181,3.19519e-06,-6.98075e-10,4.56575e-14,-20240.1,-82.6372], Tmin=(911.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=COJ)"""),
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
    label = '[O]C=[C]OC=C=O(11379)',
    structure = SMILES('[O]C=[C]OC=C=O'),
    E0 = (114.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10413,'amu*angstrom^2'), symmetry=1, barrier=(25.3861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10495,'amu*angstrom^2'), symmetry=1, barrier=(25.4049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735088,0.058878,-6.71021e-05,3.53139e-08,-6.86317e-12,13949.3,27.8211], Tmin=(100,'K'), Tmax=(1447.26,'K')), NASAPolynomial(coeffs=[18.0882,0.00202857,1.03083e-06,-3.1427e-10,2.42441e-14,9857.18,-59.0971], Tmin=(1447.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'O=[C]CO[C]=C=O(11809)',
    structure = SMILES('O=[C]CO[C]=C=O'),
    E0 = (115.002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1685,370,2120,512.5,787.5,344.99,345.016,345.124],'cm^-1')),
        HinderedRotor(inertia=(0.216148,'amu*angstrom^2'), symmetry=1, barrier=(18.2791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216327,'amu*angstrom^2'), symmetry=1, barrier=(18.2765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216477,'amu*angstrom^2'), symmetry=1, barrier=(18.2791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23019,0.0566721,-6.70204e-05,3.79281e-08,-8.22628e-12,13935,25.306], Tmin=(100,'K'), Tmax=(1141.43,'K')), NASAPolynomial(coeffs=[15.1303,0.0079606,-3.00647e-06,5.39827e-10,-3.73544e-14,10761.8,-43.593], Tmin=(1141.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsCJ=O) + radical(C=CJO)"""),
)

species(
    label = 'C#COC[C]=O(2656)',
    structure = SMILES('C#COC[C]=O'),
    E0 = (111.902,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.26218,'amu*angstrom^2'), symmetry=1, barrier=(29.0201,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26003,'amu*angstrom^2'), symmetry=1, barrier=(28.9706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26381,'amu*angstrom^2'), symmetry=1, barrier=(29.0576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40692,0.0474914,-3.31409e-05,-2.19472e-09,7.10298e-12,13560.7,18.3356], Tmin=(100,'K'), Tmax=(956.019,'K')), NASAPolynomial(coeffs=[17.1789,0.00421919,-8.91369e-07,1.72896e-10,-1.61321e-14,9506.83,-62.4756], Tmin=(956.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.902,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][C]=O(1376)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,623.763,623.847],'cm^-1')),
        HinderedRotor(inertia=(0.000434039,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44947,0.00889162,4.91274e-06,-1.14995e-08,4.62261e-12,19378.1,9.99239], Tmin=(100,'K'), Tmax=(1030.77,'K')), NASAPolynomial(coeffs=[6.10651,0.00625992,-2.43247e-06,4.78594e-10,-3.54632e-14,18422.4,-4.88573], Tmin=(1030.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]=C[O](9592)',
    structure = SMILES('[O][CH][C]=O'),
    E0 = (204.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,343.122],'cm^-1')),
        HinderedRotor(inertia=(0.528453,'amu*angstrom^2'), symmetry=1, barrier=(44.0956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.99411,0.0153541,4.13359e-06,-1.97935e-08,9.27968e-12,24622.8,13.8337], Tmin=(100,'K'), Tmax=(987.313,'K')), NASAPolynomial(coeffs=[10.2587,0.00223056,-7.0477e-07,2.03549e-10,-2.00526e-14,22393.5,-25.1462], Tmin=(987.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=OCOJ) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]=CO[C]=C[O](11389)',
    structure = SMILES('[O]C=[C]O[CH][C]=O'),
    E0 = (268.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07527,'amu*angstrom^2'), symmetry=1, barrier=(24.7225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07753,'amu*angstrom^2'), symmetry=1, barrier=(24.7746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07641,'amu*angstrom^2'), symmetry=1, barrier=(24.7487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.247163,0.0719586,-8.70005e-05,4.56261e-08,-8.61352e-12,32514.6,30.6943], Tmin=(100,'K'), Tmax=(1540.81,'K')), NASAPolynomial(coeffs=[23.5694,-0.00618032,4.94705e-06,-1.02716e-09,7.07614e-14,27111.3,-88.2201], Tmin=(1540.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]=[C]OC[C]=O(11810)',
    structure = SMILES('[O][C]=[C]OC[C]=O'),
    E0 = (314.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1670,1700,300,440,396.451,396.494,396.527,396.555],'cm^-1')),
        HinderedRotor(inertia=(0.177993,'amu*angstrom^2'), symmetry=1, barrier=(19.8612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177991,'amu*angstrom^2'), symmetry=1, barrier=(19.8619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178008,'amu*angstrom^2'), symmetry=1, barrier=(19.861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848984,0.0582439,-6.58828e-05,3.44247e-08,-6.71863e-12,37976.6,30.4119], Tmin=(100,'K'), Tmax=(1391.17,'K')), NASAPolynomial(coeffs=[18.1606,0.00280663,-4.12421e-09,-7.05398e-11,6.0402e-15,33707.8,-56.8528], Tmin=(1391.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CsCJ=O) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]OC=C[O](11297)',
    structure = SMILES('[O]C=[C]OC=C[O]'),
    E0 = (74.8979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.971118,0.0489915,-5.32961e-06,-5.06611e-08,2.95397e-11,9133.72,27.5802], Tmin=(100,'K'), Tmax=(907.716,'K')), NASAPolynomial(coeffs=[23.9431,-0.00489397,5.47946e-06,-1.13936e-09,7.50726e-14,3012.85,-91.7659], Tmin=(907.716,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.8979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'O=[C]CO[C]=[C]O(11811)',
    structure = SMILES('O=[C]CO[C]=[C]O'),
    E0 = (173.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1855,455,950,1670,1700,300,440,240.402,240.402,240.402],'cm^-1')),
        HinderedRotor(inertia=(0.470619,'amu*angstrom^2'), symmetry=1, barrier=(19.3006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470619,'amu*angstrom^2'), symmetry=1, barrier=(19.3006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470619,'amu*angstrom^2'), symmetry=1, barrier=(19.3006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470619,'amu*angstrom^2'), symmetry=1, barrier=(19.3006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.184815,0.0682351,-8.00796e-05,4.25954e-08,-8.30261e-12,20990.5,31.3829], Tmin=(100,'K'), Tmax=(1456.11,'K')), NASAPolynomial(coeffs=[20.9,-0.000203997,2.30349e-06,-5.62685e-10,4.1154e-14,16180.5,-72.1415], Tmin=(1456.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CsCJ=O) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'O=[C][CH]OC[C]=O(11150)',
    structure = SMILES('[O][C]=COC[C]=O'),
    E0 = (74.9892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1855,455,950,1685,370,205.135,205.135,205.135,205.135],'cm^-1')),
        HinderedRotor(inertia=(0.813479,'amu*angstrom^2'), symmetry=1, barrier=(24.2913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.81348,'amu*angstrom^2'), symmetry=1, barrier=(24.2913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813479,'amu*angstrom^2'), symmetry=1, barrier=(24.2913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995855,0.0523372,-2.4758e-05,-2.24589e-08,1.707e-11,9139.99,26.753], Tmin=(100,'K'), Tmax=(927.568,'K')), NASAPolynomial(coeffs=[21.2354,0.000482909,1.80938e-06,-3.79382e-10,2.17064e-14,3861.31,-77.5841], Tmin=(927.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.9892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CsCJ=O) + radical(C=CJO)"""),
)

species(
    label = '[O]C=CO[CH][C]=O(11291)',
    structure = SMILES('[O]C=CO[CH][C]=O'),
    E0 = (29.1724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416979,0.0599772,-2.46439e-05,-3.8972e-08,2.71516e-11,3655.28,25.1767], Tmin=(100,'K'), Tmax=(910.45,'K')), NASAPolynomial(coeffs=[27.8973,-0.0101551,7.53526e-06,-1.49074e-09,9.7606e-14,-3445.81,-116.339], Tmin=(910.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.1724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]=[C]OCC=O(11812)',
    structure = SMILES('[O][C]=[C]OCC=O'),
    E0 = (160.062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1670,1700,300,440,395.704,395.705,395.705,395.705],'cm^-1')),
        HinderedRotor(inertia=(0.17156,'amu*angstrom^2'), symmetry=1, barrier=(19.0627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17156,'amu*angstrom^2'), symmetry=1, barrier=(19.0627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17156,'amu*angstrom^2'), symmetry=1, barrier=(19.0627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621074,0.0598731,-6.28993e-05,3.08781e-08,-5.65869e-12,19385.1,30.5579], Tmin=(100,'K'), Tmax=(1510.64,'K')), NASAPolynomial(coeffs=[18.6672,0.00451916,-4.18762e-07,-1.25467e-11,2.45837e-15,14796.5,-61.0902], Tmin=(1510.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'O=[C][CH]O[C]=CO(11813)',
    structure = SMILES('O=[C][CH]O[C]=CO'),
    E0 = (127.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01434,'amu*angstrom^2'), symmetry=1, barrier=(23.3218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01491,'amu*angstrom^2'), symmetry=1, barrier=(23.3347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01293,'amu*angstrom^2'), symmetry=1, barrier=(23.2892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01372,'amu*angstrom^2'), symmetry=1, barrier=(23.3073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.913834,0.0819787,-0.000101295,5.39176e-08,-1.0245e-11,15528.6,31.6743], Tmin=(100,'K'), Tmax=(1553.21,'K')), NASAPolynomial(coeffs=[25.9597,-0.00869698,7.00654e-06,-1.46653e-09,1.0184e-13,9770.07,-101.472], Tmin=(1553.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = 'O=C1CO[C][CH]O1(11814)',
    structure = SMILES('O=C1CO[C][CH]O1'),
    E0 = (131.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62796,0.0356439,1.59871e-05,-4.95481e-08,2.10045e-11,15965.3,17.2488], Tmin=(100,'K'), Tmax=(1035.26,'K')), NASAPolynomial(coeffs=[17.7025,0.0125107,-6.96591e-06,1.59762e-09,-1.28005e-13,10548.5,-70.9457], Tmin=(1035.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(CCsJOC(O)) + radical(CH2_triplet)"""),
)

species(
    label = '[O]C1[C]OCC1=O(11815)',
    structure = SMILES('[O]C1[C]OCC1=O'),
    E0 = (238.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08401,0.0255365,3.49863e-05,-7.08719e-08,3.05665e-11,28776.4,22.5913], Tmin=(100,'K'), Tmax=(954.571,'K')), NASAPolynomial(coeffs=[16.3369,0.00736077,-1.74249e-06,3.77321e-10,-3.55924e-14,24162.3,-55.4234], Tmin=(954.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + ring(Cyclopentane) + radical(C=OCOJ) + radical(CH2_triplet)"""),
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
    E0 = (74.9892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (557.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (710.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (908.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (77.918,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (83.2735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (97.2901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (163.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (122.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (337.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (334.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (377.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (115.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (365.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (480.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (526.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (208.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (331.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (244.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (216.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (193.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (179.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (131.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (238.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]C=[C]OC[C]=O(11163)'],
    products = ['OCHCO(3676)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=C[O](6859)', '[O]C[C]=O(2320)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(2355)', '[CH2]O[C]=C[O](11426)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=[C]OC[C]=O(2658)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[O]C=[C]OC[C]=O(11163)'],
    products = ['O=[C]COC1=CO1(11807)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O]C=[C]OC[C]=O(11163)'],
    products = ['[O]C=C1OCC1=O(11808)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O]C=[C]OC[C]=O(11163)'],
    products = ['O=[C]COC=C=O(11167)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]C=[C]OC[C]=O(11163)'],
    products = ['[O]C=COC=C=O(11302)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CO(2039)', '[CH2]O[C]=C[O](11426)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(754.069,'m^3/(mol*s)'), n=1.0822, Ea=(24.7919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [COm;C_pri_rad] for rate rule [COm;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[O]C=[C]OC=C=O(11379)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ck;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', 'O=[C]CO[C]=C=O(11809)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O(4)', 'C#COC[C]=O(2656)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.30864e+06,'m^3/(mol*s)'), n=-0.19959, Ea=(22.3126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_atom_triplet] + [Ct_Ct;YJ] for rate rule [Ct_Ct;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['OCHCO(3676)', '[CH2][C]=O(1376)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(29.883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C]=C[O](9592)', '[CH2][C]=O(1376)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[O][C]=CO[C]=C[O](11389)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[O][C]=[C]OC[C]=O(11810)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O]C=[C]OC[C]=O(11163)'],
    products = ['[O]C=[C]OC=C[O](11297)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00141351,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O=[C]CO[C]=[C]O(11811)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.25466e+06,'s^-1'), n=1.80084, Ea=(158.227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;O_H_out] + [R2H_S;Cd_rad_out;XH_out] for rate rule [R2H_S;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=[C][CH]OC[C]=O(11150)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C=[C]OC[C]=O(11163)'],
    products = ['[O]C=CO[CH][C]=O(11291)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_Cd;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_Cd;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][C]=[C]OCC=O(11812)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_1;Cd_rad_out;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O=[C][CH]O[C]=CO(11813)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.10205e+06,'s^-1'), n=1.54368, Ea=(52.1315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;XH_out] for rate rule [R5HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O]C=[C]OC[C]=O(11163)'],
    products = ['O=C1CO[C][CH]O1(11814)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(56.9244,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;carbonyl_intra_H;radadd_intra_CO] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_CO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 54.7 to 56.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]C=[C]OC[C]=O(11163)'],
    products = ['[O]C1[C]OCC1=O(11815)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.76201e+09,'s^-1'), n=0.626373, Ea=(163.572,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra_H;radadd_intra_CO]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 159.7 to 163.6 kJ/mol to match endothermicity of reaction."""),
)

network(
    label = 'PDepNetwork #2831',
    isomers = [
        '[O]C=[C]OC[C]=O(11163)',
    ],
    reactants = [
        ('OCHCO(3676)', 'ketene(1375)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2831',
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


species(
    label = 'C=[C]OO[C]=C[O](11161)',
    structure = SMILES('C=[C]OO[C]=C[O]'),
    E0 = (475.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1670,1700,300,440,316.591,316.603],'cm^-1')),
        HinderedRotor(inertia=(0.205825,'amu*angstrom^2'), symmetry=1, barrier=(14.6398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.68203,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205822,'amu*angstrom^2'), symmetry=1, barrier=(14.6401,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35049,0.0561293,-6.14008e-05,3.35657e-08,-7.2047e-12,57263.1,31.0673], Tmin=(100,'K'), Tmax=(1137.73,'K')), NASAPolynomial(coeffs=[13.2873,0.0141615,-6.06888e-06,1.1428e-09,-8.01219e-14,54547,-28.0611], Tmin=(1137.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
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
    label = 'C=[C]O[O](1354)',
    structure = SMILES('C=[C]O[O]'),
    E0 = (337.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.00100972,'amu*angstrom^2'), symmetry=1, barrier=(11.4644,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.15943,0.0239683,-4.52732e-05,4.81856e-08,-1.85891e-11,40635.7,16.5558], Tmin=(100,'K'), Tmax=(873.282,'K')), NASAPolynomial(coeffs=[2.07161,0.0153864,-7.23324e-06,1.359e-09,-9.18757e-14,41343,24.6179], Tmin=(873.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]O[O](10707)',
    structure = SMILES('[O]C=[C]O[O]'),
    E0 = (282.27,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.622504,'amu*angstrom^2'), symmetry=1, barrier=(14.3126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0275,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28922,0.0329953,-4.0169e-05,2.30557e-08,-4.92424e-12,34014.9,20.5563], Tmin=(100,'K'), Tmax=(1295.31,'K')), NASAPolynomial(coeffs=[11.0046,0.00208215,2.60805e-07,-1.36399e-10,1.20047e-14,32092.6,-22.4503], Tmin=(1295.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.27,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(C=CJO)"""),
)

species(
    label = 'H2CC(T)(1341)',
    structure = SMILES('[C]=C'),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9483,-0.00181651,2.07124e-05,-2.3669e-08,8.45456e-12,72206.8,5.3155], Tmin=(100,'K'), Tmax=(953.668,'K')), NASAPolynomial(coeffs=[4.20274,0.00473419,-1.57302e-06,2.859e-10,-2.08638e-14,71811.9,2.2838], Tmin=(953.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=[C]OO[C]=C(2665)',
    structure = SMILES('[CH]=[C]OO[C]=C'),
    E0 = (789.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1670,1700,300,440,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.251919,'amu*angstrom^2'), symmetry=1, barrier=(5.79211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251914,'amu*angstrom^2'), symmetry=1, barrier=(5.792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00699,'amu*angstrom^2'), symmetry=1, barrier=(23.1526,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3211.77,'J/mol'), sigma=(5.41486,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=501.67 K, Pc=45.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88878,0.0510668,-7.80177e-05,6.80613e-08,-2.37876e-11,95054.1,28.267], Tmin=(100,'K'), Tmax=(806.744,'K')), NASAPolynomial(coeffs=[6.50612,0.021481,-1.05653e-05,2.03861e-09,-1.41573e-13,94526.9,8.3323], Tmin=(806.744,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(789.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OOC1=CO1(11844)',
    structure = SMILES('C=[C]OOC1=CO1'),
    E0 = (468.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47395,0.0599327,-8.05553e-05,5.91293e-08,-1.76581e-11,56387.8,25.7745], Tmin=(100,'K'), Tmax=(813.463,'K')), NASAPolynomial(coeffs=[9.13617,0.022253,-1.10704e-05,2.17965e-09,-1.54617e-13,55141.3,-9.60885], Tmin=(813.463,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.112,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OOC1=C[O](11845)',
    structure = SMILES('C=C1OOC1=C[O]'),
    E0 = (111.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63238,0.042926,-2.03256e-05,-7.53193e-09,6.43746e-12,13467.9,20.5978], Tmin=(100,'K'), Tmax=(1026.72,'K')), NASAPolynomial(coeffs=[13.9939,0.0138644,-5.76831e-06,1.1323e-09,-8.34927e-14,9922.95,-44.2672], Tmin=(1026.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(C=COJ)"""),
)

species(
    label = 'C=[C]OOC=C=O(11166)',
    structure = SMILES('C=[C]OOC=C=O'),
    E0 = (275.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5497,0.0565451,-6.88166e-05,4.43979e-08,-1.15429e-11,33229.7,26.6256], Tmin=(100,'K'), Tmax=(933.341,'K')), NASAPolynomial(coeffs=[10.214,0.0194126,-9.13972e-06,1.77181e-09,-1.2523e-13,31612.4,-14.5768], Tmin=(933.341,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=CJO)"""),
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
    label = 'C#COO[C]=C[O](11390)',
    structure = SMILES('C#COO[C]=C[O]'),
    E0 = (467.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,1685,370,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(1.473,'amu*angstrom^2'), symmetry=1, barrier=(33.8672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47206,'amu*angstrom^2'), symmetry=1, barrier=(33.8455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4724,'amu*angstrom^2'), symmetry=1, barrier=(33.8533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0373,0.0618195,-7.42739e-05,4.24262e-08,-9.31277e-12,56285.9,26.6709], Tmin=(100,'K'), Tmax=(1125.19,'K')), NASAPolynomial(coeffs=[15.8776,0.00906208,-3.94143e-06,7.54082e-10,-5.37352e-14,52946.3,-46.6753], Tmin=(1125.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OO[C]=C=O(11846)',
    structure = SMILES('C=[C]OO[C]=C=O'),
    E0 = (515.317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1670,1700,300,440,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.18339,'amu*angstrom^2'), symmetry=1, barrier=(4.21649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180648,'amu*angstrom^2'), symmetry=1, barrier=(4.15344,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.716261,'amu*angstrom^2'), symmetry=1, barrier=(16.4683,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57285,0.0606681,-0.000104809,9.63724e-08,-3.40935e-11,62058.7,29.6605], Tmin=(100,'K'), Tmax=(845.662,'K')), NASAPolynomial(coeffs=[6.43103,0.0229077,-1.16133e-05,2.23367e-09,-1.53397e-13,61765.6,10.1619], Tmin=(845.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C#COO[C]=C(2664)',
    structure = SMILES('C#COO[C]=C'),
    E0 = (534.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(1.34249,'amu*angstrom^2'), symmetry=1, barrier=(30.8664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35112,'amu*angstrom^2'), symmetry=1, barrier=(31.0648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34311,'amu*angstrom^2'), symmetry=1, barrier=(30.8808,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88133,0.0503015,-6.29738e-05,4.33221e-08,-1.22468e-11,64347.3,22.3829], Tmin=(100,'K'), Tmax=(853.703,'K')), NASAPolynomial(coeffs=[8.31092,0.0201755,-1.00404e-05,1.98546e-09,-1.41558e-13,63249.5,-7.61908], Tmin=(853.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsCt) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=CJO)"""),
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
    label = '[CH]=[C]OO[C]=C[O](11376)',
    structure = SMILES('[CH]=[C]OO[C]=C[O]'),
    E0 = (722.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,247.747],'cm^-1')),
        HinderedRotor(inertia=(0.356526,'amu*angstrom^2'), symmetry=1, barrier=(15.6361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.561425,'amu*angstrom^2'), symmetry=1, barrier=(24.2541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354825,'amu*angstrom^2'), symmetry=1, barrier=(15.6453,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33121,0.0589737,-7.55537e-05,4.75801e-08,-1.16579e-11,86980.6,31.5425], Tmin=(100,'K'), Tmax=(1004.35,'K')), NASAPolynomial(coeffs=[13.1186,0.0120283,-5.44042e-06,1.04022e-09,-7.32581e-14,84612.9,-25.3759], Tmin=(1004.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OO[C]=[C][O](11847)',
    structure = SMILES('C=[C]OO[C][C]=O'),
    E0 = (709.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,303.974,303.991],'cm^-1')),
        HinderedRotor(inertia=(0.251265,'amu*angstrom^2'), symmetry=1, barrier=(16.5577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.420217,'amu*angstrom^2'), symmetry=1, barrier=(27.3843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138714,'amu*angstrom^2'), symmetry=1, barrier=(9.16146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.415767,'amu*angstrom^2'), symmetry=1, barrier=(27.386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3949,0.0624431,-8.92056e-05,6.56153e-08,-1.93493e-11,85367.4,27.4617], Tmin=(100,'K'), Tmax=(826.556,'K')), NASAPolynomial(coeffs=[10.3558,0.0190807,-1.05186e-05,2.15378e-09,-1.55991e-13,83885.9,-14.0633], Tmin=(826.556,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(709.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=COO[C]=C[O](11299)',
    structure = SMILES('[CH]=COO[C]=C[O]'),
    E0 = (482.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.06022,'amu*angstrom^2'), symmetry=1, barrier=(24.3766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05775,'amu*angstrom^2'), symmetry=1, barrier=(24.3197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05705,'amu*angstrom^2'), symmetry=1, barrier=(24.3037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.607336,0.0632666,-6.9767e-05,3.6041e-08,-7.00959e-12,58181.9,31.0117], Tmin=(100,'K'), Tmax=(1372.08,'K')), NASAPolynomial(coeffs=[18.9685,0.00502928,-9.51623e-07,1.03385e-10,-5.76616e-15,53586.6,-61.7628], Tmin=(1372.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OO[C]=[C]O(11848)',
    structure = SMILES('C=[C]OO[C]=[C]O'),
    E0 = (573.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,314.673],'cm^-1')),
        HinderedRotor(inertia=(0.27651,'amu*angstrom^2'), symmetry=1, barrier=(19.3329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16369,'amu*angstrom^2'), symmetry=1, barrier=(11.4461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163647,'amu*angstrom^2'), symmetry=1, barrier=(11.4448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163071,'amu*angstrom^2'), symmetry=1, barrier=(11.449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15938,0.0646914,-9.09531e-05,6.5058e-08,-1.82079e-11,69086.7,33.475], Tmin=(100,'K'), Tmax=(880.522,'K')), NASAPolynomial(coeffs=[12.0525,0.0152042,-6.64643e-06,1.22447e-09,-8.33453e-14,67168.4,-17.6917], Tmin=(880.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(573.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OO[CH][C]=O(11148)',
    structure = SMILES('C=[C]OO[CH][C]=O'),
    E0 = (386.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70567,0.0573674,-6.86296e-05,3.63059e-08,-1.18264e-12,46565.9,27.4075], Tmin=(100,'K'), Tmax=(562.806,'K')), NASAPolynomial(coeffs=[7.38954,0.0268941,-1.38594e-05,2.7569e-09,-1.96255e-13,45768.9,1.85713], Tmin=(562.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = 'C=COO[C]=[C][O](11849)',
    structure = SMILES('C=COO[C][C]=O'),
    E0 = (469.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04052,'amu*angstrom^2'), symmetry=1, barrier=(23.9236,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0411,'amu*angstrom^2'), symmetry=1, barrier=(23.937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04112,'amu*angstrom^2'), symmetry=1, barrier=(23.9375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04035,'amu*angstrom^2'), symmetry=1, barrier=(23.9197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04193,0.0625391,-6.94796e-05,3.68215e-08,-7.5944e-12,56552.2,25.5877], Tmin=(100,'K'), Tmax=(1184.48,'K')), NASAPolynomial(coeffs=[15.8064,0.0126792,-6.33804e-06,1.28312e-09,-9.35596e-14,53054.6,-48.1421], Tmin=(1184.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=[C]OOC=C[O](11850)',
    structure = SMILES('[CH]=[C]OO[CH]C=O'),
    E0 = (478.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.0772868,'amu*angstrom^2'), symmetry=1, barrier=(1.77698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.59592,'amu*angstrom^2'), symmetry=1, barrier=(36.6934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.06212,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29353,0.0648156,-9.63724e-05,8.10379e-08,-2.78527e-11,57698.5,28.5617], Tmin=(100,'K'), Tmax=(744.305,'K')), NASAPolynomial(coeffs=[8.15805,0.0257386,-1.32151e-05,2.60878e-09,-1.84242e-13,56737.2,-2.1218], Tmin=(744.305,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OO[C]=CO(11851)',
    structure = SMILES('[CH]=[C]OO[C]=CO'),
    E0 = (580.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.81469,'amu*angstrom^2'), symmetry=1, barrier=(18.7313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815079,'amu*angstrom^2'), symmetry=1, barrier=(18.7403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.814411,'amu*angstrom^2'), symmetry=1, barrier=(18.7249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.814094,'amu*angstrom^2'), symmetry=1, barrier=(18.7176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.699037,0.0686276,-8.87104e-05,5.45155e-08,-1.27397e-11,69993,32.3955], Tmin=(100,'K'), Tmax=(1065.16,'K')), NASAPolynomial(coeffs=[17.074,0.00713607,-2.11739e-06,3.19436e-10,-1.98427e-14,66504.5,-47.6381], Tmin=(1065.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(580.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=C1O[CH][C]OO1(11852)',
    structure = SMILES('[CH2][C]1OC=[C]OO1'),
    E0 = (436.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61122,0.0352742,2.10131e-05,-7.07682e-08,3.54569e-11,52661.2,25.6986], Tmin=(100,'K'), Tmax=(899.487,'K')), NASAPolynomial(coeffs=[20.357,-0.00149457,4.62968e-06,-1.03718e-09,7.02141e-14,47404,-73.2306], Tmin=(899.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(124trioxene) + radical(Cs_P) + radical(CJCOOH) + radical(C=CJO)"""),
)

species(
    label = 'C=C1OO[C]C1[O](11853)',
    structure = SMILES('C=C1OO[C]C1[O]'),
    E0 = (498.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90089,0.0288531,3.12863e-05,-6.84469e-08,2.96013e-11,60023.4,22.1841], Tmin=(100,'K'), Tmax=(967.409,'K')), NASAPolynomial(coeffs=[17.3227,0.0078763,-2.5336e-06,5.79411e-10,-5.17716e-14,55037.3,-62.0546], Tmin=(967.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CH2_triplet)"""),
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
    E0 = (475.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (829.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (882.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1309.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (478.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (483.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (497.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (693.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (734.883,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (799.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (475.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (934.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (920.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (588.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (731.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (679.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (498.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (788.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (732.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (533.722,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (569.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OO[C]=C[O](11161)'],
    products = ['OCHCO(3676)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[C]=C[O](6859)', 'C=[C]O[O](1354)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]C=[C]O[O](10707)', 'H2CC(T)(1341)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=[C]OO[C]=C(2665)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]OO[C]=C[O](11161)'],
    products = ['C=[C]OOC1=CO1(11844)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C]OO[C]=C[O](11161)'],
    products = ['C=C1OOC1=C[O](11845)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]OO[C]=C[O](11161)'],
    products = ['C=[C]OOC=C=O(11166)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C#COO[C]=C[O](11390)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C=[C]OO[C]=C=O(11846)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['O(4)', 'C#COO[C]=C(2664)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.30864e+06,'m^3/(mol*s)'), n=-0.19959, Ea=(22.3126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_atom_triplet] + [Ct_Ct;YJ] for rate rule [Ct_Ct;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=C[O](9592)', '[CH2][C]=O(1376)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(109.999,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination
Ea raised from 105.2 to 110.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]=[C]OO[C]=C[O](11376)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', 'C=[C]OO[C]=[C][O](11847)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=COO[C]=C[O](11299)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]OO[C]=[C]O(11848)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.25466e+06,'s^-1'), n=1.80084, Ea=(158.227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;O_H_out] + [R2H_S;Cd_rad_out;XH_out] for rate rule [R2H_S;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C]OO[C]=C[O](11161)'],
    products = ['C=[C]OO[CH][C]=O(11148)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=COO[C]=[C][O](11849)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.658e+09,'s^-1'), n=0.699, Ea=(29.5516,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]OOC=C[O](11850)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.24929e+07,'s^-1'), n=1.719, Ea=(309.087,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_doubleC] for rate rule [R5HJ_1;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]OO[C]=CO(11851)'],
    products = ['C=[C]OO[C]=C[O](11161)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]OO[C]=C[O](11161)'],
    products = ['C=C1O[CH][C]OO1(11852)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.19745e+11,'s^-1'), n=0.440371, Ea=(58.418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6;multiplebond_intra;radadd_intra_cddouble] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=[C]OO[C]=C[O](11161)'],
    products = ['C=C1OO[C]C1[O](11853)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.83276e+09,'s^-1'), n=0.690186, Ea=(94.193,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_cddouble] + [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

network(
    label = 'PDepNetwork #2829',
    isomers = [
        'C=[C]OO[C]=C[O](11161)',
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
    label = 'PDepNetwork #2829',
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


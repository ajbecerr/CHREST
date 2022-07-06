species(
    label = '[CH2]C(C=C)C([O])=C[O](11262)',
    structure = SMILES('[CH2]C(C=C)C([O])=C[O]'),
    E0 = (116.179,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,386.932,386.934,386.935,386.959],'cm^-1')),
        HinderedRotor(inertia=(0.14209,'amu*angstrom^2'), symmetry=1, barrier=(15.0966,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142096,'amu*angstrom^2'), symmetry=1, barrier=(15.0974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142088,'amu*angstrom^2'), symmetry=1, barrier=(15.0967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.449233,0.0836711,-8.67254e-05,4.49871e-08,-8.88727e-12,14145.4,32.9887], Tmin=(100,'K'), Tmax=(1373.96,'K')), NASAPolynomial(coeffs=[20.4882,0.0149466,-3.21414e-06,3.50409e-10,-1.6487e-14,9125.3,-72.0051], Tmin=(1373.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl)"""),
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
    label = 'C=CC[CH]C([O])=C[O](12006)',
    structure = SMILES('C=CC[CH]C([O])=C[O]'),
    E0 = (119.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44621,0.062241,-1.65539e-05,-3.55796e-08,2.17479e-11,14517.3,33.1206], Tmin=(100,'K'), Tmax=(937.352,'K')), NASAPolynomial(coeffs=[21.185,0.0138092,-3.1694e-06,5.0365e-10,-3.844e-14,8869.15,-74.9799], Tmin=(937.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJCO)"""),
)

species(
    label = 'C=C[CH]CC([O])=C[O](11260)',
    structure = SMILES('[CH2]C=CCC([O])=C[O]'),
    E0 = (60.1396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315737,0.0648577,-1.99561e-05,-3.38401e-08,2.14782e-11,7380.62,30.9617], Tmin=(100,'K'), Tmax=(937.203,'K')), NASAPolynomial(coeffs=[21.8558,0.0136192,-3.08064e-06,4.86701e-10,-3.7325e-14,1555.93,-81.0942], Tmin=(937.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.1396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1[CH]CC([O])[C]1[O](12007)',
    structure = SMILES('[CH2]C1[CH]CC([O])[C]1[O]'),
    E0 = (567.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13279,0.0514275,-6.37411e-06,-2.91632e-08,1.5243e-11,68372.6,31.2653], Tmin=(100,'K'), Tmax=(971.794,'K')), NASAPolynomial(coeffs=[13.9706,0.0258893,-9.09909e-06,1.61767e-09,-1.13298e-13,64588.2,-36.9359], Tmin=(971.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopentane) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(Cs_S) + radical(C2CsJOH) + radical(Isobutyl)"""),
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
    label = 'C=C[CH]C([O])=C[O](12008)',
    structure = SMILES('[CH2]C=CC([O])=C[O]'),
    E0 = (32.5332,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.50515,'amu*angstrom^2'), symmetry=1, barrier=(34.6062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49601,'amu*angstrom^2'), symmetry=1, barrier=(34.3963,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.980366,0.0515229,-7.89604e-06,-4.26113e-08,2.51092e-11,4035.52,24.2124], Tmin=(100,'K'), Tmax=(906.957,'K')), NASAPolynomial(coeffs=[20.358,0.00614205,8.68621e-07,-3.26837e-10,2.21209e-14,-1127.91,-76.4691], Tmin=(906.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.5332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CC=CCJ)"""),
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
    label = '[CH]=C([O])C([CH2])C=C(5412)',
    structure = SMILES('[CH]=C([O])C([CH2])C=C'),
    E0 = (430.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,386.646,386.833],'cm^-1')),
        HinderedRotor(inertia=(0.0866711,'amu*angstrom^2'), symmetry=1, barrier=(9.19563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0865084,'amu*angstrom^2'), symmetry=1, barrier=(9.20059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868164,'amu*angstrom^2'), symmetry=1, barrier=(9.2058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3817.52,'J/mol'), sigma=(6.40826,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.29 K, Pc=32.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.649929,0.071954,-7.97361e-05,4.80514e-08,-1.15285e-11,51911.9,28.1778], Tmin=(100,'K'), Tmax=(1019.49,'K')), NASAPolynomial(coeffs=[13.0047,0.0234797,-8.41478e-06,1.41289e-09,-9.17639e-14,49392.8,-31.6653], Tmin=(1019.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = '[CH]C(C=C)C([O])=C[O](12009)',
    structure = SMILES('[CH]C(C=C)C([O])=C[O]'),
    E0 = (359.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.427692,0.0839287,-8.97136e-05,4.67562e-08,-9.27019e-12,43386.1,32.1359], Tmin=(100,'K'), Tmax=(1338.86,'K')), NASAPolynomial(coeffs=[21.7491,0.0120373,-2.85554e-06,3.62482e-10,-2.02048e-14,37952.9,-79.4394], Tmin=(1338.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=CC1COC1=C[O](12010)',
    structure = SMILES('C=CC1COC1=C[O]'),
    E0 = (-50.8657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.517907,0.0541999,2.29896e-05,-8.94937e-08,4.50536e-11,-5971.33,21.7109], Tmin=(100,'K'), Tmax=(902.77,'K')), NASAPolynomial(coeffs=[24.904,0.00627254,2.7261e-06,-7.58676e-10,5.1559e-14,-12824.3,-107.013], Tmin=(902.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-50.8657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C=C)C1=COO1(12011)',
    structure = SMILES('[CH2]C(C=C)C1=COO1'),
    E0 = (356.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.68312,0.062749,-3.52205e-05,-1.36755e-09,5.60772e-12,42976.5,29.2242], Tmin=(100,'K'), Tmax=(1012.82,'K')), NASAPolynomial(coeffs=[15.5904,0.0245419,-9.24391e-06,1.68116e-09,-1.17819e-13,38896.8,-48.1181], Tmin=(1012.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.26,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1COC=C1[O](12012)',
    structure = SMILES('C=CC1COC=C1[O]'),
    E0 = (-133.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.89216,0.0485951,2.58872e-05,-8.54955e-08,4.23357e-11,-15934.9,20.6944], Tmin=(100,'K'), Tmax=(896.451,'K')), NASAPolynomial(coeffs=[21.6341,0.00965302,1.34544e-06,-5.35467e-10,3.86289e-14,-21807.8,-89.1204], Tmin=(896.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-133.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(=C)C(O)=C[O](12013)',
    structure = SMILES('[CH2]C(C=C)=C(O)C=O'),
    E0 = (-122.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.343363,0.0825959,-6.72997e-05,1.48903e-08,3.54327e-12,-14531.2,26.21], Tmin=(100,'K'), Tmax=(985.852,'K')), NASAPolynomial(coeffs=[23.2716,0.0141955,-4.93949e-06,9.27807e-10,-6.91289e-14,-20519.6,-94.1388], Tmin=(985.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)C(O)=C=O(12014)',
    structure = SMILES('[CH2]C(C=C)C(O)=C=O'),
    E0 = (18.6703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.628862,0.0973998,-0.000121294,7.49496e-08,-1.7821e-11,2416.05,29.3277], Tmin=(100,'K'), Tmax=(1042.61,'K')), NASAPolynomial(coeffs=[20.3664,0.0168511,-5.40837e-06,8.50225e-10,-5.32808e-14,-1961.92,-72.8382], Tmin=(1042.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.6703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-CsCd(CCO)H) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O)"""),
)

species(
    label = 'C=CC(=C)C([O])=CO(12015)',
    structure = SMILES('C=CC(=C)C([O])=CO'),
    E0 = (-120.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.101517,0.0693496,-2.64561e-05,-3.40603e-08,2.39388e-11,-14287.5,25.6298], Tmin=(100,'K'), Tmax=(910.153,'K')), NASAPolynomial(coeffs=[23.8116,0.00923744,-5.10312e-08,-1.76529e-10,1.19172e-14,-20429.6,-96.5577], Tmin=(910.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-120.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC(C)C([O])=C=O(12016)',
    structure = SMILES('C=CC(C)C(=O)[C]=O'),
    E0 = (-78.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.851285,0.0744684,-8.24865e-05,5.41235e-08,-1.50987e-11,-9275.77,28.1204], Tmin=(100,'K'), Tmax=(853.599,'K')), NASAPolynomial(coeffs=[8.74031,0.0374998,-1.75226e-05,3.38597e-09,-2.38724e-13,-10622.6,-8.69078], Tmin=(853.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C[C]([CH2])C([O])=C[O](12017)',
    structure = SMILES('[CH2]C[C]([CH2])C([O])=C[O]'),
    E0 = (346.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,645.598,645.601,645.611,645.618],'cm^-1')),
        HinderedRotor(inertia=(0.0416004,'amu*angstrom^2'), symmetry=1, barrier=(12.3052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00888251,'amu*angstrom^2'), symmetry=1, barrier=(2.62716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0416054,'amu*angstrom^2'), symmetry=1, barrier=(12.3051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284114,'amu*angstrom^2'), symmetry=1, barrier=(84.0293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.335668,0.0675638,-3.23682e-05,-1.95351e-08,1.66567e-11,41850.2,34.7593], Tmin=(100,'K'), Tmax=(916.097,'K')), NASAPolynomial(coeffs=[20.2471,0.0151862,-3.19814e-06,4.20624e-10,-2.80125e-14,36751.7,-67.4726], Tmin=(916.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ(C)CO) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C=C)C([O])[CH][O](12018)',
    structure = SMILES('[CH2]C=C([CH2])C([O])[CH][O]'),
    E0 = (515.126,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,679.854,688.286],'cm^-1')),
        HinderedRotor(inertia=(0.0867984,'amu*angstrom^2'), symmetry=1, barrier=(29.1251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00120903,'amu*angstrom^2'), symmetry=1, barrier=(4.65964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.853069,'amu*angstrom^2'), symmetry=1, barrier=(19.6137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161155,'amu*angstrom^2'), symmetry=1, barrier=(53.149,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.462894,0.0814346,-8.79417e-05,5.17674e-08,-1.2538e-11,62079.7,32.1001], Tmin=(100,'K'), Tmax=(988.799,'K')), NASAPolynomial(coeffs=[12.2447,0.0337738,-1.56409e-05,3.02104e-09,-2.1345e-13,59749.7,-24.6076], Tmin=(988.799,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.126,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=C)C([O])[CH][O](12019)',
    structure = SMILES('[CH2]C([C]=C)C([O])[CH][O]'),
    E0 = (674.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.154701,0.0926866,-0.000140325,1.19816e-07,-4.02344e-11,81306.7,35.6595], Tmin=(100,'K'), Tmax=(860.459,'K')), NASAPolynomial(coeffs=[8.68938,0.0374306,-1.68382e-05,3.09663e-09,-2.08151e-13,80414.8,-0.881006], Tmin=(860.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(674.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(Isobutyl) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([C]=C)[C]([O])C[O](12020)',
    structure = SMILES('[CH2]C([C]=C)[C]([O])C[O]'),
    E0 = (671.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.231127,0.0931563,-0.000147934,1.32218e-07,-4.57971e-11,80860.5,35.2595], Tmin=(100,'K'), Tmax=(867.621,'K')), NASAPolynomial(coeffs=[6.89297,0.0403602,-1.8478e-05,3.41033e-09,-2.2901e-13,80535.7,8.8558], Tmin=(867.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C([O])[CH][O](12021)',
    structure = SMILES('[CH]=CC([CH2])C([O])[CH][O]'),
    E0 = (684.19,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3120,650,792.5,1650,266.325,266.432,266.443,2791.76],'cm^-1')),
        HinderedRotor(inertia=(0.0100162,'amu*angstrom^2'), symmetry=1, barrier=(55.4238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00237338,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0992,'amu*angstrom^2'), symmetry=1, barrier=(55.4342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250799,'amu*angstrom^2'), symmetry=1, barrier=(12.6493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.123068,0.0914845,-0.000131666,1.0709e-07,-3.47491e-11,82422.7,35.7036], Tmin=(100,'K'), Tmax=(846.805,'K')), NASAPolynomial(coeffs=[10.1114,0.0351147,-1.55379e-05,2.85184e-09,-1.92222e-13,81060.5,-8.87835], Tmin=(846.805,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(684.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(Isobutyl) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]CC([CH2])C([O])=[C][O](12022)',
    structure = SMILES('[CH2]CC([CH2])C([O])=[C][O]'),
    E0 = (433.935,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,375.706,375.733,375.739,1079.14],'cm^-1')),
        HinderedRotor(inertia=(0.00119437,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0011938,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00119487,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120654,'amu*angstrom^2'), symmetry=1, barrier=(12.0814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00305616,0.0792655,-8.22456e-05,4.42508e-08,-9.2721e-12,52341.9,38.2125], Tmin=(100,'K'), Tmax=(1214.71,'K')), NASAPolynomial(coeffs=[17.5661,0.0193472,-5.70616e-06,8.45221e-10,-5.0912e-14,48225.9,-49.3389], Tmin=(1214.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.935,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=CC([CH2])[C]([O])C[O](12023)',
    structure = SMILES('[CH]=CC([CH2])[C]([O])C[O]'),
    E0 = (680.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,180,1077.95,1091.52],'cm^-1')),
        HinderedRotor(inertia=(0.173781,'amu*angstrom^2'), symmetry=1, barrier=(3.99557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.174168,'amu*angstrom^2'), symmetry=1, barrier=(4.00447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169009,'amu*angstrom^2'), symmetry=1, barrier=(3.88584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171924,'amu*angstrom^2'), symmetry=1, barrier=(3.95288,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.192969,0.092037,-0.000139592,1.19945e-07,-4.05237e-11,81976.8,35.3265], Tmin=(100,'K'), Tmax=(859.758,'K')), NASAPolynomial(coeffs=[8.33467,0.0380093,-1.7157e-05,3.16053e-09,-2.12657e-13,81173.7,0.748754], Tmin=(859.758,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([CH]C)C([O])=[C][O](12024)',
    structure = SMILES('[CH2]C([CH]C)C([O])=[C][O]'),
    E0 = (423.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,377.526,380.357,382.589,2283.33],'cm^-1')),
        HinderedRotor(inertia=(0.0011703,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00115778,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00115327,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150137,'amu*angstrom^2'), symmetry=1, barrier=(15.6293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.327466,0.0752186,-7.58224e-05,4.05788e-08,-8.60737e-12,51039.9,37.7968], Tmin=(100,'K'), Tmax=(1151.31,'K')), NASAPolynomial(coeffs=[15.2969,0.02321,-8.06167e-06,1.34158e-09,-8.71902e-14,47593,-36.531], Tmin=(1151.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(423.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cs_S) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(C=C)[C]1OC1[O](12025)',
    structure = SMILES('[CH2]C(C=C)[C]1OC1[O]'),
    E0 = (355.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09963,0.0657735,-5.20835e-05,1.00174e-08,8.32405e-12,42883.8,30.3917], Tmin=(100,'K'), Tmax=(695.381,'K')), NASAPolynomial(coeffs=[9.85848,0.0308163,-9.9526e-06,1.52712e-09,-9.24435e-14,41292.7,-11.3644], Tmin=(695.381,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C1([O])[CH]O1(12026)',
    structure = SMILES('[CH2]C(C=C)C1([O])[CH]O1'),
    E0 = (353.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.20508,0.0844403,-8.29663e-05,4.08051e-08,-7.36607e-12,42749.7,35.4973], Tmin=(100,'K'), Tmax=(1656.49,'K')), NASAPolynomial(coeffs=[19.0186,0.0141348,1.39467e-07,-4.65952e-10,4.35165e-14,38995.3,-63.3859], Tmin=(1656.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CCsJO) + radical(Isobutyl)"""),
)

species(
    label = '[O]C=C([O])C1[CH]CC1(11282)',
    structure = SMILES('[O]C=C([O])C1[CH]CC1'),
    E0 = (129.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08513,0.0422197,4.12842e-05,-9.34125e-08,4.1485e-11,15745.7,29.5264], Tmin=(100,'K'), Tmax=(941.321,'K')), NASAPolynomial(coeffs=[20.5192,0.0149255,-3.32374e-06,5.7544e-10,-4.80034e-14,9637.43,-76.0673], Tmin=(941.321,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(C=C(C)OJ) + radical(C=COJ) + radical(cyclobutane)"""),
)

species(
    label = '[CH2][CH]C1CC([O])C1=O(11360)',
    structure = SMILES('[CH2][CH]C1CC([O])C1=O'),
    E0 = (309.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57861,0.0345669,4.48956e-05,-8.41619e-08,3.4772e-11,37358.3,32.9902], Tmin=(100,'K'), Tmax=(967.606,'K')), NASAPolynomial(coeffs=[15.74,0.0222613,-7.70471e-06,1.46308e-09,-1.10222e-13,32453.3,-46.0487], Tmin=(967.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ) + radical(CCJCC=O) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1[CH]COC1=C[O](12027)',
    structure = SMILES('[CH2]C1[CH]COC1=C[O]'),
    E0 = (156.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17312,0.0287775,0.000102878,-1.76272e-07,7.64403e-11,18913.3,28.8113], Tmin=(100,'K'), Tmax=(914.036,'K')), NASAPolynomial(coeffs=[27.423,0.000935477,5.74196e-06,-1.25196e-09,7.75437e-14,10479,-115.357], Tmin=(914.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[CH]COC=C1[O](12028)',
    structure = SMILES('[CH2]C1[CH]COC=C1[O]'),
    E0 = (142.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19268,0.0369951,6.42694e-05,-1.29291e-07,5.91988e-11,17214.2,25.416], Tmin=(100,'K'), Tmax=(894.064,'K')), NASAPolynomial(coeffs=[22.8853,0.00646072,3.89958e-06,-1.06142e-09,7.44067e-14,10676.7,-91.6769], Tmin=(894.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(C=C(C)OJ) + radical(CCJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC1C([O])=C[O](12029)',
    structure = SMILES('[CH2]C1CC1C([O])=C[O]'),
    E0 = (142.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.70062,0.0493285,3.28543e-05,-9.7047e-08,4.68174e-11,17279.2,29.3561], Tmin=(100,'K'), Tmax=(910.734,'K')), NASAPolynomial(coeffs=[24.3905,0.00662977,2.13704e-06,-5.97249e-10,3.80188e-14,10419.9,-96.6871], Tmin=(910.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1CC1([O])[CH][O](12030)',
    structure = SMILES('C=CC1CC1([O])[CH][O]'),
    E0 = (375.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.561836,0.0703478,-6.01599e-05,2.66527e-08,-4.76778e-12,45251.8,28.8258], Tmin=(100,'K'), Tmax=(1330.58,'K')), NASAPolynomial(coeffs=[15.0726,0.0267245,-1.09813e-05,2.01205e-09,-1.38011e-13,41390.3,-45.3245], Tmin=(1330.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1OC(=C[O])C1[CH2](11986)',
    structure = SMILES('[CH2]C1OC(=C[O])C1[CH2]'),
    E0 = (224.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.119191,0.0625568,9.84698e-06,-8.68034e-08,4.77833e-11,27204.6,25.5537], Tmin=(100,'K'), Tmax=(881.831,'K')), NASAPolynomial(coeffs=[28.2327,-0.000402983,7.11922e-06,-1.71501e-09,1.22683e-13,19736,-120.776], Tmin=(881.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ) + radical(Isobutyl) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1OC=C([O])C1[CH2](12031)',
    structure = SMILES('[CH2]C1OC=C([O])C1[CH2]'),
    E0 = (142.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.495483,0.0569175,1.29255e-05,-8.31532e-08,4.52784e-11,17240.9,24.5304], Tmin=(100,'K'), Tmax=(873.654,'K')), NASAPolynomial(coeffs=[24.9871,0.00293656,5.76195e-06,-1.49729e-09,1.10205e-13,10742.2,-103.021], Tmin=(873.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(CJC(C)OC)"""),
)

species(
    label = 'C=CC(=C)C([O])=C[O](12032)',
    structure = SMILES('C=CC(=C)C([O])=C[O]'),
    E0 = (21.3755,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.21693,'amu*angstrom^2'), symmetry=1, barrier=(27.9796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.2097,'amu*angstrom^2'), symmetry=1, barrier=(27.8133,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.512014,0.0624097,-2.32645e-05,-2.7372e-08,1.88839e-11,2709.67,25.5654], Tmin=(100,'K'), Tmax=(928.767,'K')), NASAPolynomial(coeffs=[20.4855,0.0130458,-2.74343e-06,3.94341e-10,-2.90052e-14,-2581.53,-77.8307], Tmin=(928.767,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.3755,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)(Cds-Cds)) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C=C)C([O])=C=O(12033)',
    structure = SMILES('[CH2]C(C=C)C(=O)[C]=O'),
    E0 = (132.483,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,269.074,269.613],'cm^-1')),
        HinderedRotor(inertia=(0.152349,'amu*angstrom^2'), symmetry=1, barrier=(7.81835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151801,'amu*angstrom^2'), symmetry=1, barrier=(7.82342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.456052,'amu*angstrom^2'), symmetry=1, barrier=(23.5369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152056,'amu*angstrom^2'), symmetry=1, barrier=(7.81787,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331767,0.087491,-0.000130018,1.09137e-07,-3.68368e-11,16059.6,30.3193], Tmin=(100,'K'), Tmax=(804.397,'K')), NASAPolynomial(coeffs=[9.51117,0.0341958,-1.63726e-05,3.12871e-09,-2.1633e-13,14830.3,-10.4298], Tmin=(804.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(CJC(C)C=O) + radical(CCCJ=O)"""),
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
    label = 'C=CC([O])=C[O](11103)',
    structure = SMILES('C=CC([O])=C[O]'),
    E0 = (-49.4969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20625,'amu*angstrom^2'), symmetry=1, barrier=(27.7341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36073,0.044686,-1.0018e-05,-3.58118e-08,2.21372e-11,-5845.51,19.4833], Tmin=(100,'K'), Tmax=(902.976,'K')), NASAPolynomial(coeffs=[19.5497,0.000521135,2.86666e-06,-6.71619e-10,4.54736e-14,-10614.7,-74.6305], Tmin=(902.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.4969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
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
    label = '[CH2][CH]C([O])=C[O](11447)',
    structure = SMILES('[CH2]C=C([O])[CH][O]'),
    E0 = (234.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,504.235,507.339,507.755],'cm^-1')),
        HinderedRotor(inertia=(0.000655989,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000659839,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91251,0.03846,-1.68356e-05,-7.51854e-09,6.05239e-12,28310.6,23.5168], Tmin=(100,'K'), Tmax=(1002.24,'K')), NASAPolynomial(coeffs=[11.8359,0.0149245,-5.66133e-06,1.0463e-09,-7.45272e-14,25514.4,-28.4065], Tmin=(1002.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(Allyl_P) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][C](C=C)C([O])=C[O](12034)',
    structure = SMILES('[CH2]C=C([CH2])C([O])=C[O]'),
    E0 = (146.442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,260.94,260.941,260.949],'cm^-1')),
        HinderedRotor(inertia=(0.744973,'amu*angstrom^2'), symmetry=1, barrier=(36.0021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.745138,'amu*angstrom^2'), symmetry=1, barrier=(36.0031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.745139,'amu*angstrom^2'), symmetry=1, barrier=(36.0028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.15519,0.0688071,-2.90918e-05,-2.99837e-08,2.21968e-11,17766,28.0554], Tmin=(100,'K'), Tmax=(910.909,'K')), NASAPolynomial(coeffs=[23.4272,0.00893261,-1.81253e-07,-1.42041e-10,9.54401e-15,11770.6,-91.6838], Tmin=(910.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C([C]=C)C([O])=C[O](12035)',
    structure = SMILES('[CH2]C([C]=C)C([O])=C[O]'),
    E0 = (354.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,389.761,390.732,392.946,393.242],'cm^-1')),
        HinderedRotor(inertia=(0.127038,'amu*angstrom^2'), symmetry=1, barrier=(13.7805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128662,'amu*angstrom^2'), symmetry=1, barrier=(13.8153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125961,'amu*angstrom^2'), symmetry=1, barrier=(13.7889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.08553,0.0815506,-8.78768e-05,4.33668e-08,-6.54574e-12,42724.1,31.5457], Tmin=(100,'K'), Tmax=(891.329,'K')), NASAPolynomial(coeffs=[17.9784,0.0167234,-4.81488e-06,7.13359e-10,-4.39263e-14,38919.8,-56.1659], Tmin=(891.329,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C([O])=C[O](12036)',
    structure = SMILES('[CH]=CC([CH2])C([O])=C[O]'),
    E0 = (363.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,3120,650,792.5,1650,333.817,334.357,335.67],'cm^-1')),
        HinderedRotor(inertia=(0.195473,'amu*angstrom^2'), symmetry=1, barrier=(15.5424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196153,'amu*angstrom^2'), symmetry=1, barrier=(15.5433,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195704,'amu*angstrom^2'), symmetry=1, barrier=(15.5328,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.354504,0.0852766,-9.6995e-05,5.45212e-08,-1.16387e-11,43857.7,33.0474], Tmin=(100,'K'), Tmax=(1257.04,'K')), NASAPolynomial(coeffs=[20.3511,0.0127612,-2.55412e-06,2.39894e-10,-8.92421e-15,39175.8,-69.4988], Tmin=(1257.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(C=C)C([O])=[C][O](12037)',
    structure = SMILES('[CH2]C(C=C)C([O])=[C][O]'),
    E0 = (355.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,458.983,459.299,460.182,460.286],'cm^-1')),
        HinderedRotor(inertia=(0.0624447,'amu*angstrom^2'), symmetry=1, barrier=(9.29159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619516,'amu*angstrom^2'), symmetry=1, barrier=(9.29595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0622596,'amu*angstrom^2'), symmetry=1, barrier=(9.29631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.265542,0.079494,-9.29299e-05,5.70372e-08,-1.37249e-11,42944.4,33.5514], Tmin=(100,'K'), Tmax=(1023.22,'K')), NASAPolynomial(coeffs=[15.2284,0.0210004,-7.17983e-06,1.16726e-09,-7.42402e-14,39882.4,-38.9789], Tmin=(1023.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = 'C2H2(1342)',
    structure = SMILES('C#C'),
    E0 = (218.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100,647.284,647.285,3592.16],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08745,0.00578572,8.56332e-06,-1.72824e-08,7.83594e-12,26273.8,4.4608], Tmin=(100,'K'), Tmax=(897.945,'K')), NASAPolynomial(coeffs=[5.89068,0.00209055,4.88462e-08,-5.66903e-11,4.15085e-15,25415.9,-10.7352], Tmin=(897.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.303,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(87.302,'J/(mol*K)'), label="""C2H2""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]C(=O)C[O](11458)',
    structure = SMILES('[CH2]C=C([O])C[O]'),
    E0 = (117.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,180,1986.31,3394.07],'cm^-1')),
        HinderedRotor(inertia=(0.613218,'amu*angstrom^2'), symmetry=1, barrier=(14.0991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.254474,'amu*angstrom^2'), symmetry=1, barrier=(60.136,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04046,0.0454536,-4.21713e-05,2.35123e-08,-5.69947e-12,14190.1,22.4681], Tmin=(100,'K'), Tmax=(963.158,'K')), NASAPolynomial(coeffs=[6.81765,0.0256136,-1.12726e-05,2.12492e-09,-1.48006e-13,13269.9,-0.399616], Tmin=(963.158,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'OCCO(T)(10910)',
    structure = SMILES('[O]C#C[O]'),
    E0 = (6.35398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2100,2250,500,550,344.171,344.172,344.187],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0202,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09814,0.022834,-4.17189e-05,3.85177e-08,-1.34153e-11,793.817,7.87311], Tmin=(100,'K'), Tmax=(871.13,'K')), NASAPolynomial(coeffs=[4.99347,0.00747035,-3.79496e-06,7.17517e-10,-4.83561e-14,716.336,0.441352], Tmin=(871.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.35398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(78.9875,'J/(mol*K)'), label="""OCCO(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'm1_allyl(186)',
    structure = SMILES('C=C[CH]C'),
    E0 = (121.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0800698,'amu*angstrom^2'), symmetry=1, barrier=(1.84096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.062381,'amu*angstrom^2'), symmetry=1, barrier=(19.3985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2897.21,'J/mol'), sigma=(5.21305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.54 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68514,0.0207746,2.19954e-05,-3.85566e-08,1.4972e-11,14712.8,14.0724], Tmin=(100,'K'), Tmax=(995.014,'K')), NASAPolynomial(coeffs=[7.60275,0.0207981,-7.87776e-06,1.45013e-09,-1.02662e-13,12754.4,-14.5511], Tmin=(995.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""m1_allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C[C](C)C([O])=C[O](12038)',
    structure = SMILES('[CH2]C=C(C)C([O])=C[O]'),
    E0 = (-5.05738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.148581,0.0714021,-3.81217e-05,-1.5762e-08,1.56839e-11,-457.239,27.7977], Tmin=(100,'K'), Tmax=(918.72,'K')), NASAPolynomial(coeffs=[21.1815,0.0149049,-3.14976e-06,4.19611e-10,-2.83237e-14,-5802.27,-79.9476], Tmin=(918.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.05738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C(C=C)C([O])=[C]O(12039)',
    structure = SMILES('[CH2]C(C=C)C([O])=[C]O'),
    E0 = (214.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,358.578,358.578,358.578],'cm^-1')),
        HinderedRotor(inertia=(0.129594,'amu*angstrom^2'), symmetry=1, barrier=(11.8244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129594,'amu*angstrom^2'), symmetry=1, barrier=(11.8244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129594,'amu*angstrom^2'), symmetry=1, barrier=(11.8244,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129594,'amu*angstrom^2'), symmetry=1, barrier=(11.8244,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0636832,0.0854014,-9.21387e-05,4.45919e-08,-5.91897e-12,25943.8,33.3289], Tmin=(100,'K'), Tmax=(862.335,'K')), NASAPolynomial(coeffs=[18.3601,0.0175347,-4.69029e-06,6.4522e-10,-3.74249e-14,22112.2,-56.6194], Tmin=(862.335,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(214.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C(C)C([O])=C[O](12040)',
    structure = SMILES('C=[C]C(C)C([O])=C[O]'),
    E0 = (148.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,265.654,265.655,265.656,265.656],'cm^-1')),
        HinderedRotor(inertia=(0.289796,'amu*angstrom^2'), symmetry=1, barrier=(14.5132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289799,'amu*angstrom^2'), symmetry=1, barrier=(14.5132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289796,'amu*angstrom^2'), symmetry=1, barrier=(14.5132,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0142076,0.0808844,-7.92229e-05,3.6327e-08,-5.50422e-12,18064.2,30.1184], Tmin=(100,'K'), Tmax=(988.108,'K')), NASAPolynomial(coeffs=[18.4103,0.0197095,-6.71337e-06,1.14044e-09,-7.67309e-14,13768.5,-61.8615], Tmin=(988.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C](C=C)C(O)=C[O](12041)',
    structure = SMILES('[CH2]C=C([CH2])C(O)=C[O]'),
    E0 = (8.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.208182,0.0727372,-2.13556e-05,-4.74971e-08,3.04864e-11,1208.84,27.859], Tmin=(100,'K'), Tmax=(910.718,'K')), NASAPolynomial(coeffs=[27.1615,0.00558035,1.87158e-06,-5.33086e-10,3.49217e-14,-5976.57,-113.704], Tmin=(910.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C(C=C)C(O)=[C][O](12042)',
    structure = SMILES('[CH2]C(C=C)[C](O)[C]=O'),
    E0 = (216.936,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,180,3309.33],'cm^-1')),
        HinderedRotor(inertia=(0.735419,'amu*angstrom^2'), symmetry=1, barrier=(16.9087,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.735644,'amu*angstrom^2'), symmetry=1, barrier=(16.9139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134637,'amu*angstrom^2'), symmetry=1, barrier=(82.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00217474,'amu*angstrom^2'), symmetry=1, barrier=(16.901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.73451,'amu*angstrom^2'), symmetry=1, barrier=(16.8878,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.320503,0.0765402,-7.96966e-05,4.36436e-08,-9.45858e-12,26227.6,36.3984], Tmin=(100,'K'), Tmax=(1127.11,'K')), NASAPolynomial(coeffs=[15.4321,0.02291,-8.32277e-06,1.42673e-09,-9.45163e-14,22821.1,-38.3146], Tmin=(1127.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([C]=C)C(O)=C[O](12043)',
    structure = SMILES('[CH2]C([C]=C)C(O)=C[O]'),
    E0 = (216.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,242.344,242.582,242.595],'cm^-1')),
        HinderedRotor(inertia=(0.430893,'amu*angstrom^2'), symmetry=1, barrier=(17.8973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429142,'amu*angstrom^2'), symmetry=1, barrier=(17.897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.42997,'amu*angstrom^2'), symmetry=1, barrier=(17.8979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.429201,'amu*angstrom^2'), symmetry=1, barrier=(17.8956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.816038,0.091916,-0.000103088,5.63164e-08,-1.15966e-11,26190.1,33.2742], Tmin=(100,'K'), Tmax=(1325.22,'K')), NASAPolynomial(coeffs=[22.781,0.0115002,-1.66333e-06,5.96861e-11,3.4126e-15,20743,-84.1672], Tmin=(1325.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC(C)C([O])=C[O](12044)',
    structure = SMILES('[CH]=CC(C)C([O])=C[O]'),
    E0 = (158.192,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.706668,'amu*angstrom^2'), symmetry=1, barrier=(16.2477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.707017,'amu*angstrom^2'), symmetry=1, barrier=(16.2557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706608,'amu*angstrom^2'), symmetry=1, barrier=(16.2463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.452636,0.0844158,-8.68151e-05,4.4244e-08,-8.64233e-12,19198,31.6259], Tmin=(100,'K'), Tmax=(1343.76,'K')), NASAPolynomial(coeffs=[21.3167,0.0149276,-4.01556e-06,5.69835e-10,-3.40495e-14,13770.6,-78.2562], Tmin=(1343.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.192,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CC(C)C([O])=[C][O](12045)',
    structure = SMILES('C=CC(C)C([O])=[C][O]'),
    E0 = (150.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,370.509,370.512,370.516,370.53],'cm^-1')),
        HinderedRotor(inertia=(0.107812,'amu*angstrom^2'), symmetry=1, barrier=(10.5024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107802,'amu*angstrom^2'), symmetry=1, barrier=(10.5025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107808,'amu*angstrom^2'), symmetry=1, barrier=(10.5024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.294286,0.07724,-7.83397e-05,4.1658e-08,-8.80108e-12,18279,31.6676], Tmin=(100,'K'), Tmax=(1151.65,'K')), NASAPolynomial(coeffs=[15.5837,0.0241353,-9.1715e-06,1.61758e-09,-1.09047e-13,14757.4,-44.2535], Tmin=(1151.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](C=C)C([O])=CO(12046)',
    structure = SMILES('[CH2]C=C([CH2])C([O])=CO'),
    E0 = (4.97919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.92578,0.095087,-9.85696e-05,4.79475e-08,-8.4894e-12,842.672,34.1406], Tmin=(100,'K'), Tmax=(1671.38,'K')), NASAPolynomial(coeffs=[25.2683,0.00638059,2.24324e-06,-7.20831e-10,5.52243e-14,-4947.84,-101.152], Tmin=(1671.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(4.97919,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH]=CC([CH2])C(O)=C[O](12047)',
    structure = SMILES('[CH]=CC([CH2])C(O)=C[O]'),
    E0 = (225.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05697,0.0931627,-0.00010293,5.4577e-08,-1.08005e-11,27315.3,34.0707], Tmin=(100,'K'), Tmax=(1410.88,'K')), NASAPolynomial(coeffs=[24.1445,0.00918416,-3.26384e-07,-1.99581e-10,2.08824e-14,21451.1,-91.7672], Tmin=(1410.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([C]=C)C([O])=CO(12048)',
    structure = SMILES('[CH2]C([C]=C)C([O])=CO'),
    E0 = (212.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,239.344,239.46,239.476],'cm^-1')),
        HinderedRotor(inertia=(0.369182,'amu*angstrom^2'), symmetry=1, barrier=(15.0287,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369289,'amu*angstrom^2'), symmetry=1, barrier=(15.0298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369024,'amu*angstrom^2'), symmetry=1, barrier=(15.0301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369384,'amu*angstrom^2'), symmetry=1, barrier=(15.0302,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.313445,0.0883218,-9.0296e-05,3.53493e-08,-7.43673e-13,25726.5,31.5706], Tmin=(100,'K'), Tmax=(866.618,'K')), NASAPolynomial(coeffs=[21.3529,0.0128362,-2.07847e-06,1.32346e-10,-2.17872e-15,21050.5,-75.1671], Tmin=(866.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C([O])=CO(12049)',
    structure = SMILES('[CH]=CC([CH2])C([O])=CO'),
    E0 = (221.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.720543,'amu*angstrom^2'), symmetry=1, barrier=(16.5667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.719911,'amu*angstrom^2'), symmetry=1, barrier=(16.5522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.722311,'amu*angstrom^2'), symmetry=1, barrier=(16.6073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.721506,'amu*angstrom^2'), symmetry=1, barrier=(16.5889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.999215,0.095057,-0.000110534,6.19402e-08,-1.29429e-11,26870.7,33.9472], Tmin=(100,'K'), Tmax=(1336.77,'K')), NASAPolynomial(coeffs=[23.4056,0.00929867,-1.68928e-08,-3.01548e-10,2.99886e-14,21483.6,-86.6203], Tmin=(1336.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(C=C)C(=O)C=O(11256)',
    structure = SMILES('[CH2]C(C=C)C(=O)C=O'),
    E0 = (-27.4778,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.876249,0.0776589,-8.29239e-05,3.63317e-08,3.97781e-12,-3200.88,27.6197], Tmin=(100,'K'), Tmax=(566.528,'K')), NASAPolynomial(coeffs=[8.10592,0.0401466,-1.94342e-05,3.78511e-09,-2.66874e-13,-4237.22,-5.06789], Tmin=(566.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.4778,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C=C)C([O])[C]=O(11251)',
    structure = SMILES('[CH2]C(C=C)C([O])[C]=O'),
    E0 = (284.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,180,180,2944.55],'cm^-1')),
        HinderedRotor(inertia=(0.033103,'amu*angstrom^2'), symmetry=1, barrier=(14.5134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0044906,'amu*angstrom^2'), symmetry=1, barrier=(1.96924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0330754,'amu*angstrom^2'), symmetry=1, barrier=(14.5148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.631373,'amu*angstrom^2'), symmetry=1, barrier=(14.5165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676589,0.0653481,-4.81298e-05,1.18476e-08,1.54306e-12,34288.8,36.476], Tmin=(100,'K'), Tmax=(978.579,'K')), NASAPolynomial(coeffs=[14.9937,0.0229897,-7.97767e-06,1.37265e-09,-9.30334e-14,30712.8,-36.2401], Tmin=(978.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1[CH]CO[C]1[CH][O](12050)',
    structure = SMILES('[CH2]C1[CH]CO[C]1[CH][O]'),
    E0 = (564.057,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.939371,0.0672499,-7.04912e-05,4.67461e-08,-1.25733e-11,67950.5,29.6994], Tmin=(100,'K'), Tmax=(1040.03,'K')), NASAPolynomial(coeffs=[7.92248,0.033069,-1.0631e-05,1.60474e-09,-9.48972e-14,66894.1,-2.36001], Tmin=(1040.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.057,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCOJ) + radical(C2CsJOCs) + radical(CCJCO) + radical(Isobutyl) + radical(CCsJOH)"""),
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
    label = '[CH2]C([C]=O)C=C(5422)',
    structure = SMILES('[CH2]C([C]=O)C=C'),
    E0 = (239.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,345.202],'cm^-1')),
        HinderedRotor(inertia=(0.115878,'amu*angstrom^2'), symmetry=1, barrier=(9.7983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115866,'amu*angstrom^2'), symmetry=1, barrier=(9.79811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115894,'amu*angstrom^2'), symmetry=1, barrier=(9.79834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69767,0.0540463,-5.7706e-05,3.67599e-08,-9.98018e-12,28849.4,23.3324], Tmin=(100,'K'), Tmax=(875.306,'K')), NASAPolynomial(coeffs=[7.4568,0.0277279,-1.26042e-05,2.40832e-09,-1.6882e-13,27841.2,-3.68502], Tmin=(875.306,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=CC1CC([O])C1=O(11264)',
    structure = SMILES('C=CC1CC([O])C1=O'),
    E0 = (31.2485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60019,0.0348504,4.24613e-05,-7.96061e-08,3.25519e-11,3860.79,28.4363], Tmin=(100,'K'), Tmax=(976.022,'K')), NASAPolynomial(coeffs=[14.9878,0.0240832,-8.76401e-06,1.67505e-09,-1.24941e-13,-852.976,-46.5863], Tmin=(976.022,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(31.2485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CC(=C)C(=O)C[O](12051)',
    structure = SMILES('C=CC(=C)C(=O)C[O]'),
    E0 = (33.2922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.32113,0.0735707,-6.64793e-05,3.03779e-08,-5.50844e-12,4142.66,28.1619], Tmin=(100,'K'), Tmax=(1329.82,'K')), NASAPolynomial(coeffs=[17.1818,0.0228546,-9.27256e-06,1.69876e-09,-1.16849e-13,-341.636,-57.9869], Tmin=(1329.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.2922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-Cds(Cds-O2d)(Cds-Cds)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C([C]=C)[C](O)[CH][O](12052)',
    structure = SMILES('[CH2]C([C]=C)[C](O)[CH][O]'),
    E0 = (621.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,360,370,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.182224,0.105659,-0.00018417,1.68249e-07,-5.79373e-11,74851.1,36.4844], Tmin=(100,'K'), Tmax=(889.586,'K')), NASAPolynomial(coeffs=[7.51742,0.0390578,-1.79444e-05,3.26557e-09,-2.15285e-13,74746.6,7.35118], Tmin=(889.586,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(C2CsJOH) + radical(Isobutyl) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])[C](O)[CH][O](12053)',
    structure = SMILES('[CH]=CC([CH2])[C](O)[CH][O]'),
    E0 = (630.457,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.226494,0.104616,-0.000176111,1.56362e-07,-5.28352e-11,75967.6,36.5731], Tmin=(100,'K'), Tmax=(887.789,'K')), NASAPolynomial(coeffs=[8.98232,0.0366661,-1.65993e-05,3.00997e-09,-1.98446e-13,75375.2,-0.885633], Tmin=(887.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(630.457,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(C2CsJOH) + radical(Isobutyl) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1CO[C]1[CH][O](12054)',
    structure = SMILES('C=CC1CO[C]1[CH][O]'),
    E0 = (366.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.777308,0.0733741,-8.65151e-05,6.30721e-08,-1.85596e-11,44179.8,28.7154], Tmin=(100,'K'), Tmax=(954.566,'K')), NASAPolynomial(coeffs=[7.94966,0.0348009,-1.25157e-05,2.0425e-09,-1.27693e-13,43198.6,-3.52069], Tmin=(954.566,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CCOJ) + radical(C2CsJOCs) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C(C=C)[C]1[CH]OO1(12055)',
    structure = SMILES('[CH2]C(C=C)[C]1[CH]OO1'),
    E0 = (575.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38999,0.0643739,-4.93505e-05,5.06804e-09,1.43434e-11,69302.6,26.1864], Tmin=(100,'K'), Tmax=(565.422,'K')), NASAPolynomial(coeffs=[6.13714,0.0418242,-1.8799e-05,3.55724e-09,-2.47629e-13,68589.4,4.43129], Tmin=(565.422,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C2CsJOO) + radical(CCsJOO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1[CH]CC([O])C1=O(12056)',
    structure = SMILES('[CH2]C1[CH]CC([O])C1=O'),
    E0 = (240.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37753,0.0383067,4.12571e-05,-8.45473e-08,3.59825e-11,28977.8,30.2976], Tmin=(100,'K'), Tmax=(957.923,'K')), NASAPolynomial(coeffs=[17.108,0.0207534,-6.62621e-06,1.23069e-09,-9.33725e-14,23755.7,-56.4432], Tmin=(957.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentanone) + radical(C=OCOJ) + radical(CCJCC=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C1C(=O)C([O])C1[CH2](12057)',
    structure = SMILES('[CH2]C1C(=O)C([O])C1[CH2]'),
    E0 = (311.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15087,0.047509,1.11827e-05,-5.27523e-08,2.49738e-11,37539.4,31.4024], Tmin=(100,'K'), Tmax=(950.928,'K')), NASAPolynomial(coeffs=[16.2506,0.0214199,-6.70108e-06,1.17452e-09,-8.48017e-14,32975.4,-49.5827], Tmin=(950.928,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ) + radical(Isobutyl) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2][CH]C(O)=C[O](11451)',
    structure = SMILES('[CH2]C=C(O)[CH][O]'),
    E0 = (96.9032,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,223.599,785.767],'cm^-1')),
        HinderedRotor(inertia=(0.0486319,'amu*angstrom^2'), symmetry=1, barrier=(21.3285,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0487672,'amu*angstrom^2'), symmetry=1, barrier=(21.3794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0712099,'amu*angstrom^2'), symmetry=1, barrier=(31.5916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56507,0.0422161,-8.57272e-06,-2.55601e-08,1.44823e-11,11752.7,23.2625], Tmin=(100,'K'), Tmax=(960.567,'K')), NASAPolynomial(coeffs=[15.4223,0.0118195,-3.74938e-06,6.88239e-10,-5.18665e-14,7830.7,-49.5909], Tmin=(960.567,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.9032,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)C(=O)C[O](12058)',
    structure = SMILES('[CH2]C(C=C)=C([O])C[O]'),
    E0 = (167.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,447.267,448.66,449.222,450.424],'cm^-1')),
        HinderedRotor(inertia=(0.0138969,'amu*angstrom^2'), symmetry=1, barrier=(2.00514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.357481,'amu*angstrom^2'), symmetry=1, barrier=(8.21918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.126264,'amu*angstrom^2'), symmetry=1, barrier=(17.8402,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.458752,0.0755005,-7.4335e-05,3.90153e-08,-8.25259e-12,20289.6,29.5685], Tmin=(100,'K'), Tmax=(1140.2,'K')), NASAPolynomial(coeffs=[14.1754,0.0273803,-1.10299e-05,2.00122e-09,-1.3687e-13,17161.6,-38.4059], Tmin=(1140.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=C)C(=O)C[O](12059)',
    structure = SMILES('[CH2]C([C]=C)C(=O)C[O]'),
    E0 = (361.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,347.491,347.523,347.526,3145.66],'cm^-1')),
        HinderedRotor(inertia=(0.0995989,'amu*angstrom^2'), symmetry=1, barrier=(8.53838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0996384,'amu*angstrom^2'), symmetry=1, barrier=(8.5385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0996791,'amu*angstrom^2'), symmetry=1, barrier=(8.53946,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0996284,'amu*angstrom^2'), symmetry=1, barrier=(8.53856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0964871,0.0955319,-0.000152377,1.34772e-07,-4.65453e-11,43665.5,32.7221], Tmin=(100,'K'), Tmax=(849.103,'K')), NASAPolynomial(coeffs=[8.29608,0.0382092,-1.80854e-05,3.40272e-09,-2.31588e-13,42947,-1.52635], Tmin=(849.103,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CJC(C)C=O) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C(=O)C[O](12060)',
    structure = SMILES('[CH]=CC([CH2])C(=O)C[O]'),
    E0 = (371.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3120,650,792.5,1650,277.853,277.97,3337],'cm^-1')),
        HinderedRotor(inertia=(0.304764,'amu*angstrom^2'), symmetry=1, barrier=(16.6718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143,'amu*angstrom^2'), symmetry=1, barrier=(7.82441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142856,'amu*angstrom^2'), symmetry=1, barrier=(7.82315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142606,'amu*angstrom^2'), symmetry=1, barrier=(7.82733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0655225,0.0943231,-0.000143703,1.22048e-07,-4.10752e-11,44781.5,32.7637], Tmin=(100,'K'), Tmax=(834.54,'K')), NASAPolynomial(coeffs=[9.70861,0.0359099,-1.6795e-05,3.1603e-09,-2.15857e-13,43596.6,-9.47051], Tmin=(834.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CJC(C)C=O) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=C([CH][O])OC=C(12061)',
    structure = SMILES('[CH2][CH]OC(C=C)=C[O]'),
    E0 = (199.351,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.92477,0.104467,-0.000117873,6.09635e-08,-1.15875e-11,24211.9,34.0454], Tmin=(100,'K'), Tmax=(1494.25,'K')), NASAPolynomial(coeffs=[29.6961,0.00215223,2.57161e-06,-6.86393e-10,5.09158e-14,16734.5,-124.606], Tmin=(1494.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(199.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]C([O])(C=C)C=O(12062)',
    structure = SMILES('[CH2][CH]C([O])(C=C)C=O'),
    E0 = (321.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180,180,180,1600,1840.33,2680.89,3200],'cm^-1')),
        HinderedRotor(inertia=(0.161511,'amu*angstrom^2'), symmetry=1, barrier=(3.71346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161511,'amu*angstrom^2'), symmetry=1, barrier=(3.71346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161511,'amu*angstrom^2'), symmetry=1, barrier=(3.71346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161511,'amu*angstrom^2'), symmetry=1, barrier=(3.71346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.48117,0.0832914,-0.000115053,9.37414e-08,-3.12349e-11,38777.7,33.6357], Tmin=(100,'K'), Tmax=(799.09,'K')), NASAPolynomial(coeffs=[8.75996,0.0363839,-1.67399e-05,3.15982e-09,-2.17608e-13,37629.1,-3.35587], Tmin=(799.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = 'HCO(1372)',
    structure = SMILES('[CH]=O'),
    E0 = (32.5964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1116.11,1837.47,2716.03],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06159,-0.0018195,9.37895e-06,-8.24742e-09,2.33851e-12,3918.59,3.98389], Tmin=(100,'K'), Tmax=(1105.12,'K')), NASAPolynomial(coeffs=[3.05335,0.00410464,-1.74962e-06,3.28532e-10,-2.28985e-14,4002.52,8.32038], Tmin=(1105.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.5964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]C([CH2])[C]=O(9176)',
    structure = SMILES('[CH2][CH]C([CH2])[C]=O'),
    E0 = (517.704,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.0825441,'amu*angstrom^2'), symmetry=1, barrier=(1.89785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0028141,'amu*angstrom^2'), symmetry=1, barrier=(1.89786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0825441,'amu*angstrom^2'), symmetry=1, barrier=(1.89785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0825456,'amu*angstrom^2'), symmetry=1, barrier=(1.89789,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65821,0.0540087,-5.62824e-05,3.36872e-08,-8.44873e-12,62347.6,27.948], Tmin=(100,'K'), Tmax=(950.253,'K')), NASAPolynomial(coeffs=[8.4254,0.0255231,-1.13173e-05,2.14139e-09,-1.49465e-13,61061.5,-4.35455], Tmin=(950.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.704,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CJC(C)C=O) + radical(RCCJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=CC1CC1([O])C=O(12063)',
    structure = SMILES('C=CC1CC1([O])C=O'),
    E0 = (69.7165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932348,0.0564248,-2.11934e-05,-1.48406e-08,1.02814e-11,8505.18,28.4383], Tmin=(100,'K'), Tmax=(992.258,'K')), NASAPolynomial(coeffs=[15.1182,0.0236883,-8.66647e-06,1.57585e-09,-1.11412e-13,4486.33,-45.9553], Tmin=(992.258,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.7165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'C=CC(=C)C([O])C=O(12064)',
    structure = SMILES('C=CC(=C)C([O])C=O'),
    E0 = (23.6025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.567731,0.0651787,-4.15601e-05,3.6499e-09,4.17413e-12,2971.31,30.786], Tmin=(100,'K'), Tmax=(1018.59,'K')), NASAPolynomial(coeffs=[16.6102,0.0224578,-8.50959e-06,1.56245e-09,-1.10341e-13,-1348.77,-52.0689], Tmin=(1018.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.6025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CC(C)=C([O])C=O(12065)',
    structure = SMILES('C=C[C](C)C(=O)C=O'),
    E0 = (-118.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848418,0.061237,-4.0486e-05,1.29618e-08,-1.66227e-12,-14106.7,26.7792], Tmin=(100,'K'), Tmax=(1798.25,'K')), NASAPolynomial(coeffs=[16.8939,0.0255456,-1.07143e-05,1.92456e-09,-1.27829e-13,-19877.4,-60.0467], Tmin=(1798.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-118.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(C=CCJ(C)C=O)"""),
)

species(
    label = '[CH2]C1[CH]CC1([O])C=O(12066)',
    structure = SMILES('[CH2]C1[CH]CC1([O])C=O'),
    E0 = (329.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.921614,0.0596145,-4.245e-05,1.60747e-08,-2.49558e-12,39731.1,32.122], Tmin=(100,'K'), Tmax=(1506.89,'K')), NASAPolynomial(coeffs=[13.2154,0.0269809,-9.96572e-06,1.70326e-09,-1.11303e-13,36026,-32.2298], Tmin=(1506.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(CC(C)(C=O)OJ) + radical(cyclobutane) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1CO[CH][C]1[O](12067)',
    structure = SMILES('C=CC1CO[CH][C]1[O]'),
    E0 = (288.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.882961,0.0560388,-4.92127e-06,-4.57994e-08,2.67167e-11,34842.2,24.464], Tmin=(100,'K'), Tmax=(875.979,'K')), NASAPolynomial(coeffs=[17.1939,0.0181383,-2.66121e-06,1.52566e-10,-3.09059e-15,30581.1,-60.0781], Tmin=(875.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOCs)"""),
)

species(
    label = '[CH2]C1C([CH2])C1([O])C=O(12068)',
    structure = SMILES('[CH2]C1C([CH2])C1([O])C=O'),
    E0 = (341.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.249621,0.0702877,-6.42627e-05,3.13062e-08,-5.95017e-12,41277,32.9723], Tmin=(100,'K'), Tmax=(1411.66,'K')), NASAPolynomial(coeffs=[15.8502,0.0205118,-5.45218e-06,7.36907e-10,-4.13498e-14,37427.6,-45.704], Tmin=(1411.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)(C=O)OJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC=C([O])C=O(12069)',
    structure = SMILES('[CH2]C=CC(=O)C=O'),
    E0 = (-67.0863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,227.519],'cm^-1')),
        HinderedRotor(inertia=(0.00325301,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01103,'amu*angstrom^2'), symmetry=1, barrier=(37.1205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01086,'amu*angstrom^2'), symmetry=1, barrier=(37.1205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81222,0.0495499,-3.62055e-05,1.33984e-08,-2.07814e-12,-7991.54,22.3076], Tmin=(100,'K'), Tmax=(1441.68,'K')), NASAPolynomial(coeffs=[10.0755,0.0266234,-1.23517e-05,2.36793e-09,-1.65368e-13,-10374.1,-20.5805], Tmin=(1441.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-67.0863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)C([O])C=O(12070)',
    structure = SMILES('[CH2]C=C([CH2])C([O])C=O'),
    E0 = (201.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1200.74,1200.74],'cm^-1')),
        HinderedRotor(inertia=(0.22531,'amu*angstrom^2'), symmetry=1, barrier=(5.18033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.872468,'amu*angstrom^2'), symmetry=1, barrier=(20.0598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0257713,'amu*angstrom^2'), symmetry=1, barrier=(26.3671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.872466,'amu*angstrom^2'), symmetry=1, barrier=(20.0597,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483313,0.0678006,-5.06166e-05,1.59654e-08,-1.04298e-12,24347.1,32.7872], Tmin=(100,'K'), Tmax=(1122.94,'K')), NASAPolynomial(coeffs=[16.1352,0.0246883,-9.91395e-06,1.82591e-09,-1.26955e-13,20034.8,-48.0879], Tmin=(1122.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([C]=C)C([O])C=O(12071)',
    structure = SMILES('[CH2]C([C]=C)C([O])C=O'),
    E0 = (361.922,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,180,657.068,3263.99],'cm^-1')),
        HinderedRotor(inertia=(0.802425,'amu*angstrom^2'), symmetry=1, barrier=(18.4493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0137219,'amu*angstrom^2'), symmetry=1, barrier=(18.4439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00243913,'amu*angstrom^2'), symmetry=1, barrier=(18.4375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.64237,'amu*angstrom^2'), symmetry=1, barrier=(83.7452,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.694604,0.0672031,-5.84184e-05,2.71385e-08,-5.10229e-12,43652.7,36.0682], Tmin=(100,'K'), Tmax=(1273.63,'K')), NASAPolynomial(coeffs=[13.8142,0.0259997,-9.89186e-06,1.73804e-09,-1.16473e-13,40310.8,-30.3994], Tmin=(1273.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CC([CH2])C([O])C=O(12072)',
    structure = SMILES('[CH]=CC([CH2])C([O])C=O'),
    E0 = (371.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,932.478,4000],'cm^-1')),
        HinderedRotor(inertia=(4.17957,'amu*angstrom^2'), symmetry=1, barrier=(96.0965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.697907,'amu*angstrom^2'), symmetry=1, barrier=(16.0463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135647,'amu*angstrom^2'), symmetry=1, barrier=(3.11879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0260202,'amu*angstrom^2'), symmetry=1, barrier=(16.0514,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621458,0.0665305,-5.18449e-05,1.75404e-08,-1.16832e-12,44770.6,36.2591], Tmin=(100,'K'), Tmax=(1039,'K')), NASAPolynomial(coeffs=[14.8162,0.0243547,-8.96151e-06,1.57779e-09,-1.07382e-13,41147.7,-36.005], Tmin=(1039,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Isobutyl) + radical(Cds_P)"""),
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
    E0 = (116.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (276.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (276.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (624.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (470.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (950.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (949.958,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (571.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (124.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (356.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (123.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (179.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (179.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (141.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (141.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (369.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (537.987,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (763.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (710.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (709.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (458.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (705.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (448.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (355.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (353.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (240.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (309.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (172.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (200.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (153.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (375.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (238.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (196.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (247.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (382.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (250.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (322.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (293.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (521.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (479.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (358.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (565.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (575.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (567.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (428.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (325.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (214.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (377.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (300.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (273.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (410.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (260.524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (202.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (195.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (229.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (227.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (245.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (339.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (116.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (447.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (620.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (725.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (124.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (179.579,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (660.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (655.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (366.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (575.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (240.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (311.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (376.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (297.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (463.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (404.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (513.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (484.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (605.69,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (121.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (139.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (144.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (329.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (288.662,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (341.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (354.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (359.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (503.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (415.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['OCHCO(3676)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC[CH]C([O])=C[O](12006)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C1[CH]CC([O])[C]1[O](12007)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.43608e+09,'s^-1'), n=0.0758676, Ea=(56.4793,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', 'C=C[CH]C([O])=C[O](12008)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/TwoDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', '[CH2]C([C]=C[O])C=C(9996)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(4)', '[CH]=C([O])C([CH2])C=C(5412)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C(C=C)C([O])=C[O](12009)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC1COC1=C[O](12010)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C(C=C)C1=COO1(12011)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(240.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 237.9 to 240.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC1COC=C1[O](12012)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSDS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC(=C)C(O)=C[O](12013)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C(C=C)C(O)=C=O(12014)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC(=C)C([O])=CO(12015)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC(C)C([O])=C=O(12016)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C[C]([CH2])C([O])=C[O](12017)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][C](C=C)C([O])[CH][O](12018)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([C]=C)C([O])[CH][O](12019)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([C]=C)[C]([O])C[O](12020)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=CC([CH2])C([O])[CH][O](12021)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]CC([CH2])C([O])=[C][O](12022)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=CC([CH2])[C]([O])C[O](12023)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH]C)C([O])=[C][O](12024)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C(C=C)[C]1OC1[O](12025)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(239.522,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 238.2 to 239.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C(C=C)C1([O])[CH]O1(12026)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(237.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 236.5 to 237.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[O]C=C([O])C1[CH]CC1(11282)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2][CH]C1CC([O])C1=O(11360)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(193.571,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 187.6 to 193.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C1[CH]COC1=C[O](12027)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C1[CH]COC=C1[O](12028)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.041e+08,'s^-1'), n=0.7, Ea=(84.7093,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_SDS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C1CC1C([O])=C[O](12029)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC1CC1([O])[CH][O](12030)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(259.001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 257.9 to 259.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C1OC(=C[O])C1[CH2](11986)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C1OC=C([O])C1[CH2](12031)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.249e+08,'s^-1'), n=0.846, Ea=(80.7428,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra_2H_pri;radadd_intra] for rate rule [R6_SMS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', 'C=CC(=C)C([O])=C[O](12032)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.997e+08,'cm^3/(mol*s)'), n=1.629, Ea=(14.7235,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 64 used for Cds-CdCd_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCd_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[CH2]C(C=C)C([O])=C=O(12033)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C2H3(60)', 'C=CC([O])=C[O](11103)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(9050,'cm^3/(mol*s)'), n=2.41, Ea=(13.5143,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 828 used for Cds-CdH_Cds-HH;CdsJ-H
Exact match found for rate rule [Cds-CdH_Cds-HH;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O][C]=C[O](9592)', 'butadiene13(1350)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['OCHCO(3676)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.04e+06,'m^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CO_rad] for rate rule [Cd_R;CO_rad/OneDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C2H3(60)', '[CH2][CH]C([O])=C[O](11447)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O][C]=C[O](9592)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.206e+07,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['H(3)', '[CH2][C](C=C)C([O])=C[O](12034)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(3)', '[CH2]C([C]=C)C([O])=C[O](12035)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(9.28426e+17,'m^3/(mol*s)'), n=-3.05017, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=3.08884453463, var=70.5675449198, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R
    Total Standard Deviation in ln(k): 24.6015908735
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(3)', '[CH]=CC([CH2])C([O])=C[O](12036)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['H(3)', '[CH2]C(C=C)C([O])=[C][O](12037)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C2H2(1342)', '[CH2][CH]C(=O)C[O](11458)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(4.94538e+12,'s^-1'), n=-0.334211, Ea=(312.18,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, correlation='Root_1R!H->C_2R!H->C_N-5R!H->O_Ext-5C-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_2R!H->C_N-5R!H->O_Ext-5C-R
    Total Standard Deviation in ln(k): 11.540182761524994
Exact match found for rate rule [Root_1R!H->C_2R!H->C_N-5R!H->O_Ext-5C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Retroene"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['OCCO(T)(10910)', 'm1_allyl(186)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(5.70488e+33,'s^-1'), n=-6.74695, Ea=(209.076,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.34291747530997396, var=16.416542173868212, Tref=1000.0, N=7, correlation='Root_1R!H->C_2R!H->C_N-5R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_2R!H->C_N-5R!H->O
    Total Standard Deviation in ln(k): 8.984253428860972
Exact match found for rate rule [Root_1R!H->C_2R!H->C_N-5R!H->O]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=C[C](C)C([O])=C[O](12038)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(58.4615,'s^-1'), n=3.15787, Ea=(98.4673,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_noH]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(C=C)C([O])=[C]O(12039)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=[C]C(C)C([O])=C[O](12040)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2][C](C=C)C(O)=C[O](12041)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.19923e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C(C=C)C(O)=[C][O](12042)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2]C([C]=C)C(O)=C[O](12043)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]=CC(C)C([O])=C[O](12044)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C=CC(C)C([O])=[C][O](12045)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2][C](C=C)C([O])=CO(12046)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.11e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH]=CC([CH2])C(O)=C[O](12047)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C([C]=C)C([O])=CO(12048)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSMS;Cd_rad_out;XH_out] for rate rule [R5H_SSMS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH]=CC([CH2])C([O])=CO(12049)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(173703,'s^-1'), n=1.89007, Ea=(118.15,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Cd_rad_out_singleH;XH_out] + [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 3.1622776601683795
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C(C=C)C(=O)C=O(11256)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH2]C1[CH]CO[C]1[CH][O](12050)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(6.43608e+09,'s^-1'), n=0.0758676, Ea=(56.4793,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH][O](1548)', '[CH2]C([C]=O)C=C(5422)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC1CC([O])C1=O(11264)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC(=C)C(=O)C[O](12051)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2]C([C]=C)[C](O)[CH][O](12052)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH]=CC([CH2])[C](O)[CH][O](12053)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC1CO[C]1[CH][O](12054)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(250.208,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 249.2 to 250.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction68',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C(C=C)[C]1[CH]OO1(12055)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(459.31,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra_O] for rate rule [R4_linear;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 458.2 to 459.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C1[CH]CC([O])C1=O(12056)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.45491e+07,'s^-1'), n=1.06599, Ea=(123.825,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 118.5 to 123.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction70',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C1C(=O)C([O])C1[CH2](12057)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(999564,'s^-1'), n=1.52333, Ea=(194.974,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 190.8 to 195.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C2H2(1342)', '[CH2][CH]C(O)=C[O](11451)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(5555.64,'s^-1'), n=2.32292, Ea=(260.718,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, correlation='Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_N-7CClOSSi->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_N-7CClOSSi->O
    Total Standard Deviation in ln(k): 11.540182761524994
Exact match found for rate rule [Root_1R!H->C_2R!H->C_Ext-1C-R_N-7R!H->N_N-7CClOSSi->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Retroene"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH2][C](C=C)C(=O)C[O](12058)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(11891.5,'s^-1'), n=2.58622, Ea=(129.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH2]C([C]=C)C(=O)C[O](12059)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(2.572e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction74',
    reactants = ['[CH]=CC([CH2])C(=O)C[O](12060)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[CH2]C=C([CH][O])OC=C(12061)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[CH2][CH]C([O])(C=C)C=O(12062)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction77',
    reactants = ['HCO(1372)', '[CH2][CH]C([CH2])[C]=O(9176)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction78',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC1CC1([O])C=O(12063)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction79',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC(=C)C([O])C=O(12064)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction80',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC(C)=C([O])C=O(12065)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction81',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C1[CH]CC1([O])C=O(12066)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(6.61918e+08,'s^-1'), n=0.930343, Ea=(213.187,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_pri_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 210.0 to 213.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction82',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=CC1CO[CH][C]1[O](12067)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(1.19156e+09,'s^-1'), n=0.640131, Ea=(172.483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 169.8 to 172.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction83',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['[CH2]C1C([CH2])C1([O])C=O(12068)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(9.65211e+08,'s^-1'), n=0.77, Ea=(225.809,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 223.8 to 225.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction84',
    reactants = ['CH2(T)(20)', 'C=CC=C([O])C=O(12069)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction85',
    reactants = ['[CH2][C](C=C)C([O])C=O(12070)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction86',
    reactants = ['[CH2]C([C]=C)C([O])C=O(12071)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction87',
    reactants = ['[CH]=CC([CH2])C([O])C=O(12072)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2931',
    isomers = [
        '[CH2]C(C=C)C([O])=C[O](11262)',
    ],
    reactants = [
        ('OCHCO(3676)', 'butadiene13(1350)'),
        ('OCHCO(3676)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2931',
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


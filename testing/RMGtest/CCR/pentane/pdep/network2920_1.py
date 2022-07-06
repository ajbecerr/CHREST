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
    label = 'C=C[CH]CC([O])[C]=O(11249)',
    structure = SMILES('[CH2]C=CCC([O])[C]=O'),
    E0 = (226.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950,416.263,416.376,416.385],'cm^-1')),
        HinderedRotor(inertia=(0.000972621,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0009721,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0986824,'amu*angstrom^2'), symmetry=1, barrier=(12.1432,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.098742,'amu*angstrom^2'), symmetry=1, barrier=(12.143,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.628726,0.0649128,-4.5628e-05,1.08261e-08,9.77948e-13,27333,34.9558], Tmin=(100,'K'), Tmax=(1058.95,'K')), NASAPolynomial(coeffs=[15.6144,0.0237596,-9.22297e-06,1.68714e-09,-1.17689e-13,23292.8,-42.2905], Tmin=(1058.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC[CH]C([O])[C]=O(12204)',
    structure = SMILES('C=CC[CH]C([O])[C]=O'),
    E0 = (286.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,443.436,443.448,443.448,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000857159,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000857297,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104103,'amu*angstrom^2'), symmetry=1, barrier=(14.5268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.104087,'amu*angstrom^2'), symmetry=1, barrier=(14.5268,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.763187,0.0621367,-4.15519e-05,8.10915e-09,1.67667e-12,34622.3,37.1273], Tmin=(100,'K'), Tmax=(1055.75,'K')), NASAPolynomial(coeffs=[15.0121,0.0237587,-9.2002e-06,1.68306e-09,-1.17447e-13,30743.8,-36.508], Tmin=(1055.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCJCO) + radical(CCCJ=O)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4441.08,'J/mol'), sigma=(7.16113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=693.69 K, Pc=27.44 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.449233,0.0836711,-8.67254e-05,4.49871e-08,-8.88727e-12,14145.4,32.9887], Tmin=(100,'K'), Tmax=(1373.96,'K')), NASAPolynomial(coeffs=[20.4882,0.0149466,-3.21414e-06,3.50409e-10,-1.6487e-14,9125.3,-72.0051], Tmin=(1373.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl)"""),
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
    label = 'C=C[CH]C([O])[C]=O(12205)',
    structure = SMILES('C=C[CH]C([O])[C]=O'),
    E0 = (227.627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,180,180,433.436],'cm^-1')),
        HinderedRotor(inertia=(0.044964,'amu*angstrom^2'), symmetry=1, barrier=(19.2791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.04495,'amu*angstrom^2'), symmetry=1, barrier=(19.2786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0449673,'amu*angstrom^2'), symmetry=1, barrier=(19.2788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10305,0.0552657,-3.96596e-05,7.59267e-09,2.07463e-12,27488.9,28.3649], Tmin=(100,'K'), Tmax=(1024.32,'K')), NASAPolynomial(coeffs=[15.1228,0.0165252,-6.36892e-06,1.18158e-09,-8.40215e-14,23777,-43.7084], Tmin=(1024.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(227.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CCJCO) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([CH][O])C=C(5588)',
    structure = SMILES('[CH2]C([CH][O])C=C'),
    E0 = (402.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,299.272,300.171,1811.51],'cm^-1')),
        HinderedRotor(inertia=(0.116385,'amu*angstrom^2'), symmetry=1, barrier=(7.46329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0278062,'amu*angstrom^2'), symmetry=1, barrier=(64.9565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117446,'amu*angstrom^2'), symmetry=1, barrier=(7.46693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49244,0.0585738,-6.87688e-05,5.25411e-08,-1.7037e-11,48554.8,25.9759], Tmin=(100,'K'), Tmax=(815.37,'K')), NASAPolynomial(coeffs=[6.2553,0.0319061,-1.36341e-05,2.49444e-09,-1.69196e-13,47887.9,4.64318], Tmin=(815.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Isobutyl) + radical(CCsJOH)"""),
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
    label = '[CH]C(C=C)C([O])[C]=O(12206)',
    structure = SMILES('[CH]C(C=C)C([O])[C]=O'),
    E0 = (527.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.568642,0.0671422,-5.65318e-05,2.06837e-08,-1.86774e-12,63535.1,36.087], Tmin=(100,'K'), Tmax=(1046.82,'K')), NASAPolynomial(coeffs=[16.3584,0.0198561,-7.47179e-06,1.34704e-09,-9.34554e-14,59514.4,-44.2267], Tmin=(1046.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(C=C)C1OC1=O(12207)',
    structure = SMILES('[CH2]C(C=C)C1OC1=O'),
    E0 = (17.0367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.99633,0.0575062,-2.95706e-05,2.96165e-10,3.30525e-12,2164.46,29.8018], Tmin=(100,'K'), Tmax=(1071.65,'K')), NASAPolynomial(coeffs=[12.4366,0.0297249,-1.15683e-05,2.08845e-09,-1.43541e-13,-1144.27,-30.1796], Tmin=(1071.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.0367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(2(co)oxirane) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC1COC1[C]=O(11254)',
    structure = SMILES('C=CC1COC1[C]=O'),
    E0 = (19.2573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46055,0.0367342,4.90532e-05,-1.00759e-07,4.56688e-11,2425.53,27.8204], Tmin=(100,'K'), Tmax=(897.98,'K')), NASAPolynomial(coeffs=[17.6093,0.0158489,-1.33128e-06,-4.69634e-11,5.62396e-15,-2532.9,-59.8101], Tmin=(897.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.2573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC1CC(=O)C1[O](11265)',
    structure = SMILES('C=CC1CC(=O)C1[O]'),
    E0 = (30.1719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70782,0.0299155,5.99437e-05,-9.9982e-08,4.03323e-11,3729.86,28.6104], Tmin=(100,'K'), Tmax=(965.792,'K')), NASAPolynomial(coeffs=[15.8821,0.0223206,-7.64141e-06,1.46575e-09,-1.11875e-13,-1391.7,-51.6195], Tmin=(965.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.1719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone) + radical(C=OCOJ)"""),
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
    label = 'C=CC(=C)C(O)[C]=O(12208)',
    structure = SMILES('C=CC(=C)C(O)[C]=O'),
    E0 = (-60.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.313254,0.070385,-5.04965e-05,6.92729e-09,4.65634e-12,-7094.48,31.7368], Tmin=(100,'K'), Tmax=(978.245,'K')), NASAPolynomial(coeffs=[18.7648,0.0182324,-6.24718e-06,1.11418e-09,-7.90123e-14,-11819.1,-62.5725], Tmin=(978.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
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
    label = '[CH2]C[C]([CH2])C([O])[C]=O(12209)',
    structure = SMILES('[CH2]C[C]([CH2])C([O])[C]=O'),
    E0 = (516.356,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,229.173,865.631,3094.78],'cm^-1')),
        HinderedRotor(inertia=(0.0226472,'amu*angstrom^2'), symmetry=1, barrier=(0.786123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0226472,'amu*angstrom^2'), symmetry=1, barrier=(0.786123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0226472,'amu*angstrom^2'), symmetry=1, barrier=(0.786123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0226472,'amu*angstrom^2'), symmetry=1, barrier=(0.786123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0226472,'amu*angstrom^2'), symmetry=1, barrier=(0.786123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.614009,0.0682183,-6.06389e-05,2.88055e-08,-5.50707e-12,62230.5,38.1753], Tmin=(100,'K'), Tmax=(1259.68,'K')), NASAPolynomial(coeffs=[14.217,0.0250242,-9.20541e-06,1.58569e-09,-1.05067e-13,58803.3,-30.5919], Tmin=(1259.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ(C)CO) + radical(RCCJ) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[O]C([C]=O)C1[CH]CC1(11276)',
    structure = SMILES('[O]C([C]=O)C1[CH]CC1'),
    E0 = (296.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45325,0.0418071,1.67598e-05,-4.99133e-08,2.13728e-11,35719.5,33.3144], Tmin=(100,'K'), Tmax=(996.235,'K')), NASAPolynomial(coeffs=[13.7963,0.025906,-9.97586e-06,1.89049e-09,-1.37355e-13,31590,-34.5701], Tmin=(996.235,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(C=OCOJ) + radical(cyclobutane) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1[CH]COC1[C]=O(12210)',
    structure = SMILES('[CH2]C1[CH]COC1[C]=O'),
    E0 = (216.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71536,0.0294272,6.9675e-05,-1.23811e-07,5.49124e-11,26192.3,28.4772], Tmin=(100,'K'), Tmax=(886.903,'K')), NASAPolynomial(coeffs=[17.5268,0.0142273,4.83199e-07,-4.67004e-10,3.68796e-14,21180.8,-58.3473], Tmin=(886.903,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1[CH]CC(=O)C1[O](12164)',
    structure = SMILES('[CH2]C1[CH]CC(=O)C1[O]'),
    E0 = (231.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71618,0.0291376,6.58162e-05,-1.09936e-07,4.54342e-11,27919.8,30.3709], Tmin=(100,'K'), Tmax=(940.817,'K')), NASAPolynomial(coeffs=[16.3374,0.0204099,-5.46578e-06,9.45564e-10,-7.21233e-14,22803.7,-51.8447], Tmin=(940.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentanone) + radical(C=OCOJ) + radical(CCJCC=O) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC1C([O])[C]=O(12211)',
    structure = SMILES('[CH2]C1CC1C([O])[C]=O'),
    E0 = (308.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07977,0.0488386,8.30347e-06,-5.29802e-08,2.61767e-11,37252.6,33.1014], Tmin=(100,'K'), Tmax=(935.801,'K')), NASAPolynomial(coeffs=[17.3761,0.0180941,-4.78919e-06,7.81706e-10,-5.65806e-14,32498.7,-53.541], Tmin=(935.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.744,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(C=OCOJ) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1OC([C]=O)C1[CH2](12200)',
    structure = SMILES('[CH2]C1OC([C]=O)C1[CH2]'),
    E0 = (293.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.699574,0.0747214,-6.808e-05,3.22119e-08,-5.6449e-12,35476.5,35.5774], Tmin=(100,'K'), Tmax=(1712.31,'K')), NASAPolynomial(coeffs=[15.6404,0.0176399,-1.50985e-06,-1.56199e-10,2.25818e-14,32653,-43.9462], Tmin=(1712.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(Isobutyl) + radical(CJC(C)OC) + radical(CCCJ=O)"""),
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
    label = 'C=CC(=C)C([O])[C]=O(12212)',
    structure = SMILES('C=CC(=C)C([O])[C]=O'),
    E0 = (183.563,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1855,455,950,380.747,380.75,380.752],'cm^-1')),
        HinderedRotor(inertia=(0.115512,'amu*angstrom^2'), symmetry=1, barrier=(11.883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115516,'amu*angstrom^2'), symmetry=1, barrier=(11.883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115513,'amu*angstrom^2'), symmetry=1, barrier=(11.8829,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.521179,0.0678507,-5.56685e-05,1.68506e-08,3.29213e-13,22210.5,31.7708], Tmin=(100,'K'), Tmax=(994.06,'K')), NASAPolynomial(coeffs=[17.1594,0.0180385,-6.36481e-06,1.1287e-09,-7.86588e-14,18055.8,-52.6588], Tmin=(994.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.563,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
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
    label = '[CH2]C(C=C)C=O(4917)',
    structure = SMILES('[CH2]C(C=C)C=O'),
    E0 = (80.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,279.063],'cm^-1')),
        HinderedRotor(inertia=(0.13781,'amu*angstrom^2'), symmetry=1, barrier=(7.49795,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141552,'amu*angstrom^2'), symmetry=1, barrier=(7.50677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235918,'amu*angstrom^2'), symmetry=1, barrier=(12.7234,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68421,0.054824,-5.37094e-05,3.28461e-08,-8.87771e-12,9763.62,22.2264], Tmin=(100,'K'), Tmax=(864.468,'K')), NASAPolynomial(coeffs=[6.53098,0.0323976,-1.47964e-05,2.83722e-09,-1.99413e-13,8925.63,-0.450611], Tmin=(864.468,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CJC(C)C=O)"""),
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
    label = 'C=CC([O])[C]=O(11100)',
    structure = SMILES('C=CC([O])[C]=O'),
    E0 = (133.356,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,351.422,4000],'cm^-1')),
        HinderedRotor(inertia=(0.171009,'amu*angstrom^2'), symmetry=1, barrier=(14.9744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171002,'amu*angstrom^2'), symmetry=1, barrier=(14.9751,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96577,0.0394789,-2.95633e-05,9.34197e-09,-5.31118e-13,16116.8,25.1576], Tmin=(100,'K'), Tmax=(1090.53,'K')), NASAPolynomial(coeffs=[10.8044,0.0145668,-5.62355e-06,1.01981e-09,-7.04881e-14,13742.6,-20.2963], Tmin=(1090.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.356,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C([O])[C]=O(11469)',
    structure = SMILES('[CH2][CH]C([O])[C]=O'),
    E0 = (411.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,307.307,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0164609,'amu*angstrom^2'), symmetry=1, barrier=(21.5658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164604,'amu*angstrom^2'), symmetry=1, barrier=(21.5664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00189942,'amu*angstrom^2'), symmetry=1, barrier=(21.5662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97361,0.0388712,-2.61177e-05,3.67988e-09,2.05984e-12,49612.9,29.6046], Tmin=(100,'K'), Tmax=(1004.17,'K')), NASAPolynomial(coeffs=[11.292,0.0131789,-4.80807e-06,8.64322e-10,-6.03862e-14,47165.4,-18.2588], Tmin=(1004.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C](C=C)C([O])[C]=O(12213)',
    structure = SMILES('[CH2]C=C([CH2])C([O])[C]=O'),
    E0 = (361.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,1129.89,1133.78],'cm^-1')),
        HinderedRotor(inertia=(0.246408,'amu*angstrom^2'), symmetry=1, barrier=(5.66541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0194009,'amu*angstrom^2'), symmetry=1, barrier=(17.7504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.770864,'amu*angstrom^2'), symmetry=1, barrier=(17.7237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910826,'amu*angstrom^2'), symmetry=1, barrier=(20.9417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.382081,0.0711031,-6.68517e-05,3.17939e-08,-5.94911e-12,43588.6,33.9689], Tmin=(100,'K'), Tmax=(1300.82,'K')), NASAPolynomial(coeffs=[17.1661,0.0194921,-7.33787e-06,1.29308e-09,-8.72244e-14,39222,-51.4184], Tmin=(1300.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([C]=C)C([O])[C]=O(12214)',
    structure = SMILES('[CH2]C([C]=C)C([O])[C]=O'),
    E0 = (521.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,377.111,377.245,2697.29],'cm^-1')),
        HinderedRotor(inertia=(0.00118557,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00118416,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138409,'amu*angstrom^2'), symmetry=1, barrier=(13.9523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13836,'amu*angstrom^2'), symmetry=1, barrier=(13.9534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678381,0.0695315,-7.14236e-05,3.91037e-08,-8.52205e-12,62890.6,36.9437], Tmin=(100,'K'), Tmax=(1117.99,'K')), NASAPolynomial(coeffs=[13.8446,0.0224239,-8.21824e-06,1.41303e-09,-9.36528e-14,59946.7,-28.0437], Tmin=(1117.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Isobutyl) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=CC([CH2])C([O])[C]=O(12215)',
    structure = SMILES('[CH]=CC([CH2])C([O])[C]=O'),
    E0 = (531.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,3120,650,792.5,1650,324.435,324.502],'cm^-1')),
        HinderedRotor(inertia=(0.00160081,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00160104,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00160189,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119547,'amu*angstrom^2'), symmetry=1, barrier=(8.93078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.452853,0.0706225,-7.081e-05,3.68485e-08,-7.52434e-12,64015,37.6828], Tmin=(100,'K'), Tmax=(1201.06,'K')), NASAPolynomial(coeffs=[15.9756,0.0189255,-6.24571e-06,1.01111e-09,-6.47944e-14,60286.3,-40.0496], Tmin=(1201.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Isobutyl) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = 'C=C[C](C)C([O])[C]=O(12216)',
    structure = SMILES('[CH2]C=C(C)C([O])[C]=O'),
    E0 = (209.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,334.79,334.898],'cm^-1')),
        HinderedRotor(inertia=(0.0015039,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00150332,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132449,'amu*angstrom^2'), symmetry=1, barrier=(10.5388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132459,'amu*angstrom^2'), symmetry=1, barrier=(10.5391,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.673991,0.070224,-6.39192e-05,3.07012e-08,-5.97881e-12,25352.2,32.637], Tmin=(100,'K'), Tmax=(1226.31,'K')), NASAPolynomial(coeffs=[13.871,0.0271773,-1.12648e-05,2.0761e-09,-1.43139e-13,22115.5,-33.7233], Tmin=(1226.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(209.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Allyl_P) + radical(CCCJ=O)"""),
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
    label = 'C=[C]C(C)C([O])[C]=O(12217)',
    structure = SMILES('C=[C]C(C)C([O])[C]=O'),
    E0 = (316.8,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,180,490.012,4000],'cm^-1')),
        HinderedRotor(inertia=(0.832719,'amu*angstrom^2'), symmetry=1, barrier=(19.1458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00168625,'amu*angstrom^2'), symmetry=1, barrier=(19.1458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.832718,'amu*angstrom^2'), symmetry=1, barrier=(19.1458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112372,'amu*angstrom^2'), symmetry=1, barrier=(19.1459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.62487,0.0682156,-5.99356e-05,2.74495e-08,-5.04673e-12,38228.8,35.3567], Tmin=(100,'K'), Tmax=(1304.94,'K')), NASAPolynomial(coeffs=[14.9217,0.0243916,-9.56089e-06,1.71411e-09,-1.16332e-13,34497.5,-37.4226], Tmin=(1304.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.8,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C](C=C)C(O)[C]=O(12218)',
    structure = SMILES('[CH2]C=C([CH2])C(O)[C]=O'),
    E0 = (117.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.252707,0.0727766,-5.9004e-05,1.89273e-08,-6.08328e-13,14280.2,33.6491], Tmin=(100,'K'), Tmax=(1035.26,'K')), NASAPolynomial(coeffs=[17.8908,0.0211026,-8.00485e-06,1.45845e-09,-1.02165e-13,9745.28,-56.3201], Tmin=(1035.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([C]=C)C(O)[C]=O(12219)',
    structure = SMILES('[CH2]C([C]=C)C(O)[C]=O'),
    E0 = (278.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314828,0.0738857,-7.25225e-05,3.71429e-08,-7.50622e-12,33592.4,37.4687], Tmin=(100,'K'), Tmax=(1208.79,'K')), NASAPolynomial(coeffs=[16.1675,0.0214276,-7.42654e-06,1.24139e-09,-8.11058e-14,29759.9,-42.0172], Tmin=(1208.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    label = '[CH]=CC(C)C([O])[C]=O(12220)',
    structure = SMILES('[CH]=CC(C)C([O])[C]=O'),
    E0 = (326.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1855,455,950,3120,650,792.5,1650,314.625,314.625],'cm^-1')),
        HinderedRotor(inertia=(0.00170325,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.001703,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139798,'amu*angstrom^2'), symmetry=1, barrier=(9.81884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139786,'amu*angstrom^2'), symmetry=1, barrier=(9.81888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.527275,0.0678318,-5.43767e-05,1.91737e-08,-1.67838e-12,39347.7,35.6353], Tmin=(100,'K'), Tmax=(1069.23,'K')), NASAPolynomial(coeffs=[15.9042,0.0227673,-8.63762e-06,1.55474e-09,-1.07267e-13,35347.1,-42.9098], Tmin=(1069.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CC([CH2])C(O)[C]=O(12221)',
    structure = SMILES('[CH]=CC([CH2])C(O)[C]=O'),
    E0 = (287.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0523161,0.0753627,-7.30512e-05,3.61214e-08,-6.94562e-12,34718.6,38.3442], Tmin=(100,'K'), Tmax=(1330.45,'K')), NASAPolynomial(coeffs=[18.3585,0.0178493,-5.41707e-06,8.32326e-10,-5.17535e-14,30066.6,-54.3758], Tmin=(1330.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCCJ=O) + radical(Cds_P)"""),
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
    E0 = (284.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (443.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (532.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (447.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (665.429,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (897.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (897.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (738.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (286.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (292.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (292.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (306.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (306.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (347.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (347.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (347.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (538.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (497.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (448.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (408.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (340.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (340.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (321.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (406.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (371.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (300.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (395.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (365.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (575.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (450.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (322.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (434.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (284.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (479.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (698.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (573.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (567.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (733.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (742.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (413.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (398.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (401.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (468.679,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (359.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (400.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (429.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (463.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (370.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (395.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (404.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['OCHCO(3676)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['C=C[CH]CC([O])[C]=O(11249)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=CC[CH]C([O])[C]=O(12204)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', 'C=C[CH]C([O])[C]=O(12205)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', '[CH2]C([CH][C]=O)C=C(5308)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]=O(2355)', '[CH2]C([CH][O])C=C(5588)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C(C=C)C([O])[C]=O(12206)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C(C=C)C1OC1=O(12207)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['C=CC1COC1[C]=O(11254)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['C=CC1CC(=O)C1[O](11265)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C(C=C)C(O)=C=O(12014)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C(C=C)C(=O)C=O(11256)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['C=CC(=C)C(O)[C]=O(12208)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['C=CC(=C)C([O])C=O(12064)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['C=CC(C)C([O])=C=O(12016)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C[C]([CH2])C([O])[C]=O(12209)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC([CH2])C([O])=[C][O](12022)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH]C)C([O])=[C][O](12024)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[O]C([C]=O)C1[CH]CC1(11276)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C1[CH]COC1[C]=O(12210)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C1[CH]CC(=O)C1[O](12164)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.07178e+07,'s^-1'), n=1.01592, Ea=(56.0365,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri_2H;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri_2H;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C1CC1C([O])[C]=O(12211)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.32e+08,'s^-1'), n=0.97, Ea=(37.2376,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 335 used for R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H
Exact match found for rate rule [R4_S_D;doublebond_intra_2H_pri;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C1OC([C]=O)C1[CH2](12200)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 33 used for R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O
Exact match found for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C1C(=O)C([O])C1[CH2](12057)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.65009e+08,'s^-1'), n=1.00067, Ea=(87.7294,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_2H_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_2H_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CO(2039)', '[CH2]C([CH][O])C=C(5588)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1099.39,'m^3/(mol*s)'), n=0.635039, Ea=(16.7529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;Y_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', 'C=CC(=C)C([O])[C]=O(12212)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2]C(C=C)C([O])=C=O(12033)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[C]=O(2355)', '[CH2]C(C=C)C=O(4917)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CO_birad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O(4)', '[CH2]C(C=C)C=C=O(5306)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(11.358,'m^3/(mol*s)'), n=1.81033, Ea=(28.0982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;YJ] for rate rule [Cds-CsH_Ck;O_atom_triplet]
Euclidian distance = 2.23606797749979
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O][C]=C[O](9592)', 'butadiene13(1350)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.0136943,'m^3/(mol*s)'), n=2.49, Ea=(21.946,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CdH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C2H3(60)', 'C=CC([O])[C]=O(11100)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6870,'cm^3/(mol*s)'), n=2.41, Ea=(13.7235,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 432 used for Cds-CsH_Cds-HH;CdsJ-H
Exact match found for rate rule [Cds-CsH_Cds-HH;CdsJ-H]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['OCHCO(3676)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4679.9,'m^3/(mol*s)'), n=0.573452, Ea=(84.8735,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 79.3 to 84.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O][C]=C[O](9592)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.206e+07,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C2H3(60)', '[CH2][CH]C([O])[C]=O(11469)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH2][C](C=C)C([O])[C]=O(12213)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)', '[CH2]C(C=C)C([O])=[C][O](12037)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(3)', '[CH2]C([C]=C)C([O])[C]=O(12214)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(3)', '[CH]=CC([CH2])C([O])[C]=O(12215)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C[C](C)C([O])[C]=O(12216)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.307e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C(C=C)C(O)=[C][O](12042)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['C=CC(C)C([O])=[C][O](12045)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=[C]C(C)C([O])[C]=O(12217)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 204 used for R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H
Exact match found for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2][C](C=C)C(O)[C]=O(12218)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2][C](C=C)C([O])C=O(12070)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH2]C([C]=C)C(O)[C]=O(12219)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(274,'s^-1'), n=3.09, Ea=(145.603,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;O_rad_out;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C([C]=C)C([O])C=O(12071)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=CC(C)C([O])[C]=O(12220)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    products = ['[CH]=CC([CH2])C(O)[C]=O(12221)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.468e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=CC([CH2])C([O])C=O(12072)'],
    products = ['[CH2]C(C=C)C([O])[C]=O(11251)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSS;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 3.1622776601683795
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2920',
    isomers = [
        '[CH2]C(C=C)C([O])[C]=O(11251)',
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
    label = 'PDepNetwork #2920',
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


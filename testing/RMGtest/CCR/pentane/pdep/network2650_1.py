species(
    label = '[CH2]C(CC)C([O])[C]=O(10995)',
    structure = SMILES('[CH2]C(CC)C([O])[C]=O'),
    E0 = (158.536,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.234012,0.0749462,-6.42413e-05,2.95797e-08,-5.49434e-12,19209.7,37.5589], Tmin=(100,'K'), Tmax=(1293.72,'K')), NASAPolynomial(coeffs=[15.3054,0.0283478,-1.02131e-05,1.73854e-09,-1.1429e-13,15310,-39.0329], Tmin=(1293.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.536,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2]C(C)C([O])[C]=O(10955)',
    structure = SMILES('[CH2]C(C)C([O])[C]=O'),
    E0 = (178.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,766.553,2795.38],'cm^-1')),
        HinderedRotor(inertia=(0.518413,'amu*angstrom^2'), symmetry=1, barrier=(11.9193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.027975,'amu*angstrom^2'), symmetry=1, barrier=(11.9512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125456,'amu*angstrom^2'), symmetry=1, barrier=(2.88448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0276913,'amu*angstrom^2'), symmetry=1, barrier=(11.929,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4030.69,'J/mol'), sigma=(6.74566,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=629.58 K, Pc=29.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0883,0.0577557,-4.22145e-05,1.13117e-08,8.01732e-13,21635.4,32.2451], Tmin=(100,'K'), Tmax=(978.621,'K')), NASAPolynomial(coeffs=[12.8949,0.0223973,-7.79023e-06,1.33013e-09,-8.92289e-14,18706.9,-27.6155], Tmin=(978.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC[CH]C([O])[C]=O(11769)',
    structure = SMILES('CCC[CH]C([O])[C]=O'),
    E0 = (159.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.307579,0.0715,-5.58375e-05,2.2561e-08,-3.66961e-12,19270.3,38.2989], Tmin=(100,'K'), Tmax=(1461.31,'K')), NASAPolynomial(coeffs=[16.4877,0.0272105,-1.03752e-05,1.82046e-09,-1.21317e-13,14541.5,-45.8985], Tmin=(1461.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[CH]CC([O])[C]=O(10993)',
    structure = SMILES('CC[CH]CC([O])[C]=O'),
    E0 = (153.607,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,180,400,1637.38,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0502756,'amu*angstrom^2'), symmetry=1, barrier=(4.62374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0502756,'amu*angstrom^2'), symmetry=1, barrier=(4.62374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0502756,'amu*angstrom^2'), symmetry=1, barrier=(4.62374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0502756,'amu*angstrom^2'), symmetry=1, barrier=(4.62374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0502756,'amu*angstrom^2'), symmetry=1, barrier=(4.62374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03805,0.0692662,-6.07486e-05,3.25479e-08,-7.74687e-12,18577.8,35.7471], Tmin=(100,'K'), Tmax=(967.192,'K')), NASAPolynomial(coeffs=[7.6014,0.0421222,-1.86516e-05,3.53127e-09,-2.46646e-13,17308.2,4.30162], Tmin=(967.192,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(RCCJCC) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(CC)C([O])=C[O](11005)',
    structure = SMILES('[CH2]C(CC)C([O])=C[O]'),
    E0 = (-11.0556,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,434.928,434.937,434.941,434.943],'cm^-1')),
        HinderedRotor(inertia=(0.000891128,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108148,'amu*angstrom^2'), symmetry=1, barrier=(14.5182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108149,'amu*angstrom^2'), symmetry=1, barrier=(14.5183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108149,'amu*angstrom^2'), symmetry=1, barrier=(14.5182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4443.28,'J/mol'), sigma=(7.37584,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=694.03 K, Pc=25.13 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0224168,0.0735287,-3.34221e-05,-2.18835e-08,1.79193e-11,-1173.53,33.2088], Tmin=(100,'K'), Tmax=(918.819,'K')), NASAPolynomial(coeffs=[20.9481,0.0191389,-4.55672e-06,6.5433e-10,-4.38179e-14,-6568.44,-74.4062], Tmin=(918.819,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.0556,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl)"""),
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
    label = 'CC[CH]C([O])[C]=O(11593)',
    structure = SMILES('CC[CH]C([O])[C]=O'),
    E0 = (182.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,180,374.673,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0156426,'amu*angstrom^2'), symmetry=1, barrier=(17.7199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.015657,'amu*angstrom^2'), symmetry=1, barrier=(17.7262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00156177,'amu*angstrom^2'), symmetry=1, barrier=(17.7324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771392,'amu*angstrom^2'), symmetry=1, barrier=(17.7358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1765,0.0542492,-3.41105e-05,5.37361e-09,1.88839e-12,22097.8,32.9249], Tmin=(100,'K'), Tmax=(1056.03,'K')), NASAPolynomial(coeffs=[13.093,0.0227501,-8.74003e-06,1.58622e-09,-1.09981e-13,18820.5,-28.8154], Tmin=(1056.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(CCCJ=O)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79021,0.0664689,-5.56502e-05,2.64471e-08,-5.20762e-12,30453.4,29.7261], Tmin=(100,'K'), Tmax=(1206.53,'K')), NASAPolynomial(coeffs=[11.5566,0.0307744,-1.12727e-05,1.92592e-09,-1.26606e-13,27855.4,-24.2368], Tmin=(1206.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([CH][O])CC(1756)',
    structure = SMILES('[CH2]C([CH][O])CC'),
    E0 = (277.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,180,1187.06,1194.76],'cm^-1')),
        HinderedRotor(inertia=(0.181577,'amu*angstrom^2'), symmetry=1, barrier=(4.17481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183737,'amu*angstrom^2'), symmetry=1, barrier=(4.22447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.181979,'amu*angstrom^2'), symmetry=1, barrier=(4.18404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183609,'amu*angstrom^2'), symmetry=1, barrier=(4.22153,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19047,0.0665903,-7.97859e-05,6.4446e-08,-2.19921e-11,33469.6,26.5501], Tmin=(100,'K'), Tmax=(826.775,'K')), NASAPolynomial(coeffs=[5.3147,0.0393165,-1.70223e-05,3.12742e-09,-2.12286e-13,33037.8,8.95072], Tmin=(826.775,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(CCsJOH)"""),
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
    label = '[CH]C(CC)C([O])[C]=O(11770)',
    structure = SMILES('[CH]C(CC)C([O])[C]=O'),
    E0 = (401.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273049,0.0750207,-6.66915e-05,3.07841e-08,-5.68725e-12,48449.6,36.6416], Tmin=(100,'K'), Tmax=(1301.83,'K')), NASAPolynomial(coeffs=[16.2562,0.0259107,-1.01056e-05,1.8064e-09,-1.22423e-13,44288.1,-44.6837], Tmin=(1301.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(CC)C1OC1=O(11771)',
    structure = SMILES('[CH2]C(CC)C1OC1=O'),
    E0 = (-108.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.556184,0.0669101,-4.42597e-05,1.52001e-08,-2.14459e-12,-12914.5,30.8873], Tmin=(100,'K'), Tmax=(1630.2,'K')), NASAPolynomial(coeffs=[14.7233,0.0321485,-1.22745e-05,2.11988e-09,-1.38664e-13,-17533.5,-44.3844], Tmin=(1630.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(2(co)oxirane) + radical(Isobutyl)"""),
)

species(
    label = 'CCC1COC1[C]=O(10998)',
    structure = SMILES('CCC1COC1[C]=O'),
    E0 = (-109.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15431,0.0447972,3.79002e-05,-8.87348e-08,4.06979e-11,-13062.1,28.4102], Tmin=(100,'K'), Tmax=(903.394,'K')), NASAPolynomial(coeffs=[16.7,0.0232058,-4.68851e-06,5.78681e-10,-3.68565e-14,-17798.6,-55.6786], Tmin=(903.394,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-109.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC1CC(=O)C1[O](11008)',
    structure = SMILES('CCC1CC(=O)C1[O]'),
    E0 = (-98.6799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39138,0.0380793,4.85513e-05,-8.78597e-08,3.54374e-11,-11757.3,29.2378], Tmin=(100,'K'), Tmax=(979.116,'K')), NASAPolynomial(coeffs=[15.1124,0.0294439,-1.08654e-05,2.06019e-09,-1.51784e-13,-16717.2,-48.276], Tmin=(979.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-98.6799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C(CC)C(O)=C=O(11646)',
    structure = SMILES('[CH2]C(CC)C(O)=C=O'),
    E0 = (-108.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.674857,0.0978761,-0.000110376,6.37887e-08,-1.44049e-11,-12838,30.8694], Tmin=(100,'K'), Tmax=(1089.21,'K')), NASAPolynomial(coeffs=[19.4644,0.0239158,-8.52031e-06,1.44519e-09,-9.52636e-14,-17225.1,-68.0113], Tmin=(1089.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O)"""),
)

species(
    label = '[CH2]C(CC)C(=O)C=O(10999)',
    structure = SMILES('[CH2]C(CC)C(=O)C=O'),
    E0 = (-154.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.413257,0.0853421,-9.74583e-05,6.8176e-08,-2.04895e-11,-18413.6,29.4206], Tmin=(100,'K'), Tmax=(792.757,'K')), NASAPolynomial(coeffs=[8.39136,0.0450875,-2.12922e-05,4.12506e-09,-2.90877e-13,-19678.6,-7.2164], Tmin=(792.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJC(C)C=O)"""),
)

species(
    label = 'C=C(CC)C(O)[C]=O(11772)',
    structure = SMILES('C=C(CC)C(O)[C]=O'),
    E0 = (-172.077,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.107204,0.0756088,-6.26722e-05,2.66797e-08,-4.53618e-12,-20547.7,34.726], Tmin=(100,'K'), Tmax=(1410.19,'K')), NASAPolynomial(coeffs=[17.5806,0.0260458,-9.95309e-06,1.75695e-09,-1.17869e-13,-25475.9,-55.5793], Tmin=(1410.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-172.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C(CC)C([O])C=O(11686)',
    structure = SMILES('C=C(CC)C([O])C=O'),
    E0 = (-88.3049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.430239,0.0695239,-5.03361e-05,1.85404e-08,-2.76256e-12,-10484.8,33.5339], Tmin=(100,'K'), Tmax=(1570.49,'K')), NASAPolynomial(coeffs=[16.4021,0.0288434,-1.14809e-05,2.04623e-09,-1.36882e-13,-15501.4,-50.7304], Tmin=(1570.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-88.3049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'CCC(C)C([O])=C=O(11648)',
    structure = SMILES('CCC(C)C(=O)[C]=O'),
    E0 = (-204.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.686341,0.0776377,-7.46841e-05,4.22101e-08,-1.02905e-11,-24501.4,28.9048], Tmin=(100,'K'), Tmax=(960.195,'K')), NASAPolynomial(coeffs=[9.21424,0.0421127,-1.91886e-05,3.68007e-09,-2.58898e-13,-26139.2,-11.8911], Tmin=(960.195,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
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
    label = 'C=C(CC)C([O])[C]=O(11773)',
    structure = SMILES('C=C(CC)C([O])[C]=O'),
    E0 = (71.6557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,348.168,348.205,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00139076,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156222,'amu*angstrom^2'), symmetry=1, barrier=(13.4382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156246,'amu*angstrom^2'), symmetry=1, barrier=(13.4384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156236,'amu*angstrom^2'), symmetry=1, barrier=(13.4381,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.548842,0.0704377,-5.9158e-05,2.60539e-08,-4.65065e-12,8746.89,33.9133], Tmin=(100,'K'), Tmax=(1331.31,'K')), NASAPolynomial(coeffs=[14.7267,0.0278401,-1.11635e-05,2.02052e-09,-1.37617e-13,4971.81,-38.5441], Tmin=(1331.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(71.6557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(CC)C([O])=C=O(11659)',
    structure = SMILES('[CH2]C(CC)C(=O)[C]=O'),
    E0 = (5.83529,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,1855,455,950,325.128,325.182],'cm^-1')),
        HinderedRotor(inertia=(0.09186,'amu*angstrom^2'), symmetry=1, barrier=(6.89198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0918678,'amu*angstrom^2'), symmetry=1, barrier=(6.89194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0918634,'amu*angstrom^2'), symmetry=1, barrier=(6.89149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0918808,'amu*angstrom^2'), symmetry=1, barrier=(6.89216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.252155,'amu*angstrom^2'), symmetry=1, barrier=(18.9135,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.481039,0.085921,-9.94919e-05,5.53936e-08,-6.17982e-12,820.554,30.035], Tmin=(100,'K'), Tmax=(603.684,'K')), NASAPolynomial(coeffs=[9.52124,0.0396909,-1.85887e-05,3.56025e-09,-2.48385e-13,-520.024,-11.0792], Tmin=(603.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.83529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJC(C)C=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(C=O)CC(1767)',
    structure = SMILES('[CH2]C(C=O)CC'),
    E0 = (-46.1335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,275.596],'cm^-1')),
        HinderedRotor(inertia=(0.00221825,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155542,'amu*angstrom^2'), symmetry=1, barrier=(8.39814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155877,'amu*angstrom^2'), symmetry=1, barrier=(8.40274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155826,'amu*angstrom^2'), symmetry=1, barrier=(8.39976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.1244,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51473,0.0581097,-4.65826e-05,2.21351e-08,-4.70371e-12,-5461.99,23.0227], Tmin=(100,'K'), Tmax=(1061.47,'K')), NASAPolynomial(coeffs=[7.21266,0.0366374,-1.62387e-05,3.07695e-09,-2.1498e-13,-6671.6,-4.80635], Tmin=(1061.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.1335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O)"""),
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
    label = '[CH2][C](CC)C([O])[C]=O(11774)',
    structure = SMILES('[CH2][C](CC)C([O])[C]=O'),
    E0 = (311.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3100,440,815,1455,1000,1855,455,950,245.302,852.696,1962.41],'cm^-1')),
        HinderedRotor(inertia=(0.0954738,'amu*angstrom^2'), symmetry=1, barrier=(3.30219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954738,'amu*angstrom^2'), symmetry=1, barrier=(3.30219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954738,'amu*angstrom^2'), symmetry=1, barrier=(3.30219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954738,'amu*angstrom^2'), symmetry=1, barrier=(3.30219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954738,'amu*angstrom^2'), symmetry=1, barrier=(3.30219,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.392022,0.0695653,-5.70157e-05,2.47947e-08,-4.33447e-12,37556,37.1919], Tmin=(100,'K'), Tmax=(1373.9,'K')), NASAPolynomial(coeffs=[15.4856,0.0256219,-9.03935e-06,1.515e-09,-9.84406e-14,33408.5,-40.4207], Tmin=(1373.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(311.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([CH2])C([O])[C]=O(11037)',
    structure = SMILES('[CH2]C([CH2])C([O])[C]=O'),
    E0 = (384.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,3380.53,3389.2],'cm^-1')),
        HinderedRotor(inertia=(1.74462,'amu*angstrom^2'), symmetry=1, barrier=(78.2009,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00276181,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00275125,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.74069,'amu*angstrom^2'), symmetry=1, barrier=(78.3593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18529,0.058437,-5.08245e-05,1.81181e-08,-4.31018e-14,46295.4,32.9901], Tmin=(100,'K'), Tmax=(869.45,'K')), NASAPolynomial(coeffs=[12.5304,0.0192991,-5.82829e-06,8.88268e-10,-5.52106e-14,43829.1,-22.9945], Tmin=(869.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH]C)C([O])[C]=O(11211)',
    structure = SMILES('[CH2]C([CH]C)C([O])[C]=O'),
    E0 = (353.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,187.473,803.456,3603.26],'cm^-1')),
        HinderedRotor(inertia=(0.0668632,'amu*angstrom^2'), symmetry=1, barrier=(1.60886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0668632,'amu*angstrom^2'), symmetry=1, barrier=(1.60886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0668632,'amu*angstrom^2'), symmetry=1, barrier=(1.60886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0668632,'amu*angstrom^2'), symmetry=1, barrier=(1.60886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0668632,'amu*angstrom^2'), symmetry=1, barrier=(1.60886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752892,0.0698914,-6.23929e-05,3.08805e-08,-6.31837e-12,42583.7,38.2515], Tmin=(100,'K'), Tmax=(1159.99,'K')), NASAPolynomial(coeffs=[11.9112,0.0314146,-1.26384e-05,2.28591e-09,-1.55749e-13,39995,-17.2368], Tmin=(1159.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(CC)C([O])=[C][O](11662)',
    structure = SMILES('[CH2]C(CC)C([O])=[C][O]'),
    E0 = (228.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,435.225,435.225,435.226,2246.88],'cm^-1')),
        HinderedRotor(inertia=(0.00760102,'amu*angstrom^2'), symmetry=1, barrier=(1.02171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101431,'amu*angstrom^2'), symmetry=1, barrier=(13.634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101431,'amu*angstrom^2'), symmetry=1, barrier=(13.634,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101431,'amu*angstrom^2'), symmetry=1, barrier=(13.634,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.183871,0.0801714,-7.72818e-05,3.87755e-08,-7.58423e-12,27665.6,37.0784], Tmin=(100,'K'), Tmax=(1328.42,'K')), NASAPolynomial(coeffs=[18.2967,0.0207796,-5.9902e-06,8.75639e-10,-5.23563e-14,23086,-56.0845], Tmin=(1328.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(228.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2]CC([CH2])C([O])[C]=O(11775)',
    structure = SMILES('[CH2]CC([CH2])C([O])[C]=O'),
    E0 = (363.783,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,216.864,874.939,2060.13],'cm^-1')),
        HinderedRotor(inertia=(0.0983504,'amu*angstrom^2'), symmetry=1, barrier=(3.00774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0983504,'amu*angstrom^2'), symmetry=1, barrier=(3.00774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0983504,'amu*angstrom^2'), symmetry=1, barrier=(3.00774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0983504,'amu*angstrom^2'), symmetry=1, barrier=(3.00774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0983504,'amu*angstrom^2'), symmetry=1, barrier=(3.00774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.420745,0.0739616,-6.89109e-05,3.46818e-08,-7.03694e-12,43885.8,38.6727], Tmin=(100,'K'), Tmax=(1189.09,'K')), NASAPolynomial(coeffs=[14.2526,0.0274323,-1.02153e-05,1.77378e-09,-1.18174e-13,40596.4,-30.4534], Tmin=(1189.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[C](C)C([O])[C]=O(11776)',
    structure = SMILES('CC[C](C)C([O])[C]=O'),
    E0 = (106.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.556554,0.0657953,-3.75396e-05,3.70034e-09,2.75447e-12,12884.6,34.8155], Tmin=(100,'K'), Tmax=(1085.77,'K')), NASAPolynomial(coeffs=[14.5767,0.0307803,-1.2148e-05,2.22078e-09,-1.53926e-13,8859.5,-38.4927], Tmin=(1085.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ(C)CO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(CC)C(O)=[C][O](11666)',
    structure = SMILES('[CH2]C(CC)C(O)=[C][O]'),
    E0 = (90.8837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.951447,0.0887369,-8.52221e-05,4.09721e-08,-7.48766e-12,11126.3,38.3418], Tmin=(100,'K'), Tmax=(1511.31,'K')), NASAPolynomial(coeffs=[22.203,0.0169975,-3.64138e-06,4.07371e-10,-2.01622e-14,5321.74,-78.9764], Tmin=(1511.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.8837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]C(C)C([O])[C]=O(11086)',
    structure = SMILES('C[CH]C(C)C([O])[C]=O'),
    E0 = (147.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.677823,0.0688503,-5.18995e-05,2.0469e-08,-3.32627e-12,17922.6,36.7396], Tmin=(100,'K'), Tmax=(1427.08,'K')), NASAPolynomial(coeffs=[13.7876,0.0321044,-1.32757e-05,2.42559e-09,-1.65356e-13,14180.9,-31.1694], Tmin=(1427.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Cs_S) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC(C)C([O])=[C][O](11669)',
    structure = SMILES('CCC(C)C([O])=[C][O]'),
    E0 = (23.6062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,1685,370,432.266,432.278,432.289,4000],'cm^-1')),
        HinderedRotor(inertia=(0.112117,'amu*angstrom^2'), symmetry=1, barrier=(14.8672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112112,'amu*angstrom^2'), symmetry=1, barrier=(14.867,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11211,'amu*angstrom^2'), symmetry=1, barrier=(14.8672,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112113,'amu*angstrom^2'), symmetry=1, barrier=(14.8672,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.163579,0.0742236,-5.01106e-05,7.62772e-09,3.80584e-12,2986.29,34.0474], Tmin=(100,'K'), Tmax=(989.29,'K')), NASAPolynomial(coeffs=[17.2222,0.0262879,-9.32641e-06,1.63938e-09,-1.12897e-13,-1418.36,-53.27], Tmin=(989.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.6062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](CC)C(O)[C]=O(11777)',
    structure = SMILES('[CH2][C](CC)C(O)[C]=O'),
    E0 = (67.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.400487,0.0696754,-4.39825e-05,5.59167e-09,3.54871e-12,8241.5,36.3733], Tmin=(100,'K'), Tmax=(1005.35,'K')), NASAPolynomial(coeffs=[15.4851,0.0283956,-1.03492e-05,1.8275e-09,-1.25294e-13,4261.52,-41.1911], Tmin=(1005.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(67.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C](CC)C([O])C=O(11689)',
    structure = SMILES('[CH2][C](CC)C([O])C=O'),
    E0 = (151.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,244.6,865.688,1956.8],'cm^-1')),
        HinderedRotor(inertia=(0.0968052,'amu*angstrom^2'), symmetry=1, barrier=(3.32912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968052,'amu*angstrom^2'), symmetry=1, barrier=(3.32912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968052,'amu*angstrom^2'), symmetry=1, barrier=(3.32912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968052,'amu*angstrom^2'), symmetry=1, barrier=(3.32912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968052,'amu*angstrom^2'), symmetry=1, barrier=(3.32912,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.639831,0.064625,-3.54749e-05,2.68184e-09,2.99997e-12,18308,35.4782], Tmin=(100,'K'), Tmax=(1068.64,'K')), NASAPolynomial(coeffs=[13.5215,0.0323107,-1.24385e-05,2.23589e-09,-1.53378e-13,14646.8,-31.772], Tmin=(1068.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C)C(O)[C]=O(11778)',
    structure = SMILES('[CH2]C([CH]C)C(O)[C]=O'),
    E0 = (109.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.378158,0.0743826,-6.39763e-05,2.9525e-08,-5.54158e-12,13286,38.816], Tmin=(100,'K'), Tmax=(1270.71,'K')), NASAPolynomial(coeffs=[14.4637,0.0300437,-1.16368e-05,2.06572e-09,-1.39237e-13,9706.21,-32.5131], Tmin=(1270.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(109.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH]C)C([O])C=O(11690)',
    structure = SMILES('[CH2]C([CH]C)C([O])C=O'),
    E0 = (193.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,226.602,654.338,3305.15],'cm^-1')),
        HinderedRotor(inertia=(0.0417101,'amu*angstrom^2'), symmetry=1, barrier=(1.41979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0417101,'amu*angstrom^2'), symmetry=1, barrier=(1.41979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0417101,'amu*angstrom^2'), symmetry=1, barrier=(1.41979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0417101,'amu*angstrom^2'), symmetry=1, barrier=(1.41979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0417101,'amu*angstrom^2'), symmetry=1, barrier=(1.41979,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.753891,0.0677683,-5.01668e-05,1.99222e-08,-3.29988e-12,23346.3,37.4281], Tmin=(100,'K'), Tmax=(1392.38,'K')), NASAPolynomial(coeffs=[12.5261,0.0339491,-1.37334e-05,2.47783e-09,-1.67748e-13,20068.1,-23.2627], Tmin=(1392.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(C)C([O])[C]=O(11779)',
    structure = SMILES('[CH2]CC(C)C([O])[C]=O'),
    E0 = (158.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314337,0.0732165,-5.91763e-05,2.4972e-08,-4.25902e-12,19226.4,37.2791], Tmin=(100,'K'), Tmax=(1393.42,'K')), NASAPolynomial(coeffs=[15.8697,0.0285626,-1.1107e-05,1.97367e-09,-1.32786e-13,14891.3,-42.9273], Tmin=(1393.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]CC([CH2])C(O)[C]=O(11780)',
    structure = SMILES('[CH2]CC([CH2])C(O)[C]=O'),
    E0 = (120.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0263961,0.0786472,-7.10272e-05,3.38562e-08,-6.43328e-12,14589,39.3104], Tmin=(100,'K'), Tmax=(1275.15,'K')), NASAPolynomial(coeffs=[16.7635,0.0261448,-9.2671e-06,1.56712e-09,-1.02842e-13,10320.6,-45.5046], Tmin=(1275.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.049,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]CC([CH2])C([O])C=O(11691)',
    structure = SMILES('[CH2]CC([CH2])C([O])C=O'),
    E0 = (203.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2782.5,750,1395,475,1775,1000,261.41,916.392,2028.68],'cm^-1')),
        HinderedRotor(inertia=(0.0867738,'amu*angstrom^2'), symmetry=1, barrier=(3.4084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0867738,'amu*angstrom^2'), symmetry=1, barrier=(3.4084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0867738,'amu*angstrom^2'), symmetry=1, barrier=(3.4084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0867738,'amu*angstrom^2'), symmetry=1, barrier=(3.4084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0867738,'amu*angstrom^2'), symmetry=1, barrier=(3.4084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.391646,0.0721258,-5.74327e-05,2.44261e-08,-4.23605e-12,24650,37.9626], Tmin=(100,'K'), Tmax=(1366.02,'K')), NASAPolynomial(coeffs=[14.675,0.0303015,-1.15065e-05,2.01267e-09,-1.34108e-13,20747.7,-35.4013], Tmin=(1366.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (158.536,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (597.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (318.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (318.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (321.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (620.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (771.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (771.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (613.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (161.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (166.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (166.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (181.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (181.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (221.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (221.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (221.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (175.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (283.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (239.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (448.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (323.531,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (270.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (207.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (238.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (521.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (460.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (522.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (520.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (564.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (440.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (575.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (261.216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (272.508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (305.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (276.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (233.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (274.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (212.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (252.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (241.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (190.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (247.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['OCHCO(3676)', 'butene1(35)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2]C(C)C([O])[C]=O(10955)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['CCC[CH]C([O])[C]=O(11769)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['CC[CH]CC([O])[C]=O(10993)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['[CH2]C(CC)C([O])=C[O](11005)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', 'CC[CH]C([O])[C]=O(11593)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(4)', '[CH2]C([CH][C]=O)CC(2468)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[C]=O(2355)', '[CH2]C([CH][O])CC(1756)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C(CC)C([O])[C]=O(11770)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['[CH2]C(CC)C1OC1=O(11771)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['CCC1COC1[C]=O(10998)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['CCC1CC(=O)C1[O](11008)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['[CH2]C(CC)C(O)=C=O(11646)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['[CH2]C(CC)C(=O)C=O(10999)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['C=C(CC)C(O)[C]=O(11772)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['C=C(CC)C([O])C=O(11686)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['CCC(C)C([O])=C=O(11648)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CO(2039)', '[CH2]C([CH][O])CC(1756)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1099.39,'m^3/(mol*s)'), n=0.635039, Ea=(16.7529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;Y_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', 'C=C(CC)C([O])[C]=O(11773)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2]C(CC)C([O])=C=O(11659)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[C]=O(2355)', '[CH2]C(C=O)CC(1767)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CO_birad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['O(4)', '[CH2]C(C=C=O)CC(2465)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(11.358,'m^3/(mol*s)'), n=1.81033, Ea=(28.0982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;YJ] for rate rule [Cds-CsH_Ck;O_atom_triplet]
Euclidian distance = 2.23606797749979
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C2H5(14)', 'C=CC([O])[C]=O(11100)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1020,'cm^3/(mol*s)'), n=2.41, Ea=(27.3634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 418 used for Cds-CsH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O][C]=C[O](9592)', 'butene1(35)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['OCHCO(3676)', '[CH2][CH]CC(37)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C2H5(14)', '[CH2][CH]C([O])[C]=O(11469)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.36745e+07,'m^3/(mol*s)'), n=-0.263863, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00481396807501, var=0.0768145972539, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.567716674236
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][C]=C[O](9592)', '[CH2][CH]CC(37)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2][C](CC)C([O])[C]=O(11774)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH3(17)', '[CH2]C([CH2])C([O])[C]=O(11037)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.7349e+07,'m^3/(mol*s)'), n=-0.263863, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00481396807501, var=0.0768145972539, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.567716674236
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH2]C([CH]C)C([O])[C]=O(11211)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2]C(CC)C([O])=[C][O](11662)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2]CC([CH2])C([O])[C]=O(11775)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['CC[C](C)C([O])[C]=O(11776)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['[CH2]C(CC)C(O)=[C][O](11666)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['C[CH]C(C)C([O])[C]=O(11086)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['CCC(C)C([O])=[C][O](11669)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['[CH2][C](CC)C(O)[C]=O(11777)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['[CH2][C](CC)C([O])C=O(11689)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['[CH2]C([CH]C)C(O)[C]=O(11778)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(210000,'s^-1'), n=1.76, Ea=(53.8899,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 326 used for R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([CH]C)C([O])C=O(11690)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(0.164999,'s^-1'), n=3.35583, Ea=(59.1095,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_H/NonDeC;XH_out] for rate rule [R4H_SSS;C_rad_out_H/NonDeC;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]CC(C)C([O])[C]=O(11779)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    products = ['[CH2]CC([CH2])C(O)[C]=O(11780)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.2e+11,'s^-1'), n=0, Ea=(31.8402,'kJ/mol'), T0=(1,'K'), Tmin=(200,'K'), Tmax=(1000,'K'), comment="""From training reaction 306 used for R5H_CCC;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]CC([CH2])C([O])C=O(11691)'],
    products = ['[CH2]C(CC)C([O])[C]=O(10995)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(208392,'s^-1'), n=1.49392, Ea=(43.3506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_CCC;C_rad_out_2H;XH_out] for rate rule [R5H_CCC_O;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2650',
    isomers = [
        '[CH2]C(CC)C([O])[C]=O(10995)',
    ],
    reactants = [
        ('OCHCO(3676)', 'butene1(35)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2650',
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


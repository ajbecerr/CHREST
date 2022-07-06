species(
    label = 'C=C[CH]CC([O])=C[O](11260)',
    structure = SMILES('[CH2]C=CCC([O])=C[O]'),
    E0 = (60.1396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,371.32,371.337,371.345,371.374],'cm^-1')),
        HinderedRotor(inertia=(0.194783,'amu*angstrom^2'), symmetry=1, barrier=(19.0688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194838,'amu*angstrom^2'), symmetry=1, barrier=(19.0685,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194864,'amu*angstrom^2'), symmetry=1, barrier=(19.0684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.315737,0.0648577,-1.99561e-05,-3.38401e-08,2.14782e-11,7380.62,30.9617], Tmin=(100,'K'), Tmax=(937.203,'K')), NASAPolynomial(coeffs=[21.8558,0.0136192,-3.08064e-06,4.86701e-10,-3.7325e-14,1555.93,-81.0942], Tmin=(937.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(60.1396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C1[CH]C[C]([O])C1[O](12107)',
    structure = SMILES('[CH2]C1[CH]C[C]([O])C1[O]'),
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
    label = '[CH]=CCC([O])=C[O](12108)',
    structure = SMILES('[CH]=CCC([O])=C[O]'),
    E0 = (191.762,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,196.311,198.112,199.386],'cm^-1')),
        HinderedRotor(inertia=(0.64388,'amu*angstrom^2'), symmetry=1, barrier=(18.0167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.638576,'amu*angstrom^2'), symmetry=1, barrier=(17.9737,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907963,0.0550667,-2.37639e-05,-2.26507e-08,1.67981e-11,23187,27.1532], Tmin=(100,'K'), Tmax=(923.733,'K')), NASAPolynomial(coeffs=[19.7491,0.0071186,-5.28109e-07,2.75982e-12,-2.29247e-15,18271,-70.0182], Tmin=(923.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.762,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])CC=C[CH2](5449)',
    structure = SMILES('[CH]=C([O])CC=C[CH2]'),
    E0 = (374.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,460.616,461.315],'cm^-1')),
        HinderedRotor(inertia=(0.0930579,'amu*angstrom^2'), symmetry=1, barrier=(14.037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0929739,'amu*angstrom^2'), symmetry=1, barrier=(14.036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0927408,'amu*angstrom^2'), symmetry=1, barrier=(14.0428,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3874.43,'J/mol'), sigma=(6.44487,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.18 K, Pc=32.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.84923,0.0598081,-3.62764e-05,-5.38171e-10,5.92169e-12,45171.8,28.1803], Tmin=(100,'K'), Tmax=(978.941,'K')), NASAPolynomial(coeffs=[15.7465,0.0197907,-6.91215e-06,1.22487e-09,-8.56623e-14,41255.8,-48.4768], Tmin=(978.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_P)"""),
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
    label = '[CH]C=CCC([O])=C[O](12109)',
    structure = SMILES('[CH]C=CCC([O])=C[O]'),
    E0 = (279.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.31721,0.0672646,-2.81006e-05,-2.10477e-08,1.57009e-11,33740.1,31.7845], Tmin=(100,'K'), Tmax=(945.637,'K')), NASAPolynomial(coeffs=[19.4822,0.0198032,-6.12228e-06,1.03835e-09,-7.32635e-14,28612.9,-67.549], Tmin=(945.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=CCC1=COO1(12110)',
    structure = SMILES('[CH2]C=CCC1=COO1'),
    E0 = (300.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953082,0.0497606,1.12169e-05,-5.38413e-08,2.47164e-11,36233.4,28.9739], Tmin=(100,'K'), Tmax=(984.269,'K')), NASAPolynomial(coeffs=[18.0496,0.0213313,-8.01618e-06,1.5579e-09,-1.17078e-13,30879.4,-63.3368], Tmin=(984.269,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Allyl_P)"""),
)

species(
    label = '[O]C=C1CC=CCO1(11266)',
    structure = SMILES('[O]C=C1CC=CCO1'),
    E0 = (-143.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41647,0.0239978,0.000111288,-1.69966e-07,6.8231e-11,-17169.7,20.7867], Tmin=(100,'K'), Tmax=(956.637,'K')), NASAPolynomial(coeffs=[23.4871,0.0145895,-3.90888e-06,8.72369e-10,-8.07288e-14,-25184.6,-104.533], Tmin=(956.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-143.779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclohexane) + radical(C=COJ)"""),
)

species(
    label = '[O]C1=COCC=CC1(12111)',
    structure = SMILES('[O]C1=COCC=CC1'),
    E0 = (-117.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25634,-0.00733822,0.000270965,-4.05126e-07,1.73067e-10,-13941.4,25.1868], Tmin=(100,'K'), Tmax=(897.822,'K')), NASAPolynomial(coeffs=[49.2383,-0.0427241,3.20557e-05,-6.42918e-09,4.28529e-13,-29746.9,-241.166], Tmin=(897.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C=CCC(O)=C=O(12112)',
    structure = SMILES('[CH2]C=CCC(O)=C=O'),
    E0 = (-42.8652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.203831,0.0850441,-8.66494e-05,4.07591e-08,-6.27139e-12,-4997.49,28.0498], Tmin=(100,'K'), Tmax=(970.185,'K')), NASAPolynomial(coeffs=[19.6168,0.0180887,-5.95596e-06,9.95271e-10,-6.669e-14,-9538.25,-70.5543], Tmin=(970.185,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-42.8652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + radical(Allyl_P)"""),
)

species(
    label = 'C=C=CCC(O)=C[O](11272)',
    structure = SMILES('C=C=CCC(O)=C[O]'),
    E0 = (-52.5601,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.000295166,0.0708306,-2.95155e-05,-2.98903e-08,2.15345e-11,-6161.71,29.3663], Tmin=(100,'K'), Tmax=(931.881,'K')), NASAPolynomial(coeffs=[24.456,0.00953781,-1.16627e-06,1.28657e-10,-1.30661e-14,-12616.3,-97.0692], Tmin=(931.881,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.5601,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ)"""),
)

species(
    label = 'C=C=CCC([O])=CO(12113)',
    structure = SMILES('C=C=CCC([O])=CO'),
    E0 = (-56.2179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.29658,0.088245,-8.94641e-05,4.29257e-08,-7.59475e-12,-6546.43,34.1353], Tmin=(100,'K'), Tmax=(1632.55,'K')), NASAPolynomial(coeffs=[24.4274,0.00781651,4.22156e-07,-3.09003e-10,2.57809e-14,-12626.7,-95.4746], Tmin=(1632.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-56.2179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC=CCC([O])=C=O(12114)',
    structure = SMILES('CC=CCC(=O)[C]=O'),
    E0 = (-77.1956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24104,0.0614007,-4.17259e-05,1.3089e-08,-1.64292e-12,-9186.92,26.6207], Tmin=(100,'K'), Tmax=(1787.95,'K')), NASAPolynomial(coeffs=[16.1562,0.0280331,-1.37328e-05,2.65149e-09,-1.83527e-13,-14520.5,-54.0034], Tmin=(1787.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.1956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C[CH]C([O])=C[O](12115)',
    structure = SMILES('[CH2][CH]C[CH]C([O])=C[O]'),
    E0 = (391.429,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.396913,0.0707161,-5.52433e-05,1.39688e-08,2.01275e-12,47215.4,37.3448], Tmin=(100,'K'), Tmax=(959.689,'K')), NASAPolynomial(coeffs=[17.0496,0.0205074,-6.77662e-06,1.14742e-09,-7.79312e-14,43134.9,-46.9163], Tmin=(959.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(391.429,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJCO) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C[CH]C([O])[CH][O](12116)',
    structure = SMILES('[CH2]C=C[CH]C([O])[CH][O]'),
    E0 = (496.154,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314216,0.0867414,-0.000102668,6.82925e-08,-1.88125e-11,59801.4,30.1368], Tmin=(100,'K'), Tmax=(873.022,'K')), NASAPolynomial(coeffs=[11.1006,0.0373201,-1.7754e-05,3.44899e-09,-2.43716e-13,57918,-20.4367], Tmin=(873.022,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CCsJOH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]CC([O])[CH][O](12117)',
    structure = SMILES('[CH2]C=[C]CC([O])[CH][O]'),
    E0 = (617.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.211188,0.0910836,-0.000134157,1.14924e-07,-3.96834e-11,74346.5,33.762], Tmin=(100,'K'), Tmax=(810.825,'K')), NASAPolynomial(coeffs=[8.41372,0.0397141,-1.89523e-05,3.61568e-09,-2.49734e-13,73374.8,-1.88], Tmin=(810.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]C[C]([O])C[O](12118)',
    structure = SMILES('[CH2]C=[C]C[C]([O])C[O]'),
    E0 = (613.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271356,0.0917568,-0.000142527,1.28366e-07,-4.57008e-11,73900.9,33.4194], Tmin=(100,'K'), Tmax=(831.217,'K')), NASAPolynomial(coeffs=[6.68173,0.0425291,-2.05239e-05,3.91289e-09,-2.69201e-13,73470.2,7.49746], Tmin=(831.217,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([O])[CH][O](12119)',
    structure = SMILES('[CH2][C]=CCC([O])[CH][O]'),
    E0 = (617.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.211188,0.0910836,-0.000134157,1.14924e-07,-3.96834e-11,74346.5,33.762], Tmin=(100,'K'), Tmax=(810.825,'K')), NASAPolynomial(coeffs=[8.41372,0.0397141,-1.89523e-05,3.61568e-09,-2.49734e-13,73374.8,-1.88], Tmin=(810.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]CCC([O])=[C][O](12120)',
    structure = SMILES('[CH2][CH]CCC([O])=[C][O]'),
    E0 = (431.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.644108,0.074855,-8.03922e-05,4.81734e-08,-1.17747e-11,51989.9,37.0967], Tmin=(100,'K'), Tmax=(988.123,'K')), NASAPolynomial(coeffs=[11.724,0.0300028,-1.23055e-05,2.23681e-09,-1.52573e-13,49800.2,-16.225], Tmin=(988.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJC) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=CC[C]([O])C[O](12121)',
    structure = SMILES('[CH2][C]=CC[C]([O])C[O]'),
    E0 = (613.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271356,0.0917568,-0.000142527,1.28366e-07,-4.57008e-11,73900.9,33.4194], Tmin=(100,'K'), Tmax=(831.217,'K')), NASAPolynomial(coeffs=[6.68173,0.0425291,-2.05239e-05,3.91289e-09,-2.69201e-13,73470.2,7.49746], Tmin=(831.217,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C[CH]CC([O])=[C][O](12122)',
    structure = SMILES('[CH2]C[CH]CC([O])=[C][O]'),
    E0 = (431.283,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723838,0.0747693,-8.33369e-05,5.32432e-08,-1.40267e-11,51987,36.6427], Tmin=(100,'K'), Tmax=(915.224,'K')), NASAPolynomial(coeffs=[10.4426,0.0322924,-1.37179e-05,2.53015e-09,-1.73674e-13,50208,-9.3833], Tmin=(915.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.283,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJCC) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=CC[C]1OC1[O](12123)',
    structure = SMILES('[CH2]C=CC[C]1OC1[O]'),
    E0 = (297.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.798593,0.0690186,-6.64274e-05,3.83449e-08,-9.19279e-12,35938.8,29.7395], Tmin=(100,'K'), Tmax=(1003.83,'K')), NASAPolynomial(coeffs=[10.0075,0.032324,-1.15964e-05,1.93099e-09,-1.24214e-13,34089.9,-14.7237], Tmin=(1003.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Allyl_P)"""),
)

species(
    label = '[O]C=C([O])CC1[CH]C1(12124)',
    structure = SMILES('[O]C=C([O])CC1[CH]C1'),
    E0 = (171.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.798879,0.0515053,1.35526e-05,-6.54534e-08,3.19329e-11,20735.4,30.5052], Tmin=(100,'K'), Tmax=(938.848,'K')), NASAPolynomial(coeffs=[20.729,0.0142806,-3.16532e-06,5.21121e-10,-4.18415e-14,14891.4,-75.5818], Tmin=(938.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(C=C(C)OJ) + radical(C=COJ) + radical(cyclopropane)"""),
)

species(
    label = '[CH2]C=CCC1([O])[CH]O1(12125)',
    structure = SMILES('[CH2]C=CCC1([O])[CH]O1'),
    E0 = (295.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.337479,0.0653928,-1.57562e-05,-4.4115e-08,2.80653e-11,35724,28.2587], Tmin=(100,'K'), Tmax=(885.218,'K')), NASAPolynomial(coeffs=[21.5617,0.0129167,-4.26061e-07,-2.38567e-10,2.18201e-14,30264.9,-81.1592], Tmin=(885.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CCsJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1[CH]CC(=C[O])O1(11983)',
    structure = SMILES('[CH2]C1[CH]CC(=C[O])O1'),
    E0 = (156.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.745633,0.0410133,6.91843e-05,-1.41087e-07,6.33522e-11,19016.9,28.2153], Tmin=(100,'K'), Tmax=(921.753,'K')), NASAPolynomial(coeffs=[28.0358,0.00176706,4.19726e-06,-8.89673e-10,5.11621e-14,10622.3,-119.467], Tmin=(921.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[O][C]1CC=CCC1[O](12126)',
    structure = SMILES('[O][C]1CC=CCC1[O]'),
    E0 = (271.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06442,0.0499941,4.35081e-06,-4.13023e-08,1.90355e-11,32729.1,24.8106], Tmin=(100,'K'), Tmax=(1001.55,'K')), NASAPolynomial(coeffs=[15.518,0.0259658,-1.01288e-05,1.92777e-09,-1.40215e-13,28143.8,-53.3793], Tmin=(1001.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1[CH]CC([O])=CO1(12127)',
    structure = SMILES('[CH2]C1[CH]CC([O])=CO1'),
    E0 = (142.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.764554,0.0492483,3.04574e-05,-9.38446e-08,4.59392e-11,17317.8,24.8217], Tmin=(100,'K'), Tmax=(901.588,'K')), NASAPolynomial(coeffs=[23.4639,0.00734957,2.32219e-06,-6.9147e-10,4.73932e-14,10834.5,-95.5936], Tmin=(901.588,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.849,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(C=C(C)OJ) + radical(CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][CH]C1CC(=C[O])O1(12128)',
    structure = SMILES('[CH2][CH]C1CC(=C[O])O1'),
    E0 = (222.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627364,0.0470405,4.98183e-05,-1.22376e-07,5.7991e-11,26892.7,27.4129], Tmin=(100,'K'), Tmax=(904.304,'K')), NASAPolynomial(coeffs=[27.2237,0.00163173,5.32108e-06,-1.24015e-09,8.2387e-14,19129,-114.554], Tmin=(904.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = '[O][CH]C1([O])CC=CC1(12129)',
    structure = SMILES('[O][CH]C1([O])CC=CC1'),
    E0 = (282.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.782084,0.0634166,-4.54304e-05,1.63225e-08,-2.37084e-12,34141.3,27.081], Tmin=(100,'K'), Tmax=(1602.07,'K')), NASAPolynomial(coeffs=[15.5175,0.0266254,-1.0983e-05,1.98787e-09,-1.33938e-13,29419.9,-50.9535], Tmin=(1602.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH]C1CC([O])=CO1(12130)',
    structure = SMILES('[CH2][CH]C1CC([O])=CO1'),
    E0 = (139.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.999805,0.0414591,5.26227e-05,-1.18239e-07,5.52044e-11,16929.2,26.4027], Tmin=(100,'K'), Tmax=(899.436,'K')), NASAPolynomial(coeffs=[23.9561,0.00500787,3.94307e-06,-1.01759e-09,6.95122e-14,10144.6,-96.675], Tmin=(899.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(139.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = 'C=C=CCC([O])=C[O](12131)',
    structure = SMILES('C=C=CCC([O])=C[O]'),
    E0 = (85.2447,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,540,610,2055,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.910535,'amu*angstrom^2'), symmetry=1, barrier=(20.935,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.911346,'amu*angstrom^2'), symmetry=1, barrier=(20.9536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.35817,0.0669655,-3.74777e-05,-1.20928e-08,1.31306e-11,10395.7,29.5824], Tmin=(100,'K'), Tmax=(940.191,'K')), NASAPolynomial(coeffs=[20.7507,0.0128405,-3.1904e-06,5.12904e-10,-3.78787e-14,5118.8,-75.2125], Tmin=(940.191,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.2447,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C=CCC([O])=C=O(12132)',
    structure = SMILES('[CH2]C=CCC(=O)[C]=O'),
    E0 = (74.3036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,1855,455,950,458.995,459.025],'cm^-1')),
        HinderedRotor(inertia=(0.0974943,'amu*angstrom^2'), symmetry=1, barrier=(14.6033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0977036,'amu*angstrom^2'), symmetry=1, barrier=(14.6079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.798958,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.800149,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.851265,0.063183,-4.69141e-05,1.62101e-08,-2.21305e-12,9054.67,28.3245], Tmin=(100,'K'), Tmax=(1709.13,'K')), NASAPolynomial(coeffs=[18.5143,0.0218448,-1.0634e-05,2.05861e-09,-1.4306e-13,3016.97,-66.3568], Tmin=(1709.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.3036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C([O])[CH][O](2850)',
    structure = SMILES('[CH2]C([O])=C[O]'),
    E0 = (20.9566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,217.215,217.577],'cm^-1')),
        HinderedRotor(inertia=(0.665078,'amu*angstrom^2'), symmetry=1, barrier=(22.287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00575,0.032418,-6.0508e-07,-3.82059e-08,2.20969e-11,2603.17,17.8437], Tmin=(100,'K'), Tmax=(887.044,'K')), NASAPolynomial(coeffs=[16.9262,-0.00266727,4.27963e-06,-9.58521e-10,6.70145e-14,-1310.55,-59.4906], Tmin=(887.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.9566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=C(O)CJ)"""),
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
    label = '[CH2]C=C[CH]C([O])=C[O](12133)',
    structure = SMILES('[CH2]C=C[CH]C([O])=C[O]'),
    E0 = (177.056,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,186.932,187.479,187.527,187.739],'cm^-1')),
        HinderedRotor(inertia=(1.34917,'amu*angstrom^2'), symmetry=1, barrier=(33.3238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34202,'amu*angstrom^2'), symmetry=1, barrier=(33.3366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.34081,'amu*angstrom^2'), symmetry=1, barrier=(33.3364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.100243,0.0687758,-2.75542e-05,-3.02321e-08,2.12068e-11,21451,28.9345], Tmin=(100,'K'), Tmax=(935.722,'K')), NASAPolynomial(coeffs=[23.9776,0.00961343,-1.4976e-06,2.08833e-10,-1.91152e-14,15104,-94.711], Tmin=(935.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.056,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]CC([O])=C[O](12134)',
    structure = SMILES('[CH2]C=[C]CC([O])=C[O]'),
    E0 = (297.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,385.531,385.534,385.538,385.543],'cm^-1')),
        HinderedRotor(inertia=(0.170911,'amu*angstrom^2'), symmetry=1, barrier=(18.0274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170921,'amu*angstrom^2'), symmetry=1, barrier=(18.0274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170913,'amu*angstrom^2'), symmetry=1, barrier=(18.0273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.286021,0.0694228,-4.46494e-05,-4.63857e-09,1.051e-11,35983.8,31.5419], Tmin=(100,'K'), Tmax=(939.533,'K')), NASAPolynomial(coeffs=[20.651,0.0131447,-3.37265e-06,5.39047e-10,-3.89266e-14,30814.3,-72.583], Tmin=(939.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([O])=C[O](12135)',
    structure = SMILES('[CH2][C]=CCC([O])=C[O]'),
    E0 = (297.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,385.531,385.534,385.538,385.543],'cm^-1')),
        HinderedRotor(inertia=(0.170911,'amu*angstrom^2'), symmetry=1, barrier=(18.0274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170921,'amu*angstrom^2'), symmetry=1, barrier=(18.0274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170913,'amu*angstrom^2'), symmetry=1, barrier=(18.0273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.286021,0.0694228,-4.46494e-05,-4.63857e-09,1.051e-11,35983.8,31.5419], Tmin=(100,'K'), Tmax=(939.533,'K')), NASAPolynomial(coeffs=[20.651,0.0131447,-3.37265e-06,5.39047e-10,-3.89266e-14,30814.3,-72.583], Tmin=(939.533,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CCC([O])=[C][O](12136)',
    structure = SMILES('[CH2]C=CCC([O])=[C][O]'),
    E0 = (299.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,398.337,398.34,398.348,398.351],'cm^-1')),
        HinderedRotor(inertia=(0.00106242,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150557,'amu*angstrom^2'), symmetry=1, barrier=(16.9527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150553,'amu*angstrom^2'), symmetry=1, barrier=(16.9528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.512086,0.0667779,-4.74201e-05,5.71768e-09,4.91798e-12,36202.2,33.3853], Tmin=(100,'K'), Tmax=(966.351,'K')), NASAPolynomial(coeffs=[17.8062,0.0175913,-5.83914e-06,1.01757e-09,-7.13225e-14,31813.9,-54.8682], Tmin=(966.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'C3H3(5450)',
    structure = SMILES('[CH]=C=C'),
    E0 = (338.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,180,180,2603.58],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09172,0.0173333,-8.20209e-06,-2.63358e-09,2.66049e-12,40755.9,8.10965], Tmin=(100,'K'), Tmax=(946.054,'K')), NASAPolynomial(coeffs=[6.98214,0.0072721,-2.37773e-06,3.99152e-10,-2.69331e-14,39733.9,-11.9544], Tmin=(946.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C3H3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C=C([O])C[O](3384)',
    structure = SMILES('C=C([O])C[O]'),
    E0 = (1.93929,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,180,1873.6,2779.21],'cm^-1')),
        HinderedRotor(inertia=(0.568537,'amu*angstrom^2'), symmetry=1, barrier=(13.0718,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88356,0.0288238,-9.3742e-06,-2.8552e-08,3.06477e-11,269.127,17.4206], Tmin=(100,'K'), Tmax=(490.782,'K')), NASAPolynomial(coeffs=[4.65107,0.0216573,-9.59618e-06,1.80408e-09,-1.24992e-13,8.4506,9.26316], Tmin=(490.782,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1.93929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C)OJ)"""),
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
    label = 'buten3yl1(65)',
    structure = SMILES('[CH2]CC=C'),
    E0 = (191.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,345.872],'cm^-1')),
        HinderedRotor(inertia=(0.0750497,'amu*angstrom^2'), symmetry=1, barrier=(6.37545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0750659,'amu*angstrom^2'), symmetry=1, barrier=(6.37555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.486,0.0277404,-7.51225e-07,-1.39969e-08,6.14177e-12,23115.6,15.6915], Tmin=(100,'K'), Tmax=(1051.01,'K')), NASAPolynomial(coeffs=[7.36567,0.0210618,-8.19307e-06,1.49013e-09,-1.03098e-13,21433,-11.2176], Tmin=(1051.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""buten3yl1""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[C]=CCC([O])=C[O](12137)',
    structure = SMILES('C[C]=CCC([O])=C[O]'),
    E0 = (146.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.264425,0.0721889,-5.42381e-05,1.02236e-08,3.7743e-12,17761.2,31.3382], Tmin=(100,'K'), Tmax=(965.155,'K')), NASAPolynomial(coeffs=[18.5078,0.0189435,-6.2415e-06,1.07722e-09,-7.48506e-14,13198.1,-61.4242], Tmin=(965.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CCC([O])=[C]O(12138)',
    structure = SMILES('[CH2]C=CCC([O])=[C]O'),
    E0 = (158.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,377.931,377.963,377.981],'cm^-1')),
        HinderedRotor(inertia=(0.149185,'amu*angstrom^2'), symmetry=1, barrier=(15.1223,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149154,'amu*angstrom^2'), symmetry=1, barrier=(15.1218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149188,'amu*angstrom^2'), symmetry=1, barrier=(15.1219,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149181,'amu*angstrom^2'), symmetry=1, barrier=(15.1216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.719551,0.0832256,-8.31352e-05,4.03456e-08,-7.35237e-12,19241.2,36.4083], Tmin=(100,'K'), Tmax=(1534.2,'K')), NASAPolynomial(coeffs=[22.4144,0.0114455,-1.746e-06,1.08114e-10,-1.89841e-15,13492,-80.7039], Tmin=(1534.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(158.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=C[CH]C(O)=C[O](12139)',
    structure = SMILES('[CH2]C=C[CH]C(O)=C[O]'),
    E0 = (39.2514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.259521,0.0726635,-1.96725e-05,-4.7925e-08,2.95665e-11,4893.64,28.7252], Tmin=(100,'K'), Tmax=(930.115,'K')), NASAPolynomial(coeffs=[27.6915,0.00629597,5.35101e-07,-1.77446e-10,5.86642e-15,-2634.67,-116.616], Tmin=(930.115,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.2514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CCC(O)=[C][O](12140)',
    structure = SMILES('[CH2]C=CC[C](O)[C]=O'),
    E0 = (159.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,360,370,350,3000,3100,440,815,1455,1000,1855,455,950,180,823.963],'cm^-1')),
        HinderedRotor(inertia=(0.573989,'amu*angstrom^2'), symmetry=1, barrier=(13.1971,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574106,'amu*angstrom^2'), symmetry=1, barrier=(13.1998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121718,'amu*angstrom^2'), symmetry=1, barrier=(2.79853,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0274132,'amu*angstrom^2'), symmetry=1, barrier=(13.2005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0274147,'amu*angstrom^2'), symmetry=1, barrier=(13.2009,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406129,0.0745417,-7.17681e-05,3.56228e-08,-7.04861e-12,19266,34.3985], Tmin=(100,'K'), Tmax=(1221.15,'K')), NASAPolynomial(coeffs=[15.8037,0.0241048,-9.81313e-06,1.79905e-09,-1.23962e-13,15505.5,-42.9621], Tmin=(1221.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=[C]CC([O])=C[O](12141)',
    structure = SMILES('CC=[C]CC([O])=C[O]'),
    E0 = (146.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,307.101,307.525,307.726,307.956],'cm^-1')),
        HinderedRotor(inertia=(0.196125,'amu*angstrom^2'), symmetry=1, barrier=(13.1636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195168,'amu*angstrom^2'), symmetry=1, barrier=(13.1834,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19567,'amu*angstrom^2'), symmetry=1, barrier=(13.18,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.264425,0.0721889,-5.42381e-05,1.02236e-08,3.7743e-12,17761.2,31.3382], Tmin=(100,'K'), Tmax=(965.155,'K')), NASAPolynomial(coeffs=[18.5078,0.0189435,-6.2415e-06,1.07722e-09,-7.48506e-14,13198.1,-61.4242], Tmin=(965.155,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC(O)=C[O](12142)',
    structure = SMILES('[CH2]C=[C]CC(O)=C[O]'),
    E0 = (160.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,209.522,209.575,209.617],'cm^-1')),
        HinderedRotor(inertia=(0.664203,'amu*angstrom^2'), symmetry=1, barrier=(20.7099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664236,'amu*angstrom^2'), symmetry=1, barrier=(20.7096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664357,'amu*angstrom^2'), symmetry=1, barrier=(20.7098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664382,'amu*angstrom^2'), symmetry=1, barrier=(20.7096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.071377,0.0732818,-3.66634e-05,-2.24709e-08,1.89306e-11,19426.4,31.3242], Tmin=(100,'K'), Tmax=(930.896,'K')), NASAPolynomial(coeffs=[24.3557,0.00984314,-1.34922e-06,1.54973e-10,-1.41287e-14,13079.5,-94.4363], Tmin=(930.896,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=C[CH]C([O])=C[O](12143)',
    structure = SMILES('CC=C[CH]C([O])=C[O]'),
    E0 = (25.5571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0881544,0.0714293,-3.67511e-05,-1.58672e-08,1.46748e-11,3227.95,28.6968], Tmin=(100,'K'), Tmax=(949.595,'K')), NASAPolynomial(coeffs=[21.7858,0.0154949,-4.41407e-06,7.5825e-10,-5.59716e-14,-2491.74,-83.2785], Tmin=(949.595,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.5571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C=C[CH]C([O])=CO(12144)',
    structure = SMILES('[CH2]C=CC=C([O])[CH]O'),
    E0 = (27.3314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.578167,0.0587898,-6.08015e-06,-4.61807e-08,2.56048e-11,3425.59,31.3576], Tmin=(100,'K'), Tmax=(929.797,'K')), NASAPolynomial(coeffs=[20.5607,0.0149266,-3.23899e-06,4.82014e-10,-3.59067e-14,-2110.26,-73.3788], Tmin=(929.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.3314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C]=CCC(O)=C[O](12145)',
    structure = SMILES('[CH2][C]=CCC(O)=C[O]'),
    E0 = (160.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,209.522,209.575,209.617],'cm^-1')),
        HinderedRotor(inertia=(0.664203,'amu*angstrom^2'), symmetry=1, barrier=(20.7099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664236,'amu*angstrom^2'), symmetry=1, barrier=(20.7096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664357,'amu*angstrom^2'), symmetry=1, barrier=(20.7098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.664382,'amu*angstrom^2'), symmetry=1, barrier=(20.7096,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.071377,0.0732818,-3.66634e-05,-2.24709e-08,1.89306e-11,19426.4,31.3242], Tmin=(100,'K'), Tmax=(930.896,'K')), NASAPolynomial(coeffs=[24.3557,0.00984314,-1.34922e-06,1.54973e-10,-1.41287e-14,13079.5,-94.4363], Tmin=(930.896,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]CC([O])=CO(12146)',
    structure = SMILES('[CH2]C=[C]CC([O])=CO'),
    E0 = (156.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,225.502,225.767,225.898],'cm^-1')),
        HinderedRotor(inertia=(0.505563,'amu*angstrom^2'), symmetry=1, barrier=(18.2533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.504691,'amu*angstrom^2'), symmetry=1, barrier=(18.2541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.504113,'amu*angstrom^2'), symmetry=1, barrier=(18.2536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.505477,'amu*angstrom^2'), symmetry=1, barrier=(18.2525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23987,0.0892867,-9.21195e-05,4.50588e-08,-8.14762e-12,19035.8,35.6249], Tmin=(100,'K'), Tmax=(1585.84,'K')), NASAPolynomial(coeffs=[24.5756,0.00785444,3.3917e-07,-2.98046e-10,2.5526e-14,12899.8,-94.3543], Tmin=(1585.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CC=CCC([O])=[C][O](12147)',
    structure = SMILES('CC=CCC([O])=[C][O]'),
    E0 = (148.385,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,1685,370,299.576,299.582,299.586,299.587],'cm^-1')),
        HinderedRotor(inertia=(0.00187823,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161012,'amu*angstrom^2'), symmetry=1, barrier=(10.2548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161019,'amu*angstrom^2'), symmetry=1, barrier=(10.2546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263174,0.0721038,-6.53543e-05,3.04621e-08,-5.59976e-12,17989.6,34.0047], Tmin=(100,'K'), Tmax=(1324.18,'K')), NASAPolynomial(coeffs=[17.2318,0.0208458,-7.29014e-06,1.22922e-09,-8.06793e-14,13495.7,-52.6237], Tmin=(1324.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=CCC([O])=CO(12148)',
    structure = SMILES('[CH2][C]=CCC([O])=CO'),
    E0 = (156.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,225.502,225.767,225.898],'cm^-1')),
        HinderedRotor(inertia=(0.505563,'amu*angstrom^2'), symmetry=1, barrier=(18.2533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.504691,'amu*angstrom^2'), symmetry=1, barrier=(18.2541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.504113,'amu*angstrom^2'), symmetry=1, barrier=(18.2536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.505477,'amu*angstrom^2'), symmetry=1, barrier=(18.2525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.23987,0.0892867,-9.21195e-05,4.50588e-08,-8.14762e-12,19035.8,35.6249], Tmin=(100,'K'), Tmax=(1585.84,'K')), NASAPolynomial(coeffs=[24.5756,0.00785444,3.3917e-07,-2.98046e-10,2.5526e-14,12899.8,-94.3543], Tmin=(1585.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = '[CH]CC([O])=C[O](11435)',
    structure = SMILES('[CH]CC([O])=C[O]'),
    E0 = (287.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,496.654,496.661,496.669,496.695,496.714,496.744],'cm^-1')),
        HinderedRotor(inertia=(0.00068336,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810607,'amu*angstrom^2'), symmetry=1, barrier=(14.1871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34366,0.0467992,-1.95324e-05,-2.26557e-08,1.64173e-11,34697.9,23.4452], Tmin=(100,'K'), Tmax=(911.817,'K')), NASAPolynomial(coeffs=[18.415,0.00306204,1.17073e-06,-3.23509e-10,2.11275e-14,30289.7,-64.4395], Tmin=(911.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=CC1CC(=C[O])O1(11980)',
    structure = SMILES('C=CC1CC(=C[O])O1'),
    E0 = (-55.3249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.743993,0.0406273,7.31869e-05,-1.49094e-07,6.78479e-11,-6507.82,24.094], Tmin=(100,'K'), Tmax=(909.102,'K')), NASAPolynomial(coeffs=[28.6259,-0.000472965,6.39891e-06,-1.40936e-09,9.0964e-14,-14948.4,-126.304], Tmin=(909.102,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-55.3249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(C=COJ)"""),
)

species(
    label = 'C=CC1CC([O])=CO1(12149)',
    structure = SMILES('C=CC1CC([O])=CO1'),
    E0 = (-138.035,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11577,0.0350558,7.59462e-05,-1.44878e-07,6.5017e-11,-16471.3,23.086], Tmin=(100,'K'), Tmax=(905.279,'K')), NASAPolynomial(coeffs=[25.3546,0.00290911,5.01756e-06,-1.18603e-09,7.80259e-14,-23931.2,-108.404], Tmin=(905.279,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-138.035,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CC=CC(O)=C[O](12150)',
    structure = SMILES('[CH2]C=CC=C(O)C=O'),
    E0 = (-154.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.109043,0.0720613,-4.36917e-05,-4.88172e-09,9.55659e-12,-18381,27.0688], Tmin=(100,'K'), Tmax=(977.95,'K')), NASAPolynomial(coeffs=[20.9605,0.0171105,-5.93608e-06,1.09917e-09,-8.08441e-14,-23910,-80.479], Tmin=(977.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=CC=CC([O])=CO(12151)',
    structure = SMILES('C=CC=CC([O])=CO'),
    E0 = (-139.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-2.17546,0.0972663,-0.000103446,5.06318e-08,-8.94609e-12,-16484.7,33.4316], Tmin=(100,'K'), Tmax=(1687.41,'K')), NASAPolynomial(coeffs=[26.5363,0.00288802,3.84473e-06,-9.99533e-10,7.27613e-14,-22427.7,-109.006], Tmin=(1687.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-139.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=CCCC([O])=C=O(12152)',
    structure = SMILES('C=CCCC(=O)[C]=O'),
    E0 = (-71.1326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.960655,0.0721517,-7.58793e-05,4.73342e-08,-1.26934e-11,-8450.39,28.535], Tmin=(100,'K'), Tmax=(880.698,'K')), NASAPolynomial(coeffs=[8.44785,0.0381453,-1.79587e-05,3.48888e-09,-2.46969e-13,-9769.16,-6.63505], Tmin=(880.698,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.1326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C1CC(=O)C1[O](11361)',
    structure = SMILES('[CH2][CH]C1CC(=O)C1[O]'),
    E0 = (301.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62031,0.0375615,2.69329e-05,-5.78817e-08,2.3379e-11,36313.7,32.7603], Tmin=(100,'K'), Tmax=(1007.62,'K')), NASAPolynomial(coeffs=[12.9936,0.0277779,-1.11492e-05,2.14625e-09,-1.5655e-13,32226.3,-31.1043], Tmin=(1007.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[O]C=C1C[CH][CH]CO1(12077)',
    structure = SMILES('[O]C=C1C[CH][CH]CO1'),
    E0 = (132.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18468,0.0327995,7.97317e-05,-1.35534e-07,5.56347e-11,16011.7,26.0761], Tmin=(100,'K'), Tmax=(964.14,'K')), NASAPolynomial(coeffs=[23.2887,0.0139605,-4.32204e-06,9.72429e-10,-8.66837e-14,8362.73,-97.3188], Tmin=(964.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclohexanone) + radical(C=COJ) + radical(RCCJCC) + radical(CCJCO)"""),
)

species(
    label = '[O]C1=COC[CH][CH]C1(12153)',
    structure = SMILES('[O]C1=COC[CH][CH]C1'),
    E0 = (153.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20148,0.000459763,0.000236713,-3.61942e-07,1.55995e-10,18631.6,29.3294], Tmin=(100,'K'), Tmax=(895.704,'K')), NASAPolynomial(coeffs=[45.7783,-0.0370242,2.88858e-05,-5.85138e-09,3.92367e-13,4164.19,-217.001], Tmin=(895.704,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(CCJCO)"""),
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
    label = 'C=CC=CC([O])=C[O](12154)',
    structure = SMILES('C=CC=CC([O])=C[O]'),
    E0 = (2.27503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3005,3015,3025,975,983.333,991.667,1000,1300,1325,1350,1375,400,433.333,466.667,500,1630,1646.67,1663.33,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13612,'amu*angstrom^2'), symmetry=1, barrier=(26.1216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13721,'amu*angstrom^2'), symmetry=1, barrier=(26.1467,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0667809,0.069152,-2.78198e-05,-3.4974e-08,2.49132e-11,431.467,26.7628], Tmin=(100,'K'), Tmax=(910.876,'K')), NASAPolynomial(coeffs=[25.3998,0.00451996,1.8507e-06,-5.06797e-10,3.33173e-14,-6117.37,-103.704], Tmin=(910.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.27503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'C=CC[CH]C([O])=C[O](12006)',
    structure = SMILES('C=CC[CH]C([O])=C[O]'),
    E0 = (119.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,420.439,420.443,420.445,420.463,420.468],'cm^-1')),
        HinderedRotor(inertia=(0.000953822,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129195,'amu*angstrom^2'), symmetry=1, barrier=(16.206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129211,'amu*angstrom^2'), symmetry=1, barrier=(16.2063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44621,0.062241,-1.65539e-05,-3.55796e-08,2.17479e-11,14517.3,33.1206], Tmin=(100,'K'), Tmax=(937.352,'K')), NASAPolynomial(coeffs=[21.185,0.0138092,-3.1694e-06,5.0365e-10,-3.844e-14,8869.15,-74.9799], Tmin=(937.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJCO)"""),
)

species(
    label = 'C=[C]CCC([O])=C[O](12155)',
    structure = SMILES('C=[C]CCC([O])=C[O]'),
    E0 = (157.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,393.467,393.771,393.842,393.916,394.065],'cm^-1')),
        HinderedRotor(inertia=(0.127159,'amu*angstrom^2'), symmetry=1, barrier=(14.0137,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127345,'amu*angstrom^2'), symmetry=1, barrier=(14.0069,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127381,'amu*angstrom^2'), symmetry=1, barrier=(14.016,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.217116,0.0720796,-5.00944e-05,3.41636e-09,6.79544e-12,19084.1,31.9196], Tmin=(100,'K'), Tmax=(954.666,'K')), NASAPolynomial(coeffs=[19.4762,0.017528,-5.45793e-06,9.3052e-10,-6.53525e-14,14215.6,-66.3403], Tmin=(954.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCCC([O])=C[O](12156)',
    structure = SMILES('[CH]=CCCC([O])=C[O]'),
    E0 = (166.712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,342.133,342.211,342.217,342.223],'cm^-1')),
        HinderedRotor(inertia=(0.185452,'amu*angstrom^2'), symmetry=1, barrier=(15.4101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185449,'amu*angstrom^2'), symmetry=1, barrier=(15.4103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185388,'amu*angstrom^2'), symmetry=1, barrier=(15.4107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.18033,0.0709497,-4.17457e-05,-8.81457e-09,1.20247e-11,20200.3,31.9815], Tmin=(100,'K'), Tmax=(943.39,'K')), NASAPolynomial(coeffs=[20.8947,0.015215,-4.15815e-06,6.85621e-10,-4.93997e-14,14863.8,-74.3158], Tmin=(943.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=CCCC([O])=[C][O](12157)',
    structure = SMILES('C=CCCC([O])=[C][O]'),
    E0 = (159.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,393.65,393.656,393.656,393.659,393.662],'cm^-1')),
        HinderedRotor(inertia=(0.00108758,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103317,'amu*angstrom^2'), symmetry=1, barrier=(11.3623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103305,'amu*angstrom^2'), symmetry=1, barrier=(11.3624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.428373,0.0696004,-5.33915e-05,1.43556e-08,1.009e-12,19303.1,33.8166], Tmin=(100,'K'), Tmax=(1000.62,'K')), NASAPolynomial(coeffs=[16.7469,0.021782,-7.81475e-06,1.38339e-09,-9.56361e-14,15165.5,-49.2775], Tmin=(1000.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C[CH]CC(O)=C[O](12158)',
    structure = SMILES('[CH]C=CCC(O)=C[O]'),
    E0 = (141.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0389325,0.0711096,-2.00709e-05,-3.89264e-08,2.41353e-11,17182.6,31.5622], Tmin=(100,'K'), Tmax=(935.911,'K')), NASAPolynomial(coeffs=[23.1769,0.0165187,-4.10857e-06,6.56565e-10,-4.86546e-14,10882.4,-89.3456], Tmin=(935.911,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C[CH]CC([O])=CO(12159)',
    structure = SMILES('[CH]C=CCC([O])=CO'),
    E0 = (137.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.088429,0.0741536,-3.11534e-05,-2.784e-08,2.07568e-11,16742.8,31.8312], Tmin=(100,'K'), Tmax=(921.756,'K')), NASAPolynomial(coeffs=[22.7571,0.0160811,-3.47935e-06,4.79112e-10,-3.33022e-14,10786.6,-85.9873], Tmin=(921.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(137.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C[CH]CC(=O)C=O(11257)',
    structure = SMILES('[CH2]C=CCC(=O)C=O'),
    E0 = (-85.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.744489,0.0620956,-3.74344e-05,7.90661e-09,-4.44321e-14,-10177.2,27.9066], Tmin=(100,'K'), Tmax=(1418.61,'K')), NASAPolynomial(coeffs=[17.7459,0.0265174,-1.28842e-05,2.51117e-09,-1.75967e-13,-16244.5,-64.4436], Tmin=(1418.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-85.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(Allyl_P)"""),
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
    label = '[CH2]C1[CH]C[C]([CH][O])O1(12160)',
    structure = SMILES('[CH2]C1[CH]C[C]([CH][O])O1'),
    E0 = (565.883,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.591721,0.0781961,-9.88138e-05,7.4152e-08,-2.20331e-11,68179.7,28.8558], Tmin=(100,'K'), Tmax=(959.291,'K')), NASAPolynomial(coeffs=[8.67714,0.0335651,-1.1956e-05,1.92637e-09,-1.18819e-13,67130.7,-7.19752], Tmin=(959.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.883,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCOJ) + radical(CCJCO) + radical(C2CsJOCs) + radical(CJC(C)OC) + radical(CCsJOH)"""),
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
    label = '[CH2]C=CC[C]=O(2613)',
    structure = SMILES('[CH2]C=CC[C]=O'),
    E0 = (182.295,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.167415,'amu*angstrom^2'), symmetry=1, barrier=(17.9371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.042989,'amu*angstrom^2'), symmetry=1, barrier=(17.9464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0431769,'amu*angstrom^2'), symmetry=1, barrier=(17.9449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84068,0.0373103,-1.18252e-06,-2.21772e-08,9.72876e-12,22011.4,22.3283], Tmin=(100,'K'), Tmax=(1088.78,'K')), NASAPolynomial(coeffs=[12.2928,0.021254,-9.84393e-06,1.97423e-09,-1.4453e-13,18411.1,-35.0681], Tmin=(1088.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1CC=CCC1=O(11267)',
    structure = SMILES('[O]C1CC=CCC1=O'),
    E0 = (-60.6881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08925,0.017187,9.94429e-05,-1.32233e-07,4.76752e-11,-7208.64,22.1129], Tmin=(100,'K'), Tmax=(1024.31,'K')), NASAPolynomial(coeffs=[14.8108,0.0314275,-1.50147e-05,3.18338e-09,-2.44519e-13,-13168,-55.935], Tmin=(1024.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.6881,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclohexane) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C=CCC(=O)C[O](11273)',
    structure = SMILES('C=C=CCC(=O)C[O]'),
    E0 = (91.0486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11198,0.0640954,-4.71328e-05,1.66835e-08,-2.4049e-12,11053.3,27.223], Tmin=(100,'K'), Tmax=(1570.49,'K')), NASAPolynomial(coeffs=[14.5388,0.0298977,-1.44702e-05,2.8184e-09,-1.97772e-13,6835.9,-43.6145], Tmin=(1570.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.0486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C=[C]C[C](O)[CH][O](12161)',
    structure = SMILES('[CH2]C=[C]C[C](O)[CH][O]'),
    E0 = (563.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.15278,0.104381,-0.000179132,1.64732e-07,-5.78855e-11,67892,34.6832], Tmin=(100,'K'), Tmax=(863.208,'K')), NASAPolynomial(coeffs=[7.39985,0.0410632,-1.98943e-05,3.74513e-09,-2.53548e-13,67643.1,5.46848], Tmin=(863.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C2CsJOH) + radical(CCsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC[C](O)[CH][O](12162)',
    structure = SMILES('[CH2][C]=CC[C](O)[CH][O]'),
    E0 = (563.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.15278,0.104381,-0.000179132,1.64732e-07,-5.78855e-11,67892,34.6832], Tmin=(100,'K'), Tmax=(863.208,'K')), NASAPolynomial(coeffs=[7.39985,0.0410632,-1.98943e-05,3.74513e-09,-2.53548e-13,67643.1,5.46848], Tmin=(863.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C2CsJOH) + radical(CCsJOH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC[C]1[CH]OO1(12163)',
    structure = SMILES('[CH2]C=CC[C]1[CH]OO1'),
    E0 = (517.632,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27264,0.065753,-5.96615e-05,3.5035e-08,-9.57681e-12,62349.9,24.863], Tmin=(100,'K'), Tmax=(832.917,'K')), NASAPolynomial(coeffs=[5.74214,0.044289,-2.10073e-05,4.0966e-09,-2.90765e-13,61605.4,4.1173], Tmin=(832.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(517.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(12dioxetane) + radical(C2CsJOO) + radical(CCsJOO) + radical(Allyl_P)"""),
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
    label = 'C=C(O)[CH][O](12165)',
    structure = SMILES('[CH2]C(O)=C[O]'),
    E0 = (-116.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (72.0627,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6405,0.0363738,7.0227e-06,-5.55458e-08,3.02947e-11,-13954,17.6539], Tmin=(100,'K'), Tmax=(892.949,'K')), NASAPolynomial(coeffs=[20.6587,-0.00601666,6.33098e-06,-1.34924e-09,9.23671e-14,-19056.8,-81.4995], Tmin=(892.949,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-116.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C=C[CH]C(=O)C[O](12166)',
    structure = SMILES('[CH2]C=CC=C([O])C[O]'),
    E0 = (135.741,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3000,3100,440,815,1455,1000,296.324,296.331,296.333,3877.91],'cm^-1')),
        HinderedRotor(inertia=(0.00191984,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08733,'amu*angstrom^2'), symmetry=1, barrier=(67.753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08734,'amu*angstrom^2'), symmetry=1, barrier=(67.7542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799253,0.0662911,-5.53459e-05,2.51427e-08,-4.69439e-12,16444.5,30.8278], Tmin=(100,'K'), Tmax=(1268.42,'K')), NASAPolynomial(coeffs=[12.6003,0.0290761,-1.13363e-05,2.01173e-09,-1.35378e-13,13450.8,-28.9113], Tmin=(1268.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(135.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=[C]CC(=O)C[O](12167)',
    structure = SMILES('[CH2]C=[C]CC(=O)C[O]'),
    E0 = (303.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,180,180,643.283,644.056],'cm^-1')),
        HinderedRotor(inertia=(0.00999803,'amu*angstrom^2'), symmetry=1, barrier=(2.94098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00989943,'amu*angstrom^2'), symmetry=1, barrier=(2.91002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663441,'amu*angstrom^2'), symmetry=1, barrier=(15.2538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.051628,'amu*angstrom^2'), symmetry=1, barrier=(15.2527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15701,0.0652153,-4.98816e-05,1.88128e-08,-2.93652e-12,36636.3,28.76], Tmin=(100,'K'), Tmax=(1444.65,'K')), NASAPolynomial(coeffs=[12.9861,0.0324628,-1.58746e-05,3.11972e-09,-2.20825e-13,33218.5,-32.6603], Tmin=(1444.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=OCOJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC(=O)C[O](12168)',
    structure = SMILES('[CH2][C]=CCC(=O)C[O]'),
    E0 = (303.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,180,180,643.283,644.056],'cm^-1')),
        HinderedRotor(inertia=(0.00999803,'amu*angstrom^2'), symmetry=1, barrier=(2.94098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00989943,'amu*angstrom^2'), symmetry=1, barrier=(2.91002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663441,'amu*angstrom^2'), symmetry=1, barrier=(15.2538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.051628,'amu*angstrom^2'), symmetry=1, barrier=(15.2527,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15701,0.0652153,-4.98816e-05,1.88128e-08,-2.93652e-12,36636.3,28.76], Tmin=(100,'K'), Tmax=(1444.65,'K')), NASAPolynomial(coeffs=[12.9861,0.0324628,-1.58746e-05,3.11972e-09,-2.20825e-13,33218.5,-32.6603], Tmin=(1444.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=OCOJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=COC(=C)[CH][O](12169)',
    structure = SMILES('[CH2]C(=C[O])O[CH]C=C'),
    E0 = (84.7863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.140853,0.070121,-3.41508e-05,-1.82657e-08,1.54932e-11,10349.8,30.536], Tmin=(100,'K'), Tmax=(948.483,'K')), NASAPolynomial(coeffs=[21.6763,0.015272,-4.29675e-06,7.36038e-10,-5.4511e-14,4646.59,-80.7506], Tmin=(948.483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.7863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CCJ(O)C) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])([O])C=O(12170)',
    structure = SMILES('[CH2]C([O])([CH]C=C)C=O'),
    E0 = (243.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,180,180,180,296.043,1437.95,1600,2933.33,3200],'cm^-1')),
        HinderedRotor(inertia=(0.139013,'amu*angstrom^2'), symmetry=1, barrier=(3.19619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139013,'amu*angstrom^2'), symmetry=1, barrier=(3.19619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139013,'amu*angstrom^2'), symmetry=1, barrier=(3.19619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139013,'amu*angstrom^2'), symmetry=1, barrier=(3.19619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0886333,0.0919881,-0.00011225,7.1777e-08,-1.82534e-11,29446.2,30.3133], Tmin=(100,'K'), Tmax=(960.316,'K')), NASAPolynomial(coeffs=[15.3013,0.0278838,-1.21188e-05,2.26387e-09,-1.56791e-13,26490.4,-43.311], Tmin=(960.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)(C=O)OJ) + radical(C=CCJCO) + radical(CJC(C)(C=O)O)"""),
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
    label = '[CH2][CH][CH]C[C]=O(2615)',
    structure = SMILES('[CH2][CH][CH]C[C]=O'),
    E0 = (508.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,1692.39,1692.4],'cm^-1')),
        HinderedRotor(inertia=(0.00212538,'amu*angstrom^2'), symmetry=1, barrier=(4.31964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187872,'amu*angstrom^2'), symmetry=1, barrier=(4.31954,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187877,'amu*angstrom^2'), symmetry=1, barrier=(4.31967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187875,'amu*angstrom^2'), symmetry=1, barrier=(4.31961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86177,0.0523736,-7.03961e-05,6.25011e-08,-2.25463e-11,61250.9,29.3923], Tmin=(100,'K'), Tmax=(839.706,'K')), NASAPolynomial(coeffs=[3.94104,0.0317024,-1.42378e-05,2.64608e-09,-1.80066e-13,61281.3,21.9845], Tmin=(839.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJC) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1(C=O)CC=CC1(12171)',
    structure = SMILES('[O]C1(C=O)CC=CC1'),
    E0 = (-22.6089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27066,0.0482604,-2.89003e-06,-2.87837e-08,1.3805e-11,-2610.66,26.2594], Tmin=(100,'K'), Tmax=(1012.62,'K')), NASAPolynomial(coeffs=[13.3381,0.026988,-1.04792e-05,1.95489e-09,-1.39301e-13,-6407.89,-38.7924], Tmin=(1012.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.6089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopentene) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'CC=CC=C([O])C=O(12172)',
    structure = SMILES('C[CH]C=CC(=O)C=O'),
    E0 = (-100.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3728,0.0585838,-3.77916e-05,1.18095e-08,-1.51568e-12,-11948.5,25.4062], Tmin=(100,'K'), Tmax=(1712.18,'K')), NASAPolynomial(coeffs=[12.8239,0.0318318,-1.43549e-05,2.68409e-09,-1.83257e-13,-15869.8,-35.9968], Tmin=(1712.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-100.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(Allyl_S)"""),
)

species(
    label = 'C=C=CCC([O])C=O(11271)',
    structure = SMILES('C=C=CCC([O])C=O'),
    E0 = (91.3288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654487,0.0649196,-5.01871e-05,1.9596e-08,-3.06843e-12,11111.9,32.8297], Tmin=(100,'K'), Tmax=(1512.5,'K')), NASAPolynomial(coeffs=[16.0969,0.0240805,-9.68568e-06,1.74418e-09,-1.1774e-13,6440.55,-48.0604], Tmin=(1512.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.3288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=OCOJ)"""),
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
    label = 'C=C([O])C=O(2859)',
    structure = SMILES('[CH2]C(=O)C=O'),
    E0 = (-74.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.43166,'amu*angstrom^2'), symmetry=1, barrier=(9.92472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231641,'amu*angstrom^2'), symmetry=1, barrier=(31.0079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.689,0.0301409,-2.32871e-05,9.30138e-09,-1.56938e-12,-8891.4,16.525], Tmin=(100,'K'), Tmax=(1335.25,'K')), NASAPolynomial(coeffs=[7.34752,0.0161853,-7.60952e-06,1.4738e-09,-1.03803e-13,-10135.4,-7.29647], Tmin=(1335.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJCC=O)"""),
)

species(
    label = '[CH2]C=C[CH]C([O])C=O(12173)',
    structure = SMILES('[CH2]C=C[CH]C([O])C=O'),
    E0 = (183.14,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,462.538,462.581,462.699],'cm^-1')),
        HinderedRotor(inertia=(0.000788108,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192083,'amu*angstrom^2'), symmetry=1, barrier=(29.1972,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192347,'amu*angstrom^2'), symmetry=1, barrier=(29.1959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192262,'amu*angstrom^2'), symmetry=1, barrier=(29.1967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.483388,0.0658856,-3.81945e-05,1.00755e-10,4.99918e-12,22163.2,31.8589], Tmin=(100,'K'), Tmax=(1048.47,'K')), NASAPolynomial(coeffs=[17.0477,0.0244057,-9.91712e-06,1.87399e-09,-1.33713e-13,17496.2,-54.5301], Tmin=(1048.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=[C]CC([O])C=O(12174)',
    structure = SMILES('[CH2]C=[C]CC([O])C=O'),
    E0 = (304.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1685,370,180,885.52,885.567],'cm^-1')),
        HinderedRotor(inertia=(0.133749,'amu*angstrom^2'), symmetry=1, barrier=(3.07516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13376,'amu*angstrom^2'), symmetry=1, barrier=(3.07541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.656707,'amu*angstrom^2'), symmetry=1, barrier=(15.099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0271363,'amu*angstrom^2'), symmetry=1, barrier=(15.099,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71901,0.0658763,-5.25896e-05,2.14747e-08,-3.5432e-12,36693.8,34.291], Tmin=(100,'K'), Tmax=(1431.46,'K')), NASAPolynomial(coeffs=[15.0379,0.0258644,-1.06619e-05,1.948e-09,-1.32926e-13,32594.4,-39.9254], Tmin=(1431.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CCC([O])C=O(12175)',
    structure = SMILES('[CH2][C]=CCC([O])C=O'),
    E0 = (304.065,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1685,370,336.607,336.61,336.611],'cm^-1')),
        HinderedRotor(inertia=(0.013995,'amu*angstrom^2'), symmetry=1, barrier=(22.6015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281098,'amu*angstrom^2'), symmetry=1, barrier=(22.6015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00199061,'amu*angstrom^2'), symmetry=1, barrier=(22.6015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0641226,'amu*angstrom^2'), symmetry=1, barrier=(5.15576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71901,0.0658763,-5.25896e-05,2.14747e-08,-3.5432e-12,36693.8,34.291], Tmin=(100,'K'), Tmax=(1431.46,'K')), NASAPolynomial(coeffs=[15.0379,0.0258644,-1.06619e-05,1.948e-09,-1.32926e-13,32594.4,-39.9254], Tmin=(1431.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = 'C=CCC=C([O])C=O(12176)',
    structure = SMILES('C=CC[CH]C(=O)C=O'),
    E0 = (-31.1912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33185,0.0629816,-4.96885e-05,2.13375e-08,-3.98414e-12,-3659.13,28.8637], Tmin=(100,'K'), Tmax=(1199.44,'K')), NASAPolynomial(coeffs=[9.08665,0.0371203,-1.73468e-05,3.36158e-09,-2.37422e-13,-5519.42,-9.95903], Tmin=(1199.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-31.1912,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(CCJCC=O)"""),
)

species(
    label = 'C=CC=CC([O])C=O(12177)',
    structure = SMILES('C=CC=CC([O])C=O'),
    E0 = (25.1675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738409,0.061061,-3.26407e-05,-3.61901e-09,6.25023e-12,3153.7,31.379], Tmin=(100,'K'), Tmax=(1021.13,'K')), NASAPolynomial(coeffs=[15.8841,0.0233445,-8.98499e-06,1.66464e-09,-1.18093e-13,-1066.23,-47.5243], Tmin=(1021.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.1675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = '[CH]=C[CH]CC([O])C=O(12178)',
    structure = SMILES('[CH]C=CCC([O])C=O'),
    E0 = (285.409,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,2782.5,750,1395,475,1775,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.575433,0.0656338,-4.21128e-05,1.21278e-08,-1.04326e-12,34458,35.171], Tmin=(100,'K'), Tmax=(1300.69,'K')), NASAPolynomial(coeffs=[14.5631,0.0314197,-1.28068e-05,2.30985e-09,-1.56196e-13,30074.8,-38.851], Tmin=(1300.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.409,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(AllylJ2_triplet)"""),
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
    label = 'C=CC=CC(=O)C[O](12179)',
    structure = SMILES('C=CC=CC(=O)C[O]'),
    E0 = (23.7577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.988565,0.0669417,-5.69122e-05,2.5584e-08,-4.78305e-12,2965.06,27.6011], Tmin=(100,'K'), Tmax=(1244.5,'K')), NASAPolynomial(coeffs=[12.0025,0.0315417,-1.42449e-05,2.72779e-09,-1.91639e-13,223.674,-27.9439], Tmin=(1244.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(23.7577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = '[CH]=C[CH]CC(=O)C[O](12180)',
    structure = SMILES('[CH]C=CCC(=O)C[O]'),
    E0 = (285.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.927277,0.06593,-4.2364e-05,1.26488e-08,-1.51131e-12,34404.1,29.9511], Tmin=(100,'K'), Tmax=(1875.46,'K')), NASAPolynomial(coeffs=[17.397,0.0308029,-1.42691e-05,2.66193e-09,-1.80041e-13,28226.4,-59.8629], Tmin=(1875.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=OCOJ) + radical(AllylJ2_triplet)"""),
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
    E0 = (60.1396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (624.017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (894.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (629.564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (893.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (491.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (300.221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (67.6708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (68.2286,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (123.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (85.1128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (85.1128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (85.1128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (414.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (519.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (706.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (652.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (642.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (456.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (638.382,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (456.256,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (297.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (286.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (296.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (156.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (271.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (144.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (222.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (282.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (139.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (308.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (324.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (244.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (397.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (479.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (388.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (509.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (509.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (511.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (416.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (321.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (267.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (321.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (230.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (352.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (308.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (204.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (177.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (174.572,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (204.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (189.559,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (229.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (241.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (276.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (629.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (68.0473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (67.6708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (138.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (77.9216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (99.3646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (301.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (132.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (153.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (375.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (224.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (305.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (235.587,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (315.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (311.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (203.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (293.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (289.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (60.1396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (389.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (622.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (668.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (67.6708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS78',
    E0 = (91.1087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS79',
    E0 = (602.571,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS80',
    E0 = (588.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS81',
    E0 = (517.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS82',
    E0 = (231.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS83',
    E0 = (325.969,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS84',
    E0 = (301.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS85',
    E0 = (404.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS86',
    E0 = (348.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS87',
    E0 = (398.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS88',
    E0 = (400.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS89',
    E0 = (596.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS90',
    E0 = (67.2524,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS91',
    E0 = (77.9216,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS92',
    E0 = (99.3646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS93',
    E0 = (329.366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS94',
    E0 = (338.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS95',
    E0 = (378.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS96',
    E0 = (446.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS97',
    E0 = (348.374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS98',
    E0 = (69.7165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS99',
    E0 = (87.9924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS100',
    E0 = (87.9924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS101',
    E0 = (318.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS102',
    E0 = (68.0473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS103',
    E0 = (138.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS104',
    E0 = (436.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['OCHCO(3676)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C1[CH]C[C]([O])C1[O](12107)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.43608e+09,'s^-1'), n=0.0758676, Ea=(56.4793,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH2]C=CC[C]=C[O](10035)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(20)', '[CH]=CCC([O])=C[O](12108)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH]=C([O])CC=C[CH2](5449)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH]C=CCC([O])=C[O](12109)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C=CCC1=COO1(12110)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(240.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 237.9 to 240.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[O]C=C1CC=CCO1(11266)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[O]C1=COCC=CC1(12111)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R7;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C=CCC(O)=C=O(12112)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=C=CCC(O)=C[O](11272)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=C=CCC([O])=CO(12113)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['CC=CCC([O])=C=O(12114)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C[CH]C([O])=C[O](12115)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C=C[CH]C([O])[CH][O](12116)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C=[C]CC([O])[CH][O](12117)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C=[C]C[C]([O])C[O](12118)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.96625e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=CCC([O])[CH][O](12119)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]CCC([O])=[C][O](12120)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C]=CC[C]([O])C[O](12121)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C[CH]CC([O])=[C][O](12122)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C=CC[C]1OC1[O](12123)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(237.833,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[O]C=C([O])CC1[CH]C1(12124)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C=CCC1([O])[CH]O1(12125)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(236.077,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C1[CH]CC(=C[O])O1(11983)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(96.7661,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 94.0 to 96.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[O][C]1CC=CCC1[O](12126)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.68243e+09,'s^-1'), n=0.4695, Ea=(211,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_SDS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 208.8 to 211.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C1[CH]CC([O])=CO1(12127)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.041e+08,'s^-1'), n=0.7, Ea=(84.7093,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SDS_D;doublebond_intra_pri;radadd_intra] for rate rule [R6_SDS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2][CH]C1CC(=C[O])O1(12128)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(162.238,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 160.2 to 162.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[O][CH]C1([O])CC=CC1(12129)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(222.715,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_SMS_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 222.0 to 222.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2][CH]C1CC([O])=CO1(12130)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(79.528,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS_D;doublebond_intra;radadd_intra] for rate rule [R6_SMS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 79.3 to 79.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', 'C=C=CCC([O])=C[O](12131)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2]C=CCC([O])=C=O(12132)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['OCHCO(3676)', '[CH2][CH]C=C(3743)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.000357827,'m^3/(mol*s)'), n=2.74787, Ea=(45.8259,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cd_R;CsJ-CdHH] for rate rule [Ck_O;CsJ-CdHH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=C([O])[CH][O](2850)', '[CH]C=C(3739)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=-5.80997e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O][C]=C[O](9592)', '[CH2][CH]C=C(3743)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.62e+06,'m^3/(mol*s)'), n=-5.80997e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH2]C=C[CH]C([O])=C[O](12133)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)', '[CH2]C=[C]CC([O])=C[O](12134)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(3)', '[CH2][C]=CCC([O])=C[O](12135)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['H(3)', '[CH2]C=CCC([O])=[C][O](12136)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C3H3(5450)', 'C=C([O])C[O](3384)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.47269e+12,'s^-1'), n=-0.334211, Ea=(356.155,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, correlation='Root_1R!H->C_2R!H->C_N-5R!H->O_Ext-5C-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_2R!H->C_N-5R!H->O_Ext-5C-R
    Total Standard Deviation in ln(k): 11.540182761524994
Exact match found for rate rule [Root_1R!H->C_2R!H->C_N-5R!H->O_Ext-5C-R]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['OCCO(T)(10910)', 'buten3yl1(65)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.47269e+12,'s^-1'), n=-0.334211, Ea=(260.965,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.13686319048999, Tref=1000.0, N=1, correlation='Root_1R!H->C_2R!H->C_N-5R!H->O_Ext-5C-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_2R!H->C_N-5R!H->O_Ext-5C-R
    Total Standard Deviation in ln(k): 11.540182761524994
Exact match found for rate rule [Root_1R!H->C_2R!H->C_N-5R!H->O_Ext-5C-R]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C[C]=CCC([O])=C[O](12137)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C=CCC([O])=[C]O(12138)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C=C[CH]C(O)=C[O](12139)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(587605,'s^-1'), n=2.09, Ea=(169.87,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_2Cd;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C=CCC(O)=[C][O](12140)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CC=[C]CC([O])=C[O](12141)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C=[C]CC(O)=C[O](12142)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Cd_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['CC=C[CH]C([O])=C[O](12143)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(256000,'s^-1'), n=2, Ea=(117.57,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 95 used for R4H_SDS;C_rad_out_2H;Cs_H_out_H/Cd
Exact match found for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C=C[CH]C([O])=CO(12144)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.92799e+06,'s^-1'), n=1.7075, Ea=(114.432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4H_SDS;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2][C]=CCC(O)=C[O](12145)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C=[C]CC([O])=CO(12146)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSMS;Cd_rad_out;XH_out] for rate rule [R5H_SSMS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['CC=CCC([O])=[C][O](12147)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(37753.8,'s^-1'), n=1.925, Ea=(80.7512,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;Cs_H_out_2H] for rate rule [R6H_RSSMS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2][C]=CCC([O])=CO(12148)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(16140,'s^-1'), n=1.92259, Ea=(84.7802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6H_RSSMS;Y_rad_out;XH_out] for rate rule [R6H_RSSMS;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C=C)C([O])=C[O](11262)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C2H3(60)', '[CH]CC([O])=C[O](11435)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction56',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CC1CC(=C[O])O1(11980)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CC1CC([O])=CO1(12149)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R5_SSDS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CC=CC(O)=C[O](12150)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction59',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CC=CC([O])=CO(12151)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction60',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CCCC([O])=C=O(12152)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction61',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2][CH]C1CC(=O)C1[O](11361)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(3.78932e+07,'s^-1'), n=1.19089, Ea=(240.969,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_Cs_HH_D;doublebond_intra;radadd_intra_cs] for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 236.8 to 241.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction62',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[O]C=C1C[CH][CH]CO1(12077)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(71.9272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 67.9 to 71.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction63',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[O]C1=COC[CH][CH]C1(12153)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2.76476e+07,'s^-1'), n=0.815689, Ea=(93.4549,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R7_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 85.9 to 93.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction64',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CC1CC1([O])[CH][O](12030)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(315.455,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHCd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction65',
    reactants = ['H(3)', 'C=CC=CC([O])=C[O](12154)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(4.61283,'m^3/(mol*s)'), n=2.04274, Ea=(10.5229,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 101 used for Cds-CdH_Cds-CdH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CdH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction66',
    reactants = ['[O][C]=C[O](9592)', 'butadiene13(1350)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(0.0534234,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction67',
    reactants = ['C=CC[CH]C([O])=C[O](12006)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.169e+11,'s^-1'), n=0.707, Ea=(116.068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction68',
    reactants = ['C=[C]CCC([O])=C[O](12155)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH]=CCCC([O])=C[O](12156)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction70',
    reactants = ['C=CCCC([O])=[C][O](12157)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction71',
    reactants = ['[CH]=C[CH]CC(O)=C[O](12158)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH]=C[CH]CC([O])=CO(12159)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction73',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=C[CH]CC(=O)C=O(11257)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction74',
    reactants = ['C=C[CH]CC([O])[C]=O(11249)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction75',
    reactants = ['[CH2]C1[CH]C[C]([CH][O])O1(12160)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(6.43608e+09,'s^-1'), n=0.0758676, Ea=(56.4793,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[CH][O](1548)', '[CH2]C=CC[C]=O(2613)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction77',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[O]C1CC=CCC1=O(11267)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using an average for rate rule [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction78',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=C=CCC(=O)C[O](11273)'],
    transitionState = 'TS78',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(30.9691,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction79',
    reactants = ['[CH2]C=[C]C[C](O)[CH][O](12161)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS79',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction80',
    reactants = ['[CH2][C]=CC[C](O)[CH][O](12162)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS80',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction81',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C=CC[C]1[CH]OO1(12163)'],
    transitionState = 'TS81',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(457.838,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra_O] for rate rule [R4_linear;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction82',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C1[CH]CC(=O)C1[O](12164)'],
    transitionState = 'TS82',
    kinetics = Arrhenius(A=(8.5398e+06,'s^-1'), n=1.09629, Ea=(171.155,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 166.2 to 171.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction83',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C3H3(5450)', 'C=C(O)[CH][O](12165)'],
    transitionState = 'TS83',
    kinetics = Arrhenius(A=(3.59298e+09,'s^-1'), n=0.614601, Ea=(265.829,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.07280803979615884, var=3.9519504359180297, Tref=1000.0, N=5, correlation='Root_1R!H->C_2R!H->C_5R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_2R!H->C_5R!H->O
    Total Standard Deviation in ln(k): 4.168250509264941
Exact match found for rate rule [Root_1R!H->C_2R!H->C_5R!H->O]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction84',
    reactants = ['[CH2]C=C[CH]C(=O)C[O](12166)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS84',
    kinetics = Arrhenius(A=(1.38434e+09,'s^-1'), n=1.16972, Ea=(166.189,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_H/Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction85',
    reactants = ['[CH2]C=[C]CC(=O)C[O](12167)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS85',
    kinetics = Arrhenius(A=(2.572e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction86',
    reactants = ['[CH2][C]=CCC(=O)C[O](12168)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS86',
    kinetics = Arrhenius(A=(760143,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;XH_out]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction87',
    reactants = ['[CH2]C=COC(=C)[CH][O](12169)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS87',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction88',
    reactants = ['[CH2]C=CC([CH2])([O])C=O(12170)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS88',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction89',
    reactants = ['HCO(1372)', '[CH2][CH][CH]C[C]=O(2615)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS89',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction90',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[O]C1(C=O)CC=CC1(12171)'],
    transitionState = 'TS90',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSDS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction91',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['CC=CC=C([O])C=O(12172)'],
    transitionState = 'TS91',
    kinetics = Arrhenius(A=(5.18e+08,'s^-1'), n=0.311, Ea=(17.782,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_De] for rate rule [R4radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction92',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=C=CCC([O])C=O(11271)'],
    transitionState = 'TS92',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction93',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['[CH2]C1[CH]CC1([O])C=O(12066)'],
    transitionState = 'TS93',
    kinetics = Arrhenius(A=(1.66065e+07,'s^-1'), n=1.25778, Ea=(269.226,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra_pri;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 267.6 to 269.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction94',
    reactants = ['C=C([O])C=O(2859)', '[CH]C=C(3739)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS94',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction95',
    reactants = ['[CH2]C=C[CH]C([O])C=O(12173)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS95',
    kinetics = Arrhenius(A=(2.00613e+10,'s^-1'), n=0.845153, Ea=(195.403,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction96',
    reactants = ['[CH2]C=[C]CC([O])C=O(12174)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS96',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction97',
    reactants = ['[CH2][C]=CCC([O])C=O(12175)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS97',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction98',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CC1CC1([O])C=O(12063)'],
    transitionState = 'TS98',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(9.57689,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination
Ea raised from 7.9 to 9.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction99',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CCC=C([O])C=O(12176)'],
    transitionState = 'TS99',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction100',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CC=CC([O])C=O(12177)'],
    transitionState = 'TS100',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_De;XH_Rrad_De] + [R2radExo;Y_rad_De;XH_Rrad] for rate rule [R2radExo;Y_rad_De;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction101',
    reactants = ['[CH]=C[CH]CC([O])C=O(12178)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS101',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_2;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction102',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CC1CC(=O)C1[O](11265)'],
    transitionState = 'TS102',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction103',
    reactants = ['C=C[CH]CC([O])=C[O](11260)'],
    products = ['C=CC=CC(=O)C[O](12179)'],
    transitionState = 'TS103',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction104',
    reactants = ['[CH]=C[CH]CC(=O)C[O](12180)'],
    products = ['C=C[CH]CC([O])=C[O](11260)'],
    transitionState = 'TS104',
    kinetics = Arrhenius(A=(3.73886e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;XH_out]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2929',
    isomers = [
        'C=C[CH]CC([O])=C[O](11260)',
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
    label = 'PDepNetwork #2929',
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


species(
    label = 'C=[C][CH]C(C)[CH]C(15805)',
    structure = SMILES('[CH2][C]=CC(C)[CH]C'),
    E0 = (475.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,2393.66,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00280679,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.685324,'amu*angstrom^2'), symmetry=1, barrier=(29.2134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2605,'amu*angstrom^2'), symmetry=1, barrier=(11.0694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.686694,'amu*angstrom^2'), symmetry=1, barrier=(29.2119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.690016,'amu*angstrom^2'), symmetry=1, barrier=(29.2162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627625,0.0661439,-4.10175e-05,1.26595e-08,-1.5865e-12,57353.9,31.5667], Tmin=(100,'K'), Tmax=(1813.24,'K')), NASAPolynomial(coeffs=[16.228,0.0317298,-1.25487e-05,2.19253e-09,-1.4339e-13,51696.4,-52.9802], Tmin=(1813.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_P) + radical(Cds_S)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09172,0.0173333,-8.20209e-06,-2.63358e-09,2.66049e-12,40755.9,8.10965], Tmin=(100,'K'), Tmax=(946.054,'K')), NASAPolynomial(coeffs=[6.98214,0.0072721,-2.37773e-06,3.99152e-10,-2.69331e-14,39733.9,-11.9544], Tmin=(946.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C3H3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'butene2t(396)',
    structure = SMILES('CC=CC'),
    E0 = (-28.2774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.464594,'amu*angstrom^2'), symmetry=1, barrier=(10.6819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.460586,'amu*angstrom^2'), symmetry=1, barrier=(10.5898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52,0.0276811,9.44413e-07,-1.25016e-08,4.67173e-12,-3343.75,13.0195], Tmin=(100,'K'), Tmax=(1157.77,'K')), NASAPolynomial(coeffs=[6.04212,0.0260075,-1.04845e-05,1.90884e-09,-1.30603e-13,-4862.7,-7.52639], Tmin=(1157.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.2774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene2t""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=[C][CH]C[CH]C(15645)',
    structure = SMILES('[CH2][C]=CC[CH]C'),
    E0 = (507.46,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,387.783,1997.36],'cm^-1')),
        HinderedRotor(inertia=(0.0839302,'amu*angstrom^2'), symmetry=1, barrier=(9.00985,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0842009,'amu*angstrom^2'), symmetry=1, barrier=(9.00579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0057216,'amu*angstrom^2'), symmetry=1, barrier=(9.00722,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523224,'amu*angstrom^2'), symmetry=1, barrier=(56.0911,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3287.25,'J/mol'), sigma=(5.9252,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=513.46 K, Pc=35.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81646,0.0496648,-2.82663e-05,8.17951e-09,-1.00768e-12,61109.6,25.4407], Tmin=(100,'K'), Tmax=(1692.15,'K')), NASAPolynomial(coeffs=[8.71634,0.0333545,-1.38081e-05,2.48333e-09,-1.66121e-13,58774.5,-11.4765], Tmin=(1692.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[CH]C(C)C(17431)',
    structure = SMILES('[CH2][C]=C[CH]C(C)C'),
    E0 = (421.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,180,917.207],'cm^-1')),
        HinderedRotor(inertia=(0.188594,'amu*angstrom^2'), symmetry=1, barrier=(4.33615,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0413583,'amu*angstrom^2'), symmetry=1, barrier=(24.6796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07337,'amu*angstrom^2'), symmetry=1, barrier=(24.679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00221256,'amu*angstrom^2'), symmetry=1, barrier=(24.6759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.87987,'amu*angstrom^2'), symmetry=1, barrier=(89.2059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690704,0.0620563,-2.22794e-05,-9.74684e-09,6.63946e-12,50799.4,26.9301], Tmin=(100,'K'), Tmax=(1084.36,'K')), NASAPolynomial(coeffs=[13.1637,0.0362506,-1.45319e-05,2.6735e-09,-1.85752e-13,46906.5,-39.7323], Tmin=(1084.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'CHCH3(T)(21)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438698,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82365,-0.000909765,3.21389e-05,-3.73492e-08,1.33096e-11,41371.4,7.10941], Tmin=(100,'K'), Tmax=(960.802,'K')), NASAPolynomial(coeffs=[4.3048,0.00943081,-3.27566e-06,5.95138e-10,-4.27321e-14,40709.2,1.84242], Tmin=(960.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][C]=C[CH]C(16921)',
    structure = SMILES('[CH2][C]=C[CH]C'),
    E0 = (477.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,300.801],'cm^-1')),
        HinderedRotor(inertia=(0.00186762,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00186424,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706267,'amu*angstrom^2'), symmetry=1, barrier=(45.3693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1366,0.033084,1.9467e-06,-2.18723e-08,9.53324e-12,57552.8,18.2216], Tmin=(100,'K'), Tmax=(1034.94,'K')), NASAPolynomial(coeffs=[9.21568,0.0232651,-9.24634e-06,1.71508e-09,-1.20643e-13,55148,-20.7126], Tmin=(1034.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = '[CH]C(C)C=[C][CH2](16947)',
    structure = SMILES('[CH]C(C)C=[C][CH2]'),
    E0 = (753.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,931.535,931.551,931.635,3442.03],'cm^-1')),
        HinderedRotor(inertia=(0.155659,'amu*angstrom^2'), symmetry=1, barrier=(3.57891,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.704,'amu*angstrom^2'), symmetry=1, barrier=(85.1623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910342,'amu*angstrom^2'), symmetry=1, barrier=(20.9306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.910354,'amu*angstrom^2'), symmetry=1, barrier=(20.9308,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17185,0.054294,-3.21124e-05,4.08278e-09,2.04505e-12,90704.4,25.7945], Tmin=(100,'K'), Tmax=(1071.94,'K')), NASAPolynomial(coeffs=[12.6149,0.0249374,-9.70492e-06,1.75965e-09,-1.2147e-13,87484.5,-33.7823], Tmin=(1071.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(753.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][C]=CC(C)[C]C(17432)',
    structure = SMILES('[CH2][C]=CC(C)[C]C'),
    E0 = (729.476,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.308683,0.0714616,-5.34818e-05,2.07239e-08,-3.25117e-12,87876.7,30.0503], Tmin=(100,'K'), Tmax=(1503.27,'K')), NASAPolynomial(coeffs=[16.2489,0.0290466,-1.1159e-05,1.95462e-09,-1.29759e-13,83084.2,-53.3501], Tmin=(1503.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.476,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CCJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]=CC(C)[CH]C(17433)',
    structure = SMILES('[CH][C]=CC(C)[CH]C'),
    E0 = (694.989,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00947,0.064408,-3.59515e-05,9.66791e-09,-1.05966e-12,83695.6,30.9972], Tmin=(100,'K'), Tmax=(1974.57,'K')), NASAPolynomial(coeffs=[14.4632,0.0371541,-1.52478e-05,2.67776e-09,-1.74638e-13,78382.6,-43.0622], Tmin=(1974.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(694.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1=CC(C)C1C(17434)',
    structure = SMILES('C=C1[CH]C(C)C1C'),
    E0 = (192.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53857,0.0350802,5.84913e-05,-9.41135e-08,3.64127e-11,23264.3,21.2651], Tmin=(100,'K'), Tmax=(988.868,'K')), NASAPolynomial(coeffs=[13.0377,0.0355852,-1.35977e-05,2.57075e-09,-1.86791e-13,18691.1,-45.707], Tmin=(988.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(192.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C=CC(C)=CC(17435)',
    structure = SMILES('[CH2]C=CC(C)=CC'),
    E0 = (102.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.755848,0.0599969,-1.38677e-05,-2.11035e-08,1.15495e-11,12434.8,25.4102], Tmin=(100,'K'), Tmax=(1011.64,'K')), NASAPolynomial(coeffs=[13.2525,0.0350664,-1.32012e-05,2.37832e-09,-1.64861e-13,8653.66,-41.2151], Tmin=(1011.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2]C=CC(C)C=C(5249)',
    structure = SMILES('[CH2]C=CC(C)C=C'),
    E0 = (170.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.51279,0.067028,-3.84304e-05,7.5765e-09,5.44575e-13,20658.4,27.0237], Tmin=(100,'K'), Tmax=(1202.28,'K')), NASAPolynomial(coeffs=[13.7648,0.035208,-1.4039e-05,2.53977e-09,-1.73151e-13,16585.1,-43.0391], Tmin=(1202.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=CC(C)=CC(17436)',
    structure = SMILES('C=[C]C=C(C)[CH]C'),
    E0 = (336.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,332.276,333.152],'cm^-1')),
        HinderedRotor(inertia=(0.30156,'amu*angstrom^2'), symmetry=1, barrier=(23.5514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166088,'amu*angstrom^2'), symmetry=1, barrier=(13.345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.97897,'amu*angstrom^2'), symmetry=1, barrier=(76.17,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170318,'amu*angstrom^2'), symmetry=1, barrier=(13.2316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.561249,0.0668892,-4.03968e-05,6.81405e-09,1.87134e-12,40579.2,24.3036], Tmin=(100,'K'), Tmax=(1055.3,'K')), NASAPolynomial(coeffs=[13.5457,0.0323483,-1.216e-05,2.15372e-09,-1.46456e-13,37021.6,-42.9095], Tmin=(1055.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_S) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C][CH]C(C)C=C(15837)',
    structure = SMILES('[CH2][C]=CC(C)C=C'),
    E0 = (408.495,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,325.958,325.982],'cm^-1')),
        HinderedRotor(inertia=(0.118536,'amu*angstrom^2'), symmetry=1, barrier=(8.93746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185329,'amu*angstrom^2'), symmetry=1, barrier=(13.9745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118521,'amu*angstrom^2'), symmetry=1, barrier=(8.9374,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.356527,'amu*angstrom^2'), symmetry=1, barrier=(26.8825,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.871336,0.0671985,-4.85097e-05,1.87227e-08,-3.04534e-12,49244.3,26.1973], Tmin=(100,'K'), Tmax=(1396.69,'K')), NASAPolynomial(coeffs=[11.7422,0.0360652,-1.50735e-05,2.76297e-09,-1.88634e-13,46207.6,-29.8805], Tmin=(1396.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C#CC(C)[CH]C(17437)',
    structure = SMILES('[CH2]C#CC(C)[CH]C'),
    E0 = (401.664,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2100,2250,500,550,211.787,2159.06],'cm^-1')),
        HinderedRotor(inertia=(0.563512,'amu*angstrom^2'), symmetry=1, barrier=(17.9309,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(6.88674e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00370508,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.15936,'amu*angstrom^2'), symmetry=1, barrier=(102.135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.20112,'amu*angstrom^2'), symmetry=1, barrier=(102.123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.873726,0.0605013,-3.66543e-05,1.14582e-08,-1.47609e-12,48428.4,30.4047], Tmin=(100,'K'), Tmax=(1757.09,'K')), NASAPolynomial(coeffs=[13.7292,0.0312358,-1.16707e-05,1.97904e-09,-1.27383e-13,43910.8,-38.8616], Tmin=(1757.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.664,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(Cs_S) + radical(Propargyl)"""),
)

species(
    label = '[CH2][C]=CC=CC(16966)',
    structure = SMILES('C=[C]C=C[CH]C'),
    E0 = (375.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,284.968,287.827],'cm^-1')),
        HinderedRotor(inertia=(0.366974,'amu*angstrom^2'), symmetry=1, barrier=(21.2347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.27624,'amu*angstrom^2'), symmetry=1, barrier=(75.2771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361561,'amu*angstrom^2'), symmetry=1, barrier=(21.2466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46009,0.0462592,-7.78171e-06,-2.26947e-08,1.21937e-11,45245.2,20.3196], Tmin=(100,'K'), Tmax=(972.572,'K')), NASAPolynomial(coeffs=[12.1933,0.0244663,-8.64205e-06,1.52419e-09,-1.05628e-13,42100.4,-36.598], Tmin=(972.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(375.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(Allyl_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C][CH2](16918)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,260.76,263.05,263.653],'cm^-1')),
        HinderedRotor(inertia=(1.04394,'amu*angstrom^2'), symmetry=1, barrier=(50.8215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2703,0.0129731,8.09839e-06,-1.37802e-08,4.71189e-12,73935.1,11.2848], Tmin=(100,'K'), Tmax=(1104.68,'K')), NASAPolynomial(coeffs=[4.54119,0.0159815,-6.32008e-06,1.15742e-09,-7.99393e-14,73190,2.92522], Tmin=(1104.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C[CH][CH]C(1186)',
    structure = SMILES('C[CH][CH]C'),
    E0 = (244.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2670.04],'cm^-1')),
        HinderedRotor(inertia=(0.341695,'amu*angstrom^2'), symmetry=1, barrier=(7.85623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00155291,'amu*angstrom^2'), symmetry=1, barrier=(7.85531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136081,'amu*angstrom^2'), symmetry=1, barrier=(68.8251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62536,0.022997,-1.74547e-06,-3.25078e-09,7.61424e-13,29450.7,14.8448], Tmin=(100,'K'), Tmax=(2137.09,'K')), NASAPolynomial(coeffs=[11.4859,0.0201914,-8.13358e-06,1.34908e-09,-8.16544e-14,23371.9,-35.4092], Tmin=(2137.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC)"""),
)

species(
    label = '[CH2][C]=C[CH][CH]C(16968)',
    structure = SMILES('[CH2][C]=C[CH][CH]C'),
    E0 = (648.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,378.956,379.128],'cm^-1')),
        HinderedRotor(inertia=(0.0011767,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00117509,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00116876,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00116967,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.1277,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82592,0.0464691,-2.55931e-05,6.86574e-09,-7.52239e-13,78083.9,24.3143], Tmin=(100,'K'), Tmax=(1976.42,'K')), NASAPolynomial(coeffs=[11.4306,0.0270305,-1.08402e-05,1.88942e-09,-1.22777e-13,74287.4,-28.5663], Tmin=(1976.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_S) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C](C)[CH]C(17438)',
    structure = SMILES('[CH2][C]=C[C](C)[CH]C'),
    E0 = (607.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,402.595,2037.53],'cm^-1')),
        HinderedRotor(inertia=(0.00307187,'amu*angstrom^2'), symmetry=1, barrier=(9.04983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.078678,'amu*angstrom^2'), symmetry=1, barrier=(9.04982,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.320393,'amu*angstrom^2'), symmetry=1, barrier=(36.8526,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00749202,'amu*angstrom^2'), symmetry=1, barrier=(22.0716,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.577359,'amu*angstrom^2'), symmetry=1, barrier=(66.4112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.891948,0.0608996,-3.60866e-05,1.06265e-08,-1.27244e-12,73217.6,28.872], Tmin=(100,'K'), Tmax=(1882.68,'K')), NASAPolynomial(coeffs=[15.2286,0.0304393,-1.18177e-05,2.03271e-09,-1.31282e-13,67819.3,-49.365], Tmin=(1882.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_T) + radical(Cs_S) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH]C)[CH][C]=C(15825)',
    structure = SMILES('[CH2][C]=CC([CH2])[CH]C'),
    E0 = (680.886,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,233.894,235.721],'cm^-1')),
        HinderedRotor(inertia=(1.0536e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.627709,'amu*angstrom^2'), symmetry=1, barrier=(24.5288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303204,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.345,'amu*angstrom^2'), symmetry=1, barrier=(92.7764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388177,'amu*angstrom^2'), symmetry=1, barrier=(92.767,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.890035,0.0653111,-4.63234e-05,1.79337e-08,-2.93446e-12,82005.9,32.3814], Tmin=(100,'K'), Tmax=(1393.78,'K')), NASAPolynomial(coeffs=[11.3248,0.0353643,-1.40943e-05,2.51804e-09,-1.69369e-13,79097.1,-21.4249], Tmin=(1393.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]C(C)[CH][C]=C(16020)',
    structure = SMILES('[CH2][C]=CC(C)[CH][CH2]'),
    E0 = (681.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,938.152,939.124],'cm^-1')),
        HinderedRotor(inertia=(0.207219,'amu*angstrom^2'), symmetry=1, barrier=(4.76528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20659,'amu*angstrom^2'), symmetry=1, barrier=(4.75082,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00760191,'amu*angstrom^2'), symmetry=1, barrier=(4.74942,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0579901,'amu*angstrom^2'), symmetry=1, barrier=(1.33356,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.41124,'amu*angstrom^2'), symmetry=1, barrier=(32.461,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.957014,0.0637523,-4.18569e-05,1.4027e-08,-1.94973e-12,82023,32.1476], Tmin=(100,'K'), Tmax=(1613.82,'K')), NASAPolynomial(coeffs=[13.0172,0.03386,-1.40728e-05,2.54942e-09,-1.71711e-13,78130.4,-31.8079], Tmin=(1613.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(681.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C(C)[CH]C(17439)',
    structure = SMILES('[CH2][C]=[C]C(C)[CH]C'),
    E0 = (713.645,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1670,1700,300,440,411.335,1635.91],'cm^-1')),
        HinderedRotor(inertia=(0.0668185,'amu*angstrom^2'), symmetry=1, barrier=(8.03777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000988873,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.066752,'amu*angstrom^2'), symmetry=1, barrier=(8.03057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0670837,'amu*angstrom^2'), symmetry=1, barrier=(8.04469,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205188,'amu*angstrom^2'), symmetry=1, barrier=(24.7209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19862,0.0640619,-4.4365e-05,1.66088e-08,-2.7046e-12,85930.1,29.9609], Tmin=(100,'K'), Tmax=(1346.15,'K')), NASAPolynomial(coeffs=[9.34061,0.0398683,-1.74061e-05,3.25758e-09,-2.25053e-13,83738,-11.7396], Tmin=(1346.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713.645,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C](C)CC(17020)',
    structure = SMILES('[CH2][C]=C[C](C)CC'),
    E0 = (413.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,886.971,3953.56],'cm^-1')),
        HinderedRotor(inertia=(0.90052,'amu*angstrom^2'), symmetry=1, barrier=(20.7047,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900465,'amu*angstrom^2'), symmetry=1, barrier=(20.7035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.309538,'amu*angstrom^2'), symmetry=1, barrier=(7.11689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0624669,'amu*angstrom^2'), symmetry=1, barrier=(89.1808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.87815,'amu*angstrom^2'), symmetry=1, barrier=(89.1662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752176,0.0618213,-2.50454e-05,-5.21631e-09,4.86738e-12,49826.5,26.7967], Tmin=(100,'K'), Tmax=(1093.56,'K')), NASAPolynomial(coeffs=[12.1392,0.0371851,-1.45917e-05,2.63884e-09,-1.81051e-13,46318.6,-33.8096], Tmin=(1093.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_T) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]CC(C)[CH][C]=C(16023)',
    structure = SMILES('[CH2][C]=CC(C)C[CH2]'),
    E0 = (486.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,265.005,1696.49],'cm^-1')),
        HinderedRotor(inertia=(0.171011,'amu*angstrom^2'), symmetry=1, barrier=(8.52234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00240047,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00240047,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.492614,'amu*angstrom^2'), symmetry=1, barrier=(24.5494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.492615,'amu*angstrom^2'), symmetry=1, barrier=(24.5494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.426327,0.0688315,-4.3438e-05,1.21486e-08,-8.57321e-13,58650.1,31.5077], Tmin=(100,'K'), Tmax=(1249.83,'K')), NASAPolynomial(coeffs=[14.3724,0.0339995,-1.3397e-05,2.39908e-09,-1.62212e-13,54398.5,-41.9464], Tmin=(1249.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]C(C)[CH]C(17440)',
    structure = SMILES('[CH2]C=[C]C(C)[CH]C'),
    E0 = (475.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,2393.66,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00280679,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.685324,'amu*angstrom^2'), symmetry=1, barrier=(29.2134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2605,'amu*angstrom^2'), symmetry=1, barrier=(11.0694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.686694,'amu*angstrom^2'), symmetry=1, barrier=(29.2119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.690016,'amu*angstrom^2'), symmetry=1, barrier=(29.2162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.627625,0.0661439,-4.10175e-05,1.26595e-08,-1.5865e-12,57353.9,31.5667], Tmin=(100,'K'), Tmax=(1813.24,'K')), NASAPolynomial(coeffs=[16.228,0.0317298,-1.25487e-05,2.19253e-09,-1.4339e-13,51696.4,-52.9802], Tmin=(1813.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C([CH][C]=C)CC(15663)',
    structure = SMILES('[CH2][C]=CC([CH2])CC'),
    E0 = (486.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.14219,'amu*angstrom^2'), symmetry=1, barrier=(16.7226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.726776,'amu*angstrom^2'), symmetry=1, barrier=(16.71,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151343,'amu*angstrom^2'), symmetry=1, barrier=(3.47968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0814108,'amu*angstrom^2'), symmetry=1, barrier=(87.2086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00835095,'amu*angstrom^2'), symmetry=1, barrier=(87.4454,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3457.08,'J/mol'), sigma=(6.27317,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=539.99 K, Pc=31.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.459773,0.0693092,-4.45292e-05,1.21559e-08,-3.58507e-13,58628.3,31.3731], Tmin=(100,'K'), Tmax=(1114.95,'K')), NASAPolynomial(coeffs=[13.1848,0.0347865,-1.30573e-05,2.29079e-09,-1.53986e-13,55098.9,-34.5047], Tmin=(1114.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C(C)CC(17022)',
    structure = SMILES('[CH2][C]=[C]C(C)CC'),
    E0 = (519.103,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3100,440,815,1455,1000,1670,1700,300,440,189.119,189.143],'cm^-1')),
        HinderedRotor(inertia=(0.00471298,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00471232,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00471428,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.610253,'amu*angstrom^2'), symmetry=1, barrier=(15.4925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901979,'amu*angstrom^2'), symmetry=1, barrier=(22.8953,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.606246,0.0698813,-4.84443e-05,1.76461e-08,-2.67269e-12,62559.4,29.5392], Tmin=(100,'K'), Tmax=(1506.98,'K')), NASAPolynomial(coeffs=[13.4951,0.0356702,-1.43916e-05,2.58157e-09,-1.73564e-13,58674.8,-37.9276], Tmin=(1506.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=C[C](C)[CH]C(17441)',
    structure = SMILES('[CH2]C=C[C](C)[CH]C'),
    E0 = (369.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.993448,0.0557801,-1.07878e-05,-1.75386e-08,8.52e-12,44610.7,28.0143], Tmin=(100,'K'), Tmax=(1086.64,'K')), NASAPolynomial(coeffs=[11.3033,0.0385017,-1.54733e-05,2.84357e-09,-1.97163e-13,41149.6,-28.1971], Tmin=(1086.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_T) + radical(Cs_S) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=[C]C(C)[CH]C(17442)',
    structure = SMILES('C[C]=[C]C(C)[CH]C'),
    E0 = (562.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3025,407.5,1350,352.5,1670,1700,300,440,315.669,4000],'cm^-1')),
        HinderedRotor(inertia=(0.155305,'amu*angstrom^2'), symmetry=1, barrier=(10.9818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00294133,'amu*angstrom^2'), symmetry=1, barrier=(10.9818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.394241,'amu*angstrom^2'), symmetry=1, barrier=(27.8772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155302,'amu*angstrom^2'), symmetry=1, barrier=(10.9818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155303,'amu*angstrom^2'), symmetry=1, barrier=(10.9818,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30086,0.0650662,-4.66185e-05,2.0709e-08,-4.45558e-12,67702.9,29.3361], Tmin=(100,'K'), Tmax=(990.077,'K')), NASAPolynomial(coeffs=[5.29536,0.048928,-2.21685e-05,4.24553e-09,-2.98448e-13,66912,10.1048], Tmin=(990.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=CC([CH2])[CH]C(4147)',
    structure = SMILES('[CH2]C=CC([CH2])[CH]C'),
    E0 = (443.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,930.455],'cm^-1')),
        HinderedRotor(inertia=(0.180361,'amu*angstrom^2'), symmetry=1, barrier=(4.14686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171002,'amu*angstrom^2'), symmetry=1, barrier=(29.0556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00256024,'amu*angstrom^2'), symmetry=1, barrier=(29.0691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00789359,'amu*angstrom^2'), symmetry=1, barrier=(89.6243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.89892,'amu*angstrom^2'), symmetry=1, barrier=(89.6437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.713248,0.0631297,-2.98162e-05,-7.12472e-10,3.50664e-12,53411.9,32.5467], Tmin=(100,'K'), Tmax=(1093.34,'K')), NASAPolynomial(coeffs=[12.2655,0.0362389,-1.40149e-05,2.51304e-09,-1.71526e-13,49967,-28.419], Tmin=(1093.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = 'C[C]=C[C](C)[CH]C(17443)',
    structure = SMILES('C[C]=C[C](C)[CH]C'),
    E0 = (456.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27404,0.0592024,-3.11447e-05,7.73849e-09,-7.68896e-13,54976.2,27.1955], Tmin=(100,'K'), Tmax=(2199.81,'K')), NASAPolynomial(coeffs=[15.9653,0.0324877,-1.29279e-05,2.21754e-09,-1.41437e-13,48512.9,-55.2628], Tmin=(2199.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_T) + radical(Cs_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]C(C)C=C[CH2](4105)',
    structure = SMILES('[CH2][CH]C(C)C=C[CH2]'),
    E0 = (443.208,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1226.02],'cm^-1')),
        HinderedRotor(inertia=(0.182504,'amu*angstrom^2'), symmetry=1, barrier=(4.19614,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0218248,'amu*angstrom^2'), symmetry=1, barrier=(23.2796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00393486,'amu*angstrom^2'), symmetry=1, barrier=(4.19745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182714,'amu*angstrom^2'), symmetry=1, barrier=(4.20096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01233,'amu*angstrom^2'), symmetry=1, barrier=(23.2754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722107,0.0622069,-2.74146e-05,-2.0837e-09,3.45469e-12,53431.8,32.5258], Tmin=(100,'K'), Tmax=(1157.75,'K')), NASAPolynomial(coeffs=[12.8894,0.0363317,-1.4831e-05,2.72869e-09,-1.88342e-13,49531.2,-32.6348], Tmin=(1157.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH]C)C=[C]C(17026)',
    structure = SMILES('[CH2]C([CH]C)C=[C]C'),
    E0 = (529.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,365.193,1535.92],'cm^-1')),
        HinderedRotor(inertia=(0.00126399,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0722335,'amu*angstrom^2'), symmetry=1, barrier=(6.8359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0722293,'amu*angstrom^2'), symmetry=1, barrier=(6.8359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0722309,'amu*angstrom^2'), symmetry=1, barrier=(6.83589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0722326,'amu*angstrom^2'), symmetry=1, barrier=(6.83591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1777,0.0644516,-4.33603e-05,1.67258e-08,-2.90235e-12,63769.9,31.0674], Tmin=(100,'K'), Tmax=(1249.83,'K')), NASAPolynomial(coeffs=[7.69665,0.0435885,-1.83216e-05,3.37026e-09,-2.30919e-13,62140.3,-1.83671], Tmin=(1249.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]C(C)C=[C]C(17444)',
    structure = SMILES('[CH2][CH]C(C)C=[C]C'),
    E0 = (529.551,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1590.95,4000],'cm^-1')),
        HinderedRotor(inertia=(0.184228,'amu*angstrom^2'), symmetry=1, barrier=(9.04932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00503823,'amu*angstrom^2'), symmetry=1, barrier=(9.04945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.47588,'amu*angstrom^2'), symmetry=1, barrier=(23.379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.475948,'amu*angstrom^2'), symmetry=1, barrier=(23.379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184234,'amu*angstrom^2'), symmetry=1, barrier=(9.04939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34225,0.062001,-3.67119e-05,1.09172e-08,-1.37684e-12,63781.7,30.4622], Tmin=(100,'K'), Tmax=(1660.15,'K')), NASAPolynomial(coeffs=[10.1394,0.040805,-1.75607e-05,3.22666e-09,-2.18728e-13,60860.8,-16.4382], Tmin=(1660.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.551,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C(C)[CH][CH]C(17445)',
    structure = SMILES('C=[C]C(C)[CH][CH]C'),
    E0 = (530.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,1685,370,190.713,297.924,3001.49],'cm^-1')),
        HinderedRotor(inertia=(0.065062,'amu*angstrom^2'), symmetry=1, barrier=(1.80892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.065062,'amu*angstrom^2'), symmetry=1, barrier=(1.80892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.065062,'amu*angstrom^2'), symmetry=1, barrier=(1.80892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.065062,'amu*angstrom^2'), symmetry=1, barrier=(1.80892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.065062,'amu*angstrom^2'), symmetry=1, barrier=(1.80892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6938,0.0585628,-3.21805e-05,8.94464e-09,-1.10185e-12,63939.4,30.5816], Tmin=(100,'K'), Tmax=(1535.57,'K')), NASAPolynomial(coeffs=[6.39334,0.046321,-2.02223e-05,3.75303e-09,-2.56627e-13,62496.1,5.89338], Tmin=(1535.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1C(C)C1C(17446)',
    structure = SMILES('C=[C]C1C(C)C1C'),
    E0 = (286.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23017,0.0415452,4.5679e-05,-8.806e-08,3.64491e-11,34577.7,24.5888], Tmin=(100,'K'), Tmax=(965.343,'K')), NASAPolynomial(coeffs=[15.6415,0.0296714,-1.02086e-05,1.87377e-09,-1.36687e-13,29566.2,-55.9749], Tmin=(965.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=C(C)CC(17030)',
    structure = SMILES('C=[C]C=C(C)CC'),
    E0 = (195.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.183353,0.0740923,-5.55858e-05,2.2376e-08,-3.67841e-12,23621.6,26.7753], Tmin=(100,'K'), Tmax=(1437.23,'K')), NASAPolynomial(coeffs=[15.3306,0.0319355,-1.15877e-05,1.96724e-09,-1.28392e-13,19267.6,-51.7956], Tmin=(1437.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]CC(C)=CC(15802)',
    structure = SMILES('C=[C]CC(C)=CC'),
    E0 = (251.511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,302.516,302.516],'cm^-1')),
        HinderedRotor(inertia=(0.142307,'amu*angstrom^2'), symmetry=1, barrier=(9.24166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142307,'amu*angstrom^2'), symmetry=1, barrier=(9.24166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142307,'amu*angstrom^2'), symmetry=1, barrier=(9.24166,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142307,'amu*angstrom^2'), symmetry=1, barrier=(9.24166,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.974961,0.0651918,-4.12655e-05,1.34529e-08,-1.84235e-12,30359.1,26.7122], Tmin=(100,'K'), Tmax=(1606.67,'K')), NASAPolynomial(coeffs=[11.9746,0.0378069,-1.56989e-05,2.84442e-09,-1.91658e-13,26824.5,-31.5704], Tmin=(1606.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC(C)C=C(15528)',
    structure = SMILES('C=[C]CC(C)C=C'),
    E0 = (269.79,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,361.048,361.065,361.131],'cm^-1')),
        HinderedRotor(inertia=(0.14023,'amu*angstrom^2'), symmetry=1, barrier=(13.0002,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140479,'amu*angstrom^2'), symmetry=1, barrier=(12.9993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140643,'amu*angstrom^2'), symmetry=1, barrier=(12.9995,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00129558,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.732002,0.0617375,-2.49394e-05,-6.05057e-09,5.25243e-12,32574.4,29.2485], Tmin=(100,'K'), Tmax=(1096.86,'K')), NASAPolynomial(coeffs=[12.9071,0.0356065,-1.41874e-05,2.59899e-09,-1.79913e-13,28804.5,-35.6246], Tmin=(1096.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.79,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]C(C)[CH]C(17447)',
    structure = SMILES('C#C[CH]C(C)[CH]C'),
    E0 = (413.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2175,525,750,770,3400,2100,419.8,1584.31],'cm^-1')),
        HinderedRotor(inertia=(0.0763865,'amu*angstrom^2'), symmetry=1, barrier=(9.46028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0753501,'amu*angstrom^2'), symmetry=1, barrier=(9.46924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00532636,'amu*angstrom^2'), symmetry=1, barrier=(9.45109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56703,'amu*angstrom^2'), symmetry=1, barrier=(70.4911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.55438,'amu*angstrom^2'), symmetry=1, barrier=(70.538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738047,0.0670464,-5.09877e-05,2.21487e-08,-4.04876e-12,49829.7,28.4761], Tmin=(100,'K'), Tmax=(1278.99,'K')), NASAPolynomial(coeffs=[11.1921,0.0343514,-1.26426e-05,2.16138e-09,-1.41883e-13,47155.5,-24.5312], Tmin=(1278.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(413.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Cs_S) + radical(Sec_Propargyl)"""),
)

species(
    label = 'C=[C]C=CC(16927)',
    structure = SMILES('C=[C]C=CC'),
    E0 = (256.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.08432,'amu*angstrom^2'), symmetry=1, barrier=(24.9307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.08393,'amu*angstrom^2'), symmetry=1, barrier=(24.9218,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92118,0.0394904,-1.38895e-05,-9.15945e-09,6.40789e-12,30977.3,16.8588], Tmin=(100,'K'), Tmax=(973.439,'K')), NASAPolynomial(coeffs=[9.77405,0.021281,-7.49405e-06,1.29737e-09,-8.80426e-14,28782.3,-24.2368], Tmin=(973.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C[C](C)[CH]C(15804)',
    structure = SMILES('C=[C]C[C](C)[CH]C'),
    E0 = (524.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,400,558.192,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0511401,'amu*angstrom^2'), symmetry=1, barrier=(4.70325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0511401,'amu*angstrom^2'), symmetry=1, barrier=(4.70325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0511401,'amu*angstrom^2'), symmetry=1, barrier=(4.70325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0511401,'amu*angstrom^2'), symmetry=1, barrier=(4.70325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0511401,'amu*angstrom^2'), symmetry=1, barrier=(4.70325,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02107,0.0723672,-8.25381e-05,7.28746e-08,-2.82692e-11,63153.3,31.3883], Tmin=(100,'K'), Tmax=(777.749,'K')), NASAPolynomial(coeffs=[2.38579,0.0539607,-2.50761e-05,4.79384e-09,-3.33932e-13,63285.5,27.3615], Tmin=(777.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]C(C)[CH]C(14383)',
    structure = SMILES('[CH]C=CC(C)[CH]C'),
    E0 = (457.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.652417,0.0641573,-2.54866e-05,-1.91898e-09,2.65352e-12,55110.3,31.8255], Tmin=(100,'K'), Tmax=(1253.59,'K')), NASAPolynomial(coeffs=[11.8171,0.0430957,-1.77107e-05,3.21296e-09,-2.18067e-13,51166.8,-29.1253], Tmin=(1253.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH]C)C[C]=C(15527)',
    structure = SMILES('[CH2]C([CH]C)C[C]=C'),
    E0 = (543.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3373.82,'J/mol'), sigma=(6.19899,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=526.98 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03938,0.0656139,-4.41452e-05,1.64907e-08,-2.6901e-12,65523,31.9648], Tmin=(100,'K'), Tmax=(1350.38,'K')), NASAPolynomial(coeffs=[9.28375,0.0411929,-1.70183e-05,3.09842e-09,-2.10733e-13,63296.4,-10.286], Tmin=(1350.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]CC(C)[CH]C(15806)',
    structure = SMILES('[CH]=[C]CC(C)[CH]C'),
    E0 = (585.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,295.693,295.694],'cm^-1')),
        HinderedRotor(inertia=(0.00192832,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192781,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00192795,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129888,'amu*angstrom^2'), symmetry=1, barrier=(8.05892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12989,'amu*angstrom^2'), symmetry=1, barrier=(8.05898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07882,0.0659639,-4.32577e-05,1.48836e-08,-2.19572e-12,70573.3,30.4392], Tmin=(100,'K'), Tmax=(1474.04,'K')), NASAPolynomial(coeffs=[10.3937,0.0406866,-1.75352e-05,3.25002e-09,-2.22633e-13,67827.2,-18.114], Tmin=(1474.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(585.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C(C)C[C]=C(15531)',
    structure = SMILES('[CH2][CH]C(C)C[C]=C'),
    E0 = (544.073,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3373.82,'J/mol'), sigma=(6.19899,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=526.98 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1438,0.0637272,-3.89262e-05,1.19754e-08,-1.54563e-12,65538,31.587], Tmin=(100,'K'), Tmax=(1670.27,'K')), NASAPolynomial(coeffs=[11.5944,0.0387001,-1.64506e-05,3.00465e-09,-2.02928e-13,62046.9,-24.1923], Tmin=(1670.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C][CH]C(C)CC(17032)',
    structure = SMILES('[CH][C]=CC(C)CC'),
    E0 = (500.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3010,987.5,1337.5,450,1655,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.308465,0.0711835,-4.22495e-05,1.25277e-08,-1.51441e-12,60331.1,30.9928], Tmin=(100,'K'), Tmax=(1868.76,'K')), NASAPolynomial(coeffs=[16.9666,0.0355275,-1.36295e-05,2.31773e-09,-1.48543e-13,54105,-59.7892], Tmin=(1868.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(500.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
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
    E0 = (475.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (926.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (603.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (877.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (945.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (941.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (906.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (483.711,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (539.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (500.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (558.915,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (622.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (628.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (536.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (602.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (606.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (785.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (859.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (819.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (892.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (892.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (925.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (605.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (634.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (679.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (633.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (660.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (669.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (724.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (520.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (606.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (508.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (583.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (582.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (730.224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (483.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (498.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (498.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (500.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (640.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (600.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (682.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (670.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (661.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (731.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (675.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (645.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['C3H3(5450)', 'butene2t(396)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', 'C=[C][CH]C[CH]C(15645)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][C]=C[CH]C(C)C(17431)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.318e+10,'s^-1'), n=0.99, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 1 used for cCs(-HC)CJ;CsJ-CdH;CH3
Exact match found for rate rule [cCs(-HC)CJ;CsJ-CdH;CH3]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHCH3(T)(21)', '[CH2][C]=C[CH]C(16921)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3(17)', '[CH]C(C)C=[C][CH2](16947)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2][C]=CC(C)[C]C(17432)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH][C]=CC(C)[CH]C(17433)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['[CH2]C1=CC(C)C1C(17434)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['[CH2]C=CC(C)=CC(17435)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['[CH2]C=CC(C)C=C(5249)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][C]=CC(C)=CC(17436)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(72.3521,'m^3/(mol*s)'), n=1.66655, Ea=(10.8198,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeCs_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', 'C=[C][CH]C(C)C=C(15837)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH2]C#CC(C)[CH]C(17437)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH3(17)', '[CH2][C]=CC=CC(16966)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.0143359,'m^3/(mol*s)'), n=2.40604, Ea=(24.4895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C][CH2](16918)', 'butene2t(396)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00858789,'m^3/(mol*s)'), n=2.41179, Ea=(16.3987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C3H3(5450)', 'C[CH][CH]C(1186)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.04713,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH3(17)', '[CH2][C]=C[CH][CH]C(16968)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C][CH2](16918)', 'C[CH][CH]C(1186)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2][C]=C[C](C)[CH]C(17438)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2]C([CH]C)[CH][C]=C(15825)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2][CH]C(C)[CH][C]=C(16020)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][C]=[C]C(C)[CH]C(17439)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][C]=C[C](C)CC(17020)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.26e+10,'s^-1'), n=0.77, Ea=(192.464,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 164 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]CC(C)[CH][C]=C(16023)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C=[C]C(C)[CH]C(17440)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH][C]=C)CC(15663)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=[C]C(C)CC(17022)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.32e+09,'s^-1'), n=0.99, Ea=(141.419,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['[CH2]C=C[C](C)[CH]C(17441)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[C]=[C]C(C)[CH]C(17442)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['[CH2]C=CC([CH2])[CH]C(4147)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['C[C]=C[C](C)[CH]C(17443)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['[CH2][CH]C(C)C=C[CH2](4105)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([CH]C)C=[C]C(17026)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(91273.5,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]C(C)C=[C]C(17444)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C]C(C)[CH][CH]C(17445)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.37952e+08,'s^-1'), n=1.37167, Ea=(199.228,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['C=[C]C1C(C)C1C(17446)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['C=[C]C=C(C)CC(17030)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['C=[C]CC(C)=CC(15802)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['C=[C]CC(C)C=C(15528)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['H(3)', 'C#C[CH]C(C)[CH]C(17447)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['CHCH3(T)(21)', 'C=[C]C=CC(16927)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=[C]C[C](C)[CH]C(15804)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=[C][CH]C(C)[CH]C(15805)'],
    products = ['[CH]=C[CH]C(C)[CH]C(14383)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([CH]C)C[C]=C(15527)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=[C]CC(C)[CH]C(15806)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]C(C)C[C]=C(15531)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=[C][CH]C(C)CC(17032)'],
    products = ['C=[C][CH]C(C)[CH]C(15805)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R5Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3853',
    isomers = [
        'C=[C][CH]C(C)[CH]C(15805)',
    ],
    reactants = [
        ('C3H3(5450)', 'butene2t(396)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #3853',
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


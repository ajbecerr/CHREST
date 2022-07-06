species(
    label = '[CH2][CH][CH]CC(C)[CH]C(4502)',
    structure = SMILES('[CH2][CH][CH]CC(C)[CH]C'),
    E0 = (546.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,411.584,1157.99,3528.67,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24245,0.0657809,-2.97104e-05,5.55344e-09,-3.73506e-13,65813.6,34.0867], Tmin=(100,'K'), Tmax=(2892.29,'K')), NASAPolynomial(coeffs=[50.7768,0.00698566,-3.53661e-06,5.1585e-10,-2.41145e-14,34255.6,-257.63], Tmin=(2892.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]CC[CH]C(3964)',
    structure = SMILES('[CH2][CH][CH]CC[CH]C'),
    E0 = (572.952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,2091.51,2338.86,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3350.7,'J/mol'), sigma=(6.3658,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=523.37 K, Pc=29.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36111,0.0497551,-1.77213e-05,1.52485e-09,1.23143e-13,68896.7,28.2554], Tmin=(100,'K'), Tmax=(2624.33,'K')), NASAPolynomial(coeffs=[40.8399,0.00972263,-4.60953e-06,6.75852e-10,-3.24062e-14,43339.4,-199.933], Tmin=(2624.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH]C(C)C(9268)',
    structure = SMILES('[CH2][CH][CH]C[CH]C(C)C'),
    E0 = (543.574,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24245,0.0657809,-2.97104e-05,5.55344e-09,-3.73506e-13,65411,33.3935], Tmin=(100,'K'), Tmax=(2892.29,'K')), NASAPolynomial(coeffs=[50.7768,0.00698566,-3.53661e-06,5.1585e-10,-2.41145e-14,33853,-258.323], Tmin=(2892.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.574,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)[CH]C(4479)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)[CH]C'),
    E0 = (555.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,230.747,1116.42,3259.86],'cm^-1')),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497458,0.0790439,-4.94967e-05,1.65362e-08,-2.40218e-12,66910,39.2611], Tmin=(100,'K'), Tmax=(1473.06,'K')), NASAPolynomial(coeffs=[10.5374,0.0517809,-2.17351e-05,3.97208e-09,-2.69863e-13,63952.1,-13.0649], Tmin=(1473.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C[CH]C(556)',
    structure = SMILES('[CH2][CH][CH]C[CH]C'),
    E0 = (596.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,3630.26,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79918,0.0586097,-7.80157e-05,8.14666e-08,-3.30989e-11,71839.7,31.8124], Tmin=(100,'K'), Tmax=(852.594,'K')), NASAPolynomial(coeffs=[-2.48283,0.0533526,-2.41736e-05,4.49758e-09,-3.0559e-13,73491.1,57.1903], Tmin=(852.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH][CH2](5537)',
    structure = SMILES('[CH][CH][CH2]'),
    E0 = (727.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1609,1609.01,1609.02],'cm^-1')),
        HinderedRotor(inertia=(0.0337841,'amu*angstrom^2'), symmetry=1, barrier=(5.78124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00314687,'amu*angstrom^2'), symmetry=1, barrier=(5.78111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42424,0.0149081,-5.18362e-06,2.15834e-10,1.32851e-13,87550.9,15.303], Tmin=(100,'K'), Tmax=(1973.84,'K')), NASAPolynomial(coeffs=[8.10291,0.00913487,-3.61425e-06,6.37549e-10,-4.11102e-14,84981.5,-12.2801], Tmin=(1973.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C(C)[CH]C(114)',
    structure = SMILES('[CH2]C(C)[CH]C'),
    E0 = (225.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1290.57],'cm^-1')),
        HinderedRotor(inertia=(0.100358,'amu*angstrom^2'), symmetry=1, barrier=(2.30744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10319,'amu*angstrom^2'), symmetry=1, barrier=(2.37253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0987185,'amu*angstrom^2'), symmetry=1, barrier=(2.26973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0971084,'amu*angstrom^2'), symmetry=1, barrier=(2.23271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80047,0.0421184,-1.07295e-05,-7.76584e-09,4.08549e-12,27247.4,23.8047], Tmin=(100,'K'), Tmax=(1119,'K')), NASAPolynomial(coeffs=[7.75174,0.0320963,-1.23774e-05,2.20164e-09,-1.48914e-13,25211,-8.72346], Tmin=(1119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl)"""),
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
    label = '[CH]C(C)C[CH][CH][CH2](9044)',
    structure = SMILES('[CH]C(C)C[CH][CH][CH2]'),
    E0 = (821.027,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,195.277,925.785,1038.05,2107.77,2403.71,2755.61],'cm^-1')),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102419,'amu*angstrom^2'), symmetry=1, barrier=(2.39059,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.92656,0.073534,-8.10432e-05,6.70608e-08,-2.44693e-11,98851.7,35.1812], Tmin=(100,'K'), Tmax=(795.536,'K')), NASAPolynomial(coeffs=[3.46626,0.051787,-2.31118e-05,4.32889e-09,-2.97939e-13,98731.7,25.295], Tmin=(795.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(821.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][CH][CH]CC(C)[C]C(9269)',
    structure = SMILES('[CH2][CH][CH]CC(C)[C]C'),
    E0 = (800.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279622,0.0883078,-9.48441e-05,7.51878e-08,-2.67615e-11,96417,38.6516], Tmin=(100,'K'), Tmax=(774.019,'K')), NASAPolynomial(coeffs=[4.35469,0.0602661,-2.69695e-05,5.07232e-09,-3.50552e-13,95995.3,21.3867], Tmin=(774.019,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH][C]CC(C)[CH]C(9270)',
    structure = SMILES('[CH2][CH][C]CC(C)[CH]C'),
    E0 = (800.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.523546,0.0827138,-6.72086e-05,3.52667e-08,-8.72151e-12,96419.1,37.787], Tmin=(100,'K'), Tmax=(904.151,'K')), NASAPolynomial(coeffs=[6.11708,0.0579668,-2.61515e-05,4.99259e-09,-3.50332e-13,95407.7,11.3652], Tmin=(904.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]CC(C)[CH]C(9271)',
    structure = SMILES('[CH2][C][CH]CC(C)[CH]C'),
    E0 = (800.69,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982242,0.0761618,-3.44234e-05,-3.59119e-08,4.36442e-11,96399.9,36.0794], Tmin=(100,'K'), Tmax=(517.32,'K')), NASAPolynomial(coeffs=[4.72946,0.0605256,-2.77596e-05,5.33963e-09,-3.76375e-13,95833.7,18.746], Tmin=(517.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.69,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH][CH]CC(C)[CH]C(9272)',
    structure = SMILES('[CH][CH][CH]CC(C)[CH]C'),
    E0 = (789.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,186.243,558.205,622.409,1386.62,2153.79,2761.9,3513.81],'cm^-1')),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.410654,0.0866821,-9.72662e-05,8.36576e-08,-3.16413e-11,95123.6,40.046], Tmin=(100,'K'), Tmax=(788.149,'K')), NASAPolynomial(coeffs=[2.49917,0.0634269,-2.8921e-05,5.47343e-09,-3.78942e-13,95187.4,32.961], Tmin=(788.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(789.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C1CC(C)C1C(5502)',
    structure = SMILES('[CH2][CH]C1CC(C)C1C'),
    E0 = (291.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.73619,0.0513141,4.64833e-05,-8.84673e-08,3.49519e-11,35184.8,31.9828], Tmin=(100,'K'), Tmax=(1001.7,'K')), NASAPolynomial(coeffs=[14.815,0.0445165,-1.73441e-05,3.26548e-09,-2.34773e-13,29884.7,-48.3398], Tmin=(1001.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'C=CCC[C](C)[CH]C(5514)',
    structure = SMILES('C=CCC[C](C)[CH]C'),
    E0 = (265.975,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.827904,0.0771889,-5.07585e-05,2.11708e-08,-4.60904e-12,32097.3,33.8798], Tmin=(100,'K'), Tmax=(922.2,'K')), NASAPolynomial(coeffs=[3.87437,0.0639754,-2.92669e-05,5.63484e-09,-3.97515e-13,31535.4,19.429], Tmin=(922.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.975,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Cs_S)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)CC(5123)',
    structure = SMILES('[CH2]C=C[CH]C(C)CC'),
    E0 = (163.029,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00711407,0.0730936,-1.44657e-05,-2.66661e-08,1.35279e-11,19764.3,31.8469], Tmin=(100,'K'), Tmax=(1050.14,'K')), NASAPolynomial(coeffs=[15.5054,0.0447532,-1.78256e-05,3.29824e-09,-2.3113e-13,14816.8,-51.7394], Tmin=(1050.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C(C)CCC=C(5516)',
    structure = SMILES('[CH2][CH]C(C)CCC=C'),
    E0 = (285.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0733755,0.0804334,-4.83593e-05,1.42859e-08,-1.70654e-12,34527.7,37.7712], Tmin=(100,'K'), Tmax=(1892.87,'K')), NASAPolynomial(coeffs=[19.5851,0.0388912,-1.54391e-05,2.69141e-09,-1.75194e-13,27085.5,-69.6134], Tmin=(1892.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC[C](C)[CH]C(5504)',
    structure = SMILES('[CH2]C=CC[C](C)[CH]C'),
    E0 = (401.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07474,0.0731911,-4.32618e-05,1.3263e-08,-1.82479e-12,48432.1,32.7874], Tmin=(100,'K'), Tmax=(1428.69,'K')), NASAPolynomial(coeffs=[7.32823,0.0556828,-2.48797e-05,4.68537e-09,-3.23843e-13,46645.2,0.386898], Tmin=(1428.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)[CH]C(5505)',
    structure = SMILES('[CH2]C=C[CH]C(C)[CH]C'),
    E0 = (357.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.160595,0.0722031,-2.64464e-05,-8.57675e-09,6.06686e-12,43154.3,33.8592], Tmin=(100,'K'), Tmax=(1139.55,'K')), NASAPolynomial(coeffs=[14.3526,0.0441652,-1.82066e-05,3.37353e-09,-2.341e-13,38505.8,-42.6673], Tmin=(1139.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C(C)C[CH]C=C(4544)',
    structure = SMILES('[CH2][CH]C(C)CC=C[CH2]'),
    E0 = (421.705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3599.97,'J/mol'), sigma=(6.5858,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=562.31 K, Pc=28.6 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.194485,0.0800216,-5.04997e-05,1.58116e-08,-1.99676e-12,50880.7,38.0572], Tmin=(100,'K'), Tmax=(1820.99,'K')), NASAPolynomial(coeffs=[20.0416,0.0355715,-1.38855e-05,2.40727e-09,-1.56537e-13,43510.7,-71.6998], Tmin=(1820.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH][CH2](5531)',
    structure = SMILES('[CH2][CH][CH][CH2]'),
    E0 = (655.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1802.64],'cm^-1')),
        HinderedRotor(inertia=(0.00215831,'amu*angstrom^2'), symmetry=1, barrier=(4.96293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00215621,'amu*angstrom^2'), symmetry=1, barrier=(4.95965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00214837,'amu*angstrom^2'), symmetry=1, barrier=(4.95028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82147,0.0344596,-5.62911e-05,6.39381e-08,-2.62695e-11,78855.1,21.3291], Tmin=(100,'K'), Tmax=(865.068,'K')), NASAPolynomial(coeffs=[-1.36886,0.0321046,-1.45274e-05,2.71435e-09,-1.84229e-13,80393.2,45.6373], Tmin=(865.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC[CH][CH]C(4158)',
    structure = SMILES('[CH2]C=CC[CH][CH]C'),
    E0 = (440.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,295.384,982.151,2188.21],'cm^-1')),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0841928,'amu*angstrom^2'), symmetry=1, barrier=(2.96904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90721,0.0555768,-2.65323e-05,5.65859e-09,-4.66342e-13,53019.7,29.4014], Tmin=(100,'K'), Tmax=(2690.43,'K')), NASAPolynomial(coeffs=[20.1883,0.0283967,-1.13781e-05,1.90342e-09,-1.17395e-13,43183.1,-76.8867], Tmin=(2690.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(440.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]C[CH][CH]C(4161)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH]C'),
    E0 = (767.411,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,3248.34,3721.5,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0317161,'amu*angstrom^2'), symmetry=1, barrier=(1.72667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36377,0.0743961,-0.000117792,1.28512e-07,-5.18955e-11,92377.3,38.1829], Tmin=(100,'K'), Tmax=(865.662,'K')), NASAPolynomial(coeffs=[-5.72141,0.0655337,-3.03502e-05,5.65576e-09,-3.82639e-13,95162.7,80.3458], Tmin=(865.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(767.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C](C)[CH]C(9273)',
    structure = SMILES('[CH2][CH][CH]C[C](C)[CH]C'),
    E0 = (732.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,1830.53,2847.6,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0186759,'amu*angstrom^2'), symmetry=1, barrier=(9.9017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186759,'amu*angstrom^2'), symmetry=1, barrier=(9.9017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186759,'amu*angstrom^2'), symmetry=1, barrier=(9.9017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186759,'amu*angstrom^2'), symmetry=1, barrier=(9.9017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186759,'amu*angstrom^2'), symmetry=1, barrier=(9.9017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186759,'amu*angstrom^2'), symmetry=1, barrier=(9.9017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186759,'amu*angstrom^2'), symmetry=1, barrier=(9.9017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08413,0.0545227,-1.89188e-05,1.05858e-09,2.41156e-13,88019.5,28.4524], Tmin=(100,'K'), Tmax=(2678.97,'K')), NASAPolynomial(coeffs=[58.1994,-0.00143711,-1.49439e-06,1.83634e-10,-1.85502e-16,50111,-302.586], Tmin=(2678.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C(C)[CH]C(9274)',
    structure = SMILES('[CH2][CH][CH][CH]C(C)[CH]C'),
    E0 = (741.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,2156.61,2460.56,3867.98,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0211892,'amu*angstrom^2'), symmetry=1, barrier=(10.0162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88834,0.0594494,-2.43179e-05,3.32935e-09,-5.99644e-14,89181.6,34.3076], Tmin=(100,'K'), Tmax=(2631.6,'K')), NASAPolynomial(coeffs=[45.7976,0.0102602,-5.41866e-06,8.56632e-10,-4.4986e-14,61046.1,-224.771], Tmin=(2631.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])[CH]C(4576)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])[CH]C'),
    E0 = (752.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,820.044,1553.04,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0286782,'amu*angstrom^2'), symmetry=1, barrier=(2.19903,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449419,0.0868474,-0.000103352,9.16657e-08,-3.41896e-11,90564.5,42.4706], Tmin=(100,'K'), Tmax=(830.571,'K')), NASAPolynomial(coeffs=[1.69735,0.0629525,-2.78981e-05,5.17626e-09,-3.52932e-13,90974.1,40.3954], Tmin=(830.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC(C)[CH][CH2](4548)',
    structure = SMILES('[CH2][CH][CH]CC(C)[CH][CH2]'),
    E0 = (752.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,880.259,1174.39,3820.33,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207106,'amu*angstrom^2'), symmetry=1, barrier=(1.1803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.707,0.0622921,-2.82776e-05,5.28345e-09,-3.58346e-13,90474.2,34.1373], Tmin=(100,'K'), Tmax=(3091.93,'K')), NASAPolynomial(coeffs=[61.417,-0.0058351,9.76926e-07,-2.05723e-10,1.93015e-14,50428.2,-321.425], Tmin=(3091.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C](C)CC(9086)',
    structure = SMILES('[CH2][CH][CH]C[C](C)CC'),
    E0 = (537.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43929,0.060846,-2.42954e-05,3.2722e-09,-7.02366e-14,64651.4,28.2272], Tmin=(100,'K'), Tmax=(2819.53,'K')), NASAPolynomial(coeffs=[62.6546,-0.0041358,1.53583e-07,-1.1544e-10,1.79375e-14,23697.1,-332.244], Tmin=(2819.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C(C)[CH]C(5498)',
    structure = SMILES('[CH2][CH]C[CH]C(C)[CH]C'),
    E0 = (547.005,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,216.889,745.551,2193.82,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0278415,'amu*angstrom^2'), symmetry=1, barrier=(3.2455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0278415,'amu*angstrom^2'), symmetry=1, barrier=(3.2455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0278415,'amu*angstrom^2'), symmetry=1, barrier=(3.2455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0278415,'amu*angstrom^2'), symmetry=1, barrier=(3.2455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0278415,'amu*angstrom^2'), symmetry=1, barrier=(3.2455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0278415,'amu*angstrom^2'), symmetry=1, barrier=(3.2455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0278415,'amu*angstrom^2'), symmetry=1, barrier=(3.2455,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66458,0.069606,-3.43019e-05,7.42019e-09,-6.18196e-13,65852.9,36.1953], Tmin=(100,'K'), Tmax=(2605,'K')), NASAPolynomial(coeffs=[21.6232,0.0389594,-1.66551e-05,2.90404e-09,-1.84784e-13,55454.5,-79.2024], Tmin=(2605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC(C)C[CH2](4555)',
    structure = SMILES('[CH2][CH][CH]CC(C)C[CH2]'),
    E0 = (557.625,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,235.853,275.706,2296.16,3764.3],'cm^-1')),
        HinderedRotor(inertia=(0.0232543,'amu*angstrom^2'), symmetry=1, barrier=(3.27836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232543,'amu*angstrom^2'), symmetry=1, barrier=(3.27836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232543,'amu*angstrom^2'), symmetry=1, barrier=(3.27836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232543,'amu*angstrom^2'), symmetry=1, barrier=(3.27836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232543,'amu*angstrom^2'), symmetry=1, barrier=(3.27836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232543,'amu*angstrom^2'), symmetry=1, barrier=(3.27836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232543,'amu*angstrom^2'), symmetry=1, barrier=(3.27836,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.394182,0.085686,-8.47282e-05,6.64632e-08,-2.43354e-11,67190.6,40.1221], Tmin=(100,'K'), Tmax=(764.276,'K')), NASAPolynomial(coeffs=[3.12349,0.0642761,-2.87235e-05,5.41229e-09,-3.74931e-13,66981.5,29.05], Tmin=(764.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.625,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C(C)CC(9087)',
    structure = SMILES('[CH2][CH][CH][CH]C(C)CC'),
    E0 = (546.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,411.584,1157.99,3528.67,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24245,0.0657809,-2.97104e-05,5.55344e-09,-3.73506e-13,65813.6,34.0867], Tmin=(100,'K'), Tmax=(2892.29,'K')), NASAPolynomial(coeffs=[50.7768,0.00698566,-3.53661e-06,5.1585e-10,-2.41145e-14,34255.6,-257.63], Tmin=(2892.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])CC(4038)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])CC'),
    E0 = (557.461,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,278.872,559.302,2157.52,3748.91],'cm^-1')),
        HinderedRotor(inertia=(0.0281321,'amu*angstrom^2'), symmetry=1, barrier=(3.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281321,'amu*angstrom^2'), symmetry=1, barrier=(3.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281321,'amu*angstrom^2'), symmetry=1, barrier=(3.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281321,'amu*angstrom^2'), symmetry=1, barrier=(3.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281321,'amu*angstrom^2'), symmetry=1, barrier=(3.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281321,'amu*angstrom^2'), symmetry=1, barrier=(3.14,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0281321,'amu*angstrom^2'), symmetry=1, barrier=(3.14,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296787,0.0875608,-8.99057e-05,7.03435e-08,-2.47182e-11,67174.6,40.4657], Tmin=(100,'K'), Tmax=(806.465,'K')), NASAPolynomial(coeffs=[3.66123,0.0622571,-2.68156e-05,4.94189e-09,-3.37175e-13,66912.1,26.695], Tmin=(806.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[C](C)[CH]C(5499)',
    structure = SMILES('[CH2][CH]CC[C](C)[CH]C'),
    E0 = (537.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,225.388,1282.6,3122.05,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0545736,'amu*angstrom^2'), symmetry=1, barrier=(2.78363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0545736,'amu*angstrom^2'), symmetry=1, barrier=(2.78363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0545736,'amu*angstrom^2'), symmetry=1, barrier=(2.78363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0545736,'amu*angstrom^2'), symmetry=1, barrier=(2.78363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0545736,'amu*angstrom^2'), symmetry=1, barrier=(2.78363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0545736,'amu*angstrom^2'), symmetry=1, barrier=(2.78363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0545736,'amu*angstrom^2'), symmetry=1, barrier=(2.78363,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.85974,0.0646874,-2.89247e-05,5.16679e-09,-3.2094e-13,64690.7,30.3421], Tmin=(100,'K'), Tmax=(2990.97,'K')), NASAPolynomial(coeffs=[64.9421,-0.00676162,1.10149e-06,-2.31706e-10,2.21218e-14,22375,-345.843], Tmin=(2990.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH][CH]C(C)[CH]C(9275)',
    structure = SMILES('[CH2]C[CH][CH]C(C)[CH]C'),
    E0 = (547.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,243.914,1225.67,3355.25,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0468401,'amu*angstrom^2'), symmetry=1, barrier=(2.3414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468401,'amu*angstrom^2'), symmetry=1, barrier=(2.3414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468401,'amu*angstrom^2'), symmetry=1, barrier=(2.3414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468401,'amu*angstrom^2'), symmetry=1, barrier=(2.3414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468401,'amu*angstrom^2'), symmetry=1, barrier=(2.3414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468401,'amu*angstrom^2'), symmetry=1, barrier=(2.3414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468401,'amu*angstrom^2'), symmetry=1, barrier=(2.3414,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86629,0.0680775,-3.23496e-05,6.49395e-09,-4.83787e-13,65845.2,35.3072], Tmin=(100,'K'), Tmax=(2897.59,'K')), NASAPolynomial(coeffs=[47.2647,0.0103617,-5.03674e-06,8.0003e-10,-4.34396e-14,37455.9,-235.602], Tmin=(2897.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(547.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CCC([CH2])[CH]C(4579)',
    structure = SMILES('[CH2][CH]CCC([CH2])[CH]C'),
    E0 = (557.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,251.019,465.678,2074.47,3742.99],'cm^-1')),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463966,0.0829174,-6.58677e-05,3.56466e-08,-9.25248e-12,67180,39.8737], Tmin=(100,'K'), Tmax=(863.079,'K')), NASAPolynomial(coeffs=[5.41129,0.0599893,-2.60206e-05,4.8685e-09,-3.37539e-13,66326,16.7341], Tmin=(863.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C[C](C)[CH]C(5500)',
    structure = SMILES('[CH2]C[CH]C[C](C)[CH]C'),
    E0 = (537.897,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,193.604,1717.85,3132.14,3335.31],'cm^-1')),
        HinderedRotor(inertia=(0.0473139,'amu*angstrom^2'), symmetry=1, barrier=(3.42208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0473139,'amu*angstrom^2'), symmetry=1, barrier=(3.42208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0473139,'amu*angstrom^2'), symmetry=1, barrier=(3.42208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0473139,'amu*angstrom^2'), symmetry=1, barrier=(3.42208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0473139,'amu*angstrom^2'), symmetry=1, barrier=(3.42208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0473139,'amu*angstrom^2'), symmetry=1, barrier=(3.42208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0473139,'amu*angstrom^2'), symmetry=1, barrier=(3.42208,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.06341,0.0631406,-2.69309e-05,4.21047e-09,-1.80092e-13,64682.9,29.4466], Tmin=(100,'K'), Tmax=(2806.94,'K')), NASAPolynomial(coeffs=[59.1968,-0.000808383,-1.3306e-06,1.66472e-10,-1.27e-15,26850.1,-310.558], Tmin=(2806.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH]C(C)[CH]C(4501)',
    structure = SMILES('C[CH][CH][CH]C(C)[CH]C'),
    E0 = (536.217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,1698.15,2965.81,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.658945,0.0823981,-8.94173e-05,8.31501e-08,-3.39568e-11,64603.6,40.7044], Tmin=(100,'K'), Tmax=(790.11,'K')), NASAPolynomial(coeffs=[-0.589045,0.0706006,-3.26283e-05,6.21492e-09,-4.31746e-13,65366.3,50.0095], Tmin=(790.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][CH]CCC(C)[CH][CH2](4553)',
    structure = SMILES('[CH2][CH]CCC(C)[CH][CH2]'),
    E0 = (557.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,195.127,780.466,1810.75,3110.43],'cm^-1')),
        HinderedRotor(inertia=(0.0505855,'amu*angstrom^2'), symmetry=1, barrier=(2.21445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0505855,'amu*angstrom^2'), symmetry=1, barrier=(2.21445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0505855,'amu*angstrom^2'), symmetry=1, barrier=(2.21445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0505855,'amu*angstrom^2'), symmetry=1, barrier=(2.21445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0505855,'amu*angstrom^2'), symmetry=1, barrier=(2.21445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0505855,'amu*angstrom^2'), symmetry=1, barrier=(2.21445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0505855,'amu*angstrom^2'), symmetry=1, barrier=(2.21445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.699759,0.0794089,-5.48677e-05,2.37614e-08,-5.0856e-12,67190,39.0332], Tmin=(100,'K'), Tmax=(979.613,'K')), NASAPolynomial(coeffs=[5.04216,0.0616782,-2.77186e-05,5.28576e-09,-3.70658e-13,66339.2,18.1729], Tmin=(979.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]CC([CH2])[CH]C(4058)',
    structure = SMILES('[CH2]C[CH]CC([CH2])[CH]C'),
    E0 = (557.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,243.995,501.37,2127.74,3749.88],'cm^-1')),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34826,0.0694129,5.23395e-06,-1.2171e-07,1.10855e-10,67142.7,36.7387], Tmin=(100,'K'), Tmax=(465.227,'K')), NASAPolynomial(coeffs=[4.20272,0.0622043,-2.7413e-05,5.16142e-09,-3.58889e-13,66689.5,23.1358], Tmin=(465.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C[C](C)[CH]C(4500)',
    structure = SMILES('C[CH][CH]C[C](C)[CH]C'),
    E0 = (527.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1765.94,2581.8,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.4931,0.091778,-0.000127515,1.30183e-07,-5.22604e-11,63507.2,39.8762], Tmin=(100,'K'), Tmax=(837.616,'K')), NASAPolynomial(coeffs=[-3.40655,0.0753008,-3.51508e-05,6.64087e-09,-4.55645e-13,65391.8,65.3488], Tmin=(837.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]C(C)C[CH]C[CH2](4002)',
    structure = SMILES('[CH2][CH]C(C)C[CH]C[CH2]'),
    E0 = (557.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,193.941,246.055,2194.54,3733.77],'cm^-1')),
        HinderedRotor(inertia=(0.0255163,'amu*angstrom^2'), symmetry=1, barrier=(3.13333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255163,'amu*angstrom^2'), symmetry=1, barrier=(3.13333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255163,'amu*angstrom^2'), symmetry=1, barrier=(3.13333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255163,'amu*angstrom^2'), symmetry=1, barrier=(3.13333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255163,'amu*angstrom^2'), symmetry=1, barrier=(3.13333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255163,'amu*angstrom^2'), symmetry=1, barrier=(3.13333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255163,'amu*angstrom^2'), symmetry=1, barrier=(3.13333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44052,0.0675007,1.13951e-05,-1.29614e-07,1.15943e-10,67159,36.4181], Tmin=(100,'K'), Tmax=(451.695,'K')), NASAPolynomial(coeffs=[3.50125,0.0645361,-2.95166e-05,5.68088e-09,-4.00897e-13,66816.9,26.388], Tmin=(451.695,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)C[CH][CH]C(3834)',
    structure = SMILES('[CH2]C([CH]C)C[CH][CH]C'),
    E0 = (546.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180,1264.05,3750.06,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442845,0.0857364,-9.15222e-05,7.75899e-08,-2.89865e-11,65880.5,40.7089], Tmin=(100,'K'), Tmax=(812.499,'K')), NASAPolynomial(coeffs=[1.77095,0.0654815,-2.88057e-05,5.35245e-09,-3.6637e-13,66117.4,37.3635], Tmin=(812.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(C)C[CH][CH]C(3839)',
    structure = SMILES('[CH2][CH]C(C)C[CH][CH]C'),
    E0 = (546.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,411.584,1157.99,3528.67,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0334984,'amu*angstrom^2'), symmetry=1, barrier=(2.46322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24245,0.0657809,-2.97104e-05,5.55344e-09,-3.73506e-13,65813.6,34.0867], Tmin=(100,'K'), Tmax=(2892.29,'K')), NASAPolynomial(coeffs=[50.7768,0.00698566,-3.53661e-06,5.1585e-10,-2.41145e-14,34255.6,-257.63], Tmin=(2892.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    E0 = (546.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (991.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (741.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (712.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (996.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1009.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1012.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1012.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1012.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1012.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1001.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (555.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (610.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (610.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (571.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (622.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (577.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (636.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (643.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (605.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (556.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (900.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (903.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (944.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (953.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (963.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (963.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (701.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (684.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (664.169,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (669.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (704.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (663.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (721.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (633.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (699.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (621.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (631.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (631.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (620.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (631.369,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (599.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (610.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2][CH]C=C(3743)', 'butene2t(396)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2][CH][CH]CC[CH]C(3964)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2][CH][CH]C[CH]C(C)C(9268)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHCH3(T)(21)', '[CH2][CH][CH]C[CH]C(556)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][CH][CH2](5537)', '[CH2]C(C)[CH]C(114)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', '[CH]C(C)C[CH][CH][CH2](9044)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][CH][CH]CC(C)[C]C(9269)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][CH][C]CC(C)[CH]C(9270)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][C][CH]CC(C)[CH]C(9271)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH][CH][CH]CC(C)[CH]C(9272)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2][CH]C1CC(C)C1C(5502)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['C=CCC[C](C)[CH]C(5514)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2]C=C[CH]C(C)CC(5123)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2][CH]C(C)CCC=C(5516)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]C=CC[C](C)[CH]C(5504)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.72e+07,'cm^3/(mol*s)'), n=1.64, Ea=(9.07928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2568 used for Cds-CsCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2]C=C[CH]C(C)[CH]C(5505)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2][CH]C(C)C[CH]C=C(4544)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'butene2t(396)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0171758,'m^3/(mol*s)'), n=2.41179, Ea=(16.3987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH3(17)', '[CH2]C=CC[CH][CH]C(4158)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(10100,'cm^3/(mol*s)'), n=2.41, Ea=(28.4512,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 435 used for Cds-CsH_Cds-CsH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C=C(3743)', 'C[CH][CH]C(1186)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.54197,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'C[CH][CH]C(1186)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.3e+09,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH3(17)', '[CH2][CH][CH]C[CH][CH]C(4161)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][CH][CH]C[C](C)[CH]C(9273)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][CH][CH][CH]C(C)[CH]C(9274)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH2][CH][CH]CC([CH2])[CH]C(4576)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2][CH][CH]C[C](C)CC(9086)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 150 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]C[CH]C(C)[CH]C(5498)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2][CH][CH]CC(C)C[CH2](4555)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.04e-11,'s^-1'), n=6.833, Ea=(117.248,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 33 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH][CH][CH]C(C)CC(9087)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.3913e+06,'s^-1'), n=1.66106, Ea=(123.033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH]CC([CH2])CC(4038)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2][CH]CC[C](C)[CH]C(5499)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C[CH][CH]C(C)[CH]C(9275)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(867345,'s^-1'), n=1.96939, Ea=(174.054,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R3HJ;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH]CCC([CH2])[CH]C(4579)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(927.918,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2]C[CH]C[C](C)[CH]C(5500)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(322699,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_1;Y_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['C[CH][CH][CH]C(C)[CH]C(4501)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(754000,'s^-1'), n=1.63, Ea=(74.8936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]CCC(C)[CH][CH2](4553)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C[CH]CC([CH2])[CH]C(4058)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['C[CH][CH]C[C](C)[CH]C(4500)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]C(C)C[CH]C[CH2](4002)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2778.79,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2][CH]C(C)C[CH][CH]C(3839)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(64.2,'s^-1'), n=2.1, Ea=(63.1784,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2081',
    isomers = [
        '[CH2][CH][CH]CC(C)[CH]C(4502)',
    ],
    reactants = [
        ('[CH2][CH]C=C(3743)', 'butene2t(396)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2081',
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


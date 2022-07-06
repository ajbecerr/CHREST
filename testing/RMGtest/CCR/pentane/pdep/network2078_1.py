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
    label = '[CH2]C(C=C)C(C)[CH]C(3845)',
    structure = SMILES('[CH2]C(C=C)C(C)[CH]C'),
    E0 = (277.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3541.48,'J/mol'), sigma=(6.54527,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.17 K, Pc=28.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0939387,0.0790831,-4.21075e-05,6.26435e-09,1.44898e-12,33551.6,37.5716], Tmin=(100,'K'), Tmax=(1168.13,'K')), NASAPolynomial(coeffs=[14.1853,0.0450581,-1.75117e-05,3.12534e-09,-2.11592e-13,29201,-37.8793], Tmin=(1168.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl)"""),
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
    label = '[CH2][CH]C([CH2])C[CH]C(3910)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH]C'),
    E0 = (581.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,353.512,3195.31,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261169,'amu*angstrom^2'), symmetry=1, barrier=(3.92272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3372.64,'J/mol'), sigma=(6.39812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=526.80 K, Pc=29.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50564,0.0619335,-2.04379e-05,-3.52822e-08,3.5309e-11,70000.3,33.9727], Tmin=(100,'K'), Tmax=(536.66,'K')), NASAPolynomial(coeffs=[4.26049,0.0520362,-2.25021e-05,4.21173e-09,-2.92541e-13,69551.5,20.9698], Tmin=(536.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH]C(C)C(9276)',
    structure = SMILES('[CH2][CH]C([CH2])[CH]C(C)C'),
    E0 = (551.947,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497458,0.0790439,-4.94967e-05,1.65362e-08,-2.40218e-12,66507.4,38.568], Tmin=(100,'K'), Tmax=(1473.06,'K')), NASAPolynomial(coeffs=[10.5374,0.0517809,-2.17351e-05,3.97208e-09,-2.69863e-13,63549.5,-13.758], Tmin=(1473.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24245,0.0657809,-2.97104e-05,5.55344e-09,-3.73506e-13,65813.6,34.0867], Tmin=(100,'K'), Tmax=(2892.29,'K')), NASAPolynomial(coeffs=[50.7768,0.00698566,-3.53661e-06,5.1585e-10,-2.41145e-14,34255.6,-257.63], Tmin=(2892.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C([CH2])[CH]C(3871)',
    structure = SMILES('[CH2][CH]C([CH2])[CH]C'),
    E0 = (605.201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,246.588,3146],'cm^-1')),
        HinderedRotor(inertia=(0.00277254,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.6228e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0110378,'amu*angstrom^2'), symmetry=1, barrier=(28.5633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00283321,'amu*angstrom^2'), symmetry=1, barrier=(7.3317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07418,'amu*angstrom^2'), symmetry=1, barrier=(46.3499,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90503,0.0506182,-2.74878e-05,7.53461e-09,-8.92238e-13,72859.4,29.9885], Tmin=(100,'K'), Tmax=(1669.08,'K')), NASAPolynomial(coeffs=[7.3347,0.0376059,-1.57937e-05,2.86374e-09,-1.92623e-13,71046.9,1.01198], Tmin=(1669.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH2](502)',
    structure = SMILES('[CH][CH2]'),
    E0 = (557.957,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1433.18,1433.5],'cm^-1')),
        HinderedRotor(inertia=(0.00559429,'amu*angstrom^2'), symmetry=1, barrier=(8.15686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.7298,0.00311091,1.59755e-05,-2.26123e-08,9.32732e-12,67119.1,8.34543], Tmin=(100,'K'), Tmax=(875.199,'K')), NASAPolynomial(coeffs=[4.97018,0.00511238,-6.01253e-07,2.8794e-11,-6.01916e-16,66608.3,0.848354], Tmin=(875.199,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.957,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C(C)[CH]C(1188)',
    structure = SMILES('[CH2][CH]C(C)[CH]C'),
    E0 = (400.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1353.12,1353.19],'cm^-1')),
        HinderedRotor(inertia=(0.0941265,'amu*angstrom^2'), symmetry=1, barrier=(4.97426,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0038298,'amu*angstrom^2'), symmetry=1, barrier=(4.97721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0926457,'amu*angstrom^2'), symmetry=1, barrier=(4.97987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00223722,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0933422,'amu*angstrom^2'), symmetry=1, barrier=(4.97155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49855,0.0527562,-2.53265e-05,4.9251e-09,-2.68297e-13,48215.2,29.7231], Tmin=(100,'K'), Tmax=(1919.4,'K')), NASAPolynomial(coeffs=[16.308,0.0286312,-1.17384e-05,2.03443e-09,-1.29997e-13,41289,-54.613], Tmin=(1919.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C(C)[CH]C(1275)',
    structure = SMILES('[CH2][CH][CH]C(C)[CH]C'),
    E0 = (570.785,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,549.259,3239.33,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0234404,'amu*angstrom^2'), symmetry=1, barrier=(3.2004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04461,0.0572562,-2.65999e-05,5.15936e-09,-3.5928e-13,68705.3,32.6477], Tmin=(100,'K'), Tmax=(2641.81,'K')), NASAPolynomial(coeffs=[33.6327,0.0159202,-6.81576e-06,1.09699e-09,-6.28765e-14,49749.9,-154.721], Tmin=(2641.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(570.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH]C(C)C([CH2])[CH][CH2](9061)',
    structure = SMILES('[CH]C(C)C([CH2])[CH][CH2]'),
    E0 = (829.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,212.3,859.655,1193.05,1598.22,1982.5],'cm^-1')),
        HinderedRotor(inertia=(0.106827,'amu*angstrom^2'), symmetry=1, barrier=(3.07267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106827,'amu*angstrom^2'), symmetry=1, barrier=(3.07267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106827,'amu*angstrom^2'), symmetry=1, barrier=(3.07267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106827,'amu*angstrom^2'), symmetry=1, barrier=(3.07267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106827,'amu*angstrom^2'), symmetry=1, barrier=(3.07267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106827,'amu*angstrom^2'), symmetry=1, barrier=(3.07267,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.637124,0.0714859,-5.35689e-05,2.22446e-08,-3.90827e-12,99876.7,34.9751], Tmin=(100,'K'), Tmax=(1309.78,'K')), NASAPolynomial(coeffs=[11.594,0.038024,-1.52471e-05,2.73907e-09,-1.85204e-13,97006.5,-20.8424], Tmin=(1309.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(829.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][CH]C([CH2])C(C)[C]C(9277)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)[C]C'),
    E0 = (808.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0459656,0.0856377,-6.53086e-05,2.77219e-08,-5.005e-12,97439.5,38.2426], Tmin=(100,'K'), Tmax=(1271.69,'K')), NASAPolynomial(coeffs=[12.5189,0.0464055,-1.90334e-05,3.463e-09,-2.36032e-13,94267.1,-24.93], Tmin=(1271.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(808.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]C([CH2])C(C)[CH]C(9278)',
    structure = SMILES('[CH2][C]C([CH2])C(C)[CH]C'),
    E0 = (808.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0459656,0.0856377,-6.53086e-05,2.77219e-08,-5.005e-12,97439.5,38.2426], Tmin=(100,'K'), Tmax=(1271.69,'K')), NASAPolynomial(coeffs=[12.5189,0.0464055,-1.90334e-05,3.463e-09,-2.36032e-13,94267.1,-24.93], Tmin=(1271.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(808.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH][CH2])C(C)[CH]C(9279)',
    structure = SMILES('[CH]C([CH][CH2])C(C)[CH]C'),
    E0 = (798.427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,208.67,1035.85,1086.61,1335.05,1470.95,1841.33],'cm^-1')),
        HinderedRotor(inertia=(0.131715,'amu*angstrom^2'), symmetry=1, barrier=(3.52456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131715,'amu*angstrom^2'), symmetry=1, barrier=(3.52456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131715,'amu*angstrom^2'), symmetry=1, barrier=(3.52456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131715,'amu*angstrom^2'), symmetry=1, barrier=(3.52456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131715,'amu*angstrom^2'), symmetry=1, barrier=(3.52456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131715,'amu*angstrom^2'), symmetry=1, barrier=(3.52456,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131715,'amu*angstrom^2'), symmetry=1, barrier=(3.52456,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541527,0.0790771,-5.18643e-05,1.76882e-08,-2.58618e-12,96149.6,38.3243], Tmin=(100,'K'), Tmax=(1479.77,'K')), NASAPolynomial(coeffs=[11.6196,0.0491316,-2.15095e-05,4.01276e-09,-2.75785e-13,92871,-19.4626], Tmin=(1479.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(798.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C([CH2])C(C)[CH]C(9280)',
    structure = SMILES('[CH][CH]C([CH2])C(C)[CH]C'),
    E0 = (798.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,187.329,892.156,1108,1329.59,1561.27,1968.62],'cm^-1')),
        HinderedRotor(inertia=(0.138627,'amu*angstrom^2'), symmetry=1, barrier=(3.18731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138627,'amu*angstrom^2'), symmetry=1, barrier=(3.18731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138627,'amu*angstrom^2'), symmetry=1, barrier=(3.18731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138627,'amu*angstrom^2'), symmetry=1, barrier=(3.18731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138627,'amu*angstrom^2'), symmetry=1, barrier=(3.18731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138627,'amu*angstrom^2'), symmetry=1, barrier=(3.18731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138627,'amu*angstrom^2'), symmetry=1, barrier=(3.18731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.394134,0.0814068,-5.8342e-05,2.34663e-08,-4.13536e-12,96136.6,38.8604], Tmin=(100,'K'), Tmax=(1262.07,'K')), NASAPolynomial(coeffs=[10.193,0.05035,-2.143e-05,3.96804e-09,-2.72968e-13,93663.2,-10.6942], Tmin=(1262.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(798.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CC1C(C)[CH]C(5484)',
    structure = SMILES('[CH2]C1CC1C(C)[CH]C'),
    E0 = (299.019,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.423696,0.0613781,1.77109e-05,-6.16729e-08,2.67971e-11,36107.7,33.7769], Tmin=(100,'K'), Tmax=(986.083,'K')), NASAPolynomial(coeffs=[15.1517,0.0424551,-1.55989e-05,2.82787e-09,-1.98996e-13,31218.5,-47.1339], Tmin=(986.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.019,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1C([CH2])C(C)C1C(4509)',
    structure = SMILES('[CH2]C1C([CH2])C(C)C1C'),
    E0 = (292.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.613865,0.0519936,5.56286e-05,-1.0914e-07,4.61603e-11,35352.2,29.9745], Tmin=(100,'K'), Tmax=(943.05,'K')), NASAPolynomial(coeffs=[17.2521,0.0378984,-1.17835e-05,2.01993e-09,-1.42709e-13,29702.7,-62.6351], Tmin=(943.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]CC(=C)C(C)[CH]C(9281)',
    structure = SMILES('[CH2]CC(=C)C(C)[CH]C'),
    E0 = (270.218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.272231,0.0823665,-5.14525e-05,1.60621e-08,-2.0321e-12,32663.2,37.704], Tmin=(100,'K'), Tmax=(1808.09,'K')), NASAPolynomial(coeffs=[19.6501,0.0382917,-1.48871e-05,2.57961e-09,-1.6788e-13,25459,-70.2088], Tmin=(1808.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C(C)C(C)[CH]C(4476)',
    structure = SMILES('[CH2]C=C(C)C(C)[CH]C'),
    E0 = (203.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,200,800],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.357794,0.0833835,-5.22562e-05,1.63519e-08,-2.07056e-12,24593.5,35.4026], Tmin=(100,'K'), Tmax=(1811.49,'K')), NASAPolynomial(coeffs=[20.1567,0.0380844,-1.47461e-05,2.54728e-09,-1.65392e-13,17161.2,-75.7566], Tmin=(1811.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(203.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])C(C)=CC(9282)',
    structure = SMILES('[CH2]CC([CH2])C(C)=CC'),
    E0 = (267.378,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0687169,0.0850975,-5.94473e-05,2.27319e-08,-3.67841e-12,32307.9,35.0552], Tmin=(100,'K'), Tmax=(1406.87,'K')), NASAPolynomial(coeffs=[13.485,0.0465617,-1.83605e-05,3.26226e-09,-2.18666e-13,28494.3,-34.9602], Tmin=(1406.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.378,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C(C)=CC(4475)',
    structure = SMILES('[CH2][CH]C(C)C(C)=CC'),
    E0 = (256.838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.37313,0.0791044,-4.76596e-05,1.43938e-08,-1.80969e-12,31020.4,34.2348], Tmin=(100,'K'), Tmax=(1727.26,'K')), NASAPolynomial(coeffs=[14.3823,0.0466617,-1.94854e-05,3.51942e-09,-2.35738e-13,26180.9,-41.0078], Tmin=(1727.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C=C)C(C)CC(5096)',
    structure = SMILES('[CH2]C=C([CH2])C(C)CC'),
    E0 = (160.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.299675,0.0797247,-2.6475e-05,-1.93848e-08,1.2258e-11,19416.9,32.8801], Tmin=(100,'K'), Tmax=(1029.59,'K')), NASAPolynomial(coeffs=[17.3182,0.0419959,-1.62601e-05,2.97836e-09,-2.08323e-13,14161,-60.5359], Tmin=(1029.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])C(C)C=C(5100)',
    structure = SMILES('[CH2]CC([CH2])C(C)C=C'),
    E0 = (288.367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.172878,0.0802696,-3.90594e-05,-1.52907e-09,5.30445e-12,34842.6,37.0781], Tmin=(100,'K'), Tmax=(1053.74,'K')), NASAPolynomial(coeffs=[14.9563,0.0436385,-1.65227e-05,2.94455e-09,-2.01041e-13,30499.4,-42.1825], Tmin=(1053.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C(C)C=C(3843)',
    structure = SMILES('[CH2][CH]C(C)C(C)C=C'),
    E0 = (277.827,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0733375,0.0780012,-3.90628e-05,3.95746e-09,1.82957e-12,33571,37.5101], Tmin=(100,'K'), Tmax=(1224.48,'K')), NASAPolynomial(coeffs=[14.9719,0.0449148,-1.82074e-05,3.31515e-09,-2.26423e-13,28682.4,-43.038], Tmin=(1224.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C=C)[C](C)[CH]C(5487)',
    structure = SMILES('[CH2]C(C=C)[C](C)[CH]C'),
    E0 = (463.085,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,224.721,424.959,3068.65],'cm^-1')),
        HinderedRotor(inertia=(0.0830327,'amu*angstrom^2'), symmetry=1, barrier=(2.40272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830327,'amu*angstrom^2'), symmetry=1, barrier=(2.40272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830327,'amu*angstrom^2'), symmetry=1, barrier=(2.40272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830327,'amu*angstrom^2'), symmetry=1, barrier=(2.40272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830327,'amu*angstrom^2'), symmetry=1, barrier=(2.40272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830327,'amu*angstrom^2'), symmetry=1, barrier=(2.40272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.670003,0.0783893,-6.01822e-05,3.05907e-08,-7.46351e-12,55811.8,35.9753], Tmin=(100,'K'), Tmax=(907.397,'K')), NASAPolynomial(coeffs=[5.50907,0.0570577,-2.49195e-05,4.68323e-09,-3.25662e-13,54933.6,13.0999], Tmin=(907.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.085,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C=C)C(C)[CH]C(5488)',
    structure = SMILES('[CH2]C=C([CH2])C(C)[CH]C'),
    E0 = (354.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.157815,0.078995,-3.91273e-05,-2.71186e-10,4.29666e-12,42807.4,34.9326], Tmin=(100,'K'), Tmax=(1113.88,'K')), NASAPolynomial(coeffs=[16.0444,0.0415901,-1.67369e-05,3.07471e-09,-2.12945e-13,37909,-50.767], Tmin=(1113.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C(C)C([CH2])C=C(4528)',
    structure = SMILES('[CH2][CH]C(C)C([CH2])C=C'),
    E0 = (482.909,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3541.48,'J/mol'), sigma=(6.54527,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.17 K, Pc=28.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0956704,0.0801551,-5.3169e-05,1.85006e-08,-2.65427e-12,58236.1,39.3721], Tmin=(100,'K'), Tmax=(1601.93,'K')), NASAPolynomial(coeffs=[16.2723,0.039284,-1.4898e-05,2.57339e-09,-1.68624e-13,52992.1,-47.3066], Tmin=(1601.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.909,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC(C)[CH]C(1272)',
    structure = SMILES('[CH2]C=CC(C)[CH]C'),
    E0 = (237.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,1036.48],'cm^-1')),
        HinderedRotor(inertia=(0.134069,'amu*angstrom^2'), symmetry=1, barrier=(3.08252,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862074,'amu*angstrom^2'), symmetry=1, barrier=(19.8208,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00404403,'amu*angstrom^2'), symmetry=1, barrier=(3.08227,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134062,'amu*angstrom^2'), symmetry=1, barrier=(3.08235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.862078,'amu*angstrom^2'), symmetry=1, barrier=(19.8209,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.707441,0.0612551,-1.64324e-05,-1.46775e-08,7.88375e-12,28748.1,30.7888], Tmin=(100,'K'), Tmax=(1098.63,'K')), NASAPolynomial(coeffs=[12.6824,0.0392818,-1.59589e-05,2.95305e-09,-2.05538e-13,24811.7,-34.0495], Tmin=(1098.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(237.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=CC(C)[CH]C(1184)',
    structure = SMILES('C=CC(C)[CH]C'),
    E0 = (122.488,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,475.178,477.979],'cm^-1')),
        HinderedRotor(inertia=(0.057236,'amu*angstrom^2'), symmetry=1, barrier=(9.11712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0594084,'amu*angstrom^2'), symmetry=1, barrier=(9.15138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10235,'amu*angstrom^2'), symmetry=1, barrier=(2.35362,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0142242,'amu*angstrom^2'), symmetry=1, barrier=(2.34047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38814,0.0477905,-3.17684e-06,-2.119e-08,9.15098e-12,14834,26.2582], Tmin=(100,'K'), Tmax=(1092.58,'K')), NASAPolynomial(coeffs=[10.146,0.0359246,-1.46147e-05,2.70825e-09,-1.88672e-13,11714.8,-22.2858], Tmin=(1092.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.488,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S)"""),
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
    label = '[CH2]C([CH][CH]C)C=C(5043)',
    structure = SMILES('[CH2]C([CH][CH]C)C=C'),
    E0 = (498.236,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,402.67,1554.62,3088.64],'cm^-1')),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0461965,'amu*angstrom^2'), symmetry=1, barrier=(5.20511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34038,0.0602115,-3.50779e-05,1.1048e-08,-1.54084e-12,60017.7,33.1737], Tmin=(100,'K'), Tmax=(1488.46,'K')), NASAPolynomial(coeffs=[7.87184,0.0426592,-1.73893e-05,3.12547e-09,-2.10162e-13,58073.3,-0.93461], Tmin=(1488.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(498.236,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl)"""),
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
    label = '[CH2][CH]C([CH2])[CH][CH]C(9076)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][CH]C'),
    E0 = (775.867,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1455.87,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0156776,'amu*angstrom^2'), symmetry=1, barrier=(11.3079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10013,0.069552,-7.39774e-05,6.16456e-08,-2.28945e-11,93414.1,37.884], Tmin=(100,'K'), Tmax=(797.472,'K')), NASAPolynomial(coeffs=[2.76063,0.0518444,-2.30295e-05,4.3071e-09,-2.9626e-13,93447.5,32.1186], Tmin=(797.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(775.867,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C](C)[CH]C(9283)',
    structure = SMILES('[CH2][CH]C([CH2])[C](C)[CH]C'),
    E0 = (740.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1015.09,3573.33,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.311855,0.0895574,-0.000107256,9.3894e-08,-3.50254e-11,89212.2,41.0757], Tmin=(100,'K'), Tmax=(806.995,'K')), NASAPolynomial(coeffs=[2.76558,0.0622206,-2.82377e-05,5.31553e-09,-3.66127e-13,89310.3,32.8256], Tmin=(806.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(740.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C(C)[CH]C(9284)',
    structure = SMILES('[CH2][CH][C]([CH2])C(C)[CH]C'),
    E0 = (740.716,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1015.09,3573.33,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0372434,'amu*angstrom^2'), symmetry=1, barrier=(3.11718,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.311855,0.0895574,-0.000107256,9.3894e-08,-3.50254e-11,89212.2,41.0757], Tmin=(100,'K'), Tmax=(806.995,'K')), NASAPolynomial(coeffs=[2.76558,0.0622206,-2.82377e-05,5.31553e-09,-3.66127e-13,89310.3,32.8256], Tmin=(806.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(740.716,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])[CH]C(4531)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])[CH]C'),
    E0 = (760.376,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,347.706,3366.17,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0311786,'amu*angstrom^2'), symmetry=1, barrier=(3.96627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.418259,0.0816099,-6.41805e-05,3.1146e-08,-6.78753e-12,91578.7,41.3492], Tmin=(100,'K'), Tmax=(1042.54,'K')), NASAPolynomial(coeffs=[8.11405,0.0520824,-2.16959e-05,3.97833e-09,-2.72654e-13,89974.1,3.90099], Tmin=(1042.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(760.376,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)[CH][CH2](4513)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)[CH][CH2]'),
    E0 = (760.54,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,359.878,3268.61,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0282547,'amu*angstrom^2'), symmetry=1, barrier=(3.87933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0282547,'amu*angstrom^2'), symmetry=1, barrier=(3.87933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0282547,'amu*angstrom^2'), symmetry=1, barrier=(3.87933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0282547,'amu*angstrom^2'), symmetry=1, barrier=(3.87933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0282547,'amu*angstrom^2'), symmetry=1, barrier=(3.87933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0282547,'amu*angstrom^2'), symmetry=1, barrier=(3.87933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0282547,'amu*angstrom^2'), symmetry=1, barrier=(3.87933,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612925,0.0786769,-5.55307e-05,2.26583e-08,-4.17915e-12,91590.1,40.6489], Tmin=(100,'K'), Tmax=(1183.11,'K')), NASAPolynomial(coeffs=[8.22032,0.0529569,-2.29219e-05,4.28369e-09,-2.96467e-13,89790,2.66847], Tmin=(1183.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(760.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C](C)CC(9111)',
    structure = SMILES('[CH2][CH]C([CH2])[C](C)CC'),
    E0 = (546.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47047,0.0686612,2.53716e-05,-1.93541e-07,1.81961e-10,65766.9,34.6866], Tmin=(100,'K'), Tmax=(433.871,'K')), NASAPolynomial(coeffs=[4.62695,0.061726,-2.72829e-05,5.11347e-09,-3.5319e-13,65284.4,19.6898], Tmin=(433.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C]([CH2])C(C)[CH]C(5480)',
    structure = SMILES('[CH2]C[C]([CH2])C(C)[CH]C'),
    E0 = (546.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,230.341,2005.9,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47047,0.0686612,2.53716e-05,-1.93541e-07,1.81961e-10,65766.9,34.6866], Tmin=(100,'K'), Tmax=(433.871,'K')), NASAPolynomial(coeffs=[4.62695,0.061726,-2.72829e-05,5.11347e-09,-3.5319e-13,65284.4,19.6898], Tmin=(433.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][C](C)C(C)[CH]C(4478)',
    structure = SMILES('[CH2][CH][C](C)C(C)[CH]C'),
    E0 = (535.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,491.672,2105.78,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1535,0.0681231,-3.27785e-05,6.67117e-09,-5.08737e-13,64457.7,32.4957], Tmin=(100,'K'), Tmax=(3091.01,'K')), NASAPolynomial(coeffs=[40.2879,0.0187727,-8.82899e-06,1.50559e-09,-9.09322e-14,40883.7,-194.514], Tmin=(3091.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)C[CH2](4516)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)C[CH2]'),
    E0 = (565.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,207.626,1024.57,2155.31],'cm^-1')),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118619,0.0835433,-5.7073e-05,2.12728e-08,-3.39318e-12,68214.6,39.8597], Tmin=(100,'K'), Tmax=(1404,'K')), NASAPolynomial(coeffs=[12.3641,0.0486558,-1.97999e-05,3.57419e-09,-2.41707e-13,64776.1,-23.3727], Tmin=(1404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C(C)CC(9112)',
    structure = SMILES('[CH2][CH][C]([CH2])C(C)CC'),
    E0 = (546.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47047,0.0686612,2.53716e-05,-1.93541e-07,1.81961e-10,65766.9,34.6866], Tmin=(100,'K'), Tmax=(433.871,'K')), NASAPolynomial(coeffs=[4.62695,0.061726,-2.72829e-05,5.11347e-09,-3.5319e-13,65284.4,19.6898], Tmin=(433.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])CC(3979)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])CC'),
    E0 = (565.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,344.132,543.723,3250.89],'cm^-1')),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.78,'J/mol'), sigma=(6.74357,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.91 K, Pc=26.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00459272,0.0855863,-6.29211e-05,2.6579e-08,-4.83184e-12,68199.7,40.2672], Tmin=(100,'K'), Tmax=(1258.9,'K')), NASAPolynomial(coeffs=[11.5679,0.0488451,-1.91433e-05,3.39583e-09,-2.27965e-13,65288.3,-18.1815], Tmin=(1258.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])[C](C)[CH]C(5481)',
    structure = SMILES('[CH2]CC([CH2])[C](C)[CH]C'),
    E0 = (546.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,230.341,2005.9,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47047,0.0686612,2.53716e-05,-1.93541e-07,1.81961e-10,65766.9,34.6866], Tmin=(100,'K'), Tmax=(433.871,'K')), NASAPolynomial(coeffs=[4.62695,0.061726,-2.72829e-05,5.11347e-09,-3.5319e-13,65284.4,19.6898], Tmin=(433.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)[C](C)[CH]C(4477)',
    structure = SMILES('[CH2][CH]C(C)[C](C)[CH]C'),
    E0 = (535.634,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,491.672,2105.78,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0684061,'amu*angstrom^2'), symmetry=1, barrier=(1.71805,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1535,0.0681231,-3.27785e-05,6.67117e-09,-5.08737e-13,64457.7,32.4957], Tmin=(100,'K'), Tmax=(3091.01,'K')), NASAPolynomial(coeffs=[40.2879,0.0187727,-8.82899e-06,1.50559e-09,-9.09322e-14,40883.7,-194.514], Tmin=(3091.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.634,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH]C)C(C)[CH]C(4491)',
    structure = SMILES('[CH2][C]([CH]C)C(C)[CH]C'),
    E0 = (535.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,3507.17,3768.26],'cm^-1')),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.317977,0.0882846,-9.48083e-05,7.89595e-08,-2.94426e-11,64527.7,39.2693], Tmin=(100,'K'), Tmax=(783.57,'K')), NASAPolynomial(coeffs=[2.79009,0.0648381,-2.91986e-05,5.50471e-09,-3.80668e-13,64472.7,30.0666], Tmin=(783.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])C([CH2])[CH]C(4022)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])[CH]C'),
    E0 = (565.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,344.132,543.723,3250.89],'cm^-1')),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00459272,0.0855863,-6.29211e-05,2.6579e-08,-4.83184e-12,68199.7,40.2672], Tmin=(100,'K'), Tmax=(1258.9,'K')), NASAPolynomial(coeffs=[11.5679,0.0488451,-1.91433e-05,3.39583e-09,-2.27965e-13,65288.3,-18.1815], Tmin=(1258.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C([CH2])[CH]C(3840)',
    structure = SMILES('[CH2][CH]C(C)C([CH2])[CH]C'),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.78,'J/mol'), sigma=(6.74357,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.91 K, Pc=26.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497458,0.0790439,-4.94967e-05,1.65362e-08,-2.40218e-12,66910,39.2611], Tmin=(100,'K'), Tmax=(1473.06,'K')), NASAPolynomial(coeffs=[10.5374,0.0517809,-2.17351e-05,3.97208e-09,-2.69863e-13,63952.1,-13.0649], Tmin=(1473.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)[C](C)[CH]C(4490)',
    structure = SMILES('[CH2]C([CH]C)[C](C)[CH]C'),
    E0 = (535.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,3507.17,3768.26],'cm^-1')),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0313432,'amu*angstrom^2'), symmetry=1, barrier=(3.80382,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.317977,0.0882846,-9.48083e-05,7.89595e-08,-2.94426e-11,64527.7,39.2693], Tmin=(100,'K'), Tmax=(783.57,'K')), NASAPolynomial(coeffs=[2.79009,0.0648381,-2.91986e-05,5.50471e-09,-3.80668e-13,64472.7,30.0666], Tmin=(783.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(C)C([CH2])C[CH2](3980)',
    structure = SMILES('[CH2][CH]C(C)C([CH2])C[CH2]'),
    E0 = (565.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,207.626,1024.57,2155.31],'cm^-1')),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118619,0.0835433,-5.7073e-05,2.12728e-08,-3.39318e-12,68214.6,39.8597], Tmin=(100,'K'), Tmax=(1404,'K')), NASAPolynomial(coeffs=[12.3641,0.0486558,-1.97999e-05,3.57419e-09,-2.41707e-13,64776.1,-23.3727], Tmin=(1404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C(C)[CH][CH2](3847)',
    structure = SMILES('[CH2][CH]C(C)C(C)[CH][CH2]'),
    E0 = (555.458,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,239.568,1161.42,3298.7],'cm^-1')),
        HinderedRotor(inertia=(0.0576644,'amu*angstrom^2'), symmetry=1, barrier=(1.84081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0576644,'amu*angstrom^2'), symmetry=1, barrier=(1.84081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0576644,'amu*angstrom^2'), symmetry=1, barrier=(1.84081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0576644,'amu*angstrom^2'), symmetry=1, barrier=(1.84081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0576644,'amu*angstrom^2'), symmetry=1, barrier=(1.84081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0576644,'amu*angstrom^2'), symmetry=1, barrier=(1.84081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0576644,'amu*angstrom^2'), symmetry=1, barrier=(1.84081,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.78,'J/mol'), sigma=(6.74357,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.91 K, Pc=26.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.580834,0.0773921,-4.49893e-05,1.27549e-08,-1.49265e-12,66925.8,38.2656], Tmin=(100,'K'), Tmax=(1826.06,'K')), NASAPolynomial(coeffs=[14.4822,0.0469412,-1.99757e-05,3.62287e-09,-2.42414e-13,61848.8,-37.1713], Tmin=(1826.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)C([CH2])[CH]C(3835)',
    structure = SMILES('[CH2]C([CH]C)C([CH2])[CH]C'),
    E0 = (555.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,219.021,704.765,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00367365,'amu*angstrom^2'), symmetry=1, barrier=(0.130602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00367365,'amu*angstrom^2'), symmetry=1, barrier=(0.130602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00367365,'amu*angstrom^2'), symmetry=1, barrier=(0.130602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00367365,'amu*angstrom^2'), symmetry=1, barrier=(0.130602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00367365,'amu*angstrom^2'), symmetry=1, barrier=(0.130602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00367365,'amu*angstrom^2'), symmetry=1, barrier=(0.130602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00367365,'amu*angstrom^2'), symmetry=1, barrier=(0.130602,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.78,'J/mol'), sigma=(6.74357,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.91 K, Pc=26.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.349304,0.0813884,-5.60404e-05,2.24034e-08,-3.98656e-12,66897,39.1062], Tmin=(100,'K'), Tmax=(1242.73,'K')), NASAPolynomial(coeffs=[9.22457,0.0528211,-2.15588e-05,3.90549e-09,-2.65292e-13,64691.1,-5.64038], Tmin=(1242.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (555.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (555.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1000.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (750.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (712.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1004.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1013.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1008.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1021.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1020.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1020.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1010.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1010.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (560.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (563.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (563.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (578.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (578.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (618.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (618.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (618.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (580.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (580.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (683.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (574.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (697.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (659.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (680.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (643.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (663.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (556.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (900.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (912.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (952.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (952.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (972.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (972.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (709.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (713.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (696.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (672.542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (699.335,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (712.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (671.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (673.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (712.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (641.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (638.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (686.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (639.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (609.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (609.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH]C=C(3743)', 'butene2t(396)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]C(C=C)C(C)[CH]C(3845)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', '[CH2][CH]C([CH2])C[CH]C(3910)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH]C([CH2])[CH]C(C)C(9276)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHCH3(T)(21)', '[CH2][CH]C([CH2])[CH]C(3871)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH2](502)', '[CH2][CH]C(C)[CH]C(1188)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C(C)[CH]C(1275)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(17)', '[CH]C(C)C([CH2])[CH][CH2](9061)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C(C)[C]C(9277)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][C]C([CH2])C(C)[CH]C(9278)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]C([CH][CH2])C(C)[CH]C(9279)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH][CH]C([CH2])C(C)[CH]C(9280)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]C1CC1C(C)[CH]C(5484)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]C1C([CH2])C(C)C1C(4509)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH]C1CC(C)C1C(5502)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]CC(=C)C(C)[CH]C(9281)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]C=C(C)C(C)[CH]C(4476)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]CC([CH2])C(C)=CC(9282)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH]C(C)C(C)=CC(4475)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][C](C=C)C(C)CC(5096)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]CC([CH2])C(C)C=C(5100)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH]C(C)C(C)C=C(3843)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C(C=C)[C](C)[CH]C(5487)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.72e+07,'cm^3/(mol*s)'), n=1.64, Ea=(9.07928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2568 used for Cds-CsCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][C](C=C)C(C)[CH]C(5488)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH2][CH]C(C)C([CH2])C=C(4528)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH2(T)(20)', '[CH2]C=CC(C)[CH]C(1272)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH][CH2](502)', 'C=CC(C)[CH]C(1184)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'butene2t(396)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.0171758,'m^3/(mol*s)'), n=2.41179, Ea=(16.3987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH3(17)', '[CH2]C([CH][CH]C)C=C(5043)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(10100,'cm^3/(mol*s)'), n=2.41, Ea=(28.4512,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 435 used for Cds-CsH_Cds-CsH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]C=C(3743)', 'C[CH][CH]C(1186)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.54197,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'C[CH][CH]C(1186)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.3e+09,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH3(17)', '[CH2][CH]C([CH2])[CH][CH]C(9076)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[C](C)[CH]C(9283)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[CH2][CH][C]([CH2])C(C)[CH]C(9284)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C(C)[CH][CH2](4513)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH]C([CH2])[C](C)CC(9111)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 150 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]C[C]([CH2])C(C)[CH]C(5480)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH][C](C)C(C)[CH]C(4478)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH]C([CH2])C(C)C[CH2](4516)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(5.04e-11,'s^-1'), n=6.833, Ea=(117.248,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 33 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH][C]([CH2])C(C)CC(9112)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]CC([CH2])[C](C)[CH]C(5481)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][CH]C(C)[C](C)[CH]C(4477)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2][C]([CH]C)C(C)[CH]C(4491)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]CC([CH2])C([CH2])[CH]C(4022)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(927.918,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][CH]C(C)C([CH2])[CH]C(3840)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    products = ['[CH2]C([CH]C)[C](C)[CH]C(4490)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2][CH]C(C)C([CH2])C[CH2](3980)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]C(C)C(C)[CH][CH2](3847)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(182547,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH]C)C([CH2])[CH]C(3835)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(182547,'s^-1'), n=1.79, Ea=(54.1828,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R5HJ_3;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2078',
    isomers = [
        '[CH2][CH]C([CH2])C(C)[CH]C(4479)',
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
    label = 'PDepNetwork #2078',
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


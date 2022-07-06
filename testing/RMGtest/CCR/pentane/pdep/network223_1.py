species(
    label = '[CH2][CH]CC[CH]CC(330)',
    structure = SMILES('[CH2][CH]CC[CH]CC'),
    E0 = (378.506,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,405.986,1418.12,2619.56,3831.92],'cm^-1')),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33664,0.0584013,-2.57862e-05,4.71002e-09,-3.04688e-13,45560.5,29.2647], Tmin=(100,'K'), Tmax=(2817.44,'K')), NASAPolynomial(coeffs=[41.6441,0.0104726,-4.46283e-06,6.56805e-10,-3.30881e-14,20284.7,-206.636], Tmin=(2817.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'allyl(82)',
    structure = SMILES('[CH2]C=C'),
    E0 = (156.927,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.332071,'amu*angstrom^2'), symmetry=1, barrier=(25.4371,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.29615,0.00579223,4.3392e-05,-5.9989e-08,2.33815e-11,18908.2,9.01994], Tmin=(100,'K'), Tmax=(942.181,'K')), NASAPolynomial(coeffs=[8.06863,0.0101837,-2.84795e-06,5.0088e-10,-3.79628e-14,16914.7,-19.5272], Tmin=(942.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.927,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), label="""allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=CCC[CH]CC(325)',
    structure = SMILES('C=CCC[CH]CC'),
    E0 = (106.596,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.940976,0.0651351,-3.34616e-05,8.01494e-09,-7.57585e-13,12931.9,30.3726], Tmin=(100,'K'), Tmax=(2362.79,'K')), NASAPolynomial(coeffs=[20.5339,0.0319666,-1.24052e-05,2.07391e-09,-1.28993e-13,3672.96,-80.9992], Tmin=(2362.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.596,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC)"""),
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
    label = '[CH2][CH]CC[CH]C(161)',
    structure = SMILES('[CH2][CH]CC[CH]C'),
    E0 = (402.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,289.708,1713.12,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00228109,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00228109,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00228109,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00228109,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00228109,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3202.63,'J/mol'), sigma=(6.0508,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=500.24 K, Pc=32.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51744,0.0475615,-2.00024e-05,3.3543e-09,-1.76064e-13,48420.3,26.5952], Tmin=(100,'K'), Tmax=(2622.93,'K')), NASAPolynomial(coeffs=[28.602,0.015465,-6.04071e-06,9.22376e-10,-5.07077e-14,32094,-129.438], Tmin=(2622.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC([CH2])CC(309)',
    structure = SMILES('[CH2][CH]CC([CH2])CC'),
    E0 = (386.783,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,357.162,506.393,3183.71],'cm^-1')),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3372.64,'J/mol'), sigma=(6.39812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=526.80 K, Pc=29.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.856328,0.0704044,-4.55259e-05,1.70507e-08,-2.88132e-12,46631.1,33.6409], Tmin=(100,'K'), Tmax=(1279,'K')), NASAPolynomial(coeffs=[7.9267,0.0482922,-1.95929e-05,3.53334e-09,-2.39146e-13,44822.5,-2.20939], Tmin=(1279,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C[CH]CC(287)',
    structure = SMILES('[CH2]C([CH2])C[CH]CC'),
    E0 = (383.284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,357.527,506.612,3176.77],'cm^-1')),
        HinderedRotor(inertia=(0.02753,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.02753,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.02753,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.02753,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.02753,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.02753,'amu*angstrom^2'), symmetry=1, barrier=(2.90446,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3372.64,'J/mol'), sigma=(6.39812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=526.80 K, Pc=29.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.791642,0.0726374,-5.49104e-05,2.77143e-08,-6.5126e-12,46212.1,33.0159], Tmin=(100,'K'), Tmax=(961.348,'K')), NASAPolynomial(coeffs=[6.14618,0.050358,-2.01475e-05,3.60712e-09,-2.43475e-13,45182.6,7.39441], Tmin=(961.348,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CCC([CH2])C(612)',
    structure = SMILES('[CH2][CH]CCC([CH2])C'),
    E0 = (383.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,357.162,506.393,3183.71],'cm^-1')),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0329935,'amu*angstrom^2'), symmetry=1, barrier=(2.54994,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.856328,0.0704044,-4.55259e-05,1.70507e-08,-2.88132e-12,46228.5,33.6409], Tmin=(100,'K'), Tmax=(1279,'K')), NASAPolynomial(coeffs=[7.9267,0.0482922,-1.95929e-05,3.53334e-09,-2.39146e-13,44419.9,-2.20939], Tmin=(1279,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]CC(26)',
    structure = SMILES('[CH]CC'),
    E0 = (328.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,329.486,329.64,2894.72],'cm^-1')),
        HinderedRotor(inertia=(0.00155193,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00155082,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00369,0.0161456,1.44259e-05,-2.625e-08,1.01938e-11,39516.8,12.6847], Tmin=(100,'K'), Tmax=(1003.28,'K')), NASAPolynomial(coeffs=[6.63834,0.0156476,-5.75059e-06,1.05872e-09,-7.51273e-14,38083.2,-8.3721], Tmin=(1003.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C[CH2](66)',
    structure = SMILES('[CH2][CH]C[CH2]'),
    E0 = (460.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.000785383,'amu*angstrom^2'), symmetry=1, barrier=(2.49202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108433,'amu*angstrom^2'), symmetry=1, barrier=(2.4931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000786299,'amu*angstrom^2'), symmetry=1, barrier=(2.4954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.62775,0.0341381,-3.63084e-05,3.45494e-08,-1.41477e-11,55480.4,19.9882], Tmin=(100,'K'), Tmax=(807.221,'K')), NASAPolynomial(coeffs=[1.67633,0.030324,-1.33729e-05,2.51877e-09,-1.74066e-13,55911.9,26.0956], Tmin=(807.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2]C[CH]CC(101)',
    structure = SMILES('[CH2]C[CH]CC'),
    E0 = (231.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1730.81,1731.92],'cm^-1')),
        HinderedRotor(inertia=(0.154521,'amu*angstrom^2'), symmetry=1, barrier=(3.55274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154052,'amu*angstrom^2'), symmetry=1, barrier=(3.54195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154831,'amu*angstrom^2'), symmetry=1, barrier=(3.55987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154739,'amu*angstrom^2'), symmetry=1, barrier=(3.55776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19946,0.0413401,-1.64432e-05,2.00796e-09,7.91291e-14,27920.7,22.6765], Tmin=(100,'K'), Tmax=(2000.48,'K')), NASAPolynomial(coeffs=[13.367,0.0260959,-1.03257e-05,1.7398e-09,-1.08622e-13,22034.8,-42.4877], Tmin=(2000.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJ)"""),
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
    label = '[CH]CC[CH][CH2](504)',
    structure = SMILES('[CH]CC[CH][CH2]'),
    E0 = (679.824,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,220.978,406.057,1393.01,1671.69,3531.85],'cm^-1')),
        HinderedRotor(inertia=(0.0652034,'amu*angstrom^2'), symmetry=1, barrier=(2.249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0652034,'amu*angstrom^2'), symmetry=1, barrier=(2.249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0652034,'amu*angstrom^2'), symmetry=1, barrier=(2.249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0652034,'amu*angstrom^2'), symmetry=1, barrier=(2.249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16392,0.0433951,-2.95388e-05,1.27191e-08,-2.68212e-12,81827.4,23.9734], Tmin=(100,'K'), Tmax=(1002.93,'K')), NASAPolynomial(coeffs=[4.68692,0.0333326,-1.4489e-05,2.71513e-09,-1.8843e-13,81321.4,11.794], Tmin=(1002.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][CH]CC[C]CC(613)',
    structure = SMILES('[CH2][CH]CC[C]CC'),
    E0 = (632.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.934559,0.0721887,-5.44583e-05,2.52111e-08,-5.41294e-12,76150,31.7804], Tmin=(100,'K'), Tmax=(1027.63,'K')), NASAPolynomial(coeffs=[6.64843,0.0499477,-2.19936e-05,4.14988e-09,-2.89182e-13,74975.7,4.05856], Tmin=(1027.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]CC[CH]CC(614)',
    structure = SMILES('[CH2][C]CC[CH]CC'),
    E0 = (632.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918191,0.0731598,-6.0772e-05,3.42973e-08,-9.27863e-12,76151.4,31.6763], Tmin=(100,'K'), Tmax=(831.154,'K')), NASAPolynomial(coeffs=[5.19533,0.052575,-2.36212e-05,4.49779e-09,-3.15049e-13,75440.4,11.8328], Tmin=(831.154,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(632.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]CC[CH]CC(615)',
    structure = SMILES('[CH][CH]CC[CH]CC'),
    E0 = (621.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,210.361,1016.97,1092.11,1225.49,1382.97,1635.28,1863.41],'cm^-1')),
        HinderedRotor(inertia=(0.12444,'amu*angstrom^2'), symmetry=1, barrier=(3.3487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12444,'amu*angstrom^2'), symmetry=1, barrier=(3.3487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12444,'amu*angstrom^2'), symmetry=1, barrier=(3.3487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12444,'amu*angstrom^2'), symmetry=1, barrier=(3.3487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12444,'amu*angstrom^2'), symmetry=1, barrier=(3.3487,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12444,'amu*angstrom^2'), symmetry=1, barrier=(3.3487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.855163,0.0755928,-8.17354e-05,6.88215e-08,-2.57516e-11,74853.2,33.9304], Tmin=(100,'K'), Tmax=(792.447,'K')), NASAPolynomial(coeffs=[2.66998,0.0560698,-2.51662e-05,4.72977e-09,-3.26252e-13,74890.9,27.65], Tmin=(792.447,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CCC1CC(564)',
    structure = SMILES('[CH2]C1CCC1CC'),
    E0 = (129.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21261,0.0406224,6.1602e-05,-1.04593e-07,4.19803e-11,15698,26.4383], Tmin=(100,'K'), Tmax=(962.844,'K')), NASAPolynomial(coeffs=[14.4624,0.0369592,-1.27368e-05,2.30193e-09,-1.65182e-13,10764.9,-49.3505], Tmin=(962.844,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CCC=CCC(616)',
    structure = SMILES('[CH2]CCC=CCC'),
    E0 = (106.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533926,0.0656075,-2.40325e-05,-6.88139e-09,5.06987e-12,12915.3,31.2771], Tmin=(100,'K'), Tmax=(1142.44,'K')), NASAPolynomial(coeffs=[12.6745,0.0416027,-1.68086e-05,3.07976e-09,-2.12213e-13,8933.86,-34.1956], Tmin=(1142.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(106.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CCCCC(617)',
    structure = SMILES('[CH2]C=CCCCC'),
    E0 = (51.3913,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.469479,0.0657121,-1.82925e-05,-1.58856e-08,8.81728e-12,6318.15,29.3305], Tmin=(100,'K'), Tmax=(1079.49,'K')), NASAPolynomial(coeffs=[13.4612,0.040963,-1.64054e-05,3.02189e-09,-2.10222e-13,2150.38,-40.6533], Tmin=(1079.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(51.3913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CCCC=CC(473)',
    structure = SMILES('[CH2]CCCC=CC'),
    E0 = (105.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.307849,0.0707335,-4.1195e-05,1.17634e-08,-1.35121e-12,12787.2,31.4565], Tmin=(100,'K'), Tmax=(1980.69,'K')), NASAPolynomial(coeffs=[19.0062,0.0329714,-1.25968e-05,2.13748e-09,-1.36214e-13,5380.2,-71.5311], Tmin=(1980.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(105.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(461.453,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC=CCC(618)',
    structure = SMILES('[CH2][CH]CC=CCC'),
    E0 = (300.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.933212,0.0633457,-3.42967e-05,8.97968e-09,-9.48969e-13,36281.5,32.5856], Tmin=(100,'K'), Tmax=(2085.54,'K')), NASAPolynomial(coeffs=[16.1352,0.0341885,-1.33255e-05,2.27592e-09,-1.45359e-13,29940.8,-51.9287], Tmin=(2085.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC[CH]CC(619)',
    structure = SMILES('[CH2]C=CC[CH]CC'),
    E0 = (245.85,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83813,0.0645622,-3.52329e-05,9.25012e-09,-9.76147e-13,29686.5,30.5885], Tmin=(100,'K'), Tmax=(2104.07,'K')), NASAPolynomial(coeffs=[17.1845,0.0334865,-1.30789e-05,2.23071e-09,-1.42119e-13,22807.8,-60.4327], Tmin=(2104.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.85,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CCC=CC(620)',
    structure = SMILES('[CH2][CH]CCC=CC'),
    E0 = (299.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35629,0.0618883,-3.27093e-05,8.27384e-09,-8.54501e-13,36122,30.3539], Tmin=(100,'K'), Tmax=(2021.78,'K')), NASAPolynomial(coeffs=[12.4595,0.0399211,-1.64114e-05,2.89977e-09,-1.8998e-13,31632.3,-31.0293], Tmin=(2021.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH2](497)',
    structure = SMILES('[CH2][CH][CH2]'),
    E0 = (484.819,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.000863049,'amu*angstrom^2'), symmetry=1, barrier=(2.40754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000864751,'amu*angstrom^2'), symmetry=1, barrier=(2.41365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42318,0.0135015,1.25385e-06,-4.52647e-09,1.21417e-12,58330.2,15.4092], Tmin=(100,'K'), Tmax=(1578.37,'K')), NASAPolynomial(coeffs=[5.16145,0.0152713,-6.29635e-06,1.14122e-09,-7.61383e-14,57012.3,3.79304], Tmin=(1578.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.819,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][CH]CCC=C(487)',
    structure = SMILES('[CH2][CH]CCC=C'),
    E0 = (335.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,1197.45,1197.67],'cm^-1')),
        HinderedRotor(inertia=(0.131095,'amu*angstrom^2'), symmetry=1, barrier=(3.01413,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131148,'amu*angstrom^2'), symmetry=1, barrier=(3.01536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13113,'amu*angstrom^2'), symmetry=1, barrier=(3.01493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131128,'amu*angstrom^2'), symmetry=1, barrier=(3.01489,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55938,0.0510078,-2.66124e-05,6.6167e-09,-6.58045e-13,40454.3,27.8617], Tmin=(100,'K'), Tmax=(2219.44,'K')), NASAPolynomial(coeffs=[14.8429,0.027067,-1.04318e-05,1.75636e-09,-1.10562e-13,34558,-46.814], Tmin=(2219.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH][CH]CC(621)',
    structure = SMILES('[CH2][CH]C[CH][CH]CC'),
    E0 = (572.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,2351.73,2410.93,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15593,0.0742926,-9.87011e-05,1.0042e-07,-4.03765e-11,69002.7,36.1938], Tmin=(100,'K'), Tmax=(843.06,'K')), NASAPolynomial(coeffs=[-2.56545,0.0635643,-2.91099e-05,5.4533e-09,-3.72438e-13,70638.9,59.4945], Tmin=(843.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH]CC(622)',
    structure = SMILES('[CH2][CH][CH]C[CH]CC'),
    E0 = (572.964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,2351.73,2410.93,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0242903,'amu*angstrom^2'), symmetry=1, barrier=(11.4658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15593,0.0742926,-9.87011e-05,1.0042e-07,-4.03765e-11,69002.7,36.1938], Tmin=(100,'K'), Tmax=(843.06,'K')), NASAPolynomial(coeffs=[-2.56545,0.0635643,-2.91099e-05,5.4533e-09,-3.72438e-13,70638.9,59.4945], Tmin=(843.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[CH][CH2](483)',
    structure = SMILES('[CH2][CH]CC[CH][CH2]'),
    E0 = (607.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,406.91,2706.24,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168246,'amu*angstrom^2'), symmetry=1, barrier=(3.0235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96636,0.0441879,-1.8776e-05,3.20919e-09,-1.84564e-13,73082,26.0155], Tmin=(100,'K'), Tmax=(2870.12,'K')), NASAPolynomial(coeffs=[37.1311,0.00432169,-1.9902e-06,2.50805e-10,-8.8047e-15,50279.3,-180.391], Tmin=(2870.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[CH][CH]C(623)',
    structure = SMILES('[CH2][CH]CC[CH][CH]C'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36111,0.0497551,-1.77213e-05,1.52484e-09,1.23143e-13,68896.7,28.2554], Tmin=(100,'K'), Tmax=(2624.33,'K')), NASAPolynomial(coeffs=[40.8399,0.00972266,-4.60954e-06,6.75854e-10,-3.24064e-14,43339.4,-199.933], Tmin=(2624.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[CH]C[CH2](624)',
    structure = SMILES('[CH2][CH]CC[CH]C[CH2]'),
    E0 = (583.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000701331,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000701331,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000701331,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000701331,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000701331,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000701331,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79772,0.0549377,-2.43973e-05,4.46585e-09,-2.94266e-13,70221.3,29.3293], Tmin=(100,'K'), Tmax=(3038.12,'K')), NASAPolynomial(coeffs=[51.8439,-0.00204164,-1.49216e-08,-6.16343e-11,1.05805e-14,36914.5,-267.561], Tmin=(3038.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(583.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]CCC(576)',
    structure = SMILES('[CH2][CH]C[CH]CCC'),
    E0 = (378.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33664,0.0584013,-2.57862e-05,4.71002e-09,-3.04688e-13,45560.5,29.2647], Tmin=(100,'K'), Tmax=(2817.44,'K')), NASAPolynomial(coeffs=[41.6441,0.0104726,-4.46283e-06,6.56805e-10,-3.30881e-14,20284.7,-206.636], Tmin=(2817.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C[CH]CC(625)',
    structure = SMILES('[CH2]C[CH]C[CH]CC'),
    E0 = (378.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,454.626,2277.96,3712.22,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55545,0.056742,-2.35883e-05,3.62852e-09,-1.39781e-13,45551.6,28.3085], Tmin=(100,'K'), Tmax=(2625.11,'K')), NASAPolynomial(coeffs=[38.1521,0.0144893,-6.29469e-06,9.76222e-10,-5.28569e-14,22732.1,-185.647], Tmin=(2625.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CCC[CH]C(626)',
    structure = SMILES('[CH2][CH]CCC[CH]C'),
    E0 = (378.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14569,0.0598559,-2.76193e-05,5.57225e-09,-4.28388e-13,45567.3,30.1088], Tmin=(100,'K'), Tmax=(2926.01,'K')), NASAPolynomial(coeffs=[24.7959,0.0288922,-1.17461e-05,1.95572e-09,-1.19394e-13,32312.2,-103.484], Tmin=(2926.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.494,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CCCC(627)',
    structure = SMILES('[CH2][CH][CH]CCCC'),
    E0 = (378.506,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,405.986,1418.12,2619.56,3831.92],'cm^-1')),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33664,0.0584013,-2.57862e-05,4.71002e-09,-3.04688e-13,45560.5,29.2647], Tmin=(100,'K'), Tmax=(2817.44,'K')), NASAPolynomial(coeffs=[41.6441,0.0104726,-4.46283e-06,6.56805e-10,-3.30881e-14,20284.7,-206.636], Tmin=(2817.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC[CH][CH]CC(628)',
    structure = SMILES('[CH2]CC[CH][CH]CC'),
    E0 = (378.518,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,454.626,2277.96,3712.22,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0195964,'amu*angstrom^2'), symmetry=1, barrier=(13.4561,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55545,0.056742,-2.35883e-05,3.62852e-09,-1.39781e-13,45551.6,28.3085], Tmin=(100,'K'), Tmax=(2625.11,'K')), NASAPolynomial(coeffs=[38.1521,0.0144893,-6.29469e-06,9.76222e-10,-5.28569e-14,22732.1,-185.647], Tmin=(2625.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.518,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CCCC[CH2](629)',
    structure = SMILES('[CH2][CH]CCCC[CH2]'),
    E0 = (389.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,215.186,976.149,1365.66,1942.14],'cm^-1')),
        HinderedRotor(inertia=(0.107129,'amu*angstrom^2'), symmetry=1, barrier=(3.19012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107129,'amu*angstrom^2'), symmetry=1, barrier=(3.19012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107129,'amu*angstrom^2'), symmetry=1, barrier=(3.19012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107129,'amu*angstrom^2'), symmetry=1, barrier=(3.19012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107129,'amu*angstrom^2'), symmetry=1, barrier=(3.19012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107129,'amu*angstrom^2'), symmetry=1, barrier=(3.19012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14415,0.0686779,-4.21076e-05,1.44222e-08,-2.3004e-12,46918.8,32.892], Tmin=(100,'K'), Tmax=(1265.81,'K')), NASAPolynomial(coeffs=[6.293,0.0524074,-2.28266e-05,4.26748e-09,-2.94807e-13,45615.3,6.83819], Tmin=(1265.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C[CH]CC(327)',
    structure = SMILES('C[CH][CH]C[CH]CC'),
    E0 = (367.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1264.38,2963.5,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.052059,'amu*angstrom^2'), symmetry=1, barrier=(6.55256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.052059,'amu*angstrom^2'), symmetry=1, barrier=(6.55256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.052059,'amu*angstrom^2'), symmetry=1, barrier=(6.55256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.052059,'amu*angstrom^2'), symmetry=1, barrier=(6.55256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.052059,'amu*angstrom^2'), symmetry=1, barrier=(6.55256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.052059,'amu*angstrom^2'), symmetry=1, barrier=(6.55256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13535,0.0733609,-8.75664e-05,8.73427e-08,-3.56414e-11,44319.3,34.4814], Tmin=(100,'K'), Tmax=(831.598,'K')), NASAPolynomial(coeffs=[-2.44949,0.0660173,-2.99721e-05,5.61845e-09,-3.84941e-13,45765.7,56.2267], Tmin=(831.598,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C[CH][CH]CC(326)',
    structure = SMILES('C[CH]C[CH][CH]CC'),
    E0 = (367.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1264.43,2958.44,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0520598,'amu*angstrom^2'), symmetry=1, barrier=(6.5525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0520598,'amu*angstrom^2'), symmetry=1, barrier=(6.5525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0520598,'amu*angstrom^2'), symmetry=1, barrier=(6.5525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0520598,'amu*angstrom^2'), symmetry=1, barrier=(6.5525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0520598,'amu*angstrom^2'), symmetry=1, barrier=(6.5525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0520598,'amu*angstrom^2'), symmetry=1, barrier=(6.5525,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13522,0.0733625,-8.75729e-05,8.73526e-08,-3.56464e-11,44319.3,34.4819], Tmin=(100,'K'), Tmax=(831.575,'K')), NASAPolynomial(coeffs=[-2.44925,0.0660169,-2.99718e-05,5.61838e-09,-3.84936e-13,45765.6,56.2254], Tmin=(831.575,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]CCC[CH][CH]C(630)',
    structure = SMILES('[CH2]CCC[CH][CH]C'),
    E0 = (378.506,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,405.986,1418.12,2619.56,3831.92],'cm^-1')),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33664,0.0584013,-2.57862e-05,4.71002e-09,-3.04688e-13,45560.5,29.2647], Tmin=(100,'K'), Tmax=(2817.44,'K')), NASAPolynomial(coeffs=[41.6441,0.0104726,-4.46283e-06,6.56805e-10,-3.30881e-14,20284.7,-206.636], Tmin=(2817.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]CC[CH]C(328)',
    structure = SMILES('C[CH][CH]CC[CH]C'),
    E0 = (367.706,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,875.63,1946.17,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129737,'amu*angstrom^2'), symmetry=1, barrier=(5.05239,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09928,0.0728686,-8.23053e-05,7.88527e-08,-3.17483e-11,44320.4,34.783], Tmin=(100,'K'), Tmax=(826.938,'K')), NASAPolynomial(coeffs=[-1.37392,0.0640962,-2.87802e-05,5.37859e-09,-3.68367e-13,45438.5,50.5315], Tmin=(826.938,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.706,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]CCC[CH2](631)',
    structure = SMILES('[CH2]C[CH]CCC[CH2]'),
    E0 = (389.306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,215.112,737.697,1221.41,3419.63],'cm^-1')),
        HinderedRotor(inertia=(0.0556125,'amu*angstrom^2'), symmetry=1, barrier=(2.24095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0556125,'amu*angstrom^2'), symmetry=1, barrier=(2.24095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0556125,'amu*angstrom^2'), symmetry=1, barrier=(2.24095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0556125,'amu*angstrom^2'), symmetry=1, barrier=(2.24095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0556125,'amu*angstrom^2'), symmetry=1, barrier=(2.24095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0556125,'amu*angstrom^2'), symmetry=1, barrier=(2.24095,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06613,0.0701238,-4.90858e-05,2.32862e-08,-5.72113e-12,46923.5,33.0283], Tmin=(100,'K'), Tmax=(857.101,'K')), NASAPolynomial(coeffs=[3.97064,0.0565687,-2.53633e-05,4.83451e-09,-3.39128e-13,46425.7,19.4636], Tmin=(857.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]CC[CH]C(329)',
    structure = SMILES('[CH2]C[CH]CC[CH]C'),
    E0 = (378.506,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,405.986,1418.12,2619.56,3831.92],'cm^-1')),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0522738,'amu*angstrom^2'), symmetry=1, barrier=(5.35937,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33664,0.0584013,-2.57862e-05,4.71002e-09,-3.04688e-13,45560.5,29.2647], Tmin=(100,'K'), Tmax=(2817.44,'K')), NASAPolynomial(coeffs=[41.6441,0.0104726,-4.46283e-06,6.56805e-10,-3.30881e-14,20284.7,-206.636], Tmin=(2817.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    E0 = (378.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (378.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (821.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (546.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (540.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (578.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (844.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (844.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (844.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (844.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (844.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (833.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (386.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (441.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (441.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (403.479,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (518.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (465.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (517.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (483.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (449.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (495.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (740.487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (784.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (784.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (744.071,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (784.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (795.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (537.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (515.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (537.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (491.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (494.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (530.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (522.837,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (509.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (405.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (421.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (462.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (441.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['allyl(82)', 'butene1(35)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['C=CCC[CH]CC(325)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', '[CH2][CH]CC[CH]C(161)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(215646,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]CC([CH2])CC(309)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([CH2])C[CH]CC(287)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]CCC([CH2])C(612)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]CC(26)', '[CH2][CH]C[CH2](66)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][CH2](502)', '[CH2]C[CH]CC(101)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C2H5(14)', '[CH]CC[CH][CH2](504)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][CH]CC[C]CC(613)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][C]CC[CH]CC(614)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH][CH]CC[CH]CC(615)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['[CH2]C1CCC1CC(564)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['[CH2]CCC=CCC(616)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['[CH2]C=CCCCC(617)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['[CH2]CCCC=CC(473)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2][CH]CC=CCC(618)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2]C=CC[CH]CC(619)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2][CH]CCC=CC(620)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH][CH2](497)', 'butene1(35)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0344575,'m^3/(mol*s)'), n=2.32603, Ea=(14.6351,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['allyl(82)', '[CH2][CH]CC(37)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH3(17)', '[CH2][CH]CCC=C(487)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.105698,'m^3/(mol*s)'), n=2.13, Ea=(23.0748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH]CC(37)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][CH]C[CH][CH]CC(621)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][CH][CH]C[CH]CC(622)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH3(17)', '[CH2][CH]CC[CH][CH2](483)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2][CH]CC[CH][CH]C(623)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2][CH]CC[CH]C[CH2](624)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['[CH2][CH]C[CH]CCC(576)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C[CH]C[CH]CC(625)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['[CH2][CH]CCC[CH]C(626)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][CH]CCCC(627)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(138391,'s^-1'), n=2.1085, Ea=(113.052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]CC[CH][CH]CC(628)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]CCCC[CH2](629)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 108 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['C[CH][CH]C[CH]CC(327)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.50921e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['C[CH]C[CH][CH]CC(326)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]CCC[CH][CH]C(630)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.0756983,'s^-1'), n=3.26, Ea=(26.8822,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R5HJ_1;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]CC[CH]CC(330)'],
    products = ['C[CH][CH]CC[CH]C(328)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C[CH]CCC[CH2](631)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(2778.79,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_2;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C[CH]CC[CH]C(329)'],
    products = ['[CH2][CH]CC[CH]CC(330)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(64.2,'s^-1'), n=2.1, Ea=(63.1784,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #223',
    isomers = [
        '[CH2][CH]CC[CH]CC(330)',
    ],
    reactants = [
        ('allyl(82)', 'butene1(35)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #223',
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


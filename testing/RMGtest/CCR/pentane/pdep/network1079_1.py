species(
    label = '[CH2]C(CC)C[CH][CH]C(365)',
    structure = SMILES('[CH2]C(CC)C[CH][CH]C'),
    E0 = (352.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,215.647,509.563,1150.94,3517.27],'cm^-1')),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0409758,'amu*angstrom^2'), symmetry=1, barrier=(1.52519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72881,0.0625006,5.5945e-05,-2.48776e-07,2.23958e-10,42430,33.9051], Tmin=(100,'K'), Tmax=(416.553,'K')), NASAPolynomial(coeffs=[3.63613,0.064977,-2.78436e-05,5.14836e-09,-3.53243e-13,42090.7,24.2084], Tmin=(416.553,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
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
    label = '[CH2]C(CC)CC=CC(361)',
    structure = SMILES('[CH2]C(CC)CC=CC'),
    E0 = (75.4999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.47977,0.0857328,-5.40473e-05,1.77245e-08,-2.38864e-12,9252.49,36.2282], Tmin=(100,'K'), Tmax=(1700.09,'K')), NASAPolynomial(coeffs=[18.0544,0.0421252,-1.5572e-05,2.63695e-09,-1.70005e-13,2950.54,-63.0243], Tmin=(1700.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(75.4999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl)"""),
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
    label = '[CH2]C(C)C[CH][CH]C(286)',
    structure = SMILES('[CH2]C(C)C[CH][CH]C'),
    E0 = (372.648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,225.632,1984.27,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0751959,'amu*angstrom^2'), symmetry=1, barrier=(2.35017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0751959,'amu*angstrom^2'), symmetry=1, barrier=(2.35017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0751959,'amu*angstrom^2'), symmetry=1, barrier=(2.35017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0751959,'amu*angstrom^2'), symmetry=1, barrier=(2.35017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0751959,'amu*angstrom^2'), symmetry=1, barrier=(2.35017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0751959,'amu*angstrom^2'), symmetry=1, barrier=(2.35017,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1781,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3372.64,'J/mol'), sigma=(6.39812,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=526.80 K, Pc=29.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.963558,0.071332,-6.29124e-05,4.6136e-08,-1.62607e-11,44924.2,34.0423], Tmin=(100,'K'), Tmax=(797.304,'K')), NASAPolynomial(coeffs=[2.76979,0.0564461,-2.39496e-05,4.39516e-09,-2.99744e-13,44821.3,26.8983], Tmin=(797.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH][CH]C[CH]CCC(4025)',
    structure = SMILES('C[CH][CH]C[CH]CCC'),
    E0 = (343.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73535,0.0639386,-2.46763e-05,2.99099e-09,9.83904e-15,41374.8,30.7956], Tmin=(100,'K'), Tmax=(2603.56,'K')), NASAPolynomial(coeffs=[45.709,0.0155929,-7.00707e-06,1.07443e-09,-5.65365e-14,13006.7,-229.154], Tmin=(2603.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]CC[CH]CC(385)',
    structure = SMILES('C[CH][CH]CC[CH]CC'),
    E0 = (343.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,190.844,270.176,1391.55,3906.95,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0561603,'amu*angstrom^2'), symmetry=1, barrier=(1.37848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0561603,'amu*angstrom^2'), symmetry=1, barrier=(1.37848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0561603,'amu*angstrom^2'), symmetry=1, barrier=(1.37848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0561603,'amu*angstrom^2'), symmetry=1, barrier=(1.37848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0561603,'amu*angstrom^2'), symmetry=1, barrier=(1.37848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0561603,'amu*angstrom^2'), symmetry=1, barrier=(1.37848,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0561603,'amu*angstrom^2'), symmetry=1, barrier=(1.37848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.73535,0.0639386,-2.46763e-05,2.99099e-09,9.83904e-15,41374.8,30.7956], Tmin=(100,'K'), Tmax=(2603.56,'K')), NASAPolynomial(coeffs=[45.709,0.0155929,-7.00707e-06,1.07443e-09,-5.65365e-14,13006.7,-229.154], Tmin=(2603.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([CH]C)C([CH2])CC(344)',
    structure = SMILES('[CH2]C([CH]C)C([CH2])CC'),
    E0 = (360.588,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,243.071,886.017,1827.04],'cm^-1')),
        HinderedRotor(inertia=(0.110662,'amu*angstrom^2'), symmetry=1, barrier=(3.34055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110662,'amu*angstrom^2'), symmetry=1, barrier=(3.34055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110662,'amu*angstrom^2'), symmetry=1, barrier=(3.34055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110662,'amu*angstrom^2'), symmetry=1, barrier=(3.34055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110662,'amu*angstrom^2'), symmetry=1, barrier=(3.34055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110662,'amu*angstrom^2'), symmetry=1, barrier=(3.34055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110662,'amu*angstrom^2'), symmetry=1, barrier=(3.34055,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.78,'J/mol'), sigma=(6.74357,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.91 K, Pc=26.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.196696,0.0867511,-5.8883e-05,2.22418e-08,-3.57781e-12,43524.1,39.2044], Tmin=(100,'K'), Tmax=(1411.83,'K')), NASAPolynomial(coeffs=[13.1622,0.0489025,-1.86707e-05,3.25352e-09,-2.15446e-13,39752,-29.8519], Tmin=(1411.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(360.588,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH][CH]C(3874)',
    structure = SMILES('[CH][CH]C'),
    E0 = (522.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1852.85,1853.13,1853.16],'cm^-1')),
        HinderedRotor(inertia=(0.0369856,'amu*angstrom^2'), symmetry=1, barrier=(8.28142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0370949,'amu*angstrom^2'), symmetry=1, barrier=(8.27684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.19279,0.0163284,-1.76328e-06,-3.55309e-09,1.20622e-12,62877.3,14.3587], Tmin=(100,'K'), Tmax=(1426.65,'K')), NASAPolynomial(coeffs=[5.59458,0.0149965,-6.04303e-06,1.1011e-09,-7.44871e-14,61642.3,-0.00873183], Tmin=(1426.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])CC(117)',
    structure = SMILES('[CH2]C([CH2])CC'),
    E0 = (236.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1687.14],'cm^-1')),
        HinderedRotor(inertia=(0.0765202,'amu*angstrom^2'), symmetry=1, barrier=(8.22176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0765197,'amu*angstrom^2'), symmetry=1, barrier=(8.22183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00407041,'amu*angstrom^2'), symmetry=1, barrier=(8.22185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.648495,'amu*angstrom^2'), symmetry=1, barrier=(69.6757,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71505,0.0440447,-8.69024e-06,-1.71198e-08,9.74731e-12,28518.6,22.6405], Tmin=(100,'K'), Tmax=(931.111,'K')), NASAPolynomial(coeffs=[8.79896,0.0291761,-9.80969e-06,1.63337e-09,-1.07805e-13,26524.8,-14.6525], Tmin=(931.111,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(236.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH][CH]CC([CH2])CC(579)',
    structure = SMILES('[CH][CH]CC([CH2])CC'),
    E0 = (629.752,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.734022,0.0730277,-5.53513e-05,2.52195e-08,-5.09533e-12,75858.2,33.3043], Tmin=(100,'K'), Tmax=(1123.12,'K')), NASAPolynomial(coeffs=[8.47812,0.0454469,-1.8515e-05,3.35397e-09,-2.28146e-13,74118.7,-4.95565], Tmin=(1123.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(629.752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
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
    label = '[CH2]C(CC)C[C][CH]C(4026)',
    structure = SMILES('[CH2]C(CC)C[C][CH]C'),
    E0 = (605.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.120182,0.0874144,-6.77756e-05,3.14208e-08,-6.44263e-12,73019.5,36.6555], Tmin=(100,'K'), Tmax=(1109.15,'K')), NASAPolynomial(coeffs=[9.44953,0.05377,-2.22762e-05,4.07342e-09,-2.78694e-13,70949.9,-9.31988], Tmin=(1109.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C(CC)C[CH][C]C(4027)',
    structure = SMILES('[CH2]C(CC)C[CH][C]C'),
    E0 = (605.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.222181,0.0870119,-6.93976e-05,3.45433e-08,-7.79019e-12,73015.7,36.1257], Tmin=(100,'K'), Tmax=(1005.6,'K')), NASAPolynomial(coeffs=[7.87524,0.0565697,-2.3988e-05,4.43845e-09,-3.058e-13,71476.5,-0.838501], Tmin=(1005.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(CC)C[CH][CH]C(4028)',
    structure = SMILES('[CH]C(CC)C[CH][CH]C'),
    E0 = (595.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,224.547,823.515,988.218,1152.92,1317.62,1608.64,1730.13],'cm^-1')),
        HinderedRotor(inertia=(0.13045,'amu*angstrom^2'), symmetry=1, barrier=(3.37782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13045,'amu*angstrom^2'), symmetry=1, barrier=(3.37782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13045,'amu*angstrom^2'), symmetry=1, barrier=(3.37782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13045,'amu*angstrom^2'), symmetry=1, barrier=(3.37782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13045,'amu*angstrom^2'), symmetry=1, barrier=(3.37782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13045,'amu*angstrom^2'), symmetry=1, barrier=(3.37782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13045,'amu*angstrom^2'), symmetry=1, barrier=(3.37782,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40278,0.0689683,1.46379e-05,-1.49281e-07,1.36225e-10,71685,34.1933], Tmin=(100,'K'), Tmax=(449.471,'K')), NASAPolynomial(coeffs=[4.08657,0.063459,-2.82971e-05,5.35449e-09,-3.73219e-13,71258.1,21.327], Tmin=(449.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH]C1CC(CC)C1(4029)',
    structure = SMILES('C[CH]C1CC(CC)C1'),
    E0 = (95.2111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64285,0.053606,4.51355e-05,-8.618e-08,3.36361e-11,11589.8,31.5052], Tmin=(100,'K'), Tmax=(1012.21,'K')), NASAPolynomial(coeffs=[14.2653,0.0483837,-1.91603e-05,3.61032e-09,-2.58451e-13,6341.83,-46.6813], Tmin=(1012.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(95.2111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(540.441,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S)"""),
)

species(
    label = 'C=C(CC)CC[CH]C(4030)',
    structure = SMILES('C=C(CC)CC[CH]C'),
    E0 = (68.6635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0745388,0.0817005,-4.77667e-05,1.40468e-08,-1.70762e-12,8402.66,35.0658], Tmin=(100,'K'), Tmax=(1806.58,'K')), NASAPolynomial(coeffs=[15.9908,0.04646,-1.85067e-05,3.24928e-09,-2.13434e-13,2651.84,-51.1345], Tmin=(1806.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.6635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C=CC(C)CC(4031)',
    structure = SMILES('C[CH]C=CC(C)CC'),
    E0 = (10.3874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0561617,0.0717752,-7.84247e-06,-3.22197e-08,1.50391e-11,1404.11,32.0931], Tmin=(100,'K'), Tmax=(1052.76,'K')), NASAPolynomial(coeffs=[14.6734,0.0479564,-1.91004e-05,3.52995e-09,-2.46965e-13,-3431.34,-47.5265], Tmin=(1052.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.3874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(532.126,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S)"""),
)

species(
    label = 'C=C(CC)C[CH][CH]C(4032)',
    structure = SMILES('C=C(CC)C[CH][CH]C'),
    E0 = (263.122,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,2950,3100,1380,975,1025,1650,253.347,866.33,1238.22,1883.01],'cm^-1')),
        HinderedRotor(inertia=(0.125938,'amu*angstrom^2'), symmetry=1, barrier=(3.50194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125938,'amu*angstrom^2'), symmetry=1, barrier=(3.50194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125938,'amu*angstrom^2'), symmetry=1, barrier=(3.50194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125938,'amu*angstrom^2'), symmetry=1, barrier=(3.50194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125938,'amu*angstrom^2'), symmetry=1, barrier=(3.50194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125938,'amu*angstrom^2'), symmetry=1, barrier=(3.50194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.778648,0.0766449,-5.1708e-05,2.3941e-08,-5.87167e-12,31757.2,35.1227], Tmin=(100,'K'), Tmax=(846.509,'K')), NASAPolynomial(coeffs=[3.60249,0.0633012,-2.80627e-05,5.31888e-09,-3.71895e-13,31279.2,21.9699], Tmin=(846.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C=C[CH]C)CC(4033)',
    structure = SMILES('[CH2]C(C=C[CH]C)CC'),
    E0 = (215.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0780507,0.0734547,-2.05211e-05,-1.91977e-08,1.10745e-11,26067.3,33.7936], Tmin=(100,'K'), Tmax=(1035.09,'K')), NASAPolynomial(coeffs=[14.2254,0.0449699,-1.71903e-05,3.0982e-09,-2.13652e-13,21735.7,-41.7235], Tmin=(1035.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(CC)C[CH]C=C(4034)',
    structure = SMILES('[CH2]C=CCC([CH2])CC'),
    E0 = (226.999,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.195719,0.0802966,-3.70637e-05,-4.04333e-09,6.1853e-12,27463,35.4613], Tmin=(100,'K'), Tmax=(1053.62,'K')), NASAPolynomial(coeffs=[15.1901,0.0439633,-1.67691e-05,3.00363e-09,-2.05792e-13,22995.4,-45.3855], Tmin=(1053.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
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
    label = '[CH2][CH]CC=CC(553)',
    structure = SMILES('[CH2][CH]CC=CC'),
    E0 = (323.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,1782.86],'cm^-1')),
        HinderedRotor(inertia=(0.133732,'amu*angstrom^2'), symmetry=1, barrier=(3.07475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133557,'amu*angstrom^2'), symmetry=1, barrier=(3.07073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134216,'amu*angstrom^2'), symmetry=1, barrier=(3.0859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13443,'amu*angstrom^2'), symmetry=1, barrier=(3.09082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80558,0.0490136,-2.40298e-05,5.39389e-09,-4.72015e-13,38969.5,26.5298], Tmin=(100,'K'), Tmax=(2579.8,'K')), NASAPolynomial(coeffs=[18.9745,0.022393,-8.55148e-06,1.39402e-09,-8.43992e-14,30111.1,-72.5714], Tmin=(2579.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C(3856)',
    structure = SMILES('[CH2][CH][CH]C'),
    E0 = (450.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2810.32],'cm^-1')),
        HinderedRotor(inertia=(0.00169393,'amu*angstrom^2'), symmetry=1, barrier=(9.48541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.010072,'amu*angstrom^2'), symmetry=1, barrier=(56.4188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100766,'amu*angstrom^2'), symmetry=1, barrier=(56.4312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79546,0.0335488,-4.51222e-05,5.0661e-08,-2.13657e-11,54172.4,20.3345], Tmin=(100,'K'), Tmax=(854.972,'K')), NASAPolynomial(coeffs=[-1.26875,0.0345785,-1.53757e-05,2.86235e-09,-1.94746e-13,55524.7,43.1492], Tmin=(854.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C[CH][CH]C(555)',
    structure = SMILES('[CH2][CH]C[CH][CH]C'),
    E0 = (596.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,3607.59,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168466,'amu*angstrom^2'), symmetry=1, barrier=(10.7196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79917,0.0586099,-7.80163e-05,8.14675e-08,-3.30993e-11,71839.7,31.8124], Tmin=(100,'K'), Tmax=(852.592,'K')), NASAPolynomial(coeffs=[-2.4828,0.0533525,-2.41736e-05,4.49758e-09,-3.05589e-13,73491.1,57.1902], Tmin=(852.592,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](CC)C[CH][CH]C(4035)',
    structure = SMILES('[CH2][C](CC)C[CH][CH]C'),
    E0 = (537.637,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1622.59,2758.48,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0705152,'amu*angstrom^2'), symmetry=1, barrier=(6.5088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0705152,'amu*angstrom^2'), symmetry=1, barrier=(6.5088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0705152,'amu*angstrom^2'), symmetry=1, barrier=(6.5088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0705152,'amu*angstrom^2'), symmetry=1, barrier=(6.5088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0705152,'amu*angstrom^2'), symmetry=1, barrier=(6.5088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0705152,'amu*angstrom^2'), symmetry=1, barrier=(6.5088,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0705152,'amu*angstrom^2'), symmetry=1, barrier=(6.5088,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.283603,0.0950515,-0.000129493,1.24708e-07,-4.7503e-11,64783.8,39.8565], Tmin=(100,'K'), Tmax=(857.293,'K')), NASAPolynomial(coeffs=[-1.13829,0.0703442,-3.14247e-05,5.80168e-09,-3.92231e-13,66179.3,53.2147], Tmin=(857.293,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.637,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C[CH][CH]C(607)',
    structure = SMILES('[CH2]C([CH2])C[CH][CH]C'),
    E0 = (577.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,375.685,1769.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00303069,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.876276,0.0742873,-7.9907e-05,6.43255e-08,-2.21286e-11,69592,35.4413], Tmin=(100,'K'), Tmax=(870.302,'K')), NASAPolynomial(coeffs=[3.13779,0.05207,-2.12369e-05,3.77352e-09,-2.50661e-13,69646.1,27.4173], Tmin=(870.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442845,0.0857364,-9.15222e-05,7.75899e-08,-2.89865e-11,65880.5,40.7089], Tmin=(100,'K'), Tmax=(812.499,'K')), NASAPolynomial(coeffs=[1.77095,0.0654815,-2.88057e-05,5.35245e-09,-3.6637e-13,66117.4,37.3635], Tmin=(812.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH][CH][CH]C)CC(4036)',
    structure = SMILES('[CH2]C([CH][CH][CH]C)CC'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442845,0.0857364,-9.15222e-05,7.759e-08,-2.89865e-11,65880.5,40.7089], Tmin=(100,'K'), Tmax=(812.499,'K')), NASAPolynomial(coeffs=[1.77095,0.0654815,-2.88057e-05,5.35245e-09,-3.6637e-13,66117.4,37.3635], Tmin=(812.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])C[CH][CH]C(4037)',
    structure = SMILES('[CH2]CC([CH2])C[CH][CH]C'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296787,0.0875608,-8.99057e-05,7.03435e-08,-2.47182e-11,67174.6,40.4657], Tmin=(100,'K'), Tmax=(806.465,'K')), NASAPolynomial(coeffs=[3.66123,0.0622571,-2.68156e-05,4.94189e-09,-3.37175e-13,66912.1,26.695], Tmin=(806.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296787,0.0875608,-8.99057e-05,7.03435e-08,-2.47182e-11,67174.6,40.4657], Tmin=(100,'K'), Tmax=(806.465,'K')), NASAPolynomial(coeffs=[3.66123,0.0622571,-2.68156e-05,4.94189e-09,-3.37175e-13,66912.1,26.695], Tmin=(806.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C[C](C)CC(4039)',
    structure = SMILES('C[CH][CH]C[C](C)CC'),
    E0 = (332.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98806,0.0642371,-2.55546e-05,3.43802e-09,-6.58431e-14,39989.8,28.123], Tmin=(100,'K'), Tmax=(2720.54,'K')), NASAPolynomial(coeffs=[54.396,0.00666161,-3.73934e-06,5.25479e-10,-2.18015e-14,5353.71,-283.59], Tmin=(2720.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.555,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([CH]C[CH]C)CC(4040)',
    structure = SMILES('[CH2]C([CH]C[CH]C)CC'),
    E0 = (352.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,206.219,291.627,1201.19,3327.11],'cm^-1')),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.557675,0.080864,-5.17249e-05,1.96292e-08,-3.53158e-12,42491.2,37.735], Tmin=(100,'K'), Tmax=(1156.65,'K')), NASAPolynomial(coeffs=[6.28324,0.0610634,-2.60465e-05,4.8287e-09,-3.32569e-13,41166.7,9.2792], Tmin=(1156.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C[CH]CC([CH2])CC(367)',
    structure = SMILES('[CH2]C[CH]CC([CH2])CC'),
    E0 = (363.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,243.329,878.828,1260.12,1849.5],'cm^-1')),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375234,0.0840081,-5.82107e-05,2.46379e-08,-4.86908e-12,43787.4,37.4527], Tmin=(100,'K'), Tmax=(1096.67,'K')), NASAPolynomial(coeffs=[6.95045,0.0600255,-2.54078e-05,4.69685e-09,-3.23245e-13,42345.3,5.12431], Tmin=(1096.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]CC(C)[CH]C(1248)',
    structure = SMILES('C[CH][CH]CC(C)[CH]C'),
    E0 = (341.675,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,211.095,896.703,2754.69,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860586,'amu*angstrom^2'), symmetry=1, barrier=(2.56537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79078,0.0691751,-3.0975e-05,5.72233e-09,-3.69667e-13,41152.1,33.9842], Tmin=(100,'K'), Tmax=(2677.54,'K')), NASAPolynomial(coeffs=[41.8126,0.0185235,-7.71796e-06,1.20621e-09,-6.69998e-14,16444.7,-204.632], Tmin=(2677.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH][CH]C(C)CC(4041)',
    structure = SMILES('C[CH][CH][CH]C(C)CC'),
    E0 = (341.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79078,0.0691751,-3.0975e-05,5.72233e-09,-3.69666e-13,41152.1,33.9842], Tmin=(100,'K'), Tmax=(2677.53,'K')), NASAPolynomial(coeffs=[41.8122,0.0185239,-7.71817e-06,1.20624e-09,-6.70023e-14,16445,-204.629], Tmin=(2677.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(341.675,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][C](CC)CC[CH]C(4042)',
    structure = SMILES('[CH2][C](CC)CC[CH]C'),
    E0 = (343.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0723005,0.0940419,-0.000103557,8.58805e-08,-3.10518e-11,41409,38.0539], Tmin=(100,'K'), Tmax=(816.164,'K')), NASAPolynomial(coeffs=[2.98884,0.0666396,-2.91034e-05,5.38617e-09,-3.6772e-13,41369.5,27.2504], Tmin=(816.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH][CH]CC)CC(364)',
    structure = SMILES('[CH2]C([CH][CH]CC)CC'),
    E0 = (352.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,209.12,746.778,915.592,3452.43],'cm^-1')),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542927,0.0816895,-5.70204e-05,2.68127e-08,-6.3819e-12,42492.8,37.6346], Tmin=(100,'K'), Tmax=(894.836,'K')), NASAPolynomial(coeffs=[4.39928,0.0644514,-2.81247e-05,5.28517e-09,-3.67574e-13,41802.6,19.4584], Tmin=(894.836,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CCC([CH2])CC(4043)',
    structure = SMILES('[CH2][CH]CCC([CH2])CC'),
    E0 = (363.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,217.823,811.886,1253.47,1695.09],'cm^-1')),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24835,0.0847187,-5.76949e-05,2.29355e-08,-4.10466e-12,43792.1,38.0701], Tmin=(100,'K'), Tmax=(1223.44,'K')), NASAPolynomial(coeffs=[8.79164,0.0567867,-2.3449e-05,4.27454e-09,-2.91456e-13,41701.7,-4.86923], Tmin=(1223.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)CC[CH]C(4044)',
    structure = SMILES('[CH2]C([CH]C)CC[CH]C'),
    E0 = (352.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,206.219,291.627,1201.19,3327.11],'cm^-1')),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0733722,'amu*angstrom^2'), symmetry=1, barrier=(1.9792,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.557675,0.080864,-5.17249e-05,1.96292e-08,-3.53158e-12,42491.2,37.735], Tmin=(100,'K'), Tmax=(1156.65,'K')), NASAPolynomial(coeffs=[6.28324,0.0610634,-2.60465e-05,4.8287e-09,-3.32569e-13,41166.7,9.2792], Tmin=(1156.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(C)C[CH][CH]C(4009)',
    structure = SMILES('[CH2]CC(C)C[CH][CH]C'),
    E0 = (352.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,206.43,716.6,2010.24,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91664,0.0586225,7.61378e-05,-2.97843e-07,2.71085e-10,42442.7,33.283], Tmin=(100,'K'), Tmax=(393.52,'K')), NASAPolynomial(coeffs=[2.95224,0.0672715,-2.99223e-05,5.66129e-09,-3.94668e-13,42212.7,27.3662], Tmin=(393.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](CC)C[CH]CC(362)',
    structure = SMILES('[CH2][C](CC)C[CH]CC'),
    E0 = (343.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,1461.54,2280.4,3775.88],'cm^-1')),
        HinderedRotor(inertia=(0.187343,'amu*angstrom^2'), symmetry=1, barrier=(4.30739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187343,'amu*angstrom^2'), symmetry=1, barrier=(4.30739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187343,'amu*angstrom^2'), symmetry=1, barrier=(4.30739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187343,'amu*angstrom^2'), symmetry=1, barrier=(4.30739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187343,'amu*angstrom^2'), symmetry=1, barrier=(4.30739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187343,'amu*angstrom^2'), symmetry=1, barrier=(4.30739,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187343,'amu*angstrom^2'), symmetry=1, barrier=(4.30739,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.102408,0.0946116,-0.000109124,9.48241e-08,-3.51661e-11,41408.1,37.7733], Tmin=(100,'K'), Tmax=(823.543,'K')), NASAPolynomial(coeffs=[1.92786,0.0685344,-3.02794e-05,5.62218e-09,-3.83968e-13,41691.1,32.8644], Tmin=(823.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])CC[CH]C(4045)',
    structure = SMILES('[CH2]CC([CH2])CC[CH]C'),
    E0 = (363.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,217.823,811.886,1253.47,1695.09],'cm^-1')),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139611,'amu*angstrom^2'), symmetry=1, barrier=(3.55513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24835,0.0847187,-5.76949e-05,2.29355e-08,-4.10466e-12,43792.1,38.0701], Tmin=(100,'K'), Tmax=(1223.44,'K')), NASAPolynomial(coeffs=[8.79164,0.0567867,-2.3449e-05,4.27454e-09,-2.91456e-13,41701.7,-4.86923], Tmin=(1223.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)C[CH]CC(363)',
    structure = SMILES('[CH2]C([CH]C)C[CH]CC'),
    E0 = (352.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,209.12,746.778,915.592,3452.43],'cm^-1')),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0646036,'amu*angstrom^2'), symmetry=1, barrier=(2.09294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.542927,0.0816895,-5.70204e-05,2.68127e-08,-6.3819e-12,42492.8,37.6346], Tmin=(100,'K'), Tmax=(894.836,'K')), NASAPolynomial(coeffs=[4.39928,0.0644514,-2.81247e-05,5.28517e-09,-3.67574e-13,41802.6,19.4584], Tmin=(894.836,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])C[CH]CC(366)',
    structure = SMILES('[CH2]CC([CH2])C[CH]CC'),
    E0 = (363.015,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,243.329,878.828,1260.12,1849.5],'cm^-1')),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0847426,'amu*angstrom^2'), symmetry=1, barrier=(3.43956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.375234,0.0840081,-5.82107e-05,2.46379e-08,-4.86908e-12,43787.4,37.4527], Tmin=(100,'K'), Tmax=(1096.67,'K')), NASAPolynomial(coeffs=[6.95045,0.0600255,-2.54078e-05,4.69685e-09,-3.23245e-13,42345.3,5.12431], Tmin=(1096.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC(C)CC(4046)',
    structure = SMILES('[CH2][CH][CH]CC(C)CC'),
    E0 = (352.379,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,206.43,716.6,2010.24,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0830779,'amu*angstrom^2'), symmetry=1, barrier=(2.40415,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.205,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91664,0.0586225,7.61378e-05,-2.97843e-07,2.71085e-10,42442.7,33.283], Tmin=(100,'K'), Tmax=(393.52,'K')), NASAPolynomial(coeffs=[2.95224,0.0672715,-2.99223e-05,5.66129e-09,-3.94668e-13,42212.7,27.3662], Tmin=(393.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.379,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(527.969,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    E0 = (352.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (352.215,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (791.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (512.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (512.15,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (517.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (805.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (814.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (821.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (817.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (817.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (807.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (360.499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (415.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (415.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (474.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (435.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (446.565,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (460.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (453.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (414.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (706.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (705.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (749.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (714.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (758.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (758.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (769.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (769.253,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (454.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (489.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (504.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (499.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (470.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (468.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (498.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (504.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (411.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (435.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (470.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (406.353,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (377.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (405.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (405.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['m1_allyl(186)', 'butene1(35)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['[CH2]C(CC)CC=CC(361)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_20] for rate rule [Y_12_20b]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', '[CH2]C(C)C[CH][CH]C(286)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['C[CH][CH]C[CH]CCC(4025)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['C[CH][CH]CC[CH]CC(385)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C([CH]C)C([CH2])CC(344)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', 'C[CH][CH]C[CH]CC(327)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][CH]C(3874)', '[CH2]C([CH2])CC(117)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(17)', '[CH][CH]CC([CH2])CC(579)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2]C(CC)C[C][CH]C(4026)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2]C(CC)C[CH][C]C(4027)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]C(CC)C[CH][CH]C(4028)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['C[CH]C1CC(CC)C1(4029)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['C=C(CC)CC[CH]C(4030)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['C[CH]C=CC(C)CC(4031)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'C=C(CC)C[CH][CH]C(4032)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2]C(C=C[CH]C)CC(4033)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2]C(CC)C[CH]C=C(4034)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C2H5(14)', '[CH2][CH]CC=CC(553)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1020,'cm^3/(mol*s)'), n=2.41, Ea=(27.3634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 418 used for Cds-CsH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH][CH]C(3856)', 'butene1(35)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['m1_allyl(186)', '[CH2][CH]CC(37)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C2H5(14)', '[CH2][CH]C[CH][CH]C(555)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH][CH]C(3856)', '[CH2][CH]CC(37)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2][C](CC)C[CH][CH]C(4035)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH3(17)', '[CH2]C([CH2])C[CH][CH]C(607)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2]C([CH][CH][CH]C)CC(4036)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2]CC([CH2])C[CH][CH]C(4037)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', '[CH2][CH][CH]CC([CH2])CC(4038)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['C[CH][CH]C[C](C)CC(4039)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C([CH]C[CH]C)CC(4040)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C[CH]CC([CH2])CC(367)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['C[CH][CH]CC(C)[CH]C(1248)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['C[CH][CH][CH]C(C)CC(4041)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['[CH2][C](CC)CC[CH]C(4042)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C([CH][CH]CC)CC(364)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.27681e+06,'s^-1'), n=2.16, Ea=(146.492,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH]CCC([CH2])CC(4043)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH]C)CC[CH]C(4044)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(0.329997,'s^-1'), n=3.35583, Ea=(59.1095,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]CC(C)C[CH][CH]C(4009)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C(CC)C[CH][CH]C(365)'],
    products = ['[CH2][C](CC)C[CH]CC(362)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(30253,'s^-1'), n=2.05523, Ea=(118.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R4HJ_1;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]CC([CH2])CC[CH]C(4045)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(416785,'s^-1'), n=1.49392, Ea=(43.3506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_CCC;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([CH]C)C[CH]CC(363)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(0.0564,'s^-1'), n=3.28, Ea=(24.7274,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]CC([CH2])C[CH]CC(366)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R6HJ_4;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][CH][CH]CC(C)CC(4046)'],
    products = ['[CH2]C(CC)C[CH][CH]C(365)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1079',
    isomers = [
        '[CH2]C(CC)C[CH][CH]C(365)',
    ],
    reactants = [
        ('m1_allyl(186)', 'butene1(35)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1079',
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


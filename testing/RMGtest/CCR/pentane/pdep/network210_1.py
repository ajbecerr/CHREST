species(
    label = '[CH2]C([CH2])C([CH2])C(133)',
    structure = SMILES('[CH2]C([CH2])C([CH2])C'),
    E0 = (411.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2078.21],'cm^-1')),
        HinderedRotor(inertia=(0.0541491,'amu*angstrom^2'), symmetry=1, barrier=(2.96105,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138827,'amu*angstrom^2'), symmetry=1, barrier=(5.70692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00419206,'amu*angstrom^2'), symmetry=1, barrier=(0.168624,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.61498,'amu*angstrom^2'), symmetry=1, barrier=(70.1476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0229134,'amu*angstrom^2'), symmetry=1, barrier=(70.1842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10958,0.0578175,-2.40705e-05,-1.05746e-08,9.8225e-12,49660.9,28.5637], Tmin=(100,'K'), Tmax=(880.344,'K')), NASAPolynomial(coeffs=[10.6936,0.0319281,-1.00439e-05,1.58686e-09,-1.01178e-13,47289.3,-20.3384], Tmin=(880.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = 'C3H6(27)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36193e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35107e-06,1.16619e-09,-8.2762e-14,-487.139,-4.54469], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]CC([CH2])[CH2](112)',
    structure = SMILES('[CH2]CC([CH2])[CH2]'),
    E0 = (441.632,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,484.822],'cm^-1')),
        HinderedRotor(inertia=(3.48101,'amu*angstrom^2'), symmetry=1, barrier=(80.0352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0138261,'amu*angstrom^2'), symmetry=1, barrier=(79.9434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0137509,'amu*angstrom^2'), symmetry=1, barrier=(79.9414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252917,'amu*angstrom^2'), symmetry=1, barrier=(14.6729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3074.66,'J/mol'), sigma=(5.76682,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=480.25 K, Pc=36.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76406,0.0446331,-1.86055e-05,-5.60437e-09,5.64408e-12,53200.8,24.2516], Tmin=(100,'K'), Tmax=(919.276,'K')), NASAPolynomial(coeffs=[8.54585,0.026962,-9.0878e-06,1.50178e-09,-9.81125e-14,51453.7,-10.6164], Tmin=(919.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(441.632,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C[CH]C(148)',
    structure = SMILES('[CH2]C([CH2])C[CH]C'),
    E0 = (407.052,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1669.12],'cm^-1')),
        HinderedRotor(inertia=(0.0329702,'amu*angstrom^2'), symmetry=1, barrier=(65.1865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179265,'amu*angstrom^2'), symmetry=1, barrier=(4.12164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180839,'amu*angstrom^2'), symmetry=1, barrier=(4.15785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00210479,'amu*angstrom^2'), symmetry=1, barrier=(4.15814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179369,'amu*angstrom^2'), symmetry=1, barrier=(4.12405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3039,0.0586265,-4.06869e-05,1.80741e-08,-3.65723e-12,49054.6,29.096], Tmin=(100,'K'), Tmax=(1112.42,'K')), NASAPolynomial(coeffs=[6.64732,0.0394129,-1.47791e-05,2.54772e-09,-1.67926e-13,47865.8,2.74781], Tmin=(1112.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH]CC(513)',
    structure = SMILES('[CH2]C([CH2])[CH]CC'),
    E0 = (407.148,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17559,0.0586211,-3.75045e-05,1.37025e-08,-2.17533e-12,49073,29.3758], Tmin=(100,'K'), Tmax=(1407.27,'K')), NASAPolynomial(coeffs=[9.11307,0.0360597,-1.34564e-05,2.31018e-09,-1.51484e-13,46838.9,-11.6298], Tmin=(1407.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.148,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC([CH2])C(149)',
    structure = SMILES('[CH2][CH]CC([CH2])C'),
    E0 = (407.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,2376.69],'cm^-1')),
        HinderedRotor(inertia=(0.0188907,'amu*angstrom^2'), symmetry=1, barrier=(75.7394,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0020477,'amu*angstrom^2'), symmetry=1, barrier=(8.20866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14901,'amu*angstrom^2'), symmetry=1, barrier=(3.42603,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0299108,'amu*angstrom^2'), symmetry=1, barrier=(8.2028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278496,'amu*angstrom^2'), symmetry=1, barrier=(75.7433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45716,0.0561691,-3.36023e-05,1.14364e-08,-1.75321e-12,49067.8,29.2376], Tmin=(100,'K'), Tmax=(1381.33,'K')), NASAPolynomial(coeffs=[7.16532,0.0396397,-1.56529e-05,2.77354e-09,-1.85368e-13,47490.8,-0.145002], Tmin=(1381.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2]C([CH2])[CH]C(500)',
    structure = SMILES('[CH2]C([CH2])[CH]C'),
    E0 = (430.928,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1611.15],'cm^-1')),
        HinderedRotor(inertia=(0.00291815,'amu*angstrom^2'), symmetry=1, barrier=(5.37538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87676,'amu*angstrom^2'), symmetry=1, barrier=(66.1424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00291914,'amu*angstrom^2'), symmetry=1, barrier=(5.37711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.233802,'amu*angstrom^2'), symmetry=1, barrier=(5.37557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73582,0.0447701,-2.6567e-05,8.99144e-09,-1.3162e-12,51914.3,25.1256], Tmin=(100,'K'), Tmax=(1524.32,'K')), NASAPolynomial(coeffs=[8.30204,0.0275397,-9.61162e-06,1.576e-09,-1.00019e-13,49912.5,-9.32069], Tmin=(1524.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(111)',
    structure = SMILES('[CH2][CH]C([CH2])C'),
    E0 = (431.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1088.38],'cm^-1')),
        HinderedRotor(inertia=(0.00507348,'amu*angstrom^2'), symmetry=1, barrier=(4.26448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00507018,'amu*angstrom^2'), symmetry=1, barrier=(4.2639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185504,'amu*angstrom^2'), symmetry=1, barrier=(4.2651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185366,'amu*angstrom^2'), symmetry=1, barrier=(4.26194,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72887,0.0438928,-2.37827e-05,6.55025e-09,-7.41727e-13,51935.3,25.8655], Tmin=(100,'K'), Tmax=(1961.54,'K')), NASAPolynomial(coeffs=[11.3781,0.0242161,-8.73612e-06,1.43645e-09,-8.99776e-14,48149.8,-27.1877], Tmin=(1961.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH]C(C)C([CH2])[CH2](514)',
    structure = SMILES('[CH]C(C)C([CH2])[CH2]'),
    E0 = (655.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200.05,800.034,1200.15,1600.27],'cm^-1')),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995419,0.0597654,-3.34411e-05,2.51573e-10,5.18413e-12,78907.5,28.1922], Tmin=(100,'K'), Tmax=(930.231,'K')), NASAPolynomial(coeffs=[11.779,0.029247,-9.78972e-06,1.61911e-09,-1.06299e-13,76215.4,-26.7386], Tmin=(930.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])C([CH2])C(515)',
    structure = SMILES('[CH]C([CH2])C([CH2])C'),
    E0 = (655.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200.05,800.034,1200.15,1600.27],'cm^-1')),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156034,'amu*angstrom^2'), symmetry=1, barrier=(3.5887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995419,0.0597654,-3.34411e-05,2.51573e-10,5.18413e-12,78907.5,28.8853], Tmin=(100,'K'), Tmax=(930.231,'K')), NASAPolynomial(coeffs=[11.779,0.029247,-9.78972e-06,1.61911e-09,-1.06299e-13,76215.4,-26.0454], Tmin=(930.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C(C)C1CC1(516)',
    structure = SMILES('[CH2]C(C)C1CC1'),
    E0 = (157.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65005,0.0339174,5.25521e-05,-9.13999e-08,3.75506e-11,19071.5,23.3028], Tmin=(100,'K'), Tmax=(947.151,'K')), NASAPolynomial(coeffs=[13.6079,0.0280418,-8.81492e-06,1.5379e-09,-1.10133e-13,14804.6,-44.3044], Tmin=(947.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CCC1C(517)',
    structure = SMILES('[CH2]C1CCC1C'),
    E0 = (153.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93953,0.0248738,7.89213e-05,-1.17333e-07,4.62036e-11,18531.9,21.5848], Tmin=(100,'K'), Tmax=(950.37,'K')), NASAPolynomial(coeffs=[13.1726,0.0291785,-9.28881e-06,1.65648e-09,-1.20773e-13,14067.3,-44.2922], Tmin=(950.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.312,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C(=C)C(131)',
    structure = SMILES('[CH2]C(C)C(=C)C'),
    E0 = (121.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.701264,'amu*angstrom^2'), symmetry=1, barrier=(16.1234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701245,'amu*angstrom^2'), symmetry=1, barrier=(16.123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0477163,'amu*angstrom^2'), symmetry=1, barrier=(1.09709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0095191,'amu*angstrom^2'), symmetry=1, barrier=(83.1217,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11846,0.054805,-1.74788e-05,-1.1889e-08,7.64616e-12,14776.9,24.9582], Tmin=(100,'K'), Tmax=(1010.97,'K')), NASAPolynomial(coeffs=[11.0134,0.0334956,-1.23326e-05,2.17345e-09,-1.47981e-13,11864.5,-27.3961], Tmin=(1010.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C)C(C)C(518)',
    structure = SMILES('[CH2]C(=C)C(C)C'),
    E0 = (68.3542,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12165,0.0503641,4.53197e-06,-3.90431e-08,1.78383e-11,8336.18,22.0598], Tmin=(100,'K'), Tmax=(993.755,'K')), NASAPolynomial(coeffs=[13.3542,0.0310888,-1.15994e-05,2.11892e-09,-1.49587e-13,4425.48,-44.3224], Tmin=(993.755,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.3542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C(=C)C(519)',
    structure = SMILES('[CH2]C([CH2])C(=C)C'),
    E0 = (327.02,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,1775.29],'cm^-1')),
        HinderedRotor(inertia=(0.001614,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106581,'amu*angstrom^2'), symmetry=1, barrier=(7.82154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106087,'amu*angstrom^2'), symmetry=1, barrier=(7.83331,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.948151,'amu*angstrom^2'), symmetry=1, barrier=(69.8963,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18306,0.0559407,-2.80598e-05,-1.88254e-09,5.10389e-12,39438.2,25.815], Tmin=(100,'K'), Tmax=(943.753,'K')), NASAPolynomial(coeffs=[10.5513,0.0305482,-1.04512e-05,1.74949e-09,-1.15379e-13,37032.5,-22.216], Tmin=(943.753,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.02,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])C(520)',
    structure = SMILES('[CH2]C(=C)C([CH2])C'),
    E0 = (273.437,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.00352323,'amu*angstrom^2'), symmetry=1, barrier=(2.31509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926497,'amu*angstrom^2'), symmetry=1, barrier=(21.302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156157,'amu*angstrom^2'), symmetry=1, barrier=(11.5133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0857942,'amu*angstrom^2'), symmetry=1, barrier=(75.9546,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1564,0.0518607,-7.35032e-06,-2.73003e-08,1.4536e-11,32998.8,24.4093], Tmin=(100,'K'), Tmax=(961.973,'K')), NASAPolynomial(coeffs=[13.0073,0.0279449,-9.60444e-06,1.6681e-09,-1.14754e-13,29545.3,-38.4041], Tmin=(961.973,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C)C=C(109)',
    structure = SMILES('[CH2]C(C)C=C'),
    E0 = (156.808,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,663.668],'cm^-1')),
        HinderedRotor(inertia=(0.00204431,'amu*angstrom^2'), symmetry=1, barrier=(15.8768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.690941,'amu*angstrom^2'), symmetry=1, barrier=(15.8861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.57515,'amu*angstrom^2'), symmetry=1, barrier=(82.1997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99035,0.0344738,1.41962e-05,-4.03588e-08,1.76162e-11,18940.8,21.072], Tmin=(100,'K'), Tmax=(956.909,'K')), NASAPolynomial(coeffs=[9.88824,0.0252392,-8.60353e-06,1.49487e-09,-1.03137e-13,16340.6,-22.3716], Tmin=(956.909,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
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
    label = '[CH2]C([CH2])C=C(496)',
    structure = SMILES('[CH2]C([CH2])C=C'),
    E0 = (361.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,444.303],'cm^-1')),
        HinderedRotor(inertia=(0.0434666,'amu*angstrom^2'), symmetry=1, barrier=(6.11473,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.043652,'amu*angstrom^2'), symmetry=1, barrier=(6.1174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0394383,'amu*angstrom^2'), symmetry=1, barrier=(69.5439,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.03654,0.0357961,3.13511e-06,-3.00403e-08,1.51068e-11,43603,21.9965], Tmin=(100,'K'), Tmax=(906.7,'K')), NASAPolynomial(coeffs=[9.64737,0.0219243,-6.51387e-06,1.02236e-09,-6.65462e-14,41412.9,-18.4424], Tmin=(906.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'C3H6(T)(28)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238388,'amu*angstrom^2'), symmetry=1, barrier=(5.48101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0090964,'amu*angstrom^2'), symmetry=1, barrier=(22.1004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93779,0.019099,4.26879e-06,-1.44878e-08,5.7496e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.8,'K')), NASAPolynomial(coeffs=[5.93905,0.0171893,-6.69156e-06,1.21547e-09,-8.39803e-14,33151.2,-4.14862], Tmin=(1046.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH2](499)',
    structure = SMILES('[CH2][CH]C([CH2])[CH2]'),
    E0 = (636.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1200.87],'cm^-1')),
        HinderedRotor(inertia=(0.112344,'amu*angstrom^2'), symmetry=1, barrier=(2.58301,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112403,'amu*angstrom^2'), symmetry=1, barrier=(2.58438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252635,'amu*angstrom^2'), symmetry=1, barrier=(2.58576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00252511,'amu*angstrom^2'), symmetry=1, barrier=(2.5843,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95305,0.0434597,-3.02543e-05,1.30344e-08,-2.49244e-12,76589.1,26.1277], Tmin=(100,'K'), Tmax=(1189.32,'K')), NASAPolynomial(coeffs=[6.71646,0.027439,-1.00484e-05,1.70807e-09,-1.11583e-13,75456,2.32115], Tmin=(1189.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C)C([CH2])[CH2](521)',
    structure = SMILES('[CH2][C](C)C([CH2])[CH2]'),
    E0 = (597.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1862.37],'cm^-1')),
        HinderedRotor(inertia=(0.00179192,'amu*angstrom^2'), symmetry=1, barrier=(4.41175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191665,'amu*angstrom^2'), symmetry=1, barrier=(4.40675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00179277,'amu*angstrom^2'), symmetry=1, barrier=(4.41185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.79903,'amu*angstrom^2'), symmetry=1, barrier=(64.3551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192148,'amu*angstrom^2'), symmetry=1, barrier=(4.41786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15771,0.0657478,-7.30623e-05,5.47332e-08,-1.69059e-11,71951.9,29.522], Tmin=(100,'K'), Tmax=(936.636,'K')), NASAPolynomial(coeffs=[5.51368,0.038027,-1.40655e-05,2.34753e-09,-1.49281e-13,71535.8,10.9272], Tmin=(936.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C([CH2])C(522)',
    structure = SMILES('[CH2][C]([CH2])C([CH2])C'),
    E0 = (597.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1862.37],'cm^-1')),
        HinderedRotor(inertia=(0.00179192,'amu*angstrom^2'), symmetry=1, barrier=(4.41175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191665,'amu*angstrom^2'), symmetry=1, barrier=(4.40675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00179277,'amu*angstrom^2'), symmetry=1, barrier=(4.41185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.79903,'amu*angstrom^2'), symmetry=1, barrier=(64.3551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192148,'amu*angstrom^2'), symmetry=1, barrier=(4.41786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15771,0.0657478,-7.30623e-05,5.47332e-08,-1.69059e-11,71951.9,29.522], Tmin=(100,'K'), Tmax=(936.636,'K')), NASAPolynomial(coeffs=[5.51368,0.038027,-1.40655e-05,2.34753e-09,-1.49281e-13,71535.8,10.9272], Tmin=(936.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C([CH2])[CH2](485)',
    structure = SMILES('[CH2]C([CH2])C([CH2])[CH2]'),
    E0 = (617.076,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,711.401],'cm^-1')),
        HinderedRotor(inertia=(0.0740442,'amu*angstrom^2'), symmetry=1, barrier=(1.70242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00510071,'amu*angstrom^2'), symmetry=1, barrier=(1.83404,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0708634,'amu*angstrom^2'), symmetry=1, barrier=(1.62929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0768939,'amu*angstrom^2'), symmetry=1, barrier=(1.76794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00487409,'amu*angstrom^2'), symmetry=1, barrier=(1.74828,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.6689,0.0650717,-5.68349e-05,2.94125e-08,-6.08673e-12,74344,30.5296], Tmin=(100,'K'), Tmax=(1324.25,'K')), NASAPolynomial(coeffs=[11.0931,0.0274373,-7.24249e-06,9.40708e-10,-4.98543e-14,72122.2,-20.6535], Tmin=(1324.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[C](C)C(523)',
    structure = SMILES('[CH2]C([CH2])[C](C)C'),
    E0 = (392.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2431,0.0627963,-5.59495e-05,3.60841e-08,-1.06395e-11,47284.1,26.7444], Tmin=(100,'K'), Tmax=(840.57,'K')), NASAPolynomial(coeffs=[5.28963,0.0421513,-1.66299e-05,2.9336e-09,-1.95379e-13,46652.9,8.21683], Tmin=(840.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C)C([CH2])C(132)',
    structure = SMILES('[CH2][C](C)C([CH2])C'),
    E0 = (392.334,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2510.14],'cm^-1')),
        HinderedRotor(inertia=(0.289067,'amu*angstrom^2'), symmetry=1, barrier=(6.64621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.05837,'amu*angstrom^2'), symmetry=1, barrier=(70.3179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.015728,'amu*angstrom^2'), symmetry=1, barrier=(70.3291,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0014863,'amu*angstrom^2'), symmetry=1, barrier=(6.64357,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289138,'amu*angstrom^2'), symmetry=1, barrier=(6.64786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2431,0.0627963,-5.59495e-05,3.60841e-08,-1.06395e-11,47284.1,28.1307], Tmin=(100,'K'), Tmax=(840.57,'K')), NASAPolynomial(coeffs=[5.28963,0.0421513,-1.66299e-05,2.9336e-09,-1.95379e-13,46652.9,9.60313], Tmin=(840.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C(C)C(524)',
    structure = SMILES('[CH2][C]([CH2])C(C)C'),
    E0 = (392.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2431,0.0627963,-5.59495e-05,3.60841e-08,-1.06395e-11,47284.1,26.7444], Tmin=(100,'K'), Tmax=(840.57,'K')), NASAPolynomial(coeffs=[5.28963,0.0421513,-1.66299e-05,2.9336e-09,-1.95379e-13,46652.9,8.21683], Tmin=(840.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (411.994,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (860.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (571.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (606.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (569.312,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (868.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (868.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (866.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (866.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (420.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (420.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (434.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (475.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (538.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (496.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (539.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (510.456,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (528.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (462.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (769.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (772.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (809.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (809.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (828.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (514.673,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (553.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (529.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['allyl(82)', 'C3H6(27)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2]CC([CH2])[CH2](112)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2]C([CH2])C[CH]C(148)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2]C([CH2])[CH]CC(513)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2][CH]CC([CH2])C(149)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2]C([CH2])[CH]C(500)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[CH2][CH]C([CH2])C(111)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C(C)C([CH2])[CH2](514)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([CH2])C([CH2])C(515)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2]C(C)C1CC1(516)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2]C1CCC1C(517)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2]C(C)C(=C)C(131)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2]C(=C)C(C)C(518)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH2]C([CH2])C(=C)C(519)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2]C(=C)C([CH2])C(520)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['CH2(T)(20)', '[CH2]C(C)C=C(109)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;Y_1centerbirad] for rate rule [Cds-CsH_Cds-HH;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH][CH2](497)', 'C3H6(27)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(17)', '[CH2]C([CH2])C=C(496)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(10000,'cm^3/(mol*s)'), n=2.41, Ea=(29.7482,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 417 used for Cds-CsH_Cds-HH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['allyl(82)', 'C3H6(T)(28)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.0018001,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH][CH2](497)', 'C3H6(T)(28)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH3(17)', '[CH2][CH]C([CH2])[CH2](499)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.49263e+08,'m^3/(mol*s)'), n=-0.519467, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00186367249776, var=0.658318197044, Tref=1000.0, N=7, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O
    Total Standard Deviation in ln(k): 1.6312606891
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][C](C)C([CH2])[CH2](521)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][C]([CH2])C([CH2])C(522)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C([CH2])C([CH2])[CH2](485)'],
    products = ['[CH2]C([CH2])C([CH2])C(133)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.67e+08,'m^3/(mol*s)'), n=0.115342, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_N-Sp-4R!H=3C_Sp-4R!H-3C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2]C([CH2])[C](C)C(523)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2][C](C)C([CH2])C(132)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C([CH2])C([CH2])C(133)'],
    products = ['[CH2][C]([CH2])C(C)C(524)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #210',
    isomers = [
        '[CH2]C([CH2])C([CH2])C(133)',
    ],
    reactants = [
        ('allyl(82)', 'C3H6(27)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #210',
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


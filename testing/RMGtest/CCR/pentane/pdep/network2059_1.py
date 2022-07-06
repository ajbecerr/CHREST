species(
    label = '[CH2][CH][CH]CC([CH2])[CH2](4136)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])[CH2]'),
    E0 = (782.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1213.89,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0145796,'amu*angstrom^2'), symmetry=1, barrier=(11.2803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.898412,0.0752073,-9.10451e-05,7.75073e-08,-2.6973e-11,94275.3,37.1477], Tmin=(100,'K'), Tmax=(883.59,'K')), NASAPolynomial(coeffs=[2.98448,0.0496822,-2.0413e-05,3.61754e-09,-2.38926e-13,94534.4,30.8942], Tmin=(883.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(782.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]CC[CH][CH2](4163)',
    structure = SMILES('[CH2][CH][CH]CC[CH][CH2]'),
    E0 = (778.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,3516.19,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0114336,'amu*angstrom^2'), symmetry=1, barrier=(13.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114336,'amu*angstrom^2'), symmetry=1, barrier=(13.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114336,'amu*angstrom^2'), symmetry=1, barrier=(13.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114336,'amu*angstrom^2'), symmetry=1, barrier=(13.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114336,'amu*angstrom^2'), symmetry=1, barrier=(13.5764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0114336,'amu*angstrom^2'), symmetry=1, barrier=(13.5764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.80518,0.046418,-1.65621e-05,1.42179e-09,1.06398e-13,93558.7,28.3883], Tmin=(100,'K'), Tmax=(2731,'K')), NASAPolynomial(coeffs=[48.0514,-0.000118739,-1.03615e-06,8.12319e-11,4.88421e-15,62578.6,-242], Tmin=(2731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(778.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])[CH2](4093)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])[CH2]'),
    E0 = (791.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,935.446,4000],'cm^-1')),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016801,'amu*angstrom^2'), symmetry=1, barrier=(7.40964,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.594522,0.0733338,-6.44468e-05,3.47395e-08,-7.93788e-12,95301.3,36.9958], Tmin=(100,'K'), Tmax=(1040.15,'K')), NASAPolynomial(coeffs=[9.70761,0.0382888,-1.39089e-05,2.34833e-09,-1.52743e-13,93405.4,-7.32843], Tmin=(1040.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(791.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C[CH][CH2](653)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH2]'),
    E0 = (801.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2055.29,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000899643,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82635,0.0594615,-8.88659e-05,9.41844e-08,-3.76943e-11,96522.8,33.5013], Tmin=(100,'K'), Tmax=(862.635,'K')), NASAPolynomial(coeffs=[-2.63202,0.0509582,-2.33462e-05,4.3408e-09,-2.93791e-13,98377.6,60.6438], Tmin=(862.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(801.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2]C([CH2])[CH2](489)',
    structure = SMILES('[CH2]C([CH2])[CH2]'),
    E0 = (463.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00218901,'amu*angstrom^2'), symmetry=1, barrier=(6.67823,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0226841,'amu*angstrom^2'), symmetry=1, barrier=(69.1489,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.290286,'amu*angstrom^2'), symmetry=1, barrier=(6.67424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89364,0.0397791,-3.48191e-05,1.81739e-08,-3.6853e-12,55795.5,20.1671], Tmin=(100,'K'), Tmax=(1417.63,'K')), NASAPolynomial(coeffs=[8.35096,0.0156406,-3.01573e-06,2.72841e-10,-9.09425e-15,54559.4,-11.1417], Tmin=(1417.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2][CH][C]CC([CH2])[CH2](9231)',
    structure = SMILES('[CH2][CH][C]CC([CH2])[CH2]'),
    E0 = (1036.73,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3025,407.5,1350,352.5,202.424,801.617,1207.27,1612.93],'cm^-1')),
        HinderedRotor(inertia=(0.153398,'amu*angstrom^2'), symmetry=1, barrier=(3.58418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153398,'amu*angstrom^2'), symmetry=1, barrier=(3.58418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153398,'amu*angstrom^2'), symmetry=1, barrier=(3.58418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153398,'amu*angstrom^2'), symmetry=1, barrier=(3.58418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153398,'amu*angstrom^2'), symmetry=1, barrier=(3.58418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153398,'amu*angstrom^2'), symmetry=1, barrier=(3.58418,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.574837,0.077983,-8.72514e-05,6.17718e-08,-1.83031e-11,124811,35.6515], Tmin=(100,'K'), Tmax=(875.576,'K')), NASAPolynomial(coeffs=[8.27638,0.0394277,-1.54243e-05,2.68464e-09,-1.76502e-13,123592,0.257293], Tmin=(875.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1036.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]CC([CH2])[CH2](9232)',
    structure = SMILES('[CH2][C][CH]CC([CH2])[CH2]'),
    E0 = (1036.75,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,245.162,856.429,1335.49,1840.76],'cm^-1')),
        HinderedRotor(inertia=(0.11375,'amu*angstrom^2'), symmetry=1, barrier=(3.46165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11375,'amu*angstrom^2'), symmetry=1, barrier=(3.46165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11375,'amu*angstrom^2'), symmetry=1, barrier=(3.46165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11375,'amu*angstrom^2'), symmetry=1, barrier=(3.46165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11375,'amu*angstrom^2'), symmetry=1, barrier=(3.46165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11375,'amu*angstrom^2'), symmetry=1, barrier=(3.46165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.620212,0.0783569,-9.20573e-05,6.96108e-08,-2.18912e-11,124810,35.3171], Tmin=(100,'K'), Tmax=(873.094,'K')), NASAPolynomial(coeffs=[7.17415,0.0413963,-1.66445e-05,2.93135e-09,-1.93655e-13,123929,6.10121], Tmin=(873.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1036.75,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])C[CH][CH][CH2](9233)',
    structure = SMILES('[CH]C([CH2])C[CH][CH][CH2]'),
    E0 = (1026.11,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,184.568,406.574,1301.58,2145.54,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.103745,'amu*angstrom^2'), symmetry=1, barrier=(2.39221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103745,'amu*angstrom^2'), symmetry=1, barrier=(2.39221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103745,'amu*angstrom^2'), symmetry=1, barrier=(2.39221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103745,'amu*angstrom^2'), symmetry=1, barrier=(2.39221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103745,'amu*angstrom^2'), symmetry=1, barrier=(2.39221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103745,'amu*angstrom^2'), symmetry=1, barrier=(2.39221,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.85956,0.0762193,-9.69245e-05,8.3507e-08,-2.94272e-11,123519,37.2024], Tmin=(100,'K'), Tmax=(858.139,'K')), NASAPolynomial(coeffs=[3.81611,0.0474448,-2.042e-05,3.71239e-09,-2.49297e-13,123563,26.6076], Tmin=(858.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1026.11,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH][CH]CC([CH2])[CH2](9234)',
    structure = SMILES('[CH][CH][CH]CC([CH2])[CH2]'),
    E0 = (1025.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,186.52,646.734,991.963,2236.15,3093.44,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0716296,'amu*angstrom^2'), symmetry=1, barrier=(1.66445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716296,'amu*angstrom^2'), symmetry=1, barrier=(1.66445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716296,'amu*angstrom^2'), symmetry=1, barrier=(1.66445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716296,'amu*angstrom^2'), symmetry=1, barrier=(1.66445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716296,'amu*angstrom^2'), symmetry=1, barrier=(1.66445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0716296,'amu*angstrom^2'), symmetry=1, barrier=(1.66445,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.795835,0.077683,-0.000100644,8.56007e-08,-2.91762e-11,123501,36.7333], Tmin=(100,'K'), Tmax=(893.346,'K')), NASAPolynomial(coeffs=[4.13935,0.0458079,-1.87399e-05,3.29714e-09,-2.162e-13,123578,24.7552], Tmin=(893.346,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1025.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH][CH]CC1CC1(9235)',
    structure = SMILES('[CH2][CH][CH]CC1CC1'),
    E0 = (532.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30237,0.0525423,-1.68519e-05,-3.15506e-09,2.10547e-12,64094.7,32.4002], Tmin=(100,'K'), Tmax=(1349.47,'K')), NASAPolynomial(coeffs=[9.30347,0.0404055,-1.62326e-05,2.89762e-09,-1.93826e-13,60880.9,-12.505], Tmin=(1349.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C1CC([CH2])C1(5228)',
    structure = SMILES('[CH2][CH]C1CC([CH2])C1'),
    E0 = (529.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44015,0.0401354,3.98002e-05,-7.42515e-08,2.97041e-11,63769,29.9649], Tmin=(100,'K'), Tmax=(982.689,'K')), NASAPolynomial(coeffs=[12.2781,0.0353991,-1.30795e-05,2.40153e-09,-1.71048e-13,59737.6,-31.8069], Tmin=(982.689,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(529.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC[C]([CH2])C(5011)',
    structure = SMILES('[CH2]C=CC[C]([CH2])C'),
    E0 = (432.855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.102837,'amu*angstrom^2'), symmetry=1, barrier=(2.36441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0995862,'amu*angstrom^2'), symmetry=1, barrier=(2.28968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0984102,'amu*angstrom^2'), symmetry=1, barrier=(2.26264,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0982723,'amu*angstrom^2'), symmetry=1, barrier=(2.25947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00727311,'amu*angstrom^2'), symmetry=1, barrier=(2.31243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08981,0.0663958,-4.71361e-05,1.99482e-08,-3.82871e-12,52163.1,29.7387], Tmin=(100,'K'), Tmax=(1149.28,'K')), NASAPolynomial(coeffs=[7.35999,0.0445731,-1.86543e-05,3.42689e-09,-2.3492e-13,50721.8,-1.38374], Tmin=(1149.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]([CH2])CCC=C(5242)',
    structure = SMILES('[CH2][C]([CH2])CCC=C'),
    E0 = (502.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,249.82,893.784,1978.92],'cm^-1')),
        HinderedRotor(inertia=(0.0947348,'amu*angstrom^2'), symmetry=1, barrier=(3.39845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947348,'amu*angstrom^2'), symmetry=1, barrier=(3.39845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947348,'amu*angstrom^2'), symmetry=1, barrier=(3.39845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947348,'amu*angstrom^2'), symmetry=1, barrier=(3.39845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0947348,'amu*angstrom^2'), symmetry=1, barrier=(3.39845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.825631,0.0731084,-7.31838e-05,5.09823e-08,-1.57304e-11,60491.5,31.9353], Tmin=(100,'K'), Tmax=(821.812,'K')), NASAPolynomial(coeffs=[6.05842,0.0453994,-1.85209e-05,3.32291e-09,-2.23386e-13,59707.1,8.17712], Tmin=(821.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])C(5012)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])C'),
    E0 = (388.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1451.4,2646.44],'cm^-1')),
        HinderedRotor(inertia=(0.828594,'amu*angstrom^2'), symmetry=1, barrier=(19.051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0371505,'amu*angstrom^2'), symmetry=1, barrier=(19.0516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.48473,'amu*angstrom^2'), symmetry=1, barrier=(80.1207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.828569,'amu*angstrom^2'), symmetry=1, barrier=(19.0504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0535999,'amu*angstrom^2'), symmetry=1, barrier=(80.1131,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.787236,0.0586763,-8.72898e-06,-2.75874e-08,1.41857e-11,46857.4,28.5804], Tmin=(100,'K'), Tmax=(994.507,'K')), NASAPolynomial(coeffs=[13.5045,0.0344226,-1.27148e-05,2.27869e-09,-1.58143e-13,42997.8,-39.3902], Tmin=(994.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(388.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[C]([CH2])[CH2](5231)',
    structure = SMILES('[CH2]C=CC[C]([CH2])[CH2]'),
    E0 = (637.937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,1638.84],'cm^-1')),
        HinderedRotor(inertia=(0.0015895,'amu*angstrom^2'), symmetry=1, barrier=(3.02811,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131541,'amu*angstrom^2'), symmetry=1, barrier=(3.02439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131866,'amu*angstrom^2'), symmetry=1, barrier=(3.03185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.73634,'amu*angstrom^2'), symmetry=1, barrier=(62.9138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13172,'amu*angstrom^2'), symmetry=1, barrier=(3.0285,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.903897,0.0704266,-6.76191e-05,4.26068e-08,-1.17423e-11,76835.5,31.4986], Tmin=(100,'K'), Tmax=(858.754,'K')), NASAPolynomial(coeffs=[7.18433,0.0411753,-1.65296e-05,2.94847e-09,-1.97994e-13,75756.7,2.15494], Tmin=(858.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C([CH2])[CH2](5232)',
    structure = SMILES('[CH2]C=C[CH]C([CH2])[CH2]'),
    E0 = (593.627,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,984.58],'cm^-1')),
        HinderedRotor(inertia=(0.289119,'amu*angstrom^2'), symmetry=1, barrier=(81.6058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00300946,'amu*angstrom^2'), symmetry=1, barrier=(20.7869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118639,'amu*angstrom^2'), symmetry=1, barrier=(81.6085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0118144,'amu*angstrom^2'), symmetry=1, barrier=(81.6044,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.54948,'amu*angstrom^2'), symmetry=1, barrier=(81.6095,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.831027,0.0600621,-2.02034e-05,-1.64021e-08,1.11339e-11,71519.6,29.5115], Tmin=(100,'K'), Tmax=(956.339,'K')), NASAPolynomial(coeffs=[13.1306,0.0313254,-1.0747e-05,1.83436e-09,-1.23852e-13,68128.7,-34.707], Tmin=(956.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.627,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Isobutyl) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC=C[CH2](652)',
    structure = SMILES('[CH2][CH]CC=C[CH2]'),
    E0 = (474.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,348.938,2083.29],'cm^-1')),
        HinderedRotor(inertia=(0.107987,'amu*angstrom^2'), symmetry=1, barrier=(9.35741,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0235326,'amu*angstrom^2'), symmetry=1, barrier=(56.8317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108208,'amu*angstrom^2'), symmetry=1, barrier=(9.35905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0127607,'amu*angstrom^2'), symmetry=1, barrier=(30.8466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44643,0.0505193,-2.85643e-05,7.98419e-09,-9.07072e-13,57209.5,28.1169], Tmin=(100,'K'), Tmax=(1965.55,'K')), NASAPolynomial(coeffs=[13.3496,0.0262965,-1.00793e-05,1.71471e-09,-1.09676e-13,52530.1,-37.3531], Tmin=(1965.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]C[C]([CH2])[CH2](9236)',
    structure = SMILES('[CH2][CH][CH]C[C]([CH2])[CH2]'),
    E0 = (968.399,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,193.981,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00664571,'amu*angstrom^2'), symmetry=1, barrier=(13.0476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00664571,'amu*angstrom^2'), symmetry=1, barrier=(13.0476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00664571,'amu*angstrom^2'), symmetry=1, barrier=(13.0476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00664571,'amu*angstrom^2'), symmetry=1, barrier=(13.0476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00664571,'amu*angstrom^2'), symmetry=1, barrier=(13.0476,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00664571,'amu*angstrom^2'), symmetry=1, barrier=(13.0476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986664,0.0826242,-0.000138036,1.39883e-07,-5.22799e-11,116565,37.9646], Tmin=(100,'K'), Tmax=(897.467,'K')), NASAPolynomial(coeffs=[-2.26177,0.0559003,-2.4506e-05,4.39554e-09,-2.88496e-13,118807,62.5293], Tmin=(897.467,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(968.399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C([CH2])[CH2](9237)',
    structure = SMILES('[CH2][CH][CH][CH]C([CH2])[CH2]'),
    E0 = (977.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1832.14,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0115254,'amu*angstrom^2'), symmetry=1, barrier=(8.80513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0115254,'amu*angstrom^2'), symmetry=1, barrier=(8.80513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0115254,'amu*angstrom^2'), symmetry=1, barrier=(8.80513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0115254,'amu*angstrom^2'), symmetry=1, barrier=(8.80513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0115254,'amu*angstrom^2'), symmetry=1, barrier=(8.80513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0115254,'amu*angstrom^2'), symmetry=1, barrier=(8.80513,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09315,0.0739565,-0.000102416,9.58527e-08,-3.50468e-11,117663,39.0042], Tmin=(100,'K'), Tmax=(883.954,'K')), NASAPolynomial(coeffs=[0.893914,0.0506037,-2.16305e-05,3.88461e-09,-2.57446e-13,118646,45.3008], Tmin=(883.954,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(977.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C]([CH2])C(9046)',
    structure = SMILES('[CH2][CH][CH]C[C]([CH2])C'),
    E0 = (763.317,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1314.93,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136489,'amu*angstrom^2'), symmetry=1, barrier=(6.1637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998698,0.0806099,-0.000124564,1.26537e-07,-4.85632e-11,91899.9,36.8318], Tmin=(100,'K'), Tmax=(874.732,'K')), NASAPolynomial(coeffs=[-2.3313,0.0597485,-2.69057e-05,4.94172e-09,-3.31223e-13,93863.1,60.3435], Tmin=(874.732,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(763.317,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C([CH2])[CH2](5223)',
    structure = SMILES('[CH2][CH]C[CH]C([CH2])[CH2]'),
    E0 = (783.06,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,398.105,3369.16,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0219343,'amu*angstrom^2'), symmetry=1, barrier=(3.5154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0219343,'amu*angstrom^2'), symmetry=1, barrier=(3.5154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0219343,'amu*angstrom^2'), symmetry=1, barrier=(3.5154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0219343,'amu*angstrom^2'), symmetry=1, barrier=(3.5154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0219343,'amu*angstrom^2'), symmetry=1, barrier=(3.5154,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0219343,'amu*angstrom^2'), symmetry=1, barrier=(3.5154,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.860882,0.0732015,-7.73827e-05,5.81589e-08,-1.90296e-11,94289.4,37.2761], Tmin=(100,'K'), Tmax=(844.216,'K')), NASAPolynomial(coeffs=[5.13718,0.0466943,-1.91878e-05,3.43989e-09,-2.30478e-13,93790,18.6878], Tmin=(844.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(783.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C([CH2])C(9047)',
    structure = SMILES('[CH2][CH][CH][CH]C([CH2])C'),
    E0 = (772.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1692.17,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0133765,'amu*angstrom^2'), symmetry=1, barrier=(11.0799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12398,0.0717174,-8.81604e-05,8.1535e-08,-3.09609e-11,92998,37.8045], Tmin=(100,'K'), Tmax=(845.322,'K')), NASAPolynomial(coeffs=[0.71883,0.0546362,-2.41384e-05,4.4567e-09,-3.02345e-13,93745.3,43.7059], Tmin=(845.322,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(772.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[C]([CH2])[CH2](5224)',
    structure = SMILES('[CH2][CH]CC[C]([CH2])[CH2]'),
    E0 = (773.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,381.668,2909.72,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0176443,'amu*angstrom^2'), symmetry=1, barrier=(3.02532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0176443,'amu*angstrom^2'), symmetry=1, barrier=(3.02532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0176443,'amu*angstrom^2'), symmetry=1, barrier=(3.02532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0176443,'amu*angstrom^2'), symmetry=1, barrier=(3.02532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0176443,'amu*angstrom^2'), symmetry=1, barrier=(3.02532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0176443,'amu*angstrom^2'), symmetry=1, barrier=(3.02532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.719401,0.0823074,-0.000114653,1.04491e-07,-3.73112e-11,93192,36.3603], Tmin=(100,'K'), Tmax=(886.934,'K')), NASAPolynomial(coeffs=[2.1024,0.0517779,-2.19374e-05,3.92053e-09,-2.5898e-13,93902.2,35.2404], Tmin=(886.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(773.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH][CH]C([CH2])[CH2](9238)',
    structure = SMILES('[CH2]C[CH][CH]C([CH2])[CH2]'),
    E0 = (783.072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,933.234,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0132402,'amu*angstrom^2'), symmetry=1, barrier=(11.6134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0132402,'amu*angstrom^2'), symmetry=1, barrier=(11.6134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0132402,'amu*angstrom^2'), symmetry=1, barrier=(11.6134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0132402,'amu*angstrom^2'), symmetry=1, barrier=(11.6134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0132402,'amu*angstrom^2'), symmetry=1, barrier=(11.6134,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0132402,'amu*angstrom^2'), symmetry=1, barrier=(11.6134,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.897078,0.0736932,-8.26485e-05,6.66704e-08,-2.29422e-11,94288.3,36.9741], Tmin=(100,'K'), Tmax=(849.441,'K')), NASAPolynomial(coeffs=[4.05673,0.048624,-2.03848e-05,3.68098e-09,-2.47155e-13,94119.2,24.4102], Tmin=(849.441,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(783.072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C[C]([CH2])[CH2](5225)',
    structure = SMILES('[CH2]C[CH]C[C]([CH2])[CH2]'),
    E0 = (773.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,225.428,2981.42,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0104874,'amu*angstrom^2'), symmetry=1, barrier=(2.04636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104874,'amu*angstrom^2'), symmetry=1, barrier=(2.04636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104874,'amu*angstrom^2'), symmetry=1, barrier=(2.04636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104874,'amu*angstrom^2'), symmetry=1, barrier=(2.04636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104874,'amu*angstrom^2'), symmetry=1, barrier=(2.04636,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104874,'amu*angstrom^2'), symmetry=1, barrier=(2.04636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.764405,0.0826858,-0.000119474,1.12349e-07,-4.09062e-11,93190.5,36.0272], Tmin=(100,'K'), Tmax=(884.717,'K')), NASAPolynomial(coeffs=[1.00253,0.0537424,-2.31552e-05,4.16667e-09,-2.76084e-13,94239,41.0711], Tmin=(884.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(773.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH][CH][CH]C(4135)',
    structure = SMILES('[CH2]C([CH2])[CH][CH][CH]C'),
    E0 = (772.272,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,470.946,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0101099,'amu*angstrom^2'), symmetry=1, barrier=(8.80249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06545,0.0731073,-9.15492e-05,8.306e-08,-3.03861e-11,92980.3,37.3174], Tmin=(100,'K'), Tmax=(874.672,'K')), NASAPolynomial(coeffs=[1.06209,0.0529653,-2.24387e-05,4.03682e-09,-2.68864e-13,93752,41.741], Tmin=(874.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(772.272,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C[CH][CH]C(4134)',
    structure = SMILES('[CH2][C]([CH2])C[CH][CH]C'),
    E0 = (763.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1881.4,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199814,'amu*angstrom^2'), symmetry=1, barrier=(5.1281,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953113,0.0818466,-0.00012743,1.27439e-07,-4.77692e-11,91881.7,36.2985], Tmin=(100,'K'), Tmax=(892.802,'K')), NASAPolynomial(coeffs=[-2.06846,0.0582179,-2.52885e-05,4.54158e-09,-2.99396e-13,93902.5,58.8287], Tmin=(892.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(763.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (782.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (940.295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (948.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1239.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1246.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1248.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1248.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1237.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1237.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (791.066,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (791.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (805.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (846.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (846.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (861.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (813.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (857.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (832.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (796.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1140.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1180.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1189.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (924.464,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (920.064,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (900.916,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (899.112,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (957.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (935.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (858.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (856.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['allyl(82)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2][CH][CH]CC[CH][CH2](4163)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C[CH][CH2](653)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH][CH][CH2](5537)', '[CH2]C([CH2])[CH2](489)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(4.02491e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2][CH][C]CC([CH2])[CH2](9231)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2][C][CH]CC([CH2])[CH2](9232)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C([CH2])C[CH][CH][CH2](9233)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH][CH][CH]CC([CH2])[CH2](9234)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2][CH][CH]CC1CC1(9235)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2][CH]C1CC([CH2])C1(5228)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2]C=CC[C]([CH2])C(5011)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2][C]([CH2])CCC=C(5242)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2]C=C[CH]C([CH2])C(5012)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2]C=CC[C]([CH2])[CH2](5231)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]C=C[CH]C([CH2])[CH2](5232)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2(T)(20)', '[CH2][CH]CC=C[CH2](652)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;Y_1centerbirad] for rate rule [Cds-CsH_Cds-HH;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['allyl(82)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.0036002,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2][CH][CH]C[C]([CH2])[CH2](9236)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][CH][CH][CH]C([CH2])[CH2](9237)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2][CH][CH]C[C]([CH2])C(9046)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]C[CH]C([CH2])[CH2](5223)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2][CH][CH][CH]C([CH2])C(9047)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(333380,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2][CH]CC[C]([CH2])[CH2](5224)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C[CH][CH]C([CH2])[CH2](9238)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(867345,'s^-1'), n=1.96939, Ea=(174.054,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R3HJ;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2]C[CH]C[C]([CH2])[CH2](5225)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(322699,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_1;Y_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2]C([CH2])[CH][CH][CH]C(4135)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.32767e+06,'s^-1'), n=1.53625, Ea=(75.6258,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH][CH]CC([CH2])[CH2](4136)'],
    products = ['[CH2][C]([CH2])C[CH][CH]C(4134)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2059',
    isomers = [
        '[CH2][CH][CH]CC([CH2])[CH2](4136)',
    ],
    reactants = [
        ('allyl(82)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2059',
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


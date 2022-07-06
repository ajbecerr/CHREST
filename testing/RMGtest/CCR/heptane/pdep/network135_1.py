species(
    label = '[CH2][CH]CC[CH][CH2](450)',
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
    label = '[CH2]C=C(53)',
    structure = SMILES('[CH2]C=C'),
    E0 = (157.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.570286,'amu*angstrom^2'), symmetry=1, barrier=(32.8574,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2161.77,'J/mol'), sigma=(4.85,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31926,0.00566541,4.27429e-05,-5.78804e-08,2.21687e-11,18990.6,9.19662], Tmin=(100,'K'), Tmax=(952.01,'K')), NASAPolynomial(coeffs=[7.5573,0.0114809,-3.63938e-06,6.63549e-10,-4.95289e-14,17113.2,-16.6632], Tmin=(952.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CCC=C(447)',
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
    label = '[CH2][CH]CC([CH2])[CH2](426)',
    structure = SMILES('[CH2][CH]CC([CH2])[CH2]'),
    E0 = (612.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1170.47,1171.2],'cm^-1')),
        HinderedRotor(inertia=(0.169767,'amu*angstrom^2'), symmetry=1, barrier=(3.90329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00400723,'amu*angstrom^2'), symmetry=1, barrier=(3.90322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169659,'amu*angstrom^2'), symmetry=1, barrier=(3.90079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169917,'amu*angstrom^2'), symmetry=1, barrier=(3.90674,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00401056,'amu*angstrom^2'), symmetry=1, barrier=(3.90297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64314,0.0551022,-3.22884e-05,-1.61575e-09,1.00735e-11,73724.5,29.7062], Tmin=(100,'K'), Tmax=(647.823,'K')), NASAPolynomial(coeffs=[5.92384,0.0380634,-1.45841e-05,2.54566e-09,-1.69302e-13,72972.8,9.39149], Tmin=(647.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH][CH2](215)',
    structure = SMILES('[CH][CH2]'),
    E0 = (556.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1101.59,1101.66],'cm^-1')),
        HinderedRotor(inertia=(0.00420677,'amu*angstrom^2'), symmetry=1, barrier=(3.62356,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77493,-0.000462565,3.18167e-05,-4.30783e-08,1.77606e-11,66973.8,8.79001], Tmin=(100,'K'), Tmax=(870.354,'K')), NASAPolynomial(coeffs=[6.06996,0.00332438,5.85464e-07,-2.33e-10,1.82456e-14,66031.4,-5.08253], Tmin=(870.354,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C[CH2](85)',
    structure = SMILES('[CH2][CH]C[CH2]'),
    E0 = (460.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2622.36],'cm^-1')),
        HinderedRotor(inertia=(0.0897296,'amu*angstrom^2'), symmetry=1, barrier=(2.06306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00870163,'amu*angstrom^2'), symmetry=1, barrier=(42.4425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0245193,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81948,0.0272005,-9.56974e-06,5.91729e-10,1.83286e-13,55442.9,20.1592], Tmin=(100,'K'), Tmax=(1946.59,'K')), NASAPolynomial(coeffs=[9.19037,0.0193116,-7.49963e-06,1.25572e-09,-7.83181e-14,51976.9,-17.352], Tmin=(1946.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = 'H(9)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2][C]CC[CH][CH2](591)',
    structure = SMILES('[CH2][C]CC[CH][CH2]'),
    E0 = (861.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,207.395,800.047,1212.02,1599.91],'cm^-1')),
        HinderedRotor(inertia=(0.106001,'amu*angstrom^2'), symmetry=1, barrier=(3.50929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106001,'amu*angstrom^2'), symmetry=1, barrier=(3.50929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106001,'amu*angstrom^2'), symmetry=1, barrier=(3.50929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106001,'amu*angstrom^2'), symmetry=1, barrier=(3.50929,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106001,'amu*angstrom^2'), symmetry=1, barrier=(3.50929,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7295,0.0558791,-3.72741e-05,-9.53696e-10,1.37784e-11,103665,28.5203], Tmin=(100,'K'), Tmax=(564.701,'K')), NASAPolynomial(coeffs=[5.22586,0.0398281,-1.77877e-05,3.37097e-09,-2.35356e-13,103132,12.42], Tmin=(564.701,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(861.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]CC[CH][CH2](592)',
    structure = SMILES('[CH][CH]CC[CH][CH2]'),
    E0 = (850.49,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,188.915,593.676,984.496,1922.58,2850.26,4000],'cm^-1')),
        HinderedRotor(inertia=(0.119915,'amu*angstrom^2'), symmetry=1, barrier=(2.77593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119915,'amu*angstrom^2'), symmetry=1, barrier=(2.77593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119915,'amu*angstrom^2'), symmetry=1, barrier=(2.77593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119915,'amu*angstrom^2'), symmetry=1, barrier=(2.77593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119915,'amu*angstrom^2'), symmetry=1, barrier=(2.77593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49138,0.0611909,-7.35099e-05,6.47858e-08,-2.40302e-11,102375,31.3588], Tmin=(100,'K'), Tmax=(826.786,'K')), NASAPolynomial(coeffs=[2.73861,0.0432225,-1.92586e-05,3.58241e-09,-2.44659e-13,102576,28.0456], Tmin=(826.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(850.49,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1CCC1[CH2](576)',
    structure = SMILES('[CH2]C1CCC1[CH2]'),
    E0 = (358.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96626,0.0264433,6.69061e-05,-1.05634e-07,4.30348e-11,43194.9,22.5781], Tmin=(100,'K'), Tmax=(931.069,'K')), NASAPolynomial(coeffs=[12.9666,0.0258012,-7.16177e-06,1.1749e-09,-8.34162e-14,39125.9,-40.5575], Tmin=(931.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C=CCC[CH2](593)',
    structure = SMILES('[CH2]C=CCC[CH2]'),
    E0 = (280.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25947,0.0504297,-1.06899e-05,-1.69915e-08,8.70219e-12,33833.6,26.035], Tmin=(100,'K'), Tmax=(1056.37,'K')), NASAPolynomial(coeffs=[11.6223,0.0314484,-1.25029e-05,2.30642e-09,-1.61128e-13,30513.9,-29.8779], Tmin=(1056.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC=C[CH2](594)',
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
    label = '[CH2][CH][CH2](210)',
    structure = SMILES('[CH2][CH][CH2]'),
    E0 = (484.293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.00132735,'amu*angstrom^2'), symmetry=1, barrier=(2.41068,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00132737,'amu*angstrom^2'), symmetry=1, barrier=(2.41079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.34321,0.013802,2.16431e-06,-5.76334e-09,1.61333e-12,58271,14.955], Tmin=(100,'K'), Tmax=(1447.11,'K')), NASAPolynomial(coeffs=[4.39501,0.0167645,-6.99094e-06,1.25742e-09,-8.38113e-14,57352,7.36834], Tmin=(1447.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH][CH2](595)',
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
    label = '[CH2][CH]C[CH]C[CH2](596)',
    structure = SMILES('[CH2][CH]C[CH]C[CH2]'),
    E0 = (607.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18168,0.0425548,-1.66264e-05,2.15789e-09,-2.55685e-14,73073.4,25.7665], Tmin=(100,'K'), Tmax=(2657.42,'K')), NASAPolynomial(coeffs=[33.3343,0.00857218,-3.88157e-06,5.75411e-10,-2.86112e-14,53021.2,-156.748], Tmin=(2657.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC[CH2](597)',
    structure = SMILES('[CH2][CH][CH]CC[CH2]'),
    E0 = (607.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.18168,0.0425548,-1.66264e-05,2.15789e-09,-2.55685e-14,73073.4,25.7665], Tmin=(100,'K'), Tmax=(2657.42,'K')), NASAPolynomial(coeffs=[33.3343,0.00857218,-3.88157e-06,5.75411e-10,-2.86112e-14,53021.2,-156.748], Tmin=(2657.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.533,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH][CH]C(448)',
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
    label = '[CH2][CH][CH]C[CH]C(449)',
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
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
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
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (607.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (607.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (769.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1074.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1073.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1062.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (615.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (670.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (694.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (677.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (968.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1013.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (764.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (723.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (751.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (738.646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC[CH][CH2](450)'],
    products = ['[CH2]C=C(53)', '[CH2]C=C(53)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]CC[CH][CH2](450)'],
    products = ['[CH2][CH]CCC=C(447)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]CC([CH2])[CH2](426)'],
    products = ['[CH2][CH]CC[CH][CH2](450)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH][CH2](215)', '[CH2][CH]C[CH2](85)'],
    products = ['[CH2][CH]CC[CH][CH2](450)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.22604e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(9)', '[CH2][C]CC[CH][CH2](591)'],
    products = ['[CH2][CH]CC[CH][CH2](450)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(9)', '[CH][CH]CC[CH][CH2](592)'],
    products = ['[CH2][CH]CC[CH][CH2](450)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][CH]CC[CH][CH2](450)'],
    products = ['[CH2]C1CCC1[CH2](576)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH]CC[CH][CH2](450)'],
    products = ['[CH2]C=CCC[CH2](593)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(9)', '[CH2][CH]CC=C[CH2](594)'],
    products = ['[CH2][CH]CC[CH][CH2](450)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C=C(53)', '[CH2][CH][CH2](210)'],
    products = ['[CH2][CH]CC[CH][CH2](450)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(0.493875,'m^3/(mol*s)'), n=2.00579, Ea=(36.0234,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH][CH2](210)', '[CH2][CH][CH2](210)'],
    products = ['[CH2][CH]CC[CH][CH2](450)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(9)', '[CH2][CH][CH]C[CH][CH2](595)'],
    products = ['[CH2][CH]CC[CH][CH2](450)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]CC[CH][CH2](450)'],
    products = ['[CH2][CH]C[CH]C[CH2](596)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.58236e+06,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]CC[CH][CH2](450)'],
    products = ['[CH2][CH][CH]CC[CH2](597)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.1424e+06,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]CC[CH][CH2](450)'],
    products = ['[CH2][CH]C[CH][CH]C(448)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.01841e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]CC[CH][CH2](450)'],
    products = ['[CH2][CH][CH]C[CH]C(449)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1739.66,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #135',
    isomers = [
        '[CH2][CH]CC[CH][CH2](450)',
    ],
    reactants = [
        ('[CH2]C=C(53)', '[CH2]C=C(53)'),
    ],
    bathGas = {
        'Ar': 0.25,
        'He': 0.25,
        'Ne': 0.25,
        'N2': 0.25,
    },
)

pressureDependence(
    label = 'PDepNetwork #135',
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


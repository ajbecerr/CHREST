species(
    label = '[CH2][CH]CC([CH2])[CH2](484)',
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64314,0.0551022,-3.22884e-05,-1.61575e-09,1.00735e-11,73724.5,29.7062], Tmin=(100,'K'), Tmax=(647.823,'K')), NASAPolynomial(coeffs=[5.92384,0.0380634,-1.45841e-05,2.54566e-09,-1.69302e-13,72972.8,9.39149], Tmin=(647.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2]C([CH2])CC=C(486)',
    structure = SMILES('[CH2]C([CH2])CC=C'),
    E0 = (337.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,368.087,1588.75],'cm^-1')),
        HinderedRotor(inertia=(0.00124259,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0679355,'amu*angstrom^2'), symmetry=1, barrier=(6.53436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0678612,'amu*angstrom^2'), symmetry=1, barrier=(6.53507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724317,'amu*angstrom^2'), symmetry=1, barrier=(69.7204,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29622,0.0520098,-1.64684e-05,-1.38417e-08,9.25185e-12,40640.7,26.8619], Tmin=(100,'K'), Tmax=(947.419,'K')), NASAPolynomial(coeffs=[10.6916,0.0302257,-1.02919e-05,1.73514e-09,-1.15664e-13,38057.9,-22.1933], Tmin=(947.419,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3248.64,'J/mol'), sigma=(6.11673,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=507.43 K, Pc=32.21 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.6689,0.0650717,-5.68349e-05,2.94125e-08,-6.08673e-12,74344,30.5296], Tmin=(100,'K'), Tmax=(1324.25,'K')), NASAPolynomial(coeffs=[11.0931,0.0274373,-7.24249e-06,9.40708e-10,-4.98543e-14,72122.2,-20.6535], Tmin=(1324.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2][CH]C[CH][CH2](510)',
    structure = SMILES('[CH2][CH]C[CH][CH2]'),
    E0 = (631.301,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1407.99,1417.06],'cm^-1')),
        HinderedRotor(inertia=(0.00560336,'amu*angstrom^2'), symmetry=1, barrier=(7.95035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00563889,'amu*angstrom^2'), symmetry=1, barrier=(7.95031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00557143,'amu*angstrom^2'), symmetry=1, barrier=(7.9482,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00562946,'amu*angstrom^2'), symmetry=1, barrier=(7.93664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26906,0.0435841,-4.87505e-05,4.66768e-08,-1.86942e-11,75984.9,26.4118], Tmin=(100,'K'), Tmax=(836.615,'K')), NASAPolynomial(coeffs=[0.578549,0.0388266,-1.7199e-05,3.18971e-09,-2.17339e-13,76717.2,36.9515], Tmin=(836.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][C]CC([CH2])[CH2](637)',
    structure = SMILES('[CH2][C]CC([CH2])[CH2]'),
    E0 = (866.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.928003,0.064135,-5.80072e-05,3.03055e-08,-6.51295e-12,104277,29.4001], Tmin=(100,'K'), Tmax=(1118.33,'K')), NASAPolynomial(coeffs=[10.949,0.0282919,-9.9306e-06,1.64537e-09,-1.05952e-13,102036,-20.0658], Tmin=(1118.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(866.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])C[CH][CH2](638)',
    structure = SMILES('[CH]C([CH2])C[CH][CH2]'),
    E0 = (855.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,210.214,429.066,1087.85,1965.66,3038.1],'cm^-1')),
        HinderedRotor(inertia=(0.0721512,'amu*angstrom^2'), symmetry=1, barrier=(1.66657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721512,'amu*angstrom^2'), symmetry=1, barrier=(1.66657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721512,'amu*angstrom^2'), symmetry=1, barrier=(1.66657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721512,'amu*angstrom^2'), symmetry=1, barrier=(1.66657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721512,'amu*angstrom^2'), symmetry=1, barrier=(1.66657,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38542,0.0592899,-5.27549e-05,3.01052e-08,-7.57619e-12,102977,30.5121], Tmin=(100,'K'), Tmax=(929.473,'K')), NASAPolynomial(coeffs=[6.908,0.0355222,-1.43962e-05,2.59095e-09,-1.75334e-13,101951,4.27303], Tmin=(929.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(855.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]CC([CH2])[CH2](639)',
    structure = SMILES('[CH][CH]CC([CH2])[CH2]'),
    E0 = (855.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,199.748,1040.96,1367.2,1898.02,2038.31],'cm^-1')),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118398,'amu*angstrom^2'), symmetry=1, barrier=(2.82079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52842,0.0577316,-4.2467e-05,7.26258e-09,7.53295e-12,102951,29.3345], Tmin=(100,'K'), Tmax=(654.247,'K')), NASAPolynomial(coeffs=[7.16544,0.0340336,-1.28177e-05,2.20259e-09,-1.44655e-13,101983,2.76985], Tmin=(654.247,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(855.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]CC1CC1(640)',
    structure = SMILES('[CH2][CH]CC1CC1'),
    E0 = (361.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83847,0.0358897,2.43838e-05,-4.9718e-08,1.96225e-11,43552.4,25.6477], Tmin=(100,'K'), Tmax=(1004.95,'K')), NASAPolynomial(coeffs=[9.57061,0.0328878,-1.25921e-05,2.31296e-09,-1.62792e-13,40595.8,-18.6715], Tmin=(1004.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(361.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1CC([CH2])C1(641)',
    structure = SMILES('[CH2]C1CC([CH2])C1'),
    E0 = (358.394,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96626,0.0264433,6.69061e-05,-1.05634e-07,4.30348e-11,43194.9,21.8849], Tmin=(100,'K'), Tmax=(931.069,'K')), NASAPolynomial(coeffs=[12.9666,0.0258012,-7.16177e-06,1.1749e-09,-8.34162e-14,39125.9,-41.2507], Tmin=(931.069,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CC(=C)C(531)',
    structure = SMILES('[CH2][CH]CC(=C)C'),
    E0 = (320.336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180,1599.33],'cm^-1')),
        HinderedRotor(inertia=(0.0892705,'amu*angstrom^2'), symmetry=1, barrier=(2.0525,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0889334,'amu*angstrom^2'), symmetry=1, barrier=(2.04475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0891168,'amu*angstrom^2'), symmetry=1, barrier=(2.04897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0889638,'amu*angstrom^2'), symmetry=1, barrier=(2.04545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82264,0.0512649,-2.77873e-05,7.49967e-09,-8.56116e-13,38601.8,25.3925], Tmin=(100,'K'), Tmax=(1779.89,'K')), NASAPolynomial(coeffs=[8.6284,0.0359703,-1.48979e-05,2.67195e-09,-1.78032e-13,36179.1,-11.3654], Tmin=(1779.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(320.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CCC([CH2])=C(642)',
    structure = SMILES('[CH2]CCC([CH2])=C'),
    E0 = (277.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0844,0.0545993,-1.97902e-05,-9.49043e-09,6.52864e-12,33475.2,24.9206], Tmin=(100,'K'), Tmax=(1058.22,'K')), NASAPolynomial(coeffs=[12.3682,0.0305284,-1.20084e-05,2.19976e-09,-1.53007e-13,30046.7,-35.0717], Tmin=(1058.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC([CH2])C(532)',
    structure = SMILES('[CH2]C=CC([CH2])C'),
    E0 = (272.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.00511567,'amu*angstrom^2'), symmetry=1, barrier=(4.60508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.679991,'amu*angstrom^2'), symmetry=1, barrier=(15.6343,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0173188,'amu*angstrom^2'), symmetry=1, barrier=(15.6385,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.30708,'amu*angstrom^2'), symmetry=1, barrier=(76.0364,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32269,0.0477996,1.34752e-06,-3.4245e-08,1.64564e-11,32854.3,25.5549], Tmin=(100,'K'), Tmax=(967.436,'K')), NASAPolynomial(coeffs=[12.2804,0.0288314,-1.00791e-05,1.77001e-09,-1.22476e-13,29501.6,-33.3169], Tmin=(967.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC([CH2])=C(643)',
    structure = SMILES('[CH2][CH]CC([CH2])=C'),
    E0 = (471.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,246.856,1775.46],'cm^-1')),
        HinderedRotor(inertia=(0.0051858,'amu*angstrom^2'), symmetry=1, barrier=(11.5977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.34837e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.701505,'amu*angstrom^2'), symmetry=1, barrier=(30.3593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01076,'amu*angstrom^2'), symmetry=1, barrier=(87.6001,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42509,0.0531145,-3.31274e-05,1.07407e-08,-1.45659e-12,56843.8,26.4331], Tmin=(100,'K'), Tmax=(1641.34,'K')), NASAPolynomial(coeffs=[10.9943,0.029794,-1.1815e-05,2.08414e-09,-1.38063e-13,53702.6,-24.4744], Tmin=(1641.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJC) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])[CH2](644)',
    structure = SMILES('[CH2]C=CC([CH2])[CH2]'),
    E0 = (477.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,397.773],'cm^-1')),
        HinderedRotor(inertia=(0.00475792,'amu*angstrom^2'), symmetry=1, barrier=(0.534766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607471,'amu*angstrom^2'), symmetry=1, barrier=(72.6648,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0695163,'amu*angstrom^2'), symmetry=1, barrier=(7.65916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0302878,'amu*angstrom^2'), symmetry=1, barrier=(72.6694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3691,0.049131,-9.81136e-06,-2.36824e-08,1.37799e-11,57516.4,26.478], Tmin=(100,'K'), Tmax=(924.208,'K')), NASAPolynomial(coeffs=[11.9922,0.0255947,-8.03361e-06,1.30776e-09,-8.6727e-14,54594.4,-29.1203], Tmin=(924.208,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC=C(508)',
    structure = SMILES('[CH2][CH]CC=C'),
    E0 = (359.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,1243.57],'cm^-1')),
        HinderedRotor(inertia=(0.00331557,'amu*angstrom^2'), symmetry=1, barrier=(3.65555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160798,'amu*angstrom^2'), symmetry=1, barrier=(3.69707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161411,'amu*angstrom^2'), symmetry=1, barrier=(3.71116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13435,0.0369817,-1.51234e-05,1.33725e-09,3.77008e-13,43295.2,23.5599], Tmin=(100,'K'), Tmax=(1532.36,'K')), NASAPolynomial(coeffs=[8.68459,0.025936,-1.02358e-05,1.78847e-09,-1.17139e-13,40577.1,-13.1557], Tmin=(1532.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C[C]([CH2])[CH2](645)',
    structure = SMILES('[CH2][CH]C[C]([CH2])[CH2]'),
    E0 = (797.721,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,2201.03],'cm^-1')),
        HinderedRotor(inertia=(0.0014274,'amu*angstrom^2'), symmetry=1, barrier=(4.89245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188804,'amu*angstrom^2'), symmetry=1, barrier=(64.772,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.81006,'amu*angstrom^2'), symmetry=1, barrier=(64.6089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00142474,'amu*angstrom^2'), symmetry=1, barrier=(4.89662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00142217,'amu*angstrom^2'), symmetry=1, barrier=(4.89392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41035,0.066982,-9.87864e-05,9.35429e-08,-3.37862e-11,96027.4,31.6356], Tmin=(100,'K'), Tmax=(901.876,'K')), NASAPolynomial(coeffs=[1.02445,0.0436355,-1.82801e-05,3.22553e-09,-2.10453e-13,97116.1,39.1073], Tmin=(901.876,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(797.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C([CH2])[CH2](646)',
    structure = SMILES('[CH2][CH][CH]C([CH2])[CH2]'),
    E0 = (806.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,1304.76],'cm^-1')),
        HinderedRotor(inertia=(0.00319015,'amu*angstrom^2'), symmetry=1, barrier=(3.8543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318875,'amu*angstrom^2'), symmetry=1, barrier=(3.85345,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0031889,'amu*angstrom^2'), symmetry=1, barrier=(3.85323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167613,'amu*angstrom^2'), symmetry=1, barrier=(3.85376,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00319051,'amu*angstrom^2'), symmetry=1, barrier=(3.85453,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52596,0.0581986,-6.27215e-05,4.88717e-08,-1.62492e-11,97126,32.6431], Tmin=(100,'K'), Tmax=(876.179,'K')), NASAPolynomial(coeffs=[4.1578,0.0383786,-1.54281e-05,2.72027e-09,-1.79882e-13,96964.3,22.0036], Tmin=(876.179,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(806.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[C]([CH2])C(533)',
    structure = SMILES('[CH2][CH]C[C]([CH2])C'),
    E0 = (592.639,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1668.38,1668.4],'cm^-1')),
        HinderedRotor(inertia=(0.00254266,'amu*angstrom^2'), symmetry=1, barrier=(5.02143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218375,'amu*angstrom^2'), symmetry=1, barrier=(5.02086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00254274,'amu*angstrom^2'), symmetry=1, barrier=(5.02251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218363,'amu*angstrom^2'), symmetry=1, barrier=(5.0206,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218382,'amu*angstrom^2'), symmetry=1, barrier=(5.02102,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4329,0.0648358,-8.4817e-05,7.94975e-08,-2.97485e-11,71362.4,30.4657], Tmin=(100,'K'), Tmax=(864.649,'K')), NASAPolynomial(coeffs=[0.922292,0.0475413,-2.07138e-05,3.77988e-09,-2.53868e-13,72185.4,37.1039], Tmin=(864.649,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.639,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C([CH2])[CH2](647)',
    structure = SMILES('[CH2]C[CH]C([CH2])[CH2]'),
    E0 = (612.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1182.51,1182.53],'cm^-1')),
        HinderedRotor(inertia=(0.00227336,'amu*angstrom^2'), symmetry=1, barrier=(2.25593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0981085,'amu*angstrom^2'), symmetry=1, barrier=(2.25571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0981609,'amu*angstrom^2'), symmetry=1, barrier=(2.25691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0981153,'amu*angstrom^2'), symmetry=1, barrier=(2.25586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0981001,'amu*angstrom^2'), symmetry=1, barrier=(2.25551,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34718,0.0577643,-4.24404e-05,1.89839e-08,-3.752e-12,73750,30.5481], Tmin=(100,'K'), Tmax=(1153.11,'K')), NASAPolynomial(coeffs=[7.69479,0.0357454,-1.37978e-05,2.42436e-09,-1.61845e-13,72286.1,-0.979833], Tmin=(1153.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C([CH2])C(534)',
    structure = SMILES('[CH2][CH][CH]C([CH2])C'),
    E0 = (601.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1361.73,1362.64],'cm^-1')),
        HinderedRotor(inertia=(0.0028834,'amu*angstrom^2'), symmetry=1, barrier=(3.79819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165946,'amu*angstrom^2'), symmetry=1, barrier=(3.81544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166303,'amu*angstrom^2'), symmetry=1, barrier=(3.82363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165528,'amu*angstrom^2'), symmetry=1, barrier=(3.80582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00289278,'amu*angstrom^2'), symmetry=1, barrier=(3.80392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43307,0.0419375,2.62092e-05,-1.26297e-07,1.09211e-10,72423.3,28.4935], Tmin=(100,'K'), Tmax=(437.212,'K')), NASAPolynomial(coeffs=[3.59076,0.0431407,-1.83847e-05,3.4034e-09,-2.34317e-13,72209.4,22.5769], Tmin=(437.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC[C]([CH2])[CH2](648)',
    structure = SMILES('[CH2]CC[C]([CH2])[CH2]'),
    E0 = (603.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19541,0.0669506,-7.98669e-05,6.54923e-08,-2.2167e-11,72653.1,29.6725], Tmin=(100,'K'), Tmax=(879.613,'K')), NASAPolynomial(coeffs=[4.2714,0.0415085,-1.69476e-05,3.00107e-09,-1.98414e-13,72555.1,17.746], Tmin=(879.613,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])[CH][CH]C(547)',
    structure = SMILES('[CH2]C([CH2])[CH][CH]C'),
    E0 = (601.594,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1338.21,1338.68],'cm^-1')),
        HinderedRotor(inertia=(0.166056,'amu*angstrom^2'), symmetry=1, barrier=(3.81796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00300753,'amu*angstrom^2'), symmetry=1, barrier=(3.82466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166254,'amu*angstrom^2'), symmetry=1, barrier=(3.82251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.166265,'amu*angstrom^2'), symmetry=1, barrier=(3.82275,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0030015,'amu*angstrom^2'), symmetry=1, barrier=(3.8146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51846,0.0570927,-5.08681e-05,3.46656e-08,-1.09231e-11,72442,30.885], Tmin=(100,'K'), Tmax=(846.21,'K')), NASAPolynomial(coeffs=[4.26796,0.040843,-1.62975e-05,2.88724e-09,-1.92544e-13,72093.1,18.7676], Tmin=(846.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(601.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([CH2])C[CH]C(546)',
    structure = SMILES('[CH2][C]([CH2])C[CH]C'),
    E0 = (592.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2218.22,2219.39],'cm^-1')),
        HinderedRotor(inertia=(0.0022833,'amu*angstrom^2'), symmetry=1, barrier=(7.97352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186934,'amu*angstrom^2'), symmetry=1, barrier=(65.3551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0187096,'amu*angstrom^2'), symmetry=1, barrier=(65.3131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346726,'amu*angstrom^2'), symmetry=1, barrier=(7.97191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.346658,'amu*angstrom^2'), symmetry=1, barrier=(7.97034,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37826,0.0661856,-8.81064e-05,8.0988e-08,-2.92209e-11,71344.5,29.9645], Tmin=(100,'K'), Tmax=(894.757,'K')), NASAPolynomial(coeffs=[1.21591,0.0459566,-1.90646e-05,3.37207e-09,-2.21397e-13,72212.4,35.417], Tmin=(894.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(592.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (612.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (612.298,'kJ/mol'),
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
    E0 = (774.395,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1069.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1076.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1077.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1067.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1067.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (620.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (620.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (635.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (675.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (675.699,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (694.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (696.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (741.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (662.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (969.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1009.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1018.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (753.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (749.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (730.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (728.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (756.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (743.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['allyl(82)', 'allyl(82)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2]C([CH2])CC=C(486)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2][CH]CC[CH][CH2](483)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([CH2])C([CH2])[CH2](485)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.32e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', '[CH2][CH]C[CH][CH2](510)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][CH2](502)', '[CH2]C([CH2])[CH2](489)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(4.02491e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2][C]CC([CH2])[CH2](637)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C([CH2])C[CH][CH2](638)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH][CH]CC([CH2])[CH2](639)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2][CH]CC1CC1(640)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2]C1CC([CH2])C1(641)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2][CH]CC(=C)C(531)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2]CCC([CH2])=C(642)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2]C=CC([CH2])C(532)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2][CH]CC([CH2])=C(643)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]C=CC([CH2])[CH2](644)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CH2(T)(20)', '[CH2][CH]CC=C(508)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;Y_1centerbirad] for rate rule [Cds-CsH_Cds-HH;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH][CH2](497)', 'allyl(82)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.0036002,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH][CH2](497)', '[CH2][CH][CH2](497)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2][CH]C[C]([CH2])[CH2](645)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2][CH][CH]C([CH2])[CH2](646)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2][CH]C[C]([CH2])C(533)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C[CH]C([CH2])[CH2](647)'],
    products = ['[CH2][CH]CC([CH2])[CH2](484)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2][CH][CH]C([CH2])C(534)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(333380,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2]CC[C]([CH2])[CH2](648)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2]C([CH2])[CH][CH]C(547)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.50921e+06,'s^-1'), n=1.8375, Ea=(144.331,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]CC([CH2])[CH2](484)'],
    products = ['[CH2][C]([CH2])C[CH]C(546)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #241',
    isomers = [
        '[CH2][CH]CC([CH2])[CH2](484)',
    ],
    reactants = [
        ('allyl(82)', 'allyl(82)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #241',
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


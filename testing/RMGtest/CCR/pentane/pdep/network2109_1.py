species(
    label = '[CH2][CH][CH]CC[CH][O](4354)',
    structure = SMILES('[CH2][CH][CH]CC[CH][O]'),
    E0 = (648.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,476.327,863.886,2166.31,2892.04,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.944353,0.0844202,-0.000148429,1.54869e-07,-5.98592e-11,78085.5,35.2576], Tmin=(100,'K'), Tmax=(866.461,'K')), NASAPolynomial(coeffs=[-2.27782,0.0568642,-2.72688e-05,5.12854e-09,-3.47543e-13,80236.6,59.5321], Tmin=(866.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = 'vinoxy(1351)',
    structure = SMILES('[CH2]C=O'),
    E0 = (6.53237,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,1127.76,1127.76,1127.77,1945],'cm^-1')),
        HinderedRotor(inertia=(0.00342653,'amu*angstrom^2'), symmetry=1, barrier=(3.09254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.51396,0.0045106,2.69938e-05,-3.71579e-08,1.42465e-11,808.756,8.86322], Tmin=(100,'K'), Tmax=(957.493,'K')), NASAPolynomial(coeffs=[6.49564,0.00770518,-2.52913e-06,4.69083e-10,-3.5128e-14,-479.661,-9.13857], Tmin=(957.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(6.53237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""vinoxy""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][CH][CH]CCC=O(4258)',
    structure = SMILES('[CH2][CH][CH]CCC=O'),
    E0 = (319.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,180.195,1248.65,3975.57],'cm^-1')),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000782982,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3572.62,'J/mol'), sigma=(6.29932,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.04 K, Pc=32.43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.96258,0.0421737,-1.60296e-05,1.66613e-09,4.99658e-14,38383.4,23.0896], Tmin=(100,'K'), Tmax=(2776.37,'K')), NASAPolynomial(coeffs=[48.7984,-0.00559881,6.91183e-07,-1.66324e-10,1.84404e-14,7003.1,-250.678], Tmin=(2776.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH][O](4373)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH][O]'),
    E0 = (656.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,236.693,456.946,2042.85,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0325019,'amu*angstrom^2'), symmetry=1, barrier=(3.58182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.702035,0.0817735,-0.00011892,1.08066e-07,-3.90055e-11,79108.8,34.8881], Tmin=(100,'K'), Tmax=(846.966,'K')), NASAPolynomial(coeffs=[4.17012,0.0459529,-2.10488e-05,3.92751e-09,-2.67084e-13,79218.7,22.8492], Tmin=(846.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])[O](4458)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])[O]'),
    E0 = (671.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,686.59,2061.2,2278.4,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852795,0.0775059,-0.000107974,9.70252e-08,-3.50878e-11,80900.1,35.2748], Tmin=(100,'K'), Tmax=(840.862,'K')), NASAPolynomial(coeffs=[4.01348,0.0455248,-2.06945e-05,3.86014e-09,-2.62948e-13,80967.6,24.1364], Tmin=(840.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
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
    label = '[CH2]C[CH][O](1549)',
    structure = SMILES('[CH2]C[CH][O]'),
    E0 = (338.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1524.15,1525.95],'cm^-1')),
        HinderedRotor(inertia=(0.310045,'amu*angstrom^2'), symmetry=1, barrier=(7.12854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00431249,'amu*angstrom^2'), symmetry=1, barrier=(7.11449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0712,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.81289,0.029873,-3.56716e-05,3.591e-08,-1.53006e-11,40735.5,17.4078], Tmin=(100,'K'), Tmax=(773.473,'K')), NASAPolynomial(coeffs=[2.08215,0.0257045,-1.21752e-05,2.37389e-09,-1.67404e-13,41086.3,22.2823], Tmin=(773.473,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(RCCJ) + radical(CCsJOH)"""),
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
    label = '[CH2][CH][CH]C[CH2](509)',
    structure = SMILES('[CH2][CH][CH]C[CH2]'),
    E0 = (631.313,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1539.83],'cm^-1')),
        HinderedRotor(inertia=(0.00318581,'amu*angstrom^2'), symmetry=1, barrier=(5.35941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318763,'amu*angstrom^2'), symmetry=1, barrier=(5.35914,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318741,'amu*angstrom^2'), symmetry=1, barrier=(5.35718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00318288,'amu*angstrom^2'), symmetry=1, barrier=(5.35416,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30707,0.0440519,-5.39188e-05,5.50378e-08,-2.25294e-11,75983.8,26.7965], Tmin=(100,'K'), Tmax=(841.968,'K')), NASAPolynomial(coeffs=[-0.503948,0.0407601,-1.83983e-05,3.43135e-09,-2.34064e-13,77047.1,43.3784], Tmin=(841.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(631.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][CH][C]CC[CH][O](9285)',
    structure = SMILES('[CH2][CH][C]CC[CH][O]'),
    E0 = (902.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619472,0.0872153,-0.000144727,1.39297e-07,-5.12854e-11,108621,33.7658], Tmin=(100,'K'), Tmax=(858.33,'K')), NASAPolynomial(coeffs=[3.00616,0.0466232,-2.22879e-05,4.19748e-09,-2.85273e-13,109297,28.9397], Tmin=(858.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(902.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]CC[CH][O](9286)',
    structure = SMILES('[CH2][C][CH]CC[CH][O]'),
    E0 = (902.231,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.66111,0.087637,-0.000149718,1.47406e-07,-5.50033e-11,108620,33.4447], Tmin=(100,'K'), Tmax=(859.156,'K')), NASAPolynomial(coeffs=[1.91266,0.0485763,-2.34989e-05,4.44195e-09,-3.02237e-13,109631,34.7349], Tmin=(859.156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(902.231,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(CCsJOH) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH][CH]CC[C][O](9287)',
    structure = SMILES('[CH2][CH][CH]CC[C][O]'),
    E0 = (929.152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,224.131,1059.8,1219.37,1463.24,1722.45,1919.59],'cm^-1')),
        HinderedRotor(inertia=(0.108962,'amu*angstrom^2'), symmetry=1, barrier=(2.77065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108962,'amu*angstrom^2'), symmetry=1, barrier=(2.77065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108962,'amu*angstrom^2'), symmetry=1, barrier=(2.77065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108962,'amu*angstrom^2'), symmetry=1, barrier=(2.77065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108962,'amu*angstrom^2'), symmetry=1, barrier=(2.77065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.978777,0.0818921,-0.000143324,1.47093e-07,-5.64659e-11,111845,34.3166], Tmin=(100,'K'), Tmax=(859.127,'K')), NASAPolynomial(coeffs=[-0.521752,0.0515315,-2.51096e-05,4.76226e-09,-3.24602e-13,113481,49.3494], Tmin=(859.127,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(929.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH][CH][CH]CC[CH][O](9288)',
    structure = SMILES('[CH][CH][CH]CC[CH][O]'),
    E0 = (891.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,202.198,1303.51,1308.71,1339.23,1397.29,1674.04,1731.03,1973.82],'cm^-1')),
        HinderedRotor(inertia=(0.135751,'amu*angstrom^2'), symmetry=1, barrier=(3.24714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135751,'amu*angstrom^2'), symmetry=1, barrier=(3.24714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135751,'amu*angstrom^2'), symmetry=1, barrier=(3.24714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135751,'amu*angstrom^2'), symmetry=1, barrier=(3.24714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135751,'amu*angstrom^2'), symmetry=1, barrier=(3.24714,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847206,0.086824,-0.000157733,1.62503e-07,-6.18251e-11,107311,34.8242], Tmin=(100,'K'), Tmax=(870.885,'K')), NASAPolynomial(coeffs=[-1.12401,0.0529925,-2.55976e-05,4.80862e-09,-3.24861e-13,109281,53.3987], Tmin=(870.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(891.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C1CCC1[O](5547)',
    structure = SMILES('[CH2][CH]C1CCC1[O]'),
    E0 = (414.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51773,0.0394184,2.88431e-05,-6.063e-08,2.43389e-11,50009.6,28.1078], Tmin=(100,'K'), Tmax=(1007.43,'K')), NASAPolynomial(coeffs=[12.879,0.0306433,-1.21923e-05,2.32672e-09,-1.68609e-13,45876.7,-35.9393], Tmin=(1007.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(414.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C=CCC[O](9289)',
    structure = SMILES('[CH2]C=C[CH]CC[O]'),
    E0 = (282.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47171,0.0567053,-3.24312e-05,8.813e-09,-9.70955e-13,34025.1,25.9299], Tmin=(100,'K'), Tmax=(1961.57,'K')), NASAPolynomial(coeffs=[13.4071,0.032367,-1.382e-05,2.48776e-09,-1.64815e-13,29342.6,-39.693], Tmin=(1961.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]CC[CH]C=O(5553)',
    structure = SMILES('[CH2][CH]CCC=C[O]'),
    E0 = (268.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793706,0.0616451,-3.51916e-05,2.96195e-09,3.05203e-12,32389.7,31.8677], Tmin=(100,'K'), Tmax=(1050.09,'K')), NASAPolynomial(coeffs=[13.5649,0.0288505,-1.0992e-05,1.97555e-09,-1.3594e-13,28833.5,-34.5322], Tmin=(1050.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]CC=O(5550)',
    structure = SMILES('[CH2][CH][CH][CH]CC=O'),
    E0 = (519.389,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,787.57,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.000861135,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000861135,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000861135,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000861135,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000861135,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.931,0.0394388,-1.53766e-05,1.79909e-09,7.86562e-15,62430.5,25.7142], Tmin=(100,'K'), Tmax=(2816.01,'K')), NASAPolynomial(coeffs=[46.9769,-0.00677023,1.2818e-06,-2.61189e-10,2.35649e-14,32265,-237.039], Tmin=(2816.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.389,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH]C=O(5551)',
    structure = SMILES('[CH2][CH][CH]CC=C[O]'),
    E0 = (462.739,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180,991.213,991.493],'cm^-1')),
        HinderedRotor(inertia=(0.162799,'amu*angstrom^2'), symmetry=1, barrier=(3.74306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162033,'amu*angstrom^2'), symmetry=1, barrier=(3.72546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00536787,'amu*angstrom^2'), symmetry=1, barrier=(3.73664,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00535731,'amu*angstrom^2'), symmetry=1, barrier=(3.73103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41551,0.0577798,-4.39124e-05,1.92192e-08,-3.66946e-12,55746.8,32.1999], Tmin=(100,'K'), Tmax=(1188.33,'K')), NASAPolynomial(coeffs=[8.37237,0.0343626,-1.43534e-05,2.63623e-09,-1.80753e-13,54093.4,-2.56335], Tmin=(1188.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][O](1556)',
    structure = SMILES('[CH2][CH][O]'),
    E0 = (367.955,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1667.93],'cm^-1')),
        HinderedRotor(inertia=(0.00517725,'amu*angstrom^2'), symmetry=1, barrier=(10.2207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (43.0446,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.22654,0.0212781,-3.59354e-05,3.91027e-08,-1.59281e-11,44278.5,12.2199], Tmin=(100,'K'), Tmax=(830.699,'K')), NASAPolynomial(coeffs=[2.17156,0.016018,-7.76586e-06,1.51127e-09,-1.05387e-13,44810.5,19.2612], Tmin=(830.699,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.955,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO) + radical(CCsJOH)"""),
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
    label = '[CH2][CH][CH][CH]C[CH][O](9290)',
    structure = SMILES('[CH2][CH][CH][CH]C[CH][O]'),
    E0 = (842.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,1140.6,1322.73,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0375057,'amu*angstrom^2'), symmetry=1, barrier=(3.32667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0375057,'amu*angstrom^2'), symmetry=1, barrier=(3.32667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0375057,'amu*angstrom^2'), symmetry=1, barrier=(3.32667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0375057,'amu*angstrom^2'), symmetry=1, barrier=(3.32667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0375057,'amu*angstrom^2'), symmetry=1, barrier=(3.32667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20972,0.084753,-0.000171827,1.90195e-07,-7.47476e-11,101458,36.8691], Tmin=(100,'K'), Tmax=(878.965,'K')), NASAPolynomial(coeffs=[-6.60553,0.0609239,-2.9801e-05,5.59489e-09,-3.76338e-13,105126,86.6163], Tmin=(878.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(842.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH][CH][O](9291)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH][O]'),
    E0 = (848.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,190.781,863.101,2429.18,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00408961,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00408961,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00408961,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00408961,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00408961,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1251,0.079365,-0.000140325,1.45978e-07,-5.62186e-11,102123,37.1056], Tmin=(100,'K'), Tmax=(867.76,'K')), NASAPolynomial(coeffs=[-1.66481,0.0524844,-2.51644e-05,4.72844e-09,-3.20107e-13,104103,58.7909], Tmin=(867.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(848.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(CCJCO) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C[CH][O](9292)',
    structure = SMILES('[CH2][CH]C[CH]C[CH][O]'),
    E0 = (648.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,476.327,863.886,2166.31,2892.04,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0854374,'amu*angstrom^2'), symmetry=1, barrier=(3.66409,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.944353,0.0844202,-0.000148429,1.54869e-07,-5.98592e-11,78085.5,35.2576], Tmin=(100,'K'), Tmax=(866.461,'K')), NASAPolynomial(coeffs=[-2.27782,0.0568642,-2.72688e-05,5.12854e-09,-3.47543e-13,80236.6,59.5321], Tmin=(866.461,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH]C[O](5541)',
    structure = SMILES('[CH2][CH][CH]C[CH]C[O]'),
    E0 = (668.066,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,181.576,182.104,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0104895,'amu*angstrom^2'), symmetry=1, barrier=(3.54191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104895,'amu*angstrom^2'), symmetry=1, barrier=(3.54191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104895,'amu*angstrom^2'), symmetry=1, barrier=(3.54191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104895,'amu*angstrom^2'), symmetry=1, barrier=(3.54191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0104895,'amu*angstrom^2'), symmetry=1, barrier=(3.54191,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37589,0.0708599,-0.00011005,1.15138e-07,-4.58013e-11,80431.6,36.6556], Tmin=(100,'K'), Tmax=(850.856,'K')), NASAPolynomial(coeffs=[-2.24094,0.0556084,-2.63003e-05,4.9646e-09,-3.3914e-13,82214.6,60.3816], Tmin=(850.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(668.066,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[CH][CH][O](9293)',
    structure = SMILES('[CH2][CH]CC[CH][CH][O]'),
    E0 = (653.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180,728.389,1599.52,2085.23,2456],'cm^-1')),
        HinderedRotor(inertia=(0.0955195,'amu*angstrom^2'), symmetry=1, barrier=(2.19618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0955195,'amu*angstrom^2'), symmetry=1, barrier=(2.19618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0955195,'amu*angstrom^2'), symmetry=1, barrier=(2.19618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0955195,'amu*angstrom^2'), symmetry=1, barrier=(2.19618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0955195,'amu*angstrom^2'), symmetry=1, barrier=(2.19618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.875858,0.0788361,-0.000116225,1.09747e-07,-4.09641e-11,78749.4,35.4368], Tmin=(100,'K'), Tmax=(840.803,'K')), NASAPolynomial(coeffs=[2.58569,0.0485602,-2.27122e-05,4.28128e-09,-2.92924e-13,79244.5,32.1385], Tmin=(840.803,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCJCO) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]CC[O](5542)',
    structure = SMILES('[CH2][CH][CH][CH]CC[O]'),
    E0 = (662.622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180.355,203.223,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0074524,'amu*angstrom^2'), symmetry=1, barrier=(6.35724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0074524,'amu*angstrom^2'), symmetry=1, barrier=(6.35724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0074524,'amu*angstrom^2'), symmetry=1, barrier=(6.35724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0074524,'amu*angstrom^2'), symmetry=1, barrier=(6.35724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0074524,'amu*angstrom^2'), symmetry=1, barrier=(6.35724,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45305,0.076336,-0.000141852,1.5971e-07,-6.44534e-11,79767.3,36.4458], Tmin=(100,'K'), Tmax=(868.908,'K')), NASAPolynomial(coeffs=[-7.13702,0.0639699,-3.08911e-05,5.82008e-09,-3.94452e-13,83219.7,87.9571], Tmin=(868.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH][CH]C[CH][O](9294)',
    structure = SMILES('[CH2]C[CH][CH]C[CH][O]'),
    E0 = (648.474,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,501.598,1042.96,2563.55,3421.48,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0468104,'amu*angstrom^2'), symmetry=1, barrier=(2.70507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468104,'amu*angstrom^2'), symmetry=1, barrier=(2.70507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468104,'amu*angstrom^2'), symmetry=1, barrier=(2.70507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468104,'amu*angstrom^2'), symmetry=1, barrier=(2.70507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0468104,'amu*angstrom^2'), symmetry=1, barrier=(2.70507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986805,0.0848315,-0.000153381,1.6292e-07,-6.35498e-11,78084.1,34.9336], Tmin=(100,'K'), Tmax=(866.494,'K')), NASAPolynomial(coeffs=[-3.37348,0.0588211,-2.84821e-05,5.37356e-09,-3.64554e-13,80571.8,65.3394], Tmin=(866.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.474,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJCC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH][CH]O(9295)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH]O'),
    E0 = (622.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.78218,0.0818669,-0.000127084,1.19506e-07,-4.34374e-11,74993.8,37.477], Tmin=(100,'K'), Tmax=(866.297,'K')), NASAPolynomial(coeffs=[3.00158,0.0461111,-2.10052e-05,3.88319e-09,-2.61484e-13,75566.4,32.6126], Tmin=(866.297,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(622.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(CCJCO) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C[CH][CH][O](9296)',
    structure = SMILES('[CH2]C[CH]C[CH][CH][O]'),
    E0 = (653.918,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,257.603,503.029,1218.01,2854.25,3616.77],'cm^-1')),
        HinderedRotor(inertia=(0.0382721,'amu*angstrom^2'), symmetry=1, barrier=(1.63197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382721,'amu*angstrom^2'), symmetry=1, barrier=(1.63197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382721,'amu*angstrom^2'), symmetry=1, barrier=(1.63197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382721,'amu*angstrom^2'), symmetry=1, barrier=(1.63197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382721,'amu*angstrom^2'), symmetry=1, barrier=(1.63197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.915065,0.0792885,-0.000121334,1.18021e-07,-4.47577e-11,78748.2,35.1242], Tmin=(100,'K'), Tmax=(843.862,'K')), NASAPolynomial(coeffs=[1.50004,0.0504994,-2.39148e-05,4.52374e-09,-3.09718e-13,79575.7,37.8899], Tmin=(843.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(653.918,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(CCJCO) + radical(RCCJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH][CH][CH]C[CH]O(5543)',
    structure = SMILES('[CH2][CH][CH][CH]C[CH]O'),
    E0 = (617.215,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,2219.03,2276.5,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0273209,'amu*angstrom^2'), symmetry=1, barrier=(9.41522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273209,'amu*angstrom^2'), symmetry=1, barrier=(9.41522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273209,'amu*angstrom^2'), symmetry=1, barrier=(9.41522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273209,'amu*angstrom^2'), symmetry=1, barrier=(9.41522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273209,'amu*angstrom^2'), symmetry=1, barrier=(9.41522,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0273209,'amu*angstrom^2'), symmetry=1, barrier=(9.41522,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.862284,0.087313,-0.000158813,1.64059e-07,-6.21293e-11,74329.4,37.2564], Tmin=(100,'K'), Tmax=(881.131,'K')), NASAPolynomial(coeffs=[-1.92992,0.054534,-2.56319e-05,4.74725e-09,-3.17513e-13,76585.9,60.3866], Tmin=(881.131,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.215,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH]C[CH][O](4352)',
    structure = SMILES('C[CH][CH][CH]C[CH][O]'),
    E0 = (637.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,541.296,1123.53,2543.35,3501.03,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0491993,'amu*angstrom^2'), symmetry=1, barrier=(2.19155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0491993,'amu*angstrom^2'), symmetry=1, barrier=(2.19155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0491993,'amu*angstrom^2'), symmetry=1, barrier=(2.19155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0491993,'amu*angstrom^2'), symmetry=1, barrier=(2.19155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0491993,'amu*angstrom^2'), symmetry=1, barrier=(2.19155,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17648,0.0839748,-0.000161238,1.7781e-07,-7.02847e-11,76775.2,35.2018], Tmin=(100,'K'), Tmax=(874.736,'K')), NASAPolynomial(coeffs=[-6.42606,0.0632652,-3.05972e-05,5.74418e-09,-3.87509e-13,80227.6,82.9935], Tmin=(874.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(637.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH][CH]C[CH][CH][O](4353)',
    structure = SMILES('C[CH][CH]C[CH][CH][O]'),
    E0 = (643.118,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,479.669,1006.77,2649.11,3416.77,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0475632,'amu*angstrom^2'), symmetry=1, barrier=(3.03516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0475632,'amu*angstrom^2'), symmetry=1, barrier=(3.03516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0475632,'amu*angstrom^2'), symmetry=1, barrier=(3.03516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0475632,'amu*angstrom^2'), symmetry=1, barrier=(3.03516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0475632,'amu*angstrom^2'), symmetry=1, barrier=(3.03516,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09501,0.0785492,-0.000129607,1.33436e-07,-5.16979e-11,77439.7,35.427], Tmin=(100,'K'), Tmax=(860.675,'K')), NASAPolynomial(coeffs=[-1.50328,0.0548571,-2.5979e-05,4.88213e-09,-3.31646e-13,79211.7,55.2685], Tmin=(860.675,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(643.118,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCJCO) + radical(CCsJOH)"""),
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
    E0 = (648.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (648.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (814.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (829.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1121.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1117.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1114.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1114.02,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1140.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1103.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (656.746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (711.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (711.862,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (738.955,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (682.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (679.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (698.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1023.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (1054.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1060.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (785.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (805.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (770.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (778.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (822.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (803.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (806.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (771.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (724.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (722.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]CC[CH][O](4354)'],
    products = ['vinoxy(1351)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH][CH]CC[CH][O](4354)'],
    products = ['[CH2][CH][CH]CCC=O(4258)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH][CH][CH2](5537)', '[CH2]C[CH][O](1549)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][O](1548)', '[CH2][CH][CH]C[CH2](509)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2][CH][C]CC[CH][O](9285)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH2][C][CH]CC[CH][O](9286)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][CH][CH]CC[C][O](9287)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH][CH][CH]CC[CH][O](9288)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH][CH]CC[CH][O](4354)'],
    products = ['[CH2][CH]C1CCC1[O](5547)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH][CH]CC[CH][O](4354)'],
    products = ['[CH2][CH]C=CCC[O](9289)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH][CH]CC[CH][O](4354)'],
    products = ['[CH2][CH]CC[CH]C=O(5553)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH2][CH][CH][CH]CC=O(5550)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2][CH][CH]C[CH]C=O(5551)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['vinoxy(1351)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.206e+07,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2][CH][CH][CH]C[CH][O](9290)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2][CH][CH]C[CH][CH][O](9291)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C[CH]C[CH][O](9292)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH][CH]C[CH]C[O](5541)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]CC[CH][CH][O](9293)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH][CH][CH]CC[O](5542)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C[CH][CH]C[CH][O](9294)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(867345,'s^-1'), n=1.96939, Ea=(174.054,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R3HJ;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH][CH]CC[CH][O](4354)'],
    products = ['[CH2][CH][CH]C[CH][CH]O(9295)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.30814e+09,'s^-1'), n=1.19923, Ea=(155.469,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3HJ;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C[CH]C[CH][CH][O](9296)'],
    products = ['[CH2][CH][CH]CC[CH][O](4354)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(645397,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_2;Y_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH]CC[CH][O](4354)'],
    products = ['[CH2][CH][CH][CH]C[CH]O(5543)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.04154e+06,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH]CC[CH][O](4354)'],
    products = ['C[CH][CH][CH]C[CH][O](4352)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.32767e+06,'s^-1'), n=1.53625, Ea=(75.6258,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH][CH]CC[CH][O](4354)'],
    products = ['C[CH][CH]C[CH][CH][O](4353)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2109',
    isomers = [
        '[CH2][CH][CH]CC[CH][O](4354)',
    ],
    reactants = [
        ('vinoxy(1351)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2109',
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


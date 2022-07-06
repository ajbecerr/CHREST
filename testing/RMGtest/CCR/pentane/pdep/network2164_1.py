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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3350.7,'J/mol'), sigma=(6.3658,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=523.37 K, Pc=29.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.80518,0.046418,-1.65621e-05,1.42179e-09,1.06398e-13,93558.7,28.3883], Tmin=(100,'K'), Tmax=(2731,'K')), NASAPolynomial(coeffs=[48.0514,-0.000118739,-1.03615e-06,8.12319e-11,4.88421e-15,62578.6,-242], Tmin=(2731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(778.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C[CH]C([CH2])C(9388)',
    structure = SMILES('[CH2][CH][CH]C[CH]C([CH2])C'),
    E0 = (748.656,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449419,0.0868474,-0.000103352,9.16657e-08,-3.41896e-11,90161.9,42.4706], Tmin=(100,'K'), Tmax=(830.571,'K')), NASAPolynomial(coeffs=[1.69735,0.0629525,-2.78981e-05,5.17626e-09,-3.52932e-13,90571.5,40.3954], Tmin=(830.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(748.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][C]C(C)C[CH][CH][CH2](9389)',
    structure = SMILES('[CH2][C]C(C)C[CH][CH][CH2]'),
    E0 = (1005.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,255.381,869.615,1159.49,1517.94,1882.82],'cm^-1')),
        HinderedRotor(inertia=(0.121689,'amu*angstrom^2'), symmetry=1, barrier=(3.33932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121689,'amu*angstrom^2'), symmetry=1, barrier=(3.33932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121689,'amu*angstrom^2'), symmetry=1, barrier=(3.33932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121689,'amu*angstrom^2'), symmetry=1, barrier=(3.33932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121689,'amu*angstrom^2'), symmetry=1, barrier=(3.33932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121689,'amu*angstrom^2'), symmetry=1, barrier=(3.33932,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121689,'amu*angstrom^2'), symmetry=1, barrier=(3.33932,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.262892,0.0897214,-0.000107865,9.10002e-08,-3.27871e-11,121102,40.4951], Tmin=(100,'K'), Tmax=(803.701,'K')), NASAPolynomial(coeffs=[4.35326,0.0576062,-2.59829e-05,4.87685e-09,-3.35474e-13,120824,24.0177], Tmin=(803.701,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1005.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH][C]CC(C)[CH][CH2](9390)',
    structure = SMILES('[CH2][CH][C]CC(C)[CH][CH2]'),
    E0 = (1005.92,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,208.088,808.088,1077.45,1360.29,1643.99],'cm^-1')),
        HinderedRotor(inertia=(0.143306,'amu*angstrom^2'), symmetry=1, barrier=(3.56675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143306,'amu*angstrom^2'), symmetry=1, barrier=(3.56675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143306,'amu*angstrom^2'), symmetry=1, barrier=(3.56675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143306,'amu*angstrom^2'), symmetry=1, barrier=(3.56675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143306,'amu*angstrom^2'), symmetry=1, barrier=(3.56675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143306,'amu*angstrom^2'), symmetry=1, barrier=(3.56675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.143306,'amu*angstrom^2'), symmetry=1, barrier=(3.56675,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.75421,0.0800211,-5.84409e-05,6.96428e-09,1.50009e-11,121093,38.8077], Tmin=(100,'K'), Tmax=(564.922,'K')), NASAPolynomial(coeffs=[5.81553,0.0559096,-2.5554e-05,4.89663e-09,-3.44002e-13,120335,15.6241], Tmin=(564.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1005.92,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]CC(C)[CH][CH2](9391)',
    structure = SMILES('[CH2][C][CH]CC(C)[CH][CH2]'),
    E0 = (1005.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,245.851,857.419,1143.22,1486.17,1844.99],'cm^-1')),
        HinderedRotor(inertia=(0.12497,'amu*angstrom^2'), symmetry=1, barrier=(3.42922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12497,'amu*angstrom^2'), symmetry=1, barrier=(3.42922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12497,'amu*angstrom^2'), symmetry=1, barrier=(3.42922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12497,'amu*angstrom^2'), symmetry=1, barrier=(3.42922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12497,'amu*angstrom^2'), symmetry=1, barrier=(3.42922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12497,'amu*angstrom^2'), symmetry=1, barrier=(3.42922,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12497,'amu*angstrom^2'), symmetry=1, barrier=(3.42922,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42129,0.0693378,2.90918e-06,-1.44175e-07,1.4277e-10,121066,36.4325], Tmin=(100,'K'), Tmax=(444.669,'K')), NASAPolynomial(coeffs=[5.08897,0.0571774,-2.63429e-05,5.03661e-09,-3.5199e-13,120534,19.3946], Tmin=(444.669,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1005.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C(C)C[CH][CH][CH2](9392)',
    structure = SMILES('[CH][CH]C(C)C[CH][CH][CH2]'),
    E0 = (995.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180.746,372.15,707.666,1911.72,2514.28,3528.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.40976,0.087886,-0.000109436,9.81733e-08,-3.70162e-11,119808,41.8341], Tmin=(100,'K'), Tmax=(809.584,'K')), NASAPolynomial(coeffs=[2.46238,0.0608317,-2.79737e-05,5.28758e-09,-3.64684e-13,120030,35.7881], Tmin=(809.584,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(995.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH][CH]CC(C)[CH][CH2](9393)',
    structure = SMILES('[CH][CH][CH]CC(C)[CH][CH2]'),
    E0 = (995.136,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180.746,372.15,707.666,1911.72,2514.28,3528.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0739175,'amu*angstrom^2'), symmetry=1, barrier=(1.72502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.40976,0.087886,-0.000109436,9.81733e-08,-3.70162e-11,119808,41.8341], Tmin=(100,'K'), Tmax=(809.584,'K')), NASAPolynomial(coeffs=[2.46238,0.0608317,-2.79737e-05,5.28758e-09,-3.64684e-13,120030,35.7881], Tmin=(809.584,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(995.136,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C1CC(C)C1[CH2](5796)',
    structure = SMILES('[CH2][CH]C1CC(C)C1[CH2]'),
    E0 = (496.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759399,0.0529585,3.40269e-05,-7.58905e-08,3.12513e-11,59847.9,33.6798], Tmin=(100,'K'), Tmax=(985.818,'K')), NASAPolynomial(coeffs=[14.4708,0.0413645,-1.5343e-05,2.81298e-09,-1.99784e-13,55004.5,-43.1277], Tmin=(985.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C](C)CCC=C(5815)',
    structure = SMILES('[CH2][CH][C](C)CCC=C'),
    E0 = (471.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,180,400,1436.2,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0476616,'amu*angstrom^2'), symmetry=1, barrier=(4.38333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0476616,'amu*angstrom^2'), symmetry=1, barrier=(4.38333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0476616,'amu*angstrom^2'), symmetry=1, barrier=(4.38333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0476616,'amu*angstrom^2'), symmetry=1, barrier=(4.38333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0476616,'amu*angstrom^2'), symmetry=1, barrier=(4.38333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0476616,'amu*angstrom^2'), symmetry=1, barrier=(4.38333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66792,0.0632872,2.71484e-05,-1.77626e-07,1.63521e-10,56746.2,32.9186], Tmin=(100,'K'), Tmax=(425.543,'K')), NASAPolynomial(coeffs=[3.68018,0.061732,-2.8561e-05,5.5132e-09,-3.89053e-13,56417.7,23.0831], Tmin=(425.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)C[CH2](5805)',
    structure = SMILES('[CH2]C=C[CH]C(C)C[CH2]'),
    E0 = (368.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0341853,0.0739304,-2.51936e-05,-1.41905e-08,9.06e-12,44447.4,33.5372], Tmin=(100,'K'), Tmax=(1071.69,'K')), NASAPolynomial(coeffs=[15.4344,0.0422402,-1.69355e-05,3.12768e-09,-2.1825e-13,39665.6,-48.7359], Tmin=(1071.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.276,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][C](C)CC=C[CH2](5798)',
    structure = SMILES('[CH2][CH][C](C)CC=C[CH2]'),
    E0 = (607.128,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,184.098,405.743,2978.29],'cm^-1')),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0922303,'amu*angstrom^2'), symmetry=1, barrier=(2.20216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793825,0.077034,-6.20887e-05,3.37689e-08,-9.00501e-12,73130.4,35.6306], Tmin=(100,'K'), Tmax=(828.936,'K')), NASAPolynomial(coeffs=[4.87591,0.0573368,-2.64468e-05,5.10518e-09,-3.60577e-13,72453.6,16.7027], Tmin=(828.936,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(607.128,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C(C)[CH]C=C[CH2](5799)',
    structure = SMILES('[CH2][CH]C(C)[CH]C=C[CH2]'),
    E0 = (562.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.162606,0.0732669,-3.76535e-05,4.08843e-09,1.69225e-12,67838.6,35.6441], Tmin=(100,'K'), Tmax=(1226.2,'K')), NASAPolynomial(coeffs=[14.9314,0.0406437,-1.67731e-05,3.08096e-09,-2.11488e-13,63047.4,-43.3863], Tmin=(1226.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Cs_S) + radical(RCCJ) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]CC=C[CH2](5280)',
    structure = SMILES('[CH2][CH][CH]CC=C[CH2]'),
    E0 = (645.542,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,239.279,1036.11,2991.11],'cm^-1')),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0187629,'amu*angstrom^2'), symmetry=1, barrier=(0.702078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36806,0.0623164,-5.56022e-05,4.06775e-08,-1.47416e-11,77731.4,33.1934], Tmin=(100,'K'), Tmax=(737.489,'K')), NASAPolynomial(coeffs=[3.19352,0.0490934,-2.19507e-05,4.14967e-09,-2.88591e-13,77552.5,25.555], Tmin=(737.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.542,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(Allyl_P)"""),
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
    label = '[CH2][CH][CH]C[CH][CH][CH2](9244)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH][CH2]'),
    E0 = (972.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1137.5,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00835278,'amu*angstrom^2'), symmetry=1, barrier=(14.3892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00835278,'amu*angstrom^2'), symmetry=1, barrier=(14.3892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00835278,'amu*angstrom^2'), symmetry=1, barrier=(14.3892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00835278,'amu*angstrom^2'), symmetry=1, barrier=(14.3892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00835278,'amu*angstrom^2'), symmetry=1, barrier=(14.3892,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00835278,'amu*angstrom^2'), symmetry=1, barrier=(14.3892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39459,0.0752035,-0.000128484,1.41026e-07,-5.64093e-11,117060,39.1657], Tmin=(100,'K'), Tmax=(871.57,'K')), NASAPolynomial(coeffs=[-5.8884,0.0631706,-2.95412e-05,5.50341e-09,-3.71211e-13,120056,83.2057], Tmin=(871.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(972.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C](C)[CH][CH2](9394)',
    structure = SMILES('[CH2][CH][CH]C[C](C)[CH][CH2]'),
    E0 = (937.59,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3426.73,3880.96,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0139827,'amu*angstrom^2'), symmetry=1, barrier=(1.03579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139827,'amu*angstrom^2'), symmetry=1, barrier=(1.03579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139827,'amu*angstrom^2'), symmetry=1, barrier=(1.03579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139827,'amu*angstrom^2'), symmetry=1, barrier=(1.03579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139827,'amu*angstrom^2'), symmetry=1, barrier=(1.03579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139827,'amu*angstrom^2'), symmetry=1, barrier=(1.03579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139827,'amu*angstrom^2'), symmetry=1, barrier=(1.03579,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.5314,0.0511614,-1.77143e-05,9.26878e-10,2.30097e-13,112681,28.5725], Tmin=(100,'K'), Tmax=(2751.26,'K')), NASAPolynomial(coeffs=[65.645,-0.0114589,2.12509e-06,-4.15002e-10,3.71314e-14,69125.6,-346.158], Tmin=(2751.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(937.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C(C)[CH][CH2](9395)',
    structure = SMILES('[CH2][CH][CH][CH]C(C)[CH][CH2]'),
    E0 = (946.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2153.03,2663.14,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0164675,'amu*angstrom^2'), symmetry=1, barrier=(9.73051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164675,'amu*angstrom^2'), symmetry=1, barrier=(9.73051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164675,'amu*angstrom^2'), symmetry=1, barrier=(9.73051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164675,'amu*angstrom^2'), symmetry=1, barrier=(9.73051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164675,'amu*angstrom^2'), symmetry=1, barrier=(9.73051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164675,'amu*angstrom^2'), symmetry=1, barrier=(9.73051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0164675,'amu*angstrom^2'), symmetry=1, barrier=(9.73051,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.33324,0.056106,-2.31469e-05,3.21881e-09,-7.52165e-14,113844,34.4372], Tmin=(100,'K'), Tmax=(2747.09,'K')), NASAPolynomial(coeffs=[53.2504,0.00017397,-1.75317e-06,2.46763e-10,-6.75789e-15,80097.4,-268.332], Tmin=(2747.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(946.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])[CH][CH2](9025)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])[CH][CH2]'),
    E0 = (957.25,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2680.8,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0680509,'amu*angstrom^2'), symmetry=1, barrier=(3.59777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465157,0.0878423,-0.000114739,1.05124e-07,-3.91147e-11,115248,44.2], Tmin=(100,'K'), Tmax=(844.606,'K')), NASAPolynomial(coeffs=[1.5906,0.0604829,-2.70259e-05,5.00868e-09,-3.40222e-13,115844,43.612], Tmin=(844.606,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(957.25,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C](C)C[CH2](9396)',
    structure = SMILES('[CH2][CH][CH]C[C](C)C[CH2]'),
    E0 = (743.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.89682,0.0574083,-2.29518e-05,3.05491e-09,-6.47913e-14,89312.5,28.3062], Tmin=(100,'K'), Tmax=(2908.94,'K')), NASAPolynomial(coeffs=[72.0075,-0.0159037,4.36051e-06,-8.00163e-10,5.98877e-14,41078.5,-387.812], Tmin=(2908.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C(C)[CH][CH2](5788)',
    structure = SMILES('[CH2][CH]C[CH]C(C)[CH][CH2]'),
    E0 = (752.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,499.065,1113.29,3597.27,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0347676,'amu*angstrom^2'), symmetry=1, barrier=(1.59157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0347676,'amu*angstrom^2'), symmetry=1, barrier=(1.59157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0347676,'amu*angstrom^2'), symmetry=1, barrier=(1.59157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0347676,'amu*angstrom^2'), symmetry=1, barrier=(1.59157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0347676,'amu*angstrom^2'), symmetry=1, barrier=(1.59157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0347676,'amu*angstrom^2'), symmetry=1, barrier=(1.59157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0347676,'amu*angstrom^2'), symmetry=1, barrier=(1.59157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.551234,0.0824473,-8.26283e-05,6.60518e-08,-2.49723e-11,90592.9,42.107], Tmin=(100,'K'), Tmax=(734.092,'K')), NASAPolynomial(coeffs=[3.07435,0.0624105,-2.88365e-05,5.53113e-09,-3.87501e-13,90391.9,31.8684], Tmin=(734.092,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C(C)C[CH2](9397)',
    structure = SMILES('[CH2][CH][CH][CH]C(C)C[CH2]'),
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
    label = '[CH2][CH][CH]CC([CH2])C[CH2](9085)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])C[CH2]'),
    E0 = (762.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,187.146,1307.08,3863.53,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270552,'amu*angstrom^2'), symmetry=1, barrier=(4.33074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.294331,0.0887901,-0.00010221,8.5135e-08,-3.02778e-11,91858.9,42.259], Tmin=(100,'K'), Tmax=(829.918,'K')), NASAPolynomial(coeffs=[3.60694,0.0596931,-2.58869e-05,4.76055e-09,-3.23299e-13,91761.3,29.6196], Tmin=(829.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[C](C)[CH][CH2](5789)',
    structure = SMILES('[CH2][CH]CC[C](C)[CH][CH2]'),
    E0 = (743.132,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2541.55,2888.39,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0232515,'amu*angstrom^2'), symmetry=1, barrier=(11.9542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232515,'amu*angstrom^2'), symmetry=1, barrier=(11.9542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232515,'amu*angstrom^2'), symmetry=1, barrier=(11.9542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232515,'amu*angstrom^2'), symmetry=1, barrier=(11.9542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232515,'amu*angstrom^2'), symmetry=1, barrier=(11.9542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232515,'amu*angstrom^2'), symmetry=1, barrier=(11.9542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0232515,'amu*angstrom^2'), symmetry=1, barrier=(11.9542,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.308644,0.0928237,-0.000124644,1.18765e-07,-4.59374e-11,89499.7,41.5482], Tmin=(100,'K'), Tmax=(823.832,'K')), NASAPolynomial(coeffs=[0.500066,0.0666684,-3.10916e-05,5.89169e-09,-4.05836e-13,90324.2,45.8573], Tmin=(823.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.132,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Cs_S) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C(C)[CH][CH]C[CH2](9398)',
    structure = SMILES('[CH2][CH]C(C)[CH][CH]C[CH2]'),
    E0 = (752.263,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,662.479,1014.23,3766.29,4000],'cm^-1')),
        HinderedRotor(inertia=(0.012913,'amu*angstrom^2'), symmetry=1, barrier=(0.632687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012913,'amu*angstrom^2'), symmetry=1, barrier=(0.632687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012913,'amu*angstrom^2'), symmetry=1, barrier=(0.632687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012913,'amu*angstrom^2'), symmetry=1, barrier=(0.632687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012913,'amu*angstrom^2'), symmetry=1, barrier=(0.632687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012913,'amu*angstrom^2'), symmetry=1, barrier=(0.632687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.012913,'amu*angstrom^2'), symmetry=1, barrier=(0.632687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33648,0.064549,-3.08503e-05,6.18665e-09,-4.621e-13,90505.3,35.3349], Tmin=(100,'K'), Tmax=(3254.72,'K')), NASAPolynomial(coeffs=[50.3594,0.00552885,-3.6494e-06,6.15001e-10,-3.41282e-14,59245.5,-253.02], Tmin=(3254.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.263,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CCC([CH2])[CH][CH2](5791)',
    structure = SMILES('[CH2][CH]CCC([CH2])[CH][CH2]'),
    E0 = (762.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,208.296,1250.46,2018.91,4000],'cm^-1')),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112012,'amu*angstrom^2'), symmetry=1, barrier=(2.57537,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10531,0.0737795,-2.35e-05,-6.32665e-08,6.5968e-11,91836.8,39.5049], Tmin=(100,'K'), Tmax=(503.992,'K')), NASAPolynomial(coeffs=[5.23185,0.0577035,-2.5282e-05,4.73755e-09,-3.28201e-13,91209.1,20.3232], Tmin=(503.992,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][CH][C](C)C[CH]C[CH2](5790)',
    structure = SMILES('[CH2][CH][C](C)C[CH]C[CH2]'),
    E0 = (743.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,960.021,2568.19,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0328356,'amu*angstrom^2'), symmetry=1, barrier=(5.68989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0328356,'amu*angstrom^2'), symmetry=1, barrier=(5.68989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0328356,'amu*angstrom^2'), symmetry=1, barrier=(5.68989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0328356,'amu*angstrom^2'), symmetry=1, barrier=(5.68989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0328356,'amu*angstrom^2'), symmetry=1, barrier=(5.68989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0328356,'amu*angstrom^2'), symmetry=1, barrier=(5.68989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0328356,'amu*angstrom^2'), symmetry=1, barrier=(5.68989,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.52071,0.0597045,-2.55905e-05,3.99509e-09,-1.75011e-13,89344.1,29.5266], Tmin=(100,'K'), Tmax=(2913.68,'K')), NASAPolynomial(coeffs=[68.4836,-0.0125133,2.85413e-06,-5.14814e-10,4.04827e-14,44286.4,-365.712], Tmin=(2913.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)[CH][CH][CH]C(4546)',
    structure = SMILES('[CH2][CH]C(C)[CH][CH][CH]C'),
    E0 = (741.463,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,2156.61,2460.57,3867.98,4000],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.88834,0.0594493,-2.43178e-05,3.32932e-09,-5.99579e-14,89181.6,34.3076], Tmin=(100,'K'), Tmax=(2631.57,'K')), NASAPolynomial(coeffs=[45.7947,0.0102635,-5.42006e-06,8.56886e-10,-4.5003e-14,61048.2,-224.753], Tmin=(2631.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH]C[CH2](5792)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH]C[CH2]'),
    E0 = (762.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1007.55,2415.44,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0650668,'amu*angstrom^2'), symmetry=1, barrier=(6.02803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331576,0.0867925,-9.20103e-05,7.18669e-08,-2.5212e-11,91870.3,41.9491], Tmin=(100,'K'), Tmax=(780.083,'K')), NASAPolynomial(coeffs=[4.49909,0.0589563,-2.60505e-05,4.87045e-09,-3.35455e-13,91416.9,24.1394], Tmin=(780.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C](C)C[CH][CH]C(4545)',
    structure = SMILES('[CH2][CH][C](C)C[CH][CH]C'),
    E0 = (732.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,1830.52,2847.58,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0186758,'amu*angstrom^2'), symmetry=1, barrier=(9.90169,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.08413,0.0545227,-1.89188e-05,1.05859e-09,2.41154e-13,88019.5,28.4524], Tmin=(100,'K'), Tmax=(2678.98,'K')), NASAPolynomial(coeffs=[58.2003,-0.00143803,-1.49401e-06,1.83565e-10,-1.80944e-16,50110.4,-302.591], Tmin=(2678.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH][CH]C(4547)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH][CH]C'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449419,0.0868474,-0.000103352,9.16657e-08,-3.41896e-11,90564.5,42.4706], Tmin=(100,'K'), Tmax=(830.571,'K')), NASAPolynomial(coeffs=[1.69735,0.0629525,-2.78981e-05,5.17626e-09,-3.52932e-13,90974.1,40.3954], Tmin=(830.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (752.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (752.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1196.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (895.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (946.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (917.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1210.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1214.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1217.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1217.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1217.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1206.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1206.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (760.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (815.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (815.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (826.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (782.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (813.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (815.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (761.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (1105.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1109.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (1149.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (1158.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1169.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (917.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (889.255,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (868.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (880.647,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (868.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (901.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (926.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (838.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (883.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (904.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (827.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (836.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (825.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (804.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH]C=C(3743)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH]C(C)C[CH]C=C(4544)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', '[CH2][CH][CH]CC[CH][CH2](4163)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH]C(4576)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH][CH]C[CH]C([CH2])C(9388)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH][CH2](4513)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH2](502)', '[CH2][CH][CH]C[CH]C(556)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][CH][CH2](5537)', '[CH2][CH]C([CH2])C(111)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][C]C(C)C[CH][CH][CH2](9389)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][CH][C]CC(C)[CH][CH2](9390)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][C][CH]CC(C)[CH][CH2](9391)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH][CH]C(C)C[CH][CH][CH2](9392)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH][CH][CH]CC(C)[CH][CH2](9393)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH]C1CC(C)C1[CH2](5796)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH][C](C)CCC=C(5815)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2]C=C[CH]C(C)C[CH2](5805)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2][CH][C](C)CC=C[CH2](5798)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2][CH]C(C)[CH]C=C[CH2](5799)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'm1_allyl(186)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['CH3(17)', '[CH2][CH][CH]CC=C[CH2](5280)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.139979,'m^3/(mol*s)'), n=2.09962, Ea=(33.817,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C=C(3743)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH][CH][CH2](5531)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH3(17)', '[CH2][CH][CH]C[CH][CH][CH2](9244)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
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
    reactants = ['H(3)', '[CH2][CH][CH]C[C](C)[CH][CH2](9394)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][CH][CH][CH]C(C)[CH][CH2](9395)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH2][CH][CH]CC([CH2])[CH][CH2](9025)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH][CH]C[C](C)C[CH2](9396)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]C[CH]C(C)[CH][CH2](5788)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH][CH]C(C)C[CH2](9397)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH][CH]CC([CH2])C[CH2](9085)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH]CC[C](C)[CH][CH2](5789)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH][CH]C[C](C)[CH]C(9273)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.081e+09,'s^-1'), n=0.812344, Ea=(148.932,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_Cs2] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C(C)[CH][CH]C[CH2](9398)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(867345,'s^-1'), n=1.96939, Ea=(174.054,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R3HJ;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH]CCC([CH2])[CH][CH2](5791)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(927.918,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH][CH][CH]C(C)[CH]C(9274)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH][C](C)C[CH]C[CH2](5790)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(322699,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_1;Y_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH]C(C)[CH][CH][CH]C(4546)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(754000,'s^-1'), n=1.63, Ea=(74.8936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]C([CH2])C[CH]C[CH2](5792)'],
    products = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH][C](C)C[CH][CH]C(4545)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH][CH]CC(C)[CH][CH2](4548)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2164',
    isomers = [
        '[CH2][CH][CH]CC(C)[CH][CH2](4548)',
    ],
    reactants = [
        ('[CH2][CH]C=C(3743)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2164',
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


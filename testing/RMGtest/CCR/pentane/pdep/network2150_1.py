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
    label = '[CH2]C([CH]C)C([CH2])C=C(4561)',
    structure = SMILES('[CH2]C([CH]C)C([CH2])C=C'),
    E0 = (482.745,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.133767,0.0814542,-5.70154e-05,2.18831e-08,-3.49527e-12,58217.4,39.4956], Tmin=(100,'K'), Tmax=(1454.32,'K')), NASAPolynomial(coeffs=[14.7768,0.0404443,-1.47177e-05,2.49381e-09,-1.62254e-13,53880.4,-38.0237], Tmin=(1454.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.745,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.78,'J/mol'), sigma=(6.74357,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.91 K, Pc=26.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612925,0.0786769,-5.55307e-05,2.26583e-08,-4.17915e-12,91590.1,40.6489], Tmin=(100,'K'), Tmax=(1183.11,'K')), NASAPolynomial(coeffs=[8.22032,0.0529569,-2.29219e-05,4.28369e-09,-2.96467e-13,89790,2.66847], Tmin=(1183.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(760.54,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449419,0.0868474,-0.000103352,9.16657e-08,-3.41896e-11,90564.5,42.4706], Tmin=(100,'K'), Tmax=(830.571,'K')), NASAPolynomial(coeffs=[1.69735,0.0629525,-2.78981e-05,5.17626e-09,-3.52932e-13,90974.1,40.3954], Tmin=(830.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C([CH2])[CH][CH2](9040)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][CH2]'),
    E0 = (810.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1610.44,1610.45],'cm^-1')),
        HinderedRotor(inertia=(0.0851314,'amu*angstrom^2'), symmetry=1, barrier=(5.99533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0032575,'amu*angstrom^2'), symmetry=1, barrier=(5.99512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(6.49989e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0851266,'amu*angstrom^2'), symmetry=1, barrier=(5.99524,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0336839,'amu*angstrom^2'), symmetry=1, barrier=(61.9916,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.75127,0.0531098,-4.22739e-05,2.36613e-08,-6.51376e-12,97552.1,31.6731], Tmin=(100,'K'), Tmax=(809.517,'K')), NASAPolynomial(coeffs=[4.44911,0.0397791,-1.75725e-05,3.3186e-09,-2.31342e-13,97115.3,19.2277], Tmin=(809.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(810.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C([CH2])[CH]C(4150)',
    structure = SMILES('[CH2][CH][CH]C([CH2])[CH]C'),
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
    label = '[CH]C([CH2])C([CH2])[CH][CH2](9250)',
    structure = SMILES('[CH]C([CH2])C([CH2])[CH][CH2]'),
    E0 = (1034.48,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200.087,1089.31,1317.82,1952.81,2130.09],'cm^-1')),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.120637,'amu*angstrom^2'), symmetry=1, barrier=(2.92085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.632966,0.0733946,-6.67761e-05,3.57233e-08,-8.02303e-12,124541,36.7755], Tmin=(100,'K'), Tmax=(1058.34,'K')), NASAPolynomial(coeffs=[10.4918,0.0361322,-1.39622e-05,2.45405e-09,-1.6401e-13,122455,-11.3465], Tmin=(1058.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1034.48,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][CH]C([CH2])C([CH2])[C]C(9431)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])[C]C'),
    E0 = (1014.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0167561,0.0878076,-7.92842e-05,4.20226e-08,-9.40792e-12,122105,40.1353], Tmin=(100,'K'), Tmax=(1056.49,'K')), NASAPolynomial(coeffs=[11.4691,0.0444465,-1.77188e-05,3.17252e-09,-2.14481e-13,119685,-15.7446], Tmin=(1056.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1014.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]C([CH2])C([CH2])[CH]C(9432)',
    structure = SMILES('[CH2][C]C([CH2])C([CH2])[CH]C'),
    E0 = (1014.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0167561,0.0878076,-7.92842e-05,4.20226e-08,-9.40792e-12,122105,40.1353], Tmin=(100,'K'), Tmax=(1056.49,'K')), NASAPolynomial(coeffs=[11.4691,0.0444465,-1.77188e-05,3.17252e-09,-2.14481e-13,119685,-15.7446], Tmin=(1056.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1014.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH][CH2])C([CH2])[CH]C(9433)',
    structure = SMILES('[CH]C([CH][CH2])C([CH2])[CH]C'),
    E0 = (1003.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,199.93,1187.47,1248.37,1415.9,1779.76,1881.31],'cm^-1')),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465606,0.0815819,-6.62652e-05,3.18869e-08,-6.79327e-12,120818,40.4027], Tmin=(100,'K'), Tmax=(1068.08,'K')), NASAPolynomial(coeffs=[8.93792,0.0498525,-2.17044e-05,4.07301e-09,-2.82976e-13,119008,-1.02931], Tmin=(1068.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1003.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH]C)C([CH2])[CH][CH2](9434)',
    structure = SMILES('[CH]C([CH]C)C([CH2])[CH][CH2]'),
    E0 = (1003.51,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,199.93,1187.47,1248.37,1415.9,1779.76,1881.31],'cm^-1')),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127945,'amu*angstrom^2'), symmetry=1, barrier=(3.10521,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465606,0.0815819,-6.62652e-05,3.18869e-08,-6.79327e-12,120818,40.4027], Tmin=(100,'K'), Tmax=(1068.08,'K')), NASAPolynomial(coeffs=[8.93792,0.0498525,-2.17044e-05,4.07301e-09,-2.82976e-13,119008,-1.02931], Tmin=(1068.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1003.51,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C([CH2])C([CH2])[CH]C(9435)',
    structure = SMILES('[CH][CH]C([CH2])C([CH2])[CH]C'),
    E0 = (1003.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,197.403,1028.56,1191.03,1421.5,1788.08,2065.74],'cm^-1')),
        HinderedRotor(inertia=(0.13402,'amu*angstrom^2'), symmetry=1, barrier=(3.09,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13402,'amu*angstrom^2'), symmetry=1, barrier=(3.09,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13402,'amu*angstrom^2'), symmetry=1, barrier=(3.09,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13402,'amu*angstrom^2'), symmetry=1, barrier=(3.09,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13402,'amu*angstrom^2'), symmetry=1, barrier=(3.09,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13402,'amu*angstrom^2'), symmetry=1, barrier=(3.09,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13402,'amu*angstrom^2'), symmetry=1, barrier=(3.09,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.28169,0.0844634,-7.50192e-05,4.07941e-08,-9.66034e-12,120806,41.0585], Tmin=(100,'K'), Tmax=(986.136,'K')), NASAPolynomial(coeffs=[9.17874,0.0483747,-2.01247e-05,3.68305e-09,-2.5208e-13,119051,-1.74038], Tmin=(986.136,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1003.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C([CH2])C1CC1C(9436)',
    structure = SMILES('[CH2][CH]C([CH2])C1CC1C'),
    E0 = (504.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.461754,0.0621012,7.29905e-06,-4.94766e-08,2.23866e-11,60790.4,35.4269], Tmin=(100,'K'), Tmax=(990.296,'K')), NASAPolynomial(coeffs=[14.9364,0.0401759,-1.48391e-05,2.68728e-09,-1.88554e-13,56131.8,-43.3105], Tmin=(990.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)C1CC1[CH2](5830)',
    structure = SMILES('[CH2]C([CH]C)C1CC1[CH2]'),
    E0 = (504.101,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,185.928,885.219,885.219,885.219,885.219,885.219,885.219,885.219,885.219,885.219,885.219,885.219,2400.49],'cm^-1')),
        HinderedRotor(inertia=(0.0162372,'amu*angstrom^2'), symmetry=1, barrier=(0.425339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162372,'amu*angstrom^2'), symmetry=1, barrier=(0.425339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162372,'amu*angstrom^2'), symmetry=1, barrier=(0.425339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162372,'amu*angstrom^2'), symmetry=1, barrier=(0.425339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0162372,'amu*angstrom^2'), symmetry=1, barrier=(0.425339,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.457387,0.0628855,5.8009e-06,-4.99114e-08,2.34957e-11,60770.4,35.4372], Tmin=(100,'K'), Tmax=(963.119,'K')), NASAPolynomial(coeffs=[14.8192,0.0392875,-1.35906e-05,2.37394e-09,-1.63909e-13,56332,-41.9902], Tmin=(963.119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(504.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1C([CH2])C(C)C1[CH2](5769)',
    structure = SMILES('[CH2]C1C([CH2])C(C)C1[CH2]'),
    E0 = (497.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647295,0.053475,4.39686e-05,-9.79802e-08,4.32613e-11,60014.9,32.3306], Tmin=(100,'K'), Tmax=(924.482,'K')), NASAPolynomial(coeffs=[17.0475,0.0345206,-9.65694e-06,1.53858e-09,-1.05379e-13,54760.2,-57.5223], Tmin=(924.482,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]C1CC([CH]C)C1[CH2](5831)',
    structure = SMILES('[CH2]C1CC([CH]C)C1[CH2]'),
    E0 = (496.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.748867,0.0538163,3.22692e-05,-7.59884e-08,3.22177e-11,59828.2,33.7121], Tmin=(100,'K'), Tmax=(964.146,'K')), NASAPolynomial(coeffs=[14.3818,0.040428,-1.40668e-05,2.49312e-09,-1.74599e-13,55192.8,-41.9672], Tmin=(964.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.334,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1CCC1[CH]C(5891)',
    structure = SMILES('[CH2][CH]C1CCC1[CH]C'),
    E0 = (495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.865289,0.0531297,2.35527e-05,-5.65402e-08,2.18621e-11,59661.1,35.0531], Tmin=(100,'K'), Tmax=(1054.79,'K')), NASAPolynomial(coeffs=[12.3798,0.0463409,-1.92348e-05,3.64834e-09,-2.59318e-13,55180.6,-30.8362], Tmin=(1054.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](CC)C([CH2])C=C(5090)',
    structure = SMILES('[CH2][C](CC)C([CH2])C=C'),
    E0 = (473.626,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,207.095,814.191,1656.76],'cm^-1')),
        HinderedRotor(inertia=(0.148125,'amu*angstrom^2'), symmetry=1, barrier=(3.57273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148125,'amu*angstrom^2'), symmetry=1, barrier=(3.57273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148125,'amu*angstrom^2'), symmetry=1, barrier=(3.57273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148125,'amu*angstrom^2'), symmetry=1, barrier=(3.57273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148125,'amu*angstrom^2'), symmetry=1, barrier=(3.57273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148125,'amu*angstrom^2'), symmetry=1, barrier=(3.57273,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.288161,0.0838988,-7.10126e-05,3.82508e-08,-9.1416e-12,57095.7,36.561], Tmin=(100,'K'), Tmax=(970.734,'K')), NASAPolynomial(coeffs=[8.18002,0.051379,-2.07611e-05,3.73913e-09,-2.5337e-13,55563.6,-1.27801], Tmin=(970.734,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2]CC(=C)C([CH2])[CH]C(9437)',
    structure = SMILES('[CH2]CC(=C)C([CH2])[CH]C'),
    E0 = (475.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0278324,0.0811864,-5.58987e-05,2.05818e-08,-3.16625e-12,57313.2,38.3758], Tmin=(100,'K'), Tmax=(1485.58,'K')), NASAPolynomial(coeffs=[14.4987,0.042223,-1.65573e-05,2.92713e-09,-1.95248e-13,53013.6,-37.1656], Tmin=(1485.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(475.3,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C(C)C([CH2])[CH]C(4526)',
    structure = SMILES('[CH2]C=C(C)C([CH2])[CH]C'),
    E0 = (408.173,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0492966,0.0821282,-5.65246e-05,2.07235e-08,-3.16487e-12,49243.1,36.0421], Tmin=(100,'K'), Tmax=(1499.75,'K')), NASAPolynomial(coeffs=[14.9878,0.0420222,-1.64116e-05,2.89239e-09,-1.92485e-13,44732.7,-42.5974], Tmin=(1499.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.173,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])C([CH2])=CC(9438)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])=CC'),
    E0 = (418.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.412823,0.0864882,-6.37624e-05,2.5169e-08,-4.07471e-12,50546.8,36.5824], Tmin=(100,'K'), Tmax=(1450.37,'K')), NASAPolynomial(coeffs=[16.8987,0.0387441,-1.43844e-05,2.47219e-09,-1.62443e-13,45525.2,-53.3726], Tmin=(1450.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C(C)C([CH2])=CC(4527)',
    structure = SMILES('[CH2][CH]C(C)C([CH2])=CC'),
    E0 = (408.337,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0150269,0.0808769,-5.28384e-05,1.75183e-08,-2.38329e-12,49261.9,35.9319], Tmin=(100,'K'), Tmax=(1669.63,'K')), NASAPolynomial(coeffs=[17.0144,0.0400788,-1.61855e-05,2.88323e-09,-1.91926e-13,43575.3,-54.9547], Tmin=(1669.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_P) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C=C)C([CH2])CC(5091)',
    structure = SMILES('[CH2]C=C([CH2])C([CH2])CC'),
    E0 = (365.13,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.266113,0.0812565,-3.85889e-05,-7.16495e-09,8.66624e-12,44079.6,34.5395], Tmin=(100,'K'), Tmax=(1004.22,'K')), NASAPolynomial(coeffs=[16.8723,0.0390109,-1.4353e-05,2.54765e-09,-1.75118e-13,39325.5,-54.7481], Tmin=(1004.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.13,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH2]CC([CH2])C([CH2])C=C(5092)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])C=C'),
    E0 = (493.449,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.126171,0.0816538,-5.06979e-05,1.01418e-08,1.91406e-12,59504.7,38.6898], Tmin=(100,'K'), Tmax=(1013.91,'K')), NASAPolynomial(coeffs=[14.4076,0.0408222,-1.47105e-05,2.53583e-09,-1.69633e-13,55709.1,-35.8121], Tmin=(1013.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH]C)C([CH2])C=C(5834)',
    structure = SMILES('[CH2][C]([CH]C)C([CH2])C=C'),
    E0 = (668.168,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,226.888,578.462,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00284817,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284817,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284817,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284817,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284817,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00284817,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.39707,0.0836971,-8.62716e-05,6.20829e-08,-1.98235e-11,80487.5,38.7234], Tmin=(100,'K'), Tmax=(800.779,'K')), NASAPolynomial(coeffs=[6.10047,0.0522905,-2.19768e-05,4.00654e-09,-2.71969e-13,79667.6,13.0589], Tmin=(800.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(668.168,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C=C)C([CH2])[CH]C(5643)',
    structure = SMILES('[CH2]C=C([CH2])C([CH2])[CH]C'),
    E0 = (559.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.163735,0.081003,-5.29533e-05,1.42283e-08,-2.90468e-13,67471.8,36.7328], Tmin=(100,'K'), Tmax=(1122.12,'K')), NASAPolynomial(coeffs=[15.7028,0.0384242,-1.47241e-05,2.61879e-09,-1.77634e-13,63030.8,-45.5639], Tmin=(1122.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])C=C(5773)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])C=C'),
    E0 = (687.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2950,3100,1380,975,1025,1650,258.129,921.441,2070.15],'cm^-1')),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807886,'amu*angstrom^2'), symmetry=1, barrier=(3.10949,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3541.48,'J/mol'), sigma=(6.54527,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.17 K, Pc=28.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.105516,0.0799229,-6.00862e-05,2.53045e-08,-4.46691e-12,82891.1,40.4157], Tmin=(100,'K'), Tmax=(1320.6,'K')), NASAPolynomial(coeffs=[13.1623,0.0403744,-1.51647e-05,2.62703e-09,-1.73834e-13,79442.6,-26.2068], Tmin=(1320.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(687.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(457.296,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C([CH2])C=C(4950)',
    structure = SMILES('[CH2][CH]C([CH2])C=C'),
    E0 = (532.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,825.24,825.296],'cm^-1')),
        HinderedRotor(inertia=(0.112603,'amu*angstrom^2'), symmetry=1, barrier=(2.58896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112591,'amu*angstrom^2'), symmetry=1, barrier=(2.58869,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0053604,'amu*angstrom^2'), symmetry=1, barrier=(2.59117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.112621,'amu*angstrom^2'), symmetry=1, barrier=(2.58939,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36168,0.0511269,-2.91252e-05,7.14084e-09,-2.95895e-13,64183.7,29.9243], Tmin=(100,'K'), Tmax=(1256.84,'K')), NASAPolynomial(coeffs=[10.4271,0.0291586,-1.11215e-05,1.94845e-09,-1.29797e-13,61361.3,-18.046], Tmin=(1256.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])[CH]C(4147)',
    structure = SMILES('[CH2]C=CC([CH2])[CH]C'),
    E0 = (443.044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,930.455],'cm^-1')),
        HinderedRotor(inertia=(0.180361,'amu*angstrom^2'), symmetry=1, barrier=(4.14686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171002,'amu*angstrom^2'), symmetry=1, barrier=(29.0556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00256024,'amu*angstrom^2'), symmetry=1, barrier=(29.0691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00789359,'amu*angstrom^2'), symmetry=1, barrier=(89.6243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.89892,'amu*angstrom^2'), symmetry=1, barrier=(89.6437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.713248,0.0631297,-2.98162e-05,-7.12472e-10,3.50664e-12,53411.9,32.5467], Tmin=(100,'K'), Tmax=(1093.34,'K')), NASAPolynomial(coeffs=[12.2655,0.0362389,-1.40149e-05,2.51304e-09,-1.71526e-13,49967,-28.419], Tmin=(1093.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH]C)C=C(3803)',
    structure = SMILES('[CH2]C([CH]C)C=C'),
    E0 = (327.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1026.84,3543.87],'cm^-1')),
        HinderedRotor(inertia=(0.235199,'amu*angstrom^2'), symmetry=1, barrier=(5.40768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.235218,'amu*angstrom^2'), symmetry=1, barrier=(5.40812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0271982,'amu*angstrom^2'), symmetry=1, barrier=(20.3613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.76263,'amu*angstrom^2'), symmetry=1, barrier=(86.5103,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39619,0.0496388,-1.64699e-05,-7.33995e-09,4.82136e-12,39497.8,28.008], Tmin=(100,'K'), Tmax=(1082.86,'K')), NASAPolynomial(coeffs=[9.7195,0.0328979,-1.26799e-05,2.27041e-09,-1.54839e-13,36874.1,-16.601], Tmin=(1082.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl)"""),
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
    label = '[CH2][CH]C([CH2])[C]([CH2])[CH]C(9439)',
    structure = SMILES('[CH2][CH]C([CH2])[C]([CH2])[CH]C'),
    E0 = (945.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1531.5,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258229,0.0920677,-0.000122437,1.09289e-07,-3.94626e-11,113879,43.05], Tmin=(100,'K'), Tmax=(851.554,'K')), NASAPolynomial(coeffs=[3.08717,0.0579297,-2.55768e-05,4.7066e-09,-3.18129e-13,114153,34.2952], Tmin=(851.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(945.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C([CH2])[CH]C(9205)',
    structure = SMILES('[CH2][CH][C]([CH2])C([CH2])[CH]C'),
    E0 = (945.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,1531.5,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0139348,'amu*angstrom^2'), symmetry=1, barrier=(10.7468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.258229,0.0920677,-0.000122437,1.09289e-07,-3.94626e-11,113879,43.05], Tmin=(100,'K'), Tmax=(851.554,'K')), NASAPolynomial(coeffs=[3.08717,0.0579297,-2.55768e-05,4.7066e-09,-3.18129e-13,114153,34.2952], Tmin=(851.554,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(945.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C([CH2])[CH][CH2](9024)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])[CH][CH2]'),
    E0 = (965.622,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3050,390,425,1340,1360,335,370,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,1228.66,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0148958,'amu*angstrom^2'), symmetry=1, barrier=(12.8536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0148958,'amu*angstrom^2'), symmetry=1, barrier=(12.8536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0148958,'amu*angstrom^2'), symmetry=1, barrier=(12.8536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0148958,'amu*angstrom^2'), symmetry=1, barrier=(12.8536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0148958,'amu*angstrom^2'), symmetry=1, barrier=(12.8536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0148958,'amu*angstrom^2'), symmetry=1, barrier=(12.8536,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0148958,'amu*angstrom^2'), symmetry=1, barrier=(12.8536,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406373,0.0828491,-7.60898e-05,4.49397e-08,-1.17707e-11,116264,42.4903], Tmin=(100,'K'), Tmax=(891.968,'K')), NASAPolynomial(coeffs=[7.68762,0.0501963,-2.11779e-05,3.89756e-09,-2.67348e-13,114965,8.195], Tmin=(891.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(965.622,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(453.139,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C]([CH2])CC(9109)',
    structure = SMILES('[CH2][CH]C([CH2])[C]([CH2])CC'),
    E0 = (751.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,354.644,3087.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0842088,0.0930555,-0.000110058,8.95056e-08,-3.07162e-11,90489.6,41.1204], Tmin=(100,'K'), Tmax=(838.585,'K')), NASAPolynomial(coeffs=[5.11629,0.0571172,-2.44242e-05,4.4552e-09,-3.00929e-13,90065.3,20.2316], Tmin=(838.585,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(751.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2]C[C]([CH2])C([CH2])[CH]C(5826)',
    structure = SMILES('[CH2]C[C]([CH2])C([CH2])[CH]C'),
    E0 = (751.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,354.644,3087.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0842088,0.0930555,-0.000110058,8.95056e-08,-3.07162e-11,90489.6,41.1204], Tmin=(100,'K'), Tmax=(838.585,'K')), NASAPolynomial(coeffs=[5.11629,0.0571172,-2.44242e-05,4.4552e-09,-3.00929e-13,90065.3,20.2316], Tmin=(838.585,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(751.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][C](C)C([CH2])[CH]C(4529)',
    structure = SMILES('[CH2][CH][C](C)C([CH2])[CH]C'),
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
    label = '[CH2][CH]C([CH2])C([CH2])C[CH2](5764)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])C[CH2]'),
    E0 = (771.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,371.229,1501.2,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0469187,'amu*angstrom^2'), symmetry=1, barrier=(1.55592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0469187,'amu*angstrom^2'), symmetry=1, barrier=(1.55592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0469187,'amu*angstrom^2'), symmetry=1, barrier=(1.55592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0469187,'amu*angstrom^2'), symmetry=1, barrier=(1.55592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0469187,'amu*angstrom^2'), symmetry=1, barrier=(1.55592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0469187,'amu*angstrom^2'), symmetry=1, barrier=(1.55592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0469187,'amu*angstrom^2'), symmetry=1, barrier=(1.55592,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.133501,0.0851864,-6.9243e-05,3.33775e-08,-6.94687e-12,92878.6,41.5957], Tmin=(100,'K'), Tmax=(1114.51,'K')), NASAPolynomial(coeffs=[10.4914,0.0480115,-1.921e-05,3.44929e-09,-2.33547e-13,90569.8,-9.49828], Tmin=(1114.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(771.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C([CH2])CC(9110)',
    structure = SMILES('[CH2][CH][C]([CH2])C([CH2])CC'),
    E0 = (751.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,354.644,3087.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0842088,0.0930555,-0.000110058,8.95056e-08,-3.07162e-11,90489.6,41.1204], Tmin=(100,'K'), Tmax=(838.585,'K')), NASAPolynomial(coeffs=[5.11629,0.0571172,-2.44242e-05,4.4552e-09,-3.00929e-13,90065.3,20.2316], Tmin=(838.585,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(751.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2]CC([CH2])[C]([CH2])[CH]C(5827)',
    structure = SMILES('[CH2]CC([CH2])[C]([CH2])[CH]C'),
    E0 = (751.257,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,354.644,3087.96,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0252434,'amu*angstrom^2'), symmetry=1, barrier=(3.78533,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0842088,0.0930555,-0.000110058,8.95056e-08,-3.07162e-11,90489.6,41.1204], Tmin=(100,'K'), Tmax=(838.585,'K')), NASAPolynomial(coeffs=[5.11629,0.0571172,-2.44242e-05,4.4552e-09,-3.00929e-13,90065.3,20.2316], Tmin=(838.585,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(751.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(C)[C]([CH2])[CH]C(4530)',
    structure = SMILES('[CH2][CH]C(C)[C]([CH2])[CH]C'),
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
    label = '[CH2][C]([CH]C)C([CH2])[CH]C(4562)',
    structure = SMILES('[CH2][C]([CH]C)C([CH2])[CH]C'),
    E0 = (740.552,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,336.774,3065.57,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0223532,'amu*angstrom^2'), symmetry=1, barrier=(3.40049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223532,'amu*angstrom^2'), symmetry=1, barrier=(3.40049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223532,'amu*angstrom^2'), symmetry=1, barrier=(3.40049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223532,'amu*angstrom^2'), symmetry=1, barrier=(3.40049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223532,'amu*angstrom^2'), symmetry=1, barrier=(3.40049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223532,'amu*angstrom^2'), symmetry=1, barrier=(3.40049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0223532,'amu*angstrom^2'), symmetry=1, barrier=(3.40049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.242026,0.0910772,-0.000111057,9.58213e-08,-3.45208e-11,89195.1,41.3224], Tmin=(100,'K'), Tmax=(837.6,'K')), NASAPolynomial(coeffs=[3.20145,0.0603861,-2.64411e-05,4.87231e-09,-3.30681e-13,89280.1,31.0366], Tmin=(837.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(740.552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (760.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (760.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (903.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (917.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (917.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1209.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1213.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1218.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1213.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1226.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1225.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1225.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1215.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (1215.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1215.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (767.907,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (765.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (768.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (768.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (768.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (768.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (782.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (783.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (783.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (783.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (823.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (823.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (823.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (823.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (785.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (785.349,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (891.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (779.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (902.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (880.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (876.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (864.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (885.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (797.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (761.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (1105.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (1157.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (1157.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (1177.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (897.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (901.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (918.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (901.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (918.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (904.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (878.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (876.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (878.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (917.764,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C=C(3743)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]C([CH]C)C([CH2])C=C(4561)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]C([CH2])C(C)[CH][CH2](4513)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH][CH]CC([CH2])[CH]C(4576)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHCH3(T)(21)', '[CH2][CH]C([CH2])[CH][CH2](9040)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[CH2][CH]C([CH2])[CH][CH]C(9076)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][CH2](502)', '[CH2][CH]C([CH2])[CH]C(3871)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C([CH2])[CH]C(4150)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH3(17)', '[CH]C([CH2])C([CH2])[CH][CH2](9250)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C([CH2])[C]C(9431)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH2][C]C([CH2])C([CH2])[CH]C(9432)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH]C([CH][CH2])C([CH2])[CH]C(9433)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH]C([CH]C)C([CH2])[CH][CH2](9434)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH][CH]C([CH2])C([CH2])[CH]C(9435)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C([CH2])C1CC1C(9436)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]C([CH]C)C1CC1[CH2](5830)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]C1C([CH2])C(C)C1[CH2](5769)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C1CC(C)C1[CH2](5796)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]C1CC([CH]C)C1[CH2](5831)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C1CCC1[CH]C(5891)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][C](CC)C([CH2])C=C(5090)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]C(C=C)[C](C)[CH]C(5487)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]CC(=C)C([CH2])[CH]C(9437)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]C=C(C)C([CH2])[CH]C(4526)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]CC([CH2])C([CH2])=CC(9438)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C(C)C([CH2])=CC(4527)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][C](C=C)C([CH2])CC(5091)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][C](C=C)C(C)[CH]C(5488)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C(C)C([CH2])C=C(4528)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]CC([CH2])C([CH2])C=C(5092)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][C]([CH]C)C([CH2])C=C(5834)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2][C](C=C)C([CH2])[CH]C(5643)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C([CH2])C=C(5773)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CH2(T)(20)', '[CH2]C([CH][CH]C)C=C(5043)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(41.7,'m^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds;Y_1centerbirad] for rate rule [Cds-CsH_Cds-CsH;CH2_triplet]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction36',
    reactants = ['CHCH3(T)(21)', '[CH2][CH]C([CH2])C=C(4950)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction37',
    reactants = ['CH2(T)(20)', '[CH2]C=CC([CH2])[CH]C(4147)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH][CH2](502)', '[CH2]C([CH]C)C=C(3803)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'm1_allyl(186)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.0036002,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH]C=C(3743)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH][CH][CH2](5531)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[C]([CH2])[CH]C(9439)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['H(3)', '[CH2][CH][C]([CH2])C([CH2])[CH]C(9205)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C([CH2])[CH][CH2](9024)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(3.53274e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C([CH2])[C]([CH2])CC(9109)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(956916,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C([CH2])[C](C)[CH]C(9283)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]C[C]([CH2])C([CH2])[CH]C(5826)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH][C](C)C([CH2])[CH]C(4529)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][CH]C([CH2])C([CH2])C[CH2](5764)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH][C]([CH2])C([CH2])CC(9110)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH][C]([CH2])C(C)[CH]C(9284)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2]CC([CH2])[C]([CH2])[CH]C(5827)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C(C)[C]([CH2])[CH]C(4530)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][C]([CH]C)C([CH2])[CH]C(4562)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2150',
    isomers = [
        '[CH2][CH]C([CH2])C([CH2])[CH]C(4531)',
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
    label = 'PDepNetwork #2150',
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


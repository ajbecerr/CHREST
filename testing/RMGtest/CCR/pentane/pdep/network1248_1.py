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
    label = '[CH2]C([CH]C)CC=CC(3841)',
    structure = SMILES('[CH2]C([CH]C)CC=CC'),
    E0 = (270.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.200913,0.0791174,-4.80623e-05,1.50728e-08,-1.97754e-12,32618.3,36.3128], Tmin=(100,'K'), Tmax=(1681.41,'K')), NASAPolynomial(coeffs=[14.3182,0.0455335,-1.81022e-05,3.194e-09,-2.11364e-13,27870.9,-39.1308], Tmin=(1681.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(C)C[CH][CH]C(3839)',
    structure = SMILES('[CH2][CH]C(C)C[CH][CH]C'),
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
    label = 'C[CH][CH]CC[CH][CH]C(3833)',
    structure = SMILES('C[CH][CH]CC[CH][CH]C'),
    E0 = (538.384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,180,180,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648374,0.0888459,-0.000122779,1.26818e-07,-5.09293e-11,64858.7,40.5153], Tmin=(100,'K'), Tmax=(850.09,'K')), NASAPolynomial(coeffs=[-4.54426,0.0761568,-3.48852e-05,6.51953e-09,-4.43964e-13,67082.9,72.6128], Tmin=(850.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC)"""),
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
    label = 'C[CH][CH]C[CH][CH]C(3962)',
    structure = SMILES('C[CH][CH]C[CH][CH]C'),
    E0 = (562.164,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3527.57,3926.39,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0224796,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33443,0.0735713,-0.000107041,1.15927e-07,-4.73572e-11,67694.3,35.8086], Tmin=(100,'K'), Tmax=(858.862,'K')), NASAPolynomial(coeffs=[-5.5633,0.0679124,-3.11683e-05,5.81031e-09,-3.9425e-13,70272.7,76.1494], Tmin=(858.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(562.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC)"""),
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
    label = '[CH]C([CH2])C[CH][CH]C(4127)',
    structure = SMILES('[CH]C([CH2])C[CH][CH]C'),
    E0 = (820.863,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,202.301,1102.02,1132.92,1564.72,2182.71,2506.67],'cm^-1')),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0623099,'amu*angstrom^2'), symmetry=1, barrier=(2.4358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.845363,0.0752019,-8.54319e-05,6.98602e-08,-2.439e-11,98835,35.4677], Tmin=(100,'K'), Tmax=(839.521,'K')), NASAPolynomial(coeffs=[3.93026,0.0499017,-2.12846e-05,3.87819e-09,-2.61858e-13,98690.6,23.3496], Tmin=(839.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(820.863,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]CC([CH2])[CH]C(4143)',
    structure = SMILES('[CH][CH]CC([CH2])[CH]C'),
    E0 = (824.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,196.596,1122.26,1203.25,1401.78,1797.35,1932.23],'cm^-1')),
        HinderedRotor(inertia=(0.108539,'amu*angstrom^2'), symmetry=1, barrier=(2.89894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108539,'amu*angstrom^2'), symmetry=1, barrier=(2.89894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108539,'amu*angstrom^2'), symmetry=1, barrier=(2.89894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108539,'amu*angstrom^2'), symmetry=1, barrier=(2.89894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108539,'amu*angstrom^2'), symmetry=1, barrier=(2.89894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108539,'amu*angstrom^2'), symmetry=1, barrier=(2.89894,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924895,0.0716708,-6.57048e-05,4.14942e-08,-1.20315e-11,99246.8,35.1853], Tmin=(100,'K'), Tmax=(800.907,'K')), NASAPolynomial(coeffs=[5.76848,0.0474806,-2.04001e-05,3.78345e-09,-2.60442e-13,98470.9,12.8931], Tmin=(800.907,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(824.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
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
    label = '[CH2]C([C]C)C[CH][CH]C(4563)',
    structure = SMILES('[CH2]C([C]C)C[CH][CH]C'),
    E0 = (800.43,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64164,0.0655592,4.00188e-05,-2.44973e-07,2.35004e-10,96339.7,34.142], Tmin=(100,'K'), Tmax=(421.392,'K')), NASAPolynomial(coeffs=[5.23311,0.057598,-2.4656e-05,4.50041e-09,-3.04005e-13,95805,17.1662], Tmin=(421.392,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH]C)C[CH][C]C(4564)',
    structure = SMILES('[CH2]C([CH]C)C[CH][C]C'),
    E0 = (800.526,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.999454,0.0761754,-2.95228e-05,-5.39803e-08,5.97138e-11,96379,36.0384], Tmin=(100,'K'), Tmax=(510.544,'K')), NASAPolynomial(coeffs=[5.41295,0.0582324,-2.56818e-05,4.827e-09,-3.34979e-13,95711.5,15.5896], Tmin=(510.544,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.526,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([CH]C)C[C][CH]C(4565)',
    structure = SMILES('[CH2]C([CH]C)C[C][CH]C'),
    E0 = (800.514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.321937,0.0858854,-7.73235e-05,4.6364e-08,-1.26837e-11,96407.6,38.5004], Tmin=(100,'K'), Tmax=(845.465,'K')), NASAPolynomial(coeffs=[6.64236,0.0559823,-2.42697e-05,4.52951e-09,-3.13258e-13,95338.9,9.0691], Tmin=(845.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH]C)C[CH][CH]C(4566)',
    structure = SMILES('[CH]C([CH]C)C[CH][CH]C'),
    E0 = (789.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,186.243,558.205,622.409,1386.62,2153.79,2761.9,3513.81],'cm^-1')),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.110639,'amu*angstrom^2'), symmetry=1, barrier=(2.56475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.410654,0.0866821,-9.72662e-05,8.36576e-08,-3.16413e-11,95123.6,40.046], Tmin=(100,'K'), Tmax=(788.149,'K')), NASAPolynomial(coeffs=[2.49917,0.0634269,-2.8921e-05,5.47343e-09,-3.78942e-13,95187.4,32.961], Tmin=(788.149,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(789.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH][CH]CC1CC1C(4567)',
    structure = SMILES('C[CH][CH]CC1CC1C'),
    E0 = (293.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.57945,0.0648221,-1.34221e-05,-1.46666e-08,6.78349e-12,35490.8,34.4707], Tmin=(100,'K'), Tmax=(1163.87,'K')), NASAPolynomial(coeffs=[10.3739,0.0506669,-2.03188e-05,3.68415e-09,-2.51453e-13,31889.8,-19.9438], Tmin=(1163.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C1CC([CH]C)C1C(4524)',
    structure = SMILES('[CH2]C1CC([CH]C)C1C'),
    E0 = (291.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.720226,0.0522484,4.4394e-05,-8.8031e-08,3.56393e-11,35165.3,32.0338], Tmin=(100,'K'), Tmax=(982.822,'K')), NASAPolynomial(coeffs=[14.6923,0.0436331,-1.60968e-05,2.95215e-09,-2.10111e-13,30088.6,-46.9864], Tmin=(982.822,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(291.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH]C1CC([CH]C)C1(4568)',
    structure = SMILES('C[CH]C1CC([CH]C)C1'),
    E0 = (289.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.828701,0.0524005,3.39291e-05,-6.86001e-08,2.61717e-11,34978.4,32.0108], Tmin=(100,'K'), Tmax=(1042.73,'K')), NASAPolynomial(coeffs=[12.5222,0.048737,-2.00592e-05,3.80373e-09,-2.7096e-13,30300.2,-35.6314], Tmin=(1042.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Cs_S)"""),
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
    label = 'C[CH][CH]CC(C)=CC(4497)',
    structure = SMILES('C[CH][CH]CC(C)=CC'),
    E0 = (249.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3050,390,425,1340,1360,335,370,180,400,3152.31],'cm^-1')),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0388829,'amu*angstrom^2'), symmetry=1, barrier=(3.57598,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.17031,0.0665187,-3.05677e-05,6.00427e-09,-4.422e-13,30073.5,28.8102], Tmin=(100,'K'), Tmax=(3278.27,'K')), NASAPolynomial(coeffs=[48.8981,0.0095073,-4.48342e-06,7.00143e-10,-3.77353e-14,-565.886,-252.108], Tmin=(3278.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(=CC)CC[CH]C(4569)',
    structure = SMILES('[CH2]C(=CC)CC[CH]C'),
    E0 = (206.783,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301972,0.0804627,-4.9702e-05,1.59208e-08,-2.16066e-12,25003.1,33.4091], Tmin=(100,'K'), Tmax=(1602.13,'K')), NASAPolynomial(coeffs=[12.9247,0.0489479,-2.01961e-05,3.643e-09,-2.4481e-13,20958.5,-33.4378], Tmin=(1602.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_P)"""),
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
    label = 'C[CH]C=CC(C)[CH]C(4498)',
    structure = SMILES('C[CH]C=CC(C)[CH]C'),
    E0 = (204.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215087,0.0708218,-1.96097e-05,-1.43963e-08,7.68642e-12,24793.9,34.0858], Tmin=(100,'K'), Tmax=(1136.31,'K')), NASAPolynomial(coeffs=[13.4909,0.0474171,-1.95089e-05,3.61161e-09,-2.50456e-13,20270.7,-38.2858], Tmin=(1136.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_S)"""),
)

species(
    label = 'C=CC(C)C[CH][CH]C(3837)',
    structure = SMILES('C=CC(C)C[CH][CH]C'),
    E0 = (269.29,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,215.829,823.733,1222.4,1702.25],'cm^-1')),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132849,'amu*angstrom^2'), symmetry=1, barrier=(3.55709,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15014,0.0706793,-3.54891e-05,8.33376e-09,-7.87455e-13,32481.1,34.2808], Tmin=(100,'K'), Tmax=(2185.09,'K')), NASAPolynomial(coeffs=[14.1606,0.0468625,-1.91396e-05,3.34557e-09,-2.16748e-13,26795.3,-38.6569], Tmin=(2185.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C=C)CC[CH]C(3846)',
    structure = SMILES('[CH2]C(C=C)CC[CH]C'),
    E0 = (279.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.153311,0.0801971,-5.08638e-05,1.74619e-08,-2.5549e-12,33807.5,37.1332], Tmin=(100,'K'), Tmax=(1519.59,'K')), NASAPolynomial(coeffs=[12.6323,0.0473488,-1.8439e-05,3.23669e-09,-2.14609e-13,30014.9,-28.2926], Tmin=(1519.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=CC)C[CH][CH]C(4570)',
    structure = SMILES('[CH2]C(=CC)C[CH][CH]C'),
    E0 = (401.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,288.859,893.046,2080.43],'cm^-1')),
        HinderedRotor(inertia=(0.069921,'amu*angstrom^2'), symmetry=1, barrier=(3.35348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.069921,'amu*angstrom^2'), symmetry=1, barrier=(3.35348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.069921,'amu*angstrom^2'), symmetry=1, barrier=(3.35348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.069921,'amu*angstrom^2'), symmetry=1, barrier=(3.35348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.069921,'amu*angstrom^2'), symmetry=1, barrier=(3.35348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.069921,'amu*angstrom^2'), symmetry=1, barrier=(3.35348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64695,0.062967,2.38496e-05,-1.59802e-07,1.43296e-10,48331.1,31.4236], Tmin=(100,'K'), Tmax=(438.914,'K')), NASAPolynomial(coeffs=[3.69058,0.0611458,-2.73513e-05,5.19097e-09,-3.62742e-13,47989.9,21.4033], Tmin=(438.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(401.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH]C)C=C[CH]C(4571)',
    structure = SMILES('[CH2]C([CH]C)C=C[CH]C'),
    E0 = (410.012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.199639,0.0729337,-3.37582e-05,4.63087e-10,2.97047e-12,49458.7,35.9208], Tmin=(100,'K'), Tmax=(1158.43,'K')), NASAPolynomial(coeffs=[13.2457,0.0440979,-1.74118e-05,3.13651e-09,-2.13599e-13,45348.3,-33.6325], Tmin=(1158.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Allyl_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C=C)C[CH][CH]C(4572)',
    structure = SMILES('[CH2]C(C=C)C[CH][CH]C'),
    E0 = (474.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,400,1161.01,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32824,0.067241,-1.17818e-05,-6.33616e-08,5.75461e-11,57141.7,35.6254], Tmin=(100,'K'), Tmax=(513.968,'K')), NASAPolynomial(coeffs=[4.21148,0.0582682,-2.48959e-05,4.62596e-09,-3.19747e-13,56667.5,21.9041], Tmin=(513.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C)C[CH]C=C(4573)',
    structure = SMILES('[CH2]C=CCC([CH2])[CH]C'),
    E0 = (421.541,'kJ/mol'),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.181225,0.0808417,-5.31453e-05,1.81329e-08,-2.53764e-12,50859.4,37.9862], Tmin=(100,'K'), Tmax=(1645.16,'K')), NASAPolynomial(coeffs=[17.2071,0.0385645,-1.45986e-05,2.51277e-09,-1.64e-13,45138.1,-54.5588], Tmin=(1645.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cs_S) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = 'C[CH][CH]CC=CC(3959)',
    structure = SMILES('C[CH][CH]CC=CC'),
    E0 = (288.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,199.254,410.84,3100.51],'cm^-1')),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0343517,'amu*angstrom^2'), symmetry=1, barrier=(3.18799,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21183,0.0544927,-2.28084e-05,3.60194e-09,-1.42453e-13,34783.4,28.0308], Tmin=(100,'K'), Tmax=(2468.19,'K')), NASAPolynomial(coeffs=[28.511,0.0216028,-8.73419e-06,1.39784e-09,-8.1001e-14,18837,-128.613], Tmin=(2468.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC)"""),
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
    label = '[CH2][C]([CH]C)C[CH][CH]C(4574)',
    structure = SMILES('[CH2][C]([CH]C)C[CH][CH]C'),
    E0 = (732.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,4000,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0513981,'amu*angstrom^2'), symmetry=1, barrier=(1.90888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.24434,0.0541381,-1.90351e-05,1.53347e-09,1.31205e-13,87988.9,27.8064], Tmin=(100,'K'), Tmax=(2822.37,'K')), NASAPolynomial(coeffs=[67.9679,-0.0118015,3.05647e-06,-6.25078e-10,5.13848e-14,42311.4,-362.938], Tmin=(2822.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C)[CH][CH][CH]C(4575)',
    structure = SMILES('[CH2]C([CH]C)[CH][CH][CH]C'),
    E0 = (741.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,1904.01,2561.33,3951.84,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0256116,'amu*angstrom^2'), symmetry=1, barrier=(8.39617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.04732,0.0590742,-2.44524e-05,3.81605e-09,-1.72319e-13,89151.2,33.6666], Tmin=(100,'K'), Tmax=(2868.23,'K')), NASAPolynomial(coeffs=[56.0077,-0.000598494,-6.63779e-07,1.10938e-11,9.03774e-15,52935.7,-287.812], Tmin=(2868.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC) + radical(Isobutyl)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.449419,0.0868474,-0.000103352,9.16657e-08,-3.41896e-11,90564.5,42.4706], Tmin=(100,'K'), Tmax=(830.571,'K')), NASAPolynomial(coeffs=[1.69735,0.0629525,-2.78981e-05,5.17626e-09,-3.52932e-13,90974.1,40.3954], Tmin=(830.571,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(752.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = 'C[CH][CH]C[C](C)[CH]C(4500)',
    structure = SMILES('C[CH][CH]C[C](C)[CH]C'),
    E0 = (527.097,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,360,370,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1765.94,2581.8,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0199712,'amu*angstrom^2'), symmetry=1, barrier=(10.9325,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.4931,0.091778,-0.000127515,1.30183e-07,-5.22604e-11,63507.2,39.8762], Tmin=(100,'K'), Tmax=(837.616,'K')), NASAPolynomial(coeffs=[-3.40655,0.0753008,-3.51508e-05,6.64087e-09,-4.55645e-13,65391.8,65.3488], Tmin=(837.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.097,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([CH]C)[CH]C[CH]C(4577)',
    structure = SMILES('[CH2]C([CH]C)[CH]C[CH]C'),
    E0 = (546.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,231.575,268.692,2293.86,3742.17],'cm^-1')),
        HinderedRotor(inertia=(0.0236076,'amu*angstrom^2'), symmetry=1, barrier=(3.104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0236076,'amu*angstrom^2'), symmetry=1, barrier=(3.104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0236076,'amu*angstrom^2'), symmetry=1, barrier=(3.104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0236076,'amu*angstrom^2'), symmetry=1, barrier=(3.104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0236076,'amu*angstrom^2'), symmetry=1, barrier=(3.104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0236076,'amu*angstrom^2'), symmetry=1, barrier=(3.104,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0236076,'amu*angstrom^2'), symmetry=1, barrier=(3.104,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79559,0.0607139,5.02236e-05,-2.32077e-07,2.1034e-10,65836.1,36.1935], Tmin=(100,'K'), Tmax=(413.508,'K')), NASAPolynomial(coeffs=[3.49356,0.0632933,-2.80718e-05,5.29632e-09,-3.68353e-13,65533.2,27.5366], Tmin=(413.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl)"""),
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
    label = '[CH2]C[CH]CC([CH2])[CH]C(4058)',
    structure = SMILES('[CH2]C[CH]CC([CH2])[CH]C'),
    E0 = (557.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,243.995,501.37,2127.74,3749.88],'cm^-1')),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34826,0.0694129,5.23395e-06,-1.2171e-07,1.10855e-10,67142.7,36.7387], Tmin=(100,'K'), Tmax=(465.227,'K')), NASAPolynomial(coeffs=[4.20272,0.0622043,-2.7413e-05,5.16142e-09,-3.58889e-13,66689.5,23.1358], Tmin=(465.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = 'C[CH][CH][CH]C(C)[CH]C(4501)',
    structure = SMILES('C[CH][CH][CH]C(C)[CH]C'),
    E0 = (536.217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2762.5,2775,2787.5,2800,2812.5,2825,2837.5,2850,1350,1380,1410,1440,1470,1500,700,750,800,1000,1050,1100,1350,1375,1400,900,1000,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,1698.15,2965.81,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022146,'amu*angstrom^2'), symmetry=1, barrier=(11.5117,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.658945,0.0823981,-8.94173e-05,8.31501e-08,-3.39568e-11,64603.6,40.7044], Tmin=(100,'K'), Tmax=(790.11,'K')), NASAPolynomial(coeffs=[-0.589045,0.0706006,-3.26283e-05,6.21492e-09,-4.31746e-13,65366.3,50.0095], Tmin=(790.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][C]([CH]C)CC[CH]C(4578)',
    structure = SMILES('[CH2][C]([CH]C)CC[CH]C'),
    E0 = (537.721,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.234601,0.0932028,-0.000116504,1.06462e-07,-4.01621e-11,64798.4,40.0249], Tmin=(100,'K'), Tmax=(833.912,'K')), NASAPolynomial(coeffs=[1.0008,0.0673787,-3.02121e-05,5.62694e-09,-3.84021e-13,65440.8,41.0851], Tmin=(833.912,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.721,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJC) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH]C)[CH][CH]CC(4057)',
    structure = SMILES('[CH2]C([CH]C)[CH][CH]CC'),
    E0 = (546.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,227.918,1297.42,3287.59,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0335329,'amu*angstrom^2'), symmetry=1, barrier=(4.32063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0335329,'amu*angstrom^2'), symmetry=1, barrier=(4.32063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0335329,'amu*angstrom^2'), symmetry=1, barrier=(4.32063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0335329,'amu*angstrom^2'), symmetry=1, barrier=(4.32063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0335329,'amu*angstrom^2'), symmetry=1, barrier=(4.32063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0335329,'amu*angstrom^2'), symmetry=1, barrier=(4.32063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0335329,'amu*angstrom^2'), symmetry=1, barrier=(4.32063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.485204,0.0836771,-8.11142e-05,6.41038e-08,-2.38782e-11,65891.7,40.3808], Tmin=(100,'K'), Tmax=(763.602,'K')), NASAPolynomial(coeffs=[2.62216,0.0648188,-2.9014e-05,5.47325e-09,-3.79453e-13,65788.8,32.1107], Tmin=(763.602,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(RCCJCC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]CCC([CH2])[CH]C(4579)',
    structure = SMILES('[CH2][CH]CCC([CH2])[CH]C'),
    E0 = (557.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,251.019,465.678,2074.47,3742.99],'cm^-1')),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463966,0.0829174,-6.58677e-05,3.56466e-08,-9.25248e-12,67180,39.8737], Tmin=(100,'K'), Tmax=(863.079,'K')), NASAPolynomial(coeffs=[5.41129,0.0599893,-2.60206e-05,4.8685e-09,-3.37539e-13,66326,16.7341], Tmin=(863.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH]C)C[CH]CC(4056)',
    structure = SMILES('[CH2][C]([CH]C)C[CH]CC'),
    E0 = (537.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1275.3,1907.84,3539.47,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0763153,'amu*angstrom^2'), symmetry=1, barrier=(5.22255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763153,'amu*angstrom^2'), symmetry=1, barrier=(5.22255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763153,'amu*angstrom^2'), symmetry=1, barrier=(5.22255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763153,'amu*angstrom^2'), symmetry=1, barrier=(5.22255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763153,'amu*angstrom^2'), symmetry=1, barrier=(5.22255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763153,'amu*angstrom^2'), symmetry=1, barrier=(5.22255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763153,'amu*angstrom^2'), symmetry=1, barrier=(5.22255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271454,0.0936853,-0.000121729,1.14906e-07,-4.40364e-11,64797.3,39.7205], Tmin=(100,'K'), Tmax=(837.493,'K')), NASAPolynomial(coeffs=[-0.0783735,0.0693062,-3.14078e-05,5.86773e-09,-4.00673e-13,65769.5,46.8004], Tmin=(837.493,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(RCCJCC) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])CC[CH]C(4554)',
    structure = SMILES('[CH2][CH]C([CH2])CC[CH]C'),
    E0 = (557.545,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,251.019,465.678,2074.47,3742.99],'cm^-1')),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0312882,'amu*angstrom^2'), symmetry=1, barrier=(3.28442,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.463966,0.0829174,-6.58677e-05,3.56466e-08,-9.25248e-12,67180,39.8737], Tmin=(100,'K'), Tmax=(863.079,'K')), NASAPolynomial(coeffs=[5.41129,0.0599893,-2.60206e-05,4.8685e-09,-3.37539e-13,66326,16.7341], Tmin=(863.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[CH]CC(4001)',
    structure = SMILES('[CH2][CH]C([CH2])C[CH]CC'),
    E0 = (557.557,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,243.995,501.37,2127.74,3749.88],'cm^-1')),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292172,'amu*angstrom^2'), symmetry=1, barrier=(3.14171,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34826,0.0694129,5.23395e-06,-1.2171e-07,1.10855e-10,67142.7,36.7387], Tmin=(100,'K'), Tmax=(465.227,'K')), NASAPolynomial(coeffs=[4.20272,0.0622043,-2.7413e-05,5.16142e-09,-3.58889e-13,66689.5,23.1358], Tmin=(465.227,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.557,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24245,0.0657809,-2.97104e-05,5.55344e-09,-3.73506e-13,65813.6,34.0867], Tmin=(100,'K'), Tmax=(2892.29,'K')), NASAPolynomial(coeffs=[50.7768,0.00698566,-3.53661e-06,5.1585e-10,-2.41145e-14,34255.6,-257.63], Tmin=(2892.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    E0 = (546.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (546.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (689.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (704.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (712.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (996.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (999.966,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1008.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1012.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1016.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1012.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1012.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1012.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (1001.68,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (554.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (555.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (555.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (569.058,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (569.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (610.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (610.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (610.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (571.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (571.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (624.323,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (629.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (688.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (641.107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (671.209,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (667.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (592.5,'kJ/mol'),
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
    E0 = (943.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (953.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (963.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (963.795,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (683.761,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (688.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (683.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (705.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (699.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (669.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (664.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (662.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (693.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (698.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (664.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (631.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (600.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (600.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (599.639,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['m1_allyl(186)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['[CH2]C([CH]C)CC=CC(3841)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_20] for rate rule [Y_12_20b]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]C(C)C[CH][CH]C(3839)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([CH]C)C([CH2])[CH]C(3835)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(2.66e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CHCH3(T)(21)', '[CH2][CH]C[CH][CH]C(555)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', 'C[CH][CH]C[CH][CH]C(3962)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][CH]C(3874)', '[CH2]C([CH2])[CH]C(500)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH3(17)', '[CH]C([CH2])C[CH][CH]C(4127)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['CH3(17)', '[CH][CH]CC([CH2])[CH]C(4143)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2]C([C]C)C[CH][CH]C(4563)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH2]C([CH]C)C[CH][C]C(4564)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH2]C([CH]C)C[C][CH]C(4565)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH]C([CH]C)C[CH][CH]C(4566)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C[CH][CH]CC1CC1C(4567)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['[CH2]C1CC([CH]C)C1C(4524)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C[CH]C1CC([CH]C)C1(4568)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C=C(CC)C[CH][CH]C(4032)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C[CH][CH]CC(C)=CC(4497)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['[CH2]C(=CC)CC[CH]C(4569)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['[CH2]C(C=C[CH]C)CC(4033)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C[CH]C=CC(C)[CH]C(4498)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C=CC(C)C[CH][CH]C(3837)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['[CH2]C(C=C)CC[CH]C(3846)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2]C(=CC)C[CH][CH]C(4570)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH2]C([CH]C)C=C[CH]C(4571)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2]C(C=C)C[CH][CH]C(4572)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2]C([CH]C)C[CH]C=C(4573)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(T)(20)', 'C[CH][CH]CC=CC(3959)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(41.7,'m^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds;Y_1centerbirad] for rate rule [Cds-CsH_Cds-CsH;CH2_triplet]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction30',
    reactants = ['CHCH3(T)(21)', '[CH2][CH]CC=CC(553)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH][CH]C(3856)', 'm1_allyl(186)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.0018001,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][CH]C(3856)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2][C]([CH]C)C[CH][CH]C(4574)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', '[CH2]C([CH]C)[CH][CH][CH]C(4575)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C[CH][CH]C(4547)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH2][CH][CH]CC([CH2])[CH]C(4576)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['[CH2][C](CC)C[CH][CH]C(4035)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(956916,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C[CH][CH]C[C](C)[CH]C(4500)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C([CH]C)[CH]C[CH]C(4577)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]CC([CH2])C[CH][CH]C(4037)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C[CH]CC([CH2])[CH]C(4058)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([CH][CH][CH]C)CC(4036)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(6.3913e+06,'s^-1'), n=1.66106, Ea=(123.033,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C[CH][CH][CH]C(C)[CH]C(4501)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['[CH2][C]([CH]C)CC[CH]C(4578)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C([CH]C)[CH][CH]CC(4057)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.27681e+06,'s^-1'), n=2.16, Ea=(146.492,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][CH]CCC([CH2])[CH]C(4579)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['[CH2][C]([CH]C)C[CH]CC(4056)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(30253,'s^-1'), n=2.05523, Ea=(118.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R4HJ_1;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][CH]C([CH2])CC[CH]C(4554)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][CH]C([CH2])C[CH]CC(4001)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2][CH][CH]CC([CH2])CC(4038)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(25800,'s^-1'), n=1.67, Ea=(42.6768,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2][CH][CH]CC(C)[CH]C(4502)'],
    products = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1248',
    isomers = [
        '[CH2]C([CH]C)C[CH][CH]C(3834)',
    ],
    reactants = [
        ('m1_allyl(186)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1248',
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


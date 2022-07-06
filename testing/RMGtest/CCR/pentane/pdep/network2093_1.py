species(
    label = '[CH2][CH]C([CH2])C([CH2])[O](4207)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])[O]'),
    E0 = (680.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,248.528,985.032,2220.78],'cm^-1')),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883532,'amu*angstrom^2'), symmetry=1, barrier=(2.71801,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.752264,0.073098,-7.18049e-05,4.06589e-08,-9.63206e-12,81917.4,34.4019], Tmin=(100,'K'), Tmax=(1005.11,'K')), NASAPolynomial(coeffs=[10.3935,0.0347274,-1.4539e-05,2.67408e-09,-1.83712e-13,79979.3,-12.1603], Tmin=(1005.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(CJCO) + radical(RCCJ)"""),
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
    label = '[CH2]C([O])C([CH2])C=C(4445)',
    structure = SMILES('[CH2]C([O])C([CH2])C=C'),
    E0 = (402.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,604.318,3113.3],'cm^-1')),
        HinderedRotor(inertia=(0.0720765,'amu*angstrom^2'), symmetry=1, barrier=(18.6837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.52797,'amu*angstrom^2'), symmetry=1, barrier=(81.1149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0721023,'amu*angstrom^2'), symmetry=1, barrier=(18.6838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139872,'amu*angstrom^2'), symmetry=1, barrier=(3.21594,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3818.2,'J/mol'), sigma=(6.62498,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.39 K, Pc=29.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.491081,0.0696401,-5.37465e-05,1.82912e-08,-1.21122e-12,48543.4,31.4972], Tmin=(100,'K'), Tmax=(1021.56,'K')), NASAPolynomial(coeffs=[14.6069,0.0270162,-9.73121e-06,1.68648e-09,-1.13564e-13,44999.4,-40.1346], Tmin=(1021.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(CJCO)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.89,'J/mol'), sigma=(6.79962,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.22 K, Pc=27.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852795,0.0775059,-0.000107974,9.70252e-08,-3.50878e-11,80900.1,35.2748], Tmin=(100,'K'), Tmax=(840.862,'K')), NASAPolynomial(coeffs=[4.01348,0.0455248,-2.06945e-05,3.86014e-09,-2.62948e-13,80967.6,24.1364], Tmin=(840.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.89,'J/mol'), sigma=(6.79962,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.22 K, Pc=27.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.702035,0.0817735,-0.00011892,1.08066e-07,-3.90055e-11,79108.8,34.8881], Tmin=(100,'K'), Tmax=(846.966,'K')), NASAPolynomial(coeffs=[4.17012,0.0459529,-2.10488e-05,3.92751e-09,-2.67084e-13,79218.7,22.8492], Tmin=(846.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCsJOH) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C([CH2])[O](1582)',
    structure = SMILES('[CH2][CH]C([CH2])[O]'),
    E0 = (530.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1127.94],'cm^-1')),
        HinderedRotor(inertia=(0.109106,'amu*angstrom^2'), symmetry=1, barrier=(2.50857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0027686,'amu*angstrom^2'), symmetry=1, barrier=(2.49694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0182808,'amu*angstrom^2'), symmetry=1, barrier=(16.5026,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78639,0.0431654,-3.16516e-05,9.87645e-09,-5.24257e-13,63867.6,24.633], Tmin=(100,'K'), Tmax=(1082.15,'K')), NASAPolynomial(coeffs=[11.0036,0.0170489,-6.47571e-06,1.15861e-09,-7.93421e-14,61407,-22.7143], Tmin=(1082.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(CJCO) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C([CH2])[O](1906)',
    structure = SMILES('[CH2][CH][CH]C([CH2])[O]'),
    E0 = (700.99,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1166.42,1166.45],'cm^-1')),
        HinderedRotor(inertia=(0.135334,'amu*angstrom^2'), symmetry=1, barrier=(3.11159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00322295,'amu*angstrom^2'), symmetry=1, barrier=(3.11184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00322277,'amu*angstrom^2'), symmetry=1, barrier=(3.11181,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135363,'amu*angstrom^2'), symmetry=1, barrier=(3.11225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6183,0.0547808,-5.28654e-05,3.06285e-08,-7.61423e-12,84393.3,30.224], Tmin=(100,'K'), Tmax=(948.637,'K')), NASAPolynomial(coeffs=[7.65214,0.029338,-1.2634e-05,2.35478e-09,-1.62902e-13,83248.6,1.43248], Tmin=(948.637,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.99,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = 'O(4)',
    structure = SMILES('[O]'),
    E0 = (243.005,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(665.16,'J/mol'), sigma=(2.75,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,29226.7,5.11107], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,29226.7,5.11107], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(243.005,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O""", comment="""Thermo library: BurkeH2O2"""),
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
    label = '[CH2][CH]C([CH2])[CH][O](9332)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][O]'),
    E0 = (677.268,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1516.23,1516.83,1517.08],'cm^-1')),
        HinderedRotor(inertia=(0.0760246,'amu*angstrom^2'), symmetry=1, barrier=(1.74796,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0776275,'amu*angstrom^2'), symmetry=1, barrier=(1.78481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0775813,'amu*angstrom^2'), symmetry=1, barrier=(1.78375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0765518,'amu*angstrom^2'), symmetry=1, barrier=(1.76008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38495,0.0665394,-0.000103342,9.74129e-08,-3.55505e-11,81542,30.1923], Tmin=(100,'K'), Tmax=(862.591,'K')), NASAPolynomial(coeffs=[3.15746,0.0376971,-1.73251e-05,3.21662e-09,-2.17229e-13,82003.4,26.3503], Tmin=(862.591,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(677.268,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Cs_S) + radical(Isobutyl) + radical(CCsJOH) + radical(RCCJ)"""),
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
    label = '[CH2][C]C([CH2])C([CH2])[O](9333)',
    structure = SMILES('[CH2][C]C([CH2])C([CH2])[O]'),
    E0 = (933.812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.358488,0.0792115,-8.66439e-05,5.12282e-08,-1.21341e-11,112444,33.1597], Tmin=(100,'K'), Tmax=(1027.65,'K')), NASAPolynomial(coeffs=[13.7496,0.0270863,-1.05574e-05,1.86696e-09,-1.25414e-13,109691,-31.8096], Tmin=(1027.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(933.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(CJCO) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH][CH2])C([CH2])[O](9334)',
    structure = SMILES('[CH]C([CH][CH2])C([CH2])[O]'),
    E0 = (923.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.794447,0.0731154,-7.39856e-05,4.1456e-08,-9.64028e-12,111157,33.475], Tmin=(100,'K'), Tmax=(1023.66,'K')), NASAPolynomial(coeffs=[11.1648,0.0325919,-1.46041e-05,2.78254e-09,-1.95202e-13,109034,-16.7983], Tmin=(1023.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(923.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(CJCO) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([O])C([CH2])[CH][CH2](9016)',
    structure = SMILES('[CH]C([O])C([CH2])[CH][CH2]'),
    E0 = (916.765,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,212.014,830.042,1022.01,1245.06,1610.29,1763.34],'cm^-1')),
        HinderedRotor(inertia=(0.134105,'amu*angstrom^2'), symmetry=1, barrier=(3.31926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134105,'amu*angstrom^2'), symmetry=1, barrier=(3.31926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134105,'amu*angstrom^2'), symmetry=1, barrier=(3.31926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134105,'amu*angstrom^2'), symmetry=1, barrier=(3.31926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134105,'amu*angstrom^2'), symmetry=1, barrier=(3.31926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.849444,0.070993,-7.16579e-05,4.12735e-08,-9.88801e-12,110373,34.4236], Tmin=(100,'K'), Tmax=(997.153,'K')), NASAPolynomial(coeffs=[10.452,0.0324719,-1.37093e-05,2.52955e-09,-1.74068e-13,108458,-11.8754], Tmin=(997.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(916.765,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C([CH2])C([CH2])[O](9335)',
    structure = SMILES('[CH][CH]C([CH2])C([CH2])[O]'),
    E0 = (923.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.623229,0.0758848,-8.25018e-05,5.02265e-08,-1.25044e-11,111144,34.0825], Tmin=(100,'K'), Tmax=(969.528,'K')), NASAPolynomial(coeffs=[11.5287,0.0308916,-1.28903e-05,2.35987e-09,-1.61526e-13,109030,-18.1925], Tmin=(969.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(923.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(CJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])C1CC1[CH2](5598)',
    structure = SMILES('[CH2]C([O])C1CC1[CH2]'),
    E0 = (427.211,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925899,0.0527656,3.91057e-06,-4.80169e-08,2.39938e-11,51505.8,28.0087], Tmin=(100,'K'), Tmax=(942.785,'K')), NASAPolynomial(coeffs=[16.7932,0.0224492,-6.73016e-06,1.13946e-09,-8.07335e-14,46869.4,-56.3286], Tmin=(942.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(427.211,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(Isobutyl) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]C([CH2])C1CO1(9336)',
    structure = SMILES('[CH2][CH]C([CH2])C1CO1'),
    E0 = (424.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.43701,0.0663843,-5.45518e-05,2.60874e-08,-4.96192e-12,51182.7,32.1213], Tmin=(100,'K'), Tmax=(1455.76,'K')), NASAPolynomial(coeffs=[11.9844,0.0274186,-6.94499e-06,8.70952e-10,-4.50192e-14,48587.5,-25.2905], Tmin=(1455.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1OC([CH2])C1[CH2](5564)',
    structure = SMILES('[CH2]C1OC([CH2])C1[CH2]'),
    E0 = (424.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.206287,0.0786908,-7.77675e-05,4.17558e-08,-8.52453e-12,51172.6,27.5269], Tmin=(100,'K'), Tmax=(1409.64,'K')), NASAPolynomial(coeffs=[15.1644,0.0220958,-3.73358e-06,2.10787e-10,1.90161e-15,48128.7,-47.3309], Tmin=(1409.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(424.115,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Isobutyl) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1CC([O])C1[CH2](5528)',
    structure = SMILES('[CH2]C1CC([O])C1[CH2]'),
    E0 = (416.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39379,0.040108,3.80023e-05,-8.14092e-08,3.56403e-11,50177.1,26.799], Tmin=(100,'K'), Tmax=(931.218,'K')), NASAPolynomial(coeffs=[15.3623,0.0239487,-6.58837e-06,1.07106e-09,-7.57166e-14,45674.7,-49.8021], Tmin=(931.218,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1COC1[CH2](5577)',
    structure = SMILES('[CH2][CH]C1COC1[CH2]'),
    E0 = (420.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1266,0.0566226,-2.0742e-05,-1.61386e-08,1.25417e-11,50738.8,27.3004], Tmin=(100,'K'), Tmax=(873.886,'K')), NASAPolynomial(coeffs=[11.5652,0.0292122,-8.65742e-06,1.31599e-09,-8.23891e-14,48136.6,-26.1032], Tmin=(873.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(Cs_S) + radical(CJC(C)OC) + radical(RCCJ)"""),
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
    label = '[CH2]CC(=C)C([CH2])[O](9337)',
    structure = SMILES('[CH2]CC(=C)C([CH2])[O]'),
    E0 = (396.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.388201,0.0699898,-5.3545e-05,1.86792e-08,-1.94725e-12,47785.9,32.2111], Tmin=(100,'K'), Tmax=(1140,'K')), NASAPolynomial(coeffs=[16.1038,0.0259365,-1.01711e-05,1.84702e-09,-1.27196e-13,43482.2,-48.8272], Tmin=(1140,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C(C)C([CH2])[O](4404)',
    structure = SMILES('[CH2]C=C(C)C([CH2])[O]'),
    E0 = (329.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,277.926],'cm^-1')),
        HinderedRotor(inertia=(0.00401278,'amu*angstrom^2'), symmetry=1, barrier=(2.40336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768385,'amu*angstrom^2'), symmetry=1, barrier=(17.6667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.76799,'amu*angstrom^2'), symmetry=1, barrier=(17.6576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0294772,'amu*angstrom^2'), symmetry=1, barrier=(17.6664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.365304,0.0703335,-5.2269e-05,1.66162e-08,-1.11365e-12,39713.3,29.6801], Tmin=(100,'K'), Tmax=(1114.01,'K')), NASAPolynomial(coeffs=[16.1801,0.0263861,-1.03801e-05,1.89265e-09,-1.30891e-13,35393.2,-51.8999], Tmin=(1114.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(=C)O(4221)',
    structure = SMILES('[CH2][CH]C([CH2])C(=C)O'),
    E0 = (318.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,886.59],'cm^-1')),
        HinderedRotor(inertia=(0.0755603,'amu*angstrom^2'), symmetry=1, barrier=(1.73728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0755616,'amu*angstrom^2'), symmetry=1, barrier=(1.73731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00311458,'amu*angstrom^2'), symmetry=1, barrier=(1.73726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.700329,'amu*angstrom^2'), symmetry=1, barrier=(16.102,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288679,'amu*angstrom^2'), symmetry=1, barrier=(16.1019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.496011,0.0684083,-4.73716e-05,7.48712e-09,3.84798e-12,38411.6,33.0798], Tmin=(100,'K'), Tmax=(960.186,'K')), NASAPolynomial(coeffs=[15.9848,0.0235017,-7.86494e-06,1.33537e-09,-9.04256e-14,34532.8,-45.7245], Tmin=(960.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)=O(4230)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)=O'),
    E0 = (304.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.437306,'amu*angstrom^2'), symmetry=1, barrier=(10.0545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182105,'amu*angstrom^2'), symmetry=1, barrier=(10.0549,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00511079,'amu*angstrom^2'), symmetry=1, barrier=(33.0519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0499496,'amu*angstrom^2'), symmetry=1, barrier=(33.0523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43754,'amu*angstrom^2'), symmetry=1, barrier=(33.0519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02005,0.0696661,-6.71767e-05,3.87794e-08,-9.69296e-12,36683.9,30.6792], Tmin=(100,'K'), Tmax=(937.74,'K')), NASAPolynomial(coeffs=[8.33828,0.0384505,-1.72459e-05,3.28308e-09,-2.29969e-13,35311.4,-4.15681], Tmin=(937.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(CJC(C)C=O) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C=C)C([CH2])O(5601)',
    structure = SMILES('[CH2]C=C([CH2])C([CH2])O'),
    E0 = (250.178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.27762,0.0710047,-4.70311e-05,5.41096e-09,4.2941e-12,30233.1,30.4634], Tmin=(100,'K'), Tmax=(1006.48,'K')), NASAPolynomial(coeffs=[17.7482,0.0234214,-8.67838e-06,1.57599e-09,-1.10866e-13,25609.7,-59.4323], Tmin=(1006.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CJCO) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C=C)C(C)[O](5602)',
    structure = SMILES('[CH2]C=C([CH2])C(C)[O]'),
    E0 = (268.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.903922,'amu*angstrom^2'), symmetry=1, barrier=(20.7829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0277022,'amu*angstrom^2'), symmetry=1, barrier=(20.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0277042,'amu*angstrom^2'), symmetry=1, barrier=(20.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0277077,'amu*angstrom^2'), symmetry=1, barrier=(20.7832,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.600357,0.0617236,-2.09926e-05,-1.93417e-08,1.23232e-11,32481.1,28.4799], Tmin=(100,'K'), Tmax=(1000.58,'K')), NASAPolynomial(coeffs=[16.8745,0.0254026,-9.62415e-06,1.788e-09,-1.27995e-13,27785.8,-57.2317], Tmin=(1000.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])C(=C)[O](2504)',
    structure = SMILES('[CH2]CC([CH2])C(=C)[O]'),
    E0 = (261.521,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,405.999,406.002,406.003],'cm^-1')),
        HinderedRotor(inertia=(0.006215,'amu*angstrom^2'), symmetry=1, barrier=(0.726941,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00621715,'amu*angstrom^2'), symmetry=1, barrier=(0.727194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00621712,'amu*angstrom^2'), symmetry=1, barrier=(0.727284,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0582175,'amu*angstrom^2'), symmetry=1, barrier=(6.80968,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.667782,0.0656988,-4.36726e-05,6.65012e-09,3.58983e-12,31580.5,31.408], Tmin=(100,'K'), Tmax=(953.375,'K')), NASAPolynomial(coeffs=[14.2976,0.0260055,-8.74214e-06,1.46918e-09,-9.8094e-14,28186.7,-37.8661], Tmin=(953.375,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C(=C)[O](2743)',
    structure = SMILES('[CH2][CH]C(C)C(=C)[O]'),
    E0 = (250.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,997.736,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0268206,'amu*angstrom^2'), symmetry=1, barrier=(18.9471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00739858,'amu*angstrom^2'), symmetry=1, barrier=(5.22675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.227331,'amu*angstrom^2'), symmetry=1, barrier=(5.22678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.824074,'amu*angstrom^2'), symmetry=1, barrier=(18.9471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3816.58,'J/mol'), sigma=(6.61681,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.14 K, Pc=29.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.645949,0.0649023,-4.90165e-05,1.93488e-08,-3.09336e-12,30314.1,32.2727], Tmin=(100,'K'), Tmax=(1479.08,'K')), NASAPolynomial(coeffs=[14.9066,0.0263364,-9.90575e-06,1.7206e-09,-1.13823e-13,26095.5,-42.1091], Tmin=(1479.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C=C)C([CH2])[O](5398)',
    structure = SMILES('[CH2]C=C([CH2])C([CH2])[O]'),
    E0 = (480.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,298.841],'cm^-1')),
        HinderedRotor(inertia=(0.951101,'amu*angstrom^2'), symmetry=1, barrier=(21.8677,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0250222,'amu*angstrom^2'), symmetry=1, barrier=(21.8578,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0250061,'amu*angstrom^2'), symmetry=1, barrier=(21.8639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.024977,'amu*angstrom^2'), symmetry=1, barrier=(21.8587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.443203,0.0669768,-4.10006e-05,2.09499e-10,5.9719e-12,57933.5,29.6774], Tmin=(100,'K'), Tmax=(1003.93,'K')), NASAPolynomial(coeffs=[17.6682,0.0216545,-8.10777e-06,1.49215e-09,-1.06187e-13,53300.4,-59.3405], Tmin=(1003.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(480.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(=C)[O](4216)',
    structure = SMILES('[CH2][CH]C([CH2])C(=C)[O]'),
    E0 = (456.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,350.492,1313.27,1313.36],'cm^-1')),
        HinderedRotor(inertia=(0.00133196,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00122576,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0519993,'amu*angstrom^2'), symmetry=1, barrier=(4.79708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0475557,'amu*angstrom^2'), symmetry=1, barrier=(4.86684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3816.58,'J/mol'), sigma=(6.61681,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.14 K, Pc=29.89 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.796865,0.065193,-5.74532e-05,2.77245e-08,-5.42681e-12,54971.4,33.5014], Tmin=(100,'K'), Tmax=(1227.73,'K')), NASAPolynomial(coeffs=[12.8885,0.0257971,-9.31944e-06,1.58696e-09,-1.0435e-13,52002.5,-27.3141], Tmin=(1227.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=CC([CH2])[O](1904)',
    structure = SMILES('[CH2]C([O])[CH]C=C'),
    E0 = (346.094,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,407.463,407.526,407.731],'cm^-1')),
        HinderedRotor(inertia=(0.924763,'amu*angstrom^2'), symmetry=1, barrier=(109.057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223531,'amu*angstrom^2'), symmetry=1, barrier=(26.3511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223695,'amu*angstrom^2'), symmetry=1, barrier=(26.3526,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.919406,0.0595137,-4.50129e-05,1.35263e-08,-3.84499e-13,41743.5,23.3808], Tmin=(100,'K'), Tmax=(1077.08,'K')), NASAPolynomial(coeffs=[14.8409,0.0203864,-8.03235e-06,1.47502e-09,-1.02917e-13,38015.2,-48.2016], Tmin=(1077.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])C=C(1523)',
    structure = SMILES('[CH2]C([O])C=C'),
    E0 = (252.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,529.328,529.355],'cm^-1')),
        HinderedRotor(inertia=(0.0692523,'amu*angstrom^2'), symmetry=1, barrier=(13.7797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0692738,'amu*angstrom^2'), symmetry=1, barrier=(13.7796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95358,0.036167,-6.30157e-06,-1.92918e-08,1.03178e-11,30464.9,21.1322], Tmin=(100,'K'), Tmax=(988.725,'K')), NASAPolynomial(coeffs=[12.0826,0.0154802,-5.70142e-06,1.06015e-09,-7.65757e-14,27470.1,-32.6352], Tmin=(988.725,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([CH][O])C=C(5588)',
    structure = SMILES('[CH2]C([CH][O])C=C'),
    E0 = (402.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,299.272,300.171,1811.51],'cm^-1')),
        HinderedRotor(inertia=(0.116385,'amu*angstrom^2'), symmetry=1, barrier=(7.46329,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0278062,'amu*angstrom^2'), symmetry=1, barrier=(64.9565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.117446,'amu*angstrom^2'), symmetry=1, barrier=(7.46693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49244,0.0585738,-6.87688e-05,5.25411e-08,-1.7037e-11,48554.8,25.9759], Tmin=(100,'K'), Tmax=(815.37,'K')), NASAPolynomial(coeffs=[6.2553,0.0319061,-1.36341e-05,2.49444e-09,-1.69196e-13,47887.9,4.64318], Tmin=(815.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Isobutyl) + radical(CCsJOH)"""),
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
    label = '[CH2][CH][C]([CH2])C([CH2])[O](9169)',
    structure = SMILES('[CH2][CH][C]([CH2])C([CH2])[O]'),
    E0 = (832.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,441.081,1090.08,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0265502,'amu*angstrom^2'), symmetry=1, barrier=(0.781671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0265502,'amu*angstrom^2'), symmetry=1, barrier=(0.781671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0265502,'amu*angstrom^2'), symmetry=1, barrier=(0.781671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0265502,'amu*angstrom^2'), symmetry=1, barrier=(0.781671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0265502,'amu*angstrom^2'), symmetry=1, barrier=(0.781671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01591,0.0665828,-6.10947e-05,3.1963e-08,-7.03345e-12,100259,33.648], Tmin=(100,'K'), Tmax=(1072.66,'K')), NASAPolynomial(coeffs=[10.0957,0.0327232,-1.37446e-05,2.53389e-09,-1.7439e-13,98311,-10.7933], Tmin=(1072.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(832.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(Cs_S) + radical(CJCO) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C]([CH2])[O](9338)',
    structure = SMILES('[CH2][CH]C([CH2])[C]([CH2])[O]'),
    E0 = (856.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,215.622,712.454,2852.05],'cm^-1')),
        HinderedRotor(inertia=(0.023458,'amu*angstrom^2'), symmetry=1, barrier=(0.732871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023458,'amu*angstrom^2'), symmetry=1, barrier=(0.732871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023458,'amu*angstrom^2'), symmetry=1, barrier=(0.732871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023458,'amu*angstrom^2'), symmetry=1, barrier=(0.732871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.023458,'amu*angstrom^2'), symmetry=1, barrier=(0.732871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.434099,0.0838566,-0.000116418,9.35297e-08,-3.02241e-11,103168,34.9624], Tmin=(100,'K'), Tmax=(845.546,'K')), NASAPolynomial(coeffs=[9.25789,0.0343899,-1.49611e-05,2.73251e-09,-1.83966e-13,101952,-4.4941], Tmin=(845.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(856.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(C2CsJOH) + radical(Isobutyl) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2]C[C]([CH2])C([CH2])[O](5595)',
    structure = SMILES('[CH2]C[C]([CH2])C([CH2])[O]'),
    E0 = (638.17,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,206.788,807.033,1654.3],'cm^-1')),
        HinderedRotor(inertia=(0.144771,'amu*angstrom^2'), symmetry=1, barrier=(3.55835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144771,'amu*angstrom^2'), symmetry=1, barrier=(3.55835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144771,'amu*angstrom^2'), symmetry=1, barrier=(3.55835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144771,'amu*angstrom^2'), symmetry=1, barrier=(3.55835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144771,'amu*angstrom^2'), symmetry=1, barrier=(3.55835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.558957,0.0709902,-6.10261e-05,2.85937e-08,-5.47714e-12,76882,32.7273], Tmin=(100,'K'), Tmax=(1244.2,'K')), NASAPolynomial(coeffs=[13.3643,0.0298228,-1.13957e-05,2.00123e-09,-1.33942e-13,73695.4,-31.8491], Tmin=(1244.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(CJCO) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][C](C)C([CH2])[O](4206)',
    structure = SMILES('[CH2][CH][C](C)C([CH2])[O]'),
    E0 = (627.63,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,274.339,991.087,2160.91],'cm^-1')),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0866893,'amu*angstrom^2'), symmetry=1, barrier=(2.69235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02249,0.0646934,-4.81139e-05,1.89014e-08,-3.11709e-12,75594,31.8349], Tmin=(100,'K'), Tmax=(1381.09,'K')), NASAPolynomial(coeffs=[11.6828,0.0338185,-1.45808e-05,2.71467e-09,-1.87042e-13,72649.4,-23.0371], Tmin=(1381.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(Cs_S) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C]([CH2])O(4226)',
    structure = SMILES('[CH2][CH]C([CH2])[C]([CH2])O'),
    E0 = (626.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,360,370,350,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,2769.9],'cm^-1')),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0701369,'amu*angstrom^2'), symmetry=1, barrier=(1.79618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.266716,0.0879048,-0.000122509,9.87791e-08,-3.19003e-11,75468.1,35.755], Tmin=(100,'K'), Tmax=(856.653,'K')), NASAPolynomial(coeffs=[9.36335,0.0361122,-1.55054e-05,2.81002e-09,-1.88113e-13,74251.4,-4.72817], Tmin=(856.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(626.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(Cs_S) + radical(Isobutyl) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])[C](C)[O](4208)',
    structure = SMILES('[CH2][CH]C([CH2])[C](C)[O]'),
    E0 = (645.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,209.961,356.282,3404.73],'cm^-1')),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0382492,'amu*angstrom^2'), symmetry=1, barrier=(1.05112,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00308,0.0727402,-6.98797e-05,2.73547e-08,3.97349e-12,77698.5,32.3443], Tmin=(100,'K'), Tmax=(602.304,'K')), NASAPolynomial(coeffs=[7.94194,0.0390944,-1.7059e-05,3.17111e-09,-2.17961e-13,76637.1,0.513732], Tmin=(602.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])[C]([CH2])[O](5596)',
    structure = SMILES('[CH2]CC([CH2])[C]([CH2])[O]'),
    E0 = (662.224,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.30334,0.0842862,-0.000101839,7.04847e-08,-1.98758e-11,79777.6,32.8806], Tmin=(100,'K'), Tmax=(862.436,'K')), NASAPolynomial(coeffs=[11.2041,0.033726,-1.38974e-05,2.50268e-09,-1.68591e-13,77897.4,-18.0956], Tmin=(862.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(662.224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Isobutyl) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH]C(C)[C]([CH2])[O](4405)',
    structure = SMILES('[CH2][CH]C(C)[C]([CH2])[O]'),
    E0 = (651.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,218.705,936.035,1813.05],'cm^-1')),
        HinderedRotor(inertia=(0.101818,'amu*angstrom^2'), symmetry=1, barrier=(3.20484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101818,'amu*angstrom^2'), symmetry=1, barrier=(3.20484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101818,'amu*angstrom^2'), symmetry=1, barrier=(3.20484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101818,'amu*angstrom^2'), symmetry=1, barrier=(3.20484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101818,'amu*angstrom^2'), symmetry=1, barrier=(3.20484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621707,0.0795784,-9.40755e-05,6.72251e-08,-2.02901e-11,78496.4,32.5191], Tmin=(100,'K'), Tmax=(795.59,'K')), NASAPolynomial(coeffs=[8.65318,0.0391958,-1.79335e-05,3.41761e-09,-2.38446e-13,77218.5,-4.39106], Tmin=(795.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(C2CsJOH) + radical(RCCJ) + radical(CJCO)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C([CH2])O(9339)',
    structure = SMILES('[CH2][CH][C]([CH2])C([CH2])O'),
    E0 = (602.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.854157,0.0705428,-6.67961e-05,3.66398e-08,-8.45763e-12,72558.3,34.4221], Tmin=(100,'K'), Tmax=(1026.81,'K')), NASAPolynomial(coeffs=[10.0542,0.0347038,-1.44415e-05,2.64811e-09,-1.81624e-13,70669,-10.206], Tmin=(1026.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ(C)CO) + radical(Cs_S) + radical(CJCO) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C(C)[O](9340)',
    structure = SMILES('[CH2][CH][C]([CH2])C(C)[O]'),
    E0 = (621.123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08547,0.0624435,-4.52313e-05,1.79208e-08,-3.01521e-12,74810,32.7578], Tmin=(100,'K'), Tmax=(1352.42,'K')), NASAPolynomial(coeffs=[10.5952,0.034317,-1.40356e-05,2.54306e-09,-1.72586e-13,72237.7,-15.9921], Tmin=(1352.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH]C)C([CH2])[O](4292)',
    structure = SMILES('[CH2][C]([CH]C)C([CH2])[O]'),
    E0 = (627.466,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,261.475,954.116,2097.57],'cm^-1')),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0822944,'amu*angstrom^2'), symmetry=1, barrier=(3.2519,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.903515,0.0667947,-5.41555e-05,2.44318e-08,-4.63745e-12,75579.2,32.26], Tmin=(100,'K'), Tmax=(1227.12,'K')), NASAPolynomial(coeffs=[11.031,0.0337824,-1.38022e-05,2.50879e-09,-1.71098e-13,73093.7,-18.6723], Tmin=(1227.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJ(C)CO) + radical(Cs_S) + radical(CJCO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]([O])C([CH2])[CH]C(4446)',
    structure = SMILES('[CH2][C]([O])C([CH2])[CH]C'),
    E0 = (651.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,245.127,884.613,1961.02],'cm^-1')),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0992444,'amu*angstrom^2'), symmetry=1, barrier=(3.42771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.446731,0.0824915,-0.000103553,7.7858e-08,-2.42038e-11,78483.6,33.1335], Tmin=(100,'K'), Tmax=(804.105,'K')), NASAPolynomial(coeffs=[9.30196,0.0369728,-1.59014e-05,2.91671e-09,-1.98086e-13,77107,-7.36197], Tmin=(804.105,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(C2CsJOH) + radical(Isobutyl) + radical(CJCO)"""),
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
    E0 = (680.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (680.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (837.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (837.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1143.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1138.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1329.8,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1115.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1145.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1135.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1128.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1134.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (685.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (685.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (688.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (688.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (688.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (688.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (743.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (743.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (743.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (743.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (700.105,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (689.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (767.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (810.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (841.482,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (775.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (680.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (720.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (1023.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1044.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (1068.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (838.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (821.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (794.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (821.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (796.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (798.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (755.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (798.078,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (837.527,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (811.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['vinoxy(1351)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]C([O])C([CH2])C=C(4445)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH][CH]CC([CH2])[O](4458)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C([CH2])C[CH][O](4373)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH][CH2](502)', '[CH2][CH]C([CH2])[O](1582)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C([CH2])[O](1906)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(4)', '[CH2][CH]C([CH2])[CH][CH2](9040)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.15242e+13,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2(T)(20)', '[CH2][CH]C([CH2])[CH][O](9332)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][C]C([CH2])C([CH2])[O](9333)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C([CH][CH2])C([CH2])[O](9334)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH]C([O])C([CH2])[CH][CH2](9016)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH][CH]C([CH2])C([CH2])[O](9335)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]C([O])C1CC1[CH2](5598)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C([CH2])C1CO1(9336)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]C1OC([CH2])C1[CH2](5564)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]C1CC([O])C1[CH2](5528)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C1COC1[CH2](5577)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C1CCC1[O](5547)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]CC(=C)C([CH2])[O](9337)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]C=C(C)C([CH2])[O](4404)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C([CH2])C(=C)O(4221)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C([CH2])C(C)=O(4230)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][C](C=C)C([CH2])O(5601)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][C](C=C)C(C)[O](5602)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]CC([CH2])C(=C)[O](2504)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C(C)C(=C)[O](2743)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2][C](C=C)C([CH2])[O](5398)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C(=C)[O](4216)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CH2(T)(20)', '[CH2]C=CC([CH2])[O](1904)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH][CH2](502)', '[CH2]C([O])C=C(1523)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction31',
    reactants = ['CH2(T)(20)', '[CH2]C([CH][O])C=C(5588)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['O(4)', '[CH2][CH]C([CH2])C=C(4950)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2803 used for Cds-CsH_Cds-HH;O_atom_triplet
Exact match found for rate rule [Cds-CsH_Cds-HH;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(37.4696,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 36.5 to 37.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction34',
    reactants = ['vinoxy(1351)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4679.9,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][CH][O](1556)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.206e+07,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH2][CH][C]([CH2])C([CH2])[O](9169)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[C]([CH2])[O](9338)'],
    products = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]C[C]([CH2])C([CH2])[O](5595)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH][C](C)C([CH2])[O](4206)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C([CH2])[C]([CH2])O(4226)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C([CH2])[C](C)[O](4208)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2]CC([CH2])[C]([CH2])[O](5596)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH]C(C)[C]([CH2])[O](4405)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH][C]([CH2])C([CH2])O(9339)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(111914,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][CH][C]([CH2])C(C)[O](9340)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][C]([CH]C)C([CH2])[O](4292)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2][CH]C([CH2])C([CH2])[O](4207)'],
    products = ['[CH2][C]([O])C([CH2])[CH]C(4446)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2093',
    isomers = [
        '[CH2][CH]C([CH2])C([CH2])[O](4207)',
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
    label = 'PDepNetwork #2093',
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


species(
    label = '[CH2][CH]C([CH2])C([CH2])CC(3979)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])CC'),
    E0 = (565.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,344.132,543.723,3250.89],'cm^-1')),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00459272,0.0855863,-6.29211e-05,2.6579e-08,-4.83184e-12,68199.7,40.2672], Tmin=(100,'K'), Tmax=(1258.9,'K')), NASAPolynomial(coeffs=[11.5679,0.0488451,-1.91433e-05,3.39583e-09,-2.27965e-13,65288.3,-18.1815], Tmin=(1258.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2]C(C=C)C([CH2])CC(3850)',
    structure = SMILES('[CH2]C(C=C)C([CH2])CC'),
    E0 = (288.203,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3541.48,'J/mol'), sigma=(6.54527,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=553.17 K, Pc=28.66 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.149631,0.0807585,-3.96873e-05,-2.81867e-09,6.64223e-12,34821.4,36.9875], Tmin=(100,'K'), Tmax=(996.011,'K')), NASAPolynomial(coeffs=[14.551,0.0432224,-1.55395e-05,2.69267e-09,-1.81419e-13,30826.5,-39.2296], Tmin=(996.011,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsCsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH2][CH]C([CH2])C([CH2])C(3893)',
    structure = SMILES('[CH2][CH]C([CH2])C([CH2])C'),
    E0 = (586.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,341.472,3095.34],'cm^-1')),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0288517,'amu*angstrom^2'), symmetry=1, barrier=(3.63562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3395.46,'J/mol'), sigma=(6.43099,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=530.36 K, Pc=28.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.596231,0.0714284,-5.11594e-05,2.10748e-08,-3.72485e-12,70636.9,35.8994], Tmin=(100,'K'), Tmax=(1296.67,'K')), NASAPolynomial(coeffs=[10.6242,0.0404935,-1.53733e-05,2.67563e-09,-1.77439e-13,68036.3,-15.0851], Tmin=(1296.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C([CH2])[CH]CCC(9102)',
    structure = SMILES('[CH2][CH]C([CH2])[CH]CCC'),
    E0 = (557.641,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.61968,0.0798284,-5.30534e-05,2.01752e-08,-3.5271e-12,67185.5,39.1161], Tmin=(100,'K'), Tmax=(1210.21,'K')), NASAPolynomial(coeffs=[7.49244,0.0571124,-2.4898e-05,4.66517e-09,-3.23113e-13,65522.1,4.64774], Tmin=(1210.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296787,0.0875608,-8.99057e-05,7.03435e-08,-2.47182e-11,67174.6,40.4657], Tmin=(100,'K'), Tmax=(806.465,'K')), NASAPolynomial(coeffs=[3.66123,0.0622571,-2.68156e-05,4.94189e-09,-3.37175e-13,66912.1,26.695], Tmin=(806.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.461,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C([CH2])[CH]CC(9058)',
    structure = SMILES('[CH2][CH]C([CH2])[CH]CC'),
    E0 = (581.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,474.041,3168.64,4000],'cm^-1')),
        HinderedRotor(inertia=(0.034071,'amu*angstrom^2'), symmetry=1, barrier=(4.78934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.034071,'amu*angstrom^2'), symmetry=1, barrier=(4.78934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.034071,'amu*angstrom^2'), symmetry=1, barrier=(4.78934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.034071,'amu*angstrom^2'), symmetry=1, barrier=(4.78934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.034071,'amu*angstrom^2'), symmetry=1, barrier=(4.78934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.034071,'amu*angstrom^2'), symmetry=1, barrier=(4.78934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25261,0.0652817,-4.03042e-05,1.37754e-08,-2.15489e-12,70023.2,34.5925], Tmin=(100,'K'), Tmax=(1310.8,'K')), NASAPolynomial(coeffs=[6.8267,0.0482719,-2.08391e-05,3.87549e-09,-2.66745e-13,68561.9,6.19219], Tmin=(1310.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH]C([CH2])CC(231)',
    structure = SMILES('[CH2][CH]C([CH2])CC'),
    E0 = (410.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1043.58],'cm^-1')),
        HinderedRotor(inertia=(0.00458383,'amu*angstrom^2'), symmetry=1, barrier=(3.54033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153753,'amu*angstrom^2'), symmetry=1, barrier=(3.53509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153943,'amu*angstrom^2'), symmetry=1, barrier=(3.53945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153951,'amu*angstrom^2'), symmetry=1, barrier=(3.53964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00457906,'amu*angstrom^2'), symmetry=1, barrier=(3.5376,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22899,0.0571953,-3.33927e-05,1.01284e-08,-1.29003e-12,49493.3,29.8858], Tmin=(100,'K'), Tmax=(1719.56,'K')), NASAPolynomial(coeffs=[11.2019,0.0339964,-1.31557e-05,2.28254e-09,-1.49341e-13,46063.5,-23.6337], Tmin=(1719.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C([CH2])CC(587)',
    structure = SMILES('[CH2][CH][CH]C([CH2])CC'),
    E0 = (581.325,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,353.512,3195.31,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0261167,'amu*angstrom^2'), symmetry=1, barrier=(3.92271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5058,0.0619309,-2.04238e-05,-3.53119e-08,3.53303e-11,70000.3,33.9722], Tmin=(100,'K'), Tmax=(536.631,'K')), NASAPolynomial(coeffs=[4.26044,0.0520363,-2.25022e-05,4.21175e-09,-2.92542e-13,69551.5,20.9701], Tmin=(536.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.325,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][C]C([CH2])C([CH2])CC(9103)',
    structure = SMILES('[CH2][C]C([CH2])C([CH2])CC'),
    E0 = (819.507,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.419311,0.0919809,-7.8459e-05,3.77642e-08,-7.51234e-12,98727.4,39.1388], Tmin=(100,'K'), Tmax=(1196,'K')), NASAPolynomial(coeffs=[14.5902,0.0417815,-1.54994e-05,2.66937e-09,-1.76416e-13,95137.1,-35.9598], Tmin=(1196,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(819.507,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH][CH2])C([CH2])CC(9104)',
    structure = SMILES('[CH]C([CH][CH2])C([CH2])CC'),
    E0 = (808.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,208.514,805.664,1007.08,1204.53,1439.65,1630.82],'cm^-1')),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0459656,0.0856377,-6.53086e-05,2.77219e-08,-5.005e-12,97439.5,39.3412], Tmin=(100,'K'), Tmax=(1271.69,'K')), NASAPolynomial(coeffs=[12.5189,0.0464055,-1.90334e-05,3.463e-09,-2.36032e-13,94267.1,-23.8314], Tmin=(1271.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(808.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(CC)C([CH2])[CH][CH2](9105)',
    structure = SMILES('[CH]C(CC)C([CH2])[CH][CH2]'),
    E0 = (808.967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,208.514,805.664,1007.08,1204.53,1439.65,1630.82],'cm^-1')),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146644,'amu*angstrom^2'), symmetry=1, barrier=(3.55608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0459656,0.0856377,-6.53086e-05,2.77219e-08,-5.005e-12,97439.5,39.3412], Tmin=(100,'K'), Tmax=(1271.69,'K')), NASAPolynomial(coeffs=[12.5189,0.0464055,-1.90334e-05,3.463e-09,-2.36032e-13,94267.1,-23.8314], Tmin=(1271.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(808.967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH]C([CH2])C([CH2])CC(9106)',
    structure = SMILES('[CH][CH]C([CH2])C([CH2])CC'),
    E0 = (808.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0858517,0.0878966,-7.18986e-05,3.39159e-08,-6.77683e-12,97425.2,39.8114], Tmin=(100,'K'), Tmax=(1170.76,'K')), NASAPolynomial(coeffs=[12.2437,0.0457713,-1.79265e-05,3.18236e-09,-2.14057e-13,94538.2,-21.6155], Tmin=(1170.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(808.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C(CC)C1CC1[CH2](5086)',
    structure = SMILES('[CH2]C(CC)C1CC1[CH2]'),
    E0 = (309.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.275565,0.0639942,1.76036e-05,-6.86865e-08,3.16772e-11,37381.7,33.5338], Tmin=(100,'K'), Tmax=(946.11,'K')), NASAPolynomial(coeffs=[16.7987,0.0385533,-1.24806e-05,2.13218e-09,-1.47481e-13,32267.3,-55.7709], Tmin=(946.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(511.34,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC(CC)C1[CH2](5087)',
    structure = SMILES('[CH2]C1CC(CC)C1[CH2]'),
    E0 = (301.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.563402,0.0549708,4.38976e-05,-9.45158e-08,4.02833e-11,36439.7,31.8216], Tmin=(100,'K'), Tmax=(949.676,'K')), NASAPolynomial(coeffs=[16.3683,0.0396815,-1.29495e-05,2.24956e-09,-1.5802e-13,31125.3,-55.7865], Tmin=(949.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.792,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C1CCC1CC(5173)',
    structure = SMILES('[CH2][CH]C1CCC1CC'),
    E0 = (300.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.680897,0.0543254,3.47571e-05,-7.40599e-08,2.9274e-11,36272.5,33.1554], Tmin=(100,'K'), Tmax=(1018.47,'K')), NASAPolynomial(coeffs=[14.0711,0.0460705,-1.83816e-05,3.46541e-09,-2.47658e-13,31245.6,-42.978], Tmin=(1018.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.458,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC(=C)C([CH2])CC(9107)',
    structure = SMILES('[CH2]CC(=C)C([CH2])CC'),
    E0 = (280.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.331858,0.0844037,-5.16351e-05,1.19962e-08,4.3724e-13,33932.5,37.1112], Tmin=(100,'K'), Tmax=(1117.19,'K')), NASAPolynomial(coeffs=[15.3103,0.0432823,-1.64075e-05,2.89999e-09,-1.9589e-13,29508.6,-44.2437], Tmin=(1117.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C(C)C([CH2])CC(3976)',
    structure = SMILES('[CH2]C=C(C)C([CH2])CC'),
    E0 = (213.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.357276,0.0847734,-5.04329e-05,1.00047e-08,1.2503e-12,25860,34.5895], Tmin=(100,'K'), Tmax=(1100.16,'K')), NASAPolynomial(coeffs=[15.4212,0.0436775,-1.65869e-05,2.93892e-09,-1.99048e-13,21403.5,-47.514], Tmin=(1100.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])C(=C)CC(9108)',
    structure = SMILES('[CH2]CC([CH2])C(=C)CC'),
    E0 = (280.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.331858,0.0844037,-5.16351e-05,1.19962e-08,4.3724e-13,33932.5,37.1112], Tmin=(100,'K'), Tmax=(1117.19,'K')), NASAPolynomial(coeffs=[15.3103,0.0432823,-1.64075e-05,2.89999e-09,-1.9589e-13,29508.6,-44.2437], Tmin=(1117.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(280.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C(=C)CC(3975)',
    structure = SMILES('[CH2][CH]C(C)C(=C)CC'),
    E0 = (270.218,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.272231,0.0823665,-5.14525e-05,1.60621e-08,-2.0321e-12,32663.2,37.704], Tmin=(100,'K'), Tmax=(1808.09,'K')), NASAPolynomial(coeffs=[19.6501,0.0382917,-1.48871e-05,2.57961e-09,-1.6788e-13,25459,-70.2088], Tmin=(1808.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.218,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](C=C)C(C)CC(5096)',
    structure = SMILES('[CH2]C=C([CH2])C(C)CC'),
    E0 = (160.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.299675,0.0797247,-2.6475e-05,-1.93848e-08,1.2258e-11,19416.9,32.8801], Tmin=(100,'K'), Tmax=(1029.59,'K')), NASAPolynomial(coeffs=[17.3182,0.0419959,-1.62601e-05,2.97836e-09,-2.08323e-13,14161,-60.5359], Tmin=(1029.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.048,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CC([CH2])CC(584)',
    structure = SMILES('[CH2]C=CC([CH2])CC'),
    E0 = (248.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.00895585,'amu*angstrom^2'), symmetry=1, barrier=(5.37435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0247296,'amu*angstrom^2'), symmetry=1, barrier=(14.9822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0245819,'amu*angstrom^2'), symmetry=1, barrier=(15.0149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.231426,'amu*angstrom^2'), symmetry=1, barrier=(5.32094,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.40817,'amu*angstrom^2'), symmetry=1, barrier=(78.3606,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.576005,0.0637746,-1.672e-05,-2.06206e-08,1.19038e-11,30021.3,30.4796], Tmin=(100,'K'), Tmax=(999.081,'K')), NASAPolynomial(coeffs=[13.7051,0.036386,-1.33982e-05,2.38525e-09,-1.64393e-13,26141.4,-39.1368], Tmin=(999.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Isobutyl) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C)CC(228)',
    structure = SMILES('[CH2]C(C=C)CC'),
    E0 = (133.028,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,3600.44],'cm^-1')),
        HinderedRotor(inertia=(0.0052352,'amu*angstrom^2'), symmetry=1, barrier=(1.79571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.737787,'amu*angstrom^2'), symmetry=1, barrier=(16.9632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.733395,'amu*angstrom^2'), symmetry=1, barrier=(16.8622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.66606,'amu*angstrom^2'), symmetry=1, barrier=(84.29,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24814,0.0504002,-3.72534e-06,-2.68812e-08,1.31039e-11,16107.6,25.9804], Tmin=(100,'K'), Tmax=(989.688,'K')), NASAPolynomial(coeffs=[11.2704,0.0328644,-1.19626e-05,2.11942e-09,-1.45819e-13,12998.9,-27.9505], Tmin=(989.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
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
    label = '[CH2][CH]C([CH2])[C](C)CC(9111)',
    structure = SMILES('[CH2][CH]C([CH2])[C](C)CC'),
    E0 = (546.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47047,0.0686612,2.53716e-05,-1.93541e-07,1.81961e-10,65766.9,34.6866], Tmin=(100,'K'), Tmax=(433.871,'K')), NASAPolynomial(coeffs=[4.62695,0.061726,-2.72829e-05,5.11347e-09,-3.5319e-13,65284.4,19.6898], Tmin=(433.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C]([CH2])C([CH2])CC(5081)',
    structure = SMILES('[CH2]C[C]([CH2])C([CH2])CC'),
    E0 = (556.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,374.098,515.552,3083.78],'cm^-1')),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.053924,0.0935766,-9.58248e-05,6.69647e-08,-2.06169e-11,67099.1,39.0647], Tmin=(100,'K'), Tmax=(805.224,'K')), NASAPolynomial(coeffs=[7.05942,0.0564598,-2.33648e-05,4.22649e-09,-2.85654e-13,66011.3,6.64646], Tmin=(805.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH][C](C)C([CH2])CC(3978)',
    structure = SMILES('[CH2][CH][C](C)C([CH2])CC'),
    E0 = (546.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,230.341,2005.9,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47047,0.0686612,2.53716e-05,-1.93541e-07,1.81961e-10,65766.9,34.6866], Tmin=(100,'K'), Tmax=(433.871,'K')), NASAPolynomial(coeffs=[4.62695,0.061726,-2.72829e-05,5.11347e-09,-3.5319e-13,65284.4,19.6898], Tmin=(433.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][C]([CH2])C(C)CC(9112)',
    structure = SMILES('[CH2][CH][C]([CH2])C(C)CC'),
    E0 = (546.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47047,0.0686612,2.53716e-05,-1.93541e-07,1.81961e-10,65766.9,34.6866], Tmin=(100,'K'), Tmax=(433.871,'K')), NASAPolynomial(coeffs=[4.62695,0.061726,-2.72829e-05,5.11347e-09,-3.5319e-13,65284.4,19.6898], Tmin=(433.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)[CH]C(4479)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)[CH]C'),
    E0 = (555.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,230.747,1116.42,3259.86],'cm^-1')),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497458,0.0790439,-4.94967e-05,1.65362e-08,-2.40218e-12,66910,39.2611], Tmin=(100,'K'), Tmax=(1473.06,'K')), NASAPolynomial(coeffs=[10.5374,0.0517809,-2.17351e-05,3.97208e-09,-2.69863e-13,63952.1,-13.0649], Tmin=(1473.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])[C]([CH2])CC(5082)',
    structure = SMILES('[CH2]CC([CH2])[C]([CH2])CC'),
    E0 = (556.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,374.098,515.552,3083.78],'cm^-1')),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355318,'amu*angstrom^2'), symmetry=1, barrier=(2.54399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.053924,0.0935766,-9.58248e-05,6.69647e-08,-2.06169e-11,67099.1,39.0647], Tmin=(100,'K'), Tmax=(805.224,'K')), NASAPolynomial(coeffs=[7.05942,0.0564598,-2.33648e-05,4.22649e-09,-2.85654e-13,66011.3,6.64646], Tmin=(805.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Isobutyl) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C(C)[C]([CH2])CC(3977)',
    structure = SMILES('[CH2][CH]C(C)[C]([CH2])CC'),
    E0 = (546.174,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,230.341,2005.9,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0828678,'amu*angstrom^2'), symmetry=1, barrier=(2.67065,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47047,0.0686612,2.53716e-05,-1.93541e-07,1.81961e-10,65766.9,34.6866], Tmin=(100,'K'), Tmax=(433.871,'K')), NASAPolynomial(coeffs=[4.62695,0.061726,-2.72829e-05,5.11347e-09,-3.5319e-13,65284.4,19.6898], Tmin=(433.871,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.174,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C]([CH]C)C([CH2])CC(4021)',
    structure = SMILES('[CH2][C]([CH]C)C([CH2])CC'),
    E0 = (546.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,498.473,1693.92,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38054,0.0705821,1.87557e-05,-1.83314e-07,1.7374e-10,65750.5,34.9972], Tmin=(100,'K'), Tmax=(445.315,'K')), NASAPolynomial(coeffs=[5.32017,0.0594118,-2.5191e-05,4.5971e-09,-3.11458e-13,65159.5,16.4818], Tmin=(445.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])C([CH2])[CH]C(4022)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])[CH]C'),
    E0 = (565.834,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,344.132,543.723,3250.89],'cm^-1')),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308194,'amu*angstrom^2'), symmetry=1, barrier=(3.48165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.00459272,0.0855863,-6.29211e-05,2.6579e-08,-4.83184e-12,68199.7,40.2672], Tmin=(100,'K'), Tmax=(1258.9,'K')), NASAPolynomial(coeffs=[11.5679,0.0488451,-1.91433e-05,3.39583e-09,-2.27965e-13,65288.3,-18.1815], Tmin=(1258.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C([CH2])[CH]C(3840)',
    structure = SMILES('[CH2][CH]C(C)C([CH2])[CH]C'),
    E0 = (555.294,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,230.747,1116.42,3259.86],'cm^-1')),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0617571,'amu*angstrom^2'), symmetry=1, barrier=(1.97851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3539.78,'J/mol'), sigma=(6.74357,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=552.91 K, Pc=26.19 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497458,0.0790439,-4.94967e-05,1.65362e-08,-2.40218e-12,66910,39.2611], Tmin=(100,'K'), Tmax=(1473.06,'K')), NASAPolynomial(coeffs=[10.5374,0.0517809,-2.17351e-05,3.97208e-09,-2.69863e-13,63952.1,-13.0649], Tmin=(1473.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](CC)C([CH2])[CH]C(4020)',
    structure = SMILES('[CH2][C](CC)C([CH2])[CH]C'),
    E0 = (546.01,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,498.473,1693.92,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0359084,'amu*angstrom^2'), symmetry=1, barrier=(5.01785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38054,0.0705821,1.87557e-05,-1.83314e-07,1.7374e-10,65750.5,34.9972], Tmin=(100,'K'), Tmax=(445.315,'K')), NASAPolynomial(coeffs=[5.32017,0.0594118,-2.5191e-05,4.5971e-09,-3.11458e-13,65159.5,16.4818], Tmin=(445.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.01,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)C[CH2](4516)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)C[CH2]'),
    E0 = (565.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,207.626,1024.57,2155.31],'cm^-1')),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118619,0.0835433,-5.7073e-05,2.12728e-08,-3.39318e-12,68214.6,39.8597], Tmin=(100,'K'), Tmax=(1404,'K')), NASAPolynomial(coeffs=[12.3641,0.0486558,-1.97999e-05,3.57419e-09,-2.41707e-13,64776.1,-23.3727], Tmin=(1404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])C([CH2])C[CH2](5083)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])C[CH2]'),
    E0 = (576.538,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3014.29,3028.57,3042.86,3057.14,3071.43,3085.71,3100,415,431.667,448.333,465,780,803.333,826.667,850,1435,1448.33,1461.67,1475,900,966.667,1033.33,1100,243.83,855.046,1952.15],'cm^-1')),
        HinderedRotor(inertia=(0.0939366,'amu*angstrom^2'), symmetry=1, barrier=(3.21512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939366,'amu*angstrom^2'), symmetry=1, barrier=(3.21512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939366,'amu*angstrom^2'), symmetry=1, barrier=(3.21512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939366,'amu*angstrom^2'), symmetry=1, barrier=(3.21512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939366,'amu*angstrom^2'), symmetry=1, barrier=(3.21512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939366,'amu*angstrom^2'), symmetry=1, barrier=(3.21512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939366,'amu*angstrom^2'), symmetry=1, barrier=(3.21512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.344742,0.0898279,-6.99139e-05,3.08573e-08,-5.70784e-12,69502.6,40.0593], Tmin=(100,'K'), Tmax=(1267.25,'K')), NASAPolynomial(coeffs=[13.8831,0.0449191,-1.67574e-05,2.89335e-09,-1.91244e-13,65896.6,-31.9517], Tmin=(1267.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(576.538,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C([CH2])C[CH2](3980)',
    structure = SMILES('[CH2][CH]C(C)C([CH2])C[CH2]'),
    E0 = (565.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,207.626,1024.57,2155.31],'cm^-1')),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0732707,'amu*angstrom^2'), symmetry=1, barrier=(2.06625,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118619,0.0835433,-5.7073e-05,2.12728e-08,-3.39318e-12,68214.6,39.8597], Tmin=(100,'K'), Tmax=(1404,'K')), NASAPolynomial(coeffs=[12.3641,0.0486558,-1.97999e-05,3.57419e-09,-2.41707e-13,64776.1,-23.3727], Tmin=(1404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsCsH) + longDistanceInteraction_noncyclic(CsCs-TT) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(RCCJ)"""),
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
    E0 = (565.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (565.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1005.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (725.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (725.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (723.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1019.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1024.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1019.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1031.3,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1020.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1020.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (1020.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (571.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (574.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (574.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (588.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (588.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (629.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (629.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (629.234,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (685.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (584.696,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (669.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (690.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (658.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (669.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (566.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (911.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (919.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (963.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (963.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (927.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (972.168,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (982.872,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (668.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (724.033,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (707.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (683.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (712.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (681.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (683.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (723.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (624.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (652.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (696.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (648.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (619.889,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (618.716,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (612.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C=C(3743)', 'butene1(35)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2]C(C=C)C([CH2])CC(3850)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.31e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2 used for Y_12_10
Exact match found for rate rule [Y_12_10]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(S)(11)', '[CH2][CH]C([CH2])C([CH2])C(3893)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.31021e+06,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;C_pri] for rate rule [carbene;C_pri/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C([CH2])C[CH]CC(4001)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C([CH2])[CH]CCC(9102)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH][CH]CC([CH2])CC(4038)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[CH2][CH]C([CH2])[CH]CC(9058)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH][CH2](502)', '[CH2][CH]C([CH2])CC(231)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CH2(T)(20)', '[CH2][CH][CH]C([CH2])CC(587)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][C]C([CH2])C([CH2])CC(9103)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH]C([CH][CH2])C([CH2])CC(9104)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]C(CC)C([CH2])[CH][CH2](9105)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH][CH]C([CH2])C([CH2])CC(9106)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2]C(CC)C1CC1[CH2](5086)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2]C1CC(CC)C1[CH2](5087)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C1CCC1CC(5173)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2]CC(=C)C([CH2])CC(9107)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2]C=C(C)C([CH2])CC(3976)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2]CC([CH2])C(=C)CC(9108)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C(C)C(=C)CC(3975)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][C](C=C)C(C)CC(5096)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][C](CC)C([CH2])C=C(5090)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][C](C=C)C([CH2])CC(5091)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(T)(20)', '[CH2]C=CC([CH2])CC(584)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH][CH2](502)', '[CH2]C(C=C)CC(228)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'butene1(35)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00337229,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C2H5(14)', '[CH2][CH]C([CH2])C=C(4950)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1020,'cm^3/(mol*s)'), n=2.41, Ea=(27.3634,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 418 used for Cds-CsH_Cds-HH;CsJ-CsHH
Exact match found for rate rule [Cds-CsH_Cds-HH;CsJ-CsHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]C=C(3743)', '[CH2][CH]CC(37)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH][CH2](5531)', '[CH2][CH]CC(37)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.5e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C2H5(14)', '[CH2][CH]C([CH2])[CH][CH2](9040)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2][CH]C([CH2])[C]([CH2])CC(9109)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][CH][C]([CH2])C([CH2])CC(9110)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00414786628487, var=0.0, Tref=1000.0, N=2, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R
    Total Standard Deviation in ln(k): 0.0104217745851
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['CH3(17)', '[CH2][CH]C([CH2])C([CH2])[CH2](4093)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C([CH2])[CH]C(4531)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C([CH2])C[CH2](5764)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C([CH2])[C](C)CC(9111)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2]C[C]([CH2])C([CH2])CC(5081)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH][C](C)C([CH2])CC(3978)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH][C]([CH2])C(C)CC(9112)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C([CH2])C(C)[CH]C(4479)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2]CC([CH2])[C]([CH2])CC(5082)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C(C)[C]([CH2])CC(3977)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][C]([CH]C)C([CH2])CC(4021)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(9819.66,'s^-1'), n=2.51904, Ea=(157.388,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;XH_out] for rate rule [R3HJ;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]CC([CH2])C([CH2])[CH]C(4022)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(0.329997,'s^-1'), n=3.35583, Ea=(59.1095,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][CH]C(C)C([CH2])[CH]C(3840)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(6.44e+09,'s^-1'), n=0.13, Ea=(86.6088,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 131 used for R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2][C](CC)C([CH2])[CH]C(4020)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2][CH]C([CH2])C(C)C[CH2](4516)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]CC([CH2])C([CH2])C[CH2](5083)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(833570,'s^-1'), n=1.49392, Ea=(43.3506,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_CCC;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][CH]C(C)C([CH2])C[CH2](3980)'],
    products = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 111 used for R5H_CCC;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]C([CH2])C([CH2])CC(3979)'],
    products = ['[CH2]C([CH]C)C([CH2])[CH]C(3835)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(262000,'s^-1'), n=1.62, Ea=(46.4424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2030',
    isomers = [
        '[CH2][CH]C([CH2])C([CH2])CC(3979)',
    ],
    reactants = [
        ('[CH2][CH]C=C(3743)', 'butene1(35)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2030',
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


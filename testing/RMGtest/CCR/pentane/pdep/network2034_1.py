species(
    label = '[CH2][CH][CH]CC[CH]CC(4073)',
    structure = SMILES('[CH2][CH][CH]CC[CH]CC'),
    E0 = (549.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,271.45,1024.88,2442.93,3402.92,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17742,0.0606168,-2.35457e-05,2.9062e-09,-1.05547e-14,66037,30.9364], Tmin=(100,'K'), Tmax=(2709.5,'K')), NASAPolynomial(coeffs=[52.8622,0.00576982,-3.42484e-06,4.75721e-10,-1.88335e-14,32321.2,-270.821], Tmin=(2709.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]CC[CH]C(3964)',
    structure = SMILES('[CH2][CH][CH]CC[CH]C'),
    E0 = (572.952,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,2091.51,2338.86,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0255693,'amu*angstrom^2'), symmetry=1, barrier=(11.1195,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1702,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3350.7,'J/mol'), sigma=(6.3658,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=523.37 K, Pc=29.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.36111,0.0497551,-1.77213e-05,1.52485e-09,1.23143e-13,68896.7,28.2554], Tmin=(100,'K'), Tmax=(2624.33,'K')), NASAPolynomial(coeffs=[40.8399,0.00972263,-4.60953e-06,6.75852e-10,-3.24062e-14,43339.4,-199.933], Tmin=(2624.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.952,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]CCC([CH2])C(9089)',
    structure = SMILES('[CH2][CH][CH]CCC([CH2])C'),
    E0 = (554.114,'kJ/mol'),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.296787,0.0875608,-8.99057e-05,7.03435e-08,-2.47182e-11,66772,40.4657], Tmin=(100,'K'), Tmax=(806.465,'K')), NASAPolynomial(coeffs=[3.66123,0.0622571,-2.68156e-05,4.94189e-09,-3.37175e-13,66509.5,26.695], Tmin=(806.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(554.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH]CC(26)',
    structure = SMILES('[CH]CC'),
    E0 = (328.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,329.486,329.64,2894.72],'cm^-1')),
        HinderedRotor(inertia=(0.00155193,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00155082,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00369,0.0161456,1.44259e-05,-2.625e-08,1.01938e-11,39516.8,12.6847], Tmin=(100,'K'), Tmax=(1003.28,'K')), NASAPolynomial(coeffs=[6.63834,0.0156476,-5.75059e-06,1.05872e-09,-7.51273e-14,38083.2,-8.3721], Tmin=(1003.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet)"""),
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
    label = '[CH2]C[CH]CC(101)',
    structure = SMILES('[CH2]C[CH]CC'),
    E0 = (231.62,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1730.81,1731.92],'cm^-1')),
        HinderedRotor(inertia=(0.154521,'amu*angstrom^2'), symmetry=1, barrier=(3.55274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154052,'amu*angstrom^2'), symmetry=1, barrier=(3.54195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154831,'amu*angstrom^2'), symmetry=1, barrier=(3.55987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154739,'amu*angstrom^2'), symmetry=1, barrier=(3.55776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19946,0.0413401,-1.64432e-05,2.00796e-09,7.91291e-14,27920.7,22.6765], Tmin=(100,'K'), Tmax=(2000.48,'K')), NASAPolynomial(coeffs=[13.367,0.0260959,-1.03257e-05,1.7398e-09,-1.08622e-13,22034.8,-42.4877], Tmin=(2000.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.62,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJ)"""),
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
    label = '[CH]CC[CH][CH][CH2](9029)',
    structure = SMILES('[CH]CC[CH][CH][CH2]'),
    E0 = (850.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180.336,996.628,1017.97,2460.18,2622.71,3804.05],'cm^-1')),
        HinderedRotor(inertia=(0.129332,'amu*angstrom^2'), symmetry=1, barrier=(2.98256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129332,'amu*angstrom^2'), symmetry=1, barrier=(2.98256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129332,'amu*angstrom^2'), symmetry=1, barrier=(2.98256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129332,'amu*angstrom^2'), symmetry=1, barrier=(2.98256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129332,'amu*angstrom^2'), symmetry=1, barrier=(2.98256,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.1357,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52658,0.0616947,-7.88177e-05,7.3348e-08,-2.79603e-11,102374,31.0603], Tmin=(100,'K'), Tmax=(833.823,'K')), NASAPolynomial(coeffs=[1.6643,0.0451413,-2.04491e-05,3.82194e-09,-2.61205e-13,102903,33.7339], Tmin=(833.823,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(850.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][CH][CH]CC[C]CC(9090)',
    structure = SMILES('[CH2][CH][CH]CC[C]CC'),
    E0 = (802.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.191861,0.0917749,-0.000108535,9.27218e-08,-3.40468e-11,96700.9,39.243], Tmin=(100,'K'), Tmax=(800.413,'K')), NASAPolynomial(coeffs=[3.55742,0.0618882,-2.80368e-05,5.27762e-09,-3.63778e-13,96580.7,26.3703], Tmin=(800.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(802.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH][C]CC[CH]CC(9091)',
    structure = SMILES('[CH2][CH][C]CC[CH]CC'),
    E0 = (802.941,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.191861,0.0917749,-0.000108535,9.27218e-08,-3.40468e-11,96700.9,39.243], Tmin=(100,'K'), Tmax=(800.413,'K')), NASAPolynomial(coeffs=[3.55742,0.0618882,-2.80368e-05,5.27762e-09,-3.63778e-13,96580.7,26.3703], Tmin=(800.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(802.941,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]CC[CH]CC(9092)',
    structure = SMILES('[CH2][C][CH]CC[CH]CC'),
    E0 = (802.953,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.218542,0.0923883,-0.000114269,1.019e-07,-3.8266e-11,96700.2,38.9745], Tmin=(100,'K'), Tmax=(810.159,'K')), NASAPolynomial(coeffs=[2.5091,0.0637602,-2.91991e-05,5.51029e-09,-3.79742e-13,96897.4,31.9139], Tmin=(810.159,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(802.953,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH][CH]CC[CH]CC(9093)',
    structure = SMILES('[CH][CH][CH]CC[CH]CC'),
    E0 = (792.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,220.448,1278.78,1330.5,1349.26,1496.81,1597.13,1777.78,1844.12],'cm^-1')),
        HinderedRotor(inertia=(0.0761017,'amu*angstrom^2'), symmetry=1, barrier=(2.89963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761017,'amu*angstrom^2'), symmetry=1, barrier=(2.89963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761017,'amu*angstrom^2'), symmetry=1, barrier=(2.89963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761017,'amu*angstrom^2'), symmetry=1, barrier=(2.89963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761017,'amu*angstrom^2'), symmetry=1, barrier=(2.89963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761017,'amu*angstrom^2'), symmetry=1, barrier=(2.89963,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0761017,'amu*angstrom^2'), symmetry=1, barrier=(2.89963,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374293,0.0919499,-0.000123656,1.1881e-07,-4.58446e-11,95392.7,40.4616], Tmin=(100,'K'), Tmax=(838.167,'K')), NASAPolynomial(coeffs=[-0.393694,0.0679393,-3.1157e-05,5.84298e-09,-3.99502e-13,96493.6,49.8303], Tmin=(838.167,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(792.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
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
    label = 'C=CCC[CH][CH]CC(5193)',
    structure = SMILES('C=CCC[CH][CH]CC'),
    E0 = (277.274,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,201.808,838.855,1092.57,1379.19,1676.57],'cm^-1')),
        HinderedRotor(inertia=(0.14741,'amu*angstrom^2'), symmetry=1, barrier=(3.47084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14741,'amu*angstrom^2'), symmetry=1, barrier=(3.47084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14741,'amu*angstrom^2'), symmetry=1, barrier=(3.47084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14741,'amu*angstrom^2'), symmetry=1, barrier=(3.47084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14741,'amu*angstrom^2'), symmetry=1, barrier=(3.47084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14741,'amu*angstrom^2'), symmetry=1, barrier=(3.47084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78822,0.0673035,-3.1138e-05,6.16102e-09,-4.5385e-13,33408,32.0182], Tmin=(100,'K'), Tmax=(2878.12,'K')), NASAPolynomial(coeffs=[44.0155,0.013597,-5.74342e-06,8.80074e-10,-4.73636e-14,7038.06,-219.928], Tmin=(2878.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C=C[CH]CCCC(5183)',
    structure = SMILES('[CH2]C=C[CH]CCCC'),
    E0 = (168.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0730568,0.0763564,-2.68608e-05,-1.17196e-08,7.77135e-12,20450.5,32.4182], Tmin=(100,'K'), Tmax=(1102.5,'K')), NASAPolynomial(coeffs=[15.136,0.0457043,-1.85291e-05,3.42159e-09,-2.37807e-13,15606.2,-49.2016], Tmin=(1102.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(168.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=CCCC[CH][CH]C(5195)',
    structure = SMILES('C=CCCC[CH][CH]C'),
    E0 = (277.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58495,0.0688434,-3.31101e-05,7.09892e-09,-5.90393e-13,33415.8,32.9127], Tmin=(100,'K'), Tmax=(2581.35,'K')), NASAPolynomial(coeffs=[19.6917,0.0407855,-1.68059e-05,2.88812e-09,-1.82582e-13,24067.9,-71.6126], Tmin=(2581.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(277.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C=CC[CH][CH]CC(5176)',
    structure = SMILES('[CH2]C=CC[CH][CH]CC'),
    E0 = (416.528,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,247.663,889.367,1257.57,1981.3],'cm^-1')),
        HinderedRotor(inertia=(0.0985091,'amu*angstrom^2'), symmetry=1, barrier=(3.47307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985091,'amu*angstrom^2'), symmetry=1, barrier=(3.47307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985091,'amu*angstrom^2'), symmetry=1, barrier=(3.47307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985091,'amu*angstrom^2'), symmetry=1, barrier=(3.47307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985091,'amu*angstrom^2'), symmetry=1, barrier=(3.47307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0985091,'amu*angstrom^2'), symmetry=1, barrier=(3.47307,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09098,0.0536173,6.89645e-05,-2.67459e-07,2.42301e-10,50153,31.2533], Tmin=(100,'K'), Tmax=(391.839,'K')), NASAPolynomial(coeffs=[2.83877,0.0624921,-2.82053e-05,5.38426e-09,-3.77761e-13,49967.7,26.729], Tmin=(391.839,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJCC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=C[CH]C[CH]CC(5177)',
    structure = SMILES('[CH2]C=C[CH]C[CH]CC'),
    E0 = (363.182,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.338359,0.0747392,-4.22957e-05,1.15891e-08,-1.28243e-12,43816.9,33.5196], Tmin=(100,'K'), Tmax=(2008.18,'K')), NASAPolynomial(coeffs=[18.386,0.0387908,-1.54442e-05,2.67499e-09,-1.72704e-13,36568.3,-66.133], Tmin=(2008.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(Allyl_S) + radical(Allyl_P)"""),
)

species(
    label = 'C=C[CH]CC[CH][CH]C(4584)',
    structure = SMILES('[CH2]C=CCC[CH][CH]C'),
    E0 = (416.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,203.721,802.482,1211.16,1619.86],'cm^-1')),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3576.97,'J/mol'), sigma=(6.55299,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=558.71 K, Pc=28.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998225,0.0725476,-4.49469e-05,1.68524e-08,-3.19921e-12,50197.8,34.9872], Tmin=(100,'K'), Tmax=(1038.59,'K')), NASAPolynomial(coeffs=[4.31352,0.0597792,-2.65059e-05,5.0152e-09,-3.49866e-13,49509.2,18.8673], Tmin=(1038.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(Allyl_P)"""),
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
    label = '[CH2][CH]CC[CH]C=C(4160)',
    structure = SMILES('[CH2][CH]CCC=C[CH2]'),
    E0 = (451.084,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3433.3,'J/mol'), sigma=(6.23992,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.27 K, Pc=32.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.948154,0.0638066,-3.81353e-05,1.15283e-08,-1.44607e-12,54364.8,32.1312], Tmin=(100,'K'), Tmax=(1753.37,'K')), NASAPolynomial(coeffs=[13.0017,0.0363086,-1.46109e-05,2.58389e-09,-1.70758e-13,50138,-32.789], Tmin=(1753.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.084,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH]C[CH][CH]CC(9094)',
    structure = SMILES('[CH2][CH][CH]C[CH][CH]CC'),
    E0 = (743.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,203.401,263.013,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716103,0.0901321,-0.000138665,1.47698e-07,-5.92618e-11,89540.5,42.5801], Tmin=(100,'K'), Tmax=(857.975,'K')), NASAPolynomial(coeffs=[-5.78083,0.0757047,-3.52624e-05,6.60569e-09,-4.49e-13,92301.2,82.5203], Tmin=(857.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C[CH]CC(9095)',
    structure = SMILES('[CH2][CH][CH][CH]C[CH]CC'),
    E0 = (743.643,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,203.401,263.013,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00930035,'amu*angstrom^2'), symmetry=1, barrier=(4.92216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716103,0.0901321,-0.000138665,1.47698e-07,-5.92618e-11,89540.5,42.5801], Tmin=(100,'K'), Tmax=(857.975,'K')), NASAPolynomial(coeffs=[-5.78083,0.0757047,-3.52624e-05,6.60569e-09,-4.49e-13,92301.2,82.5203], Tmin=(857.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.643,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]CC[CH][CH]C(4586)',
    structure = SMILES('[CH2][CH][CH]CC[CH][CH]C'),
    E0 = (743.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,217.748,2872.67,2957.53,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674542,0.0897095,-0.00013367,1.39584e-07,-5.55415e-11,89541.8,42.901], Tmin=(100,'K'), Tmax=(857.285,'K')), NASAPolynomial(coeffs=[-4.6876,0.0737521,-3.40518e-05,6.36129e-09,-4.32043e-13,91967,76.7267], Tmin=(857.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC[CH]C[CH2](9096)',
    structure = SMILES('[CH2][CH][CH]CC[CH]C[CH2]'),
    E0 = (754.431,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,938.437,1411.63,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0207254,'amu*angstrom^2'), symmetry=1, barrier=(9.3768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62777,0.0572323,-2.22985e-05,2.74783e-09,-1.63911e-14,90698.6,31.0442], Tmin=(100,'K'), Tmax=(2806,'K')), NASAPolynomial(coeffs=[60.9604,-0.00487899,4.17435e-07,-1.57664e-10,2.0493e-14,50800.5,-318.468], Tmin=(2806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(754.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[CH]CCC(9078)',
    structure = SMILES('[CH2][CH][CH]C[CH]CCC'),
    E0 = (549.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17742,0.0606168,-2.35457e-05,2.9062e-09,-1.05547e-14,66037,30.9364], Tmin=(100,'K'), Tmax=(2709.5,'K')), NASAPolynomial(coeffs=[52.8622,0.00576982,-3.42484e-06,4.75721e-10,-1.88335e-14,32321.2,-270.821], Tmin=(2709.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C[CH]CC(5169)',
    structure = SMILES('[CH2][CH]C[CH]C[CH]CC'),
    E0 = (549.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,271.45,1024.88,2442.93,3402.92,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17742,0.0606168,-2.35457e-05,2.9062e-09,-1.05547e-14,66037,30.9364], Tmin=(100,'K'), Tmax=(2709.5,'K')), NASAPolynomial(coeffs=[52.8622,0.00576982,-3.42484e-06,4.75721e-10,-1.88335e-14,32321.2,-270.821], Tmin=(2709.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CCC[CH]C(4590)',
    structure = SMILES('[CH2][CH][CH]CCC[CH]C'),
    E0 = (549.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,308.263,1179.1,2686.12,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96347,0.0622395,-2.56759e-05,3.94527e-09,-1.67124e-13,66045.5,31.8732], Tmin=(100,'K'), Tmax=(2832.04,'K')), NASAPolynomial(coeffs=[56.691,0.00151423,-1.54204e-06,1.54282e-10,6.71353e-16,29534.4,-293.994], Tmin=(2832.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]CCCC(9097)',
    structure = SMILES('[CH2][CH][CH][CH]CCCC'),
    E0 = (549.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,271.45,1024.88,2442.93,3402.92,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17742,0.0606168,-2.35457e-05,2.9062e-09,-1.05547e-14,66037,30.9364], Tmin=(100,'K'), Tmax=(2709.5,'K')), NASAPolynomial(coeffs=[52.8622,0.00576982,-3.42484e-06,4.75721e-10,-1.88335e-14,32321.2,-270.821], Tmin=(2709.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[CH][CH]CC(5170)',
    structure = SMILES('[CH2][CH]CC[CH][CH]CC'),
    E0 = (549.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,271.45,1024.88,2442.93,3402.92,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17742,0.0606168,-2.35457e-05,2.9062e-09,-1.05547e-14,66037,30.9364], Tmin=(100,'K'), Tmax=(2709.5,'K')), NASAPolynomial(coeffs=[52.8622,0.00576982,-3.42484e-06,4.75721e-10,-1.88335e-14,32321.2,-270.821], Tmin=(2709.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CCCC[CH2](9098)',
    structure = SMILES('[CH2][CH][CH]CCCC[CH2]'),
    E0 = (559.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,665.087,1594.06,2122.5,3669.16],'cm^-1')),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301222,0.0892213,-9.86902e-05,8.43945e-08,-3.18083e-11,67474.7,40.7317], Tmin=(100,'K'), Tmax=(792.207,'K')), NASAPolynomial(coeffs=[2.34307,0.0658674,-2.9772e-05,5.61299e-09,-3.87765e-13,67560.5,33.9401], Tmin=(792.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH][CH]C[CH]CC(9099)',
    structure = SMILES('[CH2]C[CH][CH]C[CH]CC'),
    E0 = (549.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,208.378,1482.69,1596.53,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39913,0.0589355,-2.1307e-05,1.79894e-09,1.59451e-13,66027.9,29.9686], Tmin=(100,'K'), Tmax=(2625.2,'K')), NASAPolynomial(coeffs=[50.3792,0.00875924,-4.86866e-06,7.30636e-10,-3.46183e-14,33985,-256.076], Tmin=(2625.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C[CH][CH]CC(5171)',
    structure = SMILES('[CH2]C[CH]C[CH][CH]CC'),
    E0 = (549.196,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,208.378,1482.69,1596.53,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0763557,'amu*angstrom^2'), symmetry=1, barrier=(2.17749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39913,0.0589355,-2.1307e-05,1.79894e-09,1.59451e-13,66027.9,29.9686], Tmin=(100,'K'), Tmax=(2625.2,'K')), NASAPolynomial(coeffs=[50.3792,0.00875924,-4.86866e-06,7.30636e-10,-3.46183e-14,33985,-256.076], Tmin=(2625.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH]C[CH]CC(4071)',
    structure = SMILES('C[CH][CH][CH]C[CH]CC'),
    E0 = (538.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,207.65,2412.37,2756.83,3726.61,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688786,0.0892831,-0.00012783,1.35012e-07,-5.46869e-11,64857.4,40.8916], Tmin=(100,'K'), Tmax=(851.35,'K')), NASAPolynomial(coeffs=[-5.63412,0.0781034,-3.60923e-05,6.76306e-09,-4.60849e-13,67415.8,79.0808], Tmin=(851.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][CH]CCC[CH][CH]C(4589)',
    structure = SMILES('[CH2][CH]CCC[CH][CH]C'),
    E0 = (549.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,308.264,1179.08,2715.63,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96347,0.0622396,-2.56759e-05,3.94528e-09,-1.67126e-13,66045.5,31.8732], Tmin=(100,'K'), Tmax=(2832.05,'K')), NASAPolynomial(coeffs=[56.6929,0.00151205,-1.54117e-06,1.54127e-10,6.8147e-16,29533,-294.007], Tmin=(2832.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C[CH][CH]CC(4070)',
    structure = SMILES('C[CH][CH]C[CH][CH]CC'),
    E0 = (538.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,207.651,2412.36,2756.86,3726.59,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688929,0.0892813,-0.000127823,1.35001e-07,-5.46816e-11,64857.4,40.8911], Tmin=(100,'K'), Tmax=(851.367,'K')), NASAPolynomial(coeffs=[-5.63442,0.0781039,-3.60927e-05,6.76314e-09,-4.60856e-13,67415.9,79.0825], Tmin=(851.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]CC[CH][CH]C(4072)',
    structure = SMILES('[CH2]C[CH]CC[CH][CH]C'),
    E0 = (549.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,271.45,1024.88,2442.93,3402.92,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17742,0.0606168,-2.35457e-05,2.9062e-09,-1.05547e-14,66037,30.9364], Tmin=(100,'K'), Tmax=(2709.5,'K')), NASAPolynomial(coeffs=[52.8622,0.00576982,-3.42484e-06,4.75721e-10,-1.88335e-14,32321.2,-270.821], Tmin=(2709.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CCC[CH]C[CH2](9100)',
    structure = SMILES('[CH2][CH]CCC[CH]C[CH2]'),
    E0 = (559.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,665.087,1594.06,2122.5,3669.16],'cm^-1')),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10849,'amu*angstrom^2'), symmetry=1, barrier=(2.49439,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.301222,0.0892213,-9.86902e-05,8.43945e-08,-3.18083e-11,67474.7,40.7317], Tmin=(100,'K'), Tmax=(792.207,'K')), NASAPolynomial(coeffs=[2.34307,0.0658674,-2.9772e-05,5.61299e-09,-3.87765e-13,67560.5,33.9401], Tmin=(792.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]CC[CH]C[CH2](9101)',
    structure = SMILES('[CH2]C[CH]CC[CH]C[CH2]'),
    E0 = (559.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2764.29,2778.57,2792.86,2807.14,2821.43,2835.71,2850,1425,1433.33,1441.67,1450,1225,1241.67,1258.33,1275,1270,1293.33,1316.67,1340,700,733.333,766.667,800,300,333.333,366.667,400,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,525.171,886.361,2209.87,2598.85,3543.17],'cm^-1')),
        HinderedRotor(inertia=(0.10317,'amu*angstrom^2'), symmetry=1, barrier=(3.56249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10317,'amu*angstrom^2'), symmetry=1, barrier=(3.56249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10317,'amu*angstrom^2'), symmetry=1, barrier=(3.56249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10317,'amu*angstrom^2'), symmetry=1, barrier=(3.56249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10317,'amu*angstrom^2'), symmetry=1, barrier=(3.56249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10317,'amu*angstrom^2'), symmetry=1, barrier=(3.56249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10317,'amu*angstrom^2'), symmetry=1, barrier=(3.56249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.325374,0.0898673,-0.000104552,9.3754e-08,-3.61113e-11,67474.1,39.779], Tmin=(100,'K'), Tmax=(803.022,'K')), NASAPolynomial(coeffs=[1.30312,0.0677242,-3.09252e-05,5.84342e-09,-4.03539e-13,67874,38.744], Tmin=(803.022,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(559.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJ) + radical(RCCJ)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3496.45,'J/mol'), sigma=(6.67923,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=546.14 K, Pc=26.63 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648374,0.0888459,-0.000122779,1.26818e-07,-5.09293e-11,64858.7,40.5153], Tmin=(100,'K'), Tmax=(850.09,'K')), NASAPolynomial(coeffs=[-4.54426,0.0761568,-3.48852e-05,6.51953e-09,-4.43964e-13,67082.9,72.6128], Tmin=(850.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC)"""),
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
    E0 = (549.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (991.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (717.396,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (714.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (748.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1014.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1014.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (1015.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1014.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (1014.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (1014.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (1003.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (557.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (612.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (612.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (574.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (634.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (582.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (634.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (653.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (566.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (610.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (911.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (955.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (955.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (914.749,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (955.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (966.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (708.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (686.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (708.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (662.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (665.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (701.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (723.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (701.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (624.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (592.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (622.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (609.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (633.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (640.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (604.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['[CH2][CH]C=C(3743)', 'butene1(35)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2][CH][CH]CC[CH]C(3964)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(215646,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH][CH]CC([CH2])CC(4038)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2][CH]C([CH2])C[CH]CC(4001)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2][CH][CH]CCC([CH2])C(9089)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-HH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]CC(26)', '[CH2][CH][CH]C[CH2](509)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH][CH][CH2](5537)', '[CH2]C[CH]CC(101)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C2H5(14)', '[CH]CC[CH][CH][CH2](9029)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][CH][CH]CC[C]CC(9090)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH2][CH][C]CC[CH]CC(9091)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][C][CH]CC[CH]CC(9092)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH][CH][CH]CC[CH]CC(9093)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['[CH2][CH]C1CCC1CC(5173)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['C=CCC[CH][CH]CC(5193)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['[CH2]C=C[CH]CCCC(5183)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['C=CCCC[CH][CH]C(5195)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2]C=CC[CH][CH]CC(5176)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2]C=C[CH]C[CH]CC(5177)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', 'C=C[CH]CC[CH][CH]C(4584)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH][CH][CH2](5531)', 'butene1(35)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.0344575,'m^3/(mol*s)'), n=2.32603, Ea=(14.6351,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C=C(3743)', '[CH2][CH]CC(37)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['CH3(17)', '[CH2][CH]CC[CH]C=C(4160)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.105698,'m^3/(mol*s)'), n=2.13, Ea=(23.0748,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH][CH][CH2](5531)', '[CH2][CH]CC(37)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
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
    reactants = ['H(3)', '[CH2][CH][CH]C[CH][CH]CC(9094)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][CH][CH][CH]C[CH]CC(9095)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH3(17)', '[CH2][CH][CH]CC[CH][CH2](4163)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3.25e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH2][CH][CH]CC[CH][CH]C(4586)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2][CH][CH]CC[CH]C[CH2](9096)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['[CH2][CH][CH]C[CH]CCC(9078)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH]C[CH]C[CH]CC(5169)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['[CH2][CH][CH]CCC[CH]C(4590)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][CH][CH]CCCC(9097)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(138391,'s^-1'), n=2.1085, Ea=(113.052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][CH]CC[CH][CH]CC(5170)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][CH][CH]CCCC[CH2](9098)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 108 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C[CH][CH]C[CH]CC(9099)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(867345,'s^-1'), n=1.96939, Ea=(174.054,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R3HJ;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C[CH]C[CH][CH]CC(5171)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(645397,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_2;Y_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['C[CH][CH][CH]C[CH]CC(4071)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.32767e+06,'s^-1'), n=1.53625, Ea=(75.6258,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['[CH2][CH]CCC[CH][CH]C(4589)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(379583,'s^-1'), n=1.54051, Ea=(43.2505,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['C[CH][CH]C[CH][CH]CC(4070)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C[CH]CC[CH][CH]C(4072)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(92.2,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2][CH]CCC[CH]C[CH2](9100)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2778.79,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_2;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C[CH]CC[CH]C[CH2](9101)'],
    products = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(243903,'s^-1'), n=1.58648, Ea=(80.9836,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1062,'s^-1'), n=1.81, Ea=(55.2288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2034',
    isomers = [
        '[CH2][CH][CH]CC[CH]CC(4073)',
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
    label = 'PDepNetwork #2034',
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


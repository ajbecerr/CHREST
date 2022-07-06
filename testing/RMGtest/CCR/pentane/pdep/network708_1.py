species(
    label = 'C[CH]C(C)C[C]=O(2311)',
    structure = SMILES('C[CH]C(C)C[C]=O'),
    E0 = (74.1445,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,1855,455,950,280.823,3290.4],'cm^-1')),
        HinderedRotor(inertia=(0.00209491,'amu*angstrom^2'), symmetry=1, barrier=(0.119643,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141957,'amu*angstrom^2'), symmetry=1, barrier=(8.03689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149774,'amu*angstrom^2'), symmetry=1, barrier=(8.0337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270155,'amu*angstrom^2'), symmetry=1, barrier=(14.8,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00220697,'amu*angstrom^2'), symmetry=1, barrier=(0.121667,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35185,0.0632097,-4.76434e-05,2.18791e-08,-4.69206e-12,9008.7,28.7385], Tmin=(100,'K'), Tmax=(1020.67,'K')), NASAPolynomial(coeffs=[6.1533,0.0443927,-1.99893e-05,3.81618e-09,-2.67756e-13,8028.57,5.47606], Tmin=(1020.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.1445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCCJ=O)"""),
)

species(
    label = 'ketene(1375)',
    structure = SMILES('C=C=O'),
    E0 = (-60.9858,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,2120,512.5,787.5],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.48302,0.0083299,5.93417e-06,-1.33939e-08,5.73133e-12,-7313.53,6.51898], Tmin=(100,'K'), Tmax=(954.863,'K')), NASAPolynomial(coeffs=[5.88246,0.00580173,-1.91266e-06,3.35917e-10,-2.37398e-14,-8114.73,-6.74202], Tmin=(954.863,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.9858,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""ketene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'butene2t(396)',
    structure = SMILES('CC=CC'),
    E0 = (-28.2774,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2995,3025,975,1000,1300,1375,400,500,1630,1680],'cm^-1')),
        HinderedRotor(inertia=(0.464594,'amu*angstrom^2'), symmetry=1, barrier=(10.6819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.460586,'amu*angstrom^2'), symmetry=1, barrier=(10.5898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52,0.0276811,9.44413e-07,-1.25016e-08,4.67173e-12,-3343.75,13.0195], Tmin=(100,'K'), Tmax=(1157.77,'K')), NASAPolynomial(coeffs=[6.04212,0.0260075,-1.04845e-05,1.90884e-09,-1.30603e-13,-4862.7,-7.52639], Tmin=(1157.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-28.2774,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), label="""butene2t""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C[CH]CC[C]=O(2248)',
    structure = SMILES('C[CH]CC[C]=O'),
    E0 = (103.523,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,2026.63,2026.63],'cm^-1')),
        HinderedRotor(inertia=(0.157416,'amu*angstrom^2'), symmetry=1, barrier=(7.79474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157084,'amu*angstrom^2'), symmetry=1, barrier=(7.77066,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157403,'amu*angstrom^2'), symmetry=1, barrier=(7.77178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157528,'amu*angstrom^2'), symmetry=1, barrier=(7.78262,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3428.95,'J/mol'), sigma=(5.98513,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.59 K, Pc=36.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.66841,0.0563955,-6.69739e-05,5.77556e-08,-2.1192e-11,12530,25.804], Tmin=(100,'K'), Tmax=(815.57,'K')), NASAPolynomial(coeffs=[3.41889,0.03858,-1.72312e-05,3.21782e-09,-2.20556e-13,12551.4,19.5978], Tmin=(815.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.523,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCCJ=O)"""),
)

species(
    label = 'CC(C)[CH]C[C]=O(2713)',
    structure = SMILES('CC(C)[CH]C[C]=O'),
    E0 = (79.5044,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1855,455,950,270.606,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00230189,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00230208,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434407,'amu*angstrom^2'), symmetry=1, barrier=(22.5726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127205,'amu*angstrom^2'), symmetry=1, barrier=(6.60998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.250813,'amu*angstrom^2'), symmetry=1, barrier=(13.0333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.979587,0.0638264,-4.43812e-05,1.6477e-08,-2.57459e-12,9672.71,29.3065], Tmin=(100,'K'), Tmax=(1449.62,'K')), NASAPolynomial(coeffs=[11.588,0.0345543,-1.40919e-05,2.54737e-09,-1.72301e-13,6597.06,-25.812], Tmin=(1449.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.5044,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=O)C(C)[CH]C(2027)',
    structure = SMILES('C=C([O])C(C)[CH]C'),
    E0 = (45.7342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,180,779.905,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0334751,'amu*angstrom^2'), symmetry=1, barrier=(14.4511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00848413,'amu*angstrom^2'), symmetry=1, barrier=(3.66081,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.628529,'amu*angstrom^2'), symmetry=1, barrier=(14.4511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.62853,'amu*angstrom^2'), symmetry=1, barrier=(14.4511,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.722002,0.0629666,-3.49958e-05,3.37878e-09,2.54443e-12,5626.3,30.2046], Tmin=(100,'K'), Tmax=(1084.24,'K')), NASAPolynomial(coeffs=[13.4125,0.0312882,-1.2115e-05,2.18829e-09,-1.50469e-13,1984.5,-36.1499], Tmin=(1084.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(45.7342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S)"""),
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
    label = 'C[CH]C[C]=O(2361)',
    structure = SMILES('C[CH]C[C]=O'),
    E0 = (132.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,204.852],'cm^-1')),
        HinderedRotor(inertia=(0.00401758,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636302,'amu*angstrom^2'), symmetry=1, barrier=(18.9512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000722328,'amu*angstrom^2'), symmetry=1, barrier=(6.83812,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35404,0.0356784,-2.29653e-05,7.91178e-09,-1.1666e-12,16026.6,20.8557], Tmin=(100,'K'), Tmax=(1491.52,'K')), NASAPolynomial(coeffs=[7.59578,0.0216208,-8.82771e-06,1.59264e-09,-1.07418e-13,14463,-6.52828], Tmin=(1491.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[C]=O(2355)',
    structure = SMILES('[C]=O'),
    E0 = (438.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3323.79],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09152,0.00193295,-1.59153e-05,2.47563e-08,-1.11287e-11,52770.5,4.46624], Tmin=(100,'K'), Tmax=(866.029,'K')), NASAPolynomial(coeffs=[1.05092,0.00526657,-3.13864e-06,6.40624e-10,-4.47939e-14,53698.8,21.0169], Tmin=(866.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C(C)[CH]C(114)',
    structure = SMILES('[CH2]C(C)[CH]C'),
    E0 = (225.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1290.57],'cm^-1')),
        HinderedRotor(inertia=(0.100358,'amu*angstrom^2'), symmetry=1, barrier=(2.30744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10319,'amu*angstrom^2'), symmetry=1, barrier=(2.37253,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0987185,'amu*angstrom^2'), symmetry=1, barrier=(2.26973,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0971084,'amu*angstrom^2'), symmetry=1, barrier=(2.23271,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.1329,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80047,0.0421184,-1.07295e-05,-7.76584e-09,4.08549e-12,27247.4,23.8047], Tmin=(100,'K'), Tmax=(1119,'K')), NASAPolynomial(coeffs=[7.75174,0.0320963,-1.23774e-05,2.20164e-09,-1.48914e-13,25211,-8.72346], Tmin=(1119,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl)"""),
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
    label = '[CH]C(C)C[C]=O(2389)',
    structure = SMILES('[CH]C(C)C[C]=O'),
    E0 = (348.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,374.091,374.102,374.113,2336.58],'cm^-1')),
        HinderedRotor(inertia=(0.0012049,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00120459,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0918026,'amu*angstrom^2'), symmetry=1, barrier=(9.1175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0918121,'amu*angstrom^2'), symmetry=1, barrier=(9.11716,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57878,0.054869,-4.98784e-05,2.60779e-08,-5.79596e-12,41970.7,24.1195], Tmin=(100,'K'), Tmax=(1055.26,'K')), NASAPolynomial(coeffs=[8.5698,0.0283692,-1.221e-05,2.28056e-09,-1.58134e-13,40495.3,-9.98407], Tmin=(1055.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = 'C[C]C(C)C[C]=O(2714)',
    structure = SMILES('C[C]C(C)C[C]=O'),
    E0 = (327.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.961983,0.069288,-6.24035e-05,3.23947e-08,-7.18678e-12,39534.7,27.4815], Tmin=(100,'K'), Tmax=(1053.12,'K')), NASAPolynomial(coeffs=[9.54499,0.0366875,-1.5969e-05,2.99963e-09,-2.08657e-13,37726.9,-14.3707], Tmin=(1053.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'CC1CC(=O)C1C(2312)',
    structure = SMILES('CC1CC(=O)C1C'),
    E0 = (-178.292,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64567,0.035563,4.44197e-05,-7.6285e-08,2.97507e-11,-21344.5,22.5939], Tmin=(100,'K'), Tmax=(993.312,'K')), NASAPolynomial(coeffs=[11.8908,0.0339504,-1.30115e-05,2.4399e-09,-1.75595e-13,-25335.6,-36.6089], Tmin=(993.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone)"""),
)

species(
    label = 'CC=C(C)CC=O(2715)',
    structure = SMILES('CC=C(C)CC=O'),
    E0 = (-168.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05895,0.0565821,-2.25073e-05,-2.15824e-09,2.26109e-12,-20119.8,24.818], Tmin=(100,'K'), Tmax=(1373,'K')), NASAPolynomial(coeffs=[13.9466,0.0348215,-1.59793e-05,3.04577e-09,-2.11197e-13,-25146.6,-46.8612], Tmin=(1373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH)"""),
)

species(
    label = 'CCC(C)C=C=O(2463)',
    structure = SMILES('CCC(C)C=C=O'),
    E0 = (-158.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.726432,0.0655059,-4.46957e-05,1.57415e-08,-2.2782e-12,-18959.1,26.6532], Tmin=(100,'K'), Tmax=(1586.33,'K')), NASAPolynomial(coeffs=[14.212,0.0315017,-1.25423e-05,2.22893e-09,-1.48673e-13,-23237.6,-44.6294], Tmin=(1586.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d)"""),
)

species(
    label = 'C=CC(C)CC=O(2314)',
    structure = SMILES('C=CC(C)CC=O'),
    E0 = (-154.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.830208,0.0615302,-3.7351e-05,1.11924e-08,-1.3554e-12,-18503.7,27.7601], Tmin=(100,'K'), Tmax=(1878.52,'K')), NASAPolynomial(coeffs=[16.1214,0.0289701,-1.13517e-05,1.96555e-09,-1.2745e-13,-24248.6,-55.6516], Tmin=(1878.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH)"""),
)

species(
    label = 'CO(2039)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35309e-07,1.51269e-10,-9.88873e-15,-14292.7,6.51158], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'CC=C(C)C[C]=O(2716)',
    structure = SMILES('CC=C(C)C[C]=O'),
    E0 = (-8.25877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,237.634],'cm^-1')),
        HinderedRotor(inertia=(0.0028318,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.2945,'amu*angstrom^2'), symmetry=1, barrier=(11.9195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288366,'amu*angstrom^2'), symmetry=1, barrier=(11.9321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339467,'amu*angstrom^2'), symmetry=1, barrier=(11.9145,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995304,0.059308,-3.62871e-05,1.01549e-08,-1.11022e-12,-879.273,25.8764], Tmin=(100,'K'), Tmax=(2114.34,'K')), NASAPolynomial(coeffs=[20.9521,0.0215535,-9.50297e-06,1.70984e-09,-1.11695e-13,-9318.54,-85.3467], Tmin=(2114.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.25877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C(C)C=C=O(2717)',
    structure = SMILES('C[CH]C(C)C=C=O'),
    E0 = (35.8846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,308.584,308.584],'cm^-1')),
        HinderedRotor(inertia=(0.00177037,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00177033,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127888,'amu*angstrom^2'), symmetry=1, barrier=(8.64178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127889,'amu*angstrom^2'), symmetry=1, barrier=(8.64176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36562,0.0592537,-3.95434e-05,1.37604e-08,-2.04111e-12,4409.05,26.8974], Tmin=(100,'K'), Tmax=(1476.41,'K')), NASAPolynomial(coeffs=[10.128,0.0355141,-1.54246e-05,2.86966e-09,-1.96994e-13,1821.67,-18.79], Tmin=(1476.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.8846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Cs_S)"""),
)

species(
    label = 'C=CC(C)C[C]=O(2718)',
    structure = SMILES('C=CC(C)C[C]=O'),
    E0 = (5.10705,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,294.781,294.781],'cm^-1')),
        HinderedRotor(inertia=(0.00193998,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158066,'amu*angstrom^2'), symmetry=1, barrier=(9.74687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158069,'amu*angstrom^2'), symmetry=1, barrier=(9.74687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.158069,'amu*angstrom^2'), symmetry=1, barrier=(9.74686,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03891,0.0615754,-4.38931e-05,1.6588e-08,-2.61263e-12,723.49,27.8013], Tmin=(100,'K'), Tmax=(1454.03,'K')), NASAPolynomial(coeffs=[12.0242,0.0313552,-1.27177e-05,2.29427e-09,-1.55041e-13,-2471.13,-29.3089], Tmin=(1454.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.10705,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C]=O(1376)',
    structure = SMILES('[CH2][C]=O'),
    E0 = (160.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,623.763,623.847],'cm^-1')),
        HinderedRotor(inertia=(0.000434039,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.44947,0.00889162,4.91274e-06,-1.14995e-08,4.62261e-12,19378.1,9.99239], Tmin=(100,'K'), Tmax=(1030.77,'K')), NASAPolynomial(coeffs=[6.10651,0.00625992,-2.43247e-06,4.78594e-10,-3.54632e-14,18422.4,-4.88573], Tmin=(1030.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJC=O) + radical(CsCJ=O)"""),
)

species(
    label = 'CC=CC[C]=O(2431)',
    structure = SMILES('CC=CC[C]=O'),
    E0 = (30.7963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,387.023],'cm^-1')),
        HinderedRotor(inertia=(0.112079,'amu*angstrom^2'), symmetry=1, barrier=(11.8552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111689,'amu*angstrom^2'), symmetry=1, barrier=(11.8567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111223,'amu*angstrom^2'), symmetry=1, barrier=(11.856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83002,0.0398038,-9.11535e-06,-1.05206e-08,4.83654e-12,3788.49,22.0947], Tmin=(100,'K'), Tmax=(1207.67,'K')), NASAPolynomial(coeffs=[11.1225,0.0255305,-1.18873e-05,2.32622e-09,-1.6555e-13,340.438,-29.4728], Tmin=(1207.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.7963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][CH]C(1186)',
    structure = SMILES('C[CH][CH]C'),
    E0 = (244.869,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2670.04],'cm^-1')),
        HinderedRotor(inertia=(0.341695,'amu*angstrom^2'), symmetry=1, barrier=(7.85623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00155291,'amu*angstrom^2'), symmetry=1, barrier=(7.85531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0136081,'amu*angstrom^2'), symmetry=1, barrier=(68.8251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.1063,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62536,0.022997,-1.74547e-06,-3.25078e-09,7.61424e-13,29450.7,14.8448], Tmin=(100,'K'), Tmax=(2137.09,'K')), NASAPolynomial(coeffs=[11.4859,0.0201914,-8.13358e-06,1.34908e-09,-8.16544e-14,23371.9,-35.4092], Tmin=(2137.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(244.869,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]C[C]=O(2433)',
    structure = SMILES('C[CH][CH]C[C]=O'),
    E0 = (303.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,2057.55,2058.05],'cm^-1')),
        HinderedRotor(inertia=(0.00216421,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106956,'amu*angstrom^2'), symmetry=1, barrier=(5.89274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106314,'amu*angstrom^2'), symmetry=1, barrier=(5.90618,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108095,'amu*angstrom^2'), symmetry=1, barrier=(5.90437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8548,0.051267,-5.85796e-05,4.84354e-08,-1.73425e-11,36566.9,27.6321], Tmin=(100,'K'), Tmax=(806.443,'K')), NASAPolynomial(coeffs=[4.01809,0.0342252,-1.51418e-05,2.8214e-09,-1.93431e-13,36423.3,18.9333], Tmin=(806.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(303.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJC) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][C](C)C[C]=O(2719)',
    structure = SMILES('C[CH][C](C)C[C]=O'),
    E0 = (259.567,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3025,407.5,1350,352.5,1855,455,950,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.102702,'amu*angstrom^2'), symmetry=1, barrier=(2.36132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102692,'amu*angstrom^2'), symmetry=1, barrier=(2.36109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102669,'amu*angstrom^2'), symmetry=1, barrier=(2.36057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102644,'amu*angstrom^2'), symmetry=1, barrier=(2.35999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102694,'amu*angstrom^2'), symmetry=1, barrier=(2.36114,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05689,0.0753737,-0.000112452,1.09476e-07,-4.18903e-11,31314.4,30.916], Tmin=(100,'K'), Tmax=(836.988,'K')), NASAPolynomial(coeffs=[1.40105,0.049758,-2.35859e-05,4.47624e-09,-3.07506e-13,32096.4,34.3328], Tmin=(836.988,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.567,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C(C)[CH][C]=O(2720)',
    structure = SMILES('C[CH]C(C)[CH][C]=O'),
    E0 = (241.674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,1855,455,950,180,461.499],'cm^-1')),
        HinderedRotor(inertia=(0.0789565,'amu*angstrom^2'), symmetry=1, barrier=(1.81537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0789769,'amu*angstrom^2'), symmetry=1, barrier=(1.81583,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0789442,'amu*angstrom^2'), symmetry=1, barrier=(1.81508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0120114,'amu*angstrom^2'), symmetry=1, barrier=(1.81535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0120106,'amu*angstrom^2'), symmetry=1, barrier=(1.81526,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25858,0.0601113,-4.25222e-05,1.64913e-08,-2.74771e-12,29165.2,28.8169], Tmin=(100,'K'), Tmax=(1346.5,'K')), NASAPolynomial(coeffs=[9.69046,0.035063,-1.46185e-05,2.67593e-09,-1.82655e-13,26894.5,-14.3707], Tmin=(1346.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH]C)C[C]=O(2467)',
    structure = SMILES('[CH2]C([CH]C)C[C]=O'),
    E0 = (279.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,413.61,2621.97],'cm^-1')),
        HinderedRotor(inertia=(0.0627707,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627708,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627705,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0627707,'amu*angstrom^2'), symmetry=1, barrier=(7.62018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0141566,'amu*angstrom^2'), symmetry=1, barrier=(69.0623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05704,0.0686662,-7.37601e-05,5.27733e-08,-1.65495e-11,33685.5,31.5724], Tmin=(100,'K'), Tmax=(762.453,'K')), NASAPolynomial(coeffs=[6.46712,0.0401417,-1.73637e-05,3.21791e-09,-2.20793e-13,32864.7,6.96619], Tmin=(762.453,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C(C)C[C]=O(2721)',
    structure = SMILES('[CH2][CH]C(C)C[C]=O'),
    E0 = (279.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,1904.05,1904.05],'cm^-1')),
        HinderedRotor(inertia=(0.0839111,'amu*angstrom^2'), symmetry=1, barrier=(8.91661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00112577,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0839118,'amu*angstrom^2'), symmetry=1, barrier=(8.91662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0839121,'amu*angstrom^2'), symmetry=1, barrier=(8.91661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0839115,'amu*angstrom^2'), symmetry=1, barrier=(8.91662,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23781,0.0656711,-6.39225e-05,4.1548e-08,-1.23157e-11,33698,30.9379], Tmin=(100,'K'), Tmax=(784.457,'K')), NASAPolynomial(coeffs=[5.8029,0.0423938,-1.94138e-05,3.72327e-09,-2.61537e-13,32981.8,10.0222], Tmin=(784.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[C](C)C[C]=O(2470)',
    structure = SMILES('CC[C](C)C[C]=O'),
    E0 = (65.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885027,0.0763373,-0.0001,8.96318e-08,-3.31428e-11,7925.29,28.9786], Tmin=(100,'K'), Tmax=(820.06,'K')), NASAPolynomial(coeffs=[3.4108,0.0489797,-2.24536e-05,4.22973e-09,-2.90719e-13,8016.67,20.3772], Tmin=(820.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]C(C)C=C[O](1993)',
    structure = SMILES('C[CH]C(C)C=C[O]'),
    E0 = (55.1579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,332.073,332.343,332.741],'cm^-1')),
        HinderedRotor(inertia=(0.00152731,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00152477,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224775,'amu*angstrom^2'), symmetry=1, barrier=(17.6755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.225045,'amu*angstrom^2'), symmetry=1, barrier=(17.6618,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.797945,0.0562584,-3.6551e-06,-3.60961e-08,1.79525e-11,6761.83,29.6405], Tmin=(100,'K'), Tmax=(988.345,'K')), NASAPolynomial(coeffs=[15.8982,0.0278864,-1.02862e-05,1.89475e-09,-1.35493e-13,2177.86,-51.1223], Tmin=(988.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.1579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S)"""),
)

species(
    label = '[CH2]CC(C)C[C]=O(2473)',
    structure = SMILES('[CH2]CC(C)C[C]=O'),
    E0 = (84.8487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,180,2716.41],'cm^-1')),
        HinderedRotor(inertia=(0.130573,'amu*angstrom^2'), symmetry=1, barrier=(3.00213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0197619,'amu*angstrom^2'), symmetry=1, barrier=(13.0894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569572,'amu*angstrom^2'), symmetry=1, barrier=(13.0956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.569975,'amu*angstrom^2'), symmetry=1, barrier=(13.1049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13056,'amu*angstrom^2'), symmetry=1, barrier=(3.00184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09352,0.0665252,-5.19973e-05,2.34167e-08,-4.62755e-12,10307.3,28.8864], Tmin=(100,'K'), Tmax=(1145.29,'K')), NASAPolynomial(coeffs=[8.66729,0.0400732,-1.73528e-05,3.25029e-09,-2.2551e-13,8572.52,-8.68007], Tmin=(1145.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.8487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CCC(C)[CH][C]=O(2471)',
    structure = SMILES('CCC(C)[CH][C]=O'),
    E0 = (47.1318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.654174,0.0660354,-4.68376e-05,1.77137e-08,-2.76228e-12,5795.26,28.4414], Tmin=(100,'K'), Tmax=(1490.95,'K')), NASAPolynomial(coeffs=[13.6381,0.0312014,-1.17923e-05,2.04343e-09,-1.34718e-13,1923.57,-39.3846], Tmin=(1490.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.1318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(CC)C[C]=O(2266)',
    structure = SMILES('[CH2]C(CC)C[C]=O'),
    E0 = (84.6847,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,321.669,1667.43],'cm^-1')),
        HinderedRotor(inertia=(0.00162925,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11647,'amu*angstrom^2'), symmetry=1, barrier=(8.55107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116463,'amu*angstrom^2'), symmetry=1, barrier=(8.55097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116448,'amu*angstrom^2'), symmetry=1, barrier=(8.5511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116459,'amu*angstrom^2'), symmetry=1, barrier=(8.55118,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918467,0.0692796,-6.02271e-05,3.15739e-08,-7.15936e-12,10294.9,29.5123], Tmin=(100,'K'), Tmax=(1029.77,'K')), NASAPolynomial(coeffs=[8.74735,0.0388699,-1.59319e-05,2.89794e-09,-1.97739e-13,8682.5,-8.48733], Tmin=(1029.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.6847,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = 'CC=C(C)C[CH][O](1992)',
    structure = SMILES('CC=C(C)C[CH][O]'),
    E0 = (154.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,272.196,272.199,2169.83],'cm^-1')),
        HinderedRotor(inertia=(0.151655,'amu*angstrom^2'), symmetry=1, barrier=(7.97358,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151654,'amu*angstrom^2'), symmetry=1, barrier=(7.97348,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151658,'amu*angstrom^2'), symmetry=1, barrier=(7.9735,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151654,'amu*angstrom^2'), symmetry=1, barrier=(7.97349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.72487,0.079291,-0.000101422,8.83917e-08,-3.24364e-11,18701.8,27.5098], Tmin=(100,'K'), Tmax=(795.259,'K')), NASAPolynomial(coeffs=[4.43364,0.0491379,-2.28593e-05,4.35115e-09,-3.01556e-13,18475.5,12.7529], Tmin=(795.259,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(154.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([CH]C)CC=O(2472)',
    structure = SMILES('[CH2]C([CH]C)CC=O'),
    E0 = (119.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,180,1568.49],'cm^-1')),
        HinderedRotor(inertia=(0.0975018,'amu*angstrom^2'), symmetry=1, barrier=(2.24176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.097689,'amu*angstrom^2'), symmetry=1, barrier=(2.24606,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0980122,'amu*angstrom^2'), symmetry=1, barrier=(2.25349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0993873,'amu*angstrom^2'), symmetry=1, barrier=(2.28511,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0961202,'amu*angstrom^2'), symmetry=1, barrier=(2.20999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3169,0.0633593,-4.98606e-05,2.59941e-08,-6.48208e-12,14437.3,29.8303], Tmin=(100,'K'), Tmax=(892.182,'K')), NASAPolynomial(coeffs=[5.25569,0.0457003,-2.01709e-05,3.80905e-09,-2.65588e-13,13734.5,11.2772], Tmin=(892.182,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = 'C=CC(C)C[CH][O](1994)',
    structure = SMILES('C=CC(C)C[CH][O]'),
    E0 = (174.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,267.033,267.135,267.259,1282.08],'cm^-1')),
        HinderedRotor(inertia=(0.00236165,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00236229,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196822,'amu*angstrom^2'), symmetry=1, barrier=(9.96279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.196473,'amu*angstrom^2'), symmetry=1, barrier=(9.96201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14439,0.0677673,-5.56557e-05,2.82003e-08,-6.54396e-12,21040.6,27.7015], Tmin=(100,'K'), Tmax=(969.587,'K')), NASAPolynomial(coeffs=[6.68749,0.0448996,-2.02785e-05,3.87591e-09,-2.72178e-13,19965.7,1.13034], Tmin=(969.587,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCsJOH)"""),
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
    E0 = (74.1445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (522.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (274.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (319.745,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (532.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (720.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (540.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (539.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (82.4288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (137.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (137.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (99.1178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (131.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (212.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (258.669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (219.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (149.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (195.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (229.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (405.799,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (439.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (471.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (453.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (491.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (491.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (228.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (233.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (191.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (218.186,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (231.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (270.709,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (195.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (247.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['ketene(1375)', 'butene2t(396)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', 'C[CH]CC[C]=O(2248)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CC(C)[CH]C[C]=O(2713)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.11838e+10,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;CH3] for rate rule [cCs(-HC)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['[CH2]C(=O)C(C)[CH]C(2027)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHCH3(T)(21)', 'C[CH]C[C]=O(2361)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[C]=O(2355)', '[CH2]C(C)[CH]C(114)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', '[CH]C(C)C[C]=O(2389)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C[C]C(C)C[C]=O(2714)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['CC1CC(=O)C1C(2312)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['CC=C(C)CC=O(2715)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['CCC(C)C=C=O(2463)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['C=CC(C)CC=O(2314)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['CO(2039)', '[CH2]C(C)[CH]C(114)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(754.069,'m^3/(mol*s)'), n=1.0822, Ea=(24.7919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'CC=C(C)C[C]=O(2716)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.72e+07,'cm^3/(mol*s)'), n=1.64, Ea=(9.07928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2568 used for Cds-CsCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'C[CH]C(C)C=C=O(2717)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-CsH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'C=CC(C)C[C]=O(2718)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][C]=O(1376)', 'butene2t(396)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(0.00858789,'m^3/(mol*s)'), n=2.41179, Ea=(16.3987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(17)', 'CC=CC[C]=O(2431)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(10100,'cm^3/(mol*s)'), n=2.41, Ea=(28.4512,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 435 used for Cds-CsH_Cds-CsH;CsJ-HHH
Exact match found for rate rule [Cds-CsH_Cds-CsH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['ketene(1375)', 'C[CH][CH]C(1186)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.12155,'m^3/(mol*s)'), n=2.01066, Ea=(45.2043,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ck;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C]=O(1376)', 'C[CH][CH]C(1186)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.328e+08,'m^3/(mol*s)'), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CH3(17)', 'C[CH][CH]C[C]=O(2433)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.64e+07,'m^3/(mol*s)'), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', 'C[CH][C](C)C[C]=O(2719)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', 'C[CH]C(C)[CH][C]=O(2720)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH2][CH]C(C)C[C]=O(2721)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['CC[C](C)C[C]=O(2470)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 150 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['C[CH]C(C)C=C[O](1993)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(7.74568e+08,'s^-1'), n=1.384, Ea=(159.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['[CH2]CC(C)C[C]=O(2473)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(5.04e-11,'s^-1'), n=6.833, Ea=(117.248,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 33 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[CH]C(C)C[C]=O(2311)'],
    products = ['CCC(C)[CH][C]=O(2471)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.17661e+06,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(CC)C[C]=O(2266)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CC=C(C)C[CH][O](1992)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;Y_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([CH]C)CC=O(2472)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(463.959,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;XH_out] for rate rule [R4H_SSS;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=CC(C)C[CH][O](1994)'],
    products = ['C[CH]C(C)C[C]=O(2311)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #708',
    isomers = [
        'C[CH]C(C)C[C]=O(2311)',
    ],
    reactants = [
        ('ketene(1375)', 'butene2t(396)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #708',
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


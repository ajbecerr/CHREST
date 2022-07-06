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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3595.79,'J/mol'), sigma=(6.33206,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=561.65 K, Pc=32.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23781,0.0656711,-6.39225e-05,4.1548e-08,-1.23157e-11,33698,30.9379], Tmin=(100,'K'), Tmax=(784.457,'K')), NASAPolynomial(coeffs=[5.8029,0.0423938,-1.94138e-05,3.72327e-09,-2.61537e-13,32981.8,10.0222], Tmin=(784.457,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH][CH]CC[C]=O(2530)',
    structure = SMILES('C[CH][CH]CC[C]=O'),
    E0 = (274.201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,180,2760.26,4000],'cm^-1')),
        HinderedRotor(inertia=(0.122453,'amu*angstrom^2'), symmetry=1, barrier=(2.81543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122453,'amu*angstrom^2'), symmetry=1, barrier=(2.81543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122453,'amu*angstrom^2'), symmetry=1, barrier=(2.81543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122453,'amu*angstrom^2'), symmetry=1, barrier=(2.81543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.122453,'amu*angstrom^2'), symmetry=1, barrier=(2.81543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21179,0.0724473,-0.000107744,1.06164e-07,-4.05909e-11,33068.5,32.2496], Tmin=(100,'K'), Tmax=(854.784,'K')), NASAPolynomial(coeffs=[0.260936,0.0506181,-2.33227e-05,4.35547e-09,-2.95872e-13,34191.1,42.3034], Tmin=(854.784,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH]C)C(=C)[O](2503)',
    structure = SMILES('[CH2]C([CH]C)C(=C)[O]'),
    E0 = (250.817,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,180,1962.86],'cm^-1')),
        HinderedRotor(inertia=(0.0102752,'amu*angstrom^2'), symmetry=1, barrier=(2.7886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0102751,'amu*angstrom^2'), symmetry=1, barrier=(2.78855,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0373544,'amu*angstrom^2'), symmetry=1, barrier=(10.1375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0373544,'amu*angstrom^2'), symmetry=1, barrier=(10.1375,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.586634,0.0664197,-5.34853e-05,2.33575e-08,-4.13471e-12,30296.4,32.4745], Tmin=(100,'K'), Tmax=(1352.09,'K')), NASAPolynomial(coeffs=[14.0751,0.0265156,-9.2158e-06,1.52979e-09,-9.87714e-14,26648.9,-36.6682], Tmin=(1352.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cs_S) + radical(Isobutyl)"""),
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
    label = '[CH2][CH]C[C]=O(2359)',
    structure = SMILES('[CH2][CH]C[C]=O'),
    E0 = (338.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,2450.4],'cm^-1')),
        HinderedRotor(inertia=(0.000151284,'amu*angstrom^2'), symmetry=1, barrier=(1.71768,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.61568,'amu*angstrom^2'), symmetry=1, barrier=(37.1477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636846,'amu*angstrom^2'), symmetry=1, barrier=(14.6423,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46967,0.0353886,-2.95162e-05,1.48804e-08,-3.3285e-12,40706.3,22.2358], Tmin=(100,'K'), Tmax=(1023.26,'K')), NASAPolynomial(coeffs=[5.98538,0.0216456,-9.37056e-06,1.75551e-09,-1.21891e-13,39986.8,5.19366], Tmin=(1023.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJ) + radical(CCCJ=O)"""),
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
    label = '[CH]C([CH2])C[C]=O(2576)',
    structure = SMILES('[CH]C([CH2])C[C]=O'),
    E0 = (553.333,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1855,455,950,415.736,415.777,2757.17,2757.2],'cm^-1')),
        HinderedRotor(inertia=(0.606643,'amu*angstrom^2'), symmetry=1, barrier=(74.4107,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0137937,'amu*angstrom^2'), symmetry=1, barrier=(74.4115,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0909531,'amu*angstrom^2'), symmetry=1, barrier=(11.1568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0909614,'amu*angstrom^2'), symmetry=1, barrier=(11.157,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1005,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54162,0.0569827,-6.26091e-05,3.64919e-08,-7.1099e-12,66636.5,26.0472], Tmin=(100,'K'), Tmax=(707.203,'K')), NASAPolynomial(coeffs=[8.57968,0.0246536,-9.90087e-06,1.75835e-09,-1.1756e-13,65454.1,-6.79137], Tmin=(707.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.333,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([C]C)C[C]=O(4232)',
    structure = SMILES('[CH2]C([C]C)C[C]=O'),
    E0 = (532.9,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788995,0.0732268,-8.28009e-05,5.51873e-08,-1.52288e-11,64206.3,29.8836], Tmin=(100,'K'), Tmax=(875.177,'K')), NASAPolynomial(coeffs=[9.62657,0.0328358,-1.35754e-05,2.45637e-09,-1.6627e-13,62659.4,-11.5746], Tmin=(875.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.9,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C([CH]C)C[C]=O(4233)',
    structure = SMILES('[CH]C([CH]C)C[C]=O'),
    E0 = (522.36,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10512,0.0685392,-7.50787e-05,5.19402e-08,-1.56318e-11,62925.2,30.6291], Tmin=(100,'K'), Tmax=(788.628,'K')), NASAPolynomial(coeffs=[7.05989,0.0383359,-1.76304e-05,3.37605e-09,-2.36536e-13,61986,3.31487], Tmin=(788.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'CC1CC1C[C]=O(4234)',
    structure = SMILES('CC1CC1C[C]=O'),
    E0 = (29.8101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29639,0.0468232,6.06457e-06,-3.91189e-08,1.77282e-11,3693.84,24.9486], Tmin=(100,'K'), Tmax=(992.627,'K')), NASAPolynomial(coeffs=[13.1681,0.0284177,-1.06011e-05,1.94706e-09,-1.38234e-13,-113.074,-39.542], Tmin=(992.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.8101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1CC(=O)C1C(4224)',
    structure = SMILES('[CH2]C1CC(=O)C1C'),
    E0 = (26.7903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67031,0.0371845,3.20747e-05,-6.39039e-08,2.61588e-11,3318.59,24.2862], Tmin=(100,'K'), Tmax=(971.807,'K')), NASAPolynomial(coeffs=[11.5698,0.0307612,-1.09898e-05,1.98269e-09,-1.40227e-13,-226.275,-31.5292], Tmin=(971.807,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.7903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(Isobutyl)"""),
)

species(
    label = 'C[CH]C1CC(=O)C1(3816)',
    structure = SMILES('C[CH]C1CC(=O)C1'),
    E0 = (22.0107,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82639,0.0363309,2.35754e-05,-4.55775e-08,1.683e-11,2735.26,25.3375], Tmin=(100,'K'), Tmax=(1066.01,'K')), NASAPolynomial(coeffs=[9.0875,0.036534,-1.53346e-05,2.91132e-09,-2.06491e-13,-372.44,-17.4726], Tmin=(1066.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(22.0107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(374.151,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(Cs_S)"""),
)

species(
    label = 'C=C(CC)C[C]=O(2464)',
    structure = SMILES('C=C(CC)C[C]=O'),
    E0 = (5.12134,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,403.965,405.398],'cm^-1')),
        HinderedRotor(inertia=(0.00103039,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116029,'amu*angstrom^2'), symmetry=1, barrier=(13.4637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11592,'amu*angstrom^2'), symmetry=1, barrier=(13.4662,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11581,'amu*angstrom^2'), symmetry=1, barrier=(13.4677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.958161,0.0564238,-2.27169e-05,-5.86096e-09,4.51892e-12,734.097,27.0856], Tmin=(100,'K'), Tmax=(1176.87,'K')), NASAPolynomial(coeffs=[14.2068,0.0303204,-1.35698e-05,2.62263e-09,-1.86094e-13,-3694.99,-44.5577], Tmin=(1176.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.12134,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
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
    label = '[CH2]C(=CC)CC=O(4235)',
    structure = SMILES('[CH2]C(=CC)CC=O'),
    E0 = (-16.7202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03489,0.0545873,-1.66785e-05,-1.06809e-08,5.6934e-12,-1895.54,25.1695], Tmin=(100,'K'), Tmax=(1193.32,'K')), NASAPolynomial(coeffs=[13.844,0.0324107,-1.48971e-05,2.90201e-09,-2.0632e-13,-6430.69,-45.0841], Tmin=(1193.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.7202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(C=C=O)CC(2465)',
    structure = SMILES('[CH2]C(C=C=O)CC'),
    E0 = (52.4284,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2120,512.5,787.5,262.196,262.199],'cm^-1')),
        HinderedRotor(inertia=(0.222037,'amu*angstrom^2'), symmetry=1, barrier=(10.8319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00245211,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222039,'amu*angstrom^2'), symmetry=1, barrier=(10.8319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.222032,'amu*angstrom^2'), symmetry=1, barrier=(10.8319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.910917,0.0706994,-6.56844e-05,3.52187e-08,-8.04121e-12,6414.62,26.3981], Tmin=(100,'K'), Tmax=(1027.93,'K')), NASAPolynomial(coeffs=[9.65419,0.0366769,-1.60377e-05,3.02056e-09,-2.10457e-13,4617.11,-16.0241], Tmin=(1027.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.4284,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O)"""),
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
    label = '[CH2]C(C=C)CC=O(3818)',
    structure = SMILES('[CH2]C(C=C)CC=O'),
    E0 = (50.2288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11594,0.060483,-4.21294e-05,1.60085e-08,-2.57539e-12,6147.14,28.4863], Tmin=(100,'K'), Tmax=(1412.26,'K')), NASAPolynomial(coeffs=[10.7326,0.0332452,-1.31992e-05,2.35176e-09,-1.57846e-13,3430.91,-21.2281], Tmin=(1412.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.2288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
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
    collisionModel = TransportData(shapeIndex=1, epsilon=(815.652,'J/mol'), sigma=(3.65,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.95,'angstroms^3'), rotrelaxcollnum=1.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35309e-07,1.51269e-10,-9.88873e-15,-14292.7,6.51158], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2]C(=CC)C[C]=O(4236)',
    structure = SMILES('[CH2]C(=CC)C[C]=O'),
    E0 = (143.24,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,692.33],'cm^-1')),
        HinderedRotor(inertia=(0.771736,'amu*angstrom^2'), symmetry=1, barrier=(17.7437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00987755,'amu*angstrom^2'), symmetry=1, barrier=(3.36108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771813,'amu*angstrom^2'), symmetry=1, barrier=(17.7455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771914,'amu*angstrom^2'), symmetry=1, barrier=(17.7478,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.939437,0.057797,-3.24939e-05,4.50331e-09,1.09052e-12,17345.8,26.3324], Tmin=(100,'K'), Tmax=(1266.22,'K')), NASAPolynomial(coeffs=[14.9878,0.0270628,-1.22491e-05,2.35476e-09,-1.65549e-13,12694.3,-49.0782], Tmin=(1266.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH]C)C=C=O(4237)',
    structure = SMILES('[CH2]C([CH]C)C=C=O'),
    E0 = (246.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2120,512.5,787.5,222.787,222.788],'cm^-1')),
        HinderedRotor(inertia=(0.00339646,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0033964,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164013,'amu*angstrom^2'), symmetry=1, barrier=(5.77708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164022,'amu*angstrom^2'), symmetry=1, barrier=(5.77708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03554,0.0702019,-7.93901e-05,5.63167e-08,-1.72464e-11,29805.9,28.5113], Tmin=(100,'K'), Tmax=(778.101,'K')), NASAPolynomial(coeffs=[7.27229,0.0381424,-1.75904e-05,3.37074e-09,-2.36124e-13,28835.3,-0.0130208], Tmin=(778.101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Cs_S) + radical(CJC(C)C=C=O)"""),
)

species(
    label = '[CH2]C(C=C)C[C]=O(4238)',
    structure = SMILES('[CH2]C(C=C)C[C]=O'),
    E0 = (210.189,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,375.464,3213.34],'cm^-1')),
        HinderedRotor(inertia=(0.0996948,'amu*angstrom^2'), symmetry=1, barrier=(10.0089,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0996427,'amu*angstrom^2'), symmetry=1, barrier=(10.0153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.786842,'amu*angstrom^2'), symmetry=1, barrier=(78.8626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100131,'amu*angstrom^2'), symmetry=1, barrier=(10.0129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10296,0.0627311,-5.47256e-05,2.73666e-08,-5.73559e-12,25385.1,29.3539], Tmin=(100,'K'), Tmax=(1127.42,'K')), NASAPolynomial(coeffs=[10.1343,0.0306887,-1.20943e-05,2.15793e-09,-1.45717e-13,23348.6,-15.3004], Tmin=(1127.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(210.189,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = 'C=CC[C]=O(2318)',
    structure = SMILES('C=CC[C]=O'),
    E0 = (66.8219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,459.425],'cm^-1')),
        HinderedRotor(inertia=(0.100142,'amu*angstrom^2'), symmetry=1, barrier=(14.9135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.09966,'amu*angstrom^2'), symmetry=1, barrier=(14.9185,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51805,0.0238833,1.19498e-05,-2.85427e-08,1.09391e-11,8097.53,17.8097], Tmin=(100,'K'), Tmax=(1083.61,'K')), NASAPolynomial(coeffs=[9.78031,0.017858,-8.47809e-06,1.72443e-09,-1.27257e-13,5303.51,-23.4396], Tmin=(1083.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.8219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
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
    label = '[CH2][C]([CH]C)C[C]=O(4239)',
    structure = SMILES('[CH2][C]([CH]C)C[C]=O'),
    E0 = (464.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,2496.34,2499.05],'cm^-1')),
        HinderedRotor(inertia=(0.148494,'amu*angstrom^2'), symmetry=1, barrier=(3.41417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148264,'amu*angstrom^2'), symmetry=1, barrier=(3.40887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.147227,'amu*angstrom^2'), symmetry=1, barrier=(3.38504,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14886,'amu*angstrom^2'), symmetry=1, barrier=(3.42259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148958,'amu*angstrom^2'), symmetry=1, barrier=(3.42484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0349,0.0774909,-0.00012618,1.22922e-07,-4.54967e-11,55979.4,32.7783], Tmin=(100,'K'), Tmax=(870.885,'K')), NASAPolynomial(coeffs=[1.58887,0.0457046,-2.10663e-05,3.90143e-09,-2.62385e-13,56991.8,36.5488], Tmin=(870.885,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH]C)[CH][C]=O(4240)',
    structure = SMILES('[CH2]C([CH]C)[CH][C]=O'),
    E0 = (446.756,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,957.83,958.279],'cm^-1')),
        HinderedRotor(inertia=(0.0744504,'amu*angstrom^2'), symmetry=1, barrier=(1.71176,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749299,'amu*angstrom^2'), symmetry=1, barrier=(1.72279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0749341,'amu*angstrom^2'), symmetry=1, barrier=(1.72288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0744698,'amu*angstrom^2'), symmetry=1, barrier=(1.71221,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0748947,'amu*angstrom^2'), symmetry=1, barrier=(1.72198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2049,0.0625364,-5.72642e-05,3.16551e-08,-7.47846e-12,53832.1,30.8], Tmin=(100,'K'), Tmax=(999.296,'K')), NASAPolynomial(coeffs=[8.4989,0.0333397,-1.3438e-05,2.41688e-09,-1.63698e-13,52374.4,-4.38418], Tmin=(999.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.756,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C([CH2])C[C]=O(4179)',
    structure = SMILES('[CH2][CH]C([CH2])C[C]=O'),
    E0 = (484.473,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,180,1543.7],'cm^-1')),
        HinderedRotor(inertia=(0.135513,'amu*angstrom^2'), symmetry=1, barrier=(3.1157,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135936,'amu*angstrom^2'), symmetry=1, barrier=(3.12545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136415,'amu*angstrom^2'), symmetry=1, barrier=(3.13646,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00185875,'amu*angstrom^2'), symmetry=1, barrier=(3.13707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136822,'amu*angstrom^2'), symmetry=1, barrier=(3.14581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01958,0.0703602,-8.79602e-05,7.04702e-08,-2.35795e-11,58371.4,33.4881], Tmin=(100,'K'), Tmax=(820.932,'K')), NASAPolynomial(coeffs=[6.49036,0.037436,-1.6349e-05,3.01548e-09,-2.05119e-13,57684.4,9.46055], Tmin=(820.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.473,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C](CC)C[C]=O(2466)',
    structure = SMILES('[CH2][C](CC)C[C]=O'),
    E0 = (270.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3100,440,815,1455,1000,1855,455,950,180,3046.75],'cm^-1')),
        HinderedRotor(inertia=(0.151978,'amu*angstrom^2'), symmetry=1, barrier=(3.49427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152017,'amu*angstrom^2'), symmetry=1, barrier=(3.49516,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151921,'amu*angstrom^2'), symmetry=1, barrier=(3.49297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151989,'amu*angstrom^2'), symmetry=1, barrier=(3.49452,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.151961,'amu*angstrom^2'), symmetry=1, barrier=(3.49389,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.84503,0.0786801,-0.000114575,1.04244e-07,-3.72684e-11,32591.1,30.9046], Tmin=(100,'K'), Tmax=(866.786,'K')), NASAPolynomial(coeffs=[3.66509,0.0448084,-1.98638e-05,3.63798e-09,-2.44169e-13,32885.8,22.2224], Tmin=(866.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2]C([CH]C)[CH]C=O(4241)',
    structure = SMILES('[CH2]C([CH]C)C=C[O]'),
    E0 = (260.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.831755,0.057764,-1.55562e-05,-2.43502e-08,1.46598e-11,31424.5,31.3003], Tmin=(100,'K'), Tmax=(959.338,'K')), NASAPolynomial(coeffs=[15.5669,0.0247168,-8.27676e-06,1.44056e-09,-1.00385e-13,27290.8,-45.9857], Tmin=(959.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.24,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])C[C]=O(2469)',
    structure = SMILES('[CH2]CC([CH2])C[C]=O'),
    E0 = (289.931,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,229.445,1259.08],'cm^-1')),
        HinderedRotor(inertia=(0.00320221,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114217,'amu*angstrom^2'), symmetry=1, barrier=(4.26686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114217,'amu*angstrom^2'), symmetry=1, barrier=(4.26687,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114215,'amu*angstrom^2'), symmetry=1, barrier=(4.26686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114217,'amu*angstrom^2'), symmetry=1, barrier=(4.26685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.925285,0.0703288,-7.16022e-05,4.48254e-08,-1.19616e-11,34979,31.2773], Tmin=(100,'K'), Tmax=(892.658,'K')), NASAPolynomial(coeffs=[8.39977,0.0368361,-1.53227e-05,2.79457e-09,-1.90489e-13,33644.5,-3.93412], Tmin=(892.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH][C]=O)CC(2468)',
    structure = SMILES('[CH2]C([CH][C]=O)CC'),
    E0 = (252.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,327.303,2711.63],'cm^-1')),
        HinderedRotor(inertia=(1.01456,'amu*angstrom^2'), symmetry=1, barrier=(77.1271,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00157363,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168326,'amu*angstrom^2'), symmetry=1, barrier=(12.7961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00245239,'amu*angstrom^2'), symmetry=1, barrier=(12.7961,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01457,'amu*angstrom^2'), symmetry=1, barrier=(77.1272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79021,0.0664689,-5.56502e-05,2.64471e-08,-5.20762e-12,30453.4,29.7261], Tmin=(100,'K'), Tmax=(1206.53,'K')), NASAPolynomial(coeffs=[11.5566,0.0307744,-1.12727e-05,1.92592e-09,-1.26606e-13,27855.4,-24.2368], Tmin=(1206.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2][C]([CH]C)CC=O(4242)',
    structure = SMILES('[CH2][C]([CH]C)CC=O'),
    E0 = (304.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,2437.6,2438.68],'cm^-1')),
        HinderedRotor(inertia=(0.203118,'amu*angstrom^2'), symmetry=1, barrier=(4.67008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00109763,'amu*angstrom^2'), symmetry=1, barrier=(4.6399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202241,'amu*angstrom^2'), symmetry=1, barrier=(4.64991,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.72347,'amu*angstrom^2'), symmetry=1, barrier=(62.6179,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201995,'amu*angstrom^2'), symmetry=1, barrier=(4.64427,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12152,0.0743819,-0.000110771,1.085e-07,-4.14073e-11,36738.6,31.6479], Tmin=(100,'K'), Tmax=(850.151,'K')), NASAPolynomial(coeffs=[0.653027,0.0507774,-2.35864e-05,4.42357e-09,-3.01384e-13,37750.9,39.3174], Tmin=(850.151,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Tertalkyl) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][CH]C([CH2])CC=O(4184)',
    structure = SMILES('[CH2][CH]C([CH2])CC=O'),
    E0 = (324.513,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2782.5,750,1395,475,1775,1000,180,1199.2],'cm^-1')),
        HinderedRotor(inertia=(0.176573,'amu*angstrom^2'), symmetry=1, barrier=(4.05977,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176564,'amu*angstrom^2'), symmetry=1, barrier=(4.05955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00397902,'amu*angstrom^2'), symmetry=1, barrier=(4.06109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176523,'amu*angstrom^2'), symmetry=1, barrier=(4.0586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00397736,'amu*angstrom^2'), symmetry=1, barrier=(4.06151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7718,0.0570456,-2.14287e-05,-4.56858e-08,5.04468e-11,39102.1,30.0967], Tmin=(100,'K'), Tmax=(500.133,'K')), NASAPolynomial(coeffs=[5.11473,0.0433361,-1.93814e-05,3.66502e-09,-2.55097e-13,38604.8,14.6564], Tmin=(500.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (279.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (422.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (436.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (524.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (737.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (741.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (925.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (745.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (744.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (734.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (286.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (287.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (287.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (301.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (302.088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (342.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (342.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (342.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (304.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (304.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (336.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (366.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (469.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (424.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (413.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (410.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (303.319,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (434.33,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (611.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (676.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (658.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (696.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (416.231,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (420.715,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (436.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (437.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (423.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (397.166,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (420.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (398.132,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['ketene(1375)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH]C(C)C[C]=O(2721)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.29612e+11,'s^-1'), n=0.58375, Ea=(142.981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HC)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['C[CH][CH]CC[C]=O(2530)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2]C([CH]C)C(=C)[O](2503)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CHCH3(T)(21)', '[CH2][CH]C[C]=O(2359)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', 'C[CH][CH]C[C]=O(2433)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]=O(2355)', '[CH2]C([CH2])[CH]C(500)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH3(17)', '[CH]C([CH2])C[C]=O(2576)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2]C([C]C)C[C]=O(4232)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C([CH]C)C[C]=O(4233)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['CC1CC1C[C]=O(4234)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Cpri_rad_out_2H] for rate rule [R3_SS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2]C1CC(=O)C1C(4224)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['C[CH]C1CC(=O)C1(3816)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['C=C(CC)C[C]=O(2464)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['CC=C(C)C[C]=O(2716)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2]C(=CC)CC=O(4235)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2]C(C=C=O)CC(2465)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['C[CH]C(C)C=C=O(2717)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['C=CC(C)C[C]=O(2718)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2]C(C=C)CC=O(3818)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CO(2039)', '[CH2]C([CH2])[CH]C(500)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1508.14,'m^3/(mol*s)'), n=1.0822, Ea=(24.7919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;C_rad/H2/Cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2]C(=CC)C[C]=O(4236)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C([CH]C)C=C=O(4237)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-CsH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C(C=C)C[C]=O(4238)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH2(T)(20)', 'CC=CC[C]=O(2431)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(41.7,'m^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds;Y_1centerbirad] for rate rule [Cds-CsH_Cds-CsH;CH2_triplet]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction26',
    reactants = ['CHCH3(T)(21)', 'C=CC[C]=O(2318)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.17e+07,'cm^3/(mol*s)'), n=1.64, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.9 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=O(1376)', 'm1_allyl(186)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.0018001,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['ketene(1375)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(0.560775,'m^3/(mol*s)'), n=2.01066, Ea=(45.2043,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ck;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][C]=O(1376)', '[CH2][CH][CH]C(3856)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.64e+07,'m^3/(mol*s)'), n=-0.35, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_3R!H->O_Ext-1C-R_5R!H->C_Ext-1C-R_6R!H->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH2][C]([CH]C)C[C]=O(4239)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2]C([CH]C)[CH][C]=O(4240)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C[C]=O(4179)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2][C](CC)C[C]=O(2466)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(956916,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['C[CH][C](C)C[C]=O(2719)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(17481.2,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2]C([CH]C)[CH]C=O(4241)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]CC([CH2])C[C]=O(2469)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['[CH2]C([CH][C]=O)CC(2468)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.17661e+06,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C([CH]C)C[C]=O(2467)'],
    products = ['C[CH]C(C)[CH][C]=O(2720)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][C]([CH]C)CC=O(4242)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;Y_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][CH]C([CH2])CC=O(4184)'],
    products = ['[CH2]C([CH]C)C[C]=O(2467)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1197',
    isomers = [
        '[CH2]C([CH]C)C[C]=O(2467)',
    ],
    reactants = [
        ('ketene(1375)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1197',
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


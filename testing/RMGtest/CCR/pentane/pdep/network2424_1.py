species(
    label = 'C[CH]O[CH][C]=O(2799)',
    structure = SMILES('C[CH]O[CH][C]=O'),
    E0 = (175.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,332.399,332.4,332.4],'cm^-1')),
        HinderedRotor(inertia=(0.00152575,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203667,'amu*angstrom^2'), symmetry=1, barrier=(15.9691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203673,'amu*angstrom^2'), symmetry=1, barrier=(15.969,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.203674,'amu*angstrom^2'), symmetry=1, barrier=(15.9691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03458,0.0651052,-8.03133e-05,5.01695e-08,-1.22686e-11,21229.3,23.4124], Tmin=(100,'K'), Tmax=(1004.76,'K')), NASAPolynomial(coeffs=[13.4531,0.0156646,-6.50111e-06,1.19282e-09,-8.20296e-14,18733.8,-36.5585], Tmin=(1004.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCsJOCs) + radical(CsCJ=O)"""),
)

species(
    label = 'HCCO(2227)',
    structure = SMILES('[CH]=C=O'),
    E0 = (165.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,231.114,1089.61,3388.99],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47174,0.0113766,-1.1308e-05,6.13951e-09,-1.35087e-12,19887.1,6.87336], Tmin=(100,'K'), Tmax=(1096.2,'K')), NASAPolynomial(coeffs=[5.39217,0.00436893,-1.71893e-06,3.07751e-10,-2.08706e-14,19466,-2.56797], Tmin=(1096.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CH3CHO(1381)',
    structure = SMILES('CC=O'),
    E0 = (-177.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,192.062,1313.7,1313.74,1629.81,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0786982,'amu*angstrom^2'), symmetry=1, barrier=(2.09733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57993,0.00518977,2.269e-05,-2.73745e-08,9.28491e-12,-21369.7,8.9697], Tmin=(100,'K'), Tmax=(1028.8,'K')), NASAPolynomial(coeffs=[4.08562,0.0139062,-5.59372e-06,1.0461e-09,-7.38743e-14,-22039.1,3.76815], Tmin=(1028.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O][C]=C[O](9592)',
    structure = SMILES('[O][CH][C]=O'),
    E0 = (204.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,343.122],'cm^-1')),
        HinderedRotor(inertia=(0.528453,'amu*angstrom^2'), symmetry=1, barrier=(44.0956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.99411,0.0153541,4.13359e-06,-1.97935e-08,9.27968e-12,24622.8,13.8337], Tmin=(100,'K'), Tmax=(987.313,'K')), NASAPolynomial(coeffs=[10.2587,0.00223056,-7.0477e-07,2.03549e-10,-2.00526e-14,22393.5,-25.1462], Tmin=(987.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=OCOJ) + radical(OCJC=O) + radical(CsCJ=O)"""),
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
    label = '[CH]=[C][O](6861)',
    structure = SMILES('[CH][C]=O'),
    E0 = (425.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1491.4,1492.35],'cm^-1')),
        HinderedRotor(inertia=(0.0736325,'amu*angstrom^2'), symmetry=1, barrier=(5.08901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62981,0.00805288,-4.35162e-06,1.03741e-09,-8.43183e-14,51134.2,10.4631], Tmin=(100,'K'), Tmax=(1941.7,'K')), NASAPolynomial(coeffs=[6.62042,0.00292938,-1.19496e-06,2.28726e-10,-1.56225e-14,49777.3,-6.4529], Tmin=(1941.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = 'C[CH][O](1382)',
    structure = SMILES('C[CH][O]'),
    E0 = (156.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2362.36],'cm^-1')),
        HinderedRotor(inertia=(0.127216,'amu*angstrom^2'), symmetry=1, barrier=(2.92495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77075,0.0113794,7.94194e-07,-3.17974e-09,7.41421e-13,18809.5,9.6394], Tmin=(100,'K'), Tmax=(1908.67,'K')), NASAPolynomial(coeffs=[6.95159,0.0119688,-5.37106e-06,9.65305e-10,-6.23696e-14,16273.7,-11.2245], Tmin=(1908.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCsJOH)"""),
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
    label = '[CH]O[CH][C]=O(10368)',
    structure = SMILES('[CH]OC=[C][O]'),
    E0 = (485.443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1685,370,367.094,367.099,367.099,367.1,367.102,367.103],'cm^-1')),
        HinderedRotor(inertia=(0.0012509,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.205819,'amu*angstrom^2'), symmetry=1, barrier=(19.6825,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77653,0.0366362,-1.26809e-05,-2.85482e-08,1.89892e-11,58476.9,19.0646], Tmin=(100,'K'), Tmax=(903.754,'K')), NASAPolynomial(coeffs=[19.393,-0.00788021,5.68009e-06,-1.13371e-09,7.548e-14,53926.5,-71.7006], Tmin=(903.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(485.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(CH2_triplet)"""),
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
    label = 'C[C]O[CH][C]=O(10369)',
    structure = SMILES('C[C]OC=[C][O]'),
    E0 = (447.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01849,'amu*angstrom^2'), symmetry=1, barrier=(23.4171,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01878,'amu*angstrom^2'), symmetry=1, barrier=(23.4237,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01902,'amu*angstrom^2'), symmetry=1, barrier=(23.4293,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.188017,0.0668591,-7.33695e-05,3.59667e-08,-6.36241e-12,54042,26.8695], Tmin=(100,'K'), Tmax=(1666.66,'K')), NASAPolynomial(coeffs=[21.1919,-0.00130803,3.15089e-06,-7.09469e-10,4.92292e-14,49256.4,-80.1747], Tmin=(1666.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(447.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CH2_triplet) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]O[C]=[C][O](10217)',
    structure = SMILES('C[CH]O[C]=[C][O]'),
    E0 = (420.579,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1670,1700,300,440,300.333,300.535,300.776,301.559],'cm^-1')),
        HinderedRotor(inertia=(0.265318,'amu*angstrom^2'), symmetry=1, barrier=(17.2448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.271062,'amu*angstrom^2'), symmetry=1, barrier=(17.2688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.268864,'amu*angstrom^2'), symmetry=1, barrier=(17.2696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.589542,0.0639086,-7.4945e-05,4.09337e-08,-8.31241e-12,50716.3,27.9239], Tmin=(100,'K'), Tmax=(1358.65,'K')), NASAPolynomial(coeffs=[18.6799,0.00306693,5.97227e-07,-2.40979e-10,1.9776e-14,46500.4,-62.3206], Tmin=(1358.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.579,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'CC1OC1[C]=O(10370)',
    structure = SMILES('CC1OC1[C]=O'),
    E0 = (-60.1037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12412,0.0273497,3.4483e-05,-7.7296e-08,3.69651e-11,-7147.93,20.0102], Tmin=(100,'K'), Tmax=(869.832,'K')), NASAPolynomial(coeffs=[15.1344,0.00676723,2.29768e-06,-7.56671e-10,5.83477e-14,-10896,-49.4769], Tmin=(869.832,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.1037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCCJ=O)"""),
)

species(
    label = 'C=COC[C]=O(2338)',
    structure = SMILES('C=COC[C]=O'),
    E0 = (-97.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04711,'amu*angstrom^2'), symmetry=1, barrier=(24.0752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04769,'amu*angstrom^2'), symmetry=1, barrier=(24.0885,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04687,'amu*angstrom^2'), symmetry=1, barrier=(24.0697,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40898,0.0418256,3.06958e-06,-4.60594e-08,2.39348e-11,-11610.3,20.8716], Tmin=(100,'K'), Tmax=(937.289,'K')), NASAPolynomial(coeffs=[18.8888,0.00559264,-3.41686e-07,3.69732e-11,-8.34689e-15,-16572.2,-71.3162], Tmin=(937.289,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCJ=O)"""),
)

species(
    label = 'C=CO[CH][C]=O(2661)',
    structure = SMILES('C=CO[CH][C]=O'),
    E0 = (96.5024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05539,'amu*angstrom^2'), symmetry=1, barrier=(24.2654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0536,'amu*angstrom^2'), symmetry=1, barrier=(24.2243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05436,'amu*angstrom^2'), symmetry=1, barrier=(24.2419,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04513,0.0511737,-2.35928e-05,-2.37905e-08,1.77179e-11,11725.8,21.6517], Tmin=(100,'K'), Tmax=(924.083,'K')), NASAPolynomial(coeffs=[21.3685,-0.000848472,2.49654e-06,-5.13079e-10,3.10353e-14,6434.75,-83.0972], Tmin=(924.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.5024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(CsCJ=O)"""),
)

species(
    label = 'OCHCO(3676)',
    structure = SMILES('O=[C]C=O'),
    E0 = (-75.5464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,180,525.376,1512.41,1512.65,1513.37],'cm^-1')),
        HinderedRotor(inertia=(0.00619061,'amu*angstrom^2'), symmetry=1, barrier=(0.260399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25386,0.0154074,-1.14326e-05,4.29104e-09,-6.58698e-13,-9058.53,11.1539], Tmin=(100,'K'), Tmax=(1504.29,'K')), NASAPolynomial(coeffs=[6.43153,0.00695778,-3.00696e-06,5.56991e-10,-3.81281e-14,-10014.6,-5.47402], Tmin=(1504.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]O[CH][C]=O(2797)',
    structure = SMILES('[CH2][CH]O[CH][C]=O'),
    E0 = (387.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,376.221,376.221,376.222],'cm^-1')),
        HinderedRotor(inertia=(0.00119101,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162859,'amu*angstrom^2'), symmetry=1, barrier=(16.3579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162859,'amu*angstrom^2'), symmetry=1, barrier=(16.3579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162859,'amu*angstrom^2'), symmetry=1, barrier=(16.3579,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3432.58,'J/mol'), sigma=(5.86884,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.16 K, Pc=38.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.842119,0.0707708,-0.000101749,7.15637e-08,-1.94092e-11,46683.2,25.4299], Tmin=(100,'K'), Tmax=(913.168,'K')), NASAPolynomial(coeffs=[14.2955,0.0118385,-4.942e-06,8.87288e-10,-5.94444e-14,44226.2,-38.2527], Tmin=(913.168,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(387.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCsJOCs) + radical(CJCO) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]CO[CH][C]=O(2798)',
    structure = SMILES('[CH2]COC=[C][O]'),
    E0 = (198.496,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,264.998,265.638,265.684,265.957],'cm^-1')),
        HinderedRotor(inertia=(0.427169,'amu*angstrom^2'), symmetry=1, barrier=(21.5405,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.428673,'amu*angstrom^2'), symmetry=1, barrier=(21.5421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.430742,'amu*angstrom^2'), symmetry=1, barrier=(21.5346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.985411,0.0533117,-2.46087e-05,-2.27174e-08,1.74139e-11,23994.2,24.5369], Tmin=(100,'K'), Tmax=(914.987,'K')), NASAPolynomial(coeffs=[20.4274,0.00288864,1.37984e-06,-3.60174e-10,2.29684e-14,18989.3,-75.4398], Tmin=(914.987,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.496,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]OC[C]=O(2325)',
    structure = SMILES('[CH2][CH]OC[C]=O'),
    E0 = (206.754,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180,180,714.505],'cm^-1')),
        HinderedRotor(inertia=(0.0085066,'amu*angstrom^2'), symmetry=1, barrier=(3.0817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.534092,'amu*angstrom^2'), symmetry=1, barrier=(12.2798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.534092,'amu*angstrom^2'), symmetry=1, barrier=(12.2798,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0605882,'amu*angstrom^2'), symmetry=1, barrier=(21.9495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3432.58,'J/mol'), sigma=(5.86884,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.16 K, Pc=38.53 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2729,0.0621587,-7.66967e-05,4.98094e-08,-1.29233e-11,24963.2,23.7497], Tmin=(100,'K'), Tmax=(939.443,'K')), NASAPolynomial(coeffs=[11.2471,0.0196905,-8.88864e-06,1.69033e-09,-1.18188e-13,23089.1,-23.7469], Tmin=(939.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CJCO) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][CH]OC=C[O](1927)',
    structure = SMILES('[CH2][CH]OC=C[O]'),
    E0 = (152.679,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13522,'amu*angstrom^2'), symmetry=1, barrier=(26.101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13573,'amu*angstrom^2'), symmetry=1, barrier=(26.1127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13566,'amu*angstrom^2'), symmetry=1, barrier=(26.1111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.28393,0.080444,-9.0962e-05,4.51811e-08,-7.97971e-12,18584.4,29.7525], Tmin=(100,'K'), Tmax=(1707.69,'K')), NASAPolynomial(coeffs=[24.2187,-0.00453177,5.84919e-06,-1.26816e-09,8.7385e-14,13554.5,-96.155], Tmin=(1707.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.679,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(CJCO)"""),
)

species(
    label = 'C[CH]OC=C=O(3353)',
    structure = SMILES('C[CH]OC=C=O'),
    E0 = (-18.8964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.920742,'amu*angstrom^2'), symmetry=1, barrier=(21.1697,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92171,'amu*angstrom^2'), symmetry=1, barrier=(21.1919,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.921327,'amu*angstrom^2'), symmetry=1, barrier=(21.1831,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.105998,0.0681407,-7.49637e-05,3.82217e-08,-7.17648e-12,-2117.49,22.8036], Tmin=(100,'K'), Tmax=(1515.96,'K')), NASAPolynomial(coeffs=[20.421,0.00283552,1.23297e-06,-3.79019e-10,2.89686e-14,-6932.16,-79.2214], Tmin=(1515.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-18.8964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CCsJOC(O))"""),
)

species(
    label = 'CC1OC=C1[O](10371)',
    structure = SMILES('CC1OC=C1[O]'),
    E0 = (-114.251,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67667,0.0236268,7.50536e-05,-1.3622e-07,6.02268e-11,-13631.9,17.6756], Tmin=(100,'K'), Tmax=(915.531,'K')), NASAPolynomial(coeffs=[25.0183,-0.00671159,7.38139e-06,-1.47074e-09,9.20997e-14,-20908.4,-109.272], Tmin=(915.531,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=COC=C[O](2765)',
    structure = SMILES('C=COC=C[O]'),
    E0 = (-97.5162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38427,0.0384861,2.24397e-05,-7.41186e-08,3.63062e-11,-11616.6,21.6985], Tmin=(100,'K'), Tmax=(916.099,'K')), NASAPolynomial(coeffs=[21.5709,0.000258613,3.304e-06,-7.17302e-10,4.45496e-14,-17409.7,-85.3531], Tmin=(916.099,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-97.5162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]OC#C[O](10372)',
    structure = SMILES('C[CH]O[C]=C=O'),
    E0 = (220.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1685,370,2120,512.5,787.5,210.751,210.815,210.86],'cm^-1')),
        HinderedRotor(inertia=(0.506481,'amu*angstrom^2'), symmetry=1, barrier=(15.9779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.50667,'amu*angstrom^2'), symmetry=1, barrier=(15.9793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.506602,'amu*angstrom^2'), symmetry=1, barrier=(15.9768,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0181,0.0617168,-7.36422e-05,4.09135e-08,-8.17425e-12,26672.7,22.6524], Tmin=(100,'K'), Tmax=(944.399,'K')), NASAPolynomial(coeffs=[15.7192,0.00815171,-2.38473e-06,3.6785e-10,-2.36912e-14,23508,-49.4853], Tmin=(944.399,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CCsJOC(O)) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]O[C]=C[O](10373)',
    structure = SMILES('C[CH]O[C]=C[O]'),
    E0 = (180.835,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.959993,'amu*angstrom^2'), symmetry=1, barrier=(22.0721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958195,'amu*angstrom^2'), symmetry=1, barrier=(22.0308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.958025,'amu*angstrom^2'), symmetry=1, barrier=(22.0269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.372209,0.0707098,-7.67305e-05,3.78118e-08,-6.7487e-12,21928.8,28.2688], Tmin=(100,'K'), Tmax=(1648.76,'K')), NASAPolynomial(coeffs=[21.4681,0.000503878,2.80736e-06,-6.83452e-10,4.88076e-14,17067.4,-80.9208], Tmin=(1648.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(C=CJO)"""),
)

species(
    label = 'CCO[C]=[C][O](10374)',
    structure = SMILES('CCO[C]=[C][O]'),
    E0 = (226.651,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1670,1700,300,440,390.522,390.6,390.762,391.278],'cm^-1')),
        HinderedRotor(inertia=(0.149455,'amu*angstrom^2'), symmetry=1, barrier=(16.2084,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149077,'amu*angstrom^2'), symmetry=1, barrier=(16.2078,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149455,'amu*angstrom^2'), symmetry=1, barrier=(16.2078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34198,0.0499585,-3.20844e-05,-2.52079e-09,7.03461e-12,27363.3,25.7511], Tmin=(100,'K'), Tmax=(941.543,'K')), NASAPolynomial(coeffs=[15.4991,0.0107505,-2.97545e-06,4.85329e-10,-3.43828e-14,23769.4,-46.6244], Tmin=(941.543,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(226.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]O[C]=[C]O(10375)',
    structure = SMILES('C[CH]O[C]=[C]O'),
    E0 = (279.116,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.772838,'amu*angstrom^2'), symmetry=1, barrier=(17.7691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.773105,'amu*angstrom^2'), symmetry=1, barrier=(17.7752,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.772409,'amu*angstrom^2'), symmetry=1, barrier=(17.7592,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.772835,'amu*angstrom^2'), symmetry=1, barrier=(17.769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0675084,0.0738253,-8.89219e-05,4.88692e-08,-9.81492e-12,33729.8,28.8686], Tmin=(100,'K'), Tmax=(1423.37,'K')), NASAPolynomial(coeffs=[21.4411,3.00838e-05,2.91614e-06,-7.35211e-10,5.50301e-14,28959.3,-77.7403], Tmin=(1423.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.116,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCsJOC(O)) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]OC=[C]O(10376)',
    structure = SMILES('[CH2][CH]OC=[C]O'),
    E0 = (250.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.937307,'amu*angstrom^2'), symmetry=1, barrier=(21.5505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.937777,'amu*angstrom^2'), symmetry=1, barrier=(21.5613,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.936425,'amu*angstrom^2'), symmetry=1, barrier=(21.5303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.936365,'amu*angstrom^2'), symmetry=1, barrier=(21.5289,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.912758,0.0828938,-0.000101299,5.43865e-08,-1.04498e-11,30382.2,30.1049], Tmin=(100,'K'), Tmax=(1537.76,'K')), NASAPolynomial(coeffs=[25.0146,-0.00605176,6.43814e-06,-1.41461e-09,1.00402e-13,24950.7,-97.87], Tmin=(1537.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCsJOC(O)) + radical(CJCO) + radical(C=CJO)"""),
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
    E0 = (175.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (603.657,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (636.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (677.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (659.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (632.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (183.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (200.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (314.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (283.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (305.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (599.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (356.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (337.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (217.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (175.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (183.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (200.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (447.75,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (384.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (380.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (429.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (343.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]O[CH][C]=O(2799)'],
    products = ['HCCO(2227)', 'CH3CHO(1381)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O][C]=C[O](9592)', 'CHCH3(T)(21)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C][O](6861)', 'C[CH][O](1382)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH3(17)', '[CH]O[CH][C]=O(10368)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C[C]O[CH][C]=O(10369)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C[CH]O[C]=[C][O](10217)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C[CH]O[CH][C]=O(2799)'],
    products = ['CC1OC1[C]=O(10370)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[CH]O[CH][C]=O(2799)'],
    products = ['C=COC[C]=O(2338)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C=CO[CH][C]=O(2661)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.67e+12,'cm^3/(mol*s)'), n=0.1, Ea=(6.4601,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2816 used for Cds-HH_Cds-OsH;HJ
Exact match found for rate rule [Cds-HH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C][O](6861)', 'CH3CHO(1381)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4e+09,'cm^3/(mol*s)'), n=1.39, Ea=(35.8862,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['OCHCO(3676)', 'CHCH3(T)(21)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.86749e+08,'m^3/(mol*s)'), n=-0.730906, Ea=(37.2121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;Y_1centerbirad] + [Od_R;YJ] for rate rule [Od_R;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH2][CH]O[CH][C]=O(2797)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CO[CH][C]=O(2798)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.7e+13,'s^-1','+|-',2), n=-0.1, Ea=(158.364,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""From training reaction 347 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeO]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH]OC[C]=O(2325)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH]O[CH][C]=O(2799)'],
    products = ['[CH2][CH]OC=C[O](1927)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(505536,'s^-1'), n=1.7378, Ea=(41.5716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R5Hall;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH]O[CH][C]=O(2799)'],
    products = ['C[CH]OC=C=O(3353)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH]O[CH][C]=O(2799)'],
    products = ['CC1OC=C1[O](10371)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[CH]O[CH][C]=O(2799)'],
    products = ['C=COC=C[O](2765)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', 'C[CH]OC#C[O](10372)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C[CH]O[C]=C[O](10373)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CCO[C]=[C][O](10374)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.4016e+09,'s^-1'), n=1.0875, Ea=(153.657,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;Cs_H_out_H/NonDeC] + [R3H_SS_O;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R3H_SS_O;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C[CH]O[C]=[C]O(10375)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]OC=[C]O(10376)'],
    products = ['C[CH]O[CH][C]=O(2799)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.54044,'s^-1'), n=3.06876, Ea=(92.3667,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_2H;O_H_out] + [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6Hall;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2424',
    isomers = [
        'C[CH]O[CH][C]=O(2799)',
    ],
    reactants = [
        ('HCCO(2227)', 'CH3CHO(1381)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2424',
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


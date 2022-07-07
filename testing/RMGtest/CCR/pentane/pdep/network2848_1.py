species(
    label = 'C[CH]OO[CH][C]=O(11180)',
    structure = SMILES('C[CH]OO[CH][C]=O'),
    E0 = (190.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,289.862],'cm^-1')),
        HinderedRotor(inertia=(0.116785,'amu*angstrom^2'), symmetry=1, barrier=(6.96318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768479,'amu*angstrom^2'), symmetry=1, barrier=(45.8161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768467,'amu*angstrom^2'), symmetry=1, barrier=(45.8161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0258579,'amu*angstrom^2'), symmetry=1, barrier=(45.816,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.768482,'amu*angstrom^2'), symmetry=1, barrier=(45.8161,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30607,0.0665916,-7.68893e-05,4.48171e-08,-6.53795e-12,22980,25.3222], Tmin=(100,'K'), Tmax=(582.369,'K')), NASAPolynomial(coeffs=[7.30416,0.0344505,-1.74313e-05,3.45671e-09,-2.46293e-13,22127.8,-1.69082], Tmin=(582.369,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CCsJOOC) + radical(OCJC=O) + radical(CsCJ=O)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25386,0.0154074,-1.14326e-05,4.29104e-09,-6.58698e-13,-9058.53,11.1539], Tmin=(100,'K'), Tmax=(1504.29,'K')), NASAPolynomial(coeffs=[6.43153,0.00695778,-3.00696e-06,5.56991e-10,-3.81281e-14,-10014.6,-5.47402], Tmin=(1504.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[O][C]=CO[O](9542)',
    structure = SMILES('[O]O[CH][C]=O'),
    E0 = (193.505,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.0580088,'amu*angstrom^2'), symmetry=1, barrier=(10.0718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.45491,'amu*angstrom^2'), symmetry=1, barrier=(56.4433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0275,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51037,0.036662,-6.11317e-05,5.4797e-08,-1.90871e-11,23323.2,17.3301], Tmin=(100,'K'), Tmax=(841.573,'K')), NASAPolynomial(coeffs=[5.80665,0.0135881,-6.80396e-06,1.30268e-09,-8.93697e-14,23030.7,3.55435], Tmin=(841.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(193.505,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(OCJC=O) + radical(CsCJ=O)"""),
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
    label = 'C[CH]O[O](57)',
    structure = SMILES('C[CH]O[O]'),
    E0 = (163.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5],'cm^-1')),
        HinderedRotor(inertia=(0.385786,'amu*angstrom^2'), symmetry=1, barrier=(8.86998,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00119371,'amu*angstrom^2'), symmetry=1, barrier=(8.85907,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (60.052,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.06777,0.0231618,-1.22527e-05,-6.92044e-09,1.03446e-11,19636.5,14.5222], Tmin=(100,'K'), Tmax=(545.361,'K')), NASAPolynomial(coeffs=[4.49941,0.0169793,-7.12428e-06,1.32855e-09,-9.21409e-14,19416.2,7.89469], Tmin=(545.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCsJOOH)"""),
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
    label = '[CH]OO[CH][C]=O(12729)',
    structure = SMILES('[CH]OO[CH][C]=O'),
    E0 = (502.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,1855,455,950,383.068,383.129,1847.05],'cm^-1')),
        HinderedRotor(inertia=(0.352905,'amu*angstrom^2'), symmetry=1, barrier=(36.7565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0151839,'amu*angstrom^2'), symmetry=1, barrier=(36.7575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352856,'amu*angstrom^2'), symmetry=1, barrier=(36.7567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0769794,'amu*angstrom^2'), symmetry=1, barrier=(8.01921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (86.0462,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58812,0.0581888,-9.26903e-05,7.9565e-08,-2.76254e-11,60508.3,21.9696], Tmin=(100,'K'), Tmax=(718.579,'K')), NASAPolynomial(coeffs=[8.30445,0.0199585,-1.11257e-05,2.25941e-09,-1.61798e-13,59564.9,-8.0617], Tmin=(718.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(502.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(166.289,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-OsHHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CsCJ=O) + radical(CH2_triplet)"""),
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
    label = 'C[C]OO[CH][C]=O(12730)',
    structure = SMILES('C[C]OO[CH][C]=O'),
    E0 = (464.866,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,311.778,312.136],'cm^-1')),
        HinderedRotor(inertia=(0.139537,'amu*angstrom^2'), symmetry=1, barrier=(9.6295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.302897,'amu*angstrom^2'), symmetry=1, barrier=(20.9195,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.566366,'amu*angstrom^2'), symmetry=1, barrier=(39.1279,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.566072,'amu*angstrom^2'), symmetry=1, barrier=(39.1269,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56654,'amu*angstrom^2'), symmetry=1, barrier=(39.1277,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.982733,0.072434,-9.73887e-05,7.06122e-08,-2.09123e-11,56013.7,24.8913], Tmin=(100,'K'), Tmax=(817.768,'K')), NASAPolynomial(coeffs=[10.2413,0.0271458,-1.43164e-05,2.88786e-09,-2.07755e-13,54499.4,-17.9132], Tmin=(817.768,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(OCJC=O) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = 'C[CH]OO[C]=[C][O](12532)',
    structure = SMILES('C[CH]OO[C][C]=O'),
    E0 = (512.821,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,339.402,339.891],'cm^-1')),
        HinderedRotor(inertia=(0.0976028,'amu*angstrom^2'), symmetry=1, barrier=(7.98424,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.427111,'amu*angstrom^2'), symmetry=1, barrier=(34.8188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173183,'amu*angstrom^2'), symmetry=1, barrier=(14.198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.426588,'amu*angstrom^2'), symmetry=1, barrier=(34.8131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.425705,'amu*angstrom^2'), symmetry=1, barrier=(34.814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04794,0.0708775,-9.35404e-05,6.62157e-08,-1.91494e-11,61779.3,25.1964], Tmin=(100,'K'), Tmax=(836.315,'K')), NASAPolynomial(coeffs=[10.3164,0.0265484,-1.40347e-05,2.83953e-09,-2.04806e-13,60228.9,-17.8622], Tmin=(836.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(512.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CCsJOOC) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = 'CC1OOC1[C]=O(11181)',
    structure = SMILES('CC1OOC1[C]=O'),
    E0 = (-8.80742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61139,0.0435321,-1.51232e-05,-1.10312e-08,6.96168e-12,-965.527,22.0674], Tmin=(100,'K'), Tmax=(1045.04,'K')), NASAPolynomial(coeffs=[12.5366,0.0202434,-8.29064e-06,1.57568e-09,-1.12765e-13,-4260.75,-35.9624], Tmin=(1045.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.80742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(12dioxetane) + radical(CCCJ=O)"""),
)

species(
    label = 'C=COOC[C]=O(11185)',
    structure = SMILES('C=COOC[C]=O'),
    E0 = (8.30988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07373,0.0586722,-5.13011e-05,2.21318e-08,-3.79139e-12,1109.42,26.7443], Tmin=(100,'K'), Tmax=(1397.17,'K')), NASAPolynomial(coeffs=[15.1102,0.0184868,-8.15813e-06,1.54593e-09,-1.07903e-13,-2812.85,-45.6679], Tmin=(1397.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.30988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][CH]OOC=C=O(11124)',
    structure = SMILES('[CH2][CH]OOC=C=O'),
    E0 = (293.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.0608961,'amu*angstrom^2'), symmetry=1, barrier=(50.851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.51292,'amu*angstrom^2'), symmetry=1, barrier=(11.793,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.179714,'amu*angstrom^2'), symmetry=1, barrier=(11.7703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.70815,'amu*angstrom^2'), symmetry=1, barrier=(39.2738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04812,0.0687846,-9.02645e-05,6.29434e-08,-1.76281e-11,35380.5,26.9081], Tmin=(100,'K'), Tmax=(870.187,'K')), NASAPolynomial(coeffs=[11.0699,0.022716,-1.08504e-05,2.10072e-09,-1.47812e-13,33636.4,-20.0472], Tmin=(870.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CCsJOOC) + radical(CJCOOH)"""),
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
    label = '[CH2][CH]OO[CH][C]=O(11110)',
    structure = SMILES('[CH2][CH]OO[CH][C]=O'),
    E0 = (404.28,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,1191.32],'cm^-1')),
        HinderedRotor(inertia=(0.208399,'amu*angstrom^2'), symmetry=1, barrier=(4.7915,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00315,'amu*angstrom^2'), symmetry=1, barrier=(46.0563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00293,'amu*angstrom^2'), symmetry=1, barrier=(46.0513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.00295,'amu*angstrom^2'), symmetry=1, barrier=(46.0517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0457296,'amu*angstrom^2'), symmetry=1, barrier=(46.0505,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3579.64,'J/mol'), sigma=(6.06667,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.13 K, Pc=36.38 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.919297,0.073907,-0.000110967,9.44225e-08,-3.27548e-11,48728.8,28.6598], Tmin=(100,'K'), Tmax=(756.818,'K')), NASAPolynomial(coeffs=[8.5186,0.0296767,-1.52448e-05,3.00448e-09,-2.11801e-13,47695,-5.11579], Tmin=(756.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(404.28,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(236.962,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOOC) + radical(OCJC=O) + radical(CJCOOH) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]COO[CH][C]=O(12731)',
    structure = SMILES('[CH2]COO[CH][C]=O'),
    E0 = (217.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09686,0.0696132,-8.68927e-05,6.20227e-08,-1.85238e-11,26299.5,27.7442], Tmin=(100,'K'), Tmax=(804.311,'K')), NASAPolynomial(coeffs=[8.75772,0.0315174,-1.58517e-05,3.14418e-09,-2.24382e-13,25067,-7.5474], Tmin=(804.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CJCOOH) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][CH]OOC[C]=O(12732)',
    structure = SMILES('[CH2][CH]OOC[C]=O'),
    E0 = (265.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,348.447],'cm^-1')),
        HinderedRotor(inertia=(0.0888847,'amu*angstrom^2'), symmetry=1, barrier=(7.65836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.476667,'amu*angstrom^2'), symmetry=1, barrier=(41.0552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0889156,'amu*angstrom^2'), symmetry=1, barrier=(7.66132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0319769,'amu*angstrom^2'), symmetry=1, barrier=(41.0567,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.476242,'amu*angstrom^2'), symmetry=1, barrier=(41.0568,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16503,0.0680182,-8.28915e-05,5.7397e-08,-1.66466e-11,32064.9,28.0389], Tmin=(100,'K'), Tmax=(826.489,'K')), NASAPolynomial(coeffs=[8.82882,0.0309272,-1.55742e-05,3.09685e-09,-2.21518e-13,30798.1,-7.47396], Tmin=(826.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOOC) + radical(CJCOOH) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][CH]OOC=C[O](12649)',
    structure = SMILES('[CH2][CH]OO[CH]C=O'),
    E0 = (249.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,232.019],'cm^-1')),
        HinderedRotor(inertia=(1.23384,'amu*angstrom^2'), symmetry=1, barrier=(44.9953,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0308364,'amu*angstrom^2'), symmetry=1, barrier=(45.031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.242244,'amu*angstrom^2'), symmetry=1, barrier=(9.41041,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23548,'amu*angstrom^2'), symmetry=1, barrier=(45.012,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0309024,'amu*angstrom^2'), symmetry=1, barrier=(45.0214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00725,0.0717089,-9.39784e-05,7.15326e-08,-2.27592e-11,30123.5,27.678], Tmin=(100,'K'), Tmax=(757.935,'K')), NASAPolynomial(coeffs=[8.44494,0.0324598,-1.63081e-05,3.22072e-09,-2.28805e-13,28995.9,-6.14372], Tmin=(757.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCsJOOC) + radical(OCJC=O) + radical(CJCOOH)"""),
)

species(
    label = 'C[CH]OOC=C=O(11193)',
    structure = SMILES('C[CH]OOC=C=O'),
    E0 = (79.3517,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20382,0.0649859,-7.32637e-05,4.52664e-08,-1.15088e-11,9641.58,24.3552], Tmin=(100,'K'), Tmax=(944.391,'K')), NASAPolynomial(coeffs=[10.2208,0.0267942,-1.26028e-05,2.44444e-09,-1.72918e-13,7938.46,-18.6306], Tmin=(944.391,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.3517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CCsJOOC)"""),
)

species(
    label = 'CC1OOC=C1[O](12695)',
    structure = SMILES('CC1OO[CH]C1=O'),
    E0 = (-117.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08167,0.0214395,6.67688e-05,-1.03935e-07,4.12203e-11,-14003.9,21.8217], Tmin=(100,'K'), Tmax=(970.217,'K')), NASAPolynomial(coeffs=[15.8804,0.0165076,-5.93474e-06,1.21806e-09,-9.75441e-14,-19126.9,-56.9344], Tmin=(970.217,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-117.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclopentane) + radical(OCJC=O)"""),
)

species(
    label = 'O=C1[CH]O1(10908)',
    structure = SMILES('O=C1[CH]O1'),
    E0 = (-45.2073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42411,0.0069757,1.79186e-05,-2.83142e-08,1.12355e-11,-5411.34,9.61388], Tmin=(100,'K'), Tmax=(971.625,'K')), NASAPolynomial(coeffs=[7.22361,0.00557768,-1.91278e-06,3.80669e-10,-2.98294e-14,-6822.02,-12.0671], Tmin=(971.625,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-45.2073,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(OCJC=O)"""),
)

species(
    label = 'C=COOC=C[O](11195)',
    structure = SMILES('C=COO[CH]C=O'),
    E0 = (-7.87623,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14262,0.0596679,-5.28045e-05,2.34668e-08,-4.1925e-12,-841.914,25.5714], Tmin=(100,'K'), Tmax=(1330.34,'K')), NASAPolynomial(coeffs=[13.8682,0.0214049,-9.66147e-06,1.84657e-09,-1.29547e-13,-4227.76,-39.4544], Tmin=(1330.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.87623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(OCJC=O)"""),
)

species(
    label = 'C[CH]OO[C]=C=O(12646)',
    structure = SMILES('C[CH]OO[C]=C=O'),
    E0 = (319.096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1685,370,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.401666,'amu*angstrom^2'), symmetry=1, barrier=(9.2351,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00263214,'amu*angstrom^2'), symmetry=1, barrier=(9.23827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97333,'amu*angstrom^2'), symmetry=1, barrier=(45.3708,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97362,'amu*angstrom^2'), symmetry=1, barrier=(45.3775,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (100.073,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16478,0.0698935,-0.000112276,1.01645e-07,-3.62054e-11,38473.2,27.61], Tmin=(100,'K'), Tmax=(826.89,'K')), NASAPolynomial(coeffs=[6.49961,0.0301816,-1.50132e-05,2.89123e-09,-1.99827e-13,38066.3,5.76114], Tmin=(826.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CCsJOOC) + radical(C=CJO)"""),
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
    label = 'C[CH]OO[C]=C[O](11190)',
    structure = SMILES('C[CH]OO[C]=C[O]'),
    E0 = (279.083,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,642.558,793.969],'cm^-1')),
        HinderedRotor(inertia=(0.68345,'amu*angstrom^2'), symmetry=1, barrier=(15.7139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.683383,'amu*angstrom^2'), symmetry=1, barrier=(15.7123,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.3047,'amu*angstrom^2'), symmetry=1, barrier=(29.9976,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.50365,'amu*angstrom^2'), symmetry=1, barrier=(80.5559,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3860.13,'J/mol'), sigma=(6.39506,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=602.94 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01023,0.0645328,-6.58337e-05,3.45421e-08,-7.25624e-12,33674.6,28.7746], Tmin=(100,'K'), Tmax=(1148.4,'K')), NASAPolynomial(coeffs=[13.4131,0.0213323,-9.40662e-06,1.78516e-09,-1.25257e-13,30825.9,-32.7784], Tmin=(1148.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.083,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOOC) + radical(C=CJO)"""),
)

species(
    label = 'C[CH]OO[C]=[C]O(12648)',
    structure = SMILES('C[CH]OO[C]=[C]O'),
    E0 = (377.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1670,1700,300,440,343.753],'cm^-1')),
        HinderedRotor(inertia=(0.125005,'amu*angstrom^2'), symmetry=1, barrier=(10.4817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124999,'amu*angstrom^2'), symmetry=1, barrier=(10.4817,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.325099,'amu*angstrom^2'), symmetry=1, barrier=(27.2574,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125012,'amu*angstrom^2'), symmetry=1, barrier=(10.4818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.883618,'amu*angstrom^2'), symmetry=1, barrier=(74.098,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.810066,0.0731614,-9.54508e-05,6.59273e-08,-1.81495e-11,45498.7,31.2177], Tmin=(100,'K'), Tmax=(888.285,'K')), NASAPolynomial(coeffs=[12.0293,0.0226421,-1.01442e-05,1.90571e-09,-1.31774e-13,43505.4,-21.5799], Tmin=(888.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCsJOOC) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'CCOO[C]=[C][O](12650)',
    structure = SMILES('CCOO[C][C]=O'),
    E0 = (326.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,289.06,289.061],'cm^-1')),
        HinderedRotor(inertia=(0.195215,'amu*angstrom^2'), symmetry=1, barrier=(11.5749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195215,'amu*angstrom^2'), symmetry=1, barrier=(11.575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.437867,'amu*angstrom^2'), symmetry=1, barrier=(25.9626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.476484,'amu*angstrom^2'), symmetry=1, barrier=(28.2523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.476483,'amu*angstrom^2'), symmetry=1, barrier=(28.2523,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1221,0.0681143,-7.6391e-05,4.52704e-08,-1.10365e-11,39354.2,24.6323], Tmin=(100,'K'), Tmax=(979.546,'K')), NASAPolynomial(coeffs=[11.0948,0.0273905,-1.403e-05,2.82833e-09,-2.04461e-13,37400.4,-23.2743], Tmin=(979.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][CH]OOC=[C]O(12733)',
    structure = SMILES('[CH2][CH]OOC=[C]O'),
    E0 = (351.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,663.343],'cm^-1')),
        HinderedRotor(inertia=(0.755321,'amu*angstrom^2'), symmetry=1, barrier=(17.3663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133288,'amu*angstrom^2'), symmetry=1, barrier=(17.3669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755339,'amu*angstrom^2'), symmetry=1, barrier=(17.3667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.946183,'amu*angstrom^2'), symmetry=1, barrier=(21.7546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.34083,'amu*angstrom^2'), symmetry=1, barrier=(76.8124,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.081,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.302441,0.0769892,-9.23008e-05,5.41303e-08,-1.22279e-11,42422.5,31.8971], Tmin=(100,'K'), Tmax=(1094.83,'K')), NASAPolynomial(coeffs=[17.7938,0.0130854,-4.74978e-06,8.19764e-10,-5.49699e-14,38592.4,-54.0738], Tmin=(1094.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCsJOOC) + radical(CJCOOH) + radical(C=CJO)"""),
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
    E0 = (190.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (537.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (588.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (694.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (676.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (724.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (198.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (215.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (511.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (360.744,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (616.073,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (378.598,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (339.415,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (323.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (190.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (197.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (279.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (215.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (545.998,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (190.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (361.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (483.22,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (527.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (399.289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (447.617,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]OO[CH][C]=O(11180)'],
    products = ['OCHCO(3676)', 'CH3CHO(1381)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O][C]=CO[O](9542)', 'CHCH3(T)(21)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C][O](6861)', 'C[CH]O[O](57)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH3(17)', '[CH]OO[CH][C]=O(12729)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', 'C[C]OO[CH][C]=O(12730)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C[CH]OO[C]=[C][O](12532)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C[CH]OO[CH][C]=O(11180)'],
    products = ['CC1OOC1[C]=O(11181)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[CH]OO[CH][C]=O(11180)'],
    products = ['C=COOC[C]=O(11185)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH2][CH]OOC=C=O(11124)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.67e+12,'cm^3/(mol*s)'), n=0.1, Ea=(6.4601,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""From training reaction 2816 used for Cds-HH_Cds-OsH;HJ
Exact match found for rate rule [Cds-HH_Cds-OsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][C]=C[O](9592)', 'C[CH][O](1382)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][CH]OO[CH][C]=O(11110)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH]OO[CH][C]=O(11180)'],
    products = ['[CH2]COO[CH][C]=O(12731)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.41e+13,'s^-1','+|-',2), n=0, Ea=(188.28,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1500,'K'), comment="""From training reaction 346 used for R2H_S;C_rad_out_H/NonDeO;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeO;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH]OOC[C]=O(12732)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][CH]OOC=C[O](12649)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6Hall;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH]OO[CH][C]=O(11180)'],
    products = ['C[CH]OOC=C=O(11193)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH]OO[CH][C]=O(11180)'],
    products = ['CC1OOC=C1[O](12695)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.24579e+11,'s^-1'), n=0.1555, Ea=(7.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R5;C_rad_out_single;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH]OO[CH][C]=O(11180)'],
    products = ['O=C1[CH]O1(10908)', 'C[CH][O](1382)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.60331e+13,'s^-1'), n=-0.0568549, Ea=(88.8943,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO;Y_rad_intra;OO_intra] for rate rule [R2OO_D;Y_rad_intra;OO_intra]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[CH]OO[CH][C]=O(11180)'],
    products = ['C=COOC=C[O](11195)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', 'C[CH]OO[C]=C=O(12646)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O][C]=C[O](9592)', 'CH3CHO(1381)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(22803.5,'m^3/(mol*s)'), n=0.695, Ea=(163.776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Od_CO-CsH;YJ] for rate rule [Od_CO-CsH;O_rad/OneDe]
Euclidian distance = 4.0
family: R_Addition_MultipleBond
Ea raised from 161.8 to 163.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['HCCO(2227)', 'C[CH]O[O](57)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(33.0536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;O_rad/NonDe] for rate rule [Ct_Ct;OJ-O2s]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C[CH]OO[C]=C[O](11190)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[CH]OO[C]=[C]O(12648)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CCOO[C]=[C][O](12650)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(11191.6,'s^-1'), n=2.1395, Ea=(72.9083,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R4H_SSS;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]OOC=[C]O(12733)'],
    products = ['C[CH]OO[CH][C]=O(11180)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(30.079,'s^-1'), n=2.77074, Ea=(96.0343,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_2H;O_H_out] + [R7Hall;C_rad_out_2H;XH_out] for rate rule [R7Hall;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2848',
    isomers = [
        'C[CH]OO[CH][C]=O(11180)',
    ],
    reactants = [
        ('OCHCO(3676)', 'CH3CHO(1381)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2848',
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


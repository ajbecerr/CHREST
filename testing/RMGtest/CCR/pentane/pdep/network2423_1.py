species(
    label = 'CC([O])[CH][C]=O(2839)',
    structure = SMILES('CC([O])[CH][C]=O'),
    E0 = (190.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,180,3194.19],'cm^-1')),
        HinderedRotor(inertia=(0.00562155,'amu*angstrom^2'), symmetry=1, barrier=(3.00918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170101,'amu*angstrom^2'), symmetry=1, barrier=(14.8165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.643915,'amu*angstrom^2'), symmetry=1, barrier=(14.8049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.8653,0.0482555,-4.77339e-05,2.65264e-08,-6.1333e-12,23005.7,24.7847], Tmin=(100,'K'), Tmax=(1029.12,'K')), NASAPolynomial(coeffs=[8.60644,0.0220551,-9.54704e-06,1.78992e-09,-1.2444e-13,21618.1,-7.93116], Tmin=(1029.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(190.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCJCO) + radical(CCCJ=O)"""),
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
    label = '[O]C[CH][C]=O(10377)',
    structure = SMILES('[O]C[CH][C]=O'),
    E0 = (222.417,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1855,455,950,3135.94,3137.1],'cm^-1')),
        HinderedRotor(inertia=(0.151809,'amu*angstrom^2'), symmetry=1, barrier=(3.49039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(5.20297,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60088,0.0383015,-6.83039e-05,7.07067e-08,-2.7222e-11,26793.7,21.2914], Tmin=(100,'K'), Tmax=(860.701,'K')), NASAPolynomial(coeffs=[1.68396,0.0242382,-1.18593e-05,2.25073e-09,-1.53356e-13,27630.3,29.5206], Tmin=(860.701,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'CC([C]=O)[CH][O](10378)',
    structure = SMILES('CC([C]=O)[CH][O]'),
    E0 = (172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60873,0.0609547,-9.98061e-05,9.4935e-08,-3.48611e-11,20764.9,24.35], Tmin=(100,'K'), Tmax=(848.6,'K')), NASAPolynomial(coeffs=[3.92589,0.0315807,-1.52685e-05,2.89848e-09,-1.98235e-13,21036,17.466], Tmin=(848.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCsJOH) + radical(CC(C)CJ=O)"""),
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
    label = 'C[CH][CH][C]=O(9638)',
    structure = SMILES('C[CH]C=[C][O]'),
    E0 = (297.005,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,484.767,484.8],'cm^-1')),
        HinderedRotor(inertia=(0.000717709,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163781,'amu*angstrom^2'), symmetry=1, barrier=(27.3151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.43094,0.0256217,1.17561e-05,-3.33743e-08,1.45333e-11,35785.9,19.0015], Tmin=(100,'K'), Tmax=(974.935,'K')), NASAPolynomial(coeffs=[10.1488,0.0154088,-5.53624e-06,1.01957e-09,-7.35997e-14,33261.5,-23.2652], Tmin=(974.935,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.005,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(C=CJO)"""),
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
    label = 'CC([O])[C]=[C][O](10234)',
    structure = SMILES('CC([O])[C][C]=O'),
    E0 = (438.963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,265.096,265.101,265.105],'cm^-1')),
        HinderedRotor(inertia=(0.00239878,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169411,'amu*angstrom^2'), symmetry=1, barrier=(8.44833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.377135,'amu*angstrom^2'), symmetry=1, barrier=(18.8058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60813,0.0564984,-7.49717e-05,5.30539e-08,-1.4372e-11,52877.6,23.1254], Tmin=(100,'K'), Tmax=(694.715,'K')), NASAPolynomial(coeffs=[8.849,0.0204519,-9.32926e-06,1.75741e-09,-1.21207e-13,51735.3,-10.1505], Tmin=(694.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = 'CC(O)=C[C]=O(10379)',
    structure = SMILES('CC(O)=C[C]=O'),
    E0 = (-171.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20371,0.0626683,-7.36995e-05,4.36085e-08,-1.02142e-11,-20503,19.8684], Tmin=(100,'K'), Tmax=(1039.42,'K')), NASAPolynomial(coeffs=[12.99,0.0173109,-8.24343e-06,1.62585e-09,-1.16571e-13,-22953.2,-37.4493], Tmin=(1039.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-171.301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC(=O)C[C]=O(2831)',
    structure = SMILES('CC(=O)C[C]=O'),
    E0 = (-167.699,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.259701,'amu*angstrom^2'), symmetry=1, barrier=(5.97103,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259487,'amu*angstrom^2'), symmetry=1, barrier=(5.96611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.259773,'amu*angstrom^2'), symmetry=1, barrier=(5.97269,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6758,0.0552233,-7.55477e-05,6.22687e-08,-2.10704e-11,-20089.8,20.8778], Tmin=(100,'K'), Tmax=(801.644,'K')), NASAPolynomial(coeffs=[6.69139,0.0256778,-1.18079e-05,2.22924e-09,-1.53549e-13,-20748.7,-1.30503], Tmin=(801.644,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-167.699,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + radical(C=OCCJ=O)"""),
)

species(
    label = 'CC(=O)[CH][C]=O(2685)',
    structure = SMILES('CC([O])=C[C]=O'),
    E0 = (-33.4967,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.643647,'amu*angstrom^2'), symmetry=1, barrier=(14.7987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.647005,'amu*angstrom^2'), symmetry=1, barrier=(14.8759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49459,0.0595615,-8.42204e-05,6.47522e-08,-2.01404e-11,-3942.59,20.3273], Tmin=(100,'K'), Tmax=(783.703,'K')), NASAPolynomial(coeffs=[9.0362,0.0210688,-1.05448e-05,2.07825e-09,-1.47208e-13,-5124.65,-14.2184], Tmin=(783.703,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-33.4967,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CCJ=O)"""),
)

species(
    label = 'CC=C[C]=O(9636)',
    structure = SMILES('C[CH]C=C=O'),
    E0 = (39.4293,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,1069.02],'cm^-1')),
        HinderedRotor(inertia=(0.163668,'amu*angstrom^2'), symmetry=1, barrier=(3.76304,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30969,'amu*angstrom^2'), symmetry=1, barrier=(30.1123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.47551,0.0278456,-3.33541e-06,-1.12635e-08,5.08824e-12,4801.91,15.1824], Tmin=(100,'K'), Tmax=(1090.48,'K')), NASAPolynomial(coeffs=[7.95588,0.0195027,-8.03539e-06,1.49913e-09,-1.04907e-13,2907.47,-14.9378], Tmin=(1090.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.4293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Allyl_S)"""),
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
    label = '[O]C=C[C]=O(6950)',
    structure = SMILES('[O]C=C[C]=O'),
    E0 = (11.9527,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.24455,'amu*angstrom^2'), symmetry=1, barrier=(28.6147,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90726,0.0409354,-4.49734e-05,2.27291e-08,-4.37773e-12,1517.23,16.7461], Tmin=(100,'K'), Tmax=(1283.03,'K')), NASAPolynomial(coeffs=[13.6845,0.0042188,-2.0481e-06,4.25142e-10,-3.18151e-14,-1504.9,-43.0074], Tmin=(1283.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(11.9527,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=COJ) + radical(C=CCJ=O)"""),
)

species(
    label = '[O][CH][CH][C]=O(10226)',
    structure = SMILES('[O][C]=C[CH][O]'),
    E0 = (368.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,601.601,603.403,605.692],'cm^-1')),
        HinderedRotor(inertia=(0.000467773,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80501,0.0208365,-1.4527e-06,-1.24356e-08,5.89833e-12,44356.3,21.0677], Tmin=(100,'K'), Tmax=(1032.4,'K')), NASAPolynomial(coeffs=[8.46942,0.0109373,-4.57394e-06,8.83144e-10,-6.39732e-14,42544.7,-9.54985], Tmin=(1032.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(C=CCJO) + radical(C=CJO)"""),
)

species(
    label = 'C[C]([O])[CH][C]=O(10230)',
    structure = SMILES('C[C]([O])[CH][C]=O'),
    E0 = (367.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,1855,455,950,307.53,2770.09],'cm^-1')),
        HinderedRotor(inertia=(0.00178662,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138146,'amu*angstrom^2'), symmetry=1, barrier=(9.24554,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138025,'amu*angstrom^2'), symmetry=1, barrier=(9.24729,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55801,0.058877,-9.18218e-05,7.86374e-08,-2.63588e-11,44256.3,25.3068], Tmin=(100,'K'), Tmax=(856.342,'K')), NASAPolynomial(coeffs=[7.45478,0.0217449,-9.98484e-06,1.85204e-09,-1.24998e-13,43597.9,-0.174522], Tmin=(856.342,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([O])[CH][C]=O(2671)',
    structure = SMILES('[CH2]C([O])[CH][C]=O'),
    E0 = (402.239,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,851.202,2647.18],'cm^-1')),
        HinderedRotor(inertia=(0.229533,'amu*angstrom^2'), symmetry=1, barrier=(5.27742,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0274428,'amu*angstrom^2'), symmetry=1, barrier=(14.1083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.613626,'amu*angstrom^2'), symmetry=1, barrier=(14.1085,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3870.7,'J/mol'), sigma=(6.40089,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.59 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63814,0.0543197,-7.05236e-05,4.96507e-08,-1.40152e-11,48461.1,26.2344], Tmin=(100,'K'), Tmax=(865.41,'K')), NASAPolynomial(coeffs=[9.44214,0.0182491,-8.00344e-06,1.48877e-09,-1.02264e-13,47110.3,-10.2874], Tmin=(865.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(402.239,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCJCO) + radical(CJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'C[C](O)[CH][C]=O(10380)',
    structure = SMILES('C[C](O)[CH][C]=O'),
    E0 = (136.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39485,0.0628716,-9.77069e-05,8.35934e-08,-2.7898e-11,16555.8,26.0845], Tmin=(100,'K'), Tmax=(866.882,'K')), NASAPolynomial(coeffs=[7.54724,0.0234903,-1.05429e-05,1.93287e-09,-1.29426e-13,15902.1,-0.336042], Tmin=(866.882,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'C[C]([O])C[C]=O(2837)',
    structure = SMILES('C[C]([O])C[C]=O'),
    E0 = (167.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.178356,'amu*angstrom^2'), symmetry=1, barrier=(4.10077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.175464,'amu*angstrom^2'), symmetry=1, barrier=(4.03427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176952,'amu*angstrom^2'), symmetry=1, barrier=(4.06847,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37769,0.0639271,-9.99079e-05,8.7506e-08,-2.99911e-11,20219,23.4574], Tmin=(100,'K'), Tmax=(854.64,'K')), NASAPolynomial(coeffs=[6.83943,0.0261287,-1.20917e-05,2.25272e-09,-1.52484e-13,19732.3,0.57984], Tmin=(854.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(O)[CH][C]=O(2838)',
    structure = SMILES('[CH2]C(O)[CH][C]=O'),
    E0 = (171.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44932,0.0586408,-7.76733e-05,5.64463e-08,-1.6441e-11,20761.7,27.1025], Tmin=(100,'K'), Tmax=(840.049,'K')), NASAPolynomial(coeffs=[9.58368,0.0199075,-8.50979e-06,1.55714e-09,-1.05642e-13,19395,-10.7231], Tmin=(840.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCO) + radical(CJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([O])C[C]=O(2323)',
    structure = SMILES('[CH2]C([O])C[C]=O'),
    E0 = (202.337,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1855,455,950,246.733,246.907],'cm^-1')),
        HinderedRotor(inertia=(0.00276507,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130657,'amu*angstrom^2'), symmetry=1, barrier=(5.64502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130572,'amu*angstrom^2'), symmetry=1, barrier=(5.64507,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3870.7,'J/mol'), sigma=(6.40089,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=604.59 K, Pc=33.49 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44273,0.059557,-7.9314e-05,5.95249e-08,-1.81291e-11,24424.5,24.4385], Tmin=(100,'K'), Tmax=(800.528,'K')), NASAPolynomial(coeffs=[8.83515,0.0226197,-1.01031e-05,1.88784e-09,-1.29623e-13,23240.9,-9.58108], Tmin=(800.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'CC(=O)[CH][CH][O](2794)',
    structure = SMILES('CC([O])=C[CH][O]'),
    E0 = (83.2088,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74922,0.0426726,-3.04311e-05,1.1036e-08,-1.61336e-12,10094.4,23.8364], Tmin=(100,'K'), Tmax=(1607.37,'K')), NASAPolynomial(coeffs=[12.0126,0.0171319,-6.59658e-06,1.1505e-09,-7.58329e-14,6795.05,-30.5492], Tmin=(1607.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.2088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJO)"""),
)

species(
    label = '[CH2]C([O])C=C[O](1943)',
    structure = SMILES('[CH2]C([O])C=C[O]'),
    E0 = (185.291,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,333.305,336.401,337.542],'cm^-1')),
        HinderedRotor(inertia=(0.264343,'amu*angstrom^2'), symmetry=1, barrier=(20.9622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263358,'amu*angstrom^2'), symmetry=1, barrier=(20.9887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33903,0.0448345,-7.01928e-06,-3.46412e-08,1.96775e-11,22393.8,24.6071], Tmin=(100,'K'), Tmax=(938.528,'K')), NASAPolynomial(coeffs=[18.4331,0.00647177,-8.3256e-07,1.22255e-10,-1.32806e-14,17666.1,-64.8701], Tmin=(938.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(185.291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=COJ) + radical(CJCO)"""),
)

species(
    label = 'CC([O])C=C=O(2833)',
    structure = SMILES('CC([O])C=C=O'),
    E0 = (-52.1701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.595162,'amu*angstrom^2'), symmetry=1, barrier=(13.6839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.593936,'amu*angstrom^2'), symmetry=1, barrier=(13.6558,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49301,0.0581826,-7.45239e-05,5.47977e-08,-1.64644e-11,-6187.16,19.948], Tmin=(100,'K'), Tmax=(809.517,'K')), NASAPolynomial(coeffs=[8.50602,0.0235281,-1.03076e-05,1.91069e-09,-1.30716e-13,-7322.53,-12.4034], Tmin=(809.517,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-52.1701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CC(C)OJ)"""),
)

species(
    label = 'CC1C=C([O])O1(10381)',
    structure = SMILES('CC1[CH]C(=O)O1'),
    E0 = (-132.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89339,0.00689325,8.31154e-05,-1.15876e-07,4.65415e-11,-15938.4,19.1808], Tmin=(100,'K'), Tmax=(913.107,'K')), NASAPolynomial(coeffs=[11.6236,0.0136544,-1.92322e-06,1.89578e-10,-1.48538e-14,-19408.9,-32.4176], Tmin=(913.107,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-132.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(CCJCO)"""),
)

species(
    label = 'CC(=O)C=C[O](2791)',
    structure = SMILES('CC(=O)C=C[O]'),
    E0 = (-193.824,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85519,0.0389105,-1.35655e-05,-9.66105e-09,6.01806e-12,-23227.3,20.3807], Tmin=(100,'K'), Tmax=(1061.75,'K')), NASAPolynomial(coeffs=[11.8841,0.017903,-7.58626e-06,1.46535e-09,-1.05604e-13,-26302.5,-33.0568], Tmin=(1061.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.824,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'CC([O])C#C[O](10382)',
    structure = SMILES('CC([O])[C]=C=O'),
    E0 = (150.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.15465,'amu*angstrom^2'), symmetry=1, barrier=(26.5477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793483,'amu*angstrom^2'), symmetry=1, barrier=(18.2437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41805,0.0599894,-8.47153e-05,6.58687e-08,-2.04734e-11,18144.1,19.9228], Tmin=(100,'K'), Tmax=(829.572,'K')), NASAPolynomial(coeffs=[9.15419,0.020388,-8.9516e-06,1.64143e-09,-1.10841e-13,16939.7,-15.4772], Tmin=(829.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(150.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CC(C)OJ) + radical(CCCJ=C=O)"""),
)

species(
    label = 'CC([O])[C]=C[O](10383)',
    structure = SMILES('CC([O])[C]=C[O]'),
    E0 = (211.544,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,298.176,298.284,299.895],'cm^-1')),
        HinderedRotor(inertia=(0.299696,'amu*angstrom^2'), symmetry=1, barrier=(19.0706,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.300887,'amu*angstrom^2'), symmetry=1, barrier=(19.0638,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45889,0.0442413,-1.20638e-05,-2.44822e-08,1.4823e-11,25545,24.0167], Tmin=(100,'K'), Tmax=(950.831,'K')), NASAPolynomial(coeffs=[16.4485,0.00972059,-2.62622e-06,4.66893e-10,-3.63913e-14,21404.4,-54.3276], Tmin=(950.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'CC(O)[C]=[C][O](10384)',
    structure = SMILES('CC(O)[C][C]=O'),
    E0 = (208.602,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,247.777,247.787],'cm^-1')),
        HinderedRotor(inertia=(0.229403,'amu*angstrom^2'), symmetry=1, barrier=(9.99443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.229392,'amu*angstrom^2'), symmetry=1, barrier=(9.99441,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00274567,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.334588,'amu*angstrom^2'), symmetry=1, barrier=(14.577,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34382,0.061882,-8.68106e-05,6.77402e-08,-2.12297e-11,25181.4,24.2543], Tmin=(100,'K'), Tmax=(816.945,'K')), NASAPolynomial(coeffs=[9.0931,0.0219165,-9.71536e-06,1.79585e-09,-1.22005e-13,23982.7,-11.1515], Tmin=(816.945,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.602,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'CC([O])[C]=[C]O(10385)',
    structure = SMILES('CC([O])[C]=[C]O'),
    E0 = (309.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1670,1700,300,440,304.675,304.687],'cm^-1')),
        HinderedRotor(inertia=(0.216022,'amu*angstrom^2'), symmetry=1, barrier=(14.2342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216117,'amu*angstrom^2'), symmetry=1, barrier=(14.2341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216063,'amu*angstrom^2'), symmetry=1, barrier=(14.2342,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.829748,0.0579211,-5.92608e-05,2.9459e-08,-5.54003e-12,37387.6,27.9987], Tmin=(100,'K'), Tmax=(1447.12,'K')), NASAPolynomial(coeffs=[16.7851,0.00808353,-1.65737e-06,1.83361e-10,-9.34051e-15,33370.3,-52.7985], Tmin=(1447.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C[C]([O])C=[C]O(10386)',
    structure = SMILES('CC([O])=C[C]O'),
    E0 = (201.197,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,228.711,230.647,233.805],'cm^-1')),
        HinderedRotor(inertia=(0.181065,'amu*angstrom^2'), symmetry=1, barrier=(6.90141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180829,'amu*angstrom^2'), symmetry=1, barrier=(6.8921,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184491,'amu*angstrom^2'), symmetry=1, barrier=(6.8973,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35608,0.0580256,-6.88789e-05,4.26065e-08,-1.04035e-11,24293.8,22.5094], Tmin=(100,'K'), Tmax=(1003.28,'K')), NASAPolynomial(coeffs=[11.8156,0.0163243,-6.53157e-06,1.17756e-09,-8.01229e-14,22195,-27.9861], Tmin=(1003.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.197,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C([O])C=[C]O(10387)',
    structure = SMILES('[CH2]C([O])C=[C]O'),
    E0 = (283.573,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,355.155,355.344],'cm^-1')),
        HinderedRotor(inertia=(0.173742,'amu*angstrom^2'), symmetry=1, barrier=(15.5302,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172348,'amu*angstrom^2'), symmetry=1, barrier=(15.5464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171456,'amu*angstrom^2'), symmetry=1, barrier=(15.569,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.469509,0.0612828,-6.36418e-05,3.12487e-08,-5.6961e-12,34247.1,29.4566], Tmin=(100,'K'), Tmax=(1552.54,'K')), NASAPolynomial(coeffs=[18.4116,0.00525213,-3.49221e-08,-1.31748e-10,1.19534e-14,29857.5,-61.191], Tmin=(1552.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(CJCO) + radical(C=CJO)"""),
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
    E0 = (190.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (641.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (389.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (816.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (636.802,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (650.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (193.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (213.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (213.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (199.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (303.296,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (282.434,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (190.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (504.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (579.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (614.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (304.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (348.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (329.328,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (320.276,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (343.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (316.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (190.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (198.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (254.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (377.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (344.245,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (415.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (332.371,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (460.171,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (285.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (357.192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['HCCO(2227)', 'CH3CHO(1381)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[O]C[CH][C]=O(10377)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['CC([C]=O)[CH][O](10378)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.37952e+08,'s^-1'), n=1.37167, Ea=(199.228,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', 'C[CH][CH][C]=O(9638)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C][O](6861)', 'C[CH][O](1382)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'CC([O])[C]=[C][O](10234)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['CC1OC1[C]=O(10370)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['CC(O)=C[C]=O(10379)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['CC(=O)C[C]=O(2831)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'CC(=O)[CH][C]=O(2685)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C][O](6861)', 'CH3CHO(1381)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O(4)', 'CC=C[C]=O(9636)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH3(17)', '[O]C=C[C]=O(6950)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.61258,'m^3/(mol*s)'), n=1.485, Ea=(42.1469,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from 34.0 to 42.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH3(17)', '[O][CH][CH][C]=O(10226)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'C[C]([O])[CH][C]=O(10230)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH2]C([O])[CH][C]=O(2671)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['C[C](O)[CH][C]=O(10380)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['C[C]([O])C[C]=O(2837)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['[CH2]C(O)[CH][C]=O(2838)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C([O])C[C]=O(2323)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['CC(=O)[CH][CH][O](2794)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.10706e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C([O])C=C[O](1943)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_2;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['CC([O])C=C=O(2833)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['CC1C=C([O])O1(10381)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC([O])[CH][C]=O(2839)'],
    products = ['CC(=O)C=C[O](2791)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', 'CC([O])C#C[O](10382)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['HCCO(2227)', 'C[CH][O](1382)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CC([O])[C]=C[O](10383)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['CC(O)[C]=[C][O](10384)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(117344,'s^-1'), n=2.01217, Ea=(123.77,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS_Cs;Y_rad_out;O_H_out] + [R3H_SS_Cs;Cd_rad_out;XH_out] for rate rule [R3H_SS_Cs;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CC([O])[C]=[C]O(10385)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[C]([O])C=[C]O(10386)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C([O])C=[C]O(10387)'],
    products = ['CC([O])[CH][C]=O(2839)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2423',
    isomers = [
        'CC([O])[CH][C]=O(2839)',
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
    label = 'PDepNetwork #2423',
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


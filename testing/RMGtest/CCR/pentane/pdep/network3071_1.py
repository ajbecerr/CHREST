species(
    label = '[O][C]=CC([O])=C[O](11388)',
    structure = SMILES('[O]C=C([O])[CH][C]=O'),
    E0 = (122.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,380.143,380.433,380.493],'cm^-1')),
        HinderedRotor(inertia=(0.00116362,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20124,'amu*angstrom^2'), symmetry=1, barrier=(20.6549,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43853,0.0424997,-1.2256e-05,-2.78613e-08,1.68493e-11,14855,27.0016], Tmin=(100,'K'), Tmax=(957.826,'K')), NASAPolynomial(coeffs=[19.5475,0.000592525,5.68507e-07,-3.47189e-11,-6.39784e-15,9839.26,-67.6573], Tmin=(957.826,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJCO) + radical(CCCJ=O)"""),
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
    label = '[O]C=[C][CH][C]=O(10420)',
    structure = SMILES('[O]C=[C][CH][C]=O'),
    E0 = (399.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,280.37,281.055],'cm^-1')),
        HinderedRotor(inertia=(0.536924,'amu*angstrom^2'), symmetry=1, barrier=(29.2866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519015,'amu*angstrom^2'), symmetry=1, barrier=(29.3345,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96492,0.0301427,1.04884e-05,-4.48047e-08,2.11533e-11,48130.7,21.9171], Tmin=(100,'K'), Tmax=(973.854,'K')), NASAPolynomial(coeffs=[17.5085,0.00182392,-6.10857e-07,2.51427e-10,-2.89905e-14,43418.6,-61.309], Tmin=(973.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCJC=O) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    label = '[CH]=C([O])C=[C][O](9624)',
    structure = SMILES('[CH]C([O])=C[C]=O'),
    E0 = (337.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09785,'amu*angstrom^2'), symmetry=1, barrier=(48.2337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10192,'amu*angstrom^2'), symmetry=1, barrier=(48.3274,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3949.55,'J/mol'), sigma=(6.24613,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=616.91 K, Pc=36.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46965,0.0598019,-8.50278e-05,6.58157e-08,-2.06105e-11,40641.6,21.5237], Tmin=(100,'K'), Tmax=(778.943,'K')), NASAPolynomial(coeffs=[9.00533,0.0211058,-1.05129e-05,2.04269e-09,-1.43112e-13,39467.6,-12.9492], Tmin=(778.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
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
    label = '[O]C=C([O])[C][C]=O(13053)',
    structure = SMILES('[O][C]=[C]C([O])=C[O]'),
    E0 = (321.913,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.58092,'amu*angstrom^2'), symmetry=1, barrier=(36.3484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.327655,0.0683478,-9.11331e-05,5.37639e-08,-1.1437e-11,38860.2,27.1717], Tmin=(100,'K'), Tmax=(1365.3,'K')), NASAPolynomial(coeffs=[19.9274,-0.00403763,4.83328e-06,-1.12311e-09,8.32667e-14,34902.9,-68.3811], Tmin=(1365.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.913,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[O]C=C1OC1[C]=O(13039)',
    structure = SMILES('[O]C=C1OC1[C]=O'),
    E0 = (-6.40991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39848,0.0334868,3.60317e-05,-9.68102e-08,4.69034e-11,-655.049,21.5418], Tmin=(100,'K'), Tmax=(911.054,'K')), NASAPolynomial(coeffs=[26.2811,-0.0130714,9.47223e-06,-1.84721e-09,1.19443e-13,-7790.61,-110.463], Tmin=(911.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-6.40991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(methyleneoxirane) + radical(C=COJ) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][CH]C1=COO1(13054)',
    structure = SMILES('O=[C]C=C1[CH]OO1'),
    E0 = (245.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84765,0.0355543,-2.349e-06,-2.43565e-08,1.13444e-11,29589.7,23.107], Tmin=(100,'K'), Tmax=(1066.93,'K')), NASAPolynomial(coeffs=[14.9917,0.0125157,-6.84887e-06,1.50553e-09,-1.16558e-13,25291.4,-48.1562], Tmin=(1066.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutane) + radical(C=CCJO) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]C1=COC1[C]=O(12991)',
    structure = SMILES('[O]C1=COC1[C]=O'),
    E0 = (-35.2991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51869,0.0303671,4.40693e-05,-1.04444e-07,4.95248e-11,-4133.49,21.9723], Tmin=(100,'K'), Tmax=(909.9,'K')), NASAPolynomial(coeffs=[25.8299,-0.0125709,9.45341e-06,-1.85658e-09,1.20345e-13,-11204.3,-107.564], Tmin=(909.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.2991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][CH]C(O)=C=O(13055)',
    structure = SMILES('O=[C]C=C(O)[C]=O'),
    E0 = (-90.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279235,0.088552,-0.00015252,1.23731e-07,-3.88158e-11,-10713.7,22.5932], Tmin=(100,'K'), Tmax=(788.239,'K')), NASAPolynomial(coeffs=[15.2133,0.0127658,-8.29695e-06,1.74878e-09,-1.26646e-13,-13068,-45.9009], Tmin=(788.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-90.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CCJ=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]C(=C=O)C[C]=O(11828)',
    structure = SMILES('O=[C]CC(=O)[C]=O'),
    E0 = (-60.8595,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1850,1860,440,470,900,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.478851,'amu*angstrom^2'), symmetry=1, barrier=(11.0097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.478701,'amu*angstrom^2'), symmetry=1, barrier=(11.0063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.478388,'amu*angstrom^2'), symmetry=1, barrier=(10.9991,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17211,0.0701346,-0.000125684,1.12856e-07,-3.84676e-11,-7225.45,24.3442], Tmin=(100,'K'), Tmax=(867.39,'K')), NASAPolynomial(coeffs=[8.47541,0.019855,-1.0027e-05,1.89959e-09,-1.28288e-13,-7867.94,-6.25129], Tmin=(867.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-60.8595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-OdCsH) + group(Cds-O2d(Cds-O2d)H) + radical(C=OCCJ=O) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1O[C]1[CH][C]=O(13056)',
    structure = SMILES('[O]C1O[C]1[CH][C]=O'),
    E0 = (355.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72435,0.0573618,-9.71808e-05,8.57252e-08,-2.78242e-11,42823.8,26.9582], Tmin=(100,'K'), Tmax=(964.757,'K')), NASAPolynomial(coeffs=[5.55448,0.0205231,-7.31804e-06,1.11047e-09,-6.26835e-14,43060.1,13.6724], Tmin=(964.757,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(355.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1C([O])C1[C]=O(13057)',
    structure = SMILES('[O][C]1C([O])C1[C]=O'),
    E0 = (380.83,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45851,0.0513437,-5.55057e-05,2.93197e-08,-5.99708e-12,45898.6,24.395], Tmin=(100,'K'), Tmax=(1202.51,'K')), NASAPolynomial(coeffs=[13.8957,0.00997306,-3.90051e-06,7.10052e-10,-4.91854e-14,42907.4,-37.9008], Tmin=(1202.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[O]C1([CH][C]=O)[CH]O1(13058)',
    structure = SMILES('[O]C1([CH][C]=O)[CH]O1'),
    E0 = (353.398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.608553,0.0617628,-7.61508e-05,4.42377e-08,-9.27695e-12,42637.2,27.8066], Tmin=(100,'K'), Tmax=(1418.38,'K')), NASAPolynomial(coeffs=[15.906,0.00272957,3.08658e-06,-9.05173e-10,7.21962e-14,39896.3,-45.7062], Tmin=(1418.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CCsJO) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1[CH]C(=O)C1[O](13059)',
    structure = SMILES('[O][C]1C=C([O])C1[O]'),
    E0 = (359.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53846,0.0412278,-1.02956e-05,-2.91749e-08,1.75553e-11,43332.3,24.5751], Tmin=(100,'K'), Tmax=(939.002,'K')), NASAPolynomial(coeffs=[18.4414,0.00175576,7.90933e-07,-1.50351e-10,4.71258e-15,38723.7,-63.5449], Tmin=(939.002,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(359.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[O][CH]C1([O])[CH]C1=O(13060)',
    structure = SMILES('[O][CH]C1([O])C=C1[O]'),
    E0 = (443.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38221,0.0612549,-8.63292e-05,6.27147e-08,-1.80772e-11,53478.2,23.2026], Tmin=(100,'K'), Tmax=(849.731,'K')), NASAPolynomial(coeffs=[10.7465,0.0171755,-8.52014e-06,1.67094e-09,-1.18132e-13,51886.7,-20.45], Tmin=(849.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CC(C)2OJ) + radical(C=C(C)OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[O]C(=C=O)[CH][C]=O(13061)',
    structure = SMILES('[O][C]=C([O])C=C=O'),
    E0 = (116.888,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,1685,370,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.671234,'amu*angstrom^2'), symmetry=1, barrier=(15.433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32525,0.0583967,-8.38881e-05,5.80109e-08,-1.53044e-11,14155.1,24.5767], Tmin=(100,'K'), Tmax=(947.579,'K')), NASAPolynomial(coeffs=[13.4216,0.00709588,-2.3018e-06,3.4515e-10,-2.02062e-14,11873.3,-33.0732], Tmin=(947.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C=C=C[C]=O(10418)',
    structure = SMILES('[O]C=C=C[C]=O'),
    E0 = (146.581,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.34698,'amu*angstrom^2'), symmetry=1, barrier=(30.9696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7611,0.0442043,-4.36984e-05,2.00881e-08,-3.57078e-12,17714.4,19.2867], Tmin=(100,'K'), Tmax=(1370.85,'K')), NASAPolynomial(coeffs=[14.1522,0.00804831,-4.13608e-06,8.48302e-10,-6.20339e-14,14317.1,-44.4016], Tmin=(1370.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(146.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(C=COJ) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]=C([O])[CH][C]=O(13062)',
    structure = SMILES('[O][C]=C([O])[CH][C]=O'),
    E0 = (362.384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,350,440,435,1725,1855,455,950,1685,370,423.162,423.252,423.387],'cm^-1')),
        HinderedRotor(inertia=(0.000941034,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000940344,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63241,0.0444181,-3.95473e-05,1.11903e-08,6.48475e-13,43676.7,29.4359], Tmin=(100,'K'), Tmax=(1029.67,'K')), NASAPolynomial(coeffs=[15.6575,0.00430427,-2.04425e-06,4.62487e-10,-3.76504e-14,40026.7,-42.3361], Tmin=(1029.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJCO) + radical(C=CJO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C(=[C]O)[CH][C]=O(13063)',
    structure = SMILES('[O]C(=[C]O)[CH][C]=O'),
    E0 = (220.922,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,1685,370,1855,455,950,404.971,405],'cm^-1')),
        HinderedRotor(inertia=(0.00102779,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125663,'amu*angstrom^2'), symmetry=1, barrier=(14.6266,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125664,'amu*angstrom^2'), symmetry=1, barrier=(14.6266,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.667192,0.0577597,-6.45584e-05,3.21494e-08,-5.87952e-12,26704.1,31.5007], Tmin=(100,'K'), Tmax=(1536.71,'K')), NASAPolynomial(coeffs=[19.7871,-0.000980538,1.53597e-06,-3.23356e-10,2.13843e-14,21887,-65.5101], Tmin=(1536.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(220.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(CCJCO) + radical(C=CJO) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]=C(O)[CH][C]=O(13064)',
    structure = SMILES('O=[C][CH][C](O)[C]=O'),
    E0 = (216.667,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,3025,407.5,1350,352.5,1850,1860,440,470,900,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0210726,'amu*angstrom^2'), symmetry=1, barrier=(4.03332,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000959935,'amu*angstrom^2'), symmetry=1, barrier=(10.8992,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0481415,'amu*angstrom^2'), symmetry=1, barrier=(10.9014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0477034,'amu*angstrom^2'), symmetry=1, barrier=(10.9002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32224,0.0630165,-0.000103163,8.42965e-08,-2.65188e-11,26151.5,31.6515], Tmin=(100,'K'), Tmax=(870.428,'K')), NASAPolynomial(coeffs=[10.6247,0.0135164,-6.22589e-06,1.14047e-09,-7.57871e-14,24787.8,-10.4673], Tmin=(870.428,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.667,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(CCJCO) + radical(CCCJ=O) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]=C([O])C[C]=O(11829)',
    structure = SMILES('[O][C]=C([O])C[C]=O'),
    E0 = (162.482,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1855,455,950,1685,370,377.747,378.659,378.953],'cm^-1')),
        HinderedRotor(inertia=(0.0011805,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14693,'amu*angstrom^2'), symmetry=1, barrier=(14.8814,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32315,0.0508629,-5.18798e-05,2.45824e-08,-4.4553e-12,19645.2,28.057], Tmin=(100,'K'), Tmax=(1359.44,'K')), NASAPolynomial(coeffs=[16.3931,0.00652122,-2.9534e-06,5.88978e-10,-4.29284e-14,15547.9,-49.2744], Tmin=(1359.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.482,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCCJ=O) + radical(C=CJO)"""),
)

species(
    label = '[O][C]=C([O])[CH]C=O(13065)',
    structure = SMILES('[O][C]=C([O])C=C[O]'),
    E0 = (122.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23911,'amu*angstrom^2'), symmetry=1, barrier=(28.4896,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0555981,0.0669788,-7.897e-05,4.12954e-08,-7.74135e-12,14948.7,28.9685], Tmin=(100,'K'), Tmax=(1580.46,'K')), NASAPolynomial(coeffs=[20.764,-0.00339061,4.59397e-06,-1.03019e-09,7.32037e-14,10575.5,-74.0193], Tmin=(1580.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C=C([O])C=C=O(11380)',
    structure = SMILES('[O]C=C([O])C=C=O'),
    E0 = (-122.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10808,'amu*angstrom^2'), symmetry=1, barrier=(25.4769,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.729912,0.0612206,-7.32087e-05,4.05217e-08,-8.30047e-12,-14649.2,23.5818], Tmin=(100,'K'), Tmax=(1359.99,'K')), NASAPolynomial(coeffs=[18.0567,0.00205243,1.10232e-06,-3.4264e-10,2.70638e-14,-18603.2,-62.5473], Tmin=(1359.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C=C1C=C([O])O1(13066)',
    structure = SMILES('[O]C=C1[CH]C(=O)O1'),
    E0 = (-195.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48712,0.0319209,4.37687e-05,-9.98431e-08,4.60738e-11,-23365.5,19.4549], Tmin=(100,'K'), Tmax=(922.763,'K')), NASAPolynomial(coeffs=[23.924,-0.00487041,5.28076e-06,-1.02247e-09,6.12419e-14,-30080.7,-100.936], Tmin=(922.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-195.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + ring(4-Methylene-2-oxetanone) + radical(C=COJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]C1=COC([O])=C1(13067)',
    structure = SMILES('[O]C1=COC(=O)[CH]1'),
    E0 = (-240.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2368,0.0114053,9.25555e-05,-1.40772e-07,5.7231e-11,-28784,15.8486], Tmin=(100,'K'), Tmax=(950.187,'K')), NASAPolynomial(coeffs=[21.6058,-0.000645933,1.8862e-06,-1.93771e-10,-5.37519e-15,-35601.7,-93.1124], Tmin=(950.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-240.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=C(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[O]C=CC([O])=C=O(13068)',
    structure = SMILES('[O]C=CC(=O)[C]=O'),
    E0 = (-89.9293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33453,0.0531676,-5.74577e-05,2.99649e-08,-6.02708e-12,-10715.4,24.031], Tmin=(100,'K'), Tmax=(1224.77,'K')), NASAPolynomial(coeffs=[14.7926,0.0092147,-3.62778e-06,6.64225e-10,-4.6251e-14,-14012,-43.6251], Tmin=(1224.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.9293,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)H) + radical(C=COJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]=[C]C([O])[CH][O](13069)',
    structure = SMILES('[O][CH]C([O])[C][C]=O'),
    E0 = (685.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,1855,455,950,180,180,180,748.137,748.18],'cm^-1')),
        HinderedRotor(inertia=(0.18605,'amu*angstrom^2'), symmetry=1, barrier=(4.27766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186107,'amu*angstrom^2'), symmetry=1, barrier=(4.27896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18606,'amu*angstrom^2'), symmetry=1, barrier=(4.27788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.844853,0.0840088,-0.000171313,1.65706e-07,-5.87348e-11,82502.8,28.9637], Tmin=(100,'K'), Tmax=(886.09,'K')), NASAPolynomial(coeffs=[6.04544,0.0256566,-1.34943e-05,2.54857e-09,-1.69729e-13,82950.3,12.2286], Tmin=(886.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(685.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1C=[C]OC1[O](13070)',
    structure = SMILES('[O]C1=C[C]OC1[O]'),
    E0 = (323.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6827,0.0442748,-3.79367e-05,1.47256e-08,-1.98328e-12,38988.9,21.0805], Tmin=(100,'K'), Tmax=(1194.48,'K')), NASAPolynomial(coeffs=[13.4325,0.0115026,-5.03876e-06,9.72692e-10,-6.93449e-14,35712.9,-39.657], Tmin=(1194.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(323.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(25dihydrofuran) + radical(CCOJ) + radical(C=C(C)OJ) + radical(CH2_triplet)"""),
)

species(
    label = '[O][CH]C1([O])C=[C]O1(13071)',
    structure = SMILES('[O][CH]C1([O])C=[C]O1'),
    E0 = (468.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.284059,0.0657781,-7.88073e-05,4.24134e-08,-8.30288e-12,56495.5,26.4861], Tmin=(100,'K'), Tmax=(1469.25,'K')), NASAPolynomial(coeffs=[20.3173,-0.00163007,3.14913e-06,-7.35092e-10,5.33316e-14,51997.6,-73.1438], Tmin=(1469.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.504,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=CC(C)(O)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(C=CJO)"""),
)

species(
    label = '[O]C#CC([O])=C[O](13072)',
    structure = SMILES('[O]C=C([O])[C]=C=O'),
    E0 = (114.986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03298,'amu*angstrom^2'), symmetry=1, barrier=(23.7503,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.965147,0.0626763,-8.70884e-05,5.5767e-08,-1.33389e-11,13942.4,23.2105], Tmin=(100,'K'), Tmax=(1117.37,'K')), NASAPolynomial(coeffs=[16.5945,0.00207865,4.98652e-07,-2.13093e-10,1.88522e-14,10739.8,-52.6284], Tmin=(1117.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]C=[C]C([O])=C[O](13073)',
    structure = SMILES('[O]C=[C]C([O])=C[O]'),
    E0 = (82.1685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.567144,0.0744767,-9.10438e-05,4.87732e-08,-9.27396e-12,10069.4,27.2676], Tmin=(100,'K'), Tmax=(1580.19,'K')), NASAPolynomial(coeffs=[22.2694,-0.00588608,6.65205e-06,-1.47698e-09,1.05199e-13,5668.31,-84.4436], Tmin=(1580.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(82.1685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[O][C]=[C]C(O)=C[O](13074)',
    structure = SMILES('[O][C]=[C]C(O)=C[O]'),
    E0 = (184.108,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43363,'amu*angstrom^2'), symmetry=1, barrier=(32.9619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.42933,'amu*angstrom^2'), symmetry=1, barrier=(32.8631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.391654,0.0764243,-9.76974e-05,5.45817e-08,-1.08975e-11,22318.6,28.2561], Tmin=(100,'K'), Tmax=(1478.74,'K')), NASAPolynomial(coeffs=[22.9493,-0.00650112,6.49274e-06,-1.44019e-09,1.0362e-13,17579.1,-86.1658], Tmin=(1478.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[O]C=C([O])[C]=[C]O(13075)',
    structure = SMILES('[O]C=C([O])[C]=[C]O'),
    E0 = (180.45,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29866,'amu*angstrom^2'), symmetry=1, barrier=(29.8588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30121,'amu*angstrom^2'), symmetry=1, barrier=(29.9173,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.319757,0.0781656,-0.000104826,6.14076e-08,-1.28433e-11,21873.3,28.0808], Tmin=(100,'K'), Tmax=(1411.09,'K')), NASAPolynomial(coeffs=[22.5437,-0.00685589,7.03694e-06,-1.5919e-09,1.16518e-13,17433,-82.9655], Tmin=(1411.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(180.45,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[O][C]=[C]C([O])=CO(13076)',
    structure = SMILES('[O]C(=[C][C]=O)[CH]O'),
    E0 = (174.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,1685,370,1855,455,950,244.723,244.73],'cm^-1')),
        HinderedRotor(inertia=(0.370054,'amu*angstrom^2'), symmetry=1, barrier=(15.7273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.370024,'amu*angstrom^2'), symmetry=1, barrier=(15.7272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.37003,'amu*angstrom^2'), symmetry=1, barrier=(15.7272,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.997538,0.0681396,-9.53208e-05,6.28191e-08,-1.59824e-11,21108.8,27.1689], Tmin=(100,'K'), Tmax=(969.301,'K')), NASAPolynomial(coeffs=[15.0282,0.0102393,-5.71943e-06,1.19289e-09,-8.7829e-14,18388.8,-40.0836], Tmin=(969.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(174.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[O][C]=C([O])C=[C]O(13077)',
    structure = SMILES('[O][C]=C([O])C=[C]O'),
    E0 = (221.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.903286,'amu*angstrom^2'), symmetry=1, barrier=(20.7683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.897728,'amu*angstrom^2'), symmetry=1, barrier=(20.6405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.182356,0.07077,-9.30694e-05,5.42866e-08,-1.14399e-11,26753,29.8161], Tmin=(100,'K'), Tmax=(1372.69,'K')), NASAPolynomial(coeffs=[20.7263,-0.00392795,4.76531e-06,-1.10033e-09,8.11374e-14,22510.4,-70.714], Tmin=(1372.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH][O](1548)',
    structure = SMILES('[CH][O]'),
    E0 = (431.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([806.798,806.798,2696.84],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86882,-0.000695042,1.53508e-05,-2.00576e-08,7.65929e-12,51865.1,6.76618], Tmin=(100,'K'), Tmax=(973.285,'K')), NASAPolynomial(coeffs=[6.03251,-0.000301449,4.32946e-07,-3.67045e-11,-1.25373e-15,51004.1,-5.87326], Tmin=(973.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(H3COJ) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = '[O][C]=C[C]=O(9898)',
    structure = SMILES('O=[C][CH][C]=O'),
    E0 = (211.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1850,1860,440,470,900,1000],'cm^-1')),
        HinderedRotor(inertia=(1.2407,'amu*angstrom^2'), symmetry=1, barrier=(28.5261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24145,'amu*angstrom^2'), symmetry=1, barrier=(28.5433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0388,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81165,0.0557507,-0.00011026,1.00914e-07,-3.40854e-11,25521.6,14.8862], Tmin=(100,'K'), Tmax=(897.886,'K')), NASAPolynomial(coeffs=[7.74751,0.0112957,-5.90386e-06,1.09028e-09,-7.08098e-14,25181.6,-9.06853], Tmin=(897.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCJC=O) + radical(C=OCCJ=O) + radical(C=OCCJ=O)"""),
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
    label = '[CH]=C([O])[CH][O](10249)',
    structure = SMILES('[CH]C([O])=C[O]'),
    E0 = (232.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.14782,'amu*angstrom^2'), symmetry=1, barrier=(49.3827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06413,0.0315518,3.76054e-06,-3.87898e-08,2.07567e-11,28070.3,18.0523], Tmin=(100,'K'), Tmax=(902.966,'K')), NASAPolynomial(coeffs=[15.0822,0.00356738,9.38086e-07,-3.00037e-10,2.0679e-14,24509.2,-50.1243], Tmin=(902.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[O]C1OC1=C[C]=O(13078)',
    structure = SMILES('[O]C1OC1=C[C]=O'),
    E0 = (80.2008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54618,0.0571474,-6.91821e-05,4.08749e-08,-9.5972e-12,9731.72,19.1684], Tmin=(100,'K'), Tmax=(1031.13,'K')), NASAPolynomial(coeffs=[12.2402,0.0156631,-8.83454e-06,1.85805e-09,-1.37507e-13,7526.33,-32.7517], Tmin=(1031.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.2008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(methyleneoxirane) + radical(CCOJ) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]C1=CC(=O)C1[O](13079)',
    structure = SMILES('[O]C1=CC(=O)C1[O]'),
    E0 = (28.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10721,0.0332609,-9.82269e-06,-1.33629e-08,7.74431e-12,3541.96,22.3041], Tmin=(100,'K'), Tmax=(1020.72,'K')), NASAPolynomial(coeffs=[12.3235,0.0115464,-4.83633e-06,9.65591e-10,-7.21843e-14,501.947,-31.8684], Tmin=(1020.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(28.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(C=OCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O][CH]C(O)=C=C=O(13080)',
    structure = SMILES('[O]C=C(O)[C]=C=O'),
    E0 = (-22.8193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.354102,0.0695599,-8.98514e-05,5.21488e-08,-1.11123e-11,-2604.04,23.9005], Tmin=(100,'K'), Tmax=(1310.29,'K')), NASAPolynomial(coeffs=[20.4381,-0.00153039,2.72615e-06,-6.49673e-10,4.82589e-14,-7027.81,-75.2179], Tmin=(1310.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-22.8193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = '[O]CC([O])=C=C=O(13081)',
    structure = SMILES('[O]CC([O])=C=C=O'),
    E0 = (119.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.54939,0.0351237,-5.12819e-05,4.48313e-08,-1.57911e-11,14395.7,7.87742], Tmin=(100,'K'), Tmax=(810.12,'K')), NASAPolynomial(coeffs=[5.33274,0.0164747,-7.66791e-06,1.46497e-09,-1.01504e-13,14105.8,-3.97095], Tmin=(810.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cdd-CdsCds) + missing(Cdd-CddO2d) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O][CH][C]1OC1[C]=O(13082)',
    structure = SMILES('[O][CH][C]1OC1[C]=O'),
    E0 = (366.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42574,0.0600678,-9.29892e-05,7.32494e-08,-2.14494e-11,44130.5,26.7345], Tmin=(100,'K'), Tmax=(1019.09,'K')), NASAPolynomial(coeffs=[9.12848,0.0147093,-3.96372e-06,4.47346e-10,-1.713e-14,43345.9,-6.71905], Tmin=(1019.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(CCsJOH) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]C1O[CH]C1=O(12998)',
    structure = SMILES('[O][C]C1OC=C1[O]'),
    E0 = (412.615,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09598,0.0466245,-7.95036e-06,-4.43148e-08,2.60865e-11,49746.8,22.0556], Tmin=(100,'K'), Tmax=(922.975,'K')), NASAPolynomial(coeffs=[23.676,-0.00630479,5.05275e-06,-9.67019e-10,5.98157e-14,43664.9,-95.4369], Tmin=(922.975,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.615,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(CCOJ) + radical(CH2_triplet)"""),
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
    label = 'O=[C]C=C=O(9594)',
    structure = SMILES('O=[C]C=C=O'),
    E0 = (-7.87925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1855,455,950,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.841576,'amu*angstrom^2'), symmetry=1, barrier=(19.3495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (69.0388,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.98472,0.0259268,-2.97918e-05,2.05261e-08,-6.41245e-12,-914.289,13.446], Tmin=(100,'K'), Tmax=(741.746,'K')), NASAPolynomial(coeffs=[4.88246,0.0156926,-9.09529e-06,1.92406e-09,-1.42625e-13,-1195.81,4.85747], Tmin=(741.746,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.87925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cds-(Cdd-O2d)CsH) + group(Cds-O2d(Cds-Cds)H) + missing(Cdd-CdO2d) + radical(CsCJ=O)"""),
)

species(
    label = '[O]CC([O])=[C][C]=O(13083)',
    structure = SMILES('[O]CC([O])=[C][C]=O'),
    E0 = (283.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.268807,'amu*angstrom^2'), symmetry=1, barrier=(6.1804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267264,'amu*angstrom^2'), symmetry=1, barrier=(6.14493,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.986716,0.0785859,-0.000156066,1.51148e-07,-5.46759e-11,34137.8,27.4589], Tmin=(100,'K'), Tmax=(852.087,'K')), NASAPolynomial(coeffs=[6.59474,0.0251361,-1.42253e-05,2.81538e-09,-1.94768e-13,34166.7,7.07907], Tmin=(852.087,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]C1C(=O)C1[C]=O(13084)',
    structure = SMILES('[O]C1C(=O)C1[C]=O'),
    E0 = (122.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.30678,0.0205395,4.17146e-05,-7.83398e-08,3.40285e-11,14785.3,24.4269], Tmin=(100,'K'), Tmax=(933.917,'K')), NASAPolynomial(coeffs=[16.2623,0.00342724,6.82391e-07,-1.39406e-10,2.34352e-15,10318.3,-51.9061], Tmin=(933.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(122.297,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(cyclopropanone) + radical(C=OCOJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[O][CH][C]1[CH]C(=O)O1(13085)',
    structure = SMILES('[O]C=C1[CH][C]([O])O1'),
    E0 = (225.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35194,0.0391456,1.75186e-05,-7.24732e-08,3.71431e-11,27182.9,18.9853], Tmin=(100,'K'), Tmax=(903.636,'K')), NASAPolynomial(coeffs=[22.9471,-0.00520446,6.07816e-06,-1.27922e-09,8.50956e-14,21188,-94.5866], Tmin=(903.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(CCOJ) + radical(C=COJ) + radical(C=CCJCO) + radical(Cs_P)"""),
)

species(
    label = 'O=[C][CH][C]1[CH]OO1(13086)',
    structure = SMILES('O=[C][CH][C]1[CH]OO1'),
    E0 = (575.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79043,0.0592736,-0.00011029,1.10812e-07,-4.1536e-11,69252.2,23.5212], Tmin=(100,'K'), Tmax=(859.221,'K')), NASAPolynomial(coeffs=[2.63138,0.0301172,-1.5324e-05,2.93752e-09,-2.00731e-13,70039.4,25.0137], Tmin=(859.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxetane) + radical(C2CsJOO) + radical(CCsJOO) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]=CC([O])[C]=O(11384)',
    structure = SMILES('[O]C([C]=O)[CH][C]=O'),
    E0 = (283.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,1850,1860,440,470,900,1000,295.466,4000],'cm^-1')),
        HinderedRotor(inertia=(0.303233,'amu*angstrom^2'), symmetry=1, barrier=(18.7837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.303296,'amu*angstrom^2'), symmetry=1, barrier=(18.783,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30344e-05,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76831,0.0506918,-6.72649e-05,4.6272e-08,-1.25527e-11,34208.8,31.4107], Tmin=(100,'K'), Tmax=(905.05,'K')), NASAPolynomial(coeffs=[10.1401,0.013692,-5.94308e-06,1.1023e-09,-7.5698e-14,32693.5,-8.143], Tmin=(905.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(CCCJ=O) + radical(CCCJ=O)"""),
)

species(
    label = '[O][CH][C]1C=[C]OO1(13087)',
    structure = SMILES('[O]C=C1[CH][C]OO1'),
    E0 = (490.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58857,0.0310547,3.92762e-05,-8.91549e-08,4.02241e-11,59080.9,20.0845], Tmin=(100,'K'), Tmax=(941.927,'K')), NASAPolynomial(coeffs=[22.3741,-0.0015438,2.53599e-06,-4.05746e-10,1.5394e-14,52695.6,-92.0593], Tmin=(941.927,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.337,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(C=CCJCO) + radical(CH2_triplet)"""),
)

species(
    label = 'HCO(1372)',
    structure = SMILES('[CH]=O'),
    E0 = (32.5964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1116.11,1837.47,2716.03],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06159,-0.0018195,9.37895e-06,-8.24742e-09,2.33851e-12,3918.59,3.98389], Tmin=(100,'K'), Tmax=(1105.12,'K')), NASAPolynomial(coeffs=[3.05335,0.00410464,-1.74962e-06,3.28532e-10,-2.28985e-14,4002.52,8.32038], Tmin=(1105.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.5964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[O][C][CH][C]=O(13088)',
    structure = SMILES('[O][C][CH][C]=O'),
    E0 = (683.406,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,180,180,2717.48],'cm^-1')),
        HinderedRotor(inertia=(0.110545,'amu*angstrom^2'), symmetry=1, barrier=(2.54164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116183,'amu*angstrom^2'), symmetry=1, barrier=(2.67129,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (69.0388,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38861,0.0442293,-9.33008e-05,9.35554e-08,-3.41631e-11,82244.3,20.0926], Tmin=(100,'K'), Tmax=(876.273,'K')), NASAPolynomial(coeffs=[3.99478,0.015819,-8.58615e-06,1.65359e-09,-1.11825e-13,82772.1,17.1736], Tmin=(876.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.406,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(124.717,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCO) + radical(CH2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1[CH]OC(=O)[CH]1(13089)',
    structure = SMILES('[O][C]1[CH]C([O])=CO1'),
    E0 = (142.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72442,0.0335633,2.03296e-05,-6.83507e-08,3.43662e-11,17219.4,17.975], Tmin=(100,'K'), Tmax=(896.403,'K')), NASAPolynomial(coeffs=[19.6813,-0.00183145,4.70193e-06,-1.05708e-09,7.22555e-14,12202.8,-76.7177], Tmin=(896.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(142.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(Cs_P)"""),
)

species(
    label = '[O]C1=CC1([O])C=O(13090)',
    structure = SMILES('[O]C1=CC1([O])C=O'),
    E0 = (141.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57214,0.0564796,-7.90219e-05,5.56725e-08,-1.46667e-11,17054.8,22.6779], Tmin=(100,'K'), Tmax=(726.631,'K')), NASAPolynomial(coeffs=[10.2824,0.0149211,-6.42375e-06,1.16869e-09,-7.86677e-14,15620.2,-17.7236], Tmin=(726.631,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(CC(C)(C=O)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[O]C=C=C([O])C=O(13091)',
    structure = SMILES('[O]C=[C]C(=O)C=O'),
    E0 = (-5.92874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39104,0.0536392,-5.70604e-05,2.90057e-08,-5.73004e-12,-615.975,23.7152], Tmin=(100,'K'), Tmax=(1236.91,'K')), NASAPolynomial(coeffs=[14.6307,0.0108237,-5.13775e-06,1.02043e-09,-7.37188e-14,-3891.22,-42.9735], Tmin=(1236.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-5.92874,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)H) + radical(C=COJ) + radical(C=CJC=O)"""),
)

species(
    label = '[O][C]1[CH]OO[C]=C1(13092)',
    structure = SMILES('[O]C1=COO[C][CH]1'),
    E0 = (463.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28444,0.0450157,-6.42925e-06,-4.31838e-08,2.57746e-11,55842.2,16.2798], Tmin=(100,'K'), Tmax=(898.067,'K')), NASAPolynomial(coeffs=[20.8561,-0.00174107,4.16208e-06,-9.35206e-10,6.39802e-14,50697.1,-85.1119], Tmin=(898.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(34dihydro12dioxin) + radical(C=C(C)OJ) + radical(C=CCJCO) + radical(CH2_triplet)"""),
)

species(
    label = '[O][C]=[C]C([O])C=O(13093)',
    structure = SMILES('[O]C([C][C]=O)C=O'),
    E0 = (372.124,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1855,455,950,375.745,376.711,376.888],'cm^-1')),
        HinderedRotor(inertia=(0.0677721,'amu*angstrom^2'), symmetry=1, barrier=(6.80471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00119391,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0676661,'amu*angstrom^2'), symmetry=1, barrier=(6.80429,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62391,0.0556227,-7.89749e-05,5.93262e-08,-1.77749e-11,44838.7,28.5211], Tmin=(100,'K'), Tmax=(817.307,'K')), NASAPolynomial(coeffs=[9.50788,0.0170395,-8.16679e-06,1.57186e-09,-1.09701e-13,43549.9,-7.92438], Tmin=(817.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.124,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C][CH][O](10223)',
    structure = SMILES('[O][C][CH][O]'),
    E0 = (676.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,180,180,1781.81,1786.65],'cm^-1')),
        HinderedRotor(inertia=(0.371605,'amu*angstrom^2'), symmetry=1, barrier=(8.54393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61268,0.044967,-0.000114903,1.26746e-07,-4.85652e-11,81376.1,15.5028], Tmin=(100,'K'), Tmax=(883.837,'K')), NASAPolynomial(coeffs=[-0.0833967,0.0206941,-1.18058e-05,2.28838e-09,-1.54199e-13,83277.3,36.2363], Tmin=(883.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(676.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCOJ) + radical(CCsJOH) + radical(CH2_triplet)"""),
)

species(
    label = '[O][CH][C]([O])[C]=C[O](13094)',
    structure = SMILES('[O][CH][C]=C([O])[CH][O]'),
    E0 = (510.915,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3050,390,425,1340,1360,335,370,1685,370,740.673,742.946,744.218,745.723,746.117],'cm^-1')),
        HinderedRotor(inertia=(0.000304876,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000302941,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16627,0.0430697,-4.14037e-05,2.16578e-08,-4.78475e-12,61512.7,28.0195], Tmin=(100,'K'), Tmax=(1058.13,'K')), NASAPolynomial(coeffs=[7.98473,0.021075,-1.0225e-05,2.01448e-09,-1.4382e-13,60281.3,-0.380107], Tmin=(1058.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.915,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(CCOJ) + radical(C=CCJO) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[O][C]C([C]=O)C=O(13095)',
    structure = SMILES('[O][C]C([C]=O)C=O'),
    E0 = (370.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1855,455,950,322.368,322.69,322.79],'cm^-1')),
        HinderedRotor(inertia=(0.00481096,'amu*angstrom^2'), symmetry=1, barrier=(9.46765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0048124,'amu*angstrom^2'), symmetry=1, barrier=(9.46918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.375211,'amu*angstrom^2'), symmetry=1, barrier=(27.7294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23015,0.0463077,-4.41443e-05,-5.25325e-09,2.94591e-11,44579.9,25.0453], Tmin=(100,'K'), Tmax=(492.795,'K')), NASAPolynomial(coeffs=[6.54813,0.0235676,-1.2393e-05,2.47831e-09,-1.76474e-13,44004.9,5.7527], Tmin=(492.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CH2_triplet) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH][C][O](10218)',
    structure = SMILES('[CH][C][O]'),
    E0 = (885.265,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([418.405,418.464,1386.35,1386.4,1890.01],'cm^-1')),
        HinderedRotor(inertia=(0.0607177,'amu*angstrom^2'), symmetry=1, barrier=(7.54278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.35351,0.0168371,-3.18729e-05,3.43969e-08,-1.43789e-11,106493,11.3039], Tmin=(100,'K'), Tmax=(740.605,'K')), NASAPolynomial(coeffs=[3.99585,0.00842558,-4.82643e-06,1.03992e-09,-7.72054e-14,106534,9.313], Tmin=(740.605,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(885.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCJ2_triplet) + radical(CH2_triplet)"""),
)

species(
    label = '[O][C][CH]C([O])[C]=O(13096)',
    structure = SMILES('[O][C][CH]C([O])[C]=O'),
    E0 = (733.477,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,1855,455,950,299.739,299.775,1246.65,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00187659,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364485,'amu*angstrom^2'), symmetry=1, barrier=(23.2336,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.364535,'amu*angstrom^2'), symmetry=1, barrier=(23.232,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61838,0.0576514,-8.42937e-05,6.0511e-08,-1.53235e-11,88297.7,31.4121], Tmin=(100,'K'), Tmax=(640.97,'K')), NASAPolynomial(coeffs=[9.40718,0.0171927,-8.67927e-06,1.6966e-09,-1.18899e-13,87131.9,-4.00572], Tmin=(640.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.477,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCOJ) + radical(CCJCO) + radical(CCCJ=O) + radical(CH2_triplet)"""),
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
    label = '[O][C][C]=C([O])C[O](13097)',
    structure = SMILES('[O]C[C]([O])[C][C]=O'),
    E0 = (681.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,360,370,350,1855,455,950,180,180,180,180,1698.15],'cm^-1')),
        HinderedRotor(inertia=(0.154743,'amu*angstrom^2'), symmetry=1, barrier=(3.55785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154712,'amu*angstrom^2'), symmetry=1, barrier=(3.55713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154744,'amu*angstrom^2'), symmetry=1, barrier=(3.55786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.927785,0.084397,-0.000178615,1.7768e-07,-6.41032e-11,82056.4,28.5406], Tmin=(100,'K'), Tmax=(887.083,'K')), NASAPolynomial(coeffs=[4.22701,0.0286249,-1.5157e-05,2.86779e-09,-1.91052e-13,83080.1,22.0885], Tmin=(887.083,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(681.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C][C]=C([O])[CH]O(13098)',
    structure = SMILES('[O][C][C]=C([O])[CH]O'),
    E0 = (628.903,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,1685,370,386.468,386.486,386.641,386.666,3706.96],'cm^-1')),
        HinderedRotor(inertia=(0.000272812,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311372,'amu*angstrom^2'), symmetry=1, barrier=(32.959,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116177,'amu*angstrom^2'), symmetry=1, barrier=(12.304,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22893,0.0647364,-0.000101485,8.07606e-08,-2.51739e-11,75735.8,28.6519], Tmin=(100,'K'), Tmax=(803.383,'K')), NASAPolynomial(coeffs=[11.0429,0.0150852,-7.30969e-06,1.39124e-09,-9.55671e-14,74184.3,-16.3881], Tmin=(803.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(628.903,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJO) + radical(Cds_S) + radical(CH2_triplet)"""),
)

species(
    label = '[O][C](C=O)C1[C]O1(13099)',
    structure = SMILES('[O]C=C([O])C1[C]O1'),
    E0 = (317.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.804531,0.0700281,-8.28129e-05,4.35019e-08,-7.96302e-12,38347,29.2056], Tmin=(100,'K'), Tmax=(1706.09,'K')), NASAPolynomial(coeffs=[18.2254,-0.00358839,7.4081e-06,-1.71575e-09,1.22864e-13,36074.2,-60.3995], Tmin=(1706.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CH2_triplet)"""),
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
    E0 = (122.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (918.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (684.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (856.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (533.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (125.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (245.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (130.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (186.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (186.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (356.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (381.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (355.083,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (359.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (446.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (366.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (387.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (389.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (574.176,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (383.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (410.045,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (355.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (335.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (122.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (130.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (129.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (147.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (708,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (323.432,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (470.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (341.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (392.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (272.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (334.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (373.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (278.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (254.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (698.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (726.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (125.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (130.924,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (186.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (186.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (368.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (413.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (133.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (461.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (476.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (125.569,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (364.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (576.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (447.145,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (490.337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (771.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (173.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (142.772,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (150.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (464.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (554.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (896.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (533.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (492.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (865.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (756.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (338.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (554.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (744.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (653.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (318.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['HCCO(2227)', 'OCHCO(3676)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[O]C=[C][CH][C]=O(10420)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C][O](6861)', '[O][C]=C[O](9592)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=C([O])C=[C][O](9624)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[O]C=C([O])[C][C]=O(13053)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C=C1OC1[C]=O(13039)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['O=[C][CH]C1=COO1(13054)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(122.654,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 121.7 to 122.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C1=COC1[C]=O(12991)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['O=[C][CH]C(O)=C=O(13055)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C(=C=O)C[C]=O(11828)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C1O[C]1[CH][C]=O(13056)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(234.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][C]1C([O])C1[C]=O(13057)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(259.102,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C1([CH][C]=O)[CH]O1(13058)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(232.443,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][C]1[CH]C(=O)C1[O](13059)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.77275e+09,'s^-1'), n=0.708429, Ea=(236.81,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;doublebond_intra;radadd_intra] for rate rule [R4_linear;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 236.7 to 236.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][CH]C1([O])[CH]C1=O(13060)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.2661e+10,'s^-1'), n=0.520686, Ea=(323.831,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[O]C(=C=O)[CH][C]=O(13061)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C][O](6861)', 'OCHCO(3676)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Ck_O;YJ] for rate rule [Ck_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O(4)', '[O]C=C=C[C]=O(10418)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[O][C]=C([O])[CH][C]=O(13062)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O]C(=[C]O)[CH][C]=O(13063)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][C]=C(O)[CH][C]=O(13064)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][C]=C([O])C[C]=O(11829)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(70985.8,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O][C]=C([O])[CH]C=O(13065)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.82652e+08,'s^-1'), n=1.30038, Ea=(212.995,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_2;Cd_rad_out;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C=C([O])C=C=O(11380)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C=C1C=C([O])O1(13066)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C1=COC([O])=C1(13067)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C=CC([O])=C=O(13068)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O][C]=[C]C([O])[CH][O](13069)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][C]1C=[C]OC1[O](13070)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.50361e+10,'s^-1'), n=0.23641, Ea=(200.792,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 200.7 to 200.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][CH]C1([O])C=[C]O1(13071)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(348.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[O]C#CC([O])=C[O](13072)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['HCCO(2227)', '[O][C]=C[O](9592)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C=[C]C([O])=C[O](13073)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.44274e+08,'s^-1'), n=1.26608, Ea=(149.637,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O][C]=[C]C(O)=C[O](13074)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.87605e+06,'s^-1'), n=1.79171, Ea=(150.238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[O]C=C([O])[C]=[C]O(13075)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[O][C]=[C]C([O])=CO(13076)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out;XH_out] for rate rule [R4H_SDS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O][C]=C([O])C=[C]O(13077)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_3;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH][O](1548)', '[O][C]=C[C]=O(9898)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[C]=O(2355)', '[CH]=C([O])[CH][O](10249)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C1OC1=C[C]=O(13078)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C1=CC(=O)C1[O](13079)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][CH]C(O)=C=C=O(13080)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]CC([O])=C=C=O(13081)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][CH][C]1OC1[C]=O(13082)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(245.451,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][C]C1O[CH]C1=O(12998)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.04143e+10,'s^-1'), n=0.464715, Ea=(290.863,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_linear;multiplebond_intra;radadd_intra_O] + [R4_linear;doublebond_intra;radadd_intra] for rate rule [R4_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['CO(2039)', '[CH]=C([O])[CH][O](10249)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.51e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH][O](1548)', 'O=[C]C=C=O(9594)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.96577e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Ck_O;YJ] for rate rule [Ck_O;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[O]CC([O])=[C][C]=O(13083)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(70985.8,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C1C(=O)C1[C]=O(13084)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][CH][C]1[CH]C(=O)O1(13085)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra] for rate rule [R4_linear;carbonyl_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['O=[C][CH][C]1[CH]OO1(13086)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(454.203,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[O][C]=CC([O])[C]=O(11384)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][CH][C]1C=[C]OO1(13087)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(6.44463e+09,'s^-1'), n=0.470283, Ea=(367.697,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra] for rate rule [R5_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 366.7 to 367.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction54',
    reactants = ['HCO(1372)', '[O][C][CH][C]=O(13088)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][C]1[CH]OC(=O)[CH]1(13089)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(50.8858,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;carbonyl_intra_H;radadd_intra_CO] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_CO]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C1=CC1([O])C=O(13090)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(20.1315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O]C=C=C([O])C=O(13091)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][C]1[CH]OO[C]=C1(13092)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(341.545,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[O][C]=[C]C([O])C=O(13093)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(3.64487e+09,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['HCCO(2227)', '[O][C][CH][O](10223)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[O][CH][C]([O])[C]=C[O](13094)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[O][C]C([C]=O)C=O(13095)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CO]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH][C][O](10218)', 'OCHCO(3676)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[O][C][CH]C([O])[C]=O(13096)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction65',
    reactants = ['HCO(1372)', '[O][C]=C[C]=O(9898)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(1.04e+12,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd_R;CO_pri_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction66',
    reactants = ['CO(2039)', '[O][CH][CH][C]=O(10226)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[O][C][C]=C([O])C[O](13097)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction68',
    reactants = ['[O][C][C]=C([O])[CH]O(13098)'],
    products = ['[O][C]=CC([O])=C[O](11388)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[O][C]=CC([O])=C[O](11388)'],
    products = ['[O][C](C=O)C1[C]O1(13099)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.30963e+11,'s^-1'), n=0.419784, Ea=(195.827,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra;radadd_intra_O] + [R4;multiplebond_intra;radadd_intra_O] + [R4;doublebond_intra;radadd_intra] for rate rule [R4;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

network(
    label = 'PDepNetwork #3071',
    isomers = [
        '[O][C]=CC([O])=C[O](11388)',
    ],
    reactants = [
        ('HCCO(2227)', 'OCHCO(3676)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #3071',
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


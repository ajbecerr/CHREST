species(
    label = '[CH]C([CH2])=C[C]=O(16420)',
    structure = SMILES('[CH]C([CH2])=C[C]=O'),
    E0 = (526.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.10045,'amu*angstrom^2'), symmetry=1, barrier=(48.2934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09767,'amu*angstrom^2'), symmetry=1, barrier=(48.2295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10075,'amu*angstrom^2'), symmetry=1, barrier=(48.3004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52118,0.0593111,-6.47653e-05,4.05501e-08,-1.0837e-11,63394.7,20.8288], Tmin=(100,'K'), Tmax=(884.972,'K')), NASAPolynomial(coeffs=[8.03996,0.0298466,-1.48236e-05,2.92786e-09,-2.08896e-13,62241,-9.82393], Tmin=(884.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(526.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
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
    label = 'C3H3(5450)',
    structure = SMILES('[CH]=C=C'),
    E0 = (338.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,180,180,2603.58],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09172,0.0173333,-8.20209e-06,-2.63358e-09,2.66049e-12,40755.9,8.10965], Tmin=(100,'K'), Tmax=(946.054,'K')), NASAPolynomial(coeffs=[6.98214,0.0072721,-2.37773e-06,3.99152e-10,-2.69331e-14,39733.9,-11.9544], Tmin=(946.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C3H3""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH][C]=C[C]=O(10445)',
    structure = SMILES('[CH][C]=C[C]=O'),
    E0 = (651.783,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07768,'amu*angstrom^2'), symmetry=1, barrier=(47.7699,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08207,'amu*angstrom^2'), symmetry=1, barrier=(47.8709,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.0581,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99066,0.0507484,-8.43786e-05,8.05615e-08,-3.01685e-11,78457.4,18.5596], Tmin=(100,'K'), Tmax=(815.064,'K')), NASAPolynomial(coeffs=[4.49774,0.025284,-1.3295e-05,2.60942e-09,-1.82246e-13,78485.9,9.65873], Tmin=(815.064,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.783,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Cds_S) + radical(C=CCJ=O) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C(=[CH])[CH2](17858)',
    structure = SMILES('[CH]C([CH])=C'),
    E0 = (707.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,333.478,333.562,333.637,333.655,333.685,333.706],'cm^-1')),
        HinderedRotor(inertia=(0.644549,'amu*angstrom^2'), symmetry=1, barrier=(50.8939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.644751,'amu*angstrom^2'), symmetry=1, barrier=(50.893,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (52.0746,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52506,0.0264751,7.08406e-06,-2.07737e-08,8.00454e-12,85138.3,15.7761], Tmin=(100,'K'), Tmax=(1061.32,'K')), NASAPolynomial(coeffs=[6.48666,0.0249945,-9.83289e-06,1.79345e-09,-1.24009e-13,83539.9,-7.14092], Tmin=(1061.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(707.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([CH])=C[C]=O(18034)',
    structure = SMILES('[CH]C([CH])=C[C]=O'),
    E0 = (745.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.15839,'amu*angstrom^2'), symmetry=1, barrier=(49.6256,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15706,'amu*angstrom^2'), symmetry=1, barrier=(49.5951,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15697,'amu*angstrom^2'), symmetry=1, barrier=(49.5929,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30487,0.0645276,-8.40225e-05,6.99907e-08,-2.49351e-11,89763.5,22.418], Tmin=(100,'K'), Tmax=(750.224,'K')), NASAPolynomial(coeffs=[5.94986,0.035525,-1.75639e-05,3.40665e-09,-2.38686e-13,89185.8,2.13817], Tmin=(750.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(745.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]C1=CC(=O)C1(18035)',
    structure = SMILES('[CH]C1=CC(=O)C1'),
    E0 = (370.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.46206,0.0225822,3.31251e-05,-4.96439e-08,1.71068e-11,44646.1,16.6719], Tmin=(100,'K'), Tmax=(1107.24,'K')), NASAPolynomial(coeffs=[9.33443,0.0274624,-1.3731e-05,2.79944e-09,-2.05307e-13,41303.2,-25.4067], Tmin=(1107.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(C)=C=C=O(18036)',
    structure = SMILES('[CH]=C(C)[C]=C=O'),
    E0 = (421.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92681,0.0496031,-6.72872e-05,5.6613e-08,-1.97635e-11,50727.9,20.8092], Tmin=(100,'K'), Tmax=(786.548,'K')), NASAPolynomial(coeffs=[5.86426,0.0250444,-1.18038e-05,2.25601e-09,-1.56626e-13,50248.8,3.65041], Tmin=(786.548,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)CJ=C=O) + radical(Cds_P)"""),
)

species(
    label = '[CH][C]1CC1[C]=O(18037)',
    structure = SMILES('[CH][C]1CC1[C]=O'),
    E0 = (689.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.28148,0.0365844,-2.58569e-05,9.55922e-09,-1.47789e-12,82958.7,22.4789], Tmin=(100,'K'), Tmax=(1465.97,'K')), NASAPolynomial(coeffs=[8.63832,0.0192391,-8.10876e-06,1.48793e-09,-1.01436e-13,81094.9,-10.6207], Tmin=(1465.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Tertalkyl) + radical(CC(C)CJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH2])[CH]C1=O(18038)',
    structure = SMILES('[CH]C1([CH2])[CH]C1=O'),
    E0 = (756.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.443029,0.0594698,-6.0347e-05,2.89029e-08,-5.11288e-12,91159.1,23.1558], Tmin=(100,'K'), Tmax=(1627.32,'K')), NASAPolynomial(coeffs=[17.8601,0.00526344,1.21505e-07,-1.72213e-10,1.48662e-14,86999.2,-64.7168], Tmin=(1627.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(756.741,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsCs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CCJC=O) + radical(CJC(C)2C=O) + radical(CCJ2_triplet)"""),
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
    label = '[CH]C([CH2])=C=C=O(18039)',
    structure = SMILES('[CH]C(=C)[C]=C=O'),
    E0 = (544.775,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,2120,512.5,787.5,353.072,353.073,353.073,353.073],'cm^-1')),
        HinderedRotor(inertia=(0.580783,'amu*angstrom^2'), symmetry=1, barrier=(51.3771,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.580784,'amu*angstrom^2'), symmetry=1, barrier=(51.3771,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12634,0.0443601,-4.36504e-05,2.87186e-08,-8.6425e-12,65586,20.8087], Tmin=(100,'K'), Tmax=(773.251,'K')), NASAPolynomial(coeffs=[5.15341,0.028702,-1.32775e-05,2.53367e-09,-1.77086e-13,65117.8,6.98305], Tmin=(773.251,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.775,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)CJ=C=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH](6993)',
    structure = SMILES('[CH]'),
    E0 = (585.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([4000],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (13.0186,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.17282,-0.00431623,8.93613e-06,-7.67944e-09,2.45632e-12,70459.2,-0.856057], Tmin=(100,'K'), Tmax=(949.528,'K')), NASAPolynomial(coeffs=[3.66649,-0.000619273,6.25229e-07,-1.09644e-10,6.55478e-15,70484.8,1.18918], Tmin=(949.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(585.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJ3)"""),
)

species(
    label = 'C=C=C[C]=O(9908)',
    structure = SMILES('C=C=C[C]=O'),
    E0 = (213.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.13978,'amu*angstrom^2'), symmetry=1, barrier=(26.2058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00684,0.0284132,-1.9291e-05,6.01948e-09,-7.90367e-13,25757,13.5223], Tmin=(100,'K'), Tmax=(1551.73,'K')), NASAPolynomial(coeffs=[6.76337,0.0187298,-9.93048e-06,1.99793e-09,-1.42457e-13,24591.2,-6.25135], Tmin=(1551.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]C([CH2])=[C][C]=O(18040)',
    structure = SMILES('[CH]C([CH2])=[C][C]=O'),
    E0 = (770.347,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.054,'amu*angstrom^2'), symmetry=1, barrier=(47.2254,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04888,'amu*angstrom^2'), symmetry=1, barrier=(47.1079,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04959,'amu*angstrom^2'), symmetry=1, barrier=(47.1242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26392,0.0656604,-8.99913e-05,6.82959e-08,-2.13263e-11,92744.9,22.453], Tmin=(100,'K'), Tmax=(775.619,'K')), NASAPolynomial(coeffs=[8.90928,0.0262296,-1.37299e-05,2.74312e-09,-1.9587e-13,91559,-12.4884], Tmin=(775.619,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(770.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CJC=O) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]C([CH2])=[C]C=O(18041)',
    structure = SMILES('[CH]C(=C)[C]=C[O]'),
    E0 = (558.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07691,'amu*angstrom^2'), symmetry=1, barrier=(47.7523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07652,'amu*angstrom^2'), symmetry=1, barrier=(47.7433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16763,0.0517738,-2.13603e-05,-2.06327e-08,1.54167e-11,67305.6,20.7387], Tmin=(100,'K'), Tmax=(905.023,'K')), NASAPolynomial(coeffs=[16.7502,0.0115593,-2.20424e-06,2.43448e-10,-1.47857e-14,63311.5,-59.3666], Tmin=(905.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(C)=[C][C]=O(18042)',
    structure = SMILES('[CH]C(C)=[C][C]=O'),
    E0 = (618.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.01172,'amu*angstrom^2'), symmetry=1, barrier=(46.2533,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01739,'amu*angstrom^2'), symmetry=1, barrier=(46.3838,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.01698,'amu*angstrom^2'), symmetry=1, barrier=(46.3742,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89422,0.0571525,-3.42237e-05,-6.89074e-08,9.36472e-11,74494.6,20.0945], Tmin=(100,'K'), Tmax=(454.516,'K')), NASAPolynomial(coeffs=[7.14472,0.031311,-1.61529e-05,3.17023e-09,-2.22214e-13,73807,-3.41011], Tmin=(454.516,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=O) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
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
    label = '[CH]=[C][CH2](16918)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,260.76,263.05,263.653],'cm^-1')),
        HinderedRotor(inertia=(1.04394,'amu*angstrom^2'), symmetry=1, barrier=(50.8215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2703,0.0129731,8.09839e-06,-1.37802e-08,4.71189e-12,73935.1,11.2848], Tmin=(100,'K'), Tmax=(1104.68,'K')), NASAPolynomial(coeffs=[4.54119,0.0159815,-6.32008e-06,1.15742e-09,-7.99393e-14,73190,2.92522], Tmin=(1104.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]1[CH]C(=O)C1(18043)',
    structure = SMILES('[CH][C]1C=C([O])C1'),
    E0 = (616.158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31001,0.0214673,4.27477e-05,-7.86475e-08,3.41349e-11,74182,18.1817], Tmin=(100,'K'), Tmax=(923.509,'K')), NASAPolynomial(coeffs=[14.8623,0.00758027,-4.47027e-07,-1.6773e-12,-3.92575e-15,70137.4,-50.7232], Tmin=(923.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(616.158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Allyl_T) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])C[C]=O(17227)',
    structure = SMILES('[CH]C(=[CH])C[C]=O'),
    E0 = (645.547,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1855,455,950,3120,650,792.5,1650,364.026,364.151,364.75],'cm^-1')),
        HinderedRotor(inertia=(0.522897,'amu*angstrom^2'), symmetry=1, barrier=(49.2289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.522968,'amu*angstrom^2'), symmetry=1, barrier=(49.2055,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.523269,'amu*angstrom^2'), symmetry=1, barrier=(49.2022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61984,0.0468768,-3.21843e-05,1.0267e-08,-1.29676e-12,77731.5,22.624], Tmin=(100,'K'), Tmax=(1831.34,'K')), NASAPolynomial(coeffs=[15.052,0.0175388,-8.15461e-06,1.51944e-09,-1.02636e-13,72811.7,-50.3054], Tmin=(1831.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.547,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=[CH])[CH]C=O(18044)',
    structure = SMILES('[CH]C([CH])=CC=O'),
    E0 = (584.95,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,2782.5,750,1395,475,1775,1000,223.176,235.447,242.618,246.72,257.158,282.848],'cm^-1')),
        HinderedRotor(inertia=(1.28694,'amu*angstrom^2'), symmetry=1, barrier=(50.2367,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.30522,'amu*angstrom^2'), symmetry=1, barrier=(50.6063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.16032,'amu*angstrom^2'), symmetry=1, barrier=(50.15,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23275,0.0466647,-2.56969e-05,6.15114e-09,-5.75479e-13,70407.9,19.5144], Tmin=(100,'K'), Tmax=(2290.88,'K')), NASAPolynomial(coeffs=[14.1075,0.0259308,-1.2121e-05,2.20041e-09,-1.44343e-13,64967.2,-47.6177], Tmin=(2290.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(584.95,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=C([CH2])[CH][C]=O(18045)',
    structure = SMILES('[C]C([CH2])=C[C]=O'),
    E0 = (825.179,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.13635,'amu*angstrom^2'), symmetry=1, barrier=(26.127,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13788,'amu*angstrom^2'), symmetry=1, barrier=(26.162,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36668,0.062756,-8.99015e-05,6.58075e-08,-1.92081e-11,99336.7,18.7239], Tmin=(100,'K'), Tmax=(836.731,'K')), NASAPolynomial(coeffs=[10.7129,0.0180784,-9.81253e-06,1.99978e-09,-1.44488e-13,97772.5,-24.7009], Tmin=(836.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(825.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(C=CCJ=O) + radical(CJ3)"""),
)

species(
    label = '[CH]=C1CC1[C]=O(18027)',
    structure = SMILES('[CH]=C1CC1[C]=O'),
    E0 = (479.254,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9608,0.0389643,-2.46108e-05,3.69334e-09,1.43504e-12,57719.4,19.8342], Tmin=(100,'K'), Tmax=(1063.53,'K')), NASAPolynomial(coeffs=[10.8517,0.0157067,-6.16826e-06,1.13408e-09,-7.92455e-14,55252.5,-26.3137], Tmin=(1063.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(CC(C)CJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=CC1[C]=O(17994)',
    structure = SMILES('C=C1[CH]C1[C]=O'),
    E0 = (432.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.25513,0.0299459,1.5203e-06,-2.40392e-08,1.14425e-11,52035.1,20.8698], Tmin=(100,'K'), Tmax=(985.538,'K')), NASAPolynomial(coeffs=[10.8668,0.0147549,-5.43568e-06,1.01165e-09,-7.31651e-14,49378,-25.4201], Tmin=(985.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.06,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(CCJCC=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=[C][CH]C1=CC1(18046)',
    structure = SMILES('O=[C]C=C1[CH]C1'),
    E0 = (365.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77533,0.0420529,-2.94011e-05,8.45911e-09,-6.29156e-13,44008,17.6714], Tmin=(100,'K'), Tmax=(1307.58,'K')), NASAPolynomial(coeffs=[12.9541,0.0152152,-7.0561e-06,1.37073e-09,-9.6825e-14,40455.5,-41.6631], Tmin=(1307.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.19,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(Methylene_cyclopropane) + radical(Allyl_S) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]=C=C[C]=O(10456)',
    structure = SMILES('[CH]=C=C[C]=O'),
    E0 = (368.388,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([540,610,2055,1855,455,950,3120,650,792.5,1650,184.363,184.368,184.447,1693.07],'cm^-1')),
        HinderedRotor(inertia=(0.00837741,'amu*angstrom^2'), symmetry=1, barrier=(2.94921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.0581,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58161,0.0342376,-4.2892e-05,3.11333e-08,-9.52264e-12,44355.1,15.702], Tmin=(100,'K'), Tmax=(783.956,'K')), NASAPolynomial(coeffs=[6.13378,0.0161129,-8.21217e-06,1.64147e-09,-1.17678e-13,43798.2,-0.570526], Tmin=(783.956,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.388,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCJ=O) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=C([CH2])C=C=O(16875)',
    structure = SMILES('[CH]C(=C)C=C=O'),
    E0 = (357.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.17559,'amu*angstrom^2'), symmetry=1, barrier=(50.0211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.18674,'amu*angstrom^2'), symmetry=1, barrier=(50.2774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65598,0.0548559,-6.1877e-05,4.30589e-08,-1.28209e-11,43023.8,20.3781], Tmin=(100,'K'), Tmax=(802.514,'K')), NASAPolynomial(coeffs=[6.905,0.0286946,-1.29816e-05,2.443e-09,-1.69032e-13,42181.3,-3.7909], Tmin=(802.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])[C]=[C][O](18047)',
    structure = SMILES('[CH]C([CH2])[C][C]=O'),
    E0 = (1001.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30286,0.0625251,-9.24047e-05,7.23315e-08,-2.19885e-11,120552,26.792], Tmin=(100,'K'), Tmax=(910.779,'K')), NASAPolynomial(coeffs=[9.79638,0.0179949,-7.16185e-06,1.22242e-09,-7.80093e-14,119305,-11.7448], Tmin=(910.779,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(1001.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH][C]1C=[C]OC1(18048)',
    structure = SMILES('[CH][C]1C=[C]OC1'),
    E0 = (680.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61024,0.0125509,6.89148e-05,-1.08692e-07,4.63887e-11,81883,17.6567], Tmin=(100,'K'), Tmax=(899.319,'K')), NASAPolynomial(coeffs=[14.9574,0.00507204,2.26411e-06,-6.28375e-10,4.27548e-14,77743.9,-51.2667], Tmin=(899.319,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CCJ(C)CO) + radical(C=CJO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH2])C=[C]O1(18049)',
    structure = SMILES('[CH]C1([CH2])C=[C]O1'),
    E0 = (812.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26602,0.0447142,-3.91123e-06,-4.59024e-08,2.65466e-11,97792.3,22.0047], Tmin=(100,'K'), Tmax=(904.944,'K')), NASAPolynomial(coeffs=[21.1917,-0.0017122,4.00853e-06,-8.79367e-10,5.86717e-14,92480.6,-81.5574], Tmin=(904.944,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(812.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CJC(C)OC) + radical(C=CJO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)[C]=[C]O(18050)',
    structure = SMILES('[CH]C(=C)[C]=[C]O'),
    E0 = (656.962,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1670,1700,300,440,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.08581,'amu*angstrom^2'), symmetry=1, barrier=(47.9569,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08544,'amu*angstrom^2'), symmetry=1, barrier=(47.9484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08584,'amu*angstrom^2'), symmetry=1, barrier=(47.9575,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.497701,0.0661224,-7.18028e-05,3.88804e-08,-7.94893e-12,79149.8,24.8547], Tmin=(100,'K'), Tmax=(1362.67,'K')), NASAPolynomial(coeffs=[16.418,0.0109534,-1.7873e-06,8.30256e-11,2.50094e-15,75594.2,-54.0037], Tmin=(1362.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.962,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(AllylJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(=[CH])C=[C]O(18051)',
    structure = SMILES('[CH]C(=[CH])C=[C]O'),
    E0 = (705.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.02862,'amu*angstrom^2'), symmetry=1, barrier=(46.6419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02909,'amu*angstrom^2'), symmetry=1, barrier=(46.6527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02897,'amu*angstrom^2'), symmetry=1, barrier=(46.65,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.253805,0.0658817,-6.84477e-05,3.43129e-08,-6.41179e-12,84948.5,25.8533], Tmin=(100,'K'), Tmax=(1520.84,'K')), NASAPolynomial(coeffs=[18.1832,0.0079176,-6.18272e-07,-9.3113e-11,1.20882e-14,80744.8,-64.0536], Tmin=(1520.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(705.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C[CH][C]=O(16166)',
    structure = SMILES('[CH]=[C]C[CH][C]=O'),
    E0 = (689.327,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1685,370,1855,455,950,3120,650,792.5,1650,309.189],'cm^-1')),
        HinderedRotor(inertia=(0.187291,'amu*angstrom^2'), symmetry=1, barrier=(12.7217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.187692,'amu*angstrom^2'), symmetry=1, barrier=(12.7303,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0113573,'amu*angstrom^2'), symmetry=1, barrier=(74.5251,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3431.29,'J/mol'), sigma=(5.7779,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.96 K, Pc=40.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79282,0.0506005,-6.3105e-05,4.45979e-08,-1.27993e-11,82984.4,23.7372], Tmin=(100,'K'), Tmax=(848.507,'K')), NASAPolynomial(coeffs=[8.37481,0.0195721,-8.25322e-06,1.50146e-09,-1.01672e-13,81867.4,-6.93589], Tmin=(848.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCHO) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=CC([O])=C1(18052)',
    structure = SMILES('C=C1[CH]C([O])=C1'),
    E0 = (259.887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7916,0.0281381,4.38455e-05,-9.26794e-08,4.2041e-11,31355.7,14.5048], Tmin=(100,'K'), Tmax=(920.327,'K')), NASAPolynomial(coeffs=[20.4786,-3.62755e-05,3.31072e-06,-6.90594e-10,4.09312e-14,25669.6,-86.3022], Tmin=(920.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.887,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C([CH2])=[C][C]=O(16425)',
    structure = SMILES('[CH2]C([CH2])=[C][C]=O'),
    E0 = (551.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.04859,'amu*angstrom^2'), symmetry=1, barrier=(24.1091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04797,'amu*angstrom^2'), symmetry=1, barrier=(24.0948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04673,'amu*angstrom^2'), symmetry=1, barrier=(24.0664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37854,0.0618654,-7.68575e-05,4.8691e-08,-1.24091e-11,66380.3,21.2148], Tmin=(100,'K'), Tmax=(949.276,'K')), NASAPolynomial(coeffs=[11.3264,0.0199475,-1.06207e-05,2.17348e-09,-1.58306e-13,64491.7,-26.2601], Tmin=(949.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(Allyl_P) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH][C][CH2](17602)',
    structure = SMILES('[CH][C][CH2]'),
    E0 = (981.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,917.882,918.218,918.236],'cm^-1')),
        HinderedRotor(inertia=(0.00702614,'amu*angstrom^2'), symmetry=1, barrier=(4.20306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182792,'amu*angstrom^2'), symmetry=1, barrier=(4.20276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96474,0.0202917,-1.68078e-05,7.2211e-09,-1.21436e-12,118107,14.5039], Tmin=(100,'K'), Tmax=(1441.81,'K')), NASAPolynomial(coeffs=[8.09928,0.00604662,-1.9874e-06,3.68258e-10,-2.60955e-14,116627,-12.1458], Tmin=(1441.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(981.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]([CH2])[C]=C[O](18053)',
    structure = SMILES('[CH]C([CH2])=[C][CH][O]'),
    E0 = (880.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,1685,370,602.169,602.368,602.57,602.597,602.872],'cm^-1')),
        HinderedRotor(inertia=(0.201586,'amu*angstrom^2'), symmetry=1, barrier=(51.9519,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201955,'amu*angstrom^2'), symmetry=1, barrier=(51.9537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201919,'amu*angstrom^2'), symmetry=1, barrier=(51.9564,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98236,0.0446075,-2.90363e-05,9.52842e-09,-1.31884e-12,106024,24.0403], Tmin=(100,'K'), Tmax=(1575.26,'K')), NASAPolynomial(coeffs=[9.21187,0.0262499,-1.15557e-05,2.13049e-09,-1.44757e-13,103746,-14.1231], Tmin=(1575.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(880.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_P) + radical(C=CCJO) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])[C]=C=O(18054)',
    structure = SMILES('[CH]C([CH2])[C]=C=O'),
    E0 = (730.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,221.637,221.64,221.641,1530.87],'cm^-1')),
        HinderedRotor(inertia=(0.00343174,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00343175,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.8513,'amu*angstrom^2'), symmetry=1, barrier=(64.5344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23692,0.0656241,-0.000101439,8.42522e-08,-2.74951e-11,87919.5,23.8094], Tmin=(100,'K'), Tmax=(848.089,'K')), NASAPolynomial(coeffs=[9.03207,0.0213537,-9.86566e-06,1.83382e-09,-1.24022e-13,86867.2,-10.9221], Tmin=(848.089,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(730.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O) + radical(CC(C)CJ=C=O) + radical(CCJ2_triplet)"""),
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
    E0 = (526.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1089.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1201.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (957.363,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (534.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (589.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (689.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (756.741,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (608.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (773.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (814.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (982.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (721.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (774.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (1094.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (655.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (837.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (793.964,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (1036.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (531.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (529.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (534.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (789.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (769.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (526.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (1024.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (680.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (812.152,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (802.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (850.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (738.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (934.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (534.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (735.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (1202.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (903.794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (907.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['HCCO(2227)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(20)', '[CH][C]=C[C]=O(10445)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(2355)', '[CH]C(=[CH])[CH2](17858)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]C([CH])=C[C]=O(18034)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(2e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH]C1=CC(=O)C1(18035)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SDS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH]C(C)=C=C=O(18036)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH][C]1CC1[C]=O(18037)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(162.851,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 158.5 to 162.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH]C1([CH2])[CH]C1=O(18038)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(230.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 228.4 to 230.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['CO(2039)', '[CH]C(=[CH])[CH2](17858)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.02e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]C([CH2])=C=C=O(18039)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Cdd;HJ] for rate rule [Ca_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH](6993)', 'C=C=C[C]=O(9908)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.70446,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;CH_quartet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH]C([CH2])=[C][C]=O(18040)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C([CH2])=[C]C=O(18041)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C(C)=[C][C]=O(18042)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.17353e+06,'s^-1'), n=1.89718, Ea=(155.956,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C][O](6861)', '[CH]=[C][CH2](16918)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH][C]1[CH]C(=O)C1(18043)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.77275e+09,'s^-1'), n=0.708429, Ea=(129.065,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;doublebond_intra;radadd_intra] for rate rule [R4_linear;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]C(=[CH])C[C]=O(17227)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(26875.4,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]C(=[CH])[CH]C=O(18044)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.73633e+06,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_2;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[C]=C([CH2])[CH][C]=O(18045)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH]=C1CC1[C]=O(18027)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH2]C1=CC1[C]=O(17994)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['O=[C][CH]C1=CC1(18046)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['CH2(T)(20)', '[CH]=C=C[C]=O(10456)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C][O](6861)', 'C3H3(5450)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(6.8966,'m^3/(mol*s)'), n=1.84023, Ea=(5.48899,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cd_R;Y_1centerbirad] + [Cdd_Cds;YJ] for rate rule [Cdd_Cds;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH]=C([CH2])C=C=O(16875)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C([CH2])[C]=[C][O](18047)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH][C]1C=[C]OC1(18048)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(5.50361e+10,'s^-1'), n=0.23641, Ea=(153.872,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 146.7 to 153.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH]C1([CH2])C=[C]O1(18049)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(285.767,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 283.2 to 285.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['HCCO(2227)', '[CH]=[C][CH2](16918)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]C(=C)[C]=[C]O(18050)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C(=[CH])C=[C]O(18051)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]C[CH][C]=O(16166)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]C([CH2])=C[C]=O(16420)'],
    products = ['[CH2]C1=CC([O])=C1(18052)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH2])=[C][C]=O(16425)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.79573e+07,'s^-1'), n=1.84533, Ea=(184.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SD;Y_rad_out;Cd_H_out_singleH] + [R3H_SD;Cd_rad_out;Cd_H_out_single] for rate rule [R3H_SD;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['HCCO(2227)', '[CH][C][CH2](17602)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH][C]([CH2])[C]=C[O](18053)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]C([CH2])[C]=C=O(18054)'],
    products = ['[CH]C([CH2])=C[C]=O(16420)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4048',
    isomers = [
        '[CH]C([CH2])=C[C]=O(16420)',
    ],
    reactants = [
        ('HCCO(2227)', 'C3H3(5450)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4048',
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


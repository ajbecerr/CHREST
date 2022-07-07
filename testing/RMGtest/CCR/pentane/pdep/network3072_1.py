species(
    label = '[O][C]=CO[C]=C[O](11389)',
    structure = SMILES('[O]C=[C]O[CH][C]=O'),
    E0 = (268.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07527,'amu*angstrom^2'), symmetry=1, barrier=(24.7225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07753,'amu*angstrom^2'), symmetry=1, barrier=(24.7746,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.07641,'amu*angstrom^2'), symmetry=1, barrier=(24.7487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.247163,0.0719586,-8.70005e-05,4.56261e-08,-8.61352e-12,32514.6,30.6943], Tmin=(100,'K'), Tmax=(1540.81,'K')), NASAPolynomial(coeffs=[23.5694,-0.00618032,4.94705e-06,-1.02716e-09,7.07614e-14,27111.3,-88.2201], Tmin=(1540.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CsCJ=O)"""),
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
    label = '[C]=C[O](6859)',
    structure = SMILES('[C]C=O'),
    E0 = (491.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53714,0.00730776,2.32451e-06,-8.52039e-09,3.75032e-12,59149.2,7.56143], Tmin=(100,'K'), Tmax=(1018.09,'K')), NASAPolynomial(coeffs=[6.47558,0.00262438,-8.8465e-07,2.00876e-10,-1.68063e-14,58195.3,-8.41394], Tmin=(1018.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJ3)"""),
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
    label = '[CH]=[C]OC=[C][O](9625)',
    structure = SMILES('[CH]=[C]O[CH][C]=O'),
    E0 = (583.343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,1855,455,950,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.899785,'amu*angstrom^2'), symmetry=1, barrier=(20.6878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900286,'amu*angstrom^2'), symmetry=1, barrier=(20.6993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.901258,'amu*angstrom^2'), symmetry=1, barrier=(20.7217,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3435.15,'J/mol'), sigma=(5.65948,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.56 K, Pc=43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.971109,0.0589565,-7.59981e-05,4.40442e-08,-9.47793e-12,70275.6,25.447], Tmin=(100,'K'), Tmax=(1245.33,'K')), NASAPolynomial(coeffs=[18.0697,-0.000580964,1.27573e-06,-3.00035e-10,2.17823e-14,66374.9,-59.3578], Tmin=(1245.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(583.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CsCJ=O) + radical(Cds_P)"""),
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
    label = '[O]C=[C]O[C][C]=O(13038)',
    structure = SMILES('[O]C=[C]O[C][C]=O'),
    E0 = (535.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.13969,'amu*angstrom^2'), symmetry=1, barrier=(26.2036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13355,'amu*angstrom^2'), symmetry=1, barrier=(26.0626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.13288,'amu*angstrom^2'), symmetry=1, barrier=(26.0472,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0564172,0.0680436,-8.34651e-05,4.36125e-08,-8.17489e-12,64627.5,29.2704], Tmin=(100,'K'), Tmax=(1553.4,'K')), NASAPolynomial(coeffs=[23.3322,-0.0080361,5.30836e-06,-1.05627e-09,7.13784e-14,59273.9,-87.7105], Tmin=(1553.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(C=CJO) + radical(CH2_triplet) + radical(CsCJ=O)"""),
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
    label = 'O=[C][CH]OC1=CO1(13040)',
    structure = SMILES('O=C=CO[C]1[CH]O1'),
    E0 = (165.478,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.51763,0.0706258,-8.50767e-05,4.57288e-08,-8.64472e-12,20089.9,27.3452], Tmin=(100,'K'), Tmax=(1625.7,'K')), NASAPolynomial(coeffs=[19.2488,-0.00349171,6.82232e-06,-1.59901e-09,1.16078e-13,17030.4,-67.2646], Tmin=(1625.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.478,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + ring(Ethylene_oxide) + radical(Cs_P) + radical(CCsJO)"""),
)

species(
    label = 'O=[C][CH]OC=C=O(11370)',
    structure = SMILES('O=[C][CH]OC=C=O'),
    E0 = (69.1855,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1855,455,950,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.03712,'amu*angstrom^2'), symmetry=1, barrier=(23.8454,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03752,'amu*angstrom^2'), symmetry=1, barrier=(23.8546,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.04042,'amu*angstrom^2'), symmetry=1, barrier=(23.9214,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.177771,0.0699276,-8.67508e-05,4.75727e-08,-9.54369e-12,8470.89,24.7339], Tmin=(100,'K'), Tmax=(1398.1,'K')), NASAPolynomial(coeffs=[22.0931,-0.00330886,3.12753e-06,-6.74876e-10,4.75204e-14,3372.65,-84.6566], Tmin=(1398.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.1855,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CCsJOC(O)) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C#CO[CH][C]=O(12992)',
    structure = SMILES('O=[C][CH]O[C]=C=O'),
    E0 = (308.93,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1685,370,1855,455,950,2120,512.5,787.5,220.88,220.88,220.88],'cm^-1')),
        HinderedRotor(inertia=(0.551474,'amu*angstrom^2'), symmetry=1, barrier=(19.0926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551474,'amu*angstrom^2'), symmetry=1, barrier=(19.0926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.551474,'amu*angstrom^2'), symmetry=1, barrier=(19.0926,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.873407,0.0659581,-9.35916e-05,6.03118e-08,-1.46168e-11,37270.7,26.0595], Tmin=(100,'K'), Tmax=(1031.02,'K')), NASAPolynomial(coeffs=[17.389,0.00188441,-3.74386e-07,3.77442e-11,-1.90636e-15,33865.1,-54.1237], Tmin=(1031.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.93,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=[C]OC=C=O(9621)',
    structure = SMILES('C#CO[CH][C]=O'),
    E0 = (292.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2175,525,1855,455,950,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29406,'amu*angstrom^2'), symmetry=1, barrier=(29.7529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29132,'amu*angstrom^2'), symmetry=1, barrier=(29.6899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29532,'amu*angstrom^2'), symmetry=1, barrier=(29.7819,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.569776,0.060789,-7.4117e-05,3.96185e-08,-7.70243e-12,35298.6,20.788], Tmin=(100,'K'), Tmax=(1463.36,'K')), NASAPolynomial(coeffs=[20.2314,-0.00380299,3.21204e-06,-6.76284e-10,4.69603e-14,30705.7,-77.5851], Tmin=(1463.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CCsJOCs) + radical(CsCJ=O)"""),
)

species(
    label = 'O=[C][C]O[CH][C]=O(12989)',
    structure = SMILES('[O][C]=[C]O[CH][C]=O'),
    E0 = (508.661,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1670,1700,300,440,1855,455,950,302.173,302.174,302.179,302.18],'cm^-1')),
        HinderedRotor(inertia=(0.314059,'amu*angstrom^2'), symmetry=1, barrier=(20.3502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314057,'amu*angstrom^2'), symmetry=1, barrier=(20.3502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314072,'amu*angstrom^2'), symmetry=1, barrier=(20.3502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.600614,0.0663372,-8.86621e-05,5.23924e-08,-1.14272e-11,61307.5,30.7705], Tmin=(100,'K'), Tmax=(1241.4,'K')), NASAPolynomial(coeffs=[20.1366,-0.00280515,2.36801e-06,-5.1271e-10,3.66177e-14,56934.4,-65.7814], Tmin=(1241.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(508.661,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CsCJ=O) + radical(C=CJO)"""),
)

species(
    label = 'O=[C][CH]O[C]=[C]O(12993)',
    structure = SMILES('O=[C][CH]O[C]=[C]O'),
    E0 = (367.198,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1670,1700,300,440,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.866589,'amu*angstrom^2'), symmetry=1, barrier=(19.9246,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867249,'amu*angstrom^2'), symmetry=1, barrier=(19.9398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.866859,'amu*angstrom^2'), symmetry=1, barrier=(19.9308,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.867076,'amu*angstrom^2'), symmetry=1, barrier=(19.9358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0358009,0.0760287,-0.00010193,5.95114e-08,-1.26239e-11,44320.1,31.6399], Tmin=(100,'K'), Tmax=(1332.46,'K')), NASAPolynomial(coeffs=[23.1334,-0.00617703,4.85594e-06,-1.04301e-09,7.46377e-14,39268.9,-82.5733], Tmin=(1332.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(367.198,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(C=CJO) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]=CO[CH][C]=O(11385)',
    structure = SMILES('O=[C][CH]O[CH][C]=O'),
    E0 = (263.703,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,1850,1860,440,470,900,1000,411.165,709.003,709.032],'cm^-1')),
        HinderedRotor(inertia=(0.0475667,'amu*angstrom^2'), symmetry=1, barrier=(16.9686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.738035,'amu*angstrom^2'), symmetry=1, barrier=(16.9689,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.738043,'amu*angstrom^2'), symmetry=1, barrier=(16.9691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.73805,'amu*angstrom^2'), symmetry=1, barrier=(16.9692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3796.39,'J/mol'), sigma=(6.09434,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.99 K, Pc=38.06 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.931531,0.068822,-9.82834e-05,6.68224e-08,-1.74731e-11,31825.6,25.9793], Tmin=(100,'K'), Tmax=(946.421,'K')), NASAPolynomial(coeffs=[14.9006,0.00978491,-4.71839e-06,9.17215e-10,-6.4812e-14,29181.3,-40.6448], Tmin=(946.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.703,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CCsJOCs) + radical(CCsJOCs) + radical(CsCJ=O) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]=[C]OC[C]=O(11810)',
    structure = SMILES('[O][C]=[C]OC[C]=O'),
    E0 = (314.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1670,1700,300,440,396.451,396.494,396.527,396.555],'cm^-1')),
        HinderedRotor(inertia=(0.177993,'amu*angstrom^2'), symmetry=1, barrier=(19.8612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177991,'amu*angstrom^2'), symmetry=1, barrier=(19.8619,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178008,'amu*angstrom^2'), symmetry=1, barrier=(19.861,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.848984,0.0582439,-6.58828e-05,3.44247e-08,-6.71863e-12,37976.6,30.4119], Tmin=(100,'K'), Tmax=(1391.17,'K')), NASAPolynomial(coeffs=[18.1606,0.00280663,-4.12421e-09,-7.05398e-11,6.0402e-15,33707.8,-56.8528], Tmin=(1391.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CsCJ=O) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[O][C]=[C]O[CH]C=O(12994)',
    structure = SMILES('[O]C=CO[C][C]=O'),
    E0 = (296.233,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.23309,'amu*angstrom^2'), symmetry=1, barrier=(28.3512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23488,'amu*angstrom^2'), symmetry=1, barrier=(28.3924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23457,'amu*angstrom^2'), symmetry=1, barrier=(28.3851,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.583281,0.0563898,-2.24434e-05,-3.8969e-08,2.65976e-11,35769.2,23.838], Tmin=(100,'K'), Tmax=(915.028,'K')), NASAPolynomial(coeffs=[27.6182,-0.011944,7.86048e-06,-1.5118e-09,9.75886e-14,28734.8,-115.592], Tmin=(915.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.233,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C=[C]OC=C=O(11379)',
    structure = SMILES('[O]C=[C]OC=C=O'),
    E0 = (114.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2120,512.5,787.5,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.10413,'amu*angstrom^2'), symmetry=1, barrier=(25.3861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.10495,'amu*angstrom^2'), symmetry=1, barrier=(25.4049,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.735088,0.058878,-6.71021e-05,3.53139e-08,-6.86317e-12,13949.3,27.8211], Tmin=(100,'K'), Tmax=(1447.26,'K')), NASAPolynomial(coeffs=[18.0882,0.00202857,1.03083e-06,-3.1427e-10,2.42441e-14,9857.18,-59.0971], Tmin=(1447.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C=C1OC=C1[O](13041)',
    structure = SMILES('[O]C=C1OC=C1[O]'),
    E0 = (-46.1651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24126,0.0419254,8.25171e-06,-6.46175e-08,3.49209e-11,-5435.4,19.2356], Tmin=(100,'K'), Tmax=(901.507,'K')), NASAPolynomial(coeffs=[24.0065,-0.00832085,7.39123e-06,-1.51962e-09,1.01619e-13,-11602.8,-99.6738], Tmin=(901.507,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-46.1651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[O]C#CO[C]=C[O](13042)',
    structure = SMILES('[O]C=[C]O[C]=C=O'),
    E0 = (354.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1670,1700,300,440,2120,512.5,787.5,356.466,356.928,356.936,357.59],'cm^-1')),
        HinderedRotor(inertia=(0.192972,'amu*angstrom^2'), symmetry=1, barrier=(17.4568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193305,'amu*angstrom^2'), symmetry=1, barrier=(17.4535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44712,0.0547378,-7.34454e-05,4.75389e-08,-1.17694e-11,42748.4,28.3931], Tmin=(100,'K'), Tmax=(1003.32,'K')), NASAPolynomial(coeffs=[13.3473,0.00729469,-2.51637e-06,4.0951e-10,-2.61326e-14,40360.4,-29.058], Tmin=(1003.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]O[C]=C[O](13043)',
    structure = SMILES('[O]C=[C]O[C]=C[O]'),
    E0 = (314.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,215.956,215.956,215.956,215.956,215.956],'cm^-1')),
        HinderedRotor(inertia=(0.815223,'amu*angstrom^2'), symmetry=1, barrier=(26.9794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.815224,'amu*angstrom^2'), symmetry=1, barrier=(26.9794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.279177,0.0612237,-6.82464e-05,3.42816e-08,-6.23473e-12,37994.4,32.5101], Tmin=(100,'K'), Tmax=(1611.85,'K')), NASAPolynomial(coeffs=[19.3549,-0.00058158,2.73281e-06,-6.43823e-10,4.58836e-14,33724.3,-62.7961], Tmin=(1611.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[O]C=[C]O[C]=[C]O(13044)',
    structure = SMILES('[O]C=[C]O[C]=[C]O'),
    E0 = (412.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,1670,1685,1700,300,370,440,255.529,255.547,255.548,255.549],'cm^-1')),
        HinderedRotor(inertia=(0.411578,'amu*angstrom^2'), symmetry=1, barrier=(19.0695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411515,'amu*angstrom^2'), symmetry=1, barrier=(19.0696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.411494,'amu*angstrom^2'), symmetry=1, barrier=(19.0696,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.535356,0.0648339,-8.18504e-05,4.67906e-08,-9.78302e-12,49797.8,33.983], Tmin=(100,'K'), Tmax=(1357.17,'K')), NASAPolynomial(coeffs=[19.1034,-0.000792617,2.73102e-06,-6.75614e-10,5.07933e-14,45761.8,-57.569], Tmin=(1357.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(412.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[O][C]=[C]O[C]=CO(13045)',
    structure = SMILES('O=[C][C]O[C]=CO'),
    E0 = (394.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05683,'amu*angstrom^2'), symmetry=1, barrier=(24.2987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05661,'amu*angstrom^2'), symmetry=1, barrier=(24.2936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05769,'amu*angstrom^2'), symmetry=1, barrier=(24.3184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05709,'amu*angstrom^2'), symmetry=1, barrier=(24.3045,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.723179,0.0780647,-9.77637e-05,5.19083e-08,-9.80799e-12,47641.5,30.2508], Tmin=(100,'K'), Tmax=(1563.68,'K')), NASAPolynomial(coeffs=[25.7038,-0.0105269,7.35509e-06,-1.49296e-09,1.02255e-13,41942.9,-100.853], Tmin=(1563.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=CJO) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]=[C]OC=[C]O(12996)',
    structure = SMILES('O=[C][C]OC=[C]O'),
    E0 = (394.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.05683,'amu*angstrom^2'), symmetry=1, barrier=(24.2987,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05661,'amu*angstrom^2'), symmetry=1, barrier=(24.2936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05769,'amu*angstrom^2'), symmetry=1, barrier=(24.3184,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05709,'amu*angstrom^2'), symmetry=1, barrier=(24.3045,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.723179,0.0780647,-9.77637e-05,5.19083e-08,-9.80799e-12,47641.5,30.2508], Tmin=(100,'K'), Tmax=(1563.68,'K')), NASAPolynomial(coeffs=[25.7038,-0.0105269,7.35509e-06,-1.49296e-09,1.02255e-13,41942.9,-100.853], Tmin=(1563.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=CJO) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C1=CO[C][CH]O1(13046)',
    structure = SMILES('O=C1[CH]O[C][CH]O1'),
    E0 = (312.37,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32384,0.0427389,-3.67977e-06,-3.48315e-08,1.7522e-11,37679.9,17.7834], Tmin=(100,'K'), Tmax=(1016.55,'K')), NASAPolynomial(coeffs=[20.1456,0.00568466,-3.6097e-06,9.33845e-10,-8.08005e-14,31941.1,-82.7346], Tmin=(1016.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(CCsJOCs) + radical(CCsJOC(O)) + radical(CH2_triplet)"""),
)

species(
    label = '[C]1[CH]OO[C]=CO1(13047)',
    structure = SMILES('[C]1[CH]OO[C]=CO1'),
    E0 = (707.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38973,0.00305772,0.000186287,-3.01295e-07,1.32977e-10,85234.5,22.9705], Tmin=(100,'K'), Tmax=(895.848,'K')), NASAPolynomial(coeffs=[45.3757,-0.0489352,3.15511e-05,-6.20835e-09,4.15036e-13,71558.9,-216.74], Tmin=(895.848,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(707.471,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cycloheptane) + radical(CCsJOOC) + radical(CH2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[O]C1=CO[C]C1[O](13048)',
    structure = SMILES('[O]C1=CO[C]C1[O]'),
    E0 = (313.265,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47201,0.0331759,3.61126e-05,-9.92636e-08,4.95512e-11,37789.4,18.0898], Tmin=(100,'K'), Tmax=(887.187,'K')), NASAPolynomial(coeffs=[25.624,-0.014007,1.15518e-05,-2.40646e-09,1.65276e-13,31075.3,-109.226], Tmin=(887.187,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.265,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(CH2_triplet)"""),
)

species(
    label = '[O]C1[C]OC=[C]O1(13049)',
    structure = SMILES('[O]C1[C]OC=[C]O1'),
    E0 = (442.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0584504,0.0561749,-5.97681e-05,2.96464e-08,-5.19465e-12,53443.2,22.1548], Tmin=(100,'K'), Tmax=(1777.76,'K')), NASAPolynomial(coeffs=[14.0994,0.00291609,3.4505e-06,-9.16289e-10,6.73393e-14,51874.7,-44.034], Tmin=(1777.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(23dihydro14dioxin) + radical(CCOJ) + radical(CH2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[O][C]C1O[C]=CO1(12999)',
    structure = SMILES('[O][C]C1O[C]=CO1'),
    E0 = (437.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53524,0.0295711,4.78061e-05,-1.06057e-07,4.89645e-11,52672.4,22.0088], Tmin=(100,'K'), Tmax=(921.095,'K')), NASAPolynomial(coeffs=[25.3087,-0.00966279,7.46454e-06,-1.41703e-09,8.72441e-14,45577.7,-105.469], Tmin=(921.095,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.015,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(CCOJ) + radical(C=CJO) + radical(CH2_triplet)"""),
)

species(
    label = '[O]C1[C]OC1[C]=O(13050)',
    structure = SMILES('[O]C1[C]OC1[C]=O'),
    E0 = (465.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7739,0.0358001,7.47834e-06,-5.38548e-08,2.98307e-11,56126,22.8791], Tmin=(100,'K'), Tmax=(869.312,'K')), NASAPolynomial(coeffs=[18.1569,-0.0011997,5.08984e-06,-1.23063e-09,8.97467e-14,51827.2,-62.207], Tmin=(869.312,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.886,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(CC(C)OJ) + radical(CH2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[O][CH][C]O[C]=C[O](13051)',
    structure = SMILES('[O][CH][C]O[C]=C[O]'),
    E0 = (694.071,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355933,0.0790331,-0.000116066,7.82475e-08,-1.99356e-11,83609.7,27.3156], Tmin=(100,'K'), Tmax=(978.87,'K')), NASAPolynomial(coeffs=[18.6423,0.00430694,-1.55438e-06,2.56735e-10,-1.65478e-14,80029.8,-60.5144], Tmin=(978.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(694.071,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(CCsJOH) + radical(CH2_triplet) + radical(C=CJO)"""),
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
    label = '[O]C[C]O[C]=C=O(13052)',
    structure = SMILES('[O]C[C]O[C]=C=O'),
    E0 = (553.786,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,2120,512.5,787.5,313.868,313.869,313.869,313.872,313.872,313.873],'cm^-1')),
        HinderedRotor(inertia=(0.00171118,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213791,'amu*angstrom^2'), symmetry=1, barrier=(14.9458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.213791,'amu*angstrom^2'), symmetry=1, barrier=(14.9458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.868747,0.0745506,-0.000127065,1.06738e-07,-3.44944e-11,66712.5,26.0119], Tmin=(100,'K'), Tmax=(852.034,'K')), NASAPolynomial(coeffs=[11.5883,0.0155339,-7.86341e-06,1.49593e-09,-1.01562e-13,65201.3,-22.1354], Tmin=(852.034,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(553.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CCOJ) + radical(CH2_triplet) + radical(C=CJO)"""),
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
    E0 = (268.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (684.807,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (751.399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1102.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (747.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (271.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (271.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (291.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (528.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (453.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (379.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (557.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (720.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (525.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (473.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (527.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (329.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (268.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (277.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (581.557,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (518.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (563.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (427.555,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (508.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (327.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (707.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (363.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (442.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (437.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (465.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (716.372,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (841.489,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (703.858,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['HCCO(2227)', 'OCHCO(3676)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C][O](6861)', '[O][C]=C[O](9592)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=C[O](6859)', '[O][C]=C[O](9592)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=[C]OC=[C][O](9625)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[O]C=[C]O[C][C]=O(13038)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[O]C=C1OC1[C]=O(13039)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['O=[C][CH]OC1=CO1(13040)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['O=[C][CH]OC=C=O(11370)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[O]C#CO[CH][C]=O(12992)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[C]=C[O](6859)', 'OCHCO(3676)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.86749e+08,'m^3/(mol*s)'), n=-0.730906, Ea=(37.2121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;Y_1centerbirad] + [Od_R;YJ] for rate rule [Od_R;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C][O](6861)', 'OCHCO(3676)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(29.883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;Y_1centerbirad]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['O(4)', '[CH]=[C]OC=C=O(9621)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.30864e+06,'m^3/(mol*s)'), n=-0.19959, Ea=(22.3126,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_atom_triplet] + [Ct_Ct;YJ] for rate rule [Ct_Ct;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', 'O=[C][C]O[CH][C]=O(12989)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['O=[C][CH]O[C]=[C]O(12993)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.25466e+06,'s^-1'), n=1.80084, Ea=(158.227,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R2H_S;Y_rad_out;O_H_out] + [R2H_S;Cd_rad_out;XH_out] for rate rule [R2H_S;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[O][C]=CO[CH][C]=O(11385)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][C]=[C]OC[C]=O(11810)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.65304e+08,'s^-1'), n=1.30038, Ea=(212.995,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_1;Cd_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][C]=[C]O[CH]C=O(12994)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5Hall;Cd_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[O]C=[C]OC=C=O(11379)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[O]C=C1OC=C1[O](13041)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[O]C#CO[C]=C[O](13042)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]C=[C]O[C]=C[O](13043)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.86881e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O]C=[C]O[C]=[C]O(13044)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O][C]=[C]O[C]=CO(13045)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O][C]=[C]OC=[C]O(12996)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(15567.1,'s^-1'), n=2.15754, Ea=(114.223,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;O_H_out] + [RnH;Cd_rad_out;XH_out] + [R6Hall;Y_rad_out;XH_out] for rate rule [R6Hall;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[O]C1=CO[C][CH]O1(13046)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.39512e+11,'s^-1'), n=0.277081, Ea=(58.6284,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R6plus;carbonyl_intra_H;radadd_intra] + [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[C]1[CH]OO[C]=CO1(13047)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(438.555,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra_H;radadd_intra] for rate rule [R7_linear;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 431.9 to 438.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[O]C1=CO[C]C1[O](13048)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3.76201e+09,'s^-1'), n=0.626373, Ea=(94.9363,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra_H;radadd_intra]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[O]C1[C]OC=[C]O1(13049)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(174.034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;multiplebond_intra;radadd_intra_O] for rate rule [R7;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 170.7 to 174.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[O][C]C1O[C]=CO1(12999)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.25666e+09,'s^-1'), n=0.571459, Ea=(168.098,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cs] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 164.0 to 168.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O][C]=CO[C]=C[O](11389)'],
    products = ['[O]C1[C]OC1[C]=O(13050)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.32992e+07,'s^-1'), n=1.25825, Ea=(196.97,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;multiplebond_intra;radadd_intra_cs] for rate rule [R5;carbonylbond_intra_H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 193.4 to 197.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O][CH][C]O[C]=C[O](13051)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.48818e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction32',
    reactants = ['HCCO(2227)', '[O][C][CH][O](10223)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]C[C]O[C]=C=O(13052)'],
    products = ['[O][C]=CO[C]=C[O](11389)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.26683e+09,'s^-1'), n=1.04717, Ea=(150.072,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;Cd_rad_out_double;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3072',
    isomers = [
        '[O][C]=CO[C]=C[O](11389)',
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
    label = 'PDepNetwork #3072',
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


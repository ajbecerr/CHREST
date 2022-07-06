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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.971109,0.0589565,-7.59981e-05,4.40442e-08,-9.47793e-12,70275.6,25.447], Tmin=(100,'K'), Tmax=(1245.33,'K')), NASAPolynomial(coeffs=[18.0697,-0.000580964,1.27573e-06,-3.00035e-10,2.17823e-14,66374.9,-59.3578], Tmin=(1245.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(583.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CsCJ=O) + radical(Cds_P)"""),
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
    label = '[C]=[CH](9646)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([554.803,1738.79,3454.47],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97621,0.00212915,-8.08978e-08,-3.83305e-10,9.76908e-14,101881,6.00119], Tmin=(100,'K'), Tmax=(1982.31,'K')), NASAPolynomial(coeffs=[5.05695,0.00131032,-4.91873e-07,1.01502e-10,-7.16167e-15,101185,-0.627232], Tmin=(1982.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
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
    label = '[CH]=[C]O[C][C]=O(10432)',
    structure = SMILES('[CH]=[C]O[C][C]=O'),
    E0 = (850.404,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.968796,'amu*angstrom^2'), symmetry=1, barrier=(22.2745,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.969483,'amu*angstrom^2'), symmetry=1, barrier=(22.2903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.968498,'amu*angstrom^2'), symmetry=1, barrier=(22.2677,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0495,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17133,0.0549434,-7.21773e-05,4.17317e-08,-8.93828e-12,102388,23.9881], Tmin=(100,'K'), Tmax=(1238.36,'K')), NASAPolynomial(coeffs=[17.8227,-0.00241191,1.62012e-06,-3.24786e-10,2.20196e-14,98537.7,-58.8], Tmin=(1238.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(850.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CH2_triplet) + radical(C=CJO) + radical(Cds_P) + radical(CsCJ=O)"""),
)

species(
    label = '[C]=[C]O[CH][C]=O(10433)',
    structure = SMILES('[C]=[C]O[CH][C]=O'),
    E0 = (894.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,1685,370,270.4,271.087,272.154],'cm^-1')),
        HinderedRotor(inertia=(0.36062,'amu*angstrom^2'), symmetry=1, barrier=(18.9673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.361587,'amu*angstrom^2'), symmetry=1, barrier=(18.9556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.363398,'amu*angstrom^2'), symmetry=1, barrier=(18.9658,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0495,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24088,0.0590376,-8.67057e-05,5.71733e-08,-1.41632e-11,107666,24.456], Tmin=(100,'K'), Tmax=(1008.66,'K')), NASAPolynomial(coeffs=[15.898,0.000912576,-2.6665e-07,4.19822e-11,-2.96451e-15,104709,-46.3823], Tmin=(1008.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(894.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(CsCJ=O) + radical(C=CJO) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C1OC1[C]=O(10434)',
    structure = SMILES('[CH]=C1OC1[C]=O'),
    E0 = (308.016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95779,0.0281417,2.06817e-05,-6.46312e-08,3.17451e-11,37135,18.6673], Tmin=(100,'K'), Tmax=(915.893,'K')), NASAPolynomial(coeffs=[19.9948,-0.00660064,5.46889e-06,-1.06857e-09,6.77567e-14,31984.2,-76.8483], Tmin=(915.893,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCCJ=O) + radical(Cds_P)"""),
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
    label = '[CH]=C1OC=C1[O](10435)',
    structure = SMILES('[CH]=C1OC=C1[O]'),
    E0 = (268.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80121,0.0365737,-7.08101e-06,-3.24493e-08,1.97602e-11,32354.6,16.3588], Tmin=(100,'K'), Tmax=(902.425,'K')), NASAPolynomial(coeffs=[17.7131,-0.00183808,3.38094e-06,-7.39334e-10,4.97965e-14,28175,-66.0192], Tmin=(902.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.261,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC#C[O](10436)',
    structure = SMILES('C#CO[C][C]=O'),
    E0 = (572.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2175,525,1855,455,950,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.29707,'amu*angstrom^2'), symmetry=1, barrier=(29.8222,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29514,'amu*angstrom^2'), symmetry=1, barrier=(29.7778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.29758,'amu*angstrom^2'), symmetry=1, barrier=(29.8338,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0495,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.829673,0.0576808,-7.28651e-05,3.95516e-08,-7.83774e-12,69027.2,19.1422], Tmin=(100,'K'), Tmax=(1401.86,'K')), NASAPolynomial(coeffs=[20.1513,-0.00536831,3.06982e-06,-5.88827e-10,3.91586e-14,64388,-77.8255], Tmin=(1401.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=[C]O[C]=C[O](10437)',
    structure = SMILES('[CH]=[C]O[C]=C[O]'),
    E0 = (629.068,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,312.371,312.372,312.373],'cm^-1')),
        HinderedRotor(inertia=(0.289456,'amu*angstrom^2'), symmetry=1, barrier=(20.0421,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289449,'amu*angstrom^2'), symmetry=1, barrier=(20.0421,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5471,0.0477086,-5.57502e-05,3.11279e-08,-6.56324e-12,75753.1,27.7724], Tmin=(100,'K'), Tmax=(1259.92,'K')), NASAPolynomial(coeffs=[14.0769,0.00475074,-8.22687e-07,6.17237e-11,-1.63075e-15,72848.1,-34.5703], Tmin=(1259.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(629.068,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CO[C]=[C][O](10438)',
    structure = SMILES('[CH]=CO[C][C]=O'),
    E0 = (610.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3120,650,792.5,1650,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.14737,'amu*angstrom^2'), symmetry=1, barrier=(26.3803,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1453,'amu*angstrom^2'), symmetry=1, barrier=(26.3326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.14527,'amu*angstrom^2'), symmetry=1, barrier=(26.332,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.273749,0.0610795,-7.20136e-05,3.65379e-08,-6.65895e-12,73597.5,24.096], Tmin=(100,'K'), Tmax=(1603.74,'K')), NASAPolynomial(coeffs=[21.1869,-0.00574157,4.19675e-06,-8.42163e-10,5.65971e-14,68774.9,-80.7984], Tmin=(1603.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(610.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CH2_triplet) + radical(Cds_P) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=[C]O[C]=[C]O(10439)',
    structure = SMILES('[CH]=[C]O[C]=[C]O'),
    E0 = (727.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1670,1685,1700,300,370,440,3120,650,792.5,1650,309.624,309.642],'cm^-1')),
        HinderedRotor(inertia=(0.202906,'amu*angstrom^2'), symmetry=1, barrier=(13.8017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202908,'amu*angstrom^2'), symmetry=1, barrier=(13.802,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202977,'amu*angstrom^2'), symmetry=1, barrier=(13.8023,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48899,0.0547555,-8.02181e-05,5.62488e-08,-1.49098e-11,87570.8,29.6991], Tmin=(100,'K'), Tmax=(1007.59,'K')), NASAPolynomial(coeffs=[12.7887,0.00586517,-1.43262e-06,1.49371e-10,-5.23987e-15,85498.4,-23.8852], Tmin=(1007.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]O[C][C]=O(9884)',
    structure = SMILES('C=[C]O[C][C]=O'),
    E0 = (603.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1855,455,950,232.793,232.953,233.157,233.202],'cm^-1')),
        HinderedRotor(inertia=(0.566508,'amu*angstrom^2'), symmetry=1, barrier=(21.7682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.570022,'amu*angstrom^2'), symmetry=1, barrier=(21.7698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.566324,'amu*angstrom^2'), symmetry=1, barrier=(21.7679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06964,0.0534109,-6.21247e-05,3.24307e-08,-6.26796e-12,72676.1,23.9551], Tmin=(100,'K'), Tmax=(1416.95,'K')), NASAPolynomial(coeffs=[17.9725,-0.000249644,9.73841e-07,-2.17551e-10,1.47298e-14,68482.7,-61.3763], Tmin=(1416.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH][C]O[C]=C[O](10440)',
    structure = SMILES('[CH][C]O[C]=C[O]'),
    E0 = (896.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0572642,0.0666856,-8.27962e-05,4.40966e-08,-8.41233e-12,107936,26.9733], Tmin=(100,'K'), Tmax=(1533.6,'K')), NASAPolynomial(coeffs=[22.2539,-0.00735498,5.41484e-06,-1.1145e-09,7.68657e-14,103027,-83.4136], Tmin=(1533.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(896.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(CH2_triplet) + radical(CCJ2_triplet)"""),
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
    label = '[O][C]C1C=[C]O1(10427)',
    structure = SMILES('[O][C]C1C=[C]O1'),
    E0 = (729.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94669,0.0335553,-4.50039e-06,-2.77319e-08,1.52421e-11,87776.5,20.5722], Tmin=(100,'K'), Tmax=(959.672,'K')), NASAPolynomial(coeffs=[16.1492,0.00264288,-3.932e-07,1.26829e-10,-1.58712e-14,83748.1,-54.1479], Tmin=(959.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(729.113,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(CCOJ) + radical(C=CJO) + radical(CH2_triplet)"""),
)

species(
    label = '[C]1[CH]O[C]=CO1(10441)',
    structure = SMILES('[C]1[CH]O[C]=CO1'),
    E0 = (611.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.60798,0.0605469,-7.2204e-05,3.76997e-08,-6.73545e-12,73736.5,20.0054], Tmin=(100,'K'), Tmax=(1792.54,'K')), NASAPolynomial(coeffs=[14.3201,-0.00395895,7.87802e-06,-1.79192e-09,1.2631e-13,73396.4,-46.7473], Tmin=(1792.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(611.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(23dihydro14dioxin) + radical(CCsJOC(O)) + radical(CH2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C1OC1[C][O](10442)',
    structure = SMILES('[CH]=C1OC1[C][O]'),
    E0 = (755.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52879,0.0444686,-3.155e-05,-4.291e-09,8.25453e-12,91015.5,18.7734], Tmin=(100,'K'), Tmax=(946.336,'K')), NASAPolynomial(coeffs=[17.8953,-0.000426589,1.1212e-06,-1.91496e-10,8.26191e-15,86830.5,-65.0278], Tmin=(946.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(755.931,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCOJ) + radical(CH2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1O[C][CH]O1(10443)',
    structure = SMILES('[CH][C]1O[C]=CO1'),
    E0 = (563.61,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10515,0.0103022,0.000101101,-1.66456e-07,7.25226e-11,67884.3,20.0532], Tmin=(100,'K'), Tmax=(909.576,'K')), NASAPolynomial(coeffs=[27.3181,-0.0184034,1.29283e-05,-2.50864e-09,1.61999e-13,59898.5,-117.881], Tmin=(909.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(563.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(Cs_P) + radical(C=CJO) + radical(CCJ2_triplet)"""),
)

species(
    label = 'HC2(2881)',
    structure = SMILES('[C]#C'),
    E0 = (556.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([750,770,3400,2100],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01451,0.0139907,-3.08144e-05,3.10836e-08,-1.10946e-11,66983,5.75946], Tmin=(100,'K'), Tmax=(918.723,'K')), NASAPolynomial(coeffs=[3.14385,0.00498487,-2.32624e-06,4.08757e-10,-2.55787e-14,67315.6,7.08554], Tmin=(918.723,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(556.809,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(62.3585,'J/(mol*K)'), label="""HC2""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[C]=[C]OC[C]=O(9925)',
    structure = SMILES('[C]=[C]OC[C]=O'),
    E0 = (700.421,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1685,370,384.623,384.623,384.624],'cm^-1')),
        HinderedRotor(inertia=(0.172935,'amu*angstrom^2'), symmetry=1, barrier=(18.1544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172935,'amu*angstrom^2'), symmetry=1, barrier=(18.1543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172935,'amu*angstrom^2'), symmetry=1, barrier=(18.1543,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61454,0.0495602,-5.95011e-05,3.40185e-08,-7.4649e-12,84329.6,23.6415], Tmin=(100,'K'), Tmax=(1126.59,'K')), NASAPolynomial(coeffs=[13.5755,0.00709191,-2.95599e-06,5.57203e-10,-3.94779e-14,81634.7,-35.4888], Tmin=(1126.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.421,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCJ=O) + radical(C=CJO) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[C]#CO[CH][C]O(10444)',
    structure = SMILES('[C]=[C]OC=[C]O'),
    E0 = (798.612,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,1670,1700,300,440,203.22,203.229,203.682],'cm^-1')),
        HinderedRotor(inertia=(0.597261,'amu*angstrom^2'), symmetry=1, barrier=(17.5859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.595472,'amu*angstrom^2'), symmetry=1, barrier=(17.5836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.597914,'amu*angstrom^2'), symmetry=1, barrier=(17.5823,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38485,0.0549688,-7.03991e-05,3.83768e-08,-6.44359e-12,96147.5,26.9239], Tmin=(100,'K'), Tmax=(859.35,'K')), NASAPolynomial(coeffs=[15.286,0.00234043,3.83094e-07,-1.80402e-10,1.56506e-14,93312.3,-40.6289], Tmin=(859.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(798.612,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO) + radical(CdCdJ2_triplet)"""),
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
    E0 = (583.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (905.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1106.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1062.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1106.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (586.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (808.76,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (627.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (583.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (591.627,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (799.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (833.206,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (764.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (877.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (883.354,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (918.412,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (1050.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (729.113,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (614.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (755.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (709.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (761.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (1233.11,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (778.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (889.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]OC=[C][O](9625)'],
    products = ['HCCO(2227)', 'HCCO(2227)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C][O](6861)', '[CH]=[C][O](6861)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=[CH](9646)', '[O][C]=C[O](9592)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]=[C]O[C][C]=O(10432)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[C]=[C]O[CH][C]=O(10433)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]OC=[C][O](9625)'],
    products = ['[CH]=C1OC1[C]=O(10434)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]=[CH](9646)', 'OCHCO(3676)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.86749e+08,'m^3/(mol*s)'), n=-0.730906, Ea=(37.2121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;Y_1centerbirad] + [Od_R;YJ] for rate rule [Od_R;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['HCCO(2227)', '[CH]=[C][O](6861)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.86749e+08,'m^3/(mol*s)'), n=-0.730906, Ea=(37.2121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;Y_1centerbirad] + [Od_R;YJ] for rate rule [Od_R;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]OC=[C][O](9625)'],
    products = ['[CH]=[C]OC=C=O(9621)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]OC=[C][O](9625)'],
    products = ['[CH]=C1OC=C1[O](10435)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH]=[C]OC#C[O](10436)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]O[C]=C[O](10437)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=CO[C]=[C][O](10438)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.38739e+08,'s^-1'), n=1.16185, Ea=(153.379,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Cd_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]O[C]=[C]O(10439)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C]O[C][C]=O(9884)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.19427e+11,'s^-1'), n=0.834916, Ea=(280.046,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;Cd_H_out_singleH] for rate rule [R4HJ_2;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH][C]O[C]=C[O](10440)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['HCCO(2227)', '[CH][C][O](10218)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]OC=[C][O](9625)'],
    products = ['[O][C]C1C=[C]O1(10427)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.89861e+07,'s^-1'), n=1.13751, Ea=(145.771,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic
Ea raised from 142.7 to 145.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]OC=[C][O](9625)'],
    products = ['[C]1[CH]O[C]=CO1(10441)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.23539e+09,'s^-1'), n=0.445806, Ea=(31.0324,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;triplebond_intra_H;radadd_intra] for rate rule [R6_linear;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]OC=[C][O](9625)'],
    products = ['[CH]=C1OC1[C][O](10442)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.98513e+09,'s^-1'), n=0.768, Ea=(172.588,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_T;triplebond_intra_H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 171.0 to 172.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]OC=[C][O](9625)'],
    products = ['[CH]=C1O[C][CH]O1(10443)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.03297e+11,'s^-1'), n=0.45637, Ea=(125.756,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;triplebond_intra_H;radadd_intra] for rate rule [R6;triplebond_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HC2(2881)', '[O][C]=C[O](9592)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(5e+12,'cm^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1300,'K'), comment="""Estimated using template [R_R;CtJ_Ct] for rate rule [Od_R;CtJ_Ct]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['HC2(2881)', '[O][C][CH][O](10223)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[C]=[C]OC[C]=O(9925)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(20941.3,'s^-1'), n=2.29593, Ea=(78.4716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSS;Y_rad_out;XH_out] for rate rule [R4H_TSS;Ct_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[C]#CO[CH][C]O(10444)'],
    products = ['[CH]=[C]OC=[C][O](9625)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(342.007,'s^-1'), n=2.46341, Ea=(91.1577,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;O_H_out] + [R6Hall;Y_rad_out;XH_out] for rate rule [R6Hall;Ct_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2547',
    isomers = [
        '[CH]=[C]OC=[C][O](9625)',
    ],
    reactants = [
        ('HCCO(2227)', 'HCCO(2227)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2547',
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


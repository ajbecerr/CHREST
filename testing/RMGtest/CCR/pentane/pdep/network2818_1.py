species(
    label = 'O=[C][CH]OC[C]=O(11150)',
    structure = SMILES('[O][C]=COC[C]=O'),
    E0 = (74.9892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1855,455,950,1685,370,205.135,205.135,205.135,205.135],'cm^-1')),
        HinderedRotor(inertia=(0.813479,'amu*angstrom^2'), symmetry=1, barrier=(24.2913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.81348,'amu*angstrom^2'), symmetry=1, barrier=(24.2913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813479,'amu*angstrom^2'), symmetry=1, barrier=(24.2913,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995855,0.0523372,-2.4758e-05,-2.24589e-08,1.707e-11,9139.99,26.753], Tmin=(100,'K'), Tmax=(927.568,'K')), NASAPolynomial(coeffs=[21.2354,0.000482909,1.80938e-06,-3.79382e-10,2.17064e-14,3861.31,-77.5841], Tmin=(927.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.9892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CsCJ=O) + radical(C=CJO)"""),
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
    label = 'O=[C]COC=C=O(11167)',
    structure = SMILES('O=[C]COC=C=O'),
    E0 = (-124.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.336904,0.0627665,-6.67288e-05,3.25598e-08,-5.8629e-12,-14855.8,25.3974], Tmin=(100,'K'), Tmax=(1572.81,'K')), NASAPolynomial(coeffs=[19.952,0.00247767,6.90932e-07,-2.22907e-10,1.64316e-14,-19739.2,-74.0262], Tmin=(1572.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsCJ=O)"""),
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
    label = '[CH2]OC=[C][O](11465)',
    structure = SMILES('[CH2]OC=[C][O]'),
    E0 = (216.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,400.773,400.953,401.443],'cm^-1')),
        HinderedRotor(inertia=(0.188985,'amu*angstrom^2'), symmetry=1, barrier=(21.6878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189931,'amu*angstrom^2'), symmetry=1, barrier=(21.685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.09104,0.028776,1.00821e-05,-4.80262e-08,2.47341e-11,26115.5,20.59], Tmin=(100,'K'), Tmax=(909.214,'K')), NASAPolynomial(coeffs=[17.1662,-0.00241134,3.5706e-06,-7.50851e-10,4.89516e-14,21922,-58.6903], Tmin=(909.214,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-OsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COCJ) + radical(C=CJO)"""),
)

species(
    label = '[O]C1=COCC1=O(11933)',
    structure = SMILES('O=C1[CH]OCC1=O'),
    E0 = (-192.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93711,0.0293701,2.75421e-05,-6.37864e-08,2.80896e-11,-23016.6,19.5341], Tmin=(100,'K'), Tmax=(957.138,'K')), NASAPolynomial(coeffs=[16.4047,0.00869168,-2.39902e-06,4.94705e-10,-4.30968e-14,-27608.4,-59.1497], Tmin=(957.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + ring(Cyclopentane) + radical(CCsJOCs)"""),
)

species(
    label = '[O]C=COC=C=O(11302)',
    structure = SMILES('[O]C=COC=C=O'),
    E0 = (-124.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01811,0.0513189,-1.99387e-05,-2.98534e-08,2.06899e-11,-14893.2,23.6764], Tmin=(100,'K'), Tmax=(911.124,'K')), NASAPolynomial(coeffs=[21.6465,-0.000988181,3.19519e-06,-6.98075e-10,4.56575e-14,-20240.1,-82.6372], Tmin=(911.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-124.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=COJ)"""),
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
    label = 'O=[C]CO[C]=C=O(11809)',
    structure = SMILES('O=[C]CO[C]=C=O'),
    E0 = (115.002,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,1685,370,2120,512.5,787.5,344.99,345.016,345.124],'cm^-1')),
        HinderedRotor(inertia=(0.216148,'amu*angstrom^2'), symmetry=1, barrier=(18.2791,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216327,'amu*angstrom^2'), symmetry=1, barrier=(18.2765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.216477,'amu*angstrom^2'), symmetry=1, barrier=(18.2791,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23019,0.0566721,-6.70204e-05,3.79281e-08,-8.22628e-12,13935,25.306], Tmin=(100,'K'), Tmax=(1141.43,'K')), NASAPolynomial(coeffs=[15.1303,0.0079606,-3.00647e-06,5.39827e-10,-3.73544e-14,10761.8,-43.593], Tmin=(1141.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CsCJ=O) + radical(C=CJO)"""),
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
    label = '[O]C[C]=O(2320)',
    structure = SMILES('[O]C[C]=O'),
    E0 = (65.8932,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,298.422],'cm^-1')),
        HinderedRotor(inertia=(0.460283,'amu*angstrom^2'), symmetry=1, barrier=(28.9433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01454,0.0123699,2.06297e-05,-3.93769e-08,1.66626e-11,7969.02,14.0082], Tmin=(100,'K'), Tmax=(965.917,'K')), NASAPolynomial(coeffs=[10.932,0.0028129,-6.04041e-07,1.77093e-10,-1.91371e-14,5355.8,-29.5241], Tmin=(965.917,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(65.8932,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=OCOJ) + radical(CsCJ=O)"""),
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
    label = '[O]C=CO[CH][C]=O(11291)',
    structure = SMILES('[O]C=CO[CH][C]=O'),
    E0 = (29.1724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.416979,0.0599772,-2.46439e-05,-3.8972e-08,2.71516e-11,3655.28,25.1767], Tmin=(100,'K'), Tmax=(910.45,'K')), NASAPolynomial(coeffs=[27.8973,-0.0101551,7.53526e-06,-1.49074e-09,9.7606e-14,-3445.81,-116.339], Tmin=(910.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.1724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(CCsJOC(O)) + radical(CsCJ=O)"""),
)

species(
    label = '[O]C=[C]OC[C]=O(11163)',
    structure = SMILES('[O]C=[C]OC[C]=O'),
    E0 = (74.9892,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,204.976,205.081,205.127,205.223],'cm^-1')),
        HinderedRotor(inertia=(0.813285,'amu*angstrom^2'), symmetry=1, barrier=(24.2916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813077,'amu*angstrom^2'), symmetry=1, barrier=(24.2911,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.813048,'amu*angstrom^2'), symmetry=1, barrier=(24.2914,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4069.61,'J/mol'), sigma=(6.41546,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=635.66 K, Pc=34.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.995888,0.0523368,-2.47565e-05,-2.24611e-08,1.70709e-11,9139.99,26.7528], Tmin=(100,'K'), Tmax=(927.563,'K')), NASAPolynomial(coeffs=[21.2353,0.000483082,1.80928e-06,-3.79357e-10,2.17044e-14,3861.35,-77.5836], Tmin=(927.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(74.9892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = 'O=[C]CO[C]=[C]O(11811)',
    structure = SMILES('O=[C]CO[C]=[C]O'),
    E0 = (173.271,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,1855,455,950,1670,1700,300,440,240.402,240.402,240.402],'cm^-1')),
        HinderedRotor(inertia=(0.470619,'amu*angstrom^2'), symmetry=1, barrier=(19.3006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470619,'amu*angstrom^2'), symmetry=1, barrier=(19.3006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470619,'amu*angstrom^2'), symmetry=1, barrier=(19.3006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470619,'amu*angstrom^2'), symmetry=1, barrier=(19.3006,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.184815,0.0682351,-8.00796e-05,4.25954e-08,-8.30261e-12,20990.5,31.3829], Tmin=(100,'K'), Tmax=(1456.11,'K')), NASAPolynomial(coeffs=[20.9,-0.000203997,2.30349e-06,-5.62685e-10,4.1154e-14,16180.5,-72.1415], Tmin=(1456.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.271,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(CsCJ=O) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[O][C]=[C]OCC=O(11812)',
    structure = SMILES('[O][C]=[C]OCC=O'),
    E0 = (160.062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2782.5,750,1395,475,1775,1000,1670,1700,300,440,395.704,395.705,395.705,395.705],'cm^-1')),
        HinderedRotor(inertia=(0.17156,'amu*angstrom^2'), symmetry=1, barrier=(19.0627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17156,'amu*angstrom^2'), symmetry=1, barrier=(19.0627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17156,'amu*angstrom^2'), symmetry=1, barrier=(19.0627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621074,0.0598731,-6.28993e-05,3.08781e-08,-5.65869e-12,19385.1,30.5579], Tmin=(100,'K'), Tmax=(1510.64,'K')), NASAPolynomial(coeffs=[18.6672,0.00451916,-4.18762e-07,-1.25467e-11,2.45837e-15,14796.5,-61.0902], Tmin=(1510.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'O=[C][CH]OC=[C]O(11944)',
    structure = SMILES('O=[C][CH]OC=[C]O'),
    E0 = (127.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.01501,'amu*angstrom^2'), symmetry=1, barrier=(23.3371,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01418,'amu*angstrom^2'), symmetry=1, barrier=(23.318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01397,'amu*angstrom^2'), symmetry=1, barrier=(23.3131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01489,'amu*angstrom^2'), symmetry=1, barrier=(23.3343,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.913834,0.0819787,-0.000101295,5.39176e-08,-1.0245e-11,15528.6,31.6743], Tmin=(100,'K'), Tmax=(1553.21,'K')), NASAPolynomial(coeffs=[25.9597,-0.00869698,7.00654e-06,-1.46653e-09,1.0184e-13,9770.07,-101.472], Tmin=(1553.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = 'O=[C]C1OCC1=O(11152)',
    structure = SMILES('O=[C]C1OCC1=O'),
    E0 = (-146.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86718,0.0377478,-1.08622e-05,-1.41965e-08,8.0611e-12,-17516.7,20.8996], Tmin=(100,'K'), Tmax=(1035.66,'K')), NASAPolynomial(coeffs=[12.828,0.0150566,-6.44615e-06,1.2736e-09,-9.37937e-14,-20840.5,-37.4496], Tmin=(1035.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.347,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(C=OCCJ=O)"""),
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
    E0 = (74.9892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (74.9892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (710.67,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (82.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (99.9624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (122.036,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (291.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (341.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (203.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (264.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (365.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (490.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (475.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (526.525,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (208.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (244.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (323.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (290.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (179.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (83.2735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (116.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['O=[C][CH]OC[C]=O(11150)'],
    products = ['OCHCO(3676)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O=[C][CH]OC[C]=O(11150)'],
    products = ['O=[C]COC=C=O(11167)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(2355)', '[CH2]OC=[C][O](11465)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O=[C][CH]OC[C]=O(11150)'],
    products = ['[O]C1=COCC1=O(11933)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O=[C][CH]OC[C]=O(11150)'],
    products = ['[O]C=COC=C=O(11302)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CO(2039)', '[CH2]OC=[C][O](11465)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(754.069,'m^3/(mol*s)'), n=1.0822, Ea=(24.7919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [COm;C_pri_rad] for rate rule [COm;C_rad/H2/O]
Euclidian distance = 1.0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'O=[C][CH]OC=C=O(11370)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.64e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ck;HJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'O=[C]CO[C]=C=O(11809)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C]=C[O](9592)', 'ketene(1375)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.3e+11,'cm^3/(mol*s)'), n=0, Ea=(60.12,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd_R;O_rad/OneDe] for rate rule [Cds-HH_Ck;O_rad/OneDe]
Euclidian distance = 3.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['HCCO(2227)', '[O]C[C]=O(2320)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(33.0536,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;OJ-Cs]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=C[O](9592)', '[CH2][C]=O(1376)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C][O](6861)', '[O]C[C]=O(2320)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.505e+06,'m^3/(mol*s)'), n=1.79841e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_4R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[O][C]=CO[CH][C]=O(11385)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.82057e+14,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[O][C]=[C]OC[C]=O(11810)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['O=[C][CH]OC[C]=O(11150)'],
    products = ['[O]C=CO[CH][C]=O(11291)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.00141351,'s^-1'), n=4.515, Ea=(133.574,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=[C][CH]OC[C]=O(11150)'],
    products = ['[O]C=[C]OC[C]=O(11163)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.6462e+08,'s^-1'), n=1.3775, Ea=(169.747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['O=[C]CO[C]=[C]O(11811)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][C]=[C]OCC=O(11812)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.82783e+06,'s^-1'), n=1.57155, Ea=(130.124,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS_OCs;Y_rad_out;XH_out] for rate rule [R4H_SSS_OCs;Cd_rad_out;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O=[C][CH]OC=[C]O(11944)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.10205e+06,'s^-1'), n=1.54368, Ea=(52.1315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;XH_out] for rate rule [R5HJ_3;Y_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O=[C][CH]OC[C]=O(11150)'],
    products = ['O=[C]C1OCC1=O(11152)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['OCHCO(3676)', '[CH2][C]=O(1376)'],
    products = ['O=[C][CH]OC[C]=O(11150)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.41198,'m^3/(mol*s)'), n=2.03228, Ea=(31.3669,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CJ] + [Od_R;YJ] for rate rule [Od_R;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

network(
    label = 'PDepNetwork #2818',
    isomers = [
        'O=[C][CH]OC[C]=O(11150)',
    ],
    reactants = [
        ('OCHCO(3676)', 'ketene(1375)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2818',
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


species(
    label = '[CH]=[C]OO[CH][C]=O(11367)',
    structure = SMILES('[CH]=[C]OO[CH][C]=O'),
    E0 = (633.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,1685,370,1855,455,950,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(1.60186,'amu*angstrom^2'), symmetry=1, barrier=(36.83,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.248961,'amu*angstrom^2'), symmetry=1, barrier=(5.72411,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369899,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.98762,'amu*angstrom^2'), symmetry=1, barrier=(45.6993,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29054,0.0658579,-0.000108481,9.61166e-08,-3.36909e-11,76300.2,29.2478], Tmin=(100,'K'), Tmax=(798.808,'K')), NASAPolynomial(coeffs=[8.09903,0.0232021,-1.23031e-05,2.4299e-09,-1.70439e-13,75485.7,-0.359831], Tmin=(798.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(633.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(CsCJ=O) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]O[O](9590)',
    structure = SMILES('[CH]=[C]O[O]'),
    E0 = (584.749,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1685,370,3120,650,792.5,1650,1740.97,3927],'cm^-1')),
        HinderedRotor(inertia=(0.523695,'amu*angstrom^2'), symmetry=1, barrier=(12.0408,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.07482,0.0276588,-6.25889e-05,6.64728e-08,-2.48998e-11,70355.5,17.2575], Tmin=(100,'K'), Tmax=(889.202,'K')), NASAPolynomial(coeffs=[2.44582,0.0123284,-6.09377e-06,1.14894e-09,-7.66206e-14,71185.3,24.2551], Tmin=(889.202,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(584.749,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]OO[C]=[C][O](13112)',
    structure = SMILES('[CH]=[C]OO[C][C]=O'),
    E0 = (956.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1685,370,3120,650,792.5,1650,1855,455,950,221.397],'cm^-1')),
        HinderedRotor(inertia=(0.950122,'amu*angstrom^2'), symmetry=1, barrier=(27.5328,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.626497,'amu*angstrom^2'), symmetry=1, barrier=(18.1363,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.819257,'amu*angstrom^2'), symmetry=1, barrier=(27.5467,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.322803,'amu*angstrom^2'), symmetry=1, barrier=(9.55319,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15115,0.0681532,-0.000114445,9.56385e-08,-3.14292e-11,115094,28.7279], Tmin=(100,'K'), Tmax=(749.175,'K')), NASAPolynomial(coeffs=[11.011,0.0155104,-9.04578e-06,1.84896e-09,-1.32176e-13,113617,-15.9931], Tmin=(749.175,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(956.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(166.289,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CH2_triplet) + radical(C=CJO) + radical(Cds_P) + radical(CsCJ=O)"""),
)

species(
    label = '[C]=[C]OO[CH][C]=O(13190)',
    structure = SMILES('[C]=[C]OO[CH][C]=O'),
    E0 = (944.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,1855,455,950,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.836,'amu*angstrom^2'), symmetry=1, barrier=(42.2132,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.31009,'amu*angstrom^2'), symmetry=1, barrier=(7.12957,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.83602,'amu*angstrom^2'), symmetry=1, barrier=(42.2138,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0160016,'amu*angstrom^2'), symmetry=1, barrier=(42.1654,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29398,0.0691268,-0.000130532,1.24177e-07,-4.48086e-11,113702,29.2089], Tmin=(100,'K'), Tmax=(840.049,'K')), NASAPolynomial(coeffs=[6.6985,0.0233698,-1.30754e-05,2.58913e-09,-1.79973e-13,113501,8.28219], Tmin=(840.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(944.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(166.289,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(CsCJ=O) + radical(C=CJO) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C1OOC1[C]=O(13191)',
    structure = SMILES('[CH]=C1OOC1[C]=O'),
    E0 = (374.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07856,0.0288505,1.50225e-05,-4.48169e-08,1.98123e-11,45134.7,24.9371], Tmin=(100,'K'), Tmax=(988.872,'K')), NASAPolynomial(coeffs=[14.9772,0.00894203,-3.72339e-06,8.17992e-10,-6.68362e-14,41006.1,-45.1239], Tmin=(988.872,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCCJ=O) + radical(Cds_P)"""),
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
    label = '[CH]=[C]OOC=C=O(11377)',
    structure = SMILES('[CH]=[C]OOC=C=O'),
    E0 = (522.669,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,1685,370,2120,512.5,787.5,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.892692,'amu*angstrom^2'), symmetry=1, barrier=(20.5247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.669891,'amu*angstrom^2'), symmetry=1, barrier=(15.4021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.666654,'amu*angstrom^2'), symmetry=1, barrier=(15.3277,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37585,0.0612914,-8.9971e-05,6.79673e-08,-2.02797e-11,62953.8,27.6497], Tmin=(100,'K'), Tmax=(823.414,'K')), NASAPolynomial(coeffs=[10.6471,0.0162506,-7.91564e-06,1.52814e-09,-1.06639e-13,61427.1,-15.2768], Tmin=(823.414,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.669,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1OOC=C1[O](13163)',
    structure = SMILES('[CH]=C1OOC=C1[O]'),
    E0 = (265.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88612,0.0305624,2.16946e-05,-6.36647e-08,3.05537e-11,32071.8,20.5037], Tmin=(100,'K'), Tmax=(919.404,'K')), NASAPolynomial(coeffs=[18.2247,0.000972463,2.27412e-06,-4.95876e-10,2.96184e-14,27313.7,-66.4848], Tmin=(919.404,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.904,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(207.862,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=C(C)OJ) + radical(Cds_P)"""),
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
    label = '[CH]=[C]OO[C]=C=O(13111)',
    structure = SMILES('[CH]=[C]OO[C]=C=O'),
    E0 = (762.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,1670,1700,300,440,3120,650,792.5,1650,2120,512.5,787.5],'cm^-1')),
        HinderedRotor(inertia=(0.277274,'amu*angstrom^2'), symmetry=1, barrier=(6.37508,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.278784,'amu*angstrom^2'), symmetry=1, barrier=(6.40979,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.759169,'amu*angstrom^2'), symmetry=1, barrier=(17.4548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.0489,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48369,0.0643668,-0.000122035,1.1437e-07,-4.01895e-11,91779.1,30.3838], Tmin=(100,'K'), Tmax=(867.455,'K')), NASAPolynomial(coeffs=[6.79049,0.0198693,-1.04596e-05,2.00642e-09,-1.36154e-13,91611.9,9.87921], Tmin=(867.455,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(762.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + missing(O2d-Cdd) + group(Cds-CdsOsH) + group(Cds-(Cdd-O2d)OsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OO[C]=C[O](11376)',
    structure = SMILES('[CH]=[C]OO[C]=C[O]'),
    E0 = (722.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,247.747],'cm^-1')),
        HinderedRotor(inertia=(0.356526,'amu*angstrom^2'), symmetry=1, barrier=(15.6361,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.561425,'amu*angstrom^2'), symmetry=1, barrier=(24.2541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354825,'amu*angstrom^2'), symmetry=1, barrier=(15.6453,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3860.85,'J/mol'), sigma=(6.18146,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=603.06 K, Pc=37.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33121,0.0589737,-7.55537e-05,4.75801e-08,-1.16579e-11,86980.6,31.5425], Tmin=(100,'K'), Tmax=(1004.35,'K')), NASAPolynomial(coeffs=[13.1186,0.0120283,-5.44042e-06,1.04022e-09,-7.32581e-14,84612.9,-25.3759], Tmin=(1004.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OO[C]=[C]O(13113)',
    structure = SMILES('[CH]=[C]OO[C]=[C]O'),
    E0 = (820.681,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3615,1277.5,1000,1670,1685,1700,300,370,440,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.432652,'amu*angstrom^2'), symmetry=1, barrier=(11.6723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.432652,'amu*angstrom^2'), symmetry=1, barrier=(11.6723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.432647,'amu*angstrom^2'), symmetry=1, barrier=(11.6723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771712,'amu*angstrom^2'), symmetry=1, barrier=(20.8196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00797,0.0691888,-0.000111303,8.76466e-08,-2.65374e-11,98809.7,34.4176], Tmin=(100,'K'), Tmax=(873.974,'K')), NASAPolynomial(coeffs=[12.536,0.0119446,-5.36068e-06,9.65258e-10,-6.3402e-14,96965.9,-18.6661], Tmin=(873.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(820.681,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(C=CJO) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=COO[C]=[C][O](13114)',
    structure = SMILES('[CH]=COO[C][C]=O'),
    E0 = (716.394,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.06076,'amu*angstrom^2'), symmetry=1, barrier=(24.3889,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06071,'amu*angstrom^2'), symmetry=1, barrier=(24.3877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06045,'amu*angstrom^2'), symmetry=1, barrier=(24.3818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06082,'amu*angstrom^2'), symmetry=1, barrier=(24.3904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0536,0.0650183,-8.23705e-05,4.92497e-08,-1.14021e-11,86268.4,25.9521], Tmin=(100,'K'), Tmax=(1061.31,'K')), NASAPolynomial(coeffs=[15.4035,0.0109335,-5.92861e-06,1.23153e-09,-9.0879e-14,83222.5,-44.1316], Tmin=(1061.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(716.394,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CH2_triplet) + radical(Cds_P) + radical(CsCJ=O)"""),
)

species(
    label = 'C=[C]OO[C]=[C][O](11847)',
    structure = SMILES('C=[C]OO[C][C]=O'),
    E0 = (709.042,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,303.974,303.991],'cm^-1')),
        HinderedRotor(inertia=(0.251265,'amu*angstrom^2'), symmetry=1, barrier=(16.5577,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.420217,'amu*angstrom^2'), symmetry=1, barrier=(27.3843,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138714,'amu*angstrom^2'), symmetry=1, barrier=(9.16146,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.415767,'amu*angstrom^2'), symmetry=1, barrier=(27.386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3949,0.0624431,-8.92056e-05,6.56153e-08,-1.93493e-11,85367.4,27.4617], Tmin=(100,'K'), Tmax=(826.556,'K')), NASAPolynomial(coeffs=[10.3558,0.0190807,-1.05186e-05,2.15378e-09,-1.55991e-13,83885.9,-14.0633], Tmin=(826.556,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(709.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CH2_triplet) + radical(CsCJ=O)"""),
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
    E0 = (633.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1009.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1040.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1167.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1156.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (641.919,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (641.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (633.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (640.748,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (722.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (989.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (633.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (782.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (926.537,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (971.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (784.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (818.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]OO[CH][C]=O(11367)'],
    products = ['HCCO(2227)', 'OCHCO(3676)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C][O](6861)', '[CH]=[C]O[O](9590)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=[CH](9646)', '[O][C]=CO[O](9542)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [O_rad/NonDe;Birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]=[C]OO[C]=[C][O](13112)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[C]=[C]OO[CH][C]=O(13190)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=[C]OO[CH][C]=O(11367)'],
    products = ['[CH]=C1OOC1[C]=O(13191)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C][O](6861)', '[O][C]=C[O](9592)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(12.5775,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]OO[CH][C]=O(11367)'],
    products = ['[CH]=[C]OOC=C=O(11377)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]OO[CH][C]=O(11367)'],
    products = ['[CH]=C1OOC=C1[O](13163)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]OO[CH][C]=O(11367)'],
    products = ['[CH]=[C][O](6861)', 'O=C1[CH]O1(10908)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.60331e+13,'s^-1'), n=-0.0568549, Ea=(88.8943,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO;Y_rad_intra;OO_intra] for rate rule [R2OO_D;Y_rad_intra;OO_intra]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH]=[C]OO[C]=C=O(13111)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['HCCO(2227)', '[O][C]=C[O](9592)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.3e+11,'cm^3/(mol*s)'), n=0, Ea=(264.07,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R_R;O_rad/OneDe] for rate rule [Od_R;O_rad/OneDe]
Euclidian distance = 1.0
family: R_Addition_MultipleBond
Ea raised from 260.1 to 264.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['HCCO(2227)', '[CH]=[C]O[O](9590)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(33.0536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;O_rad/NonDe] for rate rule [Ct_Ct;OJ-O2s]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C]OO[C]=C[O](11376)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]OO[C]=[C]O(13113)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=COO[C]=[C][O](13114)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(13083.4,'s^-1'), n=2.28331, Ea=(68.3874,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_RSS;Cd_rad_out;XH_out] + [R4H_SSS;Y_rad_out;XH_out] for rate rule [R4H_SSS;Cd_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=[C]OO[C]=[C][O](11847)'],
    products = ['[CH]=[C]OO[CH][C]=O(11367)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(8820.35,'s^-1'), n=2.27943, Ea=(109.217,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cd_H_out_singleH] for rate rule [R5HJ_3;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3050',
    isomers = [
        '[CH]=[C]OO[CH][C]=O(11367)',
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
    label = 'PDepNetwork #3050',
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


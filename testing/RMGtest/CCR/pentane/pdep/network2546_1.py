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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46965,0.0598019,-8.50278e-05,6.58157e-08,-2.06105e-11,40641.6,21.5237], Tmin=(100,'K'), Tmax=(778.943,'K')), NASAPolynomial(coeffs=[9.00533,0.0211058,-1.05129e-05,2.04269e-09,-1.43112e-13,39467.6,-12.9492], Tmin=(778.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
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
    label = '[CH]C(=[CH])[O](10446)',
    structure = SMILES('[CH]C(=[CH])[O]'),
    E0 = (569.267,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3120,650,792.5,1650,319.325,319.431,319.703],'cm^-1')),
        HinderedRotor(inertia=(0.709667,'amu*angstrom^2'), symmetry=1, barrier=(51.1441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (54.0474,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64262,0.0274443,-2.15887e-05,9.33787e-09,-1.65053e-12,68517.9,15.3589], Tmin=(100,'K'), Tmax=(1341.73,'K')), NASAPolynomial(coeffs=[7.83006,0.0119795,-4.29974e-06,7.47549e-10,-4.99365e-14,67125.8,-11.1924], Tmin=(1341.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(569.267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C1=CC(=O)O1(10447)',
    structure = SMILES('[CH]=C1[CH]C(=O)O1'),
    E0 = (119.222,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.04502,0.0265888,2.83952e-05,-6.76756e-08,3.09442e-11,14424.6,16.5856], Tmin=(100,'K'), Tmax=(933.444,'K')), NASAPolynomial(coeffs=[17.6604,0.00156214,1.29929e-06,-2.48969e-10,9.98031e-15,9684.48,-67.4489], Tmin=(933.444,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(4-Methylene-2-oxetanone) + radical(C=CCJCO) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(O)=C=C=O(10448)',
    structure = SMILES('[CH]C(O)=C=C=O'),
    E0 = (279.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.36592,0.0272406,-2.64473e-07,-2.65854e-08,1.43913e-11,33694.4,3.52878], Tmin=(100,'K'), Tmax=(922.645,'K')), NASAPolynomial(coeffs=[12.773,0.00470649,-3.46158e-07,3.67962e-12,-1.86301e-15,30812.7,-51.0508], Tmin=(922.645,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cdd-CdsCds) + missing(Cdd-CddO2d) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C1([O])[CH]C1=O(10449)',
    structure = SMILES('[CH]C1([O])C=C1[O]'),
    E0 = (645.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72863,0.0414805,-2.78285e-05,-3.29999e-09,6.738e-12,77776.5,19.8398], Tmin=(100,'K'), Tmax=(954.123,'K')), NASAPolynomial(coeffs=[15.5122,0.00400667,-8.47257e-07,1.59411e-10,-1.4581e-14,74221.7,-50.8556], Tmin=(954.123,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(645.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopropene) + radical(CC(C)2OJ) + radical(C=C(C)OJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH]C([O])=C=C=O(10450)',
    structure = SMILES('[CH]C([O])=C=C=O'),
    E0 = (417.4,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,540,610,2055,2120,512.5,787.5,857.887,863.091,2481.31,2488.53],'cm^-1')),
        HinderedRotor(inertia=(0.00556569,'amu*angstrom^2'), symmetry=1, barrier=(2.82909,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0495,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72534,0.023358,-8.16998e-06,-8.85134e-09,6.00857e-12,50251.8,3.73928], Tmin=(100,'K'), Tmax=(936.351,'K')), NASAPolynomial(coeffs=[9.05658,0.00802796,-2.3811e-06,3.90475e-10,-2.68868e-14,48552.5,-29.1316], Tmin=(936.351,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.4,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cdd-CdsCds) + missing(Cdd-CddO2d) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([O])=[C][C]=O(10451)',
    structure = SMILES('[CH]C([O])=[C][C]=O'),
    E0 = (581.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,1685,370,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.03167,'amu*angstrom^2'), symmetry=1, barrier=(46.7121,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.04948,'amu*angstrom^2'), symmetry=1, barrier=(47.1216,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0495,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24231,0.0658409,-0.000109374,9.26905e-08,-3.08325e-11,69990.4,23.038], Tmin=(100,'K'), Tmax=(803.639,'K')), NASAPolynomial(coeffs=[9.98161,0.0172803,-9.28701e-06,1.82456e-09,-1.27176e-13,68749.2,-16.1966], Tmin=(803.639,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(581.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CJC=O) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]C([O])=[C]C=O(10452)',
    structure = SMILES('[CH]C([O])=C=C[O]'),
    E0 = (373.303,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,540,610,2055,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.08969,'amu*angstrom^2'), symmetry=1, barrier=(48.046,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.683756,0.055731,-5.65795e-05,2.72969e-08,-4.87559e-12,45032.2,23.8771], Tmin=(100,'K'), Tmax=(1613.4,'K')), NASAPolynomial(coeffs=[16.5017,0.00578825,-1.74595e-07,-1.3057e-10,1.28797e-14,41324.1,-55.6751], Tmin=(1613.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.303,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(O)=[C][C]=O(10453)',
    structure = SMILES('[CH]C(O)=[C][C]=O'),
    E0 = (443.344,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,1685,370,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.91799,'amu*angstrom^2'), symmetry=1, barrier=(44.0984,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91686,'amu*angstrom^2'), symmetry=1, barrier=(44.0724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.91311,'amu*angstrom^2'), symmetry=1, barrier=(43.9862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05642,0.0675266,-9.30036e-05,6.26691e-08,-1.65226e-11,53425.7,22.2144], Tmin=(100,'K'), Tmax=(931.479,'K')), NASAPolynomial(coeffs=[13.4148,0.0144581,-7.54715e-06,1.5087e-09,-1.08112e-13,51123.3,-36.5306], Tmin=(931.479,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.344,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=O) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
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
    label = '[C]=C([O])[CH][C]=O(10454)',
    structure = SMILES('[C]C([O])=C[C]=O'),
    E0 = (635.981,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.893684,'amu*angstrom^2'), symmetry=1, barrier=(20.5476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (81.0495,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.40607,0.0620355,-0.000105097,8.28276e-08,-2.4409e-11,76579.6,19.1005], Tmin=(100,'K'), Tmax=(689.01,'K')), NASAPolynomial(coeffs=[11.6847,0.0093233,-5.4919e-06,1.11195e-09,-7.84623e-14,74998,-27.8577], Tmin=(689.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(635.981,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CCJ=O) + radical(CJ3)"""),
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
    label = 'O=[C][CH]C1=CO1(10455)',
    structure = SMILES('O=[C]C=C1[CH]O1'),
    E0 = (178.713,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72713,0.0349558,2.32107e-06,-4.18692e-08,2.17799e-11,21589.8,16.8296], Tmin=(100,'K'), Tmax=(951.488,'K')), NASAPolynomial(coeffs=[19.4161,-0.00206219,1.8049e-06,-2.56879e-10,7.92609e-15,16533.2,-76.5133], Tmin=(951.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + ring(methyleneoxirane) + radical(C=CCJO) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]C1=CC1[C]=O(10417)',
    structure = SMILES('O=[C]C1[CH]C1=O'),
    E0 = (248.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75952,0.0147239,3.65219e-05,-6.26321e-08,2.63027e-11,29925.5,20.3976], Tmin=(100,'K'), Tmax=(939.774,'K')), NASAPolynomial(coeffs=[12.2009,0.00643157,-1.14847e-06,2.03071e-10,-1.93705e-14,26742.6,-32.0581], Tmin=(939.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(cyclopropanone) + radical(CCJCC=O) + radical(CC(C)CJ=O)"""),
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
    label = '[CH]=C([O])C=C=O(9622)',
    structure = SMILES('[CH]=C([O])C=C=O'),
    E0 = (191.57,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.76064,'amu*angstrom^2'), symmetry=1, barrier=(17.4886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69693,0.0510066,-7.12192e-05,4.97097e-08,-1.34031e-11,23123.1,19.2487], Tmin=(100,'K'), Tmax=(918.747,'K')), NASAPolynomial(coeffs=[11.2242,0.00953169,-3.51219e-06,5.85078e-10,-3.72633e-14,21372.3,-25.9085], Tmin=(918.747,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.57,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[O]C1=CC([O])=C1(10457)',
    structure = SMILES('O=C1[CH]C(=O)[CH]1'),
    E0 = (138.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26614,0.0270565,1.04886e-05,-3.49854e-08,1.52517e-11,16689.5,13.8288], Tmin=(100,'K'), Tmax=(1007.68,'K')), NASAPolynomial(coeffs=[12.8171,0.0113344,-5.04905e-06,1.0568e-09,-8.16885e-14,13234.9,-43.7445], Tmin=(1007.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cs-(Cds-O2d)(Cds-O2d)HH) + group(Cds-OdCsCs) + group(Cds-OdCsCs) + ring(Cyclobutane) + radical(C=OCJC=O) + radical(C=OCJC=O)"""),
)

species(
    label = '[CH2]C([O])=[C][C]=O(9909)',
    structure = SMILES('[CH2]C([O])=[C][C]=O'),
    E0 = (369.381,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.680006,'amu*angstrom^2'), symmetry=1, barrier=(15.6347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.680123,'amu*angstrom^2'), symmetry=1, barrier=(15.6374,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15508,0.0670571,-0.000114989,9.48772e-08,-3.01365e-11,44524.5,22.9321], Tmin=(100,'K'), Tmax=(829.692,'K')), NASAPolynomial(coeffs=[11.9656,0.0108003,-5.80085e-06,1.13136e-09,-7.79261e-14,42873,-26.3458], Tmin=(829.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.381,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]=C([O])[C]=[C]O(10458)',
    structure = SMILES('[CH]C([O])=C=[C]O'),
    E0 = (471.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,540,610,2055,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.10266,'amu*angstrom^2'), symmetry=1, barrier=(48.3443,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.10264,'amu*angstrom^2'), symmetry=1, barrier=(48.3439,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.922479,0.0595296,-7.07663e-05,4.0462e-08,-8.66302e-12,56836.4,24.7207], Tmin=(100,'K'), Tmax=(1293.9,'K')), NASAPolynomial(coeffs=[15.8384,0.00615121,-4.61265e-07,-1.0236e-10,1.32376e-14,53584.7,-48.7328], Tmin=(1293.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=CJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]1[CH]C(=O)O1(10459)',
    structure = SMILES('[CH][C]1[CH]C(=O)O1'),
    E0 = (515.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.75756,0.0182381,2.35483e-05,-4.835e-08,2.18528e-11,62047,22.1707], Tmin=(100,'K'), Tmax=(897.827,'K')), NASAPolynomial(coeffs=[10.3204,0.00878545,-1.15937e-06,6.9026e-11,-3.2927e-15,59711.9,-18.9419], Tmin=(897.827,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(515.445,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(C2CsJOC(O)) + radical(CCJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([O])[C]=[C][O](10460)',
    structure = SMILES('[CH]C([O])[C][C]=O'),
    E0 = (887.178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1855,455,950,383.815,383.852,383.857,383.861,383.865,383.881],'cm^-1')),
        HinderedRotor(inertia=(0.010658,'amu*angstrom^2'), symmetry=1, barrier=(1.11438,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0106634,'amu*angstrom^2'), symmetry=1, barrier=(1.11495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0106602,'amu*angstrom^2'), symmetry=1, barrier=(1.11482,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4077,0.0615184,-0.000102541,8.52578e-08,-2.72651e-11,106792,24.8354], Tmin=(100,'K'), Tmax=(867.628,'K')), NASAPolynomial(coeffs=[10.0836,0.0137652,-6.57515e-06,1.2181e-09,-8.13857e-14,105578,-14.1068], Tmin=(867.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(887.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)OJ) + radical(CCJ2_triplet) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH][C]1C=[C]OO1(10461)',
    structure = SMILES('[CH]C1=C[C]OO1'),
    E0 = (795.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13807,0.0330236,-6.76079e-06,-1.41014e-08,7.19338e-12,95754.3,20.035], Tmin=(100,'K'), Tmax=(1052.17,'K')), NASAPolynomial(coeffs=[10.8916,0.0165625,-7.26837e-06,1.41099e-09,-1.01807e-13,92981.4,-27.0647], Tmin=(1052.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(795.532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(12dioxolene) + radical(CH2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1([O])C=[C]O1(10462)',
    structure = SMILES('[CH]C1([O])C=[C]O1'),
    E0 = (670.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49546,0.0356061,1.7056e-05,-7.36351e-08,3.8619e-11,80756.3,20.0326], Tmin=(100,'K'), Tmax=(895.737,'K')), NASAPolynomial(coeffs=[24.3231,-0.0131893,9.77411e-06,-1.97974e-09,1.33601e-13,74534.8,-99.4849], Tmin=(895.737,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(670.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=CC(C)(O)OJ) + radical(C=CJO) + radical(CCJ2_triplet)"""),
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
    label = '[CH][C]([O])[C]=C[O](10463)',
    structure = SMILES('[CH]C([O])=[C][CH][O]'),
    E0 = (691.735,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,1685,370,500.078,500.095,500.096,500.098,500.106,500.109],'cm^-1')),
        HinderedRotor(inertia=(0.000674005,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.291395,'amu*angstrom^2'), symmetry=1, barrier=(51.7217,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.12248,0.0425601,-3.93424e-05,2.06697e-08,-4.62146e-12,83263.1,24.0695], Tmin=(100,'K'), Tmax=(1048.52,'K')), NASAPolynomial(coeffs=[7.55668,0.0218294,-9.68559e-06,1.81361e-09,-1.25621e-13,82123.5,-2.40485], Tmin=(1048.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(691.735,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJO) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C([O])[C]=C=O(10464)',
    structure = SMILES('[CH]C([O])[C]=C=O'),
    E0 = (598.326,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,1685,370,2120,512.5,787.5,262.074,268.193,268.472,270.255,272.1],'cm^-1')),
        HinderedRotor(inertia=(0.000156898,'amu*angstrom^2'), symmetry=1, barrier=(0.120248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.26824,'amu*angstrom^2'), symmetry=1, barrier=(67.3047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32442,0.0635505,-0.000106075,8.80064e-08,-2.79087e-11,72053.9,21.2614], Tmin=(100,'K'), Tmax=(888.495,'K')), NASAPolynomial(coeffs=[10.2075,0.0140344,-6.40037e-06,1.15197e-09,-7.52784e-14,70851.3,-18.4287], Tmin=(888.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.326,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CC(C)OJ) + radical(CCCJ=C=O) + radical(CCJ2_triplet)"""),
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
    E0 = (337.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1171.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1063.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (345.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (400.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (755.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (645.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (470.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (645.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (616.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (792.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (536.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (636.723,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (905.475,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (847.773,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (340.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (340.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (340.117,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (590.233,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (611.392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (337.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (345.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (553.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (664.663,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (579.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (910.039,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (795.532,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (670.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (1105.85,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (714.596,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (776.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['HCCO(2227)', 'HCCO(2227)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH][C]=C[C]=O(10445)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=O(2355)', '[CH]C(=[CH])[O](10446)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[CH]C1=CC(=O)O1(10447)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[CH]C(O)=C=C=O(10448)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[CH]=C1OC1[C][O](10442)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(418.743,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 416.2 to 418.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[CH]C1([O])[CH]C1=O(10449)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(308.738,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 305.8 to 308.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['CO(2039)', '[CH]C(=[CH])[O](10446)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(3.02e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C([O])=C=C=O(10450)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Cdd;HJ] for rate rule [Ca_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH](6993)', 'O=[C]C=C=O(9594)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.96577e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Ck_O;YJ] for rate rule [Ck_O;CH_quartet]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH]C([O])=[C][C]=O(10451)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]C([O])=[C]C=O(10452)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(O)=[C][C]=O(10453)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C][O](6861)', '[CH]=[C][O](6861)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[C]=C([O])[CH][C]=O(10454)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[CH]=C1OC1[C]=O(10434)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['O=[C][CH]C1=CO1(10455)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[O]C1=CC1[C]=O(10417)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HCCO(2227)', '[CH]=[C][O](6861)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction20',
    reactants = ['O(4)', '[CH]=C=C[C]=O(10456)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[CH]=C([O])C=C=O(9622)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[O]C1=CC([O])=C1(10457)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleH_rad_out;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C([O])=[C][C]=O(9909)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.97864e+06,'s^-1'), n=1.84533, Ea=(184.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SD;Y_rad_out;Cd_H_out_singleH] + [R3H_SD;Cd_rad_out;Cd_H_out_single] for rate rule [R3H_SD;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=C([O])[C]=[C]O(10458)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[CH][C]1[CH]C(=O)O1(10459)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra] for rate rule [R4_linear;carbonyl_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]C([O])[C]=[C][O](10460)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[CH][C]1C=[C]OO1(10461)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.44463e+09,'s^-1'), n=0.470283, Ea=(458.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra] for rate rule [R5_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 454.5 to 458.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C([O])C=[C][O](9624)'],
    products = ['[CH]C1([O])C=[C]O1(10462)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1'), n=0.478, Ea=(333.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;multiplebond_intra;radadd_intra_O] for rate rule [R5;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic
Ea raised from 330.1 to 333.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['HCCO(2227)', '[CH][C][O](10218)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH][C]([O])[C]=C[O](10463)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C([O])[C]=C=O(10464)'],
    products = ['[CH]=C([O])C=[C][O](9624)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(3.9024e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2546',
    isomers = [
        '[CH]=C([O])C=[C][O](9624)',
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
    label = 'PDepNetwork #2546',
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


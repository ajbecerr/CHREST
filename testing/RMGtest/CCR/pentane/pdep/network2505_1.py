species(
    label = '[CH]=C([O])C[CH]C(C)O[O](8113)',
    structure = SMILES('[CH]=C([O])C[CH]C(C)O[O]'),
    E0 = (357.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.347069,0.0995292,-0.000123981,8.4085e-08,-2.29492e-11,43176.2,38.3271], Tmin=(100,'K'), Tmax=(892.869,'K')), NASAPolynomial(coeffs=[14.113,0.0347468,-1.51444e-05,2.81852e-09,-1.94095e-13,40594.1,-29.7953], Tmin=(892.869,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJCOOH) + radical(Cds_P)"""),
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
    label = 'C=CC(C)O[O](3748)',
    structure = SMILES('C=CC(C)O[O]'),
    E0 = (34.1812,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,288.144],'cm^-1')),
        HinderedRotor(inertia=(0.201344,'amu*angstrom^2'), symmetry=1, barrier=(11.8628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201345,'amu*angstrom^2'), symmetry=1, barrier=(11.8628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201345,'amu*angstrom^2'), symmetry=1, barrier=(11.8628,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3658.84,'J/mol'), sigma=(6.16258,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.50 K, Pc=35.47 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42682,0.0494385,-3.17382e-05,7.83709e-09,1.46698e-14,4209.65,24.1077], Tmin=(100,'K'), Tmax=(1166.77,'K')), NASAPolynomial(coeffs=[11.9674,0.0225403,-9.03365e-06,1.64983e-09,-1.13538e-13,1121.19,-31.0646], Tmin=(1166.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(34.1812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
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
    label = '[CH]=C([O])C[CH]CO[O](3041)',
    structure = SMILES('[CH]=C([O])C[CH]CO[O]'),
    E0 = (394.138,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,360.205,360.318,2703.31],'cm^-1')),
        HinderedRotor(inertia=(0.111631,'amu*angstrom^2'), symmetry=1, barrier=(10.2895,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111632,'amu*angstrom^2'), symmetry=1, barrier=(10.288,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111668,'amu*angstrom^2'), symmetry=1, barrier=(10.2881,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111596,'amu*angstrom^2'), symmetry=1, barrier=(10.2866,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4356.36,'J/mol'), sigma=(7.16244,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.45 K, Pc=26.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.201024,0.0882312,-0.000126201,9.81731e-08,-3.04163e-11,47536.2,34.4513], Tmin=(100,'K'), Tmax=(827.563,'K')), NASAPolynomial(coeffs=[11.9577,0.0284615,-1.25286e-05,2.30182e-09,-1.55588e-13,45691.1,-19.4333], Tmin=(827.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(394.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])CC(C)[CH]O[O](10879)',
    structure = SMILES('[CH]=C([O])CC(C)[CH]O[O]'),
    E0 = (349.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.350652,0.0998667,-0.000125644,8.7018e-08,-2.43122e-11,42186.1,36.2489], Tmin=(100,'K'), Tmax=(872.209,'K')), NASAPolynomial(coeffs=[13.608,0.0358546,-1.55628e-05,2.88224e-09,-1.97624e-13,39751,-29.1857], Tmin=(872.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.483,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCsJOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C([CH2])C(C)O[O](8069)',
    structure = SMILES('[CH]=C([O])C([CH2])C(C)O[O]'),
    E0 = (354.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,364.051,371.232],'cm^-1')),
        HinderedRotor(inertia=(0.00119057,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0972628,'amu*angstrom^2'), symmetry=1, barrier=(9.6692,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954843,'amu*angstrom^2'), symmetry=1, barrier=(9.64801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0944899,'amu*angstrom^2'), symmetry=1, barrier=(9.6555,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0966947,'amu*angstrom^2'), symmetry=1, barrier=(9.68754,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4529.55,'J/mol'), sigma=(7.53448,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=707.51 K, Pc=24.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.514858,0.0985484,-0.000117065,7.39011e-08,-1.84444e-11,42788.9,38.1515], Tmin=(100,'K'), Tmax=(982.9,'K')), NASAPolynomial(coeffs=[16.5847,0.0289598,-1.08659e-05,1.86956e-09,-1.2313e-13,39427.5,-44.0492], Tmin=(982.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'O2(S)(666)',
    structure = SMILES('O=O'),
    E0 = (85.6848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53732,-0.00121572,5.3162e-06,-4.89446e-09,1.45846e-12,10304.5,4.68368], Tmin=(100,'K'), Tmax=(1074.55,'K')), NASAPolynomial(coeffs=[3.15382,0.00167804,-7.69974e-07,1.51275e-10,-1.08782e-14,10302.3,6.16756], Tmin=(1074.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(85.6848,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]=C([O])CC=CC(5457)',
    structure = SMILES('[CH]=C([O])CC=CC'),
    E0 = (223.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3120,650,792.5,1650,277.716,278.855],'cm^-1')),
        HinderedRotor(inertia=(0.168289,'amu*angstrom^2'), symmetry=1, barrier=(9.22014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.167602,'amu*angstrom^2'), symmetry=1, barrier=(9.21523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.168656,'amu*angstrom^2'), symmetry=1, barrier=(9.22786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.772519,0.0631786,-4.77308e-05,1.63194e-08,-1.45174e-12,26951.5,28.177], Tmin=(100,'K'), Tmax=(1090.16,'K')), NASAPolynomial(coeffs=[14.0962,0.0247788,-9.32459e-06,1.65714e-09,-1.1292e-13,23423.4,-40.1105], Tmin=(1090.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
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
    label = '[CH]=C([O])C[CH]C(C)[O](10880)',
    structure = SMILES('[CH]=C([O])C[CH]C(C)[O]'),
    E0 = (364.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,180,180,557.576,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00964257,'amu*angstrom^2'), symmetry=1, barrier=(2.12763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0707347,'amu*angstrom^2'), symmetry=1, barrier=(15.6052,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0707299,'amu*angstrom^2'), symmetry=1, barrier=(15.6051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0707339,'amu*angstrom^2'), symmetry=1, barrier=(15.6053,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0330346,0.0801481,-7.63091e-05,3.73772e-08,-7.22425e-12,43937.8,35.6691], Tmin=(100,'K'), Tmax=(1260.07,'K')), NASAPolynomial(coeffs=[17.9352,0.0231094,-8.40985e-06,1.45368e-09,-9.69846e-14,39409.6,-55.1709], Tmin=(1260.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(CCJCO) + radical(Cds_P)"""),
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
    label = '[CH]CC(=[CH])[O](9655)',
    structure = SMILES('[CH]CC(=[CH])[O]'),
    E0 = (602.036,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3120,650,792.5,1650,485.679,486.003,486.553,3552.34],'cm^-1')),
        HinderedRotor(inertia=(0.0198899,'amu*angstrom^2'), symmetry=1, barrier=(3.32905,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.119761,'amu*angstrom^2'), symmetry=1, barrier=(19.9624,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63798,0.0445574,-4.56301e-05,2.33077e-08,-4.55552e-12,72499.5,21.5227], Tmin=(100,'K'), Tmax=(1363.28,'K')), NASAPolynomial(coeffs=[12.9679,0.00812049,-2.02506e-06,2.65656e-10,-1.49276e-14,69707.1,-35.5602], Tmin=(1363.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(602.036,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([CH2])[O](9653)',
    structure = SMILES('[CH]C(=C)[O]'),
    E0 = (322.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,387.097,387.511,388.12,388.534],'cm^-1')),
        HinderedRotor(inertia=(0.483854,'amu*angstrom^2'), symmetry=1, barrier=(51.3097,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.72438,0.0238014,-4.51958e-06,-8.49326e-09,4.40325e-12,38798.3,14.6668], Tmin=(100,'K'), Tmax=(1019.83,'K')), NASAPolynomial(coeffs=[7.15237,0.0155204,-5.70458e-06,1.01811e-09,-7.00768e-14,37422.6,-9.09941], Tmin=(1019.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(322.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(C)O[O](668)',
    structure = SMILES('[CH]C(C)O[O]'),
    E0 = (380.16,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.155785,'amu*angstrom^2'), symmetry=1, barrier=(3.5818,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155733,'amu*angstrom^2'), symmetry=1, barrier=(3.58061,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155864,'amu*angstrom^2'), symmetry=1, barrier=(3.58362,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0706,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89717,0.0487248,-6.20882e-05,4.51124e-08,-1.32772e-11,45796.1,20.4395], Tmin=(100,'K'), Tmax=(828.346,'K')), NASAPolynomial(coeffs=[8.10163,0.0187643,-7.83515e-06,1.449e-09,-9.94636e-14,44768.2,-8.32502], Tmin=(828.346,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.16,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[CH]C(C)O[O](10881)',
    structure = SMILES('[CH]=[C]C[CH]C(C)O[O]'),
    E0 = (672.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,309.537,311.294],'cm^-1')),
        HinderedRotor(inertia=(0.130486,'amu*angstrom^2'), symmetry=1, barrier=(8.74419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128008,'amu*angstrom^2'), symmetry=1, barrier=(8.74471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127709,'amu*angstrom^2'), symmetry=1, barrier=(8.74464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00174065,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127766,'amu*angstrom^2'), symmetry=1, barrier=(8.74552,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.159487,0.0906289,-0.000123796,9.93742e-08,-3.27404e-11,80992.6,35.4159], Tmin=(100,'K'), Tmax=(787.108,'K')), NASAPolynomial(coeffs=[9.52753,0.0390723,-1.80181e-05,3.4081e-09,-2.35203e-13,79640.2,-6.75981], Tmin=(787.108,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJCOOH) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=C([O])C[C]C(C)O[O](10882)',
    structure = SMILES('[CH]=C([O])C[C]C(C)O[O]'),
    E0 = (605.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (127.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.460673,0.101015,-0.000129331,8.70441e-08,-2.32358e-11,72984.7,36.3331], Tmin=(100,'K'), Tmax=(917.524,'K')), NASAPolynomial(coeffs=[15.8918,0.0297252,-1.27833e-05,2.36128e-09,-1.61932e-13,69983.9,-41.1504], Tmin=(917.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(605.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[C]=C([O])C[CH]C(C)O[O](10883)',
    structure = SMILES('[C]C(=O)C[CH]C(C)O[O]'),
    E0 = (615.178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,375,552.5,462.5,1710,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (127.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.367474,0.104584,-0.000158949,1.34941e-07,-4.56937e-11,74137.9,35.786], Tmin=(100,'K'), Tmax=(816.938,'K')), NASAPolynomial(coeffs=[10.4944,0.0399211,-1.91406e-05,3.64825e-09,-2.51411e-13,72746.2,-12.0755], Tmin=(816.938,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(615.178,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(ROOJ) + radical(CCJCOOH) + radical(CJ3)"""),
)

species(
    label = 'CC([CH]CC1=CO1)O[O](10884)',
    structure = SMILES('CC([CH]CC1=CO1)O[O]'),
    E0 = (295.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.228627,0.0910502,-9.31057e-05,4.97522e-08,-1.06544e-11,35655,35.7407], Tmin=(100,'K'), Tmax=(1128.55,'K')), NASAPolynomial(coeffs=[16.7831,0.0307547,-1.29652e-05,2.41123e-09,-1.67362e-13,31815.2,-48.3885], Tmin=(1128.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(295.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(ROOJ) + radical(CCJCOOH)"""),
)

species(
    label = '[CH]=C([O])CC1OOC1C(10885)',
    structure = SMILES('[CH]=C([O])CC1OOC1C'),
    E0 = (162.007,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0321527,0.0837477,-7.57178e-05,3.57889e-08,-6.8213e-12,19634.1,29.5839], Tmin=(100,'K'), Tmax=(1256.06,'K')), NASAPolynomial(coeffs=[16.5516,0.0309363,-1.26506e-05,2.31581e-09,-1.59065e-13,15468,-54.2042], Tmin=(1256.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(162.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC(C(C)O[O])O1(10797)',
    structure = SMILES('[CH]=C1CC(C(C)O[O])O1'),
    E0 = (182.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.206511,0.0773519,-3.50799e-05,-2.49571e-08,2.01363e-11,22138.7,29.1341], Tmin=(100,'K'), Tmax=(914.353,'K')), NASAPolynomial(coeffs=[22.6435,0.0177975,-3.66865e-06,4.71955e-10,-3.11443e-14,16271,-88.2948], Tmin=(914.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.695,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'CC(O[O])C1C=C([O])C1(10748)',
    structure = SMILES('CC(O[O])C1[CH]C(=O)C1'),
    E0 = (40.6366,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.550564,0.0618699,-1.2549e-05,-2.63324e-08,1.37664e-11,5023.8,31.6212], Tmin=(100,'K'), Tmax=(1028.42,'K')), NASAPolynomial(coeffs=[16.064,0.0317767,-1.27716e-05,2.40897e-09,-1.72075e-13,233.449,-51.4335], Tmin=(1028.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(40.6366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(ROOJ) + radical(CCJC=O)"""),
)

species(
    label = '[CH]=C([O])CC1OC1C(10886)',
    structure = SMILES('[CH]=C([O])CC1OC1C'),
    E0 = (110.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.76597,0.0823128,-7.82183e-05,3.84781e-08,-7.10401e-12,13506.3,31.9988], Tmin=(100,'K'), Tmax=(1573.43,'K')), NASAPolynomial(coeffs=[18.2137,0.0184088,-2.37372e-06,1.97097e-11,1.12155e-14,9471.33,-62.0127], Tmin=(1573.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'HO2(9)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26737e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.914,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11273e-07,6.34955e-11,-4.86383e-15,83.4204,3.0934], Tmin=(923.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH]=C([O])CC=C[CH2](5449)',
    structure = SMILES('[CH]=C([O])CC=C[CH2]'),
    E0 = (374.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,460.616,461.315],'cm^-1')),
        HinderedRotor(inertia=(0.0930579,'amu*angstrom^2'), symmetry=1, barrier=(14.037,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0929739,'amu*angstrom^2'), symmetry=1, barrier=(14.036,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0927408,'amu*angstrom^2'), symmetry=1, barrier=(14.0428,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3874.43,'J/mol'), sigma=(6.44487,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.18 K, Pc=32.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.84923,0.0598081,-3.62764e-05,-5.38171e-10,5.92169e-12,45171.8,28.1803], Tmin=(100,'K'), Tmax=(978.941,'K')), NASAPolynomial(coeffs=[15.7465,0.0197907,-6.91215e-06,1.22487e-09,-8.56623e-14,41255.8,-48.4768], Tmin=(978.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])CC=C(C)OO(10887)',
    structure = SMILES('[CH]=C([O])CC=C(C)OO'),
    E0 = (133.003,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.10356,0.092222,-0.000101517,5.97685e-08,-1.42266e-11,16142.7,35.4301], Tmin=(100,'K'), Tmax=(1015.15,'K')), NASAPolynomial(coeffs=[14.7982,0.0335053,-1.47573e-05,2.79254e-09,-1.95298e-13,13117.2,-36.6864], Tmin=(1015.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)C=CC(C)O[O](10888)',
    structure = SMILES('[CH]=C(O)C=CC(C)O[O]'),
    E0 = (119.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.970054,0.100111,-0.000104755,5.15518e-08,-8.96851e-12,14518.3,34.9852], Tmin=(100,'K'), Tmax=(1004.01,'K')), NASAPolynomial(coeffs=[22.9358,0.0196826,-6.7245e-06,1.15411e-09,-7.85104e-14,8971.33,-84.1604], Tmin=(1004.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(119.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])C=CC(C)O[O](8109)',
    structure = SMILES('C=C([O])C=CC(C)O[O]'),
    E0 = (9.86503,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.580268,'amu*angstrom^2'), symmetry=1, barrier=(13.3415,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.580247,'amu*angstrom^2'), symmetry=1, barrier=(13.341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.579982,'amu*angstrom^2'), symmetry=1, barrier=(13.3349,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.580276,'amu*angstrom^2'), symmetry=1, barrier=(13.3417,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.509993,0.0923529,-9.45531e-05,4.96226e-08,-1.02496e-11,1354.69,34.4333], Tmin=(100,'K'), Tmax=(1182.67,'K')), NASAPolynomial(coeffs=[19.3218,0.0252779,-9.48041e-06,1.66728e-09,-1.12491e-13,-3336.19,-64.571], Tmin=(1182.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.86503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C([O])C=CC(C)OO(10889)',
    structure = SMILES('[CH]=C([O])C=CC(C)OO'),
    E0 = (104.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.827178,0.099211,-0.000105528,5.62901e-08,-1.17372e-11,12803.1,35.2168], Tmin=(100,'K'), Tmax=(1174.88,'K')), NASAPolynomial(coeffs=[21.3177,0.0238158,-9.26808e-06,1.66881e-09,-1.1437e-13,7599.61,-75.1885], Tmin=(1174.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)CC=C(C)O[O](9614)',
    structure = SMILES('[CH]=C(O)CC=C(C)O[O]'),
    E0 = (147.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.355342,0.0944483,-0.000105513,6.12759e-08,-1.41017e-11,17862.5,35.5857], Tmin=(100,'K'), Tmax=(1060.2,'K')), NASAPolynomial(coeffs=[17.2806,0.0279104,-1.13736e-05,2.07979e-09,-1.43033e-13,14122.9,-50.5283], Tmin=(1060.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])CC=C(C)O[O](8108)',
    structure = SMILES('C=C([O])CC=C(C)O[O]'),
    E0 = (37.9114,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,210.548,211.124,4000],'cm^-1')),
        HinderedRotor(inertia=(0.29326,'amu*angstrom^2'), symmetry=1, barrier=(9.27851,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295845,'amu*angstrom^2'), symmetry=1, barrier=(9.2778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294921,'amu*angstrom^2'), symmetry=1, barrier=(9.27678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.294163,'amu*angstrom^2'), symmetry=1, barrier=(9.27655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.208086,0.0854287,-9.07646e-05,5.33821e-08,-1.28557e-11,4694.57,34.6665], Tmin=(100,'K'), Tmax=(998.836,'K')), NASAPolynomial(coeffs=[12.8143,0.0349476,-1.49586e-05,2.78848e-09,-1.93212e-13,2176.13,-26.1372], Tmin=(998.836,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(37.9114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = '[CH]=C([O])CC=C(C)O[O](10890)',
    structure = SMILES('[CH]=C([O])CC=C(C)O[O]'),
    E0 = (285.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.240091,'amu*angstrom^2'), symmetry=1, barrier=(5.52016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240176,'amu*angstrom^2'), symmetry=1, barrier=(5.52213,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240102,'amu*angstrom^2'), symmetry=1, barrier=(5.52042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.240042,'amu*angstrom^2'), symmetry=1, barrier=(5.51904,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (127.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0884402,0.0894648,-0.000109102,7.28224e-08,-1.96171e-11,34416.4,35.5009], Tmin=(100,'K'), Tmax=(903.16,'K')), NASAPolynomial(coeffs=[13.0205,0.032189,-1.39747e-05,2.60287e-09,-1.79522e-13,32080.5,-25.5713], Tmin=(903.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C=CC(C)O[O](10891)',
    structure = SMILES('[CH]=C([O])C=CC(C)O[O]'),
    E0 = (256.961,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.618365,'amu*angstrom^2'), symmetry=1, barrier=(14.2174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618444,'amu*angstrom^2'), symmetry=1, barrier=(14.2192,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618715,'amu*angstrom^2'), symmetry=1, barrier=(14.2255,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.618544,'amu*angstrom^2'), symmetry=1, barrier=(14.2215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (127.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.471809,0.0945554,-0.000106626,6.11635e-08,-1.37391e-11,31069.6,34.7], Tmin=(100,'K'), Tmax=(1092.4,'K')), NASAPolynomial(coeffs=[18.9356,0.0234895,-9.04036e-06,1.60739e-09,-1.09055e-13,26829.6,-60.6443], Tmin=(1092.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C#CC[CH]C(C)O[O](10892)',
    structure = SMILES('C#CC[CH]C(C)O[O]'),
    E0 = (353.475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,192.863,193.005],'cm^-1')),
        HinderedRotor(inertia=(0.185774,'amu*angstrom^2'), symmetry=1, barrier=(4.90412,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.670415,'amu*angstrom^2'), symmetry=1, barrier=(17.7149,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.185677,'amu*angstrom^2'), symmetry=1, barrier=(4.90434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00453206,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.51886,'amu*angstrom^2'), symmetry=1, barrier=(66.4373,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.127,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.268267,0.0856173,-0.000104428,7.37272e-08,-2.13074e-11,42644.4,31.995], Tmin=(100,'K'), Tmax=(840.576,'K')), NASAPolynomial(coeffs=[10.8142,0.0354325,-1.4873e-05,2.69988e-09,-1.82605e-13,40871.5,-17.0515], Tmin=(840.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(CCJCOOH)"""),
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
    label = '[CH]=C([O])CC=CO[O](10694)',
    structure = SMILES('[CH]=C([O])CC=CO[O]'),
    E0 = (326.799,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.422444,'amu*angstrom^2'), symmetry=1, barrier=(9.71281,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.423053,'amu*angstrom^2'), symmetry=1, barrier=(9.72682,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.420301,'amu*angstrom^2'), symmetry=1, barrier=(9.66356,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.355939,0.0751348,-8.59647e-05,4.94984e-08,-1.10426e-11,39440.4,32.5241], Tmin=(100,'K'), Tmax=(1106.9,'K')), NASAPolynomial(coeffs=[16.8498,0.0155305,-5.19222e-06,8.50216e-10,-5.50446e-14,35789,-48.7245], Tmin=(1106.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.799,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_P)"""),
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
    label = '[CH2][CH]C(C)O[O](3749)',
    structure = SMILES('[CH2][CH]C(C)O[O]'),
    E0 = (312.397,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180],'cm^-1')),
        HinderedRotor(inertia=(0.0973544,'amu*angstrom^2'), symmetry=1, barrier=(2.23837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0900935,'amu*angstrom^2'), symmetry=1, barrier=(2.07143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0964816,'amu*angstrom^2'), symmetry=1, barrier=(2.2183,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.230235,'amu*angstrom^2'), symmetry=1, barrier=(5.29355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (87.0972,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49109,0.0589098,-6.48732e-05,4.40834e-08,-1.28398e-11,37659.6,26.5908], Tmin=(100,'K'), Tmax=(818.389,'K')), NASAPolynomial(coeffs=[7.16902,0.0311595,-1.40131e-05,2.65446e-09,-1.84819e-13,36730.2,0.335789], Tmin=(818.389,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCJCOOH) + radical(RCCJ)"""),
)

species(
    label = 'O2(2)',
    structure = SMILES('[O][O]'),
    E0 = (-8.62178,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1483.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (31.9988,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.6,'angstroms^3'), rotrelaxcollnum=3.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53764,-0.00122828,5.36761e-06,-4.93131e-09,1.45956e-12,-1037.99,4.6718], Tmin=(100,'K'), Tmax=(1087.7,'K')), NASAPolynomial(coeffs=[3.16427,0.00169454,-8.00338e-07,1.59031e-10,-1.14892e-14,-1048.44,6.08305], Tmin=(1087.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.62178,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""O2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH]=C([O])C[CH][CH]C(4320)',
    structure = SMILES('[CH]=C([O])C[CH][CH]C'),
    E0 = (495.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,3120,650,792.5,1650,180,2310.2,2310.53],'cm^-1')),
        HinderedRotor(inertia=(0.0972106,'amu*angstrom^2'), symmetry=1, barrier=(2.23506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0973647,'amu*angstrom^2'), symmetry=1, barrier=(2.23861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0973159,'amu*angstrom^2'), symmetry=1, barrier=(2.23748,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0972963,'amu*angstrom^2'), symmetry=1, barrier=(2.23703,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.8,'J/mol'), sigma=(6.58551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.90 K, Pc=30.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12164,0.0685933,-8.64924e-05,7.13449e-08,-2.43724e-11,59653.1,32.1607], Tmin=(100,'K'), Tmax=(839.736,'K')), NASAPolynomial(coeffs=[5.50699,0.038617,-1.67145e-05,3.06153e-09,-2.07026e-13,59237,13.6776], Tmin=(839.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C[CH][CH]O[O](10697)',
    structure = SMILES('[CH]=C([O])C[CH][CH]O[O]'),
    E0 = (582.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,3120,650,792.5,1650,287.451,287.454,1137.49],'cm^-1')),
        HinderedRotor(inertia=(0.00203874,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130408,'amu*angstrom^2'), symmetry=1, barrier=(7.64235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130256,'amu*angstrom^2'), symmetry=1, barrier=(7.6418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.130305,'amu*angstrom^2'), symmetry=1, barrier=(7.64207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.154687,0.0915889,-0.000145325,1.213e-07,-3.93276e-11,70216.8,34.7649], Tmin=(100,'K'), Tmax=(867.471,'K')), NASAPolynomial(coeffs=[11.2931,0.0274962,-1.24822e-05,2.28786e-09,-1.5282e-13,68763.4,-14.6268], Tmin=(867.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJCOOH) + radical(CCsJOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C[CH][C](C)O[O](10893)',
    structure = SMILES('[CH]=C([O])C[CH][C](C)O[O]'),
    E0 = (544.613,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (127.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.382958,0.104246,-0.000157508,1.31104e-07,-4.31472e-11,65652.2,38.7391], Tmin=(100,'K'), Tmax=(843.997,'K')), NASAPolynomial(coeffs=[11.1753,0.0373511,-1.70851e-05,3.17572e-09,-2.15209e-13,64132.7,-12.5059], Tmin=(843.997,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(544.613,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C2CsJOOH) + radical(CCJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])[CH][CH]C(C)O[O](10894)',
    structure = SMILES('[CH]C([O])=C[CH]C(C)O[O]'),
    E0 = (385.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,350,440,435,1725,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (127.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.449162,0.097912,-0.000103117,5.82285e-08,-1.33128e-11,46532.3,33.8494], Tmin=(100,'K'), Tmax=(1055.45,'K')), NASAPolynomial(coeffs=[15.8262,0.0362301,-1.54543e-05,2.85603e-09,-1.96828e-13,43096.8,-45.5479], Tmin=(1055.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJCO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])C[CH]C([CH2])O[O](10895)',
    structure = SMILES('[CH]=C([O])C[CH]C([CH2])O[O]'),
    E0 = (571.677,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (127.118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.503068,0.103367,-0.00014128,1.02369e-07,-2.94185e-11,68915.1,40.8787], Tmin=(100,'K'), Tmax=(853.979,'K')), NASAPolynomial(coeffs=[15.0535,0.0304992,-1.32882e-05,2.44919e-09,-1.66794e-13,66258.1,-31.7168], Tmin=(853.979,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.677,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJCOOH) + radical(CJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])CC[C](C)O[O](10896)',
    structure = SMILES('[CH]=C([O])CC[C](C)O[O]'),
    E0 = (344.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.34202,0.101382,-0.000138396,1.07393e-07,-3.38296e-11,41548.3,36.2912], Tmin=(100,'K'), Tmax=(796.023,'K')), NASAPolynomial(coeffs=[11.7021,0.0390395,-1.74887e-05,3.25979e-09,-2.23e-13,39688.5,-18.7046], Tmin=(796.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C2CsJOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])[CH]CC(C)O[O](10817)',
    structure = SMILES('[CH]C([O])=CCC(C)O[O]'),
    E0 = (268.642,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.202363,0.0936044,-9.40694e-05,5.26309e-08,-1.21505e-11,32460.6,35.7657], Tmin=(100,'K'), Tmax=(1036.26,'K')), NASAPolynomial(coeffs=[13.5585,0.0404882,-1.71847e-05,3.16904e-09,-2.17974e-13,29608.6,-31.1128], Tmin=(1036.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])C[CH][C](C)OO(10897)',
    structure = SMILES('[CH]=C([O])C[CH][C](C)OO'),
    E0 = (392.609,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.37008,0.104073,-0.000136872,9.61074e-08,-2.55078e-11,47370,37.9641], Tmin=(100,'K'), Tmax=(673.78,'K')), NASAPolynomial(coeffs=[12.3261,0.0398164,-1.85667e-05,3.53705e-09,-2.45641e-13,45406.8,-20.1471], Tmin=(673.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.609,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C2CsJOOH) + radical(CCJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])CCC([CH2])O[O](10898)',
    structure = SMILES('[CH]=C([O])CCC([CH2])O[O]'),
    E0 = (371.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.469028,0.100628,-0.000122805,7.9725e-08,-2.06388e-11,44811.5,38.4526], Tmin=(100,'K'), Tmax=(943.373,'K')), NASAPolynomial(coeffs=[15.742,0.0318915,-1.35118e-05,2.48905e-09,-1.70811e-13,41752.8,-38.8114], Tmin=(943.373,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(371.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)[CH][CH]C(C)O[O](10899)',
    structure = SMILES('[CH]C(O)=C[CH]C(C)O[O]'),
    E0 = (247.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.07504,0.105004,-0.000106698,5.56801e-08,-1.1508e-11,29986.3,34.5896], Tmin=(100,'K'), Tmax=(1177.12,'K')), NASAPolynomial(coeffs=[20.7339,0.0308952,-1.22631e-05,2.19694e-09,-1.49271e-13,24851.9,-74.1826], Tmin=(1177.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=CCJCO) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C([O])[CH][CH]C(C)O[O](8111)',
    structure = SMILES('[CH2]C([O])=C[CH]C(C)O[O]'),
    E0 = (173.791,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.62803,0.100257,-0.000112958,6.6413e-08,-1.54721e-11,21070.3,34.0694], Tmin=(100,'K'), Tmax=(1048.34,'K')), NASAPolynomial(coeffs=[17.8823,0.0296294,-1.18996e-05,2.14685e-09,-1.46269e-13,17189.3,-56.1059], Tmin=(1048.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(173.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJCO) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C([O])C[CH]C([CH2])OO(10900)',
    structure = SMILES('[CH]=C([O])C[CH]C([CH2])OO'),
    E0 = (419.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.646896,0.10548,-0.000131097,8.54812e-08,-2.21978e-11,50639.5,40.6403], Tmin=(100,'K'), Tmax=(940.209,'K')), NASAPolynomial(coeffs=[16.5498,0.032316,-1.43689e-05,2.71097e-09,-1.88657e-13,47405.9,-41.263], Tmin=(940.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJCOOH) + radical(CJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])[CH][CH]C(C)OO(10901)',
    structure = SMILES('[CH]C([O])=C[CH]C(C)OO'),
    E0 = (233.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.765683,0.102157,-0.000100784,5.19912e-08,-1.08142e-11,28264,34.2234], Tmin=(100,'K'), Tmax=(1156.13,'K')), NASAPolynomial(coeffs=[18.1941,0.0365594,-1.56748e-05,2.91416e-09,-2.01772e-13,23880,-59.9973], Tmin=(1156.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(233.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CCJCO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(O)C[CH][C](C)O[O](10902)',
    structure = SMILES('[CH]=C(O)C[CH][C](C)O[O]'),
    E0 = (406.809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.605046,0.106364,-0.000142704,1.03277e-07,-2.98939e-11,49089,38.0452], Tmin=(100,'K'), Tmax=(845.221,'K')), NASAPolynomial(coeffs=[14.551,0.0346383,-1.54136e-05,2.87699e-09,-1.97629e-13,46527,-32.5254], Tmin=(845.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ) + radical(C2CsJOOH) + radical(CCJCOOH) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])C[CH][C](C)O[O](8110)',
    structure = SMILES('C=C([O])C[CH][C](C)O[O]'),
    E0 = (297.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.235032,0.0927651,-0.000103499,4.51379e-08,5.71402e-12,35909.2,36.2041], Tmin=(100,'K'), Tmax=(582.58,'K')), NASAPolynomial(coeffs=[10.5276,0.0409405,-1.85835e-05,3.48925e-09,-2.39921e-13,34390.1,-10.6354], Tmin=(582.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(297.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C2CsJOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH]=C(O)C[CH]C([CH2])O[O](10903)',
    structure = SMILES('[CH]=C(O)C[CH]C([CH2])O[O]'),
    E0 = (433.872,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.855259,0.107198,-0.000133342,8.47518e-08,-2.11321e-11,52357.3,40.64], Tmin=(100,'K'), Tmax=(985.125,'K')), NASAPolynomial(coeffs=[18.9206,0.0269039,-1.10875e-05,2.02184e-09,-1.38316e-13,48460.8,-54.4717], Tmin=(985.125,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.872,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJCOOH) + radical(CJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([CH]CC(=C)[O])O[O](8112)',
    structure = SMILES('[CH2]C([CH]CC(=C)[O])O[O]'),
    E0 = (324.581,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.346884,0.0988345,-0.000120911,7.9887e-08,-2.1187e-11,39191.8,39.9177], Tmin=(100,'K'), Tmax=(918.932,'K')), NASAPolynomial(coeffs=[14.6277,0.0336507,-1.4507e-05,2.69178e-09,-1.85306e-13,36439.7,-31.0598], Tmin=(918.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.581,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJCOOH) + radical(CJCOOH)"""),
)

species(
    label = '[CH]C([O])[CH][CH]C(C)O[O](10904)',
    structure = SMILES('[CH]C([O])[CH][CH]C(C)O[O]'),
    E0 = (765.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.289555,0.0994573,-0.000122543,8.31352e-08,-2.29297e-11,92230.1,41.3557], Tmin=(100,'K'), Tmax=(879.568,'K')), NASAPolynomial(coeffs=[13.296,0.0376724,-1.71726e-05,3.26687e-09,-2.27915e-13,89840.3,-22.4428], Tmin=(879.568,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(765.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CCJCOOH) + radical(CCJCO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]1CC(C(C)O[O])O1(10905)',
    structure = SMILES('[CH][C]1CC(C(C)O[O])O1'),
    E0 = (487.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0707171,0.0875644,-8.71641e-05,4.07492e-08,-3.75035e-12,58714.8,32.4709], Tmin=(100,'K'), Tmax=(756.647,'K')), NASAPolynomial(coeffs=[13.3699,0.0332587,-1.12263e-05,1.78949e-09,-1.11824e-13,56244.2,-31.0082], Tmin=(756.647,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(ROOJ) + radical(C2CsJOCs) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]1C[CH]C(C)OOO1(10906)',
    structure = SMILES('[CH][C]1C[CH]C(C)OOO1'),
    E0 = (655.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.379285,0.0287097,0.000162451,-2.71581e-07,1.18139e-10,79059.8,35.5547], Tmin=(100,'K'), Tmax=(910.113,'K')), NASAPolynomial(coeffs=[40.7772,-0.0155479,1.57066e-05,-3.16585e-09,2.03359e-13,66186,-185.865], Tmin=(910.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.86,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(C2CsJOO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([O])CC1C(C)O[O](10854)',
    structure = SMILES('[CH]C1([O])CC1C(C)O[O]'),
    E0 = (499.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.623598,0.0903518,-8.50863e-05,4.04916e-08,-7.56885e-12,60240.8,34.3448], Tmin=(100,'K'), Tmax=(1303.9,'K')), NASAPolynomial(coeffs=[20.9521,0.0241643,-8.94516e-06,1.56194e-09,-1.0484e-13,54614.2,-75.4708], Tmin=(1303.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(ROOJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([O])C[CH]C(C)OO1(10907)',
    structure = SMILES('[CH]C1([O])C[CH]C(C)OO1'),
    E0 = (439.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (128.126,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0901546,0.0652319,9.92164e-06,-7.6238e-08,3.97086e-11,52977.9,30.1561], Tmin=(100,'K'), Tmax=(906.162,'K')), NASAPolynomial(coeffs=[23.6381,0.0168214,-1.87283e-06,7.24107e-11,-3.81242e-15,46430.2,-93.71], Tmin=(906.162,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(399.095,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxane) + radical(CC(C)(O)OJ) + radical(CCJCOOH) + radical(CCJ2_triplet)"""),
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
    E0 = (357.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (812.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (552.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (603.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (357.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (607.057,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (820.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (757.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (1191.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (817.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (826.97,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (360.643,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (365.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (365.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (365.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (487.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (500.389,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (421.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (421.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (421.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (382.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (382.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (382.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (508.089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (476.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (612.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (490.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (483.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (514.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (357.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (486.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (719.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (756.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (737.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (597.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (783.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (512.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (482.816,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (527.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (518.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (514.825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (549.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (500.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (417.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (592.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (493.694,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (507.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (489.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (788.459,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (487.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (655.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (499.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (479.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['HCCO(2227)', 'C=CC(C)O[O](3748)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH]=C([O])C[CH]CO[O](3041)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C([O])CC(C)[CH]O[O](10879)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;CH3] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C([O])C([CH2])C(C)O[O](8069)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['O2(S)(666)', '[CH]=C([O])CC=CC(5457)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', '[CH]=C([O])C[CH]C(C)[O](10880)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C[CH]O[O](57)', '[CH]CC(=[CH])[O](9655)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C([CH2])[O](9653)', '[CH]C(C)O[O](668)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O(4)', '[CH]=[C]C[CH]C(C)O[O](10881)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]=C([O])C[C]C(C)O[O](10882)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[C]=C([O])C[CH]C(C)O[O](10883)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['CC([CH]CC1=CO1)O[O](10884)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C([O])CC1OOC1C(10885)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C1CC(C(C)O[O])O1(10797)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['CC(O[O])C1C=C([O])C1(10748)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R4;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/NonDeC;CdsinglepriH_rad_out]
Euclidian distance = 3.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['O(4)', '[CH]=C([O])CC1OC1C(10886)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.27e+09,'s^-1'), n=1.06, Ea=(129.704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 36 used for R2OO_S;C_rad/H/NonDeC_intra;OOJ
Exact match found for rate rule [R2OO_S;C_rad/H/NonDeC_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['HO2(9)', '[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.6e+09,'s^-1','*|/',2.51189), n=1.2, Ea=(142.674,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2OO_2H
Exact match found for rate rule [R2OO_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C([O])CC=C(C)OO(10887)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C(O)C=CC(C)O[O](10888)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['C=C([O])C=CC(C)O[O](8109)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C([O])C=CC(C)OO(10889)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C(O)CC=C(C)O[O](9614)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['C=C([O])CC=C(C)O[O](8108)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH]=C([O])CC=C(C)O[O](10890)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[CH]=C([O])C=CC(C)O[O](10891)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O(4)', 'C#CC[CH]C(C)O[O](10892)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH3(17)', '[CH]=C([O])CC=CO[O](10694)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00632077,'m^3/(mol*s)'), n=2.46811, Ea=(26.7526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CsJ-HHH] for rate rule [Cds-OsH_Cds-CsH;CsJ-HHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C][O](6861)', 'C=CC(C)O[O](3748)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(532000,'cm^3/(mol*s)'), n=1.85, Ea=(23.8906,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\O2s/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['HCCO(2227)', '[CH2][CH]C(C)O[O](3749)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O2(2)', '[CH]=C([O])CC=CC(5457)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(9.03779e-05,'m^3/(mol*s)'), n=3.085, Ea=(143.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds;O2b] + [Cds-CsH_Cds-CsH;OJ] for rate rule [Cds-CsH_Cds-CsH;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 141.0 to 143.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['O2(2)', '[CH]=C([O])C[CH][CH]C(4320)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CH3(17)', '[CH]=C([O])C[CH][CH]O[O](10697)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH]=C([O])C[CH][C](C)O[O](10893)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH]C(C)O[O](3749)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', '[CH]=C([O])[CH][CH]C(C)O[O](10894)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH]=C([O])C[CH]C([CH2])O[O](10895)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C([O])CC[C](C)O[O](10896)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C([O])[CH]CC(C)O[O](10817)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C([O])C[CH][C](C)OO(10897)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C([O])CCC([CH2])O[O](10898)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C(O)[CH][CH]C(C)O[O](10899)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.39846e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['C=C([O])[CH][CH]C(C)O[O](8111)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C([O])C[CH]C([CH2])OO(10900)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.989e+11,'s^-1'), n=0.15, Ea=(143.135,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 246 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]=C([O])[CH][CH]C(C)OO(10901)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(574267,'s^-1'), n=1.61427, Ea=(59.4591,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R5Hall;Y_rad_out;Cs_H_out_H/OneDe] + [R5Hall;O_rad_out;Cs_H_out_1H] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/OneDe]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=C(O)C[CH][C](C)O[O](10902)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.46277e+08,'s^-1'), n=1.29592, Ea=(185.417,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_NonDe;XH_out] for rate rule [R5HJ_1;C_rad_out_NDMustO;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['C=C([O])C[CH][C](C)O[O](8110)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(6.04e+10,'s^-1'), n=0.59, Ea=(135.98,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_NonDe] for rate rule [R5HJ_3;Cd_rad_out_singleH;Cs_H_out_NDMustO]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=C(O)C[CH]C([CH2])O[O](10903)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_2;C_rad_out_2H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH2]C([CH]CC(=C)[O])O[O](8112)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(22.7193,'s^-1'), n=3.21897, Ea=(132.277,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_2H] for rate rule [R6HJ_3;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]C([O])[CH][CH]C(C)O[O](10904)'],
    products = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH][C]1CC(C(C)O[O])O1(10905)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(4.7342e+08,'s^-1'), n=0.920995, Ea=(129.298,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 127.9 to 129.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH][C]1C[CH]C(C)OOO1(10906)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(298.146,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra;radadd_intra] for rate rule [R7_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 289.7 to 298.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]C1([O])CC1C(C)O[O](10854)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(6.44049e+09,'s^-1'), n=0.679905, Ea=(141.692,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHNd] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_csHNd]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 140.0 to 141.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['[CH]C1([O])C[CH]C(C)OO1(10907)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(9.291e+11,'s^-1'), n=0.234, Ea=(121.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7;multiplebond_intra;radadd_intra_O] for rate rule [R7;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

network(
    label = 'PDepNetwork #2505',
    isomers = [
        '[CH]=C([O])C[CH]C(C)O[O](8113)',
    ],
    reactants = [
        ('HCCO(2227)', 'C=CC(C)O[O](3748)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2505',
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


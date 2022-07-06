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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.84923,0.0598081,-3.62764e-05,-5.38171e-10,5.92169e-12,45171.8,28.1803], Tmin=(100,'K'), Tmax=(978.941,'K')), NASAPolynomial(coeffs=[15.7465,0.0197907,-6.91215e-06,1.22487e-09,-8.56623e-14,41255.8,-48.4768], Tmin=(978.941,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_P)"""),
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
    label = 'butadiene13(1350)',
    structure = SMILES('C=CC=C'),
    E0 = (96.4553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180],'cm^-1')),
        HinderedRotor(inertia=(1.30711,'amu*angstrom^2'), symmetry=1, barrier=(30.0531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2968.28,'J/mol'), sigma=(5.18,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.80606,0.0102576,6.17291e-05,-9.01684e-08,3.59136e-11,11658.5,12.0619], Tmin=(100,'K'), Tmax=(946.037,'K')), NASAPolynomial(coeffs=[12.4692,0.0100558,-2.41229e-06,4.5713e-10,-3.93205e-14,8010.87,-43.6362], Tmin=(946.037,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.4553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), label="""butadiene13""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH]C=C(3743)',
    structure = SMILES('[CH2]C=C[CH2]'),
    E0 = (274.714,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100],'cm^-1')),
        HinderedRotor(inertia=(0.210055,'amu*angstrom^2'), symmetry=1, barrier=(25.2323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779472,'amu*angstrom^2'), symmetry=1, barrier=(93.4341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2985.34,'J/mol'), sigma=(5.28927,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=466.30 K, Pc=45.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.56316,0.022343,1.87062e-05,-3.93092e-08,1.63979e-11,33100.5,13.4098], Tmin=(100,'K'), Tmax=(974.267,'K')), NASAPolynomial(coeffs=[9.83,0.0151965,-5.22268e-06,9.67646e-10,-7.07852e-14,30607.7,-26.9852], Tmin=(974.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH]=CCC(=[CH])[O](10112)',
    structure = SMILES('[CH]=CCC(=[CH])[O]'),
    E0 = (506.188,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700,248.54],'cm^-1')),
        HinderedRotor(inertia=(0.294024,'amu*angstrom^2'), symmetry=1, barrier=(11.5673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.279141,'amu*angstrom^2'), symmetry=1, barrier=(11.6205,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44632,0.0499693,-3.99721e-05,1.06062e-08,1.20641e-12,60977.9,24.3537], Tmin=(100,'K'), Tmax=(965.042,'K')), NASAPolynomial(coeffs=[13.5736,0.0134005,-4.42235e-06,7.55587e-10,-5.18352e-14,57999.4,-37.0263], Tmin=(965.042,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(506.188,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC=C[CH2](10113)',
    structure = SMILES('[CH]=[C]CC=C[CH2]'),
    E0 = (689.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.9008,'amu*angstrom^2'), symmetry=1, barrier=(20.7112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0225369,'amu*angstrom^2'), symmetry=1, barrier=(7.12406,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.900655,'amu*angstrom^2'), symmetry=1, barrier=(20.7078,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38903,0.0502897,-3.27191e-05,8.29644e-09,-2.30726e-14,82986.9,25.1635], Tmin=(100,'K'), Tmax=(1150.67,'K')), NASAPolynomial(coeffs=[11.9984,0.022766,-9.03699e-06,1.6426e-09,-1.12829e-13,79925.9,-30.201], Tmin=(1150.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]C=CCC(=[CH])[O](10114)',
    structure = SMILES('[CH]C=CCC(=[CH])[O]'),
    E0 = (593.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3120,650,792.5,1650,373.022,373.029,373.061,373.225,373.278],'cm^-1')),
        HinderedRotor(inertia=(0.51927,'amu*angstrom^2'), symmetry=1, barrier=(51.2589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519496,'amu*angstrom^2'), symmetry=1, barrier=(51.2557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.518521,'amu*angstrom^2'), symmetry=1, barrier=(51.2563,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.802791,0.0627587,-4.61968e-05,1.43259e-08,-6.16529e-13,71533.4,29.1761], Tmin=(100,'K'), Tmax=(1053.57,'K')), NASAPolynomial(coeffs=[13.7149,0.0254061,-9.63111e-06,1.7012e-09,-1.15409e-13,68164.9,-36.8647], Tmin=(1053.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(593.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[C]=C([O])CC=C[CH2](10115)',
    structure = SMILES('[C]C(=O)CC=C[CH2]'),
    E0 = (636.942,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,180,615.859,615.867],'cm^-1')),
        HinderedRotor(inertia=(0.0728206,'amu*angstrom^2'), symmetry=1, barrier=(19.5988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0728168,'amu*angstrom^2'), symmetry=1, barrier=(19.5988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0728178,'amu*angstrom^2'), symmetry=1, barrier=(19.5988,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11978,0.0526703,-2.55487e-05,-3.87104e-09,4.3689e-12,76719,24.1471], Tmin=(100,'K'), Tmax=(1146.71,'K')), NASAPolynomial(coeffs=[15.5809,0.0216552,-1.03927e-05,2.09297e-09,-1.52555e-13,72125,-53.1692], Tmin=(1146.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(CJ3)"""),
)

species(
    label = '[CH2]C=CCC1=CO1(10116)',
    structure = SMILES('[CH2]C=CCC1=CO1'),
    E0 = (312.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03921,0.0503637,-1.53167e-06,-4.04194e-08,2.07442e-11,37647.7,25.3462], Tmin=(100,'K'), Tmax=(961.437,'K')), NASAPolynomial(coeffs=[17.5778,0.0172614,-5.59272e-06,1.02361e-09,-7.62066e-14,32817.3,-62.3744], Tmin=(961.437,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(312.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclopropene) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C1CC=CCO1(9596)',
    structure = SMILES('[CH]=C1CC=CCO1'),
    E0 = (170.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97453,0.018639,9.61451e-05,-1.38341e-07,5.34621e-11,20620.4,17.9183], Tmin=(100,'K'), Tmax=(973.958,'K')), NASAPolynomial(coeffs=[17.3355,0.020837,-7.78587e-06,1.62157e-09,-1.29999e-13,14531.7,-71.6805], Tmin=(973.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclohexane) + radical(Cds_P)"""),
)

species(
    label = '[O]C1=CCC=CC1(9599)',
    structure = SMILES('[O]C1=CCC=CC1'),
    E0 = (17.0638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86576,0.027947,5.57732e-05,-9.46644e-08,3.88509e-11,2146.41,21.8715], Tmin=(100,'K'), Tmax=(952.465,'K')), NASAPolynomial(coeffs=[15.2611,0.0193966,-5.88984e-06,1.0815e-09,-8.26406e-14,-2569.19,-53.4603], Tmin=(952.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(17.0638,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(1,4-Cyclohexadiene) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C(O)CC=C=C(9598)',
    structure = SMILES('[CH]=C(O)CC=C=C'),
    E0 = (261.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.540538,0.0657064,-4.56057e-05,3.17124e-09,6.04827e-12,31629.1,26.5604], Tmin=(100,'K'), Tmax=(962.114,'K')), NASAPolynomial(coeffs=[18.2894,0.0158054,-5.05255e-06,8.79666e-10,-6.24617e-14,27108.2,-64.1282], Tmin=(962.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.866,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Cds_P)"""),
)

species(
    label = 'C=C=CCC(=C)[O](5445)',
    structure = SMILES('C=C=CCC(=C)[O]'),
    E0 = (152.575,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055,235.182,235.325,235.397],'cm^-1')),
        HinderedRotor(inertia=(0.357469,'amu*angstrom^2'), symmetry=1, barrier=(14.1808,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.360245,'amu*angstrom^2'), symmetry=1, barrier=(14.1708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949278,0.0585744,-3.77193e-05,4.47375e-09,3.27005e-12,18467.8,26.1914], Tmin=(100,'K'), Tmax=(1008.29,'K')), NASAPolynomial(coeffs=[14.523,0.0216417,-7.94047e-06,1.42286e-09,-9.8863e-14,14870.7,-43.6699], Tmin=(1008.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C([O])[CH]C[CH][CH2](10117)',
    structure = SMILES('[CH]C([O])=CC[CH][CH2]'),
    E0 = (617.296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.181,0.0635775,-5.35687e-05,2.71687e-08,-6.03398e-12,74343.7,31.4418], Tmin=(100,'K'), Tmax=(1041.28,'K')), NASAPolynomial(coeffs=[8.04543,0.0372077,-1.55813e-05,2.84725e-09,-1.9453e-13,72914.2,-1.95265], Tmin=(1041.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(617.296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(RCCJC) + radical(RCCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])CC1[CH]C1(10118)',
    structure = SMILES('[CH]=C([O])CC1[CH]C1'),
    E0 = (485.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3434,0.0463232,-2.29489e-06,-3.27729e-08,1.66432e-11,58526.1,27.6843], Tmin=(100,'K'), Tmax=(968.33,'K')), NASAPolynomial(coeffs=[14.5721,0.0205333,-7.04376e-06,1.27039e-09,-9.11003e-14,54611.3,-42.6961], Tmin=(968.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(485.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(C=C(C)OJ) + radical(cyclopropane) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C[CH]C([CH2])O1(10052)',
    structure = SMILES('[CH]=C1C[CH]C([CH2])O1'),
    E0 = (471.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30304,0.0356884,5.37783e-05,-1.08864e-07,4.81921e-11,56807.1,25.3478], Tmin=(100,'K'), Tmax=(928.838,'K')), NASAPolynomial(coeffs=[21.7693,0.00820421,2.1322e-07,-1.15581e-10,-1.47576e-16,50388.7,-85.9637], Tmin=(928.838,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCJCO) + radical(CJC(C)OC) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1[CH]CC(=O)[CH]1(9974)',
    structure = SMILES('[CH2]C1[CH]CC([O])=C1'),
    E0 = (305.575,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83143,0.0327085,3.50915e-05,-7.08122e-08,3.03492e-11,36844,25.4419], Tmin=(100,'K'), Tmax=(941.017,'K')), NASAPolynomial(coeffs=[13.5268,0.0207751,-6.10933e-06,1.04159e-09,-7.4947e-14,32970.1,-39.1586], Tmin=(941.017,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(305.575,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(C=C(C)OJ) + radical(cyclopentene-4) + radical(Isobutyl)"""),
)

species(
    label = '[CH]=C1CC([CH][CH2])O1(10119)',
    structure = SMILES('[CH]=C1CC([CH][CH2])O1'),
    E0 = (536.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18651,0.0416987,3.44491e-05,-9.0158e-08,4.2808e-11,64682.8,24.539], Tmin=(100,'K'), Tmax=(905.914,'K')), NASAPolynomial(coeffs=[20.9332,0.00810955,1.31369e-06,-4.60562e-10,3.06226e-14,58905.6,-80.9157], Tmin=(905.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.804,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CCJCO) + radical(RCCJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2][CH]C1C=C([O])C1(10120)',
    structure = SMILES('[CH2][CH]C1[CH]C(=O)C1'),
    E0 = (389.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98901,0.0296535,3.8771e-05,-6.58992e-08,2.53878e-11,46917.5,26.9218], Tmin=(100,'K'), Tmax=(1008.78,'K')), NASAPolynomial(coeffs=[11.4653,0.0276979,-1.1286e-05,2.18449e-09,-1.59627e-13,43193.2,-27.8618], Tmin=(1008.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH]=C([O])CC=C=C(10121)',
    structure = SMILES('[CH]=C([O])CC=C=C'),
    E0 = (399.671,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,2950,3100,1380,975,1025,1650,540,610,2055,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.670823,'amu*angstrom^2'), symmetry=1, barrier=(15.4235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.671092,'amu*angstrom^2'), symmetry=1, barrier=(15.4297,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.71577,0.0639151,-6.04025e-05,2.91557e-08,-5.52464e-12,48194.6,27.4365], Tmin=(100,'K'), Tmax=(1289.24,'K')), NASAPolynomial(coeffs=[15.7953,0.0171291,-5.96801e-06,1.00756e-09,-6.63522e-14,44306.3,-49.1447], Tmin=(1289.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C#CCC=C[CH2](10122)',
    structure = SMILES('C#CC[CH]C=C'),
    E0 = (370.916,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,2175,525,750,770,3400,2100,180,778.188],'cm^-1')),
        HinderedRotor(inertia=(0.710671,'amu*angstrom^2'), symmetry=1, barrier=(16.3397,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208282,'amu*angstrom^2'), symmetry=1, barrier=(87.2726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.413451,'amu*angstrom^2'), symmetry=1, barrier=(87.1105,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (79.1198,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73011,0.0416497,-7.0109e-06,-1.98497e-08,1.06713e-11,44700,19.9241], Tmin=(100,'K'), Tmax=(969.88,'K')), NASAPolynomial(coeffs=[10.9115,0.0230749,-8.11932e-06,1.4205e-09,-9.76942e-14,42011.7,-28.7678], Tmin=(969.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(370.916,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CtCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Allyl_S)"""),
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
    label = '[CH]C=C(3739)',
    structure = SMILES('[CH]C=C'),
    E0 = (376.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,230.296,230.449,230.657],'cm^-1')),
        HinderedRotor(inertia=(1.34138,'amu*angstrom^2'), symmetry=1, barrier=(50.5187,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36193e-08,1.58213e-11,45331.5,10.6389], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35107e-06,1.16619e-09,-8.2762e-14,44095,-3.44607], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(376.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
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
    label = '[CH]=C([O])[CH]C=C[CH2](10123)',
    structure = SMILES('[CH]C([O])=CC=C[CH2]'),
    E0 = (433.856,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3100,440,815,1455,1000,324.534,324.584,324.585,324.598,324.609],'cm^-1')),
        HinderedRotor(inertia=(0.679395,'amu*angstrom^2'), symmetry=1, barrier=(50.7939,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.679424,'amu*angstrom^2'), symmetry=1, barrier=(50.794,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.679391,'amu*angstrom^2'), symmetry=1, barrier=(50.794,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00649,0.0541885,-1.13512e-05,-2.89508e-08,1.68159e-11,52299.3,26.5768], Tmin=(100,'K'), Tmax=(936.668,'K')), NASAPolynomial(coeffs=[15.2672,0.0223985,-7.05936e-06,1.17411e-09,-8.03082e-14,48350.8,-48.1066], Tmin=(936.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(433.856,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])C[C]=C[CH2](10124)',
    structure = SMILES('[CH]=C([O])C[C]=C[CH2]'),
    E0 = (612.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180,547.924],'cm^-1')),
        HinderedRotor(inertia=(0.0616746,'amu*angstrom^2'), symmetry=1, barrier=(12.962,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0859878,'amu*angstrom^2'), symmetry=1, barrier=(12.9543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0622176,'amu*angstrom^2'), symmetry=1, barrier=(12.9587,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738044,0.0653007,-6.40193e-05,3.22676e-08,-6.40407e-12,73778.5,29.0545], Tmin=(100,'K'), Tmax=(1229.97,'K')), NASAPolynomial(coeffs=[15.2223,0.0181961,-6.57321e-06,1.13072e-09,-7.52686e-14,70215.5,-43.8219], Tmin=(1229.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])CC=[C][CH2](10125)',
    structure = SMILES('[CH]=C([O])CC=[C][CH2]'),
    E0 = (612.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,543.639,545.112],'cm^-1')),
        HinderedRotor(inertia=(0.0615582,'amu*angstrom^2'), symmetry=1, barrier=(12.9568,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.563643,'amu*angstrom^2'), symmetry=1, barrier=(12.9593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0858142,'amu*angstrom^2'), symmetry=1, barrier=(12.9584,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.738044,0.0653007,-6.40193e-05,3.22676e-08,-6.40407e-12,73778.5,29.0545], Tmin=(100,'K'), Tmax=(1229.97,'K')), NASAPolynomial(coeffs=[15.2223,0.0181961,-6.57321e-06,1.13072e-09,-7.52686e-14,70215.5,-43.8219], Tmin=(1229.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(612.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])CC=[C]C(10126)',
    structure = SMILES('[CH]=C([O])CC=[C]C'),
    E0 = (460.908,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.989276,0.0648441,-6.23252e-05,3.2493e-08,-6.89231e-12,55544.1,27.8727], Tmin=(100,'K'), Tmax=(1131.19,'K')), NASAPolynomial(coeffs=[12.062,0.0256898,-1.04051e-05,1.89387e-09,-1.29697e-13,53039,-26.9117], Tmin=(1131.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])C[C]=CC(10127)',
    structure = SMILES('[CH]=C([O])C[C]=CC'),
    E0 = (460.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,1685,370,3120,650,792.5,1650,241.839,241.849],'cm^-1')),
        HinderedRotor(inertia=(0.193316,'amu*angstrom^2'), symmetry=1, barrier=(8.02505,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193327,'amu*angstrom^2'), symmetry=1, barrier=(8.0247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193316,'amu*angstrom^2'), symmetry=1, barrier=(8.02488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.989276,0.0648441,-6.23252e-05,3.2493e-08,-6.89231e-12,55544.1,27.8727], Tmin=(100,'K'), Tmax=(1131.19,'K')), NASAPolynomial(coeffs=[12.062,0.0256898,-1.04051e-05,1.89387e-09,-1.29697e-13,53039,-26.9117], Tmin=(1131.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C(O)[CH]C=C[CH2](10128)',
    structure = SMILES('[CH]C(O)=CC=C[CH2]'),
    E0 = (296.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648907,0.0580491,-3.368e-06,-4.67849e-08,2.52403e-11,35741.9,26.3598], Tmin=(100,'K'), Tmax=(928.548,'K')), NASAPolynomial(coeffs=[18.9749,0.019092,-5.03308e-06,7.89364e-10,-5.54553e-14,30614.7,-69.9768], Tmin=(928.548,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C=C[CH]C(=C)[O](5446)',
    structure = SMILES('[CH2]C=CC=C([CH2])[O]'),
    E0 = (222.089,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,445.341,448.892],'cm^-1')),
        HinderedRotor(inertia=(0.121161,'amu*angstrom^2'), symmetry=1, barrier=(16.8916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.118125,'amu*angstrom^2'), symmetry=1, barrier=(16.8475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.549124,'amu*angstrom^2'), symmetry=1, barrier=(77.0637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.959077,0.0549326,-1.53391e-05,-2.87516e-08,1.82559e-11,26831.7,26.3286], Tmin=(100,'K'), Tmax=(915.398,'K')), NASAPolynomial(coeffs=[17.0252,0.0163105,-3.80249e-06,5.35638e-10,-3.5634e-14,22567.1,-56.9883], Tmin=(915.398,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.089,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C(O)C[C]=C[CH2](10129)',
    structure = SMILES('[CH]=C(O)C[C]=C[CH2]'),
    E0 = (474.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,283.274],'cm^-1')),
        HinderedRotor(inertia=(0.298874,'amu*angstrom^2'), symmetry=1, barrier=(17.0178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298892,'amu*angstrom^2'), symmetry=1, barrier=(17.0177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298915,'amu*angstrom^2'), symmetry=1, barrier=(17.018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298913,'amu*angstrom^2'), symmetry=1, barrier=(17.0179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465919,0.0681937,-5.28863e-05,1.07718e-08,3.36321e-12,57217.4,28.5287], Tmin=(100,'K'), Tmax=(964.238,'K')), NASAPolynomial(coeffs=[18.1989,0.0160937,-5.22563e-06,9.0363e-10,-6.33282e-14,52799.9,-61.5504], Tmin=(964.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C=[C]CC(=C)[O](5447)',
    structure = SMILES('[CH2]C=[C]CC(=C)[O]'),
    E0 = (365.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,379.523,379.53,379.532],'cm^-1')),
        HinderedRotor(inertia=(0.00117034,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144568,'amu*angstrom^2'), symmetry=1, barrier=(14.777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144569,'amu*angstrom^2'), symmetry=1, barrier=(14.777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.863873,0.0611859,-4.54152e-05,1.25772e-08,3.89483e-13,44056.5,28.1986], Tmin=(100,'K'), Tmax=(1025.46,'K')), NASAPolynomial(coeffs=[14.5012,0.0218156,-8.04848e-06,1.43162e-09,-9.8479e-14,40532.7,-41.4807], Tmin=(1025.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])[CH]C=CC(10130)',
    structure = SMILES('[CH]C([O])=CC=CC'),
    E0 = (315.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.696147,0.0634623,-3.59876e-05,-1.0545e-09,5.88503e-12,38109.2,26.132], Tmin=(100,'K'), Tmax=(979.239,'K')), NASAPolynomial(coeffs=[14.7069,0.0261705,-9.40687e-06,1.64291e-09,-1.12242e-13,34409.2,-46.0489], Tmin=(979.239,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(315.801,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C(O)CC=[C][CH2](10131)',
    structure = SMILES('[CH]=C(O)CC=[C][CH2]'),
    E0 = (474.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,283.274],'cm^-1')),
        HinderedRotor(inertia=(0.298874,'amu*angstrom^2'), symmetry=1, barrier=(17.0178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298892,'amu*angstrom^2'), symmetry=1, barrier=(17.0177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298915,'amu*angstrom^2'), symmetry=1, barrier=(17.018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298913,'amu*angstrom^2'), symmetry=1, barrier=(17.0179,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.465919,0.0681937,-5.28863e-05,1.07718e-08,3.36321e-12,57217.4,28.5287], Tmin=(100,'K'), Tmax=(964.238,'K')), NASAPolynomial(coeffs=[18.1989,0.0160937,-5.22563e-06,9.0363e-10,-6.33282e-14,52799.9,-61.5504], Tmin=(964.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]=CCC(=C)[O](5448)',
    structure = SMILES('[CH2][C]=CCC(=C)[O]'),
    E0 = (365.311,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,1685,370,379.523,379.53,379.532],'cm^-1')),
        HinderedRotor(inertia=(0.00117034,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144568,'amu*angstrom^2'), symmetry=1, barrier=(14.777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.144569,'amu*angstrom^2'), symmetry=1, barrier=(14.777,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.863873,0.0611859,-4.54152e-05,1.25772e-08,3.89483e-13,44056.5,28.1986], Tmin=(100,'K'), Tmax=(1025.46,'K')), NASAPolynomial(coeffs=[14.5012,0.0218156,-8.04848e-06,1.43162e-09,-9.8479e-14,40532.7,-41.4807], Tmin=(1025.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(365.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])C([CH2])C=C(5412)',
    structure = SMILES('[CH]=C([O])C([CH2])C=C'),
    E0 = (430.605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,386.646,386.833],'cm^-1')),
        HinderedRotor(inertia=(0.0866711,'amu*angstrom^2'), symmetry=1, barrier=(9.19563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0865084,'amu*angstrom^2'), symmetry=1, barrier=(9.20059,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868164,'amu*angstrom^2'), symmetry=1, barrier=(9.2058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3817.52,'J/mol'), sigma=(6.40826,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=596.29 K, Pc=32.92 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.649929,0.071954,-7.97361e-05,4.80514e-08,-1.15285e-11,51911.9,28.1778], Tmin=(100,'K'), Tmax=(1019.49,'K')), NASAPolynomial(coeffs=[13.0047,0.0234797,-8.41478e-06,1.41289e-09,-9.17639e-14,49392.8,-31.6653], Tmin=(1019.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)CsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(Cds_P)"""),
)

species(
    label = 'C2H3(60)',
    structure = SMILES('[CH]=C'),
    E0 = (286.96,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,728.586,728.586,3521.64],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (27.0452,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1737.73,'J/mol'), sigma=(4.1,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.8307,-0.00147628,3.09002e-05,-3.80476e-08,1.43171e-11,34524.2,5.61949], Tmin=(100,'K'), Tmax=(933.662,'K')), NASAPolynomial(coeffs=[5.36086,0.00527345,-1.3196e-06,2.21564e-10,-1.68768e-14,33658.6,-4.76326], Tmin=(933.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.96,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), label="""C2H3""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]=C1CC(C=C)O1(10049)',
    structure = SMILES('[CH]=C1CC(C=C)O1'),
    E0 = (259.101,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30301,0.0352865,5.78179e-05,-1.16883e-07,5.26719e-11,31282.2,21.2205], Tmin=(100,'K'), Tmax=(911.922,'K')), NASAPolynomial(coeffs=[22.3383,0.00599984,2.3944e-06,-6.30451e-10,3.92556e-14,24826.9,-92.6817], Tmin=(911.922,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(Cds_P)"""),
)

species(
    label = 'C=CC1C=C([O])C1(9989)',
    structure = SMILES('C=CC1[CH]C(=O)C1'),
    E0 = (118.449,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07338,0.0220468,7.16335e-05,-1.07791e-07,4.22444e-11,14333.8,22.7831], Tmin=(100,'K'), Tmax=(966.974,'K')), NASAPolynomial(coeffs=[14.3586,0.022232,-7.77305e-06,1.50271e-09,-1.14843e-13,9573.39,-48.4031], Tmin=(966.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.449,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutanone) + radical(CCJC=O)"""),
)

species(
    label = '[CH]=C(O)C=CC=C(10132)',
    structure = SMILES('[CH]=C(O)C=CC=C'),
    E0 = (178.896,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.259972,0.0677715,-3.5563e-05,-2.01268e-08,1.79606e-11,21664.4,23.7016], Tmin=(100,'K'), Tmax=(916.46,'K')), NASAPolynomial(coeffs=[22.8556,0.00762484,-9.18129e-08,-1.2112e-10,7.16954e-15,15907.1,-92.1531], Tmin=(916.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_P)"""),
)

species(
    label = 'C=CC=CC(=C)[O](5460)',
    structure = SMILES('C=CC=CC(=C)[O]'),
    E0 = (69.6051,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.963786,'amu*angstrom^2'), symmetry=1, barrier=(22.1593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.970382,'amu*angstrom^2'), symmetry=1, barrier=(22.311,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689408,0.0604165,-2.70195e-05,-1.94496e-08,1.53259e-11,8502.2,23.2573], Tmin=(100,'K'), Tmax=(927.424,'K')), NASAPolynomial(coeffs=[18.89,0.0137934,-3.16863e-06,4.66234e-10,-3.28662e-14,3755.39,-70.5695], Tmin=(927.424,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.6051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH]=C1C[CH][CH]CO1(10095)',
    structure = SMILES('[CH]=C1C[CH][CH]CO1'),
    E0 = (446.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74144,0.0274467,6.46187e-05,-1.04035e-07,4.0965e-11,53801.8,23.213], Tmin=(100,'K'), Tmax=(987.675,'K')), NASAPolynomial(coeffs=[17.1905,0.0201202,-8.14968e-06,1.71019e-09,-1.35018e-13,48055.7,-64.7683], Tmin=(987.675,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Cyclohexanone) + radical(RCCJCC) + radical(CCJCO) + radical(Cds_P)"""),
)

species(
    label = 'O=C1[CH]C[CH][CH]C1(10026)',
    structure = SMILES('[O]C1=CC[CH][CH]C1'),
    E0 = (290.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.26562,0.0292527,2.27953e-05,-4.01958e-08,1.47788e-11,35056.5,24.1237], Tmin=(100,'K'), Tmax=(1041.2,'K')), NASAPolynomial(coeffs=[6.87779,0.032835,-1.30527e-05,2.40578e-09,-1.67926e-13,32941.5,-3.85819], Tmin=(1041.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(290.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclohexene) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(RCCJCC)"""),
)

species(
    label = '[CH]=C([O])C=CC=C(10133)',
    structure = SMILES('[CH]=C([O])C=CC=C'),
    E0 = (316.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00922,'amu*angstrom^2'), symmetry=1, barrier=(23.204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00976,'amu*angstrom^2'), symmetry=1, barrier=(23.2165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.139481,0.0726991,-7.37085e-05,3.61499e-08,-6.62889e-12,38255.1,26.6447], Tmin=(100,'K'), Tmax=(1536,'K')), NASAPolynomial(coeffs=[20.1187,0.00908879,-9.88883e-07,-1.31863e-11,5.83973e-15,33312.2,-75.6158], Tmin=(1536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(316.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C([O])[CH]CC=C(10064)',
    structure = SMILES('[CH]C([O])=CCC=C'),
    E0 = (346.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.902958,0.0589305,-2.85067e-05,-4.32973e-09,5.80809e-12,41812.5,28.4132], Tmin=(100,'K'), Tmax=(1015.76,'K')), NASAPolynomial(coeffs=[13.3108,0.0285114,-1.08197e-05,1.93605e-09,-1.33242e-13,38340.4,-36.3245], Tmin=(1015.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])CC[C]=C(10134)',
    structure = SMILES('[CH]=C([O])CC[C]=C'),
    E0 = (471.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,3120,650,792.5,1650,395.017,395.498,395.603],'cm^-1')),
        HinderedRotor(inertia=(0.0794573,'amu*angstrom^2'), symmetry=1, barrier=(8.80295,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0791258,'amu*angstrom^2'), symmetry=1, barrier=(8.80553,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0793304,'amu*angstrom^2'), symmetry=1, barrier=(8.79976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83364,0.0660244,-6.27298e-05,3.1608e-08,-6.39846e-12,56871.6,28.8416], Tmin=(100,'K'), Tmax=(1191.67,'K')), NASAPolynomial(coeffs=[13.5123,0.0234666,-9.16052e-06,1.6391e-09,-1.11287e-13,53849.9,-34.5489], Tmin=(1191.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=CCCC(=[CH])[O](10135)',
    structure = SMILES('[CH]=CCCC(=[CH])[O]'),
    E0 = (481.139,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3010,987.5,1337.5,450,1655,3115,3125,620,680,785,800,1600,1700,368.324,368.346],'cm^-1')),
        HinderedRotor(inertia=(0.108935,'amu*angstrom^2'), symmetry=1, barrier=(10.4543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108291,'amu*angstrom^2'), symmetry=1, barrier=(10.4562,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108047,'amu*angstrom^2'), symmetry=1, barrier=(10.457,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.577155,0.0674479,-6.31354e-05,3.0488e-08,-5.81241e-12,57997.4,29.694], Tmin=(100,'K'), Tmax=(1277.1,'K')), NASAPolynomial(coeffs=[15.8267,0.0196848,-7.03594e-06,1.2032e-09,-7.97291e-14,54102.4,-47.6063], Tmin=(1277.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(481.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_P) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C[CH]CC(=[CH])O(10136)',
    structure = SMILES('[CH]C=CCC(=[CH])O'),
    E0 = (455.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,3120,650,792.5,1650,216.875,216.878,216.882,216.916],'cm^-1')),
        HinderedRotor(inertia=(1.46077,'amu*angstrom^2'), symmetry=1, barrier=(48.7709,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46041,'amu*angstrom^2'), symmetry=1, barrier=(48.7696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46053,'amu*angstrom^2'), symmetry=1, barrier=(48.7704,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46158,'amu*angstrom^2'), symmetry=1, barrier=(48.7708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.496724,0.0660374,-3.63288e-05,-5.67522e-09,8.58427e-12,54973.7,28.7729], Tmin=(100,'K'), Tmax=(969.291,'K')), NASAPolynomial(coeffs=[17.0437,0.0227299,-7.96264e-06,1.4e-09,-9.74249e-14,50592.6,-56.5932], Tmin=(969.291,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(455.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C=CCC(=C)[O](5432)',
    structure = SMILES('[CH]C=CCC(=C)[O]'),
    E0 = (346.655,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,407.297,408.552,408.907,409.003,409.347,410.514],'cm^-1')),
        HinderedRotor(inertia=(0.43289,'amu*angstrom^2'), symmetry=1, barrier=(51.3916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.437734,'amu*angstrom^2'), symmetry=1, barrier=(51.4327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.434184,'amu*angstrom^2'), symmetry=1, barrier=(51.33,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.902958,0.0589305,-2.85067e-05,-4.32973e-09,5.80809e-12,41812.5,28.4132], Tmin=(100,'K'), Tmax=(1015.76,'K')), NASAPolynomial(coeffs=[13.3108,0.0285114,-1.08197e-05,1.93605e-09,-1.33242e-13,38340.4,-36.3245], Tmin=(1015.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(346.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]1C[CH]C([CH2])O1(10137)',
    structure = SMILES('[CH][C]1C[CH]C([CH2])O1'),
    E0 = (767.923,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.554759,0.0628229,-5.52312e-05,2.69925e-08,-5.11859e-12,92494.9,26.8772], Tmin=(100,'K'), Tmax=(1488.87,'K')), NASAPolynomial(coeffs=[12.8426,0.0212436,-4.7102e-06,5.06268e-10,-2.2313e-14,89785.4,-34.1067], Tmin=(1488.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(767.923,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CCJCO) + radical(C2CsJOCs) + radical(CJC(C)OC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([O])[CH]C=C[CH2](7208)',
    structure = SMILES('[CH]C([O])[CH]C=C[CH2]'),
    E0 = (698.194,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.34813,0.0707347,-5.76647e-05,2.01017e-08,-1.6901e-12,84113,27.8888], Tmin=(100,'K'), Tmax=(1083.68,'K')), NASAPolynomial(coeffs=[17.3773,0.0215826,-8.59941e-06,1.58735e-09,-1.11117e-13,79617.4,-59.3485], Tmin=(1083.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(698.194,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([O])C[C]=C[CH2](7210)',
    structure = SMILES('[CH]C([O])C[C]=C[CH2]'),
    E0 = (819.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.809651,0.0681515,-6.35438e-05,3.12201e-08,-6.22252e-12,98633.8,29.5054], Tmin=(100,'K'), Tmax=(1198.81,'K')), NASAPolynomial(coeffs=[13.3652,0.026258,-1.11248e-05,2.06946e-09,-1.43415e-13,95623.4,-33.3445], Tmin=(1198.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(819.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C](O)C[C]=C[CH2](10138)',
    structure = SMILES('[CH][C](O)C[C]=C[CH2]'),
    E0 = (765.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.546593,0.0801363,-0.000103507,7.42769e-08,-2.15988e-11,92175.2,30.0742], Tmin=(100,'K'), Tmax=(837.338,'K')), NASAPolynomial(coeffs=[11.0841,0.029804,-1.33526e-05,2.5061e-09,-1.72934e-13,90410.3,-18.8935], Tmin=(837.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(765.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C2CsJOH) + radical(Allyl_P) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([O])CC=[C][CH2](7211)',
    structure = SMILES('[CH]C([O])CC=[C][CH2]'),
    E0 = (819.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.809651,0.0681515,-6.35438e-05,3.12201e-08,-6.22252e-12,98633.8,29.5054], Tmin=(100,'K'), Tmax=(1198.81,'K')), NASAPolynomial(coeffs=[13.3652,0.026258,-1.11248e-05,2.06946e-09,-1.43415e-13,95623.4,-33.3445], Tmin=(1198.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(819.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(Allyl_P) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C](O)CC=[C][CH2](10139)',
    structure = SMILES('[CH][C](O)CC=[C][CH2]'),
    E0 = (765.386,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.546593,0.0801363,-0.000103507,7.42769e-08,-2.15988e-11,92175.2,30.0742], Tmin=(100,'K'), Tmax=(837.338,'K')), NASAPolynomial(coeffs=[11.0841,0.029804,-1.33526e-05,2.5061e-09,-1.72934e-13,90410.3,-18.8935], Tmin=(837.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(765.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C2CsJOH) + radical(Allyl_P) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([O])CC=CC1(10140)',
    structure = SMILES('[CH]C1([O])CC=CC1'),
    E0 = (484.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41732,0.040335,2.40284e-05,-6.27639e-08,2.74391e-11,58427.1,22.6772], Tmin=(100,'K'), Tmax=(969.764,'K')), NASAPolynomial(coeffs=[16.0285,0.0202504,-7.05826e-06,1.33383e-09,-9.9835e-14,53703.7,-57.1067], Tmin=(969.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CC(C)2OJ) + radical(CCJ2_triplet)"""),
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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09172,0.0173333,-8.20209e-06,-2.63358e-09,2.66049e-12,40755.9,8.10965], Tmin=(100,'K'), Tmax=(946.054,'K')), NASAPolynomial(coeffs=[6.98214,0.0072721,-2.37773e-06,3.99152e-10,-2.69331e-14,39733.9,-11.9544], Tmin=(946.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C3H3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]C(=C)O(10141)',
    structure = SMILES('[CH]C(=C)O'),
    E0 = (184.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.38644,0.027409,4.34284e-06,-2.73997e-08,1.32371e-11,22240.4,14.3828], Tmin=(100,'K'), Tmax=(947.673,'K')), NASAPolynomial(coeffs=[10.6449,0.012572,-3.86221e-06,6.65059e-10,-4.73709e-14,19776.1,-29.7587], Tmin=(947.673,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.373,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)OC=C[CH2](10142)',
    structure = SMILES('[CH]C(=C)O[CH]C=C'),
    E0 = (363.884,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.763411,0.0611043,-3.04454e-05,-2.47598e-09,4.75808e-12,43890.3,27.4548], Tmin=(100,'K'), Tmax=(1065.58,'K')), NASAPolynomial(coeffs=[14.1047,0.0294859,-1.19259e-05,2.1972e-09,-1.53056e-13,39998.9,-42.6749], Tmin=(1065.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(363.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CCJ(O)C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1([O])CC1C=C(10089)',
    structure = SMILES('[CH]C1([O])CC1C=C'),
    E0 = (577.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07815,0.0484951,5.81631e-06,-4.90675e-08,2.40852e-11,69542.9,24.8601], Tmin=(100,'K'), Tmax=(955.717,'K')), NASAPolynomial(coeffs=[17.884,0.0168281,-5.177e-06,9.39001e-10,-7.0659e-14,64564.5,-64.6971], Tmin=(955.717,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.219,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CCJ2_triplet)"""),
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
    E0 = (374.566,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (943.99,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1208.51,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (805.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (848.734,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (377.495,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (382.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (382.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (399.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (399.539,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (640.157,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (600.502,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (471.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (412.64,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (536.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (428.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (623.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (630.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (485.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (698.828,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (699.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (645.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (824.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (824.2,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (582.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (622.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (544.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (519.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (575.78,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (418.874,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (492.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (518.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (407.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (590.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (944.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (382.474,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (382.097,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (452.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (452.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (446.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (386.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (536.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (526.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (490.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (629.563,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (626.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (607.467,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (565.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (824.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (721.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (908.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (804.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (844.092,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (790.359,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (484.894,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (631.975,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (677.684,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (577.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['HCCO(2227)', 'butadiene13(1350)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(20)', '[CH]=CCC(=[CH])[O](10112)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH]=[C]CC=C[CH2](10113)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]C=CCC(=[CH])[O](10114)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[C]=C([O])CC=C[CH2](10115)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH2]C=CCC1=CO1(10116)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;O_rad;CdsinglepriH_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C1CC=CCO1(9596)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), Tmin=(550,'K'), Tmax=(650,'K'), comment="""Estimated using template [R6_SSSDS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6_SSSDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[O]C1=CCC=CC1(9599)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.53377e+11,'s^-1'), n=0.0685, Ea=(8.20064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad_out;Cpri_rad_out_2H] for rate rule [R6;CdsingleH_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C(O)CC=C=C(9598)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['C=C=CCC(=C)[O](5445)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C([O])[CH]C[CH][CH2](10117)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C([O])CC1[CH]C1(10118)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.05e+08,'s^-1'), n=1.192, Ea=(225.936,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra_cs2H] for rate rule [R3_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C1C[CH]C([CH2])O1(10052)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.66591e+07,'s^-1'), n=1.01661, Ea=(96.7661,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra_pri;radadd_intra] for rate rule [R5_SS_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 94.0 to 96.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH2]C1[CH]CC(=O)[CH]1(9974)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.29e+09,'s^-1'), n=0.62, Ea=(38.0744,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra_pri;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C1CC([CH][CH2])O1(10119)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(162.238,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 160.2 to 162.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH2][CH]C1C=C([O])C1(10120)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(53.5552,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R5_DS_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH]=C([O])CC=C=C(10121)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.42e+08,'cm^3/(mol*s)'), n=1.64, Ea=(11.7989,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2713 used for Ca_Cds-HH;HJ
Exact match found for rate rule [Ca_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['O(4)', 'C#CCC=C[CH2](10122)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(18.899,'m^3/(mol*s)'), n=1.76329, Ea=(16.1554,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cs_Ct-H;YJ] for rate rule [Ct-Cs_Ct-H;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['HCCO(2227)', '[CH2][CH]C=C(3743)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.000357827,'m^3/(mol*s)'), n=2.74787, Ea=(45.8259,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CsJ-CdHH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=C([CH2])[O](9653)', '[CH]C=C(3739)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.36745e+07,'m^3/(mol*s)'), n=-0.263863, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00481396807501, var=0.0768145972539, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.567716674236
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH]C=C(3743)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.7349e+07,'m^3/(mol*s)'), n=-0.263863, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00481396807501, var=0.0768145972539, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.567716674236
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH]=C([O])[CH]C=C[CH2](10123)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH]=C([O])C[C]=C[CH2](10124)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH]=C([O])CC=[C][CH2](10125)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C([O])CC=[C]C(10126)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.63e+08,'s^-1'), n=1.73, Ea=(207.731,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 123 used for R2H_S;C_rad_out_2H;Cd_H_out_doubleC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cd_H_out_doubleC]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH]=C([O])C[C]=CC(10127)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 198 used for R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C(O)[CH]C=C[CH2](10128)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(587605,'s^-1'), n=2.09, Ea=(169.87,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_H/Cd] + [R3H_SS_2Cd;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R3H_SS_2Cd;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH2]C=C[CH]C(=C)[O](5446)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C(O)C[C]=C[CH2](10129)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSS;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH2]C=[C]CC(=C)[O](5447)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cd_H_out_doubleC]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C([O])[CH]C=CC(10130)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(256000,'s^-1'), n=2, Ea=(117.57,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;C_rad_out_2H;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C(O)CC=[C][CH2](10131)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(380071,'s^-1'), n=1.62386, Ea=(44.2978,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_RSSR;Y_rad_out;XH_out] for rate rule [R5H_DSSS;Cd_rad_out;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH2][C]=CCC(=C)[O](5448)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSSD;Cd_rad_out_singleH;XH_out]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=C([O])C([CH2])C=C(5412)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C2H3(60)', '[CH]CC(=[CH])[O](9655)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C1CC(C=C)O1(10049)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.8e+12,'s^-1'), n=-0.1525, Ea=(7.90776,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_H/OneDe;Ypri_rad_out] + [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/OneDe;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['C=CC1C=C([O])C1(9989)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_H/OneDe;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_H/OneDe;CdsinglepriH_rad_out]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C(O)C=CC=C(10132)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['C=CC=CC(=C)[O](5460)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C1C[CH][CH]CO1(10095)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(9.91671e+09,'s^-1'), n=0.30082, Ea=(71.9272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri_2H;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 67.9 to 71.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['O=C1[CH]C[CH][CH]C1(10026)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.16959e+10,'s^-1'), n=0.31, Ea=(12.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6_linear;doublebond_intra_pri_2H;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['H(3)', '[CH]=C([O])C=CC=C(10133)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(105.863,'m^3/(mol*s)'), n=1.66278, Ea=(8.08939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=[C][O](6861)', 'butadiene13(1350)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.0534234,'m^3/(mol*s)'), n=2.459, Ea=(4.91982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-CdH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]=C([O])[CH]CC=C(10064)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.169e+11,'s^-1'), n=0.707, Ea=(116.068,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/OneDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=C([O])CC[C]=C(10134)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=CCCC(=[CH])[O](10135)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 194 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=C[CH]CC(=[CH])O(10136)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_2;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]C=CCC(=C)[O](5432)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.456e+11,'s^-1'), n=0.86, Ea=(191.209,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cd_H_out_singleH] for rate rule [R6HJ_3;Cd_rad_out_singleH;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH][C]1C[CH]C([CH2])O1(10137)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(6.43608e+09,'s^-1'), n=0.0758676, Ea=(56.4793,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 1 used for R5JJ
Exact match found for rate rule [R5JJ]
Euclidian distance = 0
family: 1,4_Cyclic_birad_scission"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]C([O])[CH]C=C[CH2](7208)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH]C([O])C[C]=C[CH2](7210)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2.6374e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH][C](O)C[C]=C[CH2](10138)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]C([O])CC=[C][CH2](7211)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH][C](O)CC=[C][CH2](10139)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]C1([O])CC=CC1(10140)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(3.49749e+08,'s^-1'), n=0.656505, Ea=(110.328,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_SMS;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_SMS_CO;carbonylbond_intra;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 107.9 to 110.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['C3H3(5450)', '[CH]C(=C)O(10141)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(3.59298e+09,'s^-1'), n=0.614601, Ea=(257.409,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.07280803979615884, var=3.9519504359180297, Tref=1000.0, N=5, correlation='Root_1R!H->C_2R!H->C_5R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R!H->C_2R!H->C_5R!H->O
    Total Standard Deviation in ln(k): 4.168250509264941
Exact match found for rate rule [Root_1R!H->C_2R!H->C_5R!H->O]
Euclidian distance = 0
family: Retroene"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH]C(=C)OC=C[CH2](10142)'],
    products = ['[CH]=C([O])CC=C[CH2](5449)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH]=C([O])CC=C[CH2](5449)'],
    products = ['[CH]C1([O])CC1C=C(10089)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(202.654,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCd] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_csHCd]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 201.4 to 202.7 kJ/mol to match endothermicity of reaction."""),
)

network(
    label = 'PDepNetwork #2462',
    isomers = [
        '[CH]=C([O])CC=C[CH2](5449)',
    ],
    reactants = [
        ('HCCO(2227)', 'butadiene13(1350)'),
        ('HCCO(2227)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2462',
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


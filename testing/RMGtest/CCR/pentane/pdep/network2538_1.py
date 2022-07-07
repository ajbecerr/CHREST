species(
    label = '[CH2][CH][CH]C[CH][C]=O(9128)',
    structure = SMILES('[CH2][CH][CH]C[CH][C]=O'),
    E0 = (646.977,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,1855,455,950,401.252,3128.77,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0188754,'amu*angstrom^2'), symmetry=1, barrier=(3.54249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38728,0.0671154,-0.000101705,9.70352e-08,-3.57488e-11,77898.2,33.1678], Tmin=(100,'K'), Tmax=(870.5,'K')), NASAPolynomial(coeffs=[1.93353,0.0418045,-1.88017e-05,3.45454e-09,-2.31828e-13,78667,35.5702], Tmin=(870.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(646.977,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(CCJCHO) + radical(RCCJ) + radical(CCCJ=O)"""),
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
    label = '[CH2][CH][C]=O(9632)',
    structure = SMILES('[CH2][CH][C]=O'),
    E0 = (334.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.0670133,'amu*angstrom^2'), symmetry=1, barrier=(7.51319,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00540518,'amu*angstrom^2'), symmetry=1, barrier=(7.51488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.78792,0.0268854,-2.93724e-05,1.83573e-08,-4.61895e-12,40303.1,15.8795], Tmin=(100,'K'), Tmax=(968.231,'K')), NASAPolynomial(coeffs=[6.80976,0.01027,-3.63122e-06,6.33225e-10,-4.24838e-14,39524.3,-3.39371], Tmin=(968.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(334.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJCHO) + radical(CJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[CH][CH][CH2](5537)',
    structure = SMILES('[CH][CH][CH2]'),
    E0 = (727.784,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1609,1609.01,1609.02],'cm^-1')),
        HinderedRotor(inertia=(0.0337841,'amu*angstrom^2'), symmetry=1, barrier=(5.78124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00314687,'amu*angstrom^2'), symmetry=1, barrier=(5.78111,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.42424,0.0149081,-5.18362e-06,2.15834e-10,1.32851e-13,87550.9,15.303], Tmin=(100,'K'), Tmax=(1973.84,'K')), NASAPolynomial(coeffs=[8.10291,0.00913487,-3.61425e-06,6.37549e-10,-4.11102e-14,84981.5,-12.2801], Tmin=(1973.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(727.784,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
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
    label = '[CH2][CH][CH][CH2](5531)',
    structure = SMILES('[CH2][CH][CH][CH2]'),
    E0 = (655.354,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1802.64],'cm^-1')),
        HinderedRotor(inertia=(0.00215831,'amu*angstrom^2'), symmetry=1, barrier=(4.96293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00215621,'amu*angstrom^2'), symmetry=1, barrier=(4.95965,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00214837,'amu*angstrom^2'), symmetry=1, barrier=(4.95028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.82147,0.0344596,-5.62911e-05,6.39381e-08,-2.62695e-11,78855.1,21.3291], Tmin=(100,'K'), Tmax=(865.068,'K')), NASAPolynomial(coeffs=[-1.36886,0.0321046,-1.45274e-05,2.71435e-09,-1.84229e-13,80393.2,45.6373], Tmin=(865.068,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.354,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
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
    label = '[CH2][CH][C]C[CH][C]=O(13229)',
    structure = SMILES('[CH2][CH][C]C[CH][C]=O'),
    E0 = (900.734,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06839,0.0698332,-9.76984e-05,8.10129e-08,-2.6955e-11,108434,31.6549], Tmin=(100,'K'), Tmax=(854.211,'K')), NASAPolynomial(coeffs=[7.20516,0.0315857,-1.38341e-05,2.5267e-09,-1.6983e-13,107732,5.04669], Tmin=(854.211,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(900.734,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH][CH]C[C][C]=O(13230)',
    structure = SMILES('[CH2][CH][CH]C[C][C]=O'),
    E0 = (927.663,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,210.722,884.344,1031.34,3525.62],'cm^-1')),
        HinderedRotor(inertia=(0.0555416,'amu*angstrom^2'), symmetry=1, barrier=(1.64557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555416,'amu*angstrom^2'), symmetry=1, barrier=(1.64557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555416,'amu*angstrom^2'), symmetry=1, barrier=(1.64557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555416,'amu*angstrom^2'), symmetry=1, barrier=(1.64557,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0555416,'amu*angstrom^2'), symmetry=1, barrier=(1.64557,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17188,0.0765211,-0.000138641,1.39121e-07,-5.17185e-11,111660,34.0849], Tmin=(100,'K'), Tmax=(877.607,'K')), NASAPolynomial(coeffs=[1.09496,0.041995,-2.00185e-05,3.72724e-09,-2.50031e-13,113017,42.098], Tmin=(877.607,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(927.663,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C][CH]C[CH][C]=O(13231)',
    structure = SMILES('[CH2][C][CH]C[CH][C]=O'),
    E0 = (900.746,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1855,455,950,200.06,800.04,1200.17,1599.86],'cm^-1')),
        HinderedRotor(inertia=(0.156085,'amu*angstrom^2'), symmetry=1, barrier=(3.58871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156085,'amu*angstrom^2'), symmetry=1, barrier=(3.58871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156085,'amu*angstrom^2'), symmetry=1, barrier=(3.58871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156085,'amu*angstrom^2'), symmetry=1, barrier=(3.58871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156085,'amu*angstrom^2'), symmetry=1, barrier=(3.58871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10853,0.0702741,-0.000102765,8.9232e-08,-3.0726e-11,108432,31.339], Tmin=(100,'K'), Tmax=(856.857,'K')), NASAPolynomial(coeffs=[6.11521,0.0335324,-1.50412e-05,2.77025e-09,-1.86716e-13,108065,10.8221], Tmin=(856.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(900.746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(CCJCHO) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH][CH][CH]C[CH][C]=O(13232)',
    structure = SMILES('[CH][CH][CH]C[CH][C]=O'),
    E0 = (889.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,1855,455,950,199.144,775.293,1034.98,1780.83,2436.87,4000],'cm^-1')),
        HinderedRotor(inertia=(0.116749,'amu*angstrom^2'), symmetry=1, barrier=(2.70552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116749,'amu*angstrom^2'), symmetry=1, barrier=(2.70552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116749,'amu*angstrom^2'), symmetry=1, barrier=(2.70552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116749,'amu*angstrom^2'), symmetry=1, barrier=(2.70552,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.116749,'amu*angstrom^2'), symmetry=1, barrier=(2.70552,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28821,0.0695447,-0.000111114,1.04832e-07,-3.77985e-11,107124,32.7411], Tmin=(100,'K'), Tmax=(878.267,'K')), NASAPolynomial(coeffs=[3.08827,0.037931,-1.71292e-05,3.13432e-09,-2.09119e-13,107711,29.4317], Tmin=(878.267,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(889.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(CCJCHO) + radical(RCCJC) + radical(CCCJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C1CC1[C]=O(10011)',
    structure = SMILES('[CH2][CH]C1CC1[C]=O'),
    E0 = (431.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5886,0.0432148,-4.93183e-06,-2.08235e-08,9.88487e-12,52005.7,28.7249], Tmin=(100,'K'), Tmax=(1048.08,'K')), NASAPolynomial(coeffs=[11.6117,0.0260685,-1.06e-05,1.99666e-09,-1.41734e-13,48745.4,-25.6319], Tmin=(1048.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.608,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(Cs_S) + radical(RCCJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=CC[CH][CH][C]=O(9967)',
    structure = SMILES('C=CC[CH]C=[C][O]'),
    E0 = (377.227,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,626.551,626.552,626.559,626.569],'cm^-1')),
        HinderedRotor(inertia=(0.154469,'amu*angstrom^2'), symmetry=1, barrier=(3.55155,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0750306,'amu*angstrom^2'), symmetry=1, barrier=(20.9024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.909118,'amu*angstrom^2'), symmetry=1, barrier=(20.9024,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27843,0.0493894,-1.33895e-05,-1.7535e-08,1.00144e-11,45477,28.1022], Tmin=(100,'K'), Tmax=(1016.14,'K')), NASAPolynomial(coeffs=[13.486,0.0239851,-9.32375e-06,1.73372e-09,-1.23164e-13,41826.8,-36.7416], Tmin=(1016.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.227,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=C[CH]C[C]=O(5328)',
    structure = SMILES('[CH2]C=C[CH]C[C]=O'),
    E0 = (352.235,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950,354.573,1044.56],'cm^-1')),
        HinderedRotor(inertia=(0.0840196,'amu*angstrom^2'), symmetry=1, barrier=(7.49586,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0840197,'amu*angstrom^2'), symmetry=1, barrier=(7.49587,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0840179,'amu*angstrom^2'), symmetry=1, barrier=(7.49591,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0264725,'amu*angstrom^2'), symmetry=1, barrier=(20.497,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22306,0.0574204,-4.48442e-05,1.84806e-08,-3.13231e-12,42466.9,29.7248], Tmin=(100,'K'), Tmax=(1378.77,'K')), NASAPolynomial(coeffs=[12.0671,0.0259601,-1.06173e-05,1.93092e-09,-1.31482e-13,39476.7,-26.0742], Tmin=(1378.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.235,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C=C[CH][CH][C]=O(10015)',
    structure = SMILES('[CH2]C=C[CH]C=[C][O]'),
    E0 = (478.362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,602.609,602.83,602.961],'cm^-1')),
        HinderedRotor(inertia=(0.291655,'amu*angstrom^2'), symmetry=1, barrier=(75.1758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106574,'amu*angstrom^2'), symmetry=1, barrier=(27.4797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106612,'amu*angstrom^2'), symmetry=1, barrier=(27.4803,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49026,0.039342,2.13071e-05,-6.21794e-08,2.85287e-11,57638.6,27.9866], Tmin=(100,'K'), Tmax=(942.386,'K')), NASAPolynomial(coeffs=[16.4744,0.0155637,-4.23045e-06,7.26991e-10,-5.48589e-14,53046.1,-52.7961], Tmin=(942.386,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(478.362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJC=C) + radical(Allyl_P) + radical(C=CJO)"""),
)

species(
    label = 'CH2CHCO(3668)',
    structure = SMILES('C=C[C]=O'),
    E0 = (83.3963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.27992,'amu*angstrom^2'), symmetry=1, barrier=(29.4278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0553,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.08334,0.0153205,6.54259e-06,-1.77535e-08,7.39369e-12,10067.5,11.8951], Tmin=(100,'K'), Tmax=(1002.99,'K')), NASAPolynomial(coeffs=[6.9784,0.0111885,-4.32948e-06,8.06749e-10,-5.75264e-14,8712.67,-9.76694], Tmin=(1002.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(83.3963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH2CHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2][CH][CH][CH][CH][C]=O(13233)',
    structure = SMILES('[CH2][CH][CH][CH]C=[C][O]'),
    E0 = (843.595,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1501.37,1501.37,1501.42,1501.43],'cm^-1')),
        HinderedRotor(inertia=(2.62164,'amu*angstrom^2'), symmetry=1, barrier=(60.2766,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00162778,'amu*angstrom^2'), symmetry=1, barrier=(2.60414,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113218,'amu*angstrom^2'), symmetry=1, barrier=(2.6031,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113238,'amu*angstrom^2'), symmetry=1, barrier=(2.60357,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67837,0.0560464,-6.81821e-05,5.85494e-08,-2.13289e-11,101540,33.9716], Tmin=(100,'K'), Tmax=(810.31,'K')), NASAPolynomial(coeffs=[3.86227,0.0367268,-1.66118e-05,3.116e-09,-2.14077e-13,101466,25.6247], Tmin=(810.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(843.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJCC) + radical(Allyl_S) + radical(RCCJC) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH]C[CH][CH][C]=O(10008)',
    structure = SMILES('[CH2][CH]C[CH]C=[C][O]'),
    E0 = (649.137,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,320.938,1161.84,1161.86,1161.87],'cm^-1')),
        HinderedRotor(inertia=(0.00535115,'amu*angstrom^2'), symmetry=1, barrier=(48.6604,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.28324,'amu*angstrom^2'), symmetry=1, barrier=(93.7698,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00163705,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0699164,'amu*angstrom^2'), symmetry=1, barrier=(5.10982,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52939,0.054353,-3.98316e-05,1.59854e-08,-2.74662e-12,78161.9,31.246], Tmin=(100,'K'), Tmax=(1315.24,'K')), NASAPolynomial(coeffs=[9.25591,0.0308549,-1.30328e-05,2.40188e-09,-1.64692e-13,76129.4,-8.14725], Tmin=(1315.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJC) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH][CH][CH]C[C]=O(9127)',
    structure = SMILES('[CH2][CH][CH][CH]C[C]=O'),
    E0 = (679.35,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,1855,455,950,408.216,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00941144,'amu*angstrom^2'), symmetry=1, barrier=(1.73502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00941144,'amu*angstrom^2'), symmetry=1, barrier=(1.73502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00941144,'amu*angstrom^2'), symmetry=1, barrier=(1.73502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00941144,'amu*angstrom^2'), symmetry=1, barrier=(1.73502,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00941144,'amu*angstrom^2'), symmetry=1, barrier=(1.73502,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41942,0.0682477,-0.000110508,1.10021e-07,-4.15626e-11,81788.8,35.7874], Tmin=(100,'K'), Tmax=(865.841,'K')), NASAPolynomial(coeffs=[0.724014,0.0438452,-2.03916e-05,3.79876e-09,-2.56651e-13,82944.3,45.0197], Tmin=(865.841,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(679.35,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C[CH][CH][CH][C]=O(13234)',
    structure = SMILES('[CH2]C[CH][CH]C=[C][O]'),
    E0 = (649.149,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,180,846.957,847.001,847.023],'cm^-1')),
        HinderedRotor(inertia=(0.00547899,'amu*angstrom^2'), symmetry=1, barrier=(2.78928,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00547945,'amu*angstrom^2'), symmetry=1, barrier=(2.78952,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121308,'amu*angstrom^2'), symmetry=1, barrier=(2.78912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121312,'amu*angstrom^2'), symmetry=1, barrier=(2.7892,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70795,0.0531128,-3.88276e-05,1.61061e-08,-2.97124e-12,78154.7,30.4374], Tmin=(100,'K'), Tmax=(1197.26,'K')), NASAPolynomial(coeffs=[7.41184,0.0340563,-1.49525e-05,2.81186e-09,-1.95265e-13,76788.9,1.89252], Tmin=(1197.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJCC) + radical(Allyl_S) + radical(RCCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][CH][CH][CH][CH]C=O(13235)',
    structure = SMILES('[CH2][CH][CH][CH]C=C[O]'),
    E0 = (603.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53737,0.0534881,-3.83162e-05,1.51275e-08,-2.56009e-12,72715.5,31.3459], Tmin=(100,'K'), Tmax=(1334.16,'K')), NASAPolynomial(coeffs=[9.15259,0.0306564,-1.26462e-05,2.30031e-09,-1.56466e-13,70683.6,-7.58859], Tmin=(1334.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(603.851,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_S) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH][CH][C]=O(10317)',
    structure = SMILES('C[CH][CH][CH]C=[C][O]'),
    E0 = (638.349,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3010,987.5,1337.5,450,1655,1685,370,180,1115.72,1115.78,1115.85],'cm^-1')),
        HinderedRotor(inertia=(0.181327,'amu*angstrom^2'), symmetry=1, barrier=(4.16907,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00472112,'amu*angstrom^2'), symmetry=1, barrier=(4.17203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18076,'amu*angstrom^2'), symmetry=1, barrier=(4.15602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00473836,'amu*angstrom^2'), symmetry=1, barrier=(4.16996,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69464,0.0546385,-5.51821e-05,4.2767e-08,-1.53176e-11,76854.9,31.4366], Tmin=(100,'K'), Tmax=(764.27,'K')), NASAPolynomial(coeffs=[3.86453,0.0393853,-1.75976e-05,3.31129e-09,-2.2914e-13,76637.1,22.2959], Tmin=(764.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(638.349,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(RCCJCC) + radical(Allyl_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=CC[CH][C]=O(5329)',
    structure = SMILES('[CH2]C=CC[CH][C]=O'),
    E0 = (319.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.778956,'amu*angstrom^2'), symmetry=1, barrier=(17.9097,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.403324,'amu*angstrom^2'), symmetry=1, barrier=(50.6156,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.779054,'amu*angstrom^2'), symmetry=1, barrier=(17.912,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0968321,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3657.02,'J/mol'), sigma=(6.16401,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=571.22 K, Pc=35.43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.907921,0.0595326,-4.70532e-05,1.94918e-08,-3.25348e-12,38588.9,28.1278], Tmin=(100,'K'), Tmax=(1429.35,'K')), NASAPolynomial(coeffs=[14.097,0.0226237,-8.32024e-06,1.42646e-09,-9.37883e-14,34818.5,-40.2132], Tmin=(1429.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Allyl_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C([CH2])[CH][C]=O(9158)',
    structure = SMILES('[CH2][CH]C([CH2])[CH][C]=O'),
    E0 = (652.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1855,455,950,944.543,944.577],'cm^-1')),
        HinderedRotor(inertia=(0.00506013,'amu*angstrom^2'), symmetry=1, barrier=(3.2039,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139316,'amu*angstrom^2'), symmetry=1, barrier=(3.20315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00506097,'amu*angstrom^2'), symmetry=1, barrier=(3.20391,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139369,'amu*angstrom^2'), symmetry=1, barrier=(3.20436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.139367,'amu*angstrom^2'), symmetry=1, barrier=(3.20431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20527,0.0636919,-6.9144e-05,4.57284e-08,-1.26958e-11,78516.4,32.586], Tmin=(100,'K'), Tmax=(865.982,'K')), NASAPolynomial(coeffs=[8.25811,0.0311141,-1.27137e-05,2.28552e-09,-1.54076e-13,77294.9,-0.42492], Tmin=(865.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(652.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(Cs_S) + radical(CCJCHO) + radical(Isobutyl) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH]C1C[CH]C1=O(10012)',
    structure = SMILES('[CH2][CH]C1C[CH]C1=O'),
    E0 = (398.027,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94415,0.0266983,5.65855e-05,-9.19716e-08,3.66849e-11,47962.2,27.1629], Tmin=(100,'K'), Tmax=(968.923,'K')), NASAPolynomial(coeffs=[14.2169,0.022172,-7.83596e-06,1.49995e-09,-1.13182e-13,43418.2,-42.8348], Tmin=(968.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.027,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(CCJC=O) + radical(CCJCC=O) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C[CH][CH]C=O(10020)',
    structure = SMILES('[CH2]C=C[CH]C=C[O]'),
    E0 = (238.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27925,0.0375948,4.81722e-05,-1.00984e-07,4.47846e-11,28817.6,25.6155], Tmin=(100,'K'), Tmax=(930.859,'K')), NASAPolynomial(coeffs=[20.6022,0.0114585,-1.3953e-06,1.7802e-10,-1.93605e-14,22755.2,-79.4629], Tmin=(930.859,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(238.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CCJC=C) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C=CC[C][C]=O(10004)',
    structure = SMILES('[CH2]C=CC[C][C]=O'),
    E0 = (600.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1855,455,950,201.083,893.444,893.971],'cm^-1')),
        HinderedRotor(inertia=(0.558185,'amu*angstrom^2'), symmetry=1, barrier=(12.8338,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558438,'amu*angstrom^2'), symmetry=1, barrier=(12.8396,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.558281,'amu*angstrom^2'), symmetry=1, barrier=(12.836,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.91041,'amu*angstrom^2'), symmetry=1, barrier=(89.908,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21223,0.0628526,-6.29322e-05,3.46186e-08,-7.86405e-12,72328.4,27.1781], Tmin=(100,'K'), Tmax=(1049.51,'K')), NASAPolynomial(coeffs=[10.5612,0.0272212,-1.20072e-05,2.27064e-09,-1.58656e-13,70366,-18.3773], Tmin=(1049.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(Allyl_P) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH][CH]C[C]=C[O](13236)',
    structure = SMILES('[CH2][CH][CH]C[C]=C[O]'),
    E0 = (700.58,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,180,1218.13,1218.16,1218.17],'cm^-1')),
        HinderedRotor(inertia=(0.114122,'amu*angstrom^2'), symmetry=1, barrier=(2.62389,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114122,'amu*angstrom^2'), symmetry=1, barrier=(2.62388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114088,'amu*angstrom^2'), symmetry=1, barrier=(2.62311,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11409,'amu*angstrom^2'), symmetry=1, barrier=(2.62316,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27977,0.063543,-7.26815e-05,5.38463e-08,-1.71568e-11,84354.8,33.1652], Tmin=(100,'K'), Tmax=(774.46,'K')), NASAPolynomial(coeffs=[6.67714,0.0347554,-1.51608e-05,2.8131e-09,-1.92864e-13,83546.1,8.68182], Tmin=(774.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]CC[C]=[C][O](10030)',
    structure = SMILES('[CH2][CH]CC[C][C]=O'),
    E0 = (733.204,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918776,0.0760353,-0.000114671,1.03004e-07,-3.64756e-11,88287.1,32.4301], Tmin=(100,'K'), Tmax=(850.283,'K')), NASAPolynomial(coeffs=[5.379,0.0380123,-1.75319e-05,3.27185e-09,-2.22158e-13,88144.6,15.2578], Tmin=(850.283,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH][CH]C[C]=[C]O(13237)',
    structure = SMILES('[CH2][CH][CH]C[C]=[C]O'),
    E0 = (798.862,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,1670,1700,300,440,203.999,329.52,3277.71],'cm^-1')),
        HinderedRotor(inertia=(0.0426138,'amu*angstrom^2'), symmetry=1, barrier=(1.21072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0426138,'amu*angstrom^2'), symmetry=1, barrier=(1.21072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0426138,'amu*angstrom^2'), symmetry=1, barrier=(1.21072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0426138,'amu*angstrom^2'), symmetry=1, barrier=(1.21072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0426138,'amu*angstrom^2'), symmetry=1, barrier=(1.21072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.972021,0.07359,-0.000107873,9.31239e-08,-3.16112e-11,96183.3,35.9834], Tmin=(100,'K'), Tmax=(875.153,'K')), NASAPolynomial(coeffs=[6.27319,0.0343495,-1.48871e-05,2.69087e-09,-1.78993e-13,95830.2,14.3994], Tmin=(875.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(798.862,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C[CH]C[C]=[C][O](10031)',
    structure = SMILES('[CH2]C[CH]C[C][C]=O'),
    E0 = (733.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.958964,0.0764755,-0.000119734,1.11216e-07,-4.02421e-11,88285.8,32.114], Tmin=(100,'K'), Tmax=(852.672,'K')), NASAPolynomial(coeffs=[4.28949,0.0399583,-1.87386e-05,3.5153e-09,-2.39036e-13,88477.4,21.0308], Tmin=(852.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][CH][CH][CH]C=[C]O(13238)',
    structure = SMILES('[CH2][CH]C=C[CH][C]O'),
    E0 = (689.244,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520323,0.0718154,-6.79091e-05,3.24967e-08,-6.18949e-12,83026.1,27.2957], Tmin=(100,'K'), Tmax=(1265.15,'K')), NASAPolynomial(coeffs=[16.055,0.0227006,-9.67812e-06,1.81253e-09,-1.26246e-13,79095.3,-51.3038], Tmin=(1265.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(689.244,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(C=CCJCO) + radical(RCCJ) + radical(CH2_triplet)"""),
)

species(
    label = 'C[CH][CH]C[C][C]=O(10315)',
    structure = SMILES('C[CH][CH]C[C][C]=O'),
    E0 = (722.416,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,247.496,859.297,1246.71,1853],'cm^-1')),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111507,'amu*angstrom^2'), symmetry=1, barrier=(3.35981,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14077,0.0757163,-0.000127953,1.26598e-07,-4.71935e-11,86977.3,32.4101], Tmin=(100,'K'), Tmax=(870.25,'K')), NASAPolynomial(coeffs=[1.2671,0.0443492,-2.08223e-05,3.87836e-09,-2.61356e-13,88121.1,38.5161], Tmin=(870.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.416,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    E0 = (646.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1117.91,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1135.79,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1112.53,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1139.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1112.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1101.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (649.906,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (669.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (669.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (697.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (699.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (811.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (1055.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (807.336,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (837.549,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (802.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (800.143,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (778.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (646.977,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (809.321,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (655.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (710.377,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (827.45,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (843.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (904.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (875.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (949.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (946.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (773.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (755.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['HCCO(2227)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][CH][C]=O(9632)', '[CH][CH][CH2](5537)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH2][CH][C]C[CH][C]=O(13229)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2][CH][CH]C[C][C]=O(13230)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2][C][CH]C[CH][C]=O(13231)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH][CH][CH]C[CH][C]=O(13232)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['[CH2][CH]C1CC1[C]=O(10011)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['C=CC[CH][CH][C]=O(9967)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['[CH2]C=C[CH]C[C]=O(5328)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2]C=C[CH][CH][C]=O(10015)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C][O](6861)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(106.851,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction13',
    reactants = ['CH2CHCO(3668)', '[CH][CH][CH2](5537)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[CH2][CH][CH][CH][CH][C]=O(13233)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2][CH]C[CH][CH][C]=O(10008)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][CH][CH][CH]C[C]=O(9127)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C[CH][CH][CH][C]=O(13234)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.21412e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;Y_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['[CH2][CH][CH][CH][CH]C=O(13235)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.21412e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['C[CH][CH][CH][CH][C]=O(10317)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['[CH2]C=CC[CH][C]=O(5329)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2][CH]C([CH2])[CH][C]=O(9158)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['[CH2][CH]C1C[CH]C1=O(10012)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    products = ['[CH2]C=C[CH][CH]C=O(10020)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C=CC[C][C]=O(10004)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HCCO(2227)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.04713,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH][CH]C[C]=C[O](13236)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][CH]CC[C]=[C][O](10030)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH][CH]C[C]=[C]O(13237)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2]C[CH]C[C]=[C][O](10031)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.65304e+08,'s^-1'), n=1.30038, Ea=(212.995,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_2;Cd_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH][CH][CH]C=[C]O(13238)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;O_H_out] for rate rule [R4HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[CH][CH]C[C][C]=O(10315)'],
    products = ['[CH2][CH][CH]C[CH][C]=O(9128)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2538',
    isomers = [
        '[CH2][CH][CH]C[CH][C]=O(9128)',
    ],
    reactants = [
        ('HCCO(2227)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2538',
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


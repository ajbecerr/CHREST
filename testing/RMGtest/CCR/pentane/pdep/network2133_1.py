species(
    label = '[CH2][CH][CH]CC(C)[O](4467)',
    structure = SMILES('[CH2][CH][CH]CC(C)[O]'),
    E0 = (460.177,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,429.164,522.343,2086.8,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0254566,'amu*angstrom^2'), symmetry=1, barrier=(2.93611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0254566,'amu*angstrom^2'), symmetry=1, barrier=(2.93611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0254566,'amu*angstrom^2'), symmetry=1, barrier=(2.93611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0254566,'amu*angstrom^2'), symmetry=1, barrier=(2.93611,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0254566,'amu*angstrom^2'), symmetry=1, barrier=(2.93611,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04314,0.0718628,-8.66497e-05,7.5936e-08,-2.82114e-11,55446.3,33.9585], Tmin=(100,'K'), Tmax=(812.912,'K')), NASAPolynomial(coeffs=[2.99347,0.0496657,-2.24407e-05,4.21089e-09,-2.89349e-13,55545.6,27.514], Tmin=(812.912,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.177,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CH3CHO(1381)',
    structure = SMILES('CC=O'),
    E0 = (-177.833,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,192.062,1313.7,1313.74,1629.81,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0786982,'amu*angstrom^2'), symmetry=1, barrier=(2.09733,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3625.12,'J/mol'), sigma=(3.97,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.57993,0.00518977,2.269e-05,-2.73745e-08,9.28491e-12,-21369.7,8.9697], Tmin=(100,'K'), Tmax=(1028.8,'K')), NASAPolynomial(coeffs=[4.08562,0.0139062,-5.59372e-06,1.0461e-09,-7.38743e-14,-22039.1,3.76815], Tmin=(1028.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-177.833,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), label="""CH3CHO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2][CH][CH]CC[O](1855)',
    structure = SMILES('[CH2][CH][CH]CC[O]'),
    E0 = (491.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,1757.19,1757.25,1757.25],'cm^-1')),
        HinderedRotor(inertia=(0.00194198,'amu*angstrom^2'), symmetry=1, barrier=(4.25539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184803,'amu*angstrom^2'), symmetry=1, barrier=(4.24899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.184859,'amu*angstrom^2'), symmetry=1, barrier=(4.25028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00194215,'amu*angstrom^2'), symmetry=1, barrier=(4.25555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.1164,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.88361,0.0606102,-0.000102302,1.12976e-07,-4.5796e-11,59229.9,30.0924], Tmin=(100,'K'), Tmax=(860.273,'K')), NASAPolynomial(coeffs=[-3.8809,0.0517578,-2.46962e-05,4.65748e-09,-3.17032e-13,61541.1,64.7036], Tmin=(860.273,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)[O](4413)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)[O]'),
    E0 = (468.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,280.061,952.884,2101.54],'cm^-1')),
        HinderedRotor(inertia=(0.0752685,'amu*angstrom^2'), symmetry=1, barrier=(3.39341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0752685,'amu*angstrom^2'), symmetry=1, barrier=(3.39341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0752685,'amu*angstrom^2'), symmetry=1, barrier=(3.39341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0752685,'amu*angstrom^2'), symmetry=1, barrier=(3.39341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0752685,'amu*angstrom^2'), symmetry=1, barrier=(3.39341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884591,0.0682369,-5.35431e-05,2.37339e-08,-4.49341e-12,56465.8,33.2858], Tmin=(100,'K'), Tmax=(1216.44,'K')), NASAPolynomial(coeffs=[10.304,0.0372631,-1.53491e-05,2.80179e-09,-1.91489e-13,54174.1,-14.0031], Tmin=(1216.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2][CH][CH]C[CH]C(556)',
    structure = SMILES('[CH2][CH][CH]C[CH]C'),
    E0 = (596.733,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,3630.26,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0168484,'amu*angstrom^2'), symmetry=1, barrier=(10.7202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79918,0.0586097,-7.80157e-05,8.14666e-08,-3.30989e-11,71839.7,31.8124], Tmin=(100,'K'), Tmax=(852.594,'K')), NASAPolynomial(coeffs=[-2.48283,0.0533526,-2.41736e-05,4.49758e-09,-3.0559e-13,73491.1,57.1903], Tmin=(852.594,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH2]C(C)[O](667)',
    structure = SMILES('[CH2]C(C)[O]'),
    E0 = (151.12,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,348.049],'cm^-1')),
        HinderedRotor(inertia=(0.00139827,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0820826,'amu*angstrom^2'), symmetry=1, barrier=(6.94834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0791,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.37036,0.0302848,-1.00565e-05,-6.85406e-09,4.36724e-12,18239,16.3236], Tmin=(100,'K'), Tmax=(1040.13,'K')), NASAPolynomial(coeffs=[9.01571,0.0161203,-6.05708e-06,1.11169e-09,-7.80913e-14,16240.4,-18.9599], Tmin=(1040.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.12,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CC(C)OJ) + radical(CJCO)"""),
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
    label = '[CH2][CH][C]CC(C)[O](9348)',
    structure = SMILES('[CH2][CH][C]CC(C)[O]'),
    E0 = (713.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.919792,0.072062,-7.27334e-05,4.52385e-08,-1.22109e-11,85973.5,31.7575], Tmin=(100,'K'), Tmax=(873.934,'K')), NASAPolynomial(coeffs=[7.90448,0.0400918,-1.78584e-05,3.3764e-09,-2.35241e-13,84752.8,-0.998143], Tmin=(873.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]CC(C)[O](9349)',
    structure = SMILES('[CH2][C][CH]CC(C)[O]'),
    E0 = (713.946,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.891624,0.0733482,-8.09737e-05,5.7982e-08,-1.81492e-11,85975.2,31.6841], Tmin=(100,'K'), Tmax=(759.486,'K')), NASAPolynomial(coeffs=[6.84447,0.0419935,-1.90425e-05,3.61506e-09,-2.51758e-13,85071,4.60326], Tmin=(759.486,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(713.946,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][CH][CH]CC(C)[O](9350)',
    structure = SMILES('[CH][CH][CH]CC(C)[O]'),
    E0 = (703.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,209.424,1011.85,1089.22,1223.37,1380.05,1644.86,1864.43],'cm^-1')),
        HinderedRotor(inertia=(0.124065,'amu*angstrom^2'), symmetry=1, barrier=(3.34398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124065,'amu*angstrom^2'), symmetry=1, barrier=(3.34398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124065,'amu*angstrom^2'), symmetry=1, barrier=(3.34398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124065,'amu*angstrom^2'), symmetry=1, barrier=(3.34398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124065,'amu*angstrom^2'), symmetry=1, barrier=(3.34398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.936375,0.0743842,-9.63772e-05,8.41148e-08,-3.03942e-11,84672.4,33.5592], Tmin=(100,'K'), Tmax=(823.056,'K')), NASAPolynomial(coeffs=[4.1943,0.0457109,-2.07201e-05,3.87908e-09,-2.65664e-13,84571,21.1181], Tmin=(823.056,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(703.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][CH]C1CC(C)O1(5971)',
    structure = SMILES('[CH2][CH]C1CC(C)O1'),
    E0 = (212.202,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39484,0.0410925,4.03608e-05,-8.63422e-08,3.87803e-11,25631,26.4584], Tmin=(100,'K'), Tmax=(901.679,'K')), NASAPolynomial(coeffs=[14.5126,0.025915,-5.94906e-06,8.05216e-10,-5.14526e-14,21516.8,-45.1662], Tmin=(901.679,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(212.202,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCJCO) + radical(RCCJ)"""),
)

species(
    label = 'C=CCC[C](C)[O](5982)',
    structure = SMILES('C=CCC[C](C)[O]'),
    E0 = (170.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,268.021,268.021,268.021,4000],'cm^-1')),
        HinderedRotor(inertia=(0.2271,'amu*angstrom^2'), symmetry=1, barrier=(11.5767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402417,'amu*angstrom^2'), symmetry=1, barrier=(20.5135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.402416,'amu*angstrom^2'), symmetry=1, barrier=(20.5135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00234673,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.82463,0.0699936,-5.50198e-05,2.36543e-08,-4.30661e-12,20612.7,27.5539], Tmin=(100,'K'), Tmax=(1261.77,'K')), NASAPolynomial(coeffs=[11.2242,0.0370251,-1.58263e-05,2.94591e-09,-2.0353e-13,17988.4,-25.0361], Tmin=(1261.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(170.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)O(5974)',
    structure = SMILES('[CH2]C=C[CH]C(C)O'),
    E0 = (19.6179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.247625,0.0715602,-4.36774e-05,4.99067e-09,3.27417e-12,2503.99,27.4371], Tmin=(100,'K'), Tmax=(1058.22,'K')), NASAPolynomial(coeffs=[16.5138,0.0295025,-1.15998e-05,2.13096e-09,-1.48834e-13,-2026.42,-57.098], Tmin=(1058.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(19.6179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH]CC(C)=O(4331)',
    structure = SMILES('[CH2][CH][CH]CC(C)=O'),
    E0 = (293.84,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00648823,'amu*angstrom^2'), symmetry=1, barrier=(1.81616,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3139,0.0642699,-7.04167e-05,5.80813e-08,-2.14531e-11,35432.6,32.134], Tmin=(100,'K'), Tmax=(769.892,'K')), NASAPolynomial(coeffs=[3.70612,0.0452824,-2.06448e-05,3.91782e-09,-2.72106e-13,35258.6,22.4809], Tmin=(769.892,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(CCJCC=O) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C=C[CH]C(C)[O](5621)',
    structure = SMILES('[CH2]C=C[CH]C(C)[O]'),
    E0 = (249.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,421.193,421.653,422.395],'cm^-1')),
        HinderedRotor(inertia=(0.203033,'amu*angstrom^2'), symmetry=1, barrier=(25.5492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.704557,'amu*angstrom^2'), symmetry=1, barrier=(89.151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.202815,'amu*angstrom^2'), symmetry=1, barrier=(25.56,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.201783,'amu*angstrom^2'), symmetry=1, barrier=(25.5538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.420245,0.0674528,-3.73873e-05,-5.16328e-10,5.06711e-12,30204.1,26.6255], Tmin=(100,'K'), Tmax=(1048.99,'K')), NASAPolynomial(coeffs=[16.3814,0.0278217,-1.10777e-05,2.05836e-09,-1.45074e-13,25687.3,-56.7091], Tmin=(1048.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(C=CCJCO) + radical(Allyl_P)"""),
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
    label = '[CH2]C=CC[CH][O](1848)',
    structure = SMILES('[CH2]C=CC[CH][O]'),
    E0 = (345.127,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,368.754,368.859,1940.62],'cm^-1')),
        HinderedRotor(inertia=(0.073507,'amu*angstrom^2'), symmetry=1, barrier=(7.09342,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0734922,'amu*angstrom^2'), symmetry=1, barrier=(7.09448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.524237,'amu*angstrom^2'), symmetry=1, barrier=(50.5887,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68117,0.0552208,-5.54951e-05,3.68237e-08,-1.10436e-11,41588.9,23.6157], Tmin=(100,'K'), Tmax=(778.818,'K')), NASAPolynomial(coeffs=[5.66422,0.0347641,-1.60961e-05,3.09862e-09,-2.18027e-13,40968.5,5.39543], Tmin=(778.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(345.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P)"""),
)

species(
    label = 'C[CH][O](1382)',
    structure = SMILES('C[CH][O]'),
    E0 = (156.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2362.36],'cm^-1')),
        HinderedRotor(inertia=(0.127216,'amu*angstrom^2'), symmetry=1, barrier=(2.92495,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (44.0526,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.77075,0.0113794,7.94194e-07,-3.17974e-09,7.41421e-13,18809.5,9.6394], Tmin=(100,'K'), Tmax=(1908.67,'K')), NASAPolynomial(coeffs=[6.95159,0.0119688,-5.37106e-06,9.65305e-10,-6.23696e-14,16273.7,-11.2245], Tmin=(1908.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2][CH][CH]C[CH][O](1850)',
    structure = SMILES('[CH2][CH][CH]C[CH][O]'),
    E0 = (672.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180,1701.38,1701.66,1701.95],'cm^-1')),
        HinderedRotor(inertia=(0.195647,'amu*angstrom^2'), symmetry=1, barrier=(4.49832,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00218892,'amu*angstrom^2'), symmetry=1, barrier=(4.49774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00218902,'amu*angstrom^2'), symmetry=1, barrier=(4.49819,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195437,'amu*angstrom^2'), symmetry=1, barrier=(4.49349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63772,0.0690567,-0.000132374,1.43569e-07,-5.61224e-11,80920.8,30.5248], Tmin=(100,'K'), Tmax=(875.114,'K')), NASAPolynomial(coeffs=[-3.33206,0.0486816,-2.35884e-05,4.42807e-09,-2.98564e-13,83440.6,63.2656], Tmin=(875.114,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(672.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJCC) + radical(RCCJC) + radical(CCsJOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C](C)[O](9321)',
    structure = SMILES('[CH2][CH][CH]C[C](C)[O]'),
    E0 = (636.805,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,769.75,886.25,2154.04,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0330264,'amu*angstrom^2'), symmetry=1, barrier=(3.71431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0330264,'amu*angstrom^2'), symmetry=1, barrier=(3.71431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0330264,'amu*angstrom^2'), symmetry=1, barrier=(3.71431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0330264,'amu*angstrom^2'), symmetry=1, barrier=(3.71431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0330264,'amu*angstrom^2'), symmetry=1, barrier=(3.71431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.840219,0.0811923,-0.000125845,1.20945e-07,-4.49509e-11,76692.4,34.1096], Tmin=(100,'K'), Tmax=(858.481,'K')), NASAPolynomial(coeffs=[1.89118,0.0492622,-2.28204e-05,4.25847e-09,-2.88647e-13,77508.1,35.0015], Tmin=(858.481,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(636.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C(C)[O](9323)',
    structure = SMILES('[CH2][CH][CH][CH]C(C)[O]'),
    E0 = (660.079,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,567.092,739.234,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0109053,'amu*angstrom^2'), symmetry=1, barrier=(28.7086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0109053,'amu*angstrom^2'), symmetry=1, barrier=(28.7086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0109053,'amu*angstrom^2'), symmetry=1, barrier=(28.7086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0109053,'amu*angstrom^2'), symmetry=1, barrier=(28.7086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0109053,'amu*angstrom^2'), symmetry=1, barrier=(28.7086,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22907,0.0667403,-7.82808e-05,6.66565e-08,-2.43837e-11,79483.3,35.7882], Tmin=(100,'K'), Tmax=(806.029,'K')), NASAPolynomial(coeffs=[3.593,0.0453103,-2.0351e-05,3.81436e-09,-2.62215e-13,79417.3,26.8478], Tmin=(806.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(660.079,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])[O](4458)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])[O]'),
    E0 = (671.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3000,3050,390,425,1340,1360,335,370,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,686.59,2061.2,2278.4,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0620368,'amu*angstrom^2'), symmetry=1, barrier=(6.42054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.89,'J/mol'), sigma=(6.79962,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.22 K, Pc=27.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.852795,0.0775059,-0.000107974,9.70252e-08,-3.50878e-11,80900.1,35.2748], Tmin=(100,'K'), Tmax=(840.862,'K')), NASAPolynomial(coeffs=[4.01348,0.0455248,-2.06945e-05,3.86014e-09,-2.62948e-13,80967.6,24.1364], Tmin=(840.862,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]C[C](C)O(9351)',
    structure = SMILES('[CH2][CH][CH]C[C](C)O'),
    E0 = (406.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.679388,0.0851568,-0.000131612,1.25725e-07,-4.64043e-11,48991.8,34.879], Tmin=(100,'K'), Tmax=(862.856,'K')), NASAPolynomial(coeffs=[1.97889,0.0510162,-2.33836e-05,4.34053e-09,-2.9318e-13,49814.2,34.8665], Tmin=(862.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH]C(C)[O](5967)',
    structure = SMILES('[CH2][CH]C[CH]C(C)[O]'),
    E0 = (465.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,229.252,1134.11,1472.18,2064.58],'cm^-1')),
        HinderedRotor(inertia=(0.0791247,'amu*angstrom^2'), symmetry=1, barrier=(2.91545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0791247,'amu*angstrom^2'), symmetry=1, barrier=(2.91545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0791247,'amu*angstrom^2'), symmetry=1, barrier=(2.91545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0791247,'amu*angstrom^2'), symmetry=1, barrier=(2.91545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0791247,'amu*angstrom^2'), symmetry=1, barrier=(2.91545,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18611,0.0637959,-4.55255e-05,1.82625e-08,-3.21916e-12,56100.7,33.3754], Tmin=(100,'K'), Tmax=(1258.09,'K')), NASAPolynomial(coeffs=[8.70635,0.0398859,-1.7018e-05,3.15623e-09,-2.17341e-13,54208.5,-4.63201], Tmin=(1258.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH][CH]C(C)O(9352)',
    structure = SMILES('[CH2][CH][CH][CH]C(C)O'),
    E0 = (429.718,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06203,0.0707807,-8.43192e-05,7.17811e-08,-2.59692e-11,51783,36.5797], Tmin=(100,'K'), Tmax=(817.557,'K')), NASAPolynomial(coeffs=[3.7133,0.0470064,-2.08796e-05,3.88811e-09,-2.66046e-13,51710.5,26.5309], Tmin=(817.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJCO) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC([CH2])O(4466)',
    structure = SMILES('[CH2][CH][CH]CC([CH2])O'),
    E0 = (441.405,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.689946,0.0814947,-0.000113826,1.01912e-07,-3.65812e-11,53199.6,36.0514], Tmin=(100,'K'), Tmax=(847.751,'K')), NASAPolynomial(coeffs=[4.1122,0.0472593,-2.12461e-05,3.93945e-09,-2.67249e-13,53269.3,23.9398], Tmin=(847.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(441.405,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(382.466,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC[C](C)[O](5968)',
    structure = SMILES('[CH2][CH]CC[C](C)[O]'),
    E0 = (442.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,201.384,815.196,1215.1,1680.55],'cm^-1')),
        HinderedRotor(inertia=(0.148401,'amu*angstrom^2'), symmetry=1, barrier=(3.49846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148401,'amu*angstrom^2'), symmetry=1, barrier=(3.49846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148401,'amu*angstrom^2'), symmetry=1, barrier=(3.49846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148401,'amu*angstrom^2'), symmetry=1, barrier=(3.49846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148401,'amu*angstrom^2'), symmetry=1, barrier=(3.49846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65441,0.0632846,-6.88011e-06,-1.24554e-07,1.31062e-10,53274.1,28.878], Tmin=(100,'K'), Tmax=(450.841,'K')), NASAPolynomial(coeffs=[6.16864,0.0452799,-2.03292e-05,3.80107e-09,-2.6055e-13,52643,8.21079], Tmin=(450.841,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH][CH]C(C)[O](9353)',
    structure = SMILES('[CH2]C[CH][CH]C(C)[O]'),
    E0 = (465.633,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,465.883,2087.07,3447.57],'cm^-1')),
        HinderedRotor(inertia=(0.0231111,'amu*angstrom^2'), symmetry=1, barrier=(2.79028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231111,'amu*angstrom^2'), symmetry=1, barrier=(2.79028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231111,'amu*angstrom^2'), symmetry=1, barrier=(2.79028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231111,'amu*angstrom^2'), symmetry=1, barrier=(2.79028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231111,'amu*angstrom^2'), symmetry=1, barrier=(2.79028,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32384,0.0629688,-4.5686e-05,1.95719e-08,-3.84181e-12,56095.5,32.7183], Tmin=(100,'K'), Tmax=(1112.7,'K')), NASAPolynomial(coeffs=[6.82981,0.043175,-1.90018e-05,3.58372e-09,-2.49515e-13,54870.2,5.56733], Tmin=(1112.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.633,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CCC([CH2])[O](4465)',
    structure = SMILES('[CH2][CH]CCC([CH2])[O]'),
    E0 = (477.308,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.828132,0.0740968,-7.26127e-05,4.42468e-08,-1.17812e-11,57517.3,32.8149], Tmin=(100,'K'), Tmax=(882.472,'K')), NASAPolynomial(coeffs=[7.82444,0.0423848,-1.871e-05,3.52619e-09,-2.45326e-13,56282.5,-0.0635175], Tmin=(882.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]C[C](C)[O](5969)',
    structure = SMILES('[CH2]C[CH]C[C](C)[O]'),
    E0 = (442.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,360,370,350,3000,3100,440,815,1455,1000,250.077,862.445,1261.43,1866.43],'cm^-1')),
        HinderedRotor(inertia=(0.105579,'amu*angstrom^2'), symmetry=1, barrier=(3.30322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105579,'amu*angstrom^2'), symmetry=1, barrier=(3.30322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105579,'amu*angstrom^2'), symmetry=1, barrier=(3.30322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105579,'amu*angstrom^2'), symmetry=1, barrier=(3.30322,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105579,'amu*angstrom^2'), symmetry=1, barrier=(3.30322,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.657982,0.0807676,-0.000105548,9.11905e-08,-3.26929e-11,53316.7,32.0299], Tmin=(100,'K'), Tmax=(813.967,'K')), NASAPolynomial(coeffs=[4.95195,0.0474619,-2.16808e-05,4.08032e-09,-2.80498e-13,53022,14.6813], Tmin=(813.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(442.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(C2CsJOH) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH]C(C)[O](4464)',
    structure = SMILES('C[CH][CH][CH]C(C)[O]'),
    E0 = (454.832,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,186.488,484.231,2188.02,4000],'cm^-1')),
        HinderedRotor(inertia=(0.019148,'amu*angstrom^2'), symmetry=1, barrier=(3.05425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.019148,'amu*angstrom^2'), symmetry=1, barrier=(3.05425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.019148,'amu*angstrom^2'), symmetry=1, barrier=(3.05425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.019148,'amu*angstrom^2'), symmetry=1, barrier=(3.05425,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.019148,'amu*angstrom^2'), symmetry=1, barrier=(3.05425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24597,0.0653245,-6.52512e-05,5.0834e-08,-1.83553e-11,54798.4,33.9441], Tmin=(100,'K'), Tmax=(767.74,'K')), NASAPolynomial(coeffs=[3.59289,0.047973,-2.13393e-05,4.01027e-09,-2.7733e-13,54589,24.2254], Tmin=(767.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2]C[CH]CC([CH2])[O](1817)',
    structure = SMILES('[CH2]C[CH]CC([CH2])[O]'),
    E0 = (477.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,212.706,892.407,1338.61,1994.27],'cm^-1')),
        HinderedRotor(inertia=(0.109982,'amu*angstrom^2'), symmetry=1, barrier=(3.14663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109982,'amu*angstrom^2'), symmetry=1, barrier=(3.14663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109982,'amu*angstrom^2'), symmetry=1, barrier=(3.14663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109982,'amu*angstrom^2'), symmetry=1, barrier=(3.14663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109982,'amu*angstrom^2'), symmetry=1, barrier=(3.14663,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.803974,0.0753228,-8.05758e-05,5.65245e-08,-1.7467e-11,57518.8,32.728], Tmin=(100,'K'), Tmax=(766.598,'K')), NASAPolynomial(coeffs=[6.74472,0.0443234,-1.99168e-05,3.77047e-09,-2.62326e-13,56608,5.64656], Tmin=(766.598,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(CJCO) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C[C](C)[O](4461)',
    structure = SMILES('C[CH][CH]C[C](C)[O]'),
    E0 = (431.558,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,3000,3050,390,425,1340,1360,335,370,233.359,505.304,1466.41,3371.95],'cm^-1')),
        HinderedRotor(inertia=(0.00513925,'amu*angstrom^2'), symmetry=1, barrier=(0.183471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00513925,'amu*angstrom^2'), symmetry=1, barrier=(0.183471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00513925,'amu*angstrom^2'), symmetry=1, barrier=(0.183471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00513925,'amu*angstrom^2'), symmetry=1, barrier=(0.183471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00513925,'amu*angstrom^2'), symmetry=1, barrier=(0.183471,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.815924,0.0803047,-0.000114861,1.08041e-07,-4.02724e-11,52009.2,32.4105], Tmin=(100,'K'), Tmax=(847.506,'K')), NASAPolynomial(coeffs=[2.02999,0.051675,-2.36588e-05,4.41787e-09,-3.00668e-13,52625.8,31.6061], Tmin=(847.506,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.558,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([O])C[CH][CH]C(1815)',
    structure = SMILES('[CH2]C([O])C[CH][CH]C'),
    E0 = (466.52,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,359.638,441.319,1889.65,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0447988,'amu*angstrom^2'), symmetry=1, barrier=(4.33679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.143,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3797.89,'J/mol'), sigma=(6.79962,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.22 K, Pc=27.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.83996,0.0764746,-9.64501e-05,8.3381e-08,-3.00819e-11,56216.4,33.5352], Tmin=(100,'K'), Tmax=(821.465,'K')), NASAPolynomial(coeffs=[4.10855,0.0480153,-2.15791e-05,4.03073e-09,-2.75913e-13,56102.6,20.9849], Tmin=(821.465,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.52,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(RCCJCC) + radical(RCCJC) + radical(CJCO)"""),
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
    E0 = (460.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (910.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (625.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1116.09,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (934.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (925.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (925.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (914.938,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (468.461,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (523.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (523.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (577.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (469.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (533.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (513.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (467.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (811.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (808.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (848.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (871.871,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (883.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (625.847,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (602.624,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (535.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (598.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (576.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (639.687,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (553.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (612.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (535.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (550.939,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (533.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (519.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['CH3CHO(1381)', '[CH2][CH]C=C(3743)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2][CH][CH]CC[O](1855)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2][CH]C([CH2])C(C)[O](4413)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH2][CH][CH]C[CH]C(556)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/NonDeC;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH][CH][CH2](5537)', '[CH2]C(C)[O](667)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[CH2][CH][C]CC(C)[O](9348)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH2][C][CH]CC(C)[O](9349)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH][CH][CH]CC(C)[O](9350)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['[CH2][CH]C1CC(C)O1(5971)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['C=CCC[C](C)[O](5982)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['[CH2]C=C[CH]C(C)O(5974)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH2][CH][CH]CC(C)=O(4331)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.0366254,'m^3/(mol*s)'), n=1.743, Ea=(71.4418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsCs_O;YJ] for rate rule [CO-CsCs_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH2]C=C[CH]C(C)[O](5621)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CH3CHO(1381)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.0403742,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['CH3(17)', '[CH2]C=CC[CH][O](1848)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.61258,'m^3/(mol*s)'), n=1.485, Ea=(32.0285,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO_O;CsJ-HHH] for rate rule [CO-CsH_O;CsJ-HHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH][O](1382)', '[CH2][CH]C=C(3743)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH][O](1382)', '[CH2][CH][CH][CH2](5531)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.206e+07,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['CH3(17)', '[CH2][CH][CH]C[CH][O](1850)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH2][CH][CH]C[C](C)[O](9321)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2][CH][CH][CH]C(C)[O](9323)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2][CH][CH]CC([CH2])[O](4458)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['[CH2][CH][CH]C[C](C)O(9351)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(4.56178e+08,'s^-1'), n=1.25272, Ea=(165.67,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_Cs2] for rate rule [R2H_S;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH]C[CH]C(C)[O](5967)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['[CH2][CH][CH][CH]C(C)O(9352)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['[CH2][CH][CH]CC([CH2])O(4466)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(62433.6,'s^-1'), n=2.54422, Ea=(138.678,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;Y_rad_out;Cs_H_out_2H] + [R3H_SS_Cs;O_rad_out;Cs_H_out] for rate rule [R3H_SS_Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['[CH2][CH]CC[C](C)[O](5968)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C[CH][CH]C(C)[O](9353)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(867345,'s^-1'), n=1.96939, Ea=(174.054,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R3HJ;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][CH]CCC([CH2])[O](4465)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(927.918,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['[CH2]C[CH]C[C](C)[O](5969)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(322699,'s^-1'), n=1.96429, Ea=(152.474,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;XH_out] for rate rule [R4HJ_1;Y_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['C[CH][CH][CH]C(C)[O](4464)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(754000,'s^-1'), n=1.63, Ea=(74.8936,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C[CH]CC([CH2])[O](1817)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_3;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][CH][CH]CC(C)[O](4467)'],
    products = ['C[CH][CH]C[C](C)[O](4461)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([O])C[CH][CH]C(1815)'],
    products = ['[CH2][CH][CH]CC(C)[O](4467)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2133',
    isomers = [
        '[CH2][CH][CH]CC(C)[O](4467)',
    ],
    reactants = [
        ('CH3CHO(1381)', '[CH2][CH]C=C(3743)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2133',
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


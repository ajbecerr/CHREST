species(
    label = '[CH2]CC([CH2])C(C)O[O](7403)',
    structure = SMILES('[CH2]CC([CH2])C(C)O[O]'),
    E0 = (267.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.148086,0.0863418,-8.19528e-05,4.52741e-08,-1.05303e-11,32268.8,35.0936], Tmin=(100,'K'), Tmax=(1019.71,'K')), NASAPolynomial(coeffs=[11.2946,0.0426173,-1.76332e-05,3.22287e-09,-2.20593e-13,29995.6,-18.8992], Tmin=(1019.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(267.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'C2H4(22)',
    structure = SMILES('C=C'),
    E0 = (41.9072,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.98876,-0.00674712,5.04398e-05,-5.70744e-08,2.04945e-11,5047.05,3.80494], Tmin=(100,'K'), Tmax=(946.013,'K')), NASAPolynomial(coeffs=[4.59022,0.00872722,-2.66493e-06,4.8171e-10,-3.6069e-14,4127.02,-3.3246], Tmin=(946.013,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(41.9072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), label="""C2H4""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]CC([CH2])CO[O](696)',
    structure = SMILES('[CH2]CC([CH2])CO[O]'),
    E0 = (300.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,695.965],'cm^-1')),
        HinderedRotor(inertia=(0.140985,'amu*angstrom^2'), symmetry=1, barrier=(3.24153,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138214,'amu*angstrom^2'), symmetry=1, barrier=(3.17782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135751,'amu*angstrom^2'), symmetry=1, barrier=(3.12119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137936,'amu*angstrom^2'), symmetry=1, barrier=(3.17142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.140424,'amu*angstrom^2'), symmetry=1, barrier=(3.22862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.124,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3800.62,'J/mol'), sigma=(6.68442,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=593.65 K, Pc=28.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.774703,0.0739675,-7.96967e-05,5.2536e-08,-1.46262e-11,36223,30.9461], Tmin=(100,'K'), Tmax=(861.063,'K')), NASAPolynomial(coeffs=[8.70994,0.0371054,-1.54823e-05,2.81932e-09,-1.91686e-13,34856.5,-6.14985], Tmin=(861.063,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[CH]CC(C)O[O](7602)',
    structure = SMILES('[CH2]C[CH]CC(C)O[O]'),
    E0 = (258.879,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.49405,0.0685093,5.81017e-06,-1.70129e-07,1.73606e-10,31212,31.6329], Tmin=(100,'K'), Tmax=(437.128,'K')), NASAPolynomial(coeffs=[5.95979,0.0523875,-2.37715e-05,4.47369e-09,-3.08086e-13,30585.2,11.0799], Tmin=(437.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(258.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC[CH]C(C)O[O](7402)',
    structure = SMILES('[CH2]CC[CH]C(C)O[O]'),
    E0 = (264.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.27434,0.0874905,-8.87594e-05,5.485e-08,-1.46737e-11,31981.8,35.4556], Tmin=(100,'K'), Tmax=(881.446,'K')), NASAPolynomial(coeffs=[8.95655,0.0480898,-2.17079e-05,4.13564e-09,-2.89519e-13,30451.2,-5.33523], Tmin=(881.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(264.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCJCOOH) + radical(RCCJ)"""),
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
    label = '[CH2]CC([CH2])C(C)[O](1782)',
    structure = SMILES('[CH2]CC([CH2])C(C)[O]'),
    E0 = (274.007,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.338425,0.0735846,-5.62828e-05,2.34097e-08,-4.01482e-12,33093,32.693], Tmin=(100,'K'), Tmax=(1370.38,'K')), NASAPolynomial(coeffs=[13.9828,0.0337579,-1.26889e-05,2.20196e-09,-1.45863e-13,29353.4,-37.4325], Tmin=(1370.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.007,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(Isobutyl) + radical(RCCJ)"""),
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
    label = '[CH2]C[CH]C(C)O[O](7603)',
    structure = SMILES('[CH2]C[CH]C(C)O[O]'),
    E0 = (288.616,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,324.798,3547.11],'cm^-1')),
        HinderedRotor(inertia=(0.106997,'amu*angstrom^2'), symmetry=1, barrier=(8.01178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00159754,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107014,'amu*angstrom^2'), symmetry=1, barrier=(8.01199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107001,'amu*angstrom^2'), symmetry=1, barrier=(8.01165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.653403,'amu*angstrom^2'), symmetry=1, barrier=(48.9184,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.124,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884212,0.0731761,-7.67009e-05,4.92706e-08,-1.36514e-11,34820.6,31.0183], Tmin=(100,'K'), Tmax=(854.978,'K')), NASAPolynomial(coeffs=[8.05032,0.0396478,-1.78747e-05,3.39856e-09,-2.3747e-13,33595.3,-2.43088], Tmin=(854.978,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCJCOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH2])C(C)O[O](7604)',
    structure = SMILES('[CH2]C([CH2])C(C)O[O]'),
    E0 = (287.425,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1620.7],'cm^-1')),
        HinderedRotor(inertia=(0.130008,'amu*angstrom^2'), symmetry=1, barrier=(5.46355,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129999,'amu*angstrom^2'), symmetry=1, barrier=(5.46317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129764,'amu*angstrom^2'), symmetry=1, barrier=(5.46402,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.129884,'amu*angstrom^2'), symmetry=1, barrier=(5.46366,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.57217,'amu*angstrom^2'), symmetry=1, barrier=(66.0898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (101.124,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.602469,0.074606,-7.77155e-05,4.75156e-08,-1.1941e-11,34691.7,30.5146], Tmin=(100,'K'), Tmax=(962.402,'K')), NASAPolynomial(coeffs=[10.7305,0.0325114,-1.21069e-05,2.06794e-09,-1.35264e-13,32742.2,-17.9591], Tmin=(962.402,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.425,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(361.68,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    label = '[CH]C(C[CH2])C(C)O[O](7605)',
    structure = SMILES('[CH]C(C[CH2])C(C)O[O]'),
    E0 = (510.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.190058,0.0863636,-8.41568e-05,4.61096e-08,-1.05574e-11,61508.7,34.1674], Tmin=(100,'K'), Tmax=(1035.78,'K')), NASAPolynomial(coeffs=[12.0756,0.0404647,-1.76882e-05,3.32889e-09,-2.31878e-13,59046.5,-23.5917], Tmin=(1035.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]CC([CH2])C(C)O[O](7606)',
    structure = SMILES('[CH]CC([CH2])C(C)O[O]'),
    E0 = (510.125,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0206854,0.0891056,-9.2554e-05,5.47006e-08,-1.3337e-11,61495.9,34.7687], Tmin=(100,'K'), Tmax=(985.413,'K')), NASAPolynomial(coeffs=[12.4098,0.0388165,-1.60052e-05,2.91362e-09,-1.98824e-13,59054.1,-24.8198], Tmin=(985.413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.125,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC(O[O])C1CCC1(7607)',
    structure = SMILES('CC(O[O])C1CCC1'),
    E0 = (8.3103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.667392,0.0577938,6.99607e-06,-4.594e-08,2.03446e-11,1133.04,29.2228], Tmin=(100,'K'), Tmax=(1012.63,'K')), NASAPolynomial(coeffs=[14.9573,0.0367817,-1.43686e-05,2.68243e-09,-1.91007e-13,-3577.8,-48.8678], Tmin=(1012.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(8.3103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(ROOJ)"""),
)

species(
    label = '[CH2]CC1COOC1C(7608)',
    structure = SMILES('[CH2]CC1COOC1C'),
    E0 = (-12.1992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674698,0.0551667,2.76579e-05,-8.15251e-08,3.8266e-11,-1330.77,27.3975], Tmin=(100,'K'), Tmax=(910.476,'K')), NASAPolynomial(coeffs=[17.3823,0.029049,-7.21368e-06,1.04823e-09,-6.91833e-14,-6333,-62.403], Tmin=(910.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-12.1992,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(444.824,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxolane) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1CCOOC1C(7609)',
    structure = SMILES('[CH2]C1CCOOC1C'),
    E0 = (-20.6181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.799325,0.0494361,5.01453e-05,-1.10204e-07,5.03816e-11,-2344.82,25.0528], Tmin=(100,'K'), Tmax=(896.113,'K')), NASAPolynomial(coeffs=[18.6918,0.0258085,-4.44354e-06,4.42616e-10,-2.55207e-14,-7809.63,-71.9047], Tmin=(896.113,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.6181,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(448.981,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxane) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC1COC1C(7610)',
    structure = SMILES('[CH2]CC1COC1C'),
    E0 = (15.9017,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25889,0.0454478,3.26948e-05,-7.6936e-08,3.47938e-11,2025.09,24.5149], Tmin=(100,'K'), Tmax=(904.459,'K')), NASAPolynomial(coeffs=[13.4419,0.0308402,-8.21023e-06,1.222e-09,-7.93413e-14,-1785.05,-41.9179], Tmin=(904.459,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.9017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C1CCOC1C(7611)',
    structure = SMILES('[CH2]C1CCOC1C'),
    E0 = (-62.7246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1509,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59362,0.0345067,6.63965e-05,-1.14753e-07,4.95005e-11,-7440.05,21.1105], Tmin=(100,'K'), Tmax=(892.983,'K')), NASAPolynomial(coeffs=[13.8832,0.0284988,-5.88984e-06,7.13882e-10,-4.32216e-14,-11590.2,-47.7367], Tmin=(892.983,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.7246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(Isobutyl)"""),
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
    label = '[CH2]CC([CH2])=CC(3868)',
    structure = SMILES('[CH2]CC([CH2])=CC'),
    E0 = (265.143,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.14271,'amu*angstrom^2'), symmetry=1, barrier=(3.28118,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96641,'amu*angstrom^2'), symmetry=1, barrier=(22.2197,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141835,'amu*angstrom^2'), symmetry=1, barrier=(3.26106,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0122451,'amu*angstrom^2'), symmetry=1, barrier=(22.2301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10951,0.0550342,-2.49458e-05,-1.64911e-09,3.23214e-12,32000.5,25.088], Tmin=(100,'K'), Tmax=(1122.91,'K')), NASAPolynomial(coeffs=[11.6971,0.0315481,-1.25795e-05,2.29336e-09,-1.57863e-13,28725.7,-31.2129], Tmin=(1122.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.143,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]CC([CH2])C=C(3869)',
    structure = SMILES('[CH2]CC([CH2])C=C'),
    E0 = (338.275,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,3819.65,3819.8],'cm^-1')),
        HinderedRotor(inertia=(0.918967,'amu*angstrom^2'), symmetry=1, barrier=(21.1289,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919163,'amu*angstrom^2'), symmetry=1, barrier=(21.1334,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.78395,'amu*angstrom^2'), symmetry=1, barrier=(87.0006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127208,'amu*angstrom^2'), symmetry=1, barrier=(86.9922,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.1436,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3226.55,'J/mol'), sigma=(5.8838,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.98 K, Pc=35.94 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28275,0.051164,-1.4278e-05,-1.45061e-08,8.61961e-12,40790.4,27.6427], Tmin=(100,'K'), Tmax=(1000.13,'K')), NASAPolynomial(coeffs=[11.0728,0.0305554,-1.11858e-05,1.97483e-09,-1.35046e-13,37904.6,-24.2271], Tmin=(1000.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = 'C=C(CC)C(C)O[O](7612)',
    structure = SMILES('C=C(CC)C(C)O[O]'),
    E0 = (-27.5193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0154928,0.0803301,-6.10673e-05,2.41323e-08,-3.89241e-12,-3160.52,32.843], Tmin=(100,'K'), Tmax=(1451.01,'K')), NASAPolynomial(coeffs=[16.509,0.0348624,-1.40642e-05,2.53669e-09,-1.71612e-13,-7946.96,-52.8685], Tmin=(1451.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.5193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'C=CC(C)C(C)O[O](7613)',
    structure = SMILES('C=CC(C)C(C)O[O]'),
    E0 = (-21.0149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.116256,0.0792347,-5.96925e-05,2.37041e-08,-3.87274e-12,-2382.8,32.1441], Tmin=(100,'K'), Tmax=(1425.84,'K')), NASAPolynomial(coeffs=[15.3884,0.0363908,-1.46202e-05,2.63e-09,-1.77707e-13,-6737.92,-46.9531], Tmin=(1425.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-21.0149,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(436.51,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = '[CH2]CC(=C)C(C)OO(7614)',
    structure = SMILES('[CH2]CC(=C)C(C)OO'),
    E0 = (25.7224,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.222474,0.0847527,-6.69695e-05,2.69867e-08,-4.39502e-12,3252.22,34.7268], Tmin=(100,'K'), Tmax=(1447.59,'K')), NASAPolynomial(coeffs=[18.3149,0.0335302,-1.38927e-05,2.5431e-09,-1.73608e-13,-2114.71,-61.5624], Tmin=(1447.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.7224,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C=C)C(C)OO(7615)',
    structure = SMILES('[CH2]C(C=C)C(C)OO'),
    E0 = (32.0628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.20133,0.0853695,-7.05741e-05,3.10458e-08,-5.56228e-12,4013.28,34.3061], Tmin=(100,'K'), Tmax=(1324.81,'K')), NASAPolynomial(coeffs=[16.4271,0.0351638,-1.37298e-05,2.44101e-09,-1.64407e-13,-392.641,-50.5934], Tmin=(1324.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.0628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC(=C)C(C)O[O](7616)',
    structure = SMILES('[CH2]CC(=C)C(C)O[O]'),
    E0 = (177.727,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,746.583],'cm^-1')),
        HinderedRotor(inertia=(0.48927,'amu*angstrom^2'), symmetry=1, barrier=(11.2493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.133967,'amu*angstrom^2'), symmetry=1, barrier=(3.08017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.489285,'amu*angstrom^2'), symmetry=1, barrier=(11.2496,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0284723,'amu*angstrom^2'), symmetry=1, barrier=(11.2478,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0284302,'amu*angstrom^2'), symmetry=1, barrier=(11.2487,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.247709,0.0788617,-6.42846e-05,2.76733e-08,-4.89546e-12,21513.6,33.7901], Tmin=(100,'K'), Tmax=(1325.92,'K')), NASAPolynomial(coeffs=[14.8077,0.0349377,-1.45942e-05,2.68933e-09,-1.84794e-13,17652.5,-40.561], Tmin=(1325.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.727,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(ROOJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(C=C)C(C)O[O](7475)',
    structure = SMILES('[CH2]C(C=C)C(C)O[O]'),
    E0 = (184.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,215.48,2445.14],'cm^-1')),
        HinderedRotor(inertia=(0.395801,'amu*angstrom^2'), symmetry=1, barrier=(13.0419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.39581,'amu*angstrom^2'), symmetry=1, barrier=(13.0419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00363053,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.395805,'amu*angstrom^2'), symmetry=1, barrier=(13.0419,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.26012,'amu*angstrom^2'), symmetry=1, barrier=(74.4687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.215436,0.0800314,-6.94931e-05,3.34034e-08,-6.62535e-12,22277.1,33.5665], Tmin=(100,'K'), Tmax=(1196.68,'K')), NASAPolynomial(coeffs=[13.4573,0.0357697,-1.40129e-05,2.49575e-09,-1.68449e-13,19107.8,-32.6956], Tmin=(1196.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(184.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Isobutyl)"""),
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
    label = 'buten3yl1(65)',
    structure = SMILES('[CH2]CC=C'),
    E0 = (191.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,345.872],'cm^-1')),
        HinderedRotor(inertia=(0.0750497,'amu*angstrom^2'), symmetry=1, barrier=(6.37545,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0750659,'amu*angstrom^2'), symmetry=1, barrier=(6.37555,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.486,0.0277404,-7.51225e-07,-1.39969e-08,6.14177e-12,23115.6,15.6915], Tmin=(100,'K'), Tmax=(1051.01,'K')), NASAPolynomial(coeffs=[7.36567,0.0210618,-8.19307e-06,1.49013e-09,-1.03098e-13,21433,-11.2176], Tmin=(1051.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(191.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""buten3yl1""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C2H4(T)(23)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (318.146,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,995.342,995.492,2246.13,2948.67],'cm^-1')),
        HinderedRotor(inertia=(0.00723049,'amu*angstrom^2'), symmetry=1, barrier=(5.08405,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.40735,0.0100313,6.40897e-06,-1.41287e-08,5.92653e-12,38288.2,6.11706], Tmin=(100,'K'), Tmax=(954.266,'K')), NASAPolynomial(coeffs=[5.52251,0.00856169,-2.90741e-06,5.02347e-10,-3.44567e-14,37547.7,-5.75289], Tmin=(954.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.146,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C2H4(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]CC([CH2])[CH]C(230)',
    structure = SMILES('[CH2]CC([CH2])[CH]C'),
    E0 = (410.659,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1043.58],'cm^-1')),
        HinderedRotor(inertia=(0.00458383,'amu*angstrom^2'), symmetry=1, barrier=(3.54033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153753,'amu*angstrom^2'), symmetry=1, barrier=(3.53509,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153943,'amu*angstrom^2'), symmetry=1, barrier=(3.53945,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.153951,'amu*angstrom^2'), symmetry=1, barrier=(3.53964,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00457906,'amu*angstrom^2'), symmetry=1, barrier=(3.5376,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1515,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3225.18,'J/mol'), sigma=(6.08349,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=503.77 K, Pc=32.5 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22899,0.0571953,-3.33927e-05,1.01284e-08,-1.29003e-12,49493.3,29.8858], Tmin=(100,'K'), Tmax=(1719.56,'K')), NASAPolynomial(coeffs=[11.2019,0.0339964,-1.31557e-05,2.28254e-09,-1.49341e-13,46063.5,-23.6337], Tmin=(1719.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C[CH2](66)',
    structure = SMILES('[CH2][CH]C[CH2]'),
    E0 = (460.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180],'cm^-1')),
        HinderedRotor(inertia=(0.000785383,'amu*angstrom^2'), symmetry=1, barrier=(2.49202,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108433,'amu*angstrom^2'), symmetry=1, barrier=(2.4931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000786299,'amu*angstrom^2'), symmetry=1, barrier=(2.4954,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.62775,0.0341381,-3.63084e-05,3.45494e-08,-1.41477e-11,55480.4,19.9882], Tmin=(100,'K'), Tmax=(807.221,'K')), NASAPolynomial(coeffs=[1.67633,0.030324,-1.33729e-05,2.51877e-09,-1.74066e-13,55911.9,26.0956], Tmin=(807.221,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(460.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJ) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C]([CH2])C(C)O[O](7617)',
    structure = SMILES('[CH2]C[C]([CH2])C(C)O[O]'),
    E0 = (458.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,193.822,2936.59],'cm^-1')),
        HinderedRotor(inertia=(0.0586375,'amu*angstrom^2'), symmetry=1, barrier=(1.5247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0586375,'amu*angstrom^2'), symmetry=1, barrier=(1.5247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0586375,'amu*angstrom^2'), symmetry=1, barrier=(1.5247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0586375,'amu*angstrom^2'), symmetry=1, barrier=(1.5247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0586375,'amu*angstrom^2'), symmetry=1, barrier=(1.5247,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0586375,'amu*angstrom^2'), symmetry=1, barrier=(1.5247,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.747605,0.0729423,-5.90385e-05,2.79247e-08,-5.73323e-12,55208.6,37.7827], Tmin=(100,'K'), Tmax=(1120.52,'K')), NASAPolynomial(coeffs=[9.43722,0.0419222,-1.75125e-05,3.2182e-09,-2.20893e-13,53261.3,-5.12844], Tmin=(1120.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(RCCJ) + radical(Isobutyl)"""),
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
    label = '[CH2]CC([CH2])[CH]O[O](779)',
    structure = SMILES('[CH2]CC([CH2])[CH]O[O]'),
    E0 = (488.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1273.11],'cm^-1')),
        HinderedRotor(inertia=(0.0534523,'amu*angstrom^2'), symmetry=1, barrier=(61.4499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00577282,'amu*angstrom^2'), symmetry=1, barrier=(6.6352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288535,'amu*angstrom^2'), symmetry=1, barrier=(6.63398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288786,'amu*angstrom^2'), symmetry=1, barrier=(6.63975,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.288886,'amu*angstrom^2'), symmetry=1, barrier=(6.64206,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (100.116,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.662351,0.0781752,-0.00010216,8.05905e-08,-2.59442e-11,58906.5,31.492], Tmin=(100,'K'), Tmax=(846.904,'K')), NASAPolynomial(coeffs=[8.17976,0.0358996,-1.52922e-05,2.77054e-09,-1.85974e-13,57876,-2.09255], Tmin=(846.904,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(488.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CCsJOOH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])[C](C)O[O](7618)',
    structure = SMILES('[CH2]CC([CH2])[C](C)O[O]'),
    E0 = (454.055,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1600],'cm^-1')),
        HinderedRotor(inertia=(0.149927,'amu*angstrom^2'), symmetry=1, barrier=(3.44713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149927,'amu*angstrom^2'), symmetry=1, barrier=(3.44713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149927,'amu*angstrom^2'), symmetry=1, barrier=(3.44713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149927,'amu*angstrom^2'), symmetry=1, barrier=(3.44713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149927,'amu*angstrom^2'), symmetry=1, barrier=(3.44713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149927,'amu*angstrom^2'), symmetry=1, barrier=(3.44713,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.138742,0.090654,-0.00011366,8.94297e-08,-2.93202e-11,54743.8,35.4167], Tmin=(100,'K'), Tmax=(818.951,'K')), NASAPolynomial(coeffs=[8.01287,0.045842,-1.99472e-05,3.67104e-09,-2.49431e-13,53667.1,0.301894], Tmin=(818.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.055,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CsJOOH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)O[O](7619)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)O[O]'),
    E0 = (461.698,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200.489,1750.58],'cm^-1')),
        HinderedRotor(inertia=(0.124378,'amu*angstrom^2'), symmetry=1, barrier=(3.42294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124378,'amu*angstrom^2'), symmetry=1, barrier=(3.42294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124378,'amu*angstrom^2'), symmetry=1, barrier=(3.42294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124378,'amu*angstrom^2'), symmetry=1, barrier=(3.42294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124378,'amu*angstrom^2'), symmetry=1, barrier=(3.42294,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124378,'amu*angstrom^2'), symmetry=1, barrier=(3.42294,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.343614,0.0849869,-9.25349e-05,6.20636e-08,-1.77399e-11,55657.1,36.9537], Tmin=(100,'K'), Tmax=(835.856,'K')), NASAPolynomial(coeffs=[8.8817,0.0441269,-1.92073e-05,3.57735e-09,-2.46567e-13,54229.8,-2.70652], Tmin=(835.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(461.698,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC([CH2])C([CH2])O[O](7620)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])O[O]'),
    E0 = (481.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,200.213,806.253],'cm^-1')),
        HinderedRotor(inertia=(0.155193,'amu*angstrom^2'), symmetry=1, barrier=(3.57835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155193,'amu*angstrom^2'), symmetry=1, barrier=(3.57835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155193,'amu*angstrom^2'), symmetry=1, barrier=(3.57835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155193,'amu*angstrom^2'), symmetry=1, barrier=(3.57835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155193,'amu*angstrom^2'), symmetry=1, barrier=(3.57835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155193,'amu*angstrom^2'), symmetry=1, barrier=(3.57835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.142,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.029069,0.0896728,-9.71639e-05,6.04251e-08,-1.54864e-11,58006.2,37.5174], Tmin=(100,'K'), Tmax=(940.41,'K')), NASAPolynomial(coeffs=[11.9862,0.038814,-1.60424e-05,2.91787e-09,-1.98732e-13,55757.3,-19.4343], Tmin=(940.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(481.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(CJCOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C[C](C)C(C)O[O](7621)',
    structure = SMILES('[CH2]C[C](C)C(C)O[O]'),
    E0 = (252.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.728479,0.071384,-4.72575e-05,1.63593e-08,-2.39812e-12,30544.6,36.0585], Tmin=(100,'K'), Tmax=(1508.35,'K')), NASAPolynomial(coeffs=[12.0133,0.0414577,-1.74969e-05,3.20564e-09,-2.17989e-13,27140.3,-23.0227], Tmin=(1508.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C([CH]C)C(C)O[O](7622)',
    structure = SMILES('[CH2]C([CH]C)C(C)O[O]'),
    E0 = (256.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.386213,0.0833056,-7.87259e-05,4.53423e-08,-1.13179e-11,30970.9,35.0149], Tmin=(100,'K'), Tmax=(940.67,'K')), NASAPolynomial(coeffs=[9.04613,0.0464811,-2.00051e-05,3.72593e-09,-2.57608e-13,29341.7,-6.23454], Tmin=(940.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])[C](C)OO(7623)',
    structure = SMILES('[CH2]CC([CH2])[C](C)OO'),
    E0 = (302.051,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0708499,0.091725,-9.91042e-05,6.5752e-08,-1.86646e-11,36465.1,34.9152], Tmin=(100,'K'), Tmax=(839.436,'K')), NASAPolynomial(coeffs=[9.1998,0.0482241,-2.1371e-05,4.01691e-09,-2.78464e-13,34932.5,-7.52878], Tmin=(839.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.051,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CsJOOH) + radical(Isobutyl) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC(C)[C](C)O[O](7624)',
    structure = SMILES('[CH2]CC(C)[C](C)O[O]'),
    E0 = (248.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.362838,0.0858918,-8.93469e-05,6.01117e-08,-1.78656e-11,30070.2,32.8458], Tmin=(100,'K'), Tmax=(793.546,'K')), NASAPolynomial(coeffs=[7.32341,0.0508027,-2.30134e-05,4.3791e-09,-3.05867e-13,28965.6,0.875332], Tmin=(793.546,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.973,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CsJOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]C(C)C(C)O[O](7625)',
    structure = SMILES('[CH2][CH]C(C)C(C)O[O]'),
    E0 = (256.616,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.577737,0.0803109,-6.9484e-05,3.57036e-08,-8.10769e-12,30982.8,34.3336], Tmin=(100,'K'), Tmax=(1009.36,'K')), NASAPolynomial(coeffs=[8.69096,0.0481586,-2.17022e-05,4.14401e-09,-2.90862e-13,29345,-4.88351], Tmin=(1009.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.616,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Cs_S) + radical(RCCJ)"""),
)

species(
    label = '[CH2][C](CC)C(C)O[O](7626)',
    structure = SMILES('[CH2][C](CC)C(C)O[O]'),
    E0 = (252.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.621211,0.0733348,-5.2766e-05,2.12518e-08,-3.6837e-12,30529.5,36.4433], Tmin=(100,'K'), Tmax=(1305.41,'K')), NASAPolynomial(coeffs=[10.7034,0.0424415,-1.72678e-05,3.12318e-09,-2.11895e-13,27897.2,-14.8845], Tmin=(1305.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CJCOOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(CC)[C](C)O[O](7627)',
    structure = SMILES('[CH2]C(CC)[C](C)O[O]'),
    E0 = (248.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.57821,0.0830716,-7.18868e-05,2.15639e-08,8.63019e-12,30040.8,32.1228], Tmin=(100,'K'), Tmax=(597.24,'K')), NASAPolynomial(coeffs=[7.72766,0.0490556,-2.12825e-05,3.95416e-09,-2.72119e-13,28939.5,-0.75464], Tmin=(597.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(C2CsJOOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C[C]([CH2])C(C)OO(7628)',
    structure = SMILES('[CH2]C[C]([CH2])C(C)OO'),
    E0 = (306.062,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,211.876,1184.02],'cm^-1')),
        HinderedRotor(inertia=(0.124942,'amu*angstrom^2'), symmetry=1, barrier=(3.37655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124942,'amu*angstrom^2'), symmetry=1, barrier=(3.37655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124942,'amu*angstrom^2'), symmetry=1, barrier=(3.37655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124942,'amu*angstrom^2'), symmetry=1, barrier=(3.37655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124942,'amu*angstrom^2'), symmetry=1, barrier=(3.37655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124942,'amu*angstrom^2'), symmetry=1, barrier=(3.37655,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124942,'amu*angstrom^2'), symmetry=1, barrier=(3.37655,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.390068,0.0776928,-5.85026e-05,2.39543e-08,-4.14085e-12,36941.9,38.3014], Tmin=(100,'K'), Tmax=(1325.98,'K')), NASAPolynomial(coeffs=[12.5317,0.0410661,-1.70693e-05,3.12288e-09,-2.13314e-13,33722,-23.7008], Tmin=(1325.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(306.062,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(C2CJCOOH) + radical(RCCJ) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC([CH2])C([CH2])OO(7629)',
    structure = SMILES('[CH2]CC([CH2])C([CH2])OO'),
    E0 = (329.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.201258,0.0929067,-9.13457e-05,4.97797e-08,-1.11658e-11,39734.2,37.5823], Tmin=(100,'K'), Tmax=(1066.09,'K')), NASAPolynomial(coeffs=[13.9504,0.0398099,-1.66389e-05,3.06333e-09,-2.10825e-13,36716.8,-31.5972], Tmin=(1066.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.114,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Isobutyl) + radical(CJCOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]CC(C)C([CH2])O[O](7630)',
    structure = SMILES('[CH2]CC(C)C([CH2])O[O]'),
    E0 = (276.037,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.15867,0.0862854,-7.88612e-05,4.05512e-08,-8.76617e-12,33336.4,35.2681], Tmin=(100,'K'), Tmax=(1088.75,'K')), NASAPolynomial(coeffs=[12.1445,0.0422506,-1.81937e-05,3.40332e-09,-2.3629e-13,30726.4,-23.5757], Tmin=(1088.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(276.037,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(CJCOOH) + radical(RCCJ)"""),
)

species(
    label = '[CH2]C(CC)C([CH2])O[O](7631)',
    structure = SMILES('[CH2]C(CC)C([CH2])O[O]'),
    E0 = (275.873,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00454791,0.0889371,-8.68748e-05,4.85781e-08,-1.12872e-11,33323.3,35.8487], Tmin=(100,'K'), Tmax=(1027.89,'K')), NASAPolynomial(coeffs=[12.3658,0.0407978,-1.66249e-05,3.01531e-09,-2.05515e-13,30780.2,-24.1714], Tmin=(1027.89,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.873,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(432.353,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(ROOJ) + radical(Isobutyl) + radical(CJCOOH)"""),
)

species(
    label = '[CH2][CH]C([CH2])C(C)OO(7632)',
    structure = SMILES('[CH2][CH]C([CH2])C(C)OO'),
    E0 = (309.694,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.15,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.221908,0.086901,-8.1928e-05,4.49523e-08,-1.0538e-11,37380.4,36.6318], Tmin=(100,'K'), Tmax=(1001.71,'K')), NASAPolynomial(coeffs=[10.5657,0.0455958,-2.00748e-05,3.78647e-09,-2.63917e-13,35308.2,-13.2884], Tmin=(1001.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(309.694,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(Isobutyl) + radical(RCCJ)"""),
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
    E0 = (267.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (719.023,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (427.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (427.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (517.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (726.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (725.227,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (722.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (721.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (275.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (274.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (276.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (350.054,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (329.08,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (406.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (393.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (330.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (330.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (292.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (292.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (389.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (405.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (390.873,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (365.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (375.759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (402.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (623.921,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (630.543,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (669.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (625.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (665.848,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (673.491,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (692.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (369.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (414.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (428.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (404.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (385.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (385.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (376.735,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (389.952,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (407.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (358.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (328.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (351.48,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['C2H4(22)', 'C=CC(C)O[O](3748)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2]CC([CH2])CO[O](696)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(873476,'m^3/(mol*s)'), n=0.189, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [carbene;Cs_H] for rate rule [carbene;C/H2/CsO]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -1.5 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]C[CH]CC(C)O[O](7602)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]CC[CH]C(C)O[O](7402)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH2]CC([CH2])C(C)[O](1782)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2]C[CH]C(C)O[O](7603)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/NonDeC;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH2(T)(20)', '[CH2]C([CH2])C(C)O[O](7604)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C(C[CH2])C(C)O[O](7605)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]CC([CH2])C(C)O[O](7606)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['CC(O[O])C1CCC1(7607)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]CC1COOC1C(7608)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]C1CCOOC1C(7609)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.21e+10,'s^-1'), n=0.137, Ea=(8.87008,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R6_SSSSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R6_SSSSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['O(4)', '[CH2]CC1COC1C(7610)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.47e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO_SS;C_pri_rad_intra;OO] for rate rule [R3OO_SS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['O(4)', '[CH2]C1CCOC1C(7611)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.13e+10,'s^-1'), n=0, Ea=(61.9232,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4OO_SSS;C_pri_rad_intra;OO] for rate rule [R4OO_SSS;C_pri_rad_intra;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['HO2(9)', '[CH2]CC([CH2])=CC(3868)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(9.58174e+10,'s^-1'), n=0.573333, Ea=(139.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['HO2(9)', '[CH2]CC([CH2])C=C(3869)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.937e+09,'s^-1'), n=1.17, Ea=(125.938,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 5 used for R2OO_2H_HNd
Exact match found for rate rule [R2OO_2H_HNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['C=C(CC)C(C)O[O](7612)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['C=CC(C)C(C)O[O](7613)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]CC(=C)C(C)OO(7614)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]C(C=C)C(C)OO(7615)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2]CC(=C)C(C)O[O](7616)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.0051739,'m^3/(mol*s)'), n=2.82163, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 102 used for Cds-CsCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -4.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2]C(C=C)C(C)O[O](7475)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.01e+08,'cm^3/(mol*s)'), n=1.6, Ea=(10.0416,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 10 used for Cds-CsH_Cds-HH;HJ
Exact match found for rate rule [Cds-CsH_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C[CH]O[O](57)', 'buten3yl1(65)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.00140605,'m^3/(mol*s)'), n=2.41796, Ea=(36.1629,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CsH_Cds-HH;CsJ] for rate rule [Cds-CsH_Cds-HH;CsJ-OsCsH]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C2H4(T)(23)', 'C=CC(C)O[O](3748)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.08019,'m^3/(mol*s)'), n=2.05449, Ea=(13.6169,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds-HH;CJ] + [Cds-Cs\O2s/H_Cds-HH;YJ] for rate rule [Cds-Cs\O2s/H_Cds-HH;CJ]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['C2H4(22)', '[CH2][CH]C(C)O[O](3749)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00304282,'m^3/(mol*s)'), n=2.49542, Ea=(21.4551,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O2(2)', '[CH2]CC([CH2])[CH]C(230)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.57416e+06,'m^3/(mol*s)'), n=-0.129357, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_Ext-4R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_Ext-4R!H-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_Ext-4R!H-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C[CH]O[O](57)', '[CH2][CH]C[CH2](66)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C2H4(T)(23)', '[CH2][CH]C(C)O[O](3749)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', '[CH2]C[C]([CH2])C(C)O[O](7617)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['CH3(17)', '[CH2]CC([CH2])[CH]O[O](779)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(308.407,'m^3/(mol*s)'), n=0.967216, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0049236047651, var=0.0302625193734, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.361117102774
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2]CC([CH2])[C](C)O[O](7618)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2][CH]C([CH2])C(C)O[O](7619)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2]CC([CH2])C([CH2])O[O](7620)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]C[C](C)C(C)O[O](7621)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.265e-07,'s^-1'), n=5.639, Ea=(102.68,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 38 used for R2H_S;C_rad_out_2H;Cs_H_out_Cs2
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]C([CH]C)C(C)O[O](7622)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]CC([CH2])[C](C)OO(7623)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4.97e+09,'s^-1'), n=1.01, Ea=(160.958,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 245 used for R3H_SS_O;O_rad_out;Cs_H_out_Cs2
Exact match found for rate rule [R3H_SS_O;O_rad_out;Cs_H_out_Cs2]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]CC(C)[C](C)O[O](7624)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.25e+10,'s^-1'), n=0.66, Ea=(137.654,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2][CH]C(C)C(C)O[O](7625)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2][C](CC)C(C)O[O](7626)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(83345.1,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]C(CC)[C](C)O[O](7627)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.86e+10,'s^-1'), n=0.58, Ea=(109.579,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4H_SSS;C_rad_out_2H;Cs_H_out_NonDe] for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C[C]([CH2])C(C)OO(7628)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;O_H_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    products = ['[CH2]CC([CH2])C([CH2])OO(7629)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.18e+08,'s^-1'), n=1.06, Ea=(140.206,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 247 used for R4H_SSS_O(Cs)Cs;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS_O(Cs)Cs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH2]CC(C)C([CH2])O[O](7630)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(114000,'s^-1'), n=1.74, Ea=(82.8432,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 109 used for R4H_SSS;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C(CC)C([CH2])O[O](7631)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(68850,'s^-1'), n=1.68, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 111 used for R5H_CCC;C_rad_out_2H;Cs_H_out_2H
Exact match found for rate rule [R5H_CCC;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][CH]C([CH2])C(C)OO(7632)'],
    products = ['[CH2]CC([CH2])C(C)O[O](7403)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(67170.6,'s^-1'), n=1.77845, Ea=(41.7861,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSSS;Y_rad_out;XH_out] for rate rule [R5H_SSSS;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1768',
    isomers = [
        '[CH2]CC([CH2])C(C)O[O](7403)',
    ],
    reactants = [
        ('C2H4(22)', 'C=CC(C)O[O](3748)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1768',
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


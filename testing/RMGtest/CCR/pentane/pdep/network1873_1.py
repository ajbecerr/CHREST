species(
    label = 'C=C([O])C[CH]C(C)O[O](7502)',
    structure = SMILES('C=C([O])C[CH]C(C)O[O]'),
    E0 = (110.618,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.217149,0.0953578,-0.000105108,6.38701e-08,-1.58255e-11,13454,37.4568], Tmin=(100,'K'), Tmax=(974.146,'K')), NASAPolynomial(coeffs=[13.8458,0.0376131,-1.61921e-05,3.01948e-09,-2.09078e-13,10714.1,-30.0202], Tmin=(974.146,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.618,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJCOOH)"""),
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
    label = '[CH2]C(=O)C[CH]CO[O](2201)',
    structure = SMILES('C=C([O])C[CH]CO[O]'),
    E0 = (147.041,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,192.599,192.888,192.948,3131.87],'cm^-1')),
        HinderedRotor(inertia=(0.00452505,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00451733,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.540656,'amu*angstrom^2'), symmetry=1, barrier=(14.2812,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.539218,'amu*angstrom^2'), symmetry=1, barrier=(14.2798,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4356.36,'J/mol'), sigma=(7.16244,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.45 K, Pc=26.9 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.391456,0.0832265,-0.000103837,7.25835e-08,-2.06104e-11,17811.5,33.3718], Tmin=(100,'K'), Tmax=(856.345,'K')), NASAPolynomial(coeffs=[11.385,0.0318818,-1.39108e-05,2.58448e-09,-1.77519e-13,15928.4,-17.9621], Tmin=(856.345,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(147.041,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJCOOH)"""),
)

species(
    label = 'C=C([O])CC(C)[CH]O[O](8094)',
    structure = SMILES('C=C([O])CC(C)[CH]O[O]'),
    E0 = (102.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.205813,0.0954998,-0.000105998,6.56781e-08,-1.66556e-11,12463.3,35.3264], Tmin=(100,'K'), Tmax=(951.521,'K')), NASAPolynomial(coeffs=[13.2713,0.0388441,-1.66837e-05,3.1009e-09,-2.14099e-13,9898.53,-29.0231], Tmin=(951.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(102.387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCsJOOH)"""),
)

species(
    label = '[CH2]C(C(=C)[O])C(C)O[O](7503)',
    structure = SMILES('[CH2]C(C(=C)[O])C(C)O[O]'),
    E0 = (107.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4529.55,'J/mol'), sigma=(7.53448,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=707.51 K, Pc=24.03 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.471862,0.095412,-0.000101852,5.84892e-08,-1.33972e-11,13070.4,37.5922], Tmin=(100,'K'), Tmax=(1064.42,'K')), NASAPolynomial(coeffs=[16.5537,0.0314311,-1.16888e-05,2.01802e-09,-1.33795e-13,9445.97,-45.6094], Tmin=(1064.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Isobutyl)"""),
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
    label = 'C=C([O])CC=CC(2558)',
    structure = SMILES('C=C([O])CC=CC'),
    E0 = (-24.0296,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,234.646,234.758,234.823],'cm^-1')),
        HinderedRotor(inertia=(0.00306028,'amu*angstrom^2'), symmetry=1, barrier=(0.119628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251283,'amu*angstrom^2'), symmetry=1, barrier=(9.83273,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.251422,'amu*angstrom^2'), symmetry=1, barrier=(9.83137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.885366,0.0592119,-2.96115e-05,-2.80559e-09,5.1329e-12,-2769.91,27.3679], Tmin=(100,'K'), Tmax=(1031.28,'K')), NASAPolynomial(coeffs=[13.5743,0.0280751,-1.06196e-05,1.91648e-09,-1.32745e-13,-6348.48,-38.9009], Tmin=(1031.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.0296,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(369.994,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
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
    label = 'C=C([O])C[CH]C(C)[O](8095)',
    structure = SMILES('C=C([O])C[CH]C(C)[O]'),
    E0 = (116.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.116346,0.075773,-5.68822e-05,1.6754e-08,-9.28192e-14,14214.8,34.7276], Tmin=(100,'K'), Tmax=(1038.55,'K')), NASAPolynomial(coeffs=[16.9727,0.0271161,-1.00991e-05,1.80357e-09,-1.24162e-13,9836.32,-51.4556], Tmin=(1038.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(390.78,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(CCJCO)"""),
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
    label = '[CH]CC(=C)[O](2374)',
    structure = SMILES('[CH]CC(=C)[O]'),
    E0 = (354.94,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,180,818.766,818.822,818.836,818.871],'cm^-1')),
        HinderedRotor(inertia=(0.00156519,'amu*angstrom^2'), symmetry=1, barrier=(14.8707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.646792,'amu*angstrom^2'), symmetry=1, barrier=(14.871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.0819,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.96183,0.0381144,-1.88967e-05,-6.94609e-09,6.76808e-12,42768.9,19.9558], Tmin=(100,'K'), Tmax=(946.634,'K')), NASAPolynomial(coeffs=[11.9365,0.0122825,-3.81808e-06,6.42322e-10,-4.44595e-14,40149.3,-31.4806], Tmin=(946.634,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)[O](2372)',
    structure = SMILES('[CH2]C(=C)[O]'),
    E0 = (110.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1477.83],'cm^-1')),
        HinderedRotor(inertia=(0.530917,'amu*angstrom^2'), symmetry=1, barrier=(12.2068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.69661,0.0242908,-7.63704e-06,-9.31802e-09,6.21014e-12,13330.3,14.3522], Tmin=(100,'K'), Tmax=(924.139,'K')), NASAPolynomial(coeffs=[8.66416,0.00984172,-2.65662e-06,4.14898e-10,-2.77492e-14,11741.3,-16.5963], Tmin=(924.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
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
    label = 'C=[C]C[CH]C(C)O[O](8096)',
    structure = SMILES('C=[C]C[CH]C(C)O[O]'),
    E0 = (425.214,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.387217,0.0851104,-9.92468e-05,7.0305e-08,-2.11121e-11,51266.4,34.207], Tmin=(100,'K'), Tmax=(798.292,'K')), NASAPolynomial(coeffs=[8.85606,0.042676,-1.95128e-05,3.71854e-09,-2.59514e-13,49914.2,-4.74253], Tmin=(798.292,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.214,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJCOOH) + radical(Cds_S)"""),
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
    label = 'C=C([O])C[C]C(C)O[O](8097)',
    structure = SMILES('C=C([O])C[C]C(C)O[O]'),
    E0 = (358.418,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1066.67,1333.33,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.353867,0.0971226,-0.000111464,6.8181e-08,-1.67118e-11,43263.5,35.5453], Tmin=(100,'K'), Tmax=(992.238,'K')), NASAPolynomial(coeffs=[15.6747,0.032507,-1.37826e-05,2.55091e-09,-1.75981e-13,40082.6,-41.6585], Tmin=(992.238,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(358.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(386.623,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C([O])CC1OOC1C(8098)',
    structure = SMILES('C=C([O])CC1OOC1C'),
    E0 = (-85.0895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.180214,0.0827304,-6.72892e-05,2.83033e-08,-4.79533e-12,-10075.8,29.7189], Tmin=(100,'K'), Tmax=(1404.53,'K')), NASAPolynomial(coeffs=[17.9033,0.0312288,-1.2286e-05,2.19533e-09,-1.48135e-13,-15155.5,-63.6664], Tmin=(1404.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-85.0895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1CC(C(C)O[O])O1(7994)',
    structure = SMILES('C=C1CC(C(C)O[O])O1'),
    E0 = (-64.4016,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.138691,0.0739067,-1.8732e-05,-4.18741e-08,2.58207e-11,-7580.84,27.7939], Tmin=(100,'K'), Tmax=(920.576,'K')), NASAPolynomial(coeffs=[22.4014,0.0206324,-4.70338e-06,6.70808e-10,-4.60136e-14,-13623.4,-89.363], Tmin=(920.576,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.4016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(ROOJ)"""),
)

species(
    label = 'C=C([O])CC1OC1C(8099)',
    structure = SMILES('C=C([O])CC1OC1C'),
    E0 = (-136.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.445384,0.0653234,-1.39998e-05,-4.16954e-08,2.6282e-11,-16263.1,27.2509], Tmin=(100,'K'), Tmax=(874.334,'K')), NASAPolynomial(coeffs=[18.4894,0.0209846,-3.48729e-06,2.73349e-10,-1.0308e-14,-20879,-65.7303], Tmin=(874.334,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-136.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=C(C)OJ)"""),
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
    label = 'C=C[CH]CC(=C)[O](4317)',
    structure = SMILES('[CH2]C=CCC(=C)[O]'),
    E0 = (127.47,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,3000,3100,440,815,1455,1000,367.998,368.338,368.342],'cm^-1')),
        HinderedRotor(inertia=(0.00124219,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170813,'amu*angstrom^2'), symmetry=1, barrier=(16.4576,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170662,'amu*angstrom^2'), symmetry=1, barrier=(16.4578,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.1271,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3874.43,'J/mol'), sigma=(6.44487,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.18 K, Pc=32.84 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.924688,0.0562639,-1.95361e-05,-1.80404e-08,1.18946e-11,15451.9,27.5064], Tmin=(100,'K'), Tmax=(976.244,'K')), NASAPolynomial(coeffs=[15.5,0.0226347,-7.95274e-06,1.42522e-09,-1.00663e-13,11362.9,-48.8281], Tmin=(976.244,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = 'C=C([O])CC=C(C)OO(8100)',
    structure = SMILES('C=C([O])CC=C(C)OO'),
    E0 = (-114.093,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0650243,0.0891863,-8.68419e-05,4.52419e-08,-9.60909e-12,-13575.7,34.8832], Tmin=(100,'K'), Tmax=(1125.79,'K')), NASAPolynomial(coeffs=[15.0486,0.0354867,-1.52927e-05,2.87214e-09,-2.00191e-13,-16978.6,-39.8221], Tmin=(1125.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-114.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C=CC(C)O[O](8101)',
    structure = SMILES('C=C(O)C=CC(C)O[O]'),
    E0 = (-127.94,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25529,0.100706,-0.000101922,5.12338e-08,-9.94455e-12,-15185.8,35.6126], Tmin=(100,'K'), Tmax=(1301.99,'K')), NASAPolynomial(coeffs=[24.481,0.0195941,-6.12052e-06,9.73952e-10,-6.24441e-14,-21714.2,-94.6764], Tmin=(1301.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-127.94,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = 'C=C([O])C=CC(C)OO(8102)',
    structure = SMILES('C=C([O])C=CC(C)OO'),
    E0 = (-142.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.922401,0.0976278,-9.53815e-05,4.69298e-08,-9.05148e-12,-16909.3,35.1584], Tmin=(100,'K'), Tmax=(1265.97,'K')), NASAPolynomial(coeffs=[22.0367,0.0250848,-9.42744e-06,1.66559e-09,-1.12771e-13,-22722.3,-81.0205], Tmin=(1265.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-142.14,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)CC=C(C)O[O](7517)',
    structure = SMILES('C=C(O)CC=C(C)O[O]'),
    E0 = (-99.8934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.363685,0.0919234,-9.24408e-05,4.8599e-08,-1.01881e-11,-11853.8,35.2098], Tmin=(100,'K'), Tmax=(1156.49,'K')), NASAPolynomial(coeffs=[17.5974,0.0298006,-1.18656e-05,2.15078e-09,-1.47316e-13,-16008.2,-54.0533], Tmin=(1156.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-99.8934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ)"""),
)

species(
    label = '[CH2]C([O])[CH][CH]C(C)O[O](8103)',
    structure = SMILES('[CH2]C([O])[CH][CH]C(C)O[O]'),
    E0 = (528.971,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.382853,0.10151,-0.000122482,8.22151e-08,-2.25291e-11,63773.9,41.3205], Tmin=(100,'K'), Tmax=(883.731,'K')), NASAPolynomial(coeffs=[13.2134,0.0399706,-1.80275e-05,3.41751e-09,-2.38075e-13,61370.8,-22.593], Tmin=(883.731,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.971,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(ROOJ) + radical(CCJCOOH) + radical(CCJCO) + radical(CJCO)"""),
)

species(
    label = 'CC([CH]C[C]1CO1)O[O](8104)',
    structure = SMILES('CC([CH]C[C]1CO1)O[O]'),
    E0 = (256.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.411047,0.100472,-0.000129158,9.56154e-08,-2.80433e-11,31040,36.6215], Tmin=(100,'K'), Tmax=(945.102,'K')), NASAPolynomial(coeffs=[11.5763,0.0388689,-1.41383e-05,2.31417e-09,-1.44596e-13,29259.5,-17.9661], Tmin=(945.102,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.788,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(ROOJ) + radical(C2CsJO) + radical(CCJCOOH)"""),
)

species(
    label = 'CC(O[O])C1C[C]([O])C1(8105)',
    structure = SMILES('CC(O[O])C1C[C]([O])C1'),
    E0 = (242.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.229362,0.0819823,-6.41459e-05,2.56149e-08,-4.10767e-12,29366.9,34.0456], Tmin=(100,'K'), Tmax=(1480.63,'K')), NASAPolynomial(coeffs=[18.8236,0.0305103,-1.20009e-05,2.1363e-09,-1.43415e-13,23724.8,-65.3516], Tmin=(1480.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(ROOJ) + radical(C2CsJOH)"""),
)

species(
    label = 'CC1[CH]C[C]([O])COO1(8106)',
    structure = SMILES('CC1[CH]C[C]([O])COO1'),
    E0 = (232.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0215571,0.0429937,0.000122905,-2.27248e-07,1.01288e-10,28115.2,33.7044], Tmin=(100,'K'), Tmax=(910.831,'K')), NASAPolynomial(coeffs=[38.9625,-0.00899271,1.21883e-05,-2.51045e-09,1.60821e-13,16068.5,-177.876], Tmin=(910.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.205,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(428.195,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cycloheptane) + radical(CC(C)OJ) + radical(CCJCOOH) + radical(C2CsJOH)"""),
)

species(
    label = '[CH2]C1([O])CC1C(C)O[O](8063)',
    structure = SMILES('[CH2]C1([O])CC1C(C)O[O]'),
    E0 = (263.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.601753,0.0946432,-9.129e-05,4.54172e-08,-8.98615e-12,31827.8,33.0965], Tmin=(100,'K'), Tmax=(1223.84,'K')), NASAPolynomial(coeffs=[19.2278,0.0298321,-1.18539e-05,2.1456e-09,-1.46817e-13,26974.2,-66.575], Tmin=(1223.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.208,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-CsH) + group(O2s-OsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(ROOJ) + radical(CJC(C)2O)"""),
)

species(
    label = '[CH2]C1([O])C[CH]C(C)OO1(8107)',
    structure = SMILES('[CH2]C1([O])C[CH]C(C)OO1'),
    E0 = (204.899,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0671129,0.0650163,1.5852e-05,-8.27798e-08,4.19653e-11,24805.3,31.2298], Tmin=(100,'K'), Tmax=(907.654,'K')), NASAPolynomial(coeffs=[23.4899,0.0189882,-2.60298e-06,2.00833e-10,-1.27617e-14,18197.4,-92.4804], Tmin=(907.654,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.899,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxane) + radical(CC(C)(O)OJ) + radical(CCJCOOH) + radical(CJCOOH)"""),
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
    label = 'C=C([O])CC=CO[O](3038)',
    structure = SMILES('C=C([O])CC=CO[O]'),
    E0 = (79.7029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,350,440,435,1725,2950,3100,1380,975,1025,1650,236.266,236.266,236.266],'cm^-1')),
        HinderedRotor(inertia=(0.22083,'amu*angstrom^2'), symmetry=1, barrier=(8.74759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22083,'amu*angstrom^2'), symmetry=1, barrier=(8.74759,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.22083,'amu*angstrom^2'), symmetry=1, barrier=(8.74759,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.459127,0.0712632,-6.80877e-05,3.05547e-08,-4.4767e-12,9719.38,31.7509], Tmin=(100,'K'), Tmax=(992.009,'K')), NASAPolynomial(coeffs=[16.4305,0.0186644,-6.39843e-06,1.08943e-09,-7.32521e-14,5969.95,-48.1004], Tmin=(992.009,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(79.7029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ)"""),
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
    label = 'C=C([O])C[CH][CH]C(2560)',
    structure = SMILES('C=C([O])C[CH][CH]C'),
    E0 = (248.069,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,2950,3100,1380,975,1025,1650,222.408,1946.36,1946.36,1946.36],'cm^-1')),
        HinderedRotor(inertia=(0.20821,'amu*angstrom^2'), symmetry=1, barrier=(7.30844,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00340802,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208215,'amu*angstrom^2'), symmetry=1, barrier=(7.30845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00340791,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.1351,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3795.8,'J/mol'), sigma=(6.58551,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.90 K, Pc=30.16 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02514,0.0519567,-1.23522e-06,-9.02755e-08,8.70563e-11,29898,28.692], Tmin=(100,'K'), Tmax=(467.454,'K')), NASAPolynomial(coeffs=[4.96563,0.0420012,-1.80839e-05,3.34267e-09,-2.2893e-13,29457,14.9649], Tmin=(467.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(365.837,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C=C([O])C[CH][CH]O[O](3040)',
    structure = SMILES('C=C([O])C[CH][CH]O[O]'),
    E0 = (335.624,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,2950,3100,1380,975,1025,1650,342.149,343.271,344.524,1332.43],'cm^-1')),
        HinderedRotor(inertia=(0.0939153,'amu*angstrom^2'), symmetry=1, barrier=(7.85382,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0941713,'amu*angstrom^2'), symmetry=1, barrier=(7.87926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0938533,'amu*angstrom^2'), symmetry=1, barrier=(7.86464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00142827,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (114.099,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.263024,0.0876436,-0.000127135,1.01888e-07,-3.2548e-11,40495.6,33.9741], Tmin=(100,'K'), Tmax=(835.222,'K')), NASAPolynomial(coeffs=[10.8878,0.0306165,-1.36849e-05,2.52699e-09,-1.71067e-13,38935.1,-14.0886], Tmin=(835.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(335.624,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CCJCOOH) + radical(CCsJOOH)"""),
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
    label = 'C=C([O])CC[C](C)O[O](8114)',
    structure = SMILES('C=C([O])CC[C](C)O[O]'),
    E0 = (97.1024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.126892,0.0960493,-0.000114694,7.97469e-08,-2.29768e-11,11822.6,35.1253], Tmin=(100,'K'), Tmax=(837.913,'K')), NASAPolynomial(coeffs=[11.0815,0.0425472,-1.89239e-05,3.55542e-09,-2.46035e-13,9944.1,-16.9676], Tmin=(837.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.1024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C2CsJOOH)"""),
)

species(
    label = 'C=C([O])[CH]CC(C)O[O](8044)',
    structure = SMILES('[CH2]C([O])=CCC(C)O[O]'),
    E0 = (56.8746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.379714,0.0959358,-0.000103877,6.07905e-08,-1.43047e-11,6998.54,35.98], Tmin=(100,'K'), Tmax=(1031.52,'K')), NASAPolynomial(coeffs=[15.6277,0.0338639,-1.36158e-05,2.4564e-09,-1.67119e-13,3696.08,-41.7437], Tmin=(1031.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.8746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C([O])C[CH][C](C)OO(8115)',
    structure = SMILES('C=C([O])C[CH][C](C)OO'),
    E0 = (145.513,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.312101,0.100981,-0.000123238,8.57963e-08,-2.46493e-11,17650.9,37.3397], Tmin=(100,'K'), Tmax=(840.817,'K')), NASAPolynomial(coeffs=[11.8765,0.0429973,-1.97975e-05,3.78156e-09,-2.64252e-13,15601.2,-19.3503], Tmin=(840.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(145.513,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C2CsJOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C(CCC(=C)[O])O[O](8116)',
    structure = SMILES('[CH2]C(CCC(=C)[O])O[O]'),
    E0 = (124.166,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.382153,0.0969866,-0.000105883,6.21687e-08,-1.47032e-11,15091.1,37.7351], Tmin=(100,'K'), Tmax=(1024.69,'K')), NASAPolynomial(coeffs=[15.6248,0.0345,-1.44097e-05,2.65448e-09,-1.82836e-13,11810.7,-39.8794], Tmin=(1024.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(124.166,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(CJCOOH)"""),
)

species(
    label = 'C=C(O)[CH][CH]C(C)O[O](8117)',
    structure = SMILES('[CH2]C(O)=C[CH]C(C)O[O]'),
    E0 = (35.9864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.25206,0.107311,-0.000116348,6.35681e-08,-1.3531e-11,4524.3,34.8043], Tmin=(100,'K'), Tmax=(1155.07,'K')), NASAPolynomial(coeffs=[22.6532,0.0245272,-8.84287e-06,1.51959e-09,-1.01355e-13,-998.146,-83.9712], Tmin=(1155.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(35.9864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(ROOJ) + radical(C=CCJCO) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH]=C(O)C[CH]C(C)O[O](8118)',
    structure = SMILES('[CH]=C(O)C[CH]C(C)O[O]'),
    E0 = (219.91,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,200,800],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.763542,0.104181,-0.00011919,7.09297e-08,-1.67278e-11,26621.2,38.3146], Tmin=(100,'K'), Tmax=(1035.27,'K')), NASAPolynomial(coeffs=[18.2684,0.0306462,-1.26459e-05,2.31969e-09,-1.59617e-13,22680.5,-54.1633], Tmin=(1035.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(219.91,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJCOOH) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C([CH]CC(=C)[O])OO(8119)',
    structure = SMILES('[CH2]C([CH]CC(=C)[O])OO'),
    E0 = (172.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.551299,0.10174,-0.00011385,6.75256e-08,-1.60994e-11,20918.7,39.8913], Tmin=(100,'K'), Tmax=(1016.58,'K')), NASAPolynomial(coeffs=[16.4282,0.0349294,-1.52685e-05,2.87656e-09,-2.00682e-13,17466.5,-42.3039], Tmin=(1016.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(172.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJCOOH) + radical(CJCOOH)"""),
)

species(
    label = '[CH]=C([O])CCC(C)O[O](8120)',
    structure = SMILES('[CH]=C([O])CCC(C)O[O]'),
    E0 = (157.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.353083,0.0973213,-0.000107629,6.4562e-08,-1.56585e-11,19074.2,36.0407], Tmin=(100,'K'), Tmax=(998.337,'K')), NASAPolynomial(coeffs=[15.0095,0.035769,-1.51474e-05,2.80502e-09,-1.93605e-13,16006.8,-38.0492], Tmin=(998.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(157.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])[CH][CH]C(C)OO(8121)',
    structure = SMILES('[CH2]C([O])=C[CH]C(C)OO'),
    E0 = (21.7865,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.947324,0.104514,-0.000110586,6.00506e-08,-1.29e-11,2802.17,34.455], Tmin=(100,'K'), Tmax=(1133.99,'K')), NASAPolynomial(coeffs=[20.1301,0.0301662,-1.22413e-05,2.2339e-09,-1.5363e-13,-1978.13,-69.8816], Tmin=(1133.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(21.7865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CCJCO) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C(O)C[CH][C](C)O[O](8122)',
    structure = SMILES('C=C(O)C[CH][C](C)O[O]'),
    E0 = (159.713,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.437843,0.101697,-0.000121834,8.00996e-08,-2.13408e-11,19365.3,37.0449], Tmin=(100,'K'), Tmax=(911.006,'K')), NASAPolynomial(coeffs=[14.1119,0.0378122,-1.66452e-05,3.12256e-09,-2.16385e-13,16714.3,-31.7927], Tmin=(911.006,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.713,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ) + radical(C2CsJOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C([CH]CC(=C)O)O[O](8123)',
    structure = SMILES('[CH2]C([CH]CC(=C)O)O[O]'),
    E0 = (186.776,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.803827,0.103973,-0.000117862,6.90445e-08,-1.59767e-11,22638.5,40.0497], Tmin=(100,'K'), Tmax=(1056.33,'K')), NASAPolynomial(coeffs=[18.9035,0.0293471,-1.18924e-05,2.16566e-09,-1.48574e-13,18475,-56.1067], Tmin=(1056.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.776,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJCOOH) + radical(CJCOOH)"""),
)

species(
    label = '[CH]=C([O])C[CH]C(C)OO(8124)',
    structure = SMILES('[CH]=C([O])C[CH]C(C)OO'),
    E0 = (205.71,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.52556,0.102112,-0.000115717,7.00648e-08,-1.7113e-11,24902,38.2089], Tmin=(100,'K'), Tmax=(991.601,'K')), NASAPolynomial(coeffs=[15.8148,0.0361961,-1.60054e-05,3.02701e-09,-2.11449e-13,21661.4,-40.4857], Tmin=(991.601,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJCOOH) + radical(Cds_P)"""),
)

species(
    label = 'CC([CH]CC[C]=O)O[O](7506)',
    structure = SMILES('CC([CH]CC[C]=O)O[O]'),
    E0 = (136.751,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4289.48,'J/mol'), sigma=(7.22136,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=670.01 K, Pc=25.85 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.384298,0.105016,-0.000151627,1.28349e-07,-4.38048e-11,16597,38.8648], Tmin=(100,'K'), Tmax=(821.08,'K')), NASAPolynomial(coeffs=[9.06604,0.0463933,-2.15415e-05,4.06142e-09,-2.78636e-13,15469.3,-2.28145], Tmin=(821.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.751,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJCOOH) + radical(CCCJ=O)"""),
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
    label = 'CC([CH]C[C]=O)O[O](7937)',
    structure = SMILES('CC([CH]C[C]=O)O[O]'),
    E0 = (160.531,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,279.79,2119.31],'cm^-1')),
        HinderedRotor(inertia=(0.156704,'amu*angstrom^2'), symmetry=1, barrier=(8.7193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15893,'amu*angstrom^2'), symmetry=1, barrier=(8.72451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156998,'amu*angstrom^2'), symmetry=1, barrier=(8.72245,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156015,'amu*angstrom^2'), symmetry=1, barrier=(8.71621,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155709,'amu*angstrom^2'), symmetry=1, barrier=(8.73123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (115.107,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.285769,0.0899431,-0.000136653,1.18524e-07,-4.07154e-11,19433.3,34.2145], Tmin=(100,'K'), Tmax=(839.574,'K')), NASAPolynomial(coeffs=[8.10315,0.038049,-1.7765e-05,3.33781e-09,-2.27706e-13,18636.9,0.941998], Tmin=(839.574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(160.531,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJCOOH) + radical(CCCJ=O)"""),
)

species(
    label = 'CC(O[O])C1CC(=O)C1(7509)',
    structure = SMILES('CC(O[O])C1CC(=O)C1'),
    E0 = (-121.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331041,0.070236,-4.07964e-05,8.53747e-09,2.55687e-13,-14471.5,32.0376], Tmin=(100,'K'), Tmax=(1250.98,'K')), NASAPolynomial(coeffs=[15.0006,0.035918,-1.47409e-05,2.69591e-09,-1.84426e-13,-19126.7,-45.9565], Tmin=(1250.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-121.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(419.881,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)C=CC(C)O[O](8125)',
    structure = SMILES('CC(=O)C=CC(C)O[O]'),
    E0 = (-134.462,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.840797,0.0763628,-6.18457e-05,2.79808e-08,-5.61276e-12,-16064.1,32.2669], Tmin=(100,'K'), Tmax=(1111.69,'K')), NASAPolynomial(coeffs=[8.93228,0.047248,-2.25603e-05,4.42138e-09,-3.14537e-13,-17863.1,-7.62637], Tmin=(1111.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-134.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + radical(ROOJ)"""),
)

species(
    label = 'CC(=O)CC=C(C)O[O](7518)',
    structure = SMILES('CC(=O)CC=C(C)O[O]'),
    E0 = (-122.094,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14413,0.0720395,-5.1622e-05,1.89751e-08,-3.05443e-12,-14590.1,30.164], Tmin=(100,'K'), Tmax=(1320.66,'K')), NASAPolynomial(coeffs=[9.45203,0.0468766,-2.3042e-05,4.54791e-09,-3.2336e-13,-16784.5,-12.2275], Tmin=(1320.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.094,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + radical(ROOJ)"""),
)

species(
    label = '[CH2][C]1CC(C(C)O[O])O1(8126)',
    structure = SMILES('[CH2][C]1CC(C(C)O[O])O1'),
    E0 = (249.307,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.00958996,0.0939336,-9.48224e-05,3.49941e-08,6.68848e-12,30124,31.5494], Tmin=(100,'K'), Tmax=(654.311,'K')), NASAPolynomial(coeffs=[12.9785,0.0362915,-1.25589e-05,1.99937e-09,-1.23045e-13,27958.6,-29.1608], Tmin=(654.311,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(415.724,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(ROOJ) + radical(C2CsJOCs) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]1C[CH]C(C)OOO1(8127)',
    structure = SMILES('[CH2][C]1C[CH]C(C)OOO1'),
    E0 = (421.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.357235,0.0284814,0.000168431,-2.78194e-07,1.2043e-10,50887.1,36.6249], Tmin=(100,'K'), Tmax=(910.658,'K')), NASAPolynomial(coeffs=[40.6271,-0.0133778,1.49745e-05,-3.03695e-09,1.94369e-13,37954,-184.625], Tmin=(910.658,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.607,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(424.038,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cycloheptane) + radical(CCJCOOH) + radical(C2CsJOO) + radical(CJCOOH)"""),
)

species(
    label = 'CC(=O)[CH][CH]C(C)O[O](8128)',
    structure = SMILES('CC([O])=C[CH]C(C)O[O]'),
    E0 = (14.8746,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.425456,0.0977094,-0.00010255,5.76204e-08,-1.30979e-11,1948.13,32.6558], Tmin=(100,'K'), Tmax=(1061.13,'K')), NASAPolynomial(coeffs=[15.9278,0.0360643,-1.54082e-05,2.87249e-09,-1.99297e-13,-1522.43,-47.2094], Tmin=(1061.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(14.8746,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(411.566,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(C=CCJCO)"""),
)

species(
    label = 'C[C]([CH]CC(C)=O)O[O](8129)',
    structure = SMILES('C[C]([CH]CC(C)=O)O[O]'),
    E0 = (132.599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,375,552.5,462.5,1710,3025,407.5,1350,352.5,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.204835,0.105355,-0.000165429,1.53889e-07,-5.61365e-11,16086.9,37.0403], Tmin=(100,'K'), Tmax=(838.093,'K')), NASAPolynomial(coeffs=[4.87263,0.0543347,-2.61722e-05,4.98094e-09,-3.42045e-13,16176.6,19.0537], Tmin=(838.093,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(ROOJ) + radical(C2CsJOOH) + radical(CCJCOOH)"""),
)

species(
    label = '[CH2]C([CH]CC(C)=O)O[O](8130)',
    structure = SMILES('[CH2]C([CH]CC(C)=O)O[O]'),
    E0 = (159.662,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (129.134,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.38412,0.105246,-0.000152282,1.29806e-07,-4.47417e-11,19352.3,39.3877], Tmin=(100,'K'), Tmax=(812.841,'K')), NASAPolynomial(coeffs=[8.84424,0.0473143,-2.22741e-05,4.22977e-09,-2.91542e-13,18265.7,-0.6776], Tmin=(812.841,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.662,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsCs) + radical(ROOJ) + radical(CCJCOOH) + radical(CJCOOH)"""),
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
    E0 = (110.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (565.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (305.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (356.219,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (110.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (359.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (573.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (545.96,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (944.567,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (570.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (118.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (118.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (240.322,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (253.293,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (174.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (174.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (135.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (135.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (551.833,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (341.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (242.829,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (232.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (263.208,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (228.037,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (260.993,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (242.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (243.006,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (219.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (297.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (110.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (239.447,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (472.174,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (509.309,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (473.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (385.583,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (536.373,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (569.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (265.008,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (235.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (280.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (271.024,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (267.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (412.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (253.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (201.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (171.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (345.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (260.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (357.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (382.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (598.333,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (118.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (174.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (135.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (249.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (421.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (227.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (301.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (212.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['ketene(1375)', 'C=CC(C)O[O](3748)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2]C(=O)C[CH]CO[O](2201)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C([O])CC(C)[CH]O[O](8094)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.59192e+09,'s^-1'), n=1.025, Ea=(194.765,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HR!H)CJ;CsJ;CH3] for rate rule [cCs(-HR!H)CJ;CsJ-CsH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['[CH2]C(C(=C)[O])C(C)O[O](7503)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HH)CJ;CsJ;C] for rate rule [cCs(-HH)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['O2(S)(666)', 'C=C([O])CC=CC(2558)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', 'C=C([O])C[CH]C(C)[O](8095)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C[CH]O[O](57)', '[CH]CC(=C)[O](2374)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/CsO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=C)[O](2372)', '[CH]C(C)O[O](668)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.68328e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['O(4)', 'C=[C]C[CH]C(C)O[O](8096)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'C=C([O])C[C]C(C)O[O](8097)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C([O])CC1OOC1C(8098)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Cpri_rad_out_single] for rate rule [R4_SSS;O_rad;Cpri_rad_out_H/NonDeC]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C1CC(C(C)O[O])O1(7994)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['O(4)', 'C=C([O])CC1OC1C(8099)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(3.27e+09,'s^-1'), n=1.06, Ea=(129.704,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 36 used for R2OO_S;C_rad/H/NonDeC_intra;OOJ
Exact match found for rate rule [R2OO_S;C_rad/H/NonDeC_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['HO2(9)', 'C=C[CH]CC(=C)[O](4317)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.6e+09,'s^-1','*|/',2.51189), n=1.2, Ea=(142.674,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2OO_2H
Exact match found for rate rule [R2OO_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C([O])CC=C(C)OO(8100)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C(O)C=CC(C)O[O](8101)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C([O])C=CC(C)OO(8102)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C(O)CC=C(C)O[O](7517)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([O])[CH][CH]C(C)O[O](8103)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['CC([CH]C[C]1CO1)O[O](8104)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['CC(O[O])C1C[C]([O])C1(8105)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(5.25757e+07,'s^-1'), n=1.165, Ea=(132.21,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 129.8 to 132.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['CC1[CH]C[C]([O])COO1(8106)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.62196e+06,'s^-1'), n=0.867572, Ea=(121.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra;radadd_intra] for rate rule [R7_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 115.5 to 121.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['[CH2]C1([O])CC1C(C)O[O](8063)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(152.59,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 152.3 to 152.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['[CH2]C1([O])C[CH]C(C)OO1(8107)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.31936e+10,'s^-1'), n=0.385029, Ea=(117.418,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R7;doublebond_intra;radadd_intra_O] + [R7;doublebond_intra_2H;radadd_intra] for rate rule [R7;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', 'C=C([O])CC=C(C)O[O](8108)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(170.641,'m^3/(mol*s)'), n=1.56204, Ea=(11.2897,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;HJ] for rate rule [Cds-OsCs_Cds;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', 'C=C([O])C=CC(C)O[O](8109)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.0272924,'m^3/(mol*s)'), n=2.81111, Ea=(21.1569,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 26 used for Cds-CdH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['CH3(17)', 'C=C([O])CC=CO[O](3038)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(0.00632077,'m^3/(mol*s)'), n=2.46811, Ea=(26.7526,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CsJ-HHH] for rate rule [Cds-OsH_Cds-CsH;CsJ-HHH]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2][C]=O(1376)', 'C=CC(C)O[O](3748)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(532000,'cm^3/(mol*s)'), n=1.85, Ea=(23.8906,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\O2s/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction29',
    reactants = ['ketene(1375)', '[CH2][CH]C(C)O[O](3749)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['O2(2)', 'C=C([O])CC=CC(2558)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(9.03779e-05,'m^3/(mol*s)'), n=3.085, Ea=(143.27,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cds-CsH_Cds;O2b] + [Cds-CsH_Cds-CsH;OJ] for rate rule [Cds-CsH_Cds-CsH;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 141.0 to 143.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['O2(2)', 'C=C([O])C[CH][CH]C(2560)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
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
    reactants = ['CH3(17)', 'C=C([O])C[CH][CH]O[O](3040)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', 'C=C([O])C[CH][C](C)O[O](8110)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][C]=O(1376)', '[CH2][CH]C(C)O[O](3749)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(6.08474e+07,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['H(3)', 'C=C([O])[CH][CH]C(C)O[O](8111)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', '[CH2]C([CH]CC(=C)[O])O[O](8112)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(5250.69,'m^3/(mol*s)'), n=1.27262, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)', '[CH]=C([O])C[CH]C(C)O[O](8113)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C([O])CC[C](C)O[O](8114)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NonDe] for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_NDMustO]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C([O])[CH]CC(C)O[O](8044)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.682e+10,'s^-1'), n=0.35, Ea=(125.102,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 160 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C([O])C[CH][C](C)OO(8115)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(40813.3,'s^-1'), n=2.17068, Ea=(135.172,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(CCC(=C)[O])O[O](8116)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.18e+10,'s^-1'), n=0.82, Ea=(146.858,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 186 used for R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R3H_SS_Cs;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C(O)[CH][CH]C(C)O[O](8117)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.39846e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=C(O)C[CH]C(C)O[O](8118)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['[CH2]C([CH]CC(=C)[O])OO(8119)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.989e+11,'s^-1'), n=0.15, Ea=(143.135,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 246 used for R4H_SSS_OCs;O_rad_out;Cs_H_out_2H
Exact match found for rate rule [R4H_SSS_OCs;O_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=C([O])CCC(C)O[O](8120)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['C=C([O])[CH][CH]C(C)OO(8121)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(61085.3,'s^-1'), n=1.8624, Ea=(60.8997,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=C(O)C[CH][C](C)O[O](8122)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(2.46277e+08,'s^-1'), n=1.29592, Ea=(185.417,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_NonDe;XH_out] for rate rule [R5HJ_1;C_rad_out_NDMustO;O_H_out]
Euclidian distance = 2.449489742783178
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C([CH]CC(=C)O)O[O](8123)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_2;C_rad_out_2H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=C([O])C[CH]C(C)OO(8124)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R7HJ_3;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CC([CH]CC[C]=O)O[O](7506)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction51',
    reactants = ['CH2(T)(20)', 'CC([CH]C[C]=O)O[O](7937)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['CC(O[O])C1CC(=O)C1(7509)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['CC(=O)C=CC(C)O[O](8125)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['CC(=O)CC=C(C)O[O](7518)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(2.1261e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['[CH2][C]1CC(C(C)O[O])O1(8126)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(4.7342e+08,'s^-1'), n=0.920995, Ea=(138.96,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction56',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['[CH2][C]1C[CH]C(C)OOO1(8127)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(5.649e+12,'s^-1'), n=0.287, Ea=(310.989,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;carbonyl_intra;radadd_intra] for rate rule [R7_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 303.0 to 311.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=C([O])C[CH]C(C)O[O](7502)'],
    products = ['CC(=O)[CH][CH]C(C)O[O](8128)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(4.61991,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C[C]([CH]CC(C)=O)O[O](8129)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.44378e+09,'s^-1'), n=1.25282, Ea=(169.176,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_NonDe;Cs_H_out_2H] for rate rule [R5HJ_1;C_rad_out_NDMustO;Cs_H_out_2H]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2]C([CH]CC(C)=O)O[O](8130)'],
    products = ['C=C([O])C[CH]C(C)O[O](7502)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(5436.63,'s^-1'), n=1.865, Ea=(52.7184,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;Cs_H_out_2H] for rate rule [R6HJ_2;C_rad_out_2H;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1873',
    isomers = [
        'C=C([O])C[CH]C(C)O[O](7502)',
    ],
    reactants = [
        ('ketene(1375)', 'C=CC(C)O[O](3748)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1873',
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


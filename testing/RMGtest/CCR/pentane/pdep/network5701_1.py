species(
    label = '[CH2]C(=C)C1([CH2])OC1[O](24187)',
    structure = SMILES('[CH2]C(=C)C1([CH2])OC1[O]'),
    E0 = (299.843,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.26811,0.0781585,-6.92136e-05,2.06602e-08,3.57149e-12,36201.2,27.7447], Tmin=(100,'K'), Tmax=(821.339,'K')), NASAPolynomial(coeffs=[16.2773,0.0216813,-5.31496e-06,6.49292e-10,-3.34487e-14,32846.6,-50.7517], Tmin=(821.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1OC1[O](19061)',
    structure = SMILES('C=C1OC1[O]'),
    E0 = (43.0674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2950,3100,1380,975,1025,1650,577.055,577.113,577.621,577.636,577.829,578.051,579.439,580.062,580.523,4000],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3666.67,'J/mol'), sigma=(5.87705,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=572.73 K, Pc=40.99 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50293,0.0232972,9.07826e-06,-3.28256e-08,1.52265e-11,5242.42,12.4283], Tmin=(100,'K'), Tmax=(962.443,'K')), NASAPolynomial(coeffs=[12.1509,0.00682846,-2.0815e-06,4.13797e-10,-3.38067e-14,2290.92,-39.4337], Tmin=(962.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(43.0674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(182.918,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCOJ)"""),
)

species(
    label = 'allene(458)',
    structure = SMILES('C=C=C'),
    E0 = (175.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37447,0.0070462,2.78306e-05,-3.99445e-08,1.55729e-11,21188.6,7.62048], Tmin=(100,'K'), Tmax=(949.705,'K')), NASAPolynomial(coeffs=[6.79956,0.00959979,-3.02068e-06,5.37827e-10,-3.92606e-14,19772.3,-12.7582], Tmin=(949.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""allene""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'O(S)(3112)',
    structure = SMILES('O'),
    E0 = (432.331,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (15.9994,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,51997.4,2.99252], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,51997.4,2.99252], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.331,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""O(S)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2]C(=C)C(=C)[CH][O](23578)',
    structure = SMILES('[CH2]C(=C)C([CH2])=C[O]'),
    E0 = (254.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.43333,'amu*angstrom^2'), symmetry=1, barrier=(32.9551,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.4326,'amu*angstrom^2'), symmetry=1, barrier=(32.9383,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43095,'amu*angstrom^2'), symmetry=1, barrier=(32.9003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.40155,0.0603326,-6.20328e-06,-5.23668e-08,2.94114e-11,30745.2,23.9436], Tmin=(100,'K'), Tmax=(925.932,'K')), NASAPolynomial(coeffs=[23.7407,0.00806211,-1.82938e-07,-6.85501e-11,2.11033e-16,24341.7,-98.098], Tmin=(925.932,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)[C]1CC([O])O1(24995)',
    structure = SMILES('[CH2]C([CH2])=C1CC([O])O1'),
    E0 = (198.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.669643,0.0562606,2.66258e-06,-5.44257e-08,2.82959e-11,23968,24.7067], Tmin=(100,'K'), Tmax=(929.224,'K')), NASAPolynomial(coeffs=[19.9132,0.0168706,-3.88652e-06,5.90396e-10,-4.31643e-14,18516,-76.8128], Tmin=(929.224,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(198.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(2methyleneoxetane) + radical(CCOJ) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C[C]1OC1[O](24186)',
    structure = SMILES('[CH2]C(=C)C[C]1OC1[O]'),
    E0 = (294.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.802066,0.0708621,-6.62378e-05,3.18575e-08,-4.37576e-12,35572.9,27.9983], Tmin=(100,'K'), Tmax=(784.662,'K')), NASAPolynomial(coeffs=[10.4346,0.0319795,-1.1448e-05,1.90858e-09,-1.23238e-13,33746.5,-18.1424], Tmin=(784.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1(C[C]=C)OC1[O](24189)',
    structure = SMILES('[CH2]C1(C[C]=C)OC1[O]'),
    E0 = (408.853,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4095.79,'J/mol'), sigma=(6.91196,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=639.75 K, Pc=28.14 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0402812,0.0882552,-0.000114345,8.13439e-08,-2.25853e-11,49315.1,29.4963], Tmin=(100,'K'), Tmax=(988.31,'K')), NASAPolynomial(coeffs=[12.8026,0.0277934,-9.20962e-06,1.40655e-09,-8.33133e-14,47222.7,-29.7479], Tmin=(988.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(408.853,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(Cds_S)"""),
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
    label = '[CH2]C(=C)[C]1OC1[O](16777)',
    structure = SMILES('[CH2]C([CH2])=C1OC1[O]'),
    E0 = (270.985,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00333,0.0539833,-2.5209e-05,-1.38175e-08,1.10686e-11,32710.6,21.6279], Tmin=(100,'K'), Tmax=(976.023,'K')), NASAPolynomial(coeffs=[17.6647,0.0131227,-4.55616e-06,8.61782e-10,-6.47348e-14,28152,-65.0406], Tmin=(976.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(CCOJ) + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=C)C1([CH2])[CH]O1(24996)',
    structure = SMILES('[CH2]C(=C)C1([CH2])[CH]O1'),
    E0 = (458.092,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,1000,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.796616,0.0837848,-9.09301e-05,4.87284e-08,-9.55899e-12,55287.1,28.4181], Tmin=(100,'K'), Tmax=(1510.3,'K')), NASAPolynomial(coeffs=[19.3343,0.0106556,1.37805e-06,-7.03671e-10,6.12792e-14,51466.1,-69.5211], Tmin=(1510.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.092,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCsJO) + radical(CJC(C)OC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1([C]=C)OC1[O](24636)',
    structure = SMILES('[CH2]C1([C]=C)OC1[O]'),
    E0 = (425.241,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788192,0.069004,-8.6847e-05,5.83844e-08,-1.49772e-11,51261.8,25.9874], Tmin=(100,'K'), Tmax=(1102.61,'K')), NASAPolynomial(coeffs=[12.082,0.0182244,-4.42279e-06,4.81044e-10,-1.92513e-14,49367.5,-26.8983], Tmin=(1102.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.241,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(Cds_S)"""),
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
    label = '[CH]C1(C([CH2])=C)OC1[O](24997)',
    structure = SMILES('[CH]C1(C([CH2])=C)OC1[O]'),
    E0 = (537.548,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.168626,0.0774098,-7.80757e-05,4.0818e-08,-8.10741e-12,64814.5,30.5459], Tmin=(100,'K'), Tmax=(1416.25,'K')), NASAPolynomial(coeffs=[17.3268,0.0173313,-3.14865e-06,2.30546e-10,-4.20002e-15,60928.5,-56.1719], Tmin=(1416.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.548,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(Allyl_P) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])C1([CH2])OC1[O](24197)',
    structure = SMILES('[CH]C(=C)C1([CH2])OC1[O]'),
    E0 = (519.029,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,1089.57,1089.57,1089.57,1089.57,1089.57,1089.57,1089.57,1089.57,1089.57,1089.57,1089.57,1089.57,1089.57,2236.87],'cm^-1')),
        HinderedRotor(inertia=(0.103177,'amu*angstrom^2'), symmetry=1, barrier=(2.37225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103177,'amu*angstrom^2'), symmetry=1, barrier=(2.37225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103177,'amu*angstrom^2'), symmetry=1, barrier=(2.37225,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0578905,0.0847325,-9.35282e-05,5.71003e-08,-1.36868e-11,62574.7,30.4165], Tmin=(100,'K'), Tmax=(1121.8,'K')), NASAPolynomial(coeffs=[14.2575,0.0272137,-7.96034e-06,1.10383e-09,-6.11147e-14,59770.3,-38.4766], Tmin=(1121.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(519.029,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)C12COC1O2(24998)',
    structure = SMILES('[CH2]C(=C)C12COC1O2'),
    E0 = (64.9794,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.949978,0.0505297,7.99968e-06,-4.83986e-08,2.20891e-11,7939.66,23.0076], Tmin=(100,'K'), Tmax=(999.97,'K')), NASAPolynomial(coeffs=[17.3941,0.0234191,-9.33698e-06,1.82984e-09,-1.36264e-13,2717.64,-65.9922], Tmin=(999.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(64.9794,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + polycyclic(s2_3_4_ane) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1CCC12OC2[O](24975)',
    structure = SMILES('C=C1CCC12OC2[O]'),
    E0 = (107.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35864,0.0488444,-1.12981e-05,-1.29106e-08,6.32455e-12,13086.2,24.4662], Tmin=(100,'K'), Tmax=(1135.46,'K')), NASAPolynomial(coeffs=[11.1155,0.0325893,-1.37567e-05,2.58437e-09,-1.80825e-13,9702.66,-28.9872], Tmin=(1135.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + polycyclic(s1_3_4_ane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C12OC1OCC2=C(24999)',
    structure = SMILES('[CH2]C12OC1OCC2=C'),
    E0 = (52.2414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953862,0.0547411,-1.41863e-05,-1.92526e-08,1.04694e-11,6403.4,23.0395], Tmin=(100,'K'), Tmax=(1051.28,'K')), NASAPolynomial(coeffs=[14.9993,0.0270058,-1.12909e-05,2.17047e-09,-1.56222e-13,2029.75,-52.1802], Tmin=(1051.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(52.2414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + polycyclic(s2_3_5_ane) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C](O[CH][O])C([CH2])=C(25000)',
    structure = SMILES('[CH2]C([CH2])=C([CH2])O[CH][O]'),
    E0 = (421.485,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,3025,407.5,1350,352.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0483153,0.0881937,-9.91058e-05,5.83196e-08,-1.37576e-11,50834.3,29.0396], Tmin=(100,'K'), Tmax=(1027.87,'K')), NASAPolynomial(coeffs=[15.212,0.0291835,-1.29903e-05,2.46579e-09,-1.72761e-13,47717,-44.5331], Tmin=(1027.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.485,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(OCOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(C=C(O)CJ) + radical(OCJO)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])([O])[CH][O](25001)',
    structure = SMILES('[CH2]C(=C)C([CH2])([O])[CH][O]'),
    E0 = (567.201,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0485096,0.0940914,-0.000137715,1.14908e-07,-3.8559e-11,68353.9,33.6046], Tmin=(100,'K'), Tmax=(815.343,'K')), NASAPolynomial(coeffs=[9.6972,0.0378146,-1.77322e-05,3.3537e-09,-2.30446e-13,67077.8,-9.15242], Tmin=(815.343,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.201,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + longDistanceInteraction_noncyclic(CdCs-ST) + group(Cds-CdsHH) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(C=CC(C)(O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C](C([CH2])=C)C([O])[O](25002)',
    structure = SMILES('[CH2]C([CH2])=C([CH2])C([O])[O]'),
    E0 = (437.837,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,1593.31,1593.33],'cm^-1')),
        HinderedRotor(inertia=(0.00498617,'amu*angstrom^2'), symmetry=1, barrier=(8.97694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89721,'amu*angstrom^2'), symmetry=1, barrier=(43.6207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.024228,'amu*angstrom^2'), symmetry=1, barrier=(43.6203,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.89727,'amu*angstrom^2'), symmetry=1, barrier=(43.622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28939,0.070793,-4.70678e-05,-3.845e-08,6.14524e-11,52745.9,28.1133], Tmin=(100,'K'), Tmax=(486.363,'K')), NASAPolynomial(coeffs=[6.74292,0.0442589,-2.17266e-05,4.25051e-09,-3.00183e-13,51998.8,3.50681], Tmin=(486.363,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.837,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsCs) + radical(CCOJ) + radical(CCOJ) + radical(Allyl_P) + radical(Allyl_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C1(C)OC1=O(25003)',
    structure = SMILES('[CH2]C(=C)C1(C)OC1=O'),
    E0 = (-32.1463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.651582,0.0562683,2.75169e-06,-5.2445e-08,2.65374e-11,-3729.92,23.1455], Tmin=(100,'K'), Tmax=(949.308,'K')), NASAPolynomial(coeffs=[19.791,0.0185482,-5.47371e-06,9.63995e-10,-7.18712e-14,-9297.97,-78.3831], Tmin=(949.308,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-32.1463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(cyclopropanone) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C1(C(=C)C)OC1=O(25004)',
    structure = SMILES('[CH2]C1(C(=C)C)OC1=O'),
    E0 = (26.8646,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.267806,0.0719192,-5.07447e-05,7.6561e-09,4.22606e-12,3374.53,24.6873], Tmin=(100,'K'), Tmax=(975.71,'K')), NASAPolynomial(coeffs=[17.8927,0.0217388,-7.53587e-06,1.32004e-09,-9.15432e-14,-1115.57,-65.2937], Tmin=(975.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(26.8646,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(cyclopropanone) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C([CH2])C1([CH2])O[C]1[O](25005)',
    structure = SMILES('[CH2]C([CH2])C1([CH2])O[C]1[O]'),
    E0 = (649.599,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,180,980.902,980.902,980.902,980.902,980.902,980.902,980.902,980.902,980.902,2352.93],'cm^-1')),
        HinderedRotor(inertia=(0.0748951,'amu*angstrom^2'), symmetry=1, barrier=(1.72198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0748951,'amu*angstrom^2'), symmetry=1, barrier=(1.72198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0748951,'amu*angstrom^2'), symmetry=1, barrier=(1.72198,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0748951,'amu*angstrom^2'), symmetry=1, barrier=(1.72198,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.245239,0.0942576,-0.000129257,9.32441e-08,-2.53205e-11,78280.8,33.5377], Tmin=(100,'K'), Tmax=(1072.73,'K')), NASAPolynomial(coeffs=[13.6517,0.0241574,-5.67288e-06,5.54445e-10,-1.67864e-14,76351.1,-29.5802], Tmin=(1072.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(Cs_P) + radical(CJC(C)OC) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C](C)C1([CH2])O[C]1[O](25006)',
    structure = SMILES('[CH2][C](C)C1([CH2])O[C]1[O]'),
    E0 = (597.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,980.52,980.52,980.52,980.52,980.52,980.52,980.52,980.52,980.52,2344.18],'cm^-1')),
        HinderedRotor(inertia=(0.0737556,'amu*angstrom^2'), symmetry=1, barrier=(1.69579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737556,'amu*angstrom^2'), symmetry=1, barrier=(1.69579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737556,'amu*angstrom^2'), symmetry=1, barrier=(1.69579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0737556,'amu*angstrom^2'), symmetry=1, barrier=(1.69579,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0217188,0.0857358,-0.000104611,6.97932e-08,-1.7996e-11,71958.2,31.6883], Tmin=(100,'K'), Tmax=(1080.66,'K')), NASAPolynomial(coeffs=[13.5068,0.0256477,-7.08438e-06,9.16427e-10,-4.66524e-14,69637.7,-31.6665], Tmin=(1080.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCJ(C)CO) + radical(Cs_P) + radical(CJC(C)OC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1([C]2CC2)OC1[O](25007)',
    structure = SMILES('[CH2]C1([C]2CC2)OC1[O]'),
    E0 = (349.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808848,0.0594113,-1.69917e-05,-2.97406e-08,1.9672e-11,42143.2,25.2651], Tmin=(100,'K'), Tmax=(886.181,'K')), NASAPolynomial(coeffs=[16.1522,0.0211844,-4.8084e-06,6.055e-10,-3.54582e-14,38205.4,-53.7779], Tmin=(886.181,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(349.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Cyclopropane) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCJ(C)CO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]1CCC12OC2[O](25008)',
    structure = SMILES('[CH2][C]1CCC12OC2[O]'),
    E0 = (350.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22436,0.0538127,-2.8388e-05,6.21153e-09,-3.11026e-13,42239.2,26.2767], Tmin=(100,'K'), Tmax=(1517.46,'K')), NASAPolynomial(coeffs=[12.6168,0.030711,-1.2401e-05,2.19684e-09,-1.45323e-13,37983.9,-36.0647], Tmin=(1517.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(350.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + polycyclic(s1_3_4_ane) + radical(CCOJ) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]1COC2OC12[CH2](25009)',
    structure = SMILES('[CH2][C]1COC2OC12[CH2]'),
    E0 = (289.05,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.771293,0.0606833,-3.09494e-05,-1.85645e-09,4.49829e-12,34889.7,24.6222], Tmin=(100,'K'), Tmax=(1078.73,'K')), NASAPolynomial(coeffs=[14.65,0.0279204,-1.13949e-05,2.12901e-09,-1.49742e-13,30807.4,-48.4294], Tmin=(1078.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + polycyclic(s2_3_5_ane) + radical(CCJ(C)CO) + radical(CJC(C)OC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1([CH2])CC12OC2[O](25010)',
    structure = SMILES('[CH2]C1([CH2])CC12OC2[O]'),
    E0 = (407.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16036,0.0643355,-5.01692e-05,2.05975e-08,-3.57286e-12,49086.6,24.8863], Tmin=(100,'K'), Tmax=(1308.75,'K')), NASAPolynomial(coeffs=[11.0278,0.0341778,-1.56052e-05,2.99127e-09,-2.09752e-13,46503.7,-25.3738], Tmin=(1308.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(407.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + polycyclic(s1_3_3_ane) + radical(CCOJ) + radical(Neopentyl) + radical(Neopentyl)"""),
)

species(
    label = '[CH2]C1([CH2])OC2OC21[CH2](25011)',
    structure = SMILES('[CH2]C1([CH2])OC2OC21[CH2]'),
    E0 = (409.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.378343,0.0985977,-0.00012533,8.09749e-08,-2.05964e-11,49379.5,23.7711], Tmin=(100,'K'), Tmax=(963.578,'K')), NASAPolynomial(coeffs=[17.2455,0.0254375,-1.14415e-05,2.17892e-09,-1.52707e-13,45983.1,-60.5998], Tmin=(963.578,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + polycyclic(s2_3_4_ane) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C(=C)C(=C)O[CH][O](25012)',
    structure = SMILES('[CH2]C([CH2])=C([CH2])OC=O'),
    E0 = (50.8435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3020,3040,3060,3080,3100,415,440,465,780,815,850,1435,1455,1475,900,1000,1100,2782.5,750,1395,475,1775,1000,1013.12,1297.1],'cm^-1')),
        HinderedRotor(inertia=(1.0608,'amu*angstrom^2'), symmetry=1, barrier=(24.3898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0204273,'amu*angstrom^2'), symmetry=1, barrier=(24.3897,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0608,'amu*angstrom^2'), symmetry=1, barrier=(24.3899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.06079,'amu*angstrom^2'), symmetry=1, barrier=(24.3896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.99441,'amu*angstrom^2'), symmetry=1, barrier=(91.8393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0275201,0.0788112,-7.2129e-05,3.29173e-08,-5.91475e-12,6268.21,29.4082], Tmin=(100,'K'), Tmax=(1349.17,'K')), NASAPolynomial(coeffs=[19.2107,0.0217735,-8.71421e-06,1.58183e-09,-1.08264e-13,1077.13,-69.1666], Tmin=(1349.17,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(50.8435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(336.736,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-OdOsH) + radical(Allyl_P) + radical(Allyl_P) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C(=C)C(=C)C([O])[O](23848)',
    structure = SMILES('[CH2]C(=C)C(=C)C([O])[O]'),
    E0 = (261.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,325,375,415,465,420,450,1700,1750,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,3000,3100,440,815,1455,1000,278.165,278.379,2182.91,2305.41],'cm^-1')),
        HinderedRotor(inertia=(0.278603,'amu*angstrom^2'), symmetry=1, barrier=(15.3462,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602335,'amu*angstrom^2'), symmetry=1, barrier=(33.1377,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.603958,'amu*angstrom^2'), symmetry=1, barrier=(33.1404,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.693916,0.0789506,-9.48161e-05,6.97161e-08,-2.19273e-11,31575.5,27.5582], Tmin=(100,'K'), Tmax=(760.563,'K')), NASAPolynomial(coeffs=[7.94076,0.0408379,-1.96502e-05,3.83062e-09,-2.70732e-13,30473.1,-5.42028], Tmin=(760.563,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(CCOJ) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C([CH2])([O])C=O(25013)',
    structure = SMILES('[CH2]C(=C)C([CH2])([O])C=O'),
    E0 = (265.709,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.183147,0.100269,-0.000155673,1.31978e-07,-4.39765e-11,32100.1,29.8889], Tmin=(100,'K'), Tmax=(846.83,'K')), NASAPolynomial(coeffs=[10.5635,0.0356189,-1.65573e-05,3.09261e-09,-2.09825e-13,30778,-17.2303], Tmin=(846.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.709,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=CC(C)(C=O)OJ) + radical(C=CC(O)(C=O)CJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C1([CH2])OC1=O(25014)',
    structure = SMILES('[CH2]C(=C)C1([CH2])OC1=O'),
    E0 = (178.364,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.2931,0.0691125,-4.10318e-05,-7.33363e-09,1.09978e-11,21597,24.1844], Tmin=(100,'K'), Tmax=(946.316,'K')), NASAPolynomial(coeffs=[20.0038,0.0159936,-4.69756e-06,7.89017e-10,-5.62081e-14,16514.4,-76.9645], Tmin=(946.316,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(cyclopropanone) + radical(CJC(C)OC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][C]=C(459)',
    structure = SMILES('[CH2][C]=C'),
    E0 = (395.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.110742,'amu*angstrom^2'), symmetry=1, barrier=(23.6192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28575,0.0103679,1.6796e-05,-2.70041e-08,1.05227e-11,47575.3,10.4038], Tmin=(100,'K'), Tmax=(980.899,'K')), NASAPolynomial(coeffs=[6.52804,0.0104286,-3.60844e-06,6.68464e-10,-4.85509e-14,46300.2,-8.43297], Tmin=(980.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]1OC1[O](19087)',
    structure = SMILES('[CH2][C]1OC1[O]'),
    E0 = (313.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,3000,3100,440,815,1455,1000,180,180,180,552.746,1513.57,1514.62,1514.74,1514.76,1515.06],'cm^-1')),
        HinderedRotor(inertia=(2.73416,'amu*angstrom^2'), symmetry=1, barrier=(62.8637,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.36931,0.0382955,-5.24741e-05,4.13379e-08,-1.22341e-11,37790.5,17.6825], Tmin=(100,'K'), Tmax=(1040.28,'K')), NASAPolynomial(coeffs=[5.37338,0.0158905,-4.51713e-06,5.74772e-10,-2.75886e-14,37752.7,5.89376], Tmin=(1040.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(CJCO)"""),
)

species(
    label = '[CH2]C(=C)C1([CH2])O[C]1[O](25015)',
    structure = SMILES('[CH2]C(=C)C1([CH2])O[C]1[O]'),
    E0 = (505.09,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650,180,1032.9,1032.9,1032.9,1032.9,1032.9,1032.9,1032.9,1032.9,1032.9,2295.79],'cm^-1')),
        HinderedRotor(inertia=(0.0753763,'amu*angstrom^2'), symmetry=1, barrier=(1.73305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0753763,'amu*angstrom^2'), symmetry=1, barrier=(1.73305,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0753763,'amu*angstrom^2'), symmetry=1, barrier=(1.73305,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.336355,0.0784471,-7.7655e-05,2.94488e-08,1.11297e-12,60882.6,28.8854], Tmin=(100,'K'), Tmax=(792.899,'K')), NASAPolynomial(coeffs=[16.114,0.0193169,-4.50744e-06,4.97675e-10,-2.21072e-14,57737.3,-47.6278], Tmin=(792.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(505.09,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(Cs_P) + radical(CJC(C)OC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C1([CH2])O[C]1O(25016)',
    structure = SMILES('[CH2]C(=C)C1([CH2])O[C]1O'),
    E0 = (279.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.839713,0.0912995,-0.000103409,5.83368e-08,-1.22406e-11,33789.6,32.2136], Tmin=(100,'K'), Tmax=(1367.3,'K')), NASAPolynomial(coeffs=[20.446,0.0132087,-3.84122e-07,-3.58334e-10,3.85639e-14,29447.6,-71.7287], Tmin=(1367.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(279.385,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cs_P) + radical(CJC(C)OC) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C1(C)O[C]1[O](25017)',
    structure = SMILES('[CH2]C(=C)C1(C)O[C]1[O]'),
    E0 = (294.58,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.130295,0.0756245,-7.05771e-05,3.49944e-08,-6.67035e-12,35591.5,30.7899], Tmin=(100,'K'), Tmax=(1465.53,'K')), NASAPolynomial(coeffs=[16.4895,0.0206737,-4.5195e-06,4.80443e-10,-2.11399e-14,31749.8,-52.2302], Tmin=(1465.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.58,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(Cs_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(C)C1([CH2])OC1[O](25018)',
    structure = SMILES('[CH]=C(C)C1([CH2])OC1[O]'),
    E0 = (395.44,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0698172,0.0875916,-0.000107956,7.16547e-08,-1.84157e-11,47708.7,29.7548], Tmin=(100,'K'), Tmax=(1062.9,'K')), NASAPolynomial(coeffs=[14.5226,0.0243513,-6.96106e-06,9.40995e-10,-5.03539e-14,45076.9,-39.3235], Tmin=(1062.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1(C(=C)C)O[C]1[O](25019)',
    structure = SMILES('[CH2]C1(C(=C)C)O[C]1[O]'),
    E0 = (353.591,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0424665,0.084785,-0.000101678,6.64151e-08,-1.68819e-11,42671.6,30.3299], Tmin=(100,'K'), Tmax=(1067.55,'K')), NASAPolynomial(coeffs=[14.0629,0.0249127,-7.23994e-06,1.0003e-09,-5.48742e-14,40096.3,-36.2683], Tmin=(1067.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(353.591,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(Cs_P) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH]=C([CH2])C1(C)OC1[O](25020)',
    structure = SMILES('[CH]C(=C)C1(C)OC1[O]'),
    E0 = (308.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,180,1085.12,1085.12,1085.12,1085.12,1085.12,1085.12,1085.12,1085.12,1085.12,1085.12,1085.12,1085.12,1085.12,2238.26],'cm^-1')),
        HinderedRotor(inertia=(0.0990447,'amu*angstrom^2'), symmetry=1, barrier=(2.27723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0990447,'amu*angstrom^2'), symmetry=1, barrier=(2.27723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0990447,'amu*angstrom^2'), symmetry=1, barrier=(2.27723,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0802536,0.0762572,-6.44768e-05,3.03588e-08,-5.68022e-12,37264.6,30.7523], Tmin=(100,'K'), Tmax=(1438.6,'K')), NASAPolynomial(coeffs=[14.8682,0.0283537,-7.91842e-06,1.08558e-09,-6.0771e-14,33619.7,-44.5213], Tmin=(1438.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])C1([CH2])OC1O(25021)',
    structure = SMILES('[CH]C(=C)C1([CH2])OC1O'),
    E0 = (293.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.796617,0.0920064,-9.75353e-05,5.39524e-08,-1.1341e-11,35463,32.2015], Tmin=(100,'K'), Tmax=(1339.78,'K')), NASAPolynomial(coeffs=[18.7211,0.0210339,-3.85547e-06,2.62129e-10,-2.23414e-15,31373,-63.4152], Tmin=(1339.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CJC(C)OC) + radical(AllylJ2_triplet)"""),
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
    E0 = (299.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (686.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (457.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (457.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (578.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (708.787,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (977.445,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (863.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (749.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (730.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (308.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (308.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (306.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (424.414,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (570.13,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (440.766,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (363.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (324.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (712.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (605.458,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (531.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (425.607,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (372.769,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (407.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (421.886,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (300.262,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (328.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (367.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (436.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (458.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (521.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (709.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (716.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (457.88,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (416.902,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (524.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (436.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (352.827,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (411.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['C=C1OC1[O](19061)', 'allene(458)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(S)(3112)', '[CH2]C(=C)C(=C)[CH][O](23578)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.96556e+06,'m^3/(mol*s)'), n=0, Ea=(0.08368,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [o_atom_singlet;mb_db]
Euclidian distance = 0
family: 1+2_Cycloaddition"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C(=C)[C]1CC([O])O1(24995)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C1(C[C]=C)OC1[O](24189)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', '[CH2]C(=C)[C]1OC1[O](16777)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/ODMustO;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(4)', '[CH2]C(=C)C1([CH2])[CH]O1(24996)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [C_rad/H/CsO;O_birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['CH2(T)(20)', '[CH2]C1([C]=C)OC1[O](24636)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', '[CH]C1(C([CH2])=C)OC1[O](24997)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', '[CH]=C([CH2])C1([CH2])OC1[O](24197)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C(=C)C12COC1O2(24998)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['C=C1CCC12OC2[O](24975)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for R4_SSS;C_rad_out_2H;Cpri_rad_out_2H
Exact match found for rate rule [R4_SSS;C_rad_out_2H;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C12OC1OCC2=C(24999)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.552e+10,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R5_SSSS;Y_rad_out;Cpri_rad_out_2H] for rate rule [R5_SSSS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C](O[CH][O])C([CH2])=C(25000)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(=C)C([CH2])([O])[CH][O](25001)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C](C([CH2])=C)C([O])[O](25002)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(4.36e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C(=C)C1(C)OC1=O(25003)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C1(C(=C)C)OC1=O(25004)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH2])C1([CH2])O[C]1[O](25005)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C](C)C1([CH2])O[C]1[O](25006)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.328e+09,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radEndo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C1([C]2CC2)OC1[O](25007)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2][C]1CCC12OC2[O](25008)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.02142e+08,'s^-1'), n=0.996667, Ea=(125.764,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2][C]1COC2OC12[CH2](25009)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.6996e+10,'s^-1'), n=0.261395, Ea=(72.9255,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SS_D;doublebond_intra;radadd_intra] for rate rule [R5_SS_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C1([CH2])CC12OC2[O](25010)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3.24059e+09,'s^-1'), n=0.736667, Ea=(107.451,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 105.6 to 107.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C1([CH2])OC2OC21[CH2](25011)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(2.724e+10,'s^-1','*|/',3), n=0.478, Ea=(122.043,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R5_SS_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C(=C)C(=C)O[CH][O](25012)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(7.1786e+12,'s^-1'), n=-0.0937843, Ea=(249.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H_secDe;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(=C)C(=C)C([O])[O](23848)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.43572e+13,'s^-1'), n=-0.0937843, Ea=(66.7967,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H_secDe;radadd_intra] for rate rule [R4_S_D;doublebond_intra_2H_secDe;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH2]C(=C)C([CH2])([O])C=O(25013)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.785e+11,'s^-1'), n=0.342, Ea=(102.29,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonylbond_intra_H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', '[CH2]C(=C)C1([CH2])OC1=O(25014)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.83701,'m^3/(mol*s)'), n=1.71338, Ea=(46.5052,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [CO_O;HJ] + [CO-NdNd_O;YJ] for rate rule [CO-NdNd_O;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C1OC1[O](19061)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(0.0018001,'m^3/(mol*s)'), n=2.5041, Ea=(20.5128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Cds;CJ] for rate rule [Cds-OsCs_Cds;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]1OC1[O](19087)', 'allene(458)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]1OC1[O](19087)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[CH2]C(=C)C1([CH2])O[C]1[O](25015)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.57643e+07,'m^3/(mol*s)'), n=-0.222444, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0398748526415, var=2.83447135176, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H
    Total Standard Deviation in ln(k): 3.47533765731
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C(=C)C1([CH2])O[C]1O(25016)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S;O_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C(=C)C1(C)O[C]1[O](25017)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.30996,'s^-1'), n=3.5644, Ea=(117.059,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;C_rad_out_2H;XH_out] for rate rule [R3H_SS_23cy3;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C(C)C1([CH2])OC1[O](25018)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C1(C(=C)C)O[C]1[O](25019)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(32400,'s^-1'), n=2.04, Ea=(82.4248,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C([CH2])C1(C)OC1[O](25020)'],
    products = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH]=C([CH2])C1([CH2])OC1O(25021)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.936e+06,'s^-1','*|/',3), n=1.554, Ea=(111.445,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""From training reaction 290 used for R5H_SSSD;O_rad_out;Cd_H_out_singleH
Exact match found for rate rule [R5H_SSSD;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #5701',
    isomers = [
        '[CH2]C(=C)C1([CH2])OC1[O](24187)',
    ],
    reactants = [
        ('C=C1OC1[O](19061)', 'allene(458)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #5701',
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


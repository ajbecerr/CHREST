species(
    label = '[CH2]C(=C)C[C]1OC1[O](24186)',
    structure = SMILES('[CH2]C(=C)C[C]1OC1[O]'),
    E0 = (294.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,997.405,997.405,997.405,997.405,997.405,997.405,997.405,997.405,997.405,2330.21],'cm^-1')),
        HinderedRotor(inertia=(0.0560707,'amu*angstrom^2'), symmetry=1, barrier=(1.28918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0560707,'amu*angstrom^2'), symmetry=1, barrier=(1.28918,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0560707,'amu*angstrom^2'), symmetry=1, barrier=(1.28918,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.802066,0.0708621,-6.62378e-05,3.18575e-08,-4.37576e-12,35572.9,27.9983], Tmin=(100,'K'), Tmax=(784.662,'K')), NASAPolynomial(coeffs=[10.4346,0.0319795,-1.1448e-05,1.90858e-09,-1.23238e-13,33746.5,-18.1424], Tmin=(784.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Allyl_P)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4117.5,'J/mol'), sigma=(6.93639,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=643.14 K, Pc=27.99 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.26811,0.0781585,-6.92136e-05,2.06602e-08,3.57149e-12,36201.2,27.7447], Tmin=(100,'K'), Tmax=(821.339,'K')), NASAPolynomial(coeffs=[16.2773,0.0216813,-5.31496e-06,6.49292e-10,-3.34487e-14,32846.6,-50.7517], Tmin=(821.339,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CJC(C)OC) + radical(Allyl_P)"""),
)

species(
    label = 'C=[C]CC[C]1OC1[O](24188)',
    structure = SMILES('C=[C]CC[C]1OC1[O]'),
    E0 = (396.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,1000,2950,3100,1380,975,1025,1650,1685,370,300,800,800,800,800,800,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4047.21,'J/mol'), sigma=(6.77452,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=632.17 K, Pc=29.54 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741903,0.0756115,-9.415e-05,7.22363e-08,-2.23768e-11,47793.3,30.5755], Tmin=(100,'K'), Tmax=(913.688,'K')), NASAPolynomial(coeffs=[7.55598,0.0362871,-1.40062e-05,2.38826e-09,-1.53773e-13,46944.3,0.485113], Tmin=(913.688,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Cds_S)"""),
)

species(
    label = '[O]C1[C]O1(24298)',
    structure = SMILES('[O]C1[C]O1'),
    E0 = (419.455,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,180,208.172,904.773,906.483,906.539,908.364,908.721],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.8532,0.0217083,-2.53977e-05,1.58489e-08,-3.56375e-12,50493.2,11.7167], Tmin=(100,'K'), Tmax=(1382.83,'K')), NASAPolynomial(coeffs=[6.10444,0.00516132,2.98991e-07,-2.7478e-10,2.65047e-14,50276.9,-2.5533], Tmin=(1382.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.455,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C([CH2])=C(15608)',
    structure = SMILES('[CH2]C([CH2])=C'),
    E0 = (269.033,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.0902703,'amu*angstrom^2'), symmetry=1, barrier=(41.8152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0469548,'amu*angstrom^2'), symmetry=1, barrier=(21.7476,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.55516,0.0212335,2.47955e-05,-4.79586e-08,2.00948e-11,32418.8,14.0196], Tmin=(100,'K'), Tmax=(958.92,'K')), NASAPolynomial(coeffs=[10.7099,0.0134902,-4.1901e-06,7.65438e-10,-5.71862e-14,29646.9,-31.2788], Tmin=(958.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.033,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=C)C[C]1[CH]O1(25022)',
    structure = SMILES('[CH2]C(=C)C[C]1[CH]O1'),
    E0 = (453.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.13791,0.0717128,-7.09067e-05,3.7118e-08,-7.3978e-12,54641.2,27.237], Tmin=(100,'K'), Tmax=(1422.05,'K')), NASAPolynomial(coeffs=[15.4806,0.0180835,-3.29088e-06,2.40454e-10,-4.17995e-15,51336.6,-48.4614], Tmin=(1422.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(CCsJO) + radical(Allyl_P)"""),
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
    label = 'C=[C]C[C]1OC1[O](16693)',
    structure = SMILES('C=[C]C[C]1OC1[O]'),
    E0 = (420.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,2950,3100,1380,975,1025,1650,1685,370,180,1030.05,1030.05,1030.05,1030.05,1030.05,1030.05,1030.05,1030.05,1030.05,2320.65],'cm^-1')),
        HinderedRotor(inertia=(0.037605,'amu*angstrom^2'), symmetry=1, barrier=(0.864612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.037605,'amu*angstrom^2'), symmetry=1, barrier=(0.864612,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42506,0.0603926,-7.87374e-05,6.20241e-08,-1.92465e-11,50629,25.878], Tmin=(100,'K'), Tmax=(947.706,'K')), NASAPolynomial(coeffs=[6.4518,0.028188,-1.03732e-05,1.69892e-09,-1.05711e-13,50169.7,4.50024], Tmin=(947.706,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(420.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Cds_S)"""),
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
    label = '[CH]=C([CH2])C[C]1OC1[O](24196)',
    structure = SMILES('[CH]C(=C)C[C]1OC1[O]'),
    E0 = (514,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,180,1115.54,1115.54,1115.54,1115.54,1115.54,1115.54,1115.54,1115.54,1115.54,1115.54,1115.54,1115.54,2249.74],'cm^-1')),
        HinderedRotor(inertia=(0.0873519,'amu*angstrom^2'), symmetry=1, barrier=(2.00839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0873519,'amu*angstrom^2'), symmetry=1, barrier=(2.00839,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0873519,'amu*angstrom^2'), symmetry=1, barrier=(2.00839,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.65602,0.0752082,-8.22207e-05,5.66162e-08,-1.62192e-11,61938.6,30.0312], Tmin=(100,'K'), Tmax=(924.883,'K')), NASAPolynomial(coeffs=[8.28363,0.0377568,-1.42431e-05,2.39985e-09,-1.54034e-13,60718.6,-5.14012], Tmin=(924.883,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)CC12OC1O2(25023)',
    structure = SMILES('[CH2]C(=C)CC12OC1O2'),
    E0 = (66.5779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.281968,0.0557561,2.546e-05,-9.00459e-08,4.25831e-11,8165.08,23.7381], Tmin=(100,'K'), Tmax=(950.728,'K')), NASAPolynomial(coeffs=[27.1373,0.00763554,-9.62632e-07,2.47681e-10,-3.13942e-14,127.013,-119.884], Tmin=(950.728,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(66.5779,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsOsOs) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + polycyclic(s2_3_3_ane) + radical(Allyl_P)"""),
)

species(
    label = 'C=C1CC2(C1)OC2[O](24953)',
    structure = SMILES('C=C1CC2(C1)OC2[O]'),
    E0 = (115.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19014,0.053948,-2.64572e-05,2.95219e-09,8.42013e-13,13980.1,23.6603], Tmin=(100,'K'), Tmax=(1356.7,'K')), NASAPolynomial(coeffs=[12.687,0.0310481,-1.32967e-05,2.45962e-09,-1.6812e-13,9848.52,-39.0422], Tmin=(1356.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(357.522,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + polycyclic(s1_3_4_ane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=C)CC1=C(O)O1(25024)',
    structure = SMILES('[CH2]C(=C)CC1=C(O)O1'),
    E0 = (151.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.289832,0.0832905,-8.41639e-05,4.2206e-08,-8.20399e-12,18428,28.1656], Tmin=(100,'K'), Tmax=(1265.96,'K')), NASAPolynomial(coeffs=[20.6131,0.0172447,-5.9086e-06,9.96213e-10,-6.59947e-14,13135.5,-77.6087], Tmin=(1265.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)CC1OC1=O(25025)',
    structure = SMILES('[CH2]C(=C)CC1OC1=O'),
    E0 = (-20.884,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26127,0.0357592,6.29265e-05,-1.14203e-07,4.78899e-11,-2390.89,26.6405], Tmin=(100,'K'), Tmax=(952.698,'K')), NASAPolynomial(coeffs=[20.231,0.0175025,-4.98493e-06,9.56336e-10,-7.80434e-14,-8791.33,-78.5794], Tmin=(952.698,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-20.884,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(cyclopropanone) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(C)C=C1OC1[O](25026)',
    structure = SMILES('C=C(C)C=C1OC1[O]'),
    E0 = (57.2487,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.348671,0.0679596,-3.99376e-05,-4.09664e-09,8.27042e-12,7027.87,24.0896], Tmin=(100,'K'), Tmax=(987.787,'K')), NASAPolynomial(coeffs=[19.0583,0.0191964,-6.89046e-06,1.27211e-09,-9.20323e-14,2014.38,-72.6116], Tmin=(987.787,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.2487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCOJ)"""),
)

species(
    label = '[CH2]C(=C)C=C1OC1O(25027)',
    structure = SMILES('[CH2]C(=C)C=C1OC1O'),
    E0 = (-16.9573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0434084,0.0675426,-1.67881e-05,-4.54392e-08,2.75635e-11,-1879.06,23.9119], Tmin=(100,'K'), Tmax=(931.538,'K')), NASAPolynomial(coeffs=[25.6306,0.00742202,-8.94866e-08,-5.80682e-11,-1.87724e-15,-8804.73,-109.303], Tmin=(931.538,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.9573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(C)CC1=C([O])O1(24195)',
    structure = SMILES('C=C(C)C[C]1OC1=O'),
    E0 = (7.72153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25586,0.0443703,1.96261e-05,-5.74298e-08,2.49524e-11,1041.73,28.187], Tmin=(100,'K'), Tmax=(984.164,'K')), NASAPolynomial(coeffs=[15.4045,0.0251516,-9.43608e-06,1.78541e-09,-1.30771e-13,-3597.38,-49.2663], Tmin=(984.164,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(7.72153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(cyclopropanone) + radical(C2CsJO)"""),
)

species(
    label = '[CH2]C([CH2])[CH][C]1OC1[O](25028)',
    structure = SMILES('[CH2]C([CH2])[CH][C]1OC1[O]'),
    E0 = (635.18,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,1000,3025,407.5,1350,352.5,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,901.619,901.619,901.619,901.619,901.619,901.619,901.619,901.619,901.619,2454.99],'cm^-1')),
        HinderedRotor(inertia=(0.0443136,'amu*angstrom^2'), symmetry=1, barrier=(1.01886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0443136,'amu*angstrom^2'), symmetry=1, barrier=(1.01886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0443136,'amu*angstrom^2'), symmetry=1, barrier=(1.01886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0443136,'amu*angstrom^2'), symmetry=1, barrier=(1.01886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.56901,0.0761869,-9.17484e-05,6.45472e-08,-1.76668e-11,76517.4,35.3282], Tmin=(100,'K'), Tmax=(1055.57,'K')), NASAPolynomial(coeffs=[9.49653,0.0301126,-8.8759e-06,1.21849e-09,-6.54665e-14,75314.8,-4.99363], Tmin=(1055.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(635.18,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(CCJCO) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[O]C1O[C]1C[C]1CC1(25029)',
    structure = SMILES('[O]C1O[C]1C[C]1CC1'),
    E0 = (369.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25931,0.0616959,-5.8477e-05,3.81273e-08,-1.08943e-11,44572.2,28.18], Tmin=(100,'K'), Tmax=(892.325,'K')), NASAPolynomial(coeffs=[6.1954,0.0374736,-1.4237e-05,2.44345e-09,-1.5957e-13,43774.7,5.39601], Tmin=(892.325,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(369.785,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + ring(Cyclopropane) + ring(Ethylene_oxide) + radical(CCOJ) + radical(Tertalkyl) + radical(C2CsJO)"""),
)

species(
    label = '[CH2][C]1CC2(C1)OC2[O](25030)',
    structure = SMILES('[CH2][C]1CC2(C1)OC2[O]'),
    E0 = (383.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36334,0.0644933,-7.70151e-05,6.92095e-08,-2.65891e-11,46172.7,26.6905], Tmin=(100,'K'), Tmax=(799.259,'K')), NASAPolynomial(coeffs=[2.38673,0.0472224,-2.18013e-05,4.13705e-09,-2.86337e-13,46397.1,24.4101], Tmin=(799.259,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(383.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + polycyclic(s1_3_4_ane) + radical(CCOJ) + radical(Tertalkyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2][C]1COC2O[C]2C1(25031)',
    structure = SMILES('[CH2][C]1COC2O[C]2C1'),
    E0 = (275.859,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.54944,0.0459378,-7.32242e-06,-1.26886e-08,5.19719e-12,33272.8,24.206], Tmin=(100,'K'), Tmax=(1239.51,'K')), NASAPolynomial(coeffs=[9.99929,0.0357761,-1.57269e-05,2.96617e-09,-2.06025e-13,29863.9,-23.6749], Tmin=(1239.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.859,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(353.365,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + polycyclic(s2_3_6_ane) + radical(CCJ(C)CO) + radical(C2CsJO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(=C)CC([O])=C[O](25032)',
    structure = SMILES('[CH2]C(=C)CC([O])=C[O]'),
    E0 = (57.1101,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,405.431,405.431,405.431,405.431],'cm^-1')),
        HinderedRotor(inertia=(0.169559,'amu*angstrom^2'), symmetry=1, barrier=(19.778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16956,'amu*angstrom^2'), symmetry=1, barrier=(19.778,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16956,'amu*angstrom^2'), symmetry=1, barrier=(19.778,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.149598,0.0689136,-2.86165e-05,-2.69776e-08,1.96107e-11,7021.91,29.8158], Tmin=(100,'K'), Tmax=(933.174,'K')), NASAPolynomial(coeffs=[22.5968,0.0127094,-2.59288e-06,3.81753e-10,-2.93546e-14,1090.22,-86.2616], Tmin=(933.174,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(57.1101,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P)"""),
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
    label = '[CH2]C1([CH2])C[C]2OC2O1(25033)',
    structure = SMILES('[CH2]C1([CH2])C[C]2OC2O1'),
    E0 = (318.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688675,0.0778806,-8.21077e-05,4.80112e-08,-1.17501e-11,38410.5,23.0632], Tmin=(100,'K'), Tmax=(969.545,'K')), NASAPolynomial(coeffs=[10.8588,0.0359223,-1.71931e-05,3.37536e-09,-2.40634e-13,36438.4,-25.6871], Tmin=(969.545,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(318.407,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(349.208,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + polycyclic(s2_3_5_ane) + radical(C2CsJO) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C(=C)C=C1OC1[O](25034)',
    structure = SMILES('[CH2]C(=C)C=C1OC1[O]'),
    E0 = (208.748,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.370126,0.0652074,-3.04635e-05,-1.8697e-08,1.48361e-11,25250.5,23.5999], Tmin=(100,'K'), Tmax=(960.296,'K')), NASAPolynomial(coeffs=[21.1481,0.0134848,-4.07045e-06,7.45223e-10,-5.70292e-14,19654.1,-84.161], Tmin=(960.296,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.748,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCOJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)CC1=C([O])O1(24213)',
    structure = SMILES('[CH2]C(=C)C[C]1OC1=O'),
    E0 = (159.221,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26877,0.0417223,2.87204e-05,-7.15169e-08,3.12903e-11,19264.7,27.7278], Tmin=(100,'K'), Tmax=(964.968,'K')), NASAPolynomial(coeffs=[17.5229,0.0193909,-6.58754e-06,1.25175e-09,-9.52036e-14,14030.6,-60.9762], Tmin=(964.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(159.221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(cyclopropanone) + radical(C2CsJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)CC1=CO1(22338)',
    structure = SMILES('[CH2]C(=C)CC1=CO1'),
    E0 = (308.994,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,1000,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (95.1192,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.871648,0.05444,-1.02819e-05,-3.34087e-08,1.8797e-11,37289,24.2052], Tmin=(100,'K'), Tmax=(957.54,'K')), NASAPolynomial(coeffs=[18.3115,0.0163632,-5.11135e-06,9.20115e-10,-6.83538e-14,32354.9,-67.5003], Tmin=(957.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.994,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(Allyl_P)"""),
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
    label = '[CH2]C(=C)[CH][C]1OC1[O](25035)',
    structure = SMILES('[CH2]C(=C)[CH][C]1OC1[O]'),
    E0 = (411.731,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,967.255,967.255,967.255,967.255,967.255,967.255,967.255,967.255,967.255,2335.67],'cm^-1')),
        HinderedRotor(inertia=(0.0354808,'amu*angstrom^2'), symmetry=1, barrier=(0.815774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0354808,'amu*angstrom^2'), symmetry=1, barrier=(0.815774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0354808,'amu*angstrom^2'), symmetry=1, barrier=(0.815774,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.518083,0.0756976,-7.7639e-05,4.14427e-08,-7.77254e-12,49646.1,26.9033], Tmin=(100,'K'), Tmax=(838.914,'K')), NASAPolynomial(coeffs=[12.6412,0.0278165,-9.76864e-06,1.60695e-09,-1.02994e-13,47262.9,-31.5355], Tmin=(838.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.731,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C[C]1O[C]1[O](25036)',
    structure = SMILES('[CH2]C(=C)C[C]1O[C]1[O]'),
    E0 = (500.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.656025,0.0739771,-8.61712e-05,5.83461e-08,-1.59086e-11,60263.3,29.8896], Tmin=(100,'K'), Tmax=(957.451,'K')), NASAPolynomial(coeffs=[10.5503,0.0291047,-1.03308e-05,1.68109e-09,-1.05432e-13,58530.8,-16.5677], Tmin=(957.451,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(500.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Cs_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)[CH]C1OC1[O](25037)',
    structure = SMILES('[CH2]C(=C)[CH]C1OC1[O]'),
    E0 = (231.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.123566,0.0772535,-6.99726e-05,3.36461e-08,-6.33448e-12,28018.1,28.0888], Tmin=(100,'K'), Tmax=(1413,'K')), NASAPolynomial(coeffs=[17.115,0.0229303,-6.44116e-06,9.05005e-10,-5.21829e-14,23697.9,-59.0857], Tmin=(1413,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.626,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)C[C]1O[C]1O(25038)',
    structure = SMILES('[CH2]C(=C)C[C]1O[C]1O'),
    E0 = (274.355,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.439957,0.0748406,-6.64576e-05,2.22421e-08,1.62718e-12,33129.1,29.8147], Tmin=(100,'K'), Tmax=(829.477,'K')), NASAPolynomial(coeffs=[15.0095,0.0231045,-6.39551e-06,8.90274e-10,-5.14209e-14,30074.8,-41.5918], Tmin=(829.477,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.355,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(Cs_P) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C(=C)CC1O[C]1[O](25039)',
    structure = SMILES('[CH2]C(=C)CC1O[C]1[O]'),
    E0 = (319.956,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,995.475,995.475,995.475,995.475,995.475,995.475,995.475,995.475,995.475,2320.85],'cm^-1')),
        HinderedRotor(inertia=(0.053207,'amu*angstrom^2'), symmetry=1, barrier=(1.22333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.053207,'amu*angstrom^2'), symmetry=1, barrier=(1.22333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.053207,'amu*angstrom^2'), symmetry=1, barrier=(1.22333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.331732,0.0717512,-6.50302e-05,3.2582e-08,-6.5165e-12,38621.4,29.9385], Tmin=(100,'K'), Tmax=(1280.67,'K')), NASAPolynomial(coeffs=[14.3268,0.0254238,-7.70489e-06,1.14581e-09,-6.84956e-14,35251.3,-40.2042], Tmin=(1280.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(Cs_P) + radical(Allyl_P)"""),
)

species(
    label = 'C=C(C)[CH][C]1OC1[O](25040)',
    structure = SMILES('C=C(C)[CH][C]1OC1[O]'),
    E0 = (260.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700658,0.0756967,-7.55364e-05,3.74172e-08,-4.27232e-12,31414.8,25.9883], Tmin=(100,'K'), Tmax=(711.332,'K')), NASAPolynomial(coeffs=[10.1978,0.0341577,-1.29641e-05,2.22474e-09,-1.45679e-13,29763.5,-18.7054], Tmin=(711.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(C=CCJCO)"""),
)

species(
    label = '[CH]=C(C)C[C]1OC1[O](25041)',
    structure = SMILES('[CH]=C(C)C[C]1OC1[O]'),
    E0 = (390.411,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2950,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3120,650,792.5,1650,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.571809,0.07894,-9.97966e-05,7.53951e-08,-2.28142e-11,47075.7,29.6275], Tmin=(100,'K'), Tmax=(925.237,'K')), NASAPolynomial(coeffs=[8.69671,0.0346366,-1.30926e-05,2.20092e-09,-1.40257e-13,45965,-6.81625], Tmin=(925.237,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.411,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C(=C)[CH][C]1OC1O(25042)',
    structure = SMILES('[CH2]C(=C)[CH][C]1OC1O'),
    E0 = (186.026,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.368472,0.0848272,-8.86831e-05,4.8263e-08,-1.00973e-11,22540.7,29.2112], Tmin=(100,'K'), Tmax=(1299.39,'K')), NASAPolynomial(coeffs=[18.2075,0.019842,-4.65882e-06,5.32686e-10,-2.51085e-14,18371.8,-62.7377], Tmin=(1299.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.026,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C([CH2])CC1OC1[O](25043)',
    structure = SMILES('[CH]C(=C)CC1OC1[O]'),
    E0 = (333.895,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,3150,900,1100,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,180,1062,1062,1062,1062,1062,1062,1062,1062,1062,1062,1062,1062,1062,2251.36],'cm^-1')),
        HinderedRotor(inertia=(0.0791735,'amu*angstrom^2'), symmetry=1, barrier=(1.82035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0791735,'amu*angstrom^2'), symmetry=1, barrier=(1.82035,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0791735,'amu*angstrom^2'), symmetry=1, barrier=(1.82035,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.401073,0.0721802,-5.83278e-05,2.73187e-08,-5.31919e-12,40293.7,29.8301], Tmin=(100,'K'), Tmax=(1223.18,'K')), NASAPolynomial(coeffs=[12.0446,0.0341033,-1.16328e-05,1.86819e-09,-1.17389e-13,37445.3,-28.6884], Tmin=(1223.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(333.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C(C)C[C]1O[C]1[O](25044)',
    structure = SMILES('C=C(C)C[C]1O[C]1[O]'),
    E0 = (348.561,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.69203,0.0760396,-9.31922e-05,6.97399e-08,-2.11092e-11,42038.2,30.1741], Tmin=(100,'K'), Tmax=(919.134,'K')), NASAPolynomial(coeffs=[8.19701,0.0352657,-1.34105e-05,2.26942e-09,-1.45539e-13,41001.3,-3.53609], Tmin=(919.134,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(345.051,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Cs_P)"""),
)

species(
    label = '[CH]=C([CH2])C[C]1OC1O(25045)',
    structure = SMILES('[CH]C(=C)C[C]1OC1O'),
    E0 = (288.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (111.119,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.435572,0.076117,-6.26289e-05,2.05993e-08,1.32184e-12,34804.6,29.9724], Tmin=(100,'K'), Tmax=(818.186,'K')), NASAPolynomial(coeffs=[12.749,0.0317469,-1.03026e-05,1.60787e-09,-9.99313e-14,32259.8,-30.1994], Tmin=(818.186,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(288.294,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(340.893,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(AllylJ2_triplet)"""),
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
    E0 = (294.814,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (457.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (566.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (743.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (972.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (858.014,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (725.792,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (302.675,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (303.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (317.676,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (317.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (358.214,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (329.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (329.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (658.041,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (526.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (419.079,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (358.031,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (294.943,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (407.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (389.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (426.719,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (392.581,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (551.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (452.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (521.698,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (709.048,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (623.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (711.853,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (416.444,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (408.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (477.992,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (497.39,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (519.726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (373.971,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (378.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (390.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (423.535,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['C=C1OC1[O](19061)', 'allene(458)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(=C)C1([CH2])OC1[O](24187)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-R!HR!H)CJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C=[C]CC[C]1OC1[O](24188)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]C1[C]O1(24298)', '[CH2]C([CH2])=C(15608)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(4.02491e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH2]C(=C)C[C]1[CH]O1(25022)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['CH2(T)(20)', 'C=[C]C[C]1OC1[O](16693)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', '[CH]=C([CH2])C[C]1OC1[O](24196)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2]C(=C)CC12OC1O2(25023)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.61967e+11,'s^-1'), n=0.0247333, Ea=(7.86034,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;Y_rad_out;Cpri_rad_out_single] for rate rule [R3_SS;O_rad;Cpri_rad_out_NDMustO]
Euclidian distance = 3.7416573867739413
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['C=C1CC2(C1)OC2[O](24953)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Cpri_rad_out_2H] for rate rule [R4_SSS;C_rad_out_NDMustO;Cpri_rad_out_2H]
Euclidian distance = 3.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2]C(=C)CC1=C(O)O1(25024)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2]C(=C)CC1OC1=O(25025)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['C=C(C)C=C1OC1[O](25026)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2]C(=C)C=C1OC1O(25027)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['C=C(C)CC1=C([O])O1(24195)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.42e+09,'s^-1'), n=0.137, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad_NDe] for rate rule [R5radEndo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C([CH2])[CH][C]1OC1[O](25028)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[O]C1O[C]1C[C]1CC1(25029)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2][C]1CC2(C1)OC2[O](25030)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4.04e+06,'s^-1'), n=1.58, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_csNdNd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2][C]1COC2O[C]2C1(25031)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(4.9903e+10,'s^-1'), n=0.243684, Ea=(63.2167,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(=C)CC([O])=C[O](25032)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(237.833,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2]C1([CH2])CC12OC2[O](25010)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.47e+08,'s^-1'), n=1.11, Ea=(112.48,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_csNdNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 111.8 to 112.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2]C1([CH2])C[C]2OC2O1(25033)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(2.79781e+09,'s^-1'), n=0.611011, Ea=(94.2067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H;radadd_intra_O] + [R6;doublebond_intra_2H;radadd_intra] for rate rule [R6;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2]C(=C)C=C1OC1[O](25034)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(21.7644,'m^3/(mol*s)'), n=1.8572, Ea=(6.17853,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-CdH_Cds;HJ] for rate rule [Cds-CdH_Cds-OsCs;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C(=C)CC1=C([O])O1(24213)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(4)', '[CH2]C(=C)CC1=CO1(22338)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction25',
    reactants = ['C=C1OC1[O](19061)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00431078,'m^3/(mol*s)'), n=2.47625, Ea=(14.414,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Cds;CJ] for rate rule [Cds-HH_Cds-OsCs;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]1OC1[O](19087)', 'allene(458)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00086947,'m^3/(mol*s)'), n=2.67356, Ea=(32.0272,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ca_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]1OC1[O](19087)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2]C(=C)[CH][C]1OC1[O](25035)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.57643e+07,'m^3/(mol*s)'), n=-0.222444, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0398748526415, var=2.83447135176, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H
    Total Standard Deviation in ln(k): 3.47533765731
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', '[CH2]C(=C)C[C]1O[C]1[O](25036)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.57643e+07,'m^3/(mol*s)'), n=-0.222444, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0398748526415, var=2.83447135176, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H
    Total Standard Deviation in ln(k): 3.47533765731
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2]C(=C)[CH]C1OC1[O](25037)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(20108.5,'s^-1'), n=2.606, Ea=(121.63,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;C_rad_out_NonDe;Cs_H_out_H/Cd] for rate rule [R2H_S;C_rad_out_NDMustO;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2]C(=C)C[C]1O[C]1O(25038)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C(=C)CC1O[C]1[O](25039)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.43381e+07,'s^-1'), n=1.70481, Ea=(158.036,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_NonDe] for rate rule [R2H_S_cy3;Y_rad_out;Cs_H_out_NDMustO]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['C=C(C)[CH][C]1OC1[O](25040)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(2.19788e+08,'s^-1'), n=1.58167, Ea=(202.575,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_2Cd;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C(C)C[C]1OC1[O](25041)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH2]C(=C)[CH][C]1OC1O(25042)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(608922,'s^-1'), n=1.98643, Ea=(79.1566,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;Cs_H_out_H/Cd] for rate rule [R4HJ_2;O_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=C([CH2])CC1OC1[O](25043)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;Cs_H_out_NDMustO]
Euclidian distance = 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C(C)C[C]1O[C]1[O](25044)'],
    products = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(505536,'s^-1'), n=1.7378, Ea=(41.5716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_2H] for rate rule [R5HJ_1;Y_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2]C(=C)C[C]1OC1[O](24186)'],
    products = ['[CH]=C([CH2])C[C]1OC1O(25045)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(6.854,'s^-1'), n=3.311, Ea=(128.721,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;O_rad_out;Cd_H_out_singleH] for rate rule [R6HJ_2;O_rad_out;Cd_H_out_singleH]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #5700',
    isomers = [
        '[CH2]C(=C)C[C]1OC1[O](24186)',
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
    label = 'PDepNetwork #5700',
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


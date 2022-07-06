species(
    label = 'C[CH][CH]CC1[CH]CC1(5870)',
    structure = SMILES('C[CH][CH]CC1[CH]CC1'),
    E0 = (486.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,319.028,814.89,814.89,814.89,814.89,814.89,814.89,814.89,814.89,1686.92,1686.92,1686.92,1686.92,1686.92,1686.92,1686.92,1686.92],'cm^-1')),
        HinderedRotor(inertia=(0.132309,'amu*angstrom^2'), symmetry=1, barrier=(3.44019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132309,'amu*angstrom^2'), symmetry=1, barrier=(3.44019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132309,'amu*angstrom^2'), symmetry=1, barrier=(3.44019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132309,'amu*angstrom^2'), symmetry=1, barrier=(3.44019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07303,0.0579701,-1.08142e-05,-9.6676e-09,3.68507e-12,58617.2,36.055], Tmin=(100,'K'), Tmax=(1398.5,'K')), NASAPolynomial(coeffs=[9.40129,0.0508865,-2.1168e-05,3.82557e-09,-2.56768e-13,54651.1,-12.769], Tmin=(1398.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C1=CCC1(4846)',
    structure = SMILES('C1=CCC1'),
    E0 = (144.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,204.181,1057.98,1057.98,1057.98,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,1057.99,3484.69],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3170.39,'J/mol'), sigma=(5.50295,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=495.21 K, Pc=43.17 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.39614,-0.00343641,9.09434e-05,-1.1327e-07,4.24092e-11,17412.4,10.7746], Tmin=(100,'K'), Tmax=(950.248,'K')), NASAPolynomial(coeffs=[9.42765,0.013964,-4.06869e-06,7.74979e-10,-6.20376e-14,14334.2,-28.1813], Tmin=(950.248,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene)"""),
)

species(
    label = 'm1_allyl(186)',
    structure = SMILES('C=C[CH]C'),
    E0 = (121.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0800698,'amu*angstrom^2'), symmetry=1, barrier=(1.84096,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.062381,'amu*angstrom^2'), symmetry=1, barrier=(19.3985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2897.21,'J/mol'), sigma=(5.21305,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=452.54 K, Pc=46.4 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68514,0.0207746,2.19954e-05,-3.85566e-08,1.4972e-11,14712.8,14.0724], Tmin=(100,'K'), Tmax=(995.014,'K')), NASAPolynomial(coeffs=[7.60275,0.0207981,-7.87776e-06,1.45013e-09,-1.02662e-13,12754.4,-14.5511], Tmin=(995.014,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), label="""m1_allyl""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'CC=CCC1[CH]CC1(6235)',
    structure = SMILES('CC=CCC1[CH]CC1'),
    E0 = (213.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.93471,0.0506414,3.22369e-05,-6.65933e-08,2.57477e-11,25753.2,31.2587], Tmin=(100,'K'), Tmax=(1032.66,'K')), NASAPolynomial(coeffs=[12.4278,0.0458815,-1.86007e-05,3.50981e-09,-2.49832e-13,21259.6,-34.8224], Tmin=(1032.66,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(494.711,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = 'C[CH][CH]C[CH]C1CC1(6670)',
    structure = SMILES('C[CH][CH]C[CH]C1CC1'),
    E0 = (497.577,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,180,1309.16,1309.16,1309.16,1309.16,1309.16,1309.16,1309.16,1309.16,1309.16,1309.16,1309.16,1309.16,2122.68],'cm^-1')),
        HinderedRotor(inertia=(0.195492,'amu*angstrom^2'), symmetry=1, barrier=(4.49474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195492,'amu*angstrom^2'), symmetry=1, barrier=(4.49474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195492,'amu*angstrom^2'), symmetry=1, barrier=(4.49474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195492,'amu*angstrom^2'), symmetry=1, barrier=(4.49474,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195492,'amu*angstrom^2'), symmetry=1, barrier=(4.49474,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.889184,0.0641812,-2.64213e-05,2.80196e-09,3.89593e-13,59960.1,36.9173], Tmin=(100,'K'), Tmax=(1705.31,'K')), NASAPolynomial(coeffs=[14.6061,0.0429591,-1.73881e-05,3.0368e-09,-1.96972e-13,53689.3,-41.2496], Tmin=(1705.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.577,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([CH]C)C1[CH]CC1(5828)',
    structure = SMILES('[CH2]C([CH]C)C1[CH]CC1'),
    E0 = (491.479,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3730.22,'J/mol'), sigma=(6.90681,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=582.65 K, Pc=25.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.812359,0.0560205,1.39546e-05,-4.69281e-08,1.89687e-11,59238.1,35.7196], Tmin=(100,'K'), Tmax=(1039.92,'K')), NASAPolynomial(coeffs=[11.5936,0.0465191,-1.84516e-05,3.40733e-09,-2.38525e-13,55267.3,-25.0267], Tmin=(1039.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(491.479,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Cs_S) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C1CC1C[CH][CH]C(5871)',
    structure = SMILES('[CH2]C1CC1C[CH][CH]C'),
    E0 = (499.075,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2883.33,3016.67,3150,900,966.667,1033.33,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.538544,0.0672025,-2.83655e-05,1.02794e-09,1.78971e-12,60156.8,36.3991], Tmin=(100,'K'), Tmax=(1234.06,'K')), NASAPolynomial(coeffs=[10.4709,0.046814,-1.7933e-05,3.14403e-09,-2.09396e-13,56806.4,-17.2498], Tmin=(1234.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.075,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(486.397,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH][CH]C(3874)',
    structure = SMILES('[CH][CH]C'),
    E0 = (522.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1852.85,1853.13,1853.16],'cm^-1')),
        HinderedRotor(inertia=(0.0369856,'amu*angstrom^2'), symmetry=1, barrier=(8.28142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0370949,'amu*angstrom^2'), symmetry=1, barrier=(8.27684,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0718,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.19279,0.0163284,-1.76328e-06,-3.55309e-09,1.20622e-12,62877.3,14.3587], Tmin=(100,'K'), Tmax=(1426.65,'K')), NASAPolynomial(coeffs=[5.59458,0.0149965,-6.04303e-06,1.1011e-09,-7.44871e-14,61642.3,-0.00873183], Tmin=(1426.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(522.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1[CH]CC1(6246)',
    structure = SMILES('[CH2]C1[CH]CC1'),
    E0 = (373.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (68.117,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.78372,0.012651,6.63893e-05,-8.92653e-08,3.36451e-11,45035,19.6308], Tmin=(100,'K'), Tmax=(959.508,'K')), NASAPolynomial(coeffs=[8.39752,0.0249187,-8.55245e-06,1.54895e-09,-1.11544e-13,42315.7,-15.777], Tmin=(959.508,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Isobutyl)"""),
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
    label = '[CH][CH]CC1[CH]CC1(6377)',
    structure = SMILES('[CH][CH]CC1[CH]CC1'),
    E0 = (763.991,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1543,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45762,0.0456387,6.06216e-06,-2.97497e-08,1.17353e-11,91987,30.8432], Tmin=(100,'K'), Tmax=(1087.05,'K')), NASAPolynomial(coeffs=[9.45946,0.0388964,-1.59598e-05,2.9672e-09,-2.07076e-13,88906,-14.5985], Tmin=(1087.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(763.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(394.937,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJC) + radical(CCJ2_triplet)"""),
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
    label = 'C[CH][CH]CC1[C]CC1(6671)',
    structure = SMILES('C[CH][CH]CC1[C]CC1'),
    E0 = (746.83,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.775828,0.0616028,-1.71326e-05,-8.44027e-09,4.36579e-12,89946.5,34.6645], Tmin=(100,'K'), Tmax=(1237.67,'K')), NASAPolynomial(coeffs=[10.5984,0.0465178,-1.9042e-05,3.4645e-09,-2.3582e-13,86239.1,-19.9732], Tmin=(1237.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(746.83,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(RCCJCC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH][C]CC1[CH]CC1(6672)',
    structure = SMILES('C[CH][C]CC1[CH]CC1'),
    E0 = (740.21,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.695191,0.0617209,-1.19922e-05,-1.67329e-08,7.70933e-12,89154.7,34.7308], Tmin=(100,'K'), Tmax=(1146.58,'K')), NASAPolynomial(coeffs=[11.4597,0.0455376,-1.87783e-05,3.46861e-09,-2.39824e-13,85281.6,-24.7998], Tmin=(1146.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(740.21,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C][CH]CC1[CH]CC1(6673)',
    structure = SMILES('C[C][CH]CC1[CH]CC1'),
    E0 = (740.222,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.723358,0.0621818,-1.65374e-05,-1.00818e-08,5.01076e-12,89154.1,34.4657], Tmin=(100,'K'), Tmax=(1227.34,'K')), NASAPolynomial(coeffs=[11.3151,0.0460042,-1.91823e-05,3.53087e-09,-2.42188e-13,85172.7,-24.4308], Tmin=(1227.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(740.222,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJCC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH]C1CC2CCC12(6660)',
    structure = SMILES('C[CH]C1CC2CCC12'),
    E0 = (265.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30883,0.0388722,6.65929e-05,-1.00801e-07,3.72006e-11,32021.4,25.6474], Tmin=(100,'K'), Tmax=(1020.31,'K')), NASAPolynomial(coeffs=[12.1358,0.0468357,-1.92223e-05,3.69219e-09,-2.66875e-13,27188.1,-39.6622], Tmin=(1020.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.288,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(498.868,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + polycyclic(s2_4_4_ane) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]CCC1=CCC1(6674)',
    structure = SMILES('C[CH]CCC1=CCC1'),
    E0 = (229.653,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.635137,0.0602404,3.84231e-06,-3.81681e-08,1.62576e-11,27753.8,32.2748], Tmin=(100,'K'), Tmax=(1047.52,'K')), NASAPolynomial(coeffs=[12.5101,0.0451859,-1.79751e-05,3.32173e-09,-2.32456e-13,23604,-33.4991], Tmin=(1047.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C=CC1CCC1(6675)',
    structure = SMILES('CC=C[CH]C1CCC1'),
    E0 = (166.36,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.985267,0.04467,5.97909e-05,-1.01523e-07,3.95528e-11,20136.2,27.6208], Tmin=(100,'K'), Tmax=(994.349,'K')), NASAPolynomial(coeffs=[14.7947,0.0421123,-1.62933e-05,3.08555e-09,-2.2371e-13,14770.1,-52.0971], Tmin=(994.349,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.36,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(494.711,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(Allyl_S)"""),
)

species(
    label = 'C[CH]CCC1C=CC1(6234)',
    structure = SMILES('C[CH]CCC1C=CC1'),
    E0 = (235.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.74944,0.0563316,1.62532e-05,-5.19032e-08,2.13224e-11,28492.9,31.7447], Tmin=(100,'K'), Tmax=(1024.62,'K')), NASAPolynomial(coeffs=[12.6922,0.0445461,-1.74938e-05,3.23728e-09,-2.2786e-13,24216.8,-35.0867], Tmin=(1024.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.821,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C(C=C)C[CH][CH]C(4572)',
    structure = SMILES('[CH2]C(C=C)C[CH][CH]C'),
    E0 = (474.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,180,400,1161.01,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0380157,'amu*angstrom^2'), symmetry=1, barrier=(3.49622,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3518.52,'J/mol'), sigma=(6.51243,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.58 K, Pc=28.91 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32824,0.067241,-1.17818e-05,-6.33616e-08,5.75461e-11,57141.7,35.6254], Tmin=(100,'K'), Tmax=(513.968,'K')), NASAPolynomial(coeffs=[4.21148,0.0582682,-2.48959e-05,4.62596e-09,-3.19747e-13,56667.5,21.9041], Tmin=(513.968,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(474.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]CC=CC[CH][CH]C(6676)',
    structure = SMILES('[CH2]CC=CC[CH][CH]C'),
    E0 = (471.398,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,213.58,775.834,841.885,3263.82],'cm^-1')),
        HinderedRotor(inertia=(0.0524254,'amu*angstrom^2'), symmetry=1, barrier=(1.78027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524254,'amu*angstrom^2'), symmetry=1, barrier=(1.78027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524254,'amu*angstrom^2'), symmetry=1, barrier=(1.78027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524254,'amu*angstrom^2'), symmetry=1, barrier=(1.78027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524254,'amu*angstrom^2'), symmetry=1, barrier=(1.78027,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0524254,'amu*angstrom^2'), symmetry=1, barrier=(1.78027,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.808217,0.0755754,-6.2481e-05,4.29544e-08,-1.53654e-11,56806,37.8395], Tmin=(100,'K'), Tmax=(715.774,'K')), NASAPolynomial(coeffs=[2.82797,0.0616094,-2.75995e-05,5.23717e-09,-3.65547e-13,56585.5,29.2501], Tmin=(715.774,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(471.398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'CC=CC[C]1[CH]CC1(6539)',
    structure = SMILES('CC=CC[C]1[CH]CC1'),
    E0 = (398.508,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00231,0.057909,-1.23273e-05,-1.02063e-08,4.2639e-12,48043.8,32.8684], Tmin=(100,'K'), Tmax=(1324.58,'K')), NASAPolynomial(coeffs=[10.4379,0.0470292,-1.99533e-05,3.67115e-09,-2.5012e-13,43998.9,-21.1377], Tmin=(1324.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(398.508,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(cyclobutane)"""),
)

species(
    label = 'CC=C[CH]C1[CH]CC1(6537)',
    structure = SMILES('CC=C[CH]C1[CH]CC1'),
    E0 = (354.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14562,0.0452758,4.16314e-05,-7.55609e-08,2.8879e-11,42718.2,29.3946], Tmin=(100,'K'), Tmax=(1023.41,'K')), NASAPolynomial(coeffs=[12.3883,0.0435312,-1.76604e-05,3.35252e-09,-2.40207e-13,38207.2,-35.9013], Tmin=(1023.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Allyl_S)"""),
)

species(
    label = 'CC=CCC1[CH][CH]C1(6541)',
    structure = SMILES('CC=CCC1[CH][CH]C1'),
    E0 = (400.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08692,0.0512884,1.42079e-05,-4.12263e-08,1.55243e-11,48335.6,33.0652], Tmin=(100,'K'), Tmax=(1097.6,'K')), NASAPolynomial(coeffs=[10.4083,0.0466827,-1.96277e-05,3.69918e-09,-2.60061e-13,44520.6,-20.8306], Tmin=(1097.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(400.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[CH2]C=CCC1[CH]CC1(5683)',
    structure = SMILES('[CH2]C=CCC1[CH]CC1'),
    E0 = (364.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3793.9,'J/mol'), sigma=(6.74732,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=592.60 K, Pc=28.02 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953778,0.0479476,4.13492e-05,-8.04677e-08,3.18782e-11,43975.9,31.4689], Tmin=(100,'K'), Tmax=(1004.49,'K')), NASAPolynomial(coeffs=[14.3704,0.0404067,-1.59118e-05,3.01297e-09,-2.17264e-13,38965.6,-44.8418], Tmin=(1004.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(364.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(469.768,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH][CH]C(3856)',
    structure = SMILES('[CH2][CH][CH]C'),
    E0 = (450.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,2810.32],'cm^-1')),
        HinderedRotor(inertia=(0.00169393,'amu*angstrom^2'), symmetry=1, barrier=(9.48541,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.010072,'amu*angstrom^2'), symmetry=1, barrier=(56.4188,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0100766,'amu*angstrom^2'), symmetry=1, barrier=(56.4312,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (55.0984,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.79546,0.0335488,-4.51222e-05,5.0661e-08,-2.13657e-11,54172.4,20.3345], Tmin=(100,'K'), Tmax=(854.972,'K')), NASAPolynomial(coeffs=[-1.26875,0.0345785,-1.53757e-05,2.86235e-09,-1.94746e-13,55524.7,43.1492], Tmin=(854.972,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(450.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH]1[CH]CC1(6173)',
    structure = SMILES('[CH]1[CH]CC1'),
    E0 = (389.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,180,388.387,1395.59,1395.59,1395.59,1395.6,1395.6,1395.6,1395.6,1395.6,1395.61,1802.7],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (54.0904,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.61381,-0.00133339,6.67733e-05,-7.44496e-08,2.51138e-11,46874.5,15.2546], Tmin=(100,'K'), Tmax=(1012.32,'K')), NASAPolynomial(coeffs=[3.9691,0.0238094,-9.81743e-06,1.89417e-09,-1.37335e-13,45442.3,6.81781], Tmin=(1012.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = 'C[CH][CH]C[C]1[CH]CC1(6677)',
    structure = SMILES('C[CH][CH]C[C]1[CH]CC1'),
    E0 = (671.876,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,182.861,1293.05,1293.05,1293.05,1293.05,1293.05,1293.05,1293.05,1293.05,1293.05,1293.05,1293.05,1293.05,1293.05,1293.05,2148.27],'cm^-1')),
        HinderedRotor(inertia=(0.19593,'amu*angstrom^2'), symmetry=1, barrier=(4.59691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19593,'amu*angstrom^2'), symmetry=1, barrier=(4.59691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19593,'amu*angstrom^2'), symmetry=1, barrier=(4.59691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19593,'amu*angstrom^2'), symmetry=1, barrier=(4.59691,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35764,0.0518599,-1.2849e-05,-2.82452e-09,1.07178e-12,80852.8,33.2293], Tmin=(100,'K'), Tmax=(2026.84,'K')), NASAPolynomial(coeffs=[18.3548,0.0386012,-1.65887e-05,2.86307e-09,-1.79567e-13,70606.7,-64.5283], Tmin=(2026.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.876,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(cyclobutane) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]CC1[CH]C[CH]1(6678)',
    structure = SMILES('C[CH][CH]CC1[CH]C[CH]1'),
    E0 = (674.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,180,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,2125.4],'cm^-1')),
        HinderedRotor(inertia=(0.18014,'amu*angstrom^2'), symmetry=1, barrier=(4.14178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18014,'amu*angstrom^2'), symmetry=1, barrier=(4.14178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18014,'amu*angstrom^2'), symmetry=1, barrier=(4.14178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18014,'amu*angstrom^2'), symmetry=1, barrier=(4.14178,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62682,0.0539284,-1.2769e-05,-4.50468e-09,1.72264e-12,81182.3,36.4191], Tmin=(100,'K'), Tmax=(1711.62,'K')), NASAPolynomial(coeffs=[12.0584,0.0453928,-1.91725e-05,3.3972e-09,-2.21372e-13,75290.7,-26.2926], Tmin=(1711.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(674.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(cyclobutane) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH][CH]C1[CH]CC1(6679)',
    structure = SMILES('C[CH][CH][CH]C1[CH]CC1'),
    E0 = (680.996,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,220.416,1287.21,1287.21,1287.21,1287.21,1287.21,1287.21,1287.21,1287.21,1287.21,1287.21,1287.21,1287.21,1287.21,1287.21,2219.43,2219.43],'cm^-1')),
        HinderedRotor(inertia=(0.183003,'amu*angstrom^2'), symmetry=1, barrier=(5.12486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183003,'amu*angstrom^2'), symmetry=1, barrier=(5.12486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183003,'amu*angstrom^2'), symmetry=1, barrier=(5.12486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183003,'amu*angstrom^2'), symmetry=1, barrier=(5.12486,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51885,0.0536478,-1.09701e-05,-6.47697e-09,2.31295e-12,81994.7,37.0189], Tmin=(100,'K'), Tmax=(1606.96,'K')), NASAPolynomial(coeffs=[10.791,0.0469157,-1.9946e-05,3.57756e-09,-2.36168e-13,76903.9,-18.6796], Tmin=(1606.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(680.996,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = 'C[CH][CH]CC1[CH][CH]C1(6680)',
    structure = SMILES('C[CH][CH]CC1[CH][CH]C1'),
    E0 = (674.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,180,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,1262.87,2125.4],'cm^-1')),
        HinderedRotor(inertia=(0.18014,'amu*angstrom^2'), symmetry=1, barrier=(4.14177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18014,'amu*angstrom^2'), symmetry=1, barrier=(4.14177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18014,'amu*angstrom^2'), symmetry=1, barrier=(4.14177,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18014,'amu*angstrom^2'), symmetry=1, barrier=(4.14177,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62682,0.0539285,-1.27691e-05,-4.50451e-09,1.72259e-12,81182.3,36.4191], Tmin=(100,'K'), Tmax=(1711.63,'K')), NASAPolynomial(coeffs=[12.0588,0.0453922,-1.91722e-05,3.39714e-09,-2.21368e-13,75290.5,-26.2951], Tmin=(1711.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(674.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJCC) + radical(cyclobutane) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH][CH]CC1[CH]CC1(6681)',
    structure = SMILES('[CH2][CH][CH]CC1[CH]CC1'),
    E0 = (691.7,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,183.104,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,1261.39,2135.24],'cm^-1')),
        HinderedRotor(inertia=(0.223495,'amu*angstrom^2'), symmetry=1, barrier=(5.2207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223495,'amu*angstrom^2'), symmetry=1, barrier=(5.2207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223495,'amu*angstrom^2'), symmetry=1, barrier=(5.2207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.223495,'amu*angstrom^2'), symmetry=1, barrier=(5.2207,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (108.181,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20244,0.0575542,-1.70201e-05,-3.12917e-09,1.72682e-12,83296.1,37.3823], Tmin=(100,'K'), Tmax=(1563.01,'K')), NASAPolynomial(coeffs=[11.273,0.0456694,-1.8942e-05,3.37526e-09,-2.22795e-13,78451.7,-21.1266], Tmin=(1563.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(691.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(465.61,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C[C]1CCC1(6682)',
    structure = SMILES('C[CH][CH]C[C]1CCC1'),
    E0 = (484.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52898,0.0585798,-1.80517e-05,-1.18771e-09,9.61643e-13,58301.4,33.2019], Tmin=(100,'K'), Tmax=(1856.38,'K')), NASAPolynomial(coeffs=[15.1779,0.043965,-1.81973e-05,3.15774e-09,-2.01728e-13,50684.6,-47.9561], Tmin=(1856.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(484.038,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]C[CH]C1[CH]CC1(6683)',
    structure = SMILES('C[CH]C[CH]C1[CH]CC1'),
    E0 = (486.537,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.954861,0.057115,-2.04495e-06,-2.15639e-08,8.04192e-12,58634.5,36.4586], Tmin=(100,'K'), Tmax=(1215.48,'K')), NASAPolynomial(coeffs=[9.53791,0.0508639,-2.14741e-05,3.98021e-09,-2.74112e-13,54923.2,-13.3081], Tmin=(1215.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.537,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Cs_S) + radical(RCCJC)"""),
)

species(
    label = 'C[CH][CH]CC1C[CH]C1(6684)',
    structure = SMILES('C[CH][CH]CC1C[CH]C1'),
    E0 = (486.454,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,319.028,814.89,814.89,814.89,814.89,814.89,814.89,814.89,814.89,1686.92,1686.92,1686.92,1686.92,1686.92,1686.92,1686.92,1686.92],'cm^-1')),
        HinderedRotor(inertia=(0.132309,'amu*angstrom^2'), symmetry=1, barrier=(3.44019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132309,'amu*angstrom^2'), symmetry=1, barrier=(3.44019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132309,'amu*angstrom^2'), symmetry=1, barrier=(3.44019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.132309,'amu*angstrom^2'), symmetry=1, barrier=(3.44019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07303,0.0579701,-1.08142e-05,-9.6676e-09,3.68507e-12,58617.2,36.055], Tmin=(100,'K'), Tmax=(1398.5,'K')), NASAPolynomial(coeffs=[9.40129,0.0508865,-2.1168e-05,3.82557e-09,-2.56768e-13,54651.1,-12.769], Tmin=(1398.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]CC1[CH]CC1(6346)',
    structure = SMILES('[CH2]C[CH]CC1[CH]CC1'),
    E0 = (497.254,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.826809,0.0596971,-6.89569e-06,-1.82852e-08,7.30454e-12,59928.2,35.9754], Tmin=(100,'K'), Tmax=(1207.61,'K')), NASAPolynomial(coeffs=[10.136,0.0499013,-2.08613e-05,3.85148e-09,-2.64858e-13,56145.7,-17.0439], Tmin=(1207.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.254,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJCC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH]C1CCC1(6685)',
    structure = SMILES('C[CH][CH][CH]C1CCC1'),
    E0 = (493.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,180,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,1345.43,2085.39],'cm^-1')),
        HinderedRotor(inertia=(0.266553,'amu*angstrom^2'), symmetry=1, barrier=(6.12859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266553,'amu*angstrom^2'), symmetry=1, barrier=(6.12859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266553,'amu*angstrom^2'), symmetry=1, barrier=(6.12859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266553,'amu*angstrom^2'), symmetry=1, barrier=(6.12859,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07304,0.0565483,-5.60877e-06,-1.52908e-08,5.53318e-12,59424.6,36.2583], Tmin=(100,'K'), Tmax=(1306.23,'K')), NASAPolynomial(coeffs=[8.84833,0.0516218,-2.16359e-05,3.95634e-09,-2.68722e-13,55782.4,-9.49659], Tmin=(1306.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(493.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = 'C[CH]CC[C]1[CH]CC1(6686)',
    structure = SMILES('C[CH]CC[C]1[CH]CC1'),
    E0 = (477.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.20074,0.0613725,-2.11839e-05,-1.82984e-10,9.11012e-13,57520.9,34.8669], Tmin=(100,'K'), Tmax=(1758.16,'K')), NASAPolynomial(coeffs=[14.707,0.0444202,-1.8474e-05,3.24602e-09,-2.10272e-13,50642.5,-43.9695], Tmin=(1758.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.418,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(RCCJC) + radical(cyclobutane)"""),
)

species(
    label = 'CC[CH][CH]C1[CH]CC1(6344)',
    structure = SMILES('CC[CH][CH]C1[CH]CC1'),
    E0 = (486.549,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,303.053,823.328,823.328,823.328,823.328,823.328,823.328,823.328,823.328,1619.36,1619.36,1619.36,1619.36,1619.36,1619.36,1619.36,1619.36],'cm^-1')),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.54434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.54434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.54434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150921,'amu*angstrom^2'), symmetry=1, barrier=(3.54434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0176,0.0571635,-5.14279e-06,-1.67656e-08,6.10929e-12,58632.4,36.07], Tmin=(100,'K'), Tmax=(1291.73,'K')), NASAPolynomial(coeffs=[9.52804,0.0511594,-2.18013e-05,4.02793e-09,-2.75481e-13,54736,-13.738], Tmin=(1291.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(Cs_S) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][CH]CCC1[CH]CC1(6687)',
    structure = SMILES('[CH2][CH]CCC1[CH]CC1'),
    E0 = (497.242,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2830,2910,2990,3070,3150,900,940,980,1020,1060,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,300,800,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.791978,0.0593034,-2.53476e-06,-2.47646e-08,9.95595e-12,59929.1,36.2652], Tmin=(100,'K'), Tmax=(1141.69,'K')), NASAPolynomial(coeffs=[10.4091,0.0492365,-2.0351e-05,3.76541e-09,-2.606e-13,56193.3,-18.15], Tmin=(1141.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.242,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CCC1[CH]C[CH]1(6688)',
    structure = SMILES('C[CH]CCC1[CH]C[CH]1'),
    E0 = (479.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.927992,0.0588598,-8.38587e-06,-1.45002e-08,5.60955e-12,57828.2,36.351], Tmin=(100,'K'), Tmax=(1280.68,'K')), NASAPolynomial(coeffs=[9.68491,0.0506889,-2.12804e-05,3.90636e-09,-2.66413e-13,54012.3,-14.2033], Tmin=(1280.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(cyclobutane) + radical(RCCJC)"""),
)

species(
    label = 'CC[CH]C[C]1[CH]CC1(6342)',
    structure = SMILES('CC[CH]C[C]1[CH]CC1'),
    E0 = (477.43,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42426,0.0597027,-1.90549e-05,-1.1346e-09,1.02494e-12,57511.5,33.889], Tmin=(100,'K'), Tmax=(1816.72,'K')), NASAPolynomial(coeffs=[15.4541,0.0439399,-1.85304e-05,3.25639e-09,-2.10035e-13,49917.4,-49.0436], Tmin=(1816.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.43,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Tertalkyl) + radical(RCCJCC) + radical(cyclobutane)"""),
)

species(
    label = 'C[CH]CCC1[CH][CH]C1(6689)',
    structure = SMILES('C[CH]CCC1[CH][CH]C1'),
    E0 = (479.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.927992,0.0588598,-8.38587e-06,-1.45002e-08,5.60955e-12,57828.2,36.351], Tmin=(100,'K'), Tmax=(1280.68,'K')), NASAPolynomial(coeffs=[9.68491,0.0506889,-2.12804e-05,3.90636e-09,-2.66413e-13,54012.3,-14.2033], Tmin=(1280.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.834,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(RCCJC) + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = 'CC[CH]CC1[CH]C[CH]1(6343)',
    structure = SMILES('CC[CH]CC1[CH]C[CH]1'),
    E0 = (479.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01344,0.0586352,-1.05232e-05,-1.09224e-08,4.17329e-12,57825.1,35.8815], Tmin=(100,'K'), Tmax=(1373.93,'K')), NASAPolynomial(coeffs=[9.98298,0.0505597,-2.14e-05,3.91102e-09,-2.64567e-13,53657.9,-16.4365], Tmin=(1373.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(cyclobutane) + radical(cyclobutane) + radical(RCCJCC)"""),
)

species(
    label = 'CC[CH]CC1[CH][CH]C1(6345)',
    structure = SMILES('CC[CH]CC1[CH][CH]C1'),
    E0 = (479.846,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,2950,3050,3150,900,950,1000,1050,1100,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,300,800,800,800,800,800,800,800,1600,1600,1600,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01347,0.058635,-1.05225e-05,-1.09232e-08,4.17354e-12,57825.1,35.8814], Tmin=(100,'K'), Tmax=(1373.91,'K')), NASAPolynomial(coeffs=[9.98256,0.0505603,-2.14003e-05,3.91109e-09,-2.64573e-13,53658.1,-16.434], Tmin=(1373.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(479.846,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(RCCJCC) + radical(cyclobutane) + radical(cyclobutane)"""),
)

species(
    label = '[CH2][CH][CH]CC1CCC1(6690)',
    structure = SMILES('[CH2][CH][CH]CC1CCC1'),
    E0 = (503.861,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2816.67,2883.33,2950,3016.67,3083.33,3150,900,933.333,966.667,1000,1033.33,1066.67,1100,2750,2850,1437.5,1250,1305,750,350,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,307.946,804.714,804.714,804.714,804.714,804.714,804.714,804.714,804.714,1623.3,1623.3,1623.3,1623.3,1623.3,1623.3,1623.3,1623.3,1623.3],'cm^-1')),
        HinderedRotor(inertia=(0.152795,'amu*angstrom^2'), symmetry=1, barrier=(3.56821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152795,'amu*angstrom^2'), symmetry=1, barrier=(3.56821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152795,'amu*angstrom^2'), symmetry=1, barrier=(3.56821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152795,'amu*angstrom^2'), symmetry=1, barrier=(3.56821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.881174,0.0590986,-7.43497e-06,-1.66999e-08,6.67732e-12,60720.5,36.1672], Tmin=(100,'K'), Tmax=(1213.48,'K')), NASAPolynomial(coeffs=[9.38703,0.0504645,-2.07476e-05,3.79103e-09,-2.58961e-13,57227.5,-12.4009], Tmin=(1213.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(503.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(490.554,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    E0 = (486.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (486.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (657.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (648.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (659.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (951.903,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (955.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (958.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (952.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (952.015,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (494.738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (549.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (549.854,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (511.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (598.637,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (594.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (619.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (573.765,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (618.448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (584.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (610.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (547.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (839.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (883.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (886.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (892.788,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (886.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (903.492,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (640.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (623.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (645.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (638.742,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (638.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (602.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (633.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (638.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (559.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (604.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (523.032,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (511.181,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (547.247,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (585.868,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['C1=CCC1(4846)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['CC=CCC1[CH]CC1(6235)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_20] for rate rule [Y_12_20b]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C[CH][CH]C[CH]C1CC1(6670)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.31121e+11,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-CsH;C]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH]C)C1[CH]CC1(5828)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C1CC1C[CH][CH]C(5871)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH][CH]C(3874)', '[CH2]C1[CH]CC1(6246)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', '[CH][CH]CC1[CH]CC1(6377)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C[CH][CH]CC1[C]CC1(6671)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['H(3)', 'C[CH][C]CC1[CH]CC1(6672)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'C[C][CH]CC1[CH]CC1(6673)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['C[CH]C1CC2CCC12(6660)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['C[CH]CCC1=CCC1(6674)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['C[CH]C=CC1CCC1(6675)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['C[CH]CCC1C=CC1(6234)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(C=C)C[CH][CH]C(4572)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(2.01e+08,'s^-1'), n=1.02, Ea=(124.265,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 18 used for R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_RR_D;doublebond_intra_pri_2H;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]CC=CC[CH][CH]C(6676)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.61e+08,'s^-1'), n=0.96, Ea=(123.01,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 17 used for R4_Cs_HH_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H
Exact match found for rate rule [R4_Cs_HH_D;doublebond_intra_pri_HNd_Cs;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', 'CC=CC[C]1[CH]CC1(6539)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(7.72e+07,'cm^3/(mol*s)'), n=1.64, Ea=(9.07928,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2568 used for Cds-CsCs_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsCs_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', 'CC=C[CH]C1[CH]CC1(6537)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', 'CC=CCC1[CH][CH]C1(6541)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2]C=CCC1[CH]CC1(5683)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C1=CCC1(4846)', '[CH2][CH][CH]C(3856)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.00858789,'m^3/(mol*s)'), n=2.41179, Ea=(16.3987,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-CsH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]1[CH]CC1(6173)', 'm1_allyl(186)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.270985,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]1[CH]CC1(6173)', '[CH2][CH][CH]C(3856)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1e+08,'m^3/(mol*s)'), n=-1.21605e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_1C-inRing_Ext-1C-R_Sp-3R!H-1C_2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_1C-inRing_Ext-1C-R_Sp-3R!H-1C_2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_1C-inRing_Ext-1C-R_Sp-3R!H-1C_2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', 'C[CH][CH]C[C]1[CH]CC1(6677)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.95338e+06,'m^3/(mol*s)'), n=0.346862, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', 'C[CH][CH]CC1[CH]C[CH]1(6678)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.90676e+06,'m^3/(mol*s)'), n=0.346862, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', 'C[CH][CH][CH]C1[CH]CC1(6679)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.95338e+06,'m^3/(mol*s)'), n=0.346862, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', 'C[CH][CH]CC1[CH][CH]C1(6680)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.95338e+06,'m^3/(mol*s)'), n=0.346862, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', '[CH2][CH][CH]CC1[CH]CC1(6681)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.95338e+06,'m^3/(mol*s)'), n=0.346862, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_2CNO-inRing_Ext-2CNO-R_Ext-3R!H-R_Ext-4R!H-R_Sp-4R!H-3R!H_Ext-3R!H-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['C[CH][CH]C[C]1CCC1(6682)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(7.25e+10,'s^-1'), n=0.6, Ea=(154.39,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_Cs2] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_Cs2]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C[CH]C[CH]C1[CH]CC1(6683)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.91383e+06,'s^-1'), n=2.02407, Ea=(137.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[CH][CH]CC1C[CH]C1(6684)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1.352e+10,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R2H_S_cy4;C_rad_out_H/NonDeC;Cs_H_out_H/(NonDeC/Cs/Cs)]
Euclidian distance = 2.23606797749979
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]C[CH]CC1[CH]CC1(6346)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[CH][CH][CH]C1CCC1(6685)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.93903e+07,'s^-1'), n=1.67208, Ea=(144.885,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3H_SS_23cy4;Y_rad_out;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['C[CH]CC[C]1[CH]CC1(6686)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CC[CH][CH]C1[CH]CC1(6344)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(1.27681e+06,'s^-1'), n=2.16, Ea=(146.492,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][CH]CCC1[CH]CC1(6687)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['C[CH]CCC1[CH]C[CH]1(6688)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(11191.6,'s^-1'), n=2.1395, Ea=(72.9083,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['CC[CH]C[C]1[CH]CC1(6342)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(30253,'s^-1'), n=2.05523, Ea=(118.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R4HJ_1;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['C[CH]CCC1[CH][CH]C1(6689)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(233967,'s^-1'), n=1.68042, Ea=(36.5782,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5H_CCC;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['CC[CH]CC1[CH]C[CH]1(6343)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(0.0564,'s^-1'), n=3.28, Ea=(24.7274,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R5Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R5HJ_1;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C[CH][CH]CC1[CH]CC1(5870)'],
    products = ['CC[CH]CC1[CH][CH]C1(6345)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(92.2,'s^-1'), n=3.21, Ea=(60.7935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R6HJ_1;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][CH][CH]CC1CCC1(6690)'],
    products = ['C[CH][CH]CC1[CH]CC1(5870)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(620613,'s^-1'), n=1.69, Ea=(82.0064,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] + [R6Hall;C_rad_out_2H;Cs_H_out_H/NonDeC] for rate rule [R6Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1495',
    isomers = [
        'C[CH][CH]CC1[CH]CC1(5870)',
    ],
    reactants = [
        ('C1=CCC1(4846)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1495',
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


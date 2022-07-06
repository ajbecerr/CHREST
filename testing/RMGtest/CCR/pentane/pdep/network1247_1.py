species(
    label = 'C[CH][CH]CC[CH][CH]C(3833)',
    structure = SMILES('C[CH][CH]CC[CH][CH]C'),
    E0 = (538.384,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,180,180,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00940648,'amu*angstrom^2'), symmetry=1, barrier=(12.948,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648374,0.0888459,-0.000122779,1.26818e-07,-5.09293e-11,64858.7,40.5153], Tmin=(100,'K'), Tmax=(850.09,'K')), NASAPolynomial(coeffs=[-4.54426,0.0761568,-3.48852e-05,6.51953e-09,-4.43964e-13,67082.9,72.6128], Tmin=(850.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC)"""),
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
    label = 'C[CH][CH]CCC=CC(3842)',
    structure = SMILES('C[CH][CH]CCC=CC'),
    E0 = (265.016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,180,180,986.373,3200],'cm^-1')),
        HinderedRotor(inertia=(0.0855233,'amu*angstrom^2'), symmetry=1, barrier=(1.96635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0855233,'amu*angstrom^2'), symmetry=1, barrier=(1.96635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0855233,'amu*angstrom^2'), symmetry=1, barrier=(1.96635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0855233,'amu*angstrom^2'), symmetry=1, barrier=(1.96635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0855233,'amu*angstrom^2'), symmetry=1, barrier=(1.96635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0855233,'amu*angstrom^2'), symmetry=1, barrier=(1.96635,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80334,0.0670609,-3.09235e-05,6.1269e-09,-4.54089e-13,31933,31.6916], Tmin=(100,'K'), Tmax=(2937.05,'K')), NASAPolynomial(coeffs=[44.7357,0.012789,-5.34994e-06,8.08748e-10,-4.28352e-14,4903.39,-224.77], Tmin=(2937.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(265.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C([CH]C)C[CH][CH]C(3834)',
    structure = SMILES('[CH2]C([CH]C)C[CH][CH]C'),
    E0 = (546.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,180,1264.05,3750.06,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0292152,'amu*angstrom^2'), symmetry=1, barrier=(4.39142,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3517.68,'J/mol'), sigma=(6.71113,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.45 K, Pc=26.41 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.442845,0.0857364,-9.15222e-05,7.75899e-08,-2.89865e-11,65880.5,40.7089], Tmin=(100,'K'), Tmax=(812.499,'K')), NASAPolynomial(coeffs=[1.77095,0.0654815,-2.88057e-05,5.35245e-09,-3.6637e-13,66117.4,37.3635], Tmin=(812.499,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(Cs_S) + radical(RCCJCC) + radical(RCCJC) + radical(Isobutyl)"""),
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
    label = '[CH2]C[CH][CH]C(123)',
    structure = SMILES('[CH2]C[CH][CH]C'),
    E0 = (426.067,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3000,3100,440,815,1455,1000,180,2053.01],'cm^-1')),
        HinderedRotor(inertia=(0.00174415,'amu*angstrom^2'), symmetry=1, barrier=(5.21902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.226977,'amu*angstrom^2'), symmetry=1, barrier=(5.21865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0217106,'amu*angstrom^2'), symmetry=1, barrier=(64.9248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00174461,'amu*angstrom^2'), symmetry=1, barrier=(5.21934,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (69.125,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.20251,0.0328813,-8.80397e-06,-8.44726e-10,4.254e-13,51258.1,21.7496], Tmin=(100,'K'), Tmax=(2245.78,'K')), NASAPolynomial(coeffs=[16.6245,0.0207134,-8.51707e-06,1.3975e-09,-8.32898e-14,42269.4,-60.4535], Tmin=(2245.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(426.067,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    label = '[CH][CH]CC[CH][CH]C(4156)',
    structure = SMILES('[CH][CH]CC[CH][CH]C'),
    E0 = (815.921,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,184.739,389.266,424.784,1786.44,2419.41,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0698845,'amu*angstrom^2'), symmetry=1, barrier=(1.6175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0698845,'amu*angstrom^2'), symmetry=1, barrier=(1.6175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0698845,'amu*angstrom^2'), symmetry=1, barrier=(1.6175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0698845,'amu*angstrom^2'), symmetry=1, barrier=(1.6175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0698845,'amu*angstrom^2'), symmetry=1, barrier=(1.6175,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0698845,'amu*angstrom^2'), symmetry=1, barrier=(1.6175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.1622,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.59681,0.0491516,-1.96929e-05,2.79084e-09,-8.96113e-14,98103.8,26.549], Tmin=(100,'K'), Tmax=(2910.71,'K')), NASAPolynomial(coeffs=[56.094,-0.00804573,2.08078e-06,-4.32046e-10,3.56814e-14,61211.6,-293.683], Tmin=(2910.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(815.921,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(407.409,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
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
    label = 'C[C][CH]CC[CH][CH]C(4580)',
    structure = SMILES('C[C][CH]CC[CH][CH]C'),
    E0 = (792.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,224.051,881.699,1102.12,1322.55,1629.72,1948.58],'cm^-1')),
        HinderedRotor(inertia=(0.114505,'amu*angstrom^2'), symmetry=1, barrier=(3.19786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114505,'amu*angstrom^2'), symmetry=1, barrier=(3.19786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114505,'amu*angstrom^2'), symmetry=1, barrier=(3.19786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114505,'amu*angstrom^2'), symmetry=1, barrier=(3.19786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114505,'amu*angstrom^2'), symmetry=1, barrier=(3.19786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114505,'amu*angstrom^2'), symmetry=1, barrier=(3.19786,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114505,'amu*angstrom^2'), symmetry=1, barrier=(3.19786,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.374139,0.0919519,-0.000123664,1.18822e-07,-4.58505e-11,95392.7,39.3635], Tmin=(100,'K'), Tmax=(838.14,'K')), NASAPolynomial(coeffs=[-0.393389,0.0679387,-3.11566e-05,5.8429e-09,-3.99495e-13,96493.5,48.73], Tmin=(838.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(792.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH][C]CC[CH][CH]C(4581)',
    structure = SMILES('C[CH][C]CC[CH][CH]C'),
    E0 = (792.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,229.966,819.867,1024.83,1229.8,1539.07,1758.93],'cm^-1')),
        HinderedRotor(inertia=(0.137797,'amu*angstrom^2'), symmetry=1, barrier=(3.50631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137797,'amu*angstrom^2'), symmetry=1, barrier=(3.50631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137797,'amu*angstrom^2'), symmetry=1, barrier=(3.50631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137797,'amu*angstrom^2'), symmetry=1, barrier=(3.50631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137797,'amu*angstrom^2'), symmetry=1, barrier=(3.50631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137797,'amu*angstrom^2'), symmetry=1, barrier=(3.50631,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137797,'amu*angstrom^2'), symmetry=1, barrier=(3.50631,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.336553,0.0914789,-0.000118476,1.10434e-07,-4.20035e-11,95393.9,39.6704], Tmin=(100,'K'), Tmax=(834.756,'K')), NASAPolynomial(coeffs=[0.687443,0.0660082,-2.99591e-05,5.60167e-09,-3.82806e-13,96164.2,43.0055], Tmin=(834.756,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(792.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[CH]C1CCC1[CH]C(4559)',
    structure = SMILES('C[CH]C1CCC1[CH]C'),
    E0 = (289.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.828701,0.0524005,3.39291e-05,-6.86001e-08,2.61717e-11,34978.4,32.704], Tmin=(100,'K'), Tmax=(1042.73,'K')), NASAPolynomial(coeffs=[12.5222,0.048737,-2.00592e-05,3.80373e-09,-2.7096e-13,30300.2,-34.9382], Tmin=(1042.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(289.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(515.497,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(Cs_S) + radical(Cs_S)"""),
)

species(
    label = 'C[CH]C=CCC[CH]C(4582)',
    structure = SMILES('C[CH]CC[CH]C=CC'),
    E0 = (211.671,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.461029,0.0750269,-4.02105e-05,1.02416e-08,-1.04509e-12,25587,32.9592], Tmin=(100,'K'), Tmax=(2156,'K')), NASAPolynomial(coeffs=[19.2649,0.0401401,-1.59384e-05,2.73622e-09,-1.74801e-13,17478.8,-72.2047], Tmin=(2156,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.671,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(507.183,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJC) + radical(Allyl_S)"""),
)

species(
    label = 'C[CH][CH]CC=C[CH]C(4583)',
    structure = SMILES('C[CH][CH]C[CH]C=CC'),
    E0 = (406.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,216.895,880.272,1336.05,2038],'cm^-1')),
        HinderedRotor(inertia=(0.0982859,'amu*angstrom^2'), symmetry=1, barrier=(3.00727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0982859,'amu*angstrom^2'), symmetry=1, barrier=(3.00727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0982859,'amu*angstrom^2'), symmetry=1, barrier=(3.00727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0982859,'amu*angstrom^2'), symmetry=1, barrier=(3.00727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0982859,'amu*angstrom^2'), symmetry=1, barrier=(3.00727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0982859,'amu*angstrom^2'), symmetry=1, barrier=(3.00727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73098,0.0645205,-2.95761e-05,5.73968e-09,-4.05128e-13,48912.5,30.8871], Tmin=(100,'K'), Tmax=(2668.31,'K')), NASAPolynomial(coeffs=[36.3986,0.0191695,-7.80232e-06,1.22912e-09,-6.96137e-14,28055.7,-174.803], Tmin=(2668.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(Allyl_S) + radical(RCCJC)"""),
)

species(
    label = 'C=C[CH]CC[CH][CH]C(4584)',
    structure = SMILES('[CH2]C=CCC[CH][CH]C'),
    E0 = (416.516,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,203.721,802.482,1211.16,1619.86],'cm^-1')),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152006,'amu*angstrom^2'), symmetry=1, barrier=(3.5822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.998225,0.0725476,-4.49469e-05,1.68524e-08,-3.19921e-12,50197.8,34.9872], Tmin=(100,'K'), Tmax=(1038.59,'K')), NASAPolynomial(coeffs=[4.31352,0.0597792,-2.65059e-05,5.0152e-09,-3.49866e-13,49509.2,18.8673], Tmin=(1038.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.516,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(482.239,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJCC) + radical(RCCJC) + radical(Allyl_P)"""),
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
    label = 'C[CH][CH][CH]C[CH][CH]C(4585)',
    structure = SMILES('C[CH][CH][CH]C[CH][CH]C'),
    E0 = (732.842,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3012.5,3025,3037.5,3050,390,398.75,407.5,416.25,425,1340,1345,1350,1355,1360,335,343.75,352.5,361.25,370,180,180,4000,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00527302,'amu*angstrom^2'), symmetry=1, barrier=(24.1821,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.901409,0.0893282,-0.00014671,1.62828e-07,-6.60925e-11,88231.8,42.8638], Tmin=(100,'K'), Tmax=(867.327,'K')), NASAPolynomial(coeffs=[-8.81198,0.0801113,-3.73554e-05,6.97101e-09,-4.71511e-13,91948.3,100.055], Tmin=(867.327,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.842,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][CH][CH]CC[CH][CH]C(4586)',
    structure = SMILES('[CH2][CH][CH]CC[CH][CH]C'),
    E0 = (743.631,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,3000,3100,440,815,1455,1000,217.748,2872.67,2957.53,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0300806,'amu*angstrom^2'), symmetry=1, barrier=(3.74535,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.189,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.674542,0.0897095,-0.00013367,1.39584e-07,-5.55415e-11,89541.8,42.901], Tmin=(100,'K'), Tmax=(857.285,'K')), NASAPolynomial(coeffs=[-4.6876,0.0737521,-3.40518e-05,6.36129e-09,-4.32043e-13,91967,76.7267], Tmin=(857.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(743.631,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(478.082,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH]C[CH]C[CH]C(4587)',
    structure = SMILES('C[CH][CH]C[CH]C[CH]C'),
    E0 = (538.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648424,0.0888453,-0.000122776,1.26814e-07,-5.09274e-11,64858.7,41.2083], Tmin=(100,'K'), Tmax=(850.097,'K')), NASAPolynomial(coeffs=[-4.54437,0.0761569,-3.48853e-05,6.51955e-09,-4.43966e-13,67082.9,73.3066], Tmin=(850.097,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2]C[CH]CC[CH][CH]C(4072)',
    structure = SMILES('[CH2]C[CH]CC[CH][CH]C'),
    E0 = (549.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,271.45,1024.88,2442.93,3402.92,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17742,0.0606168,-2.35457e-05,2.9062e-09,-1.05547e-14,66037,30.9364], Tmin=(100,'K'), Tmax=(2709.5,'K')), NASAPolynomial(coeffs=[52.8622,0.00576982,-3.42484e-06,4.75721e-10,-1.88335e-14,32321.2,-270.821], Tmin=(2709.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH]CC[CH]C(4588)',
    structure = SMILES('C[CH][CH][CH]CC[CH]C'),
    E0 = (538.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.648403,0.0888455,-0.000122777,1.26816e-07,-5.09282e-11,64858.7,41.2084], Tmin=(100,'K'), Tmax=(850.094,'K')), NASAPolynomial(coeffs=[-4.54432,0.0761569,-3.48853e-05,6.51954e-09,-4.43965e-13,67082.9,73.3063], Tmin=(850.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.384,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = 'C[CH][CH]C[CH][CH]CC(4070)',
    structure = SMILES('C[CH][CH]C[CH][CH]CC'),
    E0 = (538.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,207.651,2412.36,2756.86,3726.59,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619495,'amu*angstrom^2'), symmetry=1, barrier=(3.67999,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688929,0.0892813,-0.000127823,1.35001e-07,-5.46816e-11,64857.4,40.8911], Tmin=(100,'K'), Tmax=(851.367,'K')), NASAPolynomial(coeffs=[-5.63442,0.0781039,-3.60927e-05,6.76314e-09,-4.60856e-13,67415.9,79.0825], Tmin=(851.367,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CCC[CH][CH]C(4589)',
    structure = SMILES('[CH2][CH]CCC[CH][CH]C'),
    E0 = (549.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,308.264,1179.08,2715.63,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127053,'amu*angstrom^2'), symmetry=1, barrier=(3.49242,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96347,0.0622396,-2.56759e-05,3.94528e-09,-1.67126e-13,66045.5,31.8732], Tmin=(100,'K'), Tmax=(2832.05,'K')), NASAPolynomial(coeffs=[56.6929,0.00151205,-1.54117e-06,1.54127e-10,6.8147e-16,29533,-294.007], Tmin=(2832.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH][CH][CH]C[CH]CC(4071)',
    structure = SMILES('C[CH][CH][CH]C[CH]CC'),
    E0 = (538.396,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,3000,3016.67,3033.33,3050,390,401.667,413.333,425,1340,1346.67,1353.33,1360,335,346.667,358.333,370,207.65,2412.37,2756.83,3726.61,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0619442,'amu*angstrom^2'), symmetry=1, barrier=(3.68025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.688786,0.0892831,-0.00012783,1.35012e-07,-5.46869e-11,64857.4,40.8916], Tmin=(100,'K'), Tmax=(851.35,'K')), NASAPolynomial(coeffs=[-5.63412,0.0781034,-3.60923e-05,6.76306e-09,-4.60849e-13,67415.8,79.0808], Tmin=(851.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.396,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJCC)"""),
)

species(
    label = '[CH2][CH][CH]CCC[CH]C(4590)',
    structure = SMILES('[CH2][CH][CH]CCC[CH]C'),
    E0 = (549.172,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,308.263,1179.1,2686.12,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127056,'amu*angstrom^2'), symmetry=1, barrier=(3.49248,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96347,0.0622395,-2.56759e-05,3.94527e-09,-1.67124e-13,66045.5,31.8732], Tmin=(100,'K'), Tmax=(2832.04,'K')), NASAPolynomial(coeffs=[56.691,0.00151423,-1.54204e-06,1.54282e-10,6.71353e-16,29534.4,-293.994], Tmin=(2832.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH][CH]CC[CH]CC(4073)',
    structure = SMILES('[CH2][CH][CH]CC[CH]CC'),
    E0 = (549.184,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1425,1437.5,1450,1225,1250,1275,1270,1305,1340,700,750,800,300,350,400,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,3000,3100,440,815,1455,1000,271.45,1024.88,2442.93,3402.92,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0954841,'amu*angstrom^2'), symmetry=1, barrier=(2.91548,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.197,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.17742,0.0606168,-2.35457e-05,2.9062e-09,-1.05547e-14,66037,30.9364], Tmin=(100,'K'), Tmax=(2709.5,'K')), NASAPolynomial(coeffs=[52.8622,0.00576982,-3.42484e-06,4.75721e-10,-1.88335e-14,32321.2,-270.821], Tmin=(2709.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.184,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(503.026,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(RCCJCC) + radical(RCCJCC) + radical(RCCJC) + radical(RCCJ)"""),
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
    E0 = (538.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (538.384,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (704.075,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1007.86,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (1003.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (1003.93,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (546.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (601.784,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (625.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (636.082,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (608.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (900.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (944.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (955.423,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (695.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (690.672,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (654.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (684.876,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (690.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (656.433,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (622.821,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (604.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH][CH]CC[CH][CH]C(3833)'],
    products = ['m1_allyl(186)', 'm1_allyl(186)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C[CH][CH]CC[CH][CH]C(3833)'],
    products = ['C[CH][CH]CCC=CC(3842)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.98e+07,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Y_12_20] for rate rule [Y_12_20b]
Euclidian distance = 1.0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C([CH]C)C[CH][CH]C(3834)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.33e+08,'s^-1'), n=1.36, Ea=(157.318,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CsJ-HH;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH][CH]C(3874)', '[CH2]C[CH][CH]C(123)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH3(17)', '[CH][CH]CC[CH][CH]C(4156)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C[C][CH]CC[CH][CH]C(4580)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C[CH][C]CC[CH][CH]C(4581)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[CH][CH]CC[CH][CH]C(3833)'],
    products = ['C[CH]C1CCC1[CH]C(4559)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH][CH]CC[CH][CH]C(3833)'],
    products = ['C[CH]C=CCC[CH]C(4582)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.9748e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['H(3)', 'C[CH][CH]CC=C[CH]C(4583)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', 'C=C[CH]CC[CH][CH]C(4584)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2][CH][CH]C(3856)', 'm1_allyl(186)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(0.135493,'m^3/(mol*s)'), n=1.98996, Ea=(36.5104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2][CH][CH]C(3856)', '[CH2][CH][CH]C(3856)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.625e+08,'m^3/(mol*s)'), n=-0.7, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R_Ext-1C-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'C[CH][CH][CH]C[CH][CH]C(4585)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH2][CH][CH]CC[CH][CH]C(4586)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.76637e+07,'m^3/(mol*s)'), n=0.153073, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.00248871722291, var=1.13870876508, Tref=1000.0, N=5, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R
    Total Standard Deviation in ln(k): 2.14551182899
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_Ext-2CN-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH][CH]CC[CH][CH]C(3833)'],
    products = ['C[CH][CH]C[CH]C[CH]C(4587)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.58236e+06,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C[CH]CC[CH][CH]C(4072)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(34962.3,'s^-1'), n=2.56136, Ea=(141.488,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[CH][CH]CC[CH][CH]C(3833)'],
    products = ['C[CH][CH][CH]CC[CH]C(4588)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.1424e+06,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C[CH][CH]CC[CH][CH]C(3833)'],
    products = ['C[CH][CH]C[CH][CH]CC(4070)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.55363e+06,'s^-1'), n=2.16, Ea=(146.492,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC] for rate rule [R3HJ;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][CH]CCC[CH][CH]C(4589)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.064e+06,'s^-1'), n=1.93, Ea=(141.419,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)] for rate rule [R3HJ;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[CH][CH]CC[CH][CH]C(3833)'],
    products = ['C[CH][CH][CH]C[CH]CC(4071)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(121012,'s^-1'), n=2.05523, Ea=(118.049,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeC;XH_out] for rate rule [R4HJ_1;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][CH][CH]CCC[CH]C(4590)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2778.79,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R6Hall;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][CH][CH]CC[CH]CC(4073)'],
    products = ['C[CH][CH]CC[CH][CH]C(3833)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1062,'s^-1'), n=1.81, Ea=(55.2288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R7Hall;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #1247',
    isomers = [
        'C[CH][CH]CC[CH][CH]C(3833)',
    ],
    reactants = [
        ('m1_allyl(186)', 'm1_allyl(186)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #1247',
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


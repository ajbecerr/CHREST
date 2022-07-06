species(
    label = '[CH2]C[CH][C]=C(279)',
    structure = SMILES('[CH2][C]=CC[CH2]'),
    E0 = (542.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.124939,'amu*angstrom^2'), symmetry=1, barrier=(27.1579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204984,'amu*angstrom^2'), symmetry=1, barrier=(4.71298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00597409,'amu*angstrom^2'), symmetry=1, barrier=(27.137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93243,0.0395004,-1.87977e-05,1.87552e-10,1.75516e-12,65271.7,21.8507], Tmin=(100,'K'), Tmax=(1147.62,'K')), NASAPolynomial(coeffs=[9.44538,0.0226846,-9.06622e-06,1.64916e-09,-1.13131e-13,62930.3,-18.1178], Tmin=(1147.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'C=C(17)',
    structure = SMILES('C=C'),
    E0 = (42.1493,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2334.71,'J/mol'), sigma=(3.971,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.5, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97974,-0.0075756,5.52973e-05,-6.36222e-08,2.31767e-11,5077.46,4.04623], Tmin=(100,'K'), Tmax=(940.446,'K')), NASAPolynomial(coeffs=[5.20299,0.00782443,-2.12684e-06,3.79691e-10,-2.94671e-14,3936.28,-6.6241], Tmin=(940.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(42.1493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-CdsHH) + group(Cds-CdsHH)"""),
)

species(
    label = '[CH]=C=C(264)',
    structure = SMILES('C#C[CH2]'),
    E0 = (328.481,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,750,770,3400,2100,316.125],'cm^-1')),
        HinderedRotor(inertia=(0.00168831,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.32026,0.0108736,8.62047e-06,-1.82971e-08,7.68641e-12,39535.3,8.27852], Tmin=(100,'K'), Tmax=(960.557,'K')), NASAPolynomial(coeffs=[6.38512,0.00814484,-2.78733e-06,4.95345e-10,-3.50146e-14,38483.6,-8.79389], Tmin=(960.557,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(328.481,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Propargyl)"""),
)

species(
    label = 'CH2(T)(15)',
    structure = SMILES('[CH2]'),
    E0 = (381.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1066.91,2790.98,3622.38],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (14.0266,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1197.29,'J/mol'), sigma=(3.8,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.01192,-0.000154978,3.26298e-06,-2.40422e-09,5.69496e-13,45867.7,0.533201], Tmin=(100,'K'), Tmax=(1104.62,'K')), NASAPolynomial(coeffs=[3.14983,0.00296674,-9.76056e-07,1.54115e-10,-9.50339e-15,46058.1,4.77808], Tmin=(1104.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""CH2(T)""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH2][C]=C[CH2](766)',
    structure = SMILES('[CH2][C]=C[CH2]'),
    E0 = (510.938,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.997397,'amu*angstrom^2'), symmetry=1, barrier=(87.1042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.076633,'amu*angstrom^2'), symmetry=1, barrier=(26.8012,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.62091,0.0235623,5.20959e-06,-2.26283e-08,1.00948e-11,61507.3,14.9567], Tmin=(100,'K'), Tmax=(985.3,'K')), NASAPolynomial(coeffs=[8.77437,0.0145495,-5.3793e-06,9.84603e-10,-6.99429e-14,59519.6,-18.5724], Tmin=(985.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.938,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = 'H(9)',
    structure = SMILES('[H]'),
    E0 = (211.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (1.00794,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1205.6,'J/mol'), sigma=(2.05,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,9.24385e-15,-1.3678e-17,6.66185e-21,-1.00107e-24,25474.2,-0.444973], Tmin=(100,'K'), Tmax=(3459.6,'K')), NASAPolynomial(coeffs=[2.5,9.20456e-12,-3.58608e-15,6.15199e-19,-3.92042e-23,25474.2,-0.444973], Tmin=(3459.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(211.805,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""H""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = '[CH]CC=[C][CH2](767)',
    structure = SMILES('[CH]CC=[C][CH2]'),
    E0 = (785.009,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,374.388,374.417,374.454,1854.97],'cm^-1')),
        HinderedRotor(inertia=(0.00120283,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0012028,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.372822,'amu*angstrom^2'), symmetry=1, barrier=(37.0813,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.82458,0.0420401,-2.86432e-05,8.6655e-09,-6.57453e-13,94497.9,21.4551], Tmin=(100,'K'), Tmax=(1180.81,'K')), NASAPolynomial(coeffs=[10.5919,0.0188326,-7.40918e-06,1.3331e-09,-9.07987e-14,91974.8,-24.2159], Tmin=(1180.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(785.009,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=[C][CH]C[CH2](752)',
    structure = SMILES('[CH][C]=CC[CH2]'),
    E0 = (761.225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,515.587,515.6,515.61,515.611],'cm^-1')),
        HinderedRotor(inertia=(0.280388,'amu*angstrom^2'), symmetry=1, barrier=(52.9053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280431,'amu*angstrom^2'), symmetry=1, barrier=(52.9053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.280438,'amu*angstrom^2'), symmetry=1, barrier=(52.9054,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0453,0.0405274,-2.17286e-05,5.57859e-09,-5.79171e-13,91626.3,22.277], Tmin=(100,'K'), Tmax=(2104.67,'K')), NASAPolynomial(coeffs=[11.5253,0.0225101,-8.88762e-06,1.51109e-09,-9.60148e-14,87635.9,-30.5132], Tmin=(2104.67,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(761.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1=CCC1(765)',
    structure = SMILES('C=C1[CH]CC1'),
    E0 = (257.127,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.90085,0.0097037,6.85862e-05,-8.87071e-08,3.23833e-11,30977.6,13.798], Tmin=(100,'K'), Tmax=(987.001,'K')), NASAPolynomial(coeffs=[8.44602,0.0241142,-9.36782e-06,1.80785e-09,-1.33445e-13,28086.5,-21.9826], Tmin=(987.001,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.127,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(Allyl_S)"""),
)

species(
    label = '[CH2]C=CC=C(485)',
    structure = SMILES('[CH2]C=CC=C'),
    E0 = (175.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29872,0.0238033,3.93748e-05,-6.9319e-08,2.88424e-11,21235.8,16.7722], Tmin=(100,'K'), Tmax=(942.049,'K')), NASAPolynomial(coeffs=[11.9105,0.0173606,-5.09266e-06,8.77912e-10,-6.40318e-14,17899.7,-37.1199], Tmin=(942.049,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.951,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C][CH]C=C(615)',
    structure = SMILES('[CH2]C=C[C]=C'),
    E0 = (374.947,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(2.0825,'amu*angstrom^2'), symmetry=1, barrier=(47.8809,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07959,'amu*angstrom^2'), symmetry=1, barrier=(47.8139,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2251,0.0302874,1.05284e-05,-3.68284e-08,1.72731e-11,45167.7,17.3268], Tmin=(100,'K'), Tmax=(923.512,'K')), NASAPolynomial(coeffs=[10.3879,0.0174193,-5.09535e-06,8.16559e-10,-5.51188e-14,42701,-26.5962], Tmin=(923.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(374.947,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C#CC[CH2](768)',
    structure = SMILES('[CH2]C#CC[CH2]'),
    E0 = (469.054,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2100,2250,500,550,180],'cm^-1')),
        HinderedRotor(inertia=(0.0330171,'amu*angstrom^2'), symmetry=1, barrier=(85.7513,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0170365,'amu*angstrom^2'), symmetry=1, barrier=(3.59874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.415122,'amu*angstrom^2'), symmetry=1, barrier=(85.7834,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.31954,0.0319385,-7.3161e-06,-1.03319e-08,5.78688e-12,56479.1,20.2157], Tmin=(100,'K'), Tmax=(984.077,'K')), NASAPolynomial(coeffs=[7.75522,0.0209666,-7.54605e-06,1.30944e-09,-8.8396e-14,54870.7,-8.65748], Tmin=(984.077,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.054,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(RCCJ) + radical(Propargyl)"""),
)

species(
    label = '[CH]=[C][CH2](751)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,228.26,228.889,229.072],'cm^-1')),
        HinderedRotor(inertia=(1.35217,'amu*angstrom^2'), symmetry=1, barrier=(50.6528,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.27542,0.0127953,9.49555e-06,-1.56031e-08,5.42959e-12,73954,11.3502], Tmin=(100,'K'), Tmax=(1063.3,'K')), NASAPolynomial(coeffs=[4.1896,0.0168436,-6.77768e-06,1.2222e-09,-8.33566e-14,73336.4,4.8934], Tmin=(1063.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2][CH2](18)',
    structure = SMILES('[CH2][CH2]'),
    E0 = (313.796,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1608.02,1609.59,1610.81,1612.2,1612.86],'cm^-1')),
        HinderedRotor(inertia=(0.00575811,'amu*angstrom^2'), symmetry=1, barrier=(10.593,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86546,-0.00370153,4.71859e-05,-6.15022e-08,2.54427e-11,37752.3,8.63154], Tmin=(100,'K'), Tmax=(845.133,'K')), NASAPolynomial(coeffs=[5.8996,0.00441368,1.29131e-06,-4.57987e-10,3.67865e-14,36774.8,-4.58851], Tmin=(845.133,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.796,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCJ) + radical(CCJ)"""),
)

species(
    label = '[CH2][C]=C[CH][CH2](769)',
    structure = SMILES('[CH2][C]=C[CH][CH2]'),
    E0 = (683.153,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,1114.42],'cm^-1')),
        HinderedRotor(inertia=(0.00740552,'amu*angstrom^2'), symmetry=1, barrier=(6.52457,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.14018,'amu*angstrom^2'), symmetry=1, barrier=(95.191,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0340784,'amu*angstrom^2'), symmetry=1, barrier=(30.0845,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16721,0.0338824,-8.66428e-06,-9.5211e-09,5.10585e-12,82235.7,20.5922], Tmin=(100,'K'), Tmax=(1070.03,'K')), NASAPolynomial(coeffs=[9.10766,0.0208121,-8.38972e-06,1.55219e-09,-1.08387e-13,80013.4,-16.8053], Tmin=(1070.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(683.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C[CH2](770)',
    structure = SMILES('[CH2][C]=[C]C[CH2]'),
    E0 = (779.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,857.991],'cm^-1')),
        HinderedRotor(inertia=(0.000811113,'amu*angstrom^2'), symmetry=1, barrier=(9.20943,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19936,'amu*angstrom^2'), symmetry=1, barrier=(27.5757,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.27328,'amu*angstrom^2'), symmetry=1, barrier=(52.2671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22653,0.0403139,-3.06795e-05,1.31961e-08,-2.47685e-12,93860.8,21.2647], Tmin=(100,'K'), Tmax=(1200.7,'K')), NASAPolynomial(coeffs=[7.09434,0.0240975,-1.0421e-05,1.94808e-09,-1.34904e-13,92691.8,-3.11005], Tmin=(1200.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(779.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[CH]C(771)',
    structure = SMILES('[CH2][C]=C[CH]C'),
    E0 = (477.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,300.801],'cm^-1')),
        HinderedRotor(inertia=(0.00186762,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00186424,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.706267,'amu*angstrom^2'), symmetry=1, barrier=(45.3693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1366,0.033084,1.9467e-06,-2.18723e-08,9.53324e-12,57552.8,18.2216], Tmin=(100,'K'), Tmax=(1034.94,'K')), NASAPolynomial(coeffs=[9.21568,0.0232651,-9.24634e-06,1.71508e-09,-1.20643e-13,55148,-20.7126], Tmin=(1034.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]C[CH2](772)',
    structure = SMILES('[CH2]C=[C]C[CH2]'),
    E0 = (542.04,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.124939,'amu*angstrom^2'), symmetry=1, barrier=(27.1579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204984,'amu*angstrom^2'), symmetry=1, barrier=(4.71298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00597409,'amu*angstrom^2'), symmetry=1, barrier=(27.137,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93243,0.0395004,-1.87977e-05,1.87552e-10,1.75516e-12,65271.7,21.8507], Tmin=(100,'K'), Tmax=(1147.62,'K')), NASAPolynomial(coeffs=[9.44538,0.0226846,-9.06622e-06,1.64916e-09,-1.13131e-13,62930.3,-18.1178], Tmin=(1147.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.04,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]CC(773)',
    structure = SMILES('[CH2][C]=[C]CC'),
    E0 = (574.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1670,1700,300,440,199.161],'cm^-1')),
        HinderedRotor(inertia=(1.13179,'amu*angstrom^2'), symmetry=1, barrier=(31.7593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.54973,'amu*angstrom^2'), symmetry=1, barrier=(15.4267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00282898,'amu*angstrom^2'), symmetry=1, barrier=(13.2995,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0313,0.0414391,-2.66012e-05,8.87805e-09,-1.2424e-12,69184.7,20.1774], Tmin=(100,'K'), Tmax=(1590.11,'K')), NASAPolynomial(coeffs=[9.24533,0.0232918,-9.48232e-06,1.70084e-09,-1.13989e-13,66890.5,-17.972], Tmin=(1590.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(574.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH][CH]C=C(491)',
    structure = SMILES('[CH2][CH]C=C[CH2]'),
    E0 = (445.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.22073,0.0290685,1.67269e-05,-3.93351e-08,1.61886e-11,53631.5,19.2316], Tmin=(100,'K'), Tmax=(997.188,'K')), NASAPolynomial(coeffs=[10.0213,0.0217642,-8.36598e-06,1.562e-09,-1.11865e-13,50883.2,-24.3592], Tmin=(997.188,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.311,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(RCCJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C[C]=[C]C(569)',
    structure = SMILES('[CH2]C[C]=[C]C'),
    E0 = (628.383,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1670,1700,300,440,180],'cm^-1')),
        HinderedRotor(inertia=(0.171531,'amu*angstrom^2'), symmetry=1, barrier=(3.94384,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171044,'amu*angstrom^2'), symmetry=1, barrier=(3.93265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17057,'amu*angstrom^2'), symmetry=1, barrier=(3.92174,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.68332,0.034996,4.68072e-06,-7.11492e-08,6.63693e-11,75618.1,19.4749], Tmin=(100,'K'), Tmax=(449.867,'K')), NASAPolynomial(coeffs=[4.1381,0.0313544,-1.41646e-05,2.69906e-09,-1.8888e-13,75393.2,12.5733], Tmin=(449.867,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(628.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][CH]C=[C]C(568)',
    structure = SMILES('[CH2][CH]C=[C]C'),
    E0 = (531.653,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,374.37],'cm^-1')),
        HinderedRotor(inertia=(0.00197819,'amu*angstrom^2'), symmetry=1, barrier=(5.29133,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.052747,'amu*angstrom^2'), symmetry=1, barrier=(5.27016,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36375,'amu*angstrom^2'), symmetry=1, barrier=(35.9437,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1192,0.0367177,-1.73813e-05,2.65723e-09,1.62936e-13,64014.7,20.4998], Tmin=(100,'K'), Tmax=(1486.58,'K')), NASAPolynomial(coeffs=[9.41176,0.0228971,-9.29002e-06,1.65391e-09,-1.09829e-13,61205.4,-19.7301], Tmin=(1486.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(531.653,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Allyl_S) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1CC1(774)',
    structure = SMILES('C=[C]C1CC1'),
    E0 = (352.164,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.58514,0.0159726,5.69898e-05,-8.4516e-08,3.32632e-11,42420.3,17.1827], Tmin=(100,'K'), Tmax=(955.353,'K')), NASAPolynomial(coeffs=[11.3074,0.017655,-5.63269e-06,1.03945e-09,-7.82738e-14,39010.3,-33.623], Tmin=(955.353,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(352.164,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC=C(275)',
    structure = SMILES('C=[C]CC=C'),
    E0 = (326.592,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,375.103,375.104],'cm^-1')),
        HinderedRotor(inertia=(0.121697,'amu*angstrom^2'), symmetry=1, barrier=(12.151,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121699,'amu*angstrom^2'), symmetry=1, barrier=(12.151,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.18209,0.0328829,-1.51765e-06,-1.69932e-08,7.62913e-12,39351.4,19.8243], Tmin=(100,'K'), Tmax=(1050.11,'K')), NASAPolynomial(coeffs=[8.81968,0.0228909,-9.08742e-06,1.6794e-09,-1.17556e-13,37114.3,-16.5372], Tmin=(1050.11,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.592,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=CC(775)',
    structure = SMILES('C=[C]C=CC'),
    E0 = (256.891,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92118,0.0394904,-1.38895e-05,-9.15945e-09,6.40789e-12,30977.3,16.8588], Tmin=(100,'K'), Tmax=(973.439,'K')), NASAPolynomial(coeffs=[9.77405,0.021281,-7.49405e-06,1.29737e-09,-8.80426e-14,28782.3,-24.2368], Tmin=(973.439,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(256.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C=CC[CH2](750)',
    structure = SMILES('C#C[CH]C[CH2]'),
    E0 = (477.258,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2175,525,750,770,3400,2100,180],'cm^-1')),
        HinderedRotor(inertia=(0.0041749,'amu*angstrom^2'), symmetry=1, barrier=(6.96764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.94423,'amu*angstrom^2'), symmetry=1, barrier=(67.6936,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1369,'amu*angstrom^2'), symmetry=1, barrier=(67.7355,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (66.1011,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76795,0.0430915,-3.64273e-05,1.77222e-08,-3.47094e-12,57486.3,19.7996], Tmin=(100,'K'), Tmax=(1346.56,'K')), NASAPolynomial(coeffs=[9.46879,0.0175599,-5.02767e-06,7.11761e-10,-4.0866e-14,55653.2,-18.7497], Tmin=(1346.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.258,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(Sec_Propargyl) + radical(RCCJ)"""),
)

species(
    label = 'C=[C]C=C(776)',
    structure = SMILES('C=[C]C=C'),
    E0 = (292.917,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(1.58595,'amu*angstrom^2'), symmetry=1, barrier=(36.4641,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (53.0825,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.60776,0.0234236,8.49824e-06,-3.01406e-08,1.43436e-11,35286.5,12.5895], Tmin=(100,'K'), Tmax=(918.276,'K')), NASAPolynomial(coeffs=[9.56649,0.0118208,-3.11028e-06,4.74853e-10,-3.2022e-14,33219.7,-24.6842], Tmin=(918.276,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(292.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][CH]C[C]=C(278)',
    structure = SMILES('[CH2][CH]C[C]=C'),
    E0 = (597.232,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,1361.14,1361.33],'cm^-1')),
        HinderedRotor(inertia=(0.157829,'amu*angstrom^2'), symmetry=1, barrier=(3.62879,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157827,'amu*angstrom^2'), symmetry=1, barrier=(3.62875,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157809,'amu*angstrom^2'), symmetry=1, barrier=(3.62835,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.86901,0.0335619,-1.56825e-05,3.25552e-09,-2.60722e-13,71861.2,21.3135], Tmin=(100,'K'), Tmax=(2708.28,'K')), NASAPolynomial(coeffs=[12.4387,0.019428,-7.85429e-06,1.32855e-09,-8.28455e-14,66677.7,-34.3892], Tmin=(2708.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(597.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJC) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]C[CH2](490)',
    structure = SMILES('[CH]C=CC[CH2]'),
    E0 = (523.383,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.98638,0.037094,-1.47885e-06,-1.70765e-08,7.2494e-12,63027,22.0103], Tmin=(100,'K'), Tmax=(1082.65,'K')), NASAPolynomial(coeffs=[8.03056,0.0297375,-1.20336e-05,2.19827e-09,-1.51499e-13,60840.7,-11.6824], Tmin=(1082.65,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.383,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(RCCJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC[CH2](280)',
    structure = SMILES('[CH]=[C]CC[CH2]'),
    E0 = (649.882,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,1685,370,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(1.38393,'amu*angstrom^2'), symmetry=1, barrier=(31.8193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000321969,'amu*angstrom^2'), symmetry=1, barrier=(3.65565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.714807,'amu*angstrom^2'), symmetry=1, barrier=(16.4348,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3052.11,'J/mol'), sigma=(5.53315,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=476.73 K, Pc=40.88 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07344,0.0423233,-2.94338e-05,1.10648e-08,-1.78046e-12,78231.9,21.8978], Tmin=(100,'K'), Tmax=(1387.7,'K')), NASAPolynomial(coeffs=[8.19337,0.0246827,-1.03656e-05,1.90418e-09,-1.30139e-13,76533.4,-9.63256], Tmin=(1387.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.882,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(RCCJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]CC(753)',
    structure = SMILES('[CH][C]=CCC'),
    E0 = (555.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,1685,370,397.272,399.397,399.416,400.141],'cm^-1')),
        HinderedRotor(inertia=(0.463609,'amu*angstrom^2'), symmetry=1, barrier=(52.0878,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.460822,'amu*angstrom^2'), symmetry=1, barrier=(52.0283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.459987,'amu*angstrom^2'), symmetry=1, barrier=(52.1211,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.1091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85888,0.0414772,-1.68835e-05,2.79592e-10,1.01452e-12,66950.4,21.1656], Tmin=(100,'K'), Tmax=(1338.63,'K')), NASAPolynomial(coeffs=[8.63268,0.0290645,-1.17464e-05,2.08992e-09,-1.39464e-13,64435.5,-16.1094], Tmin=(1338.63,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(555.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'Ar',
    structure = SMILES('[Ar]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (39.348,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(1134.93,'J/mol'), sigma=(3.33,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,4.37967], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""Ar""", comment="""Thermo library: primaryThermoLibrary"""),
)

species(
    label = 'He',
    structure = SMILES('[He]'),
    E0 = (-6.19738,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (4.0026,'amu'),
    collisionModel = TransportData(shapeIndex=0, epsilon=(84.8076,'J/mol'), sigma=(2.576,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.5,0,0,0,0,-745.375,0.928724], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-6.19738,'kJ/mol'), Cp0=(20.7862,'J/(mol*K)'), CpInf=(20.7862,'J/(mol*K)'), label="""He""", comment="""Thermo library: primaryThermoLibrary"""),
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
    E0 = (-8.64289,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0135,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(810.913,'J/mol'), sigma=(3.621,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.76,'angstroms^3'), rotrelaxcollnum=4.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.53101,-0.000123661,-5.02999e-07,2.43531e-09,-1.40881e-12,-1046.98,2.96747], Tmin=(200,'K'), Tmax=(1000,'K')), NASAPolynomial(coeffs=[2.95258,0.0013969,-4.92632e-07,7.8601e-11,-4.60755e-15,-923.949,5.87189], Tmin=(1000,'K'), Tmax=(6000,'K'))], Tmin=(200,'K'), Tmax=(6000,'K'), E0=(-8.64289,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""N2""", comment="""Thermo library: primaryThermoLibrary"""),
)

transitionState(
    label = 'TS1',
    E0 = (542.04,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (949.116,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (996.813,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (973.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (550.324,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (605.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (594.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (695.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (678.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (664.962,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (928.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (894.957,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (991.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (677.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (746.177,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (726.514,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (735.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (790.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (673.165,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (547.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (564.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (564.901,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (704.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (713.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (755.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (737.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (795.067,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (589.019,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C[CH][C]=C(279)'],
    products = ['C=C(17)', '[CH]=C=C(264)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(T)(15)', '[CH2][C]=C[CH2](766)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.22604e+07,'m^3/(mol*s)'), n=-0.594572, Ea=(56.8079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['H(9)', '[CH]CC=[C][CH2](767)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(9)', '[CH]=[C][CH]C[CH2](752)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH2]C[CH][C]=C(279)'],
    products = ['[CH2]C1=CCC1(765)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSD;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C[CH][C]=C(279)'],
    products = ['[CH2]C=CC=C(485)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(9)', 'C=[C][CH]C=C(615)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(105.529,'m^3/(mol*s)'), n=1.6629, Ea=(8.08712,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-OneDeH_Cds;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(9)', '[CH2]C#CC[CH2](768)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C(17)', '[CH]=[C][CH2](751)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(0.00310793,'m^3/(mol*s)'), n=2.49201, Ea=(21.4614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-HH;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][CH2](18)', '[CH]=C=C(264)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.04713,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2][CH2](18)', '[CH]=[C][CH2](751)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.9668e+08,'m^3/(mol*s)'), n=-0.557189, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00108714239892, var=1.14816174436, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R
    Total Standard Deviation in ln(k): 2.15085144951
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R_N-3R!H->O_Ext-1C-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(9)', '[CH2][C]=C[CH][CH2](769)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(9)', '[CH2][C]=[C]C[CH2](770)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.62598e+07,'m^3/(mol*s)'), n=0.255122, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_Sp-4R!H=3C_N-3C-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_Sp-4R!H=3C_N-3C-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_N-2CNO->O_3R!H->C_Sp-3C-2CN_Ext-3C-R_Sp-4R!H=3C_N-3C-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2][C]=C[CH]C(771)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(6.18e+09,'s^-1'), n=1.22, Ea=(199.995,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 157 used for R2H_S;C_rad_out_H/Cd;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C=[C]C[CH2](772)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C]=[C]CC(773)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.304e+09,'s^-1'), n=1.24, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C[CH][C]=C(279)'],
    products = ['[CH2][CH][CH]C=C(491)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(70985.8,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C[C]=[C]C(569)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(7.74e+09,'s^-1'), n=1.08, Ea=(161.921,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_Cs;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C[CH][C]=C(279)'],
    products = ['[CH2][CH]C=[C]C(568)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(869.832,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C[CH][C]=C(279)'],
    products = ['C=[C]C1CC1(774)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C[CH][C]=C(279)'],
    products = ['C=[C]CC=C(275)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C[CH][C]=C(279)'],
    products = ['C=[C]C=CC(775)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(9)', '[CH]=C=CC[CH2](750)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['CH2(T)(15)', 'C=[C]C=C(776)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][CH]C[C]=C(278)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.7863e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2]C[CH][C]=C(279)'],
    products = ['[CH]=C[CH]C[CH2](490)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(9.55842e+08,'s^-1'), n=1.508, Ea=(195.122,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=[C]CC[CH2](280)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_H/NonDeC] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C][CH]CC(753)'],
    products = ['[CH2]C[CH][C]=C(279)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using an average for rate rule [R5Hall;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #183',
    isomers = [
        '[CH2]C[CH][C]=C(279)',
    ],
    reactants = [
        ('C=C(17)', '[CH]=C=C(264)'),
    ],
    bathGas = {
        'Ar': 0.25,
        'He': 0.25,
        'Ne': 0.25,
        'N2': 0.25,
    },
)

pressureDependence(
    label = 'PDepNetwork #183',
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


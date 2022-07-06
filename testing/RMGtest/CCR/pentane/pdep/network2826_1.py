species(
    label = 'C=C([O])C([O])=C[O](11158)',
    structure = SMILES('C=C([O])C([O])=C[O]'),
    E0 = (-125.585,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.04411,'amu*angstrom^2'), symmetry=1, barrier=(24.0062,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0729095,0.0724912,-8.7986e-05,4.83801e-08,-9.6461e-12,-14943,25.3649], Tmin=(100,'K'), Tmax=(1456.33,'K')), NASAPolynomial(coeffs=[21.1206,-0.00117263,3.8032e-06,-9.24506e-10,6.85246e-14,-19477.2,-79.2224], Tmin=(1456.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-125.585,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = 'OCHCO(3676)',
    structure = SMILES('O=[C]C=O'),
    E0 = (-75.5464,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,180,525.376,1512.41,1512.65,1513.37],'cm^-1')),
        HinderedRotor(inertia=(0.00619061,'amu*angstrom^2'), symmetry=1, barrier=(0.260399,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3660.03,'J/mol'), sigma=(4.01,'angstroms'), dipoleMoment=(0,'De'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=2.0, comment="""NOx2018"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.25386,0.0154074,-1.14326e-05,4.29104e-09,-6.58698e-13,-9058.53,11.1539], Tmin=(100,'K'), Tmax=(1504.29,'K')), NASAPolynomial(coeffs=[6.43153,0.00695778,-3.00696e-06,5.56991e-10,-3.81281e-14,-10014.6,-5.47402], Tmin=(1504.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-75.5464,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), label="""OCHCO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C([O])=[C]C=O(9910)',
    structure = SMILES('C=C([O])[C]=C[O]'),
    E0 = (149.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.55298,'amu*angstrom^2'), symmetry=1, barrier=(35.706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.562656,0.0598018,-6.95646e-05,3.74263e-08,-7.30439e-12,18117.9,21.9387], Tmin=(100,'K'), Tmax=(1502.94,'K')), NASAPolynomial(coeffs=[17.396,0.00121988,2.6563e-06,-7.09798e-10,5.39892e-14,14614.4,-60.9531], Tmin=(1502.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]C([O])=C[O](11900)',
    structure = SMILES('C=[C]C([O])=C[O]'),
    E0 = (149.499,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.55298,'amu*angstrom^2'), symmetry=1, barrier=(35.706,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.562656,0.0598018,-6.95646e-05,3.74263e-08,-7.30439e-12,18117.9,21.9387], Tmin=(100,'K'), Tmax=(1502.94,'K')), NASAPolynomial(coeffs=[17.396,0.00121988,2.6563e-06,-7.09798e-10,5.39892e-14,14614.4,-60.9531], Tmin=(1502.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(149.499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C([O])C(=C)[O](2705)',
    structure = SMILES('[CH]=C([O])C(=C)[O]'),
    E0 = (188.841,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.74757,'amu*angstrom^2'), symmetry=1, barrier=(17.1881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4088.81,'J/mol'), sigma=(6.45843,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=638.66 K, Pc=34.44 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26098,0.0578202,-6.96775e-05,3.50958e-08,-4.46465e-12,22813.5,19.7231], Tmin=(100,'K'), Tmax=(835.128,'K')), NASAPolynomial(coeffs=[15.2434,0.00498286,-1.61285e-07,-1.31999e-10,1.43248e-14,19985.1,-48.1658], Tmin=(835.128,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.841,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C1OOC1=C[O](11845)',
    structure = SMILES('C=C1OOC1=C[O]'),
    E0 = (111.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63238,0.042926,-2.03256e-05,-7.53193e-09,6.43746e-12,13467.9,20.5978], Tmin=(100,'K'), Tmax=(1026.72,'K')), NASAPolynomial(coeffs=[13.9939,0.0138644,-5.76831e-06,1.1323e-09,-8.34927e-14,9922.95,-44.2672], Tmin=(1026.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.203,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(C=COJ)"""),
)

species(
    label = 'C=C([O])C1=COO1(11901)',
    structure = SMILES('C=C([O])C1=COO1'),
    E0 = (114.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30029,0.0486513,-2.58897e-05,-1.23038e-08,1.12255e-11,13877.8,20.7412], Tmin=(100,'K'), Tmax=(942.821,'K')), NASAPolynomial(coeffs=[17.4613,0.00665124,-1.33184e-06,2.15498e-10,-1.82128e-14,9649.73,-62.5359], Tmin=(942.821,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C1OOC=C1[O](11902)',
    structure = SMILES('C=C1OOC=C1[O]'),
    E0 = (18.8076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.95798,0.0270659,3.82398e-05,-8.0866e-08,3.63734e-11,2352.1,19.1492], Tmin=(100,'K'), Tmax=(924.496,'K')), NASAPolynomial(coeffs=[17.9743,0.00382202,1.23067e-06,-2.94916e-10,1.45712e-14,-2577.4,-67.5071], Tmin=(924.496,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(18.8076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C([O])C(O)=C=O(11903)',
    structure = SMILES('C=C([O])C(O)=C=O'),
    E0 = (-216.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.223029,0.0742128,-9.77354e-05,5.90701e-08,-1.31553e-11,-25850.8,22.9157], Tmin=(100,'K'), Tmax=(1254.96,'K')), NASAPolynomial(coeffs=[19.9302,0.00100599,2.18797e-06,-6.10738e-10,4.83209e-14,-29978.8,-73.3751], Tmin=(1254.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-216.131,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(O)C([O])=C=O(11904)',
    structure = SMILES('C=C(O)C(=O)[C]=O'),
    E0 = (-229.156,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.841057,0.0647585,-7.86609e-05,4.58289e-08,-1.02035e-11,-27443.3,23.1978], Tmin=(100,'K'), Tmax=(1114.75,'K')), NASAPolynomial(coeffs=[16.5494,0.00839318,-2.8164e-06,4.70957e-10,-3.13422e-14,-30945.5,-54.2925], Tmin=(1114.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-229.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([O])C([O])=[C][O](11905)',
    structure = SMILES('[CH2]C([O])C([O])=[C][O]'),
    E0 = (348.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,442.418,443.009,443.454,444.264],'cm^-1')),
        HinderedRotor(inertia=(0.000854147,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0944573,'amu*angstrom^2'), symmetry=1, barrier=(13.1764,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.614758,0.0648882,-7.64177e-05,4.2237e-08,-8.75302e-12,42018.7,31.8926], Tmin=(100,'K'), Tmax=(1297.62,'K')), NASAPolynomial(coeffs=[18.5298,0.0043291,-2.46972e-07,-6.47605e-11,7.24171e-15,37818.5,-57.4741], Tmin=(1297.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(348.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([O])C([O])[CH][O](11906)',
    structure = SMILES('[CH]=C([O])C([O])[CH][O]'),
    E0 = (457.55,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,420.957,421.746,422.457,3928.67],'cm^-1')),
        HinderedRotor(inertia=(0.0597868,'amu*angstrom^2'), symmetry=1, barrier=(7.56726,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183402,'amu*angstrom^2'), symmetry=1, barrier=(22.8508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.861635,0.0738398,-0.000118138,9.64368e-08,-3.06321e-11,55139,29.4937], Tmin=(100,'K'), Tmax=(846.71,'K')), NASAPolynomial(coeffs=[11.3689,0.0180821,-8.51889e-06,1.59101e-09,-1.07596e-13,53579,-18.1541], Tmin=(846.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_P)"""),
)

species(
    label = 'C[C]([O])C([O])=[C][O](11907)',
    structure = SMILES('C[C]([O])C([O])=[C][O]'),
    E0 = (313.321,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,1685,370,373.13,373.217,373.345,373.353],'cm^-1')),
        HinderedRotor(inertia=(0.00121034,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0818712,'amu*angstrom^2'), symmetry=1, barrier=(8.08613,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.895355,0.0651017,-8.207e-05,5.02022e-08,-1.17613e-11,37798.3,29.6767], Tmin=(100,'K'), Tmax=(1059.28,'K')), NASAPolynomial(coeffs=[15.6619,0.00934104,-3.10955e-06,5.07724e-10,-3.29066e-14,34670,-42.4134], Tmin=(1059.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(313.321,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C2CsJOH) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([O])[C]([O])C[O](11908)',
    structure = SMILES('[CH]C([O])=C([O])C[O]'),
    E0 = (259.844,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,424.984,425.016,425.084,425.247,425.33,425.375,425.63],'cm^-1')),
        HinderedRotor(inertia=(0.000932563,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000932987,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22034,0.0649225,-8.89751e-05,7.01009e-08,-2.21883e-11,31348.5,27.7913], Tmin=(100,'K'), Tmax=(861.958,'K')), NASAPolynomial(coeffs=[8.29741,0.0261677,-1.12433e-05,2.02197e-09,-1.34592e-13,30348.2,-4.02601], Tmin=(861.958,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C([O])[C]1OC1[O](11909)',
    structure = SMILES('[CH2]C([O])=C1OC1[O]'),
    E0 = (89.2047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.663102,0.0602833,-6.59584e-05,3.39403e-08,-6.4997e-12,10860.2,23.7822], Tmin=(100,'K'), Tmax=(1454.55,'K')), NASAPolynomial(coeffs=[18.1012,0.00427596,1.03205e-07,-1.44e-10,1.27022e-14,6639.13,-63.9525], Tmin=(1454.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(89.2047,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[O]C=C([O])[C]1CO1(11910)',
    structure = SMILES('[O]C=C([O])[C]1CO1'),
    E0 = (36.2561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0512711,0.0640064,-7.17636e-05,3.80925e-08,-7.18671e-12,4528.17,28.8667], Tmin=(100,'K'), Tmax=(1630.94,'K')), NASAPolynomial(coeffs=[15.7243,0.00343598,4.0671e-06,-1.12985e-09,8.63905e-14,2292.31,-46.037], Tmin=(1630.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(36.2561,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C2CsJO)"""),
)

species(
    label = 'C=C([O])C1([O])[CH]O1(11911)',
    structure = SMILES('C=C([O])C1([O])[CH]O1'),
    E0 = (127.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.773624,0.0769545,-9.63689e-05,5.4214e-08,-1.06727e-11,15504.9,28.3935], Tmin=(100,'K'), Tmax=(1571.19,'K')), NASAPolynomial(coeffs=[19.447,-0.00362551,8.3427e-06,-2.00417e-09,1.48115e-13,12742.9,-66.8649], Tmin=(1571.19,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.285,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=CC(C)(O)OJ) + radical(C=C(C)OJ) + radical(CCsJO)"""),
)

species(
    label = '[O]C=C1OC[C]1[O](11912)',
    structure = SMILES('[O][CH]C1=C([O])CO1'),
    E0 = (76.1387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77032,0.0323794,2.32331e-05,-6.42191e-08,2.97517e-11,9252.94,22.8241], Tmin=(100,'K'), Tmax=(940.15,'K')), NASAPolynomial(coeffs=[18.1369,0.00488467,-1.32176e-07,2.4673e-11,-9.27611e-15,4313.22,-65.0296], Tmin=(940.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.1387,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJO)"""),
)

species(
    label = 'C=C1OC([O])[C]1[O](11913)',
    structure = SMILES('[CH2]C1=C([O])C([O])O1'),
    E0 = (63.9733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32954,0.0483141,-2.67532e-05,-1.05691e-08,1.04118e-11,7799.88,22.158], Tmin=(100,'K'), Tmax=(945.599,'K')), NASAPolynomial(coeffs=[17.208,0.00675876,-1.46305e-06,2.45065e-10,-2.03035e-14,3651.88,-59.6128], Tmin=(945.599,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.9733,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[O][C]1COC=C1[O](11914)',
    structure = SMILES('[O]C1=C([O])CO[CH]1'),
    E0 = (-122.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99428,0.0330866,6.45245e-06,-3.75194e-08,1.8393e-11,-14609.4,21.4416], Tmin=(100,'K'), Tmax=(940.475,'K')), NASAPolynomial(coeffs=[13.6171,0.011343,-3.03191e-06,5.09943e-10,-3.79558e-14,-18020.1,-40.4288], Tmin=(940.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.154,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(25dihydrofuran) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH2]C1([O])OC1=C[O](11888)',
    structure = SMILES('[CH2]C1([O])OC1=C[O]'),
    E0 = (148.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406415,0.0522865,1.4391e-05,-9.4867e-08,5.17732e-11,17985,21.4316], Tmin=(100,'K'), Tmax=(893.028,'K')), NASAPolynomial(coeffs=[33.4042,-0.0218396,1.51476e-05,-3.04839e-09,2.06515e-13,9153.6,-150.479], Tmin=(893.028,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(148.249,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(methyleneoxirane) + radical(C=CC(C)(O)OJ) + radical(C=COJ) + radical(C=CC(O)2CJ)"""),
)

species(
    label = 'C=C1OC1([O])[CH][O](11862)',
    structure = SMILES('C=C1OC1([O])[CH][O]'),
    E0 = (248.225,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.061779,0.0732378,-8.96651e-05,4.93876e-08,-9.898e-12,30014.7,23.329], Tmin=(100,'K'), Tmax=(1431.64,'K')), NASAPolynomial(coeffs=[21.8322,-0.00186746,3.6255e-06,-8.53109e-10,6.24399e-14,25173.7,-85.1656], Tmin=(1431.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=CC(C)(O)OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C1([O])OC=C1[O](11915)',
    structure = SMILES('[CH2]C1([O])OC=C1[O]'),
    E0 = (116.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.96397,0.0851106,-0.000102046,5.21516e-08,-9.27794e-12,14247.3,29.6643], Tmin=(100,'K'), Tmax=(1741.81,'K')), NASAPolynomial(coeffs=[24.4318,-0.0101288,9.78766e-06,-2.0639e-09,1.41483e-13,10304.1,-97.2512], Tmin=(1741.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(116.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=CC(C)(O)OJ) + radical(C=C(C)OJ) + radical(CJC(C)OC)"""),
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
    label = 'C=C([O])C([O])=C=O(11916)',
    structure = SMILES('C=C([O])C([O])=C=O'),
    E0 = (-78.3262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,2120,512.5,787.5,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.707535,'amu*angstrom^2'), symmetry=1, barrier=(16.2676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.789992,0.0678257,-9.66012e-05,6.46483e-08,-1.6152e-11,-9302.4,22.3854], Tmin=(100,'K'), Tmax=(1096.26,'K')), NASAPolynomial(coeffs=[16.1399,0.00451325,2.29422e-08,-1.89469e-10,2.02143e-14,-12229,-51.0778], Tmin=(1096.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-78.3262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
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
    label = '[O][C]=C[O](9592)',
    structure = SMILES('[O][CH][C]=O'),
    E0 = (204.374,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,343.122],'cm^-1')),
        HinderedRotor(inertia=(0.528453,'amu*angstrom^2'), symmetry=1, barrier=(44.0956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.99411,0.0153541,4.13359e-06,-1.97935e-08,9.27968e-12,24622.8,13.8337], Tmin=(100,'K'), Tmax=(987.313,'K')), NASAPolynomial(coeffs=[10.2587,0.00223056,-7.0477e-07,2.03549e-10,-2.00526e-14,22393.5,-25.1462], Tmin=(987.313,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(204.374,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=OCOJ) + radical(OCJC=O) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=C([O])C([O])=C[O](11373)',
    structure = SMILES('[CH]=C([O])C([O])=C[O]'),
    E0 = (121.511,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07575,'amu*angstrom^2'), symmetry=1, barrier=(24.7337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0393348,0.0739031,-9.76344e-05,5.7183e-08,-1.21177e-11,14768.5,25.3598], Tmin=(100,'K'), Tmax=(1359.98,'K')), NASAPolynomial(coeffs=[21.4734,-0.00406537,4.82431e-06,-1.11282e-09,8.20992e-14,10318.8,-79.563], Tmin=(1359.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C([O])C([O])=[C][O](11917)',
    structure = SMILES('C=C([O])C([O])=[C][O]'),
    E0 = (114.159,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.665457,'amu*angstrom^2'), symmetry=1, barrier=(15.3002,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918865,0.0648376,-8.09219e-05,4.13582e-08,-5.4014e-12,13844.1,24.9466], Tmin=(100,'K'), Tmax=(836.85,'K')), NASAPolynomial(coeffs=[17.3741,0.00266447,9.79688e-07,-3.55236e-10,2.99791e-14,10512.9,-54.9575], Tmin=(836.85,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.159,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'C=C([O])C([O])=[C]O(11918)',
    structure = SMILES('C=C([O])C([O])=[C]O'),
    E0 = (-27.3034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.763427,'amu*angstrom^2'), symmetry=1, barrier=(17.5527,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.764288,'amu*angstrom^2'), symmetry=1, barrier=(17.5725,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0547841,0.0774424,-0.000105547,6.51105e-08,-1.46528e-11,-3133.6,26.6183], Tmin=(100,'K'), Tmax=(1266.47,'K')), NASAPolynomial(coeffs=[20.5801,-0.00109063,3.70055e-06,-9.42686e-10,7.28742e-14,-7233.33,-72.9135], Tmin=(1266.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.3034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(C=CJO)"""),
)

species(
    label = 'C=C([O])C(O)=[C][O](11919)',
    structure = SMILES('[CH2]C([O])=C(O)[C]=O'),
    E0 = (-81.1371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.825019,'amu*angstrom^2'), symmetry=1, barrier=(18.9688,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.824382,'amu*angstrom^2'), symmetry=1, barrier=(18.9542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.824596,'amu*angstrom^2'), symmetry=1, barrier=(18.9591,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.127943,0.0830439,-0.000118886,7.68095e-08,-1.82681e-11,-9617.09,24.5818], Tmin=(100,'K'), Tmax=(896.579,'K')), NASAPolynomial(coeffs=[20.0712,0.00314447,-3.95645e-07,-6.11102e-12,2.95138e-15,-13558,-71.4898], Tmin=(896.579,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.1371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=C(O)CJ) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]=C(O)C([O])=C[O](11920)',
    structure = SMILES('[CH]=C(O)C([O])=C[O]'),
    E0 = (-16.2937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12851,'amu*angstrom^2'), symmetry=1, barrier=(25.9466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.12196,'amu*angstrom^2'), symmetry=1, barrier=(25.796,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.678846,0.0819667,-0.000104156,5.79493e-08,-1.15584e-11,-1773.13,26.4401], Tmin=(100,'K'), Tmax=(1470.78,'K')), NASAPolynomial(coeffs=[24.5671,-0.00663092,6.53523e-06,-1.44088e-09,1.03293e-13,-7042.95,-97.766], Tmin=(1470.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-16.2937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=C(O)C([O])=[C][O](11921)',
    structure = SMILES('C=C(O)C([O])=[C][O]'),
    E0 = (-23.6456,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.916093,'amu*angstrom^2'), symmetry=1, barrier=(21.0628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.919587,'amu*angstrom^2'), symmetry=1, barrier=(21.1431,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0375262,0.0751165,-9.66342e-05,5.63083e-08,-1.19979e-11,-2690.8,26.5933], Tmin=(100,'K'), Tmax=(1330.93,'K')), NASAPolynomial(coeffs=[21.1834,-0.000978214,3.26264e-06,-8.10908e-10,6.13327e-14,-7208.68,-77.2954], Tmin=(1330.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-23.6456,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([O])C(O)=C[O](11922)',
    structure = SMILES('[CH]C([O])=C(O)C=O'),
    E0 = (-29.9908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,2782.5,750,1395,475,1775,1000,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.03103,'amu*angstrom^2'), symmetry=1, barrier=(46.6974,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.02897,'amu*angstrom^2'), symmetry=1, barrier=(46.65,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.03089,'amu*angstrom^2'), symmetry=1, barrier=(46.6942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.437178,0.071766,-8.08292e-05,4.39992e-08,-9.22386e-12,-3473.13,24.3658], Tmin=(100,'K'), Tmax=(1177.97,'K')), NASAPolynomial(coeffs=[18.0935,0.0118106,-4.48326e-06,7.91417e-10,-5.38755e-14,-7632.85,-63.7079], Tmin=(1177.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-29.9908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])C([O])=CO(11923)',
    structure = SMILES('[CH]=C([O])C([O])=CO'),
    E0 = (-19.9515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.995919,'amu*angstrom^2'), symmetry=1, barrier=(22.8981,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.990715,'amu*angstrom^2'), symmetry=1, barrier=(22.7785,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.608192,0.0837215,-0.000111326,6.48223e-08,-1.35212e-11,-2218.38,26.2693], Tmin=(100,'K'), Tmax=(1405.59,'K')), NASAPolynomial(coeffs=[24.1215,-0.00692974,7.05163e-06,-1.58673e-09,1.15746e-13,-7167.39,-94.3314], Tmin=(1405.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-19.9515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = 'CO(2039)',
    structure = SMILES('[C-]#[O+]'),
    E0 = (-119.219,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2084.51],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    collisionModel = TransportData(shapeIndex=1, epsilon=(815.652,'J/mol'), sigma=(3.65,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(1.95,'angstroms^3'), rotrelaxcollnum=1.8, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.5971,-0.00102424,2.83336e-06,-1.75825e-09,3.42587e-13,-14343.2,3.45822], Tmin=(100,'K'), Tmax=(1669.93,'K')), NASAPolynomial(coeffs=[2.92796,0.00181931,-8.35309e-07,1.51269e-10,-9.88873e-15,-14292.7,6.51158], Tmin=(1669.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-119.219,'kJ/mol'), Cp0=(29.1007,'J/(mol*K)'), CpInf=(37.4151,'J/(mol*K)'), label="""CO""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = 'C=C([O])[CH][O](2850)',
    structure = SMILES('[CH2]C([O])=C[O]'),
    E0 = (20.9566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,217.215,217.577],'cm^-1')),
        HinderedRotor(inertia=(0.665078,'amu*angstrom^2'), symmetry=1, barrier=(22.287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00575,0.032418,-6.0508e-07,-3.82059e-08,2.20969e-11,2603.17,17.8437], Tmin=(100,'K'), Tmax=(887.044,'K')), NASAPolynomial(coeffs=[16.9262,-0.00266727,4.27963e-06,-9.58521e-10,6.70145e-14,-1310.55,-59.4906], Tmin=(887.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.9566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=C(O)CJ)"""),
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
    label = '[O][C]=C([O])C=O(11924)',
    structure = SMILES('[O]C=C([O])[C]=O'),
    E0 = (49.1287,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.18801,'amu*angstrom^2'), symmetry=1, barrier=(27.3148,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0382,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741901,0.061846,-8.35049e-05,4.84644e-08,-1.02655e-11,6034.62,20.3696], Tmin=(100,'K'), Tmax=(1305.56,'K')), NASAPolynomial(coeffs=[20.3265,-0.00611517,3.72058e-06,-7.446e-10,5.14083e-14,1599.04,-76.7393], Tmin=(1305.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(49.1287,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=OCOJ) + radical(C=COJ) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]C(C=O)=C1CO1(11925)',
    structure = SMILES('O=CC(=O)[C]1CO1'),
    E0 = (-111.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.46088,0.0487263,-4.72903e-05,2.49751e-08,-5.11733e-12,-13316.2,23.902], Tmin=(100,'K'), Tmax=(1347.78,'K')), NASAPolynomial(coeffs=[11.1438,0.0146918,-3.51654e-06,4.06794e-10,-1.92431e-14,-15445.2,-23.9177], Tmin=(1347.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.528,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsOsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + ring(Ethylene_oxide) + radical(C2CsJO)"""),
)

species(
    label = '[O]C=C1OCC1=O(11808)',
    structure = SMILES('[O]C=C1OCC1=O'),
    E0 = (-198.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62764,0.0268973,5.88805e-05,-1.12923e-07,4.94062e-11,-23817.5,18.3616], Tmin=(100,'K'), Tmax=(938.61,'K')), NASAPolynomial(coeffs=[23.732,-0.0029222,3.64726e-06,-6.14366e-10,2.81291e-14,-30802.9,-101.986], Tmin=(938.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-198.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(C=COJ)"""),
)

species(
    label = '[CH2][C]([O])C([O])[C]=O(11926)',
    structure = SMILES('[CH2][C]([O])C([O])[C]=O'),
    E0 = (452.816,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3000,3100,440,815,1455,1000,1855,455,950,483.119,483.138,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0797829,'amu*angstrom^2'), symmetry=1, barrier=(1.83437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.075417,'amu*angstrom^2'), symmetry=1, barrier=(12.492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0754188,'amu*angstrom^2'), symmetry=1, barrier=(12.4921,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04246,0.068671,-0.000101294,7.39899e-08,-2.04966e-11,54564.4,30.898], Tmin=(100,'K'), Tmax=(743.116,'K')), NASAPolynomial(coeffs=[12.2217,0.0151633,-6.74512e-06,1.2407e-09,-8.38146e-14,52718.8,-20.9546], Tmin=(743.116,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(452.816,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O][CH]C1([O])CC1=O(11827)',
    structure = SMILES('[O][CH]C1([O])CC1=O'),
    E0 = (286.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15241,0.070336,-0.000119576,1.05788e-07,-3.55655e-11,34513.1,24.3311], Tmin=(100,'K'), Tmax=(888.436,'K')), NASAPolynomial(coeffs=[7.75515,0.0231882,-1.05617e-05,1.91612e-09,-1.26198e-13,34027.4,-2.87306], Tmin=(888.436,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CC(C)(C=O)OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([O])=C1[CH]OO1(11927)',
    structure = SMILES('[CH2]C([O])=C1[CH]OO1'),
    E0 = (254.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.89196,0.0277792,3.88659e-05,-8.0234e-08,3.51493e-11,30677.9,24.3928], Tmin=(100,'K'), Tmax=(944.79,'K')), NASAPolynomial(coeffs=[18.1771,0.00572768,-5.78137e-07,1.35359e-10,-1.88467e-14,25507.7,-64.3256], Tmin=(944.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(254.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutane) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C1=C([O])[CH]OO1(11928)',
    structure = SMILES('[CH2]C1=C([O])[CH]OO1'),
    E0 = (163.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80796,0.0323991,2.0851e-05,-6.07111e-08,2.83959e-11,19736,24.2212], Tmin=(100,'K'), Tmax=(936.974,'K')), NASAPolynomial(coeffs=[17.5203,0.00545644,-2.67675e-07,3.03684e-11,-8.53963e-15,15029.9,-59.9602], Tmin=(936.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(163.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(12dioxolene) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(C=C(O)CJ)"""),
)

species(
    label = '[O]C1=C([O])C([O])C1(11929)',
    structure = SMILES('[O]C1=C([O])C([O])C1'),
    E0 = (48.5141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.583,0.0367646,1.33949e-05,-5.91345e-08,2.97822e-11,5937.18,23.4749], Tmin=(100,'K'), Tmax=(918.773,'K')), NASAPolynomial(coeffs=[19.7539,0.00050405,2.63882e-06,-5.69732e-10,3.46905e-14,789.662,-72.4919], Tmin=(918.773,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(48.5141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'O=[C]C(=O)C=O(11207)',
    structure = SMILES('O=[C]C(=O)C=O'),
    E0 = (-204.087,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(2.13548,'amu*angstrom^2'), symmetry=1, barrier=(49.0989,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13748,'amu*angstrom^2'), symmetry=1, barrier=(49.1448,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0382,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44847,0.0371926,-4.54499e-05,3.15863e-08,-9.19947e-12,-24492.8,15.9205], Tmin=(100,'K'), Tmax=(823.126,'K')), NASAPolynomial(coeffs=[6.61565,0.0169418,-8.54579e-06,1.69643e-09,-1.2117e-13,-25178.8,-3.37251], Tmin=(823.126,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-204.087,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)H) + radical(C=OC=OCJ=O)"""),
)

species(
    label = 'C=C=C([O])C=O(11930)',
    structure = SMILES('C=[C]C(=O)C=O'),
    E0 = (61.4013,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,413.184],'cm^-1')),
        HinderedRotor(inertia=(0.569432,'amu*angstrom^2'), symmetry=1, barrier=(13.0924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03685,'amu*angstrom^2'), symmetry=1, barrier=(23.8392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34268,0.0408215,-4.11136e-05,2.37985e-08,-6.04675e-12,7440.87,19.0438], Tmin=(100,'K'), Tmax=(910.792,'K')), NASAPolynomial(coeffs=[6.38375,0.0230738,-1.18843e-05,2.40342e-09,-1.74035e-13,6704.76,-0.074362], Tmin=(910.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.4013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + radical(C=CJC=O)"""),
)

species(
    label = 'HCO(1372)',
    structure = SMILES('[CH]=O'),
    E0 = (32.5964,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1116.11,1837.47,2716.03],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4140.62,'J/mol'), sigma=(3.59,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.06159,-0.0018195,9.37895e-06,-8.24742e-09,2.33851e-12,3918.59,3.98389], Tmin=(100,'K'), Tmax=(1105.12,'K')), NASAPolynomial(coeffs=[3.05335,0.00410464,-1.74962e-06,3.28532e-10,-2.28985e-14,4002.52,8.32038], Tmin=(1105.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(32.5964,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HCO""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(=O)[C]=O(2708)',
    structure = SMILES('[CH2]C([O])=C=O'),
    E0 = (61.1925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2120,512.5,787.5,268.278],'cm^-1')),
        HinderedRotor(inertia=(0.248151,'amu*angstrom^2'), symmetry=1, barrier=(12.6709,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85946,0.0408908,-5.12141e-05,3.03128e-08,-6.5709e-12,7442.22,17.0853], Tmin=(100,'K'), Tmax=(1331.22,'K')), NASAPolynomial(coeffs=[12.1004,0.00223239,1.23234e-06,-4.02519e-10,3.31529e-14,5414.45,-32.6265], Tmin=(1331.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.1925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)HHH) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
)

species(
    label = '[CH2]C([O])=[C][O](11931)',
    structure = SMILES('[CH2]C([O])=[C][O]'),
    E0 = (260.701,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,302.894,304.202],'cm^-1')),
        HinderedRotor(inertia=(0.00184505,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.80405,0.0391264,-4.53395e-05,2.46802e-08,-4.9033e-12,31442,21.6898], Tmin=(100,'K'), Tmax=(1465.05,'K')), NASAPolynomial(coeffs=[12.45,0.00181111,1.31153e-06,-3.91401e-10,3.07946e-14,29207.9,-30.7153], Tmin=(1465.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=C(O)CJ) + radical(C=CJO)"""),
)

species(
    label = 'CC([O])=C([O])[C]=O(11932)',
    structure = SMILES('CC([O])=C([O])[C]=O'),
    E0 = (3.67933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.703635,'amu*angstrom^2'), symmetry=1, barrier=(16.178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.702839,'amu*angstrom^2'), symmetry=1, barrier=(16.1596,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.473863,0.0787321,-0.000116411,8.16723e-08,-2.19065e-11,568.484,23.4335], Tmin=(100,'K'), Tmax=(925.831,'K')), NASAPolynomial(coeffs=[16.5313,0.0093554,-4.00645e-06,7.30767e-10,-4.94773e-14,-2404.73,-52.7963], Tmin=(925.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(3.67933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + radical(C=C(C)OJ) + radical(C=OCOJ) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]C=C([O])C[C]=O(11162)',
    structure = SMILES('[O]C=C([O])C[C]=O'),
    E0 = (-77.262,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1855,455,950,338.167,338.779,339.1],'cm^-1')),
        HinderedRotor(inertia=(0.217869,'amu*angstrom^2'), symmetry=1, barrier=(17.9781,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.218111,'amu*angstrom^2'), symmetry=1, barrier=(17.9671,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4491.42,'J/mol'), sigma=(6.92459,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=701.55 K, Pc=30.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24419,0.0477059,-2.08267e-05,-1.8493e-08,1.30808e-11,-9181.61,25.2029], Tmin=(100,'K'), Tmax=(971.866,'K')), NASAPolynomial(coeffs=[19.0515,0.00477417,-1.42205e-06,3.38529e-10,-3.1611e-14,-14076.6,-67.575], Tmin=(971.866,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-77.262,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C1=COCC1=O(11933)',
    structure = SMILES('O=C1[CH]OCC1=O'),
    E0 = (-192.109,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93711,0.0293701,2.75421e-05,-6.37864e-08,2.80896e-11,-23016.6,19.5341], Tmin=(100,'K'), Tmax=(957.138,'K')), NASAPolynomial(coeffs=[16.4047,0.00869168,-2.39902e-06,4.94705e-10,-4.30968e-14,-27608.4,-59.1497], Tmin=(957.138,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-192.109,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + ring(Cyclopentane) + radical(CCsJOCs)"""),
)

species(
    label = 'CC(=O)C([O])=C=O(11934)',
    structure = SMILES('CC(=O)C(=O)[C]=O'),
    E0 = (-267.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85678,0.0510343,-5.68839e-05,3.66054e-08,-9.98738e-12,-32068.1,19.8767], Tmin=(100,'K'), Tmax=(871.388,'K')), NASAPolynomial(coeffs=[7.51769,0.025049,-1.21536e-05,2.38442e-09,-1.69549e-13,-33054.7,-6.65471], Tmin=(871.388,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-267.243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)H) + radical(C=OC=OCJ=O)"""),
)

species(
    label = '[CH2][C](O)C([O])=[C][O](11935)',
    structure = SMILES('[CH2][C](O)C([O])=[C][O]'),
    E0 = (294.549,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,371.25,371.253,371.267],'cm^-1')),
        HinderedRotor(inertia=(0.0012232,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810738,'amu*angstrom^2'), symmetry=1, barrier=(7.92955,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810786,'amu*angstrom^2'), symmetry=1, barrier=(7.92961,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.658072,0.0733109,-0.000104048,6.92134e-08,-1.71133e-11,35546.7,32.0511], Tmin=(100,'K'), Tmax=(854.2,'K')), NASAPolynomial(coeffs=[16.2731,0.00781705,-2.43203e-06,3.59893e-10,-2.11514e-14,32600.7,-42.4502], Tmin=(854.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(294.549,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C2CsJOH) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]1OOC1=C[O](11936)',
    structure = SMILES('[CH2][C]1OOC1=C[O]'),
    E0 = (381.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72815,0.0335818,2.00302e-05,-5.74984e-08,2.58124e-11,46017.1,26.1248], Tmin=(100,'K'), Tmax=(971.693,'K')), NASAPolynomial(coeffs=[17.5611,0.00843288,-2.93807e-06,6.53581e-10,-5.62661e-14,41050.5,-59.5288], Tmin=(971.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(C=COJ) + radical(C2CsJOO) + radical(CJCOOH)"""),
)

species(
    label = 'C=C([O])C(=O)C=O(11155)',
    structure = SMILES('[CH2]C(=O)C(=O)C=O'),
    E0 = (-168.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2814,0.0658337,-0.000106877,9.29977e-08,-3.16583e-11,-20210.6,24.6991], Tmin=(100,'K'), Tmax=(835.489,'K')), NASAPolynomial(coeffs=[8.22035,0.0225392,-1.10623e-05,2.11271e-09,-1.45043e-13,-21018.5,-5.42605], Tmin=(835.489,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-168.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)(Cds-O2d)) + group(Cds-O2d(Cds-O2d)H) + radical(CJCC=O)"""),
)

species(
    label = 'C=C1OC1([O])C=O(11881)',
    structure = SMILES('C=C1OC1([O])C=O'),
    E0 = (-68.2729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03246,0.0395483,3.45965e-05,-1.03109e-07,5.10501e-11,-8080.33,19.1065], Tmin=(100,'K'), Tmax=(906.005,'K')), NASAPolynomial(coeffs=[28.8709,-0.0147722,1.0979e-05,-2.17595e-09,1.43027e-13,-15939.6,-127.985], Tmin=(906.005,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-68.2729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C([O])C([O])[C]=O(11145)',
    structure = SMILES('C=C([O])C([O])[C]=O'),
    E0 = (56.6024,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,225.492,226.025,4000],'cm^-1')),
        HinderedRotor(inertia=(0.3979,'amu*angstrom^2'), symmetry=1, barrier=(14.3261,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00332241,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3606,0.0539288,-5.88187e-05,3.05704e-08,-5.71278e-12,6906.41,28.8859], Tmin=(100,'K'), Tmax=(976.711,'K')), NASAPolynomial(coeffs=[13.8314,0.0103358,-3.35687e-06,5.54759e-10,-3.68304e-14,4113.57,-32.811], Tmin=(976.711,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(56.6024,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=C(C)OJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C([O])C([O])C=O(11937)',
    structure = SMILES('[CH]=C([O])C([O])C=O'),
    E0 = (143.738,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,436.159,4000],'cm^-1')),
        HinderedRotor(inertia=(0.101483,'amu*angstrom^2'), symmetry=1, barrier=(13.7723,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10171,'amu*angstrom^2'), symmetry=1, barrier=(13.7704,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36036,0.0544636,-6.02669e-05,3.3318e-08,-7.16455e-12,17385.8,28.472], Tmin=(100,'K'), Tmax=(1143.61,'K')), NASAPolynomial(coeffs=[13.523,0.0119223,-4.46785e-06,7.89805e-10,-5.36512e-14,14604,-31.8376], Tmin=(1143.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(143.738,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[O]C1(C=O)CC1=O(11840)',
    structure = SMILES('[O]C1(C=O)CC1=O'),
    E0 = (-49.6726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78025,0.0443951,-3.81929e-05,1.71885e-08,-3.11048e-12,-5890.65,23.0823], Tmin=(100,'K'), Tmax=(1324.42,'K')), NASAPolynomial(coeffs=[11.1217,0.0161818,-6.23876e-06,1.10361e-09,-7.42119e-14,-8365,-24.6092], Tmin=(1324.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.6726,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsCs) + group(Cds-OdCsH) + ring(cyclopropanone) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = '[CH][O](1548)',
    structure = SMILES('[CH][O]'),
    E0 = (431.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([806.798,806.798,2696.84],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (29.018,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.86882,-0.000695042,1.53508e-05,-2.00576e-08,7.65929e-12,51865.1,6.76618], Tmin=(100,'K'), Tmax=(973.285,'K')), NASAPolynomial(coeffs=[6.03251,-0.000301449,4.32946e-07,-3.67045e-11,-1.25373e-15,51004.1,-5.87326], Tmin=(973.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(H3COJ) + radical(OsCsJ2H_triplet)"""),
)

species(
    label = '[O]C1CC(=O)C1=O(11164)',
    structure = SMILES('[O]C1CC(=O)C1=O'),
    E0 = (-120.299,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.19252,0.0292548,1.03486e-05,-3.33706e-08,1.41007e-11,-14394.4,22.1464], Tmin=(100,'K'), Tmax=(1023.96,'K')), NASAPolynomial(coeffs=[11.6785,0.0167161,-7.19901e-06,1.43756e-09,-1.06809e-13,-17622.4,-30.1189], Tmin=(1023.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-120.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + ring(Cyclobutane) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C1([O])C(=O)C1[O](11938)',
    structure = SMILES('[CH2]C1([O])C(=O)C1[O]'),
    E0 = (325.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.458934,0.0636209,-7.27071e-05,3.87043e-08,-7.58072e-12,39282.8,28.2226], Tmin=(100,'K'), Tmax=(1454.23,'K')), NASAPolynomial(coeffs=[18.668,0.00263771,1.43578e-06,-4.38122e-10,3.42003e-14,35139.1,-62.4826], Tmin=(1454.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(CC(C)(C=O)OJ) + radical(C=OCOJ) + radical(CJC(C)(C=O)O)"""),
)

species(
    label = 'C=C1OC([O])C1=O(11939)',
    structure = SMILES('C=C1OC([O])C1=O'),
    E0 = (-102.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91481,0.0301855,2.29196e-05,-5.73038e-08,2.50655e-11,-12241.3,19.3167], Tmin=(100,'K'), Tmax=(976.727,'K')), NASAPolynomial(coeffs=[16.4426,0.0089771,-3.30936e-06,7.32487e-10,-6.18116e-14,-16905.5,-59.7785], Tmin=(976.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-102.521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=OCOJ)"""),
)

species(
    label = '[CH]=C([O])[C](O)[CH][O](11940)',
    structure = SMILES('[CH]C([O])=C(O)[CH][O]'),
    E0 = (239.335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,504.103,508.393,510.449,510.53,512.232,512.348],'cm^-1')),
        HinderedRotor(inertia=(0.279096,'amu*angstrom^2'), symmetry=1, barrier=(50.9871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273864,'amu*angstrom^2'), symmetry=1, barrier=(51.145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275241,'amu*angstrom^2'), symmetry=1, barrier=(50.8666,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.477473,0.0646316,-6.47052e-05,3.17092e-08,-5.92705e-12,28923.1,29.5591], Tmin=(100,'K'), Tmax=(1440.48,'K')), NASAPolynomial(coeffs=[17.8316,0.0108396,-2.85677e-06,3.85302e-10,-2.21214e-14,24504.7,-58.4807], Tmin=(1440.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(239.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(CCOJ) + radical(C=CCJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])C(=O)C[O](11941)',
    structure = SMILES('[CH]C(=O)C(=O)C[O]'),
    E0 = (197.64,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,365,385,505,600,445,480,1700,1720,461.958,658.057,660.204,2679.98,2681.41,2681.74],'cm^-1')),
        HinderedRotor(inertia=(1.41816,'amu*angstrom^2'), symmetry=1, barrier=(32.6062,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.51081,'amu*angstrom^2'), symmetry=1, barrier=(11.7445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.105573,'amu*angstrom^2'), symmetry=1, barrier=(32.6037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44316,0.0606351,-8.32952e-05,6.3102e-08,-1.94553e-11,23858.7,25.2602], Tmin=(100,'K'), Tmax=(788.97,'K')), NASAPolynomial(coeffs=[8.92066,0.0227238,-1.12154e-05,2.19391e-09,-1.54806e-13,22678.8,-9.04182], Tmin=(788.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(197.64,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)Cs) + radical(C=OCOJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C([O])([C]=O)C=O(11942)',
    structure = SMILES('[CH2]C([O])([C]=O)C=O'),
    E0 = (113.15,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1855,455,950,334.284,334.553,334.575,334.622,334.676,1707.72],'cm^-1')),
        HinderedRotor(inertia=(0.114914,'amu*angstrom^2'), symmetry=1, barrier=(9.10466,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114315,'amu*angstrom^2'), symmetry=1, barrier=(9.0888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114597,'amu*angstrom^2'), symmetry=1, barrier=(9.10146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01286,0.0726532,-0.000122843,1.07155e-07,-3.59956e-11,13709.7,28.3826], Tmin=(100,'K'), Tmax=(858.197,'K')), NASAPolynomial(coeffs=[9.12694,0.0216225,-1.05575e-05,1.99102e-09,-1.34881e-13,12803.5,-6.68783], Tmin=(858.197,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(113.15,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CC(C)(C=O)OJ) + radical(CJC(C)(C=O)O) + radical(CC(C)(O)CJ=O)"""),
)

species(
    label = '[CH2][C][O](2821)',
    structure = SMILES('[CH2][C][O]'),
    E0 = (648.642,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,2091.34],'cm^-1')),
        HinderedRotor(inertia=(0.0328816,'amu*angstrom^2'), symmetry=1, barrier=(10.2243,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (42.0367,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.26468,0.0187769,-3.12134e-05,3.2382e-08,-1.33361e-11,78037.6,11.2589], Tmin=(100,'K'), Tmax=(748.797,'K')), NASAPolynomial(coeffs=[3.86755,0.0107985,-5.69995e-06,1.18131e-09,-8.60475e-14,78080.7,9.41548], Tmin=(748.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.642,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C([O])C(=O)[C]=O(11943)',
    structure = SMILES('[CH2]C([O])C(=O)[C]=O'),
    E0 = (134.086,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,1855,455,950,380.342,380.342],'cm^-1')),
        HinderedRotor(inertia=(0.00116534,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00116535,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.108031,'amu*angstrom^2'), symmetry=1, barrier=(11.0898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4427,0.0574125,-6.87401e-05,4.20657e-08,-1.02085e-11,16218,28.0562], Tmin=(100,'K'), Tmax=(1004.68,'K')), NASAPolynomial(coeffs=[11.7392,0.0164187,-7.53614e-06,1.45331e-09,-1.0275e-13,14149.1,-21.6664], Tmin=(1004.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(134.086,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(C=OCOJ) + radical(CJCO) + radical(CCCJ=O)"""),
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
    E0 = (-125.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (668.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (668.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (708.195,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (111.203,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (114.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (18.8076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (-62.1849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (-100.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (411.682,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (520.95,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (338.294,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (284.818,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (89.2047,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (60.9508,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (127.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (76.1387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (63.9733,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (-10.7189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (148.249,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (248.41,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (116.34,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (171.661,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (131.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (189.578,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (365.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (333.303,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (325.951,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (135.476,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (112.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (175.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (20.6629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (14.3178,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (13.0887,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (207.044,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (486.931,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (-111.528,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (-117.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (475.677,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (286.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (254.298,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (163.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (70.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (216.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (319.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (187.72,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (293.297,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (130.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (168.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (-118.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (-100.612,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (319.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (381.805,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (-125.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (-68.2729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (214.801,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (335.867,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (-49.6726,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (547.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (-117.301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (325.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (-102.521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (264.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (241.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (235.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (628.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (263.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['OCHCO(3676)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH2]C([O])=[C]C=O(9910)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', 'C=[C]C([O])=C[O](11900)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=C([O])C(=C)[O](2705)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C1OOC1=C[O](11845)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(236.788,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination
Ea raised from 233.6 to 236.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C([O])C1=COO1(11901)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(240.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 237.9 to 240.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C1OOC=C1[O](11902)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(144.393,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 139.4 to 144.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C([O])C(O)=C=O(11903)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C(O)C([O])=C=O(11904)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C([O])C([O])=[C][O](11905)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=C([O])C([O])[CH][O](11906)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[C]([O])C([O])=[C][O](11907)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=C([O])[C]([O])C[O](11908)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C([O])[C]1OC1[O](11909)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(214.79,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 213.2 to 214.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O]C=C([O])[C]1CO1(11910)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(186.536,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe_2H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C([O])C1([O])[CH]O1(11911)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(252.87,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 252.5 to 252.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O]C=C1OC[C]1[O](11912)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(201.724,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 197.4 to 201.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C1OC([O])[C]1[O](11913)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(189.558,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 187.2 to 189.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O][C]1COC=C1[O](11914)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(8.9814e+11,'s^-1'), n=0.0209575, Ea=(114.866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[CH2]C1([O])OC1=C[O](11888)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(273.834,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 272.9 to 273.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C1OC1([O])[CH][O](11862)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(373.995,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_(Cd)_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[CH2]C1([O])OC=C1[O](11915)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(241.925,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 241.7 to 241.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', 'C=C([O])C([O])=C=O(11916)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['OCHCO(3676)', '[CH2][C]=O(1376)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][C]=C[O](9592)', 'ketene(1375)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O][C]=C[O](9592)', '[CH2][C]=O(1376)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH]=C([O])C([O])=C[O](11373)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', 'C=C([O])C([O])=[C][O](11917)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([O])C([O])=[C]O(11918)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=C([O])C(O)=[C][O](11919)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C(O)C([O])=C[O](11920)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=C(O)C([O])=[C][O](11921)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out;XH_out] for rate rule [R4H_DSS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=C([O])C(O)=C[O](11922)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(37100,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=C([O])C([O])=CO(11923)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSMS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['CO(2039)', 'C=C([O])[CH][O](2850)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CH2(T)(20)', '[O][C]=C([O])C=O(11924)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O]C(C=O)=C1CO1(11925)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(14.0568,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 2.23606797749979
family: Birad_recombination
Ea raised from 11.2 to 14.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O]C=C1OCC1=O(11808)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Cpri_rad_out_2H] for rate rule [R4_SDS;O_rad;Cpri_rad_out_2H]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][C]([O])C([O])[C]=O(11926)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O][CH]C1([O])CC1=O(11827)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.68393e+12,'s^-1'), n=-0.105173, Ea=(411.965,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 135 used for R3_D;doublebond_intra;radadd_intra_cs2H
Exact match found for rate rule [R3_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[CH2]C([O])=C1[CH]OO1(11927)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(379.883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 375.0 to 379.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[CH2]C1=C([O])[CH]OO1(11928)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(7.96106e+11,'s^-1'), n=0.00276955, Ea=(288.903,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra] for rate rule [R5_SD_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Endocyclic
Ea raised from 284.5 to 288.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O]C1=C([O])C([O])C1(11929)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(196.02,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD;multiplebond_intra;radadd_intra_cs] for rate rule [R5_SD_CO;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.449489742783178
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['CH2(T)(20)', 'O=[C]C(=O)C=O(11207)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Ck_O;YJ] for rate rule [Ck_O;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction45',
    reactants = ['O(4)', 'C=C=C([O])C=O(11930)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.70446,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;O_atom_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction46',
    reactants = ['HCO(1372)', '[CH2]C(=O)[C]=O(2708)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd_R;CO_pri_rad] for rate rule [Ck_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction47',
    reactants = ['HCO(1372)', '[CH2]C([O])=[C][O](11931)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.81e+06,'m^3/(mol*s)'), n=-5.80997e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R_N-2R->C]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['CC([O])=C([O])[C]=O(11932)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.63218e+07,'s^-1'), n=1.50075, Ea=(126.696,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_SDS;CO_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C=C([O])C[C]=O(11162)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O]C1=COCC1=O(11933)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSDS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['CC(=O)C([O])=C=O(11934)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2][C](O)C([O])=[C][O](11935)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[CH2][C]1OOC1=C[O](11936)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(507.39,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 503.3 to 507.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction54',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C([O])C(=O)C=O(11155)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C1OC1([O])C=O(11881)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(57.3121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination
Ea raised from 54.4 to 57.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction56',
    reactants = ['C=C([O])C([O])[C]=O(11145)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;XH_out] for rate rule [R2H_S;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH]=C([O])C([O])C=O(11937)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O]C1(C=O)CC1=O(11840)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(75.9124,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination
Ea raised from 72.7 to 75.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH][O](1548)', '[CH2]C(=O)[C]=O(2708)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction60',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[O]C1CC(=O)C1=O(11164)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction61',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['[CH2]C1([O])C(=O)C1[O](11938)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(2.63856e+09,'s^-1'), n=0.755479, Ea=(451.039,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 449.9 to 451.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction62',
    reactants = ['C=C([O])C([O])=C[O](11158)'],
    products = ['C=C1OC([O])C1=O(11939)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(23.0639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination
Ea raised from 18.3 to 23.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH]=C([O])[C](O)[CH][O](11940)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH]=C([O])C(=O)C[O](11941)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction65',
    reactants = ['[CH2]C([O])([C]=O)C=O(11942)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CO]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction66',
    reactants = ['OCHCO(3676)', '[CH2][C][O](2821)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH2]C([O])C(=O)[C]=O(11943)'],
    products = ['C=C([O])C([O])=C[O](11158)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(5945.77,'s^-1'), n=2.58622, Ea=(129.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;XH_out] for rate rule [R3H_SS;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2826',
    isomers = [
        'C=C([O])C([O])=C[O](11158)',
    ],
    reactants = [
        ('OCHCO(3676)', 'ketene(1375)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2826',
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


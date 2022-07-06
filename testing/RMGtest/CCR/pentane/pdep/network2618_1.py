species(
    label = 'C[CH]CC([O])=C[O](10963)',
    structure = SMILES('C[CH]CC([O])=C[O]'),
    E0 = (10.0605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3010,987.5,1337.5,450,1655,333.024,333.025,333.028,3866.16],'cm^-1')),
        HinderedRotor(inertia=(0.00152003,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266246,'amu*angstrom^2'), symmetry=1, barrier=(20.9538,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266247,'amu*angstrom^2'), symmetry=1, barrier=(20.9538,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.895905,0.0593956,-3.58801e-05,-1.87598e-09,7.00189e-12,1329.82,29.0267], Tmin=(100,'K'), Tmax=(946.472,'K')), NASAPolynomial(coeffs=[15.4649,0.0192123,-6.09342e-06,1.01928e-09,-6.9469e-14,-2386,-45.5196], Tmin=(946.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.0605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJC)"""),
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
    label = 'C3H6(27)',
    structure = SMILES('C=CC'),
    E0 = (5.9763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650],'cm^-1')),
        HinderedRotor(inertia=(0.497558,'amu*angstrom^2'), symmetry=1, barrier=(11.4398,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2218.31,'J/mol'), sigma=(4.982,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.31912,0.00817959,3.34736e-05,-4.36193e-08,1.58213e-11,749.325,9.54025], Tmin=(100,'K'), Tmax=(983.754,'K')), NASAPolynomial(coeffs=[5.36755,0.0170743,-6.35107e-06,1.16619e-09,-8.2762e-14,-487.139,-4.54469], Tmin=(983.754,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(5.9763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), label="""C3H6""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH2]C(C)C([O])=C[O](10965)',
    structure = SMILES('[CH2]C(C)C([O])=C[O]'),
    E0 = (12.7246,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,373.192,373.294,373.321],'cm^-1')),
        HinderedRotor(inertia=(0.00120893,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177071,'amu*angstrom^2'), symmetry=1, barrier=(17.5168,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1771,'amu*angstrom^2'), symmetry=1, barrier=(17.5186,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4304.58,'J/mol'), sigma=(7.06223,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=672.37 K, Pc=27.73 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.745948,0.0578002,-1.60607e-05,-3.48842e-08,2.23667e-11,1660.5,28.3685], Tmin=(100,'K'), Tmax=(903.692,'K')), NASAPolynomial(coeffs=[19.7571,0.0111936,-1.01508e-06,-1.29987e-11,2.39075e-15,-3308.52,-69.9059], Tmin=(903.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.7246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl)"""),
)

species(
    label = 'CHCH3(T)(21)',
    structure = SMILES('[CH]C'),
    E0 = (343.893,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,592.414,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00438698,'amu*angstrom^2'), symmetry=1, barrier=(26.7685,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0532,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.82365,-0.000909765,3.21389e-05,-3.73492e-08,1.33096e-11,41371.4,7.10941], Tmin=(100,'K'), Tmax=(960.802,'K')), NASAPolynomial(coeffs=[4.3048,0.00943081,-3.27566e-06,5.95138e-10,-4.27321e-14,40709.2,1.84242], Tmin=(960.802,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.893,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""CHCH3(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH]CC([O])=C[O](11435)',
    structure = SMILES('[CH]CC([O])=C[O]'),
    E0 = (287.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,496.654,496.661,496.669,496.695,496.714,496.744],'cm^-1')),
        HinderedRotor(inertia=(0.00068336,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810607,'amu*angstrom^2'), symmetry=1, barrier=(14.1871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34366,0.0467992,-1.95324e-05,-2.26557e-08,1.64173e-11,34697.9,23.4452], Tmin=(100,'K'), Tmax=(911.817,'K')), NASAPolynomial(coeffs=[18.415,0.00306204,1.17073e-06,-3.23509e-10,2.11275e-14,30289.7,-64.4395], Tmin=(911.817,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(287.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ2_triplet)"""),
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
    label = 'C[CH]C[C]=C[O](9690)',
    structure = SMILES('C[CH]C[C]=C[O]'),
    E0 = (324.656,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,207.595,786.49,1689.13],'cm^-1')),
        HinderedRotor(inertia=(0.601732,'amu*angstrom^2'), symmetry=1, barrier=(13.835,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.602133,'amu*angstrom^2'), symmetry=1, barrier=(13.8442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0315643,'amu*angstrom^2'), symmetry=1, barrier=(13.8306,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44674,0.0498142,-3.24375e-05,7.62536e-09,5.33056e-13,39144.4,25.9661], Tmin=(100,'K'), Tmax=(1065.61,'K')), NASAPolynomial(coeffs=[11.2955,0.0228573,-8.58564e-06,1.5206e-09,-1.03372e-13,36477,-24.8415], Tmin=(1065.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.656,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])C[CH]C(2453)',
    structure = SMILES('[CH]=C([O])C[CH]C'),
    E0 = (324.487,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3120,650,792.5,1650,667.004,2355.02],'cm^-1')),
        HinderedRotor(inertia=(0.00331348,'amu*angstrom^2'), symmetry=1, barrier=(13.0477,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.56849,'amu*angstrom^2'), symmetry=1, barrier=(13.0707,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.567954,'amu*angstrom^2'), symmetry=1, barrier=(13.0584,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.1,'J/mol'), sigma=(6.27192,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.76 K, Pc=33.61 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57622,0.0525886,-4.59139e-05,2.30288e-08,-4.84485e-12,39114.6,25.7201], Tmin=(100,'K'), Tmax=(1122.91,'K')), NASAPolynomial(coeffs=[9.08313,0.0258472,-1.01918e-05,1.82038e-09,-1.23017e-13,37428.7,-11.3665], Tmin=(1122.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(324.487,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(RCCJC) + radical(Cds_P)"""),
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
    label = 'C[C]CC([O])=C[O](11548)',
    structure = SMILES('C[C]CC([O])=C[O]'),
    E0 = (263.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,359.504,359.548,359.555,359.558,359.599],'cm^-1')),
        HinderedRotor(inertia=(0.192219,'amu*angstrom^2'), symmetry=1, barrier=(17.6352,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192196,'amu*angstrom^2'), symmetry=1, barrier=(17.635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.192222,'amu*angstrom^2'), symmetry=1, barrier=(17.6346,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.619815,0.0625268,-3.68645e-05,-9.74203e-09,1.20344e-11,31863.9,27.1882], Tmin=(100,'K'), Tmax=(930.572,'K')), NASAPolynomial(coeffs=[19.6248,0.010976,-2.35302e-06,3.39629e-10,-2.47358e-14,27021.8,-70.1442], Tmin=(930.572,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC1CC(=C[O])O1(11475)',
    structure = SMILES('CC1CC(=C[O])O1'),
    E0 = (-158.991,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11295,0.0358787,6.86296e-05,-1.37524e-07,6.26982e-11,-18992.3,19.1428], Tmin=(100,'K'), Tmax=(902.81,'K')), NASAPolynomial(coeffs=[25.7081,0.000218757,6.07242e-06,-1.38442e-09,9.21888e-14,-26420.9,-113.546], Tmin=(902.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-158.991,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]CC1=COO1(11549)',
    structure = SMILES('C[CH]CC1=COO1'),
    E0 = (250.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.51399,0.0444847,-5.13593e-06,-2.17019e-08,1.0359e-11,30183.4,27.1104], Tmin=(100,'K'), Tmax=(1044.3,'K')), NASAPolynomial(coeffs=[11.9662,0.0264196,-1.07452e-05,2.0247e-09,-1.43844e-13,26802.4,-29.5044], Tmin=(1044.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.142,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(RCCJC)"""),
)

species(
    label = 'CC1CC([O])=CO1(11550)',
    structure = SMILES('CC1CC([O])=CO1'),
    E0 = (-241.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48492,0.0303028,7.14157e-05,-1.33365e-07,5.99038e-11,-28955.8,18.1343], Tmin=(100,'K'), Tmax=(897.92,'K')), NASAPolynomial(coeffs=[22.4433,0.00358993,4.6973e-06,-1.16255e-09,7.93711e-14,-35406.5,-95.6831], Tmin=(897.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-241.701,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC=CC(O)=C[O](11551)',
    structure = SMILES('CC=CC(O)=C[O]'),
    E0 = (-223.327,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.334349,0.0644106,-2.37607e-05,-3.34205e-08,2.28864e-11,-26713,23.4714], Tmin=(100,'K'), Tmax=(917.759,'K')), NASAPolynomial(coeffs=[23.3431,0.00688116,3.90531e-07,-2.05884e-10,1.19903e-14,-32736.8,-95.3669], Tmin=(917.759,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-223.327,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ)"""),
)

species(
    label = 'C[CH]CC(O)=C=O(11552)',
    structure = SMILES('C[CH]CC(O)=C=O'),
    E0 = (-91.6221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.631919,0.0706078,-7.39201e-05,4.1151e-08,-9.09444e-12,-10895.3,26.6719], Tmin=(100,'K'), Tmax=(1104.64,'K')), NASAPolynomial(coeffs=[14.0383,0.0220615,-7.99773e-06,1.36526e-09,-9.00752e-14,-13857.1,-39.34], Tmin=(1104.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-91.6221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + radical(RCCJC)"""),
)

species(
    label = 'CC=CC([O])=CO(11553)',
    structure = SMILES('CC=CC([O])=CO'),
    E0 = (-226.985,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.970022,0.0820064,-8.47898e-05,4.14627e-08,-7.42667e-12,-27097.5,28.2611], Tmin=(100,'K'), Tmax=(1633.01,'K')), NASAPolynomial(coeffs=[22.4696,0.00643365,1.30677e-06,-4.94969e-10,3.91252e-14,-32331.8,-88.9035], Tmin=(1633.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-226.985,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CCCC([O])=C=O(11554)',
    structure = SMILES('CCCC(=O)[C]=O'),
    E0 = (-175.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21921,0.0661706,-7.55151e-05,5.44606e-08,-1.70505e-11,-20968.3,24.8824], Tmin=(100,'K'), Tmax=(760.209,'K')), NASAPolynomial(coeffs=[6.83678,0.0366138,-1.71978e-05,3.32121e-09,-2.33642e-13,-21822.4,-0.679191], Tmin=(760.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-175.135,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CCC(O)=C[O](10973)',
    structure = SMILES('C=CCC(O)=C[O]'),
    E0 = (-193.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.617335,0.0554899,5.54094e-07,-5.74116e-08,3.09261e-11,-23089.9,26.2923], Tmin=(100,'K'), Tmax=(924.982,'K')), NASAPolynomial(coeffs=[23.2285,0.00662382,4.76678e-07,-1.86234e-10,7.94715e-15,-29365.4,-92.3414], Tmin=(924.982,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.139,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = 'C=CCC([O])=CO(11555)',
    structure = SMILES('C=CCC([O])=CO'),
    E0 = (-196.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.564704,0.0585677,-1.06257e-05,-4.624e-08,2.7535e-11,-23529.6,26.5726], Tmin=(100,'K'), Tmax=(911.792,'K')), NASAPolynomial(coeffs=[22.8386,0.00613578,1.13492e-06,-3.70521e-10,2.3865e-14,-29473.8,-89.1514], Tmin=(911.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-196.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C[CH][CH]C([O])[CH][O](11556)',
    structure = SMILES('C[CH][CH]C([O])[CH][O]'),
    E0 = (530.33,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3025,3050,390,407.5,425,1340,1350,1360,335,352.5,370,200.527,294.191,1249.76,1406.03,3303.17],'cm^-1')),
        HinderedRotor(inertia=(0.0823966,'amu*angstrom^2'), symmetry=1, barrier=(2.27063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0823966,'amu*angstrom^2'), symmetry=1, barrier=(2.27063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0823966,'amu*angstrom^2'), symmetry=1, barrier=(2.27063,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0823966,'amu*angstrom^2'), symmetry=1, barrier=(2.27063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.960449,0.0763923,-0.000118578,1.10356e-07,-4.01528e-11,63884.2,33.2687], Tmin=(100,'K'), Tmax=(848.052,'K')), NASAPolynomial(coeffs=[4.13825,0.0408115,-1.92223e-05,3.61857e-09,-2.46772e-13,64085.7,22.8272], Tmin=(848.052,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.33,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCJCO) + radical(RCCJC) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH]C[C]1OC1[O](11557)',
    structure = SMILES('C[CH]C[C]1OC1[O]'),
    E0 = (249.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39708,0.0632152,-8.09814e-05,6.8347e-08,-2.27856e-11,30040,27.7666], Tmin=(100,'K'), Tmax=(929.495,'K')), NASAPolynomial(coeffs=[3.63754,0.0378114,-1.45485e-05,2.45485e-09,-1.56041e-13,30304.4,20.7842], Tmin=(929.495,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(RCCJC)"""),
)

species(
    label = 'C[CH]CC1([O])[CH]O1(11558)',
    structure = SMILES('C[CH]CC1([O])[CH]O1'),
    E0 = (247.001,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.112305,0.0694226,-6.54063e-05,3.27105e-08,-6.2097e-12,29861.1,29.234], Tmin=(100,'K'), Tmax=(1512.87,'K')), NASAPolynomial(coeffs=[15.1261,0.018154,-3.0997e-06,1.98078e-10,-1.51528e-15,26642.7,-45.0376], Tmin=(1512.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(247.001,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(RCCJC) + radical(CCsJO)"""),
)

species(
    label = 'CC1C[C]([O])C1[O](11559)',
    structure = SMILES('CC1C[C]([O])C1[O]'),
    E0 = (273.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27994,0.0466979,1.4171e-06,-3.63345e-08,1.73654e-11,32999.2,23.9631], Tmin=(100,'K'), Tmax=(989.408,'K')), NASAPolynomial(coeffs=[14.8038,0.0226229,-8.47484e-06,1.58915e-09,-1.15264e-13,28825.3,-48.7068], Tmin=(989.408,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.46,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + ring(Cyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH)"""),
)

species(
    label = 'CC1CC1([O])[CH][O](11504)',
    structure = SMILES('CC1CC1([O])[CH][O]'),
    E0 = (270.108,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05427,0.0618366,-5.11892e-05,2.23812e-08,-4.01278e-12,32594.9,24.3031], Tmin=(100,'K'), Tmax=(1312.33,'K')), NASAPolynomial(coeffs=[12.5366,0.0268385,-1.11863e-05,2.05971e-09,-1.41529e-13,29581.1,-34.2135], Tmin=(1312.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.108,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsOs) + group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Cyclopropane) + radical(CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH)"""),
)

species(
    label = 'CC=CC([O])=C[O](11560)',
    structure = SMILES('CC=CC([O])=C[O]'),
    E0 = (-85.5225,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.929737,'amu*angstrom^2'), symmetry=1, barrier=(21.3765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.933283,'amu*angstrom^2'), symmetry=1, barrier=(21.458,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.695856,0.0605033,-3.15808e-05,-1.57935e-08,1.45464e-11,-10155.7,23.6744], Tmin=(100,'K'), Tmax=(920.727,'K')), NASAPolynomial(coeffs=[19.6156,0.0102215,-1.65535e-06,1.83502e-10,-1.32488e-14,-14992.4,-73.3855], Tmin=(920.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-85.5225,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2][CH]CC(=O)C=O(11041)',
    structure = SMILES('[CH2][CH]CC(=O)C=O'),
    E0 = (70.0532,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,320.033,320.548],'cm^-1')),
        HinderedRotor(inertia=(0.00164131,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00164035,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00163947,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00164242,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6094,0.0576313,-5.8052e-05,3.67953e-08,-1.04584e-11,8507.11,26.952], Tmin=(100,'K'), Tmax=(818.521,'K')), NASAPolynomial(coeffs=[6.19767,0.0352092,-1.69624e-05,3.32903e-09,-2.36991e-13,7755.98,5.73504], Tmin=(818.521,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(70.0532,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCJCC=O) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CC([O])=C=O(11561)',
    structure = SMILES('C[CH]CC(=O)[C]=O'),
    E0 = (24.7674,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,375,552.5,462.5,1710,1855,455,950,324.673,3247.65],'cm^-1')),
        HinderedRotor(inertia=(0.1362,'amu*angstrom^2'), symmetry=1, barrier=(10.1804,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.131067,'amu*angstrom^2'), symmetry=1, barrier=(10.1531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00154737,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.472565,'amu*angstrom^2'), symmetry=1, barrier=(34.9788,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4514,0.0604789,-6.5022e-05,4.21597e-08,-1.17711e-11,3066.74,26.5477], Tmin=(100,'K'), Tmax=(848.222,'K')), NASAPolynomial(coeffs=[7.43776,0.0322485,-1.50988e-05,2.92186e-09,-2.06231e-13,2051.19,-1.34757], Tmin=(848.222,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(24.7674,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCJCC=O) + radical(CCCJ=O)"""),
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
    label = 'C3H6(T)(28)',
    structure = SMILES('[CH2][CH]C'),
    E0 = (284.865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000],'cm^-1')),
        HinderedRotor(inertia=(0.238388,'amu*angstrom^2'), symmetry=1, barrier=(5.48101,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0090964,'amu*angstrom^2'), symmetry=1, barrier=(22.1004,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (42.0797,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.93779,0.019099,4.26879e-06,-1.44878e-08,5.7496e-12,34303.2,12.9695], Tmin=(100,'K'), Tmax=(1046.8,'K')), NASAPolynomial(coeffs=[5.93905,0.0171893,-6.69156e-06,1.21547e-09,-8.39803e-14,33151.2,-4.14862], Tmin=(1046.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(284.865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), label="""C3H6(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C[CH][CH]C([O])=C[O](11562)',
    structure = SMILES('C[CH]C=C([O])[CH][O]'),
    E0 = (201.676,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,3010,987.5,1337.5,450,1655,350,440,435,1725,518.621,518.64,518.651,518.651],'cm^-1')),
        HinderedRotor(inertia=(0.000626679,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000626707,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000626756,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4196,0.0480737,-2.03705e-05,-6.49106e-09,5.4144e-12,24356.4,26.8131], Tmin=(100,'K'), Tmax=(1056.09,'K')), NASAPolynomial(coeffs=[12.3712,0.0234867,-9.44247e-06,1.75693e-09,-1.23608e-13,21101.1,-31.08], Tmin=(1056.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(201.676,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(Allyl_S) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][CH]CC([O])=C[O](11045)',
    structure = SMILES('[CH2][CH]CC([O])=C[O]'),
    E0 = (215.307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,388.649,388.677,388.945,388.978],'cm^-1')),
        HinderedRotor(inertia=(0.00111423,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111527,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00111552,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.533841,0.0647548,-6.21241e-05,3.03168e-08,-5.7138e-12,26030,32.1176], Tmin=(100,'K'), Tmax=(1404.45,'K')), NASAPolynomial(coeffs=[16.6344,0.0146285,-4.02632e-06,5.73747e-10,-3.39807e-14,21928.7,-49.5272], Tmin=(1404.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CC([O])=[C][O](11563)',
    structure = SMILES('C[CH]CC([O])=[C][O]'),
    E0 = (249.805,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,1685,370,320.14,320.142,2631.54,2631.56],'cm^-1')),
        HinderedRotor(inertia=(0.00164486,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244851,'amu*angstrom^2'), symmetry=1, barrier=(17.8075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.244851,'amu*angstrom^2'), symmetry=1, barrier=(17.8075,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.188,0.0601625,-5.91826e-05,3.20654e-08,-7.04898e-12,30147.3,31.1083], Tmin=(100,'K'), Tmax=(1097.56,'K')), NASAPolynomial(coeffs=[11.2626,0.023446,-9.00328e-06,1.58596e-09,-1.06439e-13,27935.8,-18.4338], Tmin=(1097.56,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(249.805,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJC) + radical(C=CJO)"""),
)

species(
    label = 'CC[CH]C([O])=C[O](11489)',
    structure = SMILES('CC[CH]C([O])=C[O]'),
    E0 = (15.5162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3010,987.5,1337.5,450,1655,448.958,448.967,449.685,450.369],'cm^-1')),
        HinderedRotor(inertia=(0.00083872,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0971749,'amu*angstrom^2'), symmetry=1, barrier=(13.8032,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0963356,'amu*angstrom^2'), symmetry=1, barrier=(13.8491,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.856231,0.0543835,-9.16952e-06,-3.83248e-08,2.20071e-11,1992.96,28.9305], Tmin=(100,'K'), Tmax=(931.397,'K')), NASAPolynomial(coeffs=[19.3241,0.0127049,-2.65536e-06,3.9431e-10,-2.995e-14,-3079.58,-67.6166], Tmin=(931.397,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(15.5162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJCO)"""),
)

species(
    label = '[CH2]CCC([O])=C[O](11564)',
    structure = SMILES('[CH2]CCC([O])=C[O]'),
    E0 = (20.8605,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,421.021,421.205,421.225,421.246],'cm^-1')),
        HinderedRotor(inertia=(0.000961598,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114949,'amu*angstrom^2'), symmetry=1, barrier=(14.3696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.115483,'amu*angstrom^2'), symmetry=1, barrier=(14.3816,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.70423,0.0602707,-2.80521e-05,-1.68042e-08,1.38009e-11,2638.8,28.7659], Tmin=(100,'K'), Tmax=(938.965,'K')), NASAPolynomial(coeffs=[18.553,0.0147068,-3.94254e-06,6.40128e-10,-4.58017e-14,-2056.35,-63.3728], Tmin=(938.965,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.8605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CC([O])=[C]O(11565)',
    structure = SMILES('C[CH]CC([O])=[C]O'),
    E0 = (108.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,1685,370,336.133,336.229,3788.56],'cm^-1')),
        HinderedRotor(inertia=(0.00149175,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177778,'amu*angstrom^2'), symmetry=1, barrier=(14.2461,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177486,'amu*angstrom^2'), symmetry=1, barrier=(14.2463,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.177445,'amu*angstrom^2'), symmetry=1, barrier=(14.2466,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.706493,0.0679694,-6.55778e-05,2.98522e-08,-4.07964e-12,13153.2,31.4258], Tmin=(100,'K'), Tmax=(921.466,'K')), NASAPolynomial(coeffs=[14.5871,0.0196379,-6.31029e-06,1.01485e-09,-6.54922e-14,10088.9,-37.1515], Tmin=(921.466,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(108.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(RCCJC) + radical(C=CJO)"""),
)

species(
    label = 'C[CH][CH]C(O)=C[O](11566)',
    structure = SMILES('C[CH]C=C(O)[CH][O]'),
    E0 = (63.8712,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,350,440,435,1725,3010,987.5,1337.5,450,1655,511.455,511.46,511.467],'cm^-1')),
        HinderedRotor(inertia=(0.000644507,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000644462,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15576,'amu*angstrom^2'), symmetry=1, barrier=(28.9119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.155757,'amu*angstrom^2'), symmetry=1, barrier=(28.9115,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08559,0.0516989,-1.1795e-05,-2.46974e-08,1.37961e-11,7797.88,26.509], Tmin=(100,'K'), Tmax=(989.6,'K')), NASAPolynomial(coeffs=[15.7461,0.0207267,-7.72346e-06,1.44342e-09,-1.0458e-13,3511.23,-51.0642], Tmin=(989.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(63.8712,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(Allyl_S) + radical(C=CCJO)"""),
)

species(
    label = 'C[CH]CC(O)=[C][O](11567)',
    structure = SMILES('C[CH]C[C](O)[C]=O'),
    E0 = (110.269,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,360,370,350,1855,455,950,256.683,2835.38],'cm^-1')),
        HinderedRotor(inertia=(0.00255744,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.19373,'amu*angstrom^2'), symmetry=1, barrier=(9.07333,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.194236,'amu*angstrom^2'), symmetry=1, barrier=(9.07681,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00255775,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00255724,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11148,0.0674108,-8.14782e-05,5.93077e-08,-1.80395e-11,13362.8,32.0478], Tmin=(100,'K'), Tmax=(793.635,'K')), NASAPolynomial(coeffs=[8.19653,0.0317012,-1.3985e-05,2.61183e-09,-1.7981e-13,12238.2,-0.495783], Tmin=(793.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.269,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(C2CsJOH) + radical(CCCJ=O)"""),
)

species(
    label = 'CCCC([O])=[C][O](11568)',
    structure = SMILES('CCCC([O])=[C][O]'),
    E0 = (55.3584,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,1685,370,384.927,385.606,386.035,386.316],'cm^-1')),
        HinderedRotor(inertia=(0.00113349,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0807777,'amu*angstrom^2'), symmetry=1, barrier=(8.52346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0808249,'amu*angstrom^2'), symmetry=1, barrier=(8.52809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.843822,0.0616853,-4.58412e-05,1.14605e-08,1.29705e-12,6778.54,29.6067], Tmin=(100,'K'), Tmax=(992.653,'K')), NASAPolynomial(coeffs=[14.8302,0.0207704,-7.3533e-06,1.28633e-09,-8.81555e-14,3240.87,-41.5993], Tmin=(992.653,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.3584,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = 'C[CH][CH]C([O])=CO(11569)',
    structure = SMILES('C[CH]C=C([O])[CH]O'),
    E0 = (-24.0291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12091,0.0501339,-6.03642e-06,-3.3564e-08,1.80145e-11,-2774.44,27.0213], Tmin=(100,'K'), Tmax=(958.332,'K')), NASAPolynomial(coeffs=[16.4374,0.0181085,-5.8468e-06,1.04302e-09,-7.57691e-14,-7175.13,-53.8635], Tmin=(958.332,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-24.0291,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(Allyl_S) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][CH]CC(O)=C[O](11570)',
    structure = SMILES('[CH2][CH]CC(O)=C[O]'),
    E0 = (77.502,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,480.401,482.251,486.56],'cm^-1')),
        HinderedRotor(inertia=(0.000708592,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0907973,'amu*angstrom^2'), symmetry=1, barrier=(15.1764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0895621,'amu*angstrom^2'), symmetry=1, barrier=(15.1359,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0920769,'amu*angstrom^2'), symmetry=1, barrier=(15.1853,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.585057,0.0638754,-3.79399e-05,-7.9953e-09,1.12206e-11,9454.65,30.4286], Tmin=(100,'K'), Tmax=(928.944,'K')), NASAPolynomial(coeffs=[18.9172,0.0136948,-3.3467e-06,5.03237e-10,-3.49441e-14,4807.96,-63.3408], Tmin=(928.944,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(77.502,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = '[CH2][CH]CC([O])=CO(11571)',
    structure = SMILES('[CH2][CH]CC([O])=CO'),
    E0 = (73.8442,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,323.462,323.475,3996.17],'cm^-1')),
        HinderedRotor(inertia=(0.00161138,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264984,'amu*angstrom^2'), symmetry=1, barrier=(19.6728,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265009,'amu*angstrom^2'), symmetry=1, barrier=(19.6718,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265134,'amu*angstrom^2'), symmetry=1, barrier=(19.6727,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.135308,0.0747968,-7.64658e-05,3.86356e-08,-7.34668e-12,9044.12,33.107], Tmin=(100,'K'), Tmax=(1470.01,'K')), NASAPolynomial(coeffs=[19.4164,0.0115526,-1.68394e-06,7.38795e-11,1.74373e-15,4380.94,-65.0608], Tmin=(1470.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.8442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(RCCJC) + radical(RCCJ)"""),
)

species(
    label = 'C[CH]CC(=O)C=O(10960)',
    structure = SMILES('C[CH]CC(=O)C=O'),
    E0 = (-135.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71562,0.0553488,-4.27763e-05,1.8916e-08,-3.79728e-12,-16182,24.7745], Tmin=(100,'K'), Tmax=(1092.69,'K')), NASAPolynomial(coeffs=[6.78972,0.036774,-1.72774e-05,3.35863e-09,-2.37845e-13,-17290.9,-0.15486], Tmin=(1092.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCJCC=O)"""),
)

species(
    label = 'C[CH]CC([O])[C]=O(10953)',
    structure = SMILES('C[CH]CC([O])[C]=O'),
    E0 = (177.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,180,180,2611.7],'cm^-1')),
        HinderedRotor(inertia=(0.137837,'amu*angstrom^2'), symmetry=1, barrier=(3.16913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0459305,'amu*angstrom^2'), symmetry=1, barrier=(16.3437,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00337659,'amu*angstrom^2'), symmetry=1, barrier=(16.3436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00337661,'amu*angstrom^2'), symmetry=1, barrier=(16.3436,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.56728,0.0550808,-4.60221e-05,2.23791e-08,-4.70764e-12,21419.5,31.7645], Tmin=(100,'K'), Tmax=(1097.54,'K')), NASAPolynomial(coeffs=[8.1522,0.0310823,-1.32239e-05,2.45704e-09,-1.69815e-13,19974.1,-0.616895], Tmin=(1097.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(RCCJC) + radical(CCCJ=O)"""),
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
    label = 'C[CH]C[C]=O(2361)',
    structure = SMILES('C[CH]C[C]=O'),
    E0 = (132.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,204.852],'cm^-1')),
        HinderedRotor(inertia=(0.00401758,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.636302,'amu*angstrom^2'), symmetry=1, barrier=(18.9512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000722328,'amu*angstrom^2'), symmetry=1, barrier=(6.83812,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.35404,0.0356784,-2.29653e-05,7.91178e-09,-1.1666e-12,16026.6,20.8557], Tmin=(100,'K'), Tmax=(1491.52,'K')), NASAPolynomial(coeffs=[7.59578,0.0216208,-8.82771e-06,1.59264e-09,-1.07418e-13,14463,-6.52828], Tmin=(1491.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O) + radical(CCCJ=O)"""),
)

species(
    label = 'CC1CC(=O)C1[O](10968)',
    structure = SMILES('CC1CC(=O)C1[O]'),
    E0 = (-74.8997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.11912,0.0223279,6.58427e-05,-1.00498e-07,3.95829e-11,-8923.46,24.381], Tmin=(100,'K'), Tmax=(964.984,'K')), NASAPolynomial(coeffs=[13.7856,0.0217242,-7.45169e-06,1.4227e-09,-1.08026e-13,-13398.5,-43.0082], Tmin=(964.984,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.8997,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ)"""),
)

species(
    label = 'CC=CC(=O)C[O](11572)',
    structure = SMILES('CC=CC(=O)C[O]'),
    E0 = (-64.0398,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.10733,0.0498707,-1.35515e-05,-5.86252e-08,6.26414e-11,-7642.49,22.8991], Tmin=(100,'K'), Tmax=(462.336,'K')), NASAPolynomial(coeffs=[4.55793,0.040187,-1.95033e-05,3.84212e-09,-2.74032e-13,-7992.19,11.6356], Tmin=(462.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-64.0398,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CCC(=O)C[O](10974)',
    structure = SMILES('C=CCC(=O)C[O]'),
    E0 = (-49.5301,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15246,0.0548976,-3.59601e-05,1.07596e-08,-1.25588e-12,-5847.91,26.268], Tmin=(100,'K'), Tmax=(1999.75,'K')), NASAPolynomial(coeffs=[19.6072,0.0179837,-8.27135e-06,1.5289e-09,-1.01902e-13,-13228.9,-75.5551], Tmin=(1999.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.5301,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'CC1C[C]([CH][O])O1(11573)',
    structure = SMILES('CC1C[C]([CH][O])O1'),
    E0 = (257.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13468,0.0665944,-8.14253e-05,6.24979e-08,-1.91394e-11,31095.6,24.6143], Tmin=(100,'K'), Tmax=(956.07,'K')), NASAPolynomial(coeffs=[6.3629,0.0335697,-1.21173e-05,1.97044e-09,-1.22365e-13,30605.5,2.29128], Tmin=(956.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.714,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCOJ) + radical(C2CsJOCs) + radical(CCsJOH)"""),
)

species(
    label = 'C[CH]C[C]1[CH]OO1(11574)',
    structure = SMILES('C[CH]C[C]1[CH]OO1'),
    E0 = (468.822,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.79919,0.0403181,-1.69136e-05,2.39101e-09,-6.59615e-14,56355.4,15.7187], Tmin=(100,'K'), Tmax=(2844.92,'K')), NASAPolynomial(coeffs=[47.4256,-0.00771341,1.39476e-06,-2.55049e-10,2.20747e-14,26147.3,-249.831], Tmin=(2844.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(468.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(12dioxetane) + radical(C2CsJOO) + radical(RCCJC) + radical(CCsJOO)"""),
)

species(
    label = 'C[CH][CH]C(=O)C[O](11575)',
    structure = SMILES('C[CH]C=C([O])C[O]'),
    E0 = (84.3809,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,350,440,435,1725,3010,987.5,1337.5,450,1655,386.495,386.617,1637.2,1637.79],'cm^-1')),
        HinderedRotor(inertia=(0.00112674,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0741804,'amu*angstrom^2'), symmetry=1, barrier=(7.86542,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0743204,'amu*angstrom^2'), symmetry=1, barrier=(7.8571,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.7275,0.0530541,-3.91309e-05,1.65782e-08,-3.14768e-12,10227.8,25.1112], Tmin=(100,'K'), Tmax=(1162.09,'K')), NASAPolynomial(coeffs=[7.10234,0.0345533,-1.52503e-05,2.87824e-09,-2.0039e-13,8978.63,-1.62664], Tmin=(1162.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(84.3809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(Allyl_S)"""),
)

species(
    label = '[CH2][CH]CC(=O)C[O](11576)',
    structure = SMILES('[CH2][CH]CC(=O)C[O]'),
    E0 = (221.654,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,375,552.5,462.5,1710,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,282.23,1457.42,1457.43,3060.36],'cm^-1')),
        HinderedRotor(inertia=(0.00211632,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.178332,'amu*angstrom^2'), symmetry=1, barrier=(10.0824,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.17833,'amu*angstrom^2'), symmetry=1, barrier=(10.0822,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00211596,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2342,0.0653215,-7.56769e-05,5.47533e-08,-1.69766e-11,26754.3,30.1951], Tmin=(100,'K'), Tmax=(771.491,'K')), NASAPolynomial(coeffs=[7.17402,0.0345251,-1.58005e-05,3.0128e-09,-2.10282e-13,25837.7,3.07986], Tmin=(771.491,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(221.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + radical(C=OCOJ) + radical(CCJCC=O) + radical(RCCJ)"""),
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
    label = 'C[CH]C[C][O](11577)',
    structure = SMILES('C[CH]C[C][O]'),
    E0 = (577.008,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,180,180,1758.67,1758.7],'cm^-1')),
        HinderedRotor(inertia=(0.189007,'amu*angstrom^2'), symmetry=1, barrier=(4.34565,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189067,'amu*angstrom^2'), symmetry=1, barrier=(4.34703,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.18923,'amu*angstrom^2'), symmetry=1, barrier=(4.35076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07622,0.0499436,-7.68464e-05,7.63257e-08,-2.94778e-11,69459.9,21.539], Tmin=(100,'K'), Tmax=(839.333,'K')), NASAPolynomial(coeffs=[1.82334,0.0335417,-1.60679e-05,3.05749e-09,-2.10129e-13,70122.6,26.4093], Tmin=(839.333,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(577.008,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CH2_triplet)"""),
)

species(
    label = 'CC1CC1([O])C=O(11526)',
    structure = SMILES('CC1CC1([O])C=O'),
    E0 = (-35.3552,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34367,0.0488368,-1.52925e-05,-1.53596e-08,9.53382e-12,-4148.14,24.2088], Tmin=(100,'K'), Tmax=(994.509,'K')), NASAPolynomial(coeffs=[13.0222,0.0230913,-8.47642e-06,1.53272e-09,-1.07557e-13,-7520.71,-37.3464], Tmin=(994.509,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-35.3552,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(CC(C)(C=O)OJ)"""),
)

species(
    label = 'CCC=C([O])C=O(11578)',
    structure = SMILES('CC[CH]C(=O)C=O'),
    E0 = (-135.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71562,0.0553488,-4.27763e-05,1.8916e-08,-3.79728e-12,-16182,24.7745], Tmin=(100,'K'), Tmax=(1092.69,'K')), NASAPolynomial(coeffs=[6.78972,0.036774,-1.72774e-05,3.35863e-09,-2.37845e-13,-17290.9,-0.15486], Tmin=(1092.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-135.193,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCJCC=O)"""),
)

species(
    label = 'CC=CC([O])C=O(11579)',
    structure = SMILES('CC=CC([O])C=O'),
    E0 = (-62.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36376,0.0522304,-3.4614e-05,1.14453e-08,-1.53567e-12,-7433.17,28.3179], Tmin=(100,'K'), Tmax=(1708.45,'K')), NASAPolynomial(coeffs=[13.4564,0.023918,-9.75626e-06,1.74543e-09,-1.16295e-13,-11565.2,-36.4991], Tmin=(1708.45,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-62.63,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ)"""),
)

species(
    label = 'C=CCC([O])C=O(10972)',
    structure = SMILES('C=CCC([O])C=O'),
    E0 = (-49.2499,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3684,0.0486167,-1.7644e-05,-9.7805e-09,6.5101e-12,-5820.72,29.3967], Tmin=(100,'K'), Tmax=(1058.2,'K')), NASAPolynomial(coeffs=[12.545,0.0248019,-1.00153e-05,1.87465e-09,-1.32419e-13,-9218.17,-30.0326], Tmin=(1058.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.2499,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'CC1C[C]([O])[CH]O1(11580)',
    structure = SMILES('CC1C[C]([O])[CH]O1'),
    E0 = (179.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24209,0.0492278,3.37333e-07,-4.67032e-08,2.63399e-11,21757.9,20.3571], Tmin=(100,'K'), Tmax=(869.315,'K')), NASAPolynomial(coeffs=[15.6347,0.0168612,-2.2367e-06,7.43758e-11,2.74107e-15,17976.2,-54.4217], Tmin=(869.315,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(179.989,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOCs)"""),
)

species(
    label = 'C=C([O])C=O(2859)',
    structure = SMILES('[CH2]C(=O)C=O'),
    E0 = (-74.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.43166,'amu*angstrom^2'), symmetry=1, barrier=(9.92472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231641,'amu*angstrom^2'), symmetry=1, barrier=(31.0079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.689,0.0301409,-2.32871e-05,9.30138e-09,-1.56938e-12,-8891.4,16.525], Tmin=(100,'K'), Tmax=(1335.25,'K')), NASAPolynomial(coeffs=[7.34752,0.0161853,-7.60952e-06,1.4738e-09,-1.03803e-13,-10135.4,-7.29647], Tmin=(1335.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJCC=O)"""),
)

species(
    label = 'C[CH][CH]C([O])C=O(11581)',
    structure = SMILES('C[CH][CH]C([O])C=O'),
    E0 = (217.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,2782.5,750,1395,475,1775,1000,444.056,2370.57,2370.57],'cm^-1')),
        HinderedRotor(inertia=(0.0027732,'amu*angstrom^2'), symmetry=1, barrier=(31.4871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0157725,'amu*angstrom^2'), symmetry=1, barrier=(11.5234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.36948,'amu*angstrom^2'), symmetry=1, barrier=(31.4871,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0027732,'amu*angstrom^2'), symmetry=1, barrier=(31.4871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59886,0.0496762,-3.18747e-05,1.04957e-08,-1.43466e-12,26225.8,33.3268], Tmin=(100,'K'), Tmax=(1639.48,'K')), NASAPolynomial(coeffs=[11.0307,0.0266644,-1.08208e-05,1.93444e-09,-1.29186e-13,23133.1,-16.8394], Tmin=(1639.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(RCCJC)"""),
)

species(
    label = '[CH2][CH]CC([O])C=O(11582)',
    structure = SMILES('[CH2][CH]CC([O])C=O'),
    E0 = (222.66,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1568.99,3336.21,3336.27],'cm^-1')),
        HinderedRotor(inertia=(0.00211946,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.741616,'amu*angstrom^2'), symmetry=1, barrier=(41.8566,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165936,'amu*angstrom^2'), symmetry=1, barrier=(9.36714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0131107,'amu*angstrom^2'), symmetry=1, barrier=(22.9003,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77769,0.0517078,-3.75802e-05,1.54767e-08,-2.83895e-12,26857.4,31.9739], Tmin=(100,'K'), Tmax=(1200.33,'K')), NASAPolynomial(coeffs=[7.2712,0.0334013,-1.47035e-05,2.77101e-09,-1.92682e-13,25538.6,4.46771], Tmin=(1200.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(222.66,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(RCCJC) + radical(RCCJ)"""),
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
    E0 = (10.0605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (172.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (420.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (479.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (844.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (843.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (475.621,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (18.3448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (250.142,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (17.3825,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (73.4606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (73.4606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (44.5785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (18.4285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (35.0337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (35.0337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (553.191,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (249.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (247.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (273.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (270.538,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (147.427,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (284.356,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (274.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (222.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (255.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (489.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (413.468,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (427.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (461.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (206.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (172.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (271.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (218.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (303.648,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (99.6669,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (123.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (151.121,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (147.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (10.0605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (340.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (619.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (18.3448,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (73.4606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (35.0337,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (258.409,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (469.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (213.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (295.273,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (664.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (17.5917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (35.8759,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (32.9218,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (49.2855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (179.989,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (269.584,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (375.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (353.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['OCHCO(3676)', 'C3H6(27)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['CHCH3(T)(21)', 'C=C([O])[CH][O](2850)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CH3(17)', '[CH]CC([O])=C[O](11435)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', 'C[CH]C[C]=C[O](9690)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', '[CH]=C([O])C[CH]C(2453)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['H(3)', 'C[C]CC([O])=C[O](11548)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC1CC(=C[O])O1(11475)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C[CH]CC1=COO1(11549)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(240.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 237.9 to 240.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC1CC([O])=CO1(11550)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.24579e+11,'s^-1'), n=0.1555, Ea=(7.322,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;C_rad_out_1H;Ypri_rad_out] + [R5;C_rad_out_single;Ypri_rad_out] for rate rule [R5_SSDS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC=CC(O)=C[O](11551)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C[CH]CC(O)=C=O(11552)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC=CC([O])=CO(11553)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.552e+09,'s^-1'), n=0.311, Ea=(34.518,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad_NDe] for rate rule [R4radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CCCC([O])=C=O(11554)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.76e+08,'s^-1'), n=0.311, Ea=(8.368,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [R4;Y_rad_NDe;XH_Rrad] for rate rule [R4radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C=CCC(O)=C[O](10973)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C=CCC([O])=CO(11555)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;Y_rad;XH_Rrad] for rate rule [R6radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH][CH]C([O])[CH][O](11556)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C[CH]C[C]1OC1[O](11557)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(239.089,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction19',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C[CH]CC1([O])[CH]O1(11558)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(237.333,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC1C[C]([O])C1[O](11559)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.25757e+07,'s^-1'), n=1.165, Ea=(263.399,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_HH_D;doublebond_intra;radadd_intra_csHCs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 261.8 to 263.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC1CC1([O])[CH][O](11504)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.73316e+11,'s^-1'), n=0.231216, Ea=(260.478,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_csHNd]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', 'CC=CC([O])=C[O](11560)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(0.0272924,'m^3/(mol*s)'), n=2.81111, Ea=(21.1569,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 26 used for Cds-CdH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CdH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][CH]CC(=O)C=O(11041)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', 'C[CH]CC([O])=C=O(11561)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][C]=C[O](9592)', 'C3H6(27)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00620445,'m^3/(mol*s)'), n=2.46568, Ea=(12.4666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['OCHCO(3676)', 'C3H6(T)(28)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][C]=C[O](9592)', 'C3H6(T)(28)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', 'C[CH][CH]C([O])=C[O](11562)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', '[CH2][CH]CC([O])=C[O](11045)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', 'C[CH]CC([O])=[C][O](11563)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CC[CH]C([O])=C[O](11489)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.9e+10,'s^-1'), n=0.75, Ea=(190.79,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 159 used for R2H_S;C_rad_out_H/Cd;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2]CCC([O])=C[O](11564)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(718000,'s^-1'), n=2.05, Ea=(151.879,'kJ/mol'), T0=(1,'K'), Tmin=(500,'K'), Tmax=(2000,'K'), comment="""From training reaction 147 used for R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[CH]CC([O])=[C]O(11565)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C[CH][CH]C(O)=C[O](11566)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(52326.9,'s^-1'), n=2.1859, Ea=(154.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C[CH]CC(O)=[C][O](11567)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['CCCC([O])=[C][O](11568)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out;Cs_H_out_1H] for rate rule [R4H_DSS;Cd_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C[CH][CH]C([O])=CO(11569)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]CC(O)=C[O](11570)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2][CH]CC([O])=CO(11571)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1389.4,'s^-1'), n=2.18252, Ea=(73.6483,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;C_rad_out_2H;XH_out] for rate rule [R6HJ_1;C_rad_out_2H;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C[CH]CC(=O)C=O(10960)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH][O](1548)', 'C[CH]C[C]=O(2361)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC1CC(=O)C1[O](10968)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC=CC(=O)C[O](11572)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C=CCC(=O)C[O](10974)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC1C[C]([CH][O])O1(11573)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.7342e+08,'s^-1'), n=0.920995, Ea=(248.348,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_csHCs] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_csHCs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction47',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C[CH]C[C]1[CH]OO1(11574)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(459.093,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra_O] for rate rule [R4_linear;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C[CH][CH]C(=O)C[O](11575)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(11891.5,'s^-1'), n=2.58622, Ea=(129.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][CH]CC(=O)C[O](11576)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(8.47295e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['HCO(1372)', 'C[CH]C[C][O](11577)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC1CC1([O])C=O(11526)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2e+12,'s^-1'), n=0, Ea=(7.5312,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_1H;Ypri_rad_out] for rate rule [R3_SS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction52',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CCC=C([O])C=O(11578)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.64935e+10,'s^-1'), n=0.2551, Ea=(25.8154,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad_De]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC=CC([O])C=O(11579)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 0 used for R2radExo;Y_rad_De;XH_Rrad_NDe
Exact match found for rate rule [R2radExo;Y_rad_De;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['C=CCC([O])C=O(10972)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(1.34494e+09,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radEndo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C[CH]CC([O])=C[O](10963)'],
    products = ['CC1C[C]([O])[CH]O1(11580)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(4.64e+06,'s^-1'), n=1.15, Ea=(169.929,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_csHCs] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_csHCs]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 169.0 to 169.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction56',
    reactants = ['CHCH3(T)(21)', 'C=C([O])C=O(2859)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction57',
    reactants = ['C[CH][CH]C([O])C=O(11581)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH2][CH]CC([O])C=O(11582)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(434.916,'s^-1'), n=2.67444, Ea=(131.125,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_2H;XH_out] for rate rule [R4HJ_1;C_rad_out_2H;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2618',
    isomers = [
        'C[CH]CC([O])=C[O](10963)',
    ],
    reactants = [
        ('OCHCO(3676)', 'C3H6(27)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #2618',
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


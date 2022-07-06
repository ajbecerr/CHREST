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
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.745948,0.0578002,-1.60607e-05,-3.48842e-08,2.23667e-11,1660.5,28.3685], Tmin=(100,'K'), Tmax=(903.692,'K')), NASAPolynomial(coeffs=[19.7571,0.0111936,-1.01508e-06,-1.29987e-11,2.39075e-15,-3308.52,-69.9059], Tmin=(903.692,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(12.7246,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl)"""),
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
    label = '[CH2]CC([O])=C[O](10944)',
    structure = SMILES('[CH2]CC([O])=C[O]'),
    E0 = (44.6408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,497.907,497.908,497.909],'cm^-1')),
        HinderedRotor(inertia=(0.0814769,'amu*angstrom^2'), symmetry=1, barrier=(14.3337,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.081477,'amu*angstrom^2'), symmetry=1, barrier=(14.3337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(4147.05,'J/mol'), sigma=(6.71796,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=647.76 K, Pc=31.04 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43016,0.0445211,-1.06591e-05,-2.97655e-08,1.81862e-11,5472.72,23.9167], Tmin=(100,'K'), Tmax=(919.676,'K')), NASAPolynomial(coeffs=[17.3218,0.0068291,-4.39556e-07,-1.81256e-11,-3.4228e-16,1220.68,-58.6458], Tmin=(919.676,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(44.6408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJ)"""),
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
    label = 'C[CH]CC([O])=C[O](10963)',
    structure = SMILES('C[CH]CC([O])=C[O]'),
    E0 = (10.0605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.895905,0.0593956,-3.58801e-05,-1.87598e-09,7.00189e-12,1329.82,29.0267], Tmin=(100,'K'), Tmax=(946.472,'K')), NASAPolynomial(coeffs=[15.4649,0.0192123,-6.09342e-06,1.01928e-09,-6.9469e-14,-2386,-45.5196], Tmin=(946.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.0605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJC)"""),
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
    label = 'C[CH]C([O])=C[O](11449)',
    structure = SMILES('C[CH]C([O])=C[O]'),
    E0 = (39.2965,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,350,440,435,1725,3010,987.5,1337.5,450,1655,378.4,378.613,378.816],'cm^-1')),
        HinderedRotor(inertia=(0.0011754,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169462,'amu*angstrom^2'), symmetry=1, barrier=(17.2798,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.57806,0.0386823,8.05629e-06,-5.10777e-08,2.63096e-11,4827.06,24.096], Tmin=(100,'K'), Tmax=(916.47,'K')), NASAPolynomial(coeffs=[18.1151,0.00478919,8.6961e-07,-2.69149e-10,1.59419e-14,188.145,-63.0149], Tmin=(916.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(39.2965,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJCO)"""),
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
    label = '[CH2]C(C)[C]=C[O](9677)',
    structure = SMILES('[CH2]C(C)[C]=C[O]'),
    E0 = (327.32,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,469.183,469.2],'cm^-1')),
        HinderedRotor(inertia=(0.0883299,'amu*angstrom^2'), symmetry=1, barrier=(13.7968,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142589,'amu*angstrom^2'), symmetry=1, barrier=(22.2731,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0883173,'amu*angstrom^2'), symmetry=1, barrier=(13.7975,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34907,0.0476504,-1.09057e-05,-2.71027e-08,1.63636e-11,39472.8,25.1177], Tmin=(100,'K'), Tmax=(920.856,'K')), NASAPolynomial(coeffs=[15.1079,0.0156353,-3.95906e-06,5.93706e-10,-4.01706e-14,35762.3,-46.5151], Tmin=(920.856,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Isobutyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])C([CH2])C(2418)',
    structure = SMILES('[CH]=C([O])C([CH2])C'),
    E0 = (327.151,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,3120,650,792.5,1650,542.106],'cm^-1')),
        HinderedRotor(inertia=(0.0412581,'amu*angstrom^2'), symmetry=1, barrier=(8.61011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0412807,'amu*angstrom^2'), symmetry=1, barrier=(8.61057,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0412764,'amu*angstrom^2'), symmetry=1, barrier=(8.60984,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3675.75,'J/mol'), sigma=(6.3037,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=574.14 K, Pc=33.3 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.76918,0.0587136,-5.30515e-05,2.51057e-08,-4.59034e-12,39474,27.422], Tmin=(100,'K'), Tmax=(1496.09,'K')), NASAPolynomial(coeffs=[14.5955,0.0157063,-3.87533e-06,4.93779e-10,-2.66565e-14,36012.9,-42.5929], Tmin=(1496.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.151,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(Cds_P)"""),
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
    label = '[CH]C(C)C([O])=C[O](11490)',
    structure = SMILES('[CH]C(C)C([O])=C[O]'),
    E0 = (255.857,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,506.634,506.636,506.636,506.638,506.639,506.639],'cm^-1')),
        HinderedRotor(inertia=(0.000656751,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0774507,'amu*angstrom^2'), symmetry=1, barrier=(14.1072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0774487,'amu*angstrom^2'), symmetry=1, barrier=(14.1072,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.66694,0.059302,-2.37131e-05,-2.65366e-08,1.89092e-11,30905.5,27.8729], Tmin=(100,'K'), Tmax=(921.294,'K')), NASAPolynomial(coeffs=[20.7631,0.00865302,-8.44357e-07,3.93809e-11,-4.42607e-15,25649.2,-75.8629], Tmin=(921.294,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(255.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'CC1COC1=C[O](11491)',
    structure = SMILES('CC1COC1=C[O]'),
    E0 = (-154.32,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16275,0.034833,7.08228e-05,-1.39587e-07,6.35053e-11,-18432.3,19.0642], Tmin=(100,'K'), Tmax=(900.921,'K')), NASAPolynomial(coeffs=[25.4547,0.000310103,6.20865e-06,-1.42655e-09,9.57536e-14,-25785.3,-112.112], Tmin=(900.921,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-154.32,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C)C1=COO1(11492)',
    structure = SMILES('[CH2]C(C)C1=COO1'),
    E0 = (252.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37418,0.0428617,1.42788e-05,-5.33174e-08,2.46789e-11,30513.6,26.4104], Tmin=(100,'K'), Tmax=(955.157,'K')), NASAPolynomial(coeffs=[15.7834,0.01918,-6.10448e-06,1.09379e-09,-8.02676e-14,26088.6,-51.199], Tmin=(955.157,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.806,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(Isobutyl)"""),
)

species(
    label = 'CC1COC=C1[O](11493)',
    structure = SMILES('CC1COC=C1[O]'),
    E0 = (-237.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53513,0.0292514,7.36344e-05,-1.35471e-07,6.07342e-11,-28395.8,18.0543], Tmin=(100,'K'), Tmax=(895.825,'K')), NASAPolynomial(coeffs=[22.1911,0.00367946,4.83453e-06,-1.2049e-09,8.29546e-14,-34771.4,-94.2557], Tmin=(895.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-237.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C(C)C(O)=C[O](11494)',
    structure = SMILES('C=C(C)C(O)=C[O]'),
    E0 = (-224.892,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.16871,0.0684575,-3.2373e-05,-2.66505e-08,2.10746e-11,-26895.6,22.861], Tmin=(100,'K'), Tmax=(913.635,'K')), NASAPolynomial(coeffs=[24.0932,0.0059563,8.86789e-07,-3.12807e-10,2.01225e-14,-33030.3,-100.049], Tmin=(913.635,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-224.892,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C)C(O)=C=O(11495)',
    structure = SMILES('[CH2]C(C)C(O)=C=O'),
    E0 = (-84.3961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0552764,0.0820876,-9.29252e-05,5.09254e-08,-1.01648e-11,-10004.2,26.0047], Tmin=(100,'K'), Tmax=(941.534,'K')), NASAPolynomial(coeffs=[18.0269,0.0163854,-5.21642e-06,8.33771e-10,-5.3678e-14,-13860.4,-62.1218], Tmin=(941.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-84.3961,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + radical(CJC(C)C=C=O)"""),
)

species(
    label = 'C=C(C)C([O])=CO(11496)',
    structure = SMILES('C=C(C)C([O])=CO'),
    E0 = (-228.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.05417,0.085205,-9.09214e-05,4.56477e-08,-8.38887e-12,-27283.9,27.3502], Tmin=(100,'K'), Tmax=(1585.54,'K')), NASAPolynomial(coeffs=[23.3544,0.00539708,1.82746e-06,-6.01711e-10,4.69125e-14,-32732.6,-94.4313], Tmin=(1585.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-228.55,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ)"""),
)

species(
    label = 'CC(C)C([O])=C=O(11497)',
    structure = SMILES('CC(C)C(=O)[C]=O'),
    E0 = (-180.895,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30224,0.0632349,-6.22378e-05,3.60209e-08,-8.96265e-12,-21662.8,23.7539], Tmin=(100,'K'), Tmax=(944.203,'K')), NASAPolynomial(coeffs=[8.26088,0.0337558,-1.54066e-05,2.95548e-09,-2.07909e-13,-22976.9,-9.41807], Tmin=(944.203,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-180.895,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C](C)C([O])[CH][O](11498)',
    structure = SMILES('[CH2][C](C)C([O])[CH][O]'),
    E0 = (487.944,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,243.714,243.931,604.154,2681.99],'cm^-1')),
        HinderedRotor(inertia=(1.28496,'amu*angstrom^2'), symmetry=1, barrier=(54.3111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00283294,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0106416,'amu*angstrom^2'), symmetry=1, barrier=(54.2988,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0106386,'amu*angstrom^2'), symmetry=1, barrier=(54.3063,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.888229,0.0736287,-9.74546e-05,7.91463e-08,-2.63883e-11,58793.1,29.9961], Tmin=(100,'K'), Tmax=(827.323,'K')), NASAPolynomial(coeffs=[7.24092,0.0357626,-1.58338e-05,2.9266e-09,-1.98903e-13,57986.7,2.03145], Tmin=(827.323,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCJ(C)CO) + radical(CCsJOH) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C([O])[CH][O](11499)',
    structure = SMILES('[CH2]C([CH2])C([O])[CH][O]'),
    E0 = (540.452,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,180,1038.73,1040.03,1043.93],'cm^-1')),
        HinderedRotor(inertia=(0.124824,'amu*angstrom^2'), symmetry=1, barrier=(2.86994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.124054,'amu*angstrom^2'), symmetry=1, barrier=(2.85225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.128012,'amu*angstrom^2'), symmetry=1, barrier=(2.94325,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125695,'amu*angstrom^2'), symmetry=1, barrier=(2.88997,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.562757,0.0828065,-0.000124159,1.0475e-07,-3.43143e-11,65118.3,32.0569], Tmin=(100,'K'), Tmax=(899.471,'K')), NASAPolynomial(coeffs=[7.94717,0.0333216,-1.38754e-05,2.4356e-09,-1.58346e-13,64463.3,0.95688], Tmin=(899.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(540.452,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(Isobutyl) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C([CH2])[C]([O])C[O](11500)',
    structure = SMILES('[CH2]C([CH2])[C]([O])C[O]'),
    E0 = (536.782,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,180,1331.13,1333.78,1335.2],'cm^-1')),
        HinderedRotor(inertia=(0.160889,'amu*angstrom^2'), symmetry=1, barrier=(3.69916,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164655,'amu*angstrom^2'), symmetry=1, barrier=(3.78575,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163103,'amu*angstrom^2'), symmetry=1, barrier=(3.75006,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162993,'amu*angstrom^2'), symmetry=1, barrier=(3.74753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647839,0.0831676,-0.000131357,1.16578e-07,-3.96139e-11,64671.7,31.6262], Tmin=(100,'K'), Tmax=(898.27,'K')), NASAPolynomial(coeffs=[6.12293,0.0363001,-1.55442e-05,2.75628e-09,-1.79791e-13,64595.3,10.8492], Tmin=(898.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(536.782,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)[C]1OC1[O](11501)',
    structure = SMILES('[CH2]C(C)[C]1OC1[O]'),
    E0 = (250.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11781,0.0635277,-6.88639e-05,4.63965e-08,-1.25677e-11,30247.2,27.5312], Tmin=(100,'K'), Tmax=(1038.48,'K')), NASAPolynomial(coeffs=[7.91944,0.0298978,-9.55413e-06,1.4307e-09,-8.39126e-14,29235.3,-3.61014], Tmin=(1038.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C1([O])[CH]O1(11502)',
    structure = SMILES('[CH2]C(C)C1([O])[CH]O1'),
    E0 = (248.595,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.711,0.0759176,-7.39805e-05,3.65325e-08,-6.61369e-12,30092.7,30.968], Tmin=(100,'K'), Tmax=(1663.13,'K')), NASAPolynomial(coeffs=[16.6348,0.0140092,6.50853e-08,-4.47922e-10,4.23809e-14,27115.3,-53.1447], Tmin=(1663.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CCsJO) + radical(Isobutyl)"""),
)

species(
    label = 'CC1CC([O])[C]1[O](11503)',
    structure = SMILES('CC1CC([O])[C]1[O]'),
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
    label = 'C=C(C)C([O])=C[O](11505)',
    structure = SMILES('C=C(C)C([O])=C[O]'),
    E0 = (-87.0875,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.978649,'amu*angstrom^2'), symmetry=1, barrier=(22.5011,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.979747,'amu*angstrom^2'), symmetry=1, barrier=(22.5263,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.387243,0.075182,-7.66172e-05,3.73454e-08,-6.75338e-12,-10297.9,26.3693], Tmin=(100,'K'), Tmax=(1580.13,'K')), NASAPolynomial(coeffs=[21.0079,0.00785331,-2.02288e-07,-1.68562e-10,1.63037e-14,-15415.4,-81.4364], Tmin=(1580.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-87.0875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH2]C(C)C([O])=C=O(11506)',
    structure = SMILES('[CH2]C(C)C(=O)[C]=O'),
    E0 = (29.6155,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,1855,455,950,260.182],'cm^-1')),
        HinderedRotor(inertia=(0.146533,'amu*angstrom^2'), symmetry=1, barrier=(6.8058,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.145881,'amu*angstrom^2'), symmetry=1, barrier=(6.81218,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.14624,'amu*angstrom^2'), symmetry=1, barrier=(6.79633,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.435211,'amu*angstrom^2'), symmetry=1, barrier=(20.3695,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.788727,0.0761095,-0.000108882,8.93539e-08,-2.97493e-11,3672.33,26.6294], Tmin=(100,'K'), Tmax=(796.573,'K')), NASAPolynomial(coeffs=[8.85953,0.030773,-1.44544e-05,2.7472e-09,-1.89723e-13,2539.11,-9.51431], Tmin=(796.573,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(29.6155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJC(C)C=O) + radical(CCCJ=O)"""),
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
    label = 'C=CC([O])=C[O](11103)',
    structure = SMILES('C=CC([O])=C[O]'),
    E0 = (-49.4969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.20625,'amu*angstrom^2'), symmetry=1, barrier=(27.7341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36073,0.044686,-1.0018e-05,-3.58118e-08,2.21372e-11,-5845.51,19.4833], Tmin=(100,'K'), Tmax=(902.976,'K')), NASAPolynomial(coeffs=[19.5497,0.000521135,2.86666e-06,-6.71619e-10,4.54736e-14,-10614.7,-74.6305], Tmin=(902.976,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-49.4969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
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
    label = '[CH2][CH]C([O])=C[O](11447)',
    structure = SMILES('[CH2]C=C([O])[CH][O]'),
    E0 = (234.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,3025,407.5,1350,352.5,504.235,507.339,507.755],'cm^-1')),
        HinderedRotor(inertia=(0.000655989,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000659839,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91251,0.03846,-1.68356e-05,-7.51854e-09,6.05239e-12,28310.6,23.5168], Tmin=(100,'K'), Tmax=(1002.24,'K')), NASAPolynomial(coeffs=[11.8359,0.0149245,-5.66133e-06,1.0463e-09,-7.45272e-14,25514.4,-28.4065], Tmin=(1002.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(Allyl_P) + radical(C=CCJO)"""),
)

species(
    label = '[CH2][C](C)C([O])=C[O](11507)',
    structure = SMILES('[CH2][C](C)C([O])=C[O]'),
    E0 = (165.299,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,528.541,528.6,528.602],'cm^-1')),
        HinderedRotor(inertia=(0.000603298,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0860247,'amu*angstrom^2'), symmetry=1, barrier=(17.0702,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.142762,'amu*angstrom^2'), symmetry=1, barrier=(28.3025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01426,0.0511983,-4.92023e-06,-4.42708e-08,2.52979e-11,20001.9,28.2935], Tmin=(100,'K'), Tmax=(905.985,'K')), NASAPolynomial(coeffs=[19.2901,0.00948768,-3.97407e-07,-1.10606e-10,8.12785e-15,15090.6,-66.9011], Tmin=(905.985,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.299,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C([O])=C[O](11047)',
    structure = SMILES('[CH2]C([CH2])C([O])=C[O]'),
    E0 = (217.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,470.994,472.956,480.09],'cm^-1')),
        HinderedRotor(inertia=(0.000730474,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0715053,'amu*angstrom^2'), symmetry=1, barrier=(11.639,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.574826,'amu*angstrom^2'), symmetry=1, barrier=(93.1123,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.363967,0.072711,-7.46279e-05,3.7386e-08,-6.86842e-12,26373.3,33.4457], Tmin=(100,'K'), Tmax=(1607.87,'K')), NASAPolynomial(coeffs=[18.4543,0.00888314,7.88781e-07,-4.64388e-10,3.99306e-14,22520.9,-59.4403], Tmin=(1607.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(217.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C([O])=[C][O](11508)',
    structure = SMILES('[CH2]C(C)C([O])=[C][O]'),
    E0 = (252.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,490.763,490.765,490.768],'cm^-1')),
        HinderedRotor(inertia=(0.000699932,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0565415,'amu*angstrom^2'), symmetry=1, barrier=(9.66368,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0565432,'amu*angstrom^2'), symmetry=1, barrier=(9.66369,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.425525,0.0658142,-6.48893e-05,3.25737e-08,-6.2357e-12,30504.6,32.6465], Tmin=(100,'K'), Tmax=(1438.98,'K')), NASAPolynomial(coeffs=[16.5408,0.0136778,-2.89087e-06,3.05647e-10,-1.37971e-14,26626.6,-48.3246], Tmin=(1438.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = 'C[C](C)C([O])=C[O](11509)',
    structure = SMILES('C[C](C)C([O])=C[O]'),
    E0 = (-39.7838,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,360,370,350,350,440,435,1725,3010,987.5,1337.5,450,1655,538.308,538.313,538.316],'cm^-1')),
        HinderedRotor(inertia=(0.0810079,'amu*angstrom^2'), symmetry=1, barrier=(16.6581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810065,'amu*angstrom^2'), symmetry=1, barrier=(16.6581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0810058,'amu*angstrom^2'), symmetry=1, barrier=(16.6581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.981877,0.0497163,6.67366e-06,-5.52199e-08,2.80382e-11,-4660.87,25.933], Tmin=(100,'K'), Tmax=(931.534,'K')), NASAPolynomial(coeffs=[19.4426,0.0129527,-2.57364e-06,3.82352e-10,-3.01581e-14,-9944.48,-71.719], Tmin=(931.534,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-39.7838,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[CH2]C(C)C([O])=[C]O(11510)',
    structure = SMILES('[CH2]C(C)C([O])=[C]O'),
    E0 = (111.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,362.238,362.789],'cm^-1')),
        HinderedRotor(inertia=(0.0012745,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135925,'amu*angstrom^2'), symmetry=1, barrier=(12.7205,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.136869,'amu*angstrom^2'), symmetry=1, barrier=(12.7243,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135793,'amu*angstrom^2'), symmetry=1, barrier=(12.7073,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.240589,0.0758265,-7.91524e-05,4.08205e-08,-7.84803e-12,13518.6,33.6246], Tmin=(100,'K'), Tmax=(1485.93,'K')), NASAPolynomial(coeffs=[19.1685,0.0108263,-6.63653e-07,-1.69316e-10,1.99977e-14,9158.36,-62.9625], Tmin=(1485.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](C)C(O)=C[O](11511)',
    structure = SMILES('[CH2][C](C)C(O)=C[O]'),
    E0 = (27.4937,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,510.804,510.805],'cm^-1')),
        HinderedRotor(inertia=(0.107719,'amu*angstrom^2'), symmetry=1, barrier=(19.9449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107719,'amu*angstrom^2'), symmetry=1, barrier=(19.9449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107719,'amu*angstrom^2'), symmetry=1, barrier=(19.9449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.10772,'amu*angstrom^2'), symmetry=1, barrier=(19.9449,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.650385,0.0551348,2.79165e-06,-6.17493e-08,3.3571e-11,3444.7,28.0989], Tmin=(100,'K'), Tmax=(906.899,'K')), NASAPolynomial(coeffs=[23.0256,0.00613336,1.65665e-06,-5.01947e-10,3.35304e-14,-2657.03,-88.9275], Tmin=(906.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(27.4937,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C(O)=[C][O](11512)',
    structure = SMILES('[CH2]C(C)[C](O)[C]=O'),
    E0 = (111.864,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3100,440,815,1455,1000,1855,455,950,298.148],'cm^-1')),
        HinderedRotor(inertia=(0.00189662,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00189691,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00189694,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141862,'amu*angstrom^2'), symmetry=1, barrier=(8.94605,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.141859,'amu*angstrom^2'), symmetry=1, barrier=(8.9459,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.775381,0.0684336,-7.19716e-05,4.07694e-08,-9.21509e-12,13572.4,32.0131], Tmin=(100,'K'), Tmax=(1078.21,'K')), NASAPolynomial(coeffs=[13.1001,0.0227099,-8.35993e-06,1.43698e-09,-9.50709e-14,10914.7,-28.3742], Tmin=(1078.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(111.864,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C2CsJOH) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C([CH2])C(O)=C[O](11513)',
    structure = SMILES('[CH2]C([CH2])C(O)=C[O]'),
    E0 = (80.0021,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,316.006,316.753],'cm^-1')),
        HinderedRotor(inertia=(0.290569,'amu*angstrom^2'), symmetry=1, barrier=(20.3713,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.43327,'amu*angstrom^2'), symmetry=1, barrier=(100.899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00170306,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.289744,'amu*angstrom^2'), symmetry=1, barrier=(20.3738,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.18733,0.0818418,-8.41761e-05,4.12369e-08,-7.32386e-12,9836.71,34.9164], Tmin=(100,'K'), Tmax=(1688.06,'K')), NASAPolynomial(coeffs=[20.8829,0.00711672,2.15342e-06,-7.27919e-10,5.67399e-14,5580.98,-73.6505], Tmin=(1688.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(80.0021,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = 'CC(C)C([O])=[C][O](11514)',
    structure = SMILES('CC(C)C([O])=[C][O]'),
    E0 = (47.3865,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2770,2790,2810,2830,2850,1350,1400,1450,1500,700,800,1000,1100,1350,1400,900,1100,350,440,435,1725,1685,370,384.352,384.352,384.373],'cm^-1')),
        HinderedRotor(inertia=(0.00114112,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107841,'amu*angstrom^2'), symmetry=1, barrier=(11.3099,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.107865,'amu*angstrom^2'), symmetry=1, barrier=(11.3106,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.912731,0.0582115,-3.18714e-05,-6.28784e-09,8.51669e-12,5819.21,28.4212], Tmin=(100,'K'), Tmax=(952.764,'K')), NASAPolynomial(coeffs=[15.8213,0.0186954,-5.9866e-06,1.01943e-09,-7.06005e-14,1931.04,-48.2786], Tmin=(952.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(47.3865,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](C)C([O])=CO(11515)',
    structure = SMILES('[CH2]C(C)=C([O])[CH]O'),
    E0 = (-30.0521,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.7507,0.0607193,-3.35696e-05,-7.13734e-09,9.22436e-12,-3487.74,27.5816], Tmin=(100,'K'), Tmax=(958.532,'K')), NASAPolynomial(coeffs=[17.1298,0.0176416,-5.70717e-06,9.91207e-10,-6.99496e-14,-7788.74,-56.8008], Tmin=(958.532,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-30.0521,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C([O])=CO(11516)',
    structure = SMILES('[CH2]C([CH2])C([O])=CO'),
    E0 = (76.3443,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,350,440,435,1725,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,471.426,471.637],'cm^-1')),
        HinderedRotor(inertia=(0.0814294,'amu*angstrom^2'), symmetry=1, barrier=(12.829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0814259,'amu*angstrom^2'), symmetry=1, barrier=(12.8293,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.154125,'amu*angstrom^2'), symmetry=1, barrier=(24.3318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0813693,'amu*angstrom^2'), symmetry=1, barrier=(12.8196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.02998,0.0827253,-8.89098e-05,4.56682e-08,-8.49813e-12,9387.27,34.4232], Tmin=(100,'K'), Tmax=(1607.09,'K')), NASAPolynomial(coeffs=[20.7614,0.00648251,2.7907e-06,-8.9164e-10,7.00899e-14,5224.77,-72.205], Tmin=(1607.09,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(76.3443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(Isobutyl) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C(C)C(=O)C=O(10959)',
    structure = SMILES('[CH2]C(C)C(=O)C=O'),
    E0 = (-130.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00525,0.0712915,-8.65946e-05,6.45853e-08,-2.05396e-11,-15574,25.0441], Tmin=(100,'K'), Tmax=(753.854,'K')), NASAPolynomial(coeffs=[7.56142,0.036503,-1.73712e-05,3.36623e-09,-2.3696e-13,-16562.5,-4.73296], Tmin=(753.854,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-130.345,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C(C)C([O])[C]=O(10955)',
    structure = SMILES('[CH2]C(C)C([O])[C]=O'),
    E0 = (178.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950,766.553,2795.38],'cm^-1')),
        HinderedRotor(inertia=(0.518413,'amu*angstrom^2'), symmetry=1, barrier=(11.9193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.027975,'amu*angstrom^2'), symmetry=1, barrier=(11.9512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.125456,'amu*angstrom^2'), symmetry=1, barrier=(2.88448,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0276913,'amu*angstrom^2'), symmetry=1, barrier=(11.929,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0883,0.0577557,-4.22145e-05,1.13117e-08,8.01732e-13,21635.4,32.2451], Tmin=(100,'K'), Tmax=(978.621,'K')), NASAPolynomial(coeffs=[12.8949,0.0223973,-7.79023e-06,1.33013e-09,-8.92289e-14,18706.9,-27.6155], Tmin=(978.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Isobutyl) + radical(CCCJ=O)"""),
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
    label = '[CH2]C(C)[C]=O(2423)',
    structure = SMILES('[CH2]C(C)[C]=O'),
    E0 = (136.335,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.0939094,'amu*angstrom^2'), symmetry=1, barrier=(8.49727,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0939151,'amu*angstrom^2'), symmetry=1, barrier=(8.49711,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.093799,'amu*angstrom^2'), symmetry=1, barrier=(8.49804,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1249,0.0431168,-3.86264e-05,2.03172e-08,-4.61453e-12,16463.3,19.7425], Tmin=(100,'K'), Tmax=(1025,'K')), NASAPolynomial(coeffs=[7.06838,0.0238253,-1.03952e-05,1.95561e-09,-1.36146e-13,15449.9,-4.22909], Tmin=(1025,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(136.335,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CJC(C)C=O) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'CC1CC([O])C1=O(10967)',
    structure = SMILES('CC1CC([O])C1=O'),
    E0 = (-71.6187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.0709,0.0233306,6.3939e-05,-9.95533e-08,3.95494e-11,-8527.03,24.6959], Tmin=(100,'K'), Tmax=(962.252,'K')), NASAPolynomial(coeffs=[14.1681,0.0209432,-7.00778e-06,1.33161e-09,-1.01479e-13,-13072.7,-44.7237], Tmin=(962.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.6187,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C(C)C(=O)C[O](11517)',
    structure = SMILES('C=C(C)C(=O)C[O]'),
    E0 = (-71.3576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9469,0.0529016,-2.72682e-05,-2.83268e-08,3.78247e-11,-8516.12,22.9883], Tmin=(100,'K'), Tmax=(492.352,'K')), NASAPolynomial(coeffs=[4.82816,0.0398923,-1.9315e-05,3.80185e-09,-2.71098e-13,-8925.88,9.8495], Tmin=(492.352,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.3576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = '[CH2]C([CH2])[C](O)[CH][O](11518)',
    structure = SMILES('[CH2]C([CH2])[C](O)[CH][O]'),
    E0 = (486.719,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,1380,1390,370,380,2900,435,360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3025,407.5,1350,352.5,241.417,835.125,1821.95],'cm^-1')),
        HinderedRotor(inertia=(0.100161,'amu*angstrom^2'), symmetry=1, barrier=(3.35544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100161,'amu*angstrom^2'), symmetry=1, barrier=(3.35544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100161,'amu*angstrom^2'), symmetry=1, barrier=(3.35544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100161,'amu*angstrom^2'), symmetry=1, barrier=(3.35544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.100161,'amu*angstrom^2'), symmetry=1, barrier=(3.35544,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.2416,0.0955916,-0.000167355,1.52385e-07,-5.17148e-11,58662,32.8254], Tmin=(100,'K'), Tmax=(913.844,'K')), NASAPolynomial(coeffs=[6.68163,0.035111,-1.50765e-05,2.62721e-09,-1.67375e-13,58833.3,9.71359], Tmin=(913.844,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(486.719,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + longDistanceInteraction_noncyclic(CsCs-ST) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(C2CsJOH) + radical(Isobutyl) + radical(Isobutyl) + radical(CCsJOH)"""),
)

species(
    label = 'CC1CO[C]1[CH][O](11519)',
    structure = SMILES('CC1CO[C]1[CH][O]'),
    E0 = (261.315,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18395,0.0658432,-8.08216e-05,6.28319e-08,-1.943e-11,31526.7,24.5025], Tmin=(100,'K'), Tmax=(958.464,'K')), NASAPolynomial(coeffs=[5.87065,0.0341745,-1.23085e-05,1.99532e-09,-1.235e-13,31184.5,4.99235], Tmin=(958.464,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(261.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Oxetane) + radical(CCOJ) + radical(C2CsJOCs) + radical(CCsJOH)"""),
)

species(
    label = '[CH2]C(C)[C]1[CH]OO1(11520)',
    structure = SMILES('[CH2]C(C)[C]1[CH]OO1'),
    E0 = (470.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.33602,0.063544,-7.48732e-05,6.16228e-08,-2.18013e-11,56669.2,23.5554], Tmin=(100,'K'), Tmax=(803.454,'K')), NASAPolynomial(coeffs=[4.58713,0.040217,-1.79904e-05,3.36152e-09,-2.307e-13,56377.3,10.0166], Tmin=(803.454,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(470.417,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + ring(12dioxetane) + radical(C2CsJOO) + radical(Isobutyl) + radical(CCsJOO)"""),
)

species(
    label = '[CH2][C](C)C(=O)C[O](11521)',
    structure = SMILES('[CH2]C(C)=C([O])C[O]'),
    E0 = (78.3578,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,279.504,279.633,1764.17],'cm^-1')),
        HinderedRotor(inertia=(0.00215064,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0946127,'amu*angstrom^2'), symmetry=1, barrier=(5.238,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.094224,'amu*angstrom^2'), symmetry=1, barrier=(5.24056,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22464,0.0650775,-7.1199e-05,4.84165e-08,-1.41452e-11,9520.58,26.1566], Tmin=(100,'K'), Tmax=(815.523,'K')), NASAPolynomial(coeffs=[7.39169,0.0348299,-1.55653e-05,2.93849e-09,-2.04167e-13,8514.69,-2.33834], Tmin=(815.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(78.3578,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(CCOJ) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH2]C([CH2])C(=O)C[O](11522)',
    structure = SMILES('[CH2]C([CH2])C(=O)C[O]'),
    E0 = (231.766,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,280.342,286.801,2896.96],'cm^-1')),
        HinderedRotor(inertia=(0.00213667,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.165971,'amu*angstrom^2'), symmetry=1, barrier=(9.18268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157977,'amu*angstrom^2'), symmetry=1, barrier=(9.19272,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164158,'amu*angstrom^2'), symmetry=1, barrier=(9.17802,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.308112,0.089469,-0.000141802,1.22791e-07,-4.14928e-11,28000,29.097], Tmin=(100,'K'), Tmax=(855.183,'K')), NASAPolynomial(coeffs=[9.0499,0.0329727,-1.53307e-05,2.85806e-09,-1.93351e-13,27075.6,-8.37269], Tmin=(855.183,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.766,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsCs) + radical(C=OCOJ) + radical(CJC(C)C=O) + radical(CJC(C)C=O)"""),
)

species(
    label = '[CH2]C=C([CH][O])OC(11523)',
    structure = SMILES('[CH2]C=C([CH][O])OC'),
    E0 = (115.587,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,180,180,881.966,882.314],'cm^-1')),
        HinderedRotor(inertia=(0.0440359,'amu*angstrom^2'), symmetry=1, barrier=(24.3185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0440226,'amu*angstrom^2'), symmetry=1, barrier=(24.3159,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00892531,'amu*angstrom^2'), symmetry=1, barrier=(4.92724,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05761,'amu*angstrom^2'), symmetry=1, barrier=(24.3165,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.996781,0.0545331,-2.01684e-05,-1.47762e-08,9.81789e-12,14020.2,26.9384], Tmin=(100,'K'), Tmax=(1015.61,'K')), NASAPolynomial(coeffs=[15.5074,0.022079,-8.71096e-06,1.64639e-09,-1.18548e-13,9799.06,-49.5624], Tmin=(1015.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.587,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-OsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=CCJO) + radical(Allyl_P)"""),
)

species(
    label = '[CH2][CH]C(C)([O])C=O(11524)',
    structure = SMILES('[CH2][CH]C(C)([O])C=O'),
    E0 = (216.238,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,180,180,180,1478.36,1600,1900.2,3200],'cm^-1')),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.157779,'amu*angstrom^2'), symmetry=1, barrier=(3.62765,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18088,0.0650514,-7.03975e-05,4.46845e-08,-1.18796e-11,26106.3,30.506], Tmin=(100,'K'), Tmax=(900.764,'K')), NASAPolynomial(coeffs=[8.88742,0.0308293,-1.3409e-05,2.50671e-09,-1.73495e-13,24718,-5.86814], Tmin=(900.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(216.238,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CC(C)(C=O)OJ) + radical(CCJCO) + radical(RCCJ)"""),
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
    label = '[CH2]C(C)[C][O](11525)',
    structure = SMILES('[CH2]C(C)[C][O]'),
    E0 = (578.603,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3100,440,815,1455,1000,228.611,228.799,2859.22],'cm^-1')),
        HinderedRotor(inertia=(0.122586,'amu*angstrom^2'), symmetry=1, barrier=(4.54399,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0032218,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.86797,'amu*angstrom^2'), symmetry=1, barrier=(69.3366,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0898,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.9221,0.0486457,-5.83975e-05,4.50023e-08,-1.46658e-11,69661.9,20.863], Tmin=(100,'K'), Tmax=(806.152,'K')), NASAPolynomial(coeffs=[6.01479,0.0258054,-1.11855e-05,2.06144e-09,-1.40435e-13,69084.3,2.51058], Tmin=(806.152,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(578.603,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(Isobutyl) + radical(CH2_triplet)"""),
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
    label = 'C=C(C)C([O])C=O(11527)',
    structure = SMILES('C=C(C)C([O])C=O'),
    E0 = (-65.6594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30452,0.0551735,-4.00281e-05,1.49313e-08,-2.28467e-12,-7796.83,27.4712], Tmin=(100,'K'), Tmax=(1508.83,'K')), NASAPolynomial(coeffs=[12.4638,0.0255894,-1.06171e-05,1.93619e-09,-1.31485e-13,-11164.3,-30.9561], Tmin=(1508.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-65.6594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ)"""),
)

species(
    label = 'CC(C)=C([O])C=O(11528)',
    structure = SMILES('C[C](C)C(=O)C=O'),
    E0 = (-184.028,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.918501,0.0724905,-9.32732e-05,7.39262e-08,-2.43364e-11,-22027,22.3971], Tmin=(100,'K'), Tmax=(806.776,'K')), NASAPolynomial(coeffs=[7.53275,0.0352203,-1.56551e-05,2.90983e-09,-1.98856e-13,-22948.5,-7.18973], Tmin=(806.776,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-184.028,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(C2CJCHO)"""),
)

species(
    label = 'CC1CO[CH][C]1[O](11529)',
    structure = SMILES('CC1CO[CH][C]1[O]'),
    E0 = (183.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2973,0.0483961,1.27963e-06,-4.69086e-08,2.63337e-11,22188.8,20.2247], Tmin=(100,'K'), Tmax=(865.501,'K')), NASAPolynomial(coeffs=[15.1487,0.0174565,-2.42299e-06,9.81848e-11,1.6898e-15,18552.2,-51.7563], Tmin=(865.501,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(183.59,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + ring(Tetrahydrofuran) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOCs)"""),
)

species(
    label = 'CC=C([O])C=O(11462)',
    structure = SMILES('C[CH]C(=O)C=O'),
    E0 = (-111.413,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,246.751],'cm^-1')),
        HinderedRotor(inertia=(0.00269295,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00334131,'amu*angstrom^2'), symmetry=1, barrier=(9.95884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.752681,'amu*angstrom^2'), symmetry=1, barrier=(32.993,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34071,0.0408507,-3.0048e-05,1.2415e-08,-2.3527e-12,-13343.9,20.2827], Tmin=(100,'K'), Tmax=(1129.61,'K')), NASAPolynomial(coeffs=[5.86847,0.0283588,-1.34602e-05,2.62542e-09,-1.86123e-13,-14140.9,2.83337], Tmin=(1129.61,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-111.413,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CCJCC=O)"""),
)

species(
    label = '[CH2][C](C)C([O])C=O(11530)',
    structure = SMILES('[CH2][C](C)C([O])C=O'),
    E0 = (171.583,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,591.499,591.521],'cm^-1')),
        HinderedRotor(inertia=(0.00808743,'amu*angstrom^2'), symmetry=1, barrier=(2.00777,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00808665,'amu*angstrom^2'), symmetry=1, barrier=(2.00753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00808383,'amu*angstrom^2'), symmetry=1, barrier=(2.00696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0873176,'amu*angstrom^2'), symmetry=1, barrier=(2.0076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4112,0.048388,-1.66447e-05,-1.16777e-08,7.74218e-12,20737.3,30.4632], Tmin=(100,'K'), Tmax=(1003.57,'K')), NASAPolynomial(coeffs=[11.8254,0.0251994,-9.36771e-06,1.67804e-09,-1.16143e-13,17724.5,-24.4128], Tmin=(1003.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.583,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ(C)CO) + radical(Isobutyl)"""),
)

species(
    label = '[CH2]C([CH2])C([O])C=O(11531)',
    structure = SMILES('[CH2]C([CH2])C([O])C=O'),
    E0 = (224.091,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1383.33,1386.67,1390,370,373.333,376.667,380,2800,3000,430,440,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,2782.5,750,1395,475,1775,1000,907.622,910.816],'cm^-1')),
        HinderedRotor(inertia=(0.00524476,'amu*angstrom^2'), symmetry=1, barrier=(3.06262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.13402,'amu*angstrom^2'), symmetry=1, barrier=(3.08139,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135111,'amu*angstrom^2'), symmetry=1, barrier=(3.10647,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135301,'amu*angstrom^2'), symmetry=1, barrier=(3.11084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19112,0.0563233,-3.90771e-05,8.6625e-09,1.85428e-12,27058,32.1466], Tmin=(100,'K'), Tmax=(941.86,'K')), NASAPolynomial(coeffs=[11.8979,0.023843,-8.03777e-06,1.33606e-09,-8.79815e-14,24464.9,-21.9249], Tmin=(941.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(224.091,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Isobutyl) + radical(Isobutyl)"""),
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
    E0 = (12.7246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (463.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (197.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (172.66,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (477.098,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (846.674,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (846.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (467.65,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (21.0089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (252.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (19.8374,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (76.1248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (76.1248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (37.6979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (37.6979,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (510.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (603.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (561.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (250.767,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (249.011,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (273.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (270.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (124.705,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (279.602,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (116.592,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (230.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (255.507,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (371.258,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (489.24,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (377.091,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (429.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (464.261,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (163.558,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (273.785,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (181.863,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (305.242,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (163.055,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (91.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (126.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (165.568,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (12.7246,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (342.342,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (622.89,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (21.0089,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (76.1248,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (511.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (261.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (470.771,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (207.815,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (307.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (429.387,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (415.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (666.589,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (18.1192,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (35.586,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (40.5775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (183.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (309.9,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (329.782,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (342.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['OCHCO(3676)', 'C3H6(27)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['CH2(S)(11)', '[CH2]CC([O])=C[O](10944)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(143764,'m^3/(mol*s)'), n=0.444, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [carbene;R_H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: 1,2_Insertion_carbene
Ea raised from -5.1 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CC[CH]C([O])=C[O](11489)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(6.59e+09,'s^-1'), n=0.99, Ea=(182.004,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCsCJ;CsJ-CdH;CH3] for rate rule [cCs(-HH)CJ;CsJ-CdH;CH3]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction5',
    reactants = ['CH2(T)(20)', 'C[CH]C([O])=C[O](11449)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H/OneDeC;Birad]
Euclidian distance = 4.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['O(4)', '[CH2]C(C)[C]=C[O](9677)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['O(4)', '[CH]=C([O])C([CH2])C(2418)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', '[CH]C(C)C([O])=C[O](11490)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC1COC1=C[O](11491)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['[CH2]C(C)C1=COO1(11492)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(240.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 237.9 to 240.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC1COC=C1[O](11493)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSDS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['C=C(C)C(O)=C[O](11494)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['[CH2]C(C)C(O)=C=O(11495)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['C=C(C)C([O])=CO(11496)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC(C)C([O])=C=O(11497)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2][C](C)C([O])[CH][O](11498)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C([CH2])C([O])[CH][O](11499)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C([CH2])[C]([O])C[O](11500)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(2.05689e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['[CH2]C(C)[C]1OC1[O](11501)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(238.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['[CH2]C(C)C1([O])[CH]O1(11502)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(236.287,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC1CC([O])[C]1[O](11503)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.51071e+08,'s^-1'), n=0.996667, Ea=(260.735,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_Cs_RR_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 259.0 to 260.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC1CC1([O])[CH][O](11504)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(257.716,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', 'C=C(C)C([O])=C[O](11505)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.88727,'m^3/(mol*s)'), n=1.9687, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 18 used for Cds-CdCs_Cds-HH;HJ
Exact match found for rate rule [Cds-CdCs_Cds-HH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -8.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[CH2]C(C)C([O])=C=O(11506)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CH3(17)', 'C=CC([O])=C[O](11103)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(13200,'cm^3/(mol*s)'), n=2.41, Ea=(29.539,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 813 used for Cds-CdH_Cds-HH;CsJ-HHH
Exact match found for rate rule [Cds-CdH_Cds-HH;CsJ-HHH]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O][C]=C[O](9592)', 'C3H6(27)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(0.00168615,'m^3/(mol*s)'), n=2.52599, Ea=(19.6608,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-CsH_Cds-HH;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['OCHCO(3676)', 'C3H6(T)(28)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction28',
    reactants = ['CH3(17)', '[CH2][CH]C([O])=C[O](11447)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O][C]=C[O](9592)', 'C3H6(T)(28)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH2][C](C)C([O])=C[O](11507)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['H(3)', '[CH2]C([CH2])C([O])=C[O](11047)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(5.6103e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['H(3)', '[CH2]C(C)C([O])=[C][O](11508)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C[C](C)C([O])=C[O](11509)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(4.614e+09,'s^-1'), n=1.31, Ea=(203.342,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 163 used for R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H
Exact match found for rate rule [R2H_S;C_rad_out_OneDe/Cs;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C(C)C([O])=[C]O(11510)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][C](C)C(O)=C[O](11511)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(52326.9,'s^-1'), n=2.1859, Ea=(154.369,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;Y_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2]C(C)C(O)=[C][O](11512)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH2]C([CH2])C(O)=C[O](11513)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(0.0192749,'s^-1'), n=3.795, Ea=(83.0524,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;C_rad_out_2H;O_H_out] + [R4H_SS(Cd)S;C_rad_out_2H;XH_out] for rate rule [R4H_SS(Cd)S;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['CC(C)C([O])=[C][O](11514)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(222600,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 6.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['[CH2][C](C)C([O])=CO(11515)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.11e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2]C([CH2])C([O])=CO(11516)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(969255,'s^-1'), n=1.705, Ea=(89.2238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5H_SSMS;C_rad_out_2H;XH_out] for rate rule [R5H_SSMS;C_rad_out_2H;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['[CH2]C(C)C(=O)C=O(10959)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(C)C([O])[C]=O(10955)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH][O](1548)', '[CH2]C(C)[C]=O(2423)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC1CC([O])C1=O(10967)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['C=C(C)C(=O)C[O](11517)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2]C([CH2])[C](O)[CH][O](11518)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC1CO[C]1[CH][O](11519)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.27074e+08,'s^-1'), n=0.924088, Ea=(249.065,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs2H] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['[CH2]C(C)[C]1[CH]OO1(11520)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(458.047,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra_O] for rate rule [R4_linear;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2][C](C)C(=O)C[O](11521)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(11891.5,'s^-1'), n=2.58622, Ea=(129.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH2]C([CH2])C(=O)C[O](11522)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1855.84,'s^-1'), n=2.50105, Ea=(76.0245,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_SSS;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2]C=C([CH][O])OC(11523)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2][CH]C(C)([O])C=O(11524)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.37952e+08,'s^-1'), n=1.37167, Ea=(199.228,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CH3]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction53',
    reactants = ['HCO(1372)', '[CH2]C(C)[C][O](11525)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC1CC1([O])C=O(11526)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(5.39457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['C=C(C)C([O])C=O(11527)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC(C)=C([O])C=O(11528)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH2]C(C)C([O])=C[O](10965)'],
    products = ['CC1CO[CH][C]1[O](11529)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.19156e+09,'s^-1'), n=0.640131, Ea=(170.866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 169.7 to 170.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction58',
    reactants = ['CH2(T)(20)', 'CC=C([O])C=O(11462)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(9.2456e+13,'m^3/(mol*s)'), n=-2.73375, Ea=(38.9008,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CH2_triplet] + [Cd_R;Y_1centerbirad] for rate rule [Cd_R;CH2_triplet]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction59',
    reactants = ['[CH2][C](C)C([O])C=O(11530)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH2]C([CH2])C([O])C=O(11531)'],
    products = ['[CH2]C(C)C([O])=C[O](10965)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(166690,'s^-1'), n=2.17519, Ea=(117.939,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_2H;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2620',
    isomers = [
        '[CH2]C(C)C([O])=C[O](10965)',
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
    label = 'PDepNetwork #2620',
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


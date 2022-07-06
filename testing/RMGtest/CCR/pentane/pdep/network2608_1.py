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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4030.69,'J/mol'), sigma=(6.74566,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=629.58 K, Pc=29.8 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0883,0.0577557,-4.22145e-05,1.13117e-08,8.01732e-13,21635.4,32.2451], Tmin=(100,'K'), Tmax=(978.621,'K')), NASAPolynomial(coeffs=[12.8949,0.0223973,-7.79023e-06,1.33013e-09,-8.92289e-14,18706.9,-27.6155], Tmin=(978.621,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Isobutyl) + radical(CCCJ=O)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4286.44,'J/mol'), sigma=(7.03247,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.53 K, Pc=27.97 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.895905,0.0593956,-3.58801e-05,-1.87598e-09,7.00189e-12,1329.82,29.0267], Tmin=(100,'K'), Tmax=(946.472,'K')), NASAPolynomial(coeffs=[15.4649,0.0192123,-6.09342e-06,1.01928e-09,-6.9469e-14,-2386,-45.5196], Tmin=(946.472,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(10.0605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(RCCJC)"""),
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
    label = '[CH2]C([O])[C]=O(2866)',
    structure = SMILES('[CH2]C([O])[C]=O'),
    E0 = (242.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1855,455,950,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0550228,'amu*angstrom^2'), symmetry=1, barrier=(13.826,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0550218,'amu*angstrom^2'), symmetry=1, barrier=(13.8262,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.32663,0.0326765,-2.63448e-05,6.89772e-09,8.77003e-13,29179.2,22.4436], Tmin=(100,'K'), Tmax=(963.067,'K')), NASAPolynomial(coeffs=[10.4241,0.00826416,-2.68163e-06,4.57445e-10,-3.15247e-14,27192,-18.5368], Tmin=(963.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CJCO) + radical(CCCJ=O)"""),
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
    label = 'C[CH]C[CH][C]=O(2434)',
    structure = SMILES('C[CH]C[CH][C]=O'),
    E0 = (271.053,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,2191.4,2191.84],'cm^-1')),
        HinderedRotor(inertia=(0.00214683,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.109702,'amu*angstrom^2'), symmetry=1, barrier=(6.30235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.21024,'amu*angstrom^2'), symmetry=1, barrier=(67.5894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.113548,'amu*angstrom^2'), symmetry=1, barrier=(6.32349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3428.95,'J/mol'), sigma=(5.98513,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.59 K, Pc=36.29 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23108,0.0440156,-2.01448e-05,-2.09479e-08,2.52039e-11,32658.9,23.6185], Tmin=(100,'K'), Tmax=(539.722,'K')), NASAPolynomial(coeffs=[4.99773,0.0326275,-1.38305e-05,2.54716e-09,-1.74683e-13,32227.4,10.7471], Tmin=(539.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(271.053,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[C]=O(2355)',
    structure = SMILES('[C]=O'),
    E0 = (438.818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3323.79],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (28.0101,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.09152,0.00193295,-1.59153e-05,2.47563e-08,-1.11287e-11,52770.5,4.46624], Tmin=(100,'K'), Tmax=(866.029,'K')), NASAPolynomial(coeffs=[1.05092,0.00526657,-3.13864e-06,6.40624e-10,-4.47939e-14,53698.8,21.0169], Tmin=(866.029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(438.818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'C[CH]C[CH][O](1560)',
    structure = SMILES('C[CH]C[CH][O]'),
    E0 = (296.318,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,180,1857.88,1857.94],'cm^-1')),
        HinderedRotor(inertia=(0.102893,'amu*angstrom^2'), symmetry=1, barrier=(2.36571,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.103017,'amu*angstrom^2'), symmetry=1, barrier=(2.36856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102914,'amu*angstrom^2'), symmetry=1, barrier=(2.36619,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.23223,0.0293516,-1.06816e-05,9.14899e-10,6.61011e-14,35593.9,14.3951], Tmin=(100,'K'), Tmax=(2829.07,'K')), NASAPolynomial(coeffs=[40.9762,-0.00885552,2.28855e-06,-4.24163e-10,3.26746e-14,9303.22,-210.806], Tmin=(2829.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + radical(CCOJ) + radical(RCCJC) + radical(CCsJOH)"""),
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
    label = '[CH]CC([O])[C]=O(11466)',
    structure = SMILES('[CH]CC([O])[C]=O'),
    E0 = (454.924,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1855,455,950,712.262,712.316,712.321,712.369,4000],'cm^-1')),
        HinderedRotor(inertia=(0.154807,'amu*angstrom^2'), symmetry=1, barrier=(3.55931,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00988277,'amu*angstrom^2'), symmetry=1, barrier=(3.55829,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.536417,'amu*angstrom^2'), symmetry=1, barrier=(12.3333,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64269,0.0469256,-4.54413e-05,2.24036e-08,-4.3449e-12,54803.7,27.5152], Tmin=(100,'K'), Tmax=(1257.38,'K')), NASAPolynomial(coeffs=[12.3632,0.0128213,-4.75647e-06,8.32441e-10,-5.59955e-14,52107.8,-26.6606], Tmin=(1257.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(454.924,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = 'C[C]CC([O])[C]=O(11611)',
    structure = SMILES('C[C]CC([O])[C]=O'),
    E0 = (431.144,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,1855,455,950,180,180,383.035,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0168718,'amu*angstrom^2'), symmetry=1, barrier=(16.1109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.700732,'amu*angstrom^2'), symmetry=1, barrier=(16.1112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.700725,'amu*angstrom^2'), symmetry=1, barrier=(16.111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.016872,'amu*angstrom^2'), symmetry=1, barrier=(16.1118,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05115,0.0610864,-5.72261e-05,2.79466e-08,-5.46925e-12,51963.9,30.7836], Tmin=(100,'K'), Tmax=(1230.74,'K')), NASAPolynomial(coeffs=[13.3569,0.0210915,-8.4808e-06,1.54209e-09,-1.05661e-13,48934.9,-31.1393], Tmin=(1230.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.144,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]CC1OC1=O(11612)',
    structure = SMILES('C[CH]CC1OC1=O'),
    E0 = (-89.6295,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00062,0.0459887,-2.32544e-05,5.27561e-09,-4.63274e-13,-10710.5,24.6728], Tmin=(100,'K'), Tmax=(2584.75,'K')), NASAPolynomial(coeffs=[19.0426,0.0196159,-7.94976e-06,1.32825e-09,-8.14862e-14,-19520.5,-73.729], Tmin=(2584.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.6295,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(2(co)oxirane) + radical(RCCJC)"""),
)

species(
    label = 'CC1CC([C]=O)O1(10957)',
    structure = SMILES('CC1CC([C]=O)O1'),
    E0 = (-89.4158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81789,0.0299505,5.41819e-05,-1.0143e-07,4.51562e-11,-10658.7,23.7197], Tmin=(100,'K'), Tmax=(894.656,'K')), NASAPolynomial(coeffs=[16.0399,0.0145886,-9.16349e-07,-1.22919e-10,1.12721e-14,-15133.4,-54.0962], Tmin=(894.656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-89.4158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(CCCJ=O)"""),
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
    label = 'CC=CC(O)[C]=O(11613)',
    structure = SMILES('CC=CC(O)[C]=O'),
    E0 = (-146.403,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.05386,0.0581792,-4.65624e-05,1.91991e-08,-3.1885e-12,-17496.7,29.4616], Tmin=(100,'K'), Tmax=(1429.88,'K')), NASAPolynomial(coeffs=[13.9396,0.0221322,-8.74795e-06,1.56857e-09,-1.05996e-13,-21181.7,-37.3127], Tmin=(1429.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-146.403,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCCJ=O)"""),
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
    label = 'C=CCC(O)[C]=O(10962)',
    structure = SMILES('C=CCC(O)[C]=O'),
    E0 = (-133.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10729,0.0539221,-2.70319e-05,-5.75477e-09,6.5971e-12,-15886.2,30.37], Tmin=(100,'K'), Tmax=(1010.2,'K')), NASAPolynomial(coeffs=[14.6189,0.0207028,-7.82128e-06,1.44179e-09,-1.02321e-13,-19651,-40.0751], Tmin=(1010.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-133.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCCJ=O)"""),
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
    label = 'CC=CC([O])[C]=O(11614)',
    structure = SMILES('CC=CC([O])[C]=O'),
    E0 = (97.3306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.161063,'amu*angstrom^2'), symmetry=1, barrier=(3.70315,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160332,'amu*angstrom^2'), symmetry=1, barrier=(3.68635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.161692,'amu*angstrom^2'), symmetry=1, barrier=(3.71761,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4906,0.0530521,-4.31524e-05,1.86604e-08,-3.32661e-12,11798.2,28.6677], Tmin=(100,'K'), Tmax=(1313.72,'K')), NASAPolynomial(coeffs=[11.0212,0.0240335,-1.0019e-05,1.84631e-09,-1.2691e-13,9294.06,-19.9124], Tmin=(1313.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(97.3306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
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
    label = 'C=CCC([O])[C]=O(11043)',
    structure = SMILES('C=CCC([O])[C]=O'),
    E0 = (110.711,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,313.654,313.793,4000],'cm^-1')),
        HinderedRotor(inertia=(0.00171319,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234198,'amu*angstrom^2'), symmetry=1, barrier=(16.3913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.234731,'amu*angstrom^2'), symmetry=1, barrier=(16.3889,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30815,0.0514594,-3.23939e-05,4.3152e-09,2.25592e-12,13419,30.43], Tmin=(100,'K'), Tmax=(1047.03,'K')), NASAPolynomial(coeffs=[13.1076,0.0203556,-7.85325e-06,1.43654e-09,-1.00357e-13,10182.2,-30.6953], Tmin=(1047.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.711,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C[CH]CC=O(1680)',
    structure = SMILES('C[CH]CC=O'),
    E0 = (-27.2013,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2782.5,750,1395,475,1775,1000,1638.36],'cm^-1')),
        HinderedRotor(inertia=(0.00206383,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134446,'amu*angstrom^2'), symmetry=1, barrier=(7.78225,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.134363,'amu*angstrom^2'), symmetry=1, barrier=(7.78254,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0978,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2085,0.0350017,-1.46906e-05,8.4427e-10,6.0281e-13,-3203.51,20.5797], Tmin=(100,'K'), Tmax=(1434.58,'K')), NASAPolynomial(coeffs=[8.72968,0.0234851,-9.61916e-06,1.72668e-09,-1.15452e-13,-5760.5,-15.6254], Tmin=(1434.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-27.2013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCC=O)"""),
)

species(
    label = 'C[CH]CC=C=O(2432)',
    structure = SMILES('C[CH]CC=C=O'),
    E0 = (68.9827,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,246.668,1625.9],'cm^-1')),
        HinderedRotor(inertia=(0.115588,'amu*angstrom^2'), symmetry=1, barrier=(5.00744,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.11461,'amu*angstrom^2'), symmetry=1, barrier=(5.00193,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114333,'amu*angstrom^2'), symmetry=1, barrier=(4.96075,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.1085,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.14409,0.0442309,-3.17307e-05,1.54891e-08,-3.86259e-12,8360.58,22.4694], Tmin=(100,'K'), Tmax=(855.785,'K')), NASAPolynomial(coeffs=[4.10428,0.0350693,-1.56735e-05,2.98095e-09,-2.088e-13,8025.06,13.3177], Tmin=(855.785,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(68.9827,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(RCCJC)"""),
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
    label = 'C[CH][CH]C([O])[C]=O(11615)',
    structure = SMILES('C[CH][CH]C([O])[C]=O'),
    E0 = (377.277,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,700.517,2932.96,2933.78],'cm^-1')),
        HinderedRotor(inertia=(0.143908,'amu*angstrom^2'), symmetry=1, barrier=(3.30874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00967865,'amu*angstrom^2'), symmetry=1, barrier=(3.34657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.146146,'amu*angstrom^2'), symmetry=1, barrier=(3.36017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.053869,'amu*angstrom^2'), symmetry=1, barrier=(18.6996,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64853,0.0512526,-4.24436e-05,1.96698e-08,-3.83957e-12,45460.9,33.9654], Tmin=(100,'K'), Tmax=(1195.15,'K')), NASAPolynomial(coeffs=[9.20933,0.0259477,-1.06843e-05,1.95415e-09,-1.33852e-13,43653.7,-3.85895], Tmin=(1195.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(377.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(RCCJC) + radical(CCCJ=O)"""),
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
    label = '[CH2][CH]CC([O])[C]=O(11035)',
    structure = SMILES('[CH2][CH]CC([O])[C]=O'),
    E0 = (382.621,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,2126.01,2126.01,2126.02],'cm^-1')),
        HinderedRotor(inertia=(0.0264061,'amu*angstrom^2'), symmetry=1, barrier=(11.6678,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.507464,'amu*angstrom^2'), symmetry=1, barrier=(11.6676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0355688,'amu*angstrom^2'), symmetry=1, barrier=(2.90401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0264116,'amu*angstrom^2'), symmetry=1, barrier=(11.6679,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0999,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58089,0.0560124,-5.68327e-05,3.47251e-08,-9.06541e-12,46103.4,33.5079], Tmin=(100,'K'), Tmax=(908.16,'K')), NASAPolynomial(coeffs=[7.63512,0.0293473,-1.27917e-05,2.39629e-09,-1.66153e-13,45003.7,4.88274], Tmin=(908.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(382.621,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(RCCJC) + radical(RCCJ) + radical(CCCJ=O)"""),
)

species(
    label = 'CC[CH]C([O])[C]=O(11593)',
    structure = SMILES('CC[CH]C([O])[C]=O'),
    E0 = (182.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,1855,455,950,180,374.673,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0156426,'amu*angstrom^2'), symmetry=1, barrier=(17.7199,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.015657,'amu*angstrom^2'), symmetry=1, barrier=(17.7262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00156177,'amu*angstrom^2'), symmetry=1, barrier=(17.7324,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.771392,'amu*angstrom^2'), symmetry=1, barrier=(17.7358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.1765,0.0542492,-3.41105e-05,5.37361e-09,1.88839e-12,22097.8,32.9249], Tmin=(100,'K'), Tmax=(1056.03,'K')), NASAPolynomial(coeffs=[13.093,0.0227501,-8.74003e-06,1.58622e-09,-1.09981e-13,18820.5,-28.8154], Tmin=(1056.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(182.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(CCCJ=O)"""),
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
    label = '[CH2]CCC([O])[C]=O(11616)',
    structure = SMILES('[CH2]CCC([O])[C]=O'),
    E0 = (188.175,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1855,455,950,237.8,4000,4000],'cm^-1')),
        HinderedRotor(inertia=(0.0027107,'amu*angstrom^2'), symmetry=1, barrier=(3.22676,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.525503,'amu*angstrom^2'), symmetry=1, barrier=(21.0886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.525528,'amu*angstrom^2'), symmetry=1, barrier=(21.0886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.525503,'amu*angstrom^2'), symmetry=1, barrier=(21.0886,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11387,0.0590621,-4.90782e-05,2.14801e-08,-3.82194e-12,22739.7,32.4404], Tmin=(100,'K'), Tmax=(1332.53,'K')), NASAPolynomial(coeffs=[12.7915,0.0240085,-9.61967e-06,1.73927e-09,-1.18349e-13,19627.5,-27.2502], Tmin=(1332.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.175,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(RCCJ) + radical(CCCJ=O)"""),
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
    label = 'C[CH][CH]C(O)[C]=O(11617)',
    structure = SMILES('C[CH][CH]C(O)[C]=O'),
    E0 = (133.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25328,0.0559644,-4.46997e-05,1.9054e-08,-3.32555e-12,16164.1,34.605], Tmin=(100,'K'), Tmax=(1351.3,'K')), NASAPolynomial(coeffs=[11.9599,0.0242714,-9.51897e-06,1.69739e-09,-1.14434e-13,13270.5,-20.2716], Tmin=(1351.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCO) + radical(RCCJC) + radical(CCCJ=O)"""),
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
    label = '[CH2][CH]CC(O)[C]=O(11618)',
    structure = SMILES('[CH2][CH]CC(O)[C]=O'),
    E0 = (138.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.1079,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.36314,0.058752,-5.27459e-05,2.66144e-08,-5.63929e-12,16798.7,33.5026], Tmin=(100,'K'), Tmax=(1110.26,'K')), NASAPolynomial(coeffs=[9.69315,0.0287405,-1.21985e-05,2.26697e-09,-1.56825e-13,14949,-7.55601], Tmin=(1110.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(138.888,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(RCCJC) + radical(RCCJ) + radical(CCCJ=O)"""),
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
    E0 = (177.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (338.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (340.747,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (641.361,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (790.406,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (790.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (646.864,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (642.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (180.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (185.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (185.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (200.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (200.236,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (240.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (240.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (240.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (202.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (202.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (193.852,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (314.855,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (258.128,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (325.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (467.703,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (340.086,'kJ/mol'),
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
    E0 = (267.542,'kJ/mol'),
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
    E0 = (589.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (461.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (594.413,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (341.823,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (291.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (335.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (321.416,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (252.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (333.452,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (237.269,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (296.28,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['OCHCO(3676)', 'C3H6(27)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[CH2]C(C)C([O])[C]=O(10955)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(6.55606e+10,'s^-1'), n=0.64, Ea=(159.935,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [cCs(-HC)CJ;CsJ;C] for rate rule [cCs(-HC)CJ;CsJ-HH;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['C[CH]CC([O])=C[O](10963)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['CHCH3(T)(21)', '[CH2]C([O])[C]=O(2866)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', 'C[CH]C[CH][C]=O(2434)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 2 used for Y_rad;O_birad
Exact match found for rate rule [Y_rad;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[C]=O(2355)', 'C[CH]C[CH][O](1560)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['CH3(17)', '[CH]CC([O])[C]=O(11466)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_methyl;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['H(3)', 'C[C]CC([O])[C]=O(11611)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['C[CH]CC1OC1=O(11612)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['CC1CC([C]=O)O1(10957)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['CC1CC([O])C1=O(10967)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_single;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_H/NonDeC;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['C[CH]CC(O)=C=O(11552)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['C[CH]CC(=O)C=O(10960)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['CC=CC(O)[C]=O(11613)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['CC=CC([O])C=O(11579)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad_NDe]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['CCCC([O])=C=O(11554)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['C=CCC(O)[C]=O(10962)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['C=CCC([O])C=O(10972)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.37831e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction19',
    reactants = ['CO(2039)', 'C[CH]C[CH][O](1560)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1099.39,'m^3/(mol*s)'), n=0.635039, Ea=(16.7529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;Y_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', 'CC=CC([O])[C]=O(11614)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.46e+08,'cm^3/(mol*s)'), n=1.64, Ea=(5.73208,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2555 used for Cds-CsH_Cds-CsH;HJ
Exact match found for rate rule [Cds-CsH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', 'C[CH]CC([O])=C=O(11561)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(92.1383,'m^3/(mol*s)'), n=1.68375, Ea=(21.5685,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', 'C=CCC([O])[C]=O(11043)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.36e+08,'cm^3/(mol*s)'), n=1.56, Ea=(2.5104,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 9 used for Cds-HH_Cds-CsH;HJ
Exact match found for rate rule [Cds-HH_Cds-CsH;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[C]=O(2355)', 'C[CH]CC=O(1680)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.0201871,'m^3/(mol*s)'), n=2.2105, Ea=(56.0866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [CO-CsH_O;YJ] for rate rule [CO-CsH_O;CO_birad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['O(4)', 'C[CH]CC=C=O(2432)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(11.358,'m^3/(mol*s)'), n=1.81033, Ea=(28.0982,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;YJ] for rate rule [Cds-CsH_Ck;O_atom_triplet]
Euclidian distance = 2.23606797749979
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][C]=C[O](9592)', 'C3H6(27)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(0.00620445,'m^3/(mol*s)'), n=2.46568, Ea=(12.4666,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Cds-Cs\H3/H;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction26',
    reactants = ['OCHCO(3676)', 'C3H6(T)(28)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2339.95,'m^3/(mol*s)'), n=0.573452, Ea=(58.2237,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO_O;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][C]=C[O](9592)', 'C3H6(T)(28)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['H(3)', 'C[CH][CH]C([O])[C]=O(11615)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['H(3)', 'C[CH]CC([O])=[C][O](11563)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['H(3)', '[CH2][CH]CC([O])[C]=O(11035)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['CC[CH]C([O])[C]=O(11593)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(6.76e+09,'s^-1'), n=0.88, Ea=(158.992,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 357 used for R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC
Exact match found for rate rule [R2H_S;C_rad_out_H/NonDeC;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['C[CH]CC(O)=[C][O](11567)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.15e+14,'s^-1','+|-',2), n=-0.27, Ea=(113.972,'kJ/mol'), T0=(1,'K'), Tmin=(700,'K'), Tmax=(1800,'K'), comment="""Estimated using an average for rate rule [R2H_S;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2]CCC([O])[C]=O(11616)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.48e+07,'s^-1'), n=1.57, Ea=(147.695,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 106 used for R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)
Exact match found for rate rule [R2H_S;C_rad_out_2H;Cs_H_out_H/(NonDeC/Cs)]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['CCCC([O])=[C][O](11568)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(588307,'s^-1'), n=1.79367, Ea=(144.041,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;C_rad_out_H/NonDeC;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['C[CH][CH]C(O)[C]=O(11617)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(223829,'s^-1'), n=2.27675, Ea=(75.2806,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C[CH][CH]C([O])C=O(11581)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(285601,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;Y_rad_out;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C[CH]CC([O])[C]=O(10953)'],
    products = ['[CH2][CH]CC(O)[C]=O(11618)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(7.8e+08,'s^-1'), n=0.775, Ea=(59.894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH2][CH]CC([O])C=O(11582)'],
    products = ['C[CH]CC([O])[C]=O(10953)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5HJ_1;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2608',
    isomers = [
        'C[CH]CC([O])[C]=O(10953)',
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
    label = 'PDepNetwork #2608',
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


species(
    label = '[O]OC=[C]CC[C]=O(21049)',
    structure = SMILES('[O]OC=[C]CC[C]=O'),
    E0 = (342.408,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.310148,'amu*angstrom^2'), symmetry=1, barrier=(7.13091,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310591,'amu*angstrom^2'), symmetry=1, barrier=(7.14109,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310173,'amu*angstrom^2'), symmetry=1, barrier=(7.13148,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.310478,'amu*angstrom^2'), symmetry=1, barrier=(7.13849,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.570522,0.0817866,-0.000125582,1.05095e-07,-3.45277e-11,41299.6,32.6177], Tmin=(100,'K'), Tmax=(855.488,'K')), NASAPolynomial(coeffs=[9.65862,0.0283992,-1.28718e-05,2.37598e-09,-1.60036e-13,40143.3,-7.47866], Tmin=(855.488,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(342.408,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    label = 'C=C=CO[O](16806)',
    structure = SMILES('C=C=CO[O]'),
    E0 = (250.435,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055],'cm^-1')),
        HinderedRotor(inertia=(0.895395,'amu*angstrom^2'), symmetry=1, barrier=(20.5869,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3584.1,'J/mol'), sigma=(5.7752,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=559.83 K, Pc=42.22 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2049,0.0363638,-3.70353e-05,1.95556e-08,-4.06022e-12,30187.6,17.323], Tmin=(100,'K'), Tmax=(1179.28,'K')), NASAPolynomial(coeffs=[9.983,0.00998085,-3.47674e-06,5.84074e-10,-3.83218e-14,28353.2,-21.4843], Tmin=(1179.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(250.435,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = 'C=C([CH]O[O])C[C]=O(19754)',
    structure = SMILES('C=C([CH]O[O])C[C]=O'),
    E0 = (215.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.731857,0.0627192,-5.10548e-05,1.96349e-08,-2.96228e-12,26036.1,30.8655], Tmin=(100,'K'), Tmax=(1582.77,'K')), NASAPolynomial(coeffs=[18.63,0.0174865,-8.18728e-06,1.57894e-09,-1.10307e-13,20370.4,-63.7011], Tmin=(1582.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C([O])C[C]=CO[O](21045)',
    structure = SMILES('C=C([O])C[C]=CO[O]'),
    E0 = (317.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.528941,0.0746253,-8.84397e-05,5.39374e-08,-1.28791e-11,38318.3,31.9761], Tmin=(100,'K'), Tmax=(1029.93,'K')), NASAPolynomial(coeffs=[14.9298,0.0186956,-6.98313e-06,1.21107e-09,-8.06107e-14,35351.9,-37.9246], Tmin=(1029.93,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(317.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(ROOJ) + radical(Cds_S)"""),
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
    label = 'C#CCC[C]=O(15744)',
    structure = SMILES('C#CCC[C]=O'),
    E0 = (202.963,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2175,525,1855,455,950,750,770,3400,2100,186.169],'cm^-1')),
        HinderedRotor(inertia=(0.307164,'amu*angstrom^2'), symmetry=1, barrier=(7.55447,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307159,'amu*angstrom^2'), symmetry=1, barrier=(7.55449,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.69484,'amu*angstrom^2'), symmetry=1, barrier=(90.8731,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68386,0.0526803,-6.42762e-05,4.63544e-08,-1.36212e-11,24492.5,21.3307], Tmin=(100,'K'), Tmax=(877.307,'K')), NASAPolynomial(coeffs=[7.93162,0.0221852,-8.70158e-06,1.513e-09,-9.92497e-14,23473.6,-7.55272], Tmin=(877.307,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(202.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(CCCJ=O)"""),
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
    label = '[O]C=[C]CC[C]=O(21732)',
    structure = SMILES('[O]C=[C]CC[C]=O'),
    E0 = (207.371,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,327.011,327.102,327.133],'cm^-1')),
        HinderedRotor(inertia=(0.163166,'amu*angstrom^2'), symmetry=1, barrier=(12.3861,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.163177,'amu*angstrom^2'), symmetry=1, barrier=(12.3857,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.16317,'amu*angstrom^2'), symmetry=1, barrier=(12.3864,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04954,0.062582,-6.63192e-05,3.6427e-08,-7.93006e-12,25049.2,27.5075], Tmin=(100,'K'), Tmax=(1118.92,'K')), NASAPolynomial(coeffs=[13.3351,0.0186628,-7.44211e-06,1.3473e-09,-9.22269e-14,22299.9,-33.1435], Tmin=(1118.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(207.371,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(C=COJ) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[C]=CO[O](21346)',
    structure = SMILES('[C]=CO[O]'),
    E0 = (656.011,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655],'cm^-1')),
        HinderedRotor(inertia=(0.353101,'amu*angstrom^2'), symmetry=1, barrier=(8.11848,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.89582,0.0288275,-5.64855e-05,5.40338e-08,-1.90645e-11,78935.3,14.7462], Tmin=(100,'K'), Tmax=(877.636,'K')), NASAPolynomial(coeffs=[5.05128,0.00860812,-4.16047e-06,7.90553e-10,-5.33044e-14,78957.3,6.90957], Tmin=(877.636,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(656.011,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH2]C[C]=O(2354)',
    structure = SMILES('[CH2]C[C]=O'),
    E0 = (167.212,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(0.182443,'amu*angstrom^2'), symmetry=1, barrier=(4.19471,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.180479,'amu*angstrom^2'), symmetry=1, barrier=(4.14956,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.63264,0.0331121,-4.63433e-05,4.02607e-08,-1.40727e-11,20157.3,15.987], Tmin=(100,'K'), Tmax=(833.295,'K')), NASAPolynomial(coeffs=[4.91248,0.0168159,-7.37391e-06,1.37543e-09,-9.40242e-14,19963.2,6.51892], Tmin=(833.295,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.212,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJCC=O) + radical(CCCJ=O)"""),
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
    label = '[CH2]C[C]=CO[O](21345)',
    structure = SMILES('[CH2]C[C]=CO[O]'),
    E0 = (494.273,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,180],'cm^-1')),
        HinderedRotor(inertia=(0.261531,'amu*angstrom^2'), symmetry=1, barrier=(6.01312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261497,'amu*angstrom^2'), symmetry=1, barrier=(6.01234,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261502,'amu*angstrom^2'), symmetry=1, barrier=(6.01244,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0813,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72566,0.0514423,-5.66939e-05,3.51403e-08,-8.91765e-12,59528,25.1679], Tmin=(100,'K'), Tmax=(951.302,'K')), NASAPolynomial(coeffs=[8.93614,0.0211229,-8.88508e-06,1.63507e-09,-1.1224e-13,58156.2,-9.25832], Tmin=(951.302,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(494.273,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(ROOJ) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = 'O=[C]CCC1=COO1(21733)',
    structure = SMILES('O=[C]CCC1=COO1'),
    E0 = (132.857,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.79537,0.0608457,-5.0646e-05,2.08239e-08,-3.37811e-12,16102.4,29.8179], Tmin=(100,'K'), Tmax=(1483,'K')), NASAPolynomial(coeffs=[16.701,0.0179439,-7.2518e-06,1.31625e-09,-8.95246e-14,11384.9,-53.1852], Tmin=(1483,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.857,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclobutene) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OC=C1CCC1=O(21734)',
    structure = SMILES('[O]OC=C1CCC1=O'),
    E0 = (55.545,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41331,0.0492758,-2.7872e-05,5.11828e-09,3.23868e-13,6779.86,26.7424], Tmin=(100,'K'), Tmax=(1321.79,'K')), NASAPolynomial(coeffs=[12.8319,0.0244284,-1.06908e-05,2.00897e-09,-1.38951e-13,2913.26,-34.7392], Tmin=(1321.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.545,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(303.478,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(ROOJ)"""),
)

species(
    label = 'O=[C]CCC1=CO1(21735)',
    structure = SMILES('O=[C]CCC1=CO1'),
    E0 = (144.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03981,0.059783,-5.8501e-05,2.92084e-08,-5.75091e-12,17509.6,25.6086], Tmin=(100,'K'), Tmax=(1236.6,'K')), NASAPolynomial(coeffs=[14.3056,0.0168723,-6.45024e-06,1.14718e-09,-7.78588e-14,14228.6,-41.2086], Tmin=(1236.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.659,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OC=C=CCC=O(19339)',
    structure = SMILES('[O]OC=C=CCC=O'),
    E0 = (115.147,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21735,0.0602105,-4.93183e-05,1.96482e-08,-3.16478e-12,13949.5,27.1726], Tmin=(100,'K'), Tmax=(1443.13,'K')), NASAPolynomial(coeffs=[14.2242,0.0241589,-1.18463e-05,2.33773e-09,-1.66028e-13,10195.4,-40.3491], Tmin=(1443.13,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(115.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ)"""),
)

species(
    label = '[O]OC=CCC=C=O(21736)',
    structure = SMILES('[O]OC=CCC=C=O'),
    E0 = (69.9729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.769605,0.0745935,-9.65362e-05,6.86239e-08,-1.9657e-11,8528.92,28.3497], Tmin=(100,'K'), Tmax=(851.559,'K')), NASAPolynomial(coeffs=[11.0295,0.0264023,-1.16527e-05,2.17355e-09,-1.49494e-13,6781.47,-19.5001], Tmin=(851.559,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(69.9729,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + radical(ROOJ)"""),
)

species(
    label = 'O=[C]CC=C=COO(19762)',
    structure = SMILES('O=[C]CC=C=COO'),
    E0 = (123.103,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.816781,0.0674814,-6.20299e-05,2.75848e-08,-4.88321e-12,14922.3,28.6733], Tmin=(100,'K'), Tmax=(1344.88,'K')), NASAPolynomial(coeffs=[16.2888,0.0214651,-1.07073e-05,2.14447e-09,-1.54228e-13,10760.6,-50.555], Tmin=(1344.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.103,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(CCCJ=O)"""),
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
    label = '[O]OC=C=CC[C]=O(19763)',
    structure = SMILES('[O]OC=C=CC[C]=O'),
    E0 = (275.107,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.663975,'amu*angstrom^2'), symmetry=1, barrier=(15.2661,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663689,'amu*angstrom^2'), symmetry=1, barrier=(15.2595,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.663987,'amu*angstrom^2'), symmetry=1, barrier=(15.2664,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2176,0.0622872,-6.12895e-05,3.02236e-08,-6.02072e-12,33187,27.9944], Tmin=(100,'K'), Tmax=(1193.91,'K')), NASAPolynomial(coeffs=[13.1214,0.0224056,-1.11831e-05,2.24468e-09,-1.62042e-13,30344.6,-31.5441], Tmin=(1193.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(275.107,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OC=[C]CC=C=O(21190)',
    structure = SMILES('[O]OC=[C]CC=C=O'),
    E0 = (307.815,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.473241,'amu*angstrom^2'), symmetry=1, barrier=(10.8807,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470317,'amu*angstrom^2'), symmetry=1, barrier=(10.8135,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.470777,'amu*angstrom^2'), symmetry=1, barrier=(10.8241,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.613409,0.0807678,-0.000127435,1.06755e-07,-3.4855e-11,37137.5,29.376], Tmin=(100,'K'), Tmax=(860.371,'K')), NASAPolynomial(coeffs=[10.1831,0.0252924,-1.15673e-05,2.1348e-09,-1.43417e-13,35897.3,-12.9905], Tmin=(860.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(307.815,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsOsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OC#CCC[C]=O(21737)',
    structure = SMILES('[O]OC#CCC[C]=O'),
    E0 = (330.08,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2100,2250,500,550,1855,455,950,180,3967.34],'cm^-1')),
        HinderedRotor(inertia=(0.564842,'amu*angstrom^2'), symmetry=1, barrier=(12.9868,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.02321,'amu*angstrom^2'), symmetry=1, barrier=(69.5095,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564919,'amu*angstrom^2'), symmetry=1, barrier=(12.9886,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.03487,'amu*angstrom^2'), symmetry=1, barrier=(69.7776,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.819138,0.076772,-0.000120285,1.02435e-07,-3.37835e-11,39807.6,30.701], Tmin=(100,'K'), Tmax=(883.721,'K')), NASAPolynomial(coeffs=[8.45877,0.0276211,-1.21242e-05,2.18233e-09,-1.44018e-13,39026.3,-1.99214], Tmin=(883.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.08,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCt) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CtCsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtOs) + radical(ROOJ) + radical(CCCJ=O)"""),
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
    label = '[CH2][C]=CO[O](16807)',
    structure = SMILES('C=[C][CH]O[O]'),
    E0 = (437.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,274.987],'cm^-1')),
        HinderedRotor(inertia=(0.170957,'amu*angstrom^2'), symmetry=1, barrier=(9.18865,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171244,'amu*angstrom^2'), symmetry=1, barrier=(9.18809,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.2197,0.041819,-5.77934e-05,4.42117e-08,-1.35872e-11,52657.9,18.5046], Tmin=(100,'K'), Tmax=(795.764,'K')), NASAPolynomial(coeffs=[7.632,0.0146118,-6.50507e-06,1.24122e-09,-8.65948e-14,51796.6,-6.36981], Tmin=(795.764,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(437.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S)"""),
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
    label = '[CH]=[C]CC[C]=O(15746)',
    structure = SMILES('[CH]=[C]CC[C]=O'),
    E0 = (521.797,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,1855,455,950,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.295067,'amu*angstrom^2'), symmetry=1, barrier=(6.78417,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295034,'amu*angstrom^2'), symmetry=1, barrier=(6.78341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.295118,'amu*angstrom^2'), symmetry=1, barrier=(6.78534,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (81.0926,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3431.29,'J/mol'), sigma=(5.7779,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=535.96 K, Pc=40.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61071,0.0572252,-8.17608e-05,6.91131e-08,-2.35744e-11,62839.2,24.6265], Tmin=(100,'K'), Tmax=(825.336,'K')), NASAPolynomial(coeffs=[6.60863,0.0258906,-1.18862e-05,2.23081e-09,-1.52662e-13,62256.5,2.9413], Tmin=(825.336,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(521.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[O]O[CH][C]=CC[C]=O(19765)',
    structure = SMILES('[O]O[CH][C]=CC[C]=O'),
    E0 = (456.307,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,1855,455,950,180,1226.48],'cm^-1')),
        HinderedRotor(inertia=(0.757327,'amu*angstrom^2'), symmetry=1, barrier=(17.4124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.755169,'amu*angstrom^2'), symmetry=1, barrier=(17.3628,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0663061,'amu*angstrom^2'), symmetry=1, barrier=(17.3434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0163394,'amu*angstrom^2'), symmetry=1, barrier=(17.3822,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37159,0.0576157,-4.88245e-05,1.99773e-08,-3.30078e-12,54975.4,30.0694], Tmin=(100,'K'), Tmax=(1410.07,'K')), NASAPolynomial(coeffs=[13.7748,0.0224311,-1.13958e-05,2.28138e-09,-1.63362e-13,51477.6,-34.0309], Tmin=(1410.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.307,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[O]O[CH][C]CC=C=O(20801)',
    structure = SMILES('[O]O[CH][C]CC=C=O'),
    E0 = (575.016,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2120,512.5,787.5,223.385,223.423,223.455,2716.9],'cm^-1')),
        HinderedRotor(inertia=(0.280985,'amu*angstrom^2'), symmetry=1, barrier=(9.95769,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.281029,'amu*angstrom^2'), symmetry=1, barrier=(9.95749,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35833,'amu*angstrom^2'), symmetry=1, barrier=(48.1307,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.35903,'amu*angstrom^2'), symmetry=1, barrier=(48.1309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.466685,0.0865729,-0.000143817,1.2741e-07,-4.37274e-11,69277.2,30.0186], Tmin=(100,'K'), Tmax=(850.868,'K')), NASAPolynomial(coeffs=[8.75721,0.0305241,-1.49079e-05,2.8237e-09,-1.92337e-13,68484.4,-5.00742], Tmin=(850.868,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(575.016,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(CCsJOOH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O]O[C]=[C]CC[C]=O(21738)',
    structure = SMILES('[O]O[C]=[C]CC[C]=O'),
    E0 = (582.152,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1670,1700,300,440,1855,455,950,180,4000],'cm^-1')),
        HinderedRotor(inertia=(0.351542,'amu*angstrom^2'), symmetry=1, barrier=(8.08265,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352124,'amu*angstrom^2'), symmetry=1, barrier=(8.09602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.352086,'amu*angstrom^2'), symmetry=1, barrier=(8.09514,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.351511,'amu*angstrom^2'), symmetry=1, barrier=(8.08193,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.764854,0.0837464,-0.000153267,1.45071e-07,-5.13311e-11,70121.2,35.0475], Tmin=(100,'K'), Tmax=(880.778,'K')), NASAPolynomial(coeffs=[5.57344,0.0324273,-1.56607e-05,2.91367e-09,-1.94578e-13,70417.7,18.9509], Tmin=(880.778,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.152,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(Cds_S) + radical(CCCJ=O) + radical(C=CJO)"""),
)

species(
    label = '[O]O[CH]C=CC[C]=O(19767)',
    structure = SMILES('[O]O[CH]C=CC[C]=O'),
    E0 = (218.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.884205,0.0588142,-4.28745e-05,1.33276e-08,-1.29494e-12,26395.5,31.3687], Tmin=(100,'K'), Tmax=(1332.36,'K')), NASAPolynomial(coeffs=[17.0801,0.0195886,-9.29361e-06,1.81885e-09,-1.28809e-13,21245.6,-54.545], Tmin=(1332.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.465,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(C=CCJO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]C=CC[C]=CO[O](21132)',
    structure = SMILES('[O]C=CC[C]=CO[O]'),
    E0 = (326.968,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.105988,0.0737044,-7.68378e-05,3.92526e-08,-7.66295e-12,39475.6,33.2075], Tmin=(100,'K'), Tmax=(1344.1,'K')), NASAPolynomial(coeffs=[19.5729,0.011712,-3.12499e-06,4.44437e-10,-2.67907e-14,34609.2,-65.1022], Tmin=(1344.1,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.968,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]O[C]=CCC[C]=O(21739)',
    structure = SMILES('[O]O[C]=CCC[C]=O'),
    E0 = (344.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.81061,0.0789901,-0.000127907,1.15041e-07,-4.00412e-11,41517.3,34.4101], Tmin=(100,'K'), Tmax=(861.282,'K')), NASAPolynomial(coeffs=[6.68898,0.0330566,-1.54591e-05,2.88291e-09,-1.9478e-13,41195.8,10.9402], Tmin=(861.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(344.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(C=CJO) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C]CC[C]=[C]OO(21740)',
    structure = SMILES('O=[C]CC[C]=[C]OO'),
    E0 = (430.147,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1670,1700,300,440,1855,455,950,236.43,236.431],'cm^-1')),
        HinderedRotor(inertia=(0.00301577,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195577,'amu*angstrom^2'), symmetry=1, barrier=(7.75801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195576,'amu*angstrom^2'), symmetry=1, barrier=(7.75801,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195577,'amu*angstrom^2'), symmetry=1, barrier=(7.758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.195577,'amu*angstrom^2'), symmetry=1, barrier=(7.75801,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.507689,0.0872631,-0.000148349,1.35634e-07,-4.76244e-11,51850.5,35.2115], Tmin=(100,'K'), Tmax=(856.662,'K')), NASAPolynomial(coeffs=[7.19986,0.0340214,-1.66127e-05,3.14504e-09,-2.1392e-13,51511,8.67137], Tmin=(856.662,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(430.147,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(Cds_S) + radical(CCCJ=O) + radical(C=CJO)"""),
)

species(
    label = '[O]OC=CC[CH][C]=O(10594)',
    structure = SMILES('[O]OC=CC[CH][C]=O'),
    E0 = (272.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1855,455,950,278.406,279.615],'cm^-1')),
        HinderedRotor(inertia=(0.00215126,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171404,'amu*angstrom^2'), symmetry=1, barrier=(9.47852,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170755,'amu*angstrom^2'), symmetry=1, barrier=(9.48814,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170473,'amu*angstrom^2'), symmetry=1, barrier=(9.50175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.806059,0.0703343,-8.14039e-05,5.04446e-08,-1.24592e-11,32840.5,31.0621], Tmin=(100,'K'), Tmax=(988.657,'K')), NASAPolynomial(coeffs=[12.6081,0.0225861,-8.96203e-06,1.59737e-09,-1.0767e-13,30506.8,-25.7416], Tmin=(988.657,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]O[CH][C]=CCC=O(19772)',
    structure = SMILES('[O]O[CH][C]=CCC=O'),
    E0 = (296.346,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,3025,407.5,1350,352.5,1685,370,180,1554.25],'cm^-1')),
        HinderedRotor(inertia=(0.118343,'amu*angstrom^2'), symmetry=1, barrier=(21.0675,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.916391,'amu*angstrom^2'), symmetry=1, barrier=(21.0696,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.916282,'amu*angstrom^2'), symmetry=1, barrier=(21.0671,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0122872,'amu*angstrom^2'), symmetry=1, barrier=(21.0678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2166,0.0570855,-4.1118e-05,1.35957e-08,-1.77454e-12,35745.4,29.8233], Tmin=(100,'K'), Tmax=(1772.53,'K')), NASAPolynomial(coeffs=[17.3666,0.0206404,-1.02764e-05,1.99581e-09,-1.38474e-13,30020.1,-57.3353], Tmin=(1772.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(296.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S)"""),
)

species(
    label = '[O]O[C]=[C]CCC=O(21741)',
    structure = SMILES('[O]O[C]=[C]CCC=O'),
    E0 = (422.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2782.5,750,1395,475,1775,1000,1670,1700,300,440,180,3521.01],'cm^-1')),
        HinderedRotor(inertia=(0.341271,'amu*angstrom^2'), symmetry=1, barrier=(7.84649,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341125,'amu*angstrom^2'), symmetry=1, barrier=(7.84313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.341111,'amu*angstrom^2'), symmetry=1, barrier=(7.8428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.340965,'amu*angstrom^2'), symmetry=1, barrier=(7.83946,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.847426,0.0806844,-0.000138012,1.3082e-07,-4.72921e-11,50880.5,33.9317], Tmin=(100,'K'), Tmax=(859.794,'K')), NASAPolynomial(coeffs=[4.6656,0.0374513,-1.8152e-05,3.42894e-09,-2.33e-13,51165.4,21.5627], Tmin=(859.794,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'O=[C]CC=[C][CH]OO(19769)',
    structure = SMILES('O=[C]CC=[C][CH]OO'),
    E0 = (304.302,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874677,0.0637767,-5.22628e-05,2.00373e-08,-3.03743e-12,36715.4,31.1053], Tmin=(100,'K'), Tmax=(1554.71,'K')), NASAPolynomial(coeffs=[17.7671,0.0203153,-1.03307e-05,2.05657e-09,-1.46097e-13,31462.8,-57.8454], Tmin=(1554.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.302,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=CCJO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][CH]C[C]=COO(21742)',
    structure = SMILES('O=[C][CH]C[C]=COO'),
    E0 = (357.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,1855,455,950,337.722,337.749],'cm^-1')),
        HinderedRotor(inertia=(0.00147785,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150453,'amu*angstrom^2'), symmetry=1, barrier=(12.1797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.15045,'amu*angstrom^2'), symmetry=1, barrier=(12.1797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.150472,'amu*angstrom^2'), symmetry=1, barrier=(12.1797,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.884549,'amu*angstrom^2'), symmetry=1, barrier=(71.5977,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (113.091,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.575121,0.0776829,-9.82454e-05,6.58148e-08,-1.7552e-11,43170.6,31.6103], Tmin=(100,'K'), Tmax=(916.943,'K')), NASAPolynomial(coeffs=[12.8928,0.0239523,-1.03543e-05,1.91716e-09,-1.31673e-13,40911.6,-26.7483], Tmin=(916.943,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(286.849,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    E0 = (342.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (512.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (588.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (342.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (450.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (878.613,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (988.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (350.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (350.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (473.995,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (420.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (431.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (367.381,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (399.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (495.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (530.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (556.982,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (427.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (421.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (342.408,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (513.175,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (598.241,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (668.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (786.808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (793.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (500.087,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (499.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (576.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (567.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (484.385,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (458.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (455.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (393.285,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (432.376,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['ketene(1375)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['C=C([CH]O[O])C[C]=O(19754)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(1.74842e+09,'s^-1'), n=1.084, Ea=(170.038,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [cCsCJ;CdsJ;C] + [cCs(-HH)CJ;CJ;C] for rate rule [cCs(-HH)CJ;CdsJ;C]
Euclidian distance = 1.0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['C=C([O])C[C]=CO[O](21045)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['O2(S)(666)', 'C#CCC[C]=O(15744)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[O]C=[C]CC[C]=O(21732)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [O_sec_rad;O_birad] for rate rule [O_rad/OneDe;O_birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction6',
    reactants = ['[C]=CO[O](21346)', '[CH2]C[C]=O(2354)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[C]=O(2355)', '[CH2]C[C]=CO[O](21345)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cs;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['O=[C]CCC1=COO1(21733)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;CdsinglepriND_rad_out]
Euclidian distance = 2.449489742783178
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['[O]OC=C1CCC1=O(21734)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['O(4)', 'O=[C]CCC1=CO1(21735)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO;Y_rad_intra;OOJ] for rate rule [R2OO_D;Cd_rad_in/NonDeC;OOJ]
Euclidian distance = 3.1622776601683795
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['[O]OC=C=CCC=O(19339)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['[O]OC=CCC=C=O(21736)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(5.2748e+09,'s^-1'), n=0.37, Ea=(88.9686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad_De;XH_Rrad] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['O=[C]CC=C=COO(19762)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['CO(2039)', '[CH2]C[C]=CO[O](21345)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(754.069,'m^3/(mol*s)'), n=1.0822, Ea=(24.7919,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [COm;C_rad/H2/Cs]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[O]OC=C=CC[C]=O(19763)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.23e+08,'cm^3/(mol*s)'), n=1.64, Ea=(8.49352,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2566 used for Cds-CsH_Ca;HJ
Exact match found for rate rule [Cds-CsH_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[O]OC=[C]CC=C=O(21190)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.82e-16,'cm^3/(molecule*s)'), n=1.61, Ea=(10.992,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ck;HJ] for rate rule [Cds-CsH_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[O]OC#CCC[C]=O(21737)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2][C]=O(1376)', 'C=C=CO[O](16806)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(0.00429749,'m^3/(mol*s)'), n=2.45395, Ea=(16.6104,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ca;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['ketene(1375)', '[CH2][C]=CO[O](16807)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.560775,'m^3/(mol*s)'), n=2.01066, Ea=(45.2043,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [Cds-HH_Ck;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['O2(2)', 'C#CCC[C]=O(15744)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.3016e+08,'cm^3/(mol*s)'), n=1.628, Ea=(148.067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-Cs;OJ] for rate rule [Ct-H_Ct-Cs;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -1.9 to 148.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['O2(2)', '[CH]=[C]CC[C]=O(15746)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][C]=O(1376)', '[CH2][C]=CO[O](16807)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[O]O[CH][C]=CC[C]=O(19765)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['H(3)', '[O]O[CH][C]CC=C=O(20801)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction25',
    reactants = ['H(3)', '[O]O[C]=[C]CC[C]=O(21738)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['[O]O[CH]C=CC[C]=O(19767)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.66329e+10,'s^-1'), n=0.993, Ea=(157.679,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeC]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['[O]C=CC[C]=CO[O](21132)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(791180,'s^-1'), n=2.19286, Ea=(156.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;Cs_H_out_H/NonDeC] for rate rule [R2H_S;CO_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['[O]O[C]=CCC[C]=O(21739)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.231e+11,'s^-1'), n=0.765, Ea=(234.057,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 82 used for R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['O=[C]CC[C]=[C]OO(21740)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.58963e+06,'s^-1'), n=1.78157, Ea=(137.44,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_O;Y_rad_out;XH_out] for rate rule [R3H_SS_O;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['[O]OC=CC[CH][C]=O(10594)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.823e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_Cs;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['[O]O[CH][C]=CCC=O(19772)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(571202,'s^-1'), n=2.01653, Ea=(116.136,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;XH_out] for rate rule [R3H_SS_Cs;CO_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[O]O[C]=[C]CCC=O(21741)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_single;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleNd;CO_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[O]OC=[C]CC[C]=O(21049)'],
    products = ['O=[C]CC=[C][CH]OO(19769)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.19599e+09,'s^-1'), n=0.63, Ea=(50.8774,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;O_rad_out;Cs_H_out_H/NonDeC] for rate rule [R5HJ_3;O_rad_out;Cs_H_out_H/NonDeC]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['O=[C][CH]C[C]=COO(21742)'],
    products = ['[O]OC=[C]CC[C]=O(21049)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_2;Y_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4730',
    isomers = [
        '[O]OC=[C]CC[C]=O(21049)',
    ],
    reactants = [
        ('ketene(1375)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4730',
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


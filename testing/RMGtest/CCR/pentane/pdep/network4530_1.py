species(
    label = '[CH2][C]=CC=[C]CO[O](19509)',
    structure = SMILES('[CH2][C]=CC=[C]CO[O]'),
    E0 = (685.983,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.2086,'amu*angstrom^2'), symmetry=1, barrier=(4.79612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.799987,'amu*angstrom^2'), symmetry=1, barrier=(18.3933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.208521,'amu*angstrom^2'), symmetry=1, barrier=(4.79431,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.87133,'amu*angstrom^2'), symmetry=1, barrier=(66.0176,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.772443,0.0746308,-9.04976e-05,6.40734e-08,-1.87126e-11,82617.6,30.9806], Tmin=(100,'K'), Tmax=(829.3,'K')), NASAPolynomial(coeffs=[9.53949,0.0323443,-1.40115e-05,2.58676e-09,-1.76837e-13,81163.5,-9.6744], Tmin=(829.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(685.983,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C3H3(5450)',
    structure = SMILES('[CH]=C=C'),
    E0 = (338.572,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,540,610,2055,180,180,2603.58],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2095.25,'J/mol'), sigma=(4.76,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.09172,0.0173333,-8.20209e-06,-2.63358e-09,2.66049e-12,40755.9,8.10965], Tmin=(100,'K'), Tmax=(946.054,'K')), NASAPolynomial(coeffs=[6.98214,0.0072721,-2.37773e-06,3.99152e-10,-2.69331e-14,39733.9,-11.9544], Tmin=(946.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(338.572,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), label="""C3H3""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = 'C#CCO[O](16808)',
    structure = SMILES('C#CCO[O]'),
    E0 = (242.448,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2175,525,750,770,3400,2100],'cm^-1')),
        HinderedRotor(inertia=(0.595127,'amu*angstrom^2'), symmetry=1, barrier=(13.6831,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.74124,'amu*angstrom^2'), symmetry=1, barrier=(63.0265,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3504.08,'J/mol'), sigma=(5.7666,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=547.33 K, Pc=41.46 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.29019,0.0394598,-5.33344e-05,4.08815e-08,-1.22872e-11,29219.5,16.1946], Tmin=(100,'K'), Tmax=(942.252,'K')), NASAPolynomial(coeffs=[6.8169,0.0148149,-5.46001e-06,8.95137e-10,-5.57485e-14,28607.4,-4.09622], Tmin=(942.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(242.448,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ)"""),
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
    label = '[CH2][C]=CC=C=C(16910)',
    structure = SMILES('C=[C]C=C[C]=C'),
    E0 = (543.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62996,'amu*angstrom^2'), symmetry=1, barrier=(37.4761,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62989,'amu*angstrom^2'), symmetry=1, barrier=(37.4743,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759224,0.0600512,-5.80613e-05,2.91517e-08,-5.61572e-12,65516.4,20.7626], Tmin=(100,'K'), Tmax=(1428.68,'K')), NASAPolynomial(coeffs=[14.7312,0.0143871,-3.24523e-06,3.6595e-10,-1.74458e-14,62192.1,-49.2908], Tmin=(1428.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
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
    label = '[CH2][C]=CC=[C]C[O](20255)',
    structure = SMILES('[CH2][C]=CC=[C]C[O]'),
    E0 = (688.179,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,180,1028.72,2915.4],'cm^-1')),
        HinderedRotor(inertia=(2.99818,'amu*angstrom^2'), symmetry=1, barrier=(68.9341,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.365166,'amu*angstrom^2'), symmetry=1, barrier=(8.39589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.99637,'amu*angstrom^2'), symmetry=1, barrier=(68.8925,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.21167,0.0476121,-7.06651e-06,-8.02133e-08,8.32448e-11,82824.5,25.3167], Tmin=(100,'K'), Tmax=(459.893,'K')), NASAPolynomial(coeffs=[5.31653,0.0354756,-1.59772e-05,3.00322e-09,-2.07161e-13,82381.7,11.0399], Tmin=(459.893,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(688.179,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(CCOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2]O[O](46)',
    structure = SMILES('[CH2]O[O]'),
    E0 = (205.836,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2028.66,2028.67],'cm^-1')),
        HinderedRotor(inertia=(0.221535,'amu*angstrom^2'), symmetry=1, barrier=(13.5753,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (46.0254,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.43223,0.0119348,-9.20111e-06,4.32805e-09,-8.54519e-13,24777.3,10.782], Tmin=(100,'K'), Tmax=(1197.14,'K')), NASAPolynomial(coeffs=[5.14053,0.00622677,-2.04895e-06,3.45073e-10,-2.27385e-14,24368.2,2.23311], Tmin=(1197.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CsJOOH)"""),
)

species(
    label = '[C]=CC=[C][CH2](17519)',
    structure = SMILES('[C]C=C[C]=C'),
    E0 = (926.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.54616,'amu*angstrom^2'), symmetry=1, barrier=(35.5492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (64.0853,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60723,0.0446451,-4.44308e-05,2.23547e-08,-4.30798e-12,111509,16.4407], Tmin=(100,'K'), Tmax=(1401.68,'K')), NASAPolynomial(coeffs=[12.6996,0.00908062,-2.18752e-06,2.72848e-10,-1.45836e-14,108784,-39.4491], Tmin=(1401.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(926.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(CJ3)"""),
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
    label = '[CH][C]=CC=[C]CO[O](20256)',
    structure = SMILES('[CH]=C=C[CH][C]CO[O]'),
    E0 = (907.763,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.79104,'amu*angstrom^2'), symmetry=1, barrier=(64.1714,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.63723,'amu*angstrom^2'), symmetry=1, barrier=(14.6512,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.79163,'amu*angstrom^2'), symmetry=1, barrier=(64.185,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0142611,'amu*angstrom^2'), symmetry=1, barrier=(7.05605,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.444679,0.081834,-0.000110952,8.2316e-08,-2.44478e-11,109303,29.5973], Tmin=(100,'K'), Tmax=(824.722,'K')), NASAPolynomial(coeffs=[11.6845,0.0273155,-1.17869e-05,2.14996e-09,-1.45162e-13,107450,-22.4617], Tmin=(824.722,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(907.763,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = '[CH2][C]=CC=C1COO1(20257)',
    structure = SMILES('[CH2][C]=CC=C1COO1'),
    E0 = (462.509,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.37367,0.0429002,1.57818e-05,-5.06485e-08,2.1956e-11,55734.8,28.0759], Tmin=(100,'K'), Tmax=(995.596,'K')), NASAPolynomial(coeffs=[14.6881,0.0245805,-9.61075e-06,1.84013e-09,-1.34736e-13,51340.4,-44.8543], Tmin=(995.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.509,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(Cyclobutane) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=CC=C1CO[O](20258)',
    structure = SMILES('C=C1[CH]C=C1CO[O]'),
    E0 = (354.232,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13487,0.0452386,2.09511e-05,-6.84201e-08,3.21473e-11,42723.8,24.343], Tmin=(100,'K'), Tmax=(935.656,'K')), NASAPolynomial(coeffs=[18.6225,0.0155379,-3.67254e-06,5.95225e-10,-4.56318e-14,37478.9,-69.4019], Tmin=(935.656,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.232,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(ROOJ) + radical(C=CCJC=C)"""),
)

species(
    label = '[CH2][C]=CC=C1CO1(20259)',
    structure = SMILES('[CH2][C]=CC=C1CO1'),
    E0 = (395.929,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (94.1112,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2731,0.0419579,2.22235e-05,-7.1401e-08,3.4243e-11,47734.1,21.7335], Tmin=(100,'K'), Tmax=(922.065,'K')), NASAPolynomial(coeffs=[19.5748,0.00925883,-5.45196e-07,-1.66157e-11,-2.60852e-15,42374,-75.8404], Tmin=(922.065,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.929,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(299.321,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + ring(methyleneoxirane) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC=C=COO(20260)',
    structure = SMILES('C=[C]C=C[C]=COO'),
    E0 = (459.386,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.128974,0.0859231,-9.33703e-05,4.63799e-08,-7.13777e-12,55404.7,28.5482], Tmin=(100,'K'), Tmax=(894.365,'K')), NASAPolynomial(coeffs=[18.9122,0.0169787,-4.93638e-06,7.33612e-10,-4.52655e-14,51350.2,-64.8145], Tmin=(894.365,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(459.386,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CC=C=CO[O](19516)',
    structure = SMILES('[CH2]C=CC=C=CO[O]'),
    E0 = (384.237,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.56041,0.0673203,-5.15357e-05,1.2123e-08,2.28717e-12,46344.4,29.2244], Tmin=(100,'K'), Tmax=(962.344,'K')), NASAPolynomial(coeffs=[16.4345,0.0199421,-6.68354e-06,1.13851e-09,-7.75053e-14,42427.7,-51.225], Tmin=(962.344,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(384.237,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C]=CC=C=CO[O](20261)',
    structure = SMILES('C=[C]C=C[C]=CO[O]'),
    E0 = (611.391,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,2950,3100,1380,975,1025,1650,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24577,'amu*angstrom^2'), symmetry=1, barrier=(28.6427,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24848,'amu*angstrom^2'), symmetry=1, barrier=(28.7049,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24455,'amu*angstrom^2'), symmetry=1, barrier=(28.6146,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.21518,0.0812822,-9.38495e-05,4.92278e-08,-7.60482e-12,73671.7,28.0783], Tmin=(100,'K'), Tmax=(826.812,'K')), NASAPolynomial(coeffs=[17.11,0.0156997,-4.17285e-06,5.47965e-10,-2.97909e-14,70325.9,-53.5549], Tmin=(826.812,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(611.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CC#CCO[O](20262)',
    structure = SMILES('[CH2][C]=CC#CCO[O]'),
    E0 = (622.133,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2100,2250,500,550,1685,370,180,2618.85],'cm^-1')),
        HinderedRotor(inertia=(0.511451,'amu*angstrom^2'), symmetry=1, barrier=(11.7593,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.04485,'amu*angstrom^2'), symmetry=1, barrier=(70.0071,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.04515,'amu*angstrom^2'), symmetry=1, barrier=(70.0141,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.04313,'amu*angstrom^2'), symmetry=1, barrier=(69.9676,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.932797,0.0705943,-9.03307e-05,6.7782e-08,-2.06524e-11,74932.7,28.5263], Tmin=(100,'K'), Tmax=(863.69,'K')), NASAPolynomial(coeffs=[8.99002,0.0293928,-1.20254e-05,2.12986e-09,-1.41039e-13,73685.9,-8.32532], Tmin=(863.69,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(622.133,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CtOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(ROOJ) + radical(CTCC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C#CC=[C]CO[O](20263)',
    structure = SMILES('[CH2]C#CC=[C]CO[O]'),
    E0 = (657.123,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,2100,2250,500,550,215.287,215.497],'cm^-1')),
        HinderedRotor(inertia=(0.00364717,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.738596,'amu*angstrom^2'), symmetry=1, barrier=(24.7242,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266871,'amu*angstrom^2'), symmetry=1, barrier=(8.9375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95849,'amu*angstrom^2'), symmetry=1, barrier=(65.9609,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07828,0.0690691,-8.26768e-05,5.77838e-08,-1.69181e-11,79134.6,29.4228], Tmin=(100,'K'), Tmax=(820.618,'K')), NASAPolynomial(coeffs=[8.65669,0.0321286,-1.51529e-05,2.92703e-09,-2.05901e-13,77890.8,-5.64038], Tmin=(820.618,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(657.123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CtHHH) + group(Cds-CdsCsH) + group(Cds-CdsCtH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(ROOJ) + radical(Propargyl) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C][CH2](16918)',
    structure = SMILES('[CH][C]=C'),
    E0 = (614.492,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,260.76,263.05,263.653],'cm^-1')),
        HinderedRotor(inertia=(1.04394,'amu*angstrom^2'), symmetry=1, barrier=(50.8215,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.2703,0.0129731,8.09839e-06,-1.37802e-08,4.71189e-12,73935.1,11.2848], Tmin=(100,'K'), Tmax=(1104.68,'K')), NASAPolynomial(coeffs=[4.54119,0.0159815,-6.32008e-06,1.15742e-09,-7.99393e-14,73190,2.92522], Tmin=(1104.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(614.492,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CO[O](16805)',
    structure = SMILES('[CH]=[C]CO[O]'),
    E0 = (567.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,1685,370,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.348758,'amu*angstrom^2'), symmetry=1, barrier=(8.01863,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.349832,'amu*angstrom^2'), symmetry=1, barrier=(8.04334,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92348,0.0566508,-0.000115722,1.14457e-07,-4.10713e-11,68272.1,19.3685], Tmin=(100,'K'), Tmax=(892.364,'K')), NASAPolynomial(coeffs=[4.05419,0.0203884,-9.86661e-06,1.83106e-09,-1.20972e-13,68955.4,15.2908], Tmin=(892.364,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(567.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH2][C]=CC=[C][CH2](16912)',
    structure = SMILES('[CH2][C]=CC=[C][CH2]'),
    E0 = (733.665,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,383.297],'cm^-1')),
        HinderedRotor(inertia=(2.93244,'amu*angstrom^2'), symmetry=1, barrier=(67.4226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191132,'amu*angstrom^2'), symmetry=1, barrier=(67.4298,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188775,'amu*angstrom^2'), symmetry=1, barrier=(67.4008,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (78.1118,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3375.74,'J/mol'), sigma=(5.79513,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=527.28 K, Pc=39.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84928,0.0398838,-8.2793e-06,-1.72322e-08,9.68702e-12,88323.6,22.6963], Tmin=(100,'K'), Tmax=(960.212,'K')), NASAPolynomial(coeffs=[10.3621,0.0219425,-7.62276e-06,1.3151e-09,-8.95683e-14,85881.1,-22.2335], Tmin=(960.212,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(733.665,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC=[C][CH]O[O](20264)',
    structure = SMILES('[CH2][C]=CC=[C][CH]O[O]'),
    E0 = (803.278,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1670,1700,300,440,420.589,422.201],'cm^-1')),
        HinderedRotor(inertia=(0.0830286,'amu*angstrom^2'), symmetry=1, barrier=(10.4401,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.274071,'amu*angstrom^2'), symmetry=1, barrier=(33.8599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.695948,'amu*angstrom^2'), symmetry=1, barrier=(87.7453,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694577,'amu*angstrom^2'), symmetry=1, barrier=(87.8019,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.651519,0.0675285,-6.46018e-05,3.19264e-08,-6.25719e-12,96737.8,32.0051], Tmin=(100,'K'), Tmax=(1239.07,'K')), NASAPolynomial(coeffs=[15.1829,0.0206172,-7.811e-06,1.37044e-09,-9.20237e-14,93136.7,-41.2154], Tmin=(1239.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(803.278,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C]=[C]CO[O](20265)',
    structure = SMILES('[CH2][C]=C[C]=[C]CO[O]'),
    E0 = (884.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.298147,'amu*angstrom^2'), symmetry=1, barrier=(6.85499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298145,'amu*angstrom^2'), symmetry=1, barrier=(6.85495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.72045,'amu*angstrom^2'), symmetry=1, barrier=(62.5484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.7234,'amu*angstrom^2'), symmetry=1, barrier=(62.6163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547046,0.0830213,-0.000126558,1.06694e-07,-3.49441e-11,106556,31.3788], Tmin=(100,'K'), Tmax=(890.966,'K')), NASAPolynomial(coeffs=[8.55262,0.0314615,-1.34583e-05,2.39182e-09,-1.567e-13,105749,-2.84035], Tmin=(890.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(884.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=[C]C=[C]CO[O](20266)',
    structure = SMILES('[CH2][C]=[C]C=[C]CO[O]'),
    E0 = (884.979,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.298147,'amu*angstrom^2'), symmetry=1, barrier=(6.85499,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.298145,'amu*angstrom^2'), symmetry=1, barrier=(6.85495,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.72045,'amu*angstrom^2'), symmetry=1, barrier=(62.5484,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.7234,'amu*angstrom^2'), symmetry=1, barrier=(62.6163,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.547046,0.0830213,-0.000126558,1.06694e-07,-3.49441e-11,106556,31.3788], Tmin=(100,'K'), Tmax=(890.966,'K')), NASAPolynomial(coeffs=[8.55262,0.0314615,-1.34583e-05,2.39182e-09,-1.567e-13,105749,-2.84035], Tmin=(890.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(884.979,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC=C[CH]O[O](20267)',
    structure = SMILES('[CH2][C]=CC=C[CH]O[O]'),
    E0 = (565.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716207,0.062584,-3.87787e-05,1.62257e-09,4.99297e-12,68133.1,31.2978], Tmin=(100,'K'), Tmax=(993.797,'K')), NASAPolynomial(coeffs=[15.797,0.0220585,-8.06155e-06,1.44374e-09,-1.0069e-13,64139.4,-46.377], Tmin=(993.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C]=CCO[O](20268)',
    structure = SMILES('[CH2][C]=C[C]=CCO[O]'),
    E0 = (647.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612839,0.0780071,-0.000100186,7.51713e-08,-2.29266e-11,77951.2,30.6708], Tmin=(100,'K'), Tmax=(853.53,'K')), NASAPolynomial(coeffs=[9.62796,0.0321628,-1.32999e-05,2.37152e-09,-1.57787e-13,76543.3,-10.627], Tmin=(853.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C=[C]C=[C]CO[O](20269)',
    structure = SMILES('[CH2]C=[C]C=[C]CO[O]'),
    E0 = (647.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612839,0.0780071,-0.000100186,7.51713e-08,-2.29266e-11,77951.2,30.6708], Tmin=(100,'K'), Tmax=(853.53,'K')), NASAPolynomial(coeffs=[9.62796,0.0321628,-1.32999e-05,2.37152e-09,-1.57787e-13,76543.3,-10.627], Tmin=(853.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CC=[C][CH]OO(20270)',
    structure = SMILES('[CH2][C]=CC=[C][CH]OO'),
    E0 = (651.274,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.202457,0.0731861,-6.65584e-05,3.04356e-08,-5.47488e-12,78475.5,32.8653], Tmin=(100,'K'), Tmax=(1351.28,'K')), NASAPolynomial(coeffs=[18.1602,0.020028,-7.54944e-06,1.32281e-09,-8.86957e-14,73622.4,-59.1765], Tmin=(1351.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(651.274,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C=CCO[O](20271)',
    structure = SMILES('C=[C][C]=C[CH]CO[O]'),
    E0 = (613.84,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.114131,0.0943583,-0.00013226,1.00193e-07,-2.99963e-11,73972.4,28.0987], Tmin=(100,'K'), Tmax=(876.445,'K')), NASAPolynomial(coeffs=[13.0693,0.0291304,-1.19653e-05,2.10414e-09,-1.38063e-13,71855.8,-32.6565], Tmin=(876.445,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(613.84,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJCO) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=C[C]=[C]CO[O](20272)',
    structure = SMILES('[CH2]C=C[C]=[C]CO[O]'),
    E0 = (647.137,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.612839,0.0780071,-0.000100186,7.51713e-08,-2.29266e-11,77951.2,30.6708], Tmin=(100,'K'), Tmax=(853.53,'K')), NASAPolynomial(coeffs=[9.62796,0.0321628,-1.32999e-05,2.37152e-09,-1.57787e-13,76543.3,-10.627], Tmin=(853.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.137,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = 'C[C]=[C]C=[C]CO[O](20273)',
    structure = SMILES('CC#C[CH][C]CO[O]'),
    E0 = (741.728,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2100,2250,500,550,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.357307,0.0839091,-0.000116412,8.91684e-08,-2.65255e-11,89336.9,30.6273], Tmin=(100,'K'), Tmax=(955.326,'K')), NASAPolynomial(coeffs=[10.638,0.0285851,-1.02668e-05,1.64218e-09,-1.00046e-13,87932.9,-15.5687], Tmin=(955.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(741.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(ROOJ) + radical(Sec_Propargyl) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C[C]=C[C]=[C]CO[O](20274)',
    structure = SMILES('C[C][CH]C#CCO[O]'),
    E0 = (748.419,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2750,2800,2850,1350,1500,750,1050,1375,1000,3025,407.5,1350,352.5,2100,2250,500,550,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.505999,0.0790499,-0.000102926,7.5399e-08,-2.16786e-11,90137.8,31.1447], Tmin=(100,'K'), Tmax=(967.253,'K')), NASAPolynomial(coeffs=[10.654,0.0280427,-9.80455e-06,1.553e-09,-9.44926e-14,88597.6,-15.2897], Tmin=(967.253,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(748.419,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsHHH) + group(Cs-CtOsHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(ROOJ) + radical(Sec_Propargyl) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]=C[C]=[C]COO(20275)',
    structure = SMILES('[CH2][C]=C[C]=[C]COO'),
    E0 = (732.974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,220.026,220.486],'cm^-1')),
        HinderedRotor(inertia=(0.184623,'amu*angstrom^2'), symmetry=1, barrier=(6.3523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28,'amu*angstrom^2'), symmetry=1, barrier=(78.3386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.27684,'amu*angstrom^2'), symmetry=1, barrier=(78.3532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.27913,'amu*angstrom^2'), symmetry=1, barrier=(78.3423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.27761,'amu*angstrom^2'), symmetry=1, barrier=(78.3441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314036,0.0862275,-0.000120423,9.54687e-08,-3.03687e-11,88284.2,31.4577], Tmin=(100,'K'), Tmax=(845.951,'K')), NASAPolynomial(coeffs=[10.1277,0.0331474,-1.44652e-05,2.6365e-09,-1.77168e-13,86862.7,-12.8339], Tmin=(845.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CC=[C][CH]O[O](20276)',
    structure = SMILES('[CH2]C=CC=[C][CH]O[O]'),
    E0 = (565.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.716207,0.062584,-3.87787e-05,1.62257e-09,4.99297e-12,68133.1,31.2978], Tmin=(100,'K'), Tmax=(993.797,'K')), NASAPolynomial(coeffs=[15.797,0.0220585,-8.06155e-06,1.44374e-09,-1.0069e-13,64139.4,-46.377], Tmin=(993.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(565.437,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=CC=[C][CH]O[O](20277)',
    structure = SMILES('C[C]=CC=[C][CH]O[O]'),
    E0 = (685.223,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2800,2850,1350,1500,750,1050,1375,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1670,1700,300,440,267.782,267.822],'cm^-1')),
        HinderedRotor(inertia=(0.265252,'amu*angstrom^2'), symmetry=1, barrier=(13.5217,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265585,'amu*angstrom^2'), symmetry=1, barrier=(13.5144,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.265284,'amu*angstrom^2'), symmetry=1, barrier=(13.5161,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.499442,'amu*angstrom^2'), symmetry=1, barrier=(25.4091,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690246,0.0726495,-7.45796e-05,4.06731e-08,-8.97042e-12,82532.6,30.3111], Tmin=(100,'K'), Tmax=(1091.78,'K')), NASAPolynomial(coeffs=[13.2224,0.0267346,-1.14965e-05,2.15277e-09,-1.49805e-13,79796.2,-31.2495], Tmin=(1091.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(685.223,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=[C]C=[C]COO(20278)',
    structure = SMILES('[CH2][C]=[C]C=[C]COO'),
    E0 = (732.974,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1685,1700,300,370,440,220.026,220.486],'cm^-1')),
        HinderedRotor(inertia=(0.184623,'amu*angstrom^2'), symmetry=1, barrier=(6.3523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.28,'amu*angstrom^2'), symmetry=1, barrier=(78.3386,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.27684,'amu*angstrom^2'), symmetry=1, barrier=(78.3532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.27913,'amu*angstrom^2'), symmetry=1, barrier=(78.3423,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.27761,'amu*angstrom^2'), symmetry=1, barrier=(78.3441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.314036,0.0862275,-0.000120423,9.54687e-08,-3.03687e-11,88284.2,31.4577], Tmin=(100,'K'), Tmax=(845.951,'K')), NASAPolynomial(coeffs=[10.1277,0.0331474,-1.44652e-05,2.6365e-09,-1.77168e-13,86862.7,-12.8339], Tmin=(845.951,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(732.974,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1C=C1CO[O](20238)',
    structure = SMILES('C=[C]C1C=C1CO[O]'),
    E0 = (620.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.64474,0.0794308,-0.00011066,9.04662e-08,-3.00364e-11,74686.3,27.4327], Tmin=(100,'K'), Tmax=(817.263,'K')), NASAPolynomial(coeffs=[8.50534,0.034289,-1.5567e-05,2.91108e-09,-1.99085e-13,73624.2,-7.5413], Tmin=(817.263,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.018,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropene) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C=C=CCO[O](20279)',
    structure = SMILES('C=[C]C=C=CCO[O]'),
    E0 = (467.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.471178,0.0811226,-0.000104553,7.58665e-08,-2.22714e-11,56392.4,28.2047], Tmin=(100,'K'), Tmax=(831.177,'K')), NASAPolynomial(coeffs=[11.0264,0.0303274,-1.28874e-05,2.34573e-09,-1.58593e-13,54637.7,-20.7669], Tmin=(831.177,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(467.844,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=[C]CC=C=CO[O](20092)',
    structure = SMILES('C=[C]CC=C=CO[O]'),
    E0 = (534.877,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.607411,'amu*angstrom^2'), symmetry=1, barrier=(13.9656,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.607405,'amu*angstrom^2'), symmetry=1, barrier=(13.9654,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.608521,'amu*angstrom^2'), symmetry=1, barrier=(13.9911,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.763484,0.0725592,-7.87045e-05,4.63031e-08,-1.10556e-11,64446.3,30.4417], Tmin=(100,'K'), Tmax=(1011.03,'K')), NASAPolynomial(coeffs=[12.1543,0.0274929,-1.18427e-05,2.21488e-09,-1.53817e-13,62143,-24.6376], Tmin=(1011.03,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]C=[C]CO[O](20280)',
    structure = SMILES('[CH]=C=CC=[C]CO[O]'),
    E0 = (661.167,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,540,610,2055,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.792744,'amu*angstrom^2'), symmetry=1, barrier=(18.2267,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.796024,'amu*angstrom^2'), symmetry=1, barrier=(18.3022,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.793139,'amu*angstrom^2'), symmetry=1, barrier=(18.2358,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (109.103,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.438553,0.0813041,-0.000112382,8.38726e-08,-2.47668e-11,79645.4,29.8211], Tmin=(100,'K'), Tmax=(870.475,'K')), NASAPolynomial(coeffs=[12.1376,0.0247253,-1.00268e-05,1.76127e-09,-1.15748e-13,77715.5,-24.3835], Tmin=(870.475,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.167,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(295.164,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(Cds_S) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]C[C]=[C]CO[O](20095)',
    structure = SMILES('C=[C]C[C]=[C]CO[O]'),
    E0 = (836.623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.237036,'amu*angstrom^2'), symmetry=1, barrier=(5.44993,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236946,'amu*angstrom^2'), symmetry=1, barrier=(5.44785,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.237041,'amu*angstrom^2'), symmetry=1, barrier=(5.45003,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.236793,'amu*angstrom^2'), symmetry=1, barrier=(5.44434,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.59633,0.0845956,-0.000135462,1.23326e-07,-4.37326e-11,100736,33.5417], Tmin=(100,'K'), Tmax=(849.178,'K')), NASAPolynomial(coeffs=[6.15037,0.0383385,-1.82557e-05,3.44404e-09,-2.34792e-13,100517,11.9208], Tmin=(849.178,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(836.623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH]C=[C]CO[O](20281)',
    structure = SMILES('[CH]C=CC=[C]CO[O]'),
    E0 = (700.77,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.650906,0.0773073,-8.28202e-05,5.22105e-08,-1.38745e-11,84400.6,30.4987], Tmin=(100,'K'), Tmax=(898.626,'K')), NASAPolynomial(coeffs=[9.55849,0.0376578,-1.6637e-05,3.11121e-09,-2.15066e-13,82799.7,-11.5232], Tmin=(898.626,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.77,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + radical(ROOJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]CC=[C]CO[O](19500)',
    structure = SMILES('[CH]=[C]CC=[C]CO[O]'),
    E0 = (845.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.307449,'amu*angstrom^2'), symmetry=1, barrier=(7.06887,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307457,'amu*angstrom^2'), symmetry=1, barrier=(7.06903,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307414,'amu*angstrom^2'), symmetry=1, barrier=(7.06805,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.307429,'amu*angstrom^2'), symmetry=1, barrier=(7.06841,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.563494,0.0834111,-0.000126883,1.10744e-07,-3.83322e-11,101852,33.5899], Tmin=(100,'K'), Tmax=(835.565,'K')), NASAPolynomial(coeffs=[7.56703,0.0360318,-1.69609e-05,3.20054e-09,-2.1897e-13,101165,3.95365], Tmin=(835.565,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(845.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]CC=[C][CH]O[O](20094)',
    structure = SMILES('C=[C]CC=[C][CH]O[O]'),
    E0 = (716.077,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,377.392,377.404,377.425],'cm^-1')),
        HinderedRotor(inertia=(0.00118376,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106928,'amu*angstrom^2'), symmetry=1, barrier=(10.8053,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106925,'amu*angstrom^2'), symmetry=1, barrier=(10.8054,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.30792,'amu*angstrom^2'), symmetry=1, barrier=(31.1202,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00777,0.0668246,-6.26284e-05,3.16609e-08,-6.6203e-12,86231,32.1941], Tmin=(100,'K'), Tmax=(1130.98,'K')), NASAPolynomial(coeffs=[11.5519,0.0295327,-1.31687e-05,2.50642e-09,-1.75776e-13,83846,-19.9728], Tmin=(1130.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(716.077,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C][CH]C=CCO[O](20282)',
    structure = SMILES('[CH]=C=C[CH][CH]CO[O]'),
    E0 = (659.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.520326,0.0808472,-0.000107607,8.23704e-08,-2.56537e-11,79496.5,31.7242], Tmin=(100,'K'), Tmax=(812.522,'K')), NASAPolynomial(coeffs=[10.0065,0.0321552,-1.40386e-05,2.58054e-09,-1.75063e-13,78020.7,-11.6668], Tmin=(812.522,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(659.963,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCJCOOH) + radical(Allyl_S) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C][CH]C=[C]COO(20283)',
    structure = SMILES('[CH]=C=C[CH][C]COO'),
    E0 = (755.759,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,540,610,2055,3120,650,792.5,1650,200,800,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.317484,0.0837264,-9.98763e-05,6.40855e-08,-1.65654e-11,91027.1,29.3008], Tmin=(100,'K'), Tmax=(940.066,'K')), NASAPolynomial(coeffs=[13.1243,0.0292344,-1.29292e-05,2.42673e-09,-1.68299e-13,88619.2,-31.6936], Tmin=(940.066,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(755.759,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(C=C=CJ)"""),
)

species(
    label = 'H2CC(T)(1341)',
    structure = SMILES('[C]=C'),
    E0 = (600.324,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650],'cm^-1')),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (26.0373,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.9483,-0.00181651,2.07124e-05,-2.3669e-08,8.45456e-12,72206.8,5.3155], Tmin=(100,'K'), Tmax=(953.668,'K')), NASAPolynomial(coeffs=[4.20274,0.00473419,-1.57302e-06,2.859e-10,-2.08638e-14,71811.9,2.2838], Tmin=(953.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.324,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""H2CC(T)""", comment="""Thermo library: DFT_QCI_thermo"""),
)

species(
    label = '[CH]=C[C]CO[O](20284)',
    structure = SMILES('[CH]C=[C]CO[O]'),
    E0 = (648.998,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1685,370,399.576,399.925,400.494,401.646],'cm^-1')),
        HinderedRotor(inertia=(0.474613,'amu*angstrom^2'), symmetry=1, barrier=(53.9083,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.473153,'amu*angstrom^2'), symmetry=1, barrier=(53.9369,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.472646,'amu*angstrom^2'), symmetry=1, barrier=(53.9076,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (84.0734,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72176,0.0556679,-7.58851e-05,6.70209e-08,-2.4098e-11,78133.1,24.007], Tmin=(100,'K'), Tmax=(835.497,'K')), NASAPolynomial(coeffs=[4.22008,0.0327548,-1.50853e-05,2.81734e-09,-1.92049e-13,78097.9,14.6908], Tmin=(835.497,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.998,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C1[CH]C1[C]CO[O](20285)',
    structure = SMILES('C=C1[CH]C1[C]CO[O]'),
    E0 = (784.422,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.455386,0.0688044,-5.99237e-05,2.65638e-08,-4.67159e-12,94479.4,27.4801], Tmin=(100,'K'), Tmax=(1372.59,'K')), NASAPolynomial(coeffs=[16.693,0.0214858,-8.21395e-06,1.44894e-09,-9.73252e-14,90021.8,-55.9996], Tmin=(1372.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(784.422,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(Allyl_S) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=[C]C1[CH][C]COO1(20286)',
    structure = SMILES('C=[C]C1[CH][C]COO1'),
    E0 = (797.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.02718,0.0494435,1.23514e-05,-6.3744e-08,3.25711e-11,96010.1,27.4919], Tmin=(100,'K'), Tmax=(899.341,'K')), NASAPolynomial(coeffs=[18.5361,0.0147316,-1.74309e-06,6.89771e-11,-2.22124e-15,91115.3,-64.8252], Tmin=(899.341,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(797.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxane) + radical(CCJCOOH) + radical(CCJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = 'C=[C][CH]C1[C]COO1(20287)',
    structure = SMILES('C=[C][CH]C1[C]COO1'),
    E0 = (722.352,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.615775,0.0609169,-1.48427e-05,-3.7786e-08,2.38036e-11,87013.4,24.8493], Tmin=(100,'K'), Tmax=(894.209,'K')), NASAPolynomial(coeffs=[19.132,0.0158232,-2.49673e-06,1.99858e-10,-9.60238e-15,82193.4,-70.846], Tmin=(894.209,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(722.352,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(324.264,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxolane) + radical(C=CCJCO) + radical(CCJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH][C]CO[O](20242)',
    structure = SMILES('[CH][C]CO[O]'),
    E0 = (847.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,404.609,404.615,404.632,3167.93,3167.98],'cm^-1')),
        HinderedRotor(inertia=(0.0766569,'amu*angstrom^2'), symmetry=1, barrier=(8.89488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.622243,'amu*angstrom^2'), symmetry=1, barrier=(72.2622,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0765029,'amu*angstrom^2'), symmetry=1, barrier=(8.89488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01219,0.0466378,-7.05672e-05,5.68391e-08,-1.80045e-11,102004,20.8821], Tmin=(100,'K'), Tmax=(805.913,'K')), NASAPolynomial(coeffs=[8.52167,0.0130199,-5.55933e-06,1.04759e-09,-7.2263e-14,100998,-8.85402], Tmin=(805.913,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(291.007,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C=[C][CH][C]CO[O](20288)',
    structure = SMILES('[CH2]C=[C][CH][C]CO[O]'),
    E0 = (966.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.543462,0.0809997,-0.000101607,7.20591e-08,-2.09935e-11,116306,30.3257], Tmin=(100,'K'), Tmax=(831.055,'K')), NASAPolynomial(coeffs=[10.4596,0.0332726,-1.54644e-05,2.95731e-09,-2.06552e-13,114658,-15.6791], Tmin=(831.055,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(966.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Allyl_S) + radical(Allyl_P) + radical(CCJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH][CH][C]CO[O](20289)',
    structure = SMILES('[CH]C=C[CH][C]CO[O]'),
    E0 = (947.367,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.676939,0.0775936,-8.04897e-05,4.93955e-08,-1.2967e-11,114058,30.2041], Tmin=(100,'K'), Tmax=(901.795,'K')), NASAPolynomial(coeffs=[9.09176,0.0402699,-1.84091e-05,3.50261e-09,-2.44697e-13,112540,-9.52295], Tmin=(901.795,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(947.367,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(ROOJ) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[O]OC[C]C1C=[C]C1(20290)',
    structure = SMILES('[O]OC[C]C1C=[C]C1'),
    E0 = (853.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.759389,0.0656084,-5.59062e-05,2.46202e-08,-4.3682e-12,102752,29.0886], Tmin=(100,'K'), Tmax=(1343.72,'K')), NASAPolynomial(coeffs=[14.5984,0.0244121,-9.91811e-06,1.80366e-09,-1.23138e-13,99032.7,-41.7649], Tmin=(1343.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(853.318,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(320.107,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(ROOJ) + radical(CCJ2_triplet) + radical(cyclobutene-vinyl)"""),
)

species(
    label = 'C=C1[CH][CH][C]COO1(20291)',
    structure = SMILES('[CH2]C1=C[CH][C]COO1'),
    E0 = (671.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2426,0.00905961,0.000189561,-2.94413e-07,1.26227e-10,80965.4,26.4713], Tmin=(100,'K'), Tmax=(904.672,'K')), NASAPolynomial(coeffs=[39.3365,-0.0239976,1.99129e-05,-3.98996e-09,2.6133e-13,68533.1,-184.111], Tmin=(904.672,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(671.956,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(328.422,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cycloheptane) + radical(Allyl_S) + radical(C=C(O)CJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[C]1[CH]C=[C]COOC1(20292)',
    structure = SMILES('[C]1[CH]C=[C]COOC1'),
    E0 = (774.944,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92643,0.0299765,4.56175e-05,-7.18691e-08,2.65659e-11,93292.5,21.0649], Tmin=(100,'K'), Tmax=(1034.64,'K')), NASAPolynomial(coeffs=[11.2585,0.032567,-1.42004e-05,2.79825e-09,-2.04611e-13,89291.7,-34.277], Tmin=(1034.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(774.944,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(332.579,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclooctane) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = 'C=C=[C]C[C]CO[O](20293)',
    structure = SMILES('[CH2]C#CC[C]CO[O]'),
    E0 = (751.925,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3000,3100,440,815,1455,1000,2100,2250,500,550,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.474392,0.0817772,-0.000110266,8.48579e-08,-2.6374e-11,90558.6,31.4289], Tmin=(100,'K'), Tmax=(837.733,'K')), NASAPolynomial(coeffs=[10.1945,0.0315598,-1.35355e-05,2.45668e-09,-1.65073e-13,89063.5,-12.947], Tmin=(837.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(751.925,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CtCsHH) + group(Cs-CsOsHH) + group(Cs-CtHHH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(ROOJ) + radical(Propargyl) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C=CC[C][CH]O[O](20294)',
    structure = SMILES('C=C=CC[C][CH]O[O]'),
    E0 = (800.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,540,610,2055,239.46,240.144,241.858,242.671],'cm^-1')),
        HinderedRotor(inertia=(1.82485,'amu*angstrom^2'), symmetry=1, barrier=(75.8944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188554,'amu*angstrom^2'), symmetry=1, barrier=(7.8712,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.694696,'amu*angstrom^2'), symmetry=1, barrier=(28.4667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.191523,'amu*angstrom^2'), symmetry=1, barrier=(7.87235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.25579,0.0886339,-0.000128569,1.05368e-07,-3.48359e-11,96437.6,30.902], Tmin=(100,'K'), Tmax=(803.661,'K')), NASAPolynomial(coeffs=[10.0366,0.0341497,-1.60459e-05,3.04157e-09,-2.09491e-13,95052.9,-12.981], Tmin=(803.661,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(800.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(315.95,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCsJOOH) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=C=CC[C]CO[O](20295)',
    structure = SMILES('C#C[CH]C[C]CO[O]'),
    E0 = (760.129,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,3025,407.5,1350,352.5,2175,525,750,770,3400,2100,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (110.111,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0818587,0.0910641,-0.000132846,1.04316e-07,-3.18507e-11,91558.8,30.4411], Tmin=(100,'K'), Tmax=(916.75,'K')), NASAPolynomial(coeffs=[11.7607,0.0283984,-1.11553e-05,1.89096e-09,-1.20144e-13,89909.5,-22.2043], Tmin=(916.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(760.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(311.793,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cs-CtCsHH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(Sec_Propargyl) + radical(CCJ2_triplet)"""),
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
    E0 = (685.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (685.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (931.183,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1187.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1119.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (694.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (694.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (817.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (764.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (710.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (835.608,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (852.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (884.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (880.317,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (928.367,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (685.983,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (725.043,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (1015.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (1181.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (1096.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (1096.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (825.338,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (882.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (835.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (850.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (836.329,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (879.362,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (934.806,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (917.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (959.379,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (750.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (827.757,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (807.418,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (688.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (708.284,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (710.956,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (888.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1019.19,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (855.638,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (991.063,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (812.727,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (867.878,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (907.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (1304.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (823.633,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (797.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (776.06,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (1241.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (988.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (970.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (853.318,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (748.154,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (813.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (929.611,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (936.737,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (1075.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['C3H3(5450)', 'C#CCO[O](16808)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['O2(S)(666)', '[CH2][C]=CC=C=C(16910)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH2][C]=CC=[C]C[O](20255)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]O[O](46)', '[C]=CC=[C][CH2](17519)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH][C]=CC=[C]CO[O](20256)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2][C]=CC=C1COO1(20257)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2]C1=CC=C1CO[O](20258)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleND_rad_out;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['O(4)', '[CH2][C]=CC=C1CO1(20259)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO_S;Y_rad_intra;OOJ] for rate rule [R2OO_S;Cd_rad_out;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2][C]=CC=C=COO(20260)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2]C=CC=C=CO[O](19516)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['H(3)', '[CH2][C]=CC=C=CO[O](20261)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction12',
    reactants = ['H(3)', '[CH2][C]=CC#CCO[O](20262)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(7.63e+08,'cm^3/(mol*s)'), n=1.64, Ea=(18.9954,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""From training reaction 2706 used for Ct-Cd_Ct-Cs;HJ
Exact match found for rate rule [Ct-Cd_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[CH2]C#CC=[C]CO[O](20263)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C][CH2](16918)', 'C#CCO[O](16808)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(0.109931,'m^3/(mol*s)'), n=2.3439, Ea=(23.3771,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct-H_Ct-Cs;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C3H3(5450)', '[CH]=[C]CO[O](16805)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O2(2)', '[CH2][C]=CC=C=C(16910)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(3.946e+06,'cm^3/(mol*s)'), n=2.037, Ea=(150.921,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Ca;OJ] for rate rule [Cds-HH_Ca;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: R_Addition_MultipleBond
Ea raised from -6.0 to 150.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', '[CH2][C]=CC=[C][CH2](16912)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(2.4339e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2][C]=CC=[C][CH]O[O](20264)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C][CH2](16918)', '[CH]=[C]CO[O](16805)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.63222,'m^3/(mol*s)'), n=1.59841, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Sp-4R!H-2R_N-Sp-5R!H-4R!H_N-4R!H-inRing]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH2][C]=C[C]=[C]CO[O](20265)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', '[CH2][C]=[C]C=[C]CO[O](20266)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.52804e+49,'m^3/(mol*s)'), n=-12.7885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O_Ext-2CNO-R_Ext-3CS-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2][C]=CC=C[CH]O[O](20267)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2][C]=C[C]=CCO[O](20268)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.75e+11,'s^-1'), n=0.633, Ea=(196.46,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2]C=[C]C=[C]CO[O](20269)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.44274e+08,'s^-1'), n=1.26608, Ea=(149.637,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2][C]=CC=[C][CH]OO(20270)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(5.66219e+08,'s^-1'), n=1.32333, Ea=(164.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2][C]=[C]C=CCO[O](20271)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2]C=C[C]=[C]CO[O](20272)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[C]=[C]C=[C]CO[O](20273)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C[C]=C[C]=[C]CO[O](20274)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(115.297,'s^-1'), n=2.99825, Ea=(169.04,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R4HJ_2;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C]=C[C]=[C]COO(20275)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.61756e+09,'s^-1'), n=1.29078, Ea=(226.405,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH2]C=CC=[C][CH]O[O](20276)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(26847.4,'s^-1'), n=2.09453, Ea=(64.7206,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5HJ_3;Cd_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C[C]=CC=[C][CH]O[O](20277)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(28.0827,'s^-1'), n=3.4, Ea=(142.535,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;C_rad_out_H/NonDeO;Cs_H_out_2H] for rate rule [R6Hall;C_rad_out_H/NonDeO;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH2][C]=[C]C=[C]COO(20278)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(19101.7,'s^-1'), n=1.79759, Ea=(74.4436,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6Hall;Y_rad_out;XH_out] for rate rule [R6HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.7320508075688772
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['C=[C]C1C=C1CO[O](20238)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriND_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['C=[C]C=C=CCO[O](20279)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['C=[C]CC=C=CO[O](20092)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction37',
    reactants = ['H(3)', 'C#C[CH]C=[C]CO[O](20280)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=[C]C[C]=[C]CO[O](20095)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(7.28974e+09,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=C[CH]C=[C]CO[O](20281)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.3587e+07,'s^-1'), n=1.74167, Ea=(154.868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC=[C]CO[O](19500)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.846e+10,'s^-1'), n=0.74, Ea=(145.185,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleH;Cs_H_out_1H] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/OneDe]
Euclidian distance = 1.4142135623730951
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=[C]CC=[C][CH]O[O](20094)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(55800,'s^-1'), n=1.97, Ea=(96.6504,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;C_rad_out_H/NonDeO;XH_out] for rate rule [R4HJ_1;C_rad_out_H/NonDeO;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[CH]=[C][CH]C=CCO[O](20282)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.74e+08,'s^-1'), n=1.713, Ea=(181.895,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_Cs;Cd_H_out_singleH] for rate rule [R5Hall;Cd_rad_out_Cs;Cd_H_out_singleH]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=[C][CH]C=[C]COO(20283)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R8Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['H2CC(T)(1341)', '[CH]=C[C]CO[O](20284)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['C=C1[CH]C1[C]CO[O](20285)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(5.605e+12,'s^-1'), n=0.275, Ea=(137.649,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra_pri;radadd_intra_cddouble]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['C=[C]C1[CH][C]COO1(20286)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(4.18543e+10,'s^-1'), n=0.209443, Ea=(111.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 108.3 to 111.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['C=[C][CH]C1[C]COO1(20287)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.34604e+09,'s^-1'), n=0.547362, Ea=(90.0766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra;radadd_intra_O] + [R6;doublebond_intra;radadd_intra] for rate rule [R6;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C3H3(5450)', '[CH][C]CO[O](20242)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH2]C=[C][CH][C]CO[O](20288)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=C[CH][CH][C]CO[O](20289)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[O]OC[C]C1C=[C]C1(20290)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(3.01156e+10,'s^-1'), n=0.428741, Ea=(167.335,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_CdCdd;radadd_intra_cs] + [R4;doublebond_intra;radadd_intra_cs] for rate rule [R4;doublebond_intra_CdCdd;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 166.0 to 167.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['C=C1[CH][CH][C]COO1(20291)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(8.62196e+06,'s^-1'), n=0.867572, Ea=(62.1704,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra;radadd_intra] for rate rule [R7_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH2][C]=CC=[C]CO[O](19509)'],
    products = ['[C]1[CH]C=[C]COOC1(20292)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(3.7039e+11,'s^-1'), n=-0.18575, Ea=(127.616,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;doublebond_intra_CdCdd;radadd_intra] for rate rule [R8;doublebond_intra_CdCdd;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction54',
    reactants = ['C=C=[C]C[C]CO[O](20293)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(7.80481e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C=C=CC[C][CH]O[O](20294)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.60984e+10,'s^-1'), n=0.8, Ea=(135.98,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_1H;Cs_H_out_H/Cd] for rate rule [R3Hall;C_rad_out_H/NonDeO;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]=C=CC[C]CO[O](20295)'],
    products = ['[CH2][C]=CC=[C]CO[O](19509)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(6586.33,'s^-1'), n=2.5965, Ea=(315.162,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_MMS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4530',
    isomers = [
        '[CH2][C]=CC=[C]CO[O](19509)',
    ],
    reactants = [
        ('C3H3(5450)', 'C#CCO[O](16808)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4530',
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


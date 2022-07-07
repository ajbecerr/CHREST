species(
    label = '[O][C]=CC=[C]CO[O](19437)',
    structure = SMILES('[O]OC[C]=C[CH][C]=O'),
    E0 = (501.141,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,1855,455,950,180,922.056],'cm^-1')),
        HinderedRotor(inertia=(0.00403648,'amu*angstrom^2'), symmetry=1, barrier=(2.44874,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.776933,'amu*angstrom^2'), symmetry=1, barrier=(17.8632,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777553,'amu*angstrom^2'), symmetry=1, barrier=(17.8775,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.777271,'amu*angstrom^2'), symmetry=1, barrier=(17.871,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52493,0.0583923,-5.28085e-05,2.42202e-08,-4.62241e-12,60359,29.311], Tmin=(100,'K'), Tmax=(1210.86,'K')), NASAPolynomial(coeffs=[11.0236,0.027014,-1.39371e-05,2.81857e-09,-2.03715e-13,58058.7,-18.3319], Tmin=(1210.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(501.141,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJC=O) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'HCCO(2227)',
    structure = SMILES('[CH]=C=O'),
    E0 = (165.191,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2120,512.5,787.5,231.114,1089.61,3388.99],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(1247.18,'J/mol'), sigma=(2.5,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.47174,0.0113766,-1.1308e-05,6.13951e-09,-1.35087e-12,19887.1,6.87336], Tmin=(100,'K'), Tmax=(1096.2,'K')), NASAPolynomial(coeffs=[5.39217,0.00436893,-1.71893e-06,3.07751e-10,-2.08706e-14,19466,-2.56797], Tmin=(1096.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.191,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(83.1447,'J/(mol*K)'), label="""HCCO""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = 'C=[C]C=C[C]=O(16163)',
    structure = SMILES('C=[C]C=C[C]=O'),
    E0 = (330.05,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(1.23543,'amu*angstrom^2'), symmetry=1, barrier=(28.4051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.23355,'amu*angstrom^2'), symmetry=1, barrier=(28.3617,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43864,0.0597968,-7.86998e-05,5.54009e-08,-1.56855e-11,39785,19.3972], Tmin=(100,'K'), Tmax=(860.663,'K')), NASAPolynomial(coeffs=[9.97258,0.0201332,-9.56998e-06,1.85128e-09,-1.30193e-13,38316.1,-20.4933], Tmin=(860.663,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.05,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CCJ=O)"""),
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
    label = '[O]C[C]=C[CH][C]=O(20774)',
    structure = SMILES('[O]C[C]=C[CH][C]=O'),
    E0 = (503.336,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,1855,455,950,180,180,1174.82],'cm^-1')),
        HinderedRotor(inertia=(0.919915,'amu*angstrom^2'), symmetry=1, barrier=(21.1507,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0018862,'amu*angstrom^2'), symmetry=1, barrier=(1.84782,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.021595,'amu*angstrom^2'), symmetry=1, barrier=(21.1492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44419,0.0409021,-2.51118e-05,6.50836e-09,-6.4754e-13,60585.9,25.3213], Tmin=(100,'K'), Tmax=(2246.95,'K')), NASAPolynomial(coeffs=[15.9029,0.0169427,-9.11698e-06,1.76264e-09,-1.19513e-13,54537.8,-50.5045], Tmin=(2246.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(503.336,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJC=O) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    label = '[C]=C[CH][C]=O(15285)',
    structure = SMILES('[C][CH]C=C=O'),
    E0 = (708.907,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,2120,512.5,787.5,461.81,462.577],'cm^-1')),
        HinderedRotor(inertia=(0.231582,'amu*angstrom^2'), symmetry=1, barrier=(35.1116,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (66.0581,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.44707,0.0297002,-2.26579e-05,6.17966e-09,1.60518e-13,85321.6,13.7358], Tmin=(100,'K'), Tmax=(1067.82,'K')), NASAPolynomial(coeffs=[9.9549,0.00886206,-3.62032e-06,6.83637e-10,-4.86811e-14,83302.8,-24.9227], Tmin=(1067.82,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(708.907,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(Allyl_S) + radical(CJ3)"""),
)

species(
    label = '[CH]=[C][O](6861)',
    structure = SMILES('[CH][C]=O'),
    E0 = (425.043,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1855,455,950,1491.4,1492.35],'cm^-1')),
        HinderedRotor(inertia=(0.0736325,'amu*angstrom^2'), symmetry=1, barrier=(5.08901,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (41.0287,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.62981,0.00805288,-4.35162e-06,1.03741e-09,-8.43183e-14,51134.2,10.4631], Tmin=(100,'K'), Tmax=(1941.7,'K')), NASAPolynomial(coeffs=[6.62042,0.00292938,-1.19496e-06,2.28726e-10,-1.56225e-14,49777.3,-6.4529], Tmin=(1941.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.043,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCJ2_triplet) + radical(CsCJ=O)"""),
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
    label = '[O]OC[C]=C[C][C]=O(20775)',
    structure = SMILES('[O]OC[C][CH][C]=C=O'),
    E0 = (736.017,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1685,370,2120,512.5,787.5,190.386,190.403,190.41,3178.68],'cm^-1')),
        HinderedRotor(inertia=(1.79904,'amu*angstrom^2'), symmetry=1, barrier=(46.2623,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.514341,'amu*angstrom^2'), symmetry=1, barrier=(13.232,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79867,'amu*angstrom^2'), symmetry=1, barrier=(46.2626,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.79844,'amu*angstrom^2'), symmetry=1, barrier=(46.2623,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.424069,0.0879481,-0.000151558,1.34934e-07,-4.61502e-11,88642.2,30.4679], Tmin=(100,'K'), Tmax=(855.515,'K')), NASAPolynomial(coeffs=[9.38069,0.0277557,-1.39077e-05,2.64524e-09,-1.79987e-13,87780,-7.42772], Tmin=(855.515,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(736.017,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(CCJC(C)=C=O) + radical(CCJ2_triplet) + radical(CCCJ=C=O)"""),
)

species(
    label = '[O]OCC1=CC1[C]=O(20756)',
    structure = SMILES('[O]OCC1=CC1[C]=O'),
    E0 = (356.988,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22362,0.0644014,-7.67695e-05,5.14376e-08,-1.41422e-11,43032.8,26.907], Tmin=(100,'K'), Tmax=(879.18,'K')), NASAPolynomial(coeffs=[9.57788,0.0263906,-1.19156e-05,2.25814e-09,-1.57257e-13,41563.9,-12.3215], Tmin=(879.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(356.988,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + ring(Cyclopropene) + radical(ROOJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'O=[C][CH]C=C1COO1(20776)',
    structure = SMILES('O=C=C[CH][C]1COO1'),
    E0 = (246.466,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.74049,0.059655,-3.8861e-05,-4.18234e-08,6.27255e-11,29714.2,18.9425], Tmin=(100,'K'), Tmax=(475.857,'K')), NASAPolynomial(coeffs=[6.53901,0.036189,-1.80683e-05,3.54652e-09,-2.50427e-13,29066.5,-2.65089], Tmin=(475.857,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(246.466,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + missing(O2d-Cdd) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + ring(12dioxetane) + radical(C2CsJOO) + radical(C=CCJCO)"""),
)

species(
    label = 'O=[C][CH]C=C1CO1(20777)',
    structure = SMILES('O=C=C[CH][C]1CO1'),
    E0 = (188.651,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04012,0.0621522,-7.25032e-05,4.59089e-08,-1.12311e-11,22798.6,21.3876], Tmin=(100,'K'), Tmax=(1130.9,'K')), NASAPolynomial(coeffs=[11.6437,0.0178189,-4.64342e-06,5.66352e-10,-2.72582e-14,20837,-29.1422], Tmin=(1130.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(188.651,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + missing(O2d-Cdd) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + ring(Ethylene_oxide) + radical(C2CsJO) + radical(C=CCJCO)"""),
)

species(
    label = '[O]OCC=C=C[C]=O(20778)',
    structure = SMILES('[O]OCC=C=C[C]=O'),
    E0 = (248.259,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.85975,0.0781439,-0.00012852,1.20377e-07,-4.46506e-11,29962.9,28.3235], Tmin=(100,'K'), Tmax=(797.715,'K')), NASAPolynomial(coeffs=[6.14376,0.0357402,-1.88727e-05,3.7435e-09,-2.63673e-13,29626.1,7.19792], Tmin=(797.715,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(248.259,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CsCJ=O)"""),
)

species(
    label = 'O=[C][CH]C=C=COO(20779)',
    structure = SMILES('O=[C]C=C[C]=COO'),
    E0 = (245.752,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0609058,0.091466,-0.000134446,9.9471e-08,-2.88578e-11,29694.7,28.2445], Tmin=(100,'K'), Tmax=(848.377,'K')), NASAPolynomial(coeffs=[14.9385,0.0213216,-1.04277e-05,2.01838e-09,-1.41169e-13,27170.3,-41.0852], Tmin=(848.377,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.752,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CCJ=O)"""),
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
    label = '[O]OC=C=C[CH][C]=O(20780)',
    structure = SMILES('[O]OC=[C]C=C[C]=O'),
    E0 = (397.757,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,1685,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.976916,'amu*angstrom^2'), symmetry=1, barrier=(22.4612,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.978679,'amu*angstrom^2'), symmetry=1, barrier=(22.5018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.982736,'amu*angstrom^2'), symmetry=1, barrier=(22.595,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.178594,0.089743,-0.000146387,1.19204e-07,-3.75451e-11,47971.3,28.5714], Tmin=(100,'K'), Tmax=(853.443,'K')), NASAPolynomial(coeffs=[13.6487,0.0191257,-9.11612e-06,1.69987e-09,-1.14459e-13,45944.7,-32.6824], Tmin=(853.443,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.757,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cd-Cd(CO)H) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CJC=C) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]OC[C]=C=C[C]=O(20781)',
    structure = SMILES('[O]OCC#C[CH][C]=O'),
    E0 = (417.022,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,2100,2250,500,550,1855,455,950,308.882,1047.26],'cm^-1')),
        HinderedRotor(inertia=(4.83258,'amu*angstrom^2'), symmetry=1, barrier=(111.111,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0220208,'amu*angstrom^2'), symmetry=1, barrier=(17.0428,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0218719,'amu*angstrom^2'), symmetry=1, barrier=(17.0268,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(4.83438,'amu*angstrom^2'), symmetry=1, barrier=(111.152,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29483,0.0524273,-4.17367e-05,1.32695e-08,-5.02398e-13,50259.6,30.1445], Tmin=(100,'K'), Tmax=(1052.15,'K')), NASAPolynomial(coeffs=[13.8577,0.0161155,-6.2908e-06,1.15225e-09,-8.05702e-14,46982.3,-34.1143], Tmin=(1052.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(417.022,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CtHH) + group(Cs-CtOsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtCs) + radical(ROOJ) + radical(CCJC=O) + radical(CCCJ=O)"""),
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
    label = 'C=[C][CH][CH][C]=O(16165)',
    structure = SMILES('[CH2][C]=C[CH][C]=O'),
    E0 = (582.266,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,1855,455,950,612.29],'cm^-1')),
        HinderedRotor(inertia=(0.095959,'amu*angstrom^2'), symmetry=1, barrier=(25.5686,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0959819,'amu*angstrom^2'), symmetry=1, barrier=(25.5409,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0961329,'amu*angstrom^2'), symmetry=1, barrier=(25.5573,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3515.61,'J/mol'), sigma=(5.8495,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=549.13 K, Pc=39.86 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92363,0.0346015,-6.31345e-07,-2.51693e-08,1.16693e-11,70114.7,22.8863], Tmin=(100,'K'), Tmax=(1046.06,'K')), NASAPolynomial(coeffs=[13.7985,0.0140255,-6.73462e-06,1.41395e-09,-1.07535e-13,66271.8,-41.4319], Tmin=(1046.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.266,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(Allyl_P) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[O]O[CH][C]=C[CH][C]=O(20782)',
    structure = SMILES('[O]O[CH][C]=C[CH][C]=O'),
    E0 = (618.436,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3010,987.5,1337.5,450,1655,3000,3050,390,425,1340,1360,335,370,1685,370,1855,455,950,462.75,462.757],'cm^-1')),
        HinderedRotor(inertia=(0.176977,'amu*angstrom^2'), symmetry=1, barrier=(26.8925,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000787162,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.176985,'amu*angstrom^2'), symmetry=1, barrier=(26.8926,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.787221,'amu*angstrom^2'), symmetry=1, barrier=(119.627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01138,0.0557137,-4.15332e-05,1.00772e-08,4.58687e-13,74496.6,31.759], Tmin=(100,'K'), Tmax=(1154.34,'K')), NASAPolynomial(coeffs=[17.2337,0.0144306,-7.28887e-06,1.50463e-09,-1.11323e-13,69756.7,-53.1413], Tmin=(1154.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(618.436,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJC=O) + radical(C=CCJO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OC[C]=[C][CH][C]=O(20783)',
    structure = SMILES('[O]OC[C]=[C][CH][C]=O'),
    E0 = (738.982,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1670,1700,300,440,1855,455,950,215.364,1197.92],'cm^-1')),
        HinderedRotor(inertia=(0.214208,'amu*angstrom^2'), symmetry=1, barrier=(7.03898,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.536377,'amu*angstrom^2'), symmetry=1, barrier=(17.625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.535014,'amu*angstrom^2'), symmetry=1, barrier=(17.6249,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0173044,'amu*angstrom^2'), symmetry=1, barrier=(17.6249,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38444,0.0639881,-7.99939e-05,5.54836e-08,-1.61352e-11,88967.6,30.3086], Tmin=(100,'K'), Tmax=(822.266,'K')), NASAPolynomial(coeffs=[8.65513,0.0286193,-1.54737e-05,3.17307e-09,-2.30905e-13,87771.9,-3.3455], Tmin=(822.266,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.982,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(CCJC=O) + radical(Cds_S) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[O]O[CH]C=C[CH][C]=O(20784)',
    structure = SMILES('[O]O[CH]C=C[CH][C]=O'),
    E0 = (380.594,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.07821,0.0508041,-1.61049e-05,-1.93414e-08,1.1177e-11,45891.7,31.0398], Tmin=(100,'K'), Tmax=(1045.88,'K')), NASAPolynomial(coeffs=[17.6713,0.016135,-7.67625e-06,1.60771e-09,-1.22302e-13,40846.1,-57.285], Tmin=(1045.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(380.594,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(C=CCJO) + radical(CCJC=O) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OC[C]=[C]C[C]=O(19766)',
    structure = SMILES('[O]OC[C]=[C]C[C]=O'),
    E0 = (576.854,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1670,1700,300,440,1855,455,950,214.701,214.843],'cm^-1')),
        HinderedRotor(inertia=(0.00365588,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228202,'amu*angstrom^2'), symmetry=1, barrier=(7.47312,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228277,'amu*angstrom^2'), symmetry=1, barrier=(7.47274,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.228215,'amu*angstrom^2'), symmetry=1, barrier=(7.47286,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.951705,0.0755195,-0.000122639,1.14401e-07,-4.25905e-11,69481,31.4491], Tmin=(100,'K'), Tmax=(785.995,'K')), NASAPolynomial(coeffs=[6.24614,0.034757,-1.8476e-05,3.68384e-09,-2.60575e-13,69075.6,9.8968], Tmin=(785.995,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(576.854,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OCC=[C][CH][C]=O(20785)',
    structure = SMILES('[O]OC[CH][C]=C[C]=O'),
    E0 = (439.052,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0256159,0.0998921,-0.000176895,1.6173e-07,-5.70262e-11,52939.8,29.0228], Tmin=(100,'K'), Tmax=(830.338,'K')), NASAPolynomial(coeffs=[9.60694,0.0325804,-1.75256e-05,3.44826e-09,-2.39829e-13,52060.9,-11.3175], Tmin=(830.338,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(439.052,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJCO) + radical(Cds_S) + radical(C=CCJ=O)"""),
)

species(
    label = 'O=[C][CH]C=[C][CH]OO(20786)',
    structure = SMILES('O=[C][CH]C=[C][CH]OO'),
    E0 = (466.431,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.744471,0.0593507,-3.70217e-05,1.05151e-09,4.09338e-12,56226.3,31.9572], Tmin=(100,'K'), Tmax=(1106.96,'K')), NASAPolynomial(coeffs=[18.8488,0.0160143,-8.22307e-06,1.72979e-09,-1.30027e-13,50865.1,-63.337], Tmin=(1106.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(466.431,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(C=CCJO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OC[C]=[C][CH]C=O(20787)',
    structure = SMILES('[O]C=C[C]=[C]CO[O]'),
    E0 = (497.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00555,'amu*angstrom^2'), symmetry=1, barrier=(23.1196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01055,'amu*angstrom^2'), symmetry=1, barrier=(23.2346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01144,'amu*angstrom^2'), symmetry=1, barrier=(23.255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.299151,0.0807938,-0.000106423,7.10295e-08,-1.83846e-11,60002.5,29.7002], Tmin=(100,'K'), Tmax=(955.668,'K')), NASAPolynomial(coeffs=[15.5823,0.0168263,-6.02239e-06,9.91362e-10,-6.30231e-14,57081.3,-43.3392], Tmin=(955.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(Cds_S) + radical(C=CJC=C)"""),
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
    label = 'O=[C][CH][C]=[C]COO(20788)',
    structure = SMILES('O=[C][CH][C]=[C]COO'),
    E0 = (586.978,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1670,1700,300,440,1855,455,950,263.274,1486.37],'cm^-1')),
        HinderedRotor(inertia=(0.0138832,'amu*angstrom^2'), symmetry=1, barrier=(21.7506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.436718,'amu*angstrom^2'), symmetry=1, barrier=(21.7529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.441718,'amu*angstrom^2'), symmetry=1, barrier=(21.7563,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.439082,'amu*angstrom^2'), symmetry=1, barrier=(21.7705,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.432421,'amu*angstrom^2'), symmetry=1, barrier=(21.7524,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32304,0.0652722,-6.75317e-05,3.63618e-08,-8.2011e-12,70688.2,29.7648], Tmin=(100,'K'), Tmax=(1038.05,'K')), NASAPolynomial(coeffs=[10.5634,0.0296654,-1.60789e-05,3.3169e-09,-2.42624e-13,68769.9,-15.1594], Tmin=(1038.05,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(586.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(CCJC=O) + radical(Cds_S) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[O]O[CH][C]=C[CH]C=O(20789)',
    structure = SMILES('[O]C=CC=[C][CH]O[O]'),
    E0 = (416.076,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.502758,0.0641661,-4.07368e-05,-8.07605e-09,1.18677e-11,50180.1,29.9695], Tmin=(100,'K'), Tmax=(937.612,'K')), NASAPolynomial(coeffs=[21.2478,0.00758513,-1.28481e-06,1.82391e-10,-1.58138e-14,44886.8,-76.2591], Tmin=(937.612,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(416.076,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S)"""),
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
    label = '[O]OC[C]C1[CH]C1=O(20790)',
    structure = SMILES('[O]OC[C]C1[CH]C1=O'),
    E0 = (622.503,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.589299,0.0705747,-7.5259e-05,4.06135e-08,-8.59768e-12,74996.4,27.6465], Tmin=(100,'K'), Tmax=(1156.02,'K')), NASAPolynomial(coeffs=[15.7915,0.0179724,-7.00403e-06,1.2511e-09,-8.51219e-14,71481.6,-47.8992], Tmin=(1156.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(622.503,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(ROOJ) + radical(CCJC=O) + radical(CCJ2_triplet)"""),
)

species(
    label = 'O=[C]C1[CH][C]COO1(20791)',
    structure = SMILES('O=[C]C1[CH][C]COO1'),
    E0 = (535.501,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.359,0.0402036,2.72937e-05,-8.10306e-08,3.99272e-11,64518,27.5072], Tmin=(100,'K'), Tmax=(890.44,'K')), NASAPolynomial(coeffs=[19.8772,0.00596465,2.51547e-06,-7.45098e-10,5.37637e-14,59279.6,-70.5799], Tmin=(890.44,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.501,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxane) + radical(CCJCOOH) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][CH]C1[C]COO1(20792)',
    structure = SMILES('O=[C][CH]C1[C]COO1'),
    E0 = (540.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.228689,0.0646118,-6.31051e-05,3.14072e-08,-5.82431e-12,65177.3,31.1932], Tmin=(100,'K'), Tmax=(1577,'K')), NASAPolynomial(coeffs=[15.4663,0.0127138,-1.13945e-06,-1.15627e-10,1.74919e-14,62018.7,-44.0371], Tmin=(1577,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(540.655,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(12dioxolane) + radical(CCJCO) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = '[O]OC[C]C=C=C=O(20793)',
    structure = SMILES('[O]OC[C]=C[C]=C=O'),
    E0 = (538.446,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1700,300,440,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.318976,'amu*angstrom^2'), symmetry=1, barrier=(7.3339,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318445,'amu*angstrom^2'), symmetry=1, barrier=(7.32167,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.319282,'amu*angstrom^2'), symmetry=1, barrier=(7.34093,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.860143,0.0837151,-0.000161916,1.58684e-07,-5.74797e-11,64859.2,29.6065], Tmin=(100,'K'), Tmax=(876.596,'K')), NASAPolynomial(coeffs=[4.13396,0.0339398,-1.71307e-05,3.23712e-09,-2.17824e-13,65623.7,21.8778], Tmin=(876.596,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(538.446,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[O]OC[C]C=[C]C=O(20794)',
    structure = SMILES('[O]C=[C]C=[C]CO[O]'),
    E0 = (497.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.00555,'amu*angstrom^2'), symmetry=1, barrier=(23.1196,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01055,'amu*angstrom^2'), symmetry=1, barrier=(23.2346,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.01144,'amu*angstrom^2'), symmetry=1, barrier=(23.255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.299151,0.0807938,-0.000106423,7.10295e-08,-1.83846e-11,60002.5,29.7002], Tmin=(100,'K'), Tmax=(955.668,'K')), NASAPolynomial(coeffs=[15.5823,0.0168263,-6.02239e-06,9.91362e-10,-6.30231e-14,57081.3,-43.3392], Tmin=(955.668,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(497.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = 'O=[C][C]=C[C]COO(20795)',
    structure = SMILES('O=C=[C][CH][C]COO'),
    E0 = (584.012,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,1685,370,2120,512.5,787.5,200,800,1200,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.194867,0.0911222,-0.000145402,1.23893e-07,-4.17935e-11,70370.3,30.5326], Tmin=(100,'K'), Tmax=(806.446,'K')), NASAPolynomial(coeffs=[10.8725,0.0295875,-1.50007e-05,2.91057e-09,-2.0219e-13,68926.9,-16.9553], Tmin=(806.446,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(584.012,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(CCJC(C)=C=O) + radical(CCJ2_triplet) + radical(CCCJ=C=O)"""),
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
    label = '[O]C=[C][CH][C]CO[O](20796)',
    structure = SMILES('[O]C=[C][CH][C]CO[O]'),
    E0 = (783.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,1685,370,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.456101,0.0781121,-9.61556e-05,5.98237e-08,-1.45996e-11,94327.1,29.8136], Tmin=(100,'K'), Tmax=(1005.47,'K')), NASAPolynomial(coeffs=[15.2681,0.0191855,-8.24484e-06,1.53436e-09,-1.06303e-13,91348.6,-41.7262], Tmin=(1005.47,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(783.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(Allyl_S) + radical(CCJ2_triplet) + radical(Cds_S)"""),
)

species(
    label = '[O]OC[C]C1C=[C]O1(20797)',
    structure = SMILES('[O]OC[C]C1C=[C]O1'),
    E0 = (690.364,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.594453,0.0626592,-4.04765e-05,-4.70148e-09,9.50783e-12,83165.4,29.5027], Tmin=(100,'K'), Tmax=(958.996,'K')), NASAPolynomial(coeffs=[20.2268,0.00966419,-2.77638e-06,5.06212e-10,-3.94965e-14,78071.3,-71.3171], Tmin=(958.996,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(690.364,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(ROOJ) + radical(CCJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = 'O=C1[CH][CH][C]COO1(20798)',
    structure = SMILES('O=C1[CH][CH][C]COO1'),
    E0 = (419.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71384,-0.00902127,0.000239861,-3.51577e-07,1.47928e-10,50579,29.3788], Tmin=(100,'K'), Tmax=(906.519,'K')), NASAPolynomial(coeffs=[42.4733,-0.0339443,2.47459e-05,-4.85065e-09,3.1543e-13,36823.4,-198.373], Tmin=(906.519,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.391,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(Cycloheptane) + radical(CCJCO) + radical(CCJCC=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[C]1[CH]C=[C]OOOC1(20799)',
    structure = SMILES('[C]1[CH]C=[C]COOO1'),
    E0 = (823.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70715,0.029998,5.09686e-05,-8.96055e-08,3.60839e-11,99118.4,21.5442], Tmin=(100,'K'), Tmax=(986.616,'K')), NASAPolynomial(coeffs=[17.0907,0.0179922,-7.34794e-06,1.53832e-09,-1.21301e-13,93631.6,-64.8874], Tmin=(986.616,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(823.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(O2s-OsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclooctane) + radical(C=CCJCO) + radical(Cds_S) + radical(CH2_triplet)"""),
)

species(
    label = '[O]OC[C]C[C]=C=O(20800)',
    structure = SMILES('[O]OC[C]C[C]=C=O'),
    E0 = (588.715,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2783.33,2816.67,2850,1425,1450,1225,1275,1270,1340,700,800,300,400,1685,370,2120,512.5,787.5,180,180,180,2388.23],'cm^-1')),
        HinderedRotor(inertia=(0.627064,'amu*angstrom^2'), symmetry=1, barrier=(14.4174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.625875,'amu*angstrom^2'), symmetry=1, barrier=(14.3901,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.54926,'amu*angstrom^2'), symmetry=1, barrier=(58.6124,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.49455,'amu*angstrom^2'), symmetry=1, barrier=(34.3626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.464673,0.0846749,-0.000133503,1.133e-07,-3.78191e-11,70926.7,29.5865], Tmin=(100,'K'), Tmax=(835.808,'K')), NASAPolynomial(coeffs=[9.99752,0.0284799,-1.36767e-05,2.5875e-09,-1.76851e-13,69702.4,-12.485], Tmin=(835.808,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(588.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-CsCsHH) + group(Cs-(Cds-Cdd-O2d)CsHH) + group(Cs-CsOsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(CCJ2_triplet) + radical(CCCJ=C=O)"""),
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
    label = '[O]OC[C]=CC=C=O(19433)',
    structure = SMILES('[O]OC[C]=CC=C=O'),
    E0 = (300.604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.344476,'amu*angstrom^2'), symmetry=1, barrier=(7.92019,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344439,'amu*angstrom^2'), symmetry=1, barrier=(7.91933,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.343794,'amu*angstrom^2'), symmetry=1, barrier=(7.9045,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.90366,0.0789889,-0.000136681,1.28851e-07,-4.62937e-11,36255.4,28.9768], Tmin=(100,'K'), Tmax=(859.452,'K')), NASAPolynomial(coeffs=[5.2494,0.0345689,-1.69289e-05,3.20628e-09,-2.18019e-13,36401.9,13.8679], Tmin=(859.452,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(300.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = '[O]OCC1=CC=C1[O](20802)',
    structure = SMILES('[O]OC[C]1C=CC1=O'),
    E0 = (218.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86048,0.0524983,-5.00172e-05,3.16135e-08,-9.50738e-12,26344.4,23.2776], Tmin=(100,'K'), Tmax=(755.847,'K')), NASAPolynomial(coeffs=[4.87765,0.0365314,-1.83312e-05,3.66652e-09,-2.63959e-13,25888.3,9.56598], Tmin=(755.847,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(218.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-Cds)Cs) + group(Cd-Cd(CO)H) + ring(Cyclobutene) + radical(ROOJ) + radical(CCJ(C)CO)"""),
)

species(
    label = '[O]C=CC=C=CO[O](19443)',
    structure = SMILES('[O]C=CC=C=CO[O]'),
    E0 = (234.877,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.459565,0.0782943,-8.58339e-05,4.37348e-08,-8.2199e-12,28426.7,30.7981], Tmin=(100,'K'), Tmax=(1510.91,'K')), NASAPolynomial(coeffs=[22.7521,0.00374958,1.17197e-06,-3.91144e-10,3.03938e-14,22907.1,-85.8185], Tmin=(1510.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(234.877,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = '[O][C]=[C]C=CCO[O](20803)',
    structure = SMILES('[O]OC[CH]C=[C][C]=O'),
    E0 = (445.172,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.112712,0.0999586,-0.000171135,1.51694e-07,-5.26969e-11,53680.7,29.5672], Tmin=(100,'K'), Tmax=(800.3,'K')), NASAPolynomial(coeffs=[11.3945,0.0299298,-1.64249e-05,3.27872e-09,-2.30813e-13,52239.6,-20.8813], Tmin=(800.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(445.172,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJCO) + radical(C=CJC=O) + radical(C=CCJ=O)"""),
)

species(
    label = '[O]OC[C]=C[C]=[C]O(20804)',
    structure = SMILES('[O]OC[C]=C[C]=[C]O'),
    E0 = (596.058,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.766844,'amu*angstrom^2'), symmetry=1, barrier=(17.6313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766866,'amu*angstrom^2'), symmetry=1, barrier=(17.6318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766372,'amu*angstrom^2'), symmetry=1, barrier=(17.6204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766341,'amu*angstrom^2'), symmetry=1, barrier=(17.6197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0615833,0.0899621,-0.000138325,1.05795e-07,-3.08432e-11,71828,32.2702], Tmin=(100,'K'), Tmax=(946.081,'K')), NASAPolynomial(coeffs=[14.8489,0.0170011,-6.09249e-06,9.51897e-10,-5.61169e-14,69497.2,-35.7811], Tmin=(946.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[O]OC[C]=[C]C=[C]O(20805)',
    structure = SMILES('[O]OC[C]=[C]C=[C]O'),
    E0 = (596.058,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.766844,'amu*angstrom^2'), symmetry=1, barrier=(17.6313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766866,'amu*angstrom^2'), symmetry=1, barrier=(17.6318,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766372,'amu*angstrom^2'), symmetry=1, barrier=(17.6204,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.766341,'amu*angstrom^2'), symmetry=1, barrier=(17.6197,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0615833,0.0899621,-0.000138325,1.05795e-07,-3.08432e-11,71828,32.2702], Tmin=(100,'K'), Tmax=(946.081,'K')), NASAPolynomial(coeffs=[14.8489,0.0170011,-6.09249e-06,9.51897e-10,-5.61169e-14,69497.2,-35.7811], Tmin=(946.081,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(596.058,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJC=C) + radical(C=CJO)"""),
)

species(
    label = '[O]O[CH][C]=CC=[C]O(20806)',
    structure = SMILES('[O]O[CH][C]=CC=[C]O'),
    E0 = (514.358,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,2995,3025,975,1000,1300,1375,400,500,1630,1680,3025,407.5,1350,352.5,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.913579,'amu*angstrom^2'), symmetry=1, barrier=(21.005,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.913507,'amu*angstrom^2'), symmetry=1, barrier=(21.0033,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.913443,'amu*angstrom^2'), symmetry=1, barrier=(21.0018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.913648,'amu*angstrom^2'), symmetry=1, barrier=(21.0066,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.142986,0.0781335,-8.93575e-05,4.82921e-08,-9.77765e-12,62023.3,34.0048], Tmin=(100,'K'), Tmax=(1345.55,'K')), NASAPolynomial(coeffs=[21.4501,0.00614725,-4.18358e-07,-7.91474e-11,9.55152e-15,56918,-73.9561], Tmin=(1345.55,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(514.358,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CCJO) + radical(Cds_S) + radical(C=CJO)"""),
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
    E0 = (501.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (501.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (746.341,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (970.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1047.54,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (947.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (504.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (509.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (632.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (523.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (579.388,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (526.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (621.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (636.588,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (679.022,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (501.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (573.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (830.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (950.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (640.496,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (759.42,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (705.278,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (665.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (647.82,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (587.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (620.018,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (565.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (1143.21,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (682.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (555.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (591.217,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (549.904,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (766.92,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (660.556,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (698.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (1068.12,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (806.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (690.519,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (563.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (823.279,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (766.401,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (710.996,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (501.141,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (509.425,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (526.114,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (754.986,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (651.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (789.136,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (744.5,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (614.151,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['HCCO(2227)', 'C#CCO[O](16808)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['O2(S)(666)', 'C=[C]C=C[C]=O(16163)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[O]C[C]=C[CH][C]=O(20774)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH2]O[O](46)', '[C]=C[CH][C]=O(15285)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/O;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]=[C][O](6861)', '[CH]=[C]CO[O](16805)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[O]OC[C]=C[C][C]=O(20775)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]OCC1=CC1[C]=O(20756)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriND_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['O=[C][CH]C=C1COO1(20776)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['O(4)', 'O=[C][CH]C=C1CO1(20777)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO_S;Y_rad_intra;OOJ] for rate rule [R2OO_S;Cd_rad_out;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]OCC=C=C[C]=O(20778)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['O=[C][CH]C=C=COO(20779)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.00798e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]OC=C=CC[C]=O(19763)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['H(3)', '[O]OC=C=C[CH][C]=O(20780)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(169.619,'m^3/(mol*s)'), n=1.605, Ea=(12.4249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Ca;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', '[O]OC[C]=C=C[C]=O(20781)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C][O](6861)', 'C#CCO[O](16808)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(96.4727,'m^3/(mol*s)'), n=1.65897, Ea=(11.5315,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-H_Ct-Cs;YJ] for rate rule [Ct-H_Ct-Cs;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['O2(2)', 'C=[C]C=C[C]=O(16163)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.973e+06,'cm^3/(mol*s)'), n=2.037, Ea=(179.712,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds-HH_Ca;OJ] for rate rule [Cds-HH_Ca;O2b]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from -6.0 to 179.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['O2(2)', 'C=[C][CH][CH][C]=O(16165)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[O]O[CH][C]=C[CH][C]=O(20782)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[O]OC[C]=[C][CH][C]=O(20783)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]O[CH]C=C[CH][C]=O(20784)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.89098e+10,'s^-1'), n=0.9884, Ea=(139.355,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;Cs_H_out_1H] for rate rule [R2H_S;Cd_rad_out_Cd;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O]OC[C]=[C]C[C]=O(19766)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(7.28974e+09,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]OCC=[C][CH][C]=O(20785)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['O=[C][CH]C=[C][CH]OO(20786)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(5.66219e+08,'s^-1'), n=1.32333, Ea=(164.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O]OC[C]=[C][CH]C=O(20787)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(781966,'s^-1'), n=1.9774, Ea=(150.043,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out;XH_out] for rate rule [R3HJ;Cd_rad_out;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]O[CH][C]=CC[C]=O(19765)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(323087,'s^-1'), n=1.92599, Ea=(85.9746,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R4HJ_2;Y_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O=[C][CH][C]=[C]COO(20788)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;XH_out] for rate rule [R5HJ_1;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]O[CH][C]=C[CH]C=O(20789)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(26847.4,'s^-1'), n=2.09453, Ea=(64.7206,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Y_rad_out;Cs_H_out_H/NonDeO] for rate rule [R5Hall;CO_rad_out;Cs_H_out_H/NonDeO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[C]=O(2355)', '[CH]=C[C]CO[O](20284)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]OC[C]C1[CH]C1=O(20790)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['O=[C]C1[CH][C]COO1(20791)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(4.18543e+10,'s^-1'), n=0.209443, Ea=(54.3406,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra_pri;radadd_intra] for rate rule [R6_linear;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['O=[C][CH]C1[C]COO1(20792)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(4.34604e+09,'s^-1'), n=0.547362, Ea=(90.0766,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra;radadd_intra_O] + [R6;doublebond_intra;radadd_intra] for rate rule [R6;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction32',
    reactants = ['CO(2039)', '[CH]=C[C]CO[O](20284)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.51e+11,'cm^3/(mol*s)','*|/',5), n=0, Ea=(20.125,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(2500,'K'), comment="""From training reaction 7 used for COm;Cd_pri_rad
Exact match found for rate rule [COm;Cd_pri_rad]
Euclidian distance = 0
family: R_Addition_COm"""),
)

reaction(
    label = 'reaction33',
    reactants = ['H(3)', '[O]OC[C]C=C=C=O(20793)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(6.23e+08,'cm^3/(mol*s)'), n=1.429, Ea=(16.6816,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Cdd;HJ] for rate rule [Ca_Ck;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[O]OC[C]C=[C]C=O(20794)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['O=[C][C]=C[C]COO(20795)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(15567.1,'s^-1'), n=2.15754, Ea=(114.223,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;Y_rad_out;O_H_out] + [RnH;Cd_rad_out;XH_out] + [R6Hall;Y_rad_out;XH_out] for rate rule [R6Hall;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['HCCO(2227)', '[CH][C]CO[O](20242)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[O]C=[C][CH][C]CO[O](20796)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]OC[C]C1C=[C]O1(20797)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(8.73504e+09,'s^-1'), n=0.685238, Ea=(189.379,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;multiplebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['O=C1[CH][CH][C]COO1(20798)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(8.62196e+06,'s^-1'), n=0.867572, Ea=(62.1704,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R7_linear;doublebond_intra;radadd_intra] for rate rule [R7_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[C]1[CH]C=[C]OOOC1(20799)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(3.50627e+10,'s^-1'), n=0.368321, Ea=(322.138,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6plus;multiplebond_intra;radadd_intra] for rate rule [R8;multiplebond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 318.1 to 322.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction41',
    reactants = ['[O]OC[C]C[C]=C=O(20800)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(7.80481e+09,'s^-1'), n=1.09833, Ea=(177.686,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_double;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[O]O[CH][C]CC=C=O(20801)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.60984e+10,'s^-1'), n=0.8, Ea=(135.98,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;C_rad_out_1H;Cs_H_out_H/Cd] for rate rule [R3Hall;C_rad_out_H/NonDeO;Cs_H_out_H/Cd]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]OC[C]=CC=C=O(19433)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]OCC1=CC=C1[O](20802)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_DSD;CdsingleND_rad_out;Ypri_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O]C=CC=C=CO[O](19443)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(4.25221e+09,'s^-1'), n=0.137, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad;XH_Rrad] for rate rule [R5radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction46',
    reactants = ['HCCO(2227)', '[CH]=[C]CO[O](16805)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[O][C]=CC=[C]CO[O](19437)'],
    products = ['[O][C]=[C]C=CCO[O](20803)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[O]OC[C]=C[C]=[C]O(20804)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.28371e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[O]OC[C]=[C]C=[C]O(20805)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1194.19,'s^-1'), n=2.56519, Ea=(148.442,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4Hall;Y_rad_out;O_H_out] + [R4Hall;Cd_rad_out;XH_out] for rate rule [R4HJ_2;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[O]O[CH][C]=CC=[C]O(20806)'],
    products = ['[O][C]=CC=[C]CO[O](19437)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(0.720455,'s^-1'), n=3.59722, Ea=(99.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [RnH;C_rad_out_1H;O_H_out] + [RnH;C_rad_out_H/NonDeO;XH_out] + [R6Hall;C_rad_out_1H;XH_out] for rate rule [R6Hall;C_rad_out_H/NonDeO;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4458',
    isomers = [
        '[O][C]=CC=[C]CO[O](19437)',
    ],
    reactants = [
        ('HCCO(2227)', 'C#CCO[O](16808)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4458',
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


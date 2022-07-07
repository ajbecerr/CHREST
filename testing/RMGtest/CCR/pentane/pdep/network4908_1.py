species(
    label = 'C=[C]C(C=[C][O])O[O](21197)',
    structure = SMILES('C=[C]C([CH][C]=O)O[O]'),
    E0 = (525.306,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,245.545,245.545],'cm^-1')),
        HinderedRotor(inertia=(0.00279621,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135223,'amu*angstrom^2'), symmetry=1, barrier=(5.78669,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135262,'amu*angstrom^2'), symmetry=1, barrier=(5.78673,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.135227,'amu*angstrom^2'), symmetry=1, barrier=(5.78673,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.772036,0.0763676,-0.00011509,9.46046e-08,-3.09951e-11,63290.8,34.5187], Tmin=(100,'K'), Tmax=(811.033,'K')), NASAPolynomial(coeffs=[9.97517,0.0258756,-1.22689e-05,2.32909e-09,-1.60239e-13,61965.8,-6.91907], Tmin=(811.033,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(525.306,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJCO) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    label = 'C=[C]C([O])[CH][C]=O(23311)',
    structure = SMILES('C=[C]C([O])[CH][C]=O'),
    E0 = (532.157,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,180,336.708,4000],'cm^-1')),
        HinderedRotor(inertia=(0.663985,'amu*angstrom^2'), symmetry=1, barrier=(15.2663,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.022572,'amu*angstrom^2'), symmetry=1, barrier=(15.2701,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0226698,'amu*angstrom^2'), symmetry=1, barrier=(15.2859,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52112,0.0568633,-6.51062e-05,4.01565e-08,-1.00597e-11,64091,30.1248], Tmin=(100,'K'), Tmax=(964.031,'K')), NASAPolynomial(coeffs=[10.0957,0.0212854,-9.74812e-06,1.87418e-09,-1.32059e-13,62437.8,-10.9284], Tmin=(964.031,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(532.157,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCJCO) + radical(Cds_S) + radical(CCCJ=O)"""),
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
    label = '[O]O[CH][CH][C]=O(23312)',
    structure = SMILES('[O]O[CH][CH][C]=O'),
    E0 = (419.093,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3000,3050,390,425,1340,1360,335,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.11437,'amu*angstrom^2'), symmetry=1, barrier=(2.62958,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1146,'amu*angstrom^2'), symmetry=1, barrier=(2.63488,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114572,'amu*angstrom^2'), symmetry=1, barrier=(2.63424,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (86.0462,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48119,0.0697409,-0.000150812,1.49657e-07,-5.33379e-11,50482.4,24.544], Tmin=(100,'K'), Tmax=(898.831,'K')), NASAPolynomial(coeffs=[4.56322,0.0208135,-1.03979e-05,1.9271e-09,-1.25857e-13,51350.8,17.916], Tmin=(898.831,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(419.093,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(ROOJ) + radical(CCsJOOH) + radical(CCJCO) + radical(CCCJ=O)"""),
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
    label = 'C=[C]C([C][C]=O)O[O](23313)',
    structure = SMILES('C=[C]C([C][C]=O)O[O]'),
    E0 = (773.619,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,267.295,267.305,267.327],'cm^-1')),
        HinderedRotor(inertia=(0.17258,'amu*angstrom^2'), symmetry=1, barrier=(8.75017,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.172608,'amu*angstrom^2'), symmetry=1, barrier=(8.75021,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.1726,'amu*angstrom^2'), symmetry=1, barrier=(8.75007,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.400375,'amu*angstrom^2'), symmetry=1, barrier=(20.2967,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.495416,0.0850056,-0.000144594,1.25597e-07,-4.19965e-11,93163.5,32.919], Tmin=(100,'K'), Tmax=(856.424,'K')), NASAPolynomial(coeffs=[10.4551,0.0238322,-1.17811e-05,2.22988e-09,-1.51286e-13,91995,-10.4491], Tmin=(856.424,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(773.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C1C([C]=O)C1O[O](23314)',
    structure = SMILES('C=C1C([C]=O)C1O[O]'),
    E0 = (282.861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0899,0.0614113,-5.98586e-05,2.98166e-08,-5.93146e-12,34127.3,25.9587], Tmin=(100,'K'), Tmax=(1211.34,'K')), NASAPolynomial(coeffs=[13.6158,0.0200492,-8.64004e-06,1.62827e-09,-1.13875e-13,31092.6,-36.8732], Tmin=(1211.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(ROOJ) + radical(CC(C)CJ=O)"""),
)

species(
    label = 'C=[C]C1OOC1[C]=O(23315)',
    structure = SMILES('C=[C]C1OOC1[C]=O'),
    E0 = (332.7,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18815,0.0530781,-3.57676e-05,6.75183e-09,1.33231e-12,40123.2,27.6903], Tmin=(100,'K'), Tmax=(1104.41,'K')), NASAPolynomial(coeffs=[14.6733,0.018388,-7.87137e-06,1.51435e-09,-1.08363e-13,36281.6,-42.6139], Tmin=(1104.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(332.7,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C1OOC1[CH][C]=O(23316)',
    structure = SMILES('C=C1OOC1[CH][C]=O'),
    E0 = (301.832,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45015,0.04359,-4.4124e-06,-2.70543e-08,1.32703e-11,36404.7,29.9605], Tmin=(100,'K'), Tmax=(1027.8,'K')), NASAPolynomial(coeffs=[15.1834,0.0180406,-7.8396e-06,1.5779e-09,-1.17896e-13,32108.2,-43.8392], Tmin=(1027.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.832,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C1OC1[C]=O(19001)',
    structure = SMILES('C=[C]C1OC1[C]=O'),
    E0 = (281.404,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72971,0.0366141,1.45106e-05,-5.98519e-08,3.12096e-11,33939.5,25.5262], Tmin=(100,'K'), Tmax=(874.285,'K')), NASAPolynomial(coeffs=[16.8232,0.00564275,2.30779e-06,-7.23463e-10,5.50373e-14,29844.8,-53.5873], Tmin=(874.285,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.404,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C1OC1[CH][C]=O(23317)',
    structure = SMILES('C=C1OC1[CH][C]=O'),
    E0 = (235.252,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34919,0.0426914,1.66653e-06,-4.69805e-08,2.50234e-11,28404.1,23.6172], Tmin=(100,'K'), Tmax=(929.974,'K')), NASAPolynomial(coeffs=[19.881,0.00302688,1.05377e-06,-2.39114e-10,1.09935e-14,23225.6,-73.7529], Tmin=(929.974,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C=C(C[C]=O)O[O](21754)',
    structure = SMILES('C=C=C(C[C]=O)O[O]'),
    E0 = (269.342,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.329395,'amu*angstrom^2'), symmetry=1, barrier=(9.39434,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.333425,'amu*angstrom^2'), symmetry=1, barrier=(9.43758,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.396256,'amu*angstrom^2'), symmetry=1, barrier=(9.41107,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41941,0.063032,-7.04381e-05,4.32701e-08,-1.12314e-11,32481.9,27.4878], Tmin=(100,'K'), Tmax=(910.693,'K')), NASAPolynomial(coeffs=[8.98279,0.0298102,-1.57161e-05,3.20937e-09,-2.33549e-13,31104.4,-8.29327], Tmin=(910.693,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(269.342,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(ROOJ) + radical(CCCJ=O)"""),
)

species(
    label = 'C=CC(=C[C]=O)O[O](23318)',
    structure = SMILES('[CH2]C=C(C=C=O)O[O]'),
    E0 = (165.199,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04949,0.0693992,-9.66332e-05,7.75823e-08,-2.5192e-11,19970.8,27.8351], Tmin=(100,'K'), Tmax=(824.564,'K')), NASAPolynomial(coeffs=[8.51554,0.0282016,-1.26305e-05,2.3417e-09,-1.59249e-13,18908.8,-5.71759], Tmin=(824.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(C=CC=CCJ)"""),
)

species(
    label = 'C=[C]C(=C[C]=O)OO(23319)',
    structure = SMILES('C=[C]C(=C[C]=O)OO'),
    E0 = (240.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0950707,0.100387,-0.00017583,1.56234e-07,-5.33364e-11,29081.3,28.1768], Tmin=(100,'K'), Tmax=(846.709,'K')), NASAPolynomial(coeffs=[11.0673,0.0285803,-1.48294e-05,2.86294e-09,-1.96442e-13,27874.8,-19.7804], Tmin=(846.709,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(240.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(C=CJC=C) + radical(C=CCJ=O)"""),
)

species(
    label = 'C=[C]C(=C[C]=O)O[O](23320)',
    structure = SMILES('C=[C]C(=C[C]=O)O[O]'),
    E0 = (392.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,180],'cm^-1')),
        HinderedRotor(inertia=(0.811893,'amu*angstrom^2'), symmetry=1, barrier=(18.667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.812655,'amu*angstrom^2'), symmetry=1, barrier=(18.6845,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.811544,'amu*angstrom^2'), symmetry=1, barrier=(18.659,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.158002,0.0969185,-0.00018091,1.65859e-07,-5.71052e-11,47352.2,28.0274], Tmin=(100,'K'), Tmax=(874.189,'K')), NASAPolynomial(coeffs=[9.46695,0.0269407,-1.38506e-05,2.62515e-09,-1.76561e-13,46771,-9.64671], Tmin=(874.189,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cd-Cd(CO)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CJC=C) + radical(C=CCJ=O)"""),
)

species(
    label = 'C#CC([CH][C]=O)O[O](23321)',
    structure = SMILES('C#CC([CH][C]=O)O[O]'),
    E0 = (453.414,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2175,525,1855,455,950,750,770,3400,2100,379.344],'cm^-1')),
        HinderedRotor(inertia=(0.0879291,'amu*angstrom^2'), symmetry=1, barrier=(8.98283,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0878806,'amu*angstrom^2'), symmetry=1, barrier=(8.9846,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0880437,'amu*angstrom^2'), symmetry=1, barrier=(8.98584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.708378,'amu*angstrom^2'), symmetry=1, barrier=(72.3365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700861,0.0755778,-0.000112436,8.59642e-08,-2.55298e-11,54649.1,31.9167], Tmin=(100,'K'), Tmax=(892.852,'K')), NASAPolynomial(coeffs=[12.4728,0.0179106,-7.27462e-06,1.26041e-09,-8.14243e-14,52743.4,-22.4418], Tmin=(892.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(453.414,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(ROOJ) + radical(CCJCO) + radical(CCCJ=O)"""),
)

species(
    label = '[O]OC=C[C]=O(23322)',
    structure = SMILES('[O]O[CH]C=C=O'),
    E0 = (103.279,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,2120,512.5,787.5,180],'cm^-1')),
        HinderedRotor(inertia=(0.860098,'amu*angstrom^2'), symmetry=1, barrier=(19.7753,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.12024,'amu*angstrom^2'), symmetry=1, barrier=(48.7484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (86.0462,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.85314,0.0471905,-5.8373e-05,3.67568e-08,-9.07544e-12,12498.9,19.7742], Tmin=(100,'K'), Tmax=(994.937,'K')), NASAPolynomial(coeffs=[10.6865,0.0116763,-4.82928e-06,8.78436e-10,-5.99804e-14,10741.2,-22.7966], Tmin=(994.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.279,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)OsHH) + group(Cds-(Cdd-O2d)CsH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(C=CCJO)"""),
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
    label = 'C=[C][C]([CH][C]=O)O[O](23323)',
    structure = SMILES('[CH2][C]=C(C=[C][O])O[O]'),
    E0 = (648.814,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,310.596,310.597],'cm^-1')),
        HinderedRotor(inertia=(0.114635,'amu*angstrom^2'), symmetry=1, barrier=(7.84763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114635,'amu*angstrom^2'), symmetry=1, barrier=(7.84763,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.05987,'amu*angstrom^2'), symmetry=1, barrier=(72.5562,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.951404,0.0710821,-9.55766e-05,6.3835e-08,-1.51188e-11,78140.5,31.2373], Tmin=(100,'K'), Tmax=(707.389,'K')), NASAPolynomial(coeffs=[11.7352,0.0200324,-8.37894e-06,1.49722e-09,-9.96495e-14,76366.4,-18.8108], Tmin=(707.389,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.814,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=COJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]C([CH][C]=O)O[O](23324)',
    structure = SMILES('[CH]=[C]C([CH][C]=O)O[O]'),
    E0 = (772.402,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,1685,370,1855,455,950,3120,650,792.5,1650,189.396],'cm^-1')),
        HinderedRotor(inertia=(0.00472425,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263573,'amu*angstrom^2'), symmetry=1, barrier=(6.69856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.264358,'amu*angstrom^2'), symmetry=1, barrier=(6.69799,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.263916,'amu*angstrom^2'), symmetry=1, barrier=(6.69632,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.659557,0.0803622,-0.000133448,1.14207e-07,-3.7832e-11,93012.2,35.3242], Tmin=(100,'K'), Tmax=(852.36,'K')), NASAPolynomial(coeffs=[10.4106,0.0227019,-1.10345e-05,2.08231e-09,-1.41346e-13,91782.2,-7.62544], Tmin=(852.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(772.402,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJCO) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][C](C[C]=O)O[O](21756)',
    structure = SMILES('[CH2][C]=C(C[C]=O)O[O]'),
    E0 = (482.078,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,1855,455,950,369.424],'cm^-1')),
        HinderedRotor(inertia=(0.00123537,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0951923,'amu*angstrom^2'), symmetry=1, barrier=(9.21913,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0951965,'amu*angstrom^2'), symmetry=1, barrier=(9.21896,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0952048,'amu*angstrom^2'), symmetry=1, barrier=(9.2192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.27814,0.0662945,-8.03886e-05,5.43234e-08,-1.54122e-11,58073.1,29.6959], Tmin=(100,'K'), Tmax=(841.435,'K')), NASAPolynomial(coeffs=[8.88667,0.0301258,-1.5913e-05,3.24045e-09,-2.35103e-13,56792.6,-5.69746], Tmin=(841.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(482.078,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-OsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(ROOJ) + radical(Allyl_P) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C[C]([CH][C]=O)O[O](23325)',
    structure = SMILES('[CH2]C=C(C=[C][O])O[O]'),
    E0 = (410.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.857711,0.0682615,-7.85994e-05,4.76112e-08,-1.14201e-11,49542.5,31.0828], Tmin=(100,'K'), Tmax=(1020.15,'K')), NASAPolynomial(coeffs=[13.1249,0.0201616,-7.8744e-06,1.39226e-09,-9.35294e-14,47039.6,-28.3441], Tmin=(1020.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(410.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH]=CC([CH][C]=O)O[O](23326)',
    structure = SMILES('[CH]=CC([CH][C]=O)O[O]'),
    E0 = (534.56,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,1855,455,950,267.139],'cm^-1')),
        HinderedRotor(inertia=(0.00236205,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162489,'amu*angstrom^2'), symmetry=1, barrier=(8.22884,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162489,'amu*angstrom^2'), symmetry=1, barrier=(8.22883,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.162495,'amu*angstrom^2'), symmetry=1, barrier=(8.22882,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.825821,0.0740512,-0.000101981,7.52e-08,-2.21937e-11,64403.2,34.2629], Tmin=(100,'K'), Tmax=(828.615,'K')), NASAPolynomial(coeffs=[11.2185,0.0238822,-1.11627e-05,2.13156e-09,-1.48319e-13,62680.9,-13.9219], Tmin=(828.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(534.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJCO) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C][C]([CH][C]=O)OO(23327)',
    structure = SMILES('[CH2][C]=C(C=[C][O])OO'),
    E0 = (496.809,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.647612,0.0753604,-9.45508e-05,6.18039e-08,-1.59999e-11,59871.7,31.5565], Tmin=(100,'K'), Tmax=(945.75,'K')), NASAPolynomial(coeffs=[13.3531,0.0216231,-9.32116e-06,1.72483e-09,-1.18587e-13,57468.5,-29.0315], Tmin=(945.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(496.809,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C][C]([CH]C=O)O[O](23328)',
    structure = SMILES('[CH2][C]=C(C=C[O])O[O]'),
    E0 = (409.07,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493168,0.0725881,-8.1934e-05,4.54898e-08,-9.48178e-12,49330,29.7328], Tmin=(100,'K'), Tmax=(976.907,'K')), NASAPolynomial(coeffs=[16.3798,0.0150112,-4.99952e-06,8.16845e-10,-5.30716e-14,45869.6,-48.3644], Tmin=(976.907,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.07,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=[C]C(C[C]=O)O[O](21757)',
    structure = SMILES('[CH]=[C]C(C[C]=O)O[O]'),
    E0 = (572.5,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2750,2850,1437.5,1250,1305,750,350,1685,370,1855,455,950,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.369884,'amu*angstrom^2'), symmetry=1, barrier=(8.50435,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369894,'amu*angstrom^2'), symmetry=1, barrier=(8.50458,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36986,'amu*angstrom^2'), symmetry=1, barrier=(8.5038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.369838,'amu*angstrom^2'), symmetry=1, barrier=(8.5033,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.479075,0.0854145,-0.000141543,1.23089e-07,-4.14714e-11,68975,33.4754], Tmin=(100,'K'), Tmax=(851.624,'K')), NASAPolynomial(coeffs=[9.79544,0.0270854,-1.31411e-05,2.48294e-09,-1.68828e-13,67916.6,-6.87214], Tmin=(851.624,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(572.5,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([CH][C]=O)OO(23329)',
    structure = SMILES('[CH]=[C]C([CH][C]=O)OO'),
    E0 = (620.397,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,3025,407.5,1350,352.5,1685,370,1855,455,950,3120,650,792.5,1650,345.9],'cm^-1')),
        HinderedRotor(inertia=(0.00140905,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.173864,'amu*angstrom^2'), symmetry=1, barrier=(14.7616,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101582,'amu*angstrom^2'), symmetry=1, barrier=(8.62445,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.101581,'amu*angstrom^2'), symmetry=1, barrier=(8.62436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.36624,'amu*angstrom^2'), symmetry=1, barrier=(31.0931,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.493986,0.082697,-0.000123897,9.80069e-08,-3.08871e-11,74737.6,35.1662], Tmin=(100,'K'), Tmax=(778.228,'K')), NASAPolynomial(coeffs=[11.7578,0.0247967,-1.22865e-05,2.38652e-09,-1.66857e-13,72984.6,-16.35], Tmin=(778.228,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(620.397,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCJCO) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]C([CH]C=O)O[O](23330)',
    structure = SMILES('[CH]=[C]C(C=C[O])O[O]'),
    E0 = (552.034,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.822402,'amu*angstrom^2'), symmetry=1, barrier=(18.9086,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.822112,'amu*angstrom^2'), symmetry=1, barrier=(18.902,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.821941,'amu*angstrom^2'), symmetry=1, barrier=(18.898,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.487112,0.0741741,-8.48962e-05,4.7892e-08,-1.05152e-11,66523.6,32.1785], Tmin=(100,'K'), Tmax=(1118.35,'K')), NASAPolynomial(coeffs=[16.8001,0.0158279,-6.63931e-06,1.2422e-09,-8.69959e-14,62874.8,-48.3472], Tmin=(1118.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(ROOJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]C(C=C=O)O[O](21189)',
    structure = SMILES('C=[C]C(C=C=O)O[O]'),
    E0 = (281.688,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.592712,'amu*angstrom^2'), symmetry=1, barrier=(13.6276,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.591616,'amu*angstrom^2'), symmetry=1, barrier=(13.6024,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.59336,'amu*angstrom^2'), symmetry=1, barrier=(13.6425,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.332455,0.0926989,-0.000166777,1.53007e-07,-5.28571e-11,33999.8,28.3448], Tmin=(100,'K'), Tmax=(880.403,'K')), NASAPolynomial(coeffs=[7.91628,0.030913,-1.49451e-05,2.77643e-09,-1.85134e-13,33723.6,-1.26146], Tmin=(880.403,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(281.688,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(Cds_S)"""),
)

species(
    label = 'C=C([CH][C]=O)[CH]O[O](20759)',
    structure = SMILES('[CH2]C([CH]O[O])=C[C]=O'),
    E0 = (343.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,350,440,435,1725,3025,407.5,1350,352.5,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1855,455,950,239.59],'cm^-1')),
        HinderedRotor(inertia=(0.696705,'amu*angstrom^2'), symmetry=1, barrier=(28.3754,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.440871,'amu*angstrom^2'), symmetry=1, barrier=(17.9492,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.696389,'amu*angstrom^2'), symmetry=1, barrier=(28.3755,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.440744,'amu*angstrom^2'), symmetry=1, barrier=(17.9492,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.72205,0.0767824,-9.37916e-05,5.90503e-08,-1.49814e-11,41412.1,28.4603], Tmin=(100,'K'), Tmax=(953.361,'K')), NASAPolynomial(coeffs=[12.9404,0.0255185,-1.31346e-05,2.6489e-09,-1.91433e-13,39082.3,-29.9028], Tmin=(953.361,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(343.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(ROOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(C=CCJ=O)"""),
)

species(
    label = 'C=C1C([O])=CC1O[O](23331)',
    structure = SMILES('C=C1C([O])=CC1O[O]'),
    E0 = (189.313,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.806228,0.0591382,-3.65995e-05,-5.57542e-09,9.34396e-12,22894.1,23.9317], Tmin=(100,'K'), Tmax=(950.898,'K')), NASAPolynomial(coeffs=[18.5675,0.0111711,-3.12538e-06,5.36505e-10,-3.98768e-14,18307.1,-67.2204], Tmin=(950.898,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(189.313,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C1C=C([O])OO1(23332)',
    structure = SMILES('C=[C]C1[CH]C(=O)OO1'),
    E0 = (120.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01485,0.017764,9.01467e-05,-1.35912e-07,5.43426e-11,14547.3,27.9093], Tmin=(100,'K'), Tmax=(956.416,'K')), NASAPolynomial(coeffs=[19.0554,0.0119128,-3.27365e-06,7.21064e-10,-6.58081e-14,8295.78,-69.1836], Tmin=(956.416,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(120.163,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-O2d)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCJCO) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1C=C([O])O1(23333)',
    structure = SMILES('C=[C]C1[CH]C(=O)O1'),
    E0 = (208.522,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.49791,0.0161629,6.31689e-05,-9.85519e-08,4.08877e-11,25149.1,24.701], Tmin=(100,'K'), Tmax=(918.8,'K')), NASAPolynomial(coeffs=[13.3444,0.0124757,-1.88191e-06,2.15428e-10,-1.75546e-14,21318.4,-36.7079], Tmin=(918.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(208.522,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(Beta-Propiolactone) + radical(CCJCO) + radical(Cds_S)"""),
)

species(
    label = 'HO2(9)',
    structure = SMILES('[O]O'),
    E0 = (2.67648,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1112.81,1388.53,3298.45],'cm^-1')),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (33.0067,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(892.977,'J/mol'), sigma=(3.458,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=1.0, comment="""GRI-Mech"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.02956,-0.00263985,1.5223e-05,-1.71671e-08,6.26737e-12,322.677,4.84428], Tmin=(100,'K'), Tmax=(923.914,'K')), NASAPolynomial(coeffs=[4.15133,0.00191146,-4.11273e-07,6.34955e-11,-4.86383e-15,83.4204,3.0934], Tmin=(923.914,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(2.67648,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(58.2013,'J/(mol*K)'), label="""HO2""", comment="""Thermo library: BurkeH2O2"""),
)

species(
    label = '[CH2][C]=CC#C[O](18002)',
    structure = SMILES('C=[C][CH][C]=C=O'),
    E0 = (543.519,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1670,1700,300,440,2120,512.5,787.5,194.187,194.192],'cm^-1')),
        HinderedRotor(inertia=(1.57828,'amu*angstrom^2'), symmetry=1, barrier=(42.235,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.5783,'amu*angstrom^2'), symmetry=1, barrier=(42.235,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (79.0767,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87342,0.0453999,-4.56675e-05,2.39369e-08,-5.02902e-12,65448,18.4218], Tmin=(100,'K'), Tmax=(1149.41,'K')), NASAPolynomial(coeffs=[10.5067,0.0153555,-6.45891e-06,1.19558e-09,-8.27037e-14,63463.4,-24.4311], Tmin=(1149.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(543.519,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=CCJC=C=O) + radical(Cds_S) + radical(CCCJ=C=O)"""),
)

species(
    label = 'C=C=C(C=C[O])O[O](23334)',
    structure = SMILES('C=C=C(C=C[O])O[O]'),
    E0 = (229.777,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.24253,0.0770346,-9.11549e-05,5.27273e-08,-1.17144e-11,27775.9,27.6545], Tmin=(100,'K'), Tmax=(1115.12,'K')), NASAPolynomial(coeffs=[18.2888,0.012302,-4.08069e-06,6.70969e-10,-4.38859e-14,23751.1,-61.3748], Tmin=(1115.12,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(229.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(ROOJ)"""),
)

species(
    label = 'C=[C]C(C#C[O])O[O](23335)',
    structure = SMILES('C=[C]C([C]=C=O)O[O]'),
    E0 = (483.969,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1670,1700,300,440,2120,512.5,787.5,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.64854,'amu*angstrom^2'), symmetry=1, barrier=(14.9112,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.652374,'amu*angstrom^2'), symmetry=1, barrier=(14.9994,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50952,'amu*angstrom^2'), symmetry=1, barrier=(34.7068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (111.076,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.305263,0.0938891,-0.000174547,1.60526e-07,-5.51511e-11,58329,28.1516], Tmin=(100,'K'), Tmax=(887.249,'K')), NASAPolynomial(coeffs=[8.43055,0.0280126,-1.37325e-05,2.54197e-09,-1.68204e-13,58038.3,-3.58931], Tmin=(887.249,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(483.969,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + missing(O2d-Cdd) + group(Cs-(Cds-Cdd-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-(Cdd-O2d)CsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(ROOJ) + radical(Cds_S) + radical(CCCJ=C=O)"""),
)

species(
    label = 'C=[C]C([C]=C[O])O[O](23336)',
    structure = SMILES('C=[C]C([C]=C[O])O[O]'),
    E0 = (542.78,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1670,1700,300,440,182.022,182.025,182.033],'cm^-1')),
        HinderedRotor(inertia=(0.712598,'amu*angstrom^2'), symmetry=1, barrier=(16.7599,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712897,'amu*angstrom^2'), symmetry=1, barrier=(16.76,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.712726,'amu*angstrom^2'), symmetry=1, barrier=(16.7599,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.655211,0.0737057,-8.74489e-05,5.2362e-08,-1.2344e-11,65401.8,31.6493], Tmin=(100,'K'), Tmax=(1037.23,'K')), NASAPolynomial(coeffs=[14.8109,0.0191154,-8.50289e-06,1.62059e-09,-1.14041e-13,62465.2,-37.1614], Tmin=(1037.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(542.78,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = 'C=CC([C]=[C][O])O[O](23337)',
    structure = SMILES('C=CC([C][C]=O)O[O]'),
    E0 = (535.777,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1855,455,950,287.711,287.972,288.281],'cm^-1')),
        HinderedRotor(inertia=(0.170151,'amu*angstrom^2'), symmetry=1, barrier=(9.9625,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169409,'amu*angstrom^2'), symmetry=1, barrier=(9.96659,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.354766,'amu*angstrom^2'), symmetry=1, barrier=(20.8764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.169475,'amu*angstrom^2'), symmetry=1, barrier=(9.9693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.622523,0.0792007,-0.000115128,8.95714e-08,-2.78323e-11,64556.2,31.9954], Tmin=(100,'K'), Tmax=(788.797,'K')), NASAPolynomial(coeffs=[11.331,0.0248908,-1.18366e-05,2.26148e-09,-1.56766e-13,62867,-17.1247], Tmin=(788.797,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.777,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C([C]=[C]O)O[O](23338)',
    structure = SMILES('C=[C]C([C]=[C]O)O[O]'),
    E0 = (641.061,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1670,1685,1700,300,370,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.579466,'amu*angstrom^2'), symmetry=1, barrier=(13.3231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.579589,'amu*angstrom^2'), symmetry=1, barrier=(13.3259,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.579603,'amu*angstrom^2'), symmetry=1, barrier=(13.3262,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.579573,'amu*angstrom^2'), symmetry=1, barrier=(13.3255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.360761,0.0835272,-0.00012158,9.00148e-08,-2.60677e-11,77229.7,34.4247], Tmin=(100,'K'), Tmax=(851.524,'K')), NASAPolynomial(coeffs=[14.0051,0.019429,-8.66006e-06,1.60298e-09,-1.09032e-13,74906.2,-29.2076], Tmin=(851.524,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(641.061,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C([C]=[C][O])OO(23339)',
    structure = SMILES('C=[C]C([C][C]=O)OO'),
    E0 = (621.614,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1310,387.5,850,1000,1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,200,800,1600],'cm^-1')),
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
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.285064,0.0879339,-0.00013746,1.131e-07,-3.69273e-11,74890.8,32.9176], Tmin=(100,'K'), Tmax=(786.523,'K')), NASAPolynomial(coeffs=[11.9002,0.0257486,-1.2925e-05,2.50763e-09,-1.74542e-13,73160,-19.717], Tmin=(786.523,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(621.614,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(261.906,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C][C](C=[C]O)O[O](23340)',
    structure = SMILES('[CH2][C]=C(C=[C]O)O[O]'),
    E0 = (507.351,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.434317,0.0794725,-0.000104939,6.72196e-08,-1.56016e-11,61147.9,31.6728], Tmin=(100,'K'), Tmax=(798.77,'K')), NASAPolynomial(coeffs=[15.259,0.0158708,-5.47554e-06,8.75191e-10,-5.43995e-14,58440.3,-38.6409], Tmin=(798.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(507.351,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(ROOJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]C(C=[C]O)O[O](23341)',
    structure = SMILES('[CH]=[C]C(C=[C]O)O[O]'),
    E0 = (650.316,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([492.5,1135,1000,3615,1277.5,1000,1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.656126,'amu*angstrom^2'), symmetry=1, barrier=(15.0856,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.655424,'amu*angstrom^2'), symmetry=1, barrier=(15.0695,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.656023,'amu*angstrom^2'), symmetry=1, barrier=(15.0833,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.655957,'amu*angstrom^2'), symmetry=1, barrier=(15.0817,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (112.083,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.335817,0.0822554,-0.000112713,7.70492e-08,-2.04817e-11,78345.4,34.4441], Tmin=(100,'K'), Tmax=(928.216,'K')), NASAPolynomial(coeffs=[15.4712,0.0170315,-7.30982e-06,1.34585e-09,-9.20431e-14,75535.7,-37.4476], Tmin=(928.216,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(650.316,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-(Cds-Cd)H) + group(O2s-OsH) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(ROOJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
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
    E0 = (525.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (525.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (775.162,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (917.743,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1074.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (985.411,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (528.235,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (533.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (533.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (656.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (656.893,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (553.159,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (548.167,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (588.706,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (612.223,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (678.357,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (703.603,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (689.446,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (525.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (573.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (860.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (984.194,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (683.505,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (688.085,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (639.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (690.1,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (678.472,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (781.515,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (653.438,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (585.074,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (525.306,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (619.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (533.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (532.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (608.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (703.963,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (603.553,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (710.87,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (625.185,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (746.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (711.645,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (791.407,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (705.504,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (648.582,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (801.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['HCCO(2227)', 'C=C=CO[O](16806)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['O2(S)(666)', 'C=[C]C=C[C]=O(16163)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', 'C=[C]C([O])[CH][C]=O(23311)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 3 used for O_rad/NonDe;O_birad
Exact match found for rate rule [O_rad/NonDe;O_birad]
Euclidian distance = 0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]=[C][O](6861)', '[CH2][C]=CO[O](16807)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H2CC(T)(1341)', '[O]O[CH][CH][C]=O(23312)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', 'C=[C]C([C][C]=O)O[O](23313)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=C1C([C]=O)C1O[O](23314)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=[C]C1OOC1[C]=O(23315)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=C1OOC1[CH][C]=O(23316)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['O(4)', 'C=[C]C1OC1[C]=O(19001)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2OO_S;Y_rad_intra;OOJ]
Euclidian distance = 0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['O(4)', 'C=C1OC1[CH][C]=O(23317)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.30161e+09,'s^-1'), n=1.08, Ea=(131.587,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2OO_S;Y_rad_intra;OOJ] for rate rule [R2OO_S;Cd_rad_out;OOJ]
Euclidian distance = 1.0
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=C=C(C[C]=O)O[O](21754)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=CC(=C[C]=O)O[O](23318)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=[C]C(=C[C]=O)OO(23319)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', 'C=[C]C(=C[C]=O)O[O](23320)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(8.22e+08,'cm^3/(mol*s)'), n=1.533, Ea=(7.77387,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 192 used for Cd_R;HJ
Exact match found for rate rule [Cd_R;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', 'C#CC([CH][C]=O)O[O](23321)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.255e+11,'cm^3/(mol*s)'), n=1.005, Ea=(13.1503,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 138 used for Ct-H_Ct-Cs;HJ
Exact match found for rate rule [Ct-H_Ct-Cs;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['H2CC(T)(1341)', '[O]OC=C[C]=O(23322)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C][O](6861)', 'C=C=CO[O](16806)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.00817,'m^3/(mol*s)'), n=1.99965, Ea=(13.9682,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cds_Ca;YJ] for rate rule [Cds_Ca;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['O2(2)', 'C=[C]C=C[C]=O(16163)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(0.000167406,'m^3/(mol*s)'), n=2.98833, Ea=(203.878,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O2b]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond
Ea raised from 197.8 to 203.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['O2(2)', 'C=[C][CH][CH][C]=O(16165)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.21695e+08,'m^3/(mol*s)'), n=-0.428622, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-4R!H-R_Ext-2R-R_N-2R-inRing]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', 'C=[C][C]([CH][C]=O)O[O](23323)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH]=[C]C([CH][C]=O)O[O](23324)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(9.10287e+13,'m^3/(mol*s)'), n=-2.74437, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=[C][C](C[C]=O)O[O](21756)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=C[C]([CH][C]=O)O[O](23325)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out_Cd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=CC([CH][C]=O)O[O](23326)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=[C][C]([CH][C]=O)OO(23327)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.83109e+08,'s^-1'), n=1.32333, Ea=(164.794,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS_O;O_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=[C][C]([CH]C=O)O[O](23328)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.10706e+06,'s^-1'), n=1.88838, Ea=(153.166,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Y_rad_out;XH_out] for rate rule [R3HJ;CO_rad_out;XH_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]C(C[C]=O)O[O](21757)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.73633e+06,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]C([CH][C]=O)OO(23329)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C]C([CH]C=O)O[O](23330)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5Hall;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=[C]C(C=C=O)O[O](21189)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=C([CH][C]=O)[CH]O[O](20759)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(8.66e+11,'s^-1'), n=0.438, Ea=(94.4747,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HR!H)CJ;CdsJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=C1C([O])=CC1O[O](23331)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=[C]C1C=C([O])OO1(23332)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSD;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['O(4)', 'C=[C]C1C=C([O])O1(23333)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.11355e+11,'s^-1'), n=0, Ea=(82.8981,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3OO;Y_rad_intra;OO] for rate rule [R3OO_DS;Y_rad_intra;OOJ]
Euclidian distance = 1.4142135623730951
family: Cyclic_Ether_Formation"""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['HO2(9)', '[CH2][C]=CC#C[O](18002)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(3.63e+09,'s^-1'), n=1.11, Ea=(178.657,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R2OO_0H]
Euclidian distance = 0
family: HO2_Elimination_from_PeroxyRadical"""),
)

reaction(
    label = 'reaction37',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=C=C(C=C[O])O[O](23334)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction38',
    reactants = ['H(3)', 'C=[C]C(C#C[O])O[O](23335)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction39',
    reactants = ['HCCO(2227)', '[CH2][C]=CO[O](16807)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C]C([C]=C[O])O[O](23336)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_Cs;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=CC([C]=[C][O])O[O](23337)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(3.36833e+07,'s^-1'), n=1.41, Ea=(175.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS_Cs;Y_rad_out;Cd_H_out_doubleC] for rate rule [R3H_SS_Cs;Cd_rad_out;Cd_H_out_doubleC]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=[C]C([C]=[C]O)O[O](23338)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_Cs;XH_out] for rate rule [R3HJ;Cd_rad_out_Cs;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=[C]C([C]=[C][O])OO(23339)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.00504539,'s^-1'), n=3.83, Ea=(83.8892,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;O_H_out] for rate rule [R4H_SSS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=[C]C(C=[C][O])O[O](21197)'],
    products = ['C=[C][C](C=[C]O)O[O](23340)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(520772,'s^-1'), n=1.9431, Ea=(123.276,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;O_rad_out;XH_out] for rate rule [R4HJ_1;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=[C]C(C=[C]O)O[O](23341)'],
    products = ['C=[C]C(C=[C][O])O[O](21197)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4908',
    isomers = [
        'C=[C]C(C=[C][O])O[O](21197)',
    ],
    reactants = [
        ('HCCO(2227)', 'C=C=CO[O](16806)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4908',
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


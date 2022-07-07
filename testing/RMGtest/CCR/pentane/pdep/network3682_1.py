species(
    label = '[CH2]C(=C)OC=[C][O](12363)',
    structure = SMILES('[CH2]C(=C)O[CH][C]=O'),
    E0 = (213.628,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3025,407.5,1350,352.5,1855,455,950,374.647,374.657,374.67],'cm^-1')),
        HinderedRotor(inertia=(0.183864,'amu*angstrom^2'), symmetry=1, barrier=(18.3129,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183872,'amu*angstrom^2'), symmetry=1, barrier=(18.313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183878,'amu*angstrom^2'), symmetry=1, barrier=(18.3131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.183861,'amu*angstrom^2'), symmetry=1, barrier=(18.3132,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0773396,0.075818,-8.46e-05,4.428e-08,-8.68842e-12,25852.1,27.684], Tmin=(100,'K'), Tmax=(1380.73,'K')), NASAPolynomial(coeffs=[21.7449,0.00552447,-5.49069e-07,-1.35531e-11,3.40056e-15,20500.3,-82.1939], Tmin=(1380.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.628,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(CCsJOC(O)) + radical(CsCJ=O)"""),
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
    label = 'allene(458)',
    structure = SMILES('C=C=C'),
    E0 = (175.934,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3000,3050,3100,1330,1430,900,1050,1000,1050,1600,1700,540,610,2055],'cm^-1')),
    ],
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(2746.46,'J/mol'), sigma=(4.78521,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=428.99 K, Pc=56.87 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.37447,0.0070462,2.78306e-05,-3.99445e-08,1.55729e-11,21188.6,7.62048], Tmin=(100,'K'), Tmax=(949.705,'K')), NASAPolynomial(coeffs=[6.79956,0.00959979,-3.02068e-06,5.37827e-10,-3.92606e-14,19772.3,-12.7582], Tmin=(949.705,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(175.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(157.975,'J/(mol*K)'), label="""allene""", comment="""Thermo library: DFT_QCI_thermo"""),
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
    label = '[CH2]C(=C)[O](2372)',
    structure = SMILES('[CH2]C(=C)[O]'),
    E0 = (110.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1477.83],'cm^-1')),
        HinderedRotor(inertia=(0.530917,'amu*angstrom^2'), symmetry=1, barrier=(12.2068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.69661,0.0242908,-7.63704e-06,-9.31802e-09,6.21014e-12,13330.3,14.3522], Tmin=(100,'K'), Tmax=(924.139,'K')), NASAPolynomial(coeffs=[8.66416,0.00984172,-2.65662e-06,4.14898e-10,-2.77492e-14,11741.3,-16.5963], Tmin=(924.139,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(110.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(C=C(C)OJ) + radical(C=C(O)CJ)"""),
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
    label = 'C=[C]O[CH][C]=O(2657)',
    structure = SMILES('C=[C]O[CH][C]=O'),
    E0 = (336.247,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,250.522,250.523,250.677],'cm^-1')),
        HinderedRotor(inertia=(0.452375,'amu*angstrom^2'), symmetry=1, barrier=(20.1537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.452325,'amu*angstrom^2'), symmetry=1, barrier=(20.1537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.452364,'amu*angstrom^2'), symmetry=1, barrier=(20.1532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3435.15,'J/mol'), sigma=(5.65948,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=536.56 K, Pc=43 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.874181,0.0573741,-6.57983e-05,3.45872e-08,-6.7543e-12,40563.4,25.3964], Tmin=(100,'K'), Tmax=(1406.02,'K')), NASAPolynomial(coeffs=[18.1708,0.00165655,5.89016e-07,-1.83742e-10,1.37715e-14,36343,-61.6553], Tmin=(1406.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(336.247,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(C=CJO) + radical(CsCJ=O)"""),
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
    label = '[CH]C(=C)O[CH][C]=O(16731)',
    structure = SMILES('[CH]C(=C)O[CH][C]=O'),
    E0 = (425.395,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1855,455,950,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.406579,0.0699742,-6.29845e-05,2.14761e-08,-5.74475e-13,51300.7,27.0561], Tmin=(100,'K'), Tmax=(988.484,'K')), NASAPolynomial(coeffs=[18.8243,0.0135836,-4.93803e-06,8.91261e-10,-6.32866e-14,46773.4,-66.0679], Tmin=(988.484,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.395,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(AllylJ2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]C(=C)O[C]=[C][O](16679)',
    structure = SMILES('[CH2]C(=C)O[C][C]=O'),
    E0 = (480.689,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,350.716,350.822,350.867,351.049],'cm^-1')),
        HinderedRotor(inertia=(0.224268,'amu*angstrom^2'), symmetry=1, barrier=(19.5967,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224201,'amu*angstrom^2'), symmetry=1, barrier=(19.5947,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224277,'amu*angstrom^2'), symmetry=1, barrier=(19.5978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.224397,'amu*angstrom^2'), symmetry=1, barrier=(19.5983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.118919,0.0718465,-8.09022e-05,4.20981e-08,-8.19341e-12,57964.8,26.2397], Tmin=(100,'K'), Tmax=(1387.6,'K')), NASAPolynomial(coeffs=[21.5449,0.00362123,-1.65947e-07,-4.69625e-11,4.32614e-15,52640.6,-81.9056], Tmin=(1387.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(480.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = 'C=C1CC([C]=O)O1(16732)',
    structure = SMILES('C=C1CC([C]=O)O1'),
    E0 = (-11.9106,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68005,0.0270834,6.46959e-05,-1.22222e-07,5.48685e-11,-1326.52,20.4503], Tmin=(100,'K'), Tmax=(904.825,'K')), NASAPolynomial(coeffs=[22.1921,-0.000337347,5.28595e-06,-1.18377e-09,7.7974e-14,-7627.97,-90.7663], Tmin=(904.825,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-11.9106,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(2methyleneoxetane) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][CH]O[C]1CC1(16733)',
    structure = SMILES('O=[C][CH]O[C]1CC1'),
    E0 = (302.042,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.21285,0.06026,-6.07138e-05,3.17869e-08,-6.68374e-12,36428.5,24.7247], Tmin=(100,'K'), Tmax=(1145.78,'K')), NASAPolynomial(coeffs=[12.5269,0.0207616,-9.0042e-06,1.69975e-09,-1.18955e-13,33835.9,-31.399], Tmin=(1145.78,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(302.042,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + ring(Cyclopropane) + radical(C2CsJOCs) + radical(CCsJOCs) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2][C]1CC([C]=O)O1(16734)',
    structure = SMILES('[CH2][C]1CC([C]=O)O1'),
    E0 = (301.798,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.90231,0.0591194,-6.1031e-05,3.38332e-08,-7.11736e-12,36417.2,27.3834], Tmin=(100,'K'), Tmax=(1361.74,'K')), NASAPolynomial(coeffs=[12.5554,0.0154669,-2.56725e-06,1.29866e-10,3.0551e-15,34117.2,-29.2265], Tmin=(1361.74,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.798,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(C2CsJOCs) + radical(CJC(C)OC) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C1C[C]([O])[CH]O1(16717)',
    structure = SMILES('[CH2][C]1CC([O])=CO1'),
    E0 = (169.025,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.403613,0.0682858,-7.07087e-05,3.49118e-08,-6.2124e-12,20512.3,26.872], Tmin=(100,'K'), Tmax=(1701.97,'K')), NASAPolynomial(coeffs=[17.5445,0.00574171,2.3593e-06,-7.38625e-10,5.67898e-14,17352.1,-60.5979], Tmin=(1701.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(169.025,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ) + radical(C2CsJOC(O)) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2]C1([CH2])OC1[C]=O(16735)',
    structure = SMILES('[CH2]C1([CH2])OC1[C]=O'),
    E0 = (319.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0952377,0.0800894,-0.000102235,6.29005e-08,-1.41589e-11,38541.4,27.4665], Tmin=(100,'K'), Tmax=(1297.06,'K')), NASAPolynomial(coeffs=[18.0782,0.00700935,1.97989e-06,-7.89878e-10,6.87084e-14,35260,-59.4127], Tmin=(1297.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(Cs-CsCsCsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1([CH2])O[CH]C1=O(16736)',
    structure = SMILES('[CH2]C1([CH2])OC=C1[O]'),
    E0 = (257.949,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.5936,0.0873954,-0.00010007,5.10138e-08,-9.25634e-12,31256.4,29.4647], Tmin=(100,'K'), Tmax=(1662.04,'K')), NASAPolynomial(coeffs=[25.0328,-0.00363412,6.40547e-06,-1.45022e-09,1.02529e-13,26127.6,-101.323], Tmin=(1662.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(257.949,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(CJC(C)OC) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]=C(459)',
    structure = SMILES('[CH2][C]=C'),
    E0 = (395.31,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370],'cm^-1')),
        HinderedRotor(inertia=(0.110742,'amu*angstrom^2'), symmetry=1, barrier=(23.6192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.28575,0.0103679,1.6796e-05,-2.70041e-08,1.05227e-11,47575.3,10.4038], Tmin=(100,'K'), Tmax=(980.899,'K')), NASAPolynomial(coeffs=[6.52804,0.0104286,-3.60844e-06,6.68464e-10,-4.85509e-14,46300.2,-8.43297], Tmin=(980.899,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(395.31,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(203.705,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Allyl_P) + radical(Cds_S)"""),
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
    label = '[CH]=C(C)O[CH][C]=O(16737)',
    structure = SMILES('[CH]=C(C)O[CH][C]=O'),
    E0 = (301.807,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,1855,455,950,211.827,211.953],'cm^-1')),
        HinderedRotor(inertia=(0.522268,'amu*angstrom^2'), symmetry=1, barrier=(17.2231,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.519487,'amu*angstrom^2'), symmetry=1, barrier=(17.2194,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.540098,'amu*angstrom^2'), symmetry=1, barrier=(17.2226,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.538187,'amu*angstrom^2'), symmetry=1, barrier=(17.2196,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.195233,0.0751681,-8.54675e-05,4.6265e-08,-9.5576e-12,36443.3,27.1106], Tmin=(100,'K'), Tmax=(1200.9,'K')), NASAPolynomial(coeffs=[19.997,0.00921047,-3.08084e-06,5.28075e-10,-3.60518e-14,31687.4,-72.0463], Tmin=(1200.9,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(301.807,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(Cds_P) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=C([CH2])OC[C]=O(16738)',
    structure = SMILES('[CH]C(=C)OC[C]=O'),
    E0 = (231.468,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.768587,0.0606288,-3.62263e-05,-1.09989e-09,5.86876e-12,27964.7,26.2838], Tmin=(100,'K'), Tmax=(1004.39,'K')), NASAPolynomial(coeffs=[16.4407,0.0198668,-7.68741e-06,1.42071e-09,-1.00984e-13,23724.3,-54.8309], Tmin=(1004.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.468,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CsCJ=O) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C([CH2])OC=C[O](16682)',
    structure = SMILES('[CH]C(=C)OC=C[O]'),
    E0 = (231.377,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.7735,0.0569686,-1.59056e-05,-3.00708e-08,1.84523e-11,27957.2,27.0028], Tmin=(100,'K'), Tmax=(942.923,'K')), NASAPolynomial(coeffs=[18.838,0.0150091,-4.31315e-06,7.2999e-10,-5.33237e-14,23009.1,-67.2595], Tmin=(942.923,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.377,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=C)OC=C=O(15582)',
    structure = SMILES('[CH2]C(=C)OC=C=O'),
    E0 = (59.6221,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.888121,0.062912,-6.52323e-05,3.45837e-08,-7.17784e-12,7287.63,24.873], Tmin=(100,'K'), Tmax=(1181.79,'K')), NASAPolynomial(coeffs=[14.7666,0.0159373,-5.60892e-06,9.4907e-10,-6.26278e-14,4007.35,-44.4009], Tmin=(1181.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(59.6221,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=C(O)CJ)"""),
)

species(
    label = 'C=C1CC([O])=CO1(16687)',
    structure = SMILES('C=C1CC([O])=CO1'),
    E0 = (-81.0267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.20536,0.0142031,9.26325e-05,-1.42491e-07,5.92793e-11,-9657.21,21.7295], Tmin=(100,'K'), Tmax=(922.337,'K')), NASAPolynomial(coeffs=[19.4372,0.004647,2.1788e-06,-4.97853e-10,2.58635e-14,-15608.2,-75.0395], Tmin=(922.337,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-81.0267,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(282.692,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C([CH2])O[C]=[C][O](12276)',
    structure = SMILES('[CH2]C([CH2])O[C]=[C][O]'),
    E0 = (611.248,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1670,1700,300,440,334.935,334.999,335.06,335.076],'cm^-1')),
        HinderedRotor(inertia=(0.189687,'amu*angstrom^2'), symmetry=1, barrier=(15.1067,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189682,'amu*angstrom^2'), symmetry=1, barrier=(15.1073,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189646,'amu*angstrom^2'), symmetry=1, barrier=(15.1072,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.189712,'amu*angstrom^2'), symmetry=1, barrier=(15.1068,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0370276,0.0883683,-0.000122794,8.11468e-08,-2.00649e-11,73662.1,31.415], Tmin=(100,'K'), Tmax=(862.789,'K')), NASAPolynomial(coeffs=[18.4163,0.0113395,-3.69396e-06,5.69485e-10,-3.46432e-14,70160.5,-56.7268], Tmin=(862.789,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(611.248,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](C)O[C]=[C][O](16739)',
    structure = SMILES('[CH2][C](C)O[C]=[C][O]'),
    E0 = (600.954,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,3000,3100,440,815,1455,1000,1670,1700,300,440,285.863,285.878,285.885,285.915],'cm^-1')),
        HinderedRotor(inertia=(0.00206319,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204698,'amu*angstrom^2'), symmetry=1, barrier=(11.8691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.20466,'amu*angstrom^2'), symmetry=1, barrier=(11.8691,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.204659,'amu*angstrom^2'), symmetry=1, barrier=(11.8693,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271081,0.0834272,-0.000117752,8.2303e-08,-2.22445e-11,72411.1,32.7485], Tmin=(100,'K'), Tmax=(915.627,'K')), NASAPolynomial(coeffs=[15.8482,0.0153709,-6.25033e-06,1.11122e-09,-7.4101e-14,69558.8,-41.0277], Tmin=(915.627,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(600.954,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C2CsJOC(O)) + radical(CJC(C)OC) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=C1CO[C][CH]O1(16668)',
    structure = SMILES('[CH2][C]1CO[C]=CO1'),
    E0 = (330.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.2493,0.0692251,-7.50986e-05,3.94518e-08,-7.45547e-12,39958.9,24.3177], Tmin=(100,'K'), Tmax=(1607.26,'K')), NASAPolynomial(coeffs=[16.6871,0.00635141,2.92007e-06,-9.31418e-10,7.37103e-14,37191.4,-57.1004], Tmin=(1607.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(330.791,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(278.535,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(23dihydro14dioxin) + radical(C2CsJOC(O)) + radical(CJC(C)OC) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1([CH2])O[C]=CO1(16677)',
    structure = SMILES('[CH2]C1([CH2])O[C]=CO1'),
    E0 = (282.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.145756,0.0447298,7.72104e-05,-1.80564e-07,8.7013e-11,34190.8,22.7772], Tmin=(100,'K'), Tmax=(896.094,'K')), NASAPolynomial(coeffs=[41.009,-0.0296675,2.09469e-05,-4.19691e-09,2.81974e-13,22530.8,-194.078], Tmin=(896.094,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(282.81,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(274.378,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsCsOsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=C)O[C]=C=O(16678)',
    structure = SMILES('[CH2]C(=C)O[C]=C=O'),
    E0 = (299.366,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,2120,512.5,787.5,180,180,3801.37],'cm^-1')),
        HinderedRotor(inertia=(0.835698,'amu*angstrom^2'), symmetry=1, barrier=(19.2143,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0267062,'amu*angstrom^2'), symmetry=1, barrier=(0.614028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.62297,'amu*angstrom^2'), symmetry=1, barrier=(60.3074,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30929,0.0621195,-8.29787e-05,6.13337e-08,-1.82353e-11,36099.6,26.4948], Tmin=(100,'K'), Tmax=(822.427,'K')), NASAPolynomial(coeffs=[9.59401,0.0218209,-9.47096e-06,1.74088e-09,-1.18287e-13,34737,-11.8536], Tmin=(822.427,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(299.366,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(C=C(O)CJ) + radical(C=CJO)"""),
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
    label = '[CH2]C(=C)O[C]=C[O](12273)',
    structure = SMILES('[CH2]C(=C)O[C]=C[O]'),
    E0 = (259.353,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,342.139,342.927,342.942,343.201],'cm^-1')),
        HinderedRotor(inertia=(0.00143456,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261575,'amu*angstrom^2'), symmetry=1, barrier=(21.8537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.26125,'amu*angstrom^2'), symmetry=1, barrier=(21.8509,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3877.8,'J/mol'), sigma=(6.32548,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=605.70 K, Pc=34.77 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.476684,0.0648013,-6.50412e-05,3.21076e-08,-6.03414e-12,31330.6,30.0902], Tmin=(100,'K'), Tmax=(1435.24,'K')), NASAPolynomial(coeffs=[17.8525,0.0106964,-2.56002e-06,3.28374e-10,-1.84203e-14,26927.9,-57.9786], Tmin=(1435.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(259.353,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=C(O)CJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=C)O[C]=[C]O(16680)',
    structure = SMILES('[CH2]C(=C)O[C]=[C]O'),
    E0 = (357.635,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1670,1700,300,440,347.719,347.935,348.115],'cm^-1')),
        HinderedRotor(inertia=(0.00139496,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170988,'amu*angstrom^2'), symmetry=1, barrier=(14.6859,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.171043,'amu*angstrom^2'), symmetry=1, barrier=(14.6847,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170871,'amu*angstrom^2'), symmetry=1, barrier=(14.6846,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.741759,0.0681106,-7.67608e-05,4.11152e-08,-7.66874e-12,43134.1,30.8524], Tmin=(100,'K'), Tmax=(911.301,'K')), NASAPolynomial(coeffs=[15.5935,0.0134452,-4.10393e-06,6.34953e-10,-4.00733e-14,39990.3,-41.8169], Tmin=(911.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(357.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(O)CJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=C(C)O[C]=[C][O](16683)',
    structure = SMILES('C=C(C)O[C][C]=O'),
    E0 = (321.772,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,280.524,284.548,284.658,286.465],'cm^-1')),
        HinderedRotor(inertia=(0.324927,'amu*angstrom^2'), symmetry=1, barrier=(17.9539,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.318903,'amu*angstrom^2'), symmetry=1, barrier=(17.9464,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.312352,'amu*angstrom^2'), symmetry=1, barrier=(17.9537,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.317367,'amu*angstrom^2'), symmetry=1, barrier=(17.9531,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.550874,0.0667057,-6.19389e-05,2.29225e-08,-1.69147e-12,38832.4,24.6887], Tmin=(100,'K'), Tmax=(1024,'K')), NASAPolynomial(coeffs=[18.8931,0.0112185,-4.33352e-06,8.32184e-10,-6.12933e-14,34228.6,-68.3741], Tmin=(1024,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(321.772,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=C([CH2])OC=[C]O(16740)',
    structure = SMILES('[CH]C(=C)OC=[C]O'),
    E0 = (329.658,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.0695705,0.0731116,-7.14676e-05,3.44154e-08,-6.29809e-12,39809.3,31.7567], Tmin=(100,'K'), Tmax=(1499.91,'K')), NASAPolynomial(coeffs=[19.5086,0.0127645,-2.98085e-06,3.58733e-10,-1.8897e-14,34851.3,-67.5832], Tmin=(1499.91,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.658,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C(=O)CC=[C][O](16741)',
    structure = SMILES('C=C([O])C[CH][C]=O'),
    E0 = (127.635,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11566,0.0590678,-5.84953e-05,2.90277e-08,-5.26602e-12,15458.6,26.8439], Tmin=(100,'K'), Tmax=(986.512,'K')), NASAPolynomial(coeffs=[13.268,0.0171181,-5.84677e-06,9.74628e-10,-6.40159e-14,12704.6,-33.4252], Tmin=(986.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.635,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-CdsCsOs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(CCJCHO) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2][C][CH2](15952)',
    structure = SMILES('[CH2][C][CH2]'),
    E0 = (738.707,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,618.688],'cm^-1')),
        HinderedRotor(inertia=(0.00838111,'amu*angstrom^2'), symmetry=1, barrier=(2.27602,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0989004,'amu*angstrom^2'), symmetry=1, barrier=(2.27392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (40.0639,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.05819,0.0178835,-7.35728e-06,-7.95768e-10,1.01791e-12,88882.1,14.263], Tmin=(100,'K'), Tmax=(1175.48,'K')), NASAPolynomial(coeffs=[6.73198,0.0102413,-3.80615e-06,7.07005e-10,-4.9648e-14,87682.7,-5.48288], Tmin=(1175.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(738.707,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(RCCJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]([CH2])O[C]=C[O](12275)',
    structure = SMILES('[CH2][C]([CH2])O[C]=C[O]'),
    E0 = (571.72,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3010,987.5,1337.5,450,1655,1685,370,243.193,243.194,243.194,243.194],'cm^-1')),
        HinderedRotor(inertia=(0.339519,'amu*angstrom^2'), symmetry=1, barrier=(14.2494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339519,'amu*angstrom^2'), symmetry=1, barrier=(14.2493,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339519,'amu*angstrom^2'), symmetry=1, barrier=(14.2494,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.339518,'amu*angstrom^2'), symmetry=1, barrier=(14.2494,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.339911,0.0950126,-0.000136347,9.06793e-08,-2.23344e-11,68918.8,31.565], Tmin=(100,'K'), Tmax=(856.72,'K')), NASAPolynomial(coeffs=[20.4029,0.00834465,-2.42675e-06,3.36963e-10,-1.87563e-14,64991.1,-67.4797], Tmin=(856.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.72,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C2CsJOC(O)) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]([CH2])[O](15953)',
    structure = SMILES('[CH2][C]([CH2])[O]'),
    E0 = (539.329,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,3014.49],'cm^-1')),
        HinderedRotor(inertia=(0.150723,'amu*angstrom^2'), symmetry=1, barrier=(4.5418,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.726258,'amu*angstrom^2'), symmetry=1, barrier=(21.9126,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (56.0633,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.16038,0.0431609,-6.34706e-05,5.06292e-08,-1.60119e-11,64930,16.4324], Tmin=(100,'K'), Tmax=(788.798,'K')), NASAPolynomial(coeffs=[7.96943,0.0132457,-5.71329e-06,1.07933e-09,-7.46744e-14,64027.8,-10.1245], Tmin=(788.798,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(539.329,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(270.22,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CJCO) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([CH2])O[C]=C=O(12274)',
    structure = SMILES('[CH2]C([CH2])O[C]=C=O'),
    E0 = (411.517,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3033.33,3066.67,3100,415,465,780,850,1435,1475,900,1100,1685,370,2120,512.5,787.5,263.846,263.849,263.869],'cm^-1')),
        HinderedRotor(inertia=(0.284397,'amu*angstrom^2'), symmetry=1, barrier=(14.0531,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284387,'amu*angstrom^2'), symmetry=1, barrier=(14.0529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284378,'amu*angstrom^2'), symmetry=1, barrier=(14.0534,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.284542,'amu*angstrom^2'), symmetry=1, barrier=(14.0532,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (97.092,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.023354,0.0911219,-0.000139062,1.04376e-07,-3.00824e-11,49636.5,27.6286], Tmin=(100,'K'), Tmax=(905.225,'K')), NASAPolynomial(coeffs=[16.1911,0.0151179,-5.90228e-06,9.92389e-10,-6.25059e-14,46879.4,-47.9964], Tmin=(905.225,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.517,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(266.063,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-CsCsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CJC(C)OC) + radical(CJC(C)OC) + radical(C=CJO)"""),
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
    E0 = (213.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (535.453,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (774.049,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (637.187,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (692.481,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (221.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (444.844,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (340.267,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (287.882,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (319.655,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (293.585,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (351.131,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (599.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (431.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (275.776,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (264.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (213.628,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (220.74,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (674.649,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (625.927,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (330.791,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (307.834,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (526.268,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (384.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (322.188,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (463.49,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (507.98,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (404.197,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (481.179,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (527.428,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (663.16,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (594.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (704.52,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (553.494,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['OCHCO(3676)', 'allene(458)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C][O](6861)', '[CH2]C(=C)[O](2372)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(20)', 'C=[C]O[CH][C]=O(2657)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]C(=C)O[CH][C]=O(16731)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(=C)O[C]=[C][O](16679)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['C=C1CC([C]=O)O1(16732)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(3.24e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['O=[C][CH]O[C]1CC1(16733)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs] for rate rule [R3_D;doublebond_intra_secNd_2H;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['[CH2][C]1CC([C]=O)O1(16734)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.66736e+08,'s^-1'), n=1.1093, Ea=(126.639,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['C=C1C[C]([O])[CH]O1(16717)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.10072e+11,'s^-1'), n=0.23641, Ea=(74.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['[CH2]C1([CH2])OC1[C]=O(16735)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.17979e+10,'s^-1'), n=0.520585, Ea=(106.028,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['[CH2]C1([CH2])O[CH]C1=O(16736)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(9.11991e+08,'s^-1'), n=0.766793, Ea=(79.957,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;doublebond_intra_2H;radadd_intra] for rate rule [R5;doublebond_intra_2H;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OCHCO(3676)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.41198,'m^3/(mol*s)'), n=2.03228, Ea=(31.3669,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CJ] + [Od_R;YJ] for rate rule [Od_R;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C]=C[O](9592)', '[CH2][C]=C(459)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=C(C)O[CH][C]=O(16737)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(720,'s^-1'), n=2.932, Ea=(129.315,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 57 used for R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H
Exact match found for rate rule [R3H_DS;Cd_rad_out_singleH;Cs_H_out_2H]
Euclidian distance = 0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=C([CH2])OC[C]=O(16738)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=C([CH2])OC=C[O](16682)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['[CH2]C(=C)OC=C=O(15582)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['C=C1CC([O])=CO1(16687)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.552e+10,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH2]C([CH2])O[C]=[C][O](12276)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2][C](C)O[C]=[C][O](16739)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['C=C1CO[C][CH]O1(16668)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.9903e+10,'s^-1'), n=0.243684, Ea=(117.164,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 115.1 to 117.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['[CH2]C1([CH2])O[C]=CO1(16677)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.79781e+09,'s^-1'), n=0.611011, Ea=(94.2067,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;doublebond_intra_2H;radadd_intra_O] + [R6;doublebond_intra_2H;radadd_intra] for rate rule [R6;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2]C(=C)O[C]=C=O(16678)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[O][C]=C[O](9592)', 'allene(458)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(0.000692465,'m^3/(mol*s)'), n=2.43792, Ea=(4.66472,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Cdd_Cds;OJ_sec] for rate rule [Ca_Cds-HH;O_rad/OneDe]
Euclidian distance = 2.8284271247461903
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['HCCO(2227)', '[CH2]C(=C)[O](2372)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(260000,'m^3/(mol*s)'), n=0, Ea=(46.5868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Ct_Ct;OJ_sec] for rate rule [Ct_Ct;O_rad/OneDe]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH2]C(=C)O[C]=C[O](12273)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2]C(=C)O[C]=[C]O(16680)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C(C)O[C]=[C][O](16683)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(32400,'s^-1'), n=2.04, Ea=(82.4248,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SS(Cd)S;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_SS(Cd)S;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=C([CH2])OC=[C]O(16740)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_4;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2]C(=C)OC=[C][O](12363)'],
    products = ['[CH2]C(=O)CC=[C][O](16741)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(14080,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: ketoenol"""),
)

reaction(
    label = 'reaction31',
    reactants = ['OCHCO(3676)', '[CH2][C][CH2](15952)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]([CH2])O[C]=C[O](12275)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['HCCO(2227)', '[CH2][C]([CH2])[O](15953)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(106477,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH2]C([CH2])O[C]=C=O(12274)'],
    products = ['[CH2]C(=C)OC=[C][O](12363)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_double;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #3682',
    isomers = [
        '[CH2]C(=C)OC=[C][O](12363)',
    ],
    reactants = [
        ('OCHCO(3676)', 'allene(458)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #3682',
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


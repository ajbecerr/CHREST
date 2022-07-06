species(
    label = 'C=C([O])O[CH][C]=O(11146)',
    structure = SMILES('[CH2]C(=O)O[CH][C]=O'),
    E0 = (30.111,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.58826,0.0542184,-5.75522e-05,3.11605e-08,-6.79731e-12,3707.42,26.1879], Tmin=(100,'K'), Tmax=(1101,'K')), NASAPolynomial(coeffs=[11.3875,0.0186169,-9.04852e-06,1.79092e-09,-1.2842e-13,1549.63,-22.0307], Tmin=(1101,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.111,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(CJCO) + radical(CsCJ=O)"""),
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
    label = 'C=C([O])[O](2321)',
    structure = SMILES('[CH2]C([O])=O'),
    E0 = (-8.20105,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,840.394,841.338,841.823,842,851.257],'cm^-1')),
        HinderedRotor(inertia=(0.00477063,'amu*angstrom^2'), symmetry=1, barrier=(2.45606,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.60873,0.0103835,4.16405e-06,-6.72688e-09,1.73623e-12,-973.928,12.8424], Tmin=(100,'K'), Tmax=(1564.28,'K')), NASAPolynomial(coeffs=[5.70755,0.0133435,-6.65901e-06,1.28861e-09,-8.86354e-14,-2649.35,-1.47867], Tmin=(1564.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-8.20105,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CCOJ) + radical(CJCO)"""),
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
    label = 'O=[C][CH]O[C]=O(11954)',
    structure = SMILES('O=[C][CH]O[C]=O'),
    E0 = (72.8243,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1850,1860,440,470,900,1000,462.641,462.651],'cm^-1')),
        HinderedRotor(inertia=(0.106215,'amu*angstrom^2'), symmetry=1, barrier=(16.1327,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106211,'amu*angstrom^2'), symmetry=1, barrier=(16.1326,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106216,'amu*angstrom^2'), symmetry=1, barrier=(16.1325,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (85.0382,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71335,0.046063,-5.675e-05,3.1349e-08,-6.33955e-12,8844.8,21.8058], Tmin=(100,'K'), Tmax=(1021.06,'K')), NASAPolynomial(coeffs=[14.3025,0.00233715,-7.29269e-07,1.36019e-10,-1.05478e-14,5982.42,-40.6196], Tmin=(1021.06,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.8243,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsH) + group(Cds-OdOsH) + radical(CCsJOC(O)H) + radical(CsCJ=O) + radical((O)CJOCC)"""),
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
    label = '[CH]=C([O])O[CH][C]=O(11365)',
    structure = SMILES('[CH]C(=O)O[CH][C]=O'),
    E0 = (266.737,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,1855,455,950,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68789,0.052082,-5.72866e-05,3.16143e-08,-6.98358e-12,32163.3,25.5079], Tmin=(100,'K'), Tmax=(1091.2,'K')), NASAPolynomial(coeffs=[11.4212,0.0164037,-8.24327e-06,1.65215e-09,-1.19254e-13,30039.1,-22.2994], Tmin=(1091.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(266.737,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CCsJOC(O)) + radical(CCJ2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH2]C(=O)O[C]=[C][O](11890)',
    structure = SMILES('[CH2]C(=O)O[C]=[C][O]'),
    E0 = (273.678,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1670,1700,300,440,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55014,0.0554133,-6.93558e-05,4.39998e-08,-1.10312e-11,33002.9,29.092], Tmin=(100,'K'), Tmax=(974.915,'K')), NASAPolynomial(coeffs=[11.4285,0.014883,-6.99579e-06,1.35652e-09,-9.60275e-14,31076.8,-18.3143], Tmin=(974.915,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(273.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'O=[C]C1CC(=O)O1(11151)',
    structure = SMILES('O=[C]C1CC(=O)O1'),
    E0 = (-253.138,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.64287,0.0120581,6.96878e-05,-1.06468e-07,4.43658e-11,-30380,22.8899], Tmin=(100,'K'), Tmax=(911.944,'K')), NASAPolynomial(coeffs=[14.079,0.00807599,2.80012e-07,-2.00502e-10,1.11623e-14,-34386.1,-41.7571], Tmin=(911.944,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-253.138,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + ring(Beta-Propiolactone) + radical(CCCJ=O)"""),
)

species(
    label = 'O=[C][CH]O[C]1CO1(11955)',
    structure = SMILES('[O][C]=CO[C]1CO1'),
    E0 = (181.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.809966,0.0696274,-7.68033e-05,3.89059e-08,-6.94012e-12,22032.4,31.9538], Tmin=(100,'K'), Tmax=(1744.21,'K')), NASAPolynomial(coeffs=[17.1914,0.00133728,5.15165e-06,-1.29619e-09,9.45244e-14,19861,-53.1297], Tmin=(1744.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Ethylene_oxide) + radical(C=COJ) + radical(Cs_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]1OC([C]=O)O1(11956)',
    structure = SMILES('[CH2][C]1OC([C]=O)O1'),
    E0 = (165.426,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.81896,0.0389293,-1.44142e-05,-8.97708e-09,5.50515e-12,19982.2,25.5955], Tmin=(100,'K'), Tmax=(1114.77,'K')), NASAPolynomial(coeffs=[13.1927,0.0161136,-7.9285e-06,1.62514e-09,-1.20036e-13,16328.2,-35.5274], Tmin=(1114.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(165.426,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Cyclobutane) + radical(Cs_P) + radical(CJCO) + radical(CsCJ=O)"""),
)

species(
    label = 'C=C1O[CH][C]([O])O1(11882)',
    structure = SMILES('[CH2][C]1O[CH]C(=O)O1'),
    E0 = (55.3158,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.05748,0.0260962,3.59202e-05,-7.12829e-08,3.03812e-11,6737.96,23.8247], Tmin=(100,'K'), Tmax=(959.733,'K')), NASAPolynomial(coeffs=[15.9773,0.00948603,-2.83313e-06,5.89357e-10,-5.04215e-14,2159.2,-52.693], Tmin=(959.733,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(55.3158,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsOsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Cyclopentane) + radical(Cs_P) + radical(CCsJOCs) + radical(CJCO)"""),
)

species(
    label = '[CH2]C1([O])OC1[C]=O(11957)',
    structure = SMILES('[CH2]C1([O])OC1[C]=O'),
    E0 = (176.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.276155,0.0610206,-6.49346e-05,3.32901e-08,-6.18997e-12,21332.3,28.0058], Tmin=(100,'K'), Tmax=(1612.88,'K')), NASAPolynomial(coeffs=[15.8842,0.00580736,1.76424e-06,-6.2391e-10,5.02168e-14,18444.2,-48.0998], Tmin=(1612.88,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.102,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CJC(O)2C) + radical(CCCJ=O)"""),
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
    label = 'C=C([O])OC=C=O(11168)',
    structure = SMILES('[CH2]C(=O)OC=C=O'),
    E0 = (-165.797,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24657,0.057708,-6.33703e-05,3.44275e-08,-7.30925e-12,-19839.1,23.3125], Tmin=(100,'K'), Tmax=(1152.81,'K')), NASAPolynomial(coeffs=[14.0227,0.0133769,-5.68714e-06,1.06908e-09,-7.49852e-14,-22784.8,-40.1418], Tmin=(1152.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-165.797,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CJCO)"""),
)

species(
    label = '[O]C1=COC(=O)C1(11818)',
    structure = SMILES('[O]C1=COC(=O)C1'),
    E0 = (-356.976,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.45301,0.00747765,0.000100194,-1.44442e-07,5.75343e-11,-42854.4,17.8733], Tmin=(100,'K'), Tmax=(951.494,'K')), NASAPolynomial(coeffs=[19.4853,0.00335783,3.04212e-07,8.38653e-11,-2.35668e-14,-49150.4,-79.5032], Tmin=(951.494,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-356.976,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH2]C([O])O[C]=[C][O](11958)',
    structure = SMILES('[CH2]C([O])O[C]=[C][O]'),
    E0 = (463.715,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1670,1700,300,440,458.517,459.376,459.408,459.479,460.037],'cm^-1')),
        HinderedRotor(inertia=(0.00079747,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111892,'amu*angstrom^2'), symmetry=1, barrier=(16.7042,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.111535,'amu*angstrom^2'), symmetry=1, barrier=(16.701,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06909,0.0649823,-8.26085e-05,5.19754e-08,-1.27579e-11,55877.2,31.3901], Tmin=(100,'K'), Tmax=(1001.32,'K')), NASAPolynomial(coeffs=[13.8087,0.0140917,-6.37423e-06,1.22015e-09,-8.59702e-14,53325.9,-30.088], Tmin=(1001.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(463.715,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(CJCO) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C](O)O[C]=[C][O](11959)',
    structure = SMILES('[CH2][C](O)O[C]=[C][O]'),
    E0 = (443.256,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,3000,3100,440,815,1455,1000,1670,1700,300,440,423.312,423.498,423.498,423.855],'cm^-1')),
        HinderedRotor(inertia=(0.000940933,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127421,'amu*angstrom^2'), symmetry=1, barrier=(16.2317,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.12763,'amu*angstrom^2'), symmetry=1, barrier=(16.2321,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.127453,'amu*angstrom^2'), symmetry=1, barrier=(16.2309,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.577534,0.0705191,-8.84148e-05,4.97747e-08,-9.97164e-12,53438.9,33.6687], Tmin=(100,'K'), Tmax=(945.622,'K')), NASAPolynomial(coeffs=[18.6497,0.00476068,-1.05769e-06,1.3934e-10,-8.96201e-15,49543.2,-55.0351], Tmin=(945.622,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.256,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cs_P) + radical(CJCO) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'C=C1O[CH][C]OO1(11852)',
    structure = SMILES('[CH2][C]1OC=[C]OO1'),
    E0 = (436.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61122,0.0352742,2.10131e-05,-7.07682e-08,3.54569e-11,52661.2,25.6986], Tmin=(100,'K'), Tmax=(899.487,'K')), NASAPolynomial(coeffs=[20.357,-0.00149457,4.62968e-06,-1.03718e-09,7.02141e-14,47404,-73.2306], Tmin=(899.487,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(436.999,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(124trioxene) + radical(Cs_P) + radical(CJCOOH) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C1([O])O[C]=CO1(11889)',
    structure = SMILES('[CH2]C1([O])O[C]=CO1'),
    E0 = (167.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.72824,0.0265067,5.61821e-05,-1.12165e-07,5.09469e-11,20197,22.3721], Tmin=(100,'K'), Tmax=(910.421,'K')), NASAPolynomial(coeffs=[22.9068,-0.00489899,6.3615e-06,-1.312e-09,8.43596e-14,13786,-91.8453], Tmin=(910.421,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(167.065,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(CCOJ) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = 'C=C([O])OC#C[O](11393)',
    structure = SMILES('[CH2]C(=O)O[C]=C=O'),
    E0 = (73.9475,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.60909,0.0574039,-8.16399e-05,5.95513e-08,-1.62335e-11,8975.46,25.1563], Tmin=(100,'K'), Tmax=(661.243,'K')), NASAPolynomial(coeffs=[9.06258,0.0189315,-9.37356e-06,1.82195e-09,-1.27578e-13,7845.12,-8.81308], Tmin=(661.243,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(73.9475,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + missing(O2d-Cdd) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CJCO) + radical(C=CJO)"""),
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
    label = 'C=C([O])O[C]=C[O](11159)',
    structure = SMILES('[CH2]C(=O)O[C]=C[O]'),
    E0 = (33.9343,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3995.21,'J/mol'), sigma=(6.24473,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=624.04 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09397,0.0566591,-5.34219e-05,2.00133e-08,-1.37476e-12,4192.42,27.5934], Tmin=(100,'K'), Tmax=(996.512,'K')), NASAPolynomial(coeffs=[16.3522,0.00940787,-3.36379e-06,6.18156e-10,-4.47478e-14,456.509,-49.4523], Tmin=(996.512,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(33.9343,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=O)O[C]=[C]O(11891)',
    structure = SMILES('[CH2]C(=O)O[C]=[C]O'),
    E0 = (132.216,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3000,3100,440,815,1455,1000,1670,1700,300,440,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.936331,0.0648634,-8.18558e-05,5.01514e-08,-1.1803e-11,16014.4,29.8782], Tmin=(100,'K'), Tmax=(1052.58,'K')), NASAPolynomial(coeffs=[15.3725,0.0100036,-3.67732e-06,6.36306e-10,-4.26269e-14,12975.4,-40.508], Tmin=(1052.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.216,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CJCO) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'CC(=O)O[C]=[C][O](11892)',
    structure = SMILES('CC(=O)O[C]=[C][O]'),
    E0 = (62.0894,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,1670,1700,300,440,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6938,0.0503819,-5.03696e-05,2.60246e-08,-5.4223e-12,7550.97,27.245], Tmin=(100,'K'), Tmax=(1151.49,'K')), NASAPolynomial(coeffs=[11.0338,0.0179369,-8.10472e-06,1.55483e-09,-1.09663e-13,5399.99,-19.1326], Tmin=(1151.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.0894,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
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
    label = 'C=C1OC([C]=O)O1(11952)',
    structure = SMILES('C=C1OC([C]=O)O1'),
    E0 = (-122.511,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65129,0.0304196,3.95948e-05,-8.67e-08,3.85242e-11,-14630.7,18.7972], Tmin=(100,'K'), Tmax=(947.918,'K')), NASAPolynomial(coeffs=[21.2346,0.00125838,1.11911e-06,-1.26701e-10,-4.10545e-15,-20745.9,-87.306], Tmin=(947.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-122.511,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(CsCJ=O)"""),
)

species(
    label = '[O][C]1CC([C]=O)O1(11960)',
    structure = SMILES('[O][C]1CC([C]=O)O1'),
    E0 = (177.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27788,0.0292294,1.09399e-05,-4.03803e-08,2.00615e-11,21445.9,25.57], Tmin=(100,'K'), Tmax=(886.303,'K')), NASAPolynomial(coeffs=[11.1451,0.0131924,-2.50626e-06,2.63194e-10,-1.42867e-14,18932.2,-21.4527], Tmin=(886.303,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.728,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsOsH) + group(Cds-OdCsH) + ring(Oxetane) + radical(CCOJ) + radical(Cs_P) + radical(CCCJ=O)"""),
)

species(
    label = '[O][C]1[CH]OC(=O)C1(11842)',
    structure = SMILES('[O][C]1CC([O])=CO1'),
    E0 = (25.443,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94221,0.0296169,2.80331e-05,-7.21042e-08,3.47037e-11,3148.94,19.9941], Tmin=(100,'K'), Tmax=(895.37,'K')), NASAPolynomial(coeffs=[17.5518,0.00218763,3.11111e-06,-7.77351e-10,5.38905e-14,-1342.13,-63.0579], Tmin=(895.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(25.443,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ) + radical(CCOJ) + radical(Cs_P)"""),
)

species(
    label = '[CH]=C(O)O[CH][C]=O(11961)',
    structure = SMILES('[CH]=C(O)O[CH][C]=O'),
    E0 = (186.462,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.893531,'amu*angstrom^2'), symmetry=1, barrier=(20.544,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.895212,'amu*angstrom^2'), symmetry=1, barrier=(20.5827,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.893917,'amu*angstrom^2'), symmetry=1, barrier=(20.5529,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.894041,'amu*angstrom^2'), symmetry=1, barrier=(20.5558,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0717704,0.080542,-0.000104065,5.68813e-08,-1.02938e-11,22573.1,25.3356], Tmin=(100,'K'), Tmax=(918.666,'K')), NASAPolynomial(coeffs=[22.3685,-0.000267692,1.3093e-06,-3.05542e-10,2.12546e-14,17789.7,-84.0794], Tmin=(918.666,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(186.462,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CCsJOC(O)) + radical(Cds_P) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=C([O])OC[C]=O(11962)',
    structure = SMILES('[CH]C(=O)OC[C]=O'),
    E0 = (72.8096,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1855,455,950,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99975,0.0434102,-3.30976e-05,1.23988e-08,-1.89963e-12,8829.07,24.9077], Tmin=(100,'K'), Tmax=(1493.18,'K')), NASAPolynomial(coeffs=[10.8139,0.0197984,-9.3778e-06,1.80853e-09,-1.26513e-13,6196.86,-21.1491], Tmin=(1493.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(72.8096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CsCJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=C([O])OC=C[O](11895)',
    structure = SMILES('[CH]C(=O)OC=C[O]'),
    E0 = (30.8162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,300,800,800,800,800,800,1600,1600,1600,1600,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.01004,0.0524625,-2.52575e-05,-1.95699e-08,1.5163e-11,3826.19,24.4432], Tmin=(100,'K'), Tmax=(940.72,'K')), NASAPolynomial(coeffs=[20.382,0.0033097,1.5119e-07,-4.02378e-11,-2.49909e-15,-1288.34,-75.6436], Tmin=(940.72,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(30.8162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C1OC=C([O])O1(11856)',
    structure = SMILES('C=C1O[CH]C(=O)O1'),
    E0 = (-193.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63788,0.0324863,3.12521e-05,-7.7117e-08,3.5191e-11,-23181.5,17.1761], Tmin=(100,'K'), Tmax=(941.477,'K')), NASAPolynomial(coeffs=[20.3633,0.00218402,1.05512e-06,-1.65076e-10,1.2665e-15,-28890.3,-83.6272], Tmin=(941.477,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-193.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(CCsJOC(O))"""),
)

species(
    label = 'C[C]([O])O[C]=[C][O](11963)',
    structure = SMILES('C[C]([O])O[C]=[C][O]'),
    E0 = (457.372,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,1670,1700,300,440,378.818,378.828,378.843,378.854,378.88],'cm^-1')),
        HinderedRotor(inertia=(0.00117464,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0011747,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.197605,'amu*angstrom^2'), symmetry=1, barrier=(20.1222,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.28911,0.0601645,-7.20271e-05,4.33654e-08,-1.02883e-11,55106.4,30.6547], Tmin=(100,'K'), Tmax=(1029.96,'K')), NASAPolynomial(coeffs=[12.7553,0.0156327,-7.17081e-06,1.38456e-09,-9.80768e-14,52744.5,-25.0014], Tmin=(1029.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(457.372,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(Cs_P) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = 'O=C1CO[C][CH]O1(11814)',
    structure = SMILES('O=C1CO[C][CH]O1'),
    E0 = (131.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.62796,0.0356439,1.59871e-05,-4.95481e-08,2.10045e-11,15965.3,17.2488], Tmin=(100,'K'), Tmax=(1035.26,'K')), NASAPolynomial(coeffs=[17.7025,0.0125107,-6.96591e-06,1.59762e-09,-1.28005e-13,10548.5,-70.9457], Tmin=(1035.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(131.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(CCsJOC(O)) + radical(CH2_triplet)"""),
)

species(
    label = 'C=C(O)O[C]=[C][O](11896)',
    structure = SMILES('C=C(O)O[C][C]=O'),
    E0 = (206.427,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.927126,'amu*angstrom^2'), symmetry=1, barrier=(21.3164,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92041,'amu*angstrom^2'), symmetry=1, barrier=(21.162,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.926978,'amu*angstrom^2'), symmetry=1, barrier=(21.313,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.92169,'amu*angstrom^2'), symmetry=1, barrier=(21.1915,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0129669,0.0769611,-9.75534e-05,5.54841e-08,-1.17262e-11,24980.2,24.4005], Tmin=(100,'K'), Tmax=(1262.48,'K')), NASAPolynomial(coeffs=[22.8009,-0.000834011,1.52496e-06,-3.45413e-10,2.4362e-14,19672.2,-89.0838], Tmin=(1262.48,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(206.427,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=C([O])OC=[C]O(11964)',
    structure = SMILES('[CH]C(=O)OC=[C]O'),
    E0 = (129.098,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,1685,370,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.284391,0.067315,-7.67237e-05,4.01604e-08,-7.78269e-12,15673,28.7691], Tmin=(100,'K'), Tmax=(1437.26,'K')), NASAPolynomial(coeffs=[20.468,0.00201337,9.56045e-07,-2.90326e-10,2.20937e-14,10814.1,-72.6465], Tmin=(1437.26,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(129.098,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(C=CJO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O][C]=CCC([O])=O(11965)',
    structure = SMILES('[O]C(=O)C[CH][C]=O'),
    E0 = (9.97861,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.92138,0.0553077,-9.47389e-05,9.75423e-08,-3.8326e-11,1265.67,25.9159], Tmin=(100,'K'), Tmax=(830.537,'K')), NASAPolynomial(coeffs=[1.42506,0.035284,-1.80936e-05,3.52571e-09,-2.45055e-13,2121.16,32.8722], Tmin=(830.537,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.97861,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cs-(Cds-O2d)CsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cds-OdCsOs) + group(Cds-OdCsH) + radical(CCOJ) + radical(CCJCC=O) + radical(CCCJ=O)"""),
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
    label = '[CH2][C]([O])O[C]=C[O](11966)',
    structure = SMILES('[CH2][C]([O])O[C]=C[O]'),
    E0 = (429.217,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,1685,370,499.651,499.656,499.666,499.669,499.672],'cm^-1')),
        HinderedRotor(inertia=(0.102148,'amu*angstrom^2'), symmetry=1, barrier=(18.098,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102154,'amu*angstrom^2'), symmetry=1, barrier=(18.0983,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.102155,'amu*angstrom^2'), symmetry=1, barrier=(18.0983,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.678041,0.0666163,-7.58503e-05,3.85089e-08,-6.77267e-12,51748.2,31.0401], Tmin=(100,'K'), Tmax=(993.937,'K')), NASAPolynomial(coeffs=[18.3711,0.00658973,-2.12969e-06,3.76278e-10,-2.71774e-14,47678.9,-56.9887], Tmin=(993.937,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.217,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + radical(CCOJ) + radical(C=COJ) + radical(Cs_P) + radical(CJCO) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]([O])[O](2822)',
    structure = SMILES('[CH2][C]([O])[O]'),
    E0 = (411.444,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([360,370,350,3000,3100,440,815,1455,1000,1740.04,1740.26],'cm^-1')),
        HinderedRotor(inertia=(0.0034914,'amu*angstrom^2'), symmetry=1, barrier=(7.50222,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (58.0361,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.01682,0.0315837,-7.13575e-05,8.2866e-08,-3.37674e-11,49511,16.8538], Tmin=(100,'K'), Tmax=(856.99,'K')), NASAPolynomial(coeffs=[-0.765111,0.0236783,-1.2787e-05,2.50404e-09,-1.72816e-13,51097.8,39.9915], Tmin=(856.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.444,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + radical(CCOJ) + radical(CCOJ) + radical(Cs_P) + radical(CJCO)"""),
)

species(
    label = '[CH2]C([O])O[C]=C=O(11967)',
    structure = SMILES('[CH2]C([O])O[C]=C=O'),
    E0 = (263.984,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,362.751,363.065,363.85,364.781],'cm^-1')),
        HinderedRotor(inertia=(0.054873,'amu*angstrom^2'), symmetry=1, barrier=(5.15403,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.170257,'amu*angstrom^2'), symmetry=1, barrier=(16.0323,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0163037,'amu*angstrom^2'), symmetry=1, barrier=(16.0301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11798,0.0672389,-9.67591e-05,7.19382e-08,-2.11621e-11,31850.2,27.4827], Tmin=(100,'K'), Tmax=(834.314,'K')), NASAPolynomial(coeffs=[11.3138,0.0183569,-8.87556e-06,1.71448e-09,-1.19899e-13,30148.9,-19.8594], Tmin=(834.314,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + missing(O2d-Cdd) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-(Cdd-O2d)OsH) + missing(Cdd-CdO2d) + radical(CCOJ) + radical(CJCO) + radical(C=CJO)"""),
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
    E0 = (30.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (416.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (510.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (478.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (485.471,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (38.3953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (181.486,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (272.281,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (88.6392,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (176.102,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (117.419,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (116.751,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (365.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (30.111,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (37.2238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (527.115,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (468.229,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (436.999,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (167.065,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (300.849,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (148.053,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (203.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (238.072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (282.561,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (164.212,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (855.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (38.3953,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (177.728,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (104.365,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (378.591,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (117.118,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (63.8564,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (37.2238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (482.345,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (131.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (295.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (280.618,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (343.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (573.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (451.518,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (576.634,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (405.961,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['OCHCO(3676)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C][O](6861)', 'C=C([O])[O](2321)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(30799.3,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['CH2(T)(20)', 'O=[C][CH]O[C]=O(11954)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH]=C([O])O[CH][C]=O(11365)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['H(3)', '[CH2]C(=O)O[C]=[C][O](11890)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['O=[C]C1CC(=O)O1(11151)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['O=[C][CH]O[C]1CO1(11955)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(3.95361e+10,'s^-1'), n=0.549916, Ea=(151.375,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R3_CO;carbonyl_intra_Nd;radadd_intra_cs2H]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Endocyclic
Ea raised from 150.8 to 151.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['[CH2][C]1OC([C]=O)O1(11956)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cs]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['C=C1O[CH][C]([O])O1(11882)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(6.44463e+09,'s^-1'), n=0.470283, Ea=(58.5282,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra] for rate rule [R5_linear;carbonyl_intra;radadd_intra_CO]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['[CH2]C1([O])OC1[C]=O(11957)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(2.63856e+09,'s^-1'), n=0.755479, Ea=(145.991,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S;multiplebond_intra;radadd_intra_cs] for rate rule [R4_S_CO;carbonylbond_intra;radadd_intra_cs]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 145.3 to 146.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['[CH2]C1([O])OC=C1[O](11915)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(4.70082e+08,'s^-1'), n=0.843229, Ea=(87.308,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;multiplebond_intra;radadd_intra] for rate rule [R5;carbonylbond_intra;radadd_intra_CO]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['OCHCO(3676)', '[CH2][C]=O(1376)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.41198,'m^3/(mol*s)'), n=2.03228, Ea=(31.3669,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;CJ] + [Od_R;YJ] for rate rule [Od_R;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[O][C]=C[O](9592)', '[CH2][C]=O(1376)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['C=C([O])OC=C=O(11168)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['[O]C1=COC(=O)C1(11818)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSSD;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH2]C([O])O[C]=[C][O](11958)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH2][C](O)O[C]=[C][O](11959)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction18',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['C=C1O[CH][C]OO1(11852)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(3.41956e+10,'s^-1'), n=0.267163, Ea=(406.888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra] for rate rule [R6_linear;carbonyl_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 404.3 to 406.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['[CH2]C1([O])O[C]=CO1(11889)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3.16887e+10,'s^-1'), n=0.492453, Ea=(136.954,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;multiplebond_intra;radadd_intra_O] + [R6;multiplebond_intra;radadd_intra] for rate rule [R6;carbonylbond_intra;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 133.4 to 137.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', 'C=C([O])OC#C[O](11393)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct_Ct;HJ] for rate rule [Ct-O_Ct;HJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][C]=C[O](9592)', 'ketene(1375)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(0.000692465,'m^3/(mol*s)'), n=2.43792, Ea=(4.66472,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Ck_Cds-HH;OJ_sec] for rate rule [Ck_Cds-HH;O_rad/OneDe]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['HCCO(2227)', 'C=C([O])[O](2321)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(520000,'m^3/(mol*s)'), n=0, Ea=(46.5868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;O_rad/OneDe] + [Ct_Ct;OJ_sec] for rate rule [Ct_Ct;O_rad/OneDe]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=C([O])O[C]=C[O](11159)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.43441e+10,'s^-1'), n=1.037, Ea=(204.137,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleNd;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH2]C(=O)O[C]=[C]O(11891)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(4.96975e+09,'s^-1'), n=0.933333, Ea=(150.345,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleNd;XH_out] for rate rule [R3HJ;Cd_rad_out_singleNd;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['CC(=O)O[C]=[C][O](11892)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.35985e+07,'s^-1'), n=1.33667, Ea=(102.123,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SSS;Y_rad_out;Cs_H_out_2H] for rate rule [R4H_SSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['O(4)', 'C=[C]O[CH][C]=O(2657)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['C=C1OC([C]=O)O1(11952)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['[O][C]1CC([C]=O)O1(11960)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(8.33681e+07,'s^-1'), n=1.1093, Ea=(147.617,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 143.9 to 147.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['[O][C]1[CH]OC(=O)C1(11842)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(5.50361e+10,'s^-1'), n=0.23641, Ea=(74.2541,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;doublebond_intra;radadd_intra] for rate rule [R5_linear;doublebond_intra;radadd_intra_CO]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=C(O)O[CH][C]=O(11961)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(6718.85,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out_singleH;XH_out] for rate rule [R3H_DS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=C([O])OC[C]=O(11962)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]=C([O])OC=C[O](11895)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_3;Cd_rad_out_singleH;CO_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['C=C1OC=C([O])O1(11856)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSSD;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C[C]([O])O[C]=[C][O](11963)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['O=C1CO[C][CH]O1(11814)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.49515e+10,'s^-1'), n=0.243684, Ea=(101.803,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;doublebond_intra;radadd_intra] for rate rule [R6_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 99.3 to 101.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=C(O)O[C]=[C][O](11896)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(26.449,'s^-1'), n=2.8625, Ea=(89.0146,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4H_SSS;Y_rad_out;O_H_out] + [R4H_SS(Cd)S;Y_rad_out;XH_out] for rate rule [R4H_SS(Cd)S;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=C([O])OC=[C]O(11964)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_4;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction38',
    reactants = ['C=C([O])O[CH][C]=O(11146)'],
    products = ['[O][C]=CCC([O])=O(11965)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond_CH2;R2_doublebond;R_O_C]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction39',
    reactants = ['OCHCO(3676)', '[CH2][C][O](2821)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(15399.6,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH2][C]([O])O[C]=C[O](11966)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(2.24409e+10,'s^-1'), n=0.34095, Ea=(22.3009,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_NDe;XH_Rrad] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad_NDe;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['HCCO(2227)', '[CH2][C]([O])[O](2822)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(212954,'m^3/(mol*s)'), n=0.348287, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.0108230153501, var=2.70964383578, Tref=1000.0, N=19, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R
    Total Standard Deviation in ln(k): 3.32718707999
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C([O])O[C]=C=O(11967)'],
    products = ['C=C([O])O[CH][C]=O(11146)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.4115e+09,'s^-1'), n=1.00333, Ea=(141.977,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Cd_rad_out_double;XH_out] for rate rule [R3H_SS_O;Cd_rad_out_double;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2814',
    isomers = [
        'C=C([O])O[CH][C]=O(11146)',
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
    label = 'PDepNetwork #2814',
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


species(
    label = 'C=[C]OC([O])=C[O](11160)',
    structure = SMILES('C=[C]OC([O])=C[O]'),
    E0 = (126.554,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,1685,370,215.876,216.738,218.675,220.075,220.616],'cm^-1')),
        HinderedRotor(inertia=(0.00354602,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.8711,'amu*angstrom^2'), symmetry=1, barrier=(29.5815,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06267,0.0596292,-6.66166e-05,3.64055e-08,-7.68211e-12,15330.9,27.1783], Tmin=(100,'K'), Tmax=(1168.97,'K')), NASAPolynomial(coeffs=[15.3167,0.0108545,-4.02971e-06,7.11915e-10,-4.85487e-14,11998.4,-43.8144], Tmin=(1168.97,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(126.554,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ) + radical(C=CJO)"""),
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
    label = '[O]C=C([O])[O](11106)',
    structure = SMILES('[O]C=C([O])[O]'),
    E0 = (-108.185,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (73.0275,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.23477,0.0289406,-7.48904e-06,-2.38398e-08,1.47722e-11,-12938.9,15.4323], Tmin=(100,'K'), Tmax=(918.148,'K')), NASAPolynomial(coeffs=[15.6525,-0.00349641,2.99678e-06,-5.89032e-10,3.7326e-14,-16499.4,-54.1269], Tmin=(918.148,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-108.185,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = 'C=[C]O[C]=C[O](9888)',
    structure = SMILES('C=[C]O[C]=C[O]'),
    E0 = (381.972,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,181.707,184.449,185.116],'cm^-1')),
        HinderedRotor(inertia=(0.00506491,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.25017,'amu*angstrom^2'), symmetry=1, barrier=(29.5113,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42265,0.0464145,-4.64049e-05,2.25878e-08,-4.15873e-12,46042.2,27.8232], Tmin=(100,'K'), Tmax=(1488.77,'K')), NASAPolynomial(coeffs=[14.2582,0.00685404,-1.43369e-06,1.60506e-10,-8.21627e-15,42782.7,-37.3203], Tmin=(1488.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.972,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJO) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(=O)O[C]=C(2687)',
    structure = SMILES('[CH]C(=O)O[C]=C'),
    E0 = (337.89,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,200,800,900,1000,1100,1200,1300,1400,1500,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3356.73,'J/mol'), sigma=(5.48768,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=524.31 K, Pc=46.09 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.08736,0.0425357,-4.06747e-05,2.02687e-08,-4.12913e-12,40707.3,22.4383], Tmin=(100,'K'), Tmax=(1164.53,'K')), NASAPolynomial(coeffs=[9.47998,0.017143,-7.96692e-06,1.5442e-09,-1.09366e-13,38985.5,-14.3529], Tmin=(1164.53,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(337.89,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cs-(Cds-O2d)HHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C1OC(=C[O])O1(11854)',
    structure = SMILES('C=C1O[C](C=O)O1'),
    E0 = (-71.936,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.63733,0.0307973,3.86403e-05,-8.58867e-08,3.82949e-11,-8547.47,18.8565], Tmin=(100,'K'), Tmax=(947.423,'K')), NASAPolynomial(coeffs=[21.2895,0.00113799,1.19275e-06,-1.43184e-10,-2.83114e-15,-14663.9,-87.5196], Tmin=(947.423,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-71.936,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(Cs-CsOsOsH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(Cs_P)"""),
)

species(
    label = 'C=[C]OC1=COO1(11855)',
    structure = SMILES('C=[C]OC1=COO1'),
    E0 = (362.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86555,0.045045,-3.76155e-05,1.57152e-08,-2.65837e-12,43734.3,24.6867], Tmin=(100,'K'), Tmax=(1387.81,'K')), NASAPolynomial(coeffs=[11.3179,0.0178008,-8.16865e-06,1.56968e-09,-1.10177e-13,41110.7,-24.0134], Tmin=(1387.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.978,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=CJO)"""),
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
    label = 'C=[C]OC(O)=C=O(11857)',
    structure = SMILES('C=[C]OC(O)=C=O'),
    E0 = (-7.3072,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.639325,0.0547771,-1.30008e-05,-5.09172e-08,3.19368e-11,-739.785,25.793], Tmin=(100,'K'), Tmax=(898.628,'K')), NASAPolynomial(coeffs=[27.2218,-0.0108777,8.67354e-06,-1.77343e-09,1.19601e-13,-7643.96,-111.443], Tmin=(898.628,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-7.3072,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cds-(Cdd-O2d)OsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(C=CJO)"""),
)

species(
    label = 'C=COC([O])=C=O(11858)',
    structure = SMILES('C=COC(=O)[C]=O'),
    E0 = (-229.758,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04751,0.0580629,-5.96775e-05,2.92857e-08,-5.55213e-12,-27521.4,22.9532], Tmin=(100,'K'), Tmax=(1294.25,'K')), NASAPolynomial(coeffs=[16.4557,0.0104427,-4.48726e-06,8.57364e-10,-6.0866e-14,-31509.8,-55.3565], Tmin=(1294.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-229.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(OC=OCJ=O)"""),
)

species(
    label = 'C=[C]O[C]1OC1[O](11859)',
    structure = SMILES('C=[C]O[C]1OC1[O]'),
    E0 = (274.29,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.883449,0.0571793,-6.13515e-05,3.29687e-08,-6.61349e-12,33111.5,27.9009], Tmin=(100,'K'), Tmax=(1428.94,'K')), NASAPolynomial(coeffs=[14.5298,0.00918245,-6.8362e-07,-1.33617e-10,1.736e-14,30211.7,-39.3058], Tmin=(1428.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(274.29,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(Cs_P) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]OC1([O])[CH]O1(11860)',
    structure = SMILES('C=[C]OC1([O])[CH]O1'),
    E0 = (272.619,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.773915,0.0621234,-7.48469e-05,4.47676e-08,-9.85106e-12,32912.4,26.4967], Tmin=(100,'K'), Tmax=(1330.95,'K')), NASAPolynomial(coeffs=[14.2719,0.00836835,5.99671e-07,-4.68472e-10,4.43634e-14,30487.5,-38.0938], Tmin=(1330.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.619,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsOs) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(CCsJO) + radical(C=CJO)"""),
)

species(
    label = 'C=C1O[C]([O])C1[O](11861)',
    structure = SMILES('[CH2][C]1OC(=O)C1[O]'),
    E0 = (151.689,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99953,0.0338261,4.56936e-06,-3.88145e-08,2.0531e-11,18325.6,26.7375], Tmin=(100,'K'), Tmax=(895.652,'K')), NASAPolynomial(coeffs=[13.749,0.00922026,-8.93269e-07,-9.1446e-12,2.71439e-15,15103.2,-34.8922], Tmin=(895.652,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(151.689,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cds-OdCsOs) + ring(Beta-Propiolactone) + radical(C=OCOJ) + radical(C2CsJOC(O)) + radical(CJC(C)OC)"""),
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
    label = 'C=[C]OC([O])=C=O(11863)',
    structure = SMILES('C=[C]OC(=O)[C]=O'),
    E0 = (9.98623,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1855,455,950,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.6697,0.0548591,-6.8831e-05,4.47116e-08,-1.17069e-11,1281.95,23.857], Tmin=(100,'K'), Tmax=(925.061,'K')), NASAPolynomial(coeffs=[10.1421,0.0182242,-9.42728e-06,1.90107e-09,-1.37251e-13,-285.551,-16.3574], Tmin=(925.061,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(9.98623,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(195.39,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(C=CJO) + radical(OC=OCJ=O)"""),
)

species(
    label = '[CH]=[C]OC(=O)C=O(11368)',
    structure = SMILES('C#COC([O])=C[O]'),
    E0 = (96.1362,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.44035,'amu*angstrom^2'), symmetry=1, barrier=(33.1165,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.44619,'amu*angstrom^2'), symmetry=1, barrier=(33.2508,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.845619,0.0636126,-7.62931e-05,4.2085e-08,-8.62204e-12,11681.1,22.2844], Tmin=(100,'K'), Tmax=(1036.01,'K')), NASAPolynomial(coeffs=[17.5525,0.00567345,-1.91106e-06,3.37473e-10,-2.39889e-14,7867.09,-60.6081], Tmin=(1036.01,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(96.1362,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Ct-CtOs) + group(Ct-CtH) + radical(C=COJ) + radical(C=COJ)"""),
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
    label = 'C=[C]OC([O])=[C][O](11864)',
    structure = SMILES('[CH2][C]OC(=O)[C]=O'),
    E0 = (329.401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34967,0.0627163,-8.64165e-05,6.02345e-08,-1.67178e-11,39709.4,24.3374], Tmin=(100,'K'), Tmax=(879.004,'K')), NASAPolynomial(coeffs=[11.2455,0.0176842,-9.57016e-06,1.95153e-09,-1.41316e-13,37969.7,-22.128], Tmin=(879.004,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(329.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(191.233,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(CJCO) + radical(CH2_triplet) + radical(OC=OCJ=O)"""),
)

species(
    label = '[CH]=[C]OC([O])=C[O](11375)',
    structure = SMILES('[CH]=[C]OC([O])=C[O]'),
    E0 = (373.65,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,324.026,324.089,324.282,325.447],'cm^-1')),
        HinderedRotor(inertia=(0.276329,'amu*angstrom^2'), symmetry=1, barrier=(20.6065,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.273654,'amu*angstrom^2'), symmetry=1, barrier=(20.5895,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (98.0569,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.08519,0.0620141,-7.93085e-05,4.87101e-08,-1.14763e-11,45046.5,27.5013], Tmin=(100,'K'), Tmax=(1051.4,'K')), NASAPolynomial(coeffs=[15.0579,0.00885665,-3.47171e-06,6.24649e-10,-4.28711e-14,42108.3,-40.6094], Tmin=(1051.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(373.65,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OC([O])=[C]O(11865)',
    structure = SMILES('C=[C]OC(=O)[C]O'),
    E0 = (205.41,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2950,3100,1380,975,1025,1650,1685,370,180,180,180,282.73,852.415,1600,1828.57,2971.43,3200],'cm^-1')),
        HinderedRotor(inertia=(0.149915,'amu*angstrom^2'), symmetry=1, barrier=(3.44683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149915,'amu*angstrom^2'), symmetry=1, barrier=(3.44683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149915,'amu*angstrom^2'), symmetry=1, barrier=(3.44683,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.149915,'amu*angstrom^2'), symmetry=1, barrier=(3.44683,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14524,0.0646251,-8.24162e-05,5.19861e-08,-1.28757e-11,24806.4,25.5202], Tmin=(100,'K'), Tmax=(988.4,'K')), NASAPolynomial(coeffs=[13.3321,0.0153073,-7.57384e-06,1.5074e-09,-1.08369e-13,22397.3,-33.1326], Tmin=(988.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.41,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=CJO) + radical(CH2_triplet)"""),
)

species(
    label = '[CH]=COC([O])=C[O](11300)',
    structure = SMILES('[CH]=COC([O])=C[O]'),
    E0 = (133.906,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3120,650,792.5,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19631,'amu*angstrom^2'), symmetry=1, barrier=(27.5056,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.20086,'amu*angstrom^2'), symmetry=1, barrier=(27.6101,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.318789,0.0667627,-7.49371e-05,3.88163e-08,-7.46793e-12,16249.9,27.1265], Tmin=(100,'K'), Tmax=(1438.77,'K')), NASAPolynomial(coeffs=[20.3199,0.00276934,5.23723e-07,-2.00802e-10,1.56835e-14,11362.6,-73.6267], Tmin=(1438.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(133.906,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OC(O)=[C][O](11866)',
    structure = SMILES('C=[C]O[C](O)[C]=O'),
    E0 = (132.914,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09645,0.0590704,-6.654e-05,3.60213e-08,-7.51018e-12,16094.6,29.6719], Tmin=(100,'K'), Tmax=(1182.8,'K')), NASAPolynomial(coeffs=[15.6941,0.00970357,-3.93367e-06,7.33954e-10,-5.17019e-14,12641.4,-43.2041], Tmin=(1182.8,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(132.914,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = 'C=COC([O])=[C][O](11867)',
    structure = SMILES('[CH2][CH]OC(=O)[C]=O'),
    E0 = (62.3399,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14643,0.0667703,-9.03975e-05,6.27582e-08,-1.7341e-11,7597.05,25.8066], Tmin=(100,'K'), Tmax=(883.835,'K')), NASAPolynomial(coeffs=[11.6435,0.0192635,-9.77125e-06,1.94271e-09,-1.38784e-13,5741.53,-23.5394], Tmin=(883.835,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(62.3399,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(CCsJOC(O)) + radical(CJCO) + radical(OC=OCJ=O)"""),
)

species(
    label = '[CH]=[C]OC(O)=C[O](11868)',
    structure = SMILES('[CH]=[C]O[C](O)C=O'),
    E0 = (225.34,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,360,370,350,2782.5,750,1395,475,1775,1000,1685,370,3120,650,792.5,1650,378.394,378.586],'cm^-1')),
        HinderedRotor(inertia=(0.00117802,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164452,'amu*angstrom^2'), symmetry=1, barrier=(16.7186,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164453,'amu*angstrom^2'), symmetry=1, barrier=(16.7174,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.164546,'amu*angstrom^2'), symmetry=1, barrier=(16.7193,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.990647,0.0620192,-7.29826e-05,4.11131e-08,-8.8957e-12,27214.1,29.7757], Tmin=(100,'K'), Tmax=(1142.22,'K')), NASAPolynomial(coeffs=[16.0494,0.00928398,-3.72878e-06,6.92434e-10,-4.87124e-14,23774,-44.8765], Tmin=(1142.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(225.34,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cs_P) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]OC([O])=CO(11869)',
    structure = SMILES('[CH]=[C]OC(=O)[CH]O'),
    E0 = (130.003,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3025,407.5,1350,352.5,1685,370,3120,650,792.5,1650,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18107,0.0661279,-9.37091e-05,7.01185e-08,-2.09528e-11,15733.5,26.8197], Tmin=(100,'K'), Tmax=(818.86,'K')), NASAPolynomial(coeffs=[10.5433,0.0203978,-9.94516e-06,1.92716e-09,-1.35052e-13,14200.1,-16.4772], Tmin=(818.86,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.003,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = 'C=[C]OC(=O)C=O(11153)',
    structure = SMILES('C=[C]OC(=O)C=O'),
    E0 = (-148.541,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77901,0.0538041,-6.07696e-05,3.79431e-08,-9.98495e-12,-17789.5,22.6938], Tmin=(100,'K'), Tmax=(901.153,'K')), NASAPolynomial(coeffs=[8.23716,0.0251373,-1.30517e-05,2.64091e-09,-1.91124e-13,-18953.5,-7.79074], Tmin=(901.153,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-148.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-CdsOsH) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + radical(C=CJO)"""),
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
    label = 'C=[C]O[C]=O(2686)',
    structure = SMILES('C=[C]O[C]=O'),
    E0 = (141.993,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2950,3100,1380,975,1025,1650,1685,370,1855,455,950,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.24763,'amu*angstrom^2'), symmetry=1, barrier=(28.6854,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.24862,'amu*angstrom^2'), symmetry=1, barrier=(28.7082,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.4468,0.0347823,-3.82777e-05,2.17539e-08,-4.964e-12,17133.3,16.4613], Tmin=(100,'K'), Tmax=(1058.68,'K')), NASAPolynomial(coeffs=[8.58938,0.011574,-5.3951e-06,1.04734e-09,-7.43263e-14,15832.7,-13.5235], Tmin=(1058.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(141.993,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(149.66,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cds-OdOsH) + radical(C=CJO) + radical((O)CJOC)"""),
)

species(
    label = 'C=C1OC(=O)C1[O](11870)',
    structure = SMILES('C=C1OC(=O)C1[O]'),
    E0 = (-178.257,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.84533,0.0338664,1.48943e-05,-5.68117e-08,2.88248e-11,-21349,23.6873], Tmin=(100,'K'), Tmax=(892.564,'K')), NASAPolynomial(coeffs=[16.6417,0.00460332,1.81368e-06,-5.39806e-10,3.87263e-14,-25466,-54.2817], Tmin=(892.564,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-178.257,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsOs) + group(Cds-OdCsOs) + group(Cds-CdsHH) + ring(4-Methylene-2-oxetanone) + radical(C=OCOJ)"""),
)

species(
    label = 'C=C1O[C]([CH][O])O1(11871)',
    structure = SMILES('[CH2][C]1OC(=C[O])O1'),
    E0 = (153.332,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.953676,0.0408779,3.4262e-05,-1.00064e-07,4.83395e-11,18575.6,20.9966], Tmin=(100,'K'), Tmax=(923.054,'K')), NASAPolynomial(coeffs=[28.5128,-0.0107462,7.97277e-06,-1.49998e-09,9.18857e-14,10599.4,-125.4], Tmin=(923.054,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(153.332,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsOsH) + group(Cs-CsHHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclobutane) + radical(C=COJ) + radical(Cs_P) + radical(CJCO)"""),
)

species(
    label = 'C=[C]O[C]1[CH]OO1(11872)',
    structure = SMILES('C=[C]O[C]1[CH]OO1'),
    E0 = (487.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41313,0.0547693,-5.56647e-05,2.83231e-08,-5.71368e-12,58734.8,23.1244], Tmin=(100,'K'), Tmax=(1199.25,'K')), NASAPolynomial(coeffs=[13.0288,0.016026,-7.20514e-06,1.38423e-09,-9.78777e-14,55948.8,-35.025], Tmin=(1199.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(487.56,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(Cs_P) + radical(CCsJOO) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]OC(=O)C[O](11873)',
    structure = SMILES('[CH]=[C]OC(=O)C[O]'),
    E0 = (235.251,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3120,650,792.5,1650,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42302,0.0604761,-8.19306e-05,5.98388e-08,-1.76069e-11,28383.5,26.9753], Tmin=(100,'K'), Tmax=(828.708,'K')), NASAPolynomial(coeffs=[9.66261,0.0207044,-9.94049e-06,1.92396e-09,-1.35069e-13,27017.9,-11.2277], Tmin=(828.708,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(235.251,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)(Cds-Cd)) + group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CJO) + radical(Cds_P)"""),
)

species(
    label = '[CH2][C]OC([O])[C]=O(11874)',
    structure = SMILES('[CH2][C]OC([O])[C]=O'),
    E0 = (516.369,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3000,3100,440,815,1455,1000,1855,455,950,200,800,1066.67,1333.33,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.904454,0.0747137,-0.000125963,1.09517e-07,-3.72038e-11,62209.9,28.9906], Tmin=(100,'K'), Tmax=(814.858,'K')), NASAPolynomial(coeffs=[9.96099,0.0213247,-1.12419e-05,2.20703e-09,-1.53796e-13,61030.5,-11.0281], Tmin=(814.858,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(516.369,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CJCO) + radical(CH2_triplet) + radical(CsCJ=O)"""),
)

species(
    label = 'O=C[C]1O[C]CO1(11875)',
    structure = SMILES('[O]C=C1O[C]CO1'),
    E0 = (156.122,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14313,0.0345911,5.33619e-05,-1.22053e-07,5.71211e-11,18906.4,18.1574], Tmin=(100,'K'), Tmax=(914.852,'K')), NASAPolynomial(coeffs=[28.9378,-0.0129097,9.87199e-06,-1.9147e-09,1.21714e-13,10723,-130.393], Tmin=(914.852,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(156.122,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(CH2_triplet)"""),
)

species(
    label = 'O=C1[CH]OC[C]O1(11876)',
    structure = SMILES('O=C1[CH]OC[C]O1'),
    E0 = (118.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69723,0.0332545,2.35862e-05,-5.80958e-08,2.42715e-11,14343.5,16.9707], Tmin=(100,'K'), Tmax=(1020.34,'K')), NASAPolynomial(coeffs=[17.7311,0.0120218,-6.39068e-06,1.47084e-09,-1.19129e-13,8904.73,-71.3239], Tmin=(1020.34,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(118.442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-O2d)OsHH) + group(Cds-OdCsOs) + ring(Cyclohexanone) + radical(CCsJOCs) + radical(CH2_triplet)"""),
)

species(
    label = '[O]C1(C=O)C[C]O1(11877)',
    structure = SMILES('[O]C1(C=O)C[C]O1'),
    E0 = (278.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34126,0.0371266,3.34805e-05,-9.9087e-08,5.06465e-11,33648.4,21.1736], Tmin=(100,'K'), Tmax=(872.766,'K')), NASAPolynomial(coeffs=[24.9528,-0.0104797,1.11344e-05,-2.45039e-09,1.73634e-13,27218.6,-102.749], Tmin=(872.766,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(278.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdCsH) + ring(Oxetane) + radical(C=OCOJ) + radical(CH2_triplet)"""),
)

species(
    label = '[O]C1C[C]OC1=O(11878)',
    structure = SMILES('[O]C1C[C]OC1=O'),
    E0 = (152.736,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.52322,0.0117532,7.88912e-05,-1.21252e-07,5.07901e-11,18442.3,22.3342], Tmin=(100,'K'), Tmax=(912.824,'K')), NASAPolynomial(coeffs=[16.418,0.00469552,2.03373e-06,-5.18935e-10,3.13892e-14,13663,-55.7166], Tmin=(912.824,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.736,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(232.805,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsOsHH) + group(Cds-OdCsOs) + ring(butyrolactone) + radical(C=OCOJ) + radical(CH2_triplet)"""),
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
    label = '[CH2][C]O[C]=O(11879)',
    structure = SMILES('[CH2][C]O[C]=O'),
    E0 = (462.641,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,1855,455,950,320.034,320.275,320.281],'cm^-1')),
        HinderedRotor(inertia=(0.00164831,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266197,'amu*angstrom^2'), symmetry=1, barrier=(19.3506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.266722,'amu*angstrom^2'), symmetry=1, barrier=(19.3488,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97738,0.0428924,-5.43378e-05,3.24351e-08,-7.41205e-12,55717.2,18.4812], Tmin=(100,'K'), Tmax=(1083.21,'K')), NASAPolynomial(coeffs=[12.1359,0.00537981,-2.39157e-06,4.64673e-10,-3.343e-14,53516.4,-31.3398], Tmin=(1083.21,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(462.641,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(145.503,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-OdOsH) + radical(CJCO) + radical(CH2_triplet) + radical((O)CJOCC)"""),
)

species(
    label = 'C[C]OC(=O)[C]=O(11880)',
    structure = SMILES('C[C]OC(=O)[C]=O'),
    E0 = (117.812,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.59675,0.0564452,-6.3028e-05,3.64863e-08,-8.62474e-12,14253.1,22.1212], Tmin=(100,'K'), Tmax=(1011.83,'K')), NASAPolynomial(coeffs=[10.4669,0.0213792,-1.10436e-05,2.235e-09,-1.61953e-13,12458.1,-20.7762], Tmin=(1011.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(117.812,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-O2d)) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(CH2_triplet) + radical(OC=OCJ=O)"""),
)

species(
    label = 'C=[C]OO[CH][C]=O(11148)',
    structure = SMILES('C=[C]OO[CH][C]=O'),
    E0 = (386.539,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,500,795,815,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,242.362],'cm^-1')),
        HinderedRotor(inertia=(0.993615,'amu*angstrom^2'), symmetry=1, barrier=(41.6506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.283731,'amu*angstrom^2'), symmetry=1, barrier=(11.8627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0151613,'amu*angstrom^2'), symmetry=1, barrier=(41.6601,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.993881,'amu*angstrom^2'), symmetry=1, barrier=(41.678,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.70567,0.0573674,-6.86296e-05,3.63059e-08,-1.18264e-12,46565.9,27.4075], Tmin=(100,'K'), Tmax=(562.806,'K')), NASAPolynomial(coeffs=[7.38954,0.0268941,-1.38594e-05,2.7569e-09,-1.96255e-13,45768.9,1.85713], Tmin=(562.806,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.539,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(216.176,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-O2d)OsHH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(OCJC=O) + radical(C=CJO) + radical(CsCJ=O)"""),
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
    label = '[O]C(=O)C=O(11104)',
    structure = SMILES('[O]C(=O)C=O'),
    E0 = (-220.22,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2782.5,750,1395,475,1775,1000,180,557.349,678.737,678.937,3886.95],'cm^-1')),
        HinderedRotor(inertia=(1.37256,'amu*angstrom^2'), symmetry=1, barrier=(31.558,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (73.0275,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.94681,0.0265975,-2.9742e-05,1.90551e-08,-5.37244e-12,-26451.4,14.3847], Tmin=(100,'K'), Tmax=(825.92,'K')), NASAPolynomial(coeffs=[5.38718,0.0147781,-8.27519e-06,1.72677e-09,-1.27055e-13,-26854.5,3.07817], Tmin=(825.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-220.22,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(128.874,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-O2d)H) + group(Cds-O2d(Cds-O2d)O2s) + group(Cds-O2d(Cds-O2d)H) + radical(C=OC=OOJ)"""),
)

species(
    label = 'C=[C]OC([O])[C]=O(11147)',
    structure = SMILES('C=[C]OC([O])[C]=O'),
    E0 = (171.401,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,420.166,420.172,420.178,420.191],'cm^-1')),
        HinderedRotor(inertia=(0.000954952,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121789,'amu*angstrom^2'), symmetry=1, barrier=(15.2582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.121793,'amu*angstrom^2'), symmetry=1, barrier=(15.2581,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.30631,0.0553425,-5.9661e-05,3.16319e-08,-6.52806e-12,20715.1,29.011], Tmin=(100,'K'), Tmax=(1187.92,'K')), NASAPolynomial(coeffs=[14.1784,0.0119992,-4.93106e-06,9.17167e-10,-6.41072e-14,17656.9,-35.306], Tmin=(1187.92,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(171.401,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CJO) + radical(CsCJ=O)"""),
)

species(
    label = '[CH]=[C]OC([O])C=O(11883)',
    structure = SMILES('[CH]=[C]OC([O])C=O'),
    E0 = (263.826,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,1685,370,3120,650,792.5,1650,370.283,370.664,370.896],'cm^-1')),
        HinderedRotor(inertia=(0.0012245,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152065,'amu*angstrom^2'), symmetry=1, barrier=(14.8446,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152336,'amu*angstrom^2'), symmetry=1, barrier=(14.8468,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (99.0648,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.19954,0.0583001,-6.61246e-05,3.67408e-08,-7.91792e-12,31834.7,29.1185], Tmin=(100,'K'), Tmax=(1141.16,'K')), NASAPolynomial(coeffs=[14.5233,0.0115976,-4.73673e-06,8.78176e-10,-6.13297e-14,28793.7,-36.9206], Tmin=(1141.16,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(263.826,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsOsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CJO) + radical(Cds_P)"""),
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
    E0 = (126.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (492.139,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (901.326,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (857.244,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (134.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (362.978,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (133.667,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (189.954,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (165.779,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (357.77,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (308.347,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (252.789,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (249.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (259.973,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (323.038,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (173.272,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (365.305,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (541.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (585.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (368.189,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (239.343,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (279.211,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (227.731,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (258.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (281.523,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (126.554,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (628.548,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (134.838,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (368.724,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (487.56,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (268.291,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (539.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (194.932,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (214.636,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (278.803,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (208.704,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (573.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (495.238,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (200.173,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (700.339,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (129.483,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (175.614,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (417.316,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (329.6,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (472.841,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['OCHCO(3676)', 'ketene(1375)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[O]C=C([O])[O](11106)', 'H2CC(T)(1341)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(30799.3,'m^3/(mol*s)'), n=0.932885, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [O_sec_rad;Birad] for rate rule [O_rad/OneDe;Birad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Birad_R_Recombination
Ea raised from -2.0 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', 'C=[C]O[C]=C[O](9888)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]C(=O)O[C]=C(2687)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=C1OC(=C[O])O1(11854)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=[C]OC1=COO1(11855)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(236.424,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 234.6 to 236.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=C1OC=C([O])O1(11856)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=[C]OC(O)=C=O(11857)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=COC([O])=C=O(11858)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(4.48312e+08,'s^-1'), n=0.311, Ea=(39.225,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_De;XH_Rrad] for rate rule [R4radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=[C]O[C]1OC1[O](11859)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(231.216,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=[C]OC1([O])[CH]O1(11860)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(181.793,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=C1O[C]([O])C1[O](11861)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(126.236,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=C1OC1([O])[CH][O](11862)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(122.839,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['H(3)', 'C=[C]OC([O])=C=O(11863)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction15',
    reactants = ['H(3)', '[CH]=[C]OC(=O)C=O(11368)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[O][C]=C[O](9592)', 'ketene(1375)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(11.6997,'m^3/(mol*s)'), n=2.021, Ea=(29.883,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Od_R;YJ] for rate rule [Od_Cdd;CJ]
Euclidian distance = 1.4142135623730951
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[O][C]=C[O](9592)', '[CH2][C]=O(1376)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(6.16957e+06,'m^3/(mol*s)'), n=-1.87134e-07, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=1.52051667475e-07, var=0.241111252513, Tref=1000.0, N=4, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O
    Total Standard Deviation in ln(k): 0.984387063055
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', 'C=[C]OC([O])=[C][O](11864)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH]=[C]OC([O])=C[O](11375)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['C=[C]OC([O])=[C]O(11865)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]=COC([O])=C[O](11300)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(1.08e+06,'s^-1'), n=1.99, Ea=(105.437,'kJ/mol'), T0=(1,'K'), comment="""From training reaction 17 used for R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd
Exact match found for rate rule [R2H_D;Cd_rad_out_singleH;Cd_H_out_singleNd]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=[C]OC(O)=[C][O](11866)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(8.12004e+07,'s^-1'), n=1.4524, Ea=(152.657,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SR;O_rad_out;XH_out] + [R3H_SD;Y_rad_out;XH_out] for rate rule [R3H_SD;O_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=COC([O])=[C][O](11867)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(1.286e+08,'s^-1'), n=1.323, Ea=(101.177,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSR;Cd_rad_out_Cd;XH_out] for rate rule [R4H_SSD;Cd_rad_out_Cd;XH_out]
Euclidian distance = 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]=[C]OC(O)=C[O](11868)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH]=[C]OC([O])=CO(11869)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6HJ_1;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 2.23606797749979
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=[C]OC(=O)C=O(11153)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH][O](1548)', 'C=[C]O[C]=O(2686)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=C1OC(=O)C1[O](11870)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=C1O[C]([CH][O])O1(11871)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra] for rate rule [R4_S_CO;carbonyl_intra;radadd_intra_cddouble]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction30',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=[C]O[C]1[CH]OO1(11872)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(361.006,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra_O] for rate rule [R4_linear;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 360.6 to 361.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=[C]OC(=O)C[O](11873)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH2][C]OC([O])[C]=O(11874)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['O=C[C]1O[C]CO1(11875)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.19156e+09,'s^-1'), n=0.640131, Ea=(68.3777,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;carbonyl_intra_De;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['O=C1[CH]OC[C]O1(11876)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.8912e+08,'s^-1'), n=0.529986, Ea=(88.0823,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['[O]C1(C=O)C[C]O1(11877)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.54938e+07,'s^-1'), n=1.23279, Ea=(152.249,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5;carbonylbond_intra_De;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 150.0 to 152.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['[O]C1C[C]OC1=O(11878)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(6.93455e+07,'s^-1'), n=0.812267, Ea=(82.1501,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6;carbonylbond_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['OCHCO(3676)', '[CH2][C][O](2821)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(1.505e+06,'m^3/(mol*s)'), n=1.79841e-08, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_4R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_4R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_4R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction38',
    reactants = ['HCO(1372)', '[CH2][C]O[C]=O(11879)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.13947e+08,'m^3/(mol*s)'), n=-0.514474, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=-0.000121059371351, var=1.09982367998, Tref=1000.0, N=21, correlation='Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R
    Total Standard Deviation in ln(k): 2.10271953737
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_N-1COS->O_1CS->C_N-1C-inRing_Ext-2R-R_Ext-2R-R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C[C]OC(=O)[C]=O(11880)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(4.23647e+06,'s^-1'), n=1.3192, Ea=(73.619,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;C_rad_out_2H;XH_out] for rate rule [R5Hall;C_rad_out_2H;CO_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['C=[C]OO[CH][C]=O(11148)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(7040,'s^-1'), n=2.66, Ea=(313.8,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R_ROR;R1_doublebond;R2_doublebond_H;R_O_R]
Euclidian distance = 0
family: ketoenol"""),
)

reaction(
    label = 'reaction41',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=C1OC1([O])C=O(11881)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction42',
    reactants = ['C=[C]OC([O])=C[O](11160)'],
    products = ['C=C1O[CH][C]([O])O1(11882)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(2.55236e+10,'s^-1'), n=0.405358, Ea=(49.0602,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cddouble] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[O]C(=O)C=O(11104)', 'H2CC(T)(1341)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(3.73498e+08,'m^3/(mol*s)'), n=-0.730906, Ea=(37.2121,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R_R;Y_1centerbirad] + [Od_R;YJ] for rate rule [Od_R;Y_1centerbirad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction44',
    reactants = ['C=[C]OC([O])[C]=O(11147)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Y_rad_out;XH_out] for rate rule [R2H_S;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]=[C]OC([O])C=O(11883)'],
    products = ['C=[C]OC([O])=C[O](11160)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(868165,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4Hall;Cd_rad_out_singleH;XH_out] for rate rule [R4HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #2828',
    isomers = [
        'C=[C]OC([O])=C[O](11160)',
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
    label = 'PDepNetwork #2828',
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


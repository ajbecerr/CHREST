species(
    label = 'C=[C][CH]C([O])=C[O](16700)',
    structure = SMILES('[CH2][C]=CC([O])=C[O]'),
    E0 = (270.375,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1685,370,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.51195,'amu*angstrom^2'), symmetry=1, barrier=(34.7627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51303,'amu*angstrom^2'), symmetry=1, barrier=(34.7875,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0219048,0.0669119,-6.99868e-05,3.47026e-08,-6.32141e-12,32679.5,28.1347], Tmin=(100,'K'), Tmax=(1596.64,'K')), NASAPolynomial(coeffs=[18.949,0.00549596,8.63194e-07,-3.71551e-10,3.02657e-14,28419.9,-66.446], Tmin=(1596.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    label = '[CH2][C]=C[C]=C[O](18003)',
    structure = SMILES('[CH2][C]=C[C]=C[O]'),
    E0 = (545.459,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,3000,3100,440,815,1455,1000,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.90296,'amu*angstrom^2'), symmetry=1, barrier=(43.7528,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.90469,'amu*angstrom^2'), symmetry=1, barrier=(43.7927,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.5791,0.0435018,-1.46016e-05,-2.37305e-08,1.61967e-11,65699.9,21.3905], Tmin=(100,'K'), Tmax=(893.788,'K')), NASAPolynomial(coeffs=[15.7012,0.00765852,-3.613e-07,-1.05605e-10,9.57888e-15,62082.7,-51.2681], Tmin=(893.788,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.459,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C([O])[CH][C]=C(16203)',
    structure = SMILES('[CH]C([O])=C[C]=C'),
    E0 = (550.822,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.16134,'amu*angstrom^2'), symmetry=1, barrier=(49.6934,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.15849,'amu*angstrom^2'), symmetry=1, barrier=(49.6278,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3734.97,'J/mol'), sigma=(6.1309,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=583.39 K, Pc=36.78 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.973876,0.0578595,-5.6507e-05,2.89348e-08,-5.76231e-12,66364.9,22.9252], Tmin=(100,'K'), Tmax=(1341.27,'K')), NASAPolynomial(coeffs=[13.7717,0.0154074,-4.23793e-06,5.72638e-10,-3.18357e-14,63317.4,-41.1376], Tmin=(1341.27,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(550.822,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
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
    label = '[CH][C]=CC([O])=C[O](18837)',
    structure = SMILES('[CH][C]=CC([O])=C[O]'),
    E0 = (523.004,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.08306,'amu*angstrom^2'), symmetry=1, barrier=(47.8937,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.08482,'amu*angstrom^2'), symmetry=1, barrier=(47.9341,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.150494,0.0710704,-7.71397e-05,4.04129e-08,-7.95212e-12,63053.1,27.2133], Tmin=(100,'K'), Tmax=(1408.18,'K')), NASAPolynomial(coeffs=[19.0925,0.00849722,-1.14736e-06,1.49245e-11,4.83977e-15,58587.6,-67.5682], Tmin=(1408.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.004,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1=CC(=C[O])O1(18838)',
    structure = SMILES('C=C1[CH]C(=C[O])O1'),
    E0 = (121.314,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69326,0.0262878,6.35769e-05,-1.14723e-07,4.91378e-11,14696.2,21.1756], Tmin=(100,'K'), Tmax=(940.304,'K')), NASAPolynomial(coeffs=[21.6984,0.00331565,1.11303e-06,-1.68948e-10,-1.12425e-15,8187.37,-88.7113], Tmin=(940.304,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(121.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=COJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2][C]=CC1=COO1(18839)',
    structure = SMILES('C=[C]C=C1[CH]OO1'),
    E0 = (458.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.86715,0.0275495,4.71851e-05,-8.7938e-08,3.72143e-11,55290.7,22.6602], Tmin=(100,'K'), Tmax=(952.721,'K')), NASAPolynomial(coeffs=[17.4295,0.0107702,-2.85033e-06,5.72603e-10,-4.97926e-14,50121.6,-63.2312], Tmin=(952.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(458.928,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=CCJO) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C1=CC([O])=CO1(18840)',
    structure = SMILES('C=C1C=C([O])[CH]O1'),
    E0 = (-10.0123,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.01202,0.0196018,7.59831e-05,-1.27305e-07,5.45536e-11,-1110.2,19.156], Tmin=(100,'K'), Tmax=(919.471,'K')), NASAPolynomial(coeffs=[20.4922,0.00168053,3.30178e-06,-7.11018e-10,4.13509e-14,-7149.43,-82.8092], Tmin=(919.471,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-10.0123,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=C(C)OJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH2][C]=CC(O)=C=O(18841)',
    structure = SMILES('C=[C]C=C(O)[C]=O'),
    E0 = (123.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0688516,0.0832448,-0.000112197,7.14008e-08,-1.72649e-11,14997.2,22.9704], Tmin=(100,'K'), Tmax=(1029.6,'K')), NASAPolynomial(coeffs=[19.4421,0.00798002,-2.5462e-06,4.01995e-10,-2.55523e-14,11007.8,-71.0594], Tmin=(1029.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(123.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + radical(C=CJC=C) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH2]C=CC([O])=C=O(18842)',
    structure = SMILES('[CH2]C=CC([O])=C=O'),
    E0 = (87.6442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.25648,0.0571152,-6.05029e-05,3.41553e-08,-7.63225e-12,10642.7,23.8102], Tmin=(100,'K'), Tmax=(1095.54,'K')), NASAPolynomial(coeffs=[12.1592,0.0173077,-5.99903e-06,9.88291e-10,-6.36161e-14,8253.84,-29.7839], Tmin=(1095.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(87.6442,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(C=CC=CCJ)"""),
)

species(
    label = '[CH2][C]=[C]C([O])[CH][O](18843)',
    structure = SMILES('[CH2][C]=[C]C([O])[CH][O]'),
    E0 = (878.365,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1670,1700,300,440,318.706,320.69,2516.05,2519.95],'cm^-1')),
        HinderedRotor(inertia=(0.448164,'amu*angstrom^2'), symmetry=1, barrier=(32.589,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.148067,'amu*angstrom^2'), symmetry=1, barrier=(10.6882,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.445885,'amu*angstrom^2'), symmetry=1, barrier=(32.5636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.875144,0.0776725,-0.000132634,1.20815e-07,-4.24452e-11,105747,30.5154], Tmin=(100,'K'), Tmax=(845.886,'K')), NASAPolynomial(coeffs=[7.41456,0.0289055,-1.4514e-05,2.78307e-09,-1.90878e-13,105279,3.83442], Tmin=(845.886,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(878.365,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=C[C]1OC1[O](18844)',
    structure = SMILES('C=[C][CH][C]1OC1[O]'),
    E0 = (537.129,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2039,0.064381,-8.65982e-05,6.59977e-08,-1.96855e-11,64699.6,23.8709], Tmin=(100,'K'), Tmax=(953.52,'K')), NASAPolynomial(coeffs=[8.58884,0.0241557,-8.7746e-06,1.41732e-09,-8.71884e-14,63711.6,-9.20175], Tmin=(953.52,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(537.129,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC1([O])[CH]O1(18845)',
    structure = SMILES('C=[C][CH]C1([O])[CH]O1'),
    E0 = (535.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.294265,0.0730837,-7.96911e-05,4.15895e-08,-7.91529e-12,64530.1,26.1056], Tmin=(100,'K'), Tmax=(1552.99,'K')), NASAPolynomial(coeffs=[18.904,0.00614826,1.85062e-06,-6.65343e-10,5.40959e-14,60675.9,-68.1767], Tmin=(1552.99,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CCsJO) + radical(C=CCJCO) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=C[C]([O])C1[O](18846)',
    structure = SMILES('[CH2][C]1[CH]C(=O)C1[O]'),
    E0 = (444.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.51557,0.0101091,9.01401e-05,-1.30665e-07,5.25587e-11,53581.5,26.7064], Tmin=(100,'K'), Tmax=(935.124,'K')), NASAPolynomial(coeffs=[16.153,0.00935429,-1.00981e-06,1.62927e-10,-2.08176e-14,48513.5,-51.6323], Tmin=(935.124,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(444.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(Cyclobutanone) + radical(C=OCOJ) + radical(CCJ(C)CO) + radical(CCJC=O) + radical(Isobutyl)"""),
)

species(
    label = '[O][C]1C=[C]CC1[O](18847)',
    structure = SMILES('[O]C1=C[C]CC1[O]'),
    E0 = (451.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.73526,0.0316212,3.32798e-05,-7.58189e-08,3.38968e-11,54342.5,22.6393], Tmin=(100,'K'), Tmax=(941.98,'K')), NASAPolynomial(coeffs=[18.3383,0.00709313,-8.71258e-07,1.63172e-10,-1.96159e-14,49174.8,-67.2949], Tmin=(941.98,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2]C1=CC1([O])[CH][O](18848)',
    structure = SMILES('C=C1[CH]C1([O])[CH][O]'),
    E0 = (557.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32388,0.0601396,-6.22956e-05,3.4206e-08,-7.65263e-12,67148.6,22.4164], Tmin=(100,'K'), Tmax=(1070.84,'K')), NASAPolynomial(coeffs=[11.1955,0.0232641,-1.06403e-05,2.04638e-09,-1.44392e-13,65034.5,-25.884], Tmin=(1070.84,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(557.512,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(C=CCJCO) + radical(CCsJOH)"""),
)

species(
    label = '[O][CH]C1([O])C=[C]C1(18849)',
    structure = SMILES('[O][CH]C1([O])C=[C]C1'),
    E0 = (650.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.91825,0.0491185,-4.57577e-05,2.46956e-08,-5.78156e-12,78321.8,24.6864], Tmin=(100,'K'), Tmax=(991.425,'K')), NASAPolynomial(coeffs=[7.33652,0.0272583,-1.26843e-05,2.45632e-09,-1.73751e-13,77247.5,-1.40708], Tmin=(991.425,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(650.605,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclobutene) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(cyclobutene-vinyl)"""),
)

species(
    label = '[CH2]C#CC([O])=C[O](18850)',
    structure = SMILES('[CH2]C#CC([O])=C[O]'),
    E0 = (241.515,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,2100,2250,500,550,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.47929,'amu*angstrom^2'), symmetry=1, barrier=(34.0119,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.47635,'amu*angstrom^2'), symmetry=1, barrier=(33.9441,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12796,0.0521538,-3.09926e-05,-1.08605e-08,1.18214e-11,29161.2,23.6885], Tmin=(100,'K'), Tmax=(927.245,'K')), NASAPolynomial(coeffs=[18.7218,0.00465428,-9.20355e-08,-4.81662e-11,1.05064e-15,24677.6,-66.4457], Tmin=(927.245,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(241.515,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Propargyl)"""),
)

species(
    label = '[CH2][C]=CC([O])=C=O(18851)',
    structure = SMILES('[CH2][C]=CC([O])=C=O'),
    E0 = (325.486,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,2120,512.5,787.5,180,745.45],'cm^-1')),
        HinderedRotor(inertia=(2.77854,'amu*angstrom^2'), symmetry=1, barrier=(63.8842,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.518683,'amu*angstrom^2'), symmetry=1, barrier=(11.9255,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44792,0.0587784,-7.34456e-05,4.52529e-08,-9.25122e-12,39236.5,23.6147], Tmin=(100,'K'), Tmax=(719.357,'K')), NASAPolynomial(coeffs=[10.491,0.0176493,-6.7739e-06,1.15686e-09,-7.49957e-14,37698.6,-18.6807], Tmin=(719.357,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(325.486,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CCO)H) + group(Cds-(Cdd-O2d)CsOs) + missing(Cdd-CdO2d) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
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
    label = '[CH2][C]=[C]C([O])=C[O](18852)',
    structure = SMILES('[CH2][C]=[C]C([O])=C[O]'),
    E0 = (469.37,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.95057,'amu*angstrom^2'), symmetry=1, barrier=(44.8475,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.95436,'amu*angstrom^2'), symmetry=1, barrier=(44.9347,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.414072,0.068197,-8.19407e-05,4.69949e-08,-9.97277e-12,56590.6,26.3043], Tmin=(100,'K'), Tmax=(1337.49,'K')), NASAPolynomial(coeffs=[17.738,0.00539304,8.23614e-07,-4.04237e-10,3.5668e-14,52939.8,-58.6355], Tmin=(1337.49,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(469.37,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC([O])=[C][O](18853)',
    structure = SMILES('[CH2][C]=CC([O])=[C][O]'),
    E0 = (510.119,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,380.01,380.01,380.01],'cm^-1')),
        HinderedRotor(inertia=(0.196516,'amu*angstrom^2'), symmetry=1, barrier=(20.1379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.717295,'amu*angstrom^2'), symmetry=1, barrier=(73.5047,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.18383,0.0575592,-5.84194e-05,2.39486e-08,-1.44933e-12,61458.7,27.0851], Tmin=(100,'K'), Tmax=(885.25,'K')), NASAPolynomial(coeffs=[14.9558,0.00989375,-2.32952e-06,2.91449e-10,-1.63992e-14,58449.8,-40.9009], Tmin=(885.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(510.119,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C=[C]C([O])=C[O](18854)',
    structure = SMILES('[CH2]C=[C]C([O])=C[O]'),
    E0 = (231.529,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0279479,0.0683747,-7.31432e-05,3.7527e-08,-7.07562e-12,28005.8,27.2275], Tmin=(100,'K'), Tmax=(1545.33,'K')), NASAPolynomial(coeffs=[18.4412,0.0063824,9.40969e-07,-4.34443e-10,3.65367e-14,24026,-64.0836], Tmin=(1545.33,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(231.529,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2][C]=CC([O])=[C]O(18855)',
    structure = SMILES('[CH2][C]=CC([O])=[C]O'),
    E0 = (368.657,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3100,440,815,1455,1000,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.966457,'amu*angstrom^2'), symmetry=1, barrier=(22.2207,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.966309,'amu*angstrom^2'), symmetry=1, barrier=(22.2173,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.96639,'amu*angstrom^2'), symmetry=1, barrier=(22.2192,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.269661,0.0706096,-8.38475e-05,4.74844e-08,-9.96357e-12,44483.3,28.9455], Tmin=(100,'K'), Tmax=(1346.37,'K')), NASAPolynomial(coeffs=[18.5421,0.00549589,7.58906e-07,-3.82121e-10,3.35876e-14,40544.4,-60.9985], Tmin=(1346.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.657,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=[C]C(O)=C[O](18856)',
    structure = SMILES('[CH2][C]=[C]C(O)=C[O]'),
    E0 = (331.566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.62619,'amu*angstrom^2'), symmetry=1, barrier=(37.3894,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62732,'amu*angstrom^2'), symmetry=1, barrier=(37.4152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.62679,'amu*angstrom^2'), symmetry=1, barrier=(37.4031,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.305587,0.0762732,-8.84873e-05,4.77718e-08,-9.41105e-12,40049,27.3903], Tmin=(100,'K'), Tmax=(1474.51,'K')), NASAPolynomial(coeffs=[20.9736,0.00262075,2.64073e-06,-7.55277e-10,5.86434e-14,35505.1,-77.6619], Tmin=(1474.51,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(331.566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]=CC(O)=[C][O](18857)',
    structure = SMILES('[CH2][C]=CC(O)=[C][O]'),
    E0 = (372.314,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.19832,'amu*angstrom^2'), symmetry=1, barrier=(27.5517,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.1967,'amu*angstrom^2'), symmetry=1, barrier=(27.5145,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.19738,'amu*angstrom^2'), symmetry=1, barrier=(27.5301,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.206757,0.06877,-7.64099e-05,4.03046e-08,-7.88695e-12,44928.2,29.088], Tmin=(100,'K'), Tmax=(1443.15,'K')), NASAPolynomial(coeffs=[19.2232,0.00546189,4.09489e-07,-2.71748e-10,2.38439e-14,40543.3,-65.8071], Tmin=(1443.15,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(372.314,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(C=CC=CCJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C[C]=[C]C([O])=C[O](18858)',
    structure = SMILES('C[C]=[C]C([O])=C[O]'),
    E0 = (351.315,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.12262,'amu*angstrom^2'), symmetry=1, barrier=(25.8113,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.11626,'amu*angstrom^2'), symmetry=1, barrier=(25.6651,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.641895,0.0708924,-8.23952e-05,4.16463e-08,-5.79471e-12,42377.3,23.9437], Tmin=(100,'K'), Tmax=(846.674,'K')), NASAPolynomial(coeffs=[16.8871,0.00981779,-1.96144e-06,1.77956e-10,-6.27262e-15,39064.6,-55.044], Tmin=(846.674,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(351.315,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S) + radical(C=CJC=C)"""),
)

species(
    label = '[CH2]C=CC([O])=[C][O](18859)',
    structure = SMILES('[CH2]C=CC([O])=[C][O]'),
    E0 = (272.277,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.541268,0.0608598,-6.10253e-05,3.00094e-08,-5.53158e-12,32885,28.9217], Tmin=(100,'K'), Tmax=(1533.02,'K')), NASAPolynomial(coeffs=[16.8554,0.00899195,-1.17449e-06,2.4553e-11,3.60726e-15,28975.9,-53.1903], Tmin=(1533.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(272.277,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CC=CCJ) + radical(C=CJO)"""),
)

species(
    label = '[CH2][C]=[C]C([O])=CO(18860)',
    structure = SMILES('[CH2][C]=[C]C([O])=CO'),
    E0 = (327.908,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1670,1700,300,440,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.50371,'amu*angstrom^2'), symmetry=1, barrier=(34.5732,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50495,'amu*angstrom^2'), symmetry=1, barrier=(34.6018,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.50383,'amu*angstrom^2'), symmetry=1, barrier=(34.576,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.233046,0.0780096,-9.56085e-05,5.45993e-08,-1.13606e-11,39603.7,27.2125], Tmin=(100,'K'), Tmax=(1396.57,'K')), NASAPolynomial(coeffs=[20.4611,0.00241961,3.10686e-06,-8.90229e-10,7.02515e-14,35414.9,-73.8395], Tmin=(1396.57,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(327.908,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=CC=CCJ) + radical(C=CJC=C) + radical(Cds_S)"""),
)

species(
    label = 'C[C]=CC([O])=[C][O](18861)',
    structure = SMILES('C[C]=CC([O])=[C][O]'),
    E0 = (392.063,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,3010,987.5,1337.5,450,1655,350,440,435,1725,1670,1700,300,440,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.611111,'amu*angstrom^2'), symmetry=1, barrier=(14.0506,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.612739,'amu*angstrom^2'), symmetry=1, barrier=(14.0881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.873592,0.0668635,-8.33342e-05,5.25036e-08,-1.28183e-11,47268.6,26.638], Tmin=(100,'K'), Tmax=(1012.73,'K')), NASAPolynomial(coeffs=[14.3053,0.0138125,-4.75898e-06,7.79141e-10,-4.987e-14,44548,-38.3323], Tmin=(1012.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(392.063,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = 'C=[C]C1OC1=C[O](18809)',
    structure = SMILES('C=[C]C1OC1=C[O]'),
    E0 = (252.725,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.11855,0.0414349,2.4012e-05,-8.15959e-08,3.99724e-11,30520,21.8457], Tmin=(100,'K'), Tmax=(923.326,'K')), NASAPolynomial(coeffs=[24.6577,-0.00362592,4.75438e-06,-9.31497e-10,5.58102e-14,23747,-102.977], Tmin=(923.326,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(252.725,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(C=COJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1OC=C1[O](18862)',
    structure = SMILES('C=[C]C1OC=C1[O]'),
    E0 = (223.836,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23893,0.0383138,3.20505e-05,-8.92232e-08,4.25866e-11,27041.5,22.2756], Tmin=(100,'K'), Tmax=(921.751,'K')), NASAPolynomial(coeffs=[24.2028,-0.00311921,4.73196e-06,-9.40019e-10,5.66417e-14,20334.9,-100.057], Tmin=(921.751,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(223.836,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]CC([O])=C=O(16698)',
    structure = SMILES('C=[C]CC(=O)[C]=O'),
    E0 = (196.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,180,723.215],'cm^-1')),
        HinderedRotor(inertia=(0.0453126,'amu*angstrom^2'), symmetry=1, barrier=(16.6581,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.096708,'amu*angstrom^2'), symmetry=1, barrier=(16.6497,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.724218,'amu*angstrom^2'), symmetry=1, barrier=(16.6512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.082,0.0478428,-3.75715e-05,1.4209e-08,-2.24845e-12,23717.9,22.2666], Tmin=(100,'K'), Tmax=(1396.29,'K')), NASAPolynomial(coeffs=[9.87225,0.0255259,-1.35971e-05,2.76228e-09,-1.9897e-13,21542.4,-17.9174], Tmin=(1396.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = 'C=[C]C1[C]([O])C1[O](18863)',
    structure = SMILES('C=[C]C1[C]([O])C1[O]'),
    E0 = (644.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26499,0.0522852,-3.87393e-05,7.27959e-09,2.40013e-12,77629.6,24.9845], Tmin=(100,'K'), Tmax=(996.837,'K')), NASAPolynomial(coeffs=[14.8741,0.0137702,-5.00118e-06,9.12529e-10,-6.4999e-14,74116.8,-44.6391], Tmin=(996.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(644.573,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsOsH) + group(Cs-CsCsOsH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Cyclopropane) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Cds_S)"""),
)

species(
    label = 'C#C[CH]C([O])=C[O](18864)',
    structure = SMILES('[CH]=C=CC([O])=C[O]'),
    E0 = (245.559,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,540,610,2055,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.446,'amu*angstrom^2'), symmetry=1, barrier=(33.2463,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.107982,0.0713042,-8.44489e-05,4.5646e-08,-8.92089e-12,29698.2,26.2344], Tmin=(100,'K'), Tmax=(1502.14,'K')), NASAPolynomial(coeffs=[20.4692,-0.000266934,3.7732e-06,-9.43033e-10,7.03025e-14,25409,-75.1108], Tmin=(1502.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(245.559,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = 'C=[C]C=C=C[O](18009)',
    structure = SMILES('C=[C]C=C=C[O]'),
    E0 = (366.165,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,1685,370,540,610,2055,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.63596,'amu*angstrom^2'), symmetry=1, barrier=(37.614,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.44293,0.0581897,-5.88991e-05,2.83074e-08,-4.99587e-12,44184.9,22.5048], Tmin=(100,'K'), Tmax=(1656.62,'K')), NASAPolynomial(coeffs=[16.9923,0.00540494,5.08161e-07,-2.72965e-10,2.24497e-14,40461.6,-60.3782], Tmin=(1656.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(366.165,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + group(Cds-CdsOsH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=C[CH]C([O])=C[O](18865)',
    structure = SMILES('[CH]C=CC([O])=C[O]'),
    E0 = (285.162,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2995,3010,3025,975,987.5,1000,1300,1337.5,1375,400,450,500,1630,1655,1680,350,440,435,1725,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.06803,'amu*angstrom^2'), symmetry=1, barrier=(47.5481,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06884,'amu*angstrom^2'), symmetry=1, barrier=(47.5668,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.710753,0.0602353,-3.04234e-05,-1.75558e-08,1.54107e-11,34426.8,24.7305], Tmin=(100,'K'), Tmax=(913.918,'K')), NASAPolynomial(coeffs=[19.4486,0.0105,-1.76698e-06,1.82603e-10,-1.20752e-14,29653.9,-71.3564], Tmin=(913.918,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.162,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]CC([O])=[C][O](16701)',
    structure = SMILES('C=[C]CC([O])=[C][O]'),
    E0 = (422.252,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1670,1700,300,440,334.596,335.889,337.277,338.369],'cm^-1')),
        HinderedRotor(inertia=(0.00146392,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.137636,'amu*angstrom^2'), symmetry=1, barrier=(11.1432,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03837,0.0593246,-6.37679e-05,3.44996e-08,-7.23847e-12,50896.8,29.8835], Tmin=(100,'K'), Tmax=(1175.46,'K')), NASAPolynomial(coeffs=[14.7738,0.012583,-4.12012e-06,6.69593e-10,-4.32821e-14,47667.7,-38.6025], Tmin=(1175.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(422.252,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S) + radical(C=CJO)"""),
)

species(
    label = '[CH]=[C]CC([O])=C[O](16702)',
    structure = SMILES('[CH]=[C]CC([O])=C[O]'),
    E0 = (429.604,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,225.224,226.157,226.719],'cm^-1')),
        HinderedRotor(inertia=(0.46898,'amu*angstrom^2'), symmetry=1, barrier=(16.8779,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.465946,'amu*angstrom^2'), symmetry=1, barrier=(16.8687,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.293993,0.0664589,-7.20796e-05,3.68989e-08,-7.0232e-12,51815.7,29.8339], Tmin=(100,'K'), Tmax=(1472.35,'K')), NASAPolynomial(coeffs=[19.3893,0.00508522,1.21664e-07,-1.73881e-10,1.54663e-14,47222.1,-66.1816], Tmin=(1472.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C][CH]C(O)=C[O](18866)',
    structure = SMILES('[CH][C]=CC(O)=C[O]'),
    E0 = (385.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.97914,'amu*angstrom^2'), symmetry=1, barrier=(45.5043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97914,'amu*angstrom^2'), symmetry=1, barrier=(45.5043,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97924,'amu*angstrom^2'), symmetry=1, barrier=(45.5067,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.618969,0.07966,-8.5174e-05,4.27427e-08,-7.91426e-12,46513.9,28.4834], Tmin=(100,'K'), Tmax=(1548.87,'K')), NASAPolynomial(coeffs=[22.364,0.005628,7.36747e-07,-3.53441e-10,2.93277e-14,41154.9,-86.7684], Tmin=(1548.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(385.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C][CH]C([O])=CO(18867)',
    structure = SMILES('[CH][C]=CC([O])=CO'),
    E0 = (381.541,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.99206,'amu*angstrom^2'), symmetry=1, barrier=(45.8015,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99216,'amu*angstrom^2'), symmetry=1, barrier=(45.8038,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.99211,'amu*angstrom^2'), symmetry=1, barrier=(45.8025,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.512187,0.0810469,-9.12972e-05,4.85467e-08,-9.5253e-12,46066.9,28.1787], Tmin=(100,'K'), Tmax=(1457.23,'K')), NASAPolynomial(coeffs=[21.7704,0.00557675,1.11373e-06,-4.67115e-10,3.9168e-14,41091.7,-82.5002], Tmin=(1457.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.541,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=C=CC([O])=C[O](16697)',
    structure = SMILES('C=C=CC([O])=C[O]'),
    E0 = (91.0818,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2995,3025,975,1000,1300,1375,400,500,1630,1680,2950,3100,1380,975,1025,1650,540,610,2055,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.32234,'amu*angstrom^2'), symmetry=1, barrier=(30.4032,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.160136,0.0705545,-7.64211e-05,3.8369e-08,-7.0525e-12,11122.4,25.8099], Tmin=(100,'K'), Tmax=(1574.22,'K')), NASAPolynomial(coeffs=[20.8277,0.00289944,1.69495e-06,-4.93392e-10,3.72574e-14,6289.61,-79.3298], Tmin=(1574.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(91.0818,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=COJ)"""),
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
    label = '[CH2][C]=C[C]=O(9897)',
    structure = SMILES('[CH2][C]=C[C]=O'),
    E0 = (432.598,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3000,3100,440,815,1455,1000,1685,370,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.00399,'amu*angstrom^2'), symmetry=1, barrier=(23.0837,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.00442,'amu*angstrom^2'), symmetry=1, barrier=(23.0936,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.33986,0.0429161,-4.87447e-05,1.23757e-08,1.41321e-11,52082.9,16.5479], Tmin=(100,'K'), Tmax=(512.282,'K')), NASAPolynomial(coeffs=[6.6704,0.019476,-1.04855e-05,2.11536e-09,-1.51263e-13,51503.1,-2.77645], Tmin=(512.282,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(432.598,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(Cds_S) + radical(C=CCJ=O)"""),
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
    label = '[CH]=C([O])[CH][O](10249)',
    structure = SMILES('[CH]C([O])=C[O]'),
    E0 = (232.724,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.14782,'amu*angstrom^2'), symmetry=1, barrier=(49.3827,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (70.0468,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.06413,0.0315518,3.76054e-06,-3.87898e-08,2.07567e-11,28070.3,18.0523], Tmin=(100,'K'), Tmax=(902.966,'K')), NASAPolynomial(coeffs=[15.0822,0.00356738,9.38086e-07,-3.00037e-10,2.0679e-14,24509.2,-50.1243], Tmin=(902.966,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(232.724,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(153.818,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C=C1OC1[O](18868)',
    structure = SMILES('C=[C]C=C1OC1[O]'),
    E0 = (293.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.12427,0.0543863,-3.80946e-05,1.2401e-09,6.09308e-12,35451.8,20.3017], Tmin=(100,'K'), Tmax=(946.536,'K')), NASAPolynomial(coeffs=[16.5032,0.0108417,-3.07377e-06,5.10886e-10,-3.65116e-14,31579.8,-58.1224], Tmin=(946.536,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(293.835,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(methyleneoxirane) + radical(CCOJ) + radical(C=CJC=C)"""),
)

species(
    label = 'C=C1C=C([O])C1[O](18869)',
    structure = SMILES('C=C1C=C([O])C1[O]'),
    E0 = (195.498,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.47046,0.0372324,2.28033e-05,-7.03217e-08,3.3518e-11,23621.1,21.3429], Tmin=(100,'K'), Tmax=(929.933,'K')), NASAPolynomial(coeffs=[20.3076,0.0034081,1.22539e-06,-2.6999e-10,1.17311e-14,18076.7,-79.1407], Tmin=(929.933,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(195.498,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)OsH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(C=C(C)OJ)"""),
)

species(
    label = 'C=C=C=C(O)[CH][O](18870)',
    structure = SMILES('[CH2]C#CC(O)=C[O]'),
    E0 = (103.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.767014,0.0560546,-2.31519e-05,-2.85104e-08,2.01687e-11,12603.8,23.4835], Tmin=(100,'K'), Tmax=(922.121,'K')), NASAPolynomial(coeffs=[22.4452,0.00132081,1.94987e-06,-4.36614e-10,2.62119e-14,6934.9,-88.4043], Tmin=(922.121,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(103.71,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-CtHHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-Ct(Cds-Cds)) + radical(C=COJ) + radical(Propargyl)"""),
)

species(
    label = 'C=C=C=C([O])C[O](18871)',
    structure = SMILES('[CH2]C#CC(=O)C[O]'),
    E0 = (251.845,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69734,0.0543608,-6.71234e-05,5.09469e-08,-1.63571e-11,30369.4,24.6566], Tmin=(100,'K'), Tmax=(749.929,'K')), NASAPolynomial(coeffs=[6.81708,0.0270511,-1.24951e-05,2.38059e-09,-1.65693e-13,29601.5,1.43048], Tmin=(749.929,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.845,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-CtHHH) + group(Cds-O2d(Cds-Cds)Cs) + group(Ct-CtCs) + group(Ct-CtCs) + radical(C=OCOJ) + radical(Propargyl)"""),
)

species(
    label = 'C=[C]C1O[C]1[CH][O](18872)',
    structure = SMILES('C=[C]C1O[C]1[CH][O]'),
    E0 = (627.934,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.13096,0.068811,-0.000105874,8.73217e-08,-2.71504e-11,75621.1,26.59], Tmin=(100,'K'), Tmax=(967.967,'K')), NASAPolynomial(coeffs=[7.88285,0.023326,-8.14062e-06,1.24298e-09,-7.16381e-14,75137.7,-1.50907], Tmin=(967.967,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(627.934,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(CCsJOH) + radical(Cds_S)"""),
)

species(
    label = 'C=[C]C1O[CH][C]1[O](18873)',
    structure = SMILES('C=[C]C1O[CH][C]1[O]'),
    E0 = (623.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.514949,0.06316,-6.89033e-05,3.76879e-08,-7.61388e-12,75155.6,27.007], Tmin=(100,'K'), Tmax=(1448.87,'K')), NASAPolynomial(coeffs=[15.1343,0.00962286,1.64588e-07,-3.69037e-10,3.58481e-14,72302.3,-44.1707], Tmin=(1448.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(623.739,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Oxetane) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(CCsJOCs) + radical(Cds_S)"""),
)

species(
    label = 'C=C=C[C]=O(9908)',
    structure = SMILES('C=C=C[C]=O'),
    E0 = (213.911,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,2950,3100,1380,975,1025,1650,540,610,2055,1855,455,950],'cm^-1')),
        HinderedRotor(inertia=(1.13978,'amu*angstrom^2'), symmetry=1, barrier=(26.2058,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.00684,0.0284132,-1.9291e-05,6.01948e-09,-7.90367e-13,25757,13.5223], Tmin=(100,'K'), Tmax=(1551.73,'K')), NASAPolynomial(coeffs=[6.76337,0.0187298,-9.93048e-06,1.99793e-09,-1.42457e-13,24591.2,-6.25135], Tmin=(1551.73,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(213.911,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-Cds)H) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(CsCJ=O)"""),
)

species(
    label = 'C=[C][C]=C([O])C[O](18874)',
    structure = SMILES('C=[C][C]=C([O])C[O]'),
    E0 = (451.702,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,2950,3100,1380,975,1025,1650,1670,1700,300,440,180,180,180,1212.22],'cm^-1')),
        HinderedRotor(inertia=(2.85507,'amu*angstrom^2'), symmetry=1, barrier=(65.6436,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.416354,'amu*angstrom^2'), symmetry=1, barrier=(9.5728,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.00553,0.072772,-0.000116633,9.93772e-08,-3.21667e-11,54428.6,25.9168], Tmin=(100,'K'), Tmax=(922.525,'K')), NASAPolynomial(coeffs=[8.03252,0.0247499,-1.00084e-05,1.69892e-09,-1.06887e-13,53879,-3.36931], Tmin=(922.525,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(451.702,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=CJC=C) + radical(C=CJC=C)"""),
)

species(
    label = '[CH]=[C]C=C([O])C[O](18875)',
    structure = SMILES('[CH]=[C]C=C([O])C[O]'),
    E0 = (499.803,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.0632163,'amu*angstrom^2'), symmetry=1, barrier=(1.45347,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.564788,'amu*angstrom^2'), symmetry=1, barrier=(12.9856,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.03418,0.0694573,-0.000103148,8.25128e-08,-2.57348e-11,60215.2,26.6209], Tmin=(100,'K'), Tmax=(895.763,'K')), NASAPolynomial(coeffs=[9.69255,0.0220366,-9.07584e-06,1.58634e-09,-1.03008e-13,59015.3,-12.2366], Tmin=(895.763,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(499.803,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsHH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = 'C=[C][CH]C([O])[C]=O(16747)',
    structure = SMILES('C=[C][CH]C([O])[C]=O'),
    E0 = (465.469,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,2950,3100,1380,975,1025,1650,1685,370,1855,455,950,411.58,411.588,411.588],'cm^-1')),
        HinderedRotor(inertia=(0.00099514,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186935,'amu*angstrom^2'), symmetry=1, barrier=(22.4719,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.186937,'amu*angstrom^2'), symmetry=1, barrier=(22.4719,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17705,0.0585956,-5.9949e-05,3.09068e-08,-6.28298e-12,56087.6,28.5735], Tmin=(100,'K'), Tmax=(1196.79,'K')), NASAPolynomial(coeffs=[13.889,0.0161085,-6.69732e-06,1.24298e-09,-8.64e-14,53044.9,-35.0377], Tmin=(1196.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(465.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(C=CCJCO) + radical(Cds_S) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C1=C[C]([CH][O])O1(18876)',
    structure = SMILES('[CH2][C]1[CH]C(=C[O])O1'),
    E0 = (368.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.21395,0.0803337,-8.9697e-05,4.50799e-08,-8.07644e-12,44556,25.9649], Tmin=(100,'K'), Tmax=(1688.36,'K')), NASAPolynomial(coeffs=[22.6316,-0.00132736,5.21297e-06,-1.22515e-09,8.73842e-14,40091.1,-90.9426], Tmin=(1688.36,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(368.652,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ) + radical(C2CsJOC(O)) + radical(C=CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH2][C]=C[C]1[CH]OO1(18877)',
    structure = SMILES('C=[C][CH][C]1[CH]OO1'),
    E0 = (756.917,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32848,0.0656205,-9.75307e-05,8.87434e-08,-3.27423e-11,91125.6,20.2237], Tmin=(100,'K'), Tmax=(796.485,'K')), NASAPolynomial(coeffs=[5.22726,0.0345098,-1.72251e-05,3.35047e-09,-2.34106e-13,90870.3,4.59753], Tmin=(796.485,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(756.917,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C2CsJOO) + radical(C=CCJCO) + radical(CCsJOO) + radical(Cds_S)"""),
)

species(
    label = '[O][CH][C]1C=[C]CO1(18878)',
    structure = SMILES('[O]C=C1[CH][C]CO1'),
    E0 = (431.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77616,0.0177078,9.96763e-05,-1.64444e-07,7.06868e-11,51950.5,23.4637], Tmin=(100,'K'), Tmax=(917.632,'K')), NASAPolynomial(coeffs=[26.0736,-0.0071828,7.91973e-06,-1.5605e-09,9.61126e-14,44080,-110.256], Tmin=(917.632,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.034,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(CCJCO) + radical(CCJ2_triplet)"""),
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
    label = 'C=[C][CH][C][O](18879)',
    structure = SMILES('[CH2][C][CH][C]=O'),
    E0 = (753.848,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,581.298,1161.95],'cm^-1')),
        HinderedRotor(inertia=(3.41898,'amu*angstrom^2'), symmetry=1, barrier=(78.609,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00208978,'amu*angstrom^2'), symmetry=1, barrier=(14.2588,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(3.42369,'amu*angstrom^2'), symmetry=1, barrier=(78.7175,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.24429,0.0390352,-4.6579e-05,3.03395e-08,-7.90429e-12,90729.8,20.9041], Tmin=(100,'K'), Tmax=(938.08,'K')), NASAPolynomial(coeffs=[8.31814,0.0131356,-5.16451e-06,9.06836e-10,-6.02854e-14,89590.3,-8.01032], Tmin=(938.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(753.848,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(170.447,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(CCJCHO) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = 'C=C1[CH][C]([O])[CH]O1(18880)',
    structure = SMILES('[CH2][C]1[CH]C([O])=CO1'),
    E0 = (285.942,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.632814,0.0724188,-7.92817e-05,4.00631e-08,-7.24795e-12,34583.2,24.89], Tmin=(100,'K'), Tmax=(1671.94,'K')), NASAPolynomial(coeffs=[19.6782,0.00175123,3.92305e-06,-1.01059e-09,7.44457e-14,30876.8,-74.3118], Tmin=(1671.94,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(285.942,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ) + radical(C2CsJOC(O)) + radical(C=CCJCO) + radical(CJC(C)OC)"""),
)

species(
    label = '[CH]=[C][CH]C([O])C=O(18881)',
    structure = SMILES('[CH][C]=CC([O])C=O'),
    E0 = (545.896,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,2782.5,750,1395,475,1775,1000,1685,370,534.126,534.593,535.504,535.566,536.948],'cm^-1')),
        HinderedRotor(inertia=(0.264325,'amu*angstrom^2'), symmetry=1, barrier=(53.4375,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270198,'amu*angstrom^2'), symmetry=1, barrier=(53.3944,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.261489,'amu*angstrom^2'), symmetry=1, barrier=(53.3837,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77019,0.0517246,-4.22009e-05,1.95119e-08,-3.95706e-12,65734.2,28.4331], Tmin=(100,'K'), Tmax=(1117.54,'K')), NASAPolynomial(coeffs=[7.65295,0.0306685,-1.39385e-05,2.65204e-09,-1.854e-13,64419.3,-0.601597], Tmin=(1117.54,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(545.896,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = 'C=[C]C1C(=O)C1[O](18882)',
    structure = SMILES('C=[C]C1C(=O)C1[O]'),
    E0 = (386.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.64095,0.0426466,-1.51532e-05,-1.55743e-08,1.03708e-11,46616.7,24.5369], Tmin=(100,'K'), Tmax=(961.514,'K')), NASAPolynomial(coeffs=[13.8181,0.0139234,-4.56349e-06,8.09653e-10,-5.81206e-14,43261,-39.0059], Tmin=(961.514,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(386.817,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-OdCsCs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(cyclopropanone) + radical(C=OCOJ) + radical(Cds_S)"""),
)

species(
    label = '[CH2]C1=CC1([O])C=O(18883)',
    structure = SMILES('C=C1[CH]C1([O])C=O'),
    E0 = (260.424,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.24977,0.0619518,-6.94373e-05,4.1212e-08,-9.86134e-12,31419.6,18.9587], Tmin=(100,'K'), Tmax=(1011.02,'K')), NASAPolynomial(coeffs=[11.4166,0.0217279,-9.75883e-06,1.8599e-09,-1.30541e-13,29363.9,-30.2016], Tmin=(1011.02,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(260.424,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CC(C)(C=O)OJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C=C=C([O])C=O(18884)',
    structure = SMILES('[CH2]C=[C]C(=O)C=O'),
    E0 = (176.875,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.77277,0.0530801,-5.05875e-05,2.59314e-08,-5.63117e-12,21349.9,23.17], Tmin=(100,'K'), Tmax=(1070.35,'K')), NASAPolynomial(coeffs=[8.91482,0.0263893,-1.31823e-05,2.63337e-09,-1.89418e-13,19821,-11.7717], Tmin=(1070.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(176.875,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(Allyl_P) + radical(C=CJC=O)"""),
)

species(
    label = '[O][C]1[CH]OC[C]=C1(18885)',
    structure = SMILES('[O]C1=C[C]CO[CH]1'),
    E0 = (362.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.94389,0.0342225,1.29626e-05,-4.62763e-08,2.22639e-11,43641.4,19.24], Tmin=(100,'K'), Tmax=(912.957,'K')), NASAPolynomial(coeffs=[12.7989,0.0155211,-3.726e-06,5.33945e-10,-3.56262e-14,40456.7,-38.7274], Tmin=(912.957,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(362.155,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + ring(36dihydro2hpyran) + radical(C=C(C)OJ) + radical(C=CCJ(O)C) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C]=[C]C([O])C=O(18886)',
    structure = SMILES('[CH2][C]=[C]C([O])C=O'),
    E0 = (564.553,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1670,1700,300,440,755.26,3435.47],'cm^-1')),
        HinderedRotor(inertia=(0.505617,'amu*angstrom^2'), symmetry=1, barrier=(11.6251,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.46769,'amu*angstrom^2'), symmetry=1, barrier=(33.7451,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0285636,'amu*angstrom^2'), symmetry=1, barrier=(11.6443,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.67119,0.054587,-6.08224e-05,3.81887e-08,-9.98894e-12,67981,28.4407], Tmin=(100,'K'), Tmax=(914.409,'K')), NASAPolynomial(coeffs=[8.54459,0.0245198,-1.15001e-05,2.22934e-09,-1.57616e-13,66724,-4.10449], Tmin=(914.409,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(564.553,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(Allyl_P) + radical(Cds_S) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C]C([C]=O)C=O(18887)',
    structure = SMILES('[CH2][C]C([C]=O)C=O'),
    E0 = (473.006,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,2782.5,750,1395,475,1775,1000,3000,3100,440,815,1455,1000,1855,455,950,336.336,336.361],'cm^-1')),
        HinderedRotor(inertia=(0.00149279,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00149062,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00149263,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.579285,'amu*angstrom^2'), symmetry=1, barrier=(46.4104,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.69601,0.0489826,-4.14847e-05,1.80211e-08,-3.1925e-12,56973.7,28.5084], Tmin=(100,'K'), Tmax=(1326.79,'K')), NASAPolynomial(coeffs=[11.219,0.0202727,-9.02696e-06,1.71221e-09,-1.19515e-13,54446.7,-20.1274], Tmin=(1326.79,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(473.006,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-CsCsCsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CC(C)CJ=O)"""),
)

species(
    label = '[CH][C][CH2](17602)',
    structure = SMILES('[CH][C][CH2]'),
    E0 = (981.672,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3100,440,815,1455,1000,180,917.882,918.218,918.236],'cm^-1')),
        HinderedRotor(inertia=(0.00702614,'amu*angstrom^2'), symmetry=1, barrier=(4.20306,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.182792,'amu*angstrom^2'), symmetry=1, barrier=(4.20276,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (39.0559,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.96474,0.0202917,-1.68078e-05,7.2211e-09,-1.21436e-12,118107,14.5039], Tmin=(100,'K'), Tmax=(1441.81,'K')), NASAPolynomial(coeffs=[8.09928,0.00604662,-1.9874e-06,3.68258e-10,-2.60955e-14,116627,-12.1458], Tmin=(1441.81,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(981.672,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][CH]C([O])[C]=O(18888)',
    structure = SMILES('[CH2][C][CH]C([O])[C]=O'),
    E0 = (836.292,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1855,455,950,180,867.033,867.036,3540.35],'cm^-1')),
        HinderedRotor(inertia=(0.0379554,'amu*angstrom^2'), symmetry=1, barrier=(20.248,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.159037,'amu*angstrom^2'), symmetry=1, barrier=(3.65657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.880682,'amu*angstrom^2'), symmetry=1, barrier=(20.2486,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00685144,'amu*angstrom^2'), symmetry=1, barrier=(3.65475,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.26372,0.0568439,-5.99674e-05,3.20757e-08,-6.74285e-12,100684,34.3033], Tmin=(100,'K'), Tmax=(1162.18,'K')), NASAPolynomial(coeffs=[13.4308,0.0149675,-5.91893e-06,1.07181e-09,-7.35605e-14,97856,-26.2245], Tmin=(1162.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(836.292,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsHH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsCsHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJCO) + radical(RCCJ) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
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
    label = 'C=[C][CH][CH][O](17610)',
    structure = SMILES('[CH2][C]=C[CH][O]'),
    E0 = (549.304,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,180,1441.42],'cm^-1')),
        HinderedRotor(inertia=(0.0637222,'amu*angstrom^2'), symmetry=1, barrier=(39.9635,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0270909,'amu*angstrom^2'), symmetry=1, barrier=(39.9484,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.50002,0.0283298,-1.09144e-05,-2.0554e-09,1.67008e-12,66123.7,20.3322], Tmin=(100,'K'), Tmax=(1235.14,'K')), NASAPolynomial(coeffs=[8.2763,0.0176346,-7.655e-06,1.43668e-09,-9.96449e-14,64085.8,-11.2287], Tmin=(1235.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(549.304,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=CCJO) + radical(Allyl_P) + radical(Cds_S)"""),
)

species(
    label = '[CH2][C][C]=C([O])C[O](18889)',
    structure = SMILES('[CH2][C][C]=C([O])C[O]'),
    E0 = (834.571,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,240.978,242.127,2859.06,2859.88,2860.64],'cm^-1')),
        HinderedRotor(inertia=(0.937036,'amu*angstrom^2'), symmetry=1, barrier=(38.2877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.311277,'amu*angstrom^2'), symmetry=1, barrier=(12.8679,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.929121,'amu*angstrom^2'), symmetry=1, barrier=(38.2876,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.15705,0.0705544,-0.000116876,1.05414e-07,-3.66339e-11,100470,31.1205], Tmin=(100,'K'), Tmax=(862.072,'K')), NASAPolynomial(coeffs=[6.67991,0.0279725,-1.32805e-05,2.48497e-09,-1.67958e-13,100148,8.94993], Tmin=(862.072,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(834.571,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(CCOJ) + radical(C=C(C)OJ) + radical(RCCJ) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH2][C][C]=C([O])[CH]O(18890)',
    structure = SMILES('[CH2][C][C]=C([O])[CH]O'),
    E0 = (726.161,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,350,440,435,1725,3025,407.5,1350,352.5,3000,3100,440,815,1455,1000,1685,370,499.193,499.193,499.2,499.201],'cm^-1')),
        HinderedRotor(inertia=(0.0134917,'amu*angstrom^2'), symmetry=1, barrier=(2.38582,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0134918,'amu*angstrom^2'), symmetry=1, barrier=(2.38584,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0826062,'amu*angstrom^2'), symmetry=1, barrier=(14.6077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.193794,'amu*angstrom^2'), symmetry=1, barrier=(34.2692,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.850566,0.0640389,-7.08076e-05,3.75641e-08,-7.38266e-12,87455,31.9565], Tmin=(100,'K'), Tmax=(1008.58,'K')), NASAPolynomial(coeffs=[15.9273,0.0116649,-3.95005e-06,6.65819e-10,-4.45892e-14,84036.4,-42.7798], Tmin=(1008.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(726.161,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(C=CCJO) + radical(RCCJ) + radical(Cds_S) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O][C](C=O)C1[C]C1(18891)',
    structure = SMILES('[O]C=C([O])C1[C]C1'),
    E0 = (418.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29917,0.0380306,3.09735e-05,-8.72036e-08,4.19528e-11,50445.6,24.0592], Tmin=(100,'K'), Tmax=(914.684,'K')), NASAPolynomial(coeffs=[23.279,-0.00200099,4.64215e-06,-9.72846e-10,6.12964e-14,44078.3,-92.8463], Tmin=(914.684,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(418.454,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopropane) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O][C][CH][O](10223)',
    structure = SMILES('[O][C][CH][O]'),
    E0 = (676.298,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3025,407.5,1350,352.5,180,180,1781.81,1786.65],'cm^-1')),
        HinderedRotor(inertia=(0.371605,'amu*angstrom^2'), symmetry=1, barrier=(8.54393,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (57.0281,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61268,0.044967,-0.000114903,1.26746e-07,-4.85652e-11,81376.1,15.5028], Tmin=(100,'K'), Tmax=(883.837,'K')), NASAPolynomial(coeffs=[-0.0833967,0.0206941,-1.18058e-05,2.28838e-09,-1.54199e-13,83277.3,36.2363], Tmin=(883.837,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(676.298,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(103.931,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsOsHH) + group(Cs-CsOsHH) + radical(CCOJ) + radical(CCOJ) + radical(CCsJOH) + radical(CH2_triplet)"""),
)

species(
    label = '[CH2]C=[C][C]([O])[CH][O](18892)',
    structure = SMILES('[CH2][CH][C]=C([O])[CH][O]'),
    E0 = (644.764,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3000,3050,390,425,1340,1360,335,370,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,656.32,656.336,656.398,656.412],'cm^-1')),
        HinderedRotor(inertia=(0.000391206,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.000391348,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00039132,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.68929,0.0503022,-4.48542e-05,2.11923e-08,-4.11815e-12,77630.8,28.7971], Tmin=(100,'K'), Tmax=(1214.37,'K')), NASAPolynomial(coeffs=[10.3725,0.0217008,-9.52538e-06,1.79739e-09,-1.25352e-13,75521.9,-14.7806], Tmin=(1214.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(644.764,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(Allyl_S) + radical(C=CCJO) + radical(RCCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=C[CH][C]([O])[CH][O](18893)',
    structure = SMILES('[CH][CH]C=C([O])[CH][O]'),
    E0 = (649.891,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3010,987.5,1337.5,450,1655,350,440,435,1725,3000,3050,390,425,1340,1360,335,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.34032,0.0514316,-4.08725e-05,1.43572e-08,-1.41417e-12,78265.6,28.1026], Tmin=(100,'K'), Tmax=(1119.3,'K')), NASAPolynomial(coeffs=[13.4504,0.0171144,-6.89087e-06,1.26917e-09,-8.82968e-14,74993.3,-34.1936], Tmin=(1119.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(649.891,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + radical(C=C(C)OJ) + radical(CCOJ) + radical(Allyl_S) + radical(C=CCJO) + radical(CCJ2_triplet)"""),
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
    E0 = (270.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1064.81,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1070.18,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (734.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (278.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (458.928,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (277.488,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (333.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (295.348,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (901.226,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (537.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (536.086,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (444.885,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (451.013,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (558.243,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (650.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (468.417,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (575.473,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (585.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (565.631,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (818.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (681.163,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (721.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (420.012,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (531.436,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (481.804,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (565.693,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (544.393,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (701.193,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (431.513,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (425.104,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (273.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (278.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (333.775,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (644.573,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (472.46,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (609.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (440.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (615.63,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (473.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (418.239,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (533.062,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (270.375,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (919.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (888.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (293.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (278.659,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (348.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (348.622,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (629.009,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (623.739,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (683.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (645.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (532.843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (628.842,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (408.43,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (757.846,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (431.034,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (841.835,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (319.435,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (754.911,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (386.817,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (273.304,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (298.228,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS65',
    E0 = (362.155,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS66',
    E0 = (747.119,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS67',
    E0 = (595.756,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS68',
    E0 = (961.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS69',
    E0 = (859.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS70',
    E0 = (559.125,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS71',
    E0 = (735.391,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS72',
    E0 = (897.972,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS73',
    E0 = (751.135,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS74',
    E0 = (418.454,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS75',
    E0 = (1070.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS76',
    E0 = (667.626,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS77',
    E0 = (672.753,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['OCHCO(3676)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH2][C]=C[C]=C[O](18003)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH]=C([O])[CH][C]=C(16203)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['H(3)', '[CH][C]=CC([O])=C[O](18837)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C1=CC(=C[O])O1(18838)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSD;O_rad;Ypri_rad_out]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2][C]=CC1=COO1(18839)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(188.553,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 184.8 to 188.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction7',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C1=CC([O])=CO1(18840)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2][C]=CC(O)=C=O(18841)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C=CC([O])=C=O(18842)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH2][C]=[C]C([O])[CH][O](18843)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2][C]=C[C]1OC1[O](18844)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(267.467,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2][C]=CC1([O])[CH]O1(18845)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(265.711,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction13',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C1=C[C]([O])C1[O](18846)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(174.51,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 168.3 to 174.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction14',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[O][C]1C=[C]CC1[O](18847)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(4.49459e+10,'s^-1'), n=0.314867, Ea=(180.638,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 177.3 to 180.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C1=CC1([O])[CH][O](18848)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1.69117e+12,'s^-1'), n=0.1675, Ea=(287.868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[O][CH]C1([O])C=[C]C1(18849)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.30019e+08,'s^-1'), n=1.00802, Ea=(380.23,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 379.2 to 380.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction17',
    reactants = ['H(3)', '[CH2]C#CC([O])=C[O](18850)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['H(3)', '[CH2][C]=CC([O])=C=O(18851)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction19',
    reactants = ['OCHCO(3676)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[O][C]=C[O](9592)', 'C3H3(5450)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(0.523563,'m^3/(mol*s)'), n=2.10494, Ea=(22.6844,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[O][C]=C[O](9592)', '[CH]=[C][CH2](16918)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH2][C]=[C]C([O])=C[O](18852)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction23',
    reactants = ['H(3)', '[CH2][C]=CC([O])=[C][O](18853)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction24',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C=[C]C([O])=C[O](18854)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(8.44274e+08,'s^-1'), n=1.26608, Ea=(149.637,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out;Cd_H_out_singleDe]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[CH2][C]=CC([O])=[C]O(18855)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction26',
    reactants = ['[CH2][C]=[C]C(O)=C[O](18856)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.87605e+06,'s^-1'), n=1.79171, Ea=(150.238,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;Y_rad_out;O_H_out] + [R3H_SS;Cd_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;Cd_rad_out;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction27',
    reactants = ['[CH2][C]=CC(O)=[C][O](18857)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction28',
    reactants = ['C[C]=[C]C([O])=C[O](18858)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(3.85113e+09,'s^-1'), n=1.0541, Ea=(193.078,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R3Hall;Cd_rad_out_singleDe_Cd;Cs_H_out_2H] for rate rule [R3HJ;Cd_rad_out_singleDe_Cd;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C=CC([O])=[C][O](18859)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(2.2886e+12,'s^-1'), n=1.11742, Ea=(430.818,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSD;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH2][C]=[C]C([O])=CO(18860)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.11e+08,'s^-1'), n=1.1915, Ea=(103.605,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_SDS;Cd_rad_out;XH_out] for rate rule [R4H_SDS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['C[C]=CC([O])=[C][O](18861)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(408000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out;Cs_H_out_2H] for rate rule [R5HJ_3;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction32',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=[C]C1OC1=C[O](18809)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=[C]C1OC=C1[O](18862)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction34',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=[C]CC([O])=C=O(16698)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=[C]C1[C]([O])C1[O](18863)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(2.61353e+11,'s^-1'), n=0.395207, Ea=(374.198,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_D;doublebond_intra;radadd_intra_cs]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 373.3 to 374.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['H(3)', 'C#C[CH]C([O])=C[O](18864)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(4278.27,'m^3/(mol*s)'), n=1.383, Ea=(15.1097,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Ct_Ct;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction37',
    reactants = ['O(4)', 'C=[C]C=C=C[O](18009)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;O_atom_triplet]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]=C[CH]C([O])=C[O](18865)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(3.3587e+07,'s^-1'), n=1.74167, Ea=(154.868,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_D;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction39',
    reactants = ['C=[C]CC([O])=[C][O](16701)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(70985.8,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=[C][CH]C(O)=C[O](18866)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=[C][CH]C([O])=CO(18867)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.86943e+06,'s^-1'), n=1.85754, Ea=(151.521,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;XH_out] for rate rule [R6Hall;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.4142135623730951
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=C=CC([O])=C[O](16697)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH][O](1548)', '[CH2][C]=C[C]=O(9897)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction45',
    reactants = ['H2CC(T)(1341)', '[CH]=C([O])[CH][O](10249)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=[C]C=C1OC1[O](18868)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(23.4597,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;O_rad;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination
Ea raised from 23.0 to 23.5 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction47',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=C1C=C([O])C1[O](18869)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction48',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=C=C=C(O)[CH][O](18870)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=C=C=C([O])C[O](18871)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(2.00399e+09,'s^-1'), n=0.37, Ea=(78.2471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3;Y_rad;XH_Rrad_De] + [R3radExo;Y_rad;XH_Rrad] for rate rule [R3radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction50',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=[C]C1O[C]1[CH][O](18872)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(358.634,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=[C]C1O[CH][C]1[O](18873)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(2.04143e+10,'s^-1'), n=0.464715, Ea=(353.364,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R4_linear;multiplebond_intra;radadd_intra_O] + [R4_linear;doublebond_intra;radadd_intra] for rate rule [R4_linear;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 353.1 to 353.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH][O](1548)', 'C=C=C[C]=O(9908)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Ck_O;YJ] for rate rule [Ck_O;Y_1centerbirad]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction53',
    reactants = ['C=[C][C]=C([O])C[O](18874)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(70985.8,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH]=[C]C=C([O])C[O](18875)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['C=[C][CH]C([O])[C]=O(16747)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction56',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C1=C[C]([CH][O])O1(18876)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction57',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2][C]=C[C]1[CH]OO1(18877)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(487.471,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_linear;carbonyl_intra_De;radadd_intra_O]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction58',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[O][CH][C]1C=[C]CO1(18878)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(1.19156e+09,'s^-1'), n=0.640131, Ea=(160.659,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R5_linear;carbonyl_intra;radadd_intra_cs2H]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 155.9 to 160.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction59',
    reactants = ['HCO(1372)', 'C=[C][CH][C][O](18879)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction60',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=C1[CH][C]([O])[CH]O1(18880)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(2.55236e+10,'s^-1'), n=0.405358, Ea=(49.0602,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cddouble] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_cddouble]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH]=[C][CH]C([O])C=O(18881)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(868165,'s^-1'), n=1.92199, Ea=(209.015,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4Hall;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['C=[C]C1C(=O)C1[O](18882)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(116.442,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination
Ea raised from 114.2 to 116.4 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction63',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C1=CC1([O])C=O(18883)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;Ypri_rad_out]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction64',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[CH2]C=C=C([O])C=O(18884)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.4733e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction65',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[O][C]1[CH]OC[C]=C1(18885)'],
    transitionState = 'TS65',
    kinetics = Arrhenius(A=(5.8912e+08,'s^-1'), n=0.529986, Ea=(91.7795,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cs2H] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic
Ea raised from 88.7 to 91.8 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction66',
    reactants = ['[CH2][C]=[C]C([O])C=O(18886)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS66',
    kinetics = Arrhenius(A=(3.64487e+09,'s^-1'), n=1.11169, Ea=(182.567,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction67',
    reactants = ['[CH2][C]C([C]=O)C=O(18887)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS67',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CO]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction68',
    reactants = ['OCHCO(3676)', '[CH][C][CH2](17602)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS68',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction69',
    reactants = ['[CH2][C][CH]C([O])[C]=O(18888)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS69',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction70',
    reactants = ['HCO(1372)', '[CH2][C]=C[C]=O(9897)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS70',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using an average for rate rule [Cd_R;CO_pri_rad]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction71',
    reactants = ['CO(2039)', 'C=[C][CH][CH][O](17610)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS71',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction72',
    reactants = ['[CH2][C][C]=C([O])C[O](18889)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS72',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction73',
    reactants = ['[CH2][C][C]=C([O])[CH]O(18890)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS73',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction74',
    reactants = ['C=[C][CH]C([O])=C[O](16700)'],
    products = ['[O][C](C=O)C1[C]C1(18891)'],
    transitionState = 'TS74',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(148.079,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 145.8 to 148.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction75',
    reactants = ['[O][C][CH][O](10223)', 'C3H3(5450)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS75',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [Cd_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction76',
    reactants = ['[CH2]C=[C][C]([O])[CH][O](18892)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS76',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction77',
    reactants = ['[CH]=C[CH][C]([O])[CH][O](18893)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS77',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

network(
    label = 'PDepNetwork #4078',
    isomers = [
        'C=[C][CH]C([O])=C[O](16700)',
    ],
    reactants = [
        ('OCHCO(3676)', 'C3H3(5450)'),
    ],
    bathGas = {
        'Ar': 0.333333,
        'Ne': 0.333333,
        'N2': 0.333333,
    },
)

pressureDependence(
    label = 'PDepNetwork #4078',
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


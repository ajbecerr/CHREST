species(
    label = '[CH]C(=C)C([O])=C[O](16771)',
    structure = SMILES('[CH]C(=C)C([O])=C[O]'),
    E0 = (283.597,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07555,'amu*angstrom^2'), symmetry=1, barrier=(47.721,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07358,'amu*angstrom^2'), symmetry=1, barrier=(47.6758,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.398941,0.0752525,-7.67689e-05,3.74801e-08,-6.79431e-12,34285.7,27.5192], Tmin=(100,'K'), Tmax=(1580.08,'K')), NASAPolynomial(coeffs=[20.6818,0.00836382,-4.33065e-07,-1.43713e-10,1.5485e-14,29311.9,-78.4855], Tmin=(1580.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(283.597,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C([CH2])=[C]C=O(18041)',
    structure = SMILES('[CH]C(=C)[C]=C[O]'),
    E0 = (558.68,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.07691,'amu*angstrom^2'), symmetry=1, barrier=(47.7523,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.07652,'amu*angstrom^2'), symmetry=1, barrier=(47.7433,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16763,0.0517738,-2.13603e-05,-2.06327e-08,1.54167e-11,67305.6,20.7387], Tmin=(100,'K'), Tmax=(905.023,'K')), NASAPolynomial(coeffs=[16.7502,0.0115593,-2.20424e-06,2.43448e-10,-1.47857e-14,63311.5,-59.3666], Tmin=(905.023,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(558.68,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(C=CJC=C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C(=[CH])[O](16209)',
    structure = SMILES('[CH]C(=C)C(=[CH])[O]'),
    E0 = (598.023,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3120,650,792.5,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.13889,'amu*angstrom^2'), symmetry=1, barrier=(49.1774,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13928,'amu*angstrom^2'), symmetry=1, barrier=(49.1862,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3675.45,'J/mol'), sigma=(6.08582,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=574.10 K, Pc=37 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.09494,0.0590546,-5.48217e-05,2.19943e-08,-1.83152e-12,72034.6,21.2792], Tmin=(100,'K'), Tmax=(925.341,'K')), NASAPolynomial(coeffs=[13.9432,0.0159923,-5.24267e-06,8.46817e-10,-5.50013e-14,69122.6,-42.5958], Tmin=(925.341,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(598.023,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C1=COO1(18952)',
    structure = SMILES('[CH]C(=C)C1=COO1'),
    E0 = (523.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.16871,0.0493845,-8.75914e-06,-2.92839e-08,1.60191e-11,63097.5,22.1793], Tmin=(100,'K'), Tmax=(971.726,'K')), NASAPolynomial(coeffs=[16.307,0.0174248,-6.28229e-06,1.16427e-09,-8.51031e-14,58722.3,-57.7941], Tmin=(971.726,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.678,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C(O)=C=O(18953)',
    structure = SMILES('[CH]C(=C)C(O)=C=O'),
    E0 = (196.439,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.241169,0.0800164,-9.95021e-05,6.23733e-08,-1.51571e-11,23763.8,24.2362], Tmin=(100,'K'), Tmax=(1017.43,'K')), NASAPolynomial(coeffs=[16.419,0.0164125,-5.72944e-06,9.28259e-10,-5.87364e-14,20471.9,-54.092], Tmin=(1017.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(196.439,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CCO)) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])C([O])[CH][O](18954)',
    structure = SMILES('[CH]C(=[CH])C([O])[CH][O]'),
    E0 = (865.933,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3025,407.5,1350,352.5,3120,650,792.5,1650,518.425,518.541,518.601,518.656,518.703,518.724],'cm^-1')),
        HinderedRotor(inertia=(0.276,'amu*angstrom^2'), symmetry=1, barrier=(52.6667,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275818,'amu*angstrom^2'), symmetry=1, barrier=(52.6651,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.276008,'amu*angstrom^2'), symmetry=1, barrier=(52.6608,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.737304,0.0780154,-0.000114801,9.68353e-08,-3.29113e-11,104259,29.554], Tmin=(100,'K'), Tmax=(814.231,'K')), NASAPolynomial(coeffs=[8.28244,0.0325598,-1.56061e-05,2.96396e-09,-2.03944e-13,103309,-3.5886], Tmin=(814.231,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(865.933,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C([CH2])C([O])=[C][O](18955)',
    structure = SMILES('[CH]C([CH2])C([O])=[C][O]'),
    E0 = (700.684,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,556.259,556.26,556.265,556.335,556.433,3423.16],'cm^-1')),
        HinderedRotor(inertia=(0.000544683,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0318015,'amu*angstrom^2'), symmetry=1, barrier=(6.98725,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.368779,'amu*angstrom^2'), symmetry=1, barrier=(80.9976,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.559384,0.0667814,-7.75958e-05,4.41052e-08,-9.42767e-12,84404.4,33.1631], Tmin=(100,'K'), Tmax=(1288.14,'K')), NASAPolynomial(coeffs=[16.9481,0.00844316,-9.9065e-07,-2.935e-11,8.92552e-15,80800,-47.6542], Tmin=(1288.14,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(700.684,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Isobutyl) + radical(CCJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(=[CH])[C]([O])C[O](18956)',
    structure = SMILES('[CH]C([CH])=C([O])C[O]'),
    E0 = (668.228,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,325,375,415,465,420,450,1700,1750,200,800,914.286,1028.57,1142.86,1257.14,1371.43,1485.71,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.10131,0.0690268,-8.5338e-05,7.00182e-08,-2.42088e-11,80468.5,28.5263], Tmin=(100,'K'), Tmax=(824.067,'K')), NASAPolynomial(coeffs=[5.21177,0.0406444,-1.833e-05,3.39484e-09,-2.30934e-13,80077.3,11.2279], Tmin=(824.067,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(668.228,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(CCOJ) + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C](C)C([O])=[C][O](18957)',
    structure = SMILES('[CH][C](C)C([O])=[C][O]'),
    E0 = (648.176,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,360,370,350,350,440,435,1725,1685,370,180,180,968.304,968.33,968.345,968.362],'cm^-1')),
        HinderedRotor(inertia=(0.0241452,'amu*angstrom^2'), symmetry=1, barrier=(16.0657,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00165342,'amu*angstrom^2'), symmetry=1, barrier=(16.0658,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.385373,'amu*angstrom^2'), symmetry=1, barrier=(85.4985,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.613259,0.0605988,-6.03617e-05,2.92418e-08,-5.35173e-12,78091.3,31.3972], Tmin=(100,'K'), Tmax=(1504.64,'K')), NASAPolynomial(coeffs=[17.4199,0.00889146,-1.80767e-06,1.93726e-10,-9.53214e-15,73829.3,-53.9077], Tmin=(1504.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(648.176,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CCJ(C)CO) + radical(CCJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(=C)[C]1OC1[O](18958)',
    structure = SMILES('[CH]C([CH2])=C1OC1[O]'),
    E0 = (490.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.989069,0.0565616,-3.38716e-05,-5.01051e-10,5.14648e-12,59070.8,22.5078], Tmin=(100,'K'), Tmax=(1010.28,'K')), NASAPolynomial(coeffs=[15.4388,0.0190612,-7.45847e-06,1.38091e-09,-9.79999e-14,55145.2,-52.33], Tmin=(1010.28,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(490.17,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-CsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(methyleneoxirane) + radical(CCOJ) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C1([O])[CH]O1(18959)',
    structure = SMILES('[CH]C(=C)C1([O])[CH]O1'),
    E0 = (535.668,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.808684,0.0799357,-8.81034e-05,4.69431e-08,-9.02444e-12,64621.4,28.8354], Tmin=(100,'K'), Tmax=(1574,'K')), NASAPolynomial(coeffs=[18.5297,0.0075647,2.99969e-06,-1.01866e-09,8.24417e-14,61410.8,-64.0939], Tmin=(1574,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(535.668,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(C=CC(C)(O)OJ) + radical(CCsJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]1COC1=C[O](18960)',
    structure = SMILES('[CH][C]1COC1=C[O]'),
    E0 = (446.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.38309,0.0312213,6.27485e-05,-1.29821e-07,6.03587e-11,53816.9,20.1662], Tmin=(100,'K'), Tmax=(896.079,'K')), NASAPolynomial(coeffs=[25.8841,-0.00747347,9.21458e-06,-1.97443e-09,1.33761e-13,46588.5,-111.182], Tmin=(896.079,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(446.469,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ) + radical(CCJ(C)CO) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH][C]1COC=C1[O](18961)',
    structure = SMILES('[CH]C1=C([O])[CH]OC1'),
    E0 = (286.23,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.79276,0.0381504,6.86915e-06,-3.37927e-08,1.49975e-11,34514.1,22.4727], Tmin=(100,'K'), Tmax=(992.87,'K')), NASAPolynomial(coeffs=[11.085,0.0248734,-9.57073e-06,1.75297e-09,-1.23503e-13,31478.2,-28.2872], Tmin=(992.87,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(286.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(25dihydrofuran) + radical(C=C(C)OJ) + radical(C=CCJ(O)C) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH2])OC1=C[O](18909)',
    structure = SMILES('[CH]C1([CH2])OC1=C[O]'),
    E0 = (524.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-1.77173,0.0847022,-9.64378e-05,4.79631e-08,-8.41494e-12,63332.6,30.5169], Tmin=(100,'K'), Tmax=(1742.29,'K')), NASAPolynomial(coeffs=[24.4843,-0.00557259,7.10715e-06,-1.53844e-09,1.05866e-13,58736.1,-97.6653], Tmin=(1742.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(524.544,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(methyleneoxirane) + radical(C=COJ) + radical(CJC(C)OC) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH2])OC=C1[O](18962)',
    structure = SMILES('[CH]C1([CH2])OC=C1[O]'),
    E0 = (495.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.412592,0.0577864,-7.20399e-06,-6.30009e-08,3.77844e-11,59762.6,23.4995], Tmin=(100,'K'), Tmax=(893.301,'K')), NASAPolynomial(coeffs=[28.8504,-0.0108806,9.58034e-06,-2.00269e-09,1.36788e-13,52341,-123.591], Tmin=(893.301,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(495.654,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclobutene) + radical(C=C(C)OJ) + radical(CJC(C)OC) + radical(CCJ2_triplet)"""),
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
    label = '[CH]C(=C)C([O])=C=O(18963)',
    structure = SMILES('[CH]C(=C)C(=O)[C]=O'),
    E0 = (304.742,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,375,552.5,462.5,1710,2950,3100,1380,975,1025,1650,1855,455,950,408.552,408.552,408.552,408.552],'cm^-1')),
        HinderedRotor(inertia=(0.436804,'amu*angstrom^2'), symmetry=1, barrier=(51.7379,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.436805,'amu*angstrom^2'), symmetry=1, barrier=(51.7378,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.436805,'amu*angstrom^2'), symmetry=1, barrier=(51.7378,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.61614,0.0566059,-6.24457e-05,4.24515e-08,-1.25417e-11,36734.1,24.2258], Tmin=(100,'K'), Tmax=(802.035,'K')), NASAPolynomial(coeffs=[6.72315,0.0311365,-1.48133e-05,2.85981e-09,-2.01079e-13,35914.8,0.713814], Tmin=(802.035,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.742,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(AllylJ2_triplet) + radical(CCCJ=O)"""),
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
    label = '[CH]C(=[CH])C([O])=C[O](18964)',
    structure = SMILES('[CH]C(=[CH])C([O])=C[O]'),
    E0 = (530.693,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.05476,'amu*angstrom^2'), symmetry=1, barrier=(47.2429,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.06032,'amu*angstrom^2'), symmetry=1, barrier=(47.3708,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.223611,0.0760184,-8.45626e-05,4.43681e-08,-8.62758e-12,63994.2,27.2805], Tmin=(100,'K'), Tmax=(1456.43,'K')), NASAPolynomial(coeffs=[20.9718,0.00562275,4.87115e-07,-3.06421e-10,2.68495e-14,59112.5,-78.5081], Tmin=(1456.43,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(530.693,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]C(=C)C([O])=[C][O](18965)',
    structure = SMILES('[CH]C(=C)C([O])=[C][O]'),
    E0 = (523.341,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,1685,370,213.58,214.919,217.421,218.569,219.369,220.332],'cm^-1')),
        HinderedRotor(inertia=(1.52998,'amu*angstrom^2'), symmetry=1, barrier=(49.8503,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.51976,'amu*angstrom^2'), symmetry=1, barrier=(49.7881,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.495334,0.0691539,-7.70586e-05,4.28318e-08,-9.13792e-12,63076.4,27.4233], Tmin=(100,'K'), Tmax=(1226.41,'K')), NASAPolynomial(coeffs=[16.8258,0.0123969,-3.36631e-06,4.50104e-10,-2.49757e-14,59333.6,-53.6234], Tmin=(1226.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(523.341,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(AllylJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH]C(=C)C([O])=[C]O(18966)',
    structure = SMILES('[CH]C(=C)C([O])=[C]O'),
    E0 = (381.878,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,1685,370,180,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.09095,'amu*angstrom^2'), symmetry=1, barrier=(48.075,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09028,'amu*angstrom^2'), symmetry=1, barrier=(48.0596,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.09219,'amu*angstrom^2'), symmetry=1, barrier=(48.1037,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.152108,0.0789622,-9.0676e-05,5.03262e-08,-1.04641e-11,46089.5,28.3332], Tmin=(100,'K'), Tmax=(1335.38,'K')), NASAPolynomial(coeffs=[19.9766,0.00878738,-7.50825e-07,-1.08759e-10,1.53181e-14,41594.7,-71.2988], Tmin=(1335.38,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(381.878,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=CJO) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C(O)=[C][O](18967)',
    structure = SMILES('[CH]C([CH2])=C(O)[C]=O'),
    E0 = (319.829,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,1855,455,950,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.97194,'amu*angstrom^2'), symmetry=1, barrier=(45.3388,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.96981,'amu*angstrom^2'), symmetry=1, barrier=(45.2899,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97082,'amu*angstrom^2'), symmetry=1, barrier=(45.3131,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97308,'amu*angstrom^2'), symmetry=1, barrier=(45.365,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.126875,0.083107,-9.97323e-05,5.8764e-08,-1.34746e-11,38607.8,24.4859], Tmin=(100,'K'), Tmax=(1071.96,'K')), NASAPolynomial(coeffs=[17.7904,0.0171969,-7.5052e-06,1.40749e-09,-9.82642e-14,34820.9,-61.958], Tmin=(1071.96,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(319.829,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(C=CCJ=O)"""),
)

species(
    label = '[CH]C(=[CH])C(O)=C[O](18968)',
    structure = SMILES('[CH]C([CH])=C(O)C=O'),
    E0 = (378.392,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,2782.5,750,1395,475,1775,1000,200,800,1000,1200,1400,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.488382,0.0738439,-6.99079e-05,3.39406e-08,-6.6373e-12,45639.3,24.4902], Tmin=(100,'K'), Tmax=(1225.2,'K')), NASAPolynomial(coeffs=[15.1115,0.0261031,-1.14596e-05,2.13734e-09,-1.47926e-13,42056.1,-49.0277], Tmin=(1225.2,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(378.392,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])C([O])=CO(18969)',
    structure = SMILES('[CH]C(=[CH])C([O])=CO'),
    E0 = (389.23,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,3010,987.5,1337.5,450,1655,3120,650,792.5,1650,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.97257,'amu*angstrom^2'), symmetry=1, barrier=(45.3532,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97006,'amu*angstrom^2'), symmetry=1, barrier=(45.2956,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.97117,'amu*angstrom^2'), symmetry=1, barrier=(45.3211,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.887945,0.0860132,-9.87795e-05,5.25718e-08,-1.02273e-11,47008.1,27.5588], Tmin=(100,'K'), Tmax=(1488.62,'K')), NASAPolynomial(coeffs=[23.5253,0.00287882,2.6593e-06,-7.69512e-10,5.97264e-14,41682.5,-93.4081], Tmin=(1488.62,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(389.23,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4284.46,'J/mol'), sigma=(6.81655,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=669.22 K, Pc=30.69 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.293993,0.0664589,-7.20796e-05,3.68989e-08,-7.0232e-12,51815.7,29.8339], Tmin=(100,'K'), Tmax=(1472.35,'K')), NASAPolynomial(coeffs=[19.3893,0.00508522,1.21664e-07,-1.73881e-10,1.54663e-14,47222.1,-66.1816], Tmin=(1472.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(429.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S) + radical(Cds_P)"""),
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
    label = '[CH]=[C]C([O])=C[O](13159)',
    structure = SMILES('[CH]=[C]C([O])=C[O]'),
    E0 = (396.595,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.55203,'amu*angstrom^2'), symmetry=1, barrier=(35.6842,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.690677,0.0610483,-7.87226e-05,4.57036e-08,-9.59389e-12,47828.7,21.8755], Tmin=(100,'K'), Tmax=(1385.71,'K')), NASAPolynomial(coeffs=[17.8401,-0.00178865,3.73025e-06,-9.08451e-10,6.82983e-14,44356,-61.8356], Tmin=(1385.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(396.595,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJC=C) + radical(Cds_P)"""),
)

species(
    label = '[C]=C([CH2])C([O])=C[O](18970)',
    structure = SMILES('[C]C(=C)C([O])=C[O]'),
    E0 = (582.39,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,2950,3100,1380,975,1025,1650,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.27005,'amu*angstrom^2'), symmetry=1, barrier=(29.201,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.252765,0.0753064,-9.07081e-05,4.90494e-08,-9.63595e-12,70214.2,24.3245], Tmin=(100,'K'), Tmax=(1467.18,'K')), NASAPolynomial(coeffs=[22.6042,-0.00206551,3.78808e-06,-8.83039e-10,6.41416e-14,65127.7,-89.1863], Tmin=(1467.18,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(582.39,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(228.648,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(CJ3)"""),
)

species(
    label = '[O]C=C([O])C1=CC1(18971)',
    structure = SMILES('O=CC(=O)[C]1[CH]C1'),
    E0 = (178.47,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44555,0.0494863,-4.17741e-05,1.83297e-08,-3.21403e-12,21562.6,23.2511], Tmin=(100,'K'), Tmax=(1372.39,'K')), NASAPolynomial(coeffs=[12.5952,0.0169897,-6.25626e-06,1.0764e-09,-7.113e-14,18502.2,-34.0689], Tmin=(1372.39,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.47,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsCsHH) + group(Cs-CsCsHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + ring(Cyclopropane) + radical(C2CJCHO) + radical(CCJCC=O)"""),
)

species(
    label = '[CH]=C1COC1=C[O](18896)',
    structure = SMILES('[CH]C1=C(C=O)OC1'),
    E0 = (215.196,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14977,0.0489404,-7.28608e-06,-2.95567e-08,1.53198e-11,25996.9,20.0278], Tmin=(100,'K'), Tmax=(1003.95,'K')), NASAPolynomial(coeffs=[16.7795,0.0178482,-7.41896e-06,1.46802e-09,-1.09629e-13,21287.2,-63.2642], Tmin=(1003.95,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(215.196,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-Cds(Cds-O2d)O2s) + group(Cds-O2d(Cds-Cds)H) + ring(Cyclobutene) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1=COC1=C[O](18820)',
    structure = SMILES('C=C1[CH]OC1=C[O]'),
    E0 = (86.5843,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.4976,0.0364381,2.69423e-05,-7.27464e-08,3.35647e-11,10520.9,19.1834], Tmin=(100,'K'), Tmax=(939.412,'K')), NASAPolynomial(coeffs=[19.5705,0.00656017,-5.18852e-07,8.62966e-11,-1.40013e-14,5048.05,-77.9343], Tmin=(939.412,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(86.5843,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(12methylenecyclobutane) + radical(C=COJ) + radical(C=CCJ(O)C)"""),
)

species(
    label = '[CH]=C1COC=C1[O](18972)',
    structure = SMILES('[CH]=C1COC=C1[O]'),
    E0 = (130.346,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.52146,0.0265869,7.37781e-05,-1.40075e-07,6.33967e-11,15792.4,18.9826], Tmin=(100,'K'), Tmax=(903.718,'K')), NASAPolynomial(coeffs=[26.1579,-0.00849898,9.25628e-06,-1.91983e-09,1.26728e-13,8319.42,-114.089], Tmin=(903.718,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(130.346,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=COC=C1[O](18973)',
    structure = SMILES('[CH2]C1=COC=C1[O]'),
    E0 = (-9.37808,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.71592,0.034562,2.03754e-05,-5.93962e-08,2.73725e-11,-1031.36,19.8155], Tmin=(100,'K'), Tmax=(945.934,'K')), NASAPolynomial(coeffs=[16.9374,0.0094847,-2.15994e-06,3.94215e-10,-3.39223e-14,-5668.81,-62.0643], Tmin=(945.934,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-9.37808,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsOsH) + ring(Furan) + radical(C=C(C)OJ) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C(C)C([O])=C=O(18974)',
    structure = SMILES('[CH]=C(C)C(=O)[C]=O'),
    E0 = (181.153,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39092,0.0621832,-8.74237e-05,7.23995e-08,-2.47188e-11,21877.1,24.3166], Tmin=(100,'K'), Tmax=(763.371,'K')), NASAPolynomial(coeffs=[7.4522,0.0274469,-1.33208e-05,2.57766e-09,-1.80244e-13,21038.4,-2.72108], Tmin=(763.371,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(181.153,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH2]C(=C)C(=O)[C]=O(16758)',
    structure = SMILES('[CH2]C(=C)C(=O)[C]=O'),
    E0 = (90.4088,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1855,455,950,522.18],'cm^-1')),
        HinderedRotor(inertia=(0.131867,'amu*angstrom^2'), symmetry=1, barrier=(3.03189,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.738432,'amu*angstrom^2'), symmetry=1, barrier=(16.978,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.243037,'amu*angstrom^2'), symmetry=1, barrier=(47.0259,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.53686,0.0553531,-5.64583e-05,3.06939e-08,-6.81885e-12,10961.4,22.1997], Tmin=(100,'K'), Tmax=(1076.3,'K')), NASAPolynomial(coeffs=[10.5056,0.0220212,-1.00047e-05,1.92017e-09,-1.35346e-13,9030.82,-21.7286], Tmin=(1076.3,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(90.4088,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(C=C(C=O)CJ) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=C1CC([O])[C]1[O](18975)',
    structure = SMILES('[CH]C1=C([O])C([O])C1'),
    E0 = (456.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.39759,0.0416616,1.42606e-05,-5.57411e-08,2.6393e-11,55060,24.4469], Tmin=(100,'K'), Tmax=(945.04,'K')), NASAPolynomial(coeffs=[17.0329,0.0143498,-4.07946e-06,7.15269e-10,-5.43126e-14,50369.2,-59.2833], Tmin=(945.04,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(456.897,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH2]C1=CC([O])[C]1[O](18976)',
    structure = SMILES('[CH2]C1=C([O])C([O])[CH]1'),
    E0 = (354.629,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.17716,0.0432157,1.46481e-05,-6.47066e-08,3.17997e-11,42771,21.6092], Tmin=(100,'K'), Tmax=(936.727,'K')), NASAPolynomial(coeffs=[21.5389,0.00414146,5.56063e-07,-1.1684e-10,5.21151e-17,36855.9,-86.5051], Tmin=(936.727,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(354.629,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutene) + radical(CC(C)OJ) + radical(C=C(C)OJ) + radical(C=CCJCO) + radical(Allyl_P)"""),
)

species(
    label = '[CH]=C1CC1([O])[CH][O](18936)',
    structure = SMILES('[CH]=C1CC1([O])[CH][O]'),
    E0 = (687.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.42235,0.0601808,-7.25758e-05,4.92309e-08,-1.37313e-11,82800,25.272], Tmin=(100,'K'), Tmax=(866.057,'K')), NASAPolynomial(coeffs=[9.06092,0.0249011,-1.14719e-05,2.19494e-09,-1.53669e-13,81476.9,-10.4813], Tmin=(866.057,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(687.692,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CC(C)2OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_P)"""),
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
    label = 'C#CC([O])=C[O](13171)',
    structure = SMILES('C#CC([O])=C[O]'),
    E0 = (127.289,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.58154,'amu*angstrom^2'), symmetry=1, barrier=(36.3626,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (82.0575,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.605988,0.0547459,-6.14813e-05,3.08789e-08,-5.58321e-12,15448.7,21.0975], Tmin=(100,'K'), Tmax=(1636.22,'K')), NASAPolynomial(coeffs=[17.7517,-0.00135695,2.95693e-06,-6.75319e-10,4.76604e-14,11737,-64.2598], Tmin=(1636.22,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(127.289,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=C(C)C([O])=[C][O](18977)',
    structure = SMILES('[CH]=C(C)C([O])=[C][O]'),
    E0 = (399.753,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2800,2850,1350,1500,750,1050,1375,1000,325,375,415,465,420,450,1700,1750,3120,650,792.5,1650,1685,370,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.735595,'amu*angstrom^2'), symmetry=1, barrier=(16.9128,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.735223,'amu*angstrom^2'), symmetry=1, barrier=(16.9042,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.684925,0.0696037,-8.29096e-05,4.60473e-08,-8.92883e-12,48201.6,26.0411], Tmin=(100,'K'), Tmax=(900.476,'K')), NASAPolynomial(coeffs=[16.427,0.010658,-3.01241e-06,4.39261e-10,-2.67919e-14,44921.3,-50.7274], Tmin=(900.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(399.753,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_P) + radical(C=CJO)"""),
)

species(
    label = '[CH2]C(=C)[C]([O])[C]=O(16759)',
    structure = SMILES('[CH2]C(=C)C([O])=[C][O]'),
    E0 = (304.156,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([325,375,415,465,420,450,1700,1750,3000,3100,440,815,1455,1000,2950,3100,1380,975,1025,1650,1685,370,266.213,266.498,266.614],'cm^-1')),
        HinderedRotor(inertia=(0.480299,'amu*angstrom^2'), symmetry=1, barrier=(24.1364,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.479592,'amu*angstrom^2'), symmetry=1, barrier=(24.1349,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.271561,0.0693497,-7.79531e-05,4.17181e-08,-8.34037e-12,36726.6,27.3996], Tmin=(100,'K'), Tmax=(1384.77,'K')), NASAPolynomial(coeffs=[19.2444,0.00606271,-2.11184e-07,-1.32601e-10,1.37851e-14,32285.3,-67.3729], Tmin=(1384.77,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(304.156,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)Cs) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Allyl_P) + radical(C=CJO)"""),
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
    label = '[CH]C(=C)[CH][O](17644)',
    structure = SMILES('[CH]C([CH2])=C[O]'),
    E0 = (421.922,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,180,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(2.1441,'amu*angstrom^2'), symmetry=1, barrier=(49.297,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(2.13703,'amu*angstrom^2'), symmetry=1, barrier=(49.1344,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (68.074,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.87061,0.034132,1.23923e-05,-4.75888e-08,2.28311e-11,50833.9,18.2246], Tmin=(100,'K'), Tmax=(931.786,'K')), NASAPolynomial(coeffs=[14.7241,0.0112454,-2.74727e-06,4.35245e-10,-3.25379e-14,47036.8,-50.4], Tmin=(931.786,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(421.922,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH][C]([CH2])C([O])[C]=O(18978)',
    structure = SMILES('[CH][C]([CH2])C([O])[C]=O'),
    E0 = (779.758,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,360,370,350,3000,3100,440,815,1455,1000,1855,455,950,217.42,926.693,1205.21,1629.54,2000.85],'cm^-1')),
        HinderedRotor(inertia=(0.0969767,'amu*angstrom^2'), symmetry=1, barrier=(2.99182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0969767,'amu*angstrom^2'), symmetry=1, barrier=(2.99182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0969767,'amu*angstrom^2'), symmetry=1, barrier=(2.99182,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0969767,'amu*angstrom^2'), symmetry=1, barrier=(2.99182,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.14234,0.0560994,-5.70999e-05,2.9878e-08,-6.09404e-12,93891.9,33.2508], Tmin=(100,'K'), Tmax=(1233.29,'K')), NASAPolynomial(coeffs=[13.9653,0.0135504,-4.18217e-06,6.41921e-10,-3.97019e-14,90802,-31.0053], Tmin=(1233.29,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(779.758,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-CsCsCsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + radical(C=OCOJ) + radical(CCJ(C)CO) + radical(Isobutyl) + radical(CCCJ=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C([CH2])=C1[CH]OO1(18979)',
    structure = SMILES('[CH]C([CH2])=C1[CH]OO1'),
    E0 = (655.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.76095,0.0293841,5.25832e-05,-9.11586e-08,3.68707e-11,78908.5,24.7624], Tmin=(100,'K'), Tmax=(975.252,'K')), NASAPolynomial(coeffs=[16.2047,0.0192598,-7.38968e-06,1.47896e-09,-1.14283e-13,73755.4,-56.5337], Tmin=(975.252,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(Cyclobutane) + radical(C=CCJO) + radical(Allyl_P) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH](6993)',
    structure = SMILES('[CH]'),
    E0 = (585.879,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([4000],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (13.0186,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[4.17282,-0.00431623,8.93613e-06,-7.67944e-09,2.45632e-12,70459.2,-0.856057], Tmin=(100,'K'), Tmax=(949.528,'K')), NASAPolynomial(coeffs=[3.66649,-0.000619273,6.25229e-07,-1.09644e-10,6.55478e-15,70484.8,1.18918], Tmin=(949.528,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(585.879,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(108.088,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(CJ3)"""),
)

species(
    label = 'C=C=C([O])C=O(11930)',
    structure = SMILES('C=[C]C(=O)C=O'),
    E0 = (61.4013,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,2782.5,750,1395,475,1775,1000,2950,3100,1380,975,1025,1650,1685,370,413.184],'cm^-1')),
        HinderedRotor(inertia=(0.569432,'amu*angstrom^2'), symmetry=1, barrier=(13.0924,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.03685,'amu*angstrom^2'), symmetry=1, barrier=(23.8392,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (83.0654,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.34268,0.0408215,-4.11136e-05,2.37985e-08,-6.04675e-12,7440.87,19.0438], Tmin=(100,'K'), Tmax=(910.792,'K')), NASAPolynomial(coeffs=[6.38375,0.0230738,-1.18843e-05,2.40342e-09,-1.74035e-13,6704.76,-0.074362], Tmin=(910.792,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(61.4013,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(199.547,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cds-O2d(Cds-O2d)Cs) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + radical(C=CJC=O)"""),
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
    label = '[CH]C(=C)[C]=O(17257)',
    structure = SMILES('[CH]C(=C)[C]=O'),
    E0 = (409.606,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,194.11,195.908,198.881],'cm^-1')),
        HinderedRotor(inertia=(1.8182,'amu*angstrom^2'), symmetry=1, barrier=(49.5014,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.77891,'amu*angstrom^2'), symmetry=1, barrier=(49.4512,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.84357,0.0324023,-1.85494e-05,4.5392e-09,-4.31472e-13,49298.1,15.5735], Tmin=(100,'K'), Tmax=(2261.7,'K')), NASAPolynomial(coeffs=[11.3354,0.0173839,-8.58903e-06,1.60324e-09,-1.06943e-13,45456.9,-32.3248], Tmin=(2261.7,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(409.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-CdsHH) + group(Cds-O2d(Cds-Cds)H) + radical(AllylJ2_triplet) + radical(C=C(C)CJ=O)"""),
)

species(
    label = '[CH]C([CH2])=[C][O](18980)',
    structure = SMILES('[CH]C([CH2])=[C][O]'),
    E0 = (661.666,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,1685,370,372.391,372.391,372.393,372.393],'cm^-1')),
        HinderedRotor(inertia=(0.514067,'amu*angstrom^2'), symmetry=1, barrier=(50.5877,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.514067,'amu*angstrom^2'), symmetry=1, barrier=(50.5877,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.07109,0.0360035,-1.49033e-05,-8.24232e-09,6.35593e-12,79655.3,20.6334], Tmin=(100,'K'), Tmax=(968.46,'K')), NASAPolynomial(coeffs=[10.6525,0.0152551,-5.5274e-06,9.71225e-10,-6.69592e-14,77304,-24.0496], Tmin=(968.46,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(661.666,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsOsH) + radical(C=COJ) + radical(Allyl_P) + radical(AllylJ2_triplet) + radical(C=CJO)"""),
)

species(
    label = '[CH]=C([CH2])C(=O)C=O(16881)',
    structure = SMILES('[CH]C(=C)C(=O)C=O'),
    E0 = (144.781,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.90553,0.0511228,-3.87171e-05,1.69496e-08,-3.45998e-12,17484.3,22.3662], Tmin=(100,'K'), Tmax=(1055.5,'K')), NASAPolynomial(coeffs=[5.90671,0.0359596,-1.71683e-05,3.33911e-09,-2.36273e-13,16639.7,2.84674], Tmin=(1055.5,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(144.781,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C(=C)C([O])[C]=O(16751)',
    structure = SMILES('[CH]C(=C)C([O])[C]=O'),
    E0 = (464.986,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,2950,3100,1380,975,1025,1650,1855,455,950,528.731,528.731,528.732,528.733,528.735],'cm^-1')),
        HinderedRotor(inertia=(0.267008,'amu*angstrom^2'), symmetry=1, barrier=(52.9694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267009,'amu*angstrom^2'), symmetry=1, barrier=(52.9694,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.267011,'amu*angstrom^2'), symmetry=1, barrier=(52.9694,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.32697,0.0569967,-5.13084e-05,2.4857e-08,-4.94326e-12,56022.6,29.3153], Tmin=(100,'K'), Tmax=(1194.6,'K')), NASAPolynomial(coeffs=[11.1383,0.0241443,-1.00573e-05,1.83608e-09,-1.25558e-13,53678.4,-19.7634], Tmin=(1194.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(464.986,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(CCCJ=O) + radical(AllylJ2_triplet)"""),
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
    label = '[CH]C(=[CH])[C](O)[CH][O](18981)',
    structure = SMILES('[CH]C([CH])=C(O)[CH][O]'),
    E0 = (647.718,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,325,375,415,465,420,450,1700,1750,3025,407.5,1350,352.5,200,800,933.333,1066.67,1200,1333.33,1466.67,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.660141,0.0652435,-4.90443e-05,1.60486e-08,-1.18266e-12,78029.7,29.2062], Tmin=(100,'K'), Tmax=(1084.07,'K')), NASAPolynomial(coeffs=[14.3836,0.025918,-1.02821e-05,1.83636e-09,-1.24826e-13,74389.6,-41.1751], Tmin=(1084.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(647.718,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + radical(CCOJ) + radical(C=CCJO) + radical(AllylJ2_triplet) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]C1([CH2])C(=O)C1[O](18982)',
    structure = SMILES('[CH]C1([CH2])C(=O)C1[O]'),
    E0 = (668.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.136951,0.0667228,-7.02407e-05,3.48232e-08,-6.3853e-12,80552.3,28.7634], Tmin=(100,'K'), Tmax=(1554.75,'K')), NASAPolynomial(coeffs=[19.6154,0.00506713,3.79625e-07,-2.33331e-10,1.95166e-14,75890.5,-69.3187], Tmin=(1554.75,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(668.463,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsCs) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsCs) + ring(cyclopropanone) + radical(C=OCOJ) + radical(CJC(C)2C=O) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]C(=[CH])C(=O)C[O](18983)',
    structure = SMILES('[CH]C(=[CH])C(=O)C[O]'),
    E0 = (546.423,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,350,440,435,1725,3120,650,792.5,1650,485.162,485.17,485.171,485.177,485.184],'cm^-1')),
        HinderedRotor(inertia=(0.00071616,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314083,'amu*angstrom^2'), symmetry=1, barrier=(52.4637,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.314081,'amu*angstrom^2'), symmetry=1, barrier=(52.4636,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.29419,0.0651612,-8.73176e-05,7.49762e-08,-2.71099e-11,65811.5,26.0415], Tmin=(100,'K'), Tmax=(769.991,'K')), NASAPolynomial(coeffs=[5.63221,0.035992,-1.75708e-05,3.39979e-09,-2.37736e-13,65340.1,7.52393], Tmin=(769.991,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(546.423,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC([O])C1=O(18984)',
    structure = SMILES('[CH]=C1CC([O])C1=O'),
    E0 = (308.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.27568,0.0246242,3.0957e-05,-5.57908e-08,2.19444e-11,37212.1,24.2737], Tmin=(100,'K'), Tmax=(1010.32,'K')), NASAPolynomial(coeffs=[12.0723,0.0182881,-7.81474e-06,1.58434e-09,-1.19624e-13,33576.4,-31.2864], Tmin=(1010.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(308.786,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=OCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH2]C1=CC([O])C1=O(18985)',
    structure = SMILES('C=C1[CH]C([O])C1=O'),
    E0 = (178.606,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13544,0.0250091,4.00147e-05,-6.9498e-08,2.75271e-11,21562.7,21.5728], Tmin=(100,'K'), Tmax=(998.736,'K')), NASAPolynomial(coeffs=[13.9047,0.0171961,-7.31126e-06,1.51581e-09,-1.17146e-13,17250.6,-45.0108], Tmin=(998.736,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(178.606,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cd-CdCs(CO)) + group(Cds-O2d(Cds-Cds)Cs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=OCOJ) + radical(C=CCJCO)"""),
)

species(
    label = '[CH2]C1=CO[C]1[CH][O](18986)',
    structure = SMILES('[CH2][C]1[CH]OC1=C[O]'),
    E0 = (397.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.0898,0.0391414,4.36367e-05,-1.16074e-07,5.78522e-11,47908.3,22.1667], Tmin=(100,'K'), Tmax=(879.908,'K')), NASAPolynomial(coeffs=[27.5,-0.0116151,1.20207e-05,-2.60955e-09,1.82668e-13,40577.8,-117.114], Tmin=(879.908,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(397.264,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CsCsH) + group(Cs-CsOsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ) + radical(CCJ(C)CO) + radical(CCsJOC(O)) + radical(Isobutyl)"""),
)

species(
    label = '[CH]C([CH2])([C]=O)C=O(18987)',
    structure = SMILES('[CH]C([CH2])([C]=O)C=O'),
    E0 = (472.432,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,400,1190,1240,450,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000,1855,455,950,298.576,298.692,298.865],'cm^-1')),
        HinderedRotor(inertia=(0.114876,'amu*angstrom^2'), symmetry=1, barrier=(7.25398,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114734,'amu*angstrom^2'), symmetry=1, barrier=(7.25556,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.114711,'amu*angstrom^2'), symmetry=1, barrier=(7.25455,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.344133,'amu*angstrom^2'), symmetry=1, barrier=(21.7942,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.505332,0.0818914,-0.000129185,1.04343e-07,-3.2882e-11,56941.5,27.8598], Tmin=(100,'K'), Tmax=(843.23,'K')), NASAPolynomial(coeffs=[12.3164,0.0202271,-9.46573e-06,1.76364e-09,-1.19242e-13,55149.9,-25.9197], Tmin=(843.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(472.432,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-O2d)CsCs) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-OdCsH) + group(Cds-OdCsH) + radical(CJC(C=O)2C) + radical(CCJ2_triplet) + radical(CC(C)(C=O)CJ=O)"""),
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
    label = '[CH]C([CH2])C(=O)[C]=O(18988)',
    structure = SMILES('[CH]C([CH2])C(=O)[C]=O'),
    E0 = (477.831,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,375,552.5,462.5,1710,3000,3100,440,815,1455,1000,1855,455,950,370.187,370.197,370.319,2888.15],'cm^-1')),
        HinderedRotor(inertia=(0.376121,'amu*angstrom^2'), symmetry=1, barrier=(36.6077,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00162552,'amu*angstrom^2'), symmetry=1, barrier=(9.62152,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0987901,'amu*angstrom^2'), symmetry=1, barrier=(9.62028,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.376155,'amu*angstrom^2'), symmetry=1, barrier=(36.6084,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.683917,0.0798047,-0.000130704,1.12033e-07,-3.73619e-11,57582.6,28.0078], Tmin=(100,'K'), Tmax=(845.801,'K')), NASAPolynomial(coeffs=[9.98974,0.0242833,-1.18224e-05,2.23825e-09,-1.52517e-13,56420.2,-12.8948], Tmin=(845.801,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(477.831,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)CsCsH) + group(Cs-CsHHH) + group(Cs-CsHHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJC(C)C=O) + radical(CCJ2_triplet) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]C(=[CH])C([O])C=O(18989)',
    structure = SMILES('[CH]C(=[CH])C([O])C=O'),
    E0 = (552.121,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,350,440,435,1725,2782.5,750,1395,475,1775,1000,3120,650,792.5,1650,522.506,522.511,522.516,522.52],'cm^-1')),
        HinderedRotor(inertia=(0.27057,'amu*angstrom^2'), symmetry=1, barrier=(52.4212,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270578,'amu*angstrom^2'), symmetry=1, barrier=(52.4211,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.270578,'amu*angstrom^2'), symmetry=1, barrier=(52.4212,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.41564,0.0564935,-4.91818e-05,2.306e-08,-4.50052e-12,66498.1,27.8888], Tmin=(100,'K'), Tmax=(1201.32,'K')), NASAPolynomial(coeffs=[10.4789,0.0263155,-1.15004e-05,2.14871e-09,-1.48759e-13,64320.6,-17.4987], Tmin=(1201.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(552.121,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(AllylJ2_triplet) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC1([O])C=O(18947)',
    structure = SMILES('[CH]=C1CC1([O])C=O'),
    E0 = (390.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.31311,0.0624708,-8.16689e-05,5.91342e-08,-1.73254e-11,47072.5,21.9359], Tmin=(100,'K'), Tmax=(831.721,'K')), NASAPolynomial(coeffs=[9.53839,0.0229115,-1.03215e-05,1.9436e-09,-1.34354e-13,45704.3,-16.2305], Tmin=(831.721,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(390.604,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsCsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-OdCsH) + group(Cds-CdsHH) + ring(Methylene_cyclopropane) + radical(C=CC(C)(C=O)OJ) + radical(Cds_P)"""),
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
    label = '[CH2]C1=CO[CH][C]1[O](18990)',
    structure = SMILES('[CH2]C1=C([O])[CH]O[CH]1'),
    E0 = (177.984,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.99175,0.0335062,9.64356e-06,-3.5412e-08,1.54694e-11,21488.2,21.8661], Tmin=(100,'K'), Tmax=(996.712,'K')), NASAPolynomial(coeffs=[11.5855,0.0200526,-7.80559e-06,1.47279e-09,-1.0643e-13,18331.6,-30.6281], Tmin=(996.712,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(177.984,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)OsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsCs) + group(Cds-CdsCsOs) + ring(25dihydrofuran) + radical(C=C(C)OJ) + radical(C=CCJ(O)C) + radical(C=CCJ(O)C) + radical(Allyl_P)"""),
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
    E0 = (283.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (1078.03,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (1117.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (523.678,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (346.997,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (929.334,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (764.084,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (693.201,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (673.149,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (490.17,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (535.851,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (446.469,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (398.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (524.544,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (495.654,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (554.729,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (585.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (818.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (742.485,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (735.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (544.658,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (513.207,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (422.701,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (422.271,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (675.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (834.397,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (794.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (291.686,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (291.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (291.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (290.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (290.71,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (308.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (308.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (456.897,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (421.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (687.692,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (557.512,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (532.59,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (572.941,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (444.061,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (734.974,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (608.01,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (802.62,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (655.264,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (661.933,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (536.133,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (694.263,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (283.597,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (628.358,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (896.161,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (672.691,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS53',
    E0 = (668.463,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS54',
    E0 = (590.732,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS55',
    E0 = (308.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS56',
    E0 = (291.881,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS57',
    E0 = (421.652,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS58',
    E0 = (595.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS59',
    E0 = (961.516,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS60',
    E0 = (607.288,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS61',
    E0 = (744.25,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS62',
    E0 = (390.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS63',
    E0 = (286.526,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS64',
    E0 = (324.688,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['OCHCO(3676)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['O(4)', '[CH]C([CH2])=[C]C=O(18041)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/OneDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction3',
    reactants = ['O(4)', '[CH]C(=C)C(=[CH])[O](16209)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]C(=C)C1=COO1(18952)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(240.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 237.9 to 240.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction5',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]C(=C)C(O)=C=O(18953)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction6',
    reactants = ['[CH]C(=[CH])C([O])[CH][O](18954)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]C([CH2])C([O])=[C][O](18955)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]C(=[CH])[C]([O])C[O](18956)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(2.05689e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 4.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH][C](C)C([O])=[C][O](18957)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(1.54267e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]C(=C)[C]1OC1[O](18958)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(5.92717e+11,'s^-1'), n=0.412677, Ea=(206.573,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra;radadd_intra] for rate rule [R3_D;doublebond_intra_secDe;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic
Ea raised from 205.2 to 206.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]C(=C)C1([O])[CH]O1(18959)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(252.254,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH][C]1COC1=C[O](18960)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(162.873,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra] for rate rule [R4_S_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 158.3 to 162.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH][C]1COC=C1[O](18961)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(8.9814e+11,'s^-1'), n=0.0209575, Ea=(114.866,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]C1([CH2])OC1=C[O](18909)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.62453e+10,'s^-1'), n=0.5217, Ea=(240.947,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra] for rate rule [R4_S_(Cd)_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 2.23606797749979
family: Intra_R_Add_Exocyclic
Ea raised from 240.3 to 240.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]C1([CH2])OC=C1[O](18962)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(212.058,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_SD_D;doublebond_intra;radadd_intra] for rate rule [R5_SD_D;doublebond_intra_2H;radadd_intra_O]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Exocyclic
Ea raised from 211.0 to 212.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction16',
    reactants = ['H(3)', '[CH]C(=C)C([O])=C=O(18963)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction17',
    reactants = ['OCHCO(3676)', '[CH]=[C][CH2](16918)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction18',
    reactants = ['[O][C]=C[O](9592)', '[CH]=[C][CH2](16918)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction19',
    reactants = ['H(3)', '[CH]C(=[CH])C([O])=C[O](18964)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(1.63133e+13,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: R_Recombination"""),
)

reaction(
    label = 'reaction20',
    reactants = ['H(3)', '[CH]C(=C)C([O])=[C][O](18965)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction21',
    reactants = ['[CH]C(=C)C([O])=[C]O(18966)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction22',
    reactants = ['[CH]C(=C)C(O)=[C][O](18967)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[CH]C(=[CH])C(O)=C[O](18968)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out_singleH;XH_out] for rate rule [R4H_DSS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction24',
    reactants = ['[CH]C(=[CH])C([O])=CO(18969)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(272000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1600,'K'), comment="""Estimated using template [R5H_DSMS;Cd_rad_out_singleH;XH_out] for rate rule [R5H_DSMS;Cd_rad_out_singleH;O_H_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction26',
    reactants = ['CH2(T)(20)', '[CH]=[C]C([O])=C[O](13159)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Y_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[C]=C([CH2])C([O])=C[O](18970)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[O]C=C([O])C1=CC1(18971)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(7.38971e+10,'s^-1'), n=0.0476667, Ea=(8.08907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Rn;C_rad_out_2H;Ypri_rad_out] for rate rule [R3_SD;C_rad_out_2H;CdsinglepriH_rad_out]
Euclidian distance = 2.8284271247461903
family: Birad_recombination"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]=C1COC1=C[O](18896)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using template [R4_SSS;C_rad_out_2H;Ypri_rad_out] for rate rule [R4_SSS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH2]C1=COC1=C[O](18820)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4;CdsingleH_rad_out;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]=C1COC=C1[O](18972)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;C_rad_out_2H;Ypri_rad_out] for rate rule [R5_SSDS;C_rad_out_2H;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction32',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH2]C1=COC=C1[O](18973)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5;CdsingleH_rad_out;Opri_rad]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]=C(C)C([O])=C=O(18974)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH2]C(=C)C(=O)[C]=O(16758)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]=C1CC([O])[C]1[O](18975)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(3.21867e+08,'s^-1'), n=0.926191, Ea=(173.3,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 170.1 to 173.3 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH2]C1=CC([O])[C]1[O](18976)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]=C1CC1([O])[CH][O](18936)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(2.50867e+09,'s^-1'), n=0.788889, Ea=(404.095,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra;radadd_intra_cs2H] for rate rule [R4_S_(Cd)_D;doublebond_intra;radadd_intra_cs2H]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 402.8 to 404.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH2]C1=CC1([O])[CH][O](18848)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(2.54e+10,'s^-1'), n=0.69, Ea=(273.915,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_D_D;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 272.9 to 273.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction39',
    reactants = ['CH2(T)(20)', 'C#CC([O])=C[O](13171)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(3.01923,'m^3/(mol*s)'), n=1.94267, Ea=(22.8894,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ct-Cd_Ct-H;YJ] for rate rule [Ct-Cd_Ct-H;CH2_triplet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[O][C]=C[O](9592)', 'C3H3(5450)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(0.00505236,'m^3/(mol*s)'), n=2.48091, Ea=(29.9954,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cdd_Cds;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=C(C)C([O])=[C][O](18977)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(111300,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_DSS;Cd_rad_out;Cs_H_out] for rate rule [R4H_DSS;Cd_rad_out;Cs_H_out_2H]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 3.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH2]C(=C)[C]([O])[C]=O(16759)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(9.15442e+12,'s^-1'), n=1.11742, Ea=(430.818,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSD;Cd_rad_out;Cd_H_out_singleH]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction43',
    reactants = ['CO(2039)', '[CH]C(=C)[CH][O](17644)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(0.0591985,'m^3/(mol*s)'), n=2.3675, Ea=(305.306,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [CO;R_H]
Euclidian distance = 0
family: 1,2_Insertion_CO"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH][C]([CH2])C([O])[C]=O(18978)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction45',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]C([CH2])=C1[CH]OO1(18979)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(371.667,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_CO;carbonyl_intra;radadd_intra_O] for rate rule [R4_S_CO;carbonyl_intra_H;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic
Ea raised from 367.0 to 371.7 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH](6993)', 'C=C=C([O])C=O(11930)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(1.70446,'m^3/(mol*s)'), n=2.07639, Ea=(14.6531,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Ca_Cds-HH;YJ] for rate rule [Ca_Cds-HH;CH_quartet]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction47',
    reactants = ['HCO(1372)', '[CH]C(=C)[C]=O(17257)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(5.2e+11,'cm^3/(mol*s)'), n=0, Ea=(93.9308,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Estimated using template [Cd_R;CO_pri_rad] for rate rule [Ck_O;CO_pri_rad]
Euclidian distance = 2.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction48',
    reactants = ['HCO(1372)', '[CH]C([CH2])=[C][O](18980)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(1.36745e+07,'m^3/(mol*s)'), n=-0.263863, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.00481396807501, var=0.0768145972539, Tref=1000.0, N=3, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R
    Total Standard Deviation in ln(k): 0.567716674236
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Sp-4R!H-2R]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]=C([CH2])C(=O)C=O(16881)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]C(=C)C([O])[C]=O(16751)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[CH][O](1548)', '[CH]C(=C)[C]=O(17257)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]C(=[CH])[C](O)[CH][O](18981)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radEndo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction53',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]C1([CH2])C(=O)C1[O](18982)'],
    transitionState = 'TS53',
    kinetics = Arrhenius(A=(1.17979e+10,'s^-1'), n=0.520585, Ea=(384.867,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_S_D;doublebond_intra_2H;radadd_intra_cs] for rate rule [R4_S_(CO)_D;doublebond_intra_2H;radadd_intra_cs]
Euclidian distance = 2.0
family: Intra_R_Add_Exocyclic
Ea raised from 383.3 to 384.9 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction54',
    reactants = ['[CH]C(=[CH])C(=O)C[O](18983)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS54',
    kinetics = Arrhenius(A=(148400,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4H_DSS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 4.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction55',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]=C1CC([O])C1=O(18984)'],
    transitionState = 'TS55',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(25.1889,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4_SSS;C_rad_out_2H;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination
Ea raised from 18.8 to 25.2 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction56',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH2]C1=CC([O])C1=O(18985)'],
    transitionState = 'TS56',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4;CdsingleH_rad_out;Ypri_rad_out]
Euclidian distance = 2.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction57',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH2]C1=CO[C]1[CH][O](18986)'],
    transitionState = 'TS57',
    kinetics = Arrhenius(A=(1.953e+11,'s^-1'), n=0.387, Ea=(138.055,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_D;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R4_D_CO;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.4142135623730951
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction58',
    reactants = ['[CH]C([CH2])([C]=O)C=O(18987)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS58',
    kinetics = Arrhenius(A=(8.889e+11,'s^-1'), n=0.232, Ea=(122.75,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;CO]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction59',
    reactants = ['OCHCO(3676)', '[CH][C][CH2](17602)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS59',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/OneDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction60',
    reactants = ['[CH]C([CH2])C(=O)[C]=O(18988)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS60',
    kinetics = Arrhenius(A=(5945.77,'s^-1'), n=2.58622, Ea=(129.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_SS;Y_rad_out;XH_out] for rate rule [R3H_SS;CO_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction61',
    reactants = ['[CH]C(=[CH])C([O])C=O(18989)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS61',
    kinetics = Arrhenius(A=(13437.7,'s^-1'), n=2.58467, Ea=(192.129,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_DS;Cd_rad_out_singleH;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction62',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH]=C1CC1([O])C=O(18947)'],
    transitionState = 'TS62',
    kinetics = Arrhenius(A=(5.94212e+13,'s^-1'), n=0.0123667, Ea=(107.007,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad_out;Cpri_rad_out_2H] + [R3_SS;Y_rad_out;Ypri_rad_out] for rate rule [R3_SS;Y_rad_out;Cpri_rad_out_2H]
Euclidian distance = 2.0
family: Birad_recombination
Ea raised from 105.8 to 107.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction63',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH2]C1=CC1([O])C=O(18883)'],
    transitionState = 'TS63',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3;Y_rad_out;Ypri_rad_out] for rate rule [R3_SD;Y_rad_out;CdsinglepriH_rad_out]
Euclidian distance = 2.23606797749979
family: Birad_recombination"""),
)

reaction(
    label = 'reaction64',
    reactants = ['[CH]C(=C)C([O])=C[O](16771)'],
    products = ['[CH2]C1=CO[CH][C]1[O](18990)'],
    transitionState = 'TS64',
    kinetics = Arrhenius(A=(1.95882e+10,'s^-1'), n=0.283532, Ea=(41.0907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

network(
    label = 'PDepNetwork #4069',
    isomers = [
        '[CH]C(=C)C([O])=C[O](16771)',
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
    label = 'PDepNetwork #4069',
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


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
    label = '[C]=[CH](9646)',
    structure = SMILES('[C]=[CH]'),
    E0 = (847.095,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([554.803,1738.79,3454.47],'cm^-1')),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (25.0293,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[3.97621,0.00212915,-8.08978e-08,-3.83305e-10,9.76908e-14,101881,6.00119], Tmin=(100,'K'), Tmax=(1982.31,'K')), NASAPolynomial(coeffs=[5.05695,0.00131032,-4.91873e-07,1.01502e-10,-7.16167e-15,101185,-0.627232], Tmin=(1982.31,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(847.095,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(133.032,'J/(mol*K)'), comment="""Thermo library: DFT_QCI_thermo + radical(Cds_P) + radical(CdCdJ2_triplet)"""),
)

species(
    label = 'C=C([O])[CH][O](2850)',
    structure = SMILES('[CH2]C([O])=C[O]'),
    E0 = (20.9566,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([350,440,435,1725,3000,3100,440,815,1455,1000,3010,987.5,1337.5,450,1655,217.215,217.577],'cm^-1')),
        HinderedRotor(inertia=(0.665078,'amu*angstrom^2'), symmetry=1, barrier=(22.287,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.00575,0.032418,-6.0508e-07,-3.82059e-08,2.20969e-11,2603.17,17.8437], Tmin=(100,'K'), Tmax=(887.044,'K')), NASAPolynomial(coeffs=[16.9262,-0.00266727,4.27963e-06,-9.58521e-10,6.70145e-14,-1310.55,-59.4906], Tmin=(887.044,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(20.9566,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(178.761,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=C(O)CJ)"""),
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
    label = '[CH]=[C]C[C]=C[O](18031)',
    structure = SMILES('[CH]=[C]C[C]=C[O]'),
    E0 = (744.199,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,3010,987.5,1337.5,450,1655,1670,1700,300,440,3120,650,792.5,1650,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.861713,'amu*angstrom^2'), symmetry=1, barrier=(19.8125,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.861266,'amu*angstrom^2'), symmetry=1, barrier=(19.8022,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22115,0.0524651,-5.33067e-05,2.66378e-08,-5.10804e-12,89613.9,25.4214], Tmin=(100,'K'), Tmax=(1353.42,'K')), NASAPolynomial(coeffs=[15.0989,0.00910259,-2.64666e-06,4.02387e-10,-2.52208e-14,86072.4,-44.9366], Tmin=(1353.42,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(744.199,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=COJ) + radical(Cds_S) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC(=[CH])[O](16204)',
    structure = SMILES('[CH]=[C]CC(=[CH])[O]'),
    E0 = (744.03,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1685,370,3115,3125,620,680,785,800,1600,1700,225.842],'cm^-1')),
        HinderedRotor(inertia=(0.280007,'amu*angstrom^2'), symmetry=1, barrier=(10.0108,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.275123,'amu*angstrom^2'), symmetry=1, barrier=(10.0352,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (80.0847,'amu'),
    collisionModel = TransportData(shapeIndex=2, epsilon=(3654.78,'J/mol'), sigma=(6.06207,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=570.87 K, Pc=37.23 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.45515,0.0540674,-6.29679e-05,3.75e-08,-8.72373e-12,89579.5,24.7964], Tmin=(100,'K'), Tmax=(1057.71,'K')), NASAPolynomial(coeffs=[12.3049,0.0130359,-4.77821e-06,8.23083e-10,-5.4697e-14,87284.3,-28.1558], Tmin=(1057.71,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(744.03,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(Cds_P) + radical(Cds_P)"""),
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
    label = '[C]=[C]CC([O])=C[O](18924)',
    structure = SMILES('[C]=[C]CC([O])=C[O]'),
    E0 = (740.61,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,1685,370,329.8,330.872,330.883,330.977],'cm^-1')),
        HinderedRotor(inertia=(0.189636,'amu*angstrom^2'), symmetry=1, barrier=(14.7579,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.188171,'amu*angstrom^2'), symmetry=1, barrier=(14.7541,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700601,0.0650435,-7.80476e-05,4.44957e-08,-9.58518e-12,89200.1,28.3435], Tmin=(100,'K'), Tmax=(1215.35,'K')), NASAPolynomial(coeffs=[17.6985,0.00588444,-1.06468e-06,9.09384e-11,-3.27445e-15,85305.8,-55.9997], Tmin=(1215.35,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(740.61,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S) + radical(CdCdJ2_triplet)"""),
)

species(
    label = '[CH]=C1CC(=C[O])O1(18919)',
    structure = SMILES('[CH]=C1CC(=C[O])O1'),
    E0 = (251.493,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.83664,0.0258717,5.45961e-05,-1.0106e-07,4.35458e-11,30345.5,23.865], Tmin=(100,'K'), Tmax=(939.036,'K')), NASAPolynomial(coeffs=[19.8178,0.00448705,5.64844e-07,-9.00358e-11,-4.45125e-15,24534.3,-74.7136], Tmin=(939.036,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(251.493,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=COJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC1=COO1(18925)',
    structure = SMILES('[CH]=[C]CC1=COO1'),
    E0 = (669.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.50473,0.0446633,-1.77352e-05,-1.28667e-08,8.81984e-12,80643.4,25.7844], Tmin=(100,'K'), Tmax=(1008.59,'K')), NASAPolynomial(coeffs=[14.7988,0.0142547,-5.69786e-06,1.11277e-09,-8.2532e-14,76826.7,-44.0924], Tmin=(1008.59,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(669.685,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-O2s(Cds-Cd)) + group(O2s-O2s(Cds-Cd)) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclobutene) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC([O])=CO1(18926)',
    structure = SMILES('[CH]=C1CC([O])=CO1'),
    E0 = (166.069,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.13292,0.0177065,7.60624e-05,-1.25258e-07,5.34461e-11,20062.5,23.0861], Tmin=(100,'K'), Tmax=(918.749,'K')), NASAPolynomial(coeffs=[19.6903,0.00179271,3.225e-06,-6.99467e-10,4.09652e-14,14281.8,-74.0327], Tmin=(918.749,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(166.069,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Cds-CdsHH) + ring(Cyclopentane) + radical(C=C(C)OJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CC(O)=C[O](18927)',
    structure = SMILES('[CH]=C=CC(O)=C[O]'),
    E0 = (107.754,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.893286,0.0800591,-9.29713e-05,4.84978e-08,-9.06337e-12,13159.7,27.563], Tmin=(100,'K'), Tmax=(1595.41,'K')), NASAPolynomial(coeffs=[23.1623,-0.0023399,5.26605e-06,-1.22963e-09,8.86262e-14,8295,-90.9187], Tmin=(1595.41,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(107.754,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=COJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=[C]CC(O)=C=O(18928)',
    structure = SMILES('[CH]=[C]CC(O)=C=O'),
    E0 = (326.599,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.467536,0.0784401,-0.000109857,7.35277e-08,-1.82879e-11,39407.4,24.4369], Tmin=(100,'K'), Tmax=(824.502,'K')), NASAPolynomial(coeffs=[16.0246,0.0116065,-3.98619e-06,6.33501e-10,-3.91226e-14,36548.4,-49.3956], Tmin=(824.502,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(326.599,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + missing(O2d-Cdd) + group(Cs-Cd(CCO)HH) + group(Cds-(Cdd-O2d)CsOs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + missing(Cdd-CdO2d) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CC([O])=CO(18929)',
    structure = SMILES('[CH]=C=CC([O])=CO'),
    E0 = (104.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[-0.772609,0.0813034,-9.86841e-05,5.38761e-08,-1.0532e-11,12712.1,27.2069], Tmin=(100,'K'), Tmax=(1521.83,'K')), NASAPolynomial(coeffs=[22.8937,-0.00282888,5.85422e-06,-1.38678e-09,1.01703e-13,8048.07,-88.5645], Tmin=(1521.83,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(104.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=C(C)OJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=CCC([O])=C=O(18930)',
    structure = SMILES('[CH]=CCC(=O)[C]=O'),
    E0 = (205.926,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.78505,0.049706,-3.9321e-05,1.45477e-08,-2.16784e-12,24845.5,23.2663], Tmin=(100,'K'), Tmax=(1536.76,'K')), NASAPolynomial(coeffs=[13.0075,0.0204957,-1.08099e-05,2.17937e-09,-1.55794e-13,21396.1,-35.6979], Tmin=(1536.76,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(205.926,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(Cds-CdsHH) + group(Cds-O2d(Cds-O2d)H) + radical(Cds_P) + radical(CCCJ=O)"""),
)

species(
    label = '[CH]=[C][CH]C([O])[CH][O](18931)',
    structure = SMILES('[CH][C]=CC([O])[CH][O]'),
    E0 = (859.708,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([1380,1390,370,380,2900,435,3010,987.5,1337.5,450,1655,3025,407.5,1350,352.5,1685,370,200,800,960,1120,1280,1440,1600],'cm^-1')),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.156089,'amu*angstrom^2'), symmetry=1, barrier=(3.5888,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 7,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.937359,0.0751437,-0.00011479,1.02787e-07,-3.66231e-11,103502,30.6472], Tmin=(100,'K'), Tmax=(828.041,'K')), NASAPolynomial(coeffs=[6.07851,0.0358543,-1.74334e-05,3.32292e-09,-2.28604e-13,103146,9.80638], Tmin=(828.041,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(859.708,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-CsOsHH) + group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + radical(CC(C)OJ) + radical(CCOJ) + radical(CCsJOH) + radical(Cds_S) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[C]1OC1[O](18932)',
    structure = SMILES('[CH]=[C]C[C]1OC1[O]'),
    E0 = (667.308,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.3477,0.06396,-9.55627e-05,7.96412e-08,-2.52702e-11,80348.9,26.5587], Tmin=(100,'K'), Tmax=(952.383,'K')), NASAPolynomial(coeffs=[6.70555,0.0253316,-9.32538e-06,1.49684e-09,-9.0566e-14,80059.7,4.81097], Tmin=(952.383,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(667.308,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-CsOsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CCOJ) + radical(C2CsJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=[C]CC1([O])[CH]O1(18933)',
    structure = SMILES('[CH]=[C]CC1([O])[CH]O1'),
    E0 = (665.275,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.111443,0.0697059,-7.88848e-05,4.33005e-08,-8.70484e-12,80167.7,27.8444], Tmin=(100,'K'), Tmax=(1459.37,'K')), NASAPolynomial(coeffs=[17.4272,0.00689385,1.45438e-06,-6.07766e-10,5.17018e-14,76748.4,-56.6388], Tmin=(1459.37,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(665.275,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsCs) + group(O2s-CsH) + group(Cs-CsCsOsOs) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(Ethylene_oxide) + radical(CC(C)(O)OJ) + radical(CCsJO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1C[C]([O])C1[O](18934)',
    structure = SMILES('[CH]=C1C[C]([O])C1[O]'),
    E0 = (655.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.2058,0.0531558,-4.50989e-05,1.87584e-08,-3.06864e-12,78907.8,25.3762], Tmin=(100,'K'), Tmax=(1473.68,'K')), NASAPolynomial(coeffs=[15.3156,0.0148578,-6.11683e-06,1.12355e-09,-7.70103e-14,74749.1,-48.1664], Tmin=(1473.68,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(655.182,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsCs) + group(Cds-CdsHH) + ring(methylenecyclobutane) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(Cds_P)"""),
)

species(
    label = '[O][C]1C[C]=CC1[O](18935)',
    structure = SMILES('[O][C]1C[C]=CC1[O]'),
    E0 = (571.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.65561,0.0408763,-8.00247e-06,-2.16473e-08,1.16366e-11,68838.6,23.5265], Tmin=(100,'K'), Tmax=(1001.23,'K')), NASAPolynomial(coeffs=[14.0982,0.015361,-6.02283e-06,1.16881e-09,-8.65702e-14,65134.3,-42.5733], Tmin=(1001.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(571.576,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-CsH) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-CdsCsH) + ring(Cyclopentene) + radical(CC(C)OJ) + radical(CC(C)OJ) + radical(C2CsJOH) + radical(cyclopentene-vinyl)"""),
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
    label = '[CH]=[C]CC([O])=C=O(18937)',
    structure = SMILES('[CH]=[C]CC(=O)[C]=O'),
    E0 = (443.768,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,375,552.5,462.5,1710,1685,370,1855,455,950,3120,650,792.5,1650,180],'cm^-1')),
        HinderedRotor(inertia=(0.0398572,'amu*angstrom^2'), symmetry=1, barrier=(16.7051,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.723915,'amu*angstrom^2'), symmetry=1, barrier=(16.6442,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.106151,'amu*angstrom^2'), symmetry=1, barrier=(16.7029,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.93287,0.051878,-5.4596e-05,3.06363e-08,-7.35063e-12,53442.1,23.2342], Tmin=(100,'K'), Tmax=(968.614,'K')), NASAPolynomial(coeffs=[8.20778,0.025965,-1.44667e-05,3.01646e-09,-2.2189e-13,52226.5,-6.83851], Tmin=(968.614,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(443.768,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(220.334,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-CdsCsH) + group(Cds-O2d(Cds-O2d)H) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC([O])=[C][O](18938)',
    structure = SMILES('[CH]=[C]CC([O])=[C][O]'),
    E0 = (669.348,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1670,1700,300,440,3120,650,792.5,1650,399.284,400.266,400.696],'cm^-1')),
        HinderedRotor(inertia=(0.0880122,'amu*angstrom^2'), symmetry=1, barrier=(9.91093,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0868656,'amu*angstrom^2'), symmetry=1, barrier=(9.8984,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (95.0761,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.06919,0.0616235,-7.62081e-05,4.65345e-08,-1.09371e-11,80612,30.1758], Tmin=(100,'K'), Tmax=(1054.07,'K')), NASAPolynomial(coeffs=[14.5233,0.0105665,-3.54937e-06,5.78967e-10,-3.73037e-14,77775.8,-35.4405], Tmin=(1054.07,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(669.348,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(224.491,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC([O])=[C]O(18939)',
    structure = SMILES('[CH]=[C]CC([O])=[C]O'),
    E0 = (527.885,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,1670,1700,300,440,3120,650,792.5,1650,205.749,209.045],'cm^-1')),
        HinderedRotor(inertia=(0.408696,'amu*angstrom^2'), symmetry=1, barrier=(12.7543,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.408411,'amu*angstrom^2'), symmetry=1, barrier=(12.7598,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.410718,'amu*angstrom^2'), symmetry=1, barrier=(12.749,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.700997,0.068037,-7.74298e-05,3.7104e-08,-4.62574e-12,63613.1,30.0903], Tmin=(100,'K'), Tmax=(881.207,'K')), NASAPolynomial(coeffs=[17.6902,0.00703318,-1.01749e-06,4.63468e-11,4.25078e-16,59993.2,-53.274], Tmin=(881.207,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(527.885,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(Cds_S) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = '[CH]=[C]CC(O)=[C][O](18940)',
    structure = SMILES('[CH]=[C]C[C](O)[C]=O'),
    E0 = (528.543,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([3615,1277.5,1000,2750,2850,1437.5,1250,1305,750,350,360,370,350,1685,370,1855,455,950,3120,650,792.5,1650,276.555],'cm^-1')),
        HinderedRotor(inertia=(0.00220124,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160277,'amu*angstrom^2'), symmetry=1, barrier=(8.69008,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160465,'amu*angstrom^2'), symmetry=1, barrier=(8.68767,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.160381,'amu*angstrom^2'), symmetry=1, barrier=(8.68974,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.04162,0.0684083,-9.69787e-05,7.17973e-08,-2.10105e-11,63672.6,30.9125], Tmin=(100,'K'), Tmax=(839.476,'K')), NASAPolynomial(coeffs=[11.4236,0.0189426,-8.59756e-06,1.61433e-09,-1.11017e-13,61929.4,-17.3587], Tmin=(839.476,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(528.543,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(241.12,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C2CsJOH) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
)

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
    collisionModel = TransportData(shapeIndex=2, epsilon=(4358.71,'J/mol'), sigma=(6.881,'angstroms'), dipoleMoment=(0,'C*m'), polarizability=(0,'angstroms^3'), rotrelaxcollnum=0, comment="""Epsilon & sigma estimated with Tc=680.82 K, Pc=30.36 bar (from Joback method)"""),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0219048,0.0669119,-6.99868e-05,3.47026e-08,-6.32141e-12,32679.5,28.1347], Tmin=(100,'K'), Tmax=(1596.64,'K')), NASAPolynomial(coeffs=[18.949,0.00549596,8.63194e-07,-3.71551e-10,3.02657e-14,28419.9,-66.446], Tmin=(1596.64,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(270.375,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)HHH) + group(Cds-Cds(Cds-Cds)O2s) + group(Cds-CdsCsH) + group(Cds-Cds(Cds-Cds)H) + group(Cds-CdsOsH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CC=CCJ) + radical(Cds_S)"""),
)

species(
    label = '[CH]=CCC([O])=[C][O](18941)',
    structure = SMILES('[CH]=CCC([O])=[C][O]'),
    E0 = (431.506,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.793804,0.0606091,-6.37106e-05,3.2836e-08,-6.45105e-12,52022.1,30.6932], Tmin=(100,'K'), Tmax=(1349.4,'K')), NASAPolynomial(coeffs=[16.734,0.00936266,-2.30353e-06,3.03835e-10,-1.73854e-14,48083.9,-49.6377], Tmin=(1349.4,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(431.506,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)(Cds-Cds)HH) + group(Cds-CdsCsOs) + group(Cds-CdsCsH) + group(Cds-CdsOsH) + group(Cds-CdsHH) + radical(C=C(C)OJ) + radical(C=COJ) + radical(C=CJO) + radical(Cds_P)"""),
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
    label = 'C#CCC([O])=C[O](16699)',
    structure = SMILES('C#CCC([O])=C[O]'),
    E0 = (114.046,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,350,440,435,1725,3010,987.5,1337.5,450,1655,2175,525,750,770,3400,2100,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(1.07761,'amu*angstrom^2'), symmetry=1, barrier=(24.7764,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(1.0801,'amu*angstrom^2'), symmetry=1, barrier=(24.8337,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 3,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.0965226,0.0651568,-6.76091e-05,3.29884e-08,-5.91402e-12,13874.9,26.843], Tmin=(100,'K'), Tmax=(1617.23,'K')), NASAPolynomial(coeffs=[19.1175,0.0047477,8.15627e-07,-3.276e-10,2.59884e-14,9470.17,-68.6632], Tmin=(1617.23,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(114.046,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-(Cds-Cd)H) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)CtHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + group(Ct-CtCs) + group(Ct-CtH) + radical(C=C(C)OJ) + radical(C=COJ)"""),
)

species(
    label = '[CH]=[C]CC([O])[C]=O(16748)',
    structure = SMILES('[CH]=[C]CC([O])[C]=O'),
    E0 = (595.649,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1380,1390,370,380,2900,435,1685,370,1855,455,950,3120,650,792.5,1650,214.043,4000],'cm^-1')),
        HinderedRotor(inertia=(0.471145,'amu*angstrom^2'), symmetry=1, barrier=(15.3178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.471194,'amu*angstrom^2'), symmetry=1, barrier=(15.3178,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.00367963,'amu*angstrom^2'), symmetry=1, barrier=(0.119627,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 5,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.44483,0.0566748,-6.35401e-05,3.74907e-08,-8.8494e-12,71731.7,30.8197], Tmin=(100,'K'), Tmax=(1029,'K')), NASAPolynomial(coeffs=[11.2492,0.018562,-7.98162e-06,1.49525e-09,-1.04075e-13,69713.9,-16.7611], Tmin=(1029,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(595.649,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(C=OCOJ) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
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
    label = '[CH]=[C]C[C]=O(9927)',
    structure = SMILES('[CH]=[C]C[C]=O'),
    E0 = (551.76,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,1855,455,950,3120,650,792.5,1650],'cm^-1')),
        HinderedRotor(inertia=(0.140741,'amu*angstrom^2'), symmetry=1, barrier=(13.9866,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.138478,'amu*angstrom^2'), symmetry=1, barrier=(13.9992,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 4,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.61483,0.0295896,-2.08603e-05,6.50482e-09,-7.91399e-13,66411.5,18.339], Tmin=(100,'K'), Tmax=(1892.32,'K')), NASAPolynomial(coeffs=[11.7879,0.0101995,-5.49011e-06,1.08988e-09,-7.60124e-14,62939.9,-31.7663], Tmin=(1892.32,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(551.76,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + radical(Cds_S) + radical(CCCJ=O) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C1CC(=O)C1[O](18942)',
    structure = SMILES('[CH]=C1CC(=O)C1[O]'),
    E0 = (314.204,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.1697,0.0229931,4.51024e-05,-7.50618e-08,2.93307e-11,37871.1,22.9924], Tmin=(100,'K'), Tmax=(1007.24,'K')), NASAPolynomial(coeffs=[14.8302,0.0153818,-7.10318e-06,1.54791e-09,-1.22642e-13,33156.4,-48.9223], Tmin=(1007.24,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(314.204,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cs-(Cds-O2d)(Cds-Cds)HH) + group(Cds-CdsCsCs) + group(Cds-OdCsCs) + group(Cds-CdsHH) + ring(Cyclobutane) + radical(C=OCOJ) + radical(Cds_P)"""),
)

species(
    label = '[CH]=C=CC(=O)C[O](18943)',
    structure = SMILES('C#CC=C([O])C[O]'),
    E0 = (230.497,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.55608,0.0562278,-6.26493e-05,3.87638e-08,-9.86278e-12,27808.3,23.6504], Tmin=(100,'K'), Tmax=(945.234,'K')), NASAPolynomial(coeffs=[9.31798,0.0233818,-1.05265e-05,2.00254e-09,-1.4014e-13,26340.9,-13.3593], Tmin=(945.234,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(230.497,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(O2s-(Cds-Cd)H) + group(Cs-(Cds-Cds)OsHH) + group(Cds-CdsCsOs) + group(Cds-CdsCtH) + group(Ct-Ct(Cds-Cds)) + group(Ct-CtH) + radical(CCOJ) + radical(C=C(C)OJ)"""),
)

species(
    label = '[CH][C]1CC(=C[O])O1(18907)',
    structure = SMILES('[CH][C]1CC(=C[O])O1'),
    E0 = (489.44,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[0.976629,0.0472253,9.08085e-06,-7.0069e-08,3.80667e-11,58993.2,22.1348], Tmin=(100,'K'), Tmax=(892.032,'K')), NASAPolynomial(coeffs=[24.4228,-0.00468038,6.85363e-06,-1.50906e-09,1.04088e-13,52692.4,-100.171], Tmin=(892.032,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(489.44,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2methyleneoxetane) + radical(C=COJ) + radical(C2CsJOC(O)) + radical(CCJ2_triplet)"""),
)

species(
    label = '[CH]=[C]C[C]1[CH]OO1(18944)',
    structure = SMILES('[CH]=[C]C[C]1[CH]OO1'),
    E0 = (887.096,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.43354,0.0656559,-0.000108038,1.04181e-07,-3.89248e-11,106777,23.0501], Tmin=(100,'K'), Tmax=(836.818,'K')), NASAPolynomial(coeffs=[3.58754,0.0352603,-1.75257e-05,3.37005e-09,-2.32458e-13,107120,17.2469], Tmin=(836.818,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(887.096,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-OsCs) + group(O2s-OsCs) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + ring(12dioxetane) + radical(C2CsJOO) + radical(CCsJOO) + radical(Cds_S) + radical(Cds_P)"""),
)

species(
    label = '[O][CH][C]1C[C]=CO1(18945)',
    structure = SMILES('[O]C=C1C[C][CH]O1'),
    E0 = (425.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.22139,0.0322271,6.45354e-05,-1.32871e-07,6.06989e-11,51249.8,22.4325], Tmin=(100,'K'), Tmax=(914.325,'K')), NASAPolynomial(coeffs=[28.0143,-0.00937078,8.72678e-06,-1.72782e-09,1.09501e-13,43189.6,-121.712], Tmin=(914.325,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(425.059,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(Cyclopentane) + radical(C=COJ) + radical(CCsJOC(O)) + radical(CCJ2_triplet)"""),
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
    label = '[CH]=[C]C[C][O](18946)',
    structure = SMILES('[CH]=[C]C[C][O]'),
    E0 = (995.282,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([2750,2850,1437.5,1250,1305,750,350,1685,370,3120,650,792.5,1650,180,180,180],'cm^-1')),
        HinderedRotor(inertia=(0.152441,'amu*angstrom^2'), symmetry=1, barrier=(3.50491,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.152395,'amu*angstrom^2'), symmetry=1, barrier=(3.50386,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 6,
    opticalIsomers = 1,
    molecularWeight = (67.066,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.02516,0.0506892,-9.13132e-05,8.72352e-08,-3.16573e-11,119769,20.3381], Tmin=(100,'K'), Tmax=(849.6,'K')), NASAPolynomial(coeffs=[4.99001,0.0208936,-1.07475e-05,2.07643e-09,-1.42739e-13,119837,9.88141], Tmin=(849.6,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(995.282,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsH) + group(Cds-CdsHH) + radical(CCOJ) + radical(Cds_S) + radical(CH2_triplet) + radical(Cds_P)"""),
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
    label = '[CH]=C=CC([O])C=O(18948)',
    structure = SMILES('[CH]=C=CC([O])C=O'),
    E0 = (268.451,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.48638,0.0522564,-5.04907e-05,2.58391e-08,-5.29695e-12,32380.3,27.5444], Tmin=(100,'K'), Tmax=(1180.08,'K')), NASAPolynomial(coeffs=[11.5981,0.0179817,-6.92391e-06,1.22664e-09,-8.27842e-14,29993.8,-22.9129], Tmin=(1180.08,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(268.451,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(249.434,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-CsH) + group(Cs-(Cds-O2d)CsOsH) + group(Cds-CdsCsH) + group(Cds-OdCsH) + group(Cds-CdsHH) + group(Cdd-CdsCds) + radical(C=OCOJ) + radical(C=C=CJ)"""),
)

species(
    label = '[CH]=CC=C([O])C=O(18949)',
    structure = SMILES('[CH]C=CC(=O)C=O'),
    E0 = (152.099,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.97756,0.0498308,-3.62267e-05,1.50282e-08,-2.93487e-12,18361.6,22.5585], Tmin=(100,'K'), Tmax=(1080.58,'K')), NASAPolynomial(coeffs=[5.65362,0.0362233,-1.73378e-05,3.37478e-09,-2.38815e-13,17567.1,4.53868], Tmin=(1080.58,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(152.099,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(245.277,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-Cds)HHH) + group(Cds-CdsCsH) + group(Cd-Cd(CO)H) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(AllylJ2_triplet)"""),
)

species(
    label = '[CH]=C1C[C]([O])[CH]O1(18950)',
    structure = SMILES('[CH][C]1CC([O])=CO1'),
    E0 = (406.73,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.35063,0.0416194,1.20058e-05,-6.61524e-08,3.54103e-11,49029.6,21.1193], Tmin=(100,'K'), Tmax=(883.435,'K')), NASAPolynomial(coeffs=[21.1682,-0.00132597,5.48797e-06,-1.28939e-09,9.14505e-14,43702.4,-82.3654], Tmin=(883.435,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(406.73,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(253.591,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsOsH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsHHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(2,3-Dihydrofuran) + radical(C=C(C)OJ) + radical(C2CsJOC(O)) + radical(CCJ2_triplet)"""),
)

species(
    label = '[O][C]1[CH]OC=[C]C1(18951)',
    structure = SMILES('[O]C1=CO[CH][C]C1'),
    E0 = (411.002,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
    molecularWeight = (96.0841,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[1.23964,0.0404592,2.58818e-05,-8.58434e-08,4.34459e-11,49550.7,19.0419], Tmin=(100,'K'), Tmax=(890.615,'K')), NASAPolynomial(coeffs=[23.4867,-0.0038628,6.89441e-06,-1.53965e-09,1.06562e-13,43383.1,-98.0887], Tmin=(890.615,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(411.002,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(257.749,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(O2s-Cs(Cds-Cd)) + group(O2s-(Cds-Cd)H) + group(Cs-CsCsHH) + group(Cs-(Cds-Cds)CsHH) + group(Cs-CsOsHH) + group(Cds-CdsCsOs) + group(Cds-CdsOsH) + ring(3,4-Dihydro-2H-pyran) + radical(C=C(C)OJ) + radical(CCsJOC(O)) + radical(CCJ2_triplet)"""),
)

species(
    label = 'C=C([O])C=O(2859)',
    structure = SMILES('[CH2]C(=O)C=O'),
    E0 = (-74.309,'kJ/mol'),
    modes = [
        HarmonicOscillator(frequencies=([375,552.5,462.5,1710,3000,3100,440,815,1455,1000,2782.5,750,1395,475,1775,1000],'cm^-1')),
        HinderedRotor(inertia=(0.43166,'amu*angstrom^2'), symmetry=1, barrier=(9.92472,'kJ/mol'), semiclassical=False),
        HinderedRotor(inertia=(0.0231641,'amu*angstrom^2'), symmetry=1, barrier=(31.0079,'kJ/mol'), semiclassical=False),
    ],
    spinMultiplicity = 2,
    opticalIsomers = 1,
    molecularWeight = (71.0547,'amu'),
    energyTransferModel = SingleExponentialDown(alpha0=(3.5886,'kJ/mol'), T0=(300,'K'), n=0.85),
    thermo = NASA(polynomials=[NASAPolynomial(coeffs=[2.689,0.0301409,-2.32871e-05,9.30138e-09,-1.56938e-12,-8891.4,16.525], Tmin=(100,'K'), Tmax=(1335.25,'K')), NASAPolynomial(coeffs=[7.34752,0.0161853,-7.60952e-06,1.4738e-09,-1.03803e-13,-10135.4,-7.29647], Tmin=(1335.25,'K'), Tmax=(5000,'K'))], Tmin=(100,'K'), Tmax=(5000,'K'), E0=(-74.309,'kJ/mol'), Cp0=(33.2579,'J/(mol*K)'), CpInf=(174.604,'J/(mol*K)'), comment="""Thermo group additivity estimation: group(Cs-(Cds-O2d)HHH) + group(Cds-O2d(Cds-O2d)Cs) + group(Cds-O2d(Cds-O2d)H) + radical(CJCC=O)"""),
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
    E0 = (429.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS2',
    E0 = (675.205,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS3',
    E0 = (923.441,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS4',
    E0 = (1263.55,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS5',
    E0 = (1263.38,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS6',
    E0 = (952.402,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS7',
    E0 = (437.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS8',
    E0 = (669.685,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS9',
    E0 = (436.717,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS10',
    E0 = (493.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS11',
    E0 = (493.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS12',
    E0 = (454.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS13',
    E0 = (454.577,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS14',
    E0 = (882.57,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS15',
    E0 = (667.437,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS16',
    E0 = (665.681,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS17',
    E0 = (655.182,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS18',
    E0 = (571.576,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS19',
    E0 = (687.839,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS20',
    E0 = (650.605,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS21',
    E0 = (469.059,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS22',
    E0 = (693.755,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS23',
    E0 = (569.959,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS24',
    E0 = (585.134,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS25',
    E0 = (818.866,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS26',
    E0 = (734.796,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS27',
    E0 = (881.14,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS28',
    E0 = (547.081,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS29',
    E0 = (690.665,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS30',
    E0 = (586.714,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS31',
    E0 = (721.922,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS32',
    E0 = (473.912,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS33',
    E0 = (823.909,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS34',
    E0 = (543.325,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS35',
    E0 = (462.644,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS36',
    E0 = (429.604,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS37',
    E0 = (759.021,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS38',
    E0 = (1038.31,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS39',
    E0 = (437.888,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS40',
    E0 = (493.004,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS41',
    E0 = (671.774,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS42',
    E0 = (887.442,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS43',
    E0 = (470.695,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS44',
    E0 = (629.26,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS45',
    E0 = (1083.27,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS46',
    E0 = (432.533,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS47',
    E0 = (452.465,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS48',
    E0 = (457.457,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS49',
    E0 = (484.311,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS50',
    E0 = (447.184,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS51',
    E0 = (772.786,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

transitionState(
    label = 'TS52',
    E0 = (704.095,'kJ/mol'),
    spinMultiplicity = 1,
    opticalIsomers = 1,
)

reaction(
    label = 'reaction1',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['OCHCO(3676)', 'C3H3(5450)'],
    transitionState = 'TS1',
    kinetics = Arrhenius(A=(5e+12,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), comment="""Exact match found for rate rule [RJJ]
Euclidian distance = 0
family: 1,4_Linear_birad_scission"""),
)

reaction(
    label = 'reaction2',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]C(=C)C([O])=C[O](16771)'],
    transitionState = 'TS2',
    kinetics = Arrhenius(A=(3.53e+06,'s^-1'), n=1.73, Ea=(245.601,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCs(-HH)CJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction3',
    reactants = ['[C]=[CH](9646)', 'C=C([O])[CH][O](2850)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS3',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [C_rad/H2/Cd;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction4',
    reactants = ['O(4)', '[CH]=[C]C[C]=C[O](18031)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS4',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_rad/NonDe;O_birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction5',
    reactants = ['O(4)', '[CH]=[C]CC(=[CH])[O](16204)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS5',
    kinetics = Arrhenius(A=(5.76208e+12,'m^3/(mol*s)'), n=-4.75325, Ea=(276.349,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""Estimated using template [Y_rad;O_birad] for rate rule [Cd_pri_rad;O_birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction6',
    reactants = ['H(3)', '[C]=[C]CC([O])=C[O](18924)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS6',
    kinetics = Arrhenius(A=(1e+07,'m^3/(mol*s)'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [H_rad;Birad]
Euclidian distance = 0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction7',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C1CC(=C[O])O1(18919)'],
    transitionState = 'TS7',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_SSS;Y_rad_out;Ypri_rad_out] for rate rule [R4_SSS;Y_rad_out;Opri_rad]
Euclidian distance = 1.0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction8',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=[C]CC1=COO1(18925)'],
    transitionState = 'TS8',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(240.082,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad_out;Ypri_rad_out] for rate rule [R4_SDS;O_rad;Opri_rad]
Euclidian distance = 1.7320508075688772
family: Birad_recombination
Ea raised from 237.9 to 240.1 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction9',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C1CC([O])=CO1(18926)'],
    transitionState = 'TS9',
    kinetics = Arrhenius(A=(7.76e+09,'s^-1'), n=0.311, Ea=(7.1128,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5;Y_rad_out;Ypri_rad_out] for rate rule [R5_SSDS;Y_rad_out;Opri_rad]
Euclidian distance = 1.4142135623730951
family: Birad_recombination"""),
)

reaction(
    label = 'reaction10',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C=CC(O)=C[O](18927)'],
    transitionState = 'TS10',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction11',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=[C]CC(O)=C=O(18928)'],
    transitionState = 'TS11',
    kinetics = Arrhenius(A=(7.437e+08,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction12',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C=CC([O])=CO(18929)'],
    transitionState = 'TS12',
    kinetics = Arrhenius(A=(1.02844e+09,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction13',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=CCC([O])=C=O(18930)'],
    transitionState = 'TS13',
    kinetics = Arrhenius(A=(5.14222e+08,'s^-1'), n=0.311, Ea=(24.9733,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4;Y_rad;XH_Rrad] for rate rule [R4radExo;Y_rad;XH_Rrad]
Euclidian distance = 1.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction14',
    reactants = ['[CH]=[C][CH]C([O])[CH][O](18931)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS14',
    kinetics = Arrhenius(A=(1.949e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction15',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=[C]C[C]1OC1[O](18932)'],
    transitionState = 'TS15',
    kinetics = Arrhenius(A=(3.473e+12,'s^-1'), n=0.247, Ea=(237.833,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_secNd;radadd_intra] for rate rule [R3_D;doublebond_intra_secNd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction16',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=[C]CC1([O])[CH]O1(18933)'],
    transitionState = 'TS16',
    kinetics = Arrhenius(A=(2.42595e+10,'s^-1'), n=0.7335, Ea=(236.077,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3_D;doublebond_intra_pri;radadd_intra] for rate rule [R3_D;doublebond_intra_pri;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction17',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C1C[C]([O])C1[O](18934)'],
    transitionState = 'TS17',
    kinetics = Arrhenius(A=(1.03419e+08,'s^-1'), n=1.06803, Ea=(225.579,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 224.3 to 225.6 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction18',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[O][C]1C[C]=CC1[O](18935)'],
    transitionState = 'TS18',
    kinetics = Arrhenius(A=(9.13977e+10,'s^-1'), n=0.0396934, Ea=(141.972,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5_linear;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic
Ea raised from 139.3 to 142.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction19',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C1CC1([O])[CH][O](18936)'],
    transitionState = 'TS19',
    kinetics = Arrhenius(A=(4.64245e+09,'s^-1'), n=0.690807, Ea=(258.235,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_D;doublebond_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic"""),
)

reaction(
    label = 'reaction20',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[O][CH]C1([O])C=[C]C1(18849)'],
    transitionState = 'TS20',
    kinetics = Arrhenius(A=(5.41e+10,'s^-1'), n=0.21, Ea=(221.001,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R5;doublebond_intra;radadd_intra_cdsingleH]
Euclidian distance = 0
family: Intra_R_Add_Exocyclic
Ea raised from 219.7 to 221.0 kJ/mol to match endothermicity of reaction."""),
)

reaction(
    label = 'reaction21',
    reactants = ['H(3)', 'C#C[CH]C([O])=C[O](18864)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS21',
    kinetics = Arrhenius(A=(197.535,'m^3/(mol*s)'), n=1.6075, Ea=(11.7084,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction22',
    reactants = ['H(3)', '[CH]=[C]CC([O])=C=O(18937)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS22',
    kinetics = Arrhenius(A=(1.98288e+06,'m^3/(mol*s)'), n=0.614653, Ea=(38.1948,'kJ/mol'), T0=(1,'K'), Tmin=(303.03,'K'), Tmax=(2000,'K'), comment="""From training reaction 188 used for Ck_O;HJ
Exact match found for rate rule [Ck_O;HJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction23',
    reactants = ['[O][C]=C[O](9592)', 'C3H3(5450)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS23',
    kinetics = Arrhenius(A=(0.0130473,'m^3/(mol*s)'), n=2.32138, Ea=(27.0133,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cds_Cdd;CJ]
Euclidian distance = 0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction24',
    reactants = ['OCHCO(3676)', '[CH]=[C][CH2](16918)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS24',
    kinetics = Arrhenius(A=(3992.79,'m^3/(mol*s)'), n=1.14733, Ea=(46.1888,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Cdd_Od;CJ] + [Ck_O;YJ] for rate rule [Ck_O;CJ]
Euclidian distance = 1.0
family: R_Addition_MultipleBond"""),
)

reaction(
    label = 'reaction25',
    reactants = ['[O][C]=C[O](9592)', '[CH]=[C][CH2](16918)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS25',
    kinetics = Arrhenius(A=(6.03e+06,'m^3/(mol*s)'), n=-5.81278e-09, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_N-1R->H_N-1CNOS->N_1COS->O_Ext-1O-R_Ext-2R-R_Ext-2R-R_N-4R!H->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction26',
    reactants = ['H(3)', '[CH][C]=CC([O])=C[O](18837)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS26',
    kinetics = Arrhenius(A=(2.80515e+06,'m^3/(mol*s)'), n=0.314888, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.0, var=33.1368631905, Tref=1000.0, N=1, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O
    Total Standard Deviation in ln(k): 11.5401827615
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_N-Sp-3R!H=2CCNNOO_2CNO->O_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction27',
    reactants = ['H(3)', '[CH]=[C]CC([O])=[C][O](18938)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS27',
    kinetics = Arrhenius(A=(8.15666e+12,'m^3/(mol*s)'), n=-1.49308, Ea=(0,'kJ/mol'), T0=(1,'K'), Tmin=(300,'K'), Tmax=(1500,'K'), uncertainty=RateUncertainty(mu=0.399348053434, var=9.35827249741, Tref=1000.0, N=6, correlation='Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O',), comment="""BM rule fitted to 2 training reactions at node Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O
    Total Standard Deviation in ln(k): 7.13613102162
Exact match found for rate rule [Root_1R->H_N-2R->S_N-2CHNO->H_N-2CNO-inRing_Ext-2CNO-R_Sp-3R!H=2CCNNOO_N-3R!H->O]
Euclidian distance = 0
family: R_Recombination"""),
)

reaction(
    label = 'reaction28',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C[CH]C([O])=C[O](18865)'],
    transitionState = 'TS28',
    kinetics = Arrhenius(A=(2.18474e+10,'s^-1'), n=0.962222, Ea=(117.477,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Cd_rad_out;Cs_H_out_H/Cd]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction29',
    reactants = ['[CH]=[C]CC([O])=[C]O(18939)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS29',
    kinetics = Arrhenius(A=(4.96519e+09,'s^-1'), n=1.05826, Ea=(162.779,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R2H_S;Cd_rad_out_Cd;XH_out] for rate rule [R2H_S;Cd_rad_out_Cd;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction30',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=[C][CH]C(O)=C[O](18866)'],
    transitionState = 'TS30',
    kinetics = Arrhenius(A=(2.39846e+08,'s^-1'), n=1.46351, Ea=(157.11,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [R3H_SS;O_rad_out;XH_out] + [R3H_SS_2Cd;Y_rad_out;XH_out] for rate rule [R3H_SS_2Cd;O_rad_out;XH_out]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction31',
    reactants = ['[CH]=[C]CC(O)=[C][O](18940)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS31',
    kinetics = Arrhenius(A=(35492.9,'s^-1'), n=2.38823, Ea=(193.378,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R3H_DS;Cd_rad_out;XH_out] for rate rule [R3H_DS;Cd_rad_out;O_H_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['C=[C][CH]C([O])=C[O](16700)'],
    transitionState = 'TS32',
    kinetics = Arrhenius(A=(74200,'s^-1'), n=2.23, Ea=(44.3086,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [RnH;Cd_rad_out_singleH;Cs_H_out_H/Cd] for rate rule [R3HJ;Cd_rad_out_singleH;Cs_H_out_H/Cd]
Euclidian distance = 2.0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction33',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=CCC([O])=[C][O](18941)'],
    transitionState = 'TS33',
    kinetics = Arrhenius(A=(9.01194e+11,'s^-1'), n=1.09397, Ea=(394.305,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4H_RSD;Cd_rad_out;XH_out] for rate rule [R4H_SSD;Cd_rad_out;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction34',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=[C][CH]C([O])=CO(18867)'],
    transitionState = 'TS34',
    kinetics = Arrhenius(A=(2.22e+06,'s^-1','*|/',3), n=1.78, Ea=(113.721,'kJ/mol'), T0=(1,'K'), Tmin=(600,'K'), Tmax=(2000,'K'), comment="""Estimated using an average for rate rule [R4H_SDS;O_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction35',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['C=[C]CC([O])=[C][O](16701)'],
    transitionState = 'TS35',
    kinetics = Arrhenius(A=(136000,'s^-1'), n=1.9199, Ea=(33.0402,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5Hall;Cd_rad_out_singleH;XH_out] for rate rule [R5HJ_1;Cd_rad_out_singleH;XH_out]
Euclidian distance = 1.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction36',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['C#CCC([O])=C[O](16699)'],
    transitionState = 'TS36',
    kinetics = Arrhenius(A=(1e+10,'s^-1'), n=0, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [NOS]
Euclidian distance = 0
family: 1,2-Birad_to_alkene"""),
)

reaction(
    label = 'reaction37',
    reactants = ['[CH]=[C]CC([O])[C]=O(16748)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS37',
    kinetics = Arrhenius(A=(3.36138e+09,'s^-1'), n=0.97875, Ea=(163.373,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [cCsCJ;CJ;C]
Euclidian distance = 0
family: 1,2_shiftC"""),
)

reaction(
    label = 'reaction38',
    reactants = ['[CH][O](1548)', '[CH]=[C]C[C]=O(9927)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS38',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_rad/NonDe;Birad]
Euclidian distance = 3.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction39',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C1CC(=O)C1[O](18942)'],
    transitionState = 'TS39',
    kinetics = Arrhenius(A=(1.62e+12,'s^-1'), n=-0.305, Ea=(8.28432,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_SSS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction40',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C=CC(=O)C[O](18943)'],
    transitionState = 'TS40',
    kinetics = Arrhenius(A=(1.4874e+09,'s^-1'), n=1.045, Ea=(63.4002,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3radExo;Y_rad;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction41',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH][C]1CC(=C[O])O1(18907)'],
    transitionState = 'TS41',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(242.17,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R4_S_CO;carbonyl_intra;radadd_intra]
Euclidian distance = 0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction42',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=[C]C[C]1[CH]OO1(18944)'],
    transitionState = 'TS42',
    kinetics = Arrhenius(A=(1.503e+11,'s^-1'), n=0.221, Ea=(457.838,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R4_linear;carbonyl_intra;radadd_intra_O] for rate rule [R4_linear;carbonyl_intra_Nd;radadd_intra_O]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction43',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[O][CH][C]1C[C]=CO1(18945)'],
    transitionState = 'TS43',
    kinetics = Arrhenius(A=(1.95882e+10,'s^-1'), n=0.283532, Ea=(41.0907,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R5_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R5_linear;carbonyl_intra;radadd_intra_cdsingleH]
Euclidian distance = 1.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction44',
    reactants = ['[CH]=[C]C=C([O])C[O](18875)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS44',
    kinetics = Arrhenius(A=(11891.5,'s^-1'), n=2.58622, Ea=(129.457,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3H_SS;Y_rad_out;XH_out]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: intra_H_migration"""),
)

reaction(
    label = 'reaction45',
    reactants = ['HCO(1372)', '[CH]=[C]C[C][O](18946)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS45',
    kinetics = Arrhenius(A=(1.34164e+08,'m^3/(mol*s)'), n=-0.791073, Ea=(55.3901,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [Y_rad;Birad] for rate rule [CO_pri_rad;Birad]
Euclidian distance = 2.0
family: Birad_R_Recombination"""),
)

reaction(
    label = 'reaction46',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C1CC1([O])C=O(18947)'],
    transitionState = 'TS46',
    kinetics = Arrhenius(A=(2.18e+16,'s^-1'), n=0, Ea=(2.9288,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R3_SS;Y_rad_out;Ypri_rad_out]
Euclidian distance = 0
family: Birad_recombination"""),
)

reaction(
    label = 'reaction47',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C=CC([O])C=O(18948)'],
    transitionState = 'TS47',
    kinetics = Arrhenius(A=(3.898e+11,'s^-1'), n=0.486, Ea=(22.8614,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2radExo;Y_rad_De;XH_Rrad]
Euclidian distance = 0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction48',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=CC=C([O])C=O(18949)'],
    transitionState = 'TS48',
    kinetics = Arrhenius(A=(2.94659e+10,'s^-1'), n=0.2847, Ea=(27.8529,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;Y_rad;XH_Rrad_De] + [R2radExo;Y_rad;XH_Rrad] for rate rule [R2radExo;Y_rad;XH_Rrad_De]
Euclidian distance = 1.0
Multiplied by reaction path degeneracy 2.0
family: Intra_Disproportionation"""),
)

reaction(
    label = 'reaction49',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[CH]=C1C[C]([O])[CH]O1(18950)'],
    transitionState = 'TS49',
    kinetics = Arrhenius(A=(1.90803e+11,'s^-1'), n=0.378641, Ea=(54.707,'kJ/mol'), T0=(1,'K'), comment="""Estimated using average of templates [Rn;carbonyl_intra_H;radadd_intra] + [R5_linear;multiplebond_intra;radadd_intra] for rate rule [R5_linear;carbonyl_intra_H;radadd_intra]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction50',
    reactants = ['[CH]=[C]CC([O])=C[O](16702)'],
    products = ['[O][C]1[CH]OC=[C]C1(18951)'],
    transitionState = 'TS50',
    kinetics = Arrhenius(A=(4.59219e+10,'s^-1'), n=0.253963, Ea=(17.5802,'kJ/mol'), T0=(1,'K'), comment="""Estimated using template [R6_linear;multiplebond_intra;radadd_intra_cdsingleH] for rate rule [R6_linear;carbonyl_intra_H;radadd_intra_cdsingleH]
Euclidian distance = 2.0
family: Intra_R_Add_Endocyclic"""),
)

reaction(
    label = 'reaction51',
    reactants = ['[C]=[CH](9646)', 'C=C([O])C=O(2859)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS51',
    kinetics = Arrhenius(A=(53.4257,'m^3/(mol*s)'), n=1.6025, Ea=(0,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [Cd_R;Y_1centerbirad]
Euclidian distance = 0
family: R_Addition_MultipleBond
Ea raised from -5.8 to 0 kJ/mol."""),
)

reaction(
    label = 'reaction52',
    reactants = ['[CH]=[C][CH]C([O])C=O(18881)'],
    products = ['[CH]=[C]CC([O])=C[O](16702)'],
    transitionState = 'TS52',
    kinetics = Arrhenius(A=(4.89315e+07,'s^-1'), n=1.5961, Ea=(158.199,'kJ/mol'), T0=(1,'K'), comment="""Estimated using an average for rate rule [R2H_S;Y_rad_out;XH_out]
Euclidian distance = 0
family: intra_H_migration"""),
)

network(
    label = 'PDepNetwork #4070',
    isomers = [
        '[CH]=[C]CC([O])=C[O](16702)',
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
    label = 'PDepNetwork #4070',
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

